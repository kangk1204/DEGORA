#!/usr/bin/env python
"""Write AURC, source-bootstrap, and DEGORA advantage metrics for a benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd
from scipy.stats import chi2

from degora.baselines import (
    BASELINE_RESULT_COLUMNS,
    _finalize,
    _study_level,
    validate_baseline_result,
    baseline_result_paths,
    rank_product_adapter,
    weighted_stouffer_adapter,
)
from degora.benchmark_stats import background_roc_pr, recall_stats
from degora.provenance import shell_command, write_source_sidecar
from degora.score_db import (
    _quality_weighted_consensus,
    _score_ready_harmonized,
    _source_unit_series,
    slice_consensus,
    study_gene_evidence,
)
from degora.stats import bh_adjust


POINT_METRIC_COLUMNS = [
    "corpus",
    "method_id",
    "setting_id",
    "run_status",
    "n_rows",
    "n_gold",
    "n_recovered_at_10",
    "recall_at_10",
    "recall_at_10_ci_low",
    "recall_at_10_ci_high",
    "precision_at_10",
    "hypergeom_p_at_10",
    "hypergeom_fdr_at_10",
    "n_recovered_at_20",
    "recall_at_20",
    "recall_at_20_ci_low",
    "recall_at_20_ci_high",
    "precision_at_20",
    "hypergeom_p_at_20",
    "hypergeom_fdr_at_20",
    "n_recovered_at_50",
    "recall_at_50",
    "recall_at_50_ci_low",
    "recall_at_50_ci_high",
    "precision_at_50",
    "hypergeom_p_at_50",
    "hypergeom_fdr_at_50",
    "n_recovered_at_100",
    "recall_at_100",
    "recall_at_100_ci_low",
    "recall_at_100_ci_high",
    "precision_at_100",
    "hypergeom_p_at_100",
    "hypergeom_fdr_at_100",
    "recall_at_max_k",
    "recall_at_max_k_ci_low",
    "recall_at_max_k_ci_high",
    "precision_at_max_k",
    "hypergeom_p_at_max_k",
    "hypergeom_fdr_at_max_k",
    "background_auroc",
    "background_auroc_se",
    "background_auroc_ci_low",
    "background_auroc_ci_high",
    "background_auprc",
    "background_prevalence",
    "background_auprc_enrichment",
    "direction_recall_at_10",
    "direction_recall_at_20",
    "direction_recall_at_50",
    "direction_recall_at_100",
    "direction_recall_at_max_k",
    "aurc_at_max_k",
    "direction_aurc_at_max_k",
    "top10",
    "failure_mode",
]

BOOTSTRAP_COLUMNS = [
    "corpus",
    "method_id",
    "setting_id",
    "metric",
    "n_bootstrap",
    "mean",
    "ci_low",
    "ci_high",
    "seed",
    "bootstrap_unit",
]

DEGORA_ADVANTAGE_COLUMNS = [
    "corpus",
    "score_variant",
    "subset",
    "n_genes",
    "median_n_source_units",
    "median_sign_concordance",
    "fraction_full_concordance",
    "median_loo_rank_stability_score",
    "median_loo_top50_fraction",
    "median_loo_top100_fraction",
    "median_source_quality_weight_sum",
]

def _vectorized_fisher_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    fisher_input = by_study.loc[:, ["gene_symbol", "study_id", "lfc", "p_from_signed_z"]].copy()
    fisher_input["p_from_signed_z"] = pd.to_numeric(fisher_input["p_from_signed_z"], errors="coerce")
    fisher_input = fisher_input.dropna(subset=["gene_symbol", "study_id", "p_from_signed_z"])
    fisher_input = fisher_input.loc[np.isfinite(fisher_input["p_from_signed_z"].to_numpy(dtype=float))].copy()
    fisher_input["p_from_signed_z"] = fisher_input["p_from_signed_z"].clip(lower=np.nextafter(0.0, 1.0), upper=1.0)
    fisher_input["_log_p"] = np.log(fisher_input["p_from_signed_z"].to_numpy(dtype=float))
    out = fisher_input.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        sum_log_p=("_log_p", "sum"),
        effect=("lfc", "mean"),
    )
    out["statistic"] = -2.0 * out["sum_log_p"]
    out["pvalue"] = chi2.sf(out["statistic"].to_numpy(dtype=float), 2 * out["n_studies"].to_numpy(dtype=float))
    out["effect"] = out["effect"].fillna(0.0)
    return _finalize(
        out,
        method_id="fisher",
        setting_id="default",
        effect="effect",
        pvalue="pvalue",
        n_studies="n_studies",
        runtime_s=time.perf_counter() - start,
        version="scipy.stats.chi2.sf",
    )


def _degora_quality_weighted_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    score_harmonized, _ = _score_ready_harmonized(harmonized)
    evidence = study_gene_evidence(score_harmonized)
    consensus = slice_consensus(score_harmonized, min_studies=min_studies)
    if consensus.empty or evidence.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    if "source_reliability_weight" not in evidence.columns:
        if "source_recommended_weight" in evidence.columns:
            evidence["source_reliability_weight"] = pd.to_numeric(evidence["source_recommended_weight"], errors="coerce")
        else:
            evidence["source_reliability_weight"] = pd.to_numeric(evidence.get("source_quality_weight", 0.65), errors="coerce")
    source_weights = evidence.loc[:, ["source_unit_id", "source_reliability_weight"]].drop_duplicates("source_unit_id")
    total_source_quality_weight = float(pd.to_numeric(source_weights["source_reliability_weight"], errors="coerce").fillna(0.65).sum())
    quality = _quality_weighted_consensus(
        evidence.loc[evidence["gene_symbol"].isin(set(consensus["gene_symbol"]))].copy(),
        total_source_quality_weight=total_source_quality_weight,
    )
    if quality.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    raw_support = score_harmonized.copy()
    raw_support["gene_symbol"] = raw_support["gene_symbol"].astype("string").str.upper().str.strip()
    raw_support["source_unit_id"] = _source_unit_series(raw_support)
    support = raw_support.dropna(subset=["gene_symbol", "source_unit_id"]).groupby("gene_symbol", as_index=False).agg(
        n_source_units=("source_unit_id", "nunique")
    )
    scores = (
        quality.merge(consensus[["gene_symbol", "weighted_lfc"]], on="gene_symbol", how="left")
        .merge(support, on="gene_symbol", how="left")
        .sort_values(
            ["quality_weighted_degora_score", "n_source_units", "quality_weighted_sign_concordance", "gene_symbol"],
            ascending=[False, False, False, True],
        )
        .reset_index(drop=True)
    )
    scores["quality_weighted_degora_rank"] = np.arange(1, len(scores) + 1, dtype=int)
    out = pd.DataFrame(
        {
            "method_id": "degora_quality_weighted_score",
            "setting_id": "quality_weighted_secondary",
            "gene_id": scores["gene_symbol"].astype("string"),
            "symbol": scores["gene_symbol"].astype("string"),
            "rank": pd.to_numeric(scores["quality_weighted_degora_rank"], errors="coerce"),
            "score": pd.to_numeric(scores["quality_weighted_degora_score"], errors="coerce"),
            "pvalue": np.nan,
            "padj": np.nan,
            "effect": pd.to_numeric(scores["weighted_lfc"], errors="coerce"),
            "direction": scores["quality_weighted_consensus_direction"].astype("string"),
            "n_studies": pd.to_numeric(scores["n_source_units"], errors="coerce"),
            "missingness": 0.0,
            "runtime_s": time.perf_counter() - start,
            "version": "degora_score_table quality_weighted_degora_rank",
            "status": "ok",
        }
    )
    return validate_baseline_result(out.sort_values(["rank", "symbol"]).reset_index(drop=True))


BOOTSTRAP_ADAPTERS: dict[str, Callable[[pd.DataFrame, int], pd.DataFrame]] = {
    "degora_quality_weighted_score": _degora_quality_weighted_adapter,
    "weighted_stouffer": weighted_stouffer_adapter,
    "fisher": _vectorized_fisher_adapter,
    "rank_product_approx": rank_product_adapter,
}

DEFAULT_BOOTSTRAP_METHODS = ("degora_quality_weighted_score", "weighted_stouffer", "fisher")


def _normalize_symbol(value: object) -> str:
    return str(value).strip().upper()


def _normalize_direction(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"up", "+", "1", "increase", "increased", "upregulated", "up-regulated"}:
        return "up"
    if text in {"down", "-", "-1", "decrease", "decreased", "downregulated", "down-regulated"}:
        return "down"
    return ""


def _read_gold(path: Path, gene_column: str) -> tuple[set[str], dict[str, str]]:
    frame = pd.read_csv(path)
    if gene_column not in frame.columns:
        raise ValueError(f"gold file {path} missing gene column {gene_column!r}")
    positives = {_normalize_symbol(value) for value in frame[gene_column].dropna()}
    positives = {gene for gene in positives if gene}
    directions: dict[str, str] = {}
    if "expected_direction" in frame.columns:
        for row in frame.to_dict(orient="records"):
            gene = _normalize_symbol(row.get(gene_column, ""))
            direction = _normalize_direction(row.get("expected_direction", ""))
            if gene in positives and direction:
                directions[gene] = direction
    if not positives:
        raise ValueError(f"gold file {path} has no positive genes")
    return positives, directions


def _ranked_records(frame: pd.DataFrame, symbol_column: str, rank_column: str, direction_column: str | None) -> list[tuple[str, str]]:
    if frame.empty:
        return []
    if symbol_column not in frame.columns or rank_column not in frame.columns:
        return []
    ranked = frame.dropna(subset=[symbol_column, rank_column]).copy()
    ranked[rank_column] = pd.to_numeric(ranked[rank_column], errors="coerce")
    ranked = ranked.dropna(subset=[rank_column])
    ranked = ranked.sort_values([rank_column, symbol_column])
    symbols = ranked[symbol_column].map(_normalize_symbol).tolist()
    if direction_column and direction_column in ranked.columns:
        directions = ranked[direction_column].map(_normalize_direction).tolist()
    else:
        directions = [""] * len(symbols)
    return list(zip(symbols, directions, strict=False))


def _degora_rankings(score_csv: Path) -> dict[tuple[str, str], list[tuple[str, str]]]:
    frame = pd.read_csv(score_csv)
    rankings = {
        ("degora_deg_score", "v1_2_source_unit_mean"): _ranked_records(
            frame,
            "gene_symbol",
            "degora_rank",
            "consensus_direction",
        )
    }
    if {"quality_weighted_degora_rank", "quality_weighted_consensus_direction"}.issubset(frame.columns):
        rankings[("degora_quality_weighted_score", "quality_weighted_secondary")] = _ranked_records(
            frame,
            "gene_symbol",
            "quality_weighted_degora_rank",
            "quality_weighted_consensus_direction",
        )
    return rankings


def _baseline_rankings(baseline_dir: Path) -> dict[tuple[str, str], list[tuple[str, str]]]:
    rankings: dict[tuple[str, str], list[tuple[str, str]]] = {}
    required = set(BASELINE_RESULT_COLUMNS)
    for path in baseline_result_paths(baseline_dir):
        frame = pd.read_csv(path, sep="\t")
        if frame.empty or not required.issubset(frame.columns):
            continue
        if set(frame["status"].dropna().astype(str)) != {"ok"}:
            continue
        key = (str(frame["method_id"].iloc[0]), str(frame["setting_id"].iloc[0]))
        rankings[key] = _ranked_records(frame, "symbol", "rank", "direction")
    return rankings


def _recall_at(ranked: list[tuple[str, str]], positives: set[str], k: int) -> float:
    return len({symbol for symbol, _ in ranked[:k]}.intersection(positives)) / len(positives)


def _direction_recall_at(ranked: list[tuple[str, str]], expected_directions: dict[str, str], k: int) -> float | str:
    if not expected_directions:
        return ""
    hits = {
        symbol
        for symbol, direction in ranked[:k]
        if expected_directions.get(symbol) and expected_directions[symbol] == direction
    }
    return len(hits) / len(expected_directions)


def _aurc(ranked: list[tuple[str, str]], positives: set[str], max_k: int) -> float:
    seen: set[str] = set()
    values: list[float] = []
    for i in range(max_k):
        if i < len(ranked) and ranked[i][0] in positives:
            seen.add(ranked[i][0])
        values.append(len(seen) / len(positives))
    return float(np.mean(values))


def _direction_aurc(ranked: list[tuple[str, str]], expected_directions: dict[str, str], max_k: int) -> float | str:
    if not expected_directions:
        return ""
    seen: set[str] = set()
    values: list[float] = []
    for i in range(max_k):
        if i < len(ranked):
            symbol, direction = ranked[i]
            if expected_directions.get(symbol) == direction:
                seen.add(symbol)
        values.append(len(seen) / len(expected_directions))
    return float(np.mean(values))


def point_metrics(
    *,
    corpus: str,
    baseline_dir: Path,
    gold_path: Path,
    degora_score_csv: Path,
    gold_gene_column: str = "gene_symbol",
    max_k: int = 1000,
) -> pd.DataFrame:
    positives, expected_directions = _read_gold(gold_path, gold_gene_column)
    rankings = {**_degora_rankings(degora_score_csv), **_baseline_rankings(baseline_dir)}
    failures = _failure_modes(baseline_dir)
    rows: list[dict[str, Any]] = []
    for key, ranked in sorted(rankings.items()):
        method_id, setting_id = key
        ranked_symbols = [symbol for symbol, _ in ranked]
        recall_detail = {k: recall_stats(ranked_symbols, positives, k, universe_size=len(ranked)) for k in (10, 20, 50, 100)}
        max_detail = recall_stats(ranked_symbols, positives, max_k, universe_size=len(ranked))
        background = background_roc_pr(ranked_symbols, positives)
        rows.append(
            {
                "corpus": corpus,
                "method_id": method_id,
                "setting_id": setting_id,
                "run_status": "ok",
                "n_rows": len(ranked),
                "n_gold": len(positives),
                "n_recovered_at_10": recall_detail[10].n_recovered,
                "recall_at_10": recall_detail[10].recall,
                "recall_at_10_ci_low": recall_detail[10].ci_low,
                "recall_at_10_ci_high": recall_detail[10].ci_high,
                "precision_at_10": recall_detail[10].precision,
                "hypergeom_p_at_10": recall_detail[10].hypergeom_pvalue,
                "hypergeom_fdr_at_10": "",
                "n_recovered_at_20": recall_detail[20].n_recovered,
                "recall_at_20": recall_detail[20].recall,
                "recall_at_20_ci_low": recall_detail[20].ci_low,
                "recall_at_20_ci_high": recall_detail[20].ci_high,
                "precision_at_20": recall_detail[20].precision,
                "hypergeom_p_at_20": recall_detail[20].hypergeom_pvalue,
                "hypergeom_fdr_at_20": "",
                "n_recovered_at_50": recall_detail[50].n_recovered,
                "recall_at_50": recall_detail[50].recall,
                "recall_at_50_ci_low": recall_detail[50].ci_low,
                "recall_at_50_ci_high": recall_detail[50].ci_high,
                "precision_at_50": recall_detail[50].precision,
                "hypergeom_p_at_50": recall_detail[50].hypergeom_pvalue,
                "hypergeom_fdr_at_50": "",
                "n_recovered_at_100": recall_detail[100].n_recovered,
                "recall_at_100": recall_detail[100].recall,
                "recall_at_100_ci_low": recall_detail[100].ci_low,
                "recall_at_100_ci_high": recall_detail[100].ci_high,
                "precision_at_100": recall_detail[100].precision,
                "hypergeom_p_at_100": recall_detail[100].hypergeom_pvalue,
                "hypergeom_fdr_at_100": "",
                "recall_at_max_k": max_detail.recall,
                "recall_at_max_k_ci_low": max_detail.ci_low,
                "recall_at_max_k_ci_high": max_detail.ci_high,
                "precision_at_max_k": max_detail.precision,
                "hypergeom_p_at_max_k": max_detail.hypergeom_pvalue,
                "hypergeom_fdr_at_max_k": "",
                **background,
                "direction_recall_at_10": _direction_recall_at(ranked, expected_directions, 10),
                "direction_recall_at_20": _direction_recall_at(ranked, expected_directions, 20),
                "direction_recall_at_50": _direction_recall_at(ranked, expected_directions, 50),
                "direction_recall_at_100": _direction_recall_at(ranked, expected_directions, 100),
                "direction_recall_at_max_k": _direction_recall_at(ranked, expected_directions, max_k),
                "aurc_at_max_k": _aurc(ranked, positives, max_k),
                "direction_aurc_at_max_k": _direction_aurc(ranked, expected_directions, max_k),
                "top10": ";".join(symbol for symbol, _ in ranked[:10]),
                "failure_mode": "",
            }
        )
    for method_id, failure_mode in failures.items():
        if any(row["method_id"] == method_id for row in rows):
            continue
        rows.append(
            {
                "corpus": corpus,
                "method_id": method_id,
                "setting_id": "default",
                "run_status": "blocked",
                "n_rows": 0,
                "n_gold": len(positives),
                "n_recovered_at_10": "",
                "recall_at_10": "",
                "recall_at_10_ci_low": "",
                "recall_at_10_ci_high": "",
                "precision_at_10": "",
                "hypergeom_p_at_10": "",
                "hypergeom_fdr_at_10": "",
                "n_recovered_at_20": "",
                "recall_at_20": "",
                "recall_at_20_ci_low": "",
                "recall_at_20_ci_high": "",
                "precision_at_20": "",
                "hypergeom_p_at_20": "",
                "hypergeom_fdr_at_20": "",
                "n_recovered_at_50": "",
                "recall_at_50": "",
                "recall_at_50_ci_low": "",
                "recall_at_50_ci_high": "",
                "precision_at_50": "",
                "hypergeom_p_at_50": "",
                "hypergeom_fdr_at_50": "",
                "n_recovered_at_100": "",
                "recall_at_100": "",
                "recall_at_100_ci_low": "",
                "recall_at_100_ci_high": "",
                "precision_at_100": "",
                "hypergeom_p_at_100": "",
                "hypergeom_fdr_at_100": "",
                "recall_at_max_k": "",
                "recall_at_max_k_ci_low": "",
                "recall_at_max_k_ci_high": "",
                "precision_at_max_k": "",
                "hypergeom_p_at_max_k": "",
                "hypergeom_fdr_at_max_k": "",
                "background_auroc": "",
                "background_auroc_se": "",
                "background_auroc_ci_low": "",
                "background_auroc_ci_high": "",
                "background_auprc": "",
                "background_prevalence": "",
                "background_auprc_enrichment": "",
                "direction_recall_at_10": "",
                "direction_recall_at_20": "",
                "direction_recall_at_50": "",
                "direction_recall_at_100": "",
                "direction_recall_at_max_k": "",
                "aurc_at_max_k": "",
                "direction_aurc_at_max_k": "",
                "top10": "",
                "failure_mode": failure_mode,
            }
        )
    out = pd.DataFrame.from_records(rows, columns=POINT_METRIC_COLUMNS)
    for label in ["10", "20", "50", "100", "max_k"]:
        p_column = f"hypergeom_p_at_{label}"
        fdr_column = f"hypergeom_fdr_at_{label}"
        out[fdr_column] = np.nan
        ok = out["run_status"].eq("ok") & pd.to_numeric(out[p_column], errors="coerce").notna()
        out.loc[ok, fdr_column] = bh_adjust(pd.to_numeric(out.loc[ok, p_column], errors="coerce").to_numpy(dtype=float))
    return out


def _failure_modes(baseline_dir: Path) -> dict[str, str]:
    path = baseline_dir / "baseline_failure_ledger.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    if frame.empty:
        return {}
    return {
        str(row.method_id): str(row.blocker_id)
        for row in frame.itertuples(index=False)
    }


def _bootstrap_sample(harmonized: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    frame = harmonized.copy()
    frame["_source_unit"] = _source_unit_series(frame)
    source_units = sorted(frame["_source_unit"].dropna().astype(str).unique())
    if len(source_units) < 2:
        raise ValueError("source-unit bootstrap requires at least two source units")
    sampled = rng.choice(source_units, size=len(source_units), replace=True)
    parts: list[pd.DataFrame] = []
    for i, source_unit in enumerate(sampled):
        part = frame.loc[frame["_source_unit"].astype(str).eq(str(source_unit))].copy()
        suffix = f"__boot{i:03d}"
        for column in ["study_id", "source_unit_id", "paper_id"]:
            if column in part.columns:
                part[column] = part[column].astype(str) + suffix
        parts.append(part.drop(columns=["_source_unit"]))
    return pd.concat(parts, ignore_index=True)


def bootstrap_metrics(
    *,
    corpus: str,
    harmonized_path: Path,
    gold_path: Path,
    gold_gene_column: str = "gene_symbol",
    min_studies: int = 2,
    n_bootstrap: int = 50,
    max_k: int = 1000,
    seed: int = 1729,
    bootstrap_methods: Iterable[str] = DEFAULT_BOOTSTRAP_METHODS,
) -> pd.DataFrame:
    positives, _ = _read_gold(gold_path, gold_gene_column)
    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    rng = np.random.default_rng(seed)
    selected_adapters: dict[str, Callable[[pd.DataFrame, int], pd.DataFrame]] = {}
    for method_id in bootstrap_methods:
        if method_id not in BOOTSTRAP_ADAPTERS:
            valid = ", ".join(sorted(BOOTSTRAP_ADAPTERS))
            raise ValueError(f"unknown bootstrap method {method_id!r}; valid choices: {valid}")
        selected_adapters[method_id] = BOOTSTRAP_ADAPTERS[method_id]
    values: dict[tuple[str, str], list[float]] = {}
    for _ in range(n_bootstrap):
        sampled = _bootstrap_sample(harmonized, rng)
        for method_id, adapter in selected_adapters.items():
            try:
                result = adapter(sampled, min_studies=min_studies)
            except Exception:
                continue
            ranked = _ranked_records(result, "symbol", "rank", "direction")
            values.setdefault((method_id, "recall_at_100"), []).append(_recall_at(ranked, positives, 100))
            values.setdefault((method_id, "aurc_at_max_k"), []).append(_aurc(ranked, positives, max_k))
    rows: list[dict[str, Any]] = []
    for (method_id, metric), metric_values in sorted(values.items()):
        numeric = np.asarray(metric_values, dtype=float)
        if numeric.size == 0:
            continue
        rows.append(
            {
                "corpus": corpus,
                "method_id": method_id,
                "setting_id": "default",
                "metric": metric,
                "n_bootstrap": int(numeric.size),
                "mean": float(np.mean(numeric)),
                "ci_low": float(np.quantile(numeric, 0.025)),
                "ci_high": float(np.quantile(numeric, 0.975)),
                "seed": seed,
                "bootstrap_unit": "source_unit_id",
            }
        )
    return pd.DataFrame.from_records(rows, columns=BOOTSTRAP_COLUMNS)


def degora_advantage_metrics(
    *,
    corpus: str,
    degora_score_csv: Path,
    gold_path: Path,
    gold_gene_column: str = "gene_symbol",
    top_n: int = 100,
) -> pd.DataFrame:
    positives, _ = _read_gold(gold_path, gold_gene_column)
    score = pd.read_csv(degora_score_csv)
    score["_symbol"] = score["gene_symbol"].map(_normalize_symbol)
    score["degora_rank"] = pd.to_numeric(score["degora_rank"], errors="coerce")
    score["quality_weighted_degora_rank"] = pd.to_numeric(
        score.get("quality_weighted_degora_rank", score["degora_rank"]),
        errors="coerce",
    )

    subsets = {
        "all_scored": score,
        f"top{top_n}_primary": score.nsmallest(top_n, "degora_rank"),
        f"top{top_n}_quality_weighted": score.nsmallest(top_n, "quality_weighted_degora_rank"),
        "locked_gold": score.loc[score["_symbol"].isin(positives)],
    }
    rows: list[dict[str, Any]] = []
    for subset_name, subset in subsets.items():
        rows.append(_degora_advantage_row(corpus, subset_name, subset))
    return pd.DataFrame.from_records(rows, columns=DEGORA_ADVANTAGE_COLUMNS)


def _degora_advantage_row(corpus: str, subset_name: str, subset: pd.DataFrame) -> dict[str, Any]:
    def median(column: str) -> float | str:
        if column not in subset.columns or subset.empty:
            return ""
        values = pd.to_numeric(subset[column], errors="coerce").dropna()
        return "" if values.empty else float(values.median())

    sign = pd.to_numeric(subset.get("sign_concordance", pd.Series(dtype=float)), errors="coerce")
    full_concordance = "" if sign.dropna().empty else float(sign.ge(1.0).mean())
    return {
        "corpus": corpus,
        "score_variant": "degora_quality_weighted_secondary",
        "subset": subset_name,
        "n_genes": int(len(subset)),
        "median_n_source_units": median("n_source_units"),
        "median_sign_concordance": median("sign_concordance"),
        "fraction_full_concordance": full_concordance,
        "median_loo_rank_stability_score": median("loo_rank_stability_score"),
        "median_loo_top50_fraction": median("loo_top50_fraction"),
        "median_loo_top100_fraction": median("loo_top100_fraction"),
        "median_source_quality_weight_sum": median("source_quality_weight_sum"),
    }


def write_markdown(point: pd.DataFrame, bootstrap: pd.DataFrame, advantage: pd.DataFrame, output: Path, *, title: str) -> None:
    lines = [
        f"# {title}",
        "",
        "## Point Metrics",
        "",
        "| method | status | recall@10 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | AUROC | AUPRC enrichment | AURC | top10 |",
        "| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in point.itertuples(index=False):
        def fmt(value: object) -> str:
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            return "" if pd.isna(numeric) else f"{float(numeric):.3f}"

        lines.append(
            f"| {row.method_id}/{row.setting_id} | {row.run_status} | "
            f"{fmt(row.recall_at_10)} ({fmt(row.recall_at_10_ci_low)}-{fmt(row.recall_at_10_ci_high)}) | "
            f"{fmt(row.recall_at_50)} ({fmt(row.recall_at_50_ci_low)}-{fmt(row.recall_at_50_ci_high)}) | "
            f"{fmt(row.recall_at_100)} ({fmt(row.recall_at_100_ci_low)}-{fmt(row.recall_at_100_ci_high)}) | "
            f"{fmt(row.precision_at_100)} | {fmt(row.background_auroc)} | {fmt(row.background_auprc_enrichment)} | "
            f"{fmt(row.aurc_at_max_k)} | {row.top10} |"
        )
    lines.extend(["", "## Source-Unit Bootstrap", ""])
    if bootstrap.empty:
        lines.append("No bootstrap rows emitted.")
    else:
        lines.extend(["| method | metric | mean | 95% CI | n |", "| --- | --- | ---: | --- | ---: |"])
        for row in bootstrap.itertuples(index=False):
            lines.append(
                f"| {row.method_id} | {row.metric} | {row.mean:.3f} | {row.ci_low:.3f}-{row.ci_high:.3f} | {row.n_bootstrap} |"
            )
    lines.extend(["", "## DEGORA Advantage Metrics", ""])
    lines.extend(["| subset | n | median source units | median concordance | median LOO stability |", "| --- | ---: | ---: | ---: | ---: |"])
    for row in advantage.itertuples(index=False):
        lines.append(
            f"| {row.subset} | {row.n_genes} | {row.median_n_source_units} | "
            f"{row.median_sign_concordance} | {row.median_loo_rank_stability_score} |"
        )
    output.write_text("\n".join(lines) + "\n")


def write_outputs(
    *,
    corpus: str,
    baseline_dir: Path,
    harmonized_path: Path,
    gold_path: Path,
    degora_score_csv: Path,
    output_dir: Path,
    title: str,
    gold_gene_column: str = "gene_symbol",
    min_studies: int = 2,
    max_k: int = 1000,
    n_bootstrap: int = 50,
    seed: int = 1729,
    bootstrap_methods: Iterable[str] = DEFAULT_BOOTSTRAP_METHODS,
    command: str,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    point = point_metrics(
        corpus=corpus,
        baseline_dir=baseline_dir,
        gold_path=gold_path,
        degora_score_csv=degora_score_csv,
        gold_gene_column=gold_gene_column,
        max_k=max_k,
    )
    bootstrap = bootstrap_metrics(
        corpus=corpus,
        harmonized_path=harmonized_path,
        gold_path=gold_path,
        gold_gene_column=gold_gene_column,
        min_studies=min_studies,
        n_bootstrap=n_bootstrap,
        max_k=max_k,
        seed=seed,
        bootstrap_methods=bootstrap_methods,
    )
    advantage = degora_advantage_metrics(
        corpus=corpus,
        degora_score_csv=degora_score_csv,
        gold_path=gold_path,
        gold_gene_column=gold_gene_column,
        top_n=100,
    )
    point_path = output_dir / f"{corpus}_deep_point_metrics.csv"
    bootstrap_path = output_dir / f"{corpus}_source_bootstrap_metrics.csv"
    advantage_path = output_dir / f"{corpus}_degora_advantage_metrics.csv"
    md_path = output_dir / f"{corpus}_deep_metrics_report.md"
    point.to_csv(point_path, index=False, quoting=csv.QUOTE_MINIMAL)
    bootstrap.to_csv(bootstrap_path, index=False, quoting=csv.QUOTE_MINIMAL)
    advantage.to_csv(advantage_path, index=False, quoting=csv.QUOTE_MINIMAL)
    write_markdown(point, bootstrap, advantage, md_path, title=title)
    inputs = [baseline_dir, harmonized_path, gold_path, degora_score_csv]
    for artifact in [point_path, bootstrap_path, advantage_path, md_path]:
        write_source_sidecar(artifact, command, inputs=inputs, metadata={"generator": "benchmark-deep-metrics", "corpus": corpus})
    return {
        "corpus": corpus,
        "point_metrics": str(point_path),
        "bootstrap_metrics": str(bootstrap_path),
        "degora_advantage_metrics": str(advantage_path),
        "report": str(md_path),
        "n_point_rows": int(len(point)),
        "n_bootstrap_rows": int(len(bootstrap)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--gold-gene-column", default="gene_symbol")
    parser.add_argument("--degora-score-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--title", default="Benchmark Deep Metrics")
    parser.add_argument("--min-studies", type=int, default=2)
    parser.add_argument("--max-k", type=int, default=1000)
    parser.add_argument("--n-bootstrap", type=int, default=50)
    parser.add_argument(
        "--bootstrap-methods",
        nargs="+",
        default=list(DEFAULT_BOOTSTRAP_METHODS),
        choices=sorted(BOOTSTRAP_ADAPTERS),
        help=(
            "Methods to rerun inside source-unit bootstrap. "
            "rank_product_approx is available but not default because it is much slower."
        ),
    )
    parser.add_argument("--seed", type=int, default=1729)
    args = parser.parse_args(argv)
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_benchmark_deep_metrics.py",
            "--corpus",
            args.corpus,
            "--baseline-dir",
            args.baseline_dir,
            "--harmonized",
            args.harmonized,
            "--gold",
            args.gold,
            "--gold-gene-column",
            args.gold_gene_column,
            "--degora-score-csv",
            args.degora_score_csv,
            "--output-dir",
            args.output_dir,
            "--title",
            args.title,
            "--min-studies",
            args.min_studies,
            "--max-k",
            args.max_k,
            "--n-bootstrap",
            args.n_bootstrap,
            "--bootstrap-methods",
            *args.bootstrap_methods,
            "--seed",
            args.seed,
        ]
    )
    summary = write_outputs(
        corpus=args.corpus,
        baseline_dir=args.baseline_dir,
        harmonized_path=args.harmonized,
        gold_path=args.gold,
        degora_score_csv=args.degora_score_csv,
        output_dir=args.output_dir,
        title=args.title,
        gold_gene_column=args.gold_gene_column,
        min_studies=args.min_studies,
        max_k=args.max_k,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        bootstrap_methods=args.bootstrap_methods,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
