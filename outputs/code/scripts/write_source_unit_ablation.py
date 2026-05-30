#!/usr/bin/env python
"""Compare old contrast-level consensus with source-unit-aware DEGORA consensus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

from degora.aggregate import SOURCE_UNIT_COLLAPSE_RULE, slice_consensus
from degora.gold import HIF1A_UP_TARGETS
from degora.metrics import recall_at_k
from degora.provenance import write_source_sidecar
from degora.slice_runner import catalog_include_mask, read_catalog
from degora.stats import bh_adjust


def _weight(frame: pd.DataFrame) -> pd.Series:
    n_ctrl = pd.to_numeric(frame.get("n_ctrl"), errors="coerce")
    n_treat = pd.to_numeric(frame.get("n_treat"), errors="coerce")
    total = n_ctrl + n_treat
    return pd.Series(np.where(np.isfinite(total) & total.gt(0), np.sqrt(total), 1.0), index=frame.index)


def contrast_level_consensus(harmonized: pd.DataFrame, min_studies: int) -> pd.DataFrame:
    """Return the pre-v1.1 contrast-level consensus for ablation only."""

    frame = harmonized.dropna(subset=["gene_symbol", "study_id", "signed_z", "normalized_rank"]).copy()
    frame["signed_z"] = pd.to_numeric(frame["signed_z"], errors="coerce")
    frame["lfc"] = pd.to_numeric(frame["lfc"], errors="coerce")
    frame["normalized_rank"] = pd.to_numeric(frame["normalized_rank"], errors="coerce")
    frame["weight"] = _weight(frame)
    by_study = frame.groupby(["gene_symbol", "study_id"], as_index=False).agg(
        signed_z=("signed_z", "mean"),
        lfc=("lfc", "mean"),
        normalized_rank=("normalized_rank", "mean"),
        weight=("weight", "first"),
    )
    n_studies = by_study.groupby("gene_symbol")["study_id"].nunique().rename("n_studies")
    by_study = by_study.merge(n_studies, on="gene_symbol", how="left")
    by_study = by_study.loc[by_study["n_studies"].ge(min_studies)].copy()
    if by_study.empty:
        return pd.DataFrame(
            columns=[
                "gene_symbol",
                "n_studies",
                "stouffer_z",
                "stouffer_p",
                "weighted_lfc",
                "sign_concordance",
                "stouffer_padj",
                "n_studies_rank",
                "rank_product",
                "rank_score",
                "slice_rank",
            ]
        )

    by_study["_wz"] = by_study["weight"] * by_study["signed_z"]
    by_study["_w2"] = by_study["weight"] ** 2
    by_study["_wlfc"] = by_study["weight"] * by_study["lfc"]
    stouffer = by_study.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        sum_wz=("_wz", "sum"),
        sum_w2=("_w2", "sum"),
        sum_wlfc=("_wlfc", "sum"),
        sum_w=("weight", "sum"),
    )
    stouffer["stouffer_z"] = stouffer["sum_wz"] / np.sqrt(stouffer["sum_w2"])
    stouffer["stouffer_p"] = 2.0 * norm.sf(stouffer["stouffer_z"].abs())
    stouffer["weighted_lfc"] = stouffer["sum_wlfc"] / stouffer["sum_w"]

    signs = by_study[["gene_symbol", "signed_z"]].merge(stouffer[["gene_symbol", "stouffer_z"]], on="gene_symbol")
    signs["_combined_sign"] = np.sign(signs["stouffer_z"])
    signs["_concordant"] = np.where(
        signs["_combined_sign"].ne(0),
        np.sign(signs["signed_z"]).eq(signs["_combined_sign"]),
        False,
    )
    sign_concordance = signs.groupby("gene_symbol")["_concordant"].mean().rename("sign_concordance")
    stouffer = stouffer.drop(columns=["sum_wz", "sum_w2", "sum_wlfc", "sum_w"]).merge(
        sign_concordance,
        on="gene_symbol",
        how="left",
    )
    stouffer["stouffer_padj"] = bh_adjust(stouffer["stouffer_p"].to_numpy(dtype=float))

    eps = np.finfo(float).tiny
    ranks = by_study.copy()
    ranks["_log_rank"] = np.log(ranks["normalized_rank"].clip(lower=eps, upper=1.0))
    rank_product = ranks.groupby("gene_symbol", as_index=False).agg(
        n_studies_rank=("study_id", "nunique"),
        mean_log_rank=("_log_rank", "mean"),
    )
    rank_product["rank_product"] = np.exp(rank_product["mean_log_rank"])
    rank_product["rank_score"] = -np.log(rank_product["rank_product"])
    out = (
        stouffer.merge(rank_product.drop(columns=["mean_log_rank"]), on="gene_symbol", how="left")
        .sort_values(["stouffer_padj", "stouffer_p", "gene_symbol"])
        .reset_index(drop=True)
    )
    out["slice_rank"] = np.arange(1, len(out) + 1)
    return out


def source_unit_best_signal_consensus(harmonized: pd.DataFrame, min_studies: int) -> pd.DataFrame:
    """Diagnostic-only source-unit best-signal variant for transient-response sensitivity."""

    frame = harmonized.dropna(subset=["gene_symbol", "study_id", "signed_z", "normalized_rank"]).copy()
    if frame.empty:
        return slice_consensus(frame, min_studies=min_studies)
    source_unit = frame.get("paper_id", frame["study_id"]).astype("string").fillna("").str.strip()
    study_id = frame["study_id"].astype("string").fillna("").str.strip()
    frame["_source_unit_id"] = source_unit.mask(source_unit.eq(""), study_id)
    frame["_abs_z"] = pd.to_numeric(frame["signed_z"], errors="coerce").abs()
    frame["_normalized_rank_numeric"] = pd.to_numeric(frame["normalized_rank"], errors="coerce")
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.upper().str.strip()
    best = (
        frame.sort_values(
            ["gene_symbol", "_source_unit_id", "_abs_z", "_normalized_rank_numeric", "study_id"],
            ascending=[True, True, False, True, True],
            kind="mergesort",
        )
        .drop_duplicates(["gene_symbol", "_source_unit_id"], keep="first")
        .drop(columns=["_source_unit_id", "_abs_z", "_normalized_rank_numeric"])
    )
    return slice_consensus(best, min_studies=min_studies)


def _indirect_study_ids(catalog_path: Path | None, explicit_exclude: list[str]) -> list[str]:
    if catalog_path is None or not catalog_path.exists():
        return sorted(set(explicit_exclude))
    catalog = read_catalog(catalog_path)
    active = catalog.loc[catalog_include_mask(catalog)].copy()
    text = (
        active["notes"].fillna("").astype(str)
        + " "
        + active["hypoxia_modality"].fillna("").astype(str)
    ).str.lower()
    pattern = r"provisional sensitivity|sensitivity evidence only|indirect"
    inferred = active.loc[text.str.contains(pattern, regex=True), "study_id"].astype(str).tolist()
    return sorted(set(explicit_exclude).union(inferred))


def _variant_summary(name: str, frame: pd.DataFrame, harmonized: pd.DataFrame, k_values: list[int]) -> dict[str, Any]:
    recalls = {str(k): recall_at_k(frame, HIF1A_UP_TARGETS, k) for k in k_values}
    return {
        "variant": name,
        "n_input_rows": int(len(harmonized)),
        "n_input_studies": int(harmonized["study_id"].nunique()) if "study_id" in harmonized.columns else 0,
        "n_input_source_units": int(harmonized.get("paper_id", harmonized["study_id"]).nunique())
        if "study_id" in harmonized.columns
        else 0,
        "n_scored_genes": int(len(frame)),
        "top20": frame.head(20)["gene_symbol"].astype(str).tolist(),
        "recall": recalls,
        "top50_recall": recalls["50"]["recall"] if "50" in recalls else None,
        "top100_recall": recalls["100"]["recall"] if "100" in recalls else None,
    }


def _jaccard(a: pd.DataFrame, b: pd.DataFrame, k: int) -> float:
    left = set(a.head(k)["gene_symbol"].astype(str).str.upper())
    right = set(b.head(k)["gene_symbol"].astype(str).str.upper())
    union = left.union(right)
    return round(len(left.intersection(right)) / len(union), 6) if union else 0.0


def _rank_table(variants: dict[str, pd.DataFrame]) -> pd.DataFrame:
    genes = sorted(set().union(*(set(frame["gene_symbol"].astype(str)) for frame in variants.values())))
    out = pd.DataFrame({"gene_symbol": genes})
    hif = {gene.upper() for gene in HIF1A_UP_TARGETS}
    out["is_hif_target"] = out["gene_symbol"].str.upper().isin(hif)
    for name, frame in variants.items():
        subset = frame[["gene_symbol", "slice_rank", "weighted_lfc", "sign_concordance", "n_studies"]].copy()
        subset = subset.rename(
            columns={
                "slice_rank": f"{name}_rank",
                "weighted_lfc": f"{name}_weighted_lfc",
                "sign_concordance": f"{name}_sign_concordance",
                "n_studies": f"{name}_n_units",
            }
        )
        out = out.merge(subset, on="gene_symbol", how="left")
    return out.sort_values(["is_hif_target", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)


def write_ablation(
    harmonized_path: Path,
    output_dir: Path,
    *,
    catalog_path: Path | None,
    min_studies: int,
    explicit_exclude: list[str],
    command: str,
) -> dict[str, Any]:
    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    excluded = _indirect_study_ids(catalog_path, explicit_exclude)
    filtered = harmonized.loc[~harmonized["study_id"].astype(str).isin(excluded)].copy()
    variants = {
        "contrast_level_legacy": contrast_level_consensus(harmonized, min_studies=min_studies),
        "source_unit_v1_2_mean": slice_consensus(harmonized, min_studies=min_studies),
        "source_unit_best_signal_sensitivity": source_unit_best_signal_consensus(harmonized, min_studies=min_studies),
        "source_unit_v1_2_mean_excluding_sensitivity": slice_consensus(filtered, min_studies=min_studies),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    k_values = [50, 100]
    summaries = [
        _variant_summary("contrast_level_legacy", variants["contrast_level_legacy"], harmonized, k_values),
        _variant_summary("source_unit_v1_2_mean", variants["source_unit_v1_2_mean"], harmonized, k_values),
        _variant_summary("source_unit_best_signal_sensitivity", variants["source_unit_best_signal_sensitivity"], harmonized, k_values),
        _variant_summary("source_unit_v1_2_mean_excluding_sensitivity", variants["source_unit_v1_2_mean_excluding_sensitivity"], filtered, k_values),
    ]
    summary_rows = []
    for item in summaries:
        row = {
            "variant": item["variant"],
            "n_input_rows": item["n_input_rows"],
            "n_input_studies": item["n_input_studies"],
            "n_input_source_units": item["n_input_source_units"],
            "n_scored_genes": item["n_scored_genes"],
            "recall_at_50": item["top50_recall"],
            "recall_at_100": item["top100_recall"],
            "top20": ";".join(item["top20"]),
        }
        summary_rows.append(row)
    summary_tsv = output_dir / "source_unit_ablation_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_tsv, sep="\t", index=False)

    rank_tsv = output_dir / "source_unit_ablation_gene_ranks.tsv"
    _rank_table(variants).fillna("NA").to_csv(rank_tsv, sep="\t", index=False)

    jaccard = {
        f"{left}_vs_{right}_top{k}": _jaccard(variants[left], variants[right], k)
        for left, right in [
            ("contrast_level_legacy", "source_unit_v1_2_mean"),
            ("source_unit_v1_2_mean", "source_unit_best_signal_sensitivity"),
            ("source_unit_v1_2_mean", "source_unit_v1_2_mean_excluding_sensitivity"),
            ("source_unit_best_signal_sensitivity", "source_unit_v1_2_mean_excluding_sensitivity"),
            ("contrast_level_legacy", "source_unit_v1_2_mean_excluding_sensitivity"),
        ]
        for k in [50, 100]
    }
    summary = {
        "score_versions_compared": [
            "contrast_level_legacy_pre_v1_1",
            "degora_score_v1_2_source_unit_mean",
            "source_unit_best_signal_sensitivity_diagnostic_only",
            "degora_score_v1_2_source_unit_mean_excluding_sensitivity",
        ],
        "harmonized_path": str(harmonized_path),
        "catalog_path": str(catalog_path) if catalog_path else "",
        "min_studies": min_studies,
        "excluded_sensitivity_study_ids": excluded,
        "variant_summaries": summaries,
        "topk_jaccard": jaccard,
        "summary_tsv": str(summary_tsv),
        "gene_rank_tsv": str(rank_tsv),
        "interpretation_guardrail": (
            "The legacy contrast-level variant is for ablation only. Manuscript-facing DEGORA scores "
            "should use source-unit-aware consensus to avoid counting related contrasts as independent evidence. "
            "The best-signal variant is a diagnostic-only sensitivity analysis for transient time-course signals "
            "and must not replace the primary mean aggregation without a permutation/null calibration. "
            f"Current source-unit rule: {SOURCE_UNIT_COLLAPSE_RULE}."
        ),
    }
    summary_json = output_dir / "source_unit_ablation_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    summary_md = output_dir / "source_unit_ablation_summary.md"
    lines = [
        "# Source-Unit Scoring Ablation",
        "",
        "This ablation compares legacy contrast-level aggregation with `degora_score_v1_2_source_unit_mean`.",
        "",
        f"Excluded sensitivity/provisional study IDs: {', '.join(excluded) if excluded else 'none'}",
        "",
        "| variant | scored genes | recall@50 | recall@100 | top 10 |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for item in summaries:
        lines.append(
            f"| {item['variant']} | {item['n_scored_genes']} | {item['top50_recall']:.2f} | "
            f"{item['top100_recall']:.2f} | {';'.join(item['top20'][:10])} |"
        )
    lines.extend(["", "## Top-K Jaccard", ""])
    for key, value in jaccard.items():
        lines.append(f"- `{key}`: {value}")
    lines.append("")
    summary_md.write_text("\n".join(lines))

    for artifact in [summary_json, summary_tsv, rank_tsv, summary_md]:
        write_source_sidecar(
            artifact,
            command,
            inputs=[harmonized_path, *([catalog_path] if catalog_path else [])],
            metadata={"generator": "source-unit-ablation", "score_version": "degora_score_v1_2_source_unit_mean"},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-studies", type=int, default=2)
    parser.add_argument("--exclude-study-id", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_source_unit_ablation.py"
        f" --harmonized {args.harmonized}"
        + (f" --catalog {args.catalog}" if args.catalog else "")
        + f" --output-dir {args.output_dir}"
        + f" --min-studies {args.min_studies}"
        + "".join(f" --exclude-study-id {study_id}" for study_id in args.exclude_study_id)
    )
    summary = write_ablation(
        args.harmonized,
        args.output_dir,
        catalog_path=args.catalog,
        min_studies=args.min_studies,
        explicit_exclude=args.exclude_study_id,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
