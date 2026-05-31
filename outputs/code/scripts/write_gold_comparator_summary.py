#!/usr/bin/env python
"""Summarize comparator recovery against an arbitrary locked gold panel."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from degora.baselines import BASELINE_RESULT_COLUMNS, baseline_result_paths
from degora.benchmark_stats import recall_stats
from degora.provenance import shell_command, write_source_sidecar
from degora.stats import bh_adjust


SUMMARY_COLUMNS = [
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
    "direction_recall_at_10",
    "direction_recall_at_20",
    "direction_recall_at_50",
    "direction_recall_at_100",
    "recovered_at_100",
    "direction_recovered_at_100",
    "direction_mismatched_at_100",
    "missing_at_100",
    "top10",
    "failure_mode",
]


def _join(values: Iterable[str]) -> str:
    return ";".join(values)


def _normalize_direction(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"up", "+", "1", "increase", "increased", "upregulated", "up-regulated"}:
        return "up"
    if text in {"down", "-", "-1", "decrease", "decreased", "downregulated", "down-regulated"}:
        return "down"
    return ""


def _read_gold(path: Path, gene_column: str) -> tuple[set[str], dict[str, str]]:
    gold = pd.read_csv(path)
    if gene_column not in gold.columns:
        raise ValueError(f"gold file {path} missing gene column {gene_column!r}")
    genes = gold[gene_column].dropna().astype(str).str.strip().str.upper()
    positives = {gene for gene in genes if gene}
    if not positives:
        raise ValueError(f"gold file {path} has no positive genes in column {gene_column!r}")
    directions: dict[str, str] = {}
    if "expected_direction" in gold.columns:
        gold = gold.copy()
        gold["_gene"] = gold[gene_column].astype(str).str.strip().str.upper()
        gold["_direction"] = gold["expected_direction"].map(_normalize_direction)
        directions = {
            str(gene): str(direction)
            for gene, direction in zip(gold["_gene"], gold["_direction"], strict=False)
            if str(gene) in positives and str(direction)
        }
    return positives, directions


def _ranked_records(frame: pd.DataFrame, symbol_column: str, rank_column: str, direction_column: str | None) -> list[tuple[str, str]]:
    if symbol_column not in frame.columns:
        raise ValueError(f"result frame missing symbol column {symbol_column!r}")
    if rank_column not in frame.columns:
        raise ValueError(f"result frame missing rank column {rank_column!r}")
    ranked = frame.dropna(subset=[symbol_column, rank_column]).copy()
    ranked[rank_column] = pd.to_numeric(ranked[rank_column], errors="coerce")
    ranked = ranked.dropna(subset=[rank_column])
    ranked = ranked.sort_values([rank_column, symbol_column])
    symbols = ranked[symbol_column].astype(str).str.strip().str.upper().tolist()
    if direction_column and direction_column in ranked.columns:
        directions = ranked[direction_column].map(_normalize_direction).tolist()
    else:
        directions = [""] * len(symbols)
    return list(zip(symbols, directions, strict=False))


def _recall(symbols: list[str], positives: set[str], k: int) -> tuple[float, list[str], list[str]]:
    top = set(symbols[:k])
    recovered = sorted(positives.intersection(top))
    missing = sorted(positives.difference(top))
    return len(recovered) / len(positives) if positives else 0.0, recovered, missing


def _direction_recall(
    ranked: list[tuple[str, str]],
    expected_directions: dict[str, str],
    k: int,
) -> tuple[float | str, list[str], list[str]]:
    if not expected_directions:
        return "", [], []
    expected = set(expected_directions)
    top = ranked[:k]
    recovered: list[str] = []
    mismatched: list[str] = []
    for symbol, direction in top:
        expected_direction = expected_directions.get(symbol)
        if not expected_direction:
            continue
        if direction == expected_direction:
            recovered.append(symbol)
        else:
            mismatched.append(f"{symbol}:{direction or 'unknown'}!={expected_direction}")
    return len(set(recovered)) / len(expected), sorted(set(recovered)), sorted(set(mismatched))


def _summary_record(
    *,
    method_id: str,
    setting_id: str,
    ranked: list[tuple[str, str]],
    positives: set[str],
    expected_directions: dict[str, str],
    n_rows: int,
    run_status: str = "ok",
    failure_mode: str = "",
) -> dict[str, object]:
    ranked_symbols = [symbol for symbol, _ in ranked]
    recalls = {k: _recall(ranked_symbols, positives, k) for k in (10, 20, 50, 100)}
    stats = {k: recall_stats(ranked_symbols, positives, k, universe_size=n_rows) for k in (10, 20, 50, 100)}
    direction_recalls = {k: _direction_recall(ranked, expected_directions, k) for k in (10, 20, 50, 100)}
    return {
        "method_id": method_id,
        "setting_id": setting_id,
        "run_status": run_status,
        "n_rows": n_rows,
        "n_gold": len(positives),
        "n_recovered_at_10": stats[10].n_recovered,
        "recall_at_10": recalls[10][0],
        "recall_at_10_ci_low": stats[10].ci_low,
        "recall_at_10_ci_high": stats[10].ci_high,
        "precision_at_10": stats[10].precision,
        "hypergeom_p_at_10": stats[10].hypergeom_pvalue,
        "hypergeom_fdr_at_10": "",
        "n_recovered_at_20": stats[20].n_recovered,
        "recall_at_20": recalls[20][0],
        "recall_at_20_ci_low": stats[20].ci_low,
        "recall_at_20_ci_high": stats[20].ci_high,
        "precision_at_20": stats[20].precision,
        "hypergeom_p_at_20": stats[20].hypergeom_pvalue,
        "hypergeom_fdr_at_20": "",
        "n_recovered_at_50": stats[50].n_recovered,
        "recall_at_50": recalls[50][0],
        "recall_at_50_ci_low": stats[50].ci_low,
        "recall_at_50_ci_high": stats[50].ci_high,
        "precision_at_50": stats[50].precision,
        "hypergeom_p_at_50": stats[50].hypergeom_pvalue,
        "hypergeom_fdr_at_50": "",
        "n_recovered_at_100": stats[100].n_recovered,
        "recall_at_100": recalls[100][0],
        "recall_at_100_ci_low": stats[100].ci_low,
        "recall_at_100_ci_high": stats[100].ci_high,
        "precision_at_100": stats[100].precision,
        "hypergeom_p_at_100": stats[100].hypergeom_pvalue,
        "hypergeom_fdr_at_100": "",
        "direction_recall_at_10": direction_recalls[10][0],
        "direction_recall_at_20": direction_recalls[20][0],
        "direction_recall_at_50": direction_recalls[50][0],
        "direction_recall_at_100": direction_recalls[100][0],
        "recovered_at_100": _join(recalls[100][1]),
        "direction_recovered_at_100": _join(direction_recalls[100][1]),
        "direction_mismatched_at_100": _join(direction_recalls[100][2]),
        "missing_at_100": _join(recalls[100][2]),
        "top10": _join(ranked_symbols[:10]),
        "failure_mode": failure_mode,
    }


def _baseline_rows(baseline_dir: Path, positives: set[str], expected_directions: dict[str, str]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    required = set(BASELINE_RESULT_COLUMNS)
    for path in baseline_result_paths(baseline_dir):
        frame = pd.read_csv(path, sep="\t")
        if not required.issubset(frame.columns) or frame.empty:
            continue
        method_id = str(frame["method_id"].iloc[0])
        setting_id = str(frame["setting_id"].iloc[0])
        status_values = set(frame["status"].dropna().astype(str))
        if status_values and status_values != {"ok"}:
            records.append(
                {
                    "method_id": method_id,
                    "setting_id": setting_id,
                    "run_status": "blocked",
                    "n_rows": 0,
                    "recall_at_10": "",
                    "recall_at_20": "",
                    "recall_at_50": "",
                    "recall_at_100": "",
                    "direction_recall_at_10": "",
                    "direction_recall_at_20": "",
                    "direction_recall_at_50": "",
                    "direction_recall_at_100": "",
                    "recovered_at_100": "",
                    "direction_recovered_at_100": "",
                    "direction_mismatched_at_100": "",
                    "missing_at_100": "",
                    "top10": "",
                    "failure_mode": "baseline_result_status_" + "_".join(sorted(status_values)),
                }
            )
            continue
        ranked = _ranked_records(frame, "symbol", "rank", "direction")
        records.append(
            _summary_record(
                method_id=method_id,
                setting_id=setting_id,
                ranked=ranked,
                positives=positives,
                expected_directions=expected_directions,
                n_rows=int(len(frame)),
            )
        )
    return records


def _degora_score_row(
    score_csv: Path | None,
    positives: set[str],
    expected_directions: dict[str, str],
    label: str,
) -> list[dict[str, object]]:
    if score_csv is None:
        return []
    frame = pd.read_csv(score_csv)
    if frame.empty:
        return []
    records: list[dict[str, object]] = []
    ranked = _ranked_records(frame, "gene_symbol", "degora_rank", "consensus_direction")
    records.append(
        _summary_record(
            method_id="degora_deg_score",
            setting_id=label,
            ranked=ranked,
            positives=positives,
            expected_directions=expected_directions,
            n_rows=int(len(frame)),
        )
    )
    quality_columns = {"quality_weighted_degora_rank", "quality_weighted_consensus_direction"}
    if quality_columns.issubset(frame.columns):
        quality_ranked = _ranked_records(
            frame,
            "gene_symbol",
            "quality_weighted_degora_rank",
            "quality_weighted_consensus_direction",
        )
        records.append(
            _summary_record(
                method_id="degora_quality_weighted_score",
                setting_id="quality_weighted_secondary",
                ranked=quality_ranked,
                positives=positives,
                expected_directions=expected_directions,
                n_rows=int(len(frame)),
            )
        )
    return records


def _blocked_rows(baseline_dir: Path, emitted: set[tuple[str, str]]) -> list[dict[str, object]]:
    failure_path = baseline_dir / "baseline_failure_ledger.csv"
    if not failure_path.exists():
        return []
    failure = pd.read_csv(failure_path)
    records: list[dict[str, object]] = []
    for row in failure.itertuples(index=False):
        key = (str(row.method_id), str(row.setting_id))
        if key in emitted:
            continue
        records.append(
            {
                "method_id": str(row.method_id),
                "setting_id": str(row.setting_id),
                "run_status": "blocked",
                "n_rows": 0,
                "recall_at_10": "",
                "recall_at_20": "",
                "recall_at_50": "",
                "recall_at_100": "",
                "direction_recall_at_10": "",
                "direction_recall_at_20": "",
                "direction_recall_at_50": "",
                "direction_recall_at_100": "",
                "recovered_at_100": "",
                "direction_recovered_at_100": "",
                "direction_mismatched_at_100": "",
                "missing_at_100": "",
                "top10": "",
                "failure_mode": str(row.blocker_id),
            }
        )
    return records


def build_summary(
    baseline_dir: Path,
    gold_path: Path,
    *,
    gold_gene_column: str = "gene_symbol",
    degora_score_csv: Path | None = None,
    degora_score_label: str = "score",
) -> pd.DataFrame:
    positives, expected_directions = _read_gold(gold_path, gold_gene_column)
    records = [
        *_degora_score_row(degora_score_csv, positives, expected_directions, degora_score_label),
        *_baseline_rows(baseline_dir, positives, expected_directions),
    ]
    emitted = {(str(row["method_id"]), str(row["setting_id"])) for row in records}
    records.extend(_blocked_rows(baseline_dir, emitted))
    frame = pd.DataFrame.from_records(records, columns=SUMMARY_COLUMNS)
    if frame.empty:
        return frame
    order = {
        "degora_deg_score": 0,
        "degora_quality_weighted_score": 1,
        "degora_slice": 2,
        "weighted_stouffer": 3,
        "unweighted_stouffer": 4,
        "fisher": 5,
        "rank_product_approx": 6,
        "maic": 6.5,
        "sign_vote": 7,
        "metavolcanor": 8,
        "robustrankaggreg": 9,
        "awfisher": 10,
        "metarnaseq_fisher": 11,
        "metarnaseq_invnorm": 12,
        "hstouffer": 13,
        "awmeta": 14,
        "rankprod_exact": 15,
        "metade": 16,
        "dexma": 17,
        "metaintegrator": 18,
    }
    frame["_order"] = frame["method_id"].map(order).fillna(999)
    for k in (10, 20, 50, 100):
        p_column = f"hypergeom_p_at_{k}"
        fdr_column = f"hypergeom_fdr_at_{k}"
        frame[fdr_column] = np.nan
        ok = frame["run_status"].eq("ok") & pd.to_numeric(frame[p_column], errors="coerce").notna()
        frame.loc[ok, fdr_column] = bh_adjust(pd.to_numeric(frame.loc[ok, p_column], errors="coerce").to_numpy(dtype=float))
    return frame.sort_values(["_order", "setting_id"]).drop(columns=["_order"]).reset_index(drop=True)


def write_markdown(summary: pd.DataFrame, output: Path, *, title: str, gold_path: Path) -> None:
    lines = [
        f"# {title}",
        "",
        f"Locked positives: `{gold_path}`.",
        "",
        "| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |",
        "| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in summary.itertuples(index=False):
        def fmt(value: object) -> str:
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            return "" if pd.isna(numeric) else f"{float(numeric):.2f}"

        lines.append(
            f"| {row.method_id}/{row.setting_id} | {row.run_status} | {row.n_rows} | "
            f"{fmt(row.recall_at_10)} ({fmt(row.recall_at_10_ci_low)}-{fmt(row.recall_at_10_ci_high)}) | "
            f"{fmt(row.recall_at_20)} ({fmt(row.recall_at_20_ci_low)}-{fmt(row.recall_at_20_ci_high)}) | "
            f"{fmt(row.recall_at_50)} ({fmt(row.recall_at_50_ci_low)}-{fmt(row.recall_at_50_ci_high)}) | "
            f"{fmt(row.recall_at_100)} ({fmt(row.recall_at_100_ci_low)}-{fmt(row.recall_at_100_ci_high)}) | "
            f"{fmt(row.precision_at_100)} | {fmt(row.hypergeom_fdr_at_100)} | {fmt(row.direction_recall_at_100)} | "
            f"{row.failure_mode} | {row.top10} |"
        )
    blocked = set(summary.loc[summary["run_status"].eq("blocked"), "method_id"].astype(str))
    if {"hstouffer", "awmeta"}.intersection(blocked):
        lines.extend(
            [
                "",
                "Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.",
            ]
        )
    output.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--gold-gene-column", default="gene_symbol")
    parser.add_argument("--degora-score-csv", type=Path)
    parser.add_argument("--degora-score-label", default="v1_2_source_unit_mean")
    parser.add_argument("--title", default="Gold-Panel Comparator Summary")
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args(argv)

    summary = build_summary(
        args.baseline_dir,
        args.gold,
        gold_gene_column=args.gold_gene_column,
        degora_score_csv=args.degora_score_csv,
        degora_score_label=args.degora_score_label,
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    write_markdown(summary, args.output_md, title=args.title, gold_path=args.gold)
    command_args = [
        "env",
        "PYTHONPATH=outputs/code",
        "python",
        "outputs/code/scripts/write_gold_comparator_summary.py",
        "--baseline-dir",
        args.baseline_dir,
        "--gold",
        args.gold,
        "--gold-gene-column",
        args.gold_gene_column,
    ]
    if args.degora_score_csv:
        command_args.extend(["--degora-score-csv", args.degora_score_csv])
    command_args.extend(
        [
            "--degora-score-label",
            args.degora_score_label,
            "--title",
            args.title,
            "--output-csv",
            args.output_csv,
            "--output-md",
            args.output_md,
        ]
    )
    command = shell_command(command_args)
    inputs = [args.baseline_dir, args.gold]
    if args.degora_score_csv:
        inputs.append(args.degora_score_csv)
    for artifact in (args.output_csv, args.output_md):
        write_source_sidecar(
            artifact,
            command,
            inputs=inputs,
            metadata={"generator": "gold-comparator-summary", "gold": str(args.gold)},
        )
    print(
        json.dumps(
            {
                "output_csv": str(args.output_csv),
                "output_md": str(args.output_md),
                "rows": int(len(summary)),
                "ok_rows": int(summary["run_status"].eq("ok").sum()) if "run_status" in summary.columns else 0,
                "blocked_rows": int(summary["run_status"].eq("blocked").sum()) if "run_status" in summary.columns else 0,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
