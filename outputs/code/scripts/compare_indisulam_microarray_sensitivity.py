#!/usr/bin/env python
"""Compare primary Indisulam DEGORA ranks with microarray sensitivity ranks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.provenance import write_source_sidecar


def _read_summary(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _read_anchor_table(path: Path, suffix: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    keep = [
        "gene_symbol",
        "degora_rank",
        "rank_label",
        "evidence_tier",
        "degora_score",
        "top_percent_label",
        "support_label",
        "direction_label",
        "n_source_units",
        "n_contrasts_observed",
        "source_units",
        "weighted_lfc",
    ]
    table = table[[column for column in keep if column in table.columns]].copy()
    rename = {column: f"{column}_{suffix}" for column in table.columns if column != "gene_symbol"}
    table = table.rename(columns=rename)
    table["gene_symbol"] = table["gene_symbol"].astype("string").str.upper()
    return table


def compare(
    primary_summary: Path,
    primary_ranks: Path,
    sensitivity_summary: Path,
    sensitivity_ranks: Path,
    output_json: Path,
    output_tsv: Path,
    *,
    command: str,
) -> dict[str, Any]:
    primary = _read_summary(primary_summary)
    sensitivity = _read_summary(sensitivity_summary)
    primary_table = _read_anchor_table(primary_ranks, "primary")
    sensitivity_table = _read_anchor_table(sensitivity_ranks, "microarray_sensitivity")
    comparison = primary_table.merge(sensitivity_table, on="gene_symbol", how="outer")
    for column in ["degora_rank_primary", "degora_rank_microarray_sensitivity"]:
        if column in comparison.columns:
            comparison[column] = pd.to_numeric(comparison[column], errors="coerce")
    comparison["rank_delta_microarray_minus_primary"] = (
        comparison["degora_rank_microarray_sensitivity"] - comparison["degora_rank_primary"]
    )
    comparison = comparison.sort_values(
        ["degora_rank_primary", "degora_rank_microarray_sensitivity", "gene_symbol"],
        na_position="last",
    )

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output_tsv, sep="\t", index=False)

    summary = {
        "comparison_scope": "primary_indisulam_vs_e7820_rbm39_axis_microarray_sensitivity",
        "primary_summary": str(primary_summary),
        "sensitivity_summary": str(sensitivity_summary),
        "microarray_direct_indisulam": bool(sensitivity.get("microarray_direct_indisulam", False)),
        "primary_indisulam_result_replaced": bool(sensitivity.get("primary_indisulam_result_replaced", False)),
        "as_published_author_deg_validation": bool(sensitivity.get("as_published_author_deg_validation", False)),
        "interpretation": (
            "The microarray analysis is a cross-platform sensitivity stress test. "
            "Because GSE93829 is E7820/RBM39-axis data rather than direct Indisulam treatment, "
            "rank changes should not replace the primary Indisulam score."
        ),
        "primary": {
            "n_scored_genes": primary.get("n_scored_genes"),
            "anchor_recovery_top_100": primary.get("anchor_recovery", {}).get("100", {}),
        },
        "microarray_sensitivity": {
            "n_scored_genes": sensitivity.get("n_scored_genes"),
            "anchor_recovery_top_100": sensitivity.get("anchor_recovery", {}).get("100", {}),
        },
        "anchor_rank_comparison_tsv": str(output_tsv),
    }
    output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    for artifact in [output_json, output_tsv]:
        write_source_sidecar(
            artifact,
            command,
            inputs=[primary_summary, primary_ranks, sensitivity_summary, sensitivity_ranks],
            metadata={"generator": "indisulam-microarray-sensitivity-comparison"},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-summary", type=Path, required=True)
    parser.add_argument("--primary-ranks", type=Path, required=True)
    parser.add_argument("--sensitivity-summary", type=Path, required=True)
    parser.add_argument("--sensitivity-ranks", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/compare_indisulam_microarray_sensitivity.py"
        f" --primary-summary {args.primary_summary}"
        f" --primary-ranks {args.primary_ranks}"
        f" --sensitivity-summary {args.sensitivity_summary}"
        f" --sensitivity-ranks {args.sensitivity_ranks}"
        f" --output-json {args.output_json}"
        f" --output-tsv {args.output_tsv}"
    )
    summary = compare(
        args.primary_summary,
        args.primary_ranks,
        args.sensitivity_summary,
        args.sensitivity_ranks,
        args.output_json,
        args.output_tsv,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
