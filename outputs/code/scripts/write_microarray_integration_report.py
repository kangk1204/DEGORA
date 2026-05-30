#!/usr/bin/env python
"""Write a manuscript-facing report for DEGORA microarray integration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.provenance import write_source_sidecar


KEY_ANCHORS = ["RBM39", "TYMS"]


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _method_row(summary: pd.DataFrame, method_id: str) -> dict[str, Any]:
    subset = summary.loc[summary["method_id"].astype(str).eq(method_id)]
    if subset.empty:
        return {}
    row = subset.iloc[0].to_dict()
    for key in ["recall_at_10", "recall_at_20", "recall_at_50", "recall_at_100"]:
        if row.get(key) != "":
            row[key] = float(row[key])
    if row.get("n_rows") != "":
        row["n_rows"] = int(row["n_rows"])
    return row


def _anchor_rows(path: Path) -> list[dict[str, Any]]:
    anchors = pd.read_csv(path, sep="\t")
    rows: list[dict[str, Any]] = []
    for gene in KEY_ANCHORS:
        subset = anchors.loc[anchors["gene_symbol"].astype(str).str.upper().eq(gene)]
        if subset.empty:
            continue
        row = subset.iloc[0].to_dict()
        for key in ["degora_rank", "degora_score", "n_source_units", "n_contrasts_observed", "weighted_lfc"]:
            if key in row and pd.notna(row[key]):
                row[key] = float(row[key]) if key in {"degora_score", "weighted_lfc"} else int(row[key])
        rows.append(row)
    return rows


def _best_non_degora_recall(rows: list[dict[str, Any]]) -> float | None:
    values = []
    for row in rows:
        if str(row.get("method_id", "")).startswith("degora"):
            continue
        value = pd.to_numeric(pd.Series([row.get("recall_at_100")]), errors="coerce").iloc[0]
        if pd.notna(value):
            values.append(float(value))
    return max(values) if values else None


def build_report(
    *,
    preparation_summary: Path,
    sensitivity_summary: Path,
    comparison_summary: Path,
    comparator_summary: Path,
    score_db_summary: Path,
    anchor_ranks: Path,
) -> dict[str, Any]:
    preparation = _read_json(preparation_summary)
    sensitivity = _read_json(sensitivity_summary)
    primary_vs_microarray = _read_json(comparison_summary)
    score_db = _read_json(score_db_summary)
    comparator = pd.read_csv(comparator_summary)
    degora_row = _method_row(comparator, "degora_deg_score")
    rra_row = _method_row(comparator, "robustrankaggreg")
    rank_product_row = _method_row(comparator, "rank_product_approx")
    best_non_degora_at_100 = _best_non_degora_recall(comparator.to_dict(orient="records"))
    degora_at_100 = degora_row.get("recall_at_100")
    degora_delta = (
        float(degora_at_100) - float(best_non_degora_at_100)
        if degora_at_100 is not None and best_non_degora_at_100 is not None
        else None
    )
    benchmark_interpretation = (
        "DEGORA matched or exceeded the best non-DEGORA recall@100 comparator in this microarray sensitivity run."
        if degora_delta is not None and degora_delta >= 0
        else "DEGORA did not exceed the best non-DEGORA recall@100 comparator in this microarray sensitivity run; interpret the advantage as provenance, assay metadata, and evidence-DB usability rather than recall superiority."
    )

    e7820 = preparation.get("e7820_summary", {})
    sirbm39 = preparation.get("sirbm39_summary", {})
    anchor_details = _anchor_rows(anchor_ranks)
    return {
        "claim": (
            "DEGORA supports cross-platform RNA-seq plus microarray integration with explicit assay, "
            "normalization, probe-collapse, and direct-vs-sensitivity metadata. This is a workflow/provenance "
            "advantage, not a broad anchor-recall superiority result on the historical RBM39-axis sensitivity case."
        ),
        "nutriomics_db_implication": (
            "For a NutriOmics-style evidence database, the primary value is not a calibrated meta-analysis p-value. "
            "The useful output is whether a nutrient exposure repeatedly moves a gene in the same direction across "
            "independent source units, while preserving RNA-seq/microarray provenance, treatment context, and probe-collapse "
            "rules. DEGORA's source-unit support, direction concordance, rank/percentile labels, and evidence rows are "
            "therefore directly aligned with building a queryable nutrient-gene response atlas."
        ),
        "nutrient_gene_query_contract": [
            "exposure_id or nutrient keyword, with synonym and dose/time annotations kept outside the gene score.",
            "gene_symbol plus stable identifier mappings where available.",
            "direction_consensus_label and direction_concordance_percent as first-class fields.",
            "n_source_units, n_contrasts, assay_types, platforms, and source_input_types to show evidence breadth.",
            "degora_rank, top_percent_label, evidence_tier, and component scores for intuitive prioritization.",
            "per-source evidence rows with study identifiers, treatment context, direct-vs-mechanism labels, and probe-collapse metadata.",
        ],
        "claim_boundary": {
            "microarray_direct_excluded_drug": bool(preparation.get("microarray_direct_excluded_drug", False)),
            "primary_excluded_drug_result_replaced": bool(preparation.get("primary_excluded_drug_result_replaced", False)),
            "as_published_author_deg_validation": bool(preparation.get("as_published_author_deg_validation", False)),
            "interpretation": preparation.get("microarray_interpretation", ""),
        },
        "microarray_preparation": {
            "platform": e7820.get("platform", ""),
            "normalization": e7820.get("normalization", ""),
            "source_input_type": e7820.get("source_input_type", ""),
            "collapse_rule": e7820.get("collapse_rule", ""),
            "e7820_gene_rows": e7820.get("n_gene_rows"),
            "e7820_collapsed_probe_rows": e7820.get("n_collapsed_probe_rows"),
            "sirbm39_gene_rows": sirbm39.get("n_gene_rows"),
            "sirbm39_collapsed_probe_rows": sirbm39.get("n_collapsed_probe_rows"),
        },
        "integrated_score_db": {
            "n_gene_scores": score_db.get("n_gene_scores"),
            "n_source_units": score_db.get("n_source_units"),
            "n_contrasts": score_db.get("n_contrasts"),
            "n_evidence_rows": score_db.get("n_evidence_rows"),
            "top_genes": score_db.get("top_genes", [])[:10],
        },
        "anchor_evidence": {
            "degora_score_top100_recovery": degora_row.get("recall_at_100"),
            "best_non_degora_top100_recovery": best_non_degora_at_100,
            "degora_vs_best_non_degora_delta_at_100": degora_delta,
            "benchmark_interpretation": benchmark_interpretation,
            "rank_product_top100_recovery": rank_product_row.get("recall_at_100"),
            "robustrankaggreg_top100_recovery": rra_row.get("recall_at_100"),
            "primary_top100_recovery": primary_vs_microarray.get("primary", {}).get("anchor_recovery_top_100", {}),
            "microarray_sensitivity_top100_recovery": primary_vs_microarray.get("microarray_sensitivity", {}).get("anchor_recovery_top_100", {}),
            "key_anchor_rows": anchor_details,
        },
        "strengths_vs_comparator_tools": [
            "DEGORA keeps assay_type, platform, normalization, source_input_type, and probe_collapse metadata in the evidence DB.",
            "DEGORA collapses microarray probes before gene ranking, reducing probe-level duplicate inflation.",
            "DEGORA can label microarray rows as sensitivity/mechanism evidence instead of silently merging them as direct drug evidence.",
            "DEGORA exposes per-gene source units and contrast counts, so users can inspect cross-platform support for genes such as RBM39 and TYMS.",
            "Standard p-value/rank comparators can run after forced harmonization, but their output tables do not by themselves preserve these guardrails.",
        ],
        "limitations": [
            "GSE93829 is E7820/RBM39-axis microarray sensitivity evidence, not direct treatment evidence for the now-excluded drug benchmark.",
            "The current microarray DEG derivation uses Welch tests on normalized matrix values; limma full tables remain preferred when available.",
            "Broad anchor recall is not superior to all rank-based comparators in this historical microarray sensitivity case.",
        ],
        "inputs": {
            "preparation_summary": str(preparation_summary),
            "sensitivity_summary": str(sensitivity_summary),
            "comparison_summary": str(comparison_summary),
            "comparator_summary": str(comparator_summary),
            "score_db_summary": str(score_db_summary),
            "anchor_ranks": str(anchor_ranks),
        },
    }


def write_markdown(report: dict[str, Any], output: Path) -> None:
    prep = report["microarray_preparation"]
    score = report["integrated_score_db"]
    anchor = report["anchor_evidence"]
    lines = [
        "# Microarray Integration Strength Report",
        "",
        report["claim"],
        "",
        "## NutriOmics DB Implication",
        "",
        report["nutriomics_db_implication"],
        "",
        "## Nutrient-Gene Query Contract",
        "",
        *[f"- {item}" for item in report["nutrient_gene_query_contract"]],
        "",
        "## Claim Boundary",
        "",
        f"- Direct excluded-drug microarray evidence: `{report['claim_boundary']['microarray_direct_excluded_drug']}`",
        f"- Replaces primary excluded-drug result: `{report['claim_boundary']['primary_excluded_drug_result_replaced']}`",
        f"- Author DEG validation claim: `{report['claim_boundary']['as_published_author_deg_validation']}`",
        f"- Interpretation: {report['claim_boundary']['interpretation']}",
        "",
        "## Microarray Preparation",
        "",
        f"- Platform: `{prep['platform']}`",
        f"- Normalization: `{prep['normalization']}`",
        f"- Source input type: `{prep['source_input_type']}`",
        f"- Probe collapse rule: `{prep['collapse_rule']}`",
        f"- E7820 gene rows / collapsed probe rows: {prep['e7820_gene_rows']} / {prep['e7820_collapsed_probe_rows']}",
        f"- siRBM39 gene rows / collapsed probe rows: {prep['sirbm39_gene_rows']} / {prep['sirbm39_collapsed_probe_rows']}",
        "",
        "## Integrated Evidence DB",
        "",
        f"- Gene scores: {score['n_gene_scores']}",
        f"- Source units: {score['n_source_units']}",
        f"- Contrasts: {score['n_contrasts']}",
        f"- Evidence rows: {score['n_evidence_rows']}",
        f"- Top genes: {';'.join(score['top_genes'])}",
        "",
        "## Anchor Evidence",
        "",
        f"- DEGORA score top100 anchor recovery: {anchor['degora_score_top100_recovery']}",
        f"- Best non-DEGORA top100 anchor recovery: {anchor['best_non_degora_top100_recovery']}",
        f"- Interpretation: {anchor['benchmark_interpretation']}",
        f"- Rank-product top100 anchor recovery: {anchor['rank_product_top100_recovery']}",
        f"- RobustRankAggreg top100 anchor recovery: {anchor['robustrankaggreg_top100_recovery']}",
        "",
        "| gene | DEGORA rank | top percent | source units | direction | weighted LFC |",
        "| --- | ---: | --- | --- | --- | ---: |",
    ]
    for row in anchor["key_anchor_rows"]:
        lines.append(
            f"| {row.get('gene_symbol', '')} | {row.get('degora_rank', '')} | {row.get('top_percent_label', '')} | "
            f"{row.get('support_label', '')} | {row.get('direction_label', '')} | {float(row.get('weighted_lfc', 0.0)):.3f} |"
        )
    lines.extend(["", "## Strengths Versus Comparator Tools", ""])
    lines.extend(f"- {item}" for item in report["strengths_vs_comparator_tools"])
    lines.extend(["", "## Limitations", ""])
    lines.extend(f"- {item}" for item in report["limitations"])
    output.write_text("\n".join(lines) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-summary", type=Path, required=True)
    parser.add_argument("--sensitivity-summary", type=Path, required=True)
    parser.add_argument("--comparison-summary", type=Path, required=True)
    parser.add_argument("--comparator-summary", type=Path, required=True)
    parser.add_argument("--score-db-summary", type=Path, required=True)
    parser.add_argument("--anchor-ranks", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = build_report(
        preparation_summary=args.preparation_summary,
        sensitivity_summary=args.sensitivity_summary,
        comparison_summary=args.comparison_summary,
        comparator_summary=args.comparator_summary,
        score_db_summary=args.score_db_summary,
        anchor_ranks=args.anchor_ranks,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    write_markdown(report, args.output_md)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_microarray_integration_report.py"
        f" --preparation-summary {args.preparation_summary}"
        f" --sensitivity-summary {args.sensitivity_summary}"
        f" --comparison-summary {args.comparison_summary}"
        f" --comparator-summary {args.comparator_summary}"
        f" --score-db-summary {args.score_db_summary}"
        f" --anchor-ranks {args.anchor_ranks}"
        f" --output-json {args.output_json}"
        f" --output-md {args.output_md}"
    )
    inputs = [
        args.preparation_summary,
        args.sensitivity_summary,
        args.comparison_summary,
        args.comparator_summary,
        args.score_db_summary,
        args.anchor_ranks,
    ]
    for artifact in [args.output_json, args.output_md]:
        write_source_sidecar(artifact, command, inputs=inputs, metadata={"generator": "microarray-integration-report"})
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
