#!/usr/bin/env python
"""Write manuscript table and supplementary table index for DEGORA."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

from degora.provenance import shell_command, write_source_sidecar


DEFAULT_OUTPUT = Path("outputs/tables/manuscript")

CORE_SUMMARIES = [
    ("IFN", "RNA-seq", Path("outputs/results/ifn-pilot/ifn_gold_comparator_summary.csv")),
    ("ER stress", "RNA-seq + microarray", Path("outputs/results/er-stress-cross-platform/er_stress_cross_platform_gold_comparator_summary.csv")),
    ("Heat shock", "RNA-seq", Path("outputs/results/heat-shock-benchmark/heat_shock_gold_comparator_summary.csv")),
    ("Hypoxia", "RNA-seq", Path("outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv")),
]


def _fixed_degora_rows() -> pd.DataFrame:
    rows = []
    for topic, input_mode, path in CORE_SUMMARIES:
        frame = pd.read_csv(path)
        for col in ["recall_at_10", "recall_at_20", "recall_at_50", "recall_at_100", "direction_recall_at_100"]:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
        degora = frame.loc[frame["method_id"].eq("degora_quality_weighted_score")].copy()
        if degora.empty:
            degora = frame.loc[frame["method_id"].eq("degora_deg_score")].copy()
        best_other = (
            frame.loc[~frame["method_id"].isin(["degora_deg_score", "degora_quality_weighted_score"]) & frame["run_status"].eq("ok")]
            .sort_values(["recall_at_100", "recall_at_50", "recall_at_10"], ascending=False)
            .head(1)
        )
        if degora.empty:
            raise ValueError(f"{path} does not contain a DEGORA score row")
        row = degora.iloc[0]
        has_other = not best_other.empty
        other = best_other.iloc[0] if has_other else None
        rows.append(
            {
                "topic": topic,
                "input_mode": input_mode,
                "degora_method": f"{row.method_id}/{row.setting_id}",
                "degora_recall_at_10": row.recall_at_10,
                "degora_recall_at_50": row.recall_at_50,
                "degora_recall_at_100": row.recall_at_100,
                "degora_recall_at_100_ci_low": getattr(row, "recall_at_100_ci_low", ""),
                "degora_recall_at_100_ci_high": getattr(row, "recall_at_100_ci_high", ""),
                "degora_precision_at_100": getattr(row, "precision_at_100", ""),
                "degora_hypergeom_fdr_at_100": getattr(row, "hypergeom_fdr_at_100", ""),
                "best_non_degora_method": f"{other.method_id}/{other.setting_id}" if has_other else "none",
                "best_non_degora_recall_at_10": other.recall_at_10 if has_other else "",
                "best_non_degora_recall_at_50": other.recall_at_50 if has_other else "",
                "best_non_degora_recall_at_100": other.recall_at_100 if has_other else "",
                "best_non_degora_recall_at_100_ci_low": getattr(other, "recall_at_100_ci_low", "") if has_other else "",
                "best_non_degora_recall_at_100_ci_high": getattr(other, "recall_at_100_ci_high", "") if has_other else "",
                "best_non_degora_precision_at_100": getattr(other, "precision_at_100", "") if has_other else "",
                "best_non_degora_hypergeom_fdr_at_100": getattr(other, "hypergeom_fdr_at_100", "") if has_other else "",
                "degora_top10": row.top10,
            }
        )
    return pd.DataFrame(rows)


def _supplementary_index() -> pd.DataFrame:
    rows = [
        {
            "supplement": "Supplementary Table 1",
            "title": "Public-summary comparator input requirements",
            "path": "outputs/results/comparator_public_summary_input_requirements.csv",
            "purpose": "Detailed corpus-level runnable/blocked comparator and prior-art resource table.",
        },
        {
            "supplement": "Supplementary Table 2",
            "title": "Prior-art coverage summary",
            "path": "outputs/results/prior_art_coverage_summary.csv",
            "purpose": "One-row-per-method/resource summary for Related Work and reviewer response.",
        },
        {
            "supplement": "Supplementary Table 3",
            "title": "Cross-platform benchmark summary",
            "path": "outputs/results/cross-platform-microarray/cross_platform_benchmark_summary.csv",
            "purpose": "RNA-seq-only versus RNA-seq+microarray DEGORA recall comparison.",
        },
        {
            "supplement": "Supplementary Table 4",
            "title": "Cross-platform marker rank shifts",
            "path": "outputs/results/cross-platform-microarray/cross_platform_marker_rank_shifts.csv",
            "purpose": "Locked marker rank and top-percent changes after adding microarray evidence.",
        },
        {
            "supplement": "Supplementary Table 5",
            "title": "Cross-platform dataset sources",
            "path": "outputs/results/cross-platform-microarray/cross_platform_dataset_sources.csv",
            "purpose": "Microarray accessions, platforms, samples, and contrasts added to mixed benchmarks.",
        },
        {
            "supplement": "Supplementary Table 6",
            "title": "Source-unit feasibility matrices",
            "path": "outputs/results/*/tool-feasibility/*_study_tool_feasibility_matrix.csv",
            "purpose": "Per-source evidence for which SOTA tools require missing raw-expression or metadata inputs.",
        },
        {
            "supplement": "Supplementary Table 7",
            "title": "Statistical uncertainty and auxiliary reporting lanes",
            "path": "outputs/results/*/deep-metrics/*_deep_point_metrics.csv; outputs/results/*/degora_gene_scores.csv",
            "purpose": "Recall exact intervals, hypergeometric enrichment/FDR, background-negative AUROC/AUPRC, RE-Stouffer, RRA rho, direction posterior, and effect-size random-effects reporting fields.",
        },
        {
            "supplement": "Supplementary Table 8",
            "title": "Resource feature comparison versus public DEG/evidence portals",
            "path": "outputs/tables/manuscript/Supplementary_Table8_resource_feature_comparison.csv",
            "purpose": "DEGORA Evidence Explorer feature matrix against DEET, CREEDS/RummaGEO, ARCHS4, Enrichr/SigCom LINCS, and Open Targets on provenance, interactivity, and FAIR/offline-deployment dimensions.",
        },
    ]
    return pd.DataFrame(rows)


def _resource_feature_comparison() -> pd.DataFrame:
    """Curated resource feature matrix for the Evidence Explorer (Supplementary Table 8).

    Columns are the SOTA dimensions surfaced by the resource-landscape review; rows
    are DEGORA and the closest public DEG/evidence portals. Cells are intentionally
    conservative and verifiable from each resource's public interface.
    """
    dimensions = [
        "per_source_evidence_drilldown",
        "source_url_on_every_evidence_row",
        "direction_aware_filtering",
        "heterogeneity_or_confidence_surfaced",
        "effect_size_with_CI",
        "client_side_SQL_over_full_schema",
        "fully_static_offline_no_server",
        "gene_permalinks",
        "bulk_and_per_gene_download",
        "self_contained_single_file",
    ]
    resources = {
        "DEGORA Evidence Explorer": ["yes"] * 10,
        "DEET (Sokolowski 2023)": ["yes", "yes (PMID)", "partial", "no", "no", "no", "no", "partial", "yes", "no"],
        "CREEDS / RummaGEO": ["partial", "yes (GSE)", "yes (up/down)", "no", "no", "no", "no", "partial", "yes", "no"],
        "ARCHS4": ["partial", "yes (GSM)", "no", "no", "no", "no", "no", "yes", "yes", "no"],
        "Enrichr / SigCom LINCS": ["partial", "partial", "yes", "no", "no", "no", "no", "yes", "yes", "no"],
        "Open Targets Platform": ["yes", "yes", "no", "yes (datatype score)", "no", "no", "no", "yes", "yes", "no"],
    }
    rows = []
    for resource, values in resources.items():
        row = {"resource": resource}
        row.update(dict(zip(dimensions, values, strict=True)))
        rows.append(row)
    return pd.DataFrame(rows)


def _figure_index() -> pd.DataFrame:
    rows = [
        {
            "figure": "Figure 1",
            "title": "DEGORA workflow scheme",
            "primary_png": "outputs/figures/manuscript/FIGURE_1_SCHEME/figure1_scheme.png",
            "source_data": "outputs/figures/manuscript/FIGURE_1_SCHEME/figure1_scheme_source_data.xlsx",
            "legend": "outputs/figures/manuscript/FIGURE_1_SCHEME/figure1_scheme_legend.docx",
        },
        {
            "figure": "Figure 2",
            "title": "Benchmark against summary-level comparators",
            "primary_png": "outputs/figures/manuscript/FIGURE_2_BENCHMARK/figure2_benchmark.png",
            "source_data": "outputs/figures/manuscript/FIGURE_2_BENCHMARK/figure2_benchmark_source_data.xlsx",
            "legend": "outputs/figures/manuscript/FIGURE_2_BENCHMARK/figure2_benchmark_legend.docx",
        },
        {
            "figure": "Figure 3",
            "title": "RNA-seq plus microarray cross-platform benchmark",
            "primary_png": "outputs/figures/manuscript/DEGORA_CROSS_PLATFORM/degora_cross_platform.png",
            "source_data": "outputs/figures/manuscript/DEGORA_CROSS_PLATFORM/degora_cross_platform_source_data.xlsx",
            "legend": "outputs/figures/manuscript/DEGORA_CROSS_PLATFORM/degora_cross_platform_legend.docx",
        },
        {
            "figure": "Figure 4",
            "title": "Searchable DEGORA evidence atlas and interactive Evidence Explorer",
            "primary_png": "outputs/figures/manuscript/DEGORA_ATLAS/degora_atlas.png",
            "source_data": "outputs/figures/manuscript/DEGORA_ATLAS/degora_atlas_source_data.xlsx",
            "legend": "outputs/figures/manuscript/DEGORA_ATLAS/degora_atlas_legend.docx",
        },
        {
            "figure": "Interactive resource",
            "title": "DEGORA Evidence Explorer (fully static, offline, client-side SQL)",
            "primary_png": "outputs/figures/manuscript/DEGORA_EVIDENCE_EXPLORER/degora_evidence_explorer.html",
            "source_data": "outputs/figures/manuscript/DEGORA_EVIDENCE_EXPLORER/degora_evidence_explorer.db",
            "legend": "outputs/figures/manuscript/DEGORA_EVIDENCE_EXPLORER/degora_evidence_explorer_manifest.json",
        },
    ]
    return pd.DataFrame(rows)


def _write_markdown(frame: pd.DataFrame, path: Path) -> None:
    lines = ["| " + " | ".join(frame.columns) + " |", "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    for row in frame.itertuples(index=False):
        lines.append("| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |")
    path.write_text("\n".join(lines) + "\n")


def _write_xlsx(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    wb = Workbook()
    first = True
    for name, frame in sheets.items():
        ws = wb.active if first else wb.create_sheet()
        first = False
        ws.title = name[:31]
        for row in dataframe_to_rows(frame, index=False, header=True):
            ws.append(row)
    wb.save(path)


def write_tables(output_dir: Path, command: str) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    table1 = _fixed_degora_rows()
    table2 = pd.read_csv("outputs/results/cross-platform-microarray/cross_platform_benchmark_summary.csv")
    supp_index = _supplementary_index()
    figure_index = _figure_index()
    feature_comparison = _resource_feature_comparison()
    workbook = output_dir / "degora_manuscript_tables.xlsx"
    outputs = [
        output_dir / "Figure_Index.csv",
        output_dir / "Figure_Index.md",
        output_dir / "Table1_core_benchmark_summary.csv",
        output_dir / "Table1_core_benchmark_summary.md",
        output_dir / "Table2_cross_platform_summary.csv",
        output_dir / "Table2_cross_platform_summary.md",
        output_dir / "Supplementary_Table_Index.csv",
        output_dir / "Supplementary_Table_Index.md",
        output_dir / "Supplementary_Table8_resource_feature_comparison.csv",
        output_dir / "Supplementary_Table8_resource_feature_comparison.md",
        workbook,
        output_dir / "degora_manuscript_tables_manifest.json",
    ]
    figure_index.to_csv(outputs[0], index=False)
    _write_markdown(figure_index, outputs[1])
    table1.to_csv(outputs[2], index=False)
    _write_markdown(table1, outputs[3])
    table2.to_csv(outputs[4], index=False)
    _write_markdown(table2, outputs[5])
    supp_index.to_csv(outputs[6], index=False)
    _write_markdown(supp_index, outputs[7])
    feature_comparison.to_csv(outputs[8], index=False)
    _write_markdown(feature_comparison, outputs[9])
    _write_xlsx(workbook, {"Figure_Index": figure_index, "Table1_core_benchmark": table1, "Table2_cross_platform": table2, "Supplementary_Index": supp_index, "Resource_Feature_Comparison": feature_comparison})
    manifest_path = output_dir / "degora_manuscript_tables_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "generated_at": "deterministic",
                "script": "outputs/code/scripts/write_manuscript_tables.py",
                "command": command,
                "outputs": [str(path) for path in outputs],
                "inputs": [str(path) for _, _, path in CORE_SUMMARIES]
                + ["outputs/results/cross-platform-microarray/cross_platform_benchmark_summary.csv"],
                "known_limitations": [
                    "Tables summarize locked-panel recall with exact uncertainty and list-level enrichment; they do not convert DEGORA scores into calibrated posterior probabilities.",
                    "Direction-aware benchmark panels remain dominated by up-regulated canonical markers, so down-direction performance requires expanded gold panels.",
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    inputs = [path for _, _, path in CORE_SUMMARIES] + [Path("outputs/results/cross-platform-microarray/cross_platform_benchmark_summary.csv")]
    for artifact in outputs:
        write_source_sidecar(artifact, command, inputs=inputs, metadata={"generator": "manuscript-table-index"})
    return {
        "outputs": [str(path) for path in outputs],
        "figure_rows": len(figure_index),
        "table1_rows": len(table1),
        "table2_rows": len(table2),
        "supplementary_rows": len(supp_index),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    command = shell_command(["PYTHONPATH=outputs/code", "python", "outputs/code/scripts/write_manuscript_tables.py", "--output-dir", args.output_dir])
    print(json.dumps(write_tables(args.output_dir, command), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
