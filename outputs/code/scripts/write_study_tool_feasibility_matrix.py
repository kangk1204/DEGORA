#!/usr/bin/env python
"""Write source-level feasibility rows for public-summary and raw-input tools."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.baselines import (
    HSTOUFFER_REQUIRED_ORIGINAL_COLUMNS,
    awmeta_source_input_audit,
    hstouffer_materializer_feasibility,
)
from degora.provenance import shell_command, write_source_sidecar
from degora.score_db import _source_unit_series


FEASIBILITY_COLUMNS = [
    "corpus",
    "source_unit_id",
    "study_id",
    "pipeline",
    "assay_type",
    "source_input_type",
    "table_scope",
    "source_path",
    "n_rows",
    "has_gene_symbol",
    "has_lfc",
    "has_pvalue",
    "has_padj",
    "has_rank",
    "has_sample_size",
    "awfisher",
    "metarnaseq_fisher",
    "robustrankaggreg",
    "metavolcanor",
    "rankprod_exact",
    "metade_pvalue_modes",
    "metade_effect_size_modes",
    "dexma",
    "metaintegrator",
    "omicc",
    "imageo",
    "networkanalyst",
    "crossmeta",
    "deet",
    "creeds",
    "generic_pvalue_combiners",
    "hstouffer",
    "awmeta",
    "missing_for_raw_expression_tools",
    "notes",
]


def _join_unique(values: pd.Series) -> str:
    out: list[str] = []
    for value in values.dropna().astype(str):
        label = value.strip()
        if label and label.lower() not in {"nan", "<na>"}:
            out.append(label)
    return ";".join(sorted(dict.fromkeys(out)))


def _has_nonmissing(frame: pd.DataFrame, column: str) -> bool:
    return column in frame.columns and frame[column].notna().any()


def _status(ok: bool, reason: str = "") -> str:
    return "compatible" if ok else f"blocked:{reason}"


def _split_joined_labels(value: str) -> list[str]:
    return [part.strip() for part in str(value).split(";") if part.strip()]


def _hstouffer_source_unit_status(study_ids: list[str], hstouffer_by_study: dict[str, str]) -> str:
    if not study_ids:
        return "blocked:requires_" + "_".join(HSTOUFFER_REQUIRED_ORIGINAL_COLUMNS)
    statuses = [hstouffer_by_study.get(study_id, "not_compatible_header_detected") for study_id in study_ids]
    if all(status == "compatible_header_detected" for status in statuses):
        return "compatible"
    if any(status == "compatible_header_detected" for status in statuses):
        return "blocked:partial_hstouffer_original_columns"
    return "blocked:requires_" + "_".join(HSTOUFFER_REQUIRED_ORIGINAL_COLUMNS)


def _awmeta_source_unit_status(study_ids: list[str], awmeta_by_study: dict[str, dict[str, Any]]) -> str:
    if not study_ids:
        return "blocked:requires_original_variance_or_se"
    has_variance = [bool(awmeta_by_study.get(study_id, {}).get("variance_or_weight_columns")) for study_id in study_ids]
    if all(has_variance):
        return "compatible"
    if any(has_variance):
        return "blocked:partial_original_variance_or_se"
    return "blocked:requires_original_variance_or_se"


def build_matrix(harmonized: pd.DataFrame, *, corpus: str) -> pd.DataFrame:
    frame = harmonized.copy()
    frame["_source_unit"] = _source_unit_series(frame)
    hstouffer_report = hstouffer_materializer_feasibility(frame)
    hstouffer_by_study = {
        str(row.get("study_id")): str(row.get("audit_status"))
        for row in hstouffer_report.get("source_input_audit", [])
    }
    awmeta_by_study = {
        str(row.get("study_id")): row
        for row in awmeta_source_input_audit(frame)
    }

    rows: list[dict[str, Any]] = []
    for source_unit, group in frame.groupby("_source_unit", sort=True):
        study_ids = _join_unique(group["study_id"]) if "study_id" in group.columns else str(source_unit)
        source_study_ids = _split_joined_labels(study_ids)
        has_gene = _has_nonmissing(group, "gene_symbol")
        has_lfc = _has_nonmissing(group, "lfc")
        has_p = _has_nonmissing(group, "pvalue") or _has_nonmissing(group, "signed_z")
        has_padj = _has_nonmissing(group, "padj")
        has_rank = _has_nonmissing(group, "normalized_rank")
        has_n = _has_nonmissing(group, "n_ctrl") and _has_nonmissing(group, "n_treat")
        source_input_type = _join_unique(group["source_input_type"]) if "source_input_type" in group.columns else ""
        table_scope = _join_unique(group["table_scope"]) if "table_scope" in group.columns else ""
        pipeline = _join_unique(group["pipeline"]) if "pipeline" in group.columns else ""
        assay_type = _join_unique(group["assay_type"]) if "assay_type" in group.columns else ""
        source_path = _join_unique(group["source_path"]) if "source_path" in group.columns else ""
        hst_status = _hstouffer_source_unit_status(source_study_ids, hstouffer_by_study)
        awmeta_status = _awmeta_source_unit_status(source_study_ids, awmeta_by_study)
        raw_missing = "sample-level expression matrix; phenotype/class labels; per-study variance/SE"
        rows.append(
            {
                "corpus": corpus,
                "source_unit_id": str(source_unit),
                "study_id": study_ids,
                "pipeline": pipeline,
                "assay_type": assay_type,
                "source_input_type": source_input_type,
                "table_scope": table_scope,
                "source_path": source_path,
                "n_rows": int(len(group)),
                "has_gene_symbol": has_gene,
                "has_lfc": has_lfc,
                "has_pvalue": has_p,
                "has_padj": has_padj,
                "has_rank": has_rank,
                "has_sample_size": has_n,
                "awfisher": _status(has_gene and has_p, "needs_gene_and_pvalue"),
                "metarnaseq_fisher": _status(has_gene and has_p, "needs_gene_and_pvalue"),
                "robustrankaggreg": _status(has_gene and has_rank, "needs_gene_and_rank"),
                "metavolcanor": _status(has_gene and has_lfc and has_p, "needs_gene_lfc_pvalue"),
                "rankprod_exact": "blocked:requires_replicate_expression_and_class_origin_labels",
                "metade_pvalue_modes": _status(has_gene and has_p, "needs_gene_and_pvalue"),
                "metade_effect_size_modes": "blocked:requires_effect_variance_or_raw_expression",
                "dexma": "blocked:requires_expression_matrix_and_phenotype_metadata",
                "metaintegrator": "blocked:requires_expression_and_phenotype_objects",
                "omicc": "blocked:requires_expression_matrix_and_sample_group_annotations",
                "imageo": "blocked:requires_geo_expression_workflow_and_sample_labels",
                "networkanalyst": "blocked:requires_expression_tables_and_sample_labels",
                "crossmeta": "blocked:requires_geo_expression_matrices_and_platform_annotation",
                "deet": "not_comparator:uniform_precomputed_signature_atlas",
                "creeds": "not_comparator:curated_expression_signature_database",
                "generic_pvalue_combiners": "covered:fisher_stouffer_awfisher_metarnaseq",
                "hstouffer": hst_status,
                "awmeta": awmeta_status,
                "missing_for_raw_expression_tools": raw_missing,
                "notes": "summary-compatible tools can run from this DEG row set; raw-expression/effect-size tools require additional acquisition",
            }
        )
    return pd.DataFrame.from_records(rows, columns=FEASIBILITY_COLUMNS)


def write_markdown(frame: pd.DataFrame, path: Path) -> None:
    display = frame.loc[
        :,
        [
            "corpus",
            "source_unit_id",
            "study_id",
            "source_input_type",
            "awfisher",
            "robustrankaggreg",
            "metavolcanor",
            "rankprod_exact",
            "dexma",
            "hstouffer",
            "awmeta",
        ],
    ]
    lines = [
        "| " + " | ".join(display.columns) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    for row in display.itertuples(index=False):
        values = [str(value).replace("|", "\\|").replace("\n", " ") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    path.write_text("\n".join(lines) + "\n")


def write_outputs(*, harmonized_path: Path, corpus: str, output_csv: Path, output_md: Path, command: str) -> dict[str, Any]:
    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    matrix = build_matrix(harmonized, corpus=corpus)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    write_markdown(matrix, output_md)
    for artifact in [output_csv, output_md]:
        write_source_sidecar(artifact, command, inputs=[harmonized_path], metadata={"generator": "study-tool-feasibility", "corpus": corpus})
    return {"corpus": corpus, "rows": int(len(matrix)), "output_csv": str(output_csv), "output_md": str(output_md)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args(argv)
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_study_tool_feasibility_matrix.py",
            "--harmonized",
            args.harmonized,
            "--corpus",
            args.corpus,
            "--output-csv",
            args.output_csv,
            "--output-md",
            args.output_md,
        ]
    )
    summary = write_outputs(
        harmonized_path=args.harmonized,
        corpus=args.corpus,
        output_csv=args.output_csv,
        output_md=args.output_md,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
