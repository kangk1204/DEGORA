#!/usr/bin/env python
"""Write conservative no-go evidence for hStouffer and AWmeta whole-corpus use."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.baselines import awmeta_deg_table_feasibility, hstouffer_materializer_feasibility
from degora.provenance import write_source_sidecar


ID_LIKE_COLUMNS = {"id", "gene", "gene_id", "gene_name", "symbol", "external_gene_name", "genesymbol"}
HSTOUFFER_NON_ID_COLUMNS = {"basemean", "log2foldchange", "lfcse", "stat", "pvalue", "padj"}


def _json_cell(value: Any) -> str:
    if isinstance(value, (list, dict, tuple)):
        return json.dumps(value, sort_keys=True)
    if value is None:
        return ""
    return str(value)


def _norm(value: object) -> str:
    return str(value).strip().lower().replace(" ", "_").replace("-", "_").replace(".", "_")


def _detected_columns(record: dict[str, Any]) -> list[str]:
    columns: list[str] = []
    for sheet in record.get("detected_sheets", []) or []:
        columns.extend(str(column) for column in sheet.get("columns", []) or [])
    if record.get("columns_inspected"):
        columns.extend(str(column) for column in record.get("columns_inspected", []) or [])
    return columns


def _hstouffer_subset_candidate(record: dict[str, Any]) -> bool:
    if record.get("audit_status") == "compatible_header_detected":
        return True
    pipeline = str(record.get("pipeline", "")).lower()
    if pipeline != "deseq2":
        return False
    normalized = {_norm(column) for column in _detected_columns(record)}
    has_id = bool(normalized.intersection(ID_LIKE_COLUMNS))
    return has_id and HSTOUFFER_NON_ID_COLUMNS.issubset(normalized)


def _hstouffer_rows(report: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for record in report.get("source_input_audit", []):
        candidate = _hstouffer_subset_candidate(record)
        missing = record.get("required_columns_missing", [])
        present = record.get("required_columns_present", [])
        if candidate and record.get("audit_status") != "compatible_header_detected":
            present = sorted(set(present).union({"id_like_column"}))
            missing = [field for field in missing if field != "id"]
        rows.append(
            {
                "study_id": record.get("study_id", ""),
                "paper_id": record.get("paper_id", ""),
                "pipeline": record.get("pipeline", ""),
                "source_path": record.get("source_path", ""),
                "required_fields_present": _json_cell(present),
                "required_fields_missing": _json_cell(missing),
                "compatible_subset_candidate": bool(candidate),
                "whole_corpus_status": "no_go",
                "blocker_reason": (
                    "source table appears compatible only as a tiny subset; full corpus lacks faithful hStouffer original inputs"
                    if candidate
                    else _json_cell(record.get("audit_status", "blocked"))
                ),
            }
        )
    return rows


def _awmeta_rows(report: dict[str, Any], harmonized: pd.DataFrame) -> list[dict[str, Any]]:
    paper_lookup = {}
    for row in harmonized.loc[:, [column for column in ["study_id", "paper_id"] if column in harmonized.columns]].drop_duplicates().itertuples(index=False):
        row_dict = row._asdict()
        paper_lookup[str(row_dict.get("study_id", ""))] = row_dict.get("paper_id", "")

    rows = []
    for record in report.get("source_input_audit", []):
        variance = record.get("variance_or_weight_columns", []) or []
        effect = record.get("effect_columns", []) or []
        rows.append(
            {
                "study_id": record.get("study_id", ""),
                "paper_id": paper_lookup.get(str(record.get("study_id", "")), ""),
                "pipeline": record.get("pipeline", ""),
                "source_path": record.get("source_path", ""),
                "required_fields_present": _json_cell({"effect_columns": effect, "variance_or_weight_columns": variance}),
                "required_fields_missing": "[]" if variance else '["variance_or_standard_error_or_weight"]',
                "compatible_subset_candidate": bool(variance),
                "whole_corpus_status": "no_go",
                "blocker_reason": (
                    "variance/SE candidate present only for subset; full corpus lacks original AWmeta variance/weight inputs"
                    if variance
                    else "missing original variance/SE or documented AWmeta-equivalent weights"
                ),
            }
        )
    return rows


def _write_report(
    path: Path,
    *,
    hstouffer_summary: dict[str, Any],
    awmeta_summary: dict[str, Any],
    hstouffer_rows: list[dict[str, Any]],
    awmeta_rows: list[dict[str, Any]],
) -> None:
    h_candidates = [row["study_id"] for row in hstouffer_rows if row["compatible_subset_candidate"]]
    aw_candidates = [row["study_id"] for row in awmeta_rows if row["compatible_subset_candidate"]]
    text = f"""# Direct Prior-Art No-Go Report

## Decision
hStouffer and AWmeta/AW-REM remain whole-corpus no-go comparators for this DEGORA corpus under faithful input requirements.

## hStouffer
- Whole-corpus status: no-go under faithful original DESeq2-like input requirements.
- Required original fields: id/baseMean/log2FoldChange/lfcSE/stat/pvalue/padj.
- Compatible-subset candidates: {", ".join(h_candidates) if h_candidates else "none"}.
- Subset interpretation: compatible-subset candidates are underpowered/no-claim and cannot support a superiority claim.
- Blocker id: {hstouffer_summary.get("blocker_id", "hstouffer_deg_table_materializer_blocked")}

## AWmeta / AW-REM
- Whole-corpus status: no-go under original per-gene variance/standard-error or documented equivalent-weight requirements.
- Compatible-subset candidates: {", ".join(aw_candidates) if aw_candidates else "none"}.
- Subset interpretation: compatible-subset candidates are underpowered/no-claim and cannot support a superiority claim.
- Blocker id: {awmeta_summary.get("blocker_id", "awmeta_variance_inputs_missing")}

## Claim Boundary
No superiority claim is allowed from these no-go artifacts. They document why faithful direct-prior-art execution is blocked on the whole corpus and what source-level acquisition would be needed to unblock it.
"""
    path.write_text(text)


def write_no_go_artifacts(harmonized_path: Path, output_dir: Path, *, corpus: str, command: str) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    hstouffer = hstouffer_materializer_feasibility(harmonized)
    awmeta = awmeta_deg_table_feasibility(harmonized)

    h_rows = _hstouffer_rows(hstouffer)
    aw_rows = _awmeta_rows(awmeta, harmonized)
    h_path = output_dir / "hstouffer_faithful_input_audit.tsv"
    aw_path = output_dir / "awmeta_variance_input_audit.tsv"
    pd.DataFrame(h_rows).to_csv(h_path, sep="\t", index=False)
    pd.DataFrame(aw_rows).to_csv(aw_path, sep="\t", index=False)

    summary = {
        "corpus": corpus,
        "claim_allowed": False,
        "direct_prior_art_superiority_claim_allowed": False,
        "hstouffer": {
            "whole_corpus_status": "no_go",
            "blocker_id": hstouffer.get("blocker_id", "hstouffer_deg_table_materializer_blocked"),
            "compatible_subset_candidates": [row["study_id"] for row in h_rows if row["compatible_subset_candidate"]],
            "subset_claim_status": "underpowered_no_claim",
            "n_audited_studies": len(h_rows),
            "decision": hstouffer.get("decision", ""),
        },
        "awmeta": {
            "whole_corpus_status": "no_go",
            "blocker_id": awmeta.get("blocker_id", "awmeta_variance_inputs_missing"),
            "compatible_subset_candidates": [row["study_id"] for row in aw_rows if row["compatible_subset_candidate"]],
            "subset_claim_status": "underpowered_no_claim",
            "n_audited_studies": len(aw_rows),
            "prohibited_approximations": awmeta.get("prohibited_approximations", []),
        },
    }
    summary_path = output_dir / "direct_prior_art_no_go_summary.json"
    report_path = output_dir / "direct_prior_art_no_go_report.md"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    _write_report(report_path, hstouffer_summary=hstouffer, awmeta_summary=awmeta, hstouffer_rows=h_rows, awmeta_rows=aw_rows)

    for artifact in (h_path, aw_path, summary_path, report_path):
        write_source_sidecar(
            artifact,
            command,
            inputs=[harmonized_path],
            metadata={"generator": "direct-prior-art-no-go", "corpus": corpus, "claim_allowed": False},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--corpus", default="hypoxia")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_direct_prior_art_no_go.py"
        f" --harmonized {args.harmonized}"
        f" --output-dir {args.output_dir}"
        f" --corpus {args.corpus}"
    )
    summary = write_no_go_artifacts(args.harmonized, args.output_dir, corpus=args.corpus, command=command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
