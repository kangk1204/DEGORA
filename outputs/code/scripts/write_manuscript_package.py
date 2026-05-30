#!/usr/bin/env python
"""Validate and index the current manuscript-facing paper package."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from degora.provenance import shell_command, write_source_sidecar


DEFAULT_MANUSCRIPT = Path("outputs/manuscript/main.md")
DEFAULT_TABLES = Path("outputs/tables/manuscript")
DEFAULT_FIGURES = Path("outputs/figures/manuscript")
DEFAULT_OUTPUT = Path("outputs/manuscript")

TABLE_INPUTS = [
    "Figure_Index.csv",
    "Table1_core_benchmark_summary.csv",
    "Table2_cross_platform_summary.csv",
    "Supplementary_Table_Index.csv",
    "degora_manuscript_tables_manifest.json",
]

FIGURE_MANIFESTS = [
    "FIGURE_1_SCHEME/figure1_scheme_manifest.json",
    "FIGURE_2_BENCHMARK/figure2_benchmark_manifest.json",
    "DEGORA_CROSS_PLATFORM/degora_cross_platform_manifest.json",
    "DEGORA_ATLAS/degora_atlas_manifest.json",
    "DEGORA_EVIDENCE_EXPLORER/degora_evidence_explorer_manifest.json",
]

HYPOXIA_STALE_PATTERN = re.compile(
    r"(hypoxia[^\n.]{0,120}recall[^\n.]{0,60}0\.80"
    r"|hypoxia[^\n.]{0,60}(?:=|:|of|was|is|are|reached|achieved|stayed|remained|equal[s]?)\s*0\.80"
    r"|recall[^\n.]{0,20}0\.80[^\n.]{0,20}hypoxia)",
    re.IGNORECASE,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_record(path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {"path": str(path), "exists": path.exists()}
    if path.is_file():
        record.update({"size_bytes": path.stat().st_size, "sha256": _sha256(path)})
    return record


def _require_files(paths: list[Path], failures: list[str]) -> None:
    for path in paths:
        if not path.exists():
            failures.append(f"missing required paper-package artifact: {path}")


def _format_recall(value: object) -> str:
    return f"{float(value):.2f}"


def _check_manuscript_text(text: str, failures: list[str], warnings: list[str]) -> None:
    if "[REFERENCE NEEDED" in text:
        failures.append("manuscript still contains [REFERENCE NEEDED] placeholders")
    if HYPOXIA_STALE_PATTERN.search(text):
        failures.append("manuscript contains a stale hypoxia 0.80-style claim")
    if "Citation Audit Still Needed" in text:
        warnings.append("dataset-level citation audit is still marked as pending")
    if "Data and Code Availability" not in text:
        failures.append("manuscript lacks a Data and Code Availability section")


def _check_table_claims(table1: pd.DataFrame, table2: pd.DataFrame, text: str, failures: list[str]) -> dict[str, Any]:
    required_topics = {"IFN", "ER stress", "Heat shock", "Hypoxia"}
    observed_topics = set(table1["topic"].astype(str))
    missing = required_topics.difference(observed_topics)
    if missing:
        failures.append(f"Table1 missing benchmark topics: {sorted(missing)}")

    headline = []
    for row in table1.itertuples(index=False):
        value = _format_recall(row.degora_recall_at_100)
        headline.append({"topic": row.topic, "degora_recall_at_100": value})
        if value not in text:
            failures.append(f"manuscript does not contain Table1 recall@100 value {value} for {row.topic}")

    table2_key = table2[["topic", "mode", "recall_at_20", "recall_at_50", "recall_at_100"]].copy()
    for column in ["recall_at_20", "recall_at_50", "recall_at_100"]:
        table2_key[column] = table2_key[column].map(_format_recall)
    return {
        "table1_headline_recall_at_100": headline,
        "table2_rows": table2_key.to_dict(orient="records"),
    }


def build_package(
    *,
    manuscript: Path,
    tables_dir: Path,
    figures_dir: Path,
    output_dir: Path,
    command: str,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    warnings: list[str] = []

    table_paths = [tables_dir / name for name in TABLE_INPUTS]
    figure_paths = [figures_dir / name for name in FIGURE_MANIFESTS]
    required = [manuscript, *table_paths, *figure_paths]
    _require_files(required, failures)

    text = manuscript.read_text() if manuscript.exists() else ""
    _check_manuscript_text(text, failures, warnings)

    table_claims: dict[str, Any] = {}
    table1_path = tables_dir / "Table1_core_benchmark_summary.csv"
    table2_path = tables_dir / "Table2_cross_platform_summary.csv"
    if table1_path.exists() and table2_path.exists():
        table_claims = _check_table_claims(pd.read_csv(table1_path), pd.read_csv(table2_path), text, failures)

    manifest_path = output_dir / "degora_paper_manifest.json"
    validation_path = output_dir / "degora_paper_validation.txt"
    index_path = output_dir / "degora_paper_artifact_index.csv"
    artifacts = [_artifact_record(path) for path in required]
    index = pd.DataFrame(artifacts)
    index.to_csv(index_path, index=False)

    validation_lines = [
        "DEGORA paper package validation",
        f"generated_at={datetime.now(UTC).isoformat()}",
        f"failures={len(failures)}",
        f"warnings={len(warnings)}",
        "",
        "Failures:",
        *[f"- {failure}" for failure in failures],
        "",
        "Warnings:",
        *[f"- {warning}" for warning in warnings],
    ]
    validation_path.write_text("\n".join(validation_lines) + "\n")

    manifest = {
        "generated_at": datetime.now(UTC).isoformat(),
        "script": "outputs/code/scripts/write_manuscript_package.py",
        "command": command,
        "manuscript": str(manuscript),
        "artifacts": artifacts,
        "table_claims": table_claims,
        "failures": failures,
        "warnings": warnings,
        "pre_submission_tasks": [
            "final source-study and GEO accession citation audit",
            "public repository release tag",
            "archive DOI for code, catalog, harmonized tables, figures, and source data",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    for artifact in [manifest_path, validation_path, index_path]:
        write_source_sidecar(
            artifact,
            command,
            inputs=required,
            metadata={"generator": "manuscript-package"},
            write_json=False,
        )

    if failures:
        raise ValueError("; ".join(failures))

    return {
        "manifest": str(manifest_path),
        "validation": str(validation_path),
        "artifact_index": str(index_path),
        "n_artifacts": len(required),
        "warnings": warnings,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=DEFAULT_MANUSCRIPT)
    parser.add_argument("--tables-dir", type=Path, default=DEFAULT_TABLES)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_manuscript_package.py",
            "--manuscript",
            args.manuscript,
            "--tables-dir",
            args.tables_dir,
            "--figures-dir",
            args.figures_dir,
            "--output-dir",
            args.output_dir,
        ]
    )
    print(json.dumps(build_package(manuscript=args.manuscript, tables_dir=args.tables_dir, figures_dir=args.figures_dir, output_dir=args.output_dir, command=command), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
