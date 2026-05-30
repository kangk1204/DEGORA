#!/usr/bin/env python
"""Refresh public-file comparator and prior-art feasibility tables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.baselines import (
    BASELINE_RESULT_COLUMNS,
    BASELINE_MANIFEST_COLUMNS,
    FAILURE_LEDGER_COLUMNS,
    PUBLIC_SUMMARY_TOOL_INPUT_COLUMNS,
    _write_markdown_table,
    baseline_result_paths,
    public_summary_tool_input_requirements,
)
from degora.provenance import shell_command, write_source_sidecar


DEFAULT_CORPORA = [
    "ifn-pilot",
    "er-stress-cross-platform",
    "heat-shock-benchmark",
    "hypoxia-hif1-benchmark",
    "ifn-cross-platform",
    "hypoxia-cross-platform",
]


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text()) if path.exists() else {}


def _load_ledger(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame(columns=FAILURE_LEDGER_COLUMNS)


def _load_outputs(baseline_dir: Path) -> dict[str, pd.DataFrame]:
    outputs: dict[str, pd.DataFrame] = {}
    required = set(BASELINE_RESULT_COLUMNS)
    for path in baseline_result_paths(baseline_dir):
        frame = pd.read_csv(path, sep="\t")
        if frame.empty or not required.issubset(frame.columns):
            continue
        outputs[str(frame["method_id"].iloc[0])] = frame
    return outputs


def _refresh_manifest(manifest_path: Path, csv_path: Path, md_path: Path, command: str, n_rows: int) -> None:
    if not manifest_path.exists():
        return
    manifest = pd.read_csv(manifest_path)
    if not set(BASELINE_MANIFEST_COLUMNS).issubset(manifest.columns):
        return
    updates = {
        "public_summary_tool_input_requirements": (csv_path, n_rows),
        "public_summary_tool_input_requirements_markdown": (md_path, n_rows),
    }
    for artifact_type, (artifact_path, rows) in updates.items():
        mask = manifest["artifact_type"].astype(str).eq(artifact_type)
        if mask.any():
            manifest.loc[mask, "artifact"] = str(artifact_path)
            manifest.loc[mask, "source_command"] = command
            manifest.loc[mask, "rows"] = int(rows)
            manifest.loc[mask, "notes"] = "supplementary reviewer table: why each comparator/resource can or cannot run from public supplementary DEG files"
    manifest.to_csv(manifest_path, index=False)
    write_source_sidecar(
        manifest_path,
        command,
        inputs=[csv_path, md_path],
        metadata={"generator": "public-summary-tool-requirements", "artifact": "baseline_manifest"},
        write_json=False,
    )


def refresh_corpus(corpus: str, *, results_root: Path, harmonized_root: Path, command: str) -> pd.DataFrame:
    baseline_dir = results_root / corpus / "baselines"
    harmonized_path = harmonized_root / f"{corpus}_harmonized.csv"
    if not harmonized_path.exists():
        raise FileNotFoundError(f"missing harmonized input for {corpus}: {harmonized_path}")
    if not baseline_dir.exists():
        raise FileNotFoundError(f"missing baseline directory for {corpus}: {baseline_dir}")

    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    table = public_summary_tool_input_requirements(
        corpus=corpus,
        harmonized=harmonized,
        outputs=_load_outputs(baseline_dir),
        preflight=_load_json(baseline_dir / "r_preflight_report.json"),
        ledger=_load_ledger(baseline_dir / "baseline_failure_ledger.csv"),
        hstouffer_report=_load_json(baseline_dir / "hstouffer_feasibility_report.json"),
        awmeta_report=_load_json(baseline_dir / "awmeta_feasibility_report.json"),
    )
    csv_path = baseline_dir / "public_summary_tool_input_requirements.csv"
    md_path = baseline_dir / "public_summary_tool_input_requirements.md"
    table.to_csv(csv_path, index=False, quoting=csv.QUOTE_MINIMAL)
    _write_markdown_table(table, md_path)
    inputs = [
        harmonized_path,
        baseline_dir / "baseline_manifest.csv",
        baseline_dir / "baseline_failure_ledger.csv",
        baseline_dir / "r_preflight_report.json",
        baseline_dir / "hstouffer_feasibility_report.json",
        baseline_dir / "awmeta_feasibility_report.json",
    ]
    for artifact in [csv_path, md_path]:
        write_source_sidecar(
            artifact,
            command,
            inputs=inputs,
            metadata={"generator": "public-summary-tool-requirements", "corpus": corpus},
            write_json=False,
        )
    _refresh_manifest(baseline_dir / "baseline_manifest.csv", csv_path, md_path, command, len(table))
    return table.assign(corpus=corpus)


def write_combined(frame: pd.DataFrame, output_csv: Path, output_md: Path, command: str, inputs: list[Path]) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    _write_markdown_table(frame, output_md)
    for artifact in [output_csv, output_md]:
        write_source_sidecar(
            artifact,
            command,
            inputs=inputs,
            metadata={"generator": "public-summary-tool-requirements", "scope": "combined"},
            write_json=False,
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", action="append", dest="corpora", help="Corpus id to refresh; repeatable.")
    parser.add_argument("--results-root", type=Path, default=Path("outputs/results"))
    parser.add_argument("--harmonized-root", type=Path, default=Path("data/deg/harmonized"))
    parser.add_argument("--output-csv", type=Path, default=Path("outputs/results/comparator_public_summary_input_requirements.csv"))
    parser.add_argument("--output-md", type=Path, default=Path("outputs/results/comparator_public_summary_input_requirements.md"))
    args = parser.parse_args(argv)

    corpora = args.corpora or DEFAULT_CORPORA
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_public_summary_tool_requirements.py",
            *sum([["--corpus", corpus] for corpus in corpora], []),
            "--results-root",
            args.results_root,
            "--harmonized-root",
            args.harmonized_root,
            "--output-csv",
            args.output_csv,
            "--output-md",
            args.output_md,
        ]
    )
    frames = [
        refresh_corpus(corpus, results_root=args.results_root, harmonized_root=args.harmonized_root, command=command)
        for corpus in corpora
    ]
    combined = pd.concat(frames, ignore_index=True)
    inputs = [args.results_root / corpus / "baselines" / "public_summary_tool_input_requirements.csv" for corpus in corpora]
    write_combined(combined, args.output_csv, args.output_md, command, inputs)
    print(
        json.dumps(
            {
                "corpora": corpora,
                "rows": int(len(combined)),
                "output_csv": str(args.output_csv),
                "output_md": str(args.output_md),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
