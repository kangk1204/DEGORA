#!/usr/bin/env python
"""Write one-row-per-method prior-art coverage summary from the public-file table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.baselines import _write_markdown_table
from degora.provenance import shell_command, write_source_sidecar


SUMMARY_COLUMNS = [
    "method_id",
    "method_family",
    "public_summary_deg_status",
    "can_run_from_current_public_files",
    "manuscript_use",
    "source_or_package",
    "faithful_adapter_decision",
    "corpora_covered",
    "n_corpora",
    "status_values",
    "missing_or_nonfaithful_inputs",
    "notes",
]


def _join_unique(values: pd.Series, *, max_chars: int | None = None) -> str:
    labels = [str(value).strip() for value in values.dropna().tolist()]
    labels = [label for label in labels if label and label.lower() != "nan"]
    text = "; ".join(dict.fromkeys(labels))
    if max_chars is not None and len(text) > max_chars:
        return text[:max_chars].rstrip() + " ... [truncated; see corpus-level public summary table]"
    return text


def build_summary(frame: pd.DataFrame) -> pd.DataFrame:
    if "corpus" not in frame.columns:
        raise ValueError("input table must contain a corpus column")
    rows: list[dict[str, Any]] = []
    for method_id, group in frame.groupby("method_id", sort=False):
        rows.append(
            {
                "method_id": method_id,
                "method_family": _join_unique(group["method_family"]),
                "public_summary_deg_status": _join_unique(group["public_summary_deg_status"]),
                "can_run_from_current_public_files": _join_unique(group["can_run_from_current_public_files"]),
                "manuscript_use": _join_unique(group["manuscript_use"]),
                "source_or_package": _join_unique(group["source_or_package"]),
                "faithful_adapter_decision": _join_unique(group["faithful_adapter_decision"]),
                "corpora_covered": _join_unique(group["corpus"]),
                "n_corpora": int(group["corpus"].nunique()),
                "status_values": _join_unique(group["current_pipeline_status"]),
                "missing_or_nonfaithful_inputs": _join_unique(group["missing_or_nonfaithful_inputs"], max_chars=600),
                "notes": _join_unique(group["notes"], max_chars=600),
            }
        )
    return pd.DataFrame.from_records(rows, columns=SUMMARY_COLUMNS)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=Path("outputs/results/comparator_public_summary_input_requirements.csv"))
    parser.add_argument("--output-csv", type=Path, default=Path("outputs/results/prior_art_coverage_summary.csv"))
    parser.add_argument("--output-md", type=Path, default=Path("outputs/results/prior_art_coverage_summary.md"))
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.input_csv)
    summary = build_summary(frame)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    _write_markdown_table(summary, args.output_md)
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_prior_art_coverage_summary.py",
            "--input-csv",
            args.input_csv,
            "--output-csv",
            args.output_csv,
            "--output-md",
            args.output_md,
        ]
    )
    for artifact in [args.output_csv, args.output_md]:
        write_source_sidecar(
            artifact,
            command,
            inputs=[args.input_csv],
            metadata={"generator": "prior-art-coverage-summary"},
            write_json=False,
        )
    print(json.dumps({"rows": int(len(summary)), "output_csv": str(args.output_csv), "output_md": str(args.output_md)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
