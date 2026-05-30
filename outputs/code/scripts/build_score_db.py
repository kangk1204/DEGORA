#!/usr/bin/env python
"""Build DEGORA score CSV/metadata and a local SQLite evidence database."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.score_db import write_score_database
from degora.provenance import shell_command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--min-studies", type=int, default=2)
    parser.add_argument(
        "--extra-metadata",
        type=Path,
        help="Optional JSON object merged into degora_score_metadata.json and the SQLite meta table.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    extra_metadata = None
    if args.extra_metadata:
        extra_metadata = json.loads(args.extra_metadata.read_text())
        if not isinstance(extra_metadata, dict):
            raise ValueError("--extra-metadata must contain a JSON object")
    command_args = [
        "make",
        "-C",
        "outputs/code",
        "score-db",
        f"SCORE_HARMONIZED={args.harmonized.resolve()}",
        f"CATALOG={args.catalog.resolve() if args.catalog else ''}",
        f"OUTDIR={args.output_dir.resolve()}",
        f"SCORE_DB={args.db.resolve()}",
        f"SCORE_MIN_STUDIES={args.min_studies}",
    ]
    if args.extra_metadata:
        command_args.append(f"EXTRA_METADATA={args.extra_metadata.resolve()}")
    command = shell_command(command_args)
    summary = write_score_database(
        args.harmonized,
        args.output_dir,
        catalog_path=args.catalog,
        db_path=args.db,
        min_studies=args.min_studies,
        command=command,
        extra_metadata=extra_metadata,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
