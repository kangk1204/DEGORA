#!/usr/bin/env python
"""Export manuscript-facing DEGORA score resource tables and API evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.provenance import shell_command
from degora.resource_package import export_score_resource_package


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-csv", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--api-url")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = shell_command(
        [
            "make",
            "-C",
            "outputs/code",
            "resource-package",
            f"RESOURCE_SCORE_CSV={args.score_csv.resolve()}",
            f"RESOURCE_SCORE_DB={args.db.resolve()}",
            f"RESOURCE_OUTDIR={args.output_dir.resolve()}",
            f"RESOURCE_TOP_N={args.top_n}",
            f"RESOURCE_API_URL={args.api_url or ''}",
        ]
    )
    summary = export_score_resource_package(
        args.score_csv,
        args.db,
        args.output_dir,
        top_n=args.top_n,
        api_url=args.api_url,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
