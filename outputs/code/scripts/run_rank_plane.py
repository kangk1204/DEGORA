#!/usr/bin/env python
"""Run rank-plane diagnostics for a harmonized DEGORA table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.rank_plane import write_rank_plane_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--joint-threshold", type=float, default=0.9)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = write_rank_plane_outputs(
        args.harmonized,
        args.output_dir,
        joint_threshold=args.joint_threshold,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
