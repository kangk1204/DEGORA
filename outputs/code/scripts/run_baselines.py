"""Run uniform-schema baseline adapters for a harmonized corpus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.baselines import write_baseline_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--min-studies", type=int, default=2)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = write_baseline_outputs(args.harmonized, args.output_dir, corpus=args.corpus, min_studies=args.min_studies)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
