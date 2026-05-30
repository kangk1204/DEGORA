from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.diagnostics import write_diagnostics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Write iteration diagnostic tables.")
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-studies", type=int, default=2)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = write_diagnostics(
        catalog_path=args.catalog,
        harmonized_path=args.harmonized,
        output_dir=args.output_dir,
        min_studies=args.min_studies,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
