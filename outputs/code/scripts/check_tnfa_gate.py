from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.tnfa_gate import evaluate_tnfa_gate_csv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check TNF-alpha active-row/source-family gate.")
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--scarcity-report", type=Path)
    parser.add_argument("--output", type=Path, help="Optional JSON path for the gate result.")
    parser.add_argument("--scarcity-candidate-threshold", type=int, default=20)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = evaluate_tnfa_gate_csv(
        args.ledger,
        scarcity_candidate_threshold=args.scarcity_candidate_threshold,
    )
    result_json = json.dumps(result.to_dict(), indent=2, sort_keys=True)
    print(result_json)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result_json + "\n")
    if result.passed:
        return 0
    if result.scarcity_triggered and args.scarcity_report and args.scarcity_report.exists():
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
