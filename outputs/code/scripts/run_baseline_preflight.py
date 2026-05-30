"""Write only R preflight and baseline blocker ledgers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.baselines import failure_ledger, r_preflight_report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--corpus", default="hypoxia")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    command = f"make -C outputs/code baseline-preflight CORPUS={args.corpus} BASELINE_OUTDIR={args.output_dir.resolve()}"
    preflight = r_preflight_report()
    preflight_path = args.output_dir / "r_preflight_report.json"
    preflight_path.write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    preflight_path.with_suffix(preflight_path.suffix + ".source").write_text(command + "\n")
    ledger = failure_ledger(args.corpus, preflight)
    ledger_path = args.output_dir / "baseline_failure_ledger.csv"
    ledger.to_csv(ledger_path, index=False)
    ledger_path.with_suffix(ledger_path.suffix + ".source").write_text(command + "\n")
    print(json.dumps({"output_dir": str(args.output_dir.resolve()), "open_s1_blockers": len(ledger)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
