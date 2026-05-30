#!/usr/bin/env python
"""Audit generated result artifacts for `.source` provenance sidecars."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

GENERATED_SUFFIXES = {".csv", ".db", ".json", ".md", ".parquet", ".tsv"}
AUDIT_FILENAME = "provenance_audit.json"
DEFAULT_ITERATIONS = ("iter-9", "iter-10", "iter-11", "iter-12", "iter-13")


def source_sidecar(path: Path) -> Path:
    """Return the expected `.source` sibling for a generated artifact."""

    return path.with_suffix(path.suffix + ".source")


def is_generated_artifact(path: Path) -> bool:
    """Return whether `path` is a generated artifact covered by the audit."""

    if not path.is_file():
        return False
    if path.name.endswith(".source") or ".provenance." in path.name:
        return False
    if path.name == AUDIT_FILENAME:
        return False
    return path.suffix in GENERATED_SUFFIXES


def audit_iteration(iteration_dir: Path) -> dict[str, Any]:
    """Audit one result iteration directory for missing or empty `.source` sidecars."""

    artifacts: list[dict[str, Any]] = []
    missing: list[str] = []
    empty: list[str] = []
    if not iteration_dir.exists():
        return {
            "iteration_dir": str(iteration_dir),
            "exists": False,
            "artifact_count": 0,
            "missing_source": [],
            "empty_source": [],
            "artifacts": [],
            "passed": True,
        }

    for artifact in sorted(iteration_dir.rglob("*")):
        if not is_generated_artifact(artifact):
            continue
        sidecar = source_sidecar(artifact)
        record: dict[str, Any] = {
            "artifact": str(artifact),
            "source": str(sidecar),
            "source_exists": sidecar.exists(),
            "source_first_line": "",
        }
        if not sidecar.exists():
            missing.append(str(artifact))
        else:
            text = sidecar.read_text().strip()
            first_line = text.splitlines()[0].strip() if text else ""
            record["source_first_line"] = first_line
            if not first_line:
                empty.append(str(sidecar))
        artifacts.append(record)

    return {
        "iteration_dir": str(iteration_dir),
        "exists": True,
        "artifact_count": len(artifacts),
        "missing_source": missing,
        "empty_source": empty,
        "artifacts": artifacts,
        "passed": not missing and not empty,
    }


def build_audit(results_root: Path, iterations: list[str]) -> dict[str, Any]:
    """Build a deterministic provenance audit record for requested iterations."""

    iteration_records = [audit_iteration(results_root / iteration) for iteration in iterations]
    return {
        "schema_version": "1.0",
        "results_root": str(results_root),
        "iterations": iterations,
        "iteration_records": iteration_records,
        "passed": all(record["passed"] for record in iteration_records),
    }


def write_audit(audit: dict[str, Any], output_path: Path, command: str) -> None:
    """Write the audit JSON and reviewer-facing `.source` command sidecar."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    source_sidecar(output_path).write_text(command + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=Path("../../outputs/results"))
    parser.add_argument(
        "--iteration",
        action="append",
        dest="iterations",
        default=[],
        help="Iteration directory name to audit, e.g. iter-9. May be repeated.",
    )
    parser.add_argument("--output", type=Path, default=Path("../../outputs/results/iter-13/provenance_audit.json"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    iterations = args.iterations or list(DEFAULT_ITERATIONS)
    command = (
        "make -C outputs/code provenance-check"
        f" RESULTS_ROOT={args.results_root.resolve()}"
        f" OUTPUT={args.output.resolve()}"
        f" ITERATIONS={','.join(iterations)}"
    )
    audit = build_audit(args.results_root, iterations)
    write_audit(audit, args.output, command)
    print(json.dumps({"passed": audit["passed"], "output": str(args.output), "iterations": iterations}, sort_keys=True))
    return 0 if audit["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
