#!/usr/bin/env python
"""Check generated outputs have reviewer-facing source sidecars."""

from __future__ import annotations

import argparse
from pathlib import Path

GENERATED_SUFFIXES = {".csv", ".json", ".parquet"}


def _is_generated_artifact(path: Path) -> bool:
    if not path.is_file():
        return False
    if path.name.endswith(".source") or ".provenance." in path.name:
        return False
    return path.suffix in GENERATED_SUFFIXES


def check_sources(artifact_dir: Path, *, command_prefix: str = "make -C outputs/code") -> list[str]:
    """Return source-sidecar validation errors for generated artifacts."""

    errors: list[str] = []
    for artifact in sorted(artifact_dir.rglob("*")):
        if not _is_generated_artifact(artifact):
            continue
        source = artifact.with_suffix(artifact.suffix + ".source")
        if not source.exists():
            errors.append(f"missing source sidecar: {source}")
            continue
        source_text = source.read_text().strip()
        first_line = source_text.splitlines()[0].strip() if source_text else ""
        if not first_line:
            errors.append(f"empty source sidecar: {source}")
        elif not first_line.startswith(command_prefix):
            errors.append(f"unexpected source command in {source}: {first_line}")
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        action="append",
        dest="artifact_dirs",
        default=[],
        help="Generated artifact directory to scan; may be passed multiple times.",
    )
    parser.add_argument("--results-dir", type=Path, help="Backward-compatible alias for --artifact-dir")
    parser.add_argument("--command-prefix", default="make -C outputs/code")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    artifact_dirs = list(args.artifact_dirs)
    if args.results_dir is not None:
        artifact_dirs.append(args.results_dir)
    if not artifact_dirs:
        artifact_dirs = [Path("../../outputs/results"), Path("../../data/deg/harmonized")]

    errors: list[str] = []
    for artifact_dir in artifact_dirs:
        errors.extend(check_sources(artifact_dir, command_prefix=args.command_prefix))
    if errors:
        print("Generated source sidecar check failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    joined = ", ".join(str(path) for path in artifact_dirs)
    print(f"Generated source sidecar check passed: {joined}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
