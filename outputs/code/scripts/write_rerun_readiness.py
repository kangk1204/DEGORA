#!/usr/bin/env python
"""Write an iteration rerun readiness checklist without changing analysis inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

try:
    from scripts.check_source_provenance import source_sidecar
except ModuleNotFoundError:  # pragma: no cover - direct CLI execution path
    from check_source_provenance import source_sidecar


RERUN_STEPS = [
    {
        "step": "slice",
        "command": "make -C outputs/code slice ITER={iteration}",
        "outputs": [
            "outputs/results/iter-{iteration}/slice_metrics.json",
            "outputs/results/iter-{iteration}/slice_consensus.csv",
            "data/deg/harmonized/iter-{iteration}_harmonized.csv",
        ],
    },
    {
        "step": "diagnose",
        "command": "make -C outputs/code diagnose ITER={iteration} DIAG_SOURCE_ITER={iteration} DIAG_HARMONIZED=../../data/deg/harmonized/iter-{iteration}_harmonized.csv",
        "outputs": [
            "outputs/results/iter-{iteration}/diagnostic_summary.json",
            "outputs/results/iter-{iteration}/orientation_audit.csv",
            "outputs/results/iter-{iteration}/sensitivity_metrics.csv",
        ],
    },
    {
        "step": "rank-plane",
        "command": "make -C outputs/code rank-plane ITER={iteration} HARMONIZED=../../data/deg/harmonized/iter-{iteration}_harmonized.csv",
        "outputs": [
            "outputs/results/iter-{iteration}/rank_plane_summary.json",
            "outputs/results/iter-{iteration}/rank_plane_gene_summary.csv",
        ],
    },
    {
        "step": "baseline",
        "command": "make -C outputs/code baseline CORPUS=hypoxia HARMONIZED=../../data/deg/harmonized/iter-{iteration}_harmonized.csv BASELINE_OUTDIR=../../outputs/results/iter-{iteration}/baselines",
        "outputs": [
            "outputs/results/iter-{iteration}/baselines/baseline_parity_matrix.csv",
            "outputs/results/iter-{iteration}/baselines/baseline_failure_ledger.csv",
            "outputs/results/iter-{iteration}/baselines/baseline_manifest.csv",
        ],
    },
    {
        "step": "comparator-summary",
        "command": "PYTHONPATH=outputs/code python outputs/code/scripts/write_comparator_summary.py --baseline-dir outputs/results/iter-{iteration}/baselines --output-csv outputs/results/iter-{iteration}/comparator_recall_summary.csv --output-md outputs/results/iter-{iteration}/comparator_recall_summary.md",
        "outputs": [
            "outputs/results/iter-{iteration}/comparator_recall_summary.csv",
            "outputs/results/iter-{iteration}/comparator_recall_summary.md",
        ],
    },
    {
        "step": "publication-gate",
        "command": "make -C outputs/code publication-gate ITER={iteration}",
        "outputs": ["outputs/results/iter-{iteration}/publication_gate_summary.md"],
    },
    {
        "step": "check",
        "command": "make check",
        "outputs": [],
    },
]


LANE_INPUTS = [
    {
        "id": "catalog",
        "path": "data/studies/hypoxia_catalog.csv",
        "required_for": "slice",
        "note": "Active-row decisions must already satisfy strict S0 eligibility and source-family rules.",
    },
    {
        "id": "iter14_candidate_ledger",
        "path": "data/studies/curation/iter14_hypoxia_candidate_ledger.csv",
        "required_for": "Lane A audit trail before activating new rows",
        "note": "Missing is acceptable only if Lane A reports no eligible new S0 candidates.",
    },
    {
        "id": "harmonized_iter14",
        "path": "data/deg/harmonized/iter-{iteration}_harmonized.csv",
        "required_for": "diagnose, rank-plane, baseline, comparator-summary",
        "note": "Produced by the slice step after catalog/raw-input integration.",
    },
    {
        "id": "hstouffer_feasibility",
        "path": "outputs/results/iter-{iteration}/baselines/hstouffer_feasibility_report.json",
        "required_for": "faithful hStouffer publication gate",
        "note": "Must document original DESeq2-like baseMean/lfcSE/stat availability or a blocker; no imputation.",
    },
    {
        "id": "awmeta_feasibility",
        "path": "outputs/results/iter-{iteration}/baselines/awmeta_feasibility_report.json",
        "required_for": "faithful AWmeta/AW-Fisher/AW-REM publication gate",
        "note": "Must document per-study variance/SE inputs or a blocker; no p-value-derived SE surrogates.",
    },
]


def _format(value: str, iteration: int) -> str:
    return value.format(iteration=iteration)


def _rel_exists(project_root: Path, rel_path: str) -> bool:
    return (project_root / rel_path).exists()


def _sidecar_status(project_root: Path, rel_path: str) -> str:
    path = project_root / rel_path
    if not path.exists():
        return "not_applicable_missing_artifact"
    sidecar = source_sidecar(path)
    if not sidecar.exists():
        return "missing"
    if not sidecar.read_text().strip():
        return "empty"
    return "present"


def build_readiness(project_root: Path, iteration: int) -> dict[str, Any]:
    """Return deterministic readiness state for an iteration integration rerun."""

    lane_inputs = []
    blockers: list[str] = []
    warnings: list[str] = []
    for item in LANE_INPUTS:
        rel_path = _format(item["path"], iteration)
        exists = _rel_exists(project_root, rel_path)
        sidecar = _sidecar_status(project_root, rel_path)
        lane_inputs.append({**item, "path": rel_path, "exists": exists, "source_sidecar": sidecar})
        if item["id"] in {"catalog"} and not exists:
            blockers.append(f"missing required input: {rel_path}")
        elif item["id"] in {"harmonized_iter14"} and not exists:
            blockers.append(f"downstream reruns are not executable until slice creates {rel_path}")
        elif not exists:
            warnings.append(f"pending lane artifact: {rel_path} ({item['required_for']})")

    steps = []
    for step in RERUN_STEPS:
        command = _format(step["command"], iteration)
        outputs = []
        for output in step["outputs"]:
            rel_path = _format(output, iteration)
            outputs.append(
                {
                    "path": rel_path,
                    "exists": _rel_exists(project_root, rel_path),
                    "source_sidecar": _sidecar_status(project_root, rel_path),
                }
            )
        steps.append({"step": step["step"], "command": command, "outputs": outputs})

    constraints = [
        "Do not tune DEGORA scoring, source weights, min_studies, p-value clipping, or locked HIF targets.",
        "Do not activate raw-count, expression-only, top-K-only, replicate-insufficient, non-human/mouse, or source-family-inflating rows as S0 evidence.",
        "Do not derive hStouffer/AWmeta variance inputs from p-values or signed-z surrogates.",
        "Do not claim SOTA superiority unless all direct-prior-art comparators are faithfully runnable under the equal-tuning gate.",
        "Every generated artifact must have a sibling .source sidecar before manuscript-facing use.",
    ]

    ready_to_run_full_chain = not blockers
    return {
        "schema_version": "1.0",
        "iteration": iteration,
        "ready_to_run_full_chain": ready_to_run_full_chain,
        "blockers": blockers,
        "warnings": warnings,
        "lane_inputs": lane_inputs,
        "rerun_steps": steps,
        "constraints": constraints,
    }


def render_markdown(readiness: dict[str, Any]) -> str:
    iteration = readiness["iteration"]
    lines = [
        f"# Iteration {iteration} Integration Rerun Readiness",
        "",
        f"Ready to run full chain: **{readiness['ready_to_run_full_chain']}**",
        "",
        "## Blockers",
        "",
    ]
    blockers = readiness["blockers"]
    lines.extend(f"- {item}" for item in blockers) if blockers else lines.append("- None")
    lines.extend(["", "## Pending Warnings", ""])
    warnings = readiness["warnings"]
    lines.extend(f"- {item}" for item in warnings) if warnings else lines.append("- None")
    lines.extend(["", "## Non-negotiable Constraints", ""])
    lines.extend(f"- {item}" for item in readiness["constraints"])
    lines.extend(["", "## Lane Input Status", ""])
    lines.append("| Input | Exists | Source sidecar | Required for | Path |")
    lines.append("| --- | --- | --- | --- | --- |")
    for item in readiness["lane_inputs"]:
        lines.append(
            f"| {item['id']} | {item['exists']} | {item['source_sidecar']} | {item['required_for']} | `{item['path']}` |"
        )
    lines.extend(["", "## Rerun Command Sequence", ""])
    for idx, step in enumerate(readiness["rerun_steps"], start=1):
        lines.extend([f"{idx}. `{step['command']}`"])
        if step["outputs"]:
            for output in step["outputs"]:
                lines.append(
                    f"   - `{output['path']}` exists={output['exists']} source={output['source_sidecar']}"
                )
    lines.extend(
        [
            "",
            "## Publication Gate Rule",
            "",
            "Keep the manuscript conclusion conservative until strict S0 inputs and faithful comparator inputs have landed, rerun artifacts carry `.source` sidecars, and `make check` passes.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(readiness: dict[str, Any], output_dir: Path, command: str) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "rerun_readiness.json"
    md_path = output_dir / "rerun_readiness.md"
    json_path.write_text(json.dumps(readiness, indent=2, sort_keys=True) + "\n")
    md_path.write_text(render_markdown(readiness))
    source_sidecar(json_path).write_text(command + "\n")
    source_sidecar(md_path).write_text(command + "\n")
    return json_path, md_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("../.."))
    parser.add_argument("--iteration", type=int, required=True)
    parser.add_argument("--output-dir", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    project_root = args.project_root.resolve()
    output_dir = args.output_dir or project_root / "outputs" / "results" / f"iter-{args.iteration}"
    command = f"make -C outputs/code rerun-readiness ITER={args.iteration}"
    readiness = build_readiness(project_root, args.iteration)
    json_path, md_path = write_outputs(readiness, output_dir, command)
    print(json.dumps({"ready": readiness["ready_to_run_full_chain"], "json": str(json_path), "markdown": str(md_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
