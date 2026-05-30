from __future__ import annotations

import json
from pathlib import Path

from scripts.check_source_provenance import audit_iteration, build_audit, is_generated_artifact, main, source_sidecar


def test_audit_iteration_reports_missing_source_sidecars(tmp_path: Path) -> None:
    iteration_dir = tmp_path / "iter-9"
    iteration_dir.mkdir()
    good = iteration_dir / "slice_metrics.json"
    missing = iteration_dir / "slice_consensus.csv"
    ignored_source = iteration_dir / "already.source"
    audit_output = iteration_dir / "provenance_audit.json"
    good.write_text("{}\n")
    source_sidecar(good).write_text("make -C outputs/code slice ITER=9\n")
    missing.write_text("gene_symbol\nVEGFA\n")
    ignored_source.write_text("not an artifact\n")
    audit_output.write_text("{}\n")

    audit = audit_iteration(iteration_dir)

    assert audit["passed"] is False
    assert audit["artifact_count"] == 2
    assert audit["missing_source"] == [str(missing)]
    assert is_generated_artifact(audit_output) is False


def test_main_writes_iter10_audit_json_and_source(tmp_path: Path) -> None:
    results_root = tmp_path / "outputs" / "results"
    iter9 = results_root / "iter-9"
    iter10 = results_root / "iter-10"
    iter9.mkdir(parents=True)
    iter10.mkdir()
    for artifact in [iter9 / "slice_metrics.json", iter10 / "tnfa_scarcity_report.md"]:
        artifact.write_text("{}\n")
        source_sidecar(artifact).write_text("make -C outputs/code generated\n")
    output = iter10 / "provenance_audit.json"

    exit_code = main([
        "--results-root",
        str(results_root),
        "--iteration",
        "iter-9",
        "--iteration",
        "iter-10",
        "--output",
        str(output),
    ])

    assert exit_code == 0
    audit = json.loads(output.read_text())
    assert audit["passed"] is True
    assert audit["iterations"] == ["iter-9", "iter-10"]
    assert source_sidecar(output).read_text().startswith("make -C outputs/code provenance-check")


def test_build_audit_treats_missing_future_iteration_as_non_failure(tmp_path: Path) -> None:
    audit = build_audit(tmp_path, ["iter-10"])

    assert audit["passed"] is True
    assert audit["iteration_records"][0]["exists"] is False


def test_tsv_outputs_are_generated_artifacts(tmp_path: Path) -> None:
    artifact = tmp_path / "hypoxia_fisher_default.tsv"
    artifact.write_text("gene_symbol\tscore\nVEGFA\t1.0\n")

    assert is_generated_artifact(artifact) is True
