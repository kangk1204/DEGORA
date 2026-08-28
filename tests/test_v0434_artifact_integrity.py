"""Adversarial regressions for process locks and artifact-set recovery."""

from __future__ import annotations

import json
import multiprocessing
import os
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from degora import provenance
from degora.excel_export import export_run_workbook, verify_run_workbook_export


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork semantics")
@pytest.mark.filterwarnings("ignore:This process .* multi-threaded:DeprecationWarning")
def test_forked_child_cannot_inherit_reentrant_output_lock(tmp_path: Path) -> None:
    context = multiprocessing.get_context("fork")
    output = tmp_path / "results"
    retry = context.Event()
    outcomes = context.Queue()

    def child() -> None:
        try:
            with provenance.output_directory_lock(output):
                outcomes.put("entered-while-parent-held")
        except provenance.OutputDirectoryBusyError:
            outcomes.put("busy")
        assert retry.wait(timeout=5)
        try:
            with provenance.output_directory_lock(output):
                outcomes.put("entered-after-release")
        except BaseException as exc:  # noqa: BLE001 - relay child failure
            outcomes.put(f"retry-error:{exc!r}")

    with provenance.output_directory_lock(output):
        process = context.Process(target=child)
        process.start()
        assert outcomes.get(timeout=5) == "busy"

    retry.set()
    assert outcomes.get(timeout=5) == "entered-after-release"
    process.join(timeout=5)
    assert not process.is_alive() and process.exitcode == 0


def test_publication_rollback_uses_copy_fallback_when_backup_rename_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    staged_one = tmp_path / "one.stage"
    staged_two = tmp_path / "two.stage"
    final_one = tmp_path / "one.txt"
    final_two = tmp_path / "two.txt"
    staged_one.write_text("new-one", encoding="utf-8")
    staged_two.write_text("new-two", encoding="utf-8")
    final_one.write_text("old-one", encoding="utf-8")
    final_two.write_text("old-two", encoding="utf-8")
    final_one.chmod(0o600)
    real_replace = provenance.os.replace

    def fail_publish_and_backup_rename(source: str | Path, destination: str | Path) -> None:
        source_path = Path(source)
        destination_path = Path(destination)
        if destination_path == final_two and source_path.suffix == ".pending":
            raise OSError("synthetic publication failure")
        if destination_path == final_two and source_path.suffix == ".backup":
            raise AssertionError("a failed pre-mutation replace must not be rolled back")
        if destination_path == final_one and source_path.suffix == ".backup":
            raise OSError("synthetic backup rename failure")
        real_replace(source, destination)

    monkeypatch.setattr(provenance.os, "replace", fail_publish_and_backup_rename)
    with pytest.raises(OSError, match="synthetic publication failure"):
        provenance.publish_staged_artifacts({staged_one: final_one, staged_two: final_two})

    assert final_one.read_text(encoding="utf-8") == "old-one"
    assert final_one.stat().st_mode & 0o777 == 0o600
    assert final_two.read_text(encoding="utf-8") == "old-two"
    assert not list(tmp_path.glob(".*.backup"))


def test_publication_rolls_back_replace_interrupted_after_target_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    staged_one = tmp_path / "one.stage"
    staged_two = tmp_path / "two.stage"
    final_one = tmp_path / "one.txt"
    final_two = tmp_path / "two.txt"
    staged_one.write_text("new-one", encoding="utf-8")
    staged_two.write_text("new-two", encoding="utf-8")
    final_one.write_text("old-one", encoding="utf-8")
    final_two.write_text("old-two", encoding="utf-8")
    real_replace = provenance.os.replace
    interrupted = False

    def replace_then_interrupt(source: str | Path, destination: str | Path) -> None:
        nonlocal interrupted
        source_path = Path(source)
        destination_path = Path(destination)
        if not interrupted and source_path.suffix == ".pending" and destination_path == final_one:
            interrupted = True
            real_replace(source, destination)
            assert not source_path.exists()
            assert destination_path.read_text(encoding="utf-8") == "new-one"
            raise KeyboardInterrupt("synthetic post-replace interruption")
        real_replace(source, destination)

    monkeypatch.setattr(provenance.os, "replace", replace_then_interrupt)
    with pytest.raises(KeyboardInterrupt, match="synthetic post-replace interruption"):
        provenance.publish_staged_artifacts({staged_one: final_one, staged_two: final_two})

    assert interrupted
    assert final_one.read_text(encoding="utf-8") == "old-one"
    assert final_two.read_text(encoding="utf-8") == "old-two"
    assert not list(tmp_path.glob(".*.pending"))
    assert not list(tmp_path.glob(".*.backup"))
    assert not staged_one.exists()
    assert not staged_two.exists()


def test_publication_preserves_recovery_backup_when_every_restore_route_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    staged_one = tmp_path / "one.stage"
    staged_two = tmp_path / "two.stage"
    final_one = tmp_path / "one.txt"
    final_two = tmp_path / "two.txt"
    staged_one.write_text("new-one", encoding="utf-8")
    staged_two.write_text("new-two", encoding="utf-8")
    final_one.write_text("old-one", encoding="utf-8")
    final_two.write_text("old-two", encoding="utf-8")
    real_replace = provenance.os.replace
    real_copy2 = provenance.shutil.copy2

    def fail_publish_and_backup_rename(source: str | Path, destination: str | Path) -> None:
        source_path = Path(source)
        destination_path = Path(destination)
        if destination_path == final_two and source_path.suffix == ".pending":
            raise OSError("synthetic publication failure")
        if destination_path == final_one and source_path.suffix == ".backup":
            raise OSError("synthetic backup rename failure")
        real_replace(source, destination)

    def fail_backup_copy(source: str | Path, destination: str | Path, *args, **kwargs):
        if Path(source).suffix == ".backup" and Path(destination) == final_one:
            raise OSError("synthetic backup copy failure")
        return real_copy2(source, destination, *args, **kwargs)

    monkeypatch.setattr(provenance.os, "replace", fail_publish_and_backup_rename)
    monkeypatch.setattr(provenance.shutil, "copy2", fail_backup_copy)
    with pytest.raises(RuntimeError, match="rollback was incomplete") as error:
        provenance.publish_staged_artifacts({staged_one: final_one, staged_two: final_two})

    backups = list(tmp_path.glob(".*.backup"))
    assert len(backups) == 1
    assert backups[0].read_text(encoding="utf-8") == "old-one"
    assert str(backups[0]) in str(error.value)
    assert final_two.read_text(encoding="utf-8") == "old-two"


def _minimal_result_dir(path: Path) -> Path:
    path.mkdir()
    with sqlite3.connect(path / "degora_scores.db") as connection:
        pd.DataFrame(
            {
                "gene_symbol": ["ISG15"],
                "quality_weighted_degora_rank": [1],
                "quality_weighted_degora_score": [1.0],
            }
        ).to_sql("genes", connection, index=False)
        pd.DataFrame({"gene_symbol": ["ISG15"], "study_id": ["S1"]}).to_sql(
            "gene_evidence", connection, index=False
        )
        pd.DataFrame({"study_id": ["S1"], "source_unit_id": ["U1"]}).to_sql(
            "studies", connection, index=False
        )
        pd.DataFrame({"key": ["degora_version"], "value": ["test"]}).to_sql(
            "meta", connection, index=False
        )
    (path / "degora_score_metadata.json").write_text("{}\n", encoding="utf-8")
    return path


def test_workbook_verifier_requires_and_links_manifest_own_sidecars(tmp_path: Path) -> None:
    result = export_run_workbook(_minimal_result_dir(tmp_path / "results"), command="pytest manifest audit")
    workbook = Path(result["output"])
    manifest = Path(result["manifest"])
    source = provenance.artifact_source_path(manifest)
    audit = provenance.artifact_provenance_path(manifest)
    verify_run_workbook_export(workbook)

    source_bytes = source.read_bytes()
    source.unlink()
    with pytest.raises(ValueError, match="artifact is missing"):
        verify_run_workbook_export(workbook)
    source.write_bytes(source_bytes)

    audit_bytes = audit.read_bytes()
    audit.unlink()
    with pytest.raises(ValueError, match="artifact is missing"):
        verify_run_workbook_export(workbook)
    audit.write_bytes(audit_bytes)

    source.write_text("different command\n", encoding="utf-8")
    with pytest.raises(ValueError, match="source sidecar does not match"):
        verify_run_workbook_export(workbook)
    source.write_bytes(source_bytes)

    audit_payload = json.loads(audit_bytes)
    audit_payload["artifact_sha256"] = "0" * 64
    audit.write_text(json.dumps(audit_payload), encoding="utf-8")
    with pytest.raises(ValueError, match="provenance sidecar does not match"):
        verify_run_workbook_export(workbook)
