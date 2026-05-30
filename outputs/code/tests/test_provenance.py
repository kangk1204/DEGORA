from __future__ import annotations

import json
from pathlib import Path

import pytest

from degora.provenance import artifact_provenance_path, artifact_source_path, shell_command, write_source_sidecar
from scripts.check_generated_sources import check_sources


def test_write_source_sidecar_preserves_command_first_line_and_json_record(tmp_path: Path) -> None:
    artifact = tmp_path / "result.csv"
    input_path = tmp_path / "input.csv"
    artifact.write_text("gene,value\nVEGFA,1\n")
    input_path.write_text("gene\nVEGFA\n")
    command = "make -C outputs/code slice ITER=99"

    write_source_sidecar(artifact, command, inputs=[input_path], metadata={"generator": "test"})

    assert artifact_source_path(artifact).read_text().splitlines()[0] == command
    record = json.loads(artifact_provenance_path(artifact).read_text())
    assert record["command"] == command
    assert record["artifact_path"] == str(artifact.resolve())
    assert record["inputs"][0]["path"] == str(input_path.resolve())
    assert record["inputs"][0]["exists"] is True
    assert record["metadata"] == {"generator": "test"}


def test_provenance_commands_are_single_line_and_shell_quoted(tmp_path: Path) -> None:
    artifact = tmp_path / "result.csv"
    artifact.write_text("x\n")

    command = shell_command(["make", "-C", "outputs/code", "slice", "CATALOG=space dir/evil;name.csv"])
    write_source_sidecar(artifact, command)

    source = artifact_source_path(artifact).read_text()
    assert source.count("\n") == 1
    assert "'CATALOG=space dir/evil;name.csv'" in source

    with pytest.raises(ValueError, match="single line"):
        write_source_sidecar(artifact, "make ok\nmalicious second command")


def test_generated_source_checker_flags_missing_or_unexpected_sidecars(tmp_path: Path) -> None:
    good = tmp_path / "iter-1" / "good.csv"
    bad = tmp_path / "iter-1" / "bad.json"
    good.parent.mkdir(parents=True)
    good.write_text("x\n")
    bad.write_text("{}\n")
    good.with_suffix(".csv.source").write_text("make -C outputs/code slice ITER=1\n")

    errors = check_sources(tmp_path)

    assert errors == [f"missing source sidecar: {bad.with_suffix('.json.source')}"]
