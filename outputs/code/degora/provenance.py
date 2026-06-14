"""Provenance sidecars for generated DEGORA artifacts."""

from __future__ import annotations

import errno
import hashlib
import json
import shlex
import time
from pathlib import Path
from typing import Any, Iterable

# WSL drvfs (/mnt/<drive>) intermittently raises EINVAL/EIO when opening a small
# sidecar file for writing right after a large write to the same directory (a known
# Windows-drive mount quirk). The payload is valid and the write succeeds on retry,
# so sidecar writes use a short bounded backoff instead of crashing the pipeline. The
# project documents WSL as the supported Windows path, so beginners running `make
# paper` over /mnt would otherwise hit a spurious mid-pipeline failure.
_SIDECAR_WRITE_RETRY_ERRNOS = frozenset({errno.EINVAL, errno.EIO})
_SIDECAR_WRITE_ATTEMPTS = 6
_REPO_ROOT = Path(__file__).resolve().parents[3]


def _resilient_write_text(path: Path, text: str) -> None:
    last: OSError | None = None
    for attempt in range(_SIDECAR_WRITE_ATTEMPTS):
        try:
            with path.open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
            return
        except OSError as exc:
            if exc.errno not in _SIDECAR_WRITE_RETRY_ERRNOS:
                raise
            last = exc
            time.sleep(0.25 * (attempt + 1))
    assert last is not None
    raise last


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(_REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def artifact_source_path(path: str | Path) -> Path:
    """Return the conventional ``<artifact>.<suffix>.source`` sidecar path."""

    artifact = Path(path)
    return artifact.with_suffix(artifact.suffix + ".source")


def artifact_provenance_path(path: str | Path) -> Path:
    """Return the JSON provenance sidecar path for an artifact."""

    artifact = Path(path)
    return artifact.with_suffix(artifact.suffix + ".provenance.json")


def provenance_record(
    artifact_path: str | Path,
    command: str,
    *,
    inputs: Iterable[str | Path] = (),
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a deterministic provenance record for one generated artifact."""

    artifact = Path(artifact_path).resolve()
    input_records = []
    for raw_input in inputs:
        input_path = Path(raw_input).resolve()
        record: dict[str, Any] = {"path": _portable_path(input_path), "exists": input_path.exists()}
        if input_path.is_file():
            record.update({"size_bytes": input_path.stat().st_size, "sha256": _sha256(input_path)})
        input_records.append(record)

    record = {
        "artifact_path": _portable_path(artifact),
        "command": command,
        "inputs": input_records,
    }
    if artifact.is_file():
        record.update({"artifact_size_bytes": artifact.stat().st_size, "artifact_sha256": _sha256(artifact)})
    if metadata:
        record["metadata"] = metadata
    return record


def _validate_single_line_command(command: str) -> str:
    if "\n" in command or "\r" in command:
        raise ValueError("provenance command must be a single line")
    return command


def shell_command(args: Iterable[str | Path | int | float]) -> str:
    """Return a shell-quoted, single-line regeneration command."""

    parts = [str(arg) for arg in args]
    for part in parts:
        _validate_single_line_command(part)
    return shlex.join(parts)


def write_source_sidecar(
    artifact_path: str | Path,
    command: str,
    *,
    inputs: Iterable[str | Path] = (),
    metadata: dict[str, Any] | None = None,
    write_json: bool = True,
) -> None:
    """Write command-first ``.source`` and optional JSON provenance sidecars.

    The ``.source`` sidecar intentionally keeps the exact regeneration command on
    the first line to preserve the repository's source-audit contract.
    """

    command = _validate_single_line_command(command)
    artifact = Path(artifact_path)
    _resilient_write_text(artifact_source_path(artifact), command + "\n")
    if write_json:
        record = provenance_record(artifact, command, inputs=inputs, metadata=metadata)
        _resilient_write_text(artifact_provenance_path(artifact), json.dumps(record, indent=2, sort_keys=True) + "\n")
