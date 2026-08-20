"""Provenance sidecars for generated DEGORA artifacts."""

from __future__ import annotations

import errno
import hashlib
import json
import os
import re
import shlex
import shutil
import tempfile
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

# WSL drvfs (/mnt/<drive>) intermittently raises EINVAL/EIO when opening a small
# sidecar file for writing right after a large write to the same directory (a known
# Windows-drive mount quirk). The payload is valid and the write succeeds on retry,
# so sidecar writes use a short bounded backoff instead of crashing the pipeline. The
# project documents WSL as the supported Windows path, so beginners running
# `degora run` over /mnt would otherwise hit a spurious mid-pipeline failure.
_SIDECAR_WRITE_RETRY_ERRNOS = frozenset({errno.EINVAL, errno.EIO})
_SIDECAR_WRITE_ATTEMPTS = 6
_FIXED_OOXML_DATETIME = datetime(2000, 1, 1)
_FIXED_ZIP_TIMESTAMP = (2000, 1, 1, 0, 0, 0)
_FIXED_ZIP_CREATE_SYSTEM = 3
_SECRET_METADATA_KEY_RE = re.compile(
    r"(?i)(^|[_-])(api[_-]?(?:key|token)|access[_-]?token|token|secret|password|credentials?)($|[_-])"
)
_SENSITIVE_RELATIVE_PATH_RE = re.compile(r"(?i)(?:^|/)(?:users|home)/[^/]+(?:/|$)")
EXTERNAL_PATH_PREFIX = "external-redacted://"


def reproducible_generated_at() -> str:
    """Return a deterministic ISO-8601 UTC timestamp for hash-pinned manifests.

    Figure manifests are sha256-pinned in the release package, so embedding a live
    wall-clock ``datetime.now()`` makes every re-render break its pinned hash even when
    no numeric content changed. This honors ``SOURCE_DATE_EPOCH`` (the reproducible-builds
    convention) so a re-render under a pinned epoch reproduces a byte-identical manifest;
    it falls back to the current UTC time only when the variable is unset (interactive or
    one-off renders).
    """

    epoch = os.environ.get("SOURCE_DATE_EPOCH", "").strip()
    if epoch:
        try:
            return datetime.fromtimestamp(int(epoch), tz=timezone.utc).isoformat()
        except (ValueError, OverflowError, OSError):
            pass
    return datetime.now(timezone.utc).isoformat()


def set_reproducible_workbook_properties(workbook: Any) -> None:
    """Pin mutable OOXML workbook metadata before saving it."""

    workbook.properties.created = _FIXED_OOXML_DATETIME
    workbook.properties.modified = _FIXED_OOXML_DATETIME
    workbook.properties.creator = "DEGORA"
    workbook.properties.lastModifiedBy = "DEGORA"


def normalize_ooxml_zip(path: str | Path) -> None:
    """Rewrite an XLSX/DOCX archive with stable metadata, order, and timestamps."""

    archive_path = Path(path)
    if not zipfile.is_zipfile(archive_path):
        raise ValueError(f"not an OOXML ZIP archive: {archive_path}")

    entries: list[tuple[str, bytes, int]] = []
    with zipfile.ZipFile(archive_path, "r") as source:
        for info in sorted(source.infolist(), key=lambda item: item.filename):
            data = source.read(info.filename)
            if info.filename == "docProps/core.xml":
                text = data.decode("utf-8")
                text = re.sub(
                    r"(<dcterms:created[^>]*>)[^<]+(</dcterms:created>)",
                    r"\g<1>2000-01-01T00:00:00Z\g<2>",
                    text,
                )
                text = re.sub(
                    r"(<dcterms:modified[^>]*>)[^<]+(</dcterms:modified>)",
                    r"\g<1>2000-01-01T00:00:00Z\g<2>",
                    text,
                )
                data = text.encode("utf-8")
            entries.append((info.filename, data, info.external_attr))

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{archive_path.name}.",
        suffix=".deterministic.tmp",
        dir=archive_path.parent,
    )
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        with zipfile.ZipFile(
            temporary,
            "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as target:
            for name, data, external_attr in entries:
                info = zipfile.ZipInfo(name, _FIXED_ZIP_TIMESTAMP)
                info.create_system = _FIXED_ZIP_CREATE_SYSTEM
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = external_attr
                target.writestr(info, data)
        os.replace(temporary, archive_path)
    finally:
        temporary.unlink(missing_ok=True)


def _resilient_write_text(path: Path, text: str) -> None:
    # Write to a sibling temp file and atomically rename it onto the target, so a
    # crash or interrupted write can never leave a truncated/invalid sidecar for the
    # provenance checker to flag. os.replace is atomic on the same filesystem.
    last: OSError | None = None
    for attempt in range(_SIDECAR_WRITE_ATTEMPTS):
        tmp: Path | None = None
        try:
            fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
            tmp = Path(tmp_name)
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
            tmp.replace(path)
            return
        except OSError as exc:
            if tmp is not None:
                tmp.unlink(missing_ok=True)
            if exc.errno not in _SIDECAR_WRITE_RETRY_ERRNOS:
                raise
            last = exc
            # No point sleeping after the final attempt -- the loop is exhausted and
            # the code re-raises immediately, so the trailing backoff is dead time.
            if attempt < _SIDECAR_WRITE_ATTEMPTS - 1:
                time.sleep(0.25 * (attempt + 1))
    if last is None:  # defensive: the loop can exit only after a retryable OSError
        raise RuntimeError("sidecar write retry loop exited without an error")
    raise last


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _external_path_reference(path: Path) -> str:
    digest = hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:12]
    return f"{EXTERNAL_PATH_PREFIX}{digest}/{path.name}"


def is_external_path_reference(value: Any) -> bool:
    return isinstance(value, str) and value.startswith(EXTERNAL_PATH_PREFIX)


def portable_path(path: str | Path, base_dir: str | Path) -> str:
    """Return a path relative to an explicitly documented artifact base.

    Generated outputs are commonly shared or moved between machines. Recording
    paths relative to a known base avoids embedding usernames and workstation
    layouts while retaining a replayable location. A path outside the local run
    bundle, or a cross-volume Windows path, is replaced by a unique, explicit
    ``external-redacted://`` reference rather than a plausible-but-false path.
    """

    resolved = Path(path).resolve()
    base = Path(base_dir).resolve()
    try:
        relative = os.path.relpath(resolved, start=base)
    except ValueError:
        return _external_path_reference(resolved)
    portable = Path(relative).as_posix()
    bundle_root = base.parent
    try:
        resolved.relative_to(bundle_root)
    except ValueError:
        return _external_path_reference(resolved)
    if _SENSITIVE_RELATIVE_PATH_RE.search(portable):
        return _external_path_reference(resolved)
    return portable


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
    artifact_content_path: str | Path | None = None,
    inputs: Iterable[str | Path] = (),
    allow_missing_inputs: Iterable[str | Path] = (),
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a deterministic provenance record for one generated artifact.

    ``allow_missing_inputs`` is reserved for large, intentionally unretained
    regeneration inputs. Their size and digest are recorded and verified whenever
    the file is present, while a clean checkout may omit the file.
    """

    artifact = Path(artifact_path).resolve()
    path_base = artifact.parent
    content = Path(artifact_content_path).resolve() if artifact_content_path is not None else artifact
    allowed_missing = {Path(raw_input).resolve() for raw_input in allow_missing_inputs}
    input_records = []
    for raw_input in inputs:
        input_path = Path(raw_input).resolve()
        record: dict[str, Any] = {
            "path": portable_path(input_path, path_base),
            "exists": input_path.exists(),
            "required_for_audit": input_path not in allowed_missing,
        }
        record["path_replayable"] = not is_external_path_reference(record["path"])
        if input_path.is_file():
            record.update({"size_bytes": input_path.stat().st_size, "sha256": _sha256(input_path)})
        input_records.append(record)

    portable_regeneration_command = portable_command(command, path_base)
    record = {
        "artifact_path": portable_path(artifact, path_base),
        "path_base": "artifact_directory",
        "command": portable_regeneration_command,
        "command_replayable": EXTERNAL_PATH_PREFIX not in portable_regeneration_command,
        "inputs": input_records,
    }
    if not record["command_replayable"]:
        record["replay_warning"] = (
            "One or more external paths were redacted; replace each external-redacted:// reference "
            "with the corresponding local file before replaying the command."
        )
    if content.is_file():
        record.update({"artifact_size_bytes": content.stat().st_size, "artifact_sha256": _sha256(content)})
    if metadata:
        record["metadata"] = sanitize_metadata(metadata, path_base)
    return record


def sanitize_metadata(value: Any, base_dir: str | Path | None = None) -> Any:
    """Redact secret fields and portableize absolute metadata path values."""

    if isinstance(value, dict):
        redacted: dict[str, Any] = {}
        for key, item in value.items():
            key_text = str(key)
            if _SECRET_METADATA_KEY_RE.search(key_text):
                redacted[key_text] = "[redacted]"
            else:
                redacted[key_text] = sanitize_metadata(item, base_dir)
        return redacted
    if isinstance(value, list):
        return [sanitize_metadata(item, base_dir) for item in value]
    if isinstance(value, tuple):
        return [sanitize_metadata(item, base_dir) for item in value]
    if base_dir is not None and isinstance(value, (str, Path)) and Path(value).is_absolute():
        return portable_path(value, base_dir)
    return value


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


def portable_command(command: str, base_dir: str | Path) -> str:
    """Rewrite absolute path arguments relative to ``base_dir``.

    DEGORA commands are emitted with :func:`shell_command`, so parsing and
    re-quoting them preserves argument boundaries. Both standalone absolute
    arguments and ``NAME=/absolute/path`` settings are handled.
    """

    command = _validate_single_line_command(command)
    try:
        parts = shlex.split(command, posix=True)
    except ValueError as exc:
        raise ValueError(f"provenance command is not valid shell syntax: {exc}") from exc

    portable_parts: list[str] = []
    for part in parts:
        whole_candidate = Path(part)
        if whole_candidate.is_absolute():
            portable_parts.append(portable_path(whole_candidate, base_dir))
            continue
        prefix = ""
        value = part
        if "=" in part:
            candidate_prefix, candidate_value = part.split("=", 1)
            if candidate_prefix and candidate_value:
                prefix = f"{candidate_prefix}="
                value = candidate_value
        candidate = Path(value)
        if candidate.is_absolute():
            value = portable_path(candidate, base_dir)
        portable_parts.append(f"{prefix}{value}")
    return shell_command(portable_parts)


def source_sidecar_payloads(
    artifact_path: str | Path,
    command: str,
    *,
    artifact_content_path: str | Path | None = None,
    inputs: Iterable[str | Path] = (),
    allow_missing_inputs: Iterable[str | Path] = (),
    metadata: dict[str, Any] | None = None,
    write_json: bool = True,
) -> tuple[str, str | None]:
    """Return ``.source`` and optional provenance JSON text for an artifact.

    ``artifact_path`` is the final published artifact path recorded in provenance.
    ``artifact_content_path`` lets transactional exporters hash an unpublished
    staging file with identical bytes before anything is moved into place.
    """

    artifact = Path(artifact_path).resolve()
    command = portable_command(command, artifact.parent)
    json_text: str | None = None
    if write_json:
        record = provenance_record(
            artifact_path,
            command,
            artifact_content_path=artifact_content_path,
            inputs=inputs,
            allow_missing_inputs=allow_missing_inputs,
            metadata=metadata,
        )
        json_text = json.dumps(record, indent=2, sort_keys=True) + "\n"
    return command + "\n", json_text


def write_source_sidecar(
    artifact_path: str | Path,
    command: str,
    *,
    artifact_content_path: str | Path | None = None,
    inputs: Iterable[str | Path] = (),
    allow_missing_inputs: Iterable[str | Path] = (),
    metadata: dict[str, Any] | None = None,
    write_json: bool = True,
) -> None:
    """Write command-first ``.source`` and optional JSON provenance sidecars.

    The ``.source`` sidecar keeps a portable regeneration command on the first
    line, with absolute paths made relative to the sidecar directory.
    ``allow_missing_inputs`` has the same narrow semantics as in
    :func:`provenance_record`.
    """

    artifact = Path(artifact_path)
    source_text, json_text = source_sidecar_payloads(
        artifact,
        command,
        artifact_content_path=artifact_content_path,
        inputs=inputs,
        allow_missing_inputs=allow_missing_inputs,
        metadata=metadata,
        write_json=write_json,
    )
    _resilient_write_text(artifact_source_path(artifact), source_text)
    if json_text is not None:
        _resilient_write_text(artifact_provenance_path(artifact), json_text)


def publish_staged_artifacts(staged_to_final: dict[Path, Path]) -> None:
    """Publish a complete artifact set, restoring all known targets on failure."""

    if len(set(staged_to_final.values())) != len(staged_to_final):
        raise ValueError("staged artifact publication target paths must be unique")
    if not staged_to_final:
        return

    prepared: list[tuple[Path, Path, Path | None]] = []
    published: list[tuple[Path, Path | None]] = []
    try:
        # Copy every staged artifact to a sibling pending file first. This keeps the
        # final os.replace on one filesystem even when a caller places the SQLite DB
        # outside output_dir, and ensures copy failures occur before any target changes.
        for staged_path, final_path in staged_to_final.items():
            final_path.parent.mkdir(parents=True, exist_ok=True)
            if final_path.exists() and not final_path.is_file():
                raise IsADirectoryError(f"artifact publication target is not a file: {final_path}")
            pending_path: Path | None = None
            backup_path: Path | None = None
            try:
                pending_fd, pending_name = tempfile.mkstemp(
                    prefix=f".{final_path.name}.", suffix=".pending", dir=final_path.parent
                )
                os.close(pending_fd)
                pending_path = Path(pending_name)
                shutil.copy2(staged_path, pending_path)
                if final_path.exists():
                    backup_fd, backup_name = tempfile.mkstemp(
                        prefix=f".{final_path.name}.", suffix=".backup", dir=final_path.parent
                    )
                    os.close(backup_fd)
                    backup_path = Path(backup_name)
                    backup_path.unlink()
                    try:
                        os.link(final_path, backup_path)
                    except OSError:
                        shutil.copy2(final_path, backup_path)
                prepared.append((pending_path, final_path, backup_path))
            except BaseException:
                if pending_path is not None:
                    pending_path.unlink(missing_ok=True)
                if backup_path is not None:
                    backup_path.unlink(missing_ok=True)
                raise

        try:
            for pending_path, final_path, backup_path in prepared:
                os.replace(pending_path, final_path)
                published.append((final_path, backup_path))
        except BaseException:
            for final_path, backup_path in reversed(published):
                if backup_path is not None and backup_path.exists():
                    os.replace(backup_path, final_path)
                else:
                    final_path.unlink(missing_ok=True)
            raise
    finally:
        for staged_path in staged_to_final:
            staged_path.unlink(missing_ok=True)
        for pending_path, _final_path, backup_path in prepared:
            pending_path.unlink(missing_ok=True)
            if backup_path is not None:
                backup_path.unlink(missing_ok=True)
