"""Making published text safe to open in a spreadsheet.

Kept free of other DEGORA imports so every writer can use it: excel_export needs
score_db, and score_db needs this, so anything living in excel_export would close
an import cycle.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Literal

import pandas as pd

EXCEL_ERROR_LITERALS = frozenset(
    {
        "#DIV/0!",
        "#N/A",
        "#NAME?",
        "#NULL!",
        "#NUM!",
        "#REF!",
        "#VALUE!",
        "#GETTING_DATA",
    }
)

CSV_FORMULA_PREFIXES = ("=", "+", "-", "@")
CSV_FORMULA_PREFIX_WHITESPACE = " \t\r\n"
FORMULA_GUARD_METADATA_KEY = "csv_formula_guard"
FORMULA_GUARD_SCHEME = "reversible_apostrophe_prefix_v1"
FormulaGuardState = Literal["valid", "absent", "unmarked", "invalid"]


def _formula_payload(value: str) -> str:
    # Include apostrophes in the removable prefix so the transform is bijective:
    # raw ``=x`` becomes ``'=x`` while raw ``'=x`` becomes ``''=x``. A reader
    # may therefore remove exactly one guard apostrophe without conflating them.
    return value.lstrip(CSV_FORMULA_PREFIX_WHITESPACE + "'")


def is_formula_like_text(value: Any) -> bool:
    """Return whether a text value needs a spreadsheet formula guard."""

    if not isinstance(value, str):
        return False
    payload = _formula_payload(value)
    return bool(
        payload
        and (
            payload.startswith(CSV_FORMULA_PREFIXES)
            or payload.upper() in EXCEL_ERROR_LITERALS
        )
    )


def neutralize_formula_cell(value: Any) -> Any:
    """Prefix one reversible apostrophe when spreadsheet software could execute text."""

    return "'" + value if is_formula_like_text(value) else value


def neutralize_formula_text(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy whose text cells cannot be read as a spreadsheet formula.

    The published CSVs carry gene identifiers and free text that came from tables
    DEGORA downloaded from public repositories, so a cell reading
    ``=HYPERLINK("http://...")`` is an input, not a hypothetical. Opening such a
    file in Excel executes it. The XLSX writer and the discovery exports already
    guard against this; the score table, the consensus table and the harmonized
    table - the files most likely to be shared - did not.

    Only text is touched. Numeric columns are left exactly as they are, because
    prefixing a negative number would corrupt the value it publishes.
    """

    guarded = frame.copy()
    guarded.columns = [neutralize_formula_cell(column) for column in guarded.columns]
    for column in guarded.columns:
        if not (guarded[column].dtype == object or pd.api.types.is_string_dtype(guarded[column])):
            continue
        guarded[column] = guarded[column].map(neutralize_formula_cell)
    return guarded


def restore_formula_text(frame: pd.DataFrame) -> pd.DataFrame:
    """Reverse one DEGORA formula guard without changing source identifiers."""

    restored = frame.copy()
    restored.columns = [_restore_formula_cell(column) for column in restored.columns]
    for column in restored.columns:
        if not (restored[column].dtype == object or pd.api.types.is_string_dtype(restored[column])):
            continue
        restored[column] = restored[column].map(_restore_formula_cell)
    return restored


def _restore_formula_cell(value: Any) -> Any:
    if _is_formula_guard_cell(value):
        return value[1:]
    return value


def _is_formula_guard_cell(value: Any) -> bool:
    return isinstance(value, str) and value.startswith("'") and is_formula_like_text(value[1:])


def contains_formula_guard(frame: pd.DataFrame, columns: Iterable[str] | None = None) -> bool:
    """Return whether a table contains text that could be a DEGORA guard.

    Without a matching provenance digest, ``'=BAD()`` is ambiguous: it may be
    protected ``=BAD()`` or a literal source identifier that already began with
    an apostrophe.  Callers must not silently choose either interpretation.

    ``columns`` limits the question to the columns that will be used. A DEGORA
    artifact being re-read is all DEGORA's, so every column counts; an author's
    table is read for its gene, effect and p-value columns, and an apostrophe in
    an unmapped notes column changes nothing that is scored.
    """

    if columns is None:
        if any(_is_formula_guard_cell(column) for column in frame.columns):
            return True
        text = frame.select_dtypes(include=["object", "string"])
    else:
        wanted = [name for name in columns if name and name in frame.columns]
        if not wanted:
            return False
        text = frame[wanted].select_dtypes(include=["object", "string"])
    return any(_is_formula_guard_cell(value) for value in text.to_numpy(dtype=object, copy=False).flat)


def formula_guard_metadata() -> dict[str, str]:
    """Return the provenance marker for a guarded CSV/TSV artifact."""

    return {FORMULA_GUARD_METADATA_KEY: FORMULA_GUARD_SCHEME}


def _artifact_sha256(path: Path) -> str:
    digest_builder = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest_builder.update(chunk)
    return digest_builder.hexdigest()


def formula_guard_state(path: str | Path) -> FormulaGuardState:
    """Classify formula-guard provenance for a CSV/TSV artifact."""

    artifact = Path(path)
    provenance = artifact.with_suffix(artifact.suffix + ".provenance.json")
    if not provenance.exists():
        return "absent"
    try:
        record = json.loads(provenance.read_text(encoding="utf-8"))
        if not isinstance(record, dict):
            return "invalid"
        expected_digest = str(record.get("artifact_sha256") or "")
        metadata = record.get("metadata") or {}
        if not isinstance(metadata, dict):
            return "invalid"
        scheme = metadata.get(FORMULA_GUARD_METADATA_KEY)
        if scheme is None:
            return "unmarked"
        if scheme != FORMULA_GUARD_SCHEME or not expected_digest:
            return "invalid"
        digest = _artifact_sha256(artifact)
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return "invalid"
    return "valid" if digest == expected_digest else "invalid"


def has_current_formula_guard(path: str | Path) -> bool:
    """Return whether a formula guard has a matching marker and digest."""

    return formula_guard_state(path) == "valid"


def restore_formula_text_if_marked(
    frame: pd.DataFrame, path: str | Path, *, columns: Iterable[str] | None = None
) -> pd.DataFrame:
    """Restore a verified guard and reject ambiguous or tampered guarded text.

    ``columns`` names the columns the caller will use; see contains_formula_guard.
    """

    state = formula_guard_state(path)
    if state == "valid":
        return restore_formula_text(frame)
    guarded = contains_formula_guard(frame, columns)
    if not guarded:
        return frame
    name = Path(path).name
    if state == "invalid":
        raise ValueError(
            f"{name} has invalid formula-guard provenance: its marker, digest, "
            "or provenance JSON does not match the artifact"
        )
    if state == "absent":
        # No sidecar at all: DEGORA did not write this file, so there is no
        # sidecar to restore, and saying so sent readers looking for one.
        raise ValueError(
            f"{name} has a value beginning with an apostrophe followed by formula-like text "
            "in a column DEGORA will use, and no matching DEGORA formula-guard provenance to say "
            "whether the apostrophe is part of the identifier. If this is an author's table, remove the "
            "apostrophe in the source (or map a different column); if it is a DEGORA output, "
            "keep its .provenance.json beside it."
        )
    raise ValueError(
        f"{name} contains apostrophe-guarded formula-like text but has no "
        "matching DEGORA formula-guard provenance; restore the matching "
        ".provenance.json sidecar or provide the original raw identifiers"
    )
