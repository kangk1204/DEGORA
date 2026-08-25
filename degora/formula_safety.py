"""Making published text safe to open in a spreadsheet.

Kept free of other DEGORA imports so every writer can use it: excel_export needs
score_db, and score_db needs this, so anything living in excel_export would close
an import cycle.
"""

from __future__ import annotations

from typing import Any

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
    for column in guarded.columns:
        if not (guarded[column].dtype == object or pd.api.types.is_string_dtype(guarded[column])):
            continue
        guarded[column] = guarded[column].map(_neutralized_cell)
    return guarded


def _neutralized_cell(value: Any) -> Any:
    # Anything that is not text cannot be mistaken for a formula, and must keep
    # its own type so the column still writes as numbers.
    if not isinstance(value, str):
        return value
    stripped = value.lstrip(CSV_FORMULA_PREFIX_WHITESPACE)
    if stripped.startswith(CSV_FORMULA_PREFIXES) or stripped.upper() in EXCEL_ERROR_LITERALS:
        return "'" + value
    return value


