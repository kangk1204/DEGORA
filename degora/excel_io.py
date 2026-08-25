"""Helpers for beginner-facing Excel config workbooks."""

from __future__ import annotations

from typing import Any

import pandas as pd


def _cell_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def leading_comment_rows(excel: Any, sheet_name: str) -> int:
    """Count leading note rows in a DEGORA config sheet.

    Template sheets may start with one or more human-readable rows whose first
    non-empty cell begins with '#'. Those rows are not data headers.
    """

    preview = pd.read_excel(excel, sheet_name=sheet_name, header=None, nrows=20)
    count = 0
    for _, row in preview.iterrows():
        texts = [_cell_text(value) for value in row.tolist()]
        nonempty = [text for text in texts if text]
        if not nonempty:
            count += 1
            continue
        if nonempty[0].startswith("#"):
            count += 1
            continue
        break
    return count


def read_config_sheet(excel: Any, sheet_name: str) -> pd.DataFrame:
    """Read a DEGORA config sheet while ignoring leading '#'-note rows."""

    return pd.read_excel(excel, sheet_name=sheet_name, skiprows=leading_comment_rows(excel, sheet_name))


# GoldPanel `locked` flags. A column typed as 1/0 with any blank cell arrives from
# pandas as float64, so "1" reads back as "1.0" - which used to match nothing in
# the truthy set and dropped exactly the rows the reader had marked, while the
# blank rows beside them were kept. Integer-valued floats are folded back first.
LOCKED_TRUE = frozenset({"1", "true", "t", "yes", "y", "locked", "lock", "include", "included"})
LOCKED_FALSE = frozenset({"0", "false", "f", "no", "n", "unlocked", "exclude", "excluded"})


def normalize_locked_flag(value: Any) -> str:
    """Return 'yes', 'no' or '' (blank, which counts as locked) for one cell."""

    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return ""
    if isinstance(value, (bool,)):
        return "yes" if value else "no"
    text = str(value).strip().lower()
    if not text:
        return ""
    if text not in LOCKED_TRUE and text not in LOCKED_FALSE:
        try:
            number = float(text)
        except ValueError:
            number = None
        if number is not None and number.is_integer():
            text = str(int(number))
    if text in LOCKED_TRUE:
        return "yes"
    if text in LOCKED_FALSE:
        return "no"
    # Anything else is not a flag DEGORA understands; treating it as locked would
    # silently include a row the reader may have meant to exclude.
    return "no"


def locked_panel_mask(values: pd.Series) -> pd.Series:
    """Rows a GoldPanel keeps: locked=yes, or a blank locked cell."""

    flags = values.map(normalize_locked_flag)
    return flags.isin({"yes", ""})
