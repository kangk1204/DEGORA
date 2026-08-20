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
