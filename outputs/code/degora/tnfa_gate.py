"""TNF-alpha corpus gate checks for the iteration-10 cross-domain lane."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class TnfaGateResult:
    n_active_rows: int
    n_independent_units: int
    n_candidates: int
    passed: bool
    scarcity_triggered: bool
    reasons: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_active_rows": self.n_active_rows,
            "n_independent_units": self.n_independent_units,
            "n_candidates": self.n_candidates,
            "passed": self.passed,
            "scarcity_triggered": self.scarcity_triggered,
            "reasons": list(self.reasons),
        }


def _is_active_status(value: Any) -> bool:
    if value is None or pd.isna(value):
        return False
    return str(value).strip().lower() == "active"


def _first_nonempty(row: pd.Series, columns: tuple[str, ...]) -> str:
    for column in columns:
        if column not in row.index:
            continue
        value = row[column]
        if value is None or pd.isna(value):
            continue
        text = str(value).strip()
        if text:
            return text
    return ""


def independent_source_units(active_rows: pd.DataFrame) -> set[str]:
    """Return deterministic paper/source-family units for active TNF rows.

    The TNF lane's planned catalog has an explicit source_family_id, while
    screening ledgers may only have paper/dataset fields.  Count independence
    by the most specific stable family key available and never by row count.
    """

    units: set[str] = set()
    for _, row in active_rows.iterrows():
        source_family = _first_nonempty(row, ("source_family_id", "candidate_family", "paper_or_dataset"))
        paper = _first_nonempty(row, ("paper_id", "doi_or_accession", "paper_or_dataset"))
        if source_family or paper:
            units.add(f"{paper}::{source_family}")
    return units


def evaluate_tnfa_gate(
    rows: pd.DataFrame,
    *,
    min_active_rows: int = 5,
    min_independent_units: int = 5,
    scarcity_candidate_threshold: int = 20,
) -> TnfaGateResult:
    """Evaluate the TNF-only activation gate or scarcity branch."""

    status_column = "status"
    if status_column not in rows.columns and "eligibility_status" in rows.columns:
        status_column = "eligibility_status"
    if status_column not in rows.columns:
        raise ValueError("TNF gate input is missing required column: status or eligibility_status")

    active = rows.loc[rows[status_column].map(_is_active_status)].copy()
    n_active = int(len(active))
    n_units = int(len(independent_source_units(active)))
    n_candidates = int(len(rows))

    reasons: list[str] = []
    if n_active < min_active_rows:
        reasons.append(f"active_rows_below_{min_active_rows}")
    if n_units < min_independent_units:
        reasons.append(f"independent_units_below_{min_independent_units}")

    passed = not reasons
    scarcity_triggered = (not passed) and n_candidates >= scarcity_candidate_threshold
    return TnfaGateResult(
        n_active_rows=n_active,
        n_independent_units=n_units,
        n_candidates=n_candidates,
        passed=passed,
        scarcity_triggered=scarcity_triggered,
        reasons=tuple(reasons),
    )


def evaluate_tnfa_gate_csv(path: str | Path, **kwargs: Any) -> TnfaGateResult:
    return evaluate_tnfa_gate(pd.read_csv(path), **kwargs)
