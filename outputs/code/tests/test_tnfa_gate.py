from __future__ import annotations

import pandas as pd

from degora.tnfa_gate import evaluate_tnfa_gate, independent_source_units


def test_tnfa_gate_passes_with_five_active_independent_units() -> None:
    rows = pd.DataFrame(
        {
            "status": ["active"] * 5,
            "paper_id": [f"P{i}" for i in range(5)],
            "source_family_id": [f"F{i}" for i in range(5)],
        }
    )

    result = evaluate_tnfa_gate(rows)

    assert result.passed is True
    assert result.scarcity_triggered is False
    assert result.n_active_rows == 5
    assert result.n_independent_units == 5


def test_tnfa_gate_counts_source_families_not_rows() -> None:
    rows = pd.DataFrame(
        {
            "status": ["active"] * 5,
            "paper_id": ["P1"] * 5,
            "source_family_id": ["F1"] * 5,
        }
    )

    result = evaluate_tnfa_gate(rows)

    assert result.passed is False
    assert result.n_active_rows == 5
    assert result.n_independent_units == 1
    assert "independent_units_below_5" in result.reasons


def test_tnfa_gate_triggers_scarcity_after_twenty_candidates() -> None:
    rows = pd.DataFrame(
        {
            "status": ["screened_not_activated"] * 20,
            "paper_or_dataset": [f"GSE{i}" for i in range(20)],
        }
    )

    result = evaluate_tnfa_gate(rows)

    assert result.passed is False
    assert result.scarcity_triggered is True
    assert result.n_candidates == 20
    assert result.n_active_rows == 0


def test_tnfa_gate_accepts_iteration11_eligibility_status_column() -> None:
    rows = pd.DataFrame(
        {
            "eligibility_status": ["deferred"] * 40,
            "paper_id": [f"P{i}" for i in range(40)],
            "source_family_id": [f"F{i}" for i in range(40)],
        }
    )

    result = evaluate_tnfa_gate(rows, scarcity_candidate_threshold=40)

    assert result.passed is False
    assert result.scarcity_triggered is True
    assert result.n_candidates == 40
    assert result.n_active_rows == 0


def test_independent_source_units_falls_back_to_ledger_family_fields() -> None:
    rows = pd.DataFrame(
        {
            "status": ["active", "active"],
            "doi_or_accession": ["GSE1", "GSE2"],
            "paper_or_dataset": ["Dataset 1", "Dataset 2"],
        }
    )

    assert independent_source_units(rows) == {"GSE1::Dataset 1", "GSE2::Dataset 2"}
