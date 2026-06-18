from __future__ import annotations

import pandas as pd

from degora.excel_export import _cap_evidence_for_sheet, _evidence_row_cap


def test_evidence_sheet_cap_keeps_top_ranked_genes_and_flags() -> None:
    # 60 evidence rows across 6 genes; cap to 20 must keep the lowest-rank (best) genes' evidence.
    evidence = pd.DataFrame(
        {
            "gene_symbol": [f"G{i % 6}" for i in range(60)],
            "study_id": [f"S{i}" for i in range(60)],
        }
    )
    genes = pd.DataFrame(
        {
            "gene_symbol": [f"G{i}" for i in range(6)],
            "quality_weighted_degora_rank": [3, 1, 5, 2, 6, 4],
        }
    )
    capped, was_capped = _cap_evidence_for_sheet(evidence, genes, cap=20)
    assert was_capped is True
    assert len(capped) == 20
    # The two best-ranked genes are G1 (rank 1) and G3 (rank 2); their evidence must survive the cap.
    assert {"G1", "G3"} <= set(capped["gene_symbol"])
    # The worst-ranked gene G4 (rank 6) must be dropped first.
    assert "G4" not in set(capped["gene_symbol"])


def test_evidence_sheet_cap_noop_when_small_or_disabled() -> None:
    evidence = pd.DataFrame({"gene_symbol": ["G0", "G1"], "study_id": ["S0", "S1"]})
    genes = pd.DataFrame({"gene_symbol": ["G0", "G1"], "quality_weighted_degora_rank": [1, 2]})
    sheet, was_capped = _cap_evidence_for_sheet(evidence, genes, cap=100)
    assert was_capped is False and len(sheet) == 2
    sheet0, capped0 = _cap_evidence_for_sheet(evidence, genes, cap=0)
    assert capped0 is False and len(sheet0) == 2


def test_evidence_row_cap_reads_env(monkeypatch) -> None:
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "5000")
    assert _evidence_row_cap() == 5000
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "0")
    assert _evidence_row_cap() == 0  # 0 disables the cap
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "not-an-int")
    assert _evidence_row_cap() > 0  # invalid value falls back to the default
