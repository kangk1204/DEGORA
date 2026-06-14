from __future__ import annotations

import pandas as pd

from degora.diagnostics import sensitivity_metrics, target_support_by_study


def test_target_support_by_study_includes_missing_targets() -> None:
    catalog = pd.DataFrame(
        {
            "study_id": ["S1"],
            "paper_id": ["P1"],
            "pipeline": ["DESeq2"],
            "cell_system": ["cells"],
            "hypoxia_modality": ["1% O2"],
            "duration_h": [8],
            "p_column": ["pvalue"],
            "padj_column": ["padj"],
            "notes": [""],
        }
    )
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1"],
            "gene_symbol": ["VEGFA"],
            "lfc": [2.0],
            "pvalue": [0.001],
            "padj": [0.01],
            "signed_z": [3.2],
            "within_study_rank": [1.0],
            "normalized_rank": [0.01],
            "n_genes_in_study": [100],
        }
    )

    out = target_support_by_study(harmonized, catalog)

    assert len(out) >= 20
    assert bool(out.loc[out["gene_symbol"].eq("VEGFA"), "target_present"].iloc[0]) is True
    assert bool(out.loc[out["gene_symbol"].eq("EPO"), "target_present"].iloc[0]) is False


def test_sensitivity_metrics_reports_predeclared_variants() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["HYP001", "HYP001", "HYP006", "HYP006", "HYP007", "HYP007", "HYP008", "HYP008"],
            "paper_id": ["P1", "P1", "Bauer", "Bauer", "Bauer", "Bauer", "Bentley", "Bentley"],
            "gene_symbol": ["VEGFA", "RPL13A", "VEGFA", "RPL13A", "VEGFA", "RPL13A", "VEGFA", "RPL13A"],
            "signed_z": [4.0, 0.1, 3.0, 0.2, 3.5, 0.3, 4.2, 0.1],
            "lfc": [1.5, 0.1, 1.3, 0.1, 1.4, 0.1, 1.6, 0.1],
            "n_ctrl": [3] * 8,
            "n_treat": [3] * 8,
            "normalized_rank": [0.1, 0.9, 0.1, 0.8, 0.1, 0.7, 0.1, 0.9],
        }
    )

    out = sensitivity_metrics(harmonized)

    assert "active_all" in set(out["variant"])
    assert "collapse_by_paper_id" not in set(out["variant"])
    assert out.loc[out["variant"].eq("active_all"), "n_recovered_at_50"].iloc[0] >= 1
