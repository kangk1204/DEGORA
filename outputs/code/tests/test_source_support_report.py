from __future__ import annotations

import pandas as pd

from scripts.write_source_support_report import SUMMARY_COLUMNS, TOP_GENE_COLUMNS, build_report


def test_source_support_report_counts_direction_and_source_support(tmp_path) -> None:
    score_path = tmp_path / "scores.csv"
    pd.DataFrame(
        {
            "degora_rank": [1, 2, 3, 4],
            "gene_symbol": ["IFIT1", "CMPK2", "MX1", "ACTB"],
            "consensus_direction": ["up", "up", "down", "up"],
            "n_source_units": [2, 2, 2, 1],
            "sign_concordance": [1.0, 1.0, 1.0, 1.0],
            "degora_score": [90.0, 88.0, 80.0, 10.0],
            "evidence_reliability_score": [95.0, 94.0, 93.0, 20.0],
        }
    ).to_csv(score_path, index=False)
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["IFIT1", "MX1"], "expected_direction": ["up", "up"]}).to_csv(
        gold_path,
        index=False,
    )
    marker_path = tmp_path / "markers.csv"
    pd.DataFrame(
        {
            "gene_symbol": ["CMPK2"],
            "expected_direction": ["up"],
            "marker_role": ["interpretive_isg"],
            "validation_role": ["post_output_interpretive_only"],
            "evidence_basis": ["known IFN marker"],
            "source_url": ["https://example.test"],
            "notes": ["not locked"],
        }
    ).to_csv(marker_path, index=False)

    summary, top_genes, metadata = build_report(
        score_path,
        gold_path,
        marker_panel=marker_path,
        top_n=3,
        cutoffs=[2, 3],
        min_source_units=2,
        min_sign_concordance=1.0,
    )

    assert summary.columns.tolist() == SUMMARY_COLUMNS
    assert top_genes.columns.tolist() == TOP_GENE_COLUMNS
    locked_k3 = summary.loc[
        summary["metric"].eq("locked_membership_recall") & summary["cutoff"].eq(3)
    ].iloc[0]
    direction_k3 = summary.loc[
        summary["metric"].eq("locked_direction_recall") & summary["cutoff"].eq(3)
    ].iloc[0]
    support_k3 = summary.loc[
        summary["metric"].eq("locked_source_supported_recall") & summary["cutoff"].eq(3)
    ].iloc[0]
    known_k3 = summary.loc[
        summary["metric"].eq("known_or_locked_marker_source_supported") & summary["cutoff"].eq(3)
    ].iloc[0]

    assert locked_k3["n_recovered"] == 2
    assert direction_k3["n_recovered"] == 1
    assert support_k3["n_recovered"] == 1
    assert known_k3["n_recovered"] == 2
    cmpk2 = top_genes.loc[top_genes["gene_symbol"].eq("CMPK2")].iloc[0]
    assert bool(cmpk2["interpretive_marker"]) is True
    assert bool(cmpk2["source_support_pass"]) is True
    assert metadata["n_interpretive_markers"] == 1


def test_source_support_report_handles_missing_optional_score_fields(tmp_path) -> None:
    score_path = tmp_path / "scores.csv"
    pd.DataFrame(
        {
            "degora_rank": [1, 2],
            "gene_symbol": ["HSPA5", "DDIT3"],
            "consensus_direction": ["up", "up"],
            "n_source_units": [2, 1],
            "sign_concordance": [1.0, 1.0],
        }
    ).to_csv(score_path, index=False)
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["HSPA5", "DDIT3"], "expected_direction": ["up", "up"]}).to_csv(
        gold_path,
        index=False,
    )

    summary, top_genes, _ = build_report(score_path, gold_path, top_n=2, cutoffs=[2], min_source_units=2)

    support = summary.loc[summary["metric"].eq("locked_source_supported_recall")].iloc[0]
    assert support["n_recovered"] == 1
    assert top_genes.loc[top_genes["gene_symbol"].eq("HSPA5"), "priority_score"].iloc[0] == ""
