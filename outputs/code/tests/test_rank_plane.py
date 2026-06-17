from __future__ import annotations

import pandas as pd

from degora.rank_plane import rank_plane_gene_summary, rank_plane_points, rank_plane_study_summary


def test_rank_plane_points_use_pvalue_and_signed_effect_ranks() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S1", "S2", "S2", "S2"],
            "gene_symbol": ["VEGFA", "HK2", "RPL13A", "VEGFA", "HK2", "RPL13A"],
            "pvalue": [1e-8, 0.02, 0.9, 1e-7, 0.03, 0.8],
            "lfc": [2.0, -3.0, 0.1, 1.5, -2.5, 0.2],
        }
    )

    points = rank_plane_points(harmonized)
    s1_vegfa = points.loc[points["study_id"].eq("S1") & points["gene_symbol"].eq("VEGFA")].iloc[0]
    s1_hk2 = points.loc[points["study_id"].eq("S1") & points["gene_symbol"].eq("HK2")].iloc[0]

    assert s1_vegfa["p_rank_strength"] == 1.0
    assert s1_hk2["effect_rank_strength"] == 1.0
    assert s1_hk2["signed_effect_rank"] == -1.0
    assert points["n_genes_in_rank_universe"].nunique() == 1


def test_rank_plane_gene_and_study_summaries_are_deterministic() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S1", "S2", "S2", "S2"],
            "gene_symbol": ["VEGFA", "HK2", "RPL13A", "VEGFA", "HK2", "RPL13A"],
            "pvalue": [1e-8, 0.02, 0.9, 1e-7, 0.03, 0.8],
            "lfc": [3.0, -2.0, 0.1, 2.5, -1.5, 0.2],
        }
    )

    points = rank_plane_points(harmonized)
    genes = rank_plane_gene_summary(points, joint_threshold=0.75)
    studies = rank_plane_study_summary(points, joint_threshold=0.75)

    assert genes.iloc[0]["gene_symbol"] in {"HK2", "VEGFA"}
    assert genes.loc[genes["gene_symbol"].eq("VEGFA"), "effect_sign_concordance"].iloc[0] == 1.0
    assert studies["study_id"].tolist() == ["S1", "S2"]
    assert studies["joint_high_count"].min() >= 1


def test_rank_plane_gene_summary_excludes_single_study_genes_by_default() -> None:
    points = pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S2"],
            "gene_symbol": ["SINGLE", "SHARED", "SHARED"],
            "p_rank_strength": [1.0, 0.9, 0.9],
            "effect_rank_strength": [1.0, 0.9, 0.9],
            "signed_effect_rank": [1.0, 0.9, 0.9],
            "rank_plane_delta": [0.0, 0.0, 0.0],
        }
    )

    genes = rank_plane_gene_summary(points, joint_threshold=0.75)

    assert genes["gene_symbol"].tolist() == ["SHARED"]


def test_rank_plane_reports_declared_rank_universe() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S1"],
            "gene_symbol": ["ISG15", "IFIT1"],
            "pvalue": [1e-8, 1e-7],
            "lfc": [4.0, 3.0],
            "n_genes_in_study": [20_000, 20_000],
        }
    )

    points = rank_plane_points(harmonized)

    assert points["n_genes_in_rank_universe"].unique().tolist() == [20_000]
    ifit1 = points.loc[points["gene_symbol"].eq("IFIT1")].iloc[0]
    assert ifit1["p_rank_strength"] == 1.0 - (1.0 / 19_999.0)
    assert ifit1["effect_rank_strength"] == 1.0 - (1.0 / 19_999.0)
