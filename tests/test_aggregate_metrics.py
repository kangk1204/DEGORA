from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from degora.aggregate import (
    collapse_gene_source_units,
    rank_product_consensus,
    slice_consensus,
    stouffer_consensus,
)
from degora.metrics import recall_at_k


def test_slice_consensus_and_recall() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "paper_id": ["P1", "P1", "P2", "P2", "P3", "P3"],
            "gene_symbol": ["VEGFA", "RPL13A", "VEGFA", "RPL13A", "HK2", "VEGFA"],
            "signed_z": [5.5, 0.1, 4.0, 0.2, 3.0, 3.5],
            "lfc": [2.0, 0.1, 1.7, 0.1, 1.2, 1.5],
            "n_ctrl": [3, 3, 4, 4, 5, 5],
            "n_treat": [3, 3, 4, 4, 5, 5],
            "normalized_rank": [0.1, 0.9, 0.1, 0.8, 0.2, 0.1],
        }
    )

    consensus = slice_consensus(harmonized, min_studies=2)
    assert consensus.iloc[0]["gene_symbol"] == "VEGFA"
    assert consensus.iloc[0]["stouffer_padj"] < 0.001

    metrics = recall_at_k(consensus, {"VEGFA", "HK2"}, 1)
    assert metrics["recall"] == 0.5


def test_slice_consensus_uses_source_units_not_contrast_rows() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["P1_T1", "P1_T2", "P1_T1", "P2_T1"],
            "paper_id": ["P1", "P1", "P1", "P2"],
            "gene_symbol": ["SAME_SOURCE_ONLY", "SAME_SOURCE_ONLY", "CROSS_SOURCE", "CROSS_SOURCE"],
            "signed_z": [8.0, 7.0, 3.0, 2.5],
            "lfc": [3.0, 2.5, 1.2, 1.0],
            "n_ctrl": [3, 3, 3, 3],
            "n_treat": [3, 3, 3, 3],
            "normalized_rank": [0.001, 0.002, 0.02, 0.03],
        }
    )

    consensus = slice_consensus(harmonized, min_studies=2)

    assert consensus["gene_symbol"].tolist() == ["CROSS_SOURCE"]
    assert consensus.iloc[0]["n_studies"] == 2
    assert consensus.iloc[0]["n_studies_rank"] == 2


def test_source_unit_collapse_averages_instead_of_selecting_max_abs_z() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": [f"P1_T{i}" for i in range(1, 5)] + ["P2_T1"],
            "paper_id": ["P1"] * 4 + ["P2"],
            "gene_symbol": ["NOISE"] * 5,
            "signed_z": [0.2, -0.1, 3.0, -2.5, 0.3],
            "lfc": [0.1, -0.1, 2.0, -1.8, 0.2],
            "n_ctrl": [3] * 5,
            "n_treat": [3] * 5,
            "normalized_rank": [0.8, 0.9, 0.01, 0.02, 0.7],
        }
    )

    collapsed = collapse_gene_source_units(harmonized)
    p1 = collapsed.loc[collapsed["source_unit_id"].eq("P1")].iloc[0]

    assert np.isclose(p1["signed_z"], 0.15)
    assert not np.isclose(abs(p1["signed_z"]), 3.0)
    assert np.isclose(p1["normalized_rank"], 0.4325)


def test_stouffer_weighted_lfc_ignores_missing_lfc_denominator() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "gene_symbol": ["GENE1", "GENE1"],
            "signed_z": [3.0, 3.0],
            "lfc": [2.0, np.nan],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "normalized_rank": [0.01, 0.02],
        }
    )

    consensus = stouffer_consensus(harmonized, min_studies=2)

    assert consensus.loc[0, "gene_symbol"] == "GENE1"
    assert consensus.loc[0, "weighted_lfc"] == 2.0


def test_stouffer_reports_source_unit_heterogeneity() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3"],
            "paper_id": ["P1", "P2", "P3"],
            "gene_symbol": ["GENE1", "GENE1", "GENE1"],
            "signed_z": [4.0, 3.5, -2.0],
            "lfc": [1.5, 1.2, -1.0],
            "n_ctrl": [3, 3, 3],
            "n_treat": [3, 3, 3],
            "normalized_rank": [0.01, 0.02, 0.05],
        }
    )

    consensus = stouffer_consensus(harmonized, min_studies=2)

    assert consensus.loc[0, "heterogeneity_q"] > 0
    assert consensus.loc[0, "heterogeneity_df"] == 2
    assert 0 <= consensus.loc[0, "heterogeneity_i2"] <= 1


def test_stouffer_log_probability_preserves_order_after_linear_underflow() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S1", "S2"],
            "paper_id": ["P1", "P2", "P1", "P2"],
            "gene_symbol": ["A_WEAKER", "A_WEAKER", "Z_STRONGER", "Z_STRONGER"],
            "signed_z": [37.0, 37.0, 38.0, 38.0],
            "lfc": [2.0, 2.0, 2.0, 2.0],
            "n_ctrl": [3, 3, 3, 3],
            "n_treat": [3, 3, 3, 3],
            "normalized_rank": [0.01, 0.01, 0.01, 0.01],
        }
    )

    consensus = stouffer_consensus(harmonized, min_studies=2)

    assert (consensus["stouffer_p"] == 0.0).all()
    assert consensus["gene_symbol"].tolist() == ["Z_STRONGER", "A_WEAKER"]
    assert consensus.loc[0, "stouffer_neglog10_p"] > consensus.loc[1, "stouffer_neglog10_p"]


def test_stouffer_exact_null_uses_positive_zero_neglog10() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "gene_symbol": ["NULL_GENE", "NULL_GENE"],
            "signed_z": [0.0, 0.0],
            "lfc": [0.0, 0.0],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "normalized_rank": [0.5, 0.5],
        }
    )

    consensus = stouffer_consensus(harmonized, min_studies=2)
    value = consensus.loc[0, "stouffer_neglog10_p"]

    assert value == 0.0
    assert not np.signbit(value)


def test_rank_product_exact_null_uses_positive_zero_score() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "gene_symbol": ["NULL_GENE", "NULL_GENE"],
            "signed_z": [0.0, 0.0],
            "lfc": [0.0, 0.0],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "normalized_rank": [1.0, 1.0],
        }
    )

    consensus = rank_product_consensus(harmonized, min_studies=2)
    value = consensus.loc[0, "rank_score"]

    assert value == 0.0
    assert not np.signbit(value)


def test_source_unit_sample_size_weight_is_capped() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S_BIG", "S_SMALL"],
            "paper_id": ["P_BIG", "P_SMALL"],
            "gene_symbol": ["GENE", "GENE"],
            "signed_z": [4.0, 4.0],
            "lfc": [1.0, 1.0],
            "n_ctrl": [60, 3],
            "n_treat": [60, 3],
            "normalized_rank": [0.01, 0.02],
        }
    )

    collapsed = collapse_gene_source_units(harmonized).set_index("source_unit_id")

    assert collapsed.loc["P_BIG", "weight"] == 4.0
    assert np.isclose(collapsed.loc["P_SMALL", "weight"], np.sqrt(6.0))


def test_source_unit_collapse_representative_label_is_deterministic() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["P1_T2", "P1_T1", "P2_T1"],
            "paper_id": ["P1", "P1", "P2"],
            "gene_symbol": ["GENE", "GENE", "GENE"],
            "signed_z": [2.0, 4.0, 3.0],
            "lfc": [1.0, 2.0, 1.5],
            "n_ctrl": [3, 3, 3],
            "n_treat": [3, 3, 3],
            "normalized_rank": [0.2, 0.1, 0.3],
        }
    )

    first = collapse_gene_source_units(harmonized)
    shuffled = collapse_gene_source_units(harmonized.sample(frac=1, random_state=7))

    assert first.sort_values("source_unit_id")["study_id"].tolist() == shuffled.sort_values("source_unit_id")[
        "study_id"
    ].tolist()
    assert first.loc[first["source_unit_id"].eq("P1"), "study_id"].iloc[0] == "P1_T1"


def test_direction_concordance_is_strength_weighted_across_source_units() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "gene_symbol": ["GENE", "GENE"],
            "signed_z": [8.0, -0.1],
            "lfc": [3.0, -0.01],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "normalized_rank": [0.001, 0.9],
        }
    )

    consensus = slice_consensus(harmonized, min_studies=2)

    assert consensus.iloc[0]["sign_concordance"] > 0.98
    assert not np.isclose(consensus.iloc[0]["sign_concordance"], 0.5)


def test_source_unit_collapse_supports_predeclared_time_course_modes() -> None:
    base = pd.DataFrame(
        {
            "study_id": ["T1", "T2", "T3"],
            "paper_id": ["P1", "P1", "P1"],
            "gene_symbol": ["GENE", "GENE", "GENE"],
            "signed_z": [1.0, 3.0, 5.0],
            "lfc": [0.5, 1.5, 2.5],
            "n_ctrl": [3, 3, 3],
            "n_treat": [3, 3, 3],
            "normalized_rank": [0.9, 0.5, 0.1],
            "duration_h": ["1", "8", "24"],
        }
    )

    early = collapse_gene_source_units(base.assign(time_course_mode="early"))
    late = collapse_gene_source_units(base.assign(time_course_mode="late"))
    peak = collapse_gene_source_units(base.assign(time_course_mode="peak_mean"))

    assert np.isclose(early.iloc[0]["signed_z"], 1.0)
    assert np.isclose(late.iloc[0]["signed_z"], 5.0)
    assert np.isclose(peak.iloc[0]["signed_z"], 4.0)
    assert peak.iloc[0]["n_contrast_rows"] == 2


def test_peak_mean_selects_on_statistical_strength_not_effect_size() -> None:
    """Pin what "peak" means, because the two readings disagree.

    peak_mean keeps the strongest half of a source unit's contrasts by |signed_z|,
    and signed_z comes from the p-value. A reader of a time course is at least as
    likely to mean the largest fold change, so the selection rule is documented in
    the collapse rule, the catalog help, the template and the workbook dictionary -
    and asserted here so it cannot drift away from all four.
    """

    from scipy.stats import norm

    from degora.aggregate import SOURCE_UNIT_COLLAPSE_RULE, collapse_gene_source_units

    def row(study_id: str, lfc: float, pvalue: float, duration: str) -> dict:
        return {
            "study_id": study_id,
            "paper_id": "U1",
            "source_unit_id": "U1",
            "gene_symbol": "G1",
            "lfc": lfc,
            "signed_z": float(np.sign(lfc) * norm.isf(pvalue / 2.0)),
            "normalized_rank": 0.01,
            "n_genes_in_study": 1000,
            "n_ctrl": 3,
            "n_treat": 3,
            "duration_h": duration,
            "time_course_mode": "peak_mean",
        }

    # The largest fold change and the strongest p-value are deliberately different rows.
    harmonized = pd.DataFrame(
        [
            row("U1_02h", 0.5, 0.05, "2"),
            row("U1_08h", 1.0, 1e-12, "8"),
            row("U1_24h", 4.0, 0.04, "24"),
        ]
    )

    collapsed = collapse_gene_source_units(harmonized)

    assert len(collapsed) == 1
    # Two of three kept: the p=1e-12 row and the stronger of the remaining two.
    assert int(collapsed["n_contrast_rows"].iloc[0]) == 2
    # Effect-size selection would have pooled 4.0 alone; strength selection pools 1.0 and 4.0.
    assert collapsed["lfc"].iloc[0] == pytest.approx(2.5)
    assert "rather than effect size" in SOURCE_UNIT_COLLAPSE_RULE
