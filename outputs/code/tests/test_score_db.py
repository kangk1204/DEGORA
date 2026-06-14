from __future__ import annotations

import json
import sqlite3

import numpy as np
import pandas as pd
import pytest

from degora.score_db import _write_sqlite, degora_score_table, write_score_database


def test_write_sqlite_preserves_existing_db_on_failed_rebuild(tmp_path) -> None:
    db = tmp_path / "degora_scores.db"
    genes = pd.DataFrame({"gene_symbol": ["A", "B"], "degora_rank": [1, 2], "degora_score": [0.9, 0.8]})
    evidence = pd.DataFrame({"gene_symbol": ["A", "B"], "study_id": ["S1", "S1"]})
    studies = pd.DataFrame({"source_unit_id": ["S1"]})

    _write_sqlite(db, genes, evidence, studies, {"corpus": "test"})
    assert db.exists()
    original = db.read_bytes()

    # Duplicate gene_symbol violates the unique index, so the rebuild fails mid-write.
    bad_genes = pd.DataFrame({"gene_symbol": ["A", "A"], "degora_rank": [1, 2], "degora_score": [0.9, 0.8]})
    with pytest.raises(Exception):
        _write_sqlite(db, bad_genes, evidence, studies, {"corpus": "test"})

    # The previous good DB must survive, with no lingering temp file.
    assert db.read_bytes() == original
    assert not (tmp_path / "degora_scores.db.tmp").exists()
    with sqlite3.connect(db) as connection:
        assert connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0] == 2


def _harmonized() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3", "S1", "S2", "S3", "S1", "S3"],
            "paper_id": ["P1", "P1", "P2", "P1", "P1", "P2", "P1", "P2"],
            "gene_symbol": ["VEGFA", "VEGFA", "VEGFA", "RPL13A", "RPL13A", "RPL13A", "HK2", "HK2"],
            "lfc": [2.4, 2.0, 1.8, 0.2, -0.2, 0.1, 2.1, 1.9],
            "signed_z": [6.0, 5.0, 4.5, 0.2, -0.2, 0.1, 4.0, 3.8],
            "pvalue": [1e-8, 1e-7, 1e-6, 0.8, 0.7, 0.9, 1e-5, 2e-5],
            "padj": [1e-6, 1e-5, 1e-4, 0.9, 0.8, 0.9, 1e-3, 1e-3],
            "normalized_rank": [0.02, 0.03, 0.04, 0.8, 0.7, 0.9, 0.06, 0.05],
            "n_ctrl": [3, 3, 4, 3, 3, 4, 3, 4],
            "n_treat": [3, 3, 4, 3, 3, 4, 3, 4],
            "n_genes_in_study": [1000] * 8,
            "pipeline": ["DESeq2", "DESeq2", "edgeR", "DESeq2", "DESeq2", "edgeR", "DESeq2", "edgeR"],
            "assay_type": ["RNA-seq", "RNA-seq", "microarray", "RNA-seq", "RNA-seq", "microarray", "RNA-seq", "microarray"],
            "source_input_type": ["author_deg_table"] * 8,
            "platform": [""] * 8,
            "normalization": [""] * 8,
            "probe_collapse": [""] * 8,
            "species": ["Homo sapiens"] * 8,
            "cell_system": ["A", "A", "B", "A", "A", "B", "A", "B"],
            "hypoxia_modality": ["1% O2"] * 8,
            "duration_h": ["24"] * 8,
            "source_path": ["source.csv"] * 8,
            "source_url": ["https://example.test"] * 8,
        }
    )


def test_degora_score_prioritizes_repeated_directional_source_unit_support() -> None:
    scores, evidence, metadata = degora_score_table(_harmonized(), min_studies=2)

    assert scores.iloc[0]["gene_symbol"] == "VEGFA"
    vegfa = scores.loc[scores["gene_symbol"].eq("VEGFA")].iloc[0]
    hk2 = scores.loc[scores["gene_symbol"].eq("HK2")].iloc[0]
    rpl13a = scores.loc[scores["gene_symbol"].eq("RPL13A")].iloc[0]

    assert vegfa["n_source_units"] == 2
    assert vegfa["n_contrasts_observed"] == 3
    assert vegfa["sign_concordance"] == 1.0
    assert vegfa["rank_label"] == "#1 / 3"
    assert vegfa["top_percent"] == 33.333333
    assert vegfa["percentile"] == 100.0
    assert vegfa["top_percent_label"] == "top 33.33%"
    assert vegfa["support_label"] == "2 / 2 source units"
    assert vegfa["direction_label"] == "100.0% up-concordant"
    assert vegfa["quality_weighted_degora_rank"] >= 1
    assert vegfa["quality_weighted_degora_score"] > 0
    assert vegfa["source_quality_support_score"] > 0
    assert vegfa["priority_score"] > 0
    assert vegfa["priority_rank"] >= 1
    assert vegfa["evidence_reliability_score"] > 0
    assert vegfa["direction_confidence_index"] > 0.5
    assert vegfa["direction_concordant_source_units"] == 2
    assert vegfa["direction_total_source_units"] == 2
    assert vegfa["direction_posterior_mean"] == 0.75
    assert vegfa["loo_rank_stability_score"] > 0
    assert "heterogeneity_i2" in scores.columns
    assert "heterogeneity_flag" in scores.columns
    assert "re_stouffer_z" in scores.columns
    assert "rra_rho" in scores.columns
    assert "effect_meta_log2fc_re" in scores.columns
    assert 0 <= vegfa["heterogeneity_i2"] <= 1
    assert vegfa["heterogeneity_flag"] in {"low_or_unestimated", "moderate_context_review", "high_context_dependent_review"}
    assert vegfa["re_stouffer_shrinkage_factor"] >= 1.0
    assert 0 <= vegfa["rra_rho"] <= 1
    assert vegfa["rra_rank"] >= 1
    assert vegfa["effect_meta_k"] == 2
    assert vegfa["effect_meta_ci_low"] < vegfa["effect_meta_log2fc_re"] < vegfa["effect_meta_ci_high"]
    assert vegfa["degora_score"] > hk2["degora_score"] > rpl13a["degora_score"]
    assert evidence.loc[evidence["gene_symbol"].eq("VEGFA"), "source_unit_id"].nunique() == 2
    assert set(evidence.loc[evidence["gene_symbol"].eq("VEGFA"), "assay_type"]) == {"RNA-seq", "microarray"}
    assert set(evidence["source_quality_label"]) == {"high"}
    assert len(evidence.loc[evidence["gene_symbol"].eq("VEGFA")]) == 2
    p1 = evidence.loc[evidence["gene_symbol"].eq("VEGFA") & evidence["source_unit_id"].eq("P1")].iloc[0]
    assert np.isclose(p1["signed_z"], 5.5)
    assert np.isnan(p1["aggregate_padj"])
    assert p1["min_source_padj"] == 1e-6
    assert p1["contributing_study_ids"] == "S1;S2"
    assert set(scores["evidence_tier"]) <= {"A", "B", "C", "D"}
    assert metadata["score_version"] == "degora_score_v1_2_source_unit_mean"
    assert metadata["independent_unit_for_consensus"].startswith("source_unit_id")
    assert metadata["n_contrasts_total"] == 3
    assert metadata["n_source_units_total"] == 2
    assert "max-|z|" in metadata["source_unit_collapse_rule"]
    assert "stouffer_padj" in metadata["high_confidence_rule"]
    assert "evidence_tier_rules" in metadata
    assert "quality_weighted_score_formula" in metadata
    assert "priority_score_weights" in metadata
    assert "evidence_reliability_score_weights" in metadata
    assert "direction_confidence_rule" in metadata
    assert "random_effects_stouffer_rule" in metadata
    assert "rra_rule" in metadata
    assert "effect_meta_rule" in metadata
    assert "heterogeneity_rule" in metadata
    assert "heterogeneity_flag_rule" in metadata
    assert "loo_stability_rule" in metadata
    assert "source_quality_diagnostics" in metadata


def test_degora_score_component_golden_values() -> None:
    scores, _, _ = degora_score_table(_harmonized(), min_studies=2)
    vegfa = scores.loc[scores["gene_symbol"].eq("VEGFA")].iloc[0]

    assert vegfa["support_score"] == pytest.approx(1.0, abs=1e-6)
    assert vegfa["direction_score"] == pytest.approx(1.0, abs=1e-6)
    assert vegfa["evidence_score"] == pytest.approx(0.583257, abs=1e-6)
    assert vegfa["rank_score_component"] == pytest.approx(0.968377, abs=1e-6)
    assert vegfa["effect_score"] == pytest.approx(0.629470, abs=1e-6)
    assert vegfa["degora_score"] == pytest.approx(85.305316, abs=1e-6)
    assert vegfa["quality_weighted_degora_score"] == pytest.approx(85.298205, abs=1e-6)
    assert vegfa["source_quality_support_score"] == pytest.approx(1.0, abs=1e-6)
    assert vegfa["quality_weighted_sign_concordance"] == pytest.approx(1.0, abs=1e-6)


def test_effect_meta_layer_skips_numerically_unusable_inverse_variance_weights() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A", "B"],
            "paper_id": ["A", "B"],
            "gene_symbol": ["WEAK", "WEAK"],
            "lfc": [0.1, -0.1],
            "signed_z": [1e-154, -1e-154],
            "pvalue": [1.0, 1.0],
            "padj": [1.0, 1.0],
            "normalized_rank": [0.5, 0.5],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    row = scores.loc[scores["gene_symbol"].eq("WEAK")].iloc[0]
    effect_columns = [
        "effect_meta_log2fc_re",
        "effect_meta_se",
        "effect_meta_ci_low",
        "effect_meta_ci_high",
        "effect_meta_tau2",
        "effect_meta_i2",
    ]

    assert row["effect_meta_k"] == 0
    assert not np.isinf(row[effect_columns].astype(float)).any()


def test_degora_score_does_not_treat_same_source_contrasts_as_independent() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["P1_T1", "P1_T2", "P1_T1", "P2_T1"],
            "paper_id": ["P1", "P1", "P1", "P2"],
            "gene_symbol": ["SAME_SOURCE_ONLY", "SAME_SOURCE_ONLY", "CROSS_SOURCE", "CROSS_SOURCE"],
            "lfc": [4.0, 3.5, 1.1, 1.0],
            "signed_z": [8.0, 7.0, 2.0, 2.2],
            "pvalue": [1e-15, 1e-12, 0.01, 0.02],
            "padj": [1e-12, 1e-10, 0.05, 0.05],
            "normalized_rank": [0.001, 0.002, 0.05, 0.06],
            "n_ctrl": [3, 3, 3, 3],
            "n_treat": [3, 3, 3, 3],
            "n_genes_in_study": [1000] * 4,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)

    assert scores["gene_symbol"].tolist() == ["CROSS_SOURCE"]


def test_study_gene_evidence_preserves_mixed_contrast_provenance() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1"],
            "lfc": [1.0, 2.0],
            "signed_z": [2.0, 4.0],
            "pvalue": [0.05, 0.001],
            "padj": [0.1, 0.01],
            "normalized_rank": [0.2, 0.01],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "pipeline": ["DESeq2", "limma_microarray"],
            "assay_type": ["RNA-seq", "microarray"],
            "source_input_type": ["author_deg_table", "normalized_expression_matrix"],
            "platform": ["", "GPLX"],
            "normalization": ["DESeq2", "RMA/log2"],
            "probe_collapse": ["", "min_pvalue_max_abs_lfc"],
            "species": ["human", "human"],
            "cell_system": ["A", "A"],
            "hypoxia_modality": ["drug 1h", "drug 24h"],
            "duration_h": ["1", "24"],
            "source_path": ["a.csv", "b.csv"],
            "source_url": ["u1", "u2"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["study_id"] == "T1"
    assert row["assay_type"] == "RNA-seq;microarray"
    assert row["platform"] == "GPLX"
    assert row["contributing_study_ids"] == "T1;T2"
    assert row["contributing_assay_types"] == "RNA-seq;microarray"
    assert row["contributing_source_paths"] == "a.csv;b.csv"
    assert row["source_quality_weight"] < 1.0
    assert row["source_quality_label"] == "low"
    assert row["min_source_pvalue"] == 0.001
    assert np.isclose(row["aggregate_pvalue"], 2.0 * 0.0013498980316300933)


def test_study_gene_evidence_metadata_follows_time_course_selection() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2", "T3"],
            "paper_id": ["P1", "P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1", "GENE1"],
            "lfc": [1.0, 8.0, 6.0],
            "signed_z": [2.0, 10.0, 9.0],
            "pvalue": [0.04, 1e-12, 1e-9],
            "padj": [0.08, 1e-10, 1e-7],
            "normalized_rank": [0.2, 0.001, 0.002],
            "n_ctrl": [3, 3, 3],
            "n_treat": [3, 3, 3],
            "n_genes_in_study": [1000, 1000, 1000],
            "duration_h": ["1", "8", "24"],
            "time_course_mode": ["early", "early", "early"],
            "source_path": ["early.csv", "mid.csv", "late.csv"],
            "source_url": ["u1", "u2", "u3"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["study_id"] == "T1"
    assert row["contributing_study_ids"] == "T1"
    assert row["contributing_duration_h"] == "1"
    assert row["contributing_source_paths"] == "early.csv"
    assert row["min_source_pvalue"] == 0.04
    assert row["min_source_padj"] == 0.08
    assert row["lfc"] == 1.0
    assert row["signed_z"] == 2.0


def test_source_quality_uses_conservative_replicate_counts_within_source_unit() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1"],
            "lfc": [1.0, 1.5],
            "signed_z": [2.0, 3.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.05, 0.06],
            "normalized_rank": [0.1, 0.2],
            "n_ctrl": [3, 1],
            "n_treat": [3, 1],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["n_ctrl"] == 1.0
    assert row["n_treat"] == 1.0
    assert row["source_quality_weight"] == 0.5
    assert row["source_quality_label"] == "low"


def test_source_quality_uses_conservative_table_scope_within_source_unit() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1"],
            "lfc": [1.0, 1.5],
            "signed_z": [2.0, 3.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.05, 0.06],
            "normalized_rank": [0.1, 0.2],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "deg_only"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["table_scope"] == "deg_only;full_results"
    assert row["source_quality_weight"] == 0.65
    assert row["source_quality_label"] == "medium"


def test_write_score_database_emits_sqlite_and_sidecars(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)

    summary = write_score_database(
        harmonized_path,
        tmp_path,
        db_path=tmp_path / "degora_scores.db",
        extra_metadata={"derived_count_pilot": True, "claim_scope": "unit_test"},
    )

    db_path = tmp_path / "degora_scores.db"
    assert db_path.exists()
    assert (tmp_path / "degora_gene_scores.csv").exists()
    assert (tmp_path / "degora_gene_scores.csv.source").exists()
    assert (tmp_path / "degora_score_metadata.json").exists()
    assert (tmp_path / "degora_source_quality_diagnostics.tsv").exists()
    assert (tmp_path / "degora_source_quality_diagnostics.tsv.source").exists()
    assert (tmp_path / "degora_source_quality_diagnostics.json").exists()
    assert (tmp_path / "degora_scores.db.source").exists()
    assert summary["n_gene_scores"] == 3

    with sqlite3.connect(db_path) as connection:
        top_gene = connection.execute("SELECT gene_symbol FROM genes ORDER BY degora_rank LIMIT 1").fetchone()[0]
        top_label = connection.execute("SELECT top_percent_label FROM genes ORDER BY degora_rank LIMIT 1").fetchone()[0]
        quality_rank = connection.execute("SELECT quality_weighted_degora_rank FROM genes ORDER BY degora_rank LIMIT 1").fetchone()[0]
        priority_score = connection.execute("SELECT priority_score FROM genes ORDER BY degora_rank LIMIT 1").fetchone()[0]
        reliability_score = connection.execute("SELECT evidence_reliability_score FROM genes ORDER BY degora_rank LIMIT 1").fetchone()[0]
        evidence_rows = connection.execute("SELECT COUNT(*) FROM gene_evidence WHERE gene_symbol = 'VEGFA'").fetchone()[0]
        quality_columns = [row[1] for row in connection.execute("PRAGMA table_info(gene_evidence)").fetchall()]
        metadata = dict(connection.execute("SELECT key, value FROM meta").fetchall())

    assert top_gene == "VEGFA"
    assert top_label.startswith("top ")
    assert quality_rank >= 1
    assert priority_score > 0
    assert reliability_score > 0
    assert evidence_rows == 2
    assert "source_quality_weight" in quality_columns
    assert "source_recommended_weight" in quality_columns
    assert "source_reliability_weight" in quality_columns
    assert json.loads(metadata["score_weights"])["support_score"] == 0.30
    assert "priority_score_weights" in metadata
    assert json.loads(metadata["source_quality_weight_rules"])["source_input_type_weights"]["normalized_expression_matrix"] == 0.35
    assert "evidence_tier_rules" in metadata
    assert metadata["derived_count_pilot"] == "true"
    assert metadata["claim_scope"] == "unit_test"


def test_degora_score_caps_nonfinite_lfc_for_browsing_layer() -> None:
    harmonized = _harmonized()
    harmonized.loc[harmonized["gene_symbol"].eq("VEGFA"), "lfc"] = np.inf

    scores, evidence, metadata = degora_score_table(harmonized, min_studies=2)

    assert metadata["n_nonfinite_lfc_capped_for_score"] == 3
    assert np.isfinite(scores.loc[scores["gene_symbol"].eq("VEGFA"), "weighted_lfc"].iloc[0])
    assert evidence.loc[evidence["gene_symbol"].eq("VEGFA"), "lfc"].max() == 10.0


def test_quality_weighted_secondary_score_downweights_low_quality_discordant_source() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A1", "B1", "C1", "A1", "B1", "C1"],
            "paper_id": ["A", "B", "C", "A", "B", "C"],
            "gene_symbol": ["UPR1", "UPR1", "UPR1", "CTRL", "CTRL", "CTRL"],
            "lfc": [2.0, 2.0, -8.0, 0.2, 0.2, -0.1],
            "signed_z": [4.0, 4.0, -2.0, 0.3, 0.3, -0.2],
            "pvalue": [1e-4, 1e-4, 0.05, 0.8, 0.8, 0.9],
            "padj": [1e-3, 1e-3, 0.2, 0.9, 0.9, 0.9],
            "normalized_rank": [0.01, 0.01, 0.8, 0.8, 0.8, 0.9],
            "n_ctrl": [3, 3, 2, 3, 3, 2],
            "n_treat": [3, 3, 2, 3, 3, 2],
            "n_genes_in_study": [1000] * 6,
            "pipeline": ["author", "derived", "matrix", "author", "derived", "matrix"],
            "assay_type": ["RNA-seq"] * 6,
            "source_input_type": [
                "author_deg_table",
                "derived_count_table",
                "normalized_expression_matrix",
                "author_deg_table",
                "derived_count_table",
                "normalized_expression_matrix",
            ],
            "platform": [""] * 6,
            "normalization": [""] * 6,
            "probe_collapse": [""] * 6,
            "species": ["human"] * 6,
            "cell_system": ["A", "B", "C", "A", "B", "C"],
            "hypoxia_modality": ["ER stress"] * 6,
            "duration_h": ["24"] * 6,
            "table_scope": ["full_results"] * 6,
            "source_path": ["source.csv"] * 6,
            "source_url": ["https://example.test"] * 6,
        }
    )

    scores, evidence, metadata = degora_score_table(harmonized, min_studies=2)
    upr1 = scores.loc[scores["gene_symbol"].eq("UPR1")].iloc[0]

    # UPR1 is the canonical z-vs-LFC discordant gene: two higher-quality sources are
    # strongly up (signed_z +4) while one low-quality matrix source is strongly down in
    # effect size (lfc -8) but only mildly down in z (-2). The weighted combined signed-z
    # is positive, so sign_concordance is measured against "up" (the two up sources are
    # the concordant ones). consensus_direction must report that same combined-z direction
    # ("up"); deriving it from weighted_lfc (negative here) made the direction_label
    # self-contradictory ("X% down-concordant"). Regression guard for the score_db P1
    # direction fix.
    assert upr1["stouffer_z"] > 0
    assert upr1["consensus_direction"] == "up"
    assert upr1["consensus_direction"] == ("up" if upr1["stouffer_z"] > 0 else "down")
    assert upr1["quality_weighted_consensus_direction"] == "up"
    assert upr1["quality_weighted_degora_score"] > 0
    low_quality = evidence.loc[evidence["source_unit_id"].eq("C")].iloc[0]
    assert low_quality["source_quality_label"] == "low"
    assert low_quality["source_quality_weight"] < 0.4
    assert low_quality["source_reliability_weight"] < 0.65
    assert metadata["n_source_quality_outliers"] == 0


def test_direction_confidence_penalizes_discordant_sources() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A", "B", "A", "B"],
            "paper_id": ["A", "B", "A", "B"],
            "gene_symbol": ["CONSISTENT", "CONSISTENT", "CONFLICT", "CONFLICT"],
            "lfc": [2.0, 2.0, 2.0, -2.0],
            "signed_z": [4.0, 4.0, 4.0, -4.0],
            "pvalue": [1e-4, 1e-4, 1e-4, 1e-4],
            "padj": [1e-3] * 4,
            "normalized_rank": [0.01, 0.01, 0.02, 0.02],
            "n_ctrl": [3] * 4,
            "n_treat": [3] * 4,
            "n_genes_in_study": [1000] * 4,
            "source_input_type": ["author_deg_table"] * 4,
            "table_scope": ["full_results"] * 4,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    consistent = scores.loc[scores["gene_symbol"].eq("CONSISTENT")].iloc[0]
    conflict = scores.loc[scores["gene_symbol"].eq("CONFLICT")].iloc[0]

    assert consistent["direction_confidence_index"] > conflict["direction_confidence_index"]
    assert consistent["evidence_reliability_score"] > conflict["evidence_reliability_score"]


def test_direction_confidence_counts_sources_against_consensus_direction() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["BIG_DOWN", "SMALL_UP1", "SMALL_UP2"],
            "paper_id": ["BIG", "UP1", "UP2"],
            "gene_symbol": ["GENE", "GENE", "GENE"],
            "lfc": [-1.0, 1.0, 1.0],
            "signed_z": [-3.0, 2.0, 2.0],
            "pvalue": [0.001, 0.01, 0.01],
            "padj": [0.01, 0.05, 0.05],
            "normalized_rank": [0.01, 0.02, 0.02],
            "n_ctrl": [60, 1, 1],
            "n_treat": [60, 1, 1],
            "n_genes_in_study": [1000, 1000, 1000],
            "source_input_type": ["author_deg_table"] * 3,
            "table_scope": ["full_results"] * 3,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    row = scores.iloc[0]

    assert row["stouffer_z"] < 0
    assert row["consensus_direction"] == "down"
    assert row["direction_concordant_source_units"] == 1
    assert row["direction_total_source_units"] == 3
    assert row["direction_confidence_index"] == 0.4


def test_direction_confidence_is_neutral_when_consensus_z_ties() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["UP", "DOWN"],
            "paper_id": ["UP", "DOWN"],
            "gene_symbol": ["GENE", "GENE"],
            "lfc": [1.0, -1.0],
            "signed_z": [2.0, -2.0],
            "pvalue": [0.01, 0.01],
            "padj": [0.05, 0.05],
            "normalized_rank": [0.02, 0.02],
            "n_ctrl": [1, 1],
            "n_treat": [1, 1],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    row = scores.iloc[0]

    assert row["stouffer_z"] == pytest.approx(0.0, abs=1e-12)
    assert row["consensus_direction"] == "flat"
    assert row["direction_concordant_source_units"] == 1
    assert row["direction_total_source_units"] == 2
    assert row["direction_confidence_index"] == 0.5
    assert row["direction_posterior_mean"] == 0.5
