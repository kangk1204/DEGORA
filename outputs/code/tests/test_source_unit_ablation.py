from __future__ import annotations

import json

import pandas as pd

from scripts.write_source_unit_ablation import write_ablation


def test_source_unit_ablation_reports_legacy_inflation_and_sensitivity_exclusion(tmp_path) -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["HYP003", "S1_T1", "S1_T2", "S2_T1", "S1_T1", "S2_T1", "S1_T1", "S1_T2"],
            "paper_id": ["P0", "P1", "P1", "P2", "P1", "P2", "P1", "P1"],
            "gene_symbol": [
                "VEGFA",
                "ALDOA",
                "ALDOA",
                "ALDOA",
                "RPL13A",
                "RPL13A",
                "SAME_SOURCE_ONLY",
                "SAME_SOURCE_ONLY",
            ],
            "signed_z": [10.0, 7.0, 6.0, 5.0, 0.1, 0.2, 8.0, 7.0],
            "lfc": [3.0, 2.0, 1.8, 1.5, 0.1, 0.1, 2.5, 2.4],
            "pvalue": [1e-20, 1e-12, 1e-10, 1e-8, 0.8, 0.7, 1e-14, 1e-13],
            "padj": [1e-18, 1e-10, 1e-8, 1e-6, 0.9, 0.8, 1e-12, 1e-11],
            "normalized_rank": [0.001, 0.01, 0.02, 0.03, 0.9, 0.8, 0.002, 0.003],
            "n_ctrl": [3] * 8,
            "n_treat": [3] * 8,
        }
    )
    harmonized_path = tmp_path / "harmonized.csv"
    harmonized.to_csv(harmonized_path, index=False)
    catalog = pd.DataFrame(
        {
            "study_id": ["HYP003", "S1_T1", "S1_T2", "S2_T1"],
            "paper_id": ["P0", "P1", "P1", "P2"],
            "source_path": ["x.csv"] * 4,
            "gene_column": ["gene"] * 4,
            "lfc_column": ["lfc"] * 4,
            "p_column": ["pvalue"] * 4,
            "include_in_analysis": [True] * 4,
            "notes": ["provisional sensitivity evidence only", "", "", ""],
        }
    )
    catalog_path = tmp_path / "catalog.csv"
    catalog.to_csv(catalog_path, index=False)

    summary = write_ablation(
        harmonized_path,
        tmp_path / "out",
        catalog_path=catalog_path,
        min_studies=2,
        explicit_exclude=[],
        command="unit-test",
    )

    assert summary["excluded_sensitivity_study_ids"] == ["HYP003"]
    assert (tmp_path / "out" / "source_unit_ablation_summary.tsv").exists()
    assert (tmp_path / "out" / "source_unit_ablation_summary.json.source").exists()
    persisted = json.loads((tmp_path / "out" / "source_unit_ablation_summary.json").read_text())
    variants = {row["variant"]: row for row in persisted["variant_summaries"]}
    assert variants["contrast_level_legacy"]["n_scored_genes"] == 3
    assert variants["source_unit_v1_2_mean"]["n_scored_genes"] == 2
    assert "source_unit_best_signal_sensitivity" in variants
    assert "source_unit_best_signal_sensitivity_diagnostic_only" in persisted["score_versions_compared"]
