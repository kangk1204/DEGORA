from __future__ import annotations

import json

import pandas as pd

from scripts.write_microarray_integration_report import build_report


def test_microarray_integration_report_frames_nutriomics_db_use_case(tmp_path) -> None:
    prep = {
        "microarray_direct_excluded_drug": False,
        "primary_excluded_drug_result_replaced": False,
        "as_published_author_deg_validation": False,
        "microarray_interpretation": "sensitivity evidence only",
        "e7820_summary": {
            "platform": "GPL",
            "normalization": "processed",
            "source_input_type": "normalized_expression_matrix",
            "collapse_rule": "min_pvalue_max_abs_lfc",
            "n_gene_rows": 10,
            "n_collapsed_probe_rows": 2,
        },
        "sirbm39_summary": {
            "n_gene_rows": 10,
            "n_collapsed_probe_rows": 2,
        },
    }
    sensitivity = {"unused": True}
    comparison = {
        "primary": {"anchor_recovery_top_100": {"anchor_recovery": 0.125}},
        "microarray_sensitivity": {"anchor_recovery_top_100": {"anchor_recovery": 0.0625}},
    }
    score_db = {
        "n_gene_scores": 10,
        "n_source_units": 3,
        "n_contrasts": 5,
        "n_evidence_rows": 20,
        "top_genes": ["A", "B"],
    }
    comparator = pd.DataFrame(
        {
            "method_id": ["degora_deg_score", "rank_product_approx", "robustrankaggreg"],
            "setting_id": ["v1", "default", "default"],
            "run_status": ["ok", "ok", "ok"],
            "n_rows": [10, 10, 10],
            "recall_at_10": [0, 0, 0],
            "recall_at_20": [0, 0, 0],
            "recall_at_50": [0, 0, 0],
            "recall_at_100": [0.0625, 0.125, 0.125],
            "recovered_at_100": ["TYMS", "RBM39;TYMS", "RBM39;TYMS"],
            "missing_at_100": ["RBM39", "", ""],
            "top10": ["A", "B", "C"],
            "failure_mode": ["", "", ""],
        }
    )
    anchors = pd.DataFrame(
        {
            "gene_symbol": ["RBM39", "TYMS"],
            "degora_rank": [766, 85],
            "degora_score": [72.9, 79.0],
            "top_percent_label": ["top 1.91%", "top 0.212%"],
            "support_label": ["3 / 3 source units", "3 / 3 source units"],
            "direction_label": ["100.0% up-concordant", "100.0% down-concordant"],
            "n_source_units": [3, 3],
            "n_contrasts_observed": [5, 5],
            "weighted_lfc": [0.66, -0.81],
        }
    )

    paths = {}
    for name, payload in [
        ("preparation_summary", prep),
        ("sensitivity_summary", sensitivity),
        ("comparison_summary", comparison),
        ("score_db_summary", score_db),
    ]:
        path = tmp_path / f"{name}.json"
        path.write_text(json.dumps(payload))
        paths[name] = path
    comparator_path = tmp_path / "comparator.csv"
    comparator.to_csv(comparator_path, index=False)
    anchors_path = tmp_path / "anchors.tsv"
    anchors.to_csv(anchors_path, sep="\t", index=False)

    report = build_report(
        preparation_summary=paths["preparation_summary"],
        sensitivity_summary=paths["sensitivity_summary"],
        comparison_summary=paths["comparison_summary"],
        comparator_summary=comparator_path,
        score_db_summary=paths["score_db_summary"],
        anchor_ranks=anchors_path,
    )

    assert "NutriOmics-style evidence database" in report["nutriomics_db_implication"]
    assert "direction_concordance_percent" in ";".join(report["nutrient_gene_query_contract"])
    assert report["microarray_preparation"]["collapse_rule"] == "min_pvalue_max_abs_lfc"
    assert report["anchor_evidence"]["degora_score_top100_recovery"] == 0.0625
