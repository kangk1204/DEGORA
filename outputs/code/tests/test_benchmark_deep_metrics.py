from __future__ import annotations

import sqlite3

import pandas as pd

from figures.make_evidence_card_figure import _source_evidence, parse_case, write_package
from scripts.write_benchmark_deep_metrics import (
    DEGORA_ADVANTAGE_COLUMNS,
    POINT_METRIC_COLUMNS,
    bootstrap_metrics,
    degora_advantage_metrics,
    point_metrics,
)
from scripts.write_study_tool_feasibility_matrix import FEASIBILITY_COLUMNS, build_matrix


def _harmonized() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S2", "S2"],
            "paper_id": ["P1", "P1", "P2", "P2"],
            "source_unit_id": ["P1", "P1", "P2", "P2"],
            "gene_symbol": ["GENEA", "GENEB", "GENEA", "GENEB"],
            "signed_z": [4.0, 0.1, 3.5, -0.2],
            "lfc": [2.0, 0.1, 1.5, -0.1],
            "pvalue": [1e-5, 0.8, 1e-4, 0.7],
            "padj": [1e-4, 0.9, 1e-3, 0.8],
            "normalized_rank": [0.01, 0.8, 0.02, 0.9],
            "n_ctrl": [3, 3, 4, 4],
            "n_treat": [3, 3, 4, 4],
            "pipeline": ["DESeq2", "DESeq2", "limma", "limma"],
            "assay_type": ["RNA-seq", "RNA-seq", "microarray", "microarray"],
            "source_input_type": ["author_deg_table"] * 4,
            "table_scope": ["full_results"] * 4,
            "source_path": ["missing.tsv"] * 4,
        }
    )


def _score(path) -> None:
    pd.DataFrame(
        {
            "degora_rank": [1, 2],
            "gene_symbol": ["GENEA", "GENEB"],
            "consensus_direction": ["up", "down"],
            "quality_weighted_degora_rank": [1, 2],
            "quality_weighted_consensus_direction": ["up", "down"],
            "quality_weighted_top_percent": [1.0, 2.0],
            "n_source_units": [2, 2],
            "sign_concordance": [1.0, 0.5],
            "quality_weighted_sign_concordance": [1.0, 0.5],
            "loo_rank_stability_score": [0.9, 0.4],
            "loo_top50_fraction": [1.0, 0.0],
            "loo_top100_fraction": [1.0, 1.0],
            "source_quality_weight_sum": [2.0, 1.2],
            "source_units": ["P1;P2", "P1;P2"],
        }
    ).to_csv(path, index=False)


def _gold(path) -> None:
    pd.DataFrame({"gene_symbol": ["GENEA"], "expected_direction": ["up"]}).to_csv(path, index=False)


def _baseline_dir(path) -> None:
    path.mkdir()
    result = pd.DataFrame(
        {
            "method_id": ["weighted_stouffer", "weighted_stouffer"],
            "setting_id": ["default", "default"],
            "gene_id": ["GENEA", "GENEB"],
            "symbol": ["GENEA", "GENEB"],
            "rank": [1, 2],
            "score": [10.0, 1.0],
            "pvalue": [1e-6, 0.5],
            "padj": [1e-5, 0.6],
            "effect": [3.0, -0.2],
            "direction": ["up", "down"],
            "n_studies": [2, 2],
            "missingness": [0.0, 0.0],
            "runtime_s": [0.1, 0.1],
            "version": ["test", "test"],
            "status": ["ok", "ok"],
        }
    )
    tsv = path / "test_weighted_stouffer_default.tsv"
    result.to_csv(tsv, sep="\t", index=False)
    pd.DataFrame(
        {
            "artifact": [tsv.name],
            "method_id": ["weighted_stouffer"],
            "setting_id": ["default"],
            "artifact_type": ["baseline_result"],
            "status": ["ok"],
        }
    ).to_csv(path / "baseline_manifest.csv", index=False)


def test_deep_metrics_compute_aurc_bootstrap_and_advantage(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    score_path = tmp_path / "score.csv"
    gold_path = tmp_path / "gold.csv"
    baseline_dir = tmp_path / "baselines"
    _harmonized().to_csv(harmonized_path, index=False)
    _score(score_path)
    _gold(gold_path)
    _baseline_dir(baseline_dir)

    point = point_metrics(
        corpus="toy",
        baseline_dir=baseline_dir,
        gold_path=gold_path,
        degora_score_csv=score_path,
        max_k=5,
    )
    bootstrap = bootstrap_metrics(
        corpus="toy",
        harmonized_path=harmonized_path,
        gold_path=gold_path,
        min_studies=1,
        n_bootstrap=3,
        max_k=5,
        seed=1,
    )
    advantage = degora_advantage_metrics(corpus="toy", degora_score_csv=score_path, gold_path=gold_path)

    assert point.columns.tolist() == POINT_METRIC_COLUMNS
    assert point.loc[point["method_id"].eq("degora_quality_weighted_score"), "aurc_at_max_k"].iloc[0] > 0
    degora_point = point.loc[point["method_id"].eq("degora_quality_weighted_score")].iloc[0]
    assert degora_point["n_recovered_at_10"] == 1
    assert degora_point["recall_at_10_ci_low"] <= degora_point["recall_at_10"] <= degora_point["recall_at_10_ci_high"]
    assert degora_point["precision_at_10"] == 0.5
    assert degora_point["background_auroc"] == 1.0
    assert degora_point["background_auprc_enrichment"] > 1.0
    assert set(bootstrap["method_id"]) == {"fisher", "degora_quality_weighted_score", "weighted_stouffer"}
    assert set(bootstrap["metric"]) == {"aurc_at_max_k", "recall_at_100"}
    assert advantage.columns.tolist() == DEGORA_ADVANTAGE_COLUMNS
    assert advantage.loc[advantage["subset"].eq("locked_gold"), "median_loo_rank_stability_score"].iloc[0] == 0.9


def test_study_tool_feasibility_matrix_marks_summary_and_raw_tool_requirements() -> None:
    matrix = build_matrix(_harmonized(), corpus="toy")

    assert matrix.columns.tolist() == FEASIBILITY_COLUMNS
    assert set(matrix["awfisher"]) == {"compatible"}
    assert matrix["rankprod_exact"].str.startswith("blocked:requires_replicate_expression").all()
    assert matrix["dexma"].str.contains("requires_expression_matrix").all()
    assert matrix["omicc"].str.contains("requires_expression_matrix").all()
    assert matrix["deet"].str.startswith("not_comparator:").all()
    assert matrix["generic_pvalue_combiners"].str.startswith("covered:").all()


def test_study_tool_feasibility_matrix_checks_all_studies_in_source_unit(tmp_path) -> None:
    complete = tmp_path / "S1.tsv"
    incomplete = tmp_path / "S2.tsv"
    pd.DataFrame(
        {
            "id": ["GENEA"],
            "baseMean": [100],
            "log2FoldChange": [1.5],
            "lfcSE": [0.2],
            "stat": [7.5],
            "pvalue": [1e-5],
            "padj": [1e-4],
        }
    ).to_csv(complete, sep="\t", index=False)
    pd.DataFrame(
        {
            "id": ["GENEA"],
            "baseMean": [100],
            "log2FoldChange": [1.5],
            "stat": [7.5],
            "pvalue": [1e-5],
            "padj": [1e-4],
        }
    ).to_csv(incomplete, sep="\t", index=False)
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENEA", "GENEA"],
            "lfc": [1.5, 1.5],
            "signed_z": [4.0, 4.0],
            "pvalue": [1e-5, 1e-5],
            "padj": [1e-4, 1e-4],
            "normalized_rank": [0.01, 0.02],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "pipeline": ["DESeq2", "DESeq2"],
            "assay_type": ["RNA-seq", "RNA-seq"],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
            "source_path": [str(complete), str(incomplete)],
            "source_url": ["file://S1.tsv", "file://S2.tsv"],
        }
    )

    row = build_matrix(harmonized, corpus="toy").iloc[0]

    assert row["study_id"] == "S1;S2"
    assert row["hstouffer"] == "blocked:partial_hstouffer_original_columns"
    assert row["awmeta"] == "blocked:partial_original_variance_or_se"


def test_evidence_card_figure_package_writes_required_outputs(tmp_path) -> None:
    score_path = tmp_path / "score.csv"
    db_path = tmp_path / "scores.db"
    _score(score_path)
    with sqlite3.connect(db_path) as connection:
        pd.DataFrame(
            {
                "gene_symbol": ["GENEA", "GENEA"],
                "source_unit_id": ["P1", "P2"],
                "lfc": [2.0, 1.5],
                "signed_z": [4.0, 3.5],
                "aggregate_pvalue": [1e-5, 1e-4],
                "source_quality_weight": [1.0, 0.9],
                "source_reliability_label": ["high", "high"],
                "direction": ["up", "up"],
            }
        ).to_sql("gene_evidence", connection, index=False)
    case = parse_case(f"toy|{score_path}|{db_path}|GENEA")
    output_dir = tmp_path / "figure"

    summary = write_package([case], output_dir, title="Toy evidence cards", command="test command")

    for key in ["figure_png", "figure_pdf", "figure_svg", "source_data", "legend", "manifest", "validation"]:
        assert summary[key]
        assert (output_dir / summary[key].split("/")[-1]).exists()
    validation = (output_dir / "degora_evidence_cards_validation.txt").read_text()
    assert "degora_evidence_cards_validation.txt\texists=True" in validation
    assert summary["n_gene_rows"] == 1
    assert summary["n_evidence_rows"] == 2


def test_evidence_card_source_evidence_handles_empty_gene_list(tmp_path) -> None:
    db_path = tmp_path / "scores.db"
    with sqlite3.connect(db_path) as connection:
        pd.DataFrame(
            {
                "gene_symbol": ["GENEA"],
                "source_unit_id": ["P1"],
                "lfc": [1.0],
                "signed_z": [2.0],
                "aggregate_pvalue": [0.05],
                "source_quality_weight": [1.0],
                "source_reliability_label": ["high"],
                "direction": ["up"],
            }
        ).to_sql("gene_evidence", connection, index=False)

    evidence = _source_evidence({"corpus": "toy", "db_path": db_path, "genes": []})

    assert evidence.empty
    assert "gene_symbol" in evidence.columns
