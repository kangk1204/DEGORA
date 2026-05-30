from __future__ import annotations

import sqlite3

import pandas as pd
import pytest

from degora.resource_package import (
    SOURCE_FAMILY_COLUMNS,
    TOP_GENE_COLUMNS,
    capture_gene_api_json,
    export_score_resource_package,
    source_family_collapse_table,
    top_gene_resource_table,
    validate_columns,
)
from degora.score_db import write_score_database


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
            "species": ["Homo sapiens"] * 8,
            "cell_system": ["A", "A", "B", "A", "A", "B", "A", "B"],
            "hypoxia_modality": ["1% O2"] * 8,
            "duration_h": ["24"] * 8,
            "source_path": ["source.csv"] * 8,
            "source_url": ["https://example.test"] * 8,
        }
    )


def test_export_score_resource_package_writes_tables_and_sidecars(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    summary = export_score_resource_package(
        tmp_path / "degora_gene_scores.csv",
        tmp_path / "degora_scores.db",
        tmp_path / "resource",
        top_n=2,
    )

    top_path = tmp_path / "resource" / "degora_score_resource_top_genes.tsv"
    collapse_path = tmp_path / "resource" / "degora_score_source_family_collapse_top_genes.tsv"
    assert top_path.exists()
    assert collapse_path.exists()
    assert top_path.with_suffix(".tsv.source").exists()
    assert collapse_path.with_suffix(".tsv.source").exists()
    assert summary["top_genes"] == ["VEGFA", "HK2"]

    top = pd.read_csv(top_path, sep="\t")
    collapse = pd.read_csv(collapse_path, sep="\t")
    assert list(top.columns) == TOP_GENE_COLUMNS
    assert list(collapse.columns) == SOURCE_FAMILY_COLUMNS

    vegfa = collapse.loc[collapse["gene_symbol"].eq("VEGFA")].iloc[0]
    assert vegfa["n_studies"] == 3
    assert vegfa["n_source_units"] == 2
    assert vegfa["raw_evidence_row_count"] == 3
    assert vegfa["support_label"] == "2 / 2 source units"


def test_source_family_collapse_uses_gene_evidence_rows(tmp_path) -> None:
    score_csv = tmp_path / "scores.csv"
    pd.DataFrame(
        [
            {
                "rank_label": "#1 / 1",
                "gene_symbol": "VEGFA",
                "evidence_tier": "A",
                "degora_score": 90.0,
                "priority_score": 88.0,
                "priority_top_percent": 1.0,
                "evidence_reliability_score": 91.0,
                "direction_confidence_index": 0.95,
                "loo_rank_stability_score": 0.90,
                "top_percent_label": "top 1.00%",
                "support_label": "2 / 2 source units",
                "direction_label": "100.0% up-concordant",
                "weighted_lfc": 2.1,
                "high_confidence": True,
            }
        ]
    ).to_csv(score_csv, index=False)
    db_path = tmp_path / "scores.db"
    with sqlite3.connect(db_path) as connection:
        pd.DataFrame(
            {
                "gene_symbol": ["VEGFA", "VEGFA", "VEGFA"],
                "study_id": ["S1", "S2", "S3"],
                "source_unit_id": ["P1", "P1", "P2"],
            }
        ).to_sql("gene_evidence", connection, index=False)

    collapse = source_family_collapse_table(score_csv, db_path, top_n=1)

    assert collapse.loc[0, "n_studies"] == 3
    assert collapse.loc[0, "n_source_units"] == 2
    assert collapse.loc[0, "raw_evidence_row_count"] == 3
    assert collapse.loc[0, "source_units"] == "P1;P2"


def test_source_family_collapse_handles_empty_top_scores(tmp_path) -> None:
    score_csv = tmp_path / "scores.csv"
    pd.DataFrame(columns=TOP_GENE_COLUMNS).to_csv(score_csv, index=False)
    db_path = tmp_path / "scores.db"
    with sqlite3.connect(db_path) as connection:
        pd.DataFrame({"gene_symbol": [], "study_id": [], "source_unit_id": []}).to_sql(
            "gene_evidence",
            connection,
            index=False,
        )

    collapse = source_family_collapse_table(score_csv, db_path, top_n=5)

    assert collapse.empty
    assert collapse.columns.tolist() == SOURCE_FAMILY_COLUMNS


def test_resource_package_column_validation() -> None:
    with pytest.raises(ValueError, match="missing required columns"):
        validate_columns(pd.DataFrame({"gene_symbol": ["VEGFA"]}), TOP_GENE_COLUMNS, label="bad table")


def test_top_gene_resource_table_requires_display_columns(tmp_path) -> None:
    path = tmp_path / "scores.csv"
    pd.DataFrame({"gene_symbol": ["VEGFA"]}).to_csv(path, index=False)

    with pytest.raises(ValueError, match="score CSV is missing"):
        top_gene_resource_table(path)


def test_resource_package_rejects_invalid_top_n_and_nonlocal_api(tmp_path) -> None:
    path = tmp_path / "scores.csv"
    pd.DataFrame(
        [
            {
                "rank_label": "#1 / 1",
                "gene_symbol": "VEGFA",
                "evidence_tier": "A",
                "degora_score": 90.0,
                "priority_score": 88.0,
                "priority_top_percent": 1.0,
                "evidence_reliability_score": 91.0,
                "direction_confidence_index": 0.95,
                "loo_rank_stability_score": 0.90,
                "top_percent_label": "top 1.00%",
                "support_label": "1 / 1 source units",
                "direction_label": "100.0% up-concordant",
                "weighted_lfc": 2.1,
                "high_confidence": True,
            }
        ]
    ).to_csv(path, index=False)

    with pytest.raises(ValueError, match="top_n"):
        top_gene_resource_table(path, top_n=0)
    with pytest.raises(ValueError, match="top_n"):
        source_family_collapse_table(path, tmp_path / "missing.db", top_n=0)
    with pytest.raises(ValueError, match="local DEGORA API"):
        capture_gene_api_json("file:///etc/passwd", tmp_path / "out.json", db_path=tmp_path / "missing.db", command="cmd")
