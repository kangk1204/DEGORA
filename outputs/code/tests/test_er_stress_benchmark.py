from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.write_er_stress_benchmark import (
    ContrastSpec,
    _clean_gse84450_cuffdiff,
    _derive_logcpm_welch_contrast,
    _derive_normalized_expression_contrast,
    _read_human_ensembl_counts,
)


def test_read_human_ensembl_counts_maps_hgnc_and_collapses_duplicates(tmp_path: Path) -> None:
    counts_path = tmp_path / "counts.csv"
    pd.DataFrame(
        {
            "ctrl1": [10, 3, 5],
            "ctrl2": [11, 3, 5],
            "trt1": [100, 7, 5],
            "trt2": [110, 8, 5],
        },
        index=["ENSG000001.1", "ENSG000001.2", "ENSG000002"],
    ).to_csv(counts_path)

    counts, summary = _read_human_ensembl_counts(counts_path, {"ENSG000001": "HSPA5"})

    assert list(counts.index) == ["HSPA5"]
    assert counts.loc["HSPA5", "ctrl1"] == 13
    assert summary["n_unmapped_genes"] == 1
    assert summary["n_gene_symbols_after_collapse"] == 1


def test_logcpm_welch_contrast_returns_positive_lfc_for_treatment_increase(tmp_path: Path) -> None:
    counts = pd.DataFrame(
        {
            "ctrl1": [10, 100, 0],
            "ctrl2": [11, 100, 1],
            "trt1": [100, 100, 0],
            "trt2": [110, 100, 1],
        },
        index=["HSPA5", "ACTB", "LOWCOUNT"],
    )
    spec = ContrastSpec(
        study_id="TEST",
        paper_id="P1",
        source_url="https://example.test",
        source_path=tmp_path / "out.csv",
        raw_source_key="test",
        species="Homo sapiens",
        cell_system="cell",
        condition="tunicamycin",
        duration_h="6",
        n_ctrl=2,
        n_treat=2,
        control_columns=("ctrl1", "ctrl2"),
        treat_columns=("trt1", "trt2"),
        output_name="out.csv",
        pipeline="logCPM_Welch_derived_from_public_counts",
        assay_type="RNA-seq",
        source_input_type="derived_count_table",
    )

    deg = _derive_logcpm_welch_contrast(counts, spec)

    hspa5 = deg.loc[deg["gene_symbol"].eq("HSPA5")].iloc[0]
    assert hspa5["log2FoldChange"] > 0
    assert 0 <= hspa5["pvalue"] <= 1
    assert 0 <= hspa5["padj"] <= 1
    assert "LOWCOUNT" not in set(deg["gene_symbol"])
    assert hspa5["n_genes_before_low_count_filter"] == 3
    assert hspa5["n_genes_after_low_count_filter"] == 2


def test_clean_gse84450_filters_non_ok_rows_and_collapses_gene_symbols(tmp_path: Path) -> None:
    path = tmp_path / "gse84450.diff.tsv"
    pd.DataFrame(
        [
            ["row1", "g1", "Ddit3", "loc", "WT_UT", "WT_DTT", "NOTEST", 1, 0, float("-inf"), 0, 1, 1, "no"],
            ["row2", "g2", "Ddit3", "loc", "WT_UT", "WT_DTT", "OK", 1, 10, 2.0, 5, 0.01, 0.02, "yes"],
            ["row3", "g3", "Ddit3", "loc", "WT_UT", "WT_DTT", "OK", 1, 3, 1.0, 3, 0.02, 0.04, "yes"],
        ],
        columns=[
            "test_id",
            "gene_id",
            "gene",
            "locus",
            "sample_1",
            "sample_2",
            "status",
            "value_1",
            "value_2",
            "log2(fold_change)",
            "test_stat",
            "p_value",
            "q_value",
            "significant",
        ],
    ).to_csv(path, sep="\t", index=False)

    clean, summary = _clean_gse84450_cuffdiff(path)

    assert clean["gene_symbol"].tolist() == ["DDIT3"]
    assert clean.iloc[0]["log2FoldChange"] == 2.0
    assert summary["raw_rows"] == 3
    assert summary["status_ok_rows"] == 2
    assert summary["n_gene_symbols_after_collapse"] == 1


def test_normalized_expression_contrast_collapses_duplicate_symbols(tmp_path: Path) -> None:
    matrix = tmp_path / "norm.tsv"
    pd.DataFrame(
        {
            "refSeqID": ["NM_1", "NM_2", "NM_3"],
            "name": ["Ddit3", "Ddit3", "Actb"],
            "length": [100, 200, 300],
            "TE.DMSO1": [1.0, 1.0, 5.0],
            "TE.DMSO2": [1.1, 1.1, 5.0],
            "TE.THAP1": [5.0, 2.0, 5.0],
            "TE.THAP2": [5.1, 2.1, 5.0],
        }
    ).to_csv(matrix, sep="\t", index=False)
    spec = ContrastSpec(
        study_id="TEST",
        paper_id="P1",
        source_url="https://example.test",
        source_path=tmp_path / "out.csv",
        raw_source_key="test",
        species="Mus musculus",
        cell_system="cell",
        condition="thapsigargin",
        duration_h="6",
        n_ctrl=2,
        n_treat=2,
        control_columns=("TE.DMSO1", "TE.DMSO2"),
        treat_columns=("TE.THAP1", "TE.THAP2"),
        output_name="out.csv",
        pipeline="welch_normalized_expression_matrix",
        assay_type="RNA-seq",
        source_input_type="normalized_expression_matrix",
    )

    deg, summary = _derive_normalized_expression_contrast(matrix, spec)

    assert deg["gene_symbol"].tolist().count("DDIT3") == 1
    assert deg.loc[deg["gene_symbol"].eq("DDIT3"), "log2FoldChange"].iloc[0] > 0
    assert summary["n_collapsed_duplicate_symbol_rows"] == 1
