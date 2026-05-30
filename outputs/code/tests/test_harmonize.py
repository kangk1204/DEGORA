from __future__ import annotations

import pandas as pd

from degora.harmonize import TableMapping, assess_table_scope, harmonize_frame
from degora.slice_runner import apply_gene_type_filter, catalog_include_mask


def test_harmonize_frame_computes_signed_z_and_ranks() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "RPL13A", "HK2"],
            "log2FoldChange": [2.0, 0.1, -1.5],
            "pvalue": [1e-8, 0.5, 1e-4],
            "padj": [1e-6, 0.8, 1e-3],
        }
    )
    meta = {
        "study_id": "TEST001",
        "paper_id": "PAPER001",
        "pipeline": "DESeq2",
        "n_ctrl": 3,
        "n_treat": 3,
    }

    out = harmonize_frame(
        frame,
        TableMapping("gene", "log2FoldChange", "pvalue", "padj"),
        meta,
    )

    assert list(out["gene_symbol"]) == ["VEGFA", "HK2", "RPL13A"]
    assert out.loc[out["gene_symbol"].eq("VEGFA"), "signed_z"].iloc[0] > 0
    assert out.loc[out["gene_symbol"].eq("HK2"), "signed_z"].iloc[0] < 0
    assert out["within_study_rank"].min() == 1.0
    assert out["n_genes_in_study"].nunique() == 1
    assert out["table_scope"].iloc[0] in {"ambiguous", "full_results"}


def test_deg_only_scope_uses_declared_rank_universe() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["ISG15", "IFIT1", "MX1"],
            "log2FoldChange": [5.0, 4.0, 3.0],
            "pvalue": [1e-20, 1e-12, 1e-8],
            "padj": [1e-18, 1e-10, 1e-6],
        }
    )
    meta = {
        "study_id": "DEG_ONLY",
        "paper_id": "PAPER001",
        "pipeline": "DESeq2",
        "table_scope": "deg_only",
        "rank_universe_size": 20_000,
    }

    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), meta)

    assert out["table_scope"].unique().tolist() == ["deg_only"]
    assert out["n_genes_in_study"].unique().tolist() == [20_000]
    assert out.loc[out["within_study_rank"].eq(1.0), "normalized_rank"].iloc[0] == 1 / 20_000
    assert "missing genes are unreported" in out["rank_universe_warning"].iloc[0]


def test_harmonize_preserves_assay_metadata_for_microarray_rows() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["RBM39", "RPL13A"],
            "logFC": [1.5, 0.0],
            "P.Value": [0.001, 0.9],
            "adj.P.Val": [0.01, 0.9],
        }
    )
    meta = {
        "study_id": "MICRO001",
        "paper_id": "GSE_MICRO",
        "pipeline": "limma_microarray",
        "assay_type": "microarray",
        "source_input_type": "limma_full_table",
        "platform": "GPL570",
        "normalization": "RMA/log2",
        "probe_collapse": "author_gene_level",
        "table_scope": "full_results",
    }

    out = harmonize_frame(frame, TableMapping("gene", "logFC", "P.Value", "adj.P.Val"), meta)

    rbm39 = out.loc[out["gene_symbol"].eq("RBM39")].iloc[0]
    assert rbm39["assay_type"] == "microarray"
    assert rbm39["source_input_type"] == "limma_full_table"
    assert rbm39["platform"] == "GPL570"
    assert rbm39["probe_collapse"] == "author_gene_level"


def test_harmonize_collapses_duplicate_gene_symbols_before_ranking() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX", "GENEX", "RPL13A"],
            "probe": ["p1", "p2", "p3", "p4"],
            "logFC": [0.2, 2.5, -3.0, 0.1],
            "P.Value": [0.2, 1e-5, 1e-4, 0.9],
            "adj.P.Val": [0.5, 0.001, 0.01, 0.9],
        }
    )
    meta = {
        "study_id": "MICRO_PROBES",
        "paper_id": "P1",
        "pipeline": "limma_microarray",
        "assay_type": "microarray",
        "source_input_type": "limma_full_table",
        "platform": "GPL570",
        "normalization": "RMA/log2",
        "table_scope": "full_results",
    }

    out = harmonize_frame(frame, TableMapping("gene", "logFC", "P.Value", "adj.P.Val"), meta)

    assert out["gene_symbol"].tolist() == ["GENEX", "RPL13A"]
    genex = out.loc[out["gene_symbol"].eq("GENEX")].iloc[0]
    assert genex["lfc"] == 2.5
    assert genex["n_source_rows_for_gene"] == 3
    assert genex["gene_symbol_collapse_rule"] == "min_pvalue_max_abs_lfc"
    assert "duplicate gene symbols were collapsed" in genex["gene_symbol_collapse_warning"]
    assert out["n_genes_in_study"].unique().tolist() == [2]


def test_scope_assessment_detects_full_result_tables() -> None:
    frame = pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(200)],
            "log2FoldChange": [0.1] * 200,
            "pvalue": [0.001] * 20 + [0.8] * 180,
            "padj": [0.01] * 20 + [1.0] * 180,
        }
    )

    scope = assess_table_scope(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"))

    assert scope["effective_scope"] == "full_results"
    assert scope["assessment"] == "full_results_likely"


def test_apply_gene_type_filter_keeps_requested_biotype() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "PSEUDO1", "HK2"],
            "Gene.type": ["protein_coding", "processed_pseudogene", "protein_coding"],
        }
    )

    filtered, summary = apply_gene_type_filter(frame, "Gene.type", "protein_coding")

    assert list(filtered["gene"]) == ["VEGFA", "HK2"]
    assert summary["applied"] is True
    assert summary["rows_before"] == 3
    assert summary["rows_after"] == 2


def test_catalog_include_mask_defaults_empty_values_to_active() -> None:
    catalog = pd.DataFrame({"include_in_analysis": ["true", "no", "", None]})

    mask = catalog_include_mask(catalog)

    assert mask.tolist() == [True, False, True, True]
