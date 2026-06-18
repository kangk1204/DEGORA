from __future__ import annotations

import pandas as pd
import pytest

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


def test_harmonize_repairs_excel_date_mangled_gene_symbols() -> None:
    frame = pd.DataFrame(
        {
            "gene": [
                "6-Mar",
                "Sep-2",
                "2025-03-10 00:00:00",
                "2007-09-01 00:00:00",
                "2024-12-01 00:00:00",
                "DEC2",
                "VEGFA",
            ],
            "log2FoldChange": [1.0, -1.0, 0.5, -0.5, 0.25, -0.25, 2.0],
            "pvalue": [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 1e-8],
        }
    )

    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), {"study_id": "DATES"})

    assert set(out["gene_symbol"]) == {
        "BHLHE40",
        "BHLHE41",
        "MARCHF6",
        "MARCHF10",
        "SEPTIN2",
        "SEPTIN7",
        "VEGFA",
    }


def test_harmonize_rejects_pvalues_outside_unit_interval() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["BAD_NEG", "BAD_GT1", "OK_ZERO"],
            "log2FoldChange": [2.0, -1.5, 1.0],
            "pvalue": [-0.3, 2.0, 0.0],
        }
    )

    with pytest.raises(ValueError, match=r"BAD_P: p-value column 'pvalue'.*outside \[0, 1\]"):
        harmonize_frame(
            frame,
            TableMapping("gene", "log2FoldChange", "pvalue"),
            {"study_id": "BAD_P", "paper_id": "PAPER001"},
        )


def test_harmonize_rejects_adjusted_pvalues_outside_unit_interval() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["BAD_NEG", "BAD_GT1", "OK_ZERO"],
            "log2FoldChange": [2.0, -1.5, 1.0],
            "pvalue": [0.03, 0.02, 0.0],
            "padj": [-0.1, 2.0, 0.0],
        }
    )

    with pytest.raises(ValueError, match=r"BAD_Q: adjusted p-value/FDR column 'padj'.*outside \[0, 1\]"):
        harmonize_frame(
            frame,
            TableMapping("gene", "log2FoldChange", "pvalue", "padj"),
            {"study_id": "BAD_Q", "paper_id": "PAPER001"},
        )


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


def test_harmonize_does_not_warn_on_gene_level_rnaseq_duplicate_collapse() -> None:
    # Regression: a blank/NaN probe_collapse cell must be treated as "not declared" (never the
    # string "nan"), and gene-level (RNA-seq) duplicate-symbol collapse must not emit a per-source
    # warning, while the rule actually applied is still recorded in metadata for audit.
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX", "RPL13A"],
            "log2FoldChange": [2.5, -3.0, 0.1],
            "pvalue": [1e-5, 1e-4, 0.9],
            "padj": [0.001, 0.01, 0.9],
        }
    )
    meta = {
        "study_id": "RNASEQ_DUP",
        "paper_id": "P1",
        "pipeline": "DESeq2",
        "assay_type": "RNA-seq",
        "probe_collapse": float("nan"),
    }
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), meta)
    genex = out.loc[out["gene_symbol"].eq("GENEX")].iloc[0]
    assert genex["gene_symbol_collapse_rule"] == "min_pvalue_max_abs_lfc"  # applied rule still recorded
    assert genex["requested_probe_collapse"] == ""  # NaN normalized, never the literal string "nan"
    assert genex["gene_symbol_collapse_warning"] == ""  # no RNA-seq probe-collapse noise


def test_harmonize_still_warns_for_microarray_undeclared_probe_collapse() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX", "RPL13A"],
            "logFC": [2.5, -3.0, 0.1],
            "P.Value": [1e-5, 1e-4, 0.9],
            "adj.P.Val": [0.001, 0.01, 0.9],
        }
    )
    meta = {
        "study_id": "MICRO_DUP",
        "paper_id": "P1",
        "pipeline": "limma_microarray",
        "assay_type": "microarray",
        "probe_collapse": float("nan"),
    }
    out = harmonize_frame(frame, TableMapping("gene", "logFC", "P.Value", "adj.P.Val"), meta)
    genex = out.loc[out["gene_symbol"].eq("GENEX")].iloc[0]
    assert "duplicate gene symbols were collapsed" in genex["gene_symbol_collapse_warning"]


def test_harmonize_emits_explicit_source_unit_id_for_scoring() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "HK2"],
            "log2FoldChange": [2.0, -1.5],
            "pvalue": [1e-8, 1e-4],
        }
    )
    meta = {
        "study_id": "S1",
        "paper_id": "PAPER",
        "source_unit_id": "UNIT_A",
        "pipeline": "DESeq2",
    }

    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), meta)

    assert (out["source_unit_id"] == "UNIT_A").all()
    assert (out["paper_id"] == "PAPER").all()


def test_harmonize_flags_probe_collapse_mismatch_without_silent_best_probe() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX", "RPL13A"],
            "logFC": [0.2, 2.5, 0.1],
            "P.Value": [0.2, 1e-5, 0.9],
            "adj.P.Val": [0.5, 0.001, 0.9],
        }
    )
    meta = {
        "study_id": "MICRO_MEDIAN",
        "paper_id": "P1",
        "assay_type": "microarray",
        "probe_collapse": "median_expression",
        "table_scope": "full_results",
    }

    out = harmonize_frame(frame, TableMapping("gene", "logFC", "P.Value", "adj.P.Val"), meta)

    genex = out.loc[out["gene_symbol"].eq("GENEX")].iloc[0]
    # The actually-applied collapse is best-probe, and the config asked for something else;
    # both must be recorded and the mismatch must be a non-empty (non-silent) warning.
    assert genex["gene_symbol_collapse_rule"] == "min_pvalue_max_abs_lfc"
    assert genex["requested_probe_collapse"] == "median_expression"
    warning = genex["gene_symbol_collapse_warning"]
    assert "median_expression" in warning
    assert "min_pvalue_max_abs_lfc" in warning


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
