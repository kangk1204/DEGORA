from __future__ import annotations

import pandas as pd
import pytest

from degora.harmonize import TableMapping, assess_table_scope, harmonize_frame, read_deg_table
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


def test_harmonize_rejects_reusing_pvalue_as_log_fold_change() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "HK2"],
            "log2FoldChange": [2.0, -1.0],
            "pvalue": [0.01, 0.02],
        }
    )

    with pytest.raises(ValueError, match="lfc_column .* p_column .* both map to source column 'pvalue'"):
        harmonize_frame(
            frame,
            TableMapping("gene", "pvalue", "pvalue"),
            {"study_id": "BAD_MAPPING"},
        )


@pytest.mark.parametrize(
    ("mapping", "left_role", "right_role", "shared_column"),
    [
        (TableMapping("gene", "gene", "pvalue", "padj"), "gene_column", "lfc_column", "gene"),
        (TableMapping("gene", "log2FoldChange", "gene", "padj"), "gene_column", "p_column", "gene"),
        (TableMapping("gene", "log2FoldChange", "pvalue", "gene"), "gene_column", "padj_column", "gene"),
        (
            TableMapping("gene", "log2FoldChange", "pvalue", "log2FoldChange"),
            "lfc_column",
            "padj_column",
            "log2FoldChange",
        ),
    ],
)
def test_harmonize_rejects_every_other_incompatible_mapping_pair(
    mapping: TableMapping,
    left_role: str,
    right_role: str,
    shared_column: str,
) -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "HK2"],
            "log2FoldChange": [2.0, -1.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.02, 0.03],
        }
    )

    with pytest.raises(ValueError) as excinfo:
        harmonize_frame(frame, mapping, {"study_id": "BAD_MAPPING"})

    message = str(excinfo.value)
    assert left_role in message
    assert right_role in message
    assert repr(shared_column) in message


def test_mapping_collision_is_detected_after_legacy_row_label_alias_resolution() -> None:
    frame = pd.DataFrame(
        {
            "row_name": ["VEGFA", "HK2"],
            "log2FoldChange": [2.0, -1.0],
            "pvalue": [0.01, 0.02],
        }
    )

    with pytest.raises(ValueError, match=r"gene_column .* lfc_column .* source column 'row_name'"):
        harmonize_frame(
            frame,
            TableMapping("Unnamed: 0", "row_name", "pvalue"),
            {"study_id": "ROW_LABEL_ALIAS_COLLISION"},
        )


def test_harmonize_allows_same_column_for_pvalue_and_padj() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "HK2"],
            "log2FoldChange": [2.0, -1.0],
            "qvalue": [0.01, 0.02],
        }
    )

    out = harmonize_frame(
        frame,
        TableMapping("gene", "log2FoldChange", "qvalue", "qvalue"),
        {"study_id": "P_EQUALS_PADJ"},
    )

    assert out["pvalue"].tolist() == out["padj"].tolist()


@pytest.mark.parametrize(
    ("mapping", "column_kind", "significant"),
    [
        (TableMapping("gene", "log2FoldChange", "significant"), "p-value", [1, 1, 0]),
        (
            TableMapping("gene", "log2FoldChange", "pvalue", "significant"),
            "adjusted p-value/FDR",
            [1, 1, 0],
        ),
        (TableMapping("gene", "log2FoldChange", "significant"), "p-value", [0, 0, 0]),
    ],
)
def test_harmonize_rejects_binary_significance_flags_as_probability_columns(
    mapping: TableMapping,
    column_kind: str,
    significant: list[int],
) -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "HK2", "RPL13A"],
            "log2FoldChange": [2.0, -1.0, 0.1],
            "pvalue": [0.001, 0.02, 0.4],
            "significant": significant,
        }
    )

    with pytest.raises(ValueError, match=rf"{column_kind} column 'significant'.*only 0 and/or 1"):
        harmonize_frame(frame, mapping, {"study_id": "BINARY_FLAG"})


def test_harmonize_allows_a_legitimate_adjusted_pvalue_column_that_is_all_one() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["VEGFA", "HK2", "RPL13A"],
            "log2FoldChange": [2.0, -1.0, 0.1],
            "pvalue": [0.001, 0.02, 0.4],
            "padj": [1.0, 1.0, 1.0],
        }
    )

    out = harmonize_frame(
        frame,
        TableMapping("gene", "log2FoldChange", "pvalue", "padj"),
        {"study_id": "ALL_ONE_PADJ"},
    )

    assert out["padj"].eq(1.0).all()


def test_harmonize_treats_pvalue_one_as_neutral_and_flags_every_floor_clip() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["NEUTRAL", "FLOORED", "INFORMATIVE"],
            "log2FoldChange": [2.0, 1.0, -1.0],
            "pvalue": [1.0, 1e-310, 0.01],
        }
    )

    out = harmonize_frame(
        frame,
        TableMapping("gene", "log2FoldChange", "pvalue"),
        {"study_id": "P_BOUNDARIES", "paper_id": "P_BOUNDARIES"},
    ).set_index("gene_symbol")

    assert pd.isna(out.loc["NEUTRAL", "signed_z"])
    assert out.loc["NEUTRAL", "pvalue_was_clipped"] == False  # noqa: E712
    assert out.loc["FLOORED", "pvalue"] == pytest.approx(1e-300)
    assert out.loc["FLOORED", "pvalue_was_clipped"] == True  # noqa: E712
    assert out.loc["INFORMATIVE", "signed_z"] < 0


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
        "sign_convention": "treated_vs_control_as_published",
        "table_scope": "full_results",
    }

    out = harmonize_frame(frame, TableMapping("gene", "logFC", "P.Value", "adj.P.Val"), meta)

    rbm39 = out.loc[out["gene_symbol"].eq("RBM39")].iloc[0]
    assert rbm39["assay_type"] == "microarray"
    assert rbm39["source_input_type"] == "limma_full_table"
    assert rbm39["platform"] == "GPL570"
    assert rbm39["probe_collapse"] == "author_gene_level"
    assert rbm39["sign_convention"] == "treated_vs_control_as_published"


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


def test_duplicate_gene_collapse_uses_declared_nominal_p_policy_when_fdr_conflicts() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX"],
            "logFC": [2.0, -4.0],
            "P.Value": [1e-6, 1e-5],
            "adj.P.Val": [0.9, 0.001],
        }
    )
    meta = {
        "study_id": "CONFLICTING_PROBES",
        "paper_id": "P1",
        "assay_type": "microarray",
        "probe_collapse": "min_pvalue_max_abs_lfc",
        "table_scope": "full_results",
    }

    out = harmonize_frame(frame, TableMapping("gene", "logFC", "P.Value", "adj.P.Val"), meta)

    assert len(out) == 1
    selected = out.iloc[0]
    assert selected["lfc"] == 2.0
    assert selected["pvalue"] == pytest.approx(1e-6)
    assert selected["padj"] == pytest.approx(0.9)
    assert selected["gene_symbol_collapse_rule"] == "min_pvalue_max_abs_lfc"


def test_harmonize_does_not_warn_on_gene_level_rnaseq_duplicate_collapse() -> None:
    # Regression: a blank/NaN probe_collapse cell must be treated as "not declared" (never the
    # string "nan"), and gene-level (RNA-seq) duplicate-symbol collapse must not emit a per-source
    # warning, while the rule actually applied is still recorded in metadata for audit.
    frame = pd.DataFrame(
        {
            # Same sign on both duplicate rows: this test is about probe-collapse
            # noise and the NaN probe_collapse cell, not about direction. Duplicate
            # rows that disagree on sign are a separate warning, covered by
            # test_harmonize_warns_when_duplicate_rows_disagree_on_direction.
            "gene": ["GENEX", "GENEX", "RPL13A"],
            "log2FoldChange": [2.5, 3.0, 0.1],
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
    assert "only 200 rows" in scope["reason"]
    assert "truncated supplementary table cannot be ruled out" in scope["reason"]
    assert "confirm table_scope explicitly" in scope["reason"]


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


def test_rows_dropped_for_missing_values_are_counted_and_explained() -> None:
    """Losing rows outright has to warn the way collapsing them already does.

    A table whose effect column exported as text, or whose gene column is half
    empty, lost those rows between the file and the ranking with nothing said on
    the console and nothing in any warning array. The only way to notice was to
    subtract two numbers held in different fields of the metrics file, and even
    then it could not say which column was responsible.
    """

    frame = pd.DataFrame(
        {
            "gene": ["A", "B", None, "D", "E", "F", "G", "H"],
            "log2FoldChange": ["UP", 1.0, 2.0, "#DIV/0!", -1.5, 0.5, 1.2, -0.8],
            "pvalue": [0.01, 0.02, 0.03, 0.04, None, 0.06, 0.07, 0.08],
        }
    )

    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), {"study_id": "LOSSY"})

    assert int(out["n_input_rows"].iloc[0]) == 8
    assert int(out["n_rows_dropped_unusable"].iloc[0]) == 8 - len(out)
    warning = str(out["unusable_row_warning"].iloc[0])
    assert "LOSSY" in warning
    assert "were dropped before ranking" in warning
    # The reason has to name the column, and show what was actually in the cells.
    assert "log2 fold change" in warning
    assert "'UP'" in warning and "'#DIV/0!'" in warning
    assert "gene identifier" in warning
    assert "p-value" in warning


def test_ordinary_missingness_does_not_raise_a_row_loss_warning() -> None:
    """One unusable row in a full table is normal and must stay quiet."""

    frame = pd.DataFrame(
        {
            "gene": [f"G{index}" for index in range(40)],
            "log2FoldChange": [1.0] * 39 + [None],
            "pvalue": [0.01] * 40,
        }
    )

    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), {"study_id": "TIDY"})

    assert int(out["n_rows_dropped_unusable"].iloc[0]) == 1
    assert str(out["unusable_row_warning"].iloc[0]) == ""


def test_row_loss_counts_distinct_rows_not_overlapping_reasons() -> None:
    """A row missing three things is one lost row, not three.

    Summing the per-reason counts reported more rows dropped than the table held
    - one real corpus produced "73,774 of 57,905 rows (127.4%)" - and pushed
    tables under the threshold over it, so a 4% loss warned as 12%.
    """

    size = 1000
    genes: list = [f"G{index}" for index in range(size)]
    effects: list = [1.0] * size
    pvalues: list = [0.01] * size
    for index in range(400):
        genes[index] = None
        effects[index] = None
        pvalues[index] = None

    out = harmonize_frame(
        pd.DataFrame({"gene": genes, "log2FoldChange": effects, "pvalue": pvalues}),
        TableMapping("gene", "log2FoldChange", "pvalue"),
        {"study_id": "OVERLAP"},
    )

    assert int(out["n_rows_dropped_unusable"].iloc[0]) == 400
    warning = str(out["unusable_row_warning"].iloc[0])
    assert "400 of 1,000 rows (40.0%)" in warning
    # The reasons stay as a breakdown, and say they may overlap.
    assert "a row can be missing more than one" in warning


def test_a_loss_under_the_threshold_stays_quiet() -> None:
    """Four rows in a hundred used to warn as twelve percent."""

    size = 100
    genes: list = [f"H{index}" for index in range(size)]
    effects: list = [1.0] * size
    pvalues: list = [0.01] * size
    for index in range(4):
        genes[index] = None
        effects[index] = None
        pvalues[index] = None

    out = harmonize_frame(
        pd.DataFrame({"gene": genes, "log2FoldChange": effects, "pvalue": pvalues}),
        TableMapping("gene", "log2FoldChange", "pvalue"),
        {"study_id": "UNDER"},
    )

    assert int(out["n_rows_dropped_unusable"].iloc[0]) == 4
    assert str(out["unusable_row_warning"].iloc[0]) == ""


def test_duplicate_collapse_is_not_reported_as_unusable_rows() -> None:
    """Probe rows merged into a gene were used, not lost.

    The count was taken from the final row count, after duplicate collapse, so an
    ordinary probe-level table looked as though its source could not supply rows
    it had supplied.
    """

    out = harmonize_frame(
        pd.DataFrame(
            {
                "gene": ["A", "A", "B"],
                "log2FoldChange": [1.0, 2.0, 3.0],
                "pvalue": [0.01, 0.02, 0.03],
            }
        ),
        TableMapping("gene", "log2FoldChange", "pvalue"),
        {"study_id": "DUPES"},
    )

    assert int(out["n_rows_dropped_unusable"].iloc[0]) == 0
    assert int(out["n_rows_merged_by_gene_collapse"].iloc[0]) == 1
    assert str(out["unusable_row_warning"].iloc[0]) == ""


def test_a_gzipped_workbook_is_read_rather_than_decoded_as_text(tmp_path) -> None:
    """Repositories serve supplementary workbooks gzipped, and pandas reads neither.

    A .xls.gz fell through to the CSV reader and failed on the workbook's first
    byte - "utf-8 codec can't decode byte 0xd0" - which says nothing about the
    real problem and nothing a reader could act on.
    """

    import gzip
    import shutil

    from openpyxl import Workbook

    book = Workbook()
    sheet = book.active
    sheet.append(["gene", "log2FoldChange", "pvalue"])
    for index in range(1, 6):
        sheet.append([f"G{index}", 2.0 - index * 0.1, 0.001 * index])
    plain = tmp_path / "table.xlsx"
    book.save(plain)
    compressed = tmp_path / "table.xlsx.gz"
    with plain.open("rb") as source, gzip.open(compressed, "wb") as target:
        shutil.copyfileobj(source, target)

    mapping = TableMapping("gene", "log2FoldChange", "pvalue")

    assert read_deg_table(plain, mapping).shape == (5, 3)
    assert read_deg_table(compressed, mapping).shape == (5, 3)
