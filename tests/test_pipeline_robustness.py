"""Pipeline-robustness regression tests.

Each test pins one (stage, failure-mode) cell from the adversarial robustness audit:
it reproduces an input pathology that previously crashed or silently corrupted output,
and asserts the hardened behavior (clear error or correct result).
"""

from __future__ import annotations

from decimal import Decimal
import json

import numpy as np
import pandas as pd
import pytest

from degora.aggregate import slice_consensus
from degora.harmonize import TableMapping, harmonize_frame, read_deg_table
from degora.metrics import recall_at_k
from degora.slice_runner import (
    DegoraConfigError,
    catalog_include_mask,
    run_slice,
)


def _write_deg(path) -> None:
    pd.DataFrame(
        {"gene": ["AAA", "BBB"], "log2FoldChange": [2.0, -1.0], "pvalue": [1e-4, 1e-3]}
    ).to_csv(path, index=False)


def _min_harmonized() -> pd.DataFrame:
    rows = []
    for unit in ("P1", "P2"):
        for gene, z in (("VEGFA", 5.0), ("HK2", 4.0)):
            rows.append(
                {
                    "study_id": f"{unit}_S", "paper_id": unit, "gene_symbol": gene,
                    "lfc": 2.0, "signed_z": z, "pvalue": 1e-6, "padj": 1e-5,
                    "normalized_rank": 0.02, "n_ctrl": 3, "n_treat": 3, "n_genes_in_study": 1000,
                    "pipeline": "DESeq2", "assay_type": "RNA-seq", "source_input_type": "author_deg_table",
                    "platform": "", "normalization": "DESeq2", "probe_collapse": "",
                    "species": "Homo sapiens", "cell_system": "A", "hypoxia_modality": "x",
                    "duration_h": "24", "source_path": "s.csv", "source_url": "u",
                }
            )
    return pd.DataFrame(rows)


# --- config parsing -------------------------------------------------------------------

def test_include_flag_normalizes_float_promoted_values() -> None:
    # Excel float-promotes an int include column with any blank cell -> 1.0/0.0/NaN.
    mask = catalog_include_mask(pd.DataFrame({"include_in_analysis": [1.0, 0.0, float("nan")]}))
    assert mask.tolist() == [True, False, True]


def test_run_slice_rejects_duplicate_active_study_id(tmp_path) -> None:
    # Two active rows sharing a study_id would silently double-count one contrast
    # inside the within-source-unit weighted mean; it must be rejected up front.
    src = tmp_path / "deg.csv"
    _write_deg(src)
    cfg = tmp_path / "c.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S1"], "paper_id": ["P1", "P2"], "source_path": [str(src), str(src)],
            "gene_column": ["gene", "gene"], "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
        }
    ).to_csv(cfg, index=False)
    with pytest.raises(DegoraConfigError, match="duplicate study_id"):
        run_slice(cfg, tmp_path / "o", tmp_path / "h", min_studies=1)


def test_run_slice_rejects_fractional_rank_universe(tmp_path) -> None:
    src = tmp_path / "deg.csv"
    _write_deg(src)
    cfg = tmp_path / "c.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"], "paper_id": ["P1"], "source_path": [str(src)],
            "gene_column": ["gene"], "lfc_column": ["log2FoldChange"], "p_column": ["pvalue"],
            "table_scope": ["deg_only"], "rank_universe_size": ["10.5"],
        }
    ).to_csv(cfg, index=False)
    with pytest.raises(DegoraConfigError, match="whole number of genes"):
        run_slice(cfg, tmp_path / "o", tmp_path / "h", min_studies=1)


def test_run_slice_metrics_json_valid_with_blank_excluded_notes(tmp_path) -> None:
    src = tmp_path / "deg.csv"
    _write_deg(src)
    cfg = tmp_path / "c.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S2"], "paper_id": ["P1", "P2"], "source_path": [str(src), str(src)],
            "gene_column": ["gene", "gene"], "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"], "include_in_analysis": ["yes", "no"],
        }
    ).to_csv(cfg, index=False)
    run_slice(cfg, tmp_path / "o", tmp_path / "h", min_studies=1)
    text = (tmp_path / "o" / "slice_metrics.json").read_text()
    assert "NaN" not in text
    json.loads(text)  # must parse as strict JSON


def test_run_slice_metrics_annotate_recall_ranking_basis(tmp_path) -> None:
    # C6: slice_metrics recall is computed over the Stouffer slice_consensus order, NOT
    # the primary quality_weighted_degora_rank. The metrics (and the persisted
    # slice_metrics.json) must carry an explicit ranking-basis annotation so the recall
    # cannot be misread as the quality-weighted DEGORA recall from score_db.
    src = tmp_path / "deg.csv"
    _write_deg(src)
    cfg = tmp_path / "c.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S2"], "paper_id": ["P1", "P2"], "source_path": [str(src), str(src)],
            "gene_column": ["gene", "gene"], "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
        }
    ).to_csv(cfg, index=False)
    metrics = run_slice(cfg, tmp_path / "o", tmp_path / "h", min_studies=1)
    assert metrics["recall_rank_source"] == "slice_consensus_order"
    assert "quality-weighted" in metrics["recall_rank_note"]
    persisted = json.loads((tmp_path / "o" / "slice_metrics.json").read_text())
    assert persisted["recall_rank_source"] == "slice_consensus_order"
    assert persisted["recall_rank_note"] == metrics["recall_rank_note"]


def test_run_slice_treats_parquet_engine_absence_as_optional_output_warning(tmp_path, monkeypatch) -> None:
    src = tmp_path / "deg.csv"
    _write_deg(src)
    cfg = tmp_path / "c.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"], "paper_id": ["P1"], "source_path": [str(src)],
            "gene_column": ["gene"], "lfc_column": ["log2FoldChange"], "p_column": ["pvalue"],
        }
    ).to_csv(cfg, index=False)

    def _missing_parquet_engine(self, *args, **kwargs):
        raise ImportError("no parquet engine for this test")

    monkeypatch.setattr(pd.DataFrame, "to_parquet", _missing_parquet_engine)

    metrics = run_slice(cfg, tmp_path / "o", tmp_path / "h", min_studies=1)

    assert (tmp_path / "h" / "o_harmonized.csv").exists()
    assert not (tmp_path / "h" / "o_harmonized.parquet").exists()
    assert not (tmp_path / "h" / "o_harmonized.parquet.source").exists()
    assert metrics["warnings"] == []
    assert "no usable parquet engine" in metrics["optional_output_warnings"][0]
    persisted = json.loads((tmp_path / "o" / "slice_metrics.json").read_text())
    assert persisted["optional_output_warnings"] == metrics["optional_output_warnings"]


def test_run_slice_output_dir_is_a_file_is_explained(tmp_path) -> None:
    src = tmp_path / "deg.csv"
    _write_deg(src)
    cfg = tmp_path / "c.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"], "paper_id": ["P1"], "source_path": [str(src)],
            "gene_column": ["gene"], "lfc_column": ["log2FoldChange"], "p_column": ["pvalue"],
        }
    ).to_csv(cfg, index=False)
    out_as_file = tmp_path / "out_is_a_file"
    out_as_file.write_text("x")
    with pytest.raises(DegoraConfigError, match="could not create the output directory"):
        run_slice(cfg, out_as_file, tmp_path / "h", min_studies=1)


# --- read_deg_table -------------------------------------------------------------------

def test_read_deg_table_detects_wrong_separator(tmp_path) -> None:
    p = tmp_path / "tab_in.csv"  # tab-delimited content in a .csv -> auto sep=','
    p.write_text("gene\tlfc\tpval\nA\t1.0\t0.01\n")
    with pytest.raises(ValueError, match="single column"):
        read_deg_table(p, TableMapping("gene", "lfc", "pval"))


# --- harmonize_frame ------------------------------------------------------------------

def test_harmonize_pvalue_zero_gets_finite_top_rank() -> None:
    frame = pd.DataFrame(
        {"gene": ["TOPHIT", "B", "C"], "log2FoldChange": [3.0, 1.0, 0.5], "pvalue": [0.0, 1e-3, 1e-2]}
    )
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), {"study_id": "S1", "paper_id": "P1"})
    top = out.loc[out["gene_symbol"].eq("TOPHIT")].iloc[0]
    assert np.isfinite(top["signed_z"])
    assert not pd.isna(top["normalized_rank"])
    assert top["within_study_rank"] == 1.0  # p == 0 is the most significant -> top rank


def test_harmonize_rejects_structurally_invalid_pvalues() -> None:
    frame = pd.DataFrame(
        {"gene": ["BAD_NEG", "BAD_GT1"], "log2FoldChange": [2.0, 1.0], "pvalue": [-0.3, 2.0]}
    )
    with pytest.raises(ValueError, match=r"outside \[0, 1\]"):
        harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), {"study_id": "BAD_P", "paper_id": "P1"})


def test_harmonize_missing_gene_column_lists_available() -> None:
    frame = pd.DataFrame({"symbol": ["A"], "log2FoldChange": [1.0], "pvalue": [0.01]})
    with pytest.raises(KeyError, match="Available columns"):
        harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue"), {"study_id": "S1", "paper_id": "P1"})


# --- aggregate ------------------------------------------------------------------------

def test_source_unit_id_nan_string_falls_back_to_paper_id() -> None:
    # An all-empty source_unit_id catalog column round-trips through CSV/parquet as
    # the literal string "nan"; it must be treated as blank and fall back to paper_id,
    # not collapse every gene into one literal-"nan" source unit (which drops every
    # gene below min_studies and silently yields zero scored genes).
    from degora.aggregate import _source_unit_series

    frame = pd.DataFrame(
        {
            "study_id": ["S1", "S2"], "paper_id": ["P1", "P2"], "source_unit_id": ["nan", "nan"],
            "gene_symbol": ["G", "G"], "signed_z": [4.0, 4.5], "lfc": [2.0, 1.8],
            "pvalue": [1e-4, 1e-4], "padj": [1e-3, 1e-3], "normalized_rank": [0.01, 0.02],
            "n_ctrl": [3, 3], "n_treat": [3, 3],
        }
    )
    assert set(_source_unit_series(frame)) == {"P1", "P2"}
    consensus = slice_consensus(frame, min_studies=2)
    assert "G" in set(consensus["gene_symbol"])  # gene present in 2 source units now scores


def test_aggregate_drops_non_finite_signed_z() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S1", "S2"], "paper_id": ["P1", "P2", "P1", "P2"],
            "gene_symbol": ["G", "G", "H", "H"],
            "signed_z": [np.inf, 4.0, 3.0, 3.5], "lfc": [2.0, 1.8, 1.5, 1.6],
            "pvalue": [0.0, 1e-4, 1e-3, 1e-3], "padj": [0.0, 1e-3, 1e-2, 1e-2],
            "normalized_rank": [0.01, 0.02, 0.03, 0.04], "n_ctrl": [3, 3, 3, 3], "n_treat": [3, 3, 3, 3],
        }
    )
    out = slice_consensus(harmonized, min_studies=2)
    assert np.isfinite(out["stouffer_z"]).all()
    # G's inf unit is dropped -> only 1 finite unit -> below min_studies; H keeps 2 units.
    assert "G" not in set(out["gene_symbol"])
    assert "H" in set(out["gene_symbol"])


@pytest.mark.parametrize(
    "bad_min_studies",
    [0, -1, "0", "-1", True, False, "two", 2.5, "2.5", Decimal("2.5"), np.nan, Decimal("NaN")],
)
def test_aggregate_min_studies_must_be_integer_at_least_one(bad_min_studies) -> None:
    with pytest.raises(ValueError, match="min_studies must be an integer >= 1"):
        slice_consensus(_min_harmonized(), min_studies=bad_min_studies)


def test_aggregate_min_studies_accepts_integer_string_for_public_api_compatibility() -> None:
    consensus = slice_consensus(_min_harmonized(), min_studies="2")

    assert set(consensus["gene_symbol"]) == {"VEGFA", "HK2"}


@pytest.mark.parametrize("bad_min_studies", [True, 2.5, "2.5", Decimal("2.5"), np.nan])
def test_run_slice_uses_strict_min_studies_validator(tmp_path, bad_min_studies) -> None:
    with pytest.raises(DegoraConfigError, match="min_studies must be at least 1"):
        run_slice(tmp_path / "missing.csv", tmp_path / "o", tmp_path / "h", min_studies=bad_min_studies)

    assert not (tmp_path / "o").exists()
    assert not (tmp_path / "h").exists()


# --- metrics --------------------------------------------------------------------------

def test_recall_at_k_whitespace_and_type_robust() -> None:
    consensus = pd.DataFrame({"gene_symbol": ["VEGFA ", "hk2"]})  # dirty case + whitespace
    res = recall_at_k(consensus, [" VEGFA", "HK2", np.nan, 5], 10)
    assert set(res["recovered"]) == {"VEGFA", "HK2"}  # matched despite case/whitespace
    assert res["n_positives"] == 3  # NaN dropped; the int 5 coerced to "5" without crashing


def test_recall_at_k_deduplicates_before_top_k_slice() -> None:
    consensus = pd.DataFrame({"gene_symbol": ["VEGFA", "VEGFA", "HK2"]})

    res = recall_at_k(consensus, ["HK2"], 2)

    assert res["n_recovered"] == 1
    assert res["recall"] == 1.0


def test_recall_at_k_zero_and_negative_k_returns_empty_top_set() -> None:
    consensus = pd.DataFrame({"gene_symbol": ["VEGFA", "HK2"]})

    for k in (0, -1):
        res = recall_at_k(consensus, ["VEGFA"], k)
        assert res["n_recovered"] == 0
        assert res["recall"] == 0.0
        assert res["missing"] == ["VEGFA"]


def test_recall_at_k_requires_gene_symbol_column() -> None:
    with pytest.raises(ValueError, match="gene_symbol"):
        recall_at_k(pd.DataFrame({"other": [1, 2]}), ["A"], 10)

