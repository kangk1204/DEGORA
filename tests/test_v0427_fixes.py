"""Regressions for the v0.4.26 options/exceptions review (M1-M3, S1-S3, L1)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora.slice_runner import DegoraConfigError, validate_catalog_inputs


def _table(seed: int, n: int = 150) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    p = np.clip(rng.uniform(1e-6, 0.2, n), 1e-9, 1)
    return pd.DataFrame({"gene": [f"G{i}" for i in range(n)], "log2FoldChange": rng.normal(0, 1, n), "pvalue": p, "padj": np.minimum(p * 3, 1.0)})


def _config(tmp_path: Path, **first_row) -> Path:
    _table(1).to_csv(tmp_path / "t1.csv", index=False)
    _table(2).to_csv(tmp_path / "t2.csv", index=False)
    base = {"study_id": "S1", "source_unit_id": "U1", "source_path": "t1.csv", "gene_column": "gene",
            "lfc_column": "log2FoldChange", "p_column": "pvalue", "padj_column": "padj", "include_in_analysis": "yes"}
    row1 = dict(base, **first_row)
    row2 = dict(base, study_id="S2", source_unit_id="U2", source_path="t2.csv")
    cfg = tmp_path / "cfg.csv"
    pd.DataFrame([row1, row2]).to_csv(cfg, index=False)
    return cfg


# ---- M1 ------------------------------------------------------------------------------

def test_the_analysis_result_carries_one_warning_list_and_the_card_renders_it() -> None:
    """Three warnings were computed and stored; the browser card showed a count and eight genes."""

    from degora.api import INDEX_HTML
    from degora.slice_runner import run_warning_messages

    assert run_warning_messages({"warnings": ["a"], "identifier_space_warnings": ["b"], "rank_universe_warnings": ["c"]}) == ["a", "b", "c"]
    assert '<ul id="analysisCompleteWarnings" class="analysis-warnings" role="status" hidden></ul>' in INDEX_HTML
    assert "const runWarnings = [...new Set((state.run.warnings || [])" in INDEX_HTML
    assert 'warningList.innerHTML = runWarnings.map((item) => `<li>${esc(item)}</li>`).join("");' in INDEX_HTML


def test_discovery_analysis_result_includes_run_warnings(tmp_path, monkeypatch) -> None:
    import inspect as _inspect

    import degora.discovery_run as dr

    source = _inspect.getsource(dr._execute_discovery_analysis)
    assert '"warnings": list(' in source and "run_warning_messages(metrics)" in source and "dict.fromkeys(" in source


# ---- M2 ------------------------------------------------------------------------------

def test_a_catalog_pointing_p_column_at_the_adjusted_column_is_warned_about(tmp_path) -> None:
    """init explains it and the browser requires a confirmation; a hand-written catalog said nothing."""

    result = validate_catalog_inputs(_config(tmp_path, p_column="padj", padj_column=""))
    assert any("p_column='padj' is an adjusted p-value/FDR column" in w for w in result["warnings"])

    # The same column named as both p and padj is the other way in.
    result = validate_catalog_inputs(_config(tmp_path, p_column="padj", padj_column="padj"))
    assert any("already adjusted" in w for w in result["warnings"])

    # An ordinary mapping says nothing.
    result = validate_catalog_inputs(_config(tmp_path))
    assert not any("already adjusted" in w for w in result["warnings"])


# ---- M3 ------------------------------------------------------------------------------

def test_matrix_candidates_carry_the_gene_identifier_space() -> None:
    from degora.discovery import _inspect_upstream_rows

    rows = [["gene", "GSM1", "GSM2", "GSM3", "GSM4"]] + [[f"ENSG{i:011d}", 10 + i, 11 + i, 20 + i, 21 + i] for i in range(40)]
    result = _inspect_upstream_rows(rows, declared_role="unknown_matrix")
    assert result["status"] == "upstream_matrix_ready_for_contrast"
    assert result["gene_identifier_space"] == "Ensembl ID"

    symbols = [["gene", "GSM1", "GSM2", "GSM3", "GSM4"]] + [[f"G{i}", 10 + i, 11 + i, 20 + i, 21 + i] for i in range(40)]
    assert _inspect_upstream_rows(symbols, declared_role="unknown_matrix")["gene_identifier_space"] == "gene symbol"


# ---- S1 ------------------------------------------------------------------------------

def test_an_invalid_time_course_mode_is_refused_under_a_header_that_names_it(tmp_path) -> None:
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(_config(tmp_path, time_course_mode="fastest"))
    message = str(excinfo.value)
    assert "optional contrast settings are not valid (table_scope, time_course_mode, rank_universe_size)" in message
    assert "time_course_mode" in message


# ---- S2 ------------------------------------------------------------------------------

def test_a_sheet_name_on_a_csv_source_is_reported_as_ignored(tmp_path) -> None:
    result = validate_catalog_inputs(_config(tmp_path, sheet_name="Nope"))
    assert any("sheet_name='Nope' is set but t1.csv is not a workbook, so it is ignored" in w for w in result["warnings"])


def test_a_stray_sheet_name_cannot_split_one_csv_into_two_source_units(tmp_path) -> None:
    """The same CSV under two units, one row with sheet_name set, passed the duplicate-table check."""

    _table(1).to_csv(tmp_path / "same.csv", index=False)
    base = {"study_id": "S1", "source_unit_id": "U1", "source_path": "same.csv", "gene_column": "gene",
            "lfc_column": "log2FoldChange", "p_column": "pvalue", "include_in_analysis": "yes"}
    cfg = tmp_path / "cfg.csv"
    pd.DataFrame([base, dict(base, study_id="S2", source_unit_id="U2", sheet_name="Nope")]).to_csv(cfg, index=False)
    with pytest.raises(DegoraConfigError, match="two independent source units"):
        validate_catalog_inputs(cfg)


# ---- S3 ------------------------------------------------------------------------------

def test_a_readiness_tier_with_nothing_behind_it_is_not_kept() -> None:
    """A 1959 publication wore likely_ready only because a tier was taken on trust."""

    from degora.discovery import normalize_species
    from degora.discovery_federated import _readiness

    record = {"canonical_id": "pmid:1", "year": 1959, "candidates": [], "data_readiness": {"tier": "likely_ready", "verification_state": "likely_ready"}}
    readiness = _readiness(record, normalize_species("human"))
    assert readiness["tier"] == "metadata_only"
    assert "publication_metadata_only" in readiness["basis"]


# ---- L1 ------------------------------------------------------------------------------

def test_the_preparation_marker_is_named_when_handed_to_the_analysis(tmp_path) -> None:
    from degora.discovery_run import DiscoveryError, run_discovery_analysis

    marker = {"artifact_type": "degora_discovery_prepared_bundle", "format_version": 1, "species": "human"}
    with pytest.raises(DiscoveryError, match="preparation folder's marker"):
        run_discovery_analysis(marker, [], tmp_path / "out", species="human")



def test_seurat_and_scanpy_adjusted_headers_are_warned_about_as_p_column(tmp_path) -> None:
    """p_val_adj and pvals_adj also match the nominal pattern; the warning missed exactly them."""

    for header in ("p_val_adj", "pvals_adj", "adj.P.Val", "FDR"):
        t1 = _table(1); t1 = t1.rename(columns={"padj": header}); t1.to_csv(tmp_path / "t1.csv", index=False)
        _table(2).to_csv(tmp_path / "t2.csv", index=False)
        base = {"study_id": "S1", "source_unit_id": "U1", "source_path": "t1.csv", "gene_column": "gene",
                "lfc_column": "log2FoldChange", "p_column": header, "padj_column": "", "include_in_analysis": "yes"}
        cfg = tmp_path / "cfg.csv"
        pd.DataFrame([base, dict(base, study_id="S2", source_unit_id="U2", source_path="t2.csv", p_column="pvalue", padj_column="padj")]).to_csv(cfg, index=False)
        result = validate_catalog_inputs(cfg)
        assert any("already adjusted" in x for x in result["warnings"]), header
    # A nominal header never triggers it.
    for header in ("pvalue", "P.Value", "pval", "p_val", "PValue"):
        t1 = _table(1).rename(columns={"pvalue": header}); t1.to_csv(tmp_path / "t1.csv", index=False)
        base = {"study_id": "S1", "source_unit_id": "U1", "source_path": "t1.csv", "gene_column": "gene",
                "lfc_column": "log2FoldChange", "p_column": header, "padj_column": "padj", "include_in_analysis": "yes"}
        cfg = tmp_path / "cfg.csv"
        pd.DataFrame([base, dict(base, study_id="S2", source_unit_id="U2", source_path="t2.csv", p_column="pvalue")]).to_csv(cfg, index=False)
        assert not any("already adjusted" in x for x in validate_catalog_inputs(cfg)["warnings"]), header


def test_a_blank_and_a_one_header_row_are_the_same_table_identity(tmp_path) -> None:
    _table(1).to_csv(tmp_path / "same.csv", index=False)
    base = {"study_id": "S1", "source_unit_id": "U1", "source_path": "same.csv", "gene_column": "gene",
            "lfc_column": "log2FoldChange", "p_column": "pvalue", "include_in_analysis": "yes", "header_row": ""}
    cfg = tmp_path / "cfg.csv"
    pd.DataFrame([base, dict(base, study_id="S2", source_unit_id="U2", header_row=1)]).to_csv(cfg, index=False)
    with pytest.raises(DegoraConfigError, match="two independent source units"):
        validate_catalog_inputs(cfg)


def test_the_candidate_row_names_the_identifier_space_and_the_cli_prints_run_warnings() -> None:
    from degora.api import INDEX_HTML
    from degora.cli import main  # noqa: F401 - import guard only

    assert "function identifierSpaceNote(candidate)" in INDEX_HTML
    assert INDEX_HTML.count("${identifierSpaceNote(candidate)}") >= 2
    source = Path("degora/cli.py").read_text(encoding="utf-8")
    assert 'result.get("warnings") or result.get("selection_warnings", [])' in source
