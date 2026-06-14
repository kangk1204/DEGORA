from __future__ import annotations

import json

import pandas as pd
import pytest

from degora.slice_runner import DegoraConfigError, catalog_include_mask, read_catalog, run_slice, validate_catalog_inputs


def test_read_catalog_accepts_beginner_excel_aliases(tmp_path) -> None:
    config_path = tmp_path / "degora_config.xlsx"
    pd.DataFrame(
        {
            "study_id": ["IFN_GSE001_4h"],
            "source_unit_id": ["GSE001"],
            "source_path": ["data/deg/raw/ifn_4h.csv"],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
            "time_h": [4],
            "condition": ["IFN-beta"],
            "include": ["yes"],
        }
    ).to_excel(config_path, sheet_name="Contrasts", index=False)

    catalog = read_catalog(config_path)

    assert catalog.loc[0, "paper_id"] == "GSE001"
    assert catalog.loc[0, "duration_h"] == 4
    assert catalog.loc[0, "hypoxia_modality"] == "IFN-beta"
    assert catalog.loc[0, "pipeline"] == "unknown_pipeline"
    assert "padj_column" in catalog.columns


def test_read_catalog_ignores_leading_excel_note_rows(tmp_path) -> None:
    config_path = tmp_path / "commented_config.xlsx"
    contrasts = pd.DataFrame(
        {
            "study_id": ["IFN_GSE001_4h"],
            "source_unit_id": ["GSE001"],
            "source_path": ["data/deg/raw/ifn_4h.csv"],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
            "include": ["yes"],
        }
    )
    with pd.ExcelWriter(config_path, engine="openpyxl") as writer:
        contrasts.to_excel(writer, sheet_name="Contrasts", index=False, startrow=1)
        writer.sheets["Contrasts"]["A1"] = "# Main input sheet; DEGORA should ignore this row."

    catalog = read_catalog(config_path)

    assert catalog.loc[0, "study_id"] == "IFN_GSE001_4h"
    assert catalog.loc[0, "paper_id"] == "GSE001"
    assert "source_path" in catalog.columns


def test_read_catalog_accepts_microarray_aliases(tmp_path) -> None:
    config_path = tmp_path / "microarray_config.csv"
    pd.DataFrame(
        {
            "study_id": ["MICRO001"],
            "source_unit_id": ["GSE_MICRO"],
            "source_path": ["microarray_deg.csv"],
            "gene_column": ["gene_symbol"],
            "lfc_column": ["logFC"],
            "p_column": ["P.Value"],
            "technology": ["microarray"],
            "input_type": ["limma_full_table"],
            "platform_id": ["GPL570"],
            "probe_column": ["probe_id"],
            "probe_collapse_rule": ["min_pvalue_max_abs_lfc"],
        }
    ).to_csv(config_path, index=False)

    catalog = read_catalog(config_path)

    assert catalog.loc[0, "assay_type"] == "microarray"
    assert catalog.loc[0, "source_input_type"] == "limma_full_table"
    assert catalog.loc[0, "platform"] == "GPL570"
    assert catalog.loc[0, "probe_id_column"] == "probe_id"
    assert catalog.loc[0, "probe_collapse"] == "min_pvalue_max_abs_lfc"


def test_read_catalog_explains_missing_required_columns(tmp_path) -> None:
    config_path = tmp_path / "bad_config.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"],
            "paper_id": ["P1"],
            "source_path": ["source.csv"],
            "gene_column": ["gene"],
            "p_column": ["pvalue"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        read_catalog(config_path)

    message = str(exc_info.value)
    assert "DEGORA config error" in message
    assert "lfc_column" in message
    assert "How to fix" in message
    assert "Contrasts sheet" in message


def test_excel_without_contrasts_sheet_gets_actionable_error(tmp_path) -> None:
    config_path = tmp_path / "bad_config.xlsx"
    pd.DataFrame({"field": ["project_name"], "value": ["IFN"]}).to_excel(
        config_path, sheet_name="Project", index=False
    )
    with pd.ExcelWriter(config_path, mode="a", engine="openpyxl") as writer:
        pd.DataFrame({"gene_symbol": ["ISG15"]}).to_excel(writer, sheet_name="GoldPanel", index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        read_catalog(config_path)

    message = str(exc_info.value)
    assert "missing a Contrasts sheet" in message
    assert "Found sheets" in message
    assert "Rename the sheet" in message


def test_catalog_include_mask_explains_allowed_values() -> None:
    catalog = pd.DataFrame({"include_in_analysis": ["maybe"]})

    with pytest.raises(DegoraConfigError) as exc_info:
        catalog_include_mask(catalog)

    message = str(exc_info.value)
    assert "include flag" in message
    assert "yes/true/include/1" in message
    assert "no/false/exclude/0" in message


def test_run_slice_explains_bad_table_scope(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["ISG15"],
            "log2FoldChange": [2.0],
            "pvalue": [0.001],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "bad_scope.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"],
            "paper_id": ["P1"],
            "source_path": [str(source_path)],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
            "table_scope": ["only hits"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", min_studies=1)

    message = str(exc_info.value)
    assert "table-scope settings are not valid" in message
    assert "full_results" in message
    assert "deg_only" in message


def test_run_slice_treats_whitespace_required_value_as_empty(tmp_path) -> None:
    config_path = tmp_path / "blank_required_value.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"],
            "paper_id": ["P1"],
            "source_path": ["   "],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", min_studies=1)

    message = str(exc_info.value)
    assert "empty required values" in message
    assert "source_path" in message
    assert "Fill the highlighted cells" in message


def test_validate_catalog_strips_paper_ids_before_counting_source_units(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["GENEA"],
            "log2FoldChange": [1.0],
            "pvalue": [0.01],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "blank_paper_ids.csv"
    pd.DataFrame(
            {
                "study_id": ["S1", "S2"],
                "paper_id": ["P1", " P1 "],
            "source_path": [str(source_path), str(source_path)],
            "gene_column": ["gene", "gene"],
            "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
        }
    ).to_csv(config_path, index=False)

    summary = validate_catalog_inputs(config_path)

    assert summary["source_units"] == 1


def test_read_catalog_preserves_explicit_source_unit_id_alongside_paper_id(tmp_path) -> None:
    config_path = tmp_path / "explicit_units.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["PAPER", "PAPER"],
            "source_unit_id": ["UNIT_A", "UNIT_B"],
            "source_path": ["x.csv", "x.csv"],
            "gene_column": ["gene", "gene"],
            "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
        }
    ).to_csv(config_path, index=False)

    catalog = read_catalog(config_path)

    assert "source_unit_id" in catalog.columns
    assert catalog["source_unit_id"].tolist() == ["UNIT_A", "UNIT_B"]
    assert catalog["paper_id"].tolist() == ["PAPER", "PAPER"]


def test_validate_catalog_counts_explicit_source_units_over_shared_paper_id(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame({"gene": ["GENEA"], "log2FoldChange": [1.0], "pvalue": [0.01]}).to_csv(source_path, index=False)
    config_path = tmp_path / "explicit_units.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["PAPER", "PAPER"],
            "source_unit_id": ["UNIT_A", "UNIT_B"],
            "source_path": [str(source_path), str(source_path)],
            "gene_column": ["gene", "gene"],
            "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
        }
    ).to_csv(config_path, index=False)

    summary = validate_catalog_inputs(config_path)

    assert summary["source_units"] == 2


def test_read_catalog_normalizes_header_case_spaces_and_hyphens(tmp_path) -> None:
    config_path = tmp_path / "messy_headers.csv"
    pd.DataFrame(
        {
            "Study ID": ["S1"],
            "Paper.ID": ["P1"],
            "Source Path": ["x.csv"],
            "Gene   Column": ["gene"],
            "LFC Column": ["log2FoldChange"],
            "P.Column": ["pvalue"],
        }
    ).to_csv(config_path, index=False)

    catalog = read_catalog(config_path)

    # Headers are matched case/separator-insensitively...
    assert catalog.loc[0, "study_id"] == "S1"
    assert catalog.loc[0, "paper_id"] == "P1"
    assert catalog.loc[0, "source_path"] == "x.csv"
    # ...but the mapping VALUES are preserved verbatim (matched against the DEG file exactly).
    assert catalog.loc[0, "gene_column"] == "gene"
    assert catalog.loc[0, "lfc_column"] == "log2FoldChange"
    assert catalog.loc[0, "p_column"] == "pvalue"


def test_read_catalog_rejects_ambiguous_duplicate_headers(tmp_path) -> None:
    config_path = tmp_path / "dup_headers.csv"
    # 'paper_id' plus a differently-cased 'Paper ID' both map to one field: rather than
    # silently keeping one and dropping the other, the catalog must be rejected.
    pd.DataFrame(
        {
            "study_id": ["S1"],
            "paper_id": ["CANONICAL"],
            "Paper ID": ["DUPLICATE"],
            "source_path": ["x.csv"],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        read_catalog(config_path)

    message = str(exc_info.value)
    assert "duplicate or ambiguous column headers" in message
    assert "paper_id" in message


def test_read_catalog_explains_missing_config_file(tmp_path) -> None:
    with pytest.raises(DegoraConfigError) as exc_info:
        read_catalog(tmp_path / "does_not_exist.xlsx")

    message = str(exc_info.value)
    assert "config file was not found" in message
    assert "degora template" in message


def test_validate_catalog_inputs_explains_missing_config_file(tmp_path) -> None:
    with pytest.raises(DegoraConfigError) as exc_info:
        validate_catalog_inputs(tmp_path / "missing_config.xlsx")

    assert "config file was not found" in str(exc_info.value)


def test_run_slice_explains_bad_time_course_mode(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["ISG15"],
            "log2FoldChange": [2.0],
            "pvalue": [0.001],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "bad_time_mode.csv"
    pd.DataFrame(
        {
            "study_id": ["S1"],
            "paper_id": ["P1"],
            "source_path": [str(source_path)],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
            "time_course_mode": ["strongest_only"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", min_studies=1)

    message = str(exc_info.value)
    assert "time_course_mode" in message
    assert "mean, early, late, peak_mean" in message


def test_run_slice_explains_wrong_source_column_and_suggests_fix(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["ISG15", "RPL13A"],
            "log2FoldChange": [2.0, 0.1],
            "pvalue": [0.001, 0.8],
            "padj": [0.01, 0.9],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "degora_config.csv"
    pd.DataFrame(
        {
            "study_id": ["IFN_GSE001_4h"],
            "paper_id": ["GSE001"],
            "source_path": [str(source_path)],
            "gene_column": ["gene"],
            "lfc_column": ["logFC"],
            "p_column": ["pvalue"],
            "padj_column": ["padj"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", min_studies=1)

    message = str(exc_info.value)
    assert "source table column mapping is wrong" in message
    assert "lfc_column='logFC'" in message
    assert "log2FoldChange" in message
    assert "Column names are case-sensitive" in message


def test_run_slice_explains_optional_padj_column_only_when_filled(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["ISG15", "RPL13A"],
            "log2FoldChange": [2.0, 0.1],
            "pvalue": [0.001, 0.8],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "degora_config.csv"
    pd.DataFrame(
        {
            "study_id": ["IFN_GSE001_4h"],
            "paper_id": ["GSE001"],
            "source_path": [str(source_path)],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
            "padj_column": ["padj"],
        }
    ).to_csv(config_path, index=False)

    with pytest.raises(DegoraConfigError) as exc_info:
        validate_catalog_inputs(config_path)

    message = str(exc_info.value)
    assert "optional padj_column" in message
    assert "Leave padj_column blank" in message
    assert "checked only when you fill them" in message


def test_run_slice_metrics_are_json_serializable_with_blank_optional_metadata(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["RBM39", "TYMS"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.001, 0.002],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "mixed_optional_metadata.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "source_path": [str(source_path), str(source_path)],
            "gene_column": ["gene", "gene"],
            "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
            "pipeline": ["", "welch_microarray_normalized_matrix"],
            "assay_type": ["", "microarray"],
            "source_input_type": ["", "normalized_expression_matrix"],
        }
    ).to_csv(config_path, index=False)

    metrics = run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", min_studies=1)

    metrics_path = tmp_path / "out" / "slice_metrics.json"
    persisted = json.loads(metrics_path.read_text())
    assert metrics["pipeline_counts"]["unknown"] == 1
    assert persisted["assay_type_counts"] == {"microarray": 1, "unknown": 1}
    assert persisted["gold_panel_status"] == "not_provided"
    assert persisted["recall_at_50"]["status"] == "not_applicable"


def test_run_slice_uses_excel_gold_panel_and_records_source_inputs(tmp_path) -> None:
    source_path = tmp_path / "deg.csv"
    pd.DataFrame(
        {
            "gene": ["GENEA", "GENEB"],
            "log2FoldChange": [2.0, 0.5],
            "pvalue": [1e-6, 0.2],
        }
    ).to_csv(source_path, index=False)
    config_path = tmp_path / "degora_config.xlsx"
    contrasts = pd.DataFrame(
        {
            "study_id": ["S1"],
            "source_unit_id": ["P1"],
            "source_path": [str(source_path)],
            "gene_column": ["gene"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
            "include": ["yes"],
        }
    )
    gold = pd.DataFrame({"gene_symbol": ["GENEA"], "locked": ["yes"]})
    with pd.ExcelWriter(config_path, engine="openpyxl") as writer:
        contrasts.to_excel(writer, sheet_name="Contrasts", index=False, startrow=1)
        writer.sheets["Contrasts"]["A1"] = "# Contrast rows."
        gold.to_excel(writer, sheet_name="GoldPanel", index=False, startrow=1)
        writer.sheets["GoldPanel"]["A1"] = "# Optional recall panel."

    metrics = run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", min_studies=1)
    provenance = json.loads((tmp_path / "out" / "slice_metrics.json.provenance.json").read_text())

    assert metrics["gold_panel_status"] == "locked"
    assert metrics["gold_panel_gene_count"] == 1
    assert metrics["recall_at_50"]["recall"] == 1.0
    assert str(source_path) in metrics["source_input_files"]
    assert any(record["path"] == str(source_path.resolve()) for record in provenance["inputs"])
    assert "SLICE_MIN_STUDIES=1" in provenance["command"]
