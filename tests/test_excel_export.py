from __future__ import annotations

import json
import sqlite3
from zipfile import ZipFile

import pandas as pd
from openpyxl import load_workbook

from degora.excel_export import (
    _cap_evidence_for_sheet,
    _evidence_row_cap,
    _force_formula_like_text,
    export_run_workbook,
)


def test_evidence_sheet_cap_keeps_top_ranked_genes_and_flags() -> None:
    # 60 evidence rows across 6 genes; cap to 20 must keep the lowest-rank (best) genes' evidence.
    evidence = pd.DataFrame(
        {
            "gene_symbol": [f"G{i % 6}" for i in range(60)],
            "study_id": [f"S{i}" for i in range(60)],
        }
    )
    genes = pd.DataFrame(
        {
            "gene_symbol": [f"G{i}" for i in range(6)],
            "quality_weighted_degora_rank": [3, 1, 5, 2, 6, 4],
        }
    )
    capped, was_capped = _cap_evidence_for_sheet(evidence, genes, cap=20)
    assert was_capped is True
    assert len(capped) == 20
    # The two best-ranked genes are G1 (rank 1) and G3 (rank 2); their evidence must survive the cap.
    assert {"G1", "G3"} <= set(capped["gene_symbol"])
    # The worst-ranked gene G4 (rank 6) must be dropped first.
    assert "G4" not in set(capped["gene_symbol"])


def test_evidence_sheet_cap_noop_when_small_or_disabled() -> None:
    evidence = pd.DataFrame({"gene_symbol": ["G0", "G1"], "study_id": ["S0", "S1"]})
    genes = pd.DataFrame({"gene_symbol": ["G0", "G1"], "quality_weighted_degora_rank": [1, 2]})
    sheet, was_capped = _cap_evidence_for_sheet(evidence, genes, cap=100)
    assert was_capped is False and len(sheet) == 2
    sheet0, capped0 = _cap_evidence_for_sheet(evidence, genes, cap=0)
    assert capped0 is False and len(sheet0) == 2


def test_evidence_row_cap_reads_env(monkeypatch) -> None:
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "5000")
    assert _evidence_row_cap() == 5000
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "0")
    assert _evidence_row_cap() == 0  # 0 disables the cap
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "not-an-int")
    assert _evidence_row_cap() > 0  # invalid value falls back to the default


def test_export_caps_large_evidence_in_sql_but_reports_uncapped_total(tmp_path, monkeypatch) -> None:
    result_dir = tmp_path / "results"
    result_dir.mkdir()
    db_path = result_dir / "degora_scores.db"
    genes = pd.DataFrame(
        {
            "gene_symbol": ["TOP", "LOW"],
            "quality_weighted_degora_rank": [1, 2],
            "quality_weighted_degora_score": [0.9, 0.1],
        }
    )
    evidence = pd.DataFrame(
        {
            "gene_symbol": ["LOW"] + ["TOP"] * 100_000,
            "study_id": ["LOW_SOURCE"] + [f"S{i % 3}" for i in range(100_000)],
        }
    )
    studies = pd.DataFrame({"study_id": ["S0"], "source_unit_id": ["P0"]})
    with sqlite3.connect(db_path) as connection:
        genes.to_sql("genes", connection, index=False)
        evidence.to_sql("gene_evidence", connection, index=False)
        studies.to_sql("studies", connection, index=False)
        pd.DataFrame({"key": ["corpus"], "value": ["large-test"]}).to_sql(
            "meta", connection, index=False
        )
    (result_dir / "degora_score_metadata.json").write_text("{}\n", encoding="utf-8")
    monkeypatch.setenv("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "100")

    exported = export_run_workbook(result_dir, command="pytest large evidence export")

    assert exported["rows_gene_evidence"] == 100_001
    workbook = load_workbook(exported["output"], read_only=True, data_only=False)
    assert workbook["Gene_evidence"].max_row == 101
    assert {
        row[0].value: row[1].value
        for row in workbook["Run_summary"].iter_rows(min_row=2, values_only=False)
    }["gene_evidence_rows_total"] == 100_001
    summary = {
        row[0].value: row[1].value
        for row in workbook["Run_summary"].iter_rows(min_row=2, values_only=False)
    }
    assert "complete evidence for all genes is in degora_scores.db" in summary["gene_evidence_note"]
    assert "not source-level evidence" in summary["gene_evidence_note"]
    assert summary["path_base"] == "workbook_directory"
    assert summary["result_dir"] == "."
    exported_symbols = {
        row[0].value for row in workbook["Gene_evidence"].iter_rows(min_row=2, values_only=False)
    }
    assert exported_symbols == {"TOP"}
    manifest_path = result_dir / "DEGORA_output.manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["path_base"] == "manifest_directory"
    assert manifest["inputs"] == ["degora_scores.db", "degora_score_metadata.json"]
    assert manifest["outputs"] == [
        "DEGORA_output.xlsx",
        "DEGORA_output.manifest.json",
        "DEGORA_output.validation.txt",
    ]
    assert str(tmp_path) not in json.dumps(manifest)
    assert str(tmp_path) not in (result_dir / "DEGORA_output.xlsx.source").read_text(encoding="utf-8")


def test_formula_like_strings_are_written_as_literal_xlsx_text(tmp_path) -> None:
    output = tmp_path / "formula_safe.xlsx"
    values = ["=2+2", "+cmd", "-2+3", "@SUM(A1:A2)"]

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        pd.DataFrame({"gene_symbol": values}).to_excel(writer, sheet_name="Gene_scores", index=False)
        _force_formula_like_text(writer)

    workbook = load_workbook(output, data_only=False)
    cells = list(workbook["Gene_scores"].iter_rows(min_row=2, max_col=1))
    assert [row[0].value for row in cells] == values
    assert all(row[0].data_type == "s" for row in cells)
    with ZipFile(output) as archive:
        sheet_xml = archive.read("xl/worksheets/sheet1.xml")
    assert b"<f" not in sheet_xml
