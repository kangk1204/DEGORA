from __future__ import annotations

import sqlite3

import pandas as pd
from openpyxl import load_workbook

from degora.cli import main
from degora.excel_template import TEMPLATE_SHEETS, write_template


def _write_source(path, genes, lfc_scale: float) -> None:
    pd.DataFrame(
        {
            "gene": genes,
            "log2FoldChange": [2.0 * lfc_scale, 1.5 * lfc_scale, 0.1],
            "pvalue": [1e-6, 1e-4, 0.8],
            "padj": [1e-5, 1e-3, 0.9],
        }
    ).to_csv(path, index=False)


def _write_config(path, source_a, source_b) -> None:
    project = pd.DataFrame(
        {
            "field": ["output_dir", "harmonized_dir", "min_studies"],
            "value": [str(path.parent / "results"), str(path.parent / "harmonized"), 1],
        }
    )
    contrasts = pd.DataFrame(
        {
            "study_id": ["S1_4h", "S2_4h"],
            "source_unit_id": ["P1", "P2"],
            "source_path": [str(source_a), str(source_b)],
            "gene_column": ["gene", "gene"],
            "lfc_column": ["log2FoldChange", "log2FoldChange"],
            "p_column": ["pvalue", "pvalue"],
            "padj_column": ["padj", "padj"],
            "include": ["yes", "yes"],
        }
    )
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        project.to_excel(writer, sheet_name="Project", index=False)
        contrasts.to_excel(writer, sheet_name="Contrasts", index=False)


def test_template_command_writes_beginner_workbook(tmp_path) -> None:
    output = tmp_path / "DEGORA_template.xlsx"

    assert main(["template", str(output)]) == 0

    workbook = load_workbook(output)
    assert workbook.sheetnames == TEMPLATE_SHEETS
    assert workbook["Contrasts"]["A1"].value == "study_id"
    assert workbook["ColumnGuide"]["A1"].value == "column"
    contrast_headers = [cell.value for cell in workbook["Contrasts"][1]]
    guide_columns = [cell.value for cell in workbook["ColumnGuide"]["A"]]
    assert "time_course_mode" in contrast_headers
    assert "time_course_mode" in guide_columns


def test_demo_command_writes_runnable_workspace(tmp_path) -> None:
    demo_dir = tmp_path / "demo"

    assert main(["demo", str(demo_dir)]) == 0

    config = demo_dir / "degora_demo_config.xlsx"
    assert config.exists()
    assert (demo_dir / "deg_tables" / "demo_ifn_a_4h.csv").exists()
    assert (demo_dir / "README.md").exists()
    assert main(["validate", str(config)]) == 0
    assert main(["run", str(config)]) == 0

    db = demo_dir / "results" / "degora_scores.db"
    assert db.exists()
    with sqlite3.connect(db) as connection:
        top_genes = [
            row[0]
            for row in connection.execute("SELECT gene_symbol FROM genes ORDER BY degora_rank LIMIT 4").fetchall()
        ]
        source_units = connection.execute("SELECT COUNT(DISTINCT source_unit_id) FROM studies").fetchone()[0]
    assert top_genes[0] == "ISG15"
    assert {"IFIT1", "MX1"}.issubset(set(top_genes))
    assert source_units == 2


def test_validate_command_accepts_excel_config(tmp_path) -> None:
    source = tmp_path / "source.csv"
    _write_source(source, ["ISG15", "IFIT1", "RPL13A"], 1.0)
    config = tmp_path / "config.xlsx"
    _write_config(config, source, source)

    assert main(["validate", str(config)]) == 0


def test_run_command_builds_score_database_from_excel_config(tmp_path) -> None:
    source_a = tmp_path / "source_a.csv"
    source_b = tmp_path / "source_b.csv"
    _write_source(source_a, ["ISG15", "IFIT1", "RPL13A"], 1.0)
    _write_source(source_b, ["ISG15", "IFIT1", "RPL13A"], 0.8)
    config = tmp_path / "config.xlsx"
    _write_config(config, source_a, source_b)
    db = tmp_path / "results" / "degora_scores.db"

    assert main(["run", str(config), "--db", str(db)]) == 0

    assert db.exists()
    assert (tmp_path / "results" / "degora_gene_scores.csv").exists()
    assert (tmp_path / "results" / "degora_score_metadata.json").exists()
    with sqlite3.connect(db) as connection:
        top_gene = connection.execute("SELECT gene_symbol FROM genes ORDER BY degora_rank LIMIT 1").fetchone()[0]
        source_units = connection.execute("SELECT COUNT(DISTINCT source_unit_id) FROM studies").fetchone()[0]
    assert top_gene == "ISG15"
    assert source_units == 2
