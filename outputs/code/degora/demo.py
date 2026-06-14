"""Generate a tiny runnable DEGORA demo workspace."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .excel_template import _autosize_workbook, write_note_sheet


DEMO_CONFIG_NAME = "degora_demo_config.xlsx"


def _demo_rows(scale: float, *, late: bool = False) -> pd.DataFrame:
    """Return one small IFN-like DEG table with stable top markers."""

    rows: list[dict[str, Any]] = [
        {"gene": "ISG15", "log2FoldChange": 3.2 * scale, "pvalue": 1e-8, "padj": 1e-6},
        {"gene": "IFIT1", "log2FoldChange": 2.9 * scale, "pvalue": 2e-8, "padj": 2e-6},
        {"gene": "MX1", "log2FoldChange": 2.4 * scale, "pvalue": 8e-7, "padj": 4e-5},
        {"gene": "OAS1", "log2FoldChange": 2.0 * scale, "pvalue": 2e-6, "padj": 8e-5},
        {"gene": "STAT1", "log2FoldChange": 1.6 * scale, "pvalue": 5e-5, "padj": 0.002},
        {"gene": "IRF7", "log2FoldChange": 1.4 * scale, "pvalue": 1e-4, "padj": 0.004},
        {"gene": "RPL13A", "log2FoldChange": 0.05, "pvalue": 0.78, "padj": 0.92},
        {"gene": "TBP", "log2FoldChange": -0.04, "pvalue": 0.82, "padj": 0.95},
    ]
    if late:
        rows.extend(
            [
                {"gene": "DDX58", "log2FoldChange": 1.7 * scale, "pvalue": 7e-5, "padj": 0.003},
                {"gene": "IFIH1", "log2FoldChange": 1.5 * scale, "pvalue": 9e-5, "padj": 0.004},
            ]
        )
    else:
        rows.extend(
            [
                {"gene": "CXCL10", "log2FoldChange": 1.1 * scale, "pvalue": 0.002, "padj": 0.04},
                {"gene": "HPRT1", "log2FoldChange": 0.02, "pvalue": 0.91, "padj": 0.97},
            ]
        )
    return pd.DataFrame(rows)


def _project_rows(output_dir: str, harmonized_dir: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"field": "project_name", "value": "degora_demo", "what_to_enter": "Demo project name."},
            {"field": "topic", "value": "interferon response demo", "what_to_enter": "Demo topic."},
            {"field": "organism", "value": "Homo sapiens", "what_to_enter": "Demo organism."},
            {"field": "output_dir", "value": output_dir, "what_to_enter": "Demo result folder."},
            {"field": "harmonized_dir", "value": harmonized_dir, "what_to_enter": "Demo harmonized folder."},
            {"field": "min_studies", "value": 2, "what_to_enter": "Require at least two independent source units per scored gene."},
        ]
    )


def _contrast_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "study_id": "DEMO_A_4h",
                "source_unit_id": "DEMO_A",
                "source_path": "deg_tables/demo_ifn_a_4h.csv",
                "condition": "IFN-beta vs untreated",
                "time_h": 4,
                "cell_system": "demo epithelial cells",
                "species": "Homo sapiens",
                "assay_type": "RNA-seq",
                "pipeline": "DESeq2",
                "source_input_type": "author_deg_table",
                "table_scope": "full_results",
                "gene_column": "gene",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
                "include": "yes",
                "notes": "Synthetic demo contrast.",
            },
            {
                "study_id": "DEMO_A_12h",
                "source_unit_id": "DEMO_A",
                "source_path": "deg_tables/demo_ifn_a_12h.csv",
                "condition": "IFN-beta vs untreated",
                "time_h": 12,
                "cell_system": "demo epithelial cells",
                "species": "Homo sapiens",
                "assay_type": "RNA-seq",
                "pipeline": "DESeq2",
                "source_input_type": "author_deg_table",
                "table_scope": "full_results",
                "gene_column": "gene",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
                "include": "yes",
                "notes": "Same source_unit_id as DEMO_A_4h to demonstrate time-course collapse.",
            },
            {
                "study_id": "DEMO_B_6h",
                "source_unit_id": "DEMO_B",
                "source_path": "deg_tables/demo_ifn_b_6h.csv",
                "condition": "IFN-alpha vs untreated",
                "time_h": 6,
                "cell_system": "demo hepatocyte cells",
                "species": "Homo sapiens",
                "assay_type": "RNA-seq",
                "pipeline": "edgeR",
                "source_input_type": "author_deg_table",
                "table_scope": "full_results",
                "gene_column": "gene",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
                "include": "yes",
                "notes": "Synthetic independent source unit.",
            },
            {
                "study_id": "DEMO_B_24h",
                "source_unit_id": "DEMO_B",
                "source_path": "deg_tables/demo_ifn_b_24h.csv",
                "condition": "IFN-alpha vs untreated",
                "time_h": 24,
                "cell_system": "demo hepatocyte cells",
                "species": "Homo sapiens",
                "assay_type": "RNA-seq",
                "pipeline": "edgeR",
                "source_input_type": "author_deg_table",
                "table_scope": "full_results",
                "gene_column": "gene",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
                "include": "yes",
                "notes": "Same source_unit_id as DEMO_B_6h.",
            },
        ]
    )


def _gold_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"gene_symbol": "ISG15", "expected_direction": "up", "role": "demo_marker", "evidence_basis": "Synthetic demo.", "locked": "yes"},
            {"gene_symbol": "IFIT1", "expected_direction": "up", "role": "demo_marker", "evidence_basis": "Synthetic demo.", "locked": "yes"},
            {"gene_symbol": "MX1", "expected_direction": "up", "role": "demo_marker", "evidence_basis": "Synthetic demo.", "locked": "yes"},
            {"gene_symbol": "OAS1", "expected_direction": "up", "role": "demo_marker", "evidence_basis": "Synthetic demo.", "locked": "yes"},
        ]
    )


def _readme_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"step": 1, "instruction": "Run: degora validate degora_demo_config.xlsx"},
            {"step": 2, "instruction": "Run: degora run degora_demo_config.xlsx"},
            {"step": 3, "instruction": "Run: degora serve results/degora_scores.db"},
            {"step": 4, "instruction": "Search ISG15, IFIT1, MX1, or OAS1 in the browser."},
            {"step": 5, "instruction": "Open the Contrasts sheet to see how time points share source_unit_id."},
        ]
    )


def _advanced_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"setting": "score_version", "value": "degora_score_v1", "what_it_does": "Transparent prioritization score."},
            {"setting": "browser_port", "value": 8765, "what_it_does": "Local browser/API port."},
            {"setting": "collapse_independence_by", "value": "source_unit_id", "what_it_does": "Avoids overcounting time points."},
        ]
    )


def _guide_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"column": "source_unit_id", "meaning": "DEMO_A has two time points but counts as one independent source unit."},
            {"column": "source_unit_id", "meaning": "DEMO_B has two time points but counts as one independent source unit."},
            {"column": "degora_score", "meaning": "Prioritization score, not a probability."},
        ]
    )


def write_demo_workspace(path: str | Path, *, force: bool = False) -> dict[str, str]:
    """Write a complete synthetic demo workspace with source DEG tables and Excel config."""

    output_dir = Path(path)
    if output_dir.exists() and any(output_dir.iterdir()) and not force:
        raise FileExistsError(f"Demo folder is not empty: {output_dir}. Use --force to overwrite demo files.")

    deg_dir = output_dir / "deg_tables"
    deg_dir.mkdir(parents=True, exist_ok=True)
    _demo_rows(1.0).to_csv(deg_dir / "demo_ifn_a_4h.csv", index=False)
    _demo_rows(1.1, late=True).to_csv(deg_dir / "demo_ifn_a_12h.csv", index=False)
    _demo_rows(0.9).to_csv(deg_dir / "demo_ifn_b_6h.csv", index=False)
    _demo_rows(1.2, late=True).to_csv(deg_dir / "demo_ifn_b_24h.csv", index=False)

    config_path = output_dir / DEMO_CONFIG_NAME
    with pd.ExcelWriter(config_path, engine="openpyxl") as writer:
        write_note_sheet(writer, "README", _readme_rows())
        write_note_sheet(writer, "Project", _project_rows("results", "harmonized"))
        write_note_sheet(writer, "Contrasts", _contrast_rows())
        write_note_sheet(writer, "GoldPanel", _gold_rows())
        write_note_sheet(writer, "AdvancedSettings", _advanced_rows())
        write_note_sheet(writer, "ColumnGuide", _guide_rows())
    _autosize_workbook(config_path)

    readme_path = output_dir / "README.rst"
    readme_path.write_text(
        "\n".join(
            [
                "DEGORA Demo",
                "===========",
                "",
                "This folder contains tiny synthetic IFN-like DEG tables.",
                "",
                "Run:",
                "",
                ".. code-block:: bash",
                "",
                f"degora validate {DEMO_CONFIG_NAME}",
                f"degora run {DEMO_CONFIG_NAME}",
                "degora serve results/degora_scores.db",
                "",
                "Search for ISG15, IFIT1, MX1, or OAS1 in the browser.",
                "The demo has four contrast rows but two independent source units.",
                "",
            ]
        )
    )

    return {
        "demo_dir": str(output_dir.resolve()),
        "config": str(config_path.resolve()),
        "readme": str(readme_path.resolve()),
        "deg_tables": str(deg_dir.resolve()),
    }
