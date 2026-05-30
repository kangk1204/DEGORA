from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.write_manuscript_package import FIGURE_MANIFESTS, TABLE_INPUTS, build_package


def _write_minimal_package(tmp_path: Path, manuscript_text: str) -> tuple[Path, Path, Path, Path]:
    manuscript_dir = tmp_path / "outputs" / "manuscript"
    tables_dir = tmp_path / "outputs" / "tables" / "manuscript"
    figures_dir = tmp_path / "outputs" / "figures" / "manuscript"
    manuscript_dir.mkdir(parents=True)
    tables_dir.mkdir(parents=True)
    figures_dir.mkdir(parents=True)
    manuscript = manuscript_dir / "main.md"
    manuscript.write_text(manuscript_text)

    pd.DataFrame(
        {
            "topic": ["IFN", "ER stress", "Heat shock", "Hypoxia"],
            "degora_method": ["degora_quality_weighted_score/quality_weighted_secondary"] * 4,
            "degora_recall_at_10": [0.35, 0.28, 0.38, 0.15],
            "degora_recall_at_50": [0.85, 0.83, 0.69, 0.65],
            "degora_recall_at_100": [0.90, 0.89, 0.75, 0.75],
            "best_non_degora_method": ["fisher/default"] * 4,
            "best_non_degora_recall_at_10": [0.1] * 4,
            "best_non_degora_recall_at_50": [0.1] * 4,
            "best_non_degora_recall_at_100": [0.1] * 4,
            "degora_top10": ["GENE"] * 4,
        }
    ).to_csv(tables_dir / "Table1_core_benchmark_summary.csv", index=False)
    pd.DataFrame(
        {
            "topic": ["IFN"],
            "mode": ["RNA-seq only"],
            "method_id": ["degora_quality_weighted_score"],
            "setting_id": ["quality_weighted_secondary"],
            "n_gene_scores": [100],
            "n_source_units_total": [2],
            "n_contrasts_total": [7],
            "recall_at_10": [0.35],
            "recall_at_20": [0.45],
            "recall_at_50": [0.85],
            "recall_at_100": [0.90],
            "direction_recall_at_100": [0.90],
            "top10": ["GENE"],
        }
    ).to_csv(tables_dir / "Table2_cross_platform_summary.csv", index=False)
    for name in TABLE_INPUTS:
        path = tables_dir / name
        if not path.exists():
            path.write_text("{}\n" if path.suffix == ".json" else "x\n")
    for name in FIGURE_MANIFESTS:
        path = figures_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}\n")
    return manuscript, tables_dir, figures_dir, manuscript_dir


def test_manuscript_package_builds_manifest_and_validation(tmp_path: Path) -> None:
    text = "## Results\nIFN 0.90 ER stress 0.89 Heat shock 0.75 Hypoxia 0.75\n\n## Data and Code Availability\n"
    manuscript, tables_dir, figures_dir, output_dir = _write_minimal_package(tmp_path, text)

    summary = build_package(
        manuscript=manuscript,
        tables_dir=tables_dir,
        figures_dir=figures_dir,
        output_dir=output_dir,
        command="make paper",
    )

    manifest = json.loads(Path(summary["manifest"]).read_text())
    assert manifest["failures"] == []
    assert manifest["table_claims"]["table1_headline_recall_at_100"][0]["degora_recall_at_100"] == "0.90"
    assert (output_dir / "degora_paper_manifest.json.source").exists()
    assert (output_dir / "degora_paper_validation.txt").exists()


def test_manuscript_package_blocks_reference_placeholders_and_stale_hypoxia(tmp_path: Path) -> None:
    text = "## Results\nHypoxia recall@100 was 0.80. [REFERENCE NEEDED]\n\n## Data and Code Availability\n"
    manuscript, tables_dir, figures_dir, output_dir = _write_minimal_package(tmp_path, text)

    with pytest.raises(ValueError, match="REFERENCE NEEDED"):
        build_package(
            manuscript=manuscript,
            tables_dir=tables_dir,
            figures_dir=figures_dir,
            output_dir=output_dir,
            command="make paper",
        )
