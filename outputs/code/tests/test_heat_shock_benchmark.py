from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.write_heat_shock_benchmark import (
    ContrastSpec,
    _derive_logcpm_welch_contrast,
    _read_ncbi_geneid_counts,
    _read_ncbi_geneid_symbol_map,
    _write_gold_panel,
)


def test_read_ncbi_geneid_counts_maps_symbols_and_collapses_duplicates(tmp_path: Path) -> None:
    annot_path = tmp_path / "annot.tsv"
    pd.DataFrame(
        {
            "GeneID": ["1", "2", "3", "4"],
            "Symbol": ["HSPA1A", "HSPA1A", "TRR-TCG5-1", ""],
            "GeneType": ["protein-coding", "protein-coding", "tRNA", "protein-coding"],
        }
    ).to_csv(annot_path, sep="\t", index=False)
    counts_path = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "GeneID": ["1", "2", "3", "4"],
            "ctrl1": [10, 3, 5, 5],
            "ctrl2": [11, 3, 5, 5],
            "trt1": [100, 7, 5, 5],
            "trt2": [110, 8, 5, 5],
        }
    ).to_csv(counts_path, sep="\t", index=False)

    mapping = _read_ncbi_geneid_symbol_map(annot_path)
    counts, summary = _read_ncbi_geneid_counts(counts_path, mapping)

    assert list(counts.index) == ["HSPA1A"]
    assert counts.loc["HSPA1A", "ctrl1"] == 13
    assert summary["n_unmapped_genes"] == 2
    assert summary["n_gene_symbols_after_collapse"] == 1


def test_logcpm_welch_contrast_returns_positive_lfc_for_heat_shock_increase(tmp_path: Path) -> None:
    counts = pd.DataFrame(
        {
            "ctrl1": [10, 100, 0],
            "ctrl2": [11, 100, 1],
            "hs1": [100, 100, 0],
            "hs2": [110, 100, 1],
        },
        index=["HSPA1A", "ACTB", "LOWCOUNT"],
    )
    spec = ContrastSpec(
        study_id="TEST",
        paper_id="P1",
        source_url="https://example.test",
        source_path=tmp_path / "out.csv",
        raw_source_key="test",
        species="Homo sapiens",
        cell_system="cell",
        condition="heat shock",
        duration_h="1",
        n_ctrl=2,
        n_treat=2,
        control_columns=("ctrl1", "ctrl2"),
        treat_columns=("hs1", "hs2"),
        output_name="out.csv",
    )

    deg = _derive_logcpm_welch_contrast(counts, spec)

    hspa1a = deg.loc[deg["gene_symbol"].eq("HSPA1A")].iloc[0]
    assert hspa1a["log2FoldChange"] > 0
    assert 0 <= hspa1a["pvalue"] <= 1
    assert 0 <= hspa1a["padj"] <= 1
    assert "LOWCOUNT" not in set(deg["gene_symbol"])
    assert hspa1a["n_genes_before_low_count_filter"] == 3
    assert hspa1a["n_genes_after_low_count_filter"] == 2


def test_gold_panel_is_direction_locked_before_scoring(tmp_path: Path) -> None:
    gold_path = tmp_path / "heat_shock_gold.csv"

    _write_gold_panel(gold_path, "test command")

    gold = pd.read_csv(gold_path)
    assert {"gene_symbol", "expected_direction", "locked"}.issubset(gold.columns)
    assert gold["expected_direction"].eq("up").all()
    assert gold["locked"].eq("yes").all()
    assert "HSPA1A" in set(gold["gene_symbol"])
