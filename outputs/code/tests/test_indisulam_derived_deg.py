from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.write_indisulam_derived_deg import (
    ContrastSpec,
    _derive_contrast,
    _read_gse223011_counts,
)


def test_gse223011_reader_maps_ensembl_symbols_and_collapses_duplicates(tmp_path: Path) -> None:
    workbook = tmp_path / "counts.xlsx"
    columns = [
        "Geneid",
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length",
        "../Human/3. Merged BAM/KURA_INDY-1_merged.bam",
        "../Human/3. Merged BAM/KURA_INDY-2_merged.bam",
        "../Human/3. Merged BAM/KURA_INDY-3_merged.bam",
        "../Human/3. Merged BAM/KURA_VC-1_merged.bam",
        "../Human/3. Merged BAM/KURA_VC-2_merged.bam",
        "../Human/3. Merged BAM/KURA_VC-3_merged.bam",
    ]
    frame = pd.DataFrame(
        [
            ["ENSG000001.1", "1", 1, 2, "+", 10, 10, 12, 11, 1, 1, 1],
            ["ENSG000001.2", "1", 1, 2, "+", 10, 3, 3, 3, 1, 1, 1],
            ["ENSG000002", "1", 1, 2, "+", 10, 7, 8, 9, 8, 8, 8],
        ],
        columns=columns,
    )
    frame.to_excel(workbook, sheet_name="humanCounts", index=False)

    counts, summary = _read_gse223011_counts(workbook, {"ENSG000001": "RBM39"})

    assert list(counts.index) == ["RBM39"]
    assert counts.loc["RBM39", "../Human/3. Merged BAM/KURA_INDY-1_merged.bam"] == 13
    assert summary["n_unmapped_genes"] == 1
    assert summary["n_gene_symbols_after_collapse"] == 1


def test_derive_contrast_returns_positive_lfc_for_treatment_increase(tmp_path: Path) -> None:
    counts = pd.DataFrame(
        {
            "ctrl1": [10, 100],
            "ctrl2": [11, 100],
            "trt1": [100, 100],
            "trt2": [110, 100],
        },
        index=["RBM39", "ACTB"],
    )
    spec = ContrastSpec(
        study_id="TEST",
        paper_id="P1",
        source_url="https://example.test",
        source_path=tmp_path / "out.csv",
        raw_source_key="test",
        species="Homo sapiens",
        cell_system="cell",
        condition="indisulam",
        duration_h="24",
        n_ctrl=2,
        n_treat=2,
        control_columns=("ctrl1", "ctrl2"),
        treat_columns=("trt1", "trt2"),
        output_name="out.csv",
    )

    deg = _derive_contrast(counts, spec)

    rbm39 = deg.loc[deg["gene_symbol"].eq("RBM39")].iloc[0]
    assert rbm39["log2FoldChange"] > 0
    assert 0 <= rbm39["pvalue"] <= 1
    assert 0 <= rbm39["padj"] <= 1
