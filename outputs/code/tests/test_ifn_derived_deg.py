from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.write_ifn_derived_deg import ContrastSpec, _derive_contrast


def test_ifn_derived_contrast_applies_low_count_filter(tmp_path: Path) -> None:
    counts = pd.DataFrame(
        {
            "ctrl1": [50, 100, 0],
            "ctrl2": [55, 100, 1],
            "trt1": [500, 100, 0],
            "trt2": [520, 100, 1],
        },
        index=["ISG15", "ACTB", "LOWCOUNT"],
    )
    spec = ContrastSpec(
        study_id="TEST",
        paper_id="P1",
        source_url="https://example.test",
        source_path=tmp_path / "out.csv",
        raw_source_key="test",
        species="Homo sapiens",
        cell_system="cell",
        condition="IFN",
        duration_h="4",
        n_ctrl=2,
        n_treat=2,
        control_columns=("ctrl1", "ctrl2"),
        treat_columns=("trt1", "trt2"),
        output_name="out.csv",
    )

    deg = _derive_contrast(counts, spec)

    assert "LOWCOUNT" not in set(deg["gene_symbol"])
    assert set(deg["gene_symbol"]) == {"ACTB", "ISG15"}
    assert deg["n_genes_before_low_count_filter"].iloc[0] == 3
    assert deg["n_genes_after_low_count_filter"].iloc[0] == 2
    assert deg["n_genes_removed_low_count_filter"].iloc[0] == 1
