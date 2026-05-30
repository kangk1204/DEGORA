from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts.derive_microarray_deg import derive_microarray_deg


def test_derive_microarray_deg_collapses_probe_rows_and_writes_provenance(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.csv"
    pd.DataFrame(
        {
            "probe_id": ["p1", "p2", "p3"],
            "gene_symbol": ["RBM39", "RBM39", "ACTB"],
            "ctrl1": [5.0, 5.0, 10.0],
            "ctrl2": [5.2, 5.1, 10.1],
            "drug1": [8.0, 5.4, 10.0],
            "drug2": [8.1, 5.5, 10.1],
        }
    ).to_csv(matrix, index=False)
    output = tmp_path / "microarray_deg.csv"
    summary_path = tmp_path / "summary.json"

    summary = derive_microarray_deg(
        matrix,
        output,
        gene_column="gene_symbol",
        probe_column="probe_id",
        control_samples=["ctrl1", "ctrl2"],
        treatment_samples=["drug1", "drug2"],
        command="test command",
        summary_path=summary_path,
        metadata={"platform": "GPL570", "normalization": "RMA/log2"},
    )

    deg = pd.read_csv(output)
    rbm39 = deg.loc[deg["gene_symbol"].eq("RBM39")].iloc[0]
    assert summary["n_gene_rows"] == 2
    assert summary["n_collapsed_probe_rows"] == 1
    assert rbm39["probe_id"] == "p1"
    assert rbm39["log2FoldChange"] > 0
    assert rbm39["assay_type"] == "microarray"
    assert rbm39["platform"] == "GPL570"
    assert output.with_suffix(".csv.source").exists()
    assert summary_path.with_suffix(".json.source").exists()


def test_derive_microarray_deg_rejects_single_replicate_welch_fallback(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.csv"
    pd.DataFrame(
        {
            "gene_symbol": ["RBM39", "ACTB"],
            "ctrl1": [5.0, 10.0],
            "drug1": [8.0, 10.1],
        }
    ).to_csv(matrix, index=False)

    with pytest.raises(ValueError, match="at least two control and two treatment"):
        derive_microarray_deg(
            matrix,
            tmp_path / "microarray_deg.csv",
            gene_column="gene_symbol",
            control_samples=["ctrl1"],
            treatment_samples=["drug1"],
            command="test command",
        )
