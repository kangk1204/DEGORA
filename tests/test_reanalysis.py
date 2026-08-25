from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from degora.formula_safety import restore_formula_text_if_marked
from degora.reanalysis import derive_welch_deg, validate_sample_groups


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene": ["TP53", "CDKN1A", "LOW"],
            "c1": [100, 50, 0],
            "c2": [110, 55, 1],
            "t1": [500, 200, 0],
            "t2": [520, 210, 1],
        }
    )


def test_count_fallback_uses_logcpm_welch_and_treatment_minus_control(tmp_path: Path) -> None:
    source = tmp_path / "counts.csv"
    output = tmp_path / "derived.csv"
    _matrix().to_csv(source, index=False)

    summary = derive_welch_deg(
        source,
        output,
        role="count_matrix",
        gene_column="gene",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
    )

    result = pd.read_csv(output)
    assert summary["pipeline"] == "logCPM_Welch_derived_from_public_counts"
    assert summary["source_input_type"] == "derived_count_table"
    assert set(result["gene_symbol"]) == {"TP53", "CDKN1A"}
    assert result.set_index("gene_symbol").loc["TP53", "log2FoldChange"] > 0
    assert result["pvalue"].between(0, 1).all()
    assert result["padj"].between(0, 1).all()
    sidecar = json.loads(Path(str(output) + ".provenance.json").read_text())
    assert sidecar["metadata"]["effect_direction"] == "treatment_minus_control"


def test_normalized_expression_fallback_collapses_duplicate_genes_by_median(tmp_path: Path) -> None:
    source = tmp_path / "matrix.tsv"
    output = tmp_path / "derived.csv"
    pd.DataFrame(
        {
            "symbol": ["A", "A", "B"],
            "c1": [1.0, 3.0, 8.0],
            "c2": [1.1, 3.1, 8.2],
            "t1": [5.0, 7.0, 8.0],
            "t2": [5.1, 7.1, 8.1],
        }
    ).to_csv(source, sep="\t", index=False)

    summary = derive_welch_deg(
        source,
        output,
        role="normalized_expression_matrix",
        gene_column="symbol",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        normalized_scale="log2",
    )
    result = pd.read_csv(output)
    assert summary["source_input_type"] == "normalized_expression_matrix"
    assert len(result) == 2
    assert result.set_index("gene_symbol").loc["A", "log2FoldChange"] == pytest.approx(4.0)


def test_fallback_csv_guards_formula_genes_and_restores_them_with_provenance(tmp_path: Path) -> None:
    source = tmp_path / "matrix.tsv"
    output = tmp_path / "derived.csv"
    pd.DataFrame(
        {
            "symbol": ["=BAD()", "'=LITERAL", "TP53"],
            "c1": [1.0, 2.0, 3.0],
            "c2": [1.2, 2.2, 3.2],
            "t1": [4.0, 5.0, 6.0],
            "t2": [4.2, 5.2, 6.2],
        }
    ).to_csv(source, sep="\t", index=False)

    derive_welch_deg(
        source,
        output,
        role="normalized_expression_matrix",
        gene_column="symbol",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        normalized_scale="log2",
        metadata={"csv_formula_guard": "caller_must_not_override_reserved_scheme"},
    )

    guarded = pd.read_csv(output)
    assert {"'=BAD()", "''=LITERAL"}.issubset(set(guarded["gene_symbol"]))
    restored = restore_formula_text_if_marked(guarded, output)
    assert {"=BAD()", "'=LITERAL"}.issubset(set(restored["gene_symbol"]))
    provenance = json.loads(Path(str(output) + ".provenance.json").read_text())
    assert provenance["metadata"]["csv_formula_guard"] == "reversible_apostrophe_prefix_v1"


def test_normalized_expression_requires_explicit_scale(tmp_path: Path) -> None:
    source = tmp_path / "matrix.tsv"
    output = tmp_path / "derived.csv"
    _matrix().to_csv(source, sep="\t", index=False)
    with pytest.raises(ValueError, match="normalized_scale"):
        derive_welch_deg(
            source,
            output,
            role="normalized_expression_matrix",
            gene_column="gene",
            control_samples=["c1", "c2"],
            treatment_samples=["t1", "t2"],
        )


def test_linear_normalized_expression_is_log2_transformed_before_effect_size(tmp_path: Path) -> None:
    source = tmp_path / "matrix.tsv"
    output = tmp_path / "derived.csv"
    pd.DataFrame(
        {
            "gene": ["A", "B"],
            "c1": [3.0, 7.0],
            "c2": [3.0, 7.0],
            "t1": [15.0, 7.0],
            "t2": [15.0, 7.0],
        }
    ).to_csv(source, sep="\t", index=False)
    summary = derive_welch_deg(
        source,
        output,
        role="normalized_expression_matrix",
        gene_column="gene",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        normalized_scale="linear",
    )
    result = pd.read_csv(output).set_index("gene_symbol")
    assert result.loc["A", "log2FoldChange"] == pytest.approx(2.0)
    assert summary["normalized_scale"] == "linear"


@pytest.mark.parametrize(
    ("control", "treatment", "message"),
    [
        (["c1"], ["t1", "t2"], "at least two"),
        (["c1", "c2"], ["c2", "t2"], "disjoint"),
        (["c1", "c1"], ["t1", "t2"], "duplicated within a group"),
        (["c1", "c2"], ["t1", "t1"], "duplicated within a group"),
    ],
)
def test_fallback_rejects_invalid_group_design(tmp_path: Path, control, treatment, message) -> None:
    source = tmp_path / "counts.csv"
    _matrix().to_csv(source, index=False)
    with pytest.raises(ValueError, match=message):
        derive_welch_deg(
            source,
            tmp_path / "derived.csv",
            role="count_matrix",
            gene_column="gene",
            control_samples=control,
            treatment_samples=treatment,
        )


def test_sample_validator_rejects_declared_count_mismatch() -> None:
    with pytest.raises(ValueError, match="declared n_ctrl=3 but selected 2"):
        validate_sample_groups(
            ["c1", "c2"],
            ["t1", "t2"],
            expected_control_count=3,
            expected_treatment_count=2,
            context="TEST_CONTRAST",
        )


def test_sample_validator_rejects_nonnumeric_declared_count() -> None:
    with pytest.raises(ValueError, match="declared n_ctrl must be a positive whole number"):
        validate_sample_groups(
            ["c1", "c2"],
            ["t1", "t2"],
            expected_control_count="two",
            expected_treatment_count=2,
            context="TEST_CONTRAST",
        )


def test_fallback_rejects_unannotated_microarray_probe_ids(tmp_path: Path) -> None:
    source = tmp_path / "matrix.csv"
    _matrix().rename(columns={"gene": "ID_REF"}).to_csv(source, index=False)
    with pytest.raises(ValueError, match="probe identifier"):
        derive_welch_deg(
            source,
            tmp_path / "derived.csv",
            role="normalized_expression_matrix",
            gene_column="ID_REF",
            control_samples=["c1", "c2"],
            treatment_samples=["t1", "t2"],
        )


def test_count_fallback_rejects_fractional_values_mislabeled_as_raw_counts(tmp_path: Path) -> None:
    source = tmp_path / "fractional_counts.csv"
    frame = _matrix().astype({"c1": float})
    frame.loc[0, "c1"] = 100.25
    frame.to_csv(source, index=False)
    with pytest.raises(ValueError, match="integer-like raw counts"):
        derive_welch_deg(
            source,
            tmp_path / "derived.csv",
            role="count_matrix",
            gene_column="gene",
            control_samples=["c1", "c2"],
            treatment_samples=["t1", "t2"],
        )


def test_fallback_requires_resolved_matrix_type(tmp_path: Path) -> None:
    source = tmp_path / "matrix.csv"
    _matrix().to_csv(source, index=False)
    with pytest.raises(ValueError, match="count_matrix or normalized_expression_matrix"):
        derive_welch_deg(
            source,
            tmp_path / "derived.csv",
            role="unknown_matrix",
            gene_column="gene",
            control_samples=["c1", "c2"],
            treatment_samples=["t1", "t2"],
        )
