"""Conservative fallback derivation for selected public expression matrices.

The functions here implement DEGORA's labeled fallback regime. They do
not infer biological groups: callers must supply disjoint control and treatment
sample columns, with treatment-minus-control as the effect direction.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats

from .derived_counts import attach_low_count_filter_metadata, low_count_filter_mask, low_count_filter_summary
from .formula_safety import formula_guard_metadata, neutralize_formula_text
from .harmonize import _clean_gene_symbol
from .provenance import apply_default_file_mode, write_source_sidecar


SUPPORTED_MATRIX_ROLES = frozenset({"count_matrix", "normalized_expression_matrix"})
MAX_GROUP_SAMPLES = 100
# Welch's t-test divides by the within-group variances. When both groups of a
# gene hold identical replicate values - rounded array intensities, a gene at a
# detection floor, two of two samples that agree to every digit - the variance is
# zero, t is infinite and scipy reports p = 0.0: a 1.07-fold change was ranked as
# the most significant gene in a corpus, above a 32-fold change, because its two
# replicates matched. Variances below a data-driven floor are raised to it. The
# floor is a low quantile of the positive within-group variances of the same
# matrix, so it scales with the data and touches only the genes whose variance
# is unbelievably small for that matrix.
VARIANCE_FLOOR_QUANTILE = 0.01
VARIANCE_FLOOR_MINIMUM = 1e-8


def _read_matrix(path: Path, *, sheet_name: str | int | None = None) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".xlsx"):
        return pd.read_excel(path, sheet_name=0 if sheet_name in (None, "") else sheet_name)
    separator = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return pd.read_csv(path, sep=separator, low_memory=False)


def _clean_gene_ids(values: pd.Series) -> pd.Series:
    # The one gene-identity rule DEGORA has (harmonize.canonical_gene_symbol):
    # a private copy here stripped `.\d+` from every label and merged NKX2.5
    # with NKX2.1 in matrices exactly as it once did in result tables.
    return _clean_gene_symbol(values)


def welch_with_variance_floor(
    treatment: np.ndarray,
    control: np.ndarray,
    *,
    floor_quantile: float = VARIANCE_FLOOR_QUANTILE,
    floor_minimum: float = VARIANCE_FLOOR_MINIMUM,
) -> dict[str, Any]:
    """Row-wise Welch t-test with a data-driven floor on the within-group variances.

    Rows whose variances are above the floor get exactly scipy's Welch result.
    Returns the p-values, the Welch t statistics, the degrees of freedom, the
    floor used and how many rows it touched.
    """

    treatment = np.asarray(treatment, dtype=float)
    control = np.asarray(control, dtype=float)
    n_t = treatment.shape[1]
    n_c = control.shape[1]
    mean_t = treatment.mean(axis=1)
    mean_c = control.mean(axis=1)
    var_t = treatment.var(axis=1, ddof=1)
    var_c = control.var(axis=1, ddof=1)
    positive = np.concatenate([var_t[np.isfinite(var_t) & (var_t > 0)], var_c[np.isfinite(var_c) & (var_c > 0)]])
    floor = float(max(np.quantile(positive, floor_quantile), floor_minimum)) if positive.size else float(floor_minimum)
    floored = (var_t < floor) | (var_c < floor)
    var_t = np.maximum(var_t, floor)
    var_c = np.maximum(var_c, floor)
    se_t = var_t / n_t
    se_c = var_c / n_c
    se2 = se_t + se_c
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = (mean_t - mean_c) / np.sqrt(se2)
        df = se2**2 / (se_t**2 / (n_t - 1) + se_c**2 / (n_c - 1))
    pvalue = 2.0 * stats.t.sf(np.abs(t_stat), df)
    return {
        "pvalue": pvalue,
        "t": t_stat,
        "df": df,
        "variance_floor": floor,
        "n_variance_floored": int(floored.sum()),
        "n_zero_variance_rows": int(((treatment.var(axis=1, ddof=1) == 0) | (control.var(axis=1, ddof=1) == 0)).sum()),
    }


def validate_sample_groups(
    control_samples: Iterable[str],
    treatment_samples: Iterable[str],
    *,
    expected_control_count: int | None = None,
    expected_treatment_count: int | None = None,
    max_group_samples: int | None = None,
    context: str = "Welch fallback",
) -> tuple[list[str], list[str]]:
    """Return canonical, disjoint sample lists and enforce declared group sizes.

    The function validates identifiers only. It cannot establish that a public
    sample or wet-lab tube was biologically labeled correctly; that requires
    sample-level identity and expression-QC evidence outside a DEG table.
    """

    control = [str(value).strip() for value in control_samples if str(value).strip()]
    treatment = [str(value).strip() for value in treatment_samples if str(value).strip()]
    if len(control) < 2 or len(treatment) < 2:
        raise ValueError(f"{context} requires at least two control and two treatment sample columns")
    if max_group_samples is not None and (len(control) > max_group_samples or len(treatment) > max_group_samples):
        raise ValueError(f"{context}: each sample group may contain at most {max_group_samples} columns")
    if len(set(control)) != len(control) or len(set(treatment)) != len(treatment):
        raise ValueError(f"{context}: sample columns must not be duplicated within a group")
    overlap = sorted(set(control).intersection(treatment))
    if overlap:
        raise ValueError(f"{context}: control and treatment groups must be disjoint: " + ", ".join(overlap))

    for label, expected, actual in (
        ("n_ctrl", expected_control_count, len(control)),
        ("n_treat", expected_treatment_count, len(treatment)),
    ):
        if expected is None:
            continue
        try:
            numeric = float(expected)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{context}: declared {label} must be a positive whole number, got {expected!r}"
            ) from exc
        if not np.isfinite(numeric) or numeric <= 0 or not numeric.is_integer():
            raise ValueError(f"{context}: declared {label} must be a positive whole number, got {expected!r}")
        if int(numeric) != actual:
            raise ValueError(f"{context}: declared {label}={int(numeric)} but selected {actual} sample columns")
    return control, treatment


def _validated_samples(control_samples: Iterable[str], treatment_samples: Iterable[str]) -> tuple[list[str], list[str]]:
    """Compatibility wrapper for the bounded discovery fallback."""

    return validate_sample_groups(
        control_samples,
        treatment_samples,
        max_group_samples=MAX_GROUP_SAMPLES,
        context="Welch fallback",
    )


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).clip(0.0, 1.0)
    n = len(values)
    if n == 0:
        return pd.Series(dtype=float, index=pvalues.index)
    order = np.argsort(values.to_numpy(dtype=float))
    ranked = values.to_numpy(dtype=float)[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty(n, dtype=float)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(output, index=pvalues.index)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        neutralize_formula_text(frame).to_csv(temporary, index=False, lineterminator="\n")
        apply_default_file_mode(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def derive_welch_deg(
    matrix_path: str | Path,
    output_path: str | Path,
    *,
    role: str,
    gene_column: str,
    control_samples: Iterable[str],
    treatment_samples: Iterable[str],
    normalized_scale: str | None = None,
    sheet_name: str | int | None = None,
    command: str = "degora discovery analyze",
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Derive a full DEG-like table using the documented Welch fallback."""

    if role not in SUPPORTED_MATRIX_ROLES:
        raise ValueError("role must be count_matrix or normalized_expression_matrix")
    source = Path(matrix_path).resolve()
    output = Path(output_path).resolve()
    if not source.exists():
        raise FileNotFoundError(f"matrix does not exist: {source}")
    control, treatment = _validated_samples(control_samples, treatment_samples)
    if not str(gene_column).strip():
        raise ValueError("gene_column is required")
    if str(gene_column).strip().upper() == "ID_REF":
        raise ValueError("ID_REF is a probe identifier; map probes to gene symbols before running the fallback")

    matrix = _read_matrix(source, sheet_name=sheet_name)
    required = [gene_column, *control, *treatment]
    missing = [column for column in required if column not in matrix.columns]
    if missing:
        raise ValueError("matrix is missing required columns: " + ", ".join(missing))
    genes = _clean_gene_ids(matrix[gene_column])
    values = matrix[control + treatment].apply(pd.to_numeric, errors="coerce")
    finite = np.isfinite(values.to_numpy(dtype=float)).all(axis=1)
    valid = genes.notna().to_numpy() & finite
    values = values.loc[valid].copy()
    genes = genes.loc[valid]
    if values.empty:
        raise ValueError("matrix has no complete finite rows after gene and sample validation")

    values.insert(0, "__gene__", genes.astype(str).to_numpy())
    if role == "count_matrix":
        sample_values = values[control + treatment]
        if (sample_values < 0).any().any():
            raise ValueError("count matrix contains negative values")
        numeric_counts = sample_values.to_numpy(dtype=float)
        integer_like = np.isclose(numeric_counts, np.rint(numeric_counts), rtol=0.0, atol=1e-6)
        if not bool(integer_like.all()):
            raise ValueError("count matrix must contain non-negative integer-like raw counts")
        collapsed = values.groupby("__gene__", sort=True)[control + treatment].sum(min_count=1)
        # Library sizes come from the full count matrix, before the expression filter
        # removes low-count genes. Summing the filtered matrix instead would make each
        # sample's CPM denominator depend on how much low-count mass that particular
        # sample lost, which is a per-sample shift in log space rather than the common
        # scaling logCPM is defined as.
        libraries = collapsed.sum(axis=0).replace(0.0, np.nan)
        expressed = low_count_filter_mask(collapsed)
        filter_summary = low_count_filter_summary(collapsed, expressed)
        collapsed = collapsed.loc[expressed]
        if collapsed.empty:
            raise ValueError("no genes remain after the predeclared low-count filter")
        transformed = np.log2(collapsed.divide(libraries, axis=1).mul(1_000_000.0).fillna(0.0) + 1.0)
        source_input_type = "derived_count_table"
        pipeline = "logCPM_Welch_derived_from_public_counts"
        normalization = "logCPM_from_raw_counts"
    else:
        collapsed = values.groupby("__gene__", sort=True)[control + treatment].median()
        if normalized_scale == "linear":
            if (collapsed < 0).any().any():
                raise ValueError("linear normalized scale was selected but the matrix contains negative values")
            transformed = np.log2(collapsed + 1.0)
            normalization = "public_normalized_matrix_log2x_plus_1"
        elif normalized_scale == "log2":
            transformed = collapsed
            normalization = "public_normalized_matrix_confirmed_log2_scale"
        else:
            raise ValueError("normalized_scale must explicitly be log2 or linear for normalized expression matrices")
        filter_summary = {}
        source_input_type = "normalized_expression_matrix"
        pipeline = "Welch_confirmed_log2_expression_matrix"

    control_values = transformed[control]
    treatment_values = transformed[treatment]
    log2fc = treatment_values.mean(axis=1) - control_values.mean(axis=1)
    test = welch_with_variance_floor(
        treatment_values.to_numpy(dtype=float),
        control_values.to_numpy(dtype=float),
    )
    pvalue = pd.Series(test["pvalue"], index=transformed.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    result = pd.DataFrame(
        {
            "gene_symbol": transformed.index.astype(str),
            "log2FoldChange": log2fc.to_numpy(dtype=float),
            "pvalue": pvalue.to_numpy(dtype=float),
        }
    )
    result["padj"] = _bh_adjust(result["pvalue"])
    result["source_input_type"] = source_input_type
    result["pipeline"] = pipeline
    result["normalization"] = normalization
    result["effect_direction"] = "treatment_minus_control"
    variance_summary = {
        "welch_variance_floor": float(test["variance_floor"]),
        "welch_variance_floor_quantile": float(VARIANCE_FLOOR_QUANTILE),
        "n_genes_variance_floored": int(test["n_variance_floored"]),
        "n_genes_zero_within_group_variance": int(test["n_zero_variance_rows"]),
        "welch_variance_floor_rule": (
            "within-group variances below the 1st percentile of positive within-group variances "
            "(or 1e-8) are raised to that floor before Welch's t; identical replicates no longer "
            "yield p = 0"
        ),
    }
    if filter_summary:
        result = attach_low_count_filter_metadata(result, filter_summary)
    result = result.sort_values(["pvalue", "gene_symbol"]).reset_index(drop=True)
    _atomic_csv(result, output)

    provenance = {
        "generator": "degora.reanalysis.derive_welch_deg",
        "role": role,
        "source_input_type": source_input_type,
        "pipeline": pipeline,
        "normalization": normalization,
        "gene_column": gene_column,
        "control_samples": control,
        "treatment_samples": treatment,
        "effect_direction": "treatment_minus_control",
        "normalized_scale": normalized_scale or "not_applicable",
        **variance_summary,
        **filter_summary,
        **(metadata or {}),
        **formula_guard_metadata(),
    }
    write_source_sidecar(output, command, inputs=[source], metadata=provenance)
    return {
        "matrix_path": str(source),
        "output_path": str(output),
        "source_input_type": source_input_type,
        "pipeline": pipeline,
        "normalization": normalization,
        "gene_column": "gene_symbol",
        "lfc_column": "log2FoldChange",
        "p_column": "pvalue",
        "padj_column": "padj",
        "n_input_rows": int(len(matrix)),
        "n_valid_rows": int(len(values)),
        "n_gene_rows": int(len(result)),
        "n_ctrl": len(control),
        "n_treat": len(treatment),
        "control_samples": control,
        "treatment_samples": treatment,
        "effect_direction": "treatment_minus_control",
        "normalized_scale": normalized_scale or "not_applicable",
        **variance_summary,
        **filter_summary,
    }


__all__ = ["SUPPORTED_MATRIX_ROLES", "derive_welch_deg", "validate_sample_groups", "welch_with_variance_floor"]
