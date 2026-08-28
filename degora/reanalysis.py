"""Conservative fallback derivation for selected public expression matrices.

The functions here implement DEGORA's labeled fallback regime. They do
not infer biological groups: callers must supply disjoint control and treatment
sample columns, with treatment-minus-control as the effect direction.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .derived_counts import (
    attach_low_count_filter_metadata,
    low_count_filter_mask,
    low_count_filter_summary,
)
from .formula_safety import formula_guard_metadata, neutralize_formula_text
from .harmonize import WORKBOOK_SUFFIXES, _clean_gene_symbol, _read_excel_any
from .provenance import apply_default_file_mode, write_source_sidecar


# Estimated counts are what Salmon, RSEM and kallisto write: fractional by
# design, counts in every other respect, normalised like counts (as tximport
# hands them to DESeq2) without the integer test that raw counts must pass.
SUPPORTED_MATRIX_ROLES = frozenset({"count_matrix", "estimated_count_matrix", "normalized_expression_matrix"})
COUNT_MATRIX_ROLES = frozenset({"count_matrix", "estimated_count_matrix"})
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


# The candidate inspector reads at most this much text; the sniffer below must
# see the same bytes it saw or the two can disagree about the delimiter.
SNIFF_TEXT_BYTES = 5 * 1024 * 1024
ENSEMBL_SYMBOL_SUFFIX_TRANSFORM = "ensembl_gene_symbol_suffix_v1"
ENSEMBL_PREFIX_TRANSFORM = "ensembl_gene_prefix_v1"
ENSEMBL_SYMBOL_COMPOSITE_VALUE_RE = re.compile(
    r"^(?P<ensembl>ENS[A-Z]*G\d+(?:\.\d+)?)_"
    r"(?P<symbol>(?=[A-Za-z0-9_.@-]*[A-Za-z])[A-Za-z0-9][A-Za-z0-9_.@-]*)$"
)
ENSEMBL_GENE_ID_RE = re.compile(r"^ENS[A-Z]*G\d+(?:\.\d+)?$")


def sniff_delimited_separator(path: Path, *, matrix_table: bool = False) -> str:
    """Pick the delimiter exactly the way candidate inspection did.

    Inspection tries tab, comma, semicolon and a quoted whitespace table, then
    prefers a parse with a safely recognizable gene header. Reading the same
    file back by its extension meant a
    comma-delimited ``.txt`` (or a semicolon ``.csv``) that inspection had
    presented as usable could never be activated: the columns inspection
    offered came back as one glued-together name.
    """

    opener = gzip.open if "".join(path.suffixes).lower().endswith(".gz") else open
    try:
        with opener(path, "rb") as handle:
            payload = handle.read(SNIFF_TEXT_BYTES)
    except (OSError, EOFError):
        # An unreadable file is reported by the pandas read itself.
        return "\t"
    lines = payload.decode("utf-8-sig", "replace").splitlines()
    if matrix_table:
        # Upstream-matrix inspection numbers headers inside the data table: it
        # drops GEO metadata/comments and blank lines, and starts after the
        # explicit series-matrix marker when one is present. Sniffing the raw
        # preamble can otherwise select a delimiter the inspected table never
        # used.
        for index, line in enumerate(lines):
            if line.strip().lower() == "!series_matrix_table_begin":
                lines = lines[index + 1 :]
                break
        lines = [line for line in lines if line.strip() and not line.startswith("!")]
    lines = lines[:80]
    from .discovery import (
        _author_delimited_parse_score,
        _name_upstream_row_label_column,
        is_gene_identifier_header,
    )

    best: tuple[str, tuple[int, ...]] = ("\t", (0, 0, 0))
    for separator in ("\t", ",", ";", " "):
        try:
            parsed = list(csv.reader(lines, delimiter=separator, skipinitialspace=separator == " "))
        except csv.Error:
            continue
        if matrix_table:
            width = max((len(row) for row in parsed[:10]), default=0)
            gene_header_width = 0
            for index, row in enumerate(parsed[:12]):
                repaired = _name_upstream_row_label_column(row, parsed[index + 1 : index + 21])
                if any(is_gene_identifier_header(value) for value in repaired):
                    gene_header_width = max(gene_header_width, len(repaired))
            score = (int(gene_header_width > 0), gene_header_width, width)
        else:
            score = _author_delimited_parse_score(parsed)
        if score > best[1]:
            best = (separator, score)
    return best[0]


def _safe_r_row_names(values: pd.Series | pd.Index) -> bool:
    labels = pd.Series(values, dtype="string").str.strip()
    nonempty = labels.notna() & labels.ne("")
    if len(labels) < 2 or float(nonempty.mean()) < 0.9:
        return False
    # Numeric row numbers are ambiguous with ranks or exported frame indices.
    numeric = pd.to_numeric(labels[nonempty], errors="coerce").notna()
    return not bool(numeric.any())


def _numeric_matrix_column_count(frame: pd.DataFrame, *, exclude: set[Any] | None = None) -> int:
    excluded = exclude or set()
    return sum(
        pd.to_numeric(frame[column], errors="coerce").notna().mean() >= 0.7
        for column in frame.columns
        if column not in excluded
    )


def _recover_r_row_name_column(frame: pd.DataFrame) -> pd.DataFrame:
    """Materialize the same safe row_name column that inspection reported."""

    if frame.empty or "row_name" in frame.columns:
        return frame
    if len(frame.columns):
        first = frame.columns[0]
        first_name = str(first).strip()
        if (
            (not first_name or re.fullmatch(r"Unnamed:\s*\d+", first_name, re.I))
            and _safe_r_row_names(frame[first])
            and _numeric_matrix_column_count(frame, exclude={first}) >= 4
        ):
            return frame.rename(columns={first: "row_name"})
    # pandas infers an index when every data row has one more field than the
    # header. This is how write.table(row.names=TRUE) appears when its leading
    # header cell was omitted entirely.
    if (
        not isinstance(frame.index, pd.RangeIndex)
        and frame.index.nlevels == 1
        and _safe_r_row_names(frame.index)
        and _numeric_matrix_column_count(frame) >= 4
    ):
        recovered = frame.copy()
        recovered.insert(0, "row_name", recovered.index.astype(str))
        return recovered.reset_index(drop=True)
    return frame


def read_matrix_frame(
    path: str | Path,
    *,
    sheet_name: str | int | None = None,
    header_row: int | None = None,
    nrows: int | None = None,
) -> pd.DataFrame:
    """Read the exact matrix table selected during discovery inspection.

    ``header_row`` is 1-based and is interpreted inside the matrix data table.
    For delimited files this is after blank lines and GEO ``!`` metadata have
    been removed, matching :func:`discovery.inspect_upstream_bytes`.
    """

    path = Path(path)
    if isinstance(header_row, bool):
        raise TypeError("header_row must be a positive 1-based row number")
    try:
        resolved_header = 1 if header_row in (None, "") else int(header_row)
    except (TypeError, ValueError) as exc:
        raise ValueError("header_row must be a positive 1-based row number") from exc
    if resolved_header < 1:
        raise ValueError("header_row must be a positive 1-based row number")
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(WORKBOOK_SUFFIXES):
        # .xls and gzipped workbooks were read as CSV here and died on a decode
        # error, after the inspector had opened the same file as a workbook.
        frame = _read_excel_any(
            path,
            0 if sheet_name in (None, "") else sheet_name,
            header_row=resolved_header,
        )
        frame = frame if nrows is None else frame.head(nrows)
        return _recover_r_row_name_column(frame)
    # Delimiter and decoding must match what inspection saw (utf-8-sig with
    # replacement). Blank lines and GEO metadata are excluded before the
    # 1-based header is applied, just as upstream inspection excluded them.
    separator = sniff_delimited_separator(path, matrix_table=True)
    frame = pd.read_csv(
        path,
        sep=separator,
        header=resolved_header - 1,
        comment="!",
        skip_blank_lines=True,
        skipinitialspace=separator == " ",
        low_memory=False,
        encoding="utf-8-sig",
        encoding_errors="replace",
        nrows=nrows,
    )
    return _recover_r_row_name_column(frame)


def _clean_gene_ids(values: pd.Series) -> pd.Series:
    # The one gene-identity rule DEGORA has (harmonize.canonical_gene_symbol):
    # a private copy here stripped `.\d+` from every label and merged NKX2.5
    # with NKX2.1 in matrices exactly as it once did in result tables.
    return _clean_gene_symbol(values)


def _apply_gene_value_transform(
    values: pd.Series,
    transform: str | None,
) -> tuple[pd.Series, dict[str, Any]]:
    name = str(transform or "").strip()
    if not name:
        return values, {}
    if name not in {ENSEMBL_SYMBOL_SUFFIX_TRANSFORM, ENSEMBL_PREFIX_TRANSFORM}:
        raise ValueError(f"unsupported gene_value_transform: {name}")
    labels = values.astype("string").str.strip()
    nonempty = labels.notna() & labels.ne("")
    extracted = labels.str.extract(ENSEMBL_SYMBOL_COMPOSITE_VALUE_RE, expand=True)
    matched = nonempty & extracted["ensembl"].notna() & extracted["symbol"].notna()
    denominator = int(nonempty.sum())
    fraction = float(matched.sum() / denominator) if denominator else 0.0
    if denominator < 2 or fraction < 0.9:
        raise ValueError(
            f"gene_value_transform={name} requires at least 90% ENS...G..._SYMBOL values"
        )
    if name == ENSEMBL_PREFIX_TRANSFORM:
        transformed = extracted["ensembl"].where(matched, pd.NA)
        output_space = "Ensembl gene ID"
        usable = matched
        ensembl_only_dropped = 0
    else:
        cleaned_suffixes = _clean_gene_ids(extracted["symbol"])
        symbol_suffix = (
            matched
            & cleaned_suffixes.notna()
            & ~extracted["symbol"].str.fullmatch(ENSEMBL_GENE_ID_RE, na=False)
        )
        symbol_fraction = float(symbol_suffix.sum() / denominator) if denominator else 0.0
        if symbol_fraction < 0.9:
            raise ValueError(
                "gene_value_transform=ensembl_gene_symbol_suffix_v1 requires at least 90% usable "
                "non-Ensembl gene-symbol suffixes; use the validated prefix transform instead"
            )
        transformed = extracted["symbol"].where(symbol_suffix, pd.NA)
        output_space = "gene symbol"
        usable = symbol_suffix
        ensembl_only_dropped = int((matched & ~symbol_suffix).sum())
    return transformed, {
        "gene_value_transform": name,
        "gene_value_transform_output_space": output_space,
        "gene_value_transform_match_fraction": round(fraction, 6),
        "gene_value_transform_input_nonempty": denominator,
        "gene_value_transform_matched": int(matched.sum()),
        "gene_value_transform_usable_identifiers": int(usable.sum()),
        "gene_value_transform_symbol_suffixes": (
            int(usable.sum()) if name == ENSEMBL_SYMBOL_SUFFIX_TRANSFORM else 0
        ),
        "gene_value_transform_ensembl_only_suffixes_dropped": ensembl_only_dropped,
    }


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
    header_row: int | None = None,
    gene_value_transform: str | None = None,
    command: str = "degora discovery analyze",
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Derive a full DEG-like table using the documented Welch fallback."""

    if role not in SUPPORTED_MATRIX_ROLES:
        raise ValueError("role must be count_matrix, estimated_count_matrix or normalized_expression_matrix")
    source = Path(matrix_path).resolve()
    output = Path(output_path).resolve()
    if not source.exists():
        raise FileNotFoundError(f"matrix does not exist: {source}")
    control, treatment = _validated_samples(control_samples, treatment_samples)
    if not str(gene_column).strip():
        raise ValueError("gene_column is required")
    if str(gene_column).strip().upper() == "ID_REF":
        raise ValueError("ID_REF is a probe identifier; map probes to gene symbols before running the fallback")

    matrix = read_matrix_frame(source, sheet_name=sheet_name, header_row=header_row)
    required = [gene_column, *control, *treatment]
    missing = [column for column in required if column not in matrix.columns]
    if missing:
        raise ValueError("matrix is missing required columns: " + ", ".join(missing))
    transformed_genes, gene_transform_summary = _apply_gene_value_transform(
        matrix[gene_column],
        gene_value_transform,
    )
    genes = _clean_gene_ids(transformed_genes)
    values = matrix[control + treatment].apply(pd.to_numeric, errors="coerce")
    finite = np.isfinite(values.to_numpy(dtype=float)).all(axis=1)
    valid = genes.notna().to_numpy() & finite
    values = values.loc[valid].copy()
    genes = genes.loc[valid]
    if values.empty:
        raise ValueError("matrix has no complete finite rows after gene and sample validation")

    values.insert(0, "__gene__", genes.astype(str).to_numpy())
    if role in COUNT_MATRIX_ROLES:
        sample_values = values[control + treatment]
        if (sample_values < 0).any().any():
            raise ValueError("count matrix contains negative values")
        if role == "count_matrix":
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
        normalization = "logCPM_from_raw_counts" if role == "count_matrix" else "logCPM_from_estimated_counts"
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
        "sheet_name": "" if sheet_name in (None, "") else sheet_name,
        "header_row": 1 if header_row in (None, "") else int(header_row),
        **variance_summary,
        **filter_summary,
        **(metadata or {}),
        # A caller may add provenance, but may not claim a transform that was
        # not actually validated and applied above.
        **gene_transform_summary,
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
        "sheet_name": "" if sheet_name in (None, "") else sheet_name,
        "header_row": 1 if header_row in (None, "") else int(header_row),
        **gene_transform_summary,
        **variance_summary,
        **filter_summary,
    }


__all__ = [
    "SUPPORTED_MATRIX_ROLES",
    "derive_welch_deg",
    "read_matrix_frame",
    "sniff_delimited_separator",
    "validate_sample_groups",
    "welch_with_variance_floor",
]
