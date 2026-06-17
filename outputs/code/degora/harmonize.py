"""Minimal DEG-table harmonization for the iteration-1 vertical slice."""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import norm

# Z-representable p-value floor: norm.isf(P_MIN/2) must stay finite. The smallest
# positive float (np.nextafter(0, 1) = 5e-324) overflows isf to +inf, which the
# non-finite guard then turns into NaN -- silently dropping the MOST significant gene
# (a reported p-value of 0). 1e-300 floors p=0 to a large finite signed-z instead.
P_MIN = 1e-300
AUTO_TABLE_SCOPE = "auto"
FULL_RESULTS_SCOPE = "full_results"
DEG_ONLY_SCOPE = "deg_only"
AMBIGUOUS_SCOPE = "ambiguous"
TABLE_SCOPE_ALIASES = {
    "": AUTO_TABLE_SCOPE,
    "auto": AUTO_TABLE_SCOPE,
    "infer": AUTO_TABLE_SCOPE,
    "guess": AUTO_TABLE_SCOPE,
    "full": FULL_RESULTS_SCOPE,
    "all": FULL_RESULTS_SCOPE,
    "all_genes": FULL_RESULTS_SCOPE,
    "all_results": FULL_RESULTS_SCOPE,
    "full_results": FULL_RESULTS_SCOPE,
    "full_ranked": FULL_RESULTS_SCOPE,
    "tested_genes": FULL_RESULTS_SCOPE,
    "de_results": FULL_RESULTS_SCOPE,
    "deg": DEG_ONLY_SCOPE,
    "degs": DEG_ONLY_SCOPE,
    "deg_only": DEG_ONLY_SCOPE,
    "significant": DEG_ONLY_SCOPE,
    "significant_only": DEG_ONLY_SCOPE,
    "reported_deg_only": DEG_ONLY_SCOPE,
    "list_only": DEG_ONLY_SCOPE,
    "hit_list": DEG_ONLY_SCOPE,
    "ambiguous": AMBIGUOUS_SCOPE,
    "unknown": AMBIGUOUS_SCOPE,
}

SOURCE_METADATA_COLUMNS = [
    "assay_type",
    "source_input_type",
    "platform",
    "normalization",
    "probe_id_column",
    "probe_collapse",
    "time_course_mode",
    "temporal_mode",
]
GENE_SYMBOL_COLLAPSE_RULE = "min_pvalue_max_abs_lfc"
# Declared probe_collapse values that are consistent with the post-test best-probe
# selection that harmonize_frame applies when it still sees duplicate gene symbols.
# Anything else (e.g. a pre-test expression-level rule such as median_expression)
# must be applied upstream; if it is not, harmonize_frame records the mismatch
# loudly instead of silently relabeling best-probe selection as that rule.
BEST_PROBE_COLLAPSE_ALIASES = frozenset(
    {"min_pvalue_max_abs_lfc", "best_probe", "min_p", "min_pvalue", "min_p_max_lfc"}
)
EXCEL_DATE_GENE_PREFIXES = {
    3: "MARCH",
    9: "SEPT",
    12: "DEC",
}
EXCEL_DATE_GENE_MONTHS = {
    "DEC": "DEC",
    "DECEMBER": "DEC",
    "MAR": "MARCH",
    "MARCH": "MARCH",
    "SEP": "SEPT",
    "SEPT": "SEPT",
    "SEPTEMBER": "SEPT",
}
EXCEL_DATE_GENE_CURRENT_SYMBOLS = {
    "DEC": {
        1: "BHLHE40",
        2: "BHLHE41",
    },
    "MARCH": {number: f"MARCHF{number}" for number in range(1, 12)},
    "SEPT": {number: f"SEPTIN{number}" for number in (*range(1, 13), 14)},
}


def _normalize_collapse_label(value: Any) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def _excel_date_gene_number_candidates(year: int, month: int, day: int) -> tuple[int, ...]:
    # Excel frequently stores "Sep-2" as 2002-09-01, but "2-Sep" as
    # <workbook-year>-09-02. Preserve the intended gene-family number in both cases.
    candidates: list[int] = []
    year_suffix = year % 100
    if day == 1 and 1 <= year_suffix <= 31 and year not in {2025, 2026}:
        candidates.append(year_suffix)
    candidates.append(day)
    return tuple(dict.fromkeys(candidates))


def _current_symbol_for_excel_date_gene(prefix: str, number: int) -> str | None:
    return EXCEL_DATE_GENE_CURRENT_SYMBOLS.get(prefix, {}).get(number)


def _repair_excel_date_gene_symbol_from_candidates(prefix: str, numbers: Iterable[int]) -> str | None:
    for number in numbers:
        symbol = _current_symbol_for_excel_date_gene(prefix, number)
        if symbol:
            return symbol
    return None


def _repair_excel_date_gene_symbol(value: Any) -> Any:
    """Undo common Excel date coercion and legacy aliases for date-like gene families."""

    if value is None or pd.isna(value):
        return value
    if isinstance(value, pd.Timestamp):
        value = value.to_pydatetime()
    if isinstance(value, datetime):
        prefix = EXCEL_DATE_GENE_PREFIXES.get(value.month)
        if prefix:
            symbol = _repair_excel_date_gene_symbol_from_candidates(
                prefix,
                _excel_date_gene_number_candidates(value.year, value.month, value.day),
            )
            if symbol:
                return symbol
        return value
    if isinstance(value, date):
        prefix = EXCEL_DATE_GENE_PREFIXES.get(value.month)
        if prefix:
            symbol = _repair_excel_date_gene_symbol_from_candidates(
                prefix,
                _excel_date_gene_number_candidates(value.year, value.month, value.day),
            )
            if symbol:
                return symbol
        return value

    text = str(value).strip()
    match = re.fullmatch(r"(DEC|MARCH|SEPT)(\d{1,2})", text, re.IGNORECASE)
    if match:
        symbol = _current_symbol_for_excel_date_gene(match.group(1).upper(), int(match.group(2)))
        if symbol:
            return symbol
    match = re.fullmatch(r"(\d{4})-(\d{1,2})-(\d{1,2})(?:\s+00:00:00)?", text)
    if match:
        year, month, day = map(int, match.groups())
        prefix = EXCEL_DATE_GENE_PREFIXES.get(month)
        if prefix:
            symbol = _repair_excel_date_gene_symbol_from_candidates(
                prefix,
                _excel_date_gene_number_candidates(year, month, day),
            )
            if symbol:
                return symbol
    match = re.fullmatch(r"(\d{1,2})[-/ ]([A-Za-z]{3,9})", text)
    if match:
        number = int(match.group(1))
        prefix = EXCEL_DATE_GENE_MONTHS.get(match.group(2).upper())
        if prefix:
            symbol = _current_symbol_for_excel_date_gene(prefix, number)
            if symbol:
                return symbol
    match = re.fullmatch(r"([A-Za-z]{3,9})[-/ ](\d{1,2})", text)
    if match:
        prefix = EXCEL_DATE_GENE_MONTHS.get(match.group(1).upper())
        if prefix:
            symbol = _current_symbol_for_excel_date_gene(prefix, int(match.group(2)))
            if symbol:
                return symbol
    return value


@dataclass(frozen=True)
class TableMapping:
    gene_column: str
    lfc_column: str
    p_column: str
    padj_column: str | None = None
    sep: str | None = None
    sheet_name: str | int | None = None


def normalize_table_scope(value: Any) -> str:
    """Normalize a user-entered DEG table scope label."""

    if value is None or pd.isna(value):
        return AUTO_TABLE_SCOPE
    text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if text in TABLE_SCOPE_ALIASES:
        return TABLE_SCOPE_ALIASES[text]
    raise ValueError(
        f"Unsupported table_scope={value!r}. Use auto, full_results, deg_only, or ambiguous."
    )


def read_deg_table(path: str | Path, mapping: TableMapping) -> pd.DataFrame:
    path = Path(path)
    suffixes = "".join(path.suffixes).lower()

    if suffixes.endswith((".xlsx", ".xls")):
        sheet_name: str | int | None = 0 if mapping.sheet_name in (None, "") else mapping.sheet_name
        return pd.read_excel(path, sheet_name=sheet_name)

    sep = mapping.sep
    auto_sep = sep in (None, "")
    if auto_sep:
        sep = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    frame = pd.read_csv(path, sep=sep)
    if auto_sep and frame.shape[1] == 1:
        header = str(frame.columns[0])
        other = "," if sep == "\t" else "\t"
        # Require >=3 fields when splitting on the other delimiter: a real DEG table has at
        # least gene/lfc/p columns, so a single column that splits into >=3 is the wrong
        # delimiter -- not a legitimate one-column file whose header merely contains a comma.
        if other in header and len([field for field in header.split(other) if field.strip()]) >= 3:
            used = "tab" if sep == "\t" else "comma"
            looks = "comma" if other == "," else "tab"
            raise ValueError(
                f"{path.name} parsed into a single column with the {used} delimiter, but the header "
                f"looks {looks}-delimited ({header!r}). Set 'sep' in the catalog to the correct "
                "delimiter (use \\t for a tab-separated file)."
            )
    return frame


def _series_as_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        raise KeyError(f"Required column {column!r} not found. Available columns: {list(frame.columns)!r}")
    return pd.to_numeric(frame[column], errors="coerce")


def _numeric_if_present(frame: pd.DataFrame, column: str | None) -> pd.Series | None:
    if not column or column not in frame.columns:
        return None
    return pd.to_numeric(frame[column], errors="coerce")


def assess_table_scope(frame: pd.DataFrame, mapping: TableMapping, declared_scope: Any = AUTO_TABLE_SCOPE) -> dict[str, Any]:
    """Classify whether a source table looks full-result or DEG-only.

    The result is deliberately conservative. A DEG-only table can be used by
    DEGORA, but absent genes are treated as unreported rather than
    non-differential, and rank denominators need an explicit
    ``rank_universe_size`` to avoid optimistic ranks.
    """

    declared = normalize_table_scope(declared_scope)
    value_column = mapping.padj_column if mapping.padj_column and mapping.padj_column in frame.columns else mapping.p_column
    values = _numeric_if_present(frame, value_column)
    numeric = values.dropna().clip(lower=0.0, upper=1.0) if values is not None else pd.Series(dtype=float)
    n_rows = int(len(frame))
    n_numeric = int(len(numeric))
    stats: dict[str, Any] = {
        "declared_scope": declared,
        "value_column": value_column,
        "n_rows": n_rows,
        "n_numeric_values": n_numeric,
        "n_le_0_05": int(numeric.le(0.05).sum()) if n_numeric else 0,
        "n_gt_0_05": int(numeric.gt(0.05).sum()) if n_numeric else 0,
        "fraction_le_0_05": float(numeric.le(0.05).mean()) if n_numeric else None,
        "max_value": float(numeric.max()) if n_numeric else None,
    }

    if declared != AUTO_TABLE_SCOPE:
        return {
            **stats,
            "effective_scope": declared,
            "assessment": "declared",
            "reason": f"table_scope was explicitly set to {declared}",
        }

    if n_numeric == 0:
        return {
            **stats,
            "effective_scope": AMBIGUOUS_SCOPE,
            "assessment": "auto_ambiguous",
            "reason": "no numeric p-value/FDR values were available for scope inference",
        }

    n_gt_005 = int(stats["n_gt_0_05"])
    fraction_sig = float(stats["fraction_le_0_05"] or 0.0)
    max_value = float(stats["max_value"] or 0.0)
    if (n_rows >= 10_000 and n_gt_005 >= 100 and max_value >= 0.5) or (
        n_gt_005 >= max(100, int(0.1 * n_rows)) and max_value >= 0.2
    ):
        return {
            **stats,
            "effective_scope": FULL_RESULTS_SCOPE,
            "assessment": "full_results_likely",
            "reason": "table contains many non-significant rows and high p/FDR values",
        }
    if fraction_sig >= 0.98 and max_value <= 0.05:
        return {
            **stats,
            "effective_scope": DEG_ONLY_SCOPE,
            "assessment": "deg_only_likely",
            "reason": "nearly all rows satisfy p/FDR <= 0.05",
        }
    if n_rows < 5_000 and fraction_sig >= 0.95 and max_value <= 0.1:
        return {
            **stats,
            "effective_scope": DEG_ONLY_SCOPE,
            "assessment": "deg_only_likely",
            "reason": "short table dominated by significant rows",
        }
    return {
        **stats,
        "effective_scope": AMBIGUOUS_SCOPE,
        "assessment": "auto_ambiguous",
        "reason": "scope could not be classified confidently; set table_scope explicitly",
    }


def _rank_universe_size(study_meta: dict[str, Any], observed_rows: int, scope: str) -> tuple[int, float | None, str]:
    raw = study_meta.get("rank_universe_size", "")
    declared = pd.to_numeric(pd.Series([raw]), errors="coerce").iloc[0]
    declared_value = float(declared) if pd.notna(declared) and np.isfinite(float(declared)) and float(declared) > 0 else None
    if declared_value is not None:
        used = max(int(round(declared_value)), observed_rows)
        if int(round(declared_value)) < observed_rows:
            # A declared universe smaller than the reported rows is impossible for a
            # real DEG list. Clamp up to observed rows (math stays correct) but surface
            # the misconfiguration instead of silently advertising the unused declared
            # value in rank_universe_size_declared or in a warning that claims it was used.
            warning = (
                f"declared rank_universe_size={int(round(declared_value))} < {observed_rows} reported "
                "rows; using observed rows -- check the catalog rank_universe_size"
            )
            return used, float(used), warning
        warning = (
            "DEG-only table; normalized ranks use declared rank_universe_size and missing genes are unreported"
            if scope == DEG_ONLY_SCOPE
            else ""
        )
        return used, declared_value, warning
    if scope == DEG_ONLY_SCOPE:
        return (
            observed_rows,
            None,
            "DEG-only table without rank_universe_size; normalized ranks use reported-list length and may be optimistic",
        )
    if scope == AMBIGUOUS_SCOPE:
        return (
            observed_rows,
            None,
            "table_scope ambiguous; normalized ranks use observed rows and the source should be reviewed",
        )
    return observed_rows, None, ""


def _clean_gene_symbol(values: pd.Series) -> pd.Series:
    return (
        values.map(_repair_excel_date_gene_symbol)
        .astype("string")
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
        .str.upper()
        .replace({"": pd.NA, "NAN": pd.NA, "NONE": pd.NA})
    )


def _collapse_duplicate_gene_symbols(out: pd.DataFrame, study_meta: dict[str, Any]) -> pd.DataFrame:
    """Collapse repeated gene symbols before rank calculation.

    Microarray full tables can be probe-level even when the user maps the gene
    symbol column directly. Keeping repeated symbols would rank probes instead
    of genes and would let downstream source-unit aggregation select among
    probes. Collapse to one row per gene before any within-study rank is made.
    """

    requested_probe_collapse = str(study_meta.get("probe_collapse", "") or "").strip()

    if out.empty:
        out["n_source_rows_for_gene"] = pd.Series(dtype=int)
        out["gene_symbol_collapse_rule"] = pd.Series(dtype=str)
        out["requested_probe_collapse"] = pd.Series(dtype=str)
        out["gene_symbol_collapse_warning"] = pd.Series(dtype=str)
        return out

    counts = out.groupby("gene_symbol", dropna=False)["gene_symbol"].transform("size")
    if counts.max() <= 1:
        out["n_source_rows_for_gene"] = 1
        out["gene_symbol_collapse_rule"] = "none"
        out["requested_probe_collapse"] = requested_probe_collapse
        out["gene_symbol_collapse_warning"] = ""
        return out

    frame = out.copy()
    frame["_source_row_order"] = range(len(frame))
    frame["_abs_lfc"] = frame["lfc"].abs()
    frame["_n_source_rows_for_gene"] = counts
    frame = frame.sort_values(
        ["gene_symbol", "pvalue", "_abs_lfc", "padj", "_source_row_order"],
        ascending=[True, True, False, True, True],
    )
    collapsed = frame.drop_duplicates("gene_symbol", keep="first").drop(columns=["_abs_lfc", "_source_row_order"])
    collapsed["n_source_rows_for_gene"] = collapsed["_n_source_rows_for_gene"].astype(int)
    collapsed = collapsed.drop(columns=["_n_source_rows_for_gene"])
    # Record what was actually applied (best-probe) and what the config declared, so the
    # two can never silently disagree downstream.
    collapsed["gene_symbol_collapse_rule"] = GENE_SYMBOL_COLLAPSE_RULE
    collapsed["requested_probe_collapse"] = requested_probe_collapse

    study_id = str(study_meta.get("study_id", ""))
    requested_norm = _normalize_collapse_label(requested_probe_collapse)
    if requested_probe_collapse == "":
        warning = (
            f"{study_id}: duplicate gene symbols were collapsed by "
            f"{GENE_SYMBOL_COLLAPSE_RULE}; set probe_collapse in the config if this is expected."
        )
    elif requested_norm in BEST_PROBE_COLLAPSE_ALIASES:
        warning = ""
    else:
        warning = (
            f"{study_id}: probe_collapse={requested_probe_collapse!r} was declared, but the source table still "
            f"had duplicate gene symbols that were collapsed by post-test {GENE_SYMBOL_COLLAPSE_RULE} "
            "(best-probe selection). Pre-collapse probes to one row per gene upstream "
            "(e.g. derive_microarray_deg.py --collapse-rule median_expression) so the declared rule is the one "
            "actually applied; otherwise the recorded collapse rule will not match probe_collapse."
        )
    collapsed["gene_symbol_collapse_warning"] = warning
    return collapsed


def _reject_invalid_unit_interval(
    values: pd.Series,
    *,
    column: str,
    column_kind: str,
    study_meta: dict[str, Any],
) -> None:
    invalid = values.notna() & ((values < 0.0) | (values > 1.0))
    if not bool(invalid.any()):
        return
    study_id = str(study_meta.get("study_id", "unknown_study"))
    examples = ", ".join(f"{float(value):g}" for value in values.loc[invalid].head(5))
    n_invalid = int(invalid.sum())
    raise ValueError(
        f"{study_id}: {column_kind} column {column!r} contains {n_invalid} value(s) outside [0, 1] "
        f"(examples: {examples}). Map the correct unit-interval column before running DEGORA; "
        "do not map -log10 values, test statistics, percentages, or fold-change columns."
    )


def harmonize_frame(frame: pd.DataFrame, mapping: TableMapping, study_meta: dict[str, Any]) -> pd.DataFrame:
    """Return canonical per-gene DEG rows for one study/contrast."""

    scope = assess_table_scope(frame, mapping, study_meta.get("table_scope", AUTO_TABLE_SCOPE))
    if mapping.gene_column not in frame.columns:
        raise KeyError(
            f"Required column {mapping.gene_column!r} not found. Available columns: {list(frame.columns)!r}"
        )
    genes = _clean_gene_symbol(frame[mapping.gene_column])
    lfc = _series_as_numeric(frame, mapping.lfc_column)
    pvalue = _series_as_numeric(frame, mapping.p_column)
    _reject_invalid_unit_interval(
        pvalue,
        column=mapping.p_column,
        column_kind="p-value",
        study_meta=study_meta,
    )
    if mapping.padj_column:
        padj = _series_as_numeric(frame, mapping.padj_column)
        _reject_invalid_unit_interval(
            padj,
            column=mapping.padj_column,
            column_kind="adjusted p-value/FDR",
            study_meta=study_meta,
        )
    else:
        padj = pd.Series(np.nan, index=frame.index, dtype=float)

    out = pd.DataFrame(
        {
            "study_id": study_meta["study_id"],
            "paper_id": study_meta.get("paper_id", study_meta["study_id"]),
            "source_unit_id": "" if pd.isna(study_meta.get("source_unit_id")) else str(study_meta.get("source_unit_id") or "").strip(),
            "gene_symbol": genes,
            "lfc": lfc,
            "pvalue": pvalue,
            "padj": padj,
            "pipeline": study_meta.get("pipeline", "unknown_pipeline"),
            "species": study_meta.get("species", ""),
            "cell_system": study_meta.get("cell_system", ""),
            "hypoxia_modality": study_meta.get("hypoxia_modality", ""),
            "duration_h": study_meta.get("duration_h", ""),
            "n_ctrl": study_meta.get("n_ctrl", np.nan),
            "n_treat": study_meta.get("n_treat", np.nan),
            "source_path": study_meta.get("source_path", ""),
            "source_url": study_meta.get("source_url", ""),
        }
    )
    for column in SOURCE_METADATA_COLUMNS:
        out[column] = study_meta.get(column, "")

    valid = out["gene_symbol"].notna() & out["lfc"].notna() & out["pvalue"].notna()
    out = out.loc[valid].copy()
    out["pvalue_was_clipped"] = out["pvalue"] <= 0.0
    out["pvalue"] = out["pvalue"].clip(lower=P_MIN, upper=1.0)
    out["padj"] = out["padj"].clip(lower=0.0, upper=1.0)
    out = _collapse_duplicate_gene_symbols(out, study_meta)

    lfc_sign = np.sign(out["lfc"].to_numpy(dtype=float))
    # A log2FC of exactly 0 carries no usable direction; np.sign -> 0 would make
    # signed_z == 0, which is finite and survives the guard below. Such a row would
    # then add 0 to the Stouffer numerator but a full weight to the denominator
    # (diluting the combined z) and inflate n_studies. Route lfc == 0 through NaN so
    # aggregation/scoring drop it like any other non-finite signed_z.
    out["signed_z"] = np.where(
        lfc_sign == 0.0,
        np.nan,
        lfc_sign * norm.isf(out["pvalue"].to_numpy(dtype=float) / 2.0),
    )
    out.loc[~np.isfinite(out["signed_z"]), "signed_z"] = np.nan
    out["abs_signed_z"] = out["signed_z"].abs()
    out["within_study_rank"] = out["abs_signed_z"].rank(method="average", ascending=False)
    rank_universe_used, rank_universe_declared, rank_warning = _rank_universe_size(
        study_meta,
        int(len(out)),
        str(scope["effective_scope"]),
    )
    out["n_genes_in_study"] = int(rank_universe_used)
    out["normalized_rank"] = out["within_study_rank"] / max(rank_universe_used, 1)
    out["table_scope"] = str(scope["effective_scope"])
    out["table_scope_assessment"] = str(scope["assessment"])
    out["table_scope_reason"] = str(scope["reason"])
    out["table_scope_value_column"] = str(scope["value_column"])
    out["n_rows_in_source_table"] = int(scope["n_rows"])
    out["n_reported_rows_after_filter"] = int(len(out))
    # Counts rows at value_column <= 0.05, where value_column is padj when mapped
    # else the raw p-value (see assess_table_scope). The neutral name avoids
    # asserting an adjusted-p/FDR threshold for raw-p-only tables; the exact column
    # used is recorded separately in table_scope_value_column.
    out["n_scope_significant_rows_le_0_05"] = int(scope["n_le_0_05"])
    out["rank_universe_size_declared"] = rank_universe_declared if rank_universe_declared is not None else np.nan
    out["rank_universe_size_used"] = int(rank_universe_used)
    out["rank_universe_warning"] = rank_warning
    return out.sort_values(["study_id", "within_study_rank", "gene_symbol"]).reset_index(drop=True)


def harmonize_path(path: str | Path, mapping: TableMapping, study_meta: dict[str, Any]) -> pd.DataFrame:
    frame = read_deg_table(path, mapping)
    return harmonize_frame(frame, mapping, study_meta)
