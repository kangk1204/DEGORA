"""Minimal DEG-table harmonization for the iteration-1 vertical slice."""

from __future__ import annotations

import gzip
import io
import re
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

from .formula_safety import restore_formula_text_if_marked

# Z-representable p-value floor: norm.isf(P_MIN/2) must stay finite. The smallest
# positive float (np.nextafter(0, 1) = 5e-324) overflows isf to +inf, which the
# non-finite guard then turns into NaN -- silently dropping the MOST significant gene
# (a reported p-value of 0). 1e-300 floors p=0 to a large finite signed-z instead.
P_MIN = 1e-300
# A gzipped workbook has to be expanded in memory before pandas can read it, so
# it needs the bound a streamed CSV does not.
MAX_DECOMPRESSED_WORKBOOK_BYTES = 256 * 1024 * 1024
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
    "sign_convention",
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


def _excel_date_gene_number_candidates(year: int, day: int) -> tuple[int, ...]:
    # Excel frequently stores "Sep-2" as 2002-09-01, but "2-Sep" as
    # <workbook-year>-09-02. Preserve the intended gene-family number in both cases.
    # A candidate is only ever used when it names a real gene in the family, so the
    # year suffix is offered first and silently ignored when it names nothing.
    candidates: list[int] = []
    year_suffix = year % 100
    if day == 1 and 1 <= year_suffix <= 31:
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
                _excel_date_gene_number_candidates(value.year, value.day),
            )
            if symbol:
                return symbol
        return value
    if isinstance(value, date):
        prefix = EXCEL_DATE_GENE_PREFIXES.get(value.month)
        if prefix:
            symbol = _repair_excel_date_gene_symbol_from_candidates(
                prefix,
                _excel_date_gene_number_candidates(value.year, value.day),
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
                _excel_date_gene_number_candidates(year, day),
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


MAPPING_ROLE_LABELS = {
    "gene_column": "gene symbols or gene IDs",
    "lfc_column": "numeric log2 fold change",
    "p_column": "numeric p-value in [0, 1]",
    "padj_column": "adjusted p-value/FDR in [0, 1]",
}


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


def _read_excel_any(path: Path, sheet_name: str | int | None) -> pd.DataFrame:
    """Read a workbook, decompressing it first when it arrived gzipped.

    Repositories serve supplementary workbooks as .xlsx.gz and .xls.gz, and pandas
    reads neither. Such a file used to fall through to the CSV reader and fail on
    the workbook's first byte, which said nothing about the real problem.
    """

    if not path.name.lower().endswith(".gz"):
        return pd.read_excel(path, sheet_name=sheet_name)
    with gzip.open(path, "rb") as handle:
        payload = handle.read(MAX_DECOMPRESSED_WORKBOOK_BYTES + 1)
    if len(payload) > MAX_DECOMPRESSED_WORKBOOK_BYTES:
        raise ValueError(
            f"{path.name} expands beyond the {MAX_DECOMPRESSED_WORKBOOK_BYTES // (1024 * 1024)} MB "
            "workbook safety cap; decompress it and check what it contains before using it"
        )
    return pd.read_excel(io.BytesIO(payload), sheet_name=sheet_name)


def read_deg_table(path: str | Path, mapping: TableMapping) -> pd.DataFrame:
    path = Path(path)
    suffixes = "".join(path.suffixes).lower()

    if suffixes.endswith((".xlsx", ".xls", ".xlsx.gz", ".xls.gz")):
        sheet_name: str | int | None = 0 if mapping.sheet_name in (None, "") else mapping.sheet_name
        return _read_excel_any(path, sheet_name)

    raw_sep = mapping.sep
    auto_sep = raw_sep in (None, "")
    if auto_sep:
        sep = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    else:
        sep = _normalize_separator(raw_sep)
    # A multi-character separator is a regex to pandas, and the C parser cannot take
    # one. Choosing the engine here keeps a plain ParserWarning off the user's screen.
    engine = "python" if len(sep) > 1 else None
    frame = pd.read_csv(path, sep=sep, engine=engine)
    frame = restore_formula_text_if_marked(frame, path)
    frame = _restore_unnamed_row_labels(frame)
    if frame.shape[1] == 1:
        header = str(frame.columns[0])
        # The recovery hint used to run only when the delimiter was auto-detected,
        # which is the case least likely to be wrong. A filled-in 'sep' that does not
        # match the file produced the same single mangled column with no explanation.
        for candidate, label in ((",", "comma"), ("\t", "tab"), (";", "semicolon"), ("|", "pipe")):
            if candidate == sep or candidate not in header:
                continue
            # A real DEG table has at least gene/lfc/p columns, so a header that splits
            # into >=3 fields on another delimiter is the wrong delimiter, not a
            # one-column file whose single header merely contains a comma.
            if len([field for field in header.split(candidate) if field.strip()]) < 3:
                continue
            used = _SEPARATOR_LABELS.get(sep, repr(sep))
            raise ValueError(
                f"{path.name} parsed into a single column with the {used} delimiter, but the header "
                f"looks {label}-delimited ({header[:120]!r}). Set 'sep' in the catalog to the correct "
                f"delimiter, or leave it blank to auto-detect. Accepted values include "
                f"{', '.join(sorted(_SEPARATOR_ALIASES))}."
            )
    return frame


ROW_LABEL_COLUMN = "row_name"


_PANDAS_PLACEHOLDER = re.compile(r"^Unnamed: \d+$")


def _restore_unnamed_row_labels(frame: pd.DataFrame) -> pd.DataFrame:
    """Expose R ``write.csv`` row labels as one predictable column name.

    ``write.csv(results, file)`` has two shapes in the wild. It can write a header
    with one fewer field than the data rows, which pandas resolves by consuming
    the gene identifiers as an unnamed index - somewhere no catalog mapping could
    reference. It can also write an empty leading header field, which pandas names
    ``Unnamed: 0``. Both become ``row_name`` so a catalog author learns one name.
    """

    def free_name(columns: Any) -> str:
        name = ROW_LABEL_COLUMN
        suffix = 2
        while name in columns:
            name = f"{ROW_LABEL_COLUMN}_{suffix}"
            suffix += 1
        return name

    if (
        not isinstance(frame.index, pd.MultiIndex)
        and frame.index.name is None
        and not isinstance(frame.index, pd.RangeIndex)
    ):
        restored = frame.reset_index()
        return restored.rename(columns={restored.columns[0]: free_name(frame.columns)})

    if len(frame.columns) and _PANDAS_PLACEHOLDER.match(str(frame.columns[0])):
        rest = frame.columns[1:]
        return frame.rename(columns={frame.columns[0]: free_name(rest)})
    return frame


def resolve_column_name(frame: pd.DataFrame, requested: str) -> str:
    """Accept the pandas placeholder spelling for restored row labels.

    Catalogs written before the rename refer to the first column as
    ``Unnamed: 0``; keep those working.
    """

    if requested in frame.columns:
        return requested
    if _PANDAS_PLACEHOLDER.match(str(requested)) and ROW_LABEL_COLUMN in frame.columns:
        return ROW_LABEL_COLUMN
    return requested


def validate_table_mapping_roles(frame: pd.DataFrame, mapping: TableMapping, *, study_id: Any = None) -> None:
    """Reject one source-table column being reused for incompatible statistical roles."""

    requested_roles = {
        "gene_column": mapping.gene_column,
        "lfc_column": mapping.lfc_column,
        "p_column": mapping.p_column,
    }
    if mapping.padj_column:
        requested_roles["padj_column"] = mapping.padj_column

    resolved_roles = {
        role: resolve_column_name(frame, str(column))
        for role, column in requested_roles.items()
        if column is not None and str(column).strip()
    }
    conflicts: list[str] = []
    role_items = list(resolved_roles.items())
    for left_index, (left_role, left_column) in enumerate(role_items):
        for right_role, right_column in role_items[left_index + 1 :]:
            if left_column != right_column:
                continue
            if {left_role, right_role} == {"p_column", "padj_column"}:
                continue
            conflicts.append(
                f"{left_role} ({MAPPING_ROLE_LABELS[left_role]}) and "
                f"{right_role} ({MAPPING_ROLE_LABELS[right_role]}) both map to source column {left_column!r}"
            )
    if not conflicts:
        return

    prefix = f"{study_id}: " if study_id not in (None, "") else ""
    raise ValueError(
        prefix
        + "source-table column mappings reuse one column for incompatible roles: "
        + "; ".join(conflicts)
        + ". gene_column, lfc_column and p_column must be distinct; padj_column may equal p_column only when the same unit-interval column is intentionally used as both."
    )


# Catalog authors reach for the word before the escape, so accept both. Without
# this, sep="tab" reached pandas as a three-character regex and collapsed every
# row into one column named "ORF\tGENENAME\t...".
_SEPARATOR_ALIASES = {
    "tab": "\t",
    "tabs": "\t",
    "\\t": "\t",
    "t": "\t",
    "tsv": "\t",
    "comma": ",",
    "csv": ",",
    "semicolon": ";",
    "semi": ";",
    "pipe": "|",
    "bar": "|",
    "space": " ",
    "whitespace": r"\s+",
    "ws": r"\s+",
}
_SEPARATOR_LABELS = {"\t": "tab", ",": "comma", ";": "semicolon", "|": "pipe", " ": "space", r"\s+": "whitespace"}


def _normalize_separator(value: str) -> str:
    """Map a catalog 'sep' entry onto the delimiter pandas should use."""

    text = str(value)
    stripped = text.strip()
    alias = _SEPARATOR_ALIASES.get(stripped.lower())
    if alias is not None:
        return alias
    if stripped == "" and text != "":
        return text  # a deliberate literal space
    return stripped or text


def _row_label_hint(frame: pd.DataFrame) -> str:
    if not any(str(name).startswith(ROW_LABEL_COLUMN) for name in frame.columns):
        return ""
    return (
        f" This file has unnamed row labels (R's write.csv default); they are available"
        f" as {ROW_LABEL_COLUMN!r}."
    )


def _series_as_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    column = resolve_column_name(frame, column)
    if column not in frame.columns:
        raise KeyError(
            f"Required column {column!r} not found. Available columns: {list(frame.columns)!r}."
            + _row_label_hint(frame)
        )
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


def canonical_gene_symbol(value: Any) -> str:
    """Return the symbol DEGORA stores for one user-supplied gene label.

    This is the single definition of "the same gene" across DEGORA. Source DEG
    tables, the optional GoldPanel and browser/API lookups all pass through it,
    so a panel written as ``SEPT9`` and a table written as ``9-Sep`` resolve to
    the one symbol that is actually scored (``SEPTIN9``). Returns "" when the
    value carries no usable identifier.
    """

    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        # Arrays and other non-scalars: fall through to the string form below.
        pass
    repaired = _repair_excel_date_gene_symbol(value)
    text = re.sub(r"\.\d+$", "", str(repaired).strip()).upper()
    return "" if text in {"", "NAN", "NONE", "NA", "<NA>"} else text


def original_gene_label(value: Any) -> str:
    """Return the gene label exactly as the source table carried it."""

    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _collapse_duplicate_gene_symbols(out: pd.DataFrame, study_meta: dict[str, Any]) -> pd.DataFrame:
    """Collapse repeated gene symbols before rank calculation.

    Microarray full tables can be probe-level even when the user maps the gene
    symbol column directly. Keeping repeated symbols would rank probes instead
    of genes and would let downstream source-unit aggregation select among
    probes. Collapse to one row per gene before any within-study rank is made.
    """

    _probe_collapse_raw = study_meta.get("probe_collapse", "")
    if isinstance(_probe_collapse_raw, str):
        requested_probe_collapse = _probe_collapse_raw.strip()
    elif _probe_collapse_raw is None or pd.isna(_probe_collapse_raw):
        requested_probe_collapse = ""
    else:
        requested_probe_collapse = str(_probe_collapse_raw).strip()
    # An empty/blank cell is read by pandas as NaN/None; normalize all of these (and the
    # stringified "nan"/"none") to "" so an undeclared probe_collapse is never mistaken for a
    # declared value that disagrees with the applied collapse rule.
    if requested_probe_collapse.lower() in {"", "nan", "none"}:
        requested_probe_collapse = ""

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
    if "input_gene_label" in frame.columns:
        # Two rows can reach one symbol from different labels (a table carrying
        # both "SEPT9" and "9-Sep"). Keeping only the surviving row's label would
        # hide the other one, so record every distinct label the gene came from.
        merged_labels = (
            frame.loc[frame["input_gene_label"].astype("string").fillna("").ne(""), ["gene_symbol", "input_gene_label"]]
            .astype({"input_gene_label": "string"})
            .drop_duplicates()
            .groupby("gene_symbol", dropna=False)["input_gene_label"]
            .agg(lambda labels: ";".join(sorted(dict.fromkeys(str(label) for label in labels))))
        )
        collapsed["input_gene_label"] = (
            collapsed["gene_symbol"].map(merged_labels).fillna(collapsed["input_gene_label"])
        )
    # Record what was actually applied (best-probe) and what the config declared, so the
    # two can never silently disagree downstream.
    collapsed["gene_symbol_collapse_rule"] = GENE_SYMBOL_COLLAPSE_RULE
    collapsed["requested_probe_collapse"] = requested_probe_collapse

    study_id = str(study_meta.get("study_id", ""))
    study_assay = str(study_meta.get("assay_type", "") or "").strip().lower()
    requested_norm = _normalize_collapse_label(requested_probe_collapse)
    sign_conflict_genes = _count_sign_conflicting_genes(out, counts)
    if requested_probe_collapse == "":
        # probe_collapse genuinely matters for microarray sources, so nudge the user to declare
        # it. For gene-level (RNA-seq) sources, duplicate gene symbols are an ordinary occurrence
        # and the rule actually applied is recorded in gene_symbol_collapse_rule, so emitting a
        # per-source warning for every such table is noise rather than an audit signal.
        if "microarray" in study_assay:
            warning = (
                f"{study_id}: duplicate gene symbols were collapsed by "
                f"{GENE_SYMBOL_COLLAPSE_RULE}; set probe_collapse in the config if this is expected."
            )
        elif _looks_stacked(counts):
            # A gene appearing twice in an RNA-seq table is ordinary. A table that
            # carries a whole block of rows per gene is usually several contrasts or
            # models stacked together, and keeping the most significant row of each
            # block selects on the outcome.
            dropped = int(len(out) - out["gene_symbol"].nunique())
            warning = (
                f"{study_id}: {dropped:,} of {len(out):,} rows shared a gene identifier "
                f"(up to {int(counts.max())} rows per gene) and were collapsed by "
                f"{GENE_SYMBOL_COLLAPSE_RULE}, which keeps the most significant row of each group "
                "and therefore selects on the outcome. If this table stacks several contrasts, "
                "models or cell types, filter it to one contrast before running."
            )
        else:
            warning = ""
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
    if sign_conflict_genes:
        # Directional consistency is the whole point of the score, so a source
        # whose own duplicate rows disagree about the sign must say so even when
        # the collapse rule itself is unremarkable.
        conflict_note = (
            f"{study_id}: {sign_conflict_genes:,} gene(s) had duplicate rows with opposite log2 "
            f"fold-change signs; {GENE_SYMBOL_COLLAPSE_RULE} kept the most significant row, so the "
            "direction reported for those genes comes from one row of a disagreeing pair. Check "
            "whether the table stacks more than one contrast."
        )
        warning = f"{warning} {conflict_note}".strip() if warning else conflict_note
    collapsed["gene_symbol_collapse_warning"] = warning
    return collapsed


def _count_sign_conflicting_genes(out: pd.DataFrame, counts: pd.Series) -> int:
    """Count genes whose duplicate rows disagree about the direction of change."""

    duplicated = counts.gt(1)
    if not bool(duplicated.any()):
        return 0
    signs = np.sign(pd.to_numeric(out["lfc"], errors="coerce")).loc[duplicated]
    grouped = signs.groupby(out.loc[duplicated, "gene_symbol"], dropna=False)
    highest = grouped.max()
    lowest = grouped.min()
    return int(((highest > 0) & (lowest < 0)).sum())


# A gene appearing twice is ordinary; a block of rows per gene is a stacked table.
STACKED_TABLE_MIN_ROWS_PER_GENE = 3
STACKED_TABLE_MIN_DUPLICATE_SHARE = 0.25


def _looks_stacked(counts: pd.Series) -> bool:
    """Report whether duplication looks like stacked contrasts rather than noise."""

    total = int(len(counts))
    if total == 0:
        return False
    duplicated_rows = int((counts > 1).sum())
    return (
        int(counts.max()) >= STACKED_TABLE_MIN_ROWS_PER_GENE
        and duplicated_rows / total >= STACKED_TABLE_MIN_DUPLICATE_SHARE
    )


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


def _is_binary_probability_indicator(values: pd.Series, *, column: str) -> bool:
    """Whether values/header describe a threshold flag rather than probabilities."""

    observed = pd.to_numeric(values, errors="coerce").dropna()
    if observed.empty or not bool(observed.isin((0.0, 1.0)).all()):
        return False
    distinct = {float(value) for value in observed.unique()}
    flag_like_name = bool(
        re.search(
            r"(^|[^a-z0-9])(sig(nificant)?|flag|indicator|pass|threshold(ed)?|is[_. -]?(de|deg))"
            r"($|[^a-z0-9])|^(de|deg)$",
            str(column),
            flags=re.IGNORECASE,
        )
    )
    # A genuine adjusted-p column can legitimately be constant 1 after multiple
    # testing. Reject a constant 0/1 vector only when its header also declares it
    # as a flag; seeing both 0 and 1 is itself enough evidence of thresholding.
    return distinct == {0.0, 1.0} or flag_like_name


def _reject_binary_probability_indicator(
    values: pd.Series,
    *,
    column: str,
    column_kind: str,
    study_meta: dict[str, Any],
) -> None:
    """Reject 0/1 flags that cannot carry gene-level significance evidence."""

    if not _is_binary_probability_indicator(values, column=column):
        return
    study_id = str(study_meta.get("study_id", "unknown_study"))
    raise ValueError(
        f"{study_id}: {column_kind} column {column!r} contains only 0 and/or 1. "
        "This looks like a binary significance flag or thresholded/rounded values, not usable "
        "gene-level probabilities. Map the unrounded p-value or adjusted p-value column instead."
    )


# A source that loses this share of its rows is reporting a mapping or export
# problem, not ordinary missingness, and the run should say so.
ROW_LOSS_WARNING_SHARE = 0.10


def _non_numeric_examples(raw: pd.Series, numeric: pd.Series, limit: int = 5) -> list[str]:
    """Return the original text of values that would not parse as numbers."""

    unparsed = numeric.isna() & raw.notna()
    if not bool(unparsed.any()):
        return []
    seen: list[str] = []
    for value in raw.loc[unparsed]:
        text = str(value).strip()
        if text and text not in seen:
            seen.append(text)
        if len(seen) >= limit:
            break
    return seen


def _unusable_row_warning(
    study_id: str,
    n_input_rows: int,
    n_dropped: int,
    reasons: dict[str, int],
    examples: dict[str, list[str]],
) -> str:
    """Describe rows dropped for missing gene/effect/p-value, by reason.

    Duplicate collapse already warns. Losing rows outright did not, so a table
    whose effect column exported as text, or whose gene column is half empty,
    lost half its rows between the file and the ranking with nothing said. The
    two counts needed to notice it were in different fields of the metrics file.
    """

    # The per-reason counts overlap: one row can be missing its effect and its
    # p-value, and summing them reported more dropped rows than the table had.
    # The denominator is the count of distinct rows the validity mask removed;
    # the reasons stay as a non-exclusive breakdown of why.
    dropped = int(n_dropped)
    if not dropped or not n_input_rows:
        return ""
    share = dropped / n_input_rows
    if share < ROW_LOSS_WARNING_SHARE:
        return ""
    parts = []
    for column, count in reasons.items():
        detail = f"{count:,} missing a {column}"
        sample = examples.get(column) or []
        if sample:
            detail += f" (e.g. {', '.join(repr(text) for text in sample[:3])})"
        parts.append(detail)
    return (
        f"{study_id}: {dropped:,} of {n_input_rows:,} rows ({share:.1%}) were dropped before ranking "
        f"- {'; '.join(parts)} (a row can be missing more than one). A gene symbol, a numeric log2 fold "
        f"change and a numeric p-value are all "
        "required. Check that the mapped columns are the intended ones and that the effect and p-value "
        "columns hold numbers rather than text such as 'NA', 'UP' or a spreadsheet error value."
    )


def harmonize_frame(frame: pd.DataFrame, mapping: TableMapping, study_meta: dict[str, Any]) -> pd.DataFrame:
    """Return canonical per-gene DEG rows for one study/contrast."""

    validate_table_mapping_roles(frame, mapping, study_id=study_meta.get("study_id"))
    scope = assess_table_scope(frame, mapping, study_meta.get("table_scope", AUTO_TABLE_SCOPE))
    # Catalogs written before restored row labels got one name refer to that column
    # by whatever pandas happened to call it.
    gene_column = resolve_column_name(frame, mapping.gene_column)
    if gene_column not in frame.columns:
        raise KeyError(
            f"Required column {mapping.gene_column!r} not found. Available columns: {list(frame.columns)!r}."
            + _row_label_hint(frame)
        )
    genes = _clean_gene_symbol(frame[gene_column])
    input_gene_labels = frame[gene_column].map(original_gene_label)
    lfc = _series_as_numeric(frame, mapping.lfc_column)
    pvalue = _series_as_numeric(frame, mapping.p_column)
    _reject_invalid_unit_interval(
        pvalue,
        column=mapping.p_column,
        column_kind="p-value",
        study_meta=study_meta,
    )
    _reject_binary_probability_indicator(
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
        _reject_binary_probability_indicator(
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
            # The label the source table actually carried. gene_symbol is the
            # canonical symbol DEGORA scores; without this column a repaired
            # "9-Sep" -> SEPTIN9 rename left no trace in any user-facing output.
            "input_gene_label": input_gene_labels,
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

    n_input_rows = int(len(out))
    unusable_reasons = {
        "gene identifier": int(out["gene_symbol"].isna().sum()),
        "log2 fold change": int(out["lfc"].isna().sum()),
        "p-value": int(out["pvalue"].isna().sum()),
    }
    unusable_reasons = {column: count for column, count in unusable_reasons.items() if count}
    unusable_examples = {
        "log2 fold change": _non_numeric_examples(frame[resolve_column_name(frame, mapping.lfc_column)], lfc),
        "p-value": _non_numeric_examples(frame[resolve_column_name(frame, mapping.p_column)], pvalue),
    }
    valid = out["gene_symbol"].notna() & out["lfc"].notna() & out["pvalue"].notna()
    # Settle this here. Taking it from the final row count instead counted the
    # rows that duplicate collapse merged -- ordinary probe-level rows that were
    # used, not lost -- as though their source could not supply them.
    n_rows_dropped_unusable = int(n_input_rows - int(valid.sum()))
    out = out.loc[valid].copy()
    out["pvalue_was_clipped"] = out["pvalue"] < P_MIN
    out["pvalue"] = out["pvalue"].clip(lower=P_MIN, upper=1.0)
    out["padj"] = out["padj"].clip(lower=0.0, upper=1.0)
    out = _collapse_duplicate_gene_symbols(out, study_meta)

    lfc_sign = np.sign(out["lfc"].to_numpy(dtype=float))
    # A log2FC of exactly 0 or a two-sided p-value of exactly 1 carries no usable
    # directional evidence. Either case would otherwise produce signed_z == 0,
    # survive the finite-value guard, add a full weight to the denominator, and
    # inflate source support. Route both through NaN so aggregation/scoring drop
    # the neutral row while the harmonized audit table still records it.
    neutral_evidence = (lfc_sign == 0.0) | out["pvalue"].ge(1.0).to_numpy(dtype=bool)
    out["signed_z"] = np.where(
        neutral_evidence,
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
    out["n_input_rows"] = n_input_rows
    out["n_rows_dropped_unusable"] = n_rows_dropped_unusable
    out["n_rows_merged_by_gene_collapse"] = int(n_input_rows - n_rows_dropped_unusable - len(out))
    out["unusable_row_warning"] = _unusable_row_warning(
        str(study_meta.get("study_id", "unknown_study")),
        n_input_rows,
        n_rows_dropped_unusable,
        unusable_reasons,
        unusable_examples,
    )
    return out.sort_values(["study_id", "within_study_rank", "gene_symbol"]).reset_index(drop=True)


def harmonize_path(path: str | Path, mapping: TableMapping, study_meta: dict[str, Any]) -> pd.DataFrame:
    frame = read_deg_table(path, mapping)
    return harmonize_frame(frame, mapping, study_meta)
