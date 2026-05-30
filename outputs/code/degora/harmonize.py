"""Minimal DEG-table harmonization for the iteration-1 vertical slice."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

P_MIN = np.nextafter(0.0, 1.0)
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
    if sep in (None, ""):
        sep = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return pd.read_csv(path, sep=sep)


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
        values.astype("string")
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

    if out.empty:
        out["n_source_rows_for_gene"] = pd.Series(dtype=int)
        out["gene_symbol_collapse_rule"] = pd.Series(dtype=str)
        out["gene_symbol_collapse_warning"] = pd.Series(dtype=str)
        return out

    counts = out.groupby("gene_symbol", dropna=False)["gene_symbol"].transform("size")
    if counts.max() <= 1:
        out["n_source_rows_for_gene"] = 1
        out["gene_symbol_collapse_rule"] = "none"
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
    collapsed["gene_symbol_collapse_rule"] = GENE_SYMBOL_COLLAPSE_RULE
    missing_rule = str(study_meta.get("probe_collapse", "") or "").strip() == ""
    warning = ""
    if missing_rule:
        warning = (
            f"{study_meta.get('study_id', '')}: duplicate gene symbols were collapsed by "
            f"{GENE_SYMBOL_COLLAPSE_RULE}; set probe_collapse in the config if this is expected."
        )
    collapsed["gene_symbol_collapse_warning"] = warning
    return collapsed


def harmonize_frame(frame: pd.DataFrame, mapping: TableMapping, study_meta: dict[str, Any]) -> pd.DataFrame:
    """Return canonical per-gene DEG rows for one study/contrast."""

    scope = assess_table_scope(frame, mapping, study_meta.get("table_scope", AUTO_TABLE_SCOPE))
    genes = _clean_gene_symbol(frame[mapping.gene_column])
    lfc = _series_as_numeric(frame, mapping.lfc_column)
    pvalue = _series_as_numeric(frame, mapping.p_column)
    if mapping.padj_column:
        padj = _series_as_numeric(frame, mapping.padj_column)
    else:
        padj = pd.Series(np.nan, index=frame.index, dtype=float)

    out = pd.DataFrame(
        {
            "study_id": study_meta["study_id"],
            "paper_id": study_meta.get("paper_id", study_meta["study_id"]),
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
    out["pvalue_was_clipped"] = (out["pvalue"] <= 0.0) | (out["pvalue"] > 1.0)
    out["pvalue"] = out["pvalue"].clip(lower=P_MIN, upper=1.0)
    out["padj"] = out["padj"].clip(lower=0.0, upper=1.0)
    out = _collapse_duplicate_gene_symbols(out, study_meta)

    lfc_sign = np.sign(out["lfc"].to_numpy(dtype=float))
    out["signed_z"] = lfc_sign * norm.isf(out["pvalue"].to_numpy(dtype=float) / 2.0)
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
    out["n_scope_significant_rows_padj05"] = int(scope["n_le_0_05"])
    out["rank_universe_size_declared"] = rank_universe_declared if rank_universe_declared is not None else np.nan
    out["rank_universe_size_used"] = int(rank_universe_used)
    out["rank_universe_warning"] = rank_warning
    out["sign_call"] = np.where(out["padj"].le(0.1), np.sign(out["lfc"]), 0.0)
    return out.sort_values(["study_id", "within_study_rank", "gene_symbol"]).reset_index(drop=True)


def harmonize_path(path: str | Path, mapping: TableMapping, study_meta: dict[str, Any]) -> pd.DataFrame:
    frame = read_deg_table(path, mapping)
    return harmonize_frame(frame, mapping, study_meta)
