"""Consensus aggregation methods for the thin DEGORA slice."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

from .stats import bh_adjust


STOUFFER_COLUMNS = [
    "gene_symbol",
    "n_studies",
    "stouffer_z",
    "stouffer_p",
    "weighted_lfc",
    "sign_concordance",
    "heterogeneity_q",
    "heterogeneity_df",
    "heterogeneity_i2",
    "stouffer_padj",
]
RANK_COLUMNS = ["gene_symbol", "n_studies_rank", "rank_product", "rank_score"]
COLLAPSED_SOURCE_UNIT_COLUMNS = [
    "gene_symbol",
    "study_id",
    "source_unit_id",
    "signed_z",
    "lfc",
    "normalized_rank",
    "weight",
    "n_contrast_rows",
    "n_studies_in_source_unit",
]
MAX_SOURCE_SAMPLE_WEIGHT = 4.0
SOURCE_UNIT_COLLAPSE_RULE = (
    "source-unit mean aggregation: sample-size-weighted mean signed_z, "
    "sample-size-weighted mean log2FC, mean normalized rank, and mean source weight "
    f"with per-contrast sample-size weights capped at {MAX_SOURCE_SAMPLE_WEIGHT:g}; "
    "no max-|z| representative selection. Optional time_course_mode can preselect "
    "mean, early, late, or peak_mean rows within a source unit before this aggregation."
)


def _source_unit_series(frame: pd.DataFrame) -> pd.Series:
    # Independent source unit precedence: an explicit source_unit_id wins, then
    # paper_id, then study_id. Honoring an explicit source_unit_id means a caller
    # who already grouped contrasts into source units (e.g. via the public API)
    # gets exactly the grouping they declared instead of a paper_id/study_id
    # re-derivation.
    study_id = frame["study_id"].astype("string").fillna("").str.strip()
    result = study_id
    if "paper_id" in frame.columns:
        paper_id = frame["paper_id"].astype("string").fillna("").str.strip()
        result = paper_id.mask(paper_id.eq(""), result)
    if "source_unit_id" in frame.columns:
        source_unit_id = frame["source_unit_id"].astype("string").fillna("").str.strip()
        result = source_unit_id.mask(source_unit_id.eq(""), result)
    return result


def _normalize_time_course_mode(value: Any) -> str:
    if value is None or pd.isna(value):
        return "mean"
    label = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if label in {"", "auto", "all", "source_mean", "average"}:
        return "mean"
    if label in {"first", "earliest", "early"}:
        return "early"
    if label in {"last", "latest", "late"}:
        return "late"
    if label in {"peak", "peak_mean", "strongest_window"}:
        return "peak_mean"
    return "mean"


def _duration_numeric(values: pd.Series) -> pd.Series:
    """Parse simple duration labels for predeclared temporal aggregation."""

    text = values.astype("string").fillna("").str.extract(r"([-+]?\d*\.?\d+)", expand=False)
    return pd.to_numeric(text, errors="coerce")


def _apply_time_course_mode(frame: pd.DataFrame) -> pd.DataFrame:
    mode_column = "time_course_mode" if "time_course_mode" in frame.columns else "temporal_mode" if "temporal_mode" in frame.columns else ""
    if not mode_column:
        frame["time_course_mode"] = "mean"
        return frame
    frame = frame.copy()
    frame["time_course_mode"] = frame[mode_column].map(_normalize_time_course_mode)
    if frame["time_course_mode"].eq("mean").all():
        return frame

    group_cols = ["gene_symbol", "source_unit_id"]
    selected = [frame.loc[frame["time_course_mode"].eq("mean")]]
    duration = _duration_numeric(frame["duration_h"]) if "duration_h" in frame.columns else pd.Series(np.nan, index=frame.index)

    for mode, reducer in [("early", "min"), ("late", "max")]:
        subset = frame.loc[frame["time_course_mode"].eq(mode)].copy()
        if subset.empty:
            continue
        subset["_duration"] = duration.loc[subset.index]
        finite = np.isfinite(subset["_duration"])
        if finite.any():
            targets = subset.loc[finite].groupby(group_cols, sort=False)["_duration"].transform(reducer)
            subset.loc[finite, "_target_duration"] = targets
            has_finite = subset.groupby(group_cols, sort=False)["_duration"].transform(lambda values: np.isfinite(values).any())
            subset = subset.loc[~has_finite | (subset["_duration"].eq(subset["_target_duration"]))]
        selected.append(subset.drop(columns=["_duration", "_target_duration"], errors="ignore"))

    peak = frame.loc[frame["time_course_mode"].eq("peak_mean")].copy()
    if not peak.empty:
        peak["_abs_for_peak"] = peak["signed_z"].abs()
        peak = peak.sort_values([*group_cols, "_abs_for_peak", "study_id"], ascending=[True, True, False, True], kind="mergesort")
        group_size = peak.groupby(group_cols, sort=False)["study_id"].transform("size")
        group_size_values = group_size.to_numpy(dtype=int)
        keep_n = np.where(group_size_values <= 2, group_size_values, np.maximum(2, np.ceil(group_size_values * 0.5))).astype(int)
        rank_from_top = peak.groupby(group_cols, sort=False).cumcount()
        selected.append(peak.loc[rank_from_top < keep_n].drop(columns=["_abs_for_peak"]))

    out = pd.concat([part for part in selected if not part.empty], ignore_index=True)
    return out if not out.empty else frame.iloc[0:0].copy()


def source_unit_rows_for_aggregation(harmonized: pd.DataFrame) -> pd.DataFrame:
    """Return the exact harmonized rows eligible for source-unit aggregation."""

    if harmonized.empty:
        frame = harmonized.copy()
        required_columns = [
            "gene_symbol",
            "study_id",
            "source_unit_id",
            "signed_z",
            "lfc",
            "normalized_rank",
            "_weight",
        ]
        for column in required_columns:
            if column not in frame.columns:
                dtype = float if column in {"signed_z", "lfc", "normalized_rank", "_weight"} else "string"
                frame[column] = pd.Series(dtype=dtype)
        return frame

    frame = harmonized.dropna(subset=["signed_z"]).copy()
    frame["source_unit_id"] = _source_unit_series(frame)
    frame["study_id"] = frame["study_id"].astype("string").fillna("").str.strip()
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.upper().str.strip()
    frame["signed_z"] = pd.to_numeric(frame["signed_z"], errors="coerce")
    frame["lfc"] = pd.to_numeric(frame["lfc"], errors="coerce")
    frame["normalized_rank"] = pd.to_numeric(frame["normalized_rank"], errors="coerce")
    n_ctrl = (
        pd.to_numeric(frame["n_ctrl"], errors="coerce")
        if "n_ctrl" in frame.columns
        else pd.Series(np.nan, index=frame.index)
    )
    n_treat = (
        pd.to_numeric(frame["n_treat"], errors="coerce")
        if "n_treat" in frame.columns
        else pd.Series(np.nan, index=frame.index)
    )
    frame["_weight"] = np.where(
        np.isfinite(n_ctrl) & np.isfinite(n_treat) & ((n_ctrl + n_treat) > 0),
        np.sqrt(n_ctrl + n_treat),
        1.0,
    )
    frame["_weight"] = np.minimum(frame["_weight"], MAX_SOURCE_SAMPLE_WEIGHT)
    frame = frame.dropna(subset=["gene_symbol", "study_id", "source_unit_id", "signed_z", "normalized_rank"])
    frame = frame.loc[frame["gene_symbol"].ne("") & frame["source_unit_id"].ne("")].copy()
    if frame.empty:
        return frame
    return _apply_time_course_mode(frame)


def collapse_gene_source_units(harmonized: pd.DataFrame) -> pd.DataFrame:
    """Collapse rows to one aggregate gene row per independent source unit.

    The biological replication unit for DEGORA is the independent source unit
    (`paper_id`/dataset family), not every contrast row. Multiple time points,
    cell lines, or technical table rows from one source can support a gene, but
    they must not increase the cross-study Stouffer or rank-product sample size.

    Within a source unit, DEGORA now aggregates rather than choosing the maximum
    absolute z value. This avoids reintroducing multiplicity through a
    winner-take-all representative contrast.
    """

    if harmonized.empty:
        return pd.DataFrame(columns=COLLAPSED_SOURCE_UNIT_COLUMNS)

    frame = source_unit_rows_for_aggregation(harmonized)
    if frame.empty:
        return pd.DataFrame(columns=COLLAPSED_SOURCE_UNIT_COLUMNS)

    frame = frame.sort_values(
        ["gene_symbol", "source_unit_id", "study_id", "normalized_rank"],
        kind="mergesort",
    ).reset_index(drop=True)
    lfc_is_valid = frame["lfc"].notna()

    frame["_wz"] = frame["_weight"] * frame["signed_z"]
    frame["_wlfc"] = np.where(lfc_is_valid, frame["_weight"] * frame["lfc"], 0.0)
    frame["_w_lfc_denominator"] = np.where(lfc_is_valid, frame["_weight"], 0.0)
    grouped = frame.groupby(["gene_symbol", "source_unit_id"], as_index=False, sort=False).agg(
        study_id=("study_id", "first"),
        sum_wz=("_wz", "sum"),
        sum_w=("_weight", "sum"),
        sum_wlfc=("_wlfc", "sum"),
        sum_w_lfc=("_w_lfc_denominator", "sum"),
        normalized_rank=("normalized_rank", "mean"),
        weight=("_weight", "mean"),
        n_contrast_rows=("study_id", "size"),
        n_studies_in_source_unit=("study_id", "nunique"),
    )
    grouped["signed_z"] = grouped["sum_wz"] / grouped["sum_w"]
    grouped["lfc"] = np.where(grouped["sum_w_lfc"].gt(0), grouped["sum_wlfc"] / grouped["sum_w_lfc"], np.nan)
    return grouped[COLLAPSED_SOURCE_UNIT_COLUMNS].reset_index(drop=True)


def _study_gene_stats(harmonized: pd.DataFrame) -> pd.DataFrame:
    return collapse_gene_source_units(harmonized)


def _eligible_study_gene_stats(harmonized: pd.DataFrame, min_studies: int) -> pd.DataFrame:
    by_study = _study_gene_stats(harmonized)
    if by_study.empty:
        return by_study

    n_studies = by_study.groupby("gene_symbol")["source_unit_id"].nunique().rename("n_studies")
    by_study = by_study.merge(n_studies, on="gene_symbol", how="left")
    return by_study.loc[by_study["n_studies"].ge(min_studies)].copy()


def stouffer_consensus(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    """Weighted Stouffer consensus over signed-z values."""

    by_study = _eligible_study_gene_stats(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=STOUFFER_COLUMNS)

    by_study["_wz"] = by_study["weight"] * by_study["signed_z"]
    by_study["_w2"] = by_study["weight"] ** 2
    by_study["_wlfc"] = by_study["weight"] * by_study["lfc"].fillna(0.0)
    by_study["_w_lfc_denominator"] = np.where(by_study["lfc"].notna(), by_study["weight"], 0.0)
    grouped = by_study.groupby("gene_symbol", as_index=False).agg(
        n_studies=("source_unit_id", "nunique"),
        sum_wz=("_wz", "sum"),
        sum_w=("weight", "sum"),
        sum_w2=("_w2", "sum"),
        sum_wlfc=("_wlfc", "sum"),
        sum_w_lfc=("_w_lfc_denominator", "sum"),
    )
    grouped["stouffer_z"] = grouped["sum_wz"] / np.sqrt(grouped["sum_w2"])
    grouped["mean_source_z"] = grouped["sum_wz"] / grouped["sum_w"]
    grouped["stouffer_p"] = 2.0 * norm.sf(np.abs(grouped["stouffer_z"]))
    grouped["weighted_lfc"] = np.where(
        grouped["sum_w_lfc"].gt(0),
        grouped["sum_wlfc"] / grouped["sum_w_lfc"],
        np.nan,
    )

    signs = by_study[["gene_symbol", "signed_z"]].merge(
        grouped[["gene_symbol", "stouffer_z"]],
        on="gene_symbol",
        how="left",
    )
    signs["_combined_sign"] = np.sign(signs["stouffer_z"])
    signs["_concordant"] = np.where(
        signs["_combined_sign"].ne(0),
        np.sign(signs["signed_z"]).eq(signs["_combined_sign"]),
        False,
    )
    signs["_direction_strength"] = signs["signed_z"].abs().clip(upper=8.0)
    signs["_concordant_strength"] = np.where(signs["_concordant"], signs["_direction_strength"], 0.0)
    direction = signs.groupby("gene_symbol", as_index=False).agg(
        total_strength=("_direction_strength", "sum"),
        concordant_strength=("_concordant_strength", "sum"),
    )
    direction["sign_concordance"] = np.where(
        direction["total_strength"].gt(0),
        direction["concordant_strength"] / direction["total_strength"],
        0.0,
    )
    sign_concordance = direction[["gene_symbol", "sign_concordance"]]

    heterogeneity = by_study.merge(
        grouped[["gene_symbol", "mean_source_z"]],
        on="gene_symbol",
        how="left",
    )
    heterogeneity["_q_component"] = heterogeneity["weight"] * (
        heterogeneity["signed_z"] - heterogeneity["mean_source_z"]
    ) ** 2
    heterogeneity = heterogeneity.groupby("gene_symbol", as_index=False).agg(
        heterogeneity_q=("_q_component", "sum"),
        heterogeneity_df=("source_unit_id", lambda values: max(int(values.nunique()) - 1, 0)),
    )
    heterogeneity["heterogeneity_i2"] = np.where(
        heterogeneity["heterogeneity_q"].gt(0),
        ((heterogeneity["heterogeneity_q"] - heterogeneity["heterogeneity_df"]) / heterogeneity["heterogeneity_q"]).clip(lower=0.0),
        0.0,
    )

    out = grouped.drop(columns=["sum_wz", "sum_w", "sum_w2", "sum_wlfc", "sum_w_lfc", "mean_source_z"]).merge(
        sign_concordance,
        on="gene_symbol",
        how="left",
    ).merge(
        heterogeneity,
        on="gene_symbol",
        how="left",
    )
    out["stouffer_padj"] = bh_adjust(out["stouffer_p"].to_numpy(dtype=float))
    return out.sort_values(["stouffer_padj", "stouffer_p", "gene_symbol"]).reset_index(drop=True)[STOUFFER_COLUMNS]


def rank_product_consensus(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    """Deterministic rank-product style approximation for the slice.

    This is intentionally labeled as an approximation; the full S1 baseline will
    use RobustRankAggreg via R.
    """

    by_study = _eligible_study_gene_stats(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=RANK_COLUMNS)

    eps = np.finfo(float).tiny
    by_study["_log_rank"] = np.log(by_study["normalized_rank"].clip(lower=eps, upper=1.0))
    out = by_study.groupby("gene_symbol", as_index=False).agg(
        n_studies_rank=("source_unit_id", "nunique"),
        mean_log_rank=("_log_rank", "mean"),
    )
    out["rank_product"] = np.exp(out["mean_log_rank"])
    out["rank_score"] = -np.log(out["rank_product"])
    return out.drop(columns=["mean_log_rank"]).sort_values(["rank_product", "gene_symbol"]).reset_index(drop=True)[RANK_COLUMNS]


def slice_consensus(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    stouffer = stouffer_consensus(harmonized, min_studies=min_studies)
    rank_product = rank_product_consensus(harmonized, min_studies=min_studies)
    if stouffer.empty:
        return stouffer
    merged = stouffer.merge(rank_product, on="gene_symbol", how="left")
    merged["slice_rank"] = np.arange(1, len(merged) + 1)
    return merged
