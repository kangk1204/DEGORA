"""Consensus aggregation methods for the thin DEGORA slice."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
import operator
import re
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
    "stouffer_neglog10_p",
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
    "n_genes_in_study",
    "weight",
    "n_contrast_rows",
    "n_studies_in_source_unit",
]
MAX_SOURCE_SAMPLE_WEIGHT = 4.0
# Per-contrast Stouffer weight. Both group sizes must be present and positive;
# anything else - missing, zero or negative - falls back to 1.0. Note that 1.0 is
# BELOW the sqrt(2) = 1.41 a two-sample contrast earns, so an unstated group size
# is a penalty in this lane, not the neutral value the quality lane gives it.
STOUFFER_WEIGHT_RULE = (
    "per-contrast weight = min(sqrt(n_ctrl + n_treat), "
    f"{MAX_SOURCE_SAMPLE_WEIGHT:g}) when both group sizes are present and positive, otherwise 1.0; "
    "the source-unit weight is the mean of its contrast weights. 1.0 is below the sqrt(2) a "
    "two-sample contrast earns, so an unstated, zero or negative group size lowers a source's "
    "weight here rather than leaving it neutral"
)
SOURCE_UNIT_COLLAPSE_RULE = (
    "source-unit mean aggregation: sample-size-weighted mean signed_z, "
    "sample-size-weighted mean log2FC, mean normalized rank, and mean source weight "
    f"with per-contrast sample-size weights capped at {MAX_SOURCE_SAMPLE_WEIGHT:g}; "
    "no max-|z| representative selection. Optional time_course_mode can preselect "
    "rows within a source unit before this aggregation: mean keeps every contrast; "
    "early and late keep all gene rows from the source unit's globally smallest and largest "
    "numeric duration_h; peak_mean keeps each gene's strongest half by |signed_z| (at least "
    "two, and all observations when that gene has two or fewer), where |signed_z| is "
    "statistical strength derived from the p-value rather than effect size, with ties broken "
    "by study_id."
)


def validate_min_studies(min_studies: Any) -> int:
    """Return a normalized source-unit threshold for consensus APIs."""

    if isinstance(min_studies, (bool, np.bool_)):
        raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}")
    if isinstance(min_studies, (float, np.floating)):
        numeric = float(min_studies)
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}")
        value = int(numeric)
    elif isinstance(min_studies, Decimal):
        if not min_studies.is_finite() or min_studies != min_studies.to_integral_value():
            raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}")
        value = int(min_studies)
    elif isinstance(min_studies, str):
        text = min_studies.strip()
        if not text:
            raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}")
        try:
            numeric = Decimal(text)
        except InvalidOperation:
            raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}") from None
        if not numeric.is_finite() or numeric != numeric.to_integral_value():
            raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}")
        value = int(numeric)
    else:
        try:
            value = operator.index(min_studies)
        except TypeError:
            raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}") from None
    if value < 1:
        raise ValueError(f"min_studies must be an integer >= 1, got {min_studies!r}")
    return value


def _source_unit_series(frame: pd.DataFrame) -> pd.Series:
    # Independent source unit precedence: an explicit source_unit_id wins, then
    # paper_id, then study_id. Honoring an explicit source_unit_id means a caller
    # who already grouped contrasts into source units (e.g. via the public API)
    # gets exactly the grouping they declared instead of a paper_id/study_id
    # re-derivation.
    def _clean(column: str) -> pd.Series:
        values = frame[column].astype("string").fillna("").str.strip()
        # An empty source_unit_id/paper_id column round-trips through CSV/parquet as
        # the literal strings "nan"/"none"/"<NA>"; treat those as blank so a catalog
        # with an all-empty source_unit_id column does not collapse every gene into a
        # single literal-"nan" source unit (which silently drops every gene below
        # min_studies and yields zero scored genes).
        blank = values.eq("") | values.str.lower().isin(["nan", "none", "<na>"])
        return values.mask(blank, "")

    study_id = _clean("study_id")
    result = study_id
    if "paper_id" in frame.columns:
        paper_id = _clean("paper_id")
        result = paper_id.mask(paper_id.eq(""), result)
    if "source_unit_id" in frame.columns:
        source_unit_id = _clean("source_unit_id")
        result = source_unit_id.mask(source_unit_id.eq(""), result)
    return result


def _normalize_time_course_mode(value: Any) -> str:
    if value is None or pd.isna(value):
        return "mean"
    label = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if label in {"", "auto", "all", "mean", "source_mean", "average"}:
        return "mean"
    if label in {"first", "earliest", "early"}:
        return "early"
    if label in {"last", "latest", "late"}:
        return "late"
    if label in {"peak", "peak_mean", "strongest_window"}:
        return "peak_mean"
    raise ValueError(
        f"unsupported time_course_mode={value!r}; expected mean, early, late, or peak_mean"
    )


# A duration is a number of hours, nothing else. "30min" used to parse as 30 and
# "4h" as 4, so `early` kept the 4 h contrast of a unit whose earliest point was
# 30 minutes - and inverted every gene whose direction changed between them.
DURATION_NUMBER_RE = r"^\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\s*$"


def duration_hours(value: Any) -> float:
    """Return a duration in hours, or NaN when the label is not a plain number."""

    if value is None:
        return float("nan")
    try:
        if pd.isna(value):
            return float("nan")
    except (TypeError, ValueError):
        return float("nan")
    text = str(value).strip()
    if not text or not re.match(DURATION_NUMBER_RE, text):
        return float("nan")
    number = pd.to_numeric(pd.Series([text]), errors="coerce").iloc[0]
    return float(number) if pd.notna(number) and np.isfinite(float(number)) else float("nan")


def _duration_numeric(values: pd.Series) -> pd.Series:
    """Parse plain numeric duration labels for configured temporal aggregation."""

    return pd.Series([duration_hours(value) for value in values.tolist()], index=values.index, dtype=float)


def _apply_time_course_mode(frame: pd.DataFrame) -> pd.DataFrame:
    """Preselect the contrasts of a source unit before collapsing it.

    ``mean`` keeps every contrast. ``early`` and ``late`` keep every gene row from
    the source unit's globally smallest and largest numeric ``duration_h``; a source
    unit whose rows carry no parsable duration keeps all of them. ``peak_mean`` is
    gene-specific and keeps the strongest half by ``|signed_z|`` - at least two,
    and all observations when the gene has two or fewer.

    "Strongest" here means statistical strength: ``signed_z`` is derived from the
    p-value, so a time point with a large fold change but a weak p-value is not the
    peak. Ties are broken by ``study_id`` so the selection is deterministic.
    """

    if "time_course_mode" in frame.columns:
        mode_column = "time_course_mode"
    elif "temporal_mode" in frame.columns:
        mode_column = "temporal_mode"
    else:
        mode_column = ""
    if not mode_column:
        frame["time_course_mode"] = "mean"
        return frame
    frame = frame.copy()
    frame["time_course_mode"] = frame[mode_column].map(_normalize_time_course_mode)
    modes_per_unit = frame.groupby("source_unit_id", sort=False)["time_course_mode"].nunique()
    conflicting_units = sorted(modes_per_unit.index[modes_per_unit.gt(1)].astype(str))
    if conflicting_units:
        raise ValueError(
            "each source_unit_id must use one time_course_mode; conflicting source unit(s): "
            + ", ".join(conflicting_units)
        )
    if frame["time_course_mode"].eq("mean").all():
        return frame

    gene_source_cols = ["gene_symbol", "source_unit_id"]
    selected = [frame.loc[frame["time_course_mode"].eq("mean")]]
    duration = _duration_numeric(frame["duration_h"]) if "duration_h" in frame.columns else pd.Series(np.nan, index=frame.index)

    for mode, reducer in [("early", "min"), ("late", "max")]:
        subset = frame.loc[frame["time_course_mode"].eq(mode)].copy()
        if subset.empty:
            continue
        subset["_duration"] = duration.loc[subset.index]
        finite = np.isfinite(subset["_duration"])
        if finite.any():
            targets = subset.loc[finite].groupby("source_unit_id", sort=False)["_duration"].transform(reducer)
            subset.loc[finite, "_target_duration"] = targets
            has_finite = subset.groupby("source_unit_id", sort=False)["_duration"].transform(
                lambda values: np.isfinite(values).any()
            )
            subset = subset.loc[~has_finite | (subset["_duration"].eq(subset["_target_duration"]))]
        selected.append(subset.drop(columns=["_duration", "_target_duration"], errors="ignore"))

    peak = frame.loc[frame["time_course_mode"].eq("peak_mean")].copy()
    if not peak.empty:
        peak["_abs_for_peak"] = peak["signed_z"].abs()
        peak = peak.sort_values(
            [*gene_source_cols, "_abs_for_peak", "study_id"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
        group_size = peak.groupby(gene_source_cols, sort=False)["study_id"].transform("size")
        group_size_values = group_size.to_numpy(dtype=int)
        keep_n = np.where(group_size_values <= 2, group_size_values, np.maximum(2, np.ceil(group_size_values * 0.5))).astype(int)
        rank_from_top = peak.groupby(gene_source_cols, sort=False).cumcount()
        selected.append(peak.loc[rank_from_top < keep_n].drop(columns=["_abs_for_peak"]))

    out = pd.concat([part for part in selected if not part.empty], ignore_index=True)
    return out if not out.empty else frame.iloc[0:0].copy()


TIME_COURSE_RETENTION_WARNING_FRACTION = 0.5


def time_course_selection_report(harmonized: pd.DataFrame) -> list[dict[str, Any]]:
    """Report what `early`/`late` preselection kept and dropped, per source unit.

    The selection is defensible - a gene not measured at the unit's earliest time
    point genuinely has no early observation there - but it was invisible. A unit
    pairing a 200-gene 24h table with a 2-gene 30-minute pilot keeps two rows and
    says nothing, and the genes it drops can fall below min_studies and leave the
    ranking with no warning, no count, and no diagnostic.
    """

    if harmonized.empty or "study_id" not in harmonized.columns:
        return []
    # The harmonized table carries source_unit_id only when the catalog named
    # one; a paper_id-only catalog leaves it blank on every row. Resolving the
    # unit the same way the scorer does keeps this report about the same units
    # the ranking used - grouping every row under "" raised a conflict between
    # papers that never shared a unit, and reported one unit instead of several.
    before = harmonized.copy()
    before["source_unit_id"] = _source_unit_series(before)
    after = _apply_time_course_mode(before.copy())
    if after.empty:
        return []
    modes = (
        after.groupby("source_unit_id", sort=True)["time_course_mode"].first()
        if "time_course_mode" in after.columns
        else pd.Series(dtype=object)
    )
    rows: list[dict[str, Any]] = []
    for unit, mode in modes.items():
        if str(mode) not in {"early", "late"}:
            continue
        kept = after.loc[after["source_unit_id"].eq(unit)]
        original = before.loc[before["source_unit_id"].eq(unit)]
        genes_before = int(original["gene_symbol"].nunique())
        genes_after = int(kept["gene_symbol"].nunique())
        rows.append(
            {
                "source_unit_id": str(unit),
                "time_course_mode": str(mode),
                "rows_before": int(len(original)),
                "rows_after": int(len(kept)),
                "genes_before": genes_before,
                "genes_after": genes_after,
                "gene_retention": (genes_after / genes_before) if genes_before else 1.0,
            }
        )
    return rows


def time_course_selection_warnings(report: list[dict[str, Any]]) -> list[str]:
    """Warn where preselection left a source unit contributing a small minority."""

    warnings: list[str] = []
    for entry in report:
        if entry["genes_before"] <= 0 or entry["gene_retention"] >= TIME_COURSE_RETENTION_WARNING_FRACTION:
            continue
        warnings.append(
            f"time_course_mode={entry['time_course_mode']} left source_unit_id="
            f"{entry['source_unit_id']!r} with {entry['genes_after']} of {entry['genes_before']} genes "
            f"({entry['gene_retention']:.0%}): only its "
            f"{'earliest' if entry['time_course_mode'] == 'early' else 'latest'} timed contrast is used, "
            "and genes measured at no other time in that unit contribute nothing. Genes that fall below "
            "min_studies as a result leave the ranking."
        )
    return warnings


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
            "n_genes_in_study",
            "_weight",
        ]
        for column in required_columns:
            if column not in frame.columns:
                dtype = float if column in {"signed_z", "lfc", "normalized_rank", "n_genes_in_study", "_weight"} else "string"
                frame[column] = pd.Series(dtype=dtype)
        return frame

    frame = harmonized.dropna(subset=["signed_z"]).copy()
    frame["source_unit_id"] = _source_unit_series(frame)
    frame["study_id"] = frame["study_id"].astype("string").fillna("").str.strip()
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.upper().str.strip()
    frame["signed_z"] = pd.to_numeric(frame["signed_z"], errors="coerce")
    frame["lfc"] = pd.to_numeric(frame["lfc"], errors="coerce")
    frame["normalized_rank"] = pd.to_numeric(frame["normalized_rank"], errors="coerce")
    if "n_genes_in_study" in frame.columns:
        frame["n_genes_in_study"] = pd.to_numeric(frame["n_genes_in_study"], errors="coerce")
    else:
        frame["n_genes_in_study"] = np.nan
    # Drop non-finite effect values (mirrors harmonize.py): an inf signed_z would make
    # the Stouffer combination and heterogeneity stats inf/NaN with no error.
    frame.loc[~np.isfinite(frame["signed_z"]), "signed_z"] = np.nan
    frame.loc[~np.isfinite(frame["lfc"]), "lfc"] = np.nan
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
    # Both groups must actually have samples. Testing only the sum let a contrast
    # with no controls at all be weighted like a study of its treatment group -
    # (0, 5) drew sqrt(5) - and let a negative count through as sqrt of whatever
    # the pair summed to. A contrast that cannot state both group sizes falls back
    # to the unweighted 1.0, which is the same thing a missing count gets.
    usable_counts = np.isfinite(n_ctrl) & np.isfinite(n_treat) & (n_ctrl > 0) & (n_treat > 0)
    total_samples = np.where(usable_counts, n_ctrl + n_treat, 1.0)
    frame["_weight"] = np.where(usable_counts, np.sqrt(total_samples), 1.0)
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
        n_genes_in_study=("n_genes_in_study", "max"),
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
    min_studies = validate_min_studies(min_studies)
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
    absolute_z = np.abs(grouped["stouffer_z"].to_numpy(dtype=float))
    log_two_sided_p = np.log(2.0) + norm.logsf(absolute_z)
    grouped["stouffer_p"] = np.exp(log_two_sided_p)
    stouffer_neglog10_p = -log_two_sided_p / np.log(10.0)
    stouffer_neglog10_p[stouffer_neglog10_p == 0.0] = 0.0
    grouped["stouffer_neglog10_p"] = stouffer_neglog10_p
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
    # Primary (unweighted) direction concordance: strength is |signed_z| alone, with
    # no source/effective weight. This is deliberate for the screening slice score.
    # The headline quality-weighted score weights this strength by the effective
    # weight (reliability x sqrt(n)); see score_db._quality_weighted_consensus, which
    # is the `direction = sum_concordant[w*min(|z|,8)] / sum[w*min(|z|,8)]` formula the
    # README documents. Keep the two paths distinct: do not add weight here.
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
    return out.sort_values(
        ["stouffer_padj", "stouffer_p", "stouffer_neglog10_p", "gene_symbol"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)[STOUFFER_COLUMNS]


def rank_product_consensus(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    """Deterministic rank-product style approximation for the slice.

    Note: ``rank_product`` here is the GEOMETRIC MEAN of per-source normalized ranks
    (exp(mean(log(normalized_rank)))), not the classical product-of-ranks statistic
    nor an RRA beta-order-statistic p-value. It is intentionally labeled as an
    approximation; the calibrated rank lane is rra_rho (see score_db._rra_beta_layer),
    and the full S1 baseline uses RobustRankAggreg via R.
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
    out.loc[out["rank_score"] == 0.0, "rank_score"] = 0.0
    return out.drop(columns=["mean_log_rank"]).sort_values(["rank_product", "gene_symbol"]).reset_index(drop=True)[RANK_COLUMNS]


def slice_consensus(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    stouffer = stouffer_consensus(harmonized, min_studies=min_studies)
    rank_product = rank_product_consensus(harmonized, min_studies=min_studies)
    if stouffer.empty:
        return stouffer
    merged = stouffer.merge(rank_product, on="gene_symbol", how="left")
    merged["slice_rank"] = np.arange(1, len(merged) + 1)
    return merged
