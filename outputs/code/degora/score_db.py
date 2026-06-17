"""DEGORA score tables and a local SQLite evidence database."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import beta, norm
from scipy.stats import t as t_dist

from . import runtime_version_info
from .aggregate import (
    SOURCE_UNIT_COLLAPSE_RULE,
    collapse_gene_source_units,
    slice_consensus,
    source_unit_rows_for_aggregation,
)
from .aggregate import _source_unit_series as _aggregate_source_unit_series
from .provenance import shell_command, write_source_sidecar
from .slice_runner import catalog_include_mask, read_catalog
from .stats import bh_adjust


REPO_ROOT = Path(__file__).resolve().parents[3]
SCORE_VERSION = "degora_score_v1_2_source_unit_mean"
SCORE_WEIGHTS = {
    "support_score": 0.30,
    "direction_score": 0.25,
    "evidence_score": 0.20,
    "rank_score_component": 0.15,
    "effect_score": 0.10,
}

PRIMARY_RANK_COLUMN = "quality_weighted_degora_rank"
PRIMARY_SCORE_COLUMN = "quality_weighted_degora_score"
PRIMARY_TOP_PERCENT_COLUMN = "quality_weighted_top_percent"
PRIMARY_DIRECTION_COLUMN = "quality_weighted_consensus_direction"
PRIMARY_CONCORDANCE_COLUMN = "quality_weighted_sign_concordance"
PRIMARY_RANK_DESCRIPTION = (
    "quality_weighted_degora_rank is the manuscript-facing primary rank. "
    "degora_rank and degora_score are retained as unweighted/reference outputs."
)
PRIORITY_SCORE_WEIGHTS = {
    "direction_score": 0.25,
    "evidence_score": 0.30,
    "rank_score_component": 0.30,
    "effect_score": 0.15,
}


def _portable_cli_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)

RELIABILITY_SCORE_WEIGHTS = {
    "support_score": 0.25,
    "source_quality_support_score": 0.25,
    "direction_confidence_index": 0.25,
    "loo_rank_stability_score": 0.25,
}
SCORE_FORMULA = (
    "100 * weighted_geometric_mean(support_score, direction_score, "
    "evidence_score, rank_score_component, effect_score); support is "
    "log-scaled by independent source units, direction is sign concordance, "
    "evidence is Stouffer-z strength, rank is inverse rank product, and "
    "effect is absolute weighted log2FC strength. The score is a transparent "
    "ranking aid, not a calibrated posterior probability. Consensus evidence "
    "is combined after aggregating related contrasts within each independent "
    "source unit without max-|z| representative selection."
)
QUALITY_WEIGHTED_SCORE_FORMULA = (
    "Fixed manuscript-facing source-quality-weighted ranking: same component "
    "formula as the unweighted DEGORA score, but source-unit evidence is weighted "
    "by predeclared source-quality features (source input type, table scope, and "
    "replicate count) together with a gold-panel-free source-coherence guardrail. "
    "It references no marker or gold-panel genes and is panel-agnostic by "
    "construction; it is reported as a relative prioritization index, not a "
    "calibrated probability."
)
SOURCE_INPUT_TYPE_QUALITY_WEIGHTS = {
    "author_deg_table": 1.00,
    "limma_full_table": 0.90,
    "microarray_author_deg_table": 0.90,
    "derived_count_table": 0.85,
    "normalized_expression_matrix": 0.35,
}
TABLE_SCOPE_QUALITY_MULTIPLIERS = {
    "full_results": 1.00,
    "auto": 0.90,
    "ambiguous": 0.75,
    "deg_only": 0.65,
}
EFFECT_META_MIN_WEIGHT_SUM = 1e-15

GENE_SCORE_COLUMNS = [
    "degora_rank",
    "rank_label",
    "gene_symbol",
    "evidence_tier",
    "degora_score",
    "top_percent",
    "percentile",
    "top_percent_label",
    "consensus_direction",
    "n_source_units",
    "n_contrasts_observed",
    "support_label",
    "source_units",
    "sign_concordance",
    "direction_label",
    "support_score",
    "direction_score",
    "evidence_score",
    "rank_score_component",
    "effect_score",
    "priority_rank",
    "priority_score",
    "priority_top_percent",
    "evidence_reliability_score",
    "direction_confidence_index",
    "quality_weighted_direction_confidence_index",
    "direction_concordant_source_units",
    "direction_total_source_units",
    "direction_posterior_mean",
    "loo_median_rank",
    "loo_rank_iqr",
    "loo_rank_stability_score",
    "loo_top50_fraction",
    "loo_top100_fraction",
    "quality_weighted_degora_rank",
    "quality_weighted_degora_score",
    "quality_weighted_top_percent",
    "quality_weighted_consensus_direction",
    "quality_weighted_sign_concordance",
    "source_quality_support_score",
    "source_quality_weight_sum",
    "stouffer_z",
    "stouffer_p",
    "stouffer_padj",
    "heterogeneity_q",
    "heterogeneity_df",
    "heterogeneity_i2",
    "heterogeneity_flag",
    "re_stouffer_z",
    "re_stouffer_p",
    "re_stouffer_padj",
    "re_stouffer_shrinkage_factor",
    "rra_rho",
    "rra_neglog10_rho",
    "rra_rank",
    "effect_meta_log2fc_re",
    "effect_meta_se",
    "effect_meta_ci_low",
    "effect_meta_ci_high",
    "effect_meta_tau2",
    "effect_meta_i2",
    "effect_meta_k",
    "effect_meta_se_source",
    "effect_meta_exact_weight_fraction",
    "weighted_lfc",
    "rank_product",
    "rank_score",
    "slice_rank",
    "high_confidence",
]


def primary_ranked_scores(frame: pd.DataFrame) -> pd.DataFrame:
    """Return scores in the user-facing primary rank order.

    The schema keeps both the original unweighted DEGORA rank and the
    quality-weighted rank. Output channels should use this helper when they need
    an ordered top-gene table, while preserving the original columns unchanged.
    """

    out = frame.copy()
    if out.empty:
        return out.reset_index(drop=True)
    if PRIMARY_RANK_COLUMN in out.columns:
        out["_primary_rank_sort"] = pd.to_numeric(out[PRIMARY_RANK_COLUMN], errors="coerce")
        out["_primary_rank_sort"] = out["_primary_rank_sort"].where(out["_primary_rank_sort"].gt(0))
        if out["_primary_rank_sort"].notna().any():
            sort_columns = ["_primary_rank_sort"]
            ascending = [True]
            if "gene_symbol" in out.columns:
                sort_columns.append("gene_symbol")
                ascending.append(True)
            return out.sort_values(sort_columns, ascending=ascending, na_position="last").drop(
                columns=["_primary_rank_sort"]
            ).reset_index(drop=True)
        out = out.drop(columns=["_primary_rank_sort"])
    if "degora_rank" in out.columns:
        out["_fallback_rank_sort"] = pd.to_numeric(out["degora_rank"], errors="coerce")
        out["_fallback_rank_sort"] = out["_fallback_rank_sort"].where(out["_fallback_rank_sort"].gt(0))
        sort_columns = ["_fallback_rank_sort"]
        ascending = [True]
        if "gene_symbol" in out.columns:
            sort_columns.append("gene_symbol")
            ascending.append(True)
        return out.sort_values(sort_columns, ascending=ascending, na_position="last").drop(
            columns=["_fallback_rank_sort"]
        ).reset_index(drop=True)
    if "gene_symbol" in out.columns:
        return out.sort_values("gene_symbol").reset_index(drop=True)
    return out.reset_index(drop=True)

GENE_EVIDENCE_COLUMNS = [
    "gene_symbol",
    "study_id",
    "source_unit_id",
    "paper_id",
    "pipeline",
    "assay_type",
    "source_input_type",
    "table_scope",
    "platform",
    "normalization",
    "probe_collapse",
    "species",
    "cell_system",
    "hypoxia_modality",
    "duration_h",
    "time_course_mode",
    "temporal_mode",
    "n_ctrl",
    "n_treat",
    "lfc",
    "signed_z",
    "aggregate_pvalue",
    "aggregate_padj",
    "min_source_pvalue",
    "min_source_padj",
    "normalized_rank",
    "n_genes_in_study",
    "weight",
    "source_quality_weight",
    "source_quality_label",
    "source_coherence_weight",
    "source_recommended_weight",
    "source_reliability_weight",
    "source_reliability_label",
    "source_outlier_flag",
    "direction",
    "source_path",
    "source_url",
    "contributing_study_ids",
    "contributing_pipelines",
    "contributing_assay_types",
    "contributing_source_input_types",
    "contributing_platforms",
    "contributing_normalizations",
    "contributing_probe_collapse",
    "contributing_duration_h",
    "contributing_time_course_modes",
    "contributing_source_paths",
    "contributing_source_urls",
    "n_contrast_rows",
    "n_studies_in_source_unit",
]
SOURCE_QUALITY_DIAGNOSTIC_COLUMNS = [
    "source_unit_id",
    "source_input_type",
    "assay_type",
    "pipeline",
    "n_genes",
    "n_pairwise_comparisons",
    "median_pairwise_lfc_spearman",
    "min_pairwise_lfc_spearman",
    "median_pairwise_sign_agreement",
    "source_quality_weight",
    "source_quality_label",
    "source_coherence_weight",
    "source_recommended_weight",
    "source_reliability_weight",
    "source_reliability_label",
    "source_outlier_flag",
    "recommended_role",
]


def _score_ready_harmonized(harmonized: pd.DataFrame, *, lfc_cap: float = 10.0) -> tuple[pd.DataFrame, int]:
    """Return a scoring copy with non-finite LFC values capped for display math."""

    frame = harmonized.copy()
    if "lfc" not in frame.columns:
        return frame, 0
    lfc = pd.to_numeric(frame["lfc"], errors="coerce")
    nonfinite = np.isinf(lfc.to_numpy(dtype=float))
    if not nonfinite.any():
        frame["lfc"] = lfc
        return frame, 0
    signs = np.sign(lfc.loc[nonfinite].to_numpy(dtype=float))
    signs = np.where(signs == 0, 1.0, signs)
    frame.loc[nonfinite, "lfc"] = signs * float(lfc_cap)
    frame["lfc"] = pd.to_numeric(frame["lfc"], errors="coerce")
    return frame, int(nonfinite.sum())


def _as_numeric(frame: pd.DataFrame, column: str, default: float = np.nan) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _string_column(frame: pd.DataFrame, column: str, default: str = "") -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype="string")
    return frame[column].astype("string").fillna(default)


def _source_unit_series(frame: pd.DataFrame) -> pd.Series:
    # Single source-unit definition program-wide: reuse the canonical aggregate
    # precedence (explicit source_unit_id wins, then paper_id, then study_id) so
    # scoring, evidence, support, and metadata describe the same object.
    return _aggregate_source_unit_series(frame)


def _join_unique(values: pd.Series) -> str:
    labels = []
    for value in values.dropna().astype(str):
        label = value.strip()
        if label:
            labels.append(label)
    return ";".join(sorted(dict.fromkeys(labels)))


def _min_numeric(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    return float(numeric.min()) if not numeric.empty else np.nan


def _clean_join_label(value: Any) -> str:
    if pd.isna(value):
        return ""
    label = str(value).strip()
    return "" if label.lower() in {"nan", "<na>"} else label


def _source_input_type_weight(value: Any) -> float:
    label = _clean_join_label(value).lower()
    if ";" in label:
        parts = [part.strip() for part in label.split(";") if part.strip()]
        return min((_source_input_type_weight(part) for part in parts), default=0.65)
    if label in SOURCE_INPUT_TYPE_QUALITY_WEIGHTS:
        return SOURCE_INPUT_TYPE_QUALITY_WEIGHTS[label]
    if "author" in label and "deg" in label and "table" in label:
        return 1.0
    if "derived" in label and "count" in label:
        return 0.85
    if "normalized" in label and "matrix" in label:
        return 0.35
    if not label:
        return 0.65
    return 0.70


def _table_scope_multiplier(value: Any) -> float:
    label = _clean_join_label(value).lower()
    if ";" in label:
        parts = [part.strip() for part in label.split(";") if part.strip()]
        return min((_table_scope_multiplier(part) for part in parts), default=0.75)
    return TABLE_SCOPE_QUALITY_MULTIPLIERS.get(label, 0.85 if label else 0.90)


def _replicate_quality_multiplier(n_ctrl: Any, n_treat: Any) -> float:
    ctrl = pd.to_numeric(pd.Series([n_ctrl]), errors="coerce").iloc[0]
    treat = pd.to_numeric(pd.Series([n_treat]), errors="coerce").iloc[0]
    if not np.isfinite(ctrl) or not np.isfinite(treat):
        return 0.75
    minimum = min(float(ctrl), float(treat))
    if minimum >= 3:
        return 1.00
    if minimum >= 2:
        return 0.85
    if minimum >= 1:
        return 0.50
    return 0.35


def _quality_label(weight: Any) -> str:
    value = pd.to_numeric(pd.Series([weight]), errors="coerce").iloc[0]
    if not np.isfinite(value):
        return "unknown"
    if value >= 0.85:
        return "high"
    if value >= 0.60:
        return "medium"
    return "low"


def _source_reliability_weight(
    source_quality_weight: float,
    *,
    source_coherence_weight: float = 1.0,
    n_pairwise_comparisons: int = 0,
    n_genes: int = 0,
    neutral_prior: float = 0.65,
    prior_strength: float = 1.0,
) -> float:
    """Shrink source quality toward a neutral prior before secondary weighting.

    This is an index for ranking and sensitivity analysis, not a probability.
    Static source-quality rules provide the observed value; pairwise coherence
    and gene coverage only control how strongly we trust that observed value.

    The shrinkage is symmetric: it pulls the observed quality toward the neutral
    prior (0.65, a deliberately conservative midpoint of the predeclared
    source_input_type weights; prior_strength=1 is a weak prior), so a very
    low-quality, low-coverage source is tempered UPWARD toward neutral rather
    than suppressed further. It is a sensitivity-review weight, reported
    alongside (not in place of) the static source-quality weight.
    """

    observed = float(np.clip(source_quality_weight * source_coherence_weight, 0.05, 1.0))
    coverage_strength = min(6.0, np.log1p(max(int(n_genes), 0)) / 2.0)
    comparison_strength = min(3.0, max(int(n_pairwise_comparisons), 0))
    observed_strength = max(1.0, coverage_strength + comparison_strength)
    shrunk = (prior_strength * neutral_prior + observed_strength * observed) / (prior_strength + observed_strength)
    return float(np.clip(shrunk, 0.05, 1.0))


def _source_quality_weight_frame(frame: pd.DataFrame) -> pd.Series:
    source_type = frame["source_input_type"] if "source_input_type" in frame.columns else pd.Series("", index=frame.index)
    table_scope = frame["table_scope"] if "table_scope" in frame.columns else pd.Series("", index=frame.index)
    n_ctrl = frame["n_ctrl"] if "n_ctrl" in frame.columns else pd.Series(np.nan, index=frame.index)
    n_treat = frame["n_treat"] if "n_treat" in frame.columns else pd.Series(np.nan, index=frame.index)
    weights = [
        _source_input_type_weight(source) * _table_scope_multiplier(scope) * _replicate_quality_multiplier(ctrl, treat)
        for source, scope, ctrl, treat in zip(source_type, table_scope, n_ctrl, n_treat, strict=False)
    ]
    return pd.Series(weights, index=frame.index, dtype=float).clip(0.05, 1.0)


def _metadata_for_study_gene_units(frame: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["gene_symbol", "source_unit_id"]
    source_unit_string_specs = [
        ("paper_id", "paper_id"),
        ("pipeline", "pipeline"),
        ("assay_type", "assay_type"),
        ("source_input_type", "source_input_type"),
        ("platform", "platform"),
        ("normalization", "normalization"),
        ("probe_collapse", "probe_collapse"),
        ("species", "species"),
        ("cell_system", "cell_system"),
        ("hypoxia_modality", "hypoxia_modality"),
        ("table_scope", "table_scope"),
        ("source_path", "source_path"),
        ("source_url", "source_url"),
        ("contributing_pipelines", "pipeline"),
        ("contributing_assay_types", "assay_type"),
        ("contributing_source_input_types", "source_input_type"),
        ("contributing_platforms", "platform"),
        ("contributing_normalizations", "normalization"),
        ("contributing_probe_collapse", "probe_collapse"),
    ]
    gene_source_string_specs = [
        ("duration_h", "duration_h"),
        ("time_course_mode", "time_course_mode"),
        ("temporal_mode", "temporal_mode"),
        ("contributing_study_ids", "study_id"),
        ("contributing_duration_h", "duration_h"),
        ("contributing_time_course_modes", "time_course_mode"),
        ("contributing_source_paths", "source_path"),
        ("contributing_source_urls", "source_url"),
    ]
    output_columns = [
        *group_cols,
        *[target for target, _source in source_unit_string_specs],
        *[target for target, _source in gene_source_string_specs],
        "n_ctrl",
        "n_treat",
        "n_genes_in_study",
        "min_source_pvalue",
        "min_source_padj",
    ]
    sorted_frame = frame.sort_values([*group_cols, "study_id"])
    base = sorted_frame.groupby(group_cols, as_index=False, sort=False).agg(
        n_ctrl=("n_ctrl", _min_numeric),
        n_treat=("n_treat", _min_numeric),
        n_genes_in_study=("n_genes_in_study", "max"),
        min_source_pvalue=("pvalue", "min"),
        min_source_padj=("padj", "min"),
    )

    source_unit_meta = sorted_frame[["source_unit_id"]].drop_duplicates().copy()
    for target, source in source_unit_string_specs:
        values = sorted_frame[["source_unit_id", source]].copy()
        values[source] = values[source].map(_clean_join_label)
        values = values.loc[values[source].ne("")].drop_duplicates().sort_values(["source_unit_id", source])
        joined = values.groupby("source_unit_id", sort=False)[source].agg(";".join).rename(target).reset_index()
        source_unit_meta = source_unit_meta.merge(joined, on="source_unit_id", how="left")

    out = base.merge(source_unit_meta, on="source_unit_id", how="left")
    for target, source in gene_source_string_specs:
        values = sorted_frame[[*group_cols, source]].copy()
        values[source] = values[source].map(_clean_join_label)
        values = values.loc[values[source].ne("")].drop_duplicates().sort_values([*group_cols, source])
        joined = values.groupby(group_cols, sort=False)[source].agg(";".join).rename(target).reset_index()
        out = out.merge(joined, on=group_cols, how="left")

    for column in output_columns:
        if column not in out.columns:
            out[column] = ""
    string_columns = [column for column in output_columns if column not in {*group_cols, "n_ctrl", "n_treat", "n_genes_in_study", "min_source_pvalue", "min_source_padj"}]
    out[string_columns] = out[string_columns].fillna("")
    return out[output_columns].reset_index(drop=True)


def study_gene_evidence(harmonized: pd.DataFrame) -> pd.DataFrame:
    """Collapse harmonized rows to one auditable row per gene and source unit."""

    required = {"study_id", "gene_symbol", "lfc", "signed_z", "pvalue", "normalized_rank"}
    missing = sorted(required.difference(harmonized.columns))
    if missing:
        raise ValueError(f"harmonized table is missing score columns: {missing}")

    frame = harmonized.copy()
    frame["gene_symbol"] = _string_column(frame, "gene_symbol").str.upper().str.strip()
    frame["study_id"] = _string_column(frame, "study_id").str.strip()
    frame["paper_id"] = _string_column(frame, "paper_id").str.strip()
    frame["source_unit_id"] = _source_unit_series(frame)
    frame["lfc"] = _as_numeric(frame, "lfc")
    frame["signed_z"] = _as_numeric(frame, "signed_z")
    frame["pvalue"] = _as_numeric(frame, "pvalue")
    frame["padj"] = _as_numeric(frame, "padj")
    frame["normalized_rank"] = _as_numeric(frame, "normalized_rank")
    frame["n_ctrl"] = _as_numeric(frame, "n_ctrl")
    frame["n_treat"] = _as_numeric(frame, "n_treat")
    frame["n_genes_in_study"] = _as_numeric(frame, "n_genes_in_study")
    frame = frame.dropna(subset=["gene_symbol", "study_id", "lfc", "signed_z", "pvalue", "normalized_rank"])
    frame = frame.loc[frame["gene_symbol"].ne("") & frame["study_id"].ne("")].copy()

    n_total = frame["n_ctrl"] + frame["n_treat"]
    frame["_weight"] = np.where(np.isfinite(n_total) & n_total.gt(0), np.sqrt(n_total), 1.0)

    for column, fallback in [
        ("pipeline", "unknown_pipeline"),
        ("assay_type", ""),
        ("source_input_type", ""),
        ("platform", ""),
        ("normalization", ""),
        ("probe_collapse", ""),
        ("species", ""),
        ("cell_system", ""),
        ("hypoxia_modality", ""),
        ("duration_h", ""),
        ("time_course_mode", "mean"),
        ("temporal_mode", ""),
        ("table_scope", ""),
        ("source_path", ""),
        ("source_url", ""),
    ]:
        if column not in frame.columns:
            frame[column] = fallback

    selected_frame = source_unit_rows_for_aggregation(frame)
    collapsed = collapse_gene_source_units(frame)
    meta = _metadata_for_study_gene_units(selected_frame)
    out = collapsed.merge(meta.drop(columns=["n_genes_in_study"], errors="ignore"), on=["gene_symbol", "source_unit_id"], how="left")
    out["aggregate_pvalue"] = 2.0 * norm.sf(np.abs(pd.to_numeric(out["signed_z"], errors="coerce")))
    out["aggregate_padj"] = np.nan
    out["source_quality_weight"] = _source_quality_weight_frame(out)
    out["source_quality_label"] = out["source_quality_weight"].map(_quality_label)
    out["source_coherence_weight"] = 1.0
    out["source_recommended_weight"] = out["source_quality_weight"]
    out["source_reliability_weight"] = out["source_quality_weight"]
    out["source_reliability_label"] = out["source_reliability_weight"].map(_quality_label)
    out["source_outlier_flag"] = False
    out["direction"] = np.select([out["lfc"].gt(0), out["lfc"].lt(0)], ["up", "down"], default="flat")
    return out.sort_values(["gene_symbol", "source_unit_id", "study_id"]).reset_index(drop=True)[GENE_EVIDENCE_COLUMNS]


def _component_strength_from_z(values: pd.Series) -> pd.Series:
    z = pd.to_numeric(values, errors="coerce").abs().fillna(0.0)
    return pd.Series(1.0 - np.exp(-z / 8.0), index=values.index).clip(0.0, 1.0)


def _component_strength_from_lfc(values: pd.Series) -> pd.Series:
    lfc = pd.to_numeric(values, errors="coerce").abs().fillna(0.0)
    return pd.Series(1.0 - np.exp(-lfc / 2.0), index=values.index).clip(0.0, 1.0)


def _weighted_geometric_score_with_weights(frame: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    eps = 1e-6
    total_weight = float(sum(weights.values()))
    score_log = np.zeros(len(frame), dtype=float)
    for column, weight in weights.items():
        component = pd.to_numeric(frame[column], errors="coerce").fillna(0.0).clip(eps, 1.0)
        score_log += float(weight) * np.log(component.to_numpy(dtype=float))
    return pd.Series(100.0 * np.exp(score_log / total_weight), index=frame.index)


def _weighted_geometric_score(frame: pd.DataFrame) -> pd.Series:
    return _weighted_geometric_score_with_weights(frame, SCORE_WEIGHTS)


def _weighted_geometric_score_from_components(frame: pd.DataFrame, prefix: str) -> pd.Series:
    components = pd.DataFrame(index=frame.index)
    for column in SCORE_WEIGHTS:
        quality_column = f"{prefix}{column}"
        if quality_column not in frame.columns:
            components[column] = 0.0
        else:
            components[column] = frame[quality_column]
    return _weighted_geometric_score(components)


def _source_quality_diagnostics_from_evidence(evidence: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame(columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS)

    frame = evidence.copy()
    frame["source_unit_id"] = _string_column(frame, "source_unit_id").str.strip()
    frame["gene_symbol"] = _string_column(frame, "gene_symbol").str.upper().str.strip()
    frame["lfc"] = _as_numeric(frame, "lfc")
    frame["source_quality_weight"] = _as_numeric(frame, "source_quality_weight", default=0.65).fillna(0.65)
    frame = frame.dropna(subset=["source_unit_id", "gene_symbol", "lfc"])
    frame = frame.loc[frame["source_unit_id"].ne("") & frame["gene_symbol"].ne("")].copy()

    if frame.empty:
        return pd.DataFrame(columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS)

    rows: list[dict[str, Any]] = []
    wide_lfc = frame.pivot_table(index="gene_symbol", columns="source_unit_id", values="lfc", aggfunc="mean")
    sign_frame = frame.assign(_sign=np.sign(frame["lfc"]))
    wide_sign = sign_frame.pivot_table(index="gene_symbol", columns="source_unit_id", values="_sign", aggfunc="mean")
    source_units = sorted(frame["source_unit_id"].dropna().astype(str).unique())
    pairwise: dict[str, list[dict[str, float]]] = {source_unit: [] for source_unit in source_units}
    for index, source_a in enumerate(source_units):
        for source_b in source_units[index + 1 :]:
            overlap = wide_lfc[[source_a, source_b]].dropna()
            sign_overlap = wide_sign[[source_a, source_b]].dropna()
            if len(overlap) >= 3:
                lfc_spearman = float(overlap[source_a].corr(overlap[source_b], method="spearman"))
            else:
                lfc_spearman = np.nan
            if len(sign_overlap) > 0:
                sign_agreement = float((sign_overlap[source_a] * sign_overlap[source_b] > 0).mean())
            else:
                sign_agreement = np.nan
            for source_unit in (source_a, source_b):
                pairwise[source_unit].append(
                    {
                        "lfc_spearman": lfc_spearman,
                        "sign_agreement": sign_agreement,
                        "overlap": float(len(overlap)),
                    }
                )

    grouped = frame.groupby("source_unit_id", sort=True)
    for source_unit, group in grouped:
        comparisons = pairwise.get(str(source_unit), [])
        lfc_corrs = [item["lfc_spearman"] for item in comparisons if np.isfinite(item["lfc_spearman"])]
        sign_agreements = [item["sign_agreement"] for item in comparisons if np.isfinite(item["sign_agreement"])]
        median_spearman = float(np.median(lfc_corrs)) if lfc_corrs else np.nan
        min_spearman = float(np.min(lfc_corrs)) if lfc_corrs else np.nan
        median_sign = float(np.median(sign_agreements)) if sign_agreements else np.nan
        source_quality = float(pd.to_numeric(group["source_quality_weight"], errors="coerce").median())
        low_static_quality = source_quality < 0.60
        outlier_flag = bool(
            len(lfc_corrs) >= 2
            and low_static_quality
            and (np.isfinite(median_spearman) and median_spearman < 0.05)
        )
        coherence_weight = 0.50 if outlier_flag else 1.00
        recommended_weight = max(0.05, min(1.0, source_quality * coherence_weight))
        reliability_weight = _source_reliability_weight(
            source_quality,
            source_coherence_weight=coherence_weight,
            n_pairwise_comparisons=len(lfc_corrs),
            n_genes=int(group["gene_symbol"].nunique()),
        )
        rows.append(
            {
                "source_unit_id": str(source_unit),
                "source_input_type": _join_unique(group.get("source_input_type", pd.Series(dtype=object))),
                "assay_type": _join_unique(group.get("assay_type", pd.Series(dtype=object))),
                "pipeline": _join_unique(group.get("pipeline", pd.Series(dtype=object))),
                "n_genes": int(group["gene_symbol"].nunique()),
                "n_pairwise_comparisons": int(len(lfc_corrs)),
                "median_pairwise_lfc_spearman": median_spearman,
                "min_pairwise_lfc_spearman": min_spearman,
                "median_pairwise_sign_agreement": median_sign,
                "source_quality_weight": source_quality,
                "source_quality_label": _quality_label(source_quality),
                "source_coherence_weight": coherence_weight,
                "source_recommended_weight": recommended_weight,
                "source_reliability_weight": reliability_weight,
                "source_reliability_label": _quality_label(reliability_weight),
                "source_outlier_flag": outlier_flag,
                "recommended_role": "sensitivity" if outlier_flag or source_quality < 0.60 else "primary",
            }
        )

    return pd.DataFrame.from_records(rows, columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS)


def _attach_source_quality_diagnostics(evidence: pd.DataFrame, diagnostics: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty or diagnostics.empty:
        return evidence
    mapping = diagnostics.set_index("source_unit_id")
    out = evidence.copy()
    for column in ["source_coherence_weight", "source_recommended_weight", "source_reliability_weight", "source_outlier_flag"]:
        out[column] = out["source_unit_id"].map(mapping[column]).fillna(out[column] if column in out.columns else 1.0)
    if "source_quality_weight" in mapping.columns:
        out["source_quality_weight"] = out["source_unit_id"].map(mapping["source_quality_weight"]).fillna(out["source_quality_weight"])
    out["source_quality_label"] = out["source_quality_weight"].map(_quality_label)
    out["source_reliability_label"] = out["source_reliability_weight"].map(_quality_label)
    return out


def _quality_weighted_consensus(
    evidence: pd.DataFrame,
    *,
    total_source_quality_weight: float,
) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame()

    frame = evidence.copy()
    if "source_reliability_weight" not in frame.columns:
        frame["source_reliability_weight"] = frame.get("source_recommended_weight", 0.65)
    for column in [
        "signed_z",
        "lfc",
        "normalized_rank",
        "weight",
        "source_recommended_weight",
        "source_reliability_weight",
        "source_quality_weight",
    ]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "signed_z", "normalized_rank", "weight", "source_reliability_weight"])
    if frame.empty:
        return pd.DataFrame()

    frame["_effective_weight"] = (frame["weight"] * frame["source_reliability_weight"]).clip(lower=0.05)
    frame["_source_quality"] = frame["source_reliability_weight"].clip(lower=0.05)
    frame["_wz"] = frame["_effective_weight"] * frame["signed_z"]
    frame["_w2"] = frame["_effective_weight"] ** 2
    frame["_wlfc"] = frame["_effective_weight"] * frame["lfc"].fillna(0.0)
    frame["_w_lfc_denominator"] = np.where(frame["lfc"].notna(), frame["_effective_weight"], 0.0)
    eps = np.finfo(float).tiny
    frame["_log_rank"] = np.log(frame["normalized_rank"].clip(lower=eps, upper=1.0))
    frame["_weighted_log_rank"] = frame["_effective_weight"] * frame["_log_rank"]

    grouped = frame.groupby("gene_symbol", as_index=False).agg(
        quality_sum_wz=("_wz", "sum"),
        quality_sum_w2=("_w2", "sum"),
        quality_sum_wlfc=("_wlfc", "sum"),
        quality_sum_w_lfc=("_w_lfc_denominator", "sum"),
        quality_sum_weighted_log_rank=("_weighted_log_rank", "sum"),
        quality_sum_effective_weight=("_effective_weight", "sum"),
        source_quality_weight_sum=("_source_quality", "sum"),
        n_quality_source_units=("source_unit_id", "nunique"),
    )
    grouped["quality_stouffer_z"] = grouped["quality_sum_wz"] / np.sqrt(grouped["quality_sum_w2"])
    grouped["quality_weighted_lfc"] = np.where(
        grouped["quality_sum_w_lfc"].gt(0),
        grouped["quality_sum_wlfc"] / grouped["quality_sum_w_lfc"],
        np.nan,
    )
    grouped["quality_rank_product"] = np.exp(
        grouped["quality_sum_weighted_log_rank"] / grouped["quality_sum_effective_weight"]
    )
    denominator = np.log1p(total_source_quality_weight) if total_source_quality_weight > 0 else 1.0
    grouped["source_quality_support_score"] = (
        np.log1p(grouped["source_quality_weight_sum"]) / denominator
    ).clip(0.0, 1.0)
    grouped["quality_evidence_score"] = _component_strength_from_z(grouped["quality_stouffer_z"])
    grouped["quality_rank_score_component"] = (1.0 - grouped["quality_rank_product"]).fillna(0.0).clip(0.0, 1.0)
    grouped["quality_effect_score"] = _component_strength_from_lfc(grouped["quality_weighted_lfc"])
    # Direction matches quality_weighted_sign_concordance, which is scored against
    # sign(quality_stouffer_z) below; deriving it from quality_weighted_lfc would
    # contradict the concordance for genes where z and LFC disagree in sign.
    grouped["quality_weighted_consensus_direction"] = np.select(
        [grouped["quality_stouffer_z"].gt(0), grouped["quality_stouffer_z"].lt(0)],
        ["up", "down"],
        default="flat",
    )

    signs = frame[["gene_symbol", "signed_z", "_effective_weight"]].merge(
        grouped[["gene_symbol", "quality_stouffer_z"]],
        on="gene_symbol",
        how="left",
    )
    signs["_combined_sign"] = np.sign(signs["quality_stouffer_z"])
    signs["_concordant"] = np.where(
        signs["_combined_sign"].ne(0),
        np.sign(signs["signed_z"]).eq(signs["_combined_sign"]),
        False,
    )
    signs["_direction_strength"] = signs["_effective_weight"] * signs["signed_z"].abs().clip(upper=8.0)
    signs["_concordant_strength"] = np.where(signs["_concordant"], signs["_direction_strength"], 0.0)
    direction = signs.groupby("gene_symbol", as_index=False).agg(
        total_strength=("_direction_strength", "sum"),
        concordant_strength=("_concordant_strength", "sum"),
    )
    direction["quality_weighted_sign_concordance"] = np.where(
        direction["total_strength"].gt(0),
        direction["concordant_strength"] / direction["total_strength"],
        0.0,
    )
    grouped = grouped.merge(direction[["gene_symbol", "quality_weighted_sign_concordance"]], on="gene_symbol", how="left")
    grouped["quality_direction_score"] = grouped["quality_weighted_sign_concordance"].fillna(0.0).clip(0.0, 1.0)
    grouped["quality_support_score"] = grouped["source_quality_support_score"]
    grouped["quality_weighted_degora_score"] = _weighted_geometric_score_from_components(grouped, "quality_")
    return grouped[
        [
            "gene_symbol",
            "quality_stouffer_z",
            "quality_weighted_degora_score",
            "quality_weighted_consensus_direction",
            "quality_weighted_sign_concordance",
            "source_quality_support_score",
            "source_quality_weight_sum",
        ]
    ]


def _direction_confidence_from_evidence(
    evidence: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    reference_z_column: str,
    use_reliability_weight: bool = False,
    output_column: str = "direction_confidence_index",
) -> pd.DataFrame:
    """Return a beta-binomial direction consistency index by gene.

    The primary index uses source-unit counts: x concordant source units out of
    k observed source units, with a Beta(1, 1) prior. This makes the denominator
    interpretable and avoids the older strength-sum denominator. The
    quality-weighted variant uses reliability-weighted pseudo-counts and remains
    labeled as an index rather than a calibrated posterior probability.
    """

    if evidence.empty:
        columns = ["gene_symbol", output_column]
        if output_column == "direction_confidence_index":
            columns.extend(["direction_concordant_source_units", "direction_total_source_units", "direction_posterior_mean"])
        return pd.DataFrame(columns=columns)
    frame = evidence.copy()
    if use_reliability_weight and "source_reliability_weight" not in frame.columns:
        frame["source_reliability_weight"] = frame.get("source_recommended_weight", 1.0)
    for column in ["signed_z", "weight", "source_reliability_weight"]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "signed_z"])
    frame = frame.loc[np.sign(frame["signed_z"]).ne(0)].copy()
    if frame.empty:
        columns = ["gene_symbol", output_column]
        if output_column == "direction_confidence_index":
            columns.extend(["direction_concordant_source_units", "direction_total_source_units", "direction_posterior_mean"])
        return pd.DataFrame(columns=columns)
    if use_reliability_weight:
        frame["_unit_weight"] = frame["source_reliability_weight"].fillna(1.0).clip(lower=0.05)
    else:
        frame["_unit_weight"] = 1.0
    if reference_z_column not in reference.columns:
        raise ValueError(f"reference is missing {reference_z_column!r}")
    reference_sign = reference[["gene_symbol", reference_z_column]].copy()
    reference_sign["_combined_sign"] = np.sign(pd.to_numeric(reference_sign[reference_z_column], errors="coerce"))
    signed = frame.merge(reference_sign[["gene_symbol", "_combined_sign"]], on="gene_symbol", how="inner")
    if signed.empty:
        columns = ["gene_symbol", output_column]
        if output_column == "direction_confidence_index":
            columns.extend(["direction_concordant_source_units", "direction_total_source_units", "direction_posterior_mean"])
        return pd.DataFrame(columns=columns)
    signed["_concordant"] = np.sign(signed["signed_z"]).eq(signed["_combined_sign"])
    signed["_success_weight"] = np.select(
        [signed["_combined_sign"].eq(0), signed["_concordant"]],
        [0.5 * signed["_unit_weight"], signed["_unit_weight"]],
        default=0.0,
    )
    grouped = signed.groupby("gene_symbol", as_index=False).agg(
        direction_concordant_source_units=("_success_weight", "sum"),
        direction_total_source_units=("_unit_weight", "sum"),
    )
    grouped[output_column] = np.where(
        grouped["direction_total_source_units"].gt(0),
        (1.0 + grouped["direction_concordant_source_units"]) / (2.0 + grouped["direction_total_source_units"]),
        0.5,
    )
    if output_column != "direction_confidence_index":
        return grouped[["gene_symbol", output_column]]
    grouped["direction_posterior_mean"] = grouped[output_column]
    return grouped[
        [
            "gene_symbol",
            output_column,
            "direction_concordant_source_units",
            "direction_total_source_units",
            "direction_posterior_mean",
        ]
    ]


def _random_effects_stouffer_layer(scores: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "gene_symbol",
        "re_stouffer_z",
        "re_stouffer_p",
        "re_stouffer_padj",
        "re_stouffer_shrinkage_factor",
    ]
    if scores.empty:
        return pd.DataFrame(columns=columns)
    out = scores[["gene_symbol"]].copy()
    z = pd.to_numeric(scores["stouffer_z"], errors="coerce").fillna(0.0)
    i2 = pd.to_numeric(scores.get("heterogeneity_i2", 0.0), errors="coerce").fillna(0.0).clip(0.0, 1.0)
    k_source = pd.to_numeric(scores.get("n_source_units", scores.get("n_studies", 1)), errors="coerce").fillna(1.0)
    factor = np.sqrt(1.0 + i2 * np.maximum(k_source - 1.0, 0.0))
    factor = pd.Series(factor, index=scores.index).replace(0.0, 1.0)
    out["re_stouffer_z"] = z / factor
    out["re_stouffer_p"] = 2.0 * norm.sf(np.abs(out["re_stouffer_z"].to_numpy(dtype=float)))
    out["re_stouffer_padj"] = bh_adjust(out["re_stouffer_p"].to_numpy(dtype=float))
    out["re_stouffer_shrinkage_factor"] = factor
    return out[columns]


def _rra_beta_layer(evidence: pd.DataFrame, *, total_source_units: int, min_studies: int) -> pd.DataFrame:
    columns = ["gene_symbol", "rra_rho", "rra_neglog10_rho", "rra_rank"]
    if evidence.empty or total_source_units <= 0:
        return pd.DataFrame(columns=columns)
    frame = evidence.copy()
    frame["normalized_rank"] = pd.to_numeric(frame["normalized_rank"], errors="coerce").clip(0.0, 1.0)
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "normalized_rank"])
    if frame.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    n_lists = max(int(total_source_units), 1)
    ln10 = float(np.log(10.0))
    rank_floor = float(np.finfo(float).tiny)
    for gene, group in frame.groupby("gene_symbol", sort=False):
        ranks = np.sort(group.drop_duplicates("source_unit_id")["normalized_rank"].to_numpy(dtype=float))
        if len(ranks) < min_studies:
            continue
        # Beta order-statistic RRA score in LOG space. The linear-space beta.cdf
        # underflows to 0.0 for genes ranked near-top across many lists, which
        # collapses the strongest genes into a single alphabetically-broken tie.
        # Tracking log(rho) keeps the top-of-list ordering and a usable magnitude.
        log_scores: list[float] = []
        for order_index, rank_value in enumerate(ranks, start=1):
            if order_index > n_lists:
                break
            log_scores.append(
                float(beta.logcdf(max(rank_value, rank_floor), order_index, n_lists - order_index + 1))
            )
        if not log_scores:
            continue
        rows.append({"gene_symbol": str(gene), "_log_rho": float(np.min(log_scores))})
    out = pd.DataFrame.from_records(rows, columns=["gene_symbol", "_log_rho"])
    if out.empty:
        return pd.DataFrame(columns=columns)
    out["_log_rho"] = pd.to_numeric(out["_log_rho"], errors="coerce")
    out = out.sort_values(["_log_rho", "gene_symbol"], ascending=[True, True]).reset_index(drop=True)
    out["rra_rank"] = np.arange(1, len(out) + 1, dtype=int)
    out["rra_rho"] = np.exp(out["_log_rho"].to_numpy(dtype=float)).clip(0.0, 1.0)
    out["rra_neglog10_rho"] = (-out["_log_rho"] / ln10).clip(lower=0.0)
    return out[columns]


def _effect_meta_layer(evidence: pd.DataFrame, *, min_studies: int) -> pd.DataFrame:
    columns = [
        "gene_symbol",
        "effect_meta_log2fc_re",
        "effect_meta_se",
        "effect_meta_ci_low",
        "effect_meta_ci_high",
        "effect_meta_tau2",
        "effect_meta_i2",
        "effect_meta_k",
        "effect_meta_se_source",
        "effect_meta_exact_weight_fraction",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=columns)
    frame = evidence.copy()
    for column in ["lfc", "signed_z"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "lfc", "signed_z"])
    frame = frame.loc[frame["signed_z"].abs().gt(0)].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame["_effect_se"] = (frame["lfc"].abs() / frame["signed_z"].abs()).replace([np.inf, -np.inf], np.nan)
    frame = frame.dropna(subset=["_effect_se"])
    frame = frame.loc[frame["_effect_se"].gt(0)].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, Any]] = []
    for gene, group in frame.groupby("gene_symbol", sort=False):
        group = group.drop_duplicates("source_unit_id")
        if len(group) < min_studies:
            continue
        y = group["lfc"].to_numpy(dtype=float)
        se = group["_effect_se"].to_numpy(dtype=float)
        vi = se**2
        valid = np.isfinite(y) & np.isfinite(vi) & (vi > 0)
        y = y[valid]
        vi = vi[valid]
        if len(y) < min_studies:
            continue
        w = np.divide(1.0, vi, out=np.zeros_like(vi, dtype=float), where=vi > 0)
        sum_w = float(np.sum(w))
        if not np.isfinite(sum_w) or sum_w <= EFFECT_META_MIN_WEIGHT_SUM:
            continue
        fixed = float(np.sum(w * y) / sum_w)
        q = float(np.sum(w * (y - fixed) ** 2))
        df = max(len(y) - 1, 0)
        c = float(sum_w - (np.sum(w**2) / sum_w)) if sum_w > 0 else 0.0
        tau2 = max(0.0, (q - df) / c) if c > 0 and df > 0 else 0.0
        w_re = np.divide(1.0, vi + tau2, out=np.zeros_like(vi, dtype=float), where=(vi + tau2) > 0)
        sum_w_re = float(np.sum(w_re))
        if not np.isfinite(sum_w_re) or sum_w_re <= EFFECT_META_MIN_WEIGHT_SUM:
            continue
        pooled = float(np.sum(w_re * y) / sum_w_re)
        pooled_se = float(np.sqrt(1.0 / sum_w_re))
        if not all(np.isfinite(value) for value in [pooled, pooled_se, tau2, q]):
            continue
        k_eff = int(len(y))
        # Hartung-Knapp-Sidik-Jonkman small-sample CI, truncated so it is never
        # narrower than the normal random-effects SE. Most genes have only k=2-3
        # source units, where the normal-approx z interval is anti-conservative.
        if k_eff >= 2 and sum_w_re > 0:
            q_hksj = float(np.sum(w_re * (y - pooled) ** 2) / ((k_eff - 1) * sum_w_re))
            se_ci = max(float(np.sqrt(q_hksj)), pooled_se)
            crit = float(t_dist.ppf(0.975, k_eff - 1))
        else:
            se_ci = pooled_se
            crit = 1.96
        ci_low = pooled - crit * se_ci
        ci_high = pooled + crit * se_ci
        i2 = max(0.0, (q - df) / q) if q > 0 and df > 0 else 0.0
        rows.append(
            {
                "gene_symbol": str(gene),
                "effect_meta_log2fc_re": pooled,
                "effect_meta_se": se_ci,
                "effect_meta_ci_low": ci_low,
                "effect_meta_ci_high": ci_high,
                "effect_meta_tau2": tau2,
                "effect_meta_i2": i2,
                "effect_meta_k": int(len(y)),
                "effect_meta_se_source": "derived_from_log2fc_and_two_sided_pvalue",
                "effect_meta_exact_weight_fraction": 0.0,
            }
        )
    return pd.DataFrame.from_records(rows, columns=columns)


def _priority_components_from_evidence(
    evidence: pd.DataFrame,
    *,
    support_denominator: float,
    use_reliability_weight: bool = False,
) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame()
    frame = evidence.copy()
    if use_reliability_weight and "source_reliability_weight" not in frame.columns:
        frame["source_reliability_weight"] = frame.get("source_recommended_weight", 1.0)
    for column in ["signed_z", "lfc", "normalized_rank", "weight", "source_reliability_weight"]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "signed_z", "normalized_rank", "weight"])
    if frame.empty:
        return pd.DataFrame()
    multiplier = frame["source_reliability_weight"].fillna(1.0) if use_reliability_weight else 1.0
    frame["_effective_weight"] = (frame["weight"] * multiplier).clip(lower=0.05)
    frame["_wz"] = frame["_effective_weight"] * frame["signed_z"]
    frame["_w2"] = frame["_effective_weight"] ** 2
    frame["_wlfc"] = frame["_effective_weight"] * frame["lfc"].fillna(0.0)
    frame["_w_lfc_denominator"] = np.where(frame["lfc"].notna(), frame["_effective_weight"], 0.0)
    eps = np.finfo(float).tiny
    frame["_weighted_log_rank"] = frame["_effective_weight"] * np.log(frame["normalized_rank"].clip(lower=eps, upper=1.0))
    grouped = frame.groupby("gene_symbol", as_index=False).agg(
        n_source_units=("source_unit_id", "nunique"),
        sum_wz=("_wz", "sum"),
        sum_w2=("_w2", "sum"),
        sum_wlfc=("_wlfc", "sum"),
        sum_w_lfc=("_w_lfc_denominator", "sum"),
        sum_weighted_log_rank=("_weighted_log_rank", "sum"),
        sum_effective_weight=("_effective_weight", "sum"),
    )
    grouped["stouffer_z"] = grouped["sum_wz"] / np.sqrt(grouped["sum_w2"])
    grouped["weighted_lfc"] = np.where(grouped["sum_w_lfc"].gt(0), grouped["sum_wlfc"] / grouped["sum_w_lfc"], np.nan)
    grouped["rank_product"] = np.exp(grouped["sum_weighted_log_rank"] / grouped["sum_effective_weight"])
    grouped["support_score"] = (np.log1p(grouped["n_source_units"]) / support_denominator).clip(0.0, 1.0)
    grouped["evidence_score"] = _component_strength_from_z(grouped["stouffer_z"])
    grouped["rank_score_component"] = (1.0 - grouped["rank_product"]).fillna(0.0).clip(0.0, 1.0)
    grouped["effect_score"] = _component_strength_from_lfc(grouped["weighted_lfc"])

    signs = frame[["gene_symbol", "signed_z", "_effective_weight"]].merge(
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
    signs["_direction_strength"] = signs["_effective_weight"] * signs["signed_z"].abs().clip(upper=8.0)
    signs["_concordant_strength"] = np.where(signs["_concordant"], signs["_direction_strength"], 0.0)
    direction = signs.groupby("gene_symbol", as_index=False).agg(
        total_strength=("_direction_strength", "sum"),
        concordant_strength=("_concordant_strength", "sum"),
    )
    direction["direction_score"] = np.where(
        direction["total_strength"].gt(0),
        direction["concordant_strength"] / direction["total_strength"],
        0.0,
    )
    grouped = grouped.merge(direction[["gene_symbol", "direction_score"]], on="gene_symbol", how="left")
    grouped["priority_score"] = _weighted_geometric_score_with_weights(grouped, PRIORITY_SCORE_WEIGHTS)
    return grouped[["gene_symbol", "priority_score"]]


def _leave_one_source_out_stability(
    evidence: pd.DataFrame,
    scores: pd.DataFrame,
    *,
    support_denominator: float,
) -> pd.DataFrame:
    columns = [
        "gene_symbol",
        "loo_median_rank",
        "loo_rank_iqr",
        "loo_rank_stability_score",
        "loo_top50_fraction",
        "loo_top100_fraction",
    ]
    if evidence.empty or scores.empty:
        return pd.DataFrame(columns=columns)
    evidence = evidence.copy()
    evidence["_loo_source_unit_id"] = evidence["source_unit_id"].astype(str)
    source_units = sorted(evidence["_loo_source_unit_id"].dropna().unique())
    if len(source_units) < 2:
        out = scores[["gene_symbol"]].copy()
        for column in columns[1:]:
            out[column] = 0.0
        return out

    total_genes = max(int(len(scores)), 1)
    penalty_rank = float(total_genes + 1)
    rank_records: dict[str, list[float]] = {str(gene): [] for gene in scores["gene_symbol"].astype(str)}
    for source_unit in source_units:
        subset = evidence.loc[evidence["_loo_source_unit_id"].ne(source_unit)].copy()
        components = _priority_components_from_evidence(
            subset,
            support_denominator=support_denominator,
            use_reliability_weight=False,
        )
        if components.empty:
            rank_map: dict[str, float] = {}
        else:
            components = components.loc[components["gene_symbol"].astype(str).isin(rank_records)].copy()
            ranked = components.sort_values(["priority_score", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
            rank_map = dict(zip(ranked["gene_symbol"].astype(str), (ranked.index + 1).astype(float), strict=False))
        for gene in rank_records:
            rank_records[gene].append(float(rank_map.get(gene, penalty_rank)))

    current_rank = scores.set_index("gene_symbol")["priority_rank"].astype(float).to_dict()
    rows: list[dict[str, float | str]] = []
    for gene, ranks in rank_records.items():
        values = np.array(ranks, dtype=float)
        median = float(np.median(values))
        q75, q25 = np.percentile(values, [75, 25])
        iqr = float(q75 - q25)
        shift = abs(median - float(current_rank.get(gene, penalty_rank)))
        stability = max(0.0, 1.0 - min(1.0, shift / max(float(total_genes), 1.0)))
        rows.append(
            {
                "gene_symbol": gene,
                "loo_median_rank": median,
                "loo_rank_iqr": iqr,
                "loo_rank_stability_score": stability,
                "loo_top50_fraction": float(np.mean(values <= 50.0)),
                "loo_top100_fraction": float(np.mean(values <= 100.0)),
            }
        )
    return pd.DataFrame.from_records(rows, columns=columns)


def _evidence_tier(
    top_percent: pd.Series,
    n_source_units: pd.Series,
    sign_concordance: pd.Series,
    *,
    total_source_units: int,
) -> pd.Series:
    """Assign an intentionally simple browsing tier from rank and support."""

    strong_support = max(1, min(3, int(total_source_units)))
    moderate_support = max(1, min(2, int(total_source_units)))
    tier_a = top_percent.le(1.0) & n_source_units.ge(strong_support) & sign_concordance.ge(0.90)
    tier_b = top_percent.le(5.0) & n_source_units.ge(moderate_support) & sign_concordance.ge(0.75)
    tier_c = top_percent.le(20.0) & n_source_units.ge(moderate_support)
    return pd.Series(np.select([tier_a, tier_b, tier_c], ["A", "B", "C"], default="D"), index=top_percent.index)


def _format_top_percent(value: float) -> str:
    if value < 0.01:
        return f"top {value:.4f}%"
    if value < 1.0:
        return f"top {value:.3f}%"
    return f"top {value:.2f}%"


def _format_percent(value: float) -> str:
    return f"{value:.1f}%"


def _heterogeneity_flags(i2: pd.Series) -> pd.Series:
    values = pd.to_numeric(i2, errors="coerce").fillna(0.0)
    return pd.Series(
        np.select(
            [values.ge(0.75), values.ge(0.50)],
            ["high_context_dependent_review", "moderate_context_review"],
            default="low_or_unestimated",
        ),
        index=i2.index,
        dtype="string",
    )


def degora_score_table(
    harmonized: pd.DataFrame,
    *,
    min_studies: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return gene scores, study-gene evidence rows, and score metadata."""

    score_harmonized, n_nonfinite_lfc_capped = _score_ready_harmonized(harmonized)
    evidence = study_gene_evidence(score_harmonized)
    source_quality_diagnostics = _source_quality_diagnostics_from_evidence(evidence)
    evidence = _attach_source_quality_diagnostics(evidence, source_quality_diagnostics)
    consensus = slice_consensus(score_harmonized, min_studies=min_studies)
    if consensus.empty:
        metadata = {
            "score_version": SCORE_VERSION,
            "score_formula": SCORE_FORMULA,
            "score_weights": SCORE_WEIGHTS,
            "primary_rank_column": PRIMARY_RANK_COLUMN,
            "primary_score_column": PRIMARY_SCORE_COLUMN,
            "primary_rank_interpretation": PRIMARY_RANK_DESCRIPTION,
            "priority_score_weights": PRIORITY_SCORE_WEIGHTS,
            "evidence_reliability_score_weights": RELIABILITY_SCORE_WEIGHTS,
            "min_studies": min_studies,
            "n_gene_scores": 0,
            "n_source_unit_gene_evidence_rows": int(len(evidence)),
            "n_source_units_total": int(evidence["source_unit_id"].nunique()) if not evidence.empty else 0,
            "n_contrasts_total": int(harmonized["study_id"].nunique()) if "study_id" in harmonized.columns else 0,
            "n_nonfinite_lfc_capped_for_score": n_nonfinite_lfc_capped,
            "independent_unit_for_consensus": "source_unit_id (paper_id when available, otherwise study_id)",
            "source_unit_collapse_rule": SOURCE_UNIT_COLLAPSE_RULE,
            "direction_concordance_rule": "evidence-strength-weighted concordance across independent source-unit representatives",
            "heterogeneity_rule": "source-unit z heterogeneity is reported as a Cochran-Q-style descriptive index over collapsed source-unit z values, with Q weighted by the predeclared sqrt-sample-size source weights (range 1-4, default 1) rather than by inverse variance; under homogeneity E[Q] therefore exceeds df, so heterogeneity_i2 = (Q-df)/Q is a positively biased weighted-dispersion index that is not comparable to a calibrated Higgins' I2. It is an audit/review-trigger field, not a calibrated random-effects model",
            "heterogeneity_flag_rule": "heterogeneity_i2 >= 0.75 is labeled high_context_dependent_review; >= 0.50 is moderate_context_review; flags are descriptive review aids, not score gates",
            "quality_weighted_score_formula": QUALITY_WEIGHTED_SCORE_FORMULA,
            "source_quality_weight_rules": {
                "source_input_type_weights": SOURCE_INPUT_TYPE_QUALITY_WEIGHTS,
                "table_scope_multipliers": TABLE_SCOPE_QUALITY_MULTIPLIERS,
                "replicate_multiplier": "1.0 if both groups have >=3 samples; 0.85 for >=2; 0.50 for >=1; 0.35 otherwise",
                "source_coherence_guardrail": "gold-panel-free source-source LFC Spearman check; low-quality sources with median pairwise Spearman < 0.05 receive source_coherence_weight=0.50 in the secondary score only",
                "source_reliability_shrinkage": "secondary-score weight shrunk toward neutral 0.65 using source gene coverage and pairwise-comparison evidence; not a calibrated probability",
            },
            "direction_confidence_rule": "Beta(1,1)-shrunk source-unit count concordance against the reported consensus signed-z direction: (1 + concordant source units) / (2 + observed source units); quality-weighted direction confidence uses reliability-weighted pseudo-counts against the quality-weighted consensus direction and is not a calibrated posterior probability",
            "random_effects_stouffer_rule": "descriptive heterogeneity-aware reporting lane only: stouffer_z / sqrt(1 + heterogeneity_i2 * (k_source_units - 1)). The divisor is a bounded ad-hoc shrinkage (capped at sqrt(2) for the dominant k=2 case) and inherits the descriptive heterogeneity_i2's small-k positive bias, so re_stouffer_p/padj are screening/triage fields, not calibrated formal random-effects inference or heterogeneity-corrected significance",
            "rra_rule": "parallel rank lane using beta order-statistic RobustRankAggreg-style rho over source-unit normalized ranks; missing source-unit lists are handled through the total source-unit universe; rho is computed in log space and rra_neglog10_rho (-log10 rho) preserves ordering for top genes whose rho underflows to 0; rho is not reported as a calibrated FDR",
            "effect_meta_rule": "parallel effect-size reporting lane only: DerSimonian-Laird random-effects inverse-variance meta-analysis of log2FC using SE derived from log2FC and two-sided p-value/signed_z when exact SE is unavailable; the 95% CI uses a truncated Hartung-Knapp-Sidik-Jonkman small-sample adjustment (t reference, df=k-1) so effect_meta_se/ci are not narrower than the normal random-effects interval; tau2/I2 are interpretable mainly for larger-k corpora",
            "loo_stability_rule": "leave-one-source-unit-out priority-rank stability; higher means the gene rank is less dependent on a single source unit",
            "source_quality_diagnostics": source_quality_diagnostics.to_dict(orient="records"),
        }
        return pd.DataFrame(columns=GENE_SCORE_COLUMNS), evidence, metadata

    # Support counts must come from the same eligible/collapsed evidence universe
    # that supplies the denominator, so the numerator can never exceed it. Counting
    # from raw rows let non-scoring rows (e.g. NaN signed_z dropped by
    # study_gene_evidence) inflate n_source_units above total_source_units.
    eligible_support = evidence.copy()
    eligible_support["gene_symbol"] = _string_column(eligible_support, "gene_symbol").str.upper().str.strip()
    eligible_support["source_unit_id"] = _string_column(eligible_support, "source_unit_id").str.strip()
    eligible_support["contributing_study_ids"] = _string_column(eligible_support, "contributing_study_ids")
    eligible_support = eligible_support.loc[
        eligible_support["gene_symbol"].ne("") & eligible_support["source_unit_id"].ne("")
    ].copy()

    def _count_contrasts(values: pd.Series) -> int:
        contrasts: set[str] = set()
        for entry in values.dropna().astype(str):
            for study_id in entry.split(";"):
                label = study_id.strip()
                if label:
                    contrasts.add(label)
        return len(contrasts)

    support = eligible_support.groupby("gene_symbol", as_index=False).agg(
        n_contrasts_observed=("contributing_study_ids", _count_contrasts),
        n_source_units=("source_unit_id", "nunique"),
        source_units=("source_unit_id", lambda values: ";".join(sorted(set(map(str, values))))),
    )
    total_source_units = int(evidence["source_unit_id"].nunique()) if not evidence.empty else 0
    denominator = np.log1p(total_source_units) if total_source_units > 1 else 1.0
    total_source_quality_weight = (
        float(source_quality_diagnostics["source_reliability_weight"].sum())
        if not source_quality_diagnostics.empty
        else 0.0
    )

    scores = consensus.merge(support, on="gene_symbol", how="left")
    scores["n_source_units"] = scores["n_source_units"].fillna(0).astype(int)
    scores["n_contrasts_observed"] = scores["n_contrasts_observed"].fillna(0).astype(int)
    scores["support_score"] = (np.log1p(scores["n_source_units"]) / denominator).clip(0.0, 1.0)
    scores["direction_score"] = pd.to_numeric(scores["sign_concordance"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    scores["evidence_score"] = _component_strength_from_z(scores["stouffer_z"])
    scores["rank_score_component"] = (1.0 - pd.to_numeric(scores["rank_product"], errors="coerce")).fillna(0.0).clip(0.0, 1.0)
    scores["effect_score"] = _component_strength_from_lfc(scores["weighted_lfc"])
    scores["degora_score"] = _weighted_geometric_score(scores)
    scores["priority_score"] = _weighted_geometric_score_with_weights(scores, PRIORITY_SCORE_WEIGHTS)
    quality_consensus = _quality_weighted_consensus(
        evidence.loc[evidence["gene_symbol"].isin(set(scores["gene_symbol"]))].copy(),
        total_source_quality_weight=total_source_quality_weight,
    )
    if quality_consensus.empty:
        for column in [
            "quality_weighted_degora_score",
            "quality_weighted_consensus_direction",
            "quality_weighted_sign_concordance",
            "source_quality_support_score",
            "source_quality_weight_sum",
            "quality_stouffer_z",
        ]:
            scores[column] = "" if column == "quality_weighted_consensus_direction" else 0.0
    else:
        scores = scores.merge(quality_consensus, on="gene_symbol", how="left")
        scores["quality_weighted_degora_score"] = pd.to_numeric(
            scores["quality_weighted_degora_score"], errors="coerce"
        ).fillna(0.0)
        scores["quality_weighted_consensus_direction"] = (
            scores["quality_weighted_consensus_direction"].astype("string").fillna("flat")
        )
        for column in [
            "quality_weighted_sign_concordance",
            "source_quality_support_score",
            "source_quality_weight_sum",
            "quality_stouffer_z",
        ]:
            scores[column] = pd.to_numeric(scores[column], errors="coerce").fillna(0.0)
    scores["consensus_direction"] = np.select(
        [scores["stouffer_z"].gt(0), scores["stouffer_z"].lt(0)],
        ["up", "down"],
        default="flat",
    )
    direction_confidence = _direction_confidence_from_evidence(
        evidence,
        scores[["gene_symbol", "stouffer_z"]],
        reference_z_column="stouffer_z",
        output_column="direction_confidence_index",
    )
    quality_direction_confidence = _direction_confidence_from_evidence(
        evidence,
        scores[["gene_symbol", "quality_stouffer_z"]],
        reference_z_column="quality_stouffer_z",
        use_reliability_weight=True,
        output_column="quality_weighted_direction_confidence_index",
    )
    scores = scores.merge(direction_confidence, on="gene_symbol", how="left").merge(
        quality_direction_confidence,
        on="gene_symbol",
        how="left",
    )
    scores["direction_confidence_index"] = pd.to_numeric(
        scores["direction_confidence_index"], errors="coerce"
    ).fillna(0.5)
    scores["quality_weighted_direction_confidence_index"] = pd.to_numeric(
        scores["quality_weighted_direction_confidence_index"], errors="coerce"
    ).fillna(0.5)
    for column in ["direction_concordant_source_units", "direction_total_source_units", "direction_posterior_mean"]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce").fillna(
            0.5 if column == "direction_posterior_mean" else 0.0
        )
    # Consensus direction must agree with the statistic that sign_concordance is
    # measured against. direction_score == sign_concordance counts each source
    # unit's signed_z against sign(stouffer_z) (see aggregate.stouffer_consensus),
    # so the reported direction is the combined signed-z direction, not the
    # weighted-LFC sign. Deriving it from weighted_lfc instead let the two
    # disagree for genes where the effect-size and significance-weighted z point
    # opposite ways, producing self-contradictory "<x>% down-concordant" labels.
    # weighted_lfc remains a separate reported column for the effect direction.
    scores = scores.merge(_random_effects_stouffer_layer(scores), on="gene_symbol", how="left")
    scores = scores.merge(
        _rra_beta_layer(evidence, total_source_units=total_source_units, min_studies=min_studies),
        on="gene_symbol",
        how="left",
    )
    scores = scores.merge(_effect_meta_layer(evidence, min_studies=min_studies), on="gene_symbol", how="left")
    high_confidence_min_units = max(1, min(2, total_source_units))
    scores["high_confidence"] = (
        scores["n_source_units"].ge(high_confidence_min_units)
        & scores["sign_concordance"].ge(0.75)
        & scores["rank_score_component"].ge(0.80)
        & scores["evidence_score"].ge(0.50)
    )
    priority_ranked = scores.sort_values(
        ["priority_score", "n_source_units", "direction_confidence_index", "gene_symbol"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    priority_rank_map = pd.Series(np.arange(1, len(priority_ranked) + 1), index=priority_ranked["gene_symbol"])
    scores["priority_rank"] = scores["gene_symbol"].map(priority_rank_map).fillna(0).astype(int)
    stability = _leave_one_source_out_stability(evidence, scores, support_denominator=denominator)
    scores = scores.merge(stability, on="gene_symbol", how="left")
    for column in ["loo_median_rank", "loo_rank_iqr", "loo_rank_stability_score", "loo_top50_fraction", "loo_top100_fraction"]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce").fillna(0.0)
    reliability_components = scores.copy()
    reliability_components["source_quality_support_score"] = pd.to_numeric(
        reliability_components["source_quality_support_score"], errors="coerce"
    ).fillna(0.0)
    scores["evidence_reliability_score"] = _weighted_geometric_score_with_weights(
        reliability_components,
        RELIABILITY_SCORE_WEIGHTS,
    )
    scores = scores.sort_values(
        ["degora_score", "n_source_units", "sign_concordance", "gene_symbol"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    quality_ranked = scores.sort_values(
        ["quality_weighted_degora_score", "n_source_units", "quality_weighted_sign_concordance", "gene_symbol"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    quality_rank_map = pd.Series(np.arange(1, len(quality_ranked) + 1), index=quality_ranked["gene_symbol"])
    scores["quality_weighted_degora_rank"] = scores["gene_symbol"].map(quality_rank_map).fillna(0).astype(int)
    scores["degora_rank"] = np.arange(1, len(scores) + 1)
    total_genes = max(int(len(scores)), 1)
    scores["rank_label"] = scores["quality_weighted_degora_rank"].map(lambda rank: f"#{int(rank):,} / {total_genes:,}")
    scores["top_percent"] = (100.0 * scores["degora_rank"] / total_genes).round(6)
    scores["priority_top_percent"] = (100.0 * scores["priority_rank"] / total_genes).round(6)
    scores["quality_weighted_top_percent"] = (100.0 * scores["quality_weighted_degora_rank"] / total_genes).round(6)
    scores["percentile"] = (100.0 * (1.0 - ((scores["degora_rank"] - 1.0) / total_genes))).round(6)
    # top_percent_label is the human-readable companion to the manuscript-primary
    # rank_label (quality_weighted_degora_rank), so derive it from the primary
    # quality_weighted_top_percent rather than the unweighted screening top_percent;
    # otherwise the label printed beside the primary rank reports the wrong percentile.
    scores["top_percent_label"] = scores["quality_weighted_top_percent"].map(_format_top_percent)
    scores["evidence_tier"] = _evidence_tier(
        scores["top_percent"],
        scores["n_source_units"],
        scores["sign_concordance"],
        total_source_units=total_source_units,
    )
    scores["support_label"] = scores["n_source_units"].map(lambda units: f"{int(units):,} / {total_source_units:,} source units")
    scores["direction_label"] = (
        (scores["sign_concordance"] * 100.0).map(_format_percent)
        + " "
        + scores["consensus_direction"].astype(str)
        + "-concordant"
    )
    scores["degora_score"] = scores["degora_score"].round(6)
    scores["priority_score"] = scores["priority_score"].round(6)
    scores["evidence_reliability_score"] = scores["evidence_reliability_score"].round(6)
    scores["direction_confidence_index"] = scores["direction_confidence_index"].round(6)
    scores["quality_weighted_direction_confidence_index"] = scores["quality_weighted_direction_confidence_index"].round(6)
    scores["direction_concordant_source_units"] = scores["direction_concordant_source_units"].round(3)
    scores["direction_total_source_units"] = scores["direction_total_source_units"].round(3)
    scores["direction_posterior_mean"] = scores["direction_posterior_mean"].round(6)
    scores["loo_median_rank"] = scores["loo_median_rank"].round(3)
    scores["loo_rank_iqr"] = scores["loo_rank_iqr"].round(3)
    scores["loo_rank_stability_score"] = scores["loo_rank_stability_score"].round(6)
    scores["loo_top50_fraction"] = scores["loo_top50_fraction"].round(6)
    scores["loo_top100_fraction"] = scores["loo_top100_fraction"].round(6)
    scores["quality_weighted_degora_score"] = scores["quality_weighted_degora_score"].round(6)
    scores["quality_weighted_sign_concordance"] = scores["quality_weighted_sign_concordance"].round(6)
    scores["source_quality_support_score"] = scores["source_quality_support_score"].round(6)
    scores["source_quality_weight_sum"] = scores["source_quality_weight_sum"].round(6)
    for column in ["heterogeneity_q", "heterogeneity_df", "heterogeneity_i2"]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce").fillna(0.0)
    scores["heterogeneity_q"] = scores["heterogeneity_q"].round(6)
    scores["heterogeneity_df"] = scores["heterogeneity_df"].astype(int)
    scores["heterogeneity_i2"] = scores["heterogeneity_i2"].round(6)
    scores["heterogeneity_flag"] = _heterogeneity_flags(scores["heterogeneity_i2"])
    for column in ["re_stouffer_z", "re_stouffer_p", "re_stouffer_padj", "re_stouffer_shrinkage_factor", "rra_rho"]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce")
    scores["re_stouffer_z"] = scores["re_stouffer_z"].fillna(0.0).round(6)
    scores["re_stouffer_p"] = scores["re_stouffer_p"].fillna(1.0).round(12)
    scores["re_stouffer_padj"] = scores["re_stouffer_padj"].fillna(1.0).round(12)
    scores["re_stouffer_shrinkage_factor"] = scores["re_stouffer_shrinkage_factor"].fillna(1.0).round(6)
    scores["rra_rho"] = scores["rra_rho"].fillna(1.0).round(12)
    scores["rra_neglog10_rho"] = pd.to_numeric(scores["rra_neglog10_rho"], errors="coerce").fillna(0.0).round(6)
    scores["rra_rank"] = pd.to_numeric(scores["rra_rank"], errors="coerce").fillna(0).astype(int)
    for column in [
        "effect_meta_log2fc_re",
        "effect_meta_se",
        "effect_meta_ci_low",
        "effect_meta_ci_high",
        "effect_meta_tau2",
        "effect_meta_i2",
        "effect_meta_exact_weight_fraction",
    ]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce").round(6)
    scores["effect_meta_k"] = pd.to_numeric(scores["effect_meta_k"], errors="coerce").fillna(0).astype(int)
    scores["effect_meta_se_source"] = scores["effect_meta_se_source"].fillna("")

    metadata = {
        "score_version": SCORE_VERSION,
        "score_formula": SCORE_FORMULA,
        "score_weights": SCORE_WEIGHTS,
        "primary_rank_column": PRIMARY_RANK_COLUMN,
        "primary_score_column": PRIMARY_SCORE_COLUMN,
        "primary_rank_interpretation": PRIMARY_RANK_DESCRIPTION,
        "priority_score_weights": PRIORITY_SCORE_WEIGHTS,
        "evidence_reliability_score_weights": RELIABILITY_SCORE_WEIGHTS,
        "quality_weighted_score_formula": QUALITY_WEIGHTED_SCORE_FORMULA,
        "quality_weighted_score_warning": "Quality-weighted ranking is the manuscript-facing default for benchmark tables; degora_score remains available as the unweighted/reference prioritization index.",
        "source_quality_weight_rules": {
            "source_input_type_weights": SOURCE_INPUT_TYPE_QUALITY_WEIGHTS,
            "table_scope_multipliers": TABLE_SCOPE_QUALITY_MULTIPLIERS,
            "replicate_multiplier": "1.0 if both groups have >=3 samples; 0.85 for >=2; 0.50 for >=1; 0.35 otherwise",
            "source_coherence_guardrail": "gold-panel-free source-source LFC Spearman check; low-quality sources with median pairwise Spearman < 0.05 receive source_coherence_weight=0.50 in the secondary score only",
            "source_reliability_shrinkage": "secondary-score weight shrunk toward neutral 0.65 using source gene coverage and pairwise-comparison evidence; not a calibrated probability",
        },
        "min_studies": min_studies,
        "n_gene_scores": int(len(scores)),
        "n_source_unit_gene_evidence_rows": int(len(evidence)),
        "n_source_units_total": total_source_units,
        "n_source_quality_outliers": int(source_quality_diagnostics["source_outlier_flag"].sum()) if not source_quality_diagnostics.empty else 0,
        "n_contrasts_total": int(harmonized["study_id"].nunique()) if "study_id" in harmonized.columns else 0,
        "n_nonfinite_lfc_capped_for_score": n_nonfinite_lfc_capped,
        "independent_unit_for_consensus": "source_unit_id (paper_id when available, otherwise study_id)",
        "source_unit_collapse_rule": SOURCE_UNIT_COLLAPSE_RULE,
        "direction_concordance_rule": "evidence-strength-weighted concordance across independent source-unit representatives",
        "heterogeneity_rule": "source-unit z heterogeneity is reported as a Cochran-Q-style descriptive index over collapsed source-unit z values; it is an audit field, not a calibrated random-effects model",
        "heterogeneity_flag_rule": "heterogeneity_i2 >= 0.75 is labeled high_context_dependent_review; >= 0.50 is moderate_context_review; flags are descriptive review aids, not score gates",
        "rank_interpretation": "degora_rank is the absolute rank among scored genes; top_percent is rank / total scored genes * 100, so smaller is more selective; percentile is 100 for the top-ranked gene and decreases with rank.",
        "high_confidence_rule": f"relative browsing flag: n_source_units >= min(2, total_source_units) = {high_confidence_min_units}, sign_concordance >= 0.75, rank_score_component >= 0.80, evidence_score >= 0.50; does not use stouffer_padj as a calibrated inferential gate",
        "direction_confidence_rule": "Beta(1,1)-shrunk source-unit count concordance against the reported consensus signed-z direction: (1 + concordant source units) / (2 + observed source units); quality-weighted direction confidence uses reliability-weighted pseudo-counts against the quality-weighted consensus direction and is not a calibrated posterior probability",
        "random_effects_stouffer_rule": "descriptive heterogeneity-aware reporting lane only: stouffer_z / sqrt(1 + heterogeneity_i2 * (k_source_units - 1)); p/padj are screening fields, not calibrated formal random-effects inference",
        "rra_rule": "parallel rank lane using beta order-statistic RobustRankAggreg-style rho over source-unit normalized ranks; missing source-unit lists are handled through the total source-unit universe; rho is computed in log space and rra_neglog10_rho (-log10 rho) preserves ordering for top genes whose rho underflows to 0; rho is not reported as a calibrated FDR",
        "effect_meta_rule": "parallel effect-size reporting lane only: DerSimonian-Laird random-effects inverse-variance meta-analysis of log2FC using SE derived from log2FC and two-sided p-value/signed_z when exact SE is unavailable; the 95% CI uses a truncated Hartung-Knapp-Sidik-Jonkman small-sample adjustment (t reference, df=k-1) so effect_meta_se/ci are not narrower than the normal random-effects interval; tau2/I2 are interpretable mainly for larger-k corpora",
        "loo_stability_rule": "leave-one-source-unit-out priority-rank stability; higher means the gene rank is less dependent on a single source unit",
        "evidence_tier_rules": {
            "A": "top_percent <= 1, n_source_units >= min(3, total_source_units), sign_concordance >= 0.90",
            "B": "top_percent <= 5, n_source_units >= min(2, total_source_units), sign_concordance >= 0.75",
            "C": "top_percent <= 20, n_source_units >= min(2, total_source_units)",
            "D": "lower-ranked or weakly supported",
        },
        "score_warning": "DEGORA score is for transparent prioritization and browsing, not a calibrated probability or a validation metric.",
        "source_quality_diagnostics": source_quality_diagnostics.to_dict(orient="records"),
    }
    return scores[GENE_SCORE_COLUMNS], evidence, metadata


def _active_study_table(catalog_path: Path | None, evidence: pd.DataFrame) -> pd.DataFrame:
    if catalog_path is not None and catalog_path.exists():
        catalog = read_catalog(catalog_path)
        active = catalog.loc[catalog_include_mask(catalog)].copy()
        active["source_unit_id"] = _source_unit_series(active)
        columns = [
            "study_id",
            "source_unit_id",
            "paper_id",
            "pipeline",
            "assay_type",
            "source_input_type",
            "platform",
            "normalization",
            "probe_collapse",
            "species",
            "cell_system",
            "hypoxia_modality",
            "duration_h",
            "n_ctrl",
            "n_treat",
            "source_path",
            "source_url",
            "notes",
        ]
        return active[[column for column in columns if column in active.columns]].sort_values("study_id").reset_index(drop=True)

    columns = [
        "study_id",
        "source_unit_id",
        "paper_id",
        "pipeline",
        "assay_type",
        "source_input_type",
        "platform",
        "normalization",
        "probe_collapse",
        "species",
        "cell_system",
        "hypoxia_modality",
        "duration_h",
        "n_ctrl",
        "n_treat",
        "source_path",
        "source_url",
    ]
    return evidence[columns].drop_duplicates("study_id").sort_values("study_id").reset_index(drop=True)


def _write_sqlite(
    db_path: Path,
    gene_scores: pd.DataFrame,
    evidence: pd.DataFrame,
    studies: pd.DataFrame,
    metadata: dict[str, Any],
) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    # Build the full database in a sibling temp file and atomically swap it in only after
    # every table and index succeeds. A failure mid-build therefore leaves the previous
    # good DB untouched instead of unlinking it up front and leaving an empty replacement.
    tmp_path = db_path.with_name(db_path.name + ".tmp")
    tmp_path.unlink(missing_ok=True)
    try:
        connection = sqlite3.connect(tmp_path)
        try:
            with connection:
                gene_scores.to_sql("genes", connection, index=False)
                evidence.to_sql("gene_evidence", connection, index=False)
                studies.to_sql("studies", connection, index=False)
                meta_rows = [{"key": key, "value": json.dumps(value, sort_keys=True) if not isinstance(value, str) else value} for key, value in metadata.items()]
                pd.DataFrame(meta_rows).to_sql("meta", connection, index=False)
                connection.execute("CREATE UNIQUE INDEX idx_genes_symbol ON genes(gene_symbol)")
                connection.execute("CREATE INDEX idx_genes_rank ON genes(degora_rank)")
                connection.execute("CREATE INDEX idx_genes_score ON genes(degora_score DESC)")
                if PRIMARY_RANK_COLUMN in gene_scores.columns:
                    connection.execute(f"CREATE INDEX idx_genes_primary_rank ON genes({PRIMARY_RANK_COLUMN})")
                if PRIMARY_SCORE_COLUMN in gene_scores.columns:
                    connection.execute(f"CREATE INDEX idx_genes_primary_score ON genes({PRIMARY_SCORE_COLUMN} DESC)")
                connection.execute("CREATE INDEX idx_evidence_gene ON gene_evidence(gene_symbol)")
                connection.execute("CREATE INDEX idx_evidence_study ON gene_evidence(study_id)")
                connection.execute("CREATE INDEX idx_studies_unit ON studies(source_unit_id)")
        finally:
            connection.close()
        tmp_path.replace(db_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def write_score_database(
    harmonized_path: Path,
    output_dir: Path,
    *,
    catalog_path: Path | None = None,
    db_path: Path | None = None,
    min_studies: int = 2,
    command: str | None = None,
    extra_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build score CSV, metadata JSON, and SQLite DB from a harmonized table."""

    output_dir.mkdir(parents=True, exist_ok=True)
    harmonized_path = harmonized_path.resolve()
    catalog_path = catalog_path.resolve() if catalog_path is not None else None
    db_path = (db_path or output_dir / "degora_scores.db").resolve()

    if harmonized_path.suffix.lower() in {".parquet", ".pq"}:
        harmonized = pd.read_parquet(harmonized_path)
    else:
        harmonized = pd.read_csv(harmonized_path, low_memory=False)
    gene_scores, evidence, metadata = degora_score_table(harmonized, min_studies=min_studies)
    gene_scores = primary_ranked_scores(gene_scores)
    version_info = runtime_version_info()
    metadata.update(
        {
            **version_info,
            "harmonized_path": str(harmonized_path),
            "catalog_path": str(catalog_path) if catalog_path else "",
            "db_path": str(db_path),
        }
    )
    if extra_metadata:
        metadata.update(extra_metadata)
    studies = _active_study_table(catalog_path, evidence)

    score_csv = output_dir / "degora_gene_scores.csv"
    metadata_json = output_dir / "degora_score_metadata.json"
    diagnostics_tsv = output_dir / "degora_source_quality_diagnostics.tsv"
    diagnostics_json = output_dir / "degora_source_quality_diagnostics.json"
    gene_scores.to_csv(score_csv, index=False)
    metadata_json.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    diagnostics = pd.DataFrame.from_records(
        metadata.get("source_quality_diagnostics", []),
        columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS,
    )
    diagnostics.to_csv(diagnostics_tsv, sep="\t", index=False)
    diagnostics_json.write_text(json.dumps(diagnostics.to_dict(orient="records"), indent=2, sort_keys=True) + "\n")
    _write_sqlite(db_path, gene_scores, evidence, studies, metadata)

    command = command or shell_command(
        [
            "make",
            "-C",
            "outputs/code",
            "score-db",
            f"HARMONIZED={_portable_cli_path(harmonized_path)}",
            f"OUTDIR={_portable_cli_path(output_dir)}",
            f"SCORE_DB={_portable_cli_path(db_path)}",
        ]
    )
    inputs: list[Path] = [harmonized_path]
    if catalog_path is not None:
        inputs.append(catalog_path)
    sidecar_metadata = {"generator": "degora-score-db", **version_info}
    for artifact in (score_csv, metadata_json, diagnostics_tsv, diagnostics_json, db_path):
        write_source_sidecar(artifact, command, inputs=inputs, metadata=sidecar_metadata)

    summary = {
        **version_info,
        "score_csv": str(score_csv.resolve()),
        "metadata_json": str(metadata_json.resolve()),
        "source_quality_diagnostics_tsv": str(diagnostics_tsv.resolve()),
        "source_quality_diagnostics_json": str(diagnostics_json.resolve()),
        "db_path": str(db_path),
        "n_gene_scores": int(len(gene_scores)),
        "n_evidence_rows": int(len(evidence)),
        "n_contrasts": int(len(studies)),
        "n_source_units": int(metadata["n_source_units_total"]),
        "n_source_quality_outliers": int(metadata.get("n_source_quality_outliers", 0)),
        "primary_rank_column": metadata.get("primary_rank_column", PRIMARY_RANK_COLUMN),
        "primary_score_column": metadata.get("primary_score_column", PRIMARY_SCORE_COLUMN),
        "top_genes": gene_scores.head(20)["gene_symbol"].tolist(),
    }
    summary_path = output_dir / "degora_score_db_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(summary_path, command, inputs=inputs, metadata=sidecar_metadata)
    return summary
