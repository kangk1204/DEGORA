"""DEGORA score tables and a local SQLite evidence database."""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Literal, Mapping, cast, overload

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.stats import beta, norm
from scipy.stats import t as t_dist

from . import SCORE_VERSION, runtime_version_info
from .aggregate import (
    SOURCE_UNIT_COLLAPSE_RULE,
    STOUFFER_WEIGHT_RULE,
    _collapse_preselected_source_unit_rows,
    _source_unit_series as _aggregate_source_unit_series,
    slice_consensus,
    source_unit_rows_for_aggregation,
    validate_min_studies,
    validate_normalized_rank,
)
from .formula_safety import (
    formula_guard_metadata,
    neutralize_formula_text,
    restore_formula_text_if_marked,
)
from .provenance import (
    apply_default_file_mode,
    artifact_output_lock,
    artifact_provenance_path,
    artifact_source_path,
    is_external_path_reference,
    output_directory_lock,
    portable_path,
    publication_target_lock_path,
    publish_staged_artifacts,
    sanitize_metadata,
    shell_command,
    source_sidecar_payloads,
)
from .slice_runner import catalog_include_mask, read_catalog
from .stats import bh_adjust


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
PRIMARY_SETTING_ID = "quality_weighted_primary"
PRIMARY_RANK_DESCRIPTION = (
    "quality_weighted_degora_rank is the primary browser and API rank. "
    "degora_rank and degora_score are retained as unweighted/reference outputs."
)
PRIORITY_SCORE_WEIGHTS = {
    "direction_score": 0.25,
    "evidence_score": 0.30,
    "rank_score_component": 0.30,
    "effect_score": 0.15,
}


def _json_safe_value(value: Any) -> Any:
    """Recursively replace non-finite numeric values with JSON null.

    Python's default encoder emits bare NaN/Infinity tokens, which are not valid
    RFC 8259 JSON. Score diagnostics legitimately have undefined pairwise metrics
    for a one-source corpus; those values are represented as ``None`` instead.
    """

    if isinstance(value, Mapping):
        return {str(key): _json_safe_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe_value(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe_value(value.item())
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _strict_json_dumps(value: Any, **kwargs: Any) -> str:
    """Serialize score output as standards-compliant JSON, failing closed."""

    return json.dumps(_json_safe_value(value), allow_nan=False, **kwargs)


def _diagnostic_records(diagnostics: pd.DataFrame) -> list[dict[str, Any]]:
    records = diagnostics.to_dict(orient="records")
    return cast(list[dict[str, Any]], _json_safe_value(records))


_AUTO_ABLATION_NAME = "\0auto"


@dataclass(frozen=True)
class ScoreAblation:
    """One defined ablation of the primary quality-weighted ranking.

    Component comparisons hold the number of contributing studies constant.
    That property is structural here rather than enforced by convention:

    * The gene universe and the eligible source-unit set are fixed by
      ``slice_consensus``, which no ablation touches. Every variant therefore
      scores exactly the same genes over exactly the same studies.
    * All five components of the primary lane are derived inside
      ``_quality_weighted_consensus`` from the per-source evidence frame, so an
      ablation is expressed as an override on that frame plus a weight vector,
      never as a different code path.

    ``component_weights`` may name any subset of :data:`SCORE_WEIGHTS`; the
    weighted geometric mean divides by the sum of the weights supplied, so a
    leave-one-out subset renormalises automatically and needs no rescaling by
    the caller.

    A default-constructed instance reproduces the shipped ranking exactly, so
    ``ScoreAblation()`` and ``None`` are interchangeable.
    """

    name: str = _AUTO_ABLATION_NAME
    component_weights: Mapping[str, float] | None = None
    disable_source_quality_weighting: bool = False
    disable_sample_size_weighting: bool = False
    notes: str = ""

    def __post_init__(self) -> None:
        explicit_name = self.name != _AUTO_ABLATION_NAME
        if explicit_name and (not isinstance(self.name, str) or not self.name.strip()):
            raise ValueError("ablation name must be a non-blank string")
        object.__setattr__(self, "name", self.name.strip() if explicit_name else "full")
        if self.component_weights is not None:
            weights = dict(self.component_weights)
            unknown = sorted(set(weights).difference(SCORE_WEIGHTS), key=str)
            if unknown:
                raise ValueError(
                    f"ablation {self.name!r} names unknown score components {unknown}; "
                    f"valid components are {sorted(SCORE_WEIGHTS)}"
                )
            if not weights:
                raise ValueError(f"ablation {self.name!r} must keep at least one score component")
            invalid: list[str] = []
            normalized: dict[str, float] = {}
            for component, raw_weight in weights.items():
                if isinstance(raw_weight, (bool, np.bool_)):
                    invalid.append(f"{component}={raw_weight!r}")
                    continue
                try:
                    weight = float(raw_weight)
                except (TypeError, ValueError):
                    invalid.append(f"{component}={raw_weight!r}")
                    continue
                if not np.isfinite(weight) or weight <= 0:
                    invalid.append(f"{component}={raw_weight!r}")
                    continue
                normalized[component] = weight
            if invalid:
                raise ValueError(
                    f"ablation {self.name!r} weights must be finite positive numbers; invalid: "
                    + ", ".join(invalid)
                )
            total = float(sum(normalized.values()))
            if not np.isfinite(total) or total <= 0:
                raise ValueError(f"ablation {self.name!r} has non-positive total weight {total!r}")
            # A frozen dataclass does not freeze a caller-owned dict. Snapshot the
            # mapping so post-construction mutation cannot bypass the validation above
            # and silently redefine an already-recorded ablation.
            object.__setattr__(self, "component_weights", MappingProxyType(normalized))
        if explicit_name and self.name == "full" and not self.is_default:
            raise ValueError(
                "ablation name 'full' is reserved for the canonical default score configuration"
            )
        if not explicit_name and not self.is_default:
            # Preserve the long-standing convenience of ScoreAblation(weights=...)
            # without mislabeling that custom configuration as canonical `full`.
            object.__setattr__(self, "name", "custom")

    @property
    def weights(self) -> dict[str, float]:
        return dict(self.component_weights) if self.component_weights is not None else dict(SCORE_WEIGHTS)

    @property
    def is_default(self) -> bool:
        return (
            self.weights == dict(SCORE_WEIGHTS)
            and not self.disable_source_quality_weighting
            and not self.disable_sample_size_weighting
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "component_weights": self.weights,
            "disable_source_quality_weighting": self.disable_source_quality_weighting,
            "disable_sample_size_weighting": self.disable_sample_size_weighting,
            "is_default": self.is_default,
            "notes": self.notes,
        }


# Columns overridden to a neutral 1.0 when source-quality weighting is ablated.
# All four are read downstream: _quality_weighted_consensus consumes the
# evidence frame (source_reliability_weight), while the reported diagnostics and
# the reliability lane consume the source-unit table, so both must be neutralised
# together or the ablation is only half applied.
_SOURCE_QUALITY_WEIGHT_COLUMNS = (
    "source_quality_weight",
    "source_coherence_weight",
    "source_recommended_weight",
    "source_reliability_weight",
)

# The five components of the primary quality-weighted lane, in SCORE_WEIGHTS order.
QUALITY_COMPONENT_COLUMNS = tuple(f"quality_{column}" for column in SCORE_WEIGHTS)
# The statistics the primary quality-weighted lane is computed from. The
# unweighted lane exports all three of its equivalents (stouffer_z, weighted_lfc,
# rank_product), so a reader can recompute degora_score from published output.
# The primary lane exported none of them: quality_weighted_consensus_direction is
# the sign of quality_stouffer_z, and a reader who saw a positive stouffer_z next
# to a "down" call had nothing in any output that could explain it.
QUALITY_STATISTIC_COLUMNS = ("quality_stouffer_z", "quality_weighted_lfc", "quality_rank_product")


def _portable_cli_path(path: Path, base_dir: Path) -> str:
    """Return a shareable path relative to an explicit output base."""

    return portable_path(path, base_dir)


def _portable_source_path_value(value: Any, source_base: Path, output_base: Path) -> Any:
    """Store source paths relative to the declared score-output path base."""

    if pd.isna(value):
        return value
    labels: list[str] = []
    for raw_label in str(value).split(";"):
        label = raw_label.strip()
        if label and not is_external_path_reference(label) and "://" not in label:
            source_path = Path(label)
            if not source_path.is_absolute():
                source_path = source_base / source_path
            label = portable_path(source_path, output_base)
        labels.append(label)
    return ";".join(labels)


def _portable_source_path_columns(
    frame: pd.DataFrame,
    source_base: Path,
    output_base: Path,
) -> pd.DataFrame:
    out = frame.copy()
    for column in ("source_path", "contributing_source_paths"):
        if column in out.columns:
            out[column] = _map_unique(
                out[column],
                lambda value: _portable_source_path_value(value, source_base, output_base),
            )
    return out

RELIABILITY_SCORE_WEIGHTS = {
    "support_score": 0.25,
    "source_quality_support_score": 0.25,
    "direction_confidence_index": 0.25,
    "loo_rank_stability_score": 0.25,
}
SCORE_FORMULA = (
    "100 * weighted_geometric_mean(support_score, direction_score, "
    "evidence_score, rank_score_component, effect_score); support is "
    "log-scaled by independent source units relative to the complete corpus "
    "(log1p(gene source units) / log1p(total source units), including 1.0 for "
    "a gene supported by the only source in a one-source corpus), direction is sign concordance, "
    "evidence is Stouffer-z strength, rank is one minus the 0-1 rank product "
    "(its complement, so higher is stronger), and "
    "effect is absolute weighted log2FC strength. The score is a transparent "
    "ranking aid, not a calibrated posterior probability. Consensus evidence "
    "is combined after aggregating related contrasts within each independent "
    "source unit without max-|z| representative selection."
)
SUPPORT_NORMALIZATION_RULE = (
    "support_score = log1p(gene source-unit count) / log1p(total corpus source-unit count), "
    "with a zero-source corpus left unscored; source_quality_support_score uses the same ratio "
    "over summed source-reliability weights. Therefore a gene supported by the only source in "
    "a one-source corpus has 1.0 in both support lanes."
)
QUALITY_WEIGHTED_SCORE_FORMULA = (
    "Fixed source-quality-weighted ranking: same component "
    "formula as the unweighted DEGORA score, but source-unit evidence is weighted "
    "by fixed heuristic source-quality features (source input type, table scope, and "
    "replicate count) together with a gold-panel-free source-coherence guardrail. "
    "It references no marker or gold-panel genes and is panel-agnostic by "
    "construction; it is reported as a relative prioritization index, not a "
    "calibrated probability."
)
REPLICATE_MULTIPLIER_RULE = (
    "1.0 if both groups have >=3 samples; 0.85 for >=2; 0.50 for >=1; "
    "0.35 if either count is zero; 0.75 if either count is missing or non-numeric"
)
HETEROGENEITY_RULE = (
    "source-unit z heterogeneity is reported as a Cochran-Q-style descriptive index over "
    "collapsed source-unit z values, with Q weighted by fixed sqrt-sample-size source weights "
    "(range 1-4, default 1) rather than inverse variance. Q is therefore not chi-square "
    "distributed under any null and its scale is arbitrary, so heterogeneity_i2 = (Q-df)/Q "
    "has no calibrated bias direction: raw values are frequently negative and are clamped to "
    "0, which makes the reported index effectively bimodal - near 0 for most genes, near 1 "
    "where the collapsed z values disagree. It is an audit/review-trigger field, not a "
    "calibrated Higgins' I2 and not a random-effects model; effect_meta_i2 is the "
    "inverse-variance-weighted estimate"
)
RANDOM_EFFECTS_STOUFFER_RULE = (
    "descriptive heterogeneity-aware reporting lane only: stouffer_z / sqrt(1 + "
    "heterogeneity_i2 * (k_source_units - 1)). The divisor is a bounded ad-hoc shrinkage "
    "based on the descriptive heterogeneity index; "
    "re_stouffer_p/padj are screening/triage fields, not calibrated formal random-effects "
    "inference or heterogeneity-corrected significance; finite tail probabilities are emitted "
    "without decimal-place rounding so small nonzero values are not printed as zero"
)
EFFECT_META_RULE = (
    "parallel effect-size reporting lane only: DerSimonian-Laird random-effects "
    "inverse-variance meta-analysis of log2FC using SE derived from log2FC and two-sided "
    "p-value/signed_z when exact SE is unavailable; the 95% CI uses a truncated "
    "Hartung-Knapp-Sidik-Jonkman small-sample adjustment (t reference, df=k-1) so "
    "effect_meta_se/ci are not narrower than the normal random-effects interval. For k=2, "
    "the t critical value is very large and the interval is often wide; tau2/I2 and interval "
    "estimates are descriptive and mainly interpretable for larger-k corpora"
)
LOO_STABILITY_RULE = (
    "leave-one-source-unit-out priority-rank stability over global corpus source-unit folds, "
    "using the same eligibility threshold and deterministic tie-break as the full priority "
    "lane; a fold in which a gene falls below min_studies is penalized as ineligible when at "
    "least one fold remains rank-evaluable. When zero folds are rank-evaluable, numeric LOO "
    "diagnostics are unavailable rather than zero"
)
EVIDENCE_RELIABILITY_VERSION = "evidence_reliability_v1_1_available_components"
EVIDENCE_RELIABILITY_RULE = (
    "conditional 0-100 weighted geometric mean over available diagnostics: support_score, "
    "source_quality_support_score, and direction_confidence_index are mandatory; "
    "loo_rank_stability_score participates only when at least one LOO fold is rank-evaluable. "
    "Available weights are renormalized row-wise, while an evaluated numeric zero remains "
    "negative evidence and is not treated as missing. This auxiliary summary is not used for "
    "the primary quality_weighted_degora_rank ordering"
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
    "evidence_reliability_components_used",
    "direction_confidence_index",
    "quality_weighted_direction_confidence_index",
    "direction_concordant_source_units",
    "direction_total_source_units",
    "direction_posterior_mean",
    "loo_total_folds",
    "loo_rank_evaluable_folds",
    "loo_penalty_folds",
    "loo_component_available",
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
    "stouffer_neglog10_p",
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
    "sign_convention",
    "platform",
    "normalization",
    "probe_collapse",
    "species",
    "cell_system",
    "hypoxia_modality",
    "condition",
    "duration_h",
    "time_course_mode",
    "temporal_mode",
    "n_ctrl",
    "n_treat",
    "lfc",
    "signed_z",
    "aggregate_pvalue",
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
    "source_direction_conflict_flag",
    "recommended_role",
]

# A source whose log2 fold changes run against every other source is, far more
# often than not, a contrast written the other way round - the one mistake the
# README says nothing downstream can catch. The coherence guardrail only ever
# down-weighted low-quality sources with a near-zero correlation, so a
# well-documented author table with its sign inverted kept full weight and was
# never mentioned. This flag is advisory: it changes no weight and no rank.
DIRECTION_CONFLICT_SPEARMAN = -0.10
DIRECTION_CONFLICT_ALPHA = 0.05
DIRECTION_CONFLICT_RULE = (
    "source_direction_conflict_flag is set when, for at least half of a source unit's pairwise "
    f"comparisons, the log2FC Spearman correlation is below {DIRECTION_CONFLICT_SPEARMAN:g} and "
    f"significantly negative (one-sided t approximation, p < {DIRECTION_CONFLICT_ALPHA:g}), AND the "
    f"unit's own median pairwise correlation is below {DIRECTION_CONFLICT_SPEARMAN:g}; two "
    "small tables of unrelated noise are therefore not flagged, and neither is a well-formed unit "
    "whose only conflicting comparison is the one against a reversed partner. It is an advisory "
    "review flag for a possibly reversed contrast direction and changes no weight or rank"
)

# These diagnostics deliberately sit outside every score component.  Their job
# is to surface two provenance mistakes that exact-path/hash checks cannot prove:
# a source table exported twice with harmless numeric reformatting, and two
# contrasts inside one declared source unit that appear to use opposite sign
# conventions.  Conservative thresholds keep the messages advisory; neither
# diagnostic changes source weights, eligibility, scores, or ranks.
NEAR_DUPLICATE_MIN_SHARED_GENES = 100
NEAR_DUPLICATE_MIN_SMALLER_COVERAGE = 0.90
NEAR_DUPLICATE_MIN_NEAR_IDENTICAL_FRACTION = 0.90
NEAR_DUPLICATE_LFC_ATOL = 1e-3
NEAR_DUPLICATE_SIGNED_Z_ATOL = 2e-2
NEAR_DUPLICATE_MAX_SOURCE_PAIRS = 500
NEAR_DUPLICATE_SOURCE_RULE = (
    f"advisory only: two declared source units are named when they share at least "
    f"{NEAR_DUPLICATE_MIN_SHARED_GENES} genes covering at least "
    f"{NEAR_DUPLICATE_MIN_SMALLER_COVERAGE:.0%} of the smaller source, and at least "
    f"{NEAR_DUPLICATE_MIN_NEAR_IDENTICAL_FRACTION:.0%} of both log2FC and signed-z values "
    f"are numerically near-identical; at most {NEAR_DUPLICATE_MAX_SOURCE_PAIRS} source pairs "
    "are checked per run, and the check changes no weight or rank"
)
WITHIN_SOURCE_DIRECTION_SPEARMAN = -0.80
WITHIN_SOURCE_DIRECTION_MAX_SIGN_AGREEMENT = 0.10
WITHIN_SOURCE_DIRECTION_MIN_SHARED_GENES = 20
WITHIN_SOURCE_DIRECTION_MAX_CONTRAST_PAIRS = 1_000
WITHIN_SOURCE_DIRECTION_RULE = (
    f"advisory only: two contrasts inside one source unit are named when at least "
    f"{WITHIN_SOURCE_DIRECTION_MIN_SHARED_GENES} genes overlap, log2FC Spearman is at most "
    f"{WITHIN_SOURCE_DIRECTION_SPEARMAN:g}, and same-sign agreement is at most "
    f"{WITHIN_SOURCE_DIRECTION_MAX_SIGN_AGREEMENT:.0%}; at most "
    f"{WITHIN_SOURCE_DIRECTION_MAX_CONTRAST_PAIRS} within-source contrast pairs are checked per run, "
    "and the check changes no value or rank"
)


def _negative_correlation_is_significant(rho: float, n_overlap: float) -> bool:
    """One-sided test that a Spearman rho is below zero, via the t approximation."""

    if not np.isfinite(rho) or rho >= 0 or n_overlap < 4:
        return False
    denominator = max(1.0 - rho * rho, 1e-12)
    t_value = rho * np.sqrt((n_overlap - 2) / denominator)
    return bool(t_dist.cdf(t_value, n_overlap - 2) < DIRECTION_CONFLICT_ALPHA)


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


def _validate_active_evidence_contract(harmonized: pd.DataFrame) -> None:
    """Reject rows that cannot satisfy the harmonized evidence contract.

    A finite signed z and normalized rank make a row eligible for the consensus
    lanes.  The evidence lane also needs the originating p-value and effect.  If
    either is missing, silently letting consensus count the source creates false
    replication and metadata that contradict the score. An inactive neutral audit
    row may omit the unscored counterpart field, but a present malformed effect or
    p-value is never treated as missing merely because the row will be dropped.
    Infinite effects remain the one documented exception:
    ``_score_ready_harmonized`` caps and reports them before either lane runs.
    """

    required = {"signed_z", "normalized_rank", "lfc", "pvalue"}
    if not required.issubset(harmonized.columns) or harmonized.empty:
        return
    raw_signed_z = harmonized["signed_z"]
    raw_lfc = harmonized["lfc"]
    raw_pvalue = harmonized["pvalue"]
    signed_z = pd.to_numeric(raw_signed_z, errors="coerce")
    normalized_rank = pd.to_numeric(harmonized["normalized_rank"], errors="coerce")
    lfc = pd.to_numeric(raw_lfc, errors="coerce")
    pvalue = pd.to_numeric(raw_pvalue, errors="coerce")
    signed_finite = np.isfinite(signed_z.to_numpy(dtype=float))
    rank_finite = np.isfinite(normalized_rank.to_numpy(dtype=float))

    mismatched_score_pair = signed_finite ^ rank_finite
    active = signed_finite & rank_finite
    signed_bool = np.fromiter(
        (isinstance(value, (bool, np.bool_)) for value in raw_signed_z.to_numpy(dtype=object)),
        dtype=bool,
        count=len(harmonized),
    )
    lfc_bool = np.fromiter(
        (isinstance(value, (bool, np.bool_)) for value in raw_lfc.to_numpy(dtype=object)),
        dtype=bool,
        count=len(harmonized),
    )
    pvalue_bool = np.fromiter(
        (isinstance(value, (bool, np.bool_)) for value in raw_pvalue.to_numpy(dtype=object)),
        dtype=bool,
        count=len(harmonized),
    )
    pvalue_values = pvalue.to_numpy(dtype=float)
    lfc_missing = raw_lfc.isna().to_numpy(dtype=bool)
    pvalue_missing = raw_pvalue.isna().to_numpy(dtype=bool)
    # Missing evidence is permitted only for an inactive neutral audit row.
    # A present-but-malformed value is never equivalent to missing and must fail
    # even if the other field proves that the row is neutral.
    malformed_present_evidence = (
        (~lfc_missing & (lfc.isna().to_numpy(dtype=bool) | lfc_bool))
        | (
            ~pvalue_missing
            & (
                ~np.isfinite(pvalue_values)
                | (pvalue_values < 0.0)
                | (pvalue_values > 1.0)
                | pvalue_bool
            )
        )
    )
    incomplete_active_evidence = active & (lfc_missing | pvalue_missing)
    invalid = (
        mismatched_score_pair
        | (active & signed_bool)
        | malformed_present_evidence
        | incomplete_active_evidence
    )
    if invalid.any():
        raise ValueError(
            "scoring evidence contract failed: active rows require a finite, non-boolean "
            "signed_z/normalized_rank pair and non-missing lfc/pvalue; every present lfc "
            "must be numeric and non-boolean, and every present pvalue must be a "
            "non-boolean finite value in [0, 1]; "
            f"inconsistent value(s) in {int(invalid.sum())} row(s)"
        )


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


def _map_unique(values: pd.Series, func: Callable[[Any], Any]) -> pd.Series:
    """Apply ``func`` once per distinct value rather than once per row.

    These frames repeat a handful of labels - source paths, pipeline names -
    across tens of thousands of rows, and pandas' arrow-backed string columns
    make per-element access expensive. Mapping the distinct values collapsed
    hundreds of thousands of calls into a few dozen.
    """

    if values.empty:
        return values
    codes, uniques = pd.factorize(values, use_na_sentinel=False)
    mapped: NDArray[np.object_] = np.empty(len(uniques), dtype=object)
    for index, unique_value in enumerate(uniques):
        mapped[index] = func(unique_value)
    return pd.Series(mapped[codes], index=values.index)


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


def _quality_label_frame(values: pd.Series) -> pd.Series:
    """Vectorized equivalent of ``_quality_label`` for evidence tables."""

    numeric = pd.to_numeric(values, errors="coerce")
    finite = pd.Series(np.isfinite(numeric.to_numpy(dtype=float)), index=values.index)
    labels = np.select(
        [finite & numeric.ge(0.85), finite & numeric.ge(0.60), finite],
        ["high", "medium", "low"],
        default="unknown",
    )
    return pd.Series(labels, index=values.index, dtype="string")


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
    prior (0.65, a deliberately conservative midpoint of the fixed heuristic
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
    source_labels = source_type.astype("string").fillna("")
    scope_labels = table_scope.astype("string").fillna("")
    source_mapping = {label: _source_input_type_weight(label) for label in source_labels.unique()}
    scope_mapping = {label: _table_scope_multiplier(label) for label in scope_labels.unique()}
    source_weights = source_labels.map(source_mapping).astype(float)
    scope_weights = scope_labels.map(scope_mapping).astype(float)

    ctrl = pd.to_numeric(n_ctrl, errors="coerce").to_numpy(dtype=float)
    treat = pd.to_numeric(n_treat, errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(ctrl) & np.isfinite(treat)
    minimum = np.minimum(ctrl, treat)
    replicate_weights = np.select(
        [~valid, minimum >= 3.0, minimum >= 2.0, minimum >= 1.0],
        [0.75, 1.00, 0.85, 0.50],
        default=0.35,
    )
    weights = source_weights.to_numpy(dtype=float) * scope_weights.to_numpy(dtype=float) * replicate_weights
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
        # The catalog's generic `condition` column is stored under a topic-specific
        # name that predates the tool being topic-neutral. Emit both so the API,
        # the database schema and the shipped workbook can carry the neutral name
        # without breaking a reader that already depends on the old one.
        ("condition", "hypoxia_modality"),
        ("table_scope", "table_scope"),
        ("sign_convention", "sign_convention"),
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
    # Coerce once instead of running a Python-level to_numeric per group.
    for count_column in ("n_ctrl", "n_treat"):
        if count_column in sorted_frame.columns:
            sorted_frame[count_column] = pd.to_numeric(sorted_frame[count_column], errors="coerce")
    base = sorted_frame.groupby(group_cols, as_index=False, sort=False).agg(
        n_ctrl=("n_ctrl", "min"),
        n_treat=("n_treat", "min"),
        n_genes_in_study=("n_genes_in_study", "max"),
        min_source_pvalue=("pvalue", "min"),
        min_source_padj=("padj", "min"),
    )

    source_unit_meta = sorted_frame[["source_unit_id"]].drop_duplicates().copy()
    for target, source in source_unit_string_specs:
        values = sorted_frame[["source_unit_id", source]].copy()
        values[source] = _map_unique(values[source], _clean_join_label)
        values = values.loc[values[source].ne("")].drop_duplicates().sort_values(["source_unit_id", source])
        joined = values.groupby("source_unit_id", sort=False)[source].agg(";".join).rename(target).reset_index()
        source_unit_meta = source_unit_meta.merge(joined, on="source_unit_id", how="left")

    out = base.merge(source_unit_meta, on="source_unit_id", how="left")
    for target, source in gene_source_string_specs:
        values = sorted_frame[[*group_cols, source]].copy()
        values[source] = _map_unique(values[source], _clean_join_label)
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
    frame["normalized_rank"] = validate_normalized_rank(
        frame,
        context="study-gene evidence scoring",
    )
    frame["n_ctrl"] = _as_numeric(frame, "n_ctrl")
    frame["n_treat"] = _as_numeric(frame, "n_treat")
    frame["n_genes_in_study"] = _as_numeric(frame, "n_genes_in_study")
    # Per-contrast source weights are derived downstream by
    # aggregate.source_unit_rows_for_aggregation, which applies the documented
    # MAX_SOURCE_SAMPLE_WEIGHT cap. Deriving them a second time here would leave an
    # uncapped copy of the rule in the tree for a later reader to follow.
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
        ("sign_convention", ""),
        ("source_path", ""),
        ("source_url", ""),
    ]:
        if column not in frame.columns:
            frame[column] = fallback

    # Select early/late/peak rows before dropping unusable score values. Otherwise
    # an unusable documented early time point is removed and a later point silently
    # takes its place. Metadata and collapsed scores both derive from this one exact
    # selected frame.
    selected_frame = source_unit_rows_for_aggregation(frame)
    pvalue = pd.to_numeric(selected_frame["pvalue"], errors="coerce")
    lfc = pd.to_numeric(selected_frame["lfc"], errors="coerce")
    eligible = np.isfinite(pvalue.to_numpy(dtype=float)) & np.isfinite(lfc.to_numpy(dtype=float))
    selected_frame = selected_frame.loc[eligible].copy()
    collapsed = _collapse_preselected_source_unit_rows(selected_frame)
    meta = _metadata_for_study_gene_units(selected_frame)
    out = collapsed.merge(meta.drop(columns=["n_genes_in_study"], errors="ignore"), on=["gene_symbol", "source_unit_id"], how="left")
    out["aggregate_pvalue"] = 2.0 * norm.sf(np.abs(pd.to_numeric(out["signed_z"], errors="coerce")))
    out["source_quality_weight"] = _source_quality_weight_frame(out)
    out["source_quality_label"] = _quality_label_frame(out["source_quality_weight"])
    out["source_coherence_weight"] = 1.0
    out["source_recommended_weight"] = out["source_quality_weight"]
    out["source_reliability_weight"] = out["source_quality_weight"]
    out["source_reliability_label"] = _quality_label_frame(out["source_reliability_weight"])
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
    total_weight = float(sum(weights.values()))
    missing = [column for column in weights if column not in frame.columns]
    if missing:
        raise ValueError(f"missing required score component column(s): {', '.join(missing)}")
    score: NDArray[np.float64] = np.ones(len(frame), dtype=float)
    for column, weight in weights.items():
        component = pd.to_numeric(frame[column], errors="coerce")
        bad = ~np.isfinite(component.to_numpy(dtype=float))
        if bad.any():
            raise ValueError(
                f"non-finite score component {column!r} in {int(bad.sum())} row(s); "
                "required score components must be finite before aggregation"
            )
        outside = component.lt(0.0) | component.gt(1.0)
        if outside.any():
            raise ValueError(
                f"score component {column!r} must satisfy 0 <= component <= 1; "
                f"out-of-range value(s) in {int(outside.sum())} row(s)"
            )
        # Components are constructed on [0, 1]. Preserve an evaluated zero as
        # actual negative evidence: in a geometric mean it is absorbing. The old
        # 1e-6 floor made a documented zero produce a positive score.
        score *= np.power(component.to_numpy(dtype=float), float(weight) / total_weight)
    return pd.Series(100.0 * score, index=frame.index)


def _weighted_geometric_score(frame: pd.DataFrame) -> pd.Series:
    return _weighted_geometric_score_with_weights(frame, SCORE_WEIGHTS)


def _evidence_reliability_score(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Score available reliability diagnostics without conflating N/A with zero."""

    availability_column = "loo_component_available"
    if availability_column not in frame.columns:
        raise ValueError(f"missing required reliability availability column: {availability_column}")
    availability = frame[availability_column]
    if availability.isna().any():
        raise ValueError("loo_component_available must be defined for every scored gene")
    availability = availability.astype(bool)

    loo_column = "loo_rank_stability_score"
    if loo_column not in frame.columns:
        raise ValueError(f"missing required reliability component column: {loo_column}")
    loo = pd.to_numeric(frame[loo_column], errors="coerce")
    finite_loo = pd.Series(np.isfinite(loo.to_numpy(dtype=float)), index=frame.index)
    if (availability & ~finite_loo).any():
        raise ValueError("LOO is marked available but loo_rank_stability_score is non-finite")
    if ((~availability) & finite_loo).any():
        raise ValueError("LOO is marked unavailable but loo_rank_stability_score is numeric")

    mandatory_weights = {
        column: weight
        for column, weight in RELIABILITY_SCORE_WEIGHTS.items()
        if column != "loo_rank_stability_score"
    }
    result = pd.Series(index=frame.index, dtype=float)
    components_used = pd.Series(3, index=frame.index, dtype=int)
    if (~availability).any():
        result.loc[~availability] = _weighted_geometric_score_with_weights(
            frame.loc[~availability],
            mandatory_weights,
        )
    if availability.any():
        result.loc[availability] = _weighted_geometric_score_with_weights(
            frame.loc[availability],
            RELIABILITY_SCORE_WEIGHTS,
        )
        components_used.loc[availability] = 4
    return result, components_used


def _weighted_geometric_score_from_components(
    frame: pd.DataFrame,
    prefix: str,
    weights: Mapping[str, float] | None = None,
) -> pd.Series:
    active = dict(weights) if weights is not None else dict(SCORE_WEIGHTS)
    components = pd.DataFrame(index=frame.index)
    for column in active:
        quality_column = f"{prefix}{column}"
        if quality_column not in frame.columns:
            raise ValueError(f"missing required score component column: {quality_column}")
        components[column] = frame[quality_column]
    return _weighted_geometric_score_with_weights(components, active)


def _stouffer_inference_warning(evidence: pd.DataFrame) -> str:
    table_scope = (
        evidence["table_scope"].astype("string").fillna("").str.strip().str.lower()
        if "table_scope" in evidence.columns
        else pd.Series(dtype="string")
    )
    if table_scope.eq("deg_only").any():
        return (
            "At least one active source is a DEG-only/significant-gene-only table. Its rows "
            "were selected conditional on source significance, so stouffer_p, stouffer_padj, "
            "re_stouffer_p, and re_stouffer_padj are descriptive screening fields and are not "
            "valid inferential p-values or FDR estimates."
        )
    return (
        "Stouffer and random-effects-Stouffer p-value fields are descriptive screening outputs; "
        "their calibration depends on the completeness and selection process of every input table."
    )


@overload
def _source_quality_diagnostics_from_evidence(
    evidence: pd.DataFrame,
    *,
    return_pairwise: Literal[False] = False,
) -> pd.DataFrame: ...


@overload
def _source_quality_diagnostics_from_evidence(
    evidence: pd.DataFrame,
    *,
    return_pairwise: Literal[True],
) -> tuple[pd.DataFrame, pd.DataFrame]: ...


def _source_quality_diagnostics_from_evidence(
    evidence: pd.DataFrame,
    *,
    return_pairwise: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    pairwise_columns = ["source_a", "source_b", "lfc_spearman", "sign_agreement", "overlap"]
    if evidence.empty:
        diagnostics = pd.DataFrame(columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS)
        if return_pairwise:
            return diagnostics, pd.DataFrame(columns=pairwise_columns)
        return diagnostics

    frame = evidence.copy()
    frame["source_unit_id"] = _string_column(frame, "source_unit_id").str.strip()
    frame["gene_symbol"] = _string_column(frame, "gene_symbol").str.upper().str.strip()
    frame["lfc"] = _as_numeric(frame, "lfc")
    frame["source_quality_weight"] = _as_numeric(frame, "source_quality_weight", default=0.65).fillna(0.65)
    frame = frame.dropna(subset=["source_unit_id", "gene_symbol", "lfc"])
    frame = frame.loc[frame["source_unit_id"].ne("") & frame["gene_symbol"].ne("")].copy()

    if frame.empty:
        diagnostics = pd.DataFrame(columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS)
        if return_pairwise:
            return diagnostics, pd.DataFrame(columns=pairwise_columns)
        return diagnostics

    rows: list[dict[str, Any]] = []
    pairwise_rows: list[dict[str, Any]] = []
    wide_lfc = frame.pivot_table(index="gene_symbol", columns="source_unit_id", values="lfc", aggfunc="mean")
    sign_frame = frame.assign(_sign=np.sign(frame["lfc"]))
    wide_sign = sign_frame.pivot_table(index="gene_symbol", columns="source_unit_id", values="_sign", aggfunc="mean")
    source_units = sorted(frame["source_unit_id"].dropna().astype(str).unique())
    pairwise: dict[str, list[dict[str, float]]] = {source_unit: [] for source_unit in source_units}
    for index, source_a in enumerate(source_units):
        for source_b in source_units[index + 1 :]:
            overlap = wide_lfc[[source_a, source_b]].dropna()
            sign_overlap = wide_sign[[source_a, source_b]].dropna()
            if len(overlap) >= 3 and overlap[source_a].nunique() > 1 and overlap[source_b].nunique() > 1:
                lfc_spearman = float(overlap[source_a].corr(overlap[source_b], method="spearman"))
            else:
                # Spearman is undefined (and scipy warns) when either side is constant.
                lfc_spearman = np.nan
            if len(sign_overlap) > 0:
                sign_agreement = float((sign_overlap[source_a] * sign_overlap[source_b] > 0).mean())
            else:
                sign_agreement = np.nan
            pairwise_rows.append(
                {
                    "source_a": source_a,
                    "source_b": source_b,
                    "lfc_spearman": lfc_spearman,
                    "sign_agreement": sign_agreement,
                    "overlap": float(len(overlap)),
                }
            )
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
        conflicting_pairs = [
            item
            for item in comparisons
            if np.isfinite(item["lfc_spearman"])
            and item["lfc_spearman"] < DIRECTION_CONFLICT_SPEARMAN
            and _negative_correlation_is_significant(float(item["lfc_spearman"]), float(item["overlap"]))
        ]
        # "At least half the comparisons conflict" is satisfied by a single
        # conflicting partner whenever a unit has two comparisons, so in a
        # three-source corpus with one reversed contrast every unit qualified -
        # the two well-formed ones with a median Spearman of 0.00, printed inside
        # a message asserting they disagreed with the rest. The unit's own median
        # must clear the same threshold before it is named.
        direction_conflict_flag = bool(
            lfc_corrs
            and len(conflicting_pairs) * 2 >= len(lfc_corrs)
            and np.isfinite(median_spearman)
            and median_spearman < DIRECTION_CONFLICT_SPEARMAN
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
                "source_direction_conflict_flag": direction_conflict_flag,
                "recommended_role": "sensitivity" if outlier_flag or source_quality < 0.60 else "primary",
            }
        )

    diagnostics = pd.DataFrame.from_records(rows, columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS)
    if return_pairwise:
        pairwise_frame = pd.DataFrame.from_records(pairwise_rows, columns=pairwise_columns)
        return diagnostics, pairwise_frame
    return diagnostics


def near_duplicate_source_unit_warnings(evidence: pd.DataFrame) -> list[str]:
    """Return advisory warnings for almost identical declared source units.

    Exact path and byte-hash checks run earlier in the pipeline.  This second
    check operates on harmonized values so a table that was merely re-saved or
    rounded cannot silently masquerade as independent replication.  Because two
    genuinely independent experiments can still agree unusually closely, the
    result is never used as a score or eligibility gate.
    """

    required = {"source_unit_id", "gene_symbol", "lfc", "signed_z"}
    if evidence.empty or not required.issubset(evidence.columns):
        return []
    frame = evidence[["source_unit_id", "gene_symbol", "lfc", "signed_z"]].copy()
    frame["source_unit_id"] = _string_column(frame, "source_unit_id").str.strip()
    frame["gene_symbol"] = _string_column(frame, "gene_symbol").str.upper().str.strip()
    frame["lfc"] = _as_numeric(frame, "lfc")
    frame["signed_z"] = _as_numeric(frame, "signed_z")
    frame = frame.loc[frame["source_unit_id"].ne("") & frame["gene_symbol"].ne("")].copy()
    frame = frame.dropna(subset=["lfc", "signed_z"])
    if frame.empty:
        return []
    frame = frame.groupby(["source_unit_id", "gene_symbol"], as_index=False)[["lfc", "signed_z"]].mean()
    source_frames = {
        str(source_unit): group.set_index("gene_symbol")[["lfc", "signed_z"]].copy()
        for source_unit, group in frame.groupby("source_unit_id", sort=True)
    }
    source_units = sorted(source_frames)
    n_pairs = len(source_units) * (len(source_units) - 1) // 2
    if n_pairs > NEAR_DUPLICATE_MAX_SOURCE_PAIRS:
        return [
            f"near-duplicate source-unit advisory was skipped: {len(source_units):,} source units "
            f"would require {n_pairs:,} pair checks, above the bounded limit of "
            f"{NEAR_DUPLICATE_MAX_SOURCE_PAIRS:,}. Review source-unit provenance manually; weights "
            "and ranks are unchanged."
        ]
    warnings: list[str] = []
    for index, source_a in enumerate(source_units):
        left = source_frames[source_a]
        for source_b in source_units[index + 1 :]:
            right = source_frames[source_b]
            shared = left.index.intersection(right.index, sort=False)
            n_overlap = int(len(shared))
            smaller_n = min(int(len(left)), int(len(right)))
            if smaller_n <= 0 or n_overlap < NEAR_DUPLICATE_MIN_SHARED_GENES:
                continue
            smaller_coverage = n_overlap / smaller_n
            if smaller_coverage < NEAR_DUPLICATE_MIN_SMALLER_COVERAGE:
                continue
            left_overlap = left.loc[shared]
            right_overlap = right.loc[shared]
            lfc_close = np.isclose(
                left_overlap["lfc"].to_numpy(dtype=float),
                right_overlap["lfc"].to_numpy(dtype=float),
                rtol=1e-3,
                atol=NEAR_DUPLICATE_LFC_ATOL,
            )
            signed_z_close = np.isclose(
                left_overlap["signed_z"].to_numpy(dtype=float),
                right_overlap["signed_z"].to_numpy(dtype=float),
                rtol=1e-3,
                atol=NEAR_DUPLICATE_SIGNED_Z_ATOL,
            )
            lfc_fraction = float(lfc_close.mean())
            signed_z_fraction = float(signed_z_close.mean())
            if (
                lfc_fraction < NEAR_DUPLICATE_MIN_NEAR_IDENTICAL_FRACTION
                or signed_z_fraction < NEAR_DUPLICATE_MIN_NEAR_IDENTICAL_FRACTION
            ):
                continue
            warnings.append(
                f"source units {source_a!r} and {source_b!r} may be duplicate exports: "
                f"{n_overlap:,} shared genes cover {smaller_coverage:.1%} of the smaller source, "
                f"with {lfc_fraction:.1%} near-identical log2FC and {signed_z_fraction:.1%} "
                "near-identical signed-z values. Verify that they are independent publications or "
                "datasets before interpreting replication; weights and ranks are unchanged."
            )
    return warnings


def within_source_direction_warnings(harmonized: pd.DataFrame) -> list[str]:
    """Return advisory warnings for strongly reversed contrasts in one source unit."""

    required = {"study_id", "gene_symbol", "lfc"}
    if harmonized.empty or not required.issubset(harmonized.columns):
        return []
    frame = harmonized.copy()
    frame["source_unit_id"] = _source_unit_series(frame)
    frame["study_id"] = _string_column(frame, "study_id").str.strip()
    frame["gene_symbol"] = _string_column(frame, "gene_symbol").str.upper().str.strip()
    frame["lfc"] = _as_numeric(frame, "lfc")
    frame = frame.loc[
        frame["source_unit_id"].ne("")
        & frame["study_id"].ne("")
        & frame["gene_symbol"].ne("")
    ].dropna(subset=["lfc"])
    if frame.empty:
        return []
    frame = frame.groupby(["source_unit_id", "study_id", "gene_symbol"], as_index=False)["lfc"].mean()
    source_groups = list(frame.groupby("source_unit_id", sort=True))
    n_pairs = sum(
        int(source_frame["study_id"].nunique())
        * (int(source_frame["study_id"].nunique()) - 1)
        // 2
        for _, source_frame in source_groups
    )
    if n_pairs > WITHIN_SOURCE_DIRECTION_MAX_CONTRAST_PAIRS:
        return [
            f"within-source direction advisory was skipped: the declared source units would require "
            f"{n_pairs:,} contrast-pair checks, above the bounded limit of "
            f"{WITHIN_SOURCE_DIRECTION_MAX_CONTRAST_PAIRS:,}. Review contrast directions manually; "
            "values and ranks are unchanged."
        ]
    warnings: list[str] = []
    for source_unit, source_frame in source_groups:
        contrast_frames = {
            str(study_id): group.set_index("gene_symbol")[["lfc"]].copy()
            for study_id, group in source_frame.groupby("study_id", sort=True)
        }
        contrasts = sorted(contrast_frames)
        for index, contrast_a in enumerate(contrasts):
            left = contrast_frames[contrast_a]
            for contrast_b in contrasts[index + 1 :]:
                right = contrast_frames[contrast_b]
                shared = left.index.intersection(right.index, sort=False)
                if len(shared) < WITHIN_SOURCE_DIRECTION_MIN_SHARED_GENES:
                    continue
                lfc_a = left.loc[shared, "lfc"].to_numpy(dtype=float)
                lfc_b = right.loc[shared, "lfc"].to_numpy(dtype=float)
                finite = np.isfinite(lfc_a) & np.isfinite(lfc_b)
                lfc_a = lfc_a[finite]
                lfc_b = lfc_b[finite]
                if len(lfc_a) < WITHIN_SOURCE_DIRECTION_MIN_SHARED_GENES:
                    continue
                if np.unique(lfc_a).size < 2 or np.unique(lfc_b).size < 2:
                    continue
                spearman = float(pd.Series(lfc_a).corr(pd.Series(lfc_b), method="spearman"))
                directional = (lfc_a != 0.0) & (lfc_b != 0.0)
                if not directional.any():
                    continue
                same_sign = float((np.sign(lfc_a[directional]) == np.sign(lfc_b[directional])).mean())
                if (
                    not np.isfinite(spearman)
                    or spearman > WITHIN_SOURCE_DIRECTION_SPEARMAN
                    or same_sign > WITHIN_SOURCE_DIRECTION_MAX_SIGN_AGREEMENT
                ):
                    continue
                warnings.append(
                    f"contrasts {contrast_a!r} and {contrast_b!r} inside source unit "
                    f"{str(source_unit)!r} have strongly reversed effect patterns across "
                    f"{len(lfc_a):,} shared genes (log2FC Spearman {spearman:.2f}; "
                    f"same-sign agreement {same_sign:.1%}). Verify that both effects use the declared "
                    "treatment-minus-control direction; values and ranks are unchanged."
                )
    return warnings


def _attach_source_quality_diagnostics(evidence: pd.DataFrame, diagnostics: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty or diagnostics.empty:
        return evidence
    mapping = diagnostics.set_index("source_unit_id")
    out = evidence.copy()
    for column in ["source_coherence_weight", "source_recommended_weight", "source_reliability_weight", "source_outlier_flag"]:
        out[column] = out["source_unit_id"].map(mapping[column]).fillna(out[column] if column in out.columns else 1.0)
    if "source_quality_weight" in mapping.columns:
        out["source_quality_weight"] = out["source_unit_id"].map(mapping["source_quality_weight"]).fillna(out["source_quality_weight"])
    out["source_quality_label"] = _quality_label_frame(out["source_quality_weight"])
    out["source_reliability_label"] = _quality_label_frame(out["source_reliability_weight"])
    return out


def _quality_weighted_consensus(
    evidence: pd.DataFrame,
    *,
    total_source_quality_weight: float,
    component_weights: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame()

    frame = evidence.copy()
    frame["normalized_rank"] = validate_normalized_rank(
        frame,
        context="quality-weighted consensus",
    )
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
    frame["_log_rank"] = np.log(frame["normalized_rank"])
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

    signs = frame[["gene_symbol", "source_unit_id", "signed_z", "_effective_weight"]].merge(
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
    grouped["quality_weighted_degora_score"] = _weighted_geometric_score_from_components(
        grouped,
        "quality_",
        component_weights,
    )
    # The five quality_* component columns are returned alongside the score. They
    # are the inputs the primary lane is actually built from, so an ablation or a
    # weight sweep can re-score from them arithmetically instead of re-running the
    # pipeline. Without them a downstream analysis would be forced to fall back on
    # the identically-named unweighted components of degora_score, which belong to
    # a different ranking.
    return grouped[
        [
            "gene_symbol",
            "quality_weighted_degora_score",
            "quality_weighted_consensus_direction",
            "quality_weighted_sign_concordance",
            "source_quality_support_score",
            "source_quality_weight_sum",
            *QUALITY_STATISTIC_COLUMNS,
            *QUALITY_COMPONENT_COLUMNS,
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
    min_studies = validate_min_studies(min_studies)
    columns = ["gene_symbol", "rra_rho", "rra_neglog10_rho", "rra_rank"]
    if evidence.empty or total_source_units <= 0:
        return pd.DataFrame(columns=columns)
    frame = evidence.copy()
    frame["normalized_rank"] = validate_normalized_rank(
        frame,
        context="RRA consensus",
    )
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "normalized_rank"])
    if frame.empty:
        return pd.DataFrame(columns=columns)
    n_lists = max(int(total_source_units), 1)
    ln10 = float(np.log(10.0))
    rank_floor = float(np.finfo(float).tiny)
    # One row per (gene, source unit), ordered by rank within the gene, so the
    # k-th smallest rank of every gene is scored in one vectorised beta.logcdf
    # call instead of a Python loop over tens of thousands of genes.
    ordered = (
        frame.drop_duplicates(["gene_symbol", "source_unit_id"])[["gene_symbol", "normalized_rank"]]
        .sort_values(["gene_symbol", "normalized_rank"], kind="mergesort")
        .reset_index(drop=True)
    )
    ordered["_order"] = ordered.groupby("gene_symbol", sort=False).cumcount() + 1
    counts = ordered.groupby("gene_symbol", sort=False)["_order"].transform("size")
    ordered = ordered.loc[counts.ge(min_studies) & ordered["_order"].le(n_lists)].copy()
    if ordered.empty:
        return pd.DataFrame(columns=columns)
    # Beta order-statistic RRA score in LOG space. The linear-space beta.cdf
    # underflows to 0.0 for genes ranked near-top across many lists, which
    # collapses the strongest genes into a single alphabetically-broken tie.
    # Tracking log(rho) keeps the top-of-list ordering and a usable magnitude.
    order = ordered["_order"].to_numpy(dtype=float)
    ordered["_log_score"] = beta.logcdf(
        np.maximum(ordered["normalized_rank"].to_numpy(dtype=float), rank_floor),
        order,
        n_lists - order + 1.0,
    )
    out = ordered.groupby("gene_symbol", sort=False, as_index=False)["_log_score"].min().rename(columns={"_log_score": "_log_rho"})
    out["gene_symbol"] = out["gene_symbol"].astype(str)
    if out.empty:
        return pd.DataFrame(columns=columns)
    out["_log_rho"] = pd.to_numeric(out["_log_rho"], errors="coerce")
    # Supported SciPy builds can disagree by a few ULPs in beta.logcdf even
    # when the underlying normalized ranks represent the same decimal value.
    # Quantize only the private ordering key: the raw log-rho still determines
    # the reported rho magnitudes, while 12-decimal ties fall through to the
    # deterministic gene-symbol key.
    out["_rra_sort_key"] = out["_log_rho"].round(12)
    out = out.sort_values(["_rra_sort_key", "gene_symbol"], ascending=[True, True]).reset_index(drop=True)
    out["rra_rank"] = np.arange(1, len(out) + 1, dtype=int)
    out["rra_rho"] = np.exp(out["_log_rho"].to_numpy(dtype=float)).clip(0.0, 1.0)
    out["rra_neglog10_rho"] = (-out["_log_rho"] / ln10).clip(lower=0.0)
    out.loc[out["rra_neglog10_rho"] == 0.0, "rra_neglog10_rho"] = 0.0
    return out[columns]


def _effect_meta_layer(evidence: pd.DataFrame, *, min_studies: int) -> pd.DataFrame:
    min_studies = validate_min_studies(min_studies)
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

    # DerSimonian-Laird random-effects pooling, vectorised over genes: every
    # per-gene sum is a groupby sum over the (gene, source unit) rows, which
    # replaced a Python loop that dominated the run time of a large corpus.
    frame = frame.drop_duplicates(["gene_symbol", "source_unit_id"]).copy()
    frame["_y"] = frame["lfc"].astype(float)
    frame["_vi"] = frame["_effect_se"].astype(float) ** 2
    frame = frame.loc[np.isfinite(frame["_y"]) & np.isfinite(frame["_vi"]) & frame["_vi"].gt(0)].copy()
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame["_w"] = 1.0 / frame["_vi"]
    frame["_wy"] = frame["_w"] * frame["_y"]
    frame["_w2"] = frame["_w"] ** 2
    grouped = frame.groupby("gene_symbol", sort=False)
    fixed = grouped[["_w", "_wy", "_w2"]].sum()
    fixed["k"] = grouped.size()
    fixed = fixed.loc[fixed["k"].ge(min_studies) & np.isfinite(fixed["_w"]) & fixed["_w"].gt(EFFECT_META_MIN_WEIGHT_SUM)]
    if fixed.empty:
        return pd.DataFrame(columns=columns)
    fixed["fixed"] = fixed["_wy"] / fixed["_w"]
    frame = frame.loc[frame["gene_symbol"].isin(fixed.index)].copy()
    frame["_fixed"] = frame["gene_symbol"].map(fixed["fixed"]).astype(float)
    frame["_q_part"] = frame["_w"] * (frame["_y"] - frame["_fixed"]) ** 2
    q = frame.groupby("gene_symbol", sort=False)["_q_part"].sum()
    fixed["q"] = q.reindex(fixed.index).astype(float)
    fixed["df"] = (fixed["k"] - 1).clip(lower=0)
    fixed["c"] = fixed["_w"] - fixed["_w2"] / fixed["_w"]
    with np.errstate(divide="ignore", invalid="ignore"):
        tau2 = (fixed["q"] - fixed["df"]) / fixed["c"]
    fixed["tau2"] = np.where((fixed["c"] > 0) & (fixed["df"] > 0), np.maximum(tau2, 0.0), 0.0)
    frame["_tau2"] = frame["gene_symbol"].map(fixed["tau2"]).astype(float)
    frame["_w_re"] = 1.0 / (frame["_vi"] + frame["_tau2"])
    frame["_w_re_y"] = frame["_w_re"] * frame["_y"]
    random = frame.groupby("gene_symbol", sort=False)[["_w_re", "_w_re_y"]].sum()
    fixed["sum_w_re"] = random["_w_re"].reindex(fixed.index).astype(float)
    fixed["pooled"] = (random["_w_re_y"] / random["_w_re"]).reindex(fixed.index).astype(float)
    fixed = fixed.loc[np.isfinite(fixed["sum_w_re"]) & fixed["sum_w_re"].gt(EFFECT_META_MIN_WEIGHT_SUM)]
    if fixed.empty:
        return pd.DataFrame(columns=columns)
    fixed["pooled_se"] = np.sqrt(1.0 / fixed["sum_w_re"])
    frame = frame.loc[frame["gene_symbol"].isin(fixed.index)].copy()
    frame["_pooled"] = frame["gene_symbol"].map(fixed["pooled"]).astype(float)
    frame["_hksj_part"] = frame["_w_re"] * (frame["_y"] - frame["_pooled"]) ** 2
    hksj = frame.groupby("gene_symbol", sort=False)["_hksj_part"].sum().reindex(fixed.index).astype(float)
    # Hartung-Knapp-Sidik-Jonkman small-sample CI, truncated so it is never
    # narrower than the normal random-effects SE. Most genes have only k=2-3
    # source units, where the normal-approx z interval is anti-conservative.
    k_eff = fixed["k"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        q_hksj = hksj.to_numpy(dtype=float) / ((k_eff - 1.0) * fixed["sum_w_re"].to_numpy(dtype=float))
    small = k_eff >= 2
    se_ci = np.where(small, np.maximum(np.sqrt(np.where(small, q_hksj, 0.0)), fixed["pooled_se"].to_numpy(dtype=float)), fixed["pooled_se"].to_numpy(dtype=float))
    crit = np.where(small, t_dist.ppf(0.975, np.maximum(k_eff - 1.0, 1.0)), 1.96)
    pooled = fixed["pooled"].to_numpy(dtype=float)
    q_values = fixed["q"].to_numpy(dtype=float)
    df_values = fixed["df"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        i2 = np.where((q_values > 0) & (df_values > 0), np.maximum((q_values - df_values) / q_values, 0.0), 0.0)
    out = pd.DataFrame(
        {
            "gene_symbol": fixed.index.astype(str),
            "effect_meta_log2fc_re": pooled,
            "effect_meta_se": se_ci,
            "effect_meta_ci_low": pooled - crit * se_ci,
            "effect_meta_ci_high": pooled + crit * se_ci,
            "effect_meta_tau2": fixed["tau2"].to_numpy(dtype=float),
            "effect_meta_i2": i2,
            "effect_meta_k": fixed["k"].to_numpy(dtype=int),
            "effect_meta_se_source": "derived_from_log2fc_and_two_sided_pvalue",
        }
    )
    finite = np.isfinite(out["effect_meta_log2fc_re"]) & np.isfinite(out["effect_meta_se"]) & np.isfinite(out["effect_meta_tau2"])
    out = out.loc[finite].reset_index(drop=True)
    return out[columns]


def _priority_components_from_evidence(
    evidence: pd.DataFrame,
    *,
    support_denominator: float,
    min_studies: int,
) -> pd.DataFrame:
    """Recompute the unweighted priority lane from collapsed source evidence.

    This is the fold scorer used by the leave-one-source-unit-out diagnostic. Its
    formulas intentionally mirror the full priority lane: sample-size weights enter
    Stouffer Z and the pooled log2 fold change, whereas normalized ranks and direction
    concordance remain source-unit-unweighted. Keeping those contracts identical
    prevents a scorer change from being misreported as leave-one-source instability.
    """

    min_studies = validate_min_studies(min_studies)
    if evidence.empty:
        return pd.DataFrame()
    frame = evidence.copy()
    frame["normalized_rank"] = validate_normalized_rank(
        frame,
        context="priority consensus",
    )
    for column in ["signed_z", "lfc", "normalized_rank", "weight"]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["gene_symbol", "source_unit_id", "signed_z", "normalized_rank", "weight"])
    if frame.empty:
        return pd.DataFrame()
    frame["_wz"] = frame["weight"] * frame["signed_z"]
    frame["_w2"] = frame["weight"] ** 2
    frame["_wlfc"] = frame["weight"] * frame["lfc"].fillna(0.0)
    frame["_w_lfc_denominator"] = np.where(frame["lfc"].notna(), frame["weight"], 0.0)
    frame["_log_rank"] = np.log(frame["normalized_rank"])
    grouped = frame.groupby("gene_symbol", as_index=False).agg(
        n_source_units=("source_unit_id", "nunique"),
        sum_wz=("_wz", "sum"),
        sum_w2=("_w2", "sum"),
        sum_wlfc=("_wlfc", "sum"),
        sum_w_lfc=("_w_lfc_denominator", "sum"),
        mean_log_rank=("_log_rank", "mean"),
    )
    grouped = grouped.loc[grouped["n_source_units"].ge(min_studies)].copy()
    if grouped.empty:
        return pd.DataFrame()
    grouped["stouffer_z"] = grouped["sum_wz"] / np.sqrt(grouped["sum_w2"])
    grouped["weighted_lfc"] = np.where(grouped["sum_w_lfc"].gt(0), grouped["sum_wlfc"] / grouped["sum_w_lfc"], np.nan)
    grouped["rank_product"] = np.exp(grouped["mean_log_rank"])
    grouped["support_score"] = (np.log1p(grouped["n_source_units"]) / support_denominator).clip(0.0, 1.0)
    grouped["evidence_score"] = _component_strength_from_z(grouped["stouffer_z"])
    grouped["rank_score_component"] = (1.0 - grouped["rank_product"]).fillna(0.0).clip(0.0, 1.0)
    grouped["effect_score"] = _component_strength_from_lfc(grouped["weighted_lfc"])

    signs = frame[["gene_symbol", "source_unit_id", "signed_z"]].merge(
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
    direction["direction_score"] = np.where(
        direction["total_strength"].gt(0),
        direction["concordant_strength"] / direction["total_strength"],
        0.0,
    )
    direction_confidence = _direction_confidence_from_evidence(
        frame,
        grouped[["gene_symbol", "stouffer_z"]],
        reference_z_column="stouffer_z",
        output_column="direction_confidence_index",
    )
    grouped = grouped.merge(
        direction[["gene_symbol", "direction_score"]],
        on="gene_symbol",
        how="left",
    ).merge(
        direction_confidence[["gene_symbol", "direction_confidence_index"]],
        on="gene_symbol",
        how="left",
    )
    grouped["direction_confidence_index"] = pd.to_numeric(
        grouped["direction_confidence_index"], errors="coerce"
    ).fillna(0.5)
    grouped["priority_score"] = _weighted_geometric_score_with_weights(grouped, PRIORITY_SCORE_WEIGHTS)
    return grouped[["gene_symbol", "priority_score", "n_source_units", "direction_confidence_index"]]


def _leave_one_source_out_stability(
    evidence: pd.DataFrame,
    scores: pd.DataFrame,
    *,
    support_denominator: float,
    min_studies: int,
) -> pd.DataFrame:
    min_studies = validate_min_studies(min_studies)
    columns = [
        "gene_symbol",
        "loo_total_folds",
        "loo_rank_evaluable_folds",
        "loo_penalty_folds",
        "loo_component_available",
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
        out["loo_total_folds"] = len(source_units)
        out["loo_rank_evaluable_folds"] = 0
        out["loo_penalty_folds"] = len(source_units)
        out["loo_component_available"] = False
        for column in columns[5:]:
            out[column] = np.nan
        return out

    total_genes = max(int(len(scores)), 1)
    penalty_rank = float(total_genes + 1)
    rank_records: dict[str, list[float]] = {str(gene): [] for gene in scores["gene_symbol"].astype(str)}
    eligible_records: dict[str, list[bool]] = {str(gene): [] for gene in scores["gene_symbol"].astype(str)}
    for source_unit in source_units:
        subset = evidence.loc[evidence["_loo_source_unit_id"].ne(source_unit)].copy()
        components = _priority_components_from_evidence(
            subset,
            support_denominator=support_denominator,
            min_studies=min_studies,
        )
        if components.empty:
            rank_map: dict[str, float] = {}
        else:
            components = components.loc[components["gene_symbol"].astype(str).isin(rank_records)].copy()
            components["priority_score"] = pd.to_numeric(
                components["priority_score"], errors="coerce"
            ).round(6)
            components["direction_confidence_index"] = pd.to_numeric(
                components["direction_confidence_index"], errors="coerce"
            ).round(6)
            ranked = components.sort_values(
                ["priority_score", "n_source_units", "direction_confidence_index", "gene_symbol"],
                ascending=[False, False, False, True],
            ).reset_index(drop=True)
            rank_map = dict(zip(ranked["gene_symbol"].astype(str), (ranked.index + 1).astype(float), strict=False))
        for gene in rank_records:
            eligible_records[gene].append(gene in rank_map)
            rank_records[gene].append(float(rank_map.get(gene, penalty_rank)))

    current_rank = scores.set_index("gene_symbol")["priority_rank"].astype(float).to_dict()
    rows: list[dict[str, Any]] = []
    for gene, ranks in rank_records.items():
        values = np.array(ranks, dtype=float)
        eligible_folds = int(sum(eligible_records[gene]))
        penalty_folds = len(source_units) - eligible_folds
        available = eligible_folds > 0
        if available:
            eligibility: NDArray[np.bool_] = np.asarray(eligible_records[gene], dtype=bool)
            evaluable_values = values[eligibility]
            median = float(np.median(evaluable_values))
            q75, q25 = np.percentile(evaluable_values, [75, 25])
            iqr = float(q75 - q25)
            penalized_median = float(np.median(values))
            shift = abs(penalized_median - float(current_rank.get(gene, penalty_rank)))
            stability = max(0.0, 1.0 - min(1.0, shift / max(float(total_genes), 1.0)))
            top50 = float(np.mean(evaluable_values <= 50.0))
            top100 = float(np.mean(evaluable_values <= 100.0))
        else:
            median = np.nan
            iqr = np.nan
            stability = np.nan
            top50 = np.nan
            top100 = np.nan
        rows.append(
            {
                "gene_symbol": gene,
                "loo_total_folds": len(source_units),
                "loo_rank_evaluable_folds": eligible_folds,
                "loo_penalty_folds": penalty_folds,
                "loo_component_available": available,
                "loo_median_rank": median,
                "loo_rank_iqr": iqr,
                "loo_rank_stability_score": stability,
                "loo_top50_fraction": top50,
                "loo_top100_fraction": top100,
            }
        )
    return pd.DataFrame.from_records(rows, columns=columns)


def _loo_stability_not_computed(scores: pd.DataFrame, total_source_units: int) -> pd.DataFrame:
    """LOO columns in their 'unavailable' shape, for callers that skip the diagnostic."""

    out = scores[["gene_symbol"]].copy()
    out["loo_total_folds"] = int(total_source_units)
    out["loo_rank_evaluable_folds"] = 0
    out["loo_penalty_folds"] = int(total_source_units)
    out["loo_component_available"] = False
    for column in ["loo_median_rank", "loo_rank_iqr", "loo_rank_stability_score", "loo_top50_fraction", "loo_top100_fraction"]:
        out[column] = np.nan
    return out


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
    ablation: ScoreAblation | None = None,
    include_loo_stability: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return gene scores, study-gene evidence rows, and score metadata.

    ``ablation`` re-scores the primary quality-weighted lane under a modified
    configuration; ``None`` reproduces the shipped ranking. See
    :class:`ScoreAblation` for why the gene universe and the source-unit set stay
    fixed across every variant.

    ``include_loo_stability=False`` skips the leave-one-source-out diagnostic,
    which re-scores the corpus once per source unit and dominates the run time
    of a large corpus. The primary rank does not depend on it; the LOO columns
    are then reported as unavailable, exactly as they are when no fold is
    rank-evaluable, and the reliability summary uses its three mandatory
    diagnostics. Ablation sweeps use this.
    """

    min_studies = validate_min_studies(min_studies)
    ablation = ablation or ScoreAblation()
    # Validate at the public Python/standalone scoring boundary before any row
    # filtering or temporal selection can hide malformed ranks.
    validate_normalized_rank(harmonized, context="DEGORA score table")
    _validate_active_evidence_contract(harmonized)

    score_harmonized, n_nonfinite_lfc_capped = _score_ready_harmonized(harmonized)
    evidence = study_gene_evidence(score_harmonized)
    # Component/weight ablations repeatedly score the same harmonized evidence.
    # These provenance advisories are invariant to every ablation and are not a
    # score input, so compute them only for the canonical run rather than paying
    # the source/contrast pair cost once per variant.
    if ablation.is_default:
        advisory_diagnostics_status = "computed"
        source_similarity_warnings = near_duplicate_source_unit_warnings(evidence)
        within_source_warnings = within_source_direction_warnings(score_harmonized)
    else:
        advisory_diagnostics_status = "not_computed_for_nondefault_ablation"
        source_similarity_warnings = []
        within_source_warnings = []
    source_quality_diagnostics = _source_quality_diagnostics_from_evidence(evidence)
    evidence = _attach_source_quality_diagnostics(evidence, source_quality_diagnostics)
    if ablation.disable_source_quality_weighting:
        # Neutralise every source-quality channel to 1.0 on both the evidence rows
        # and the source-unit diagnostics. Applied here, after the diagnostics are
        # derived, so the sample-size and coherence inputs are computed exactly as
        # in the default run and only their downstream weighting is removed.
        for column in _SOURCE_QUALITY_WEIGHT_COLUMNS:
            if column in evidence.columns:
                evidence[column] = 1.0
            if column in source_quality_diagnostics.columns:
                source_quality_diagnostics[column] = 1.0
    if ablation.disable_sample_size_weighting:
        # The primary lane multiplies this by the reliability weight to form the
        # effective per-source weight; setting it to 1.0 gives every source unit an
        # equal vote regardless of its replicate count.
        evidence["weight"] = 1.0
    consensus = slice_consensus(score_harmonized, min_studies=min_studies)
    if consensus.empty:
        metadata = {
            "score_version": SCORE_VERSION,
            "score_formula": SCORE_FORMULA,
            "support_normalization_rule": SUPPORT_NORMALIZATION_RULE,
            "score_weights": ablation.weights,
            "score_ablation": ablation.to_dict(),
            "primary_rank_column": PRIMARY_RANK_COLUMN,
            "primary_score_column": PRIMARY_SCORE_COLUMN,
            "primary_rank_interpretation": PRIMARY_RANK_DESCRIPTION,
            "priority_score_weights": PRIORITY_SCORE_WEIGHTS,
            "evidence_reliability_score_weights": RELIABILITY_SCORE_WEIGHTS,
            "evidence_reliability_score_version": EVIDENCE_RELIABILITY_VERSION,
            "evidence_reliability_score_rule": EVIDENCE_RELIABILITY_RULE,
            "min_studies": min_studies,
            "n_gene_scores": 0,
            "n_source_unit_gene_evidence_rows": int(len(evidence)),
            "n_source_units_total": int(evidence["source_unit_id"].nunique()) if not evidence.empty else 0,
            "n_contrasts_total": int(harmonized["study_id"].nunique()) if "study_id" in harmonized.columns else 0,
            "n_nonfinite_lfc_capped_for_score": n_nonfinite_lfc_capped,
            "independent_unit_for_consensus": "source_unit_id (paper_id when available, otherwise study_id)",
            "source_unit_collapse_rule": SOURCE_UNIT_COLLAPSE_RULE,
            "stouffer_weight_rule": STOUFFER_WEIGHT_RULE,
            "direction_concordance_rule": "|signed_z|-strength-weighted concordance across independent source-unit representatives; the primary path is source-weight-free, while the quality-weighted variant additionally multiplies the strength by the source reliability weight",
            "heterogeneity_rule": HETEROGENEITY_RULE,
            "heterogeneity_flag_rule": "heterogeneity_i2 >= 0.75 is labeled high_context_dependent_review; >= 0.50 is moderate_context_review; flags are descriptive review aids, not score gates",
            "quality_weighted_score_formula": QUALITY_WEIGHTED_SCORE_FORMULA,
            "source_quality_weight_rules": {
                "source_input_type_weights": SOURCE_INPUT_TYPE_QUALITY_WEIGHTS,
                "table_scope_multipliers": TABLE_SCOPE_QUALITY_MULTIPLIERS,
                "replicate_multiplier": REPLICATE_MULTIPLIER_RULE,
                "replicate_multiplier_reachability": (
                    "the CLI rejects a zero, negative or fractional group size during validation, "
                    "so the 0.35 zero-count branch is reachable only through the Python API"
                ),
                "source_coherence_guardrail": "gold-panel-free source-source LFC Spearman check; low-quality sources with median pairwise Spearman < 0.05 receive source_coherence_weight=0.50 in the default primary quality-weighted lane; high-quality incoherence and direction-conflict flags remain advisory and change no weight or rank",
                "source_reliability_shrinkage": "primary quality-weighted-lane source weight shrunk toward neutral 0.65 using source gene coverage and pairwise-comparison evidence; not a calibrated probability",
            },
            "near_duplicate_source_unit_rule": NEAR_DUPLICATE_SOURCE_RULE,
            "near_duplicate_source_unit_warnings": source_similarity_warnings,
            "within_source_direction_rule": WITHIN_SOURCE_DIRECTION_RULE,
            "within_source_direction_warnings": within_source_warnings,
            "advisory_diagnostics_status": advisory_diagnostics_status,
            "direction_confidence_rule": "Beta(1,1)-shrunk source-unit count concordance against the reported consensus signed-z direction: (1 + concordant source units) / (2 + observed source units). When the consensus z is exactly 0 the direction is a tie and every source unit is credited one half rather than zero, so the index is 0.5 rather than the 0.25 the formula alone would give, and direction_concordant_source_units carries that half-credit and is not a whole number. Quality-weighted direction confidence uses reliability-weighted pseudo-counts against the quality-weighted consensus direction and is not a calibrated posterior probability",
            "random_effects_stouffer_rule": RANDOM_EFFECTS_STOUFFER_RULE,
            "stouffer_inference_warning": _stouffer_inference_warning(evidence),
            "rra_rule": "parallel rank lane using beta order-statistic RobustRankAggreg-style rho over source-unit normalized ranks; missing source-unit lists are handled through the total source-unit universe; rho is computed in log space, and rra_rank sorts log-rho rounded to 12 decimal places then gene_symbol while reported rho magnitudes retain the unquantized calculation; rra_neglog10_rho (-log10 rho) preserves ordering for top genes whose rho underflows to 0; rho is not reported as a calibrated FDR",
            "effect_meta_rule": EFFECT_META_RULE,
            "effect_meta_small_k_warning": "For effect_meta_k = 2 the HKSJ t critical value is 12.71, so the interval is wide enough to be uninformative in practice: it will usually span zero whatever the pooled estimate is. For k = 3 it is 4.30. Read these intervals as descriptive only, and do not read an interval covering zero at small k as evidence of no effect.",
            "loo_stability_rule": LOO_STABILITY_RULE,
            "source_quality_diagnostics": _diagnostic_records(source_quality_diagnostics),
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
    denominator = np.log1p(total_source_units) if total_source_units > 0 else 1.0
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
        component_weights=ablation.component_weights,
    )
    if quality_consensus.empty:
        for column in [
            "quality_weighted_degora_score",
            "quality_weighted_consensus_direction",
            "quality_weighted_sign_concordance",
            "source_quality_support_score",
            "source_quality_weight_sum",
            *QUALITY_STATISTIC_COLUMNS,
            *QUALITY_COMPONENT_COLUMNS,
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
            *QUALITY_STATISTIC_COLUMNS,
            *QUALITY_COMPONENT_COLUMNS,
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
    # Rank on the exact values that are published.  Ranking first and rounding
    # later made a small number of rows impossible to reconstruct from the CSV:
    # two distinct internal floats could both be emitted as (for example)
    # 47.347037, while their rank still reflected the hidden sub-6-decimal
    # difference.  Rounding every displayed ordering component before rank
    # assignment makes the documented (score, suffix keys, symbol) contract
    # complete and independently reproducible.
    scores["priority_score"] = scores["priority_score"].round(6)
    scores["direction_confidence_index"] = scores["direction_confidence_index"].round(6)
    priority_ranked = scores.sort_values(
        ["priority_score", "n_source_units", "direction_confidence_index", "gene_symbol"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    priority_rank_map = pd.Series(np.arange(1, len(priority_ranked) + 1), index=priority_ranked["gene_symbol"])
    scores["priority_rank"] = scores["gene_symbol"].map(priority_rank_map).fillna(0).astype(int)
    if include_loo_stability:
        stability = _leave_one_source_out_stability(
            evidence,
            scores,
            support_denominator=denominator,
            min_studies=min_studies,
        )
    else:
        stability = _loo_stability_not_computed(scores, total_source_units)
    scores = scores.merge(stability, on="gene_symbol", how="left")
    for column in ["loo_median_rank", "loo_rank_iqr", "loo_rank_stability_score", "loo_top50_fraction", "loo_top100_fraction"]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce")
    for column in ["loo_total_folds", "loo_rank_evaluable_folds", "loo_penalty_folds"]:
        numeric = pd.to_numeric(scores[column], errors="coerce")
        if (~np.isfinite(numeric.to_numpy(dtype=float))).any():
            raise ValueError(f"non-finite LOO fold count in {column}")
        scores[column] = numeric.astype(int)
    if not (
        scores["loo_total_folds"]
        == scores["loo_rank_evaluable_folds"] + scores["loo_penalty_folds"]
    ).all():
        raise ValueError("LOO fold counts do not reconcile")
    scores["loo_component_available"] = scores["loo_component_available"].astype(bool)
    if not (
        scores["loo_component_available"]
        == scores["loo_rank_evaluable_folds"].gt(0)
    ).all():
        raise ValueError("LOO availability does not match rank-evaluable fold counts")
    reliability_components = scores.copy()
    reliability_components["source_quality_support_score"] = pd.to_numeric(
        reliability_components["source_quality_support_score"], errors="coerce"
    )
    (
        scores["evidence_reliability_score"],
        scores["evidence_reliability_components_used"],
    ) = _evidence_reliability_score(
        reliability_components
    )
    scores["degora_score"] = scores["degora_score"].round(6)
    scores["quality_weighted_degora_score"] = scores["quality_weighted_degora_score"].round(6)
    scores["quality_weighted_sign_concordance"] = scores["quality_weighted_sign_concordance"].round(6)
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
    # top_percent_label is the human-readable companion to the primary
    # rank_label (quality_weighted_degora_rank), so derive it from the primary
    # quality_weighted_top_percent rather than the unweighted screening top_percent;
    # otherwise the label printed beside the primary rank reports the wrong percentile.
    scores["top_percent_label"] = scores["quality_weighted_top_percent"].map(_format_top_percent)
    scores["evidence_tier"] = _evidence_tier(
        scores["quality_weighted_top_percent"],
        scores["n_source_units"],
        scores["quality_weighted_sign_concordance"],
        total_source_units=total_source_units,
    )
    scores["support_label"] = scores["n_source_units"].map(lambda units: f"{int(units):,} / {total_source_units:,} source units")
    scores["direction_label"] = (
        (scores["quality_weighted_sign_concordance"] * 100.0).map(_format_percent)
        + " "
        + scores["quality_weighted_consensus_direction"].astype(str)
        + "-concordant"
    )
    scores["evidence_reliability_score"] = scores["evidence_reliability_score"].round(6)
    scores["quality_weighted_direction_confidence_index"] = scores["quality_weighted_direction_confidence_index"].round(6)
    scores["direction_concordant_source_units"] = scores["direction_concordant_source_units"].round(3)
    scores["direction_total_source_units"] = scores["direction_total_source_units"].round(3)
    scores["direction_posterior_mean"] = scores["direction_posterior_mean"].round(6)
    scores["loo_median_rank"] = scores["loo_median_rank"].round(3)
    scores["loo_rank_iqr"] = scores["loo_rank_iqr"].round(3)
    scores["loo_rank_stability_score"] = scores["loo_rank_stability_score"].round(6)
    scores["loo_top50_fraction"] = scores["loo_top50_fraction"].round(6)
    scores["loo_top100_fraction"] = scores["loo_top100_fraction"].round(6)
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
    # Decimal-place rounding collapsed every finite probability below 5e-13 to
    # exactly zero in CSV/SQLite output. Preserve scipy/BH's raw float tail; the
    # values are auxiliary and do not participate in any score or rank.
    scores["re_stouffer_p"] = scores["re_stouffer_p"].fillna(1.0)
    scores["re_stouffer_padj"] = scores["re_stouffer_padj"].fillna(1.0)
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
    ]:
        scores[column] = pd.to_numeric(scores[column], errors="coerce").round(6)
    scores["effect_meta_k"] = pd.to_numeric(scores["effect_meta_k"], errors="coerce").fillna(0).astype(int)
    scores["effect_meta_se_source"] = scores["effect_meta_se_source"].fillna("")

    metadata = {
        "score_version": SCORE_VERSION,
        "score_formula": SCORE_FORMULA,
        "support_normalization_rule": SUPPORT_NORMALIZATION_RULE,
        "score_weights": ablation.weights,
        "score_ablation": ablation.to_dict(),
        "primary_rank_column": PRIMARY_RANK_COLUMN,
        "primary_score_column": PRIMARY_SCORE_COLUMN,
        "primary_rank_interpretation": PRIMARY_RANK_DESCRIPTION,
        "priority_score_weights": PRIORITY_SCORE_WEIGHTS,
        "evidence_reliability_score_weights": RELIABILITY_SCORE_WEIGHTS,
        "evidence_reliability_score_version": EVIDENCE_RELIABILITY_VERSION,
        "evidence_reliability_score_rule": EVIDENCE_RELIABILITY_RULE,
        "quality_weighted_score_formula": QUALITY_WEIGHTED_SCORE_FORMULA,
        "quality_weighted_score_warning": "Quality-weighted ranking is the default browser and API ordering; degora_score remains available as the unweighted/reference prioritization index.",
        "source_quality_weight_rules": {
            "source_input_type_weights": SOURCE_INPUT_TYPE_QUALITY_WEIGHTS,
            "table_scope_multipliers": TABLE_SCOPE_QUALITY_MULTIPLIERS,
            "replicate_multiplier": REPLICATE_MULTIPLIER_RULE,
            "replicate_multiplier_reachability": (
                "the CLI rejects a zero, negative or fractional group size during validation, "
                "so the 0.35 zero-count branch is reachable only through the Python API"
            ),
            "source_coherence_guardrail": "gold-panel-free source-source LFC Spearman check; low-quality sources with median pairwise Spearman < 0.05 receive source_coherence_weight=0.50 in the default primary quality-weighted lane; high-quality incoherence and direction-conflict flags remain advisory and change no weight or rank",
            "source_reliability_shrinkage": "primary quality-weighted-lane source weight shrunk toward neutral 0.65 using source gene coverage and pairwise-comparison evidence; not a calibrated probability",
        },
        "near_duplicate_source_unit_rule": NEAR_DUPLICATE_SOURCE_RULE,
        "near_duplicate_source_unit_warnings": source_similarity_warnings,
        "within_source_direction_rule": WITHIN_SOURCE_DIRECTION_RULE,
        "within_source_direction_warnings": within_source_warnings,
        "advisory_diagnostics_status": advisory_diagnostics_status,
        "min_studies": min_studies,
        "n_gene_scores": int(len(scores)),
        "n_source_unit_gene_evidence_rows": int(len(evidence)),
        "n_source_units_total": total_source_units,
        "n_source_quality_outliers": int(source_quality_diagnostics["source_outlier_flag"].sum()) if not source_quality_diagnostics.empty else 0,
        "n_source_direction_conflicts": int(source_quality_diagnostics["source_direction_conflict_flag"].sum()) if not source_quality_diagnostics.empty else 0,
        "source_direction_conflict_rule": DIRECTION_CONFLICT_RULE,
        "n_contrasts_total": int(harmonized["study_id"].nunique()) if "study_id" in harmonized.columns else 0,
        "n_nonfinite_lfc_capped_for_score": n_nonfinite_lfc_capped,
        "independent_unit_for_consensus": "source_unit_id (paper_id when available, otherwise study_id)",
        "source_unit_collapse_rule": SOURCE_UNIT_COLLAPSE_RULE,
        "stouffer_weight_rule": STOUFFER_WEIGHT_RULE,
        "direction_concordance_rule": "evidence-strength-weighted concordance across independent source-unit representatives",
        "heterogeneity_rule": HETEROGENEITY_RULE,
        "heterogeneity_flag_rule": "heterogeneity_i2 >= 0.75 is labeled high_context_dependent_review; >= 0.50 is moderate_context_review; flags are descriptive review aids, not score gates",
        "rank_interpretation": "degora_rank is the absolute rank among scored genes; top_percent is rank / total scored genes * 100, so smaller is more selective; percentile is 100 for the top-ranked gene and decreases with rank.",
        "high_confidence_rule": f"relative browsing flag: n_source_units >= min(2, total_source_units) = {high_confidence_min_units}, sign_concordance >= 0.75, rank_score_component >= 0.80, evidence_score >= 0.50; does not use stouffer_padj as a calibrated inferential gate",
        "direction_confidence_rule": "Beta(1,1)-shrunk source-unit count concordance against the reported consensus signed-z direction: (1 + concordant source units) / (2 + observed source units). When the consensus z is exactly 0 the direction is a tie and every source unit is credited one half rather than zero, so the index is 0.5 rather than the 0.25 the formula alone would give, and direction_concordant_source_units carries that half-credit and is not a whole number. Quality-weighted direction confidence uses reliability-weighted pseudo-counts against the quality-weighted consensus direction and is not a calibrated posterior probability",
        "random_effects_stouffer_rule": RANDOM_EFFECTS_STOUFFER_RULE,
        "stouffer_inference_warning": _stouffer_inference_warning(evidence),
        "rra_rule": "parallel rank lane using beta order-statistic RobustRankAggreg-style rho over source-unit normalized ranks; missing source-unit lists are handled through the total source-unit universe; rho is computed in log space, and rra_rank sorts log-rho rounded to 12 decimal places then gene_symbol while reported rho magnitudes retain the unquantized calculation; rra_neglog10_rho (-log10 rho) preserves ordering for top genes whose rho underflows to 0; rho is not reported as a calibrated FDR",
        "effect_meta_rule": EFFECT_META_RULE,
        "effect_meta_small_k_warning": "For effect_meta_k = 2 the HKSJ t critical value is 12.71, so the interval is wide enough to be uninformative in practice: it will usually span zero whatever the pooled estimate is. For k = 3 it is 4.30. Read these intervals as descriptive only, and do not read an interval covering zero at small k as evidence of no effect.",
        "loo_stability_rule": LOO_STABILITY_RULE,
        "evidence_tier_rules": {
            "basis": "quality_weighted_top_percent and quality_weighted_sign_concordance (the primary lane)",
            "A": "quality_weighted_top_percent <= 1, n_source_units >= min(3, total_source_units), quality_weighted_sign_concordance >= 0.90",
            "B": "quality_weighted_top_percent <= 5, n_source_units >= min(2, total_source_units), quality_weighted_sign_concordance >= 0.75",
            "C": "quality_weighted_top_percent <= 20, n_source_units >= min(2, total_source_units)",
            "D": "lower-ranked or weakly supported",
        },
        "score_warning": "DEGORA score is for transparent prioritization and browsing, not a calibrated probability or a validation metric.",
        "source_quality_diagnostics": _diagnostic_records(source_quality_diagnostics),
    }
    # The five quality_* components are appended rather than folded into
    # GENE_SCORE_COLUMNS so the established column order is untouched and existing
    # readers keep working. They are exported because the components already in
    # GENE_SCORE_COLUMNS (support_score, direction_score, …) belong to the
    # unweighted degora_score; an ablation or weight sweep built on those would be
    # describing a different ranking from the primary one.
    output_columns = [
        *GENE_SCORE_COLUMNS,
        *(c for c in QUALITY_COMPONENT_COLUMNS if c in scores.columns),
        *(c for c in QUALITY_STATISTIC_COLUMNS if c in scores.columns),
    ]
    return scores[output_columns], evidence, metadata


STUDY_TABLE_COLUMNS = [
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
    "condition",
    "duration_h",
    "n_ctrl",
    "n_treat",
    "source_path",
    "source_url",
    "notes",
]


def _active_study_table(catalog_path: Path | None, harmonized: pd.DataFrame) -> pd.DataFrame:
    """Return one row per active contrast for the SQLite ``studies`` table.

    Both branches keep the same grain -- one row per contrast -- and the same
    column set, so the table does not change shape or meaning depending on
    whether a catalog was supplied.
    """

    if catalog_path is not None and catalog_path.exists():
        catalog = read_catalog(catalog_path)
        active = catalog.loc[catalog_include_mask(catalog)].copy()
        active["source_unit_id"] = _source_unit_series(active)
        if "condition" not in active.columns and "hypoxia_modality" in active.columns:
            active["condition"] = active["hypoxia_modality"]
        return (
            active[[column for column in STUDY_TABLE_COLUMNS if column in active.columns]]
            .sort_values("study_id")
            .reset_index(drop=True)
        )

    # Derive from the harmonized contrast rows, not from the collapsed evidence
    # table. Evidence carries one row per (gene, source unit) and its study_id is
    # the first contributing contrast for that gene, so a contrast that is never
    # any gene's first -- for example a follow-up time point covering fewer genes
    # than its sibling in the same source unit -- would never appear here, and the
    # reported contrast count would disagree with n_contrasts_total.
    frame = harmonized.copy()
    if frame.empty:
        return pd.DataFrame(columns=STUDY_TABLE_COLUMNS)
    frame["source_unit_id"] = _source_unit_series(frame)
    if "condition" not in frame.columns and "hypoxia_modality" in frame.columns:
        frame["condition"] = frame["hypoxia_modality"]
    for column in STUDY_TABLE_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    return (
        frame[STUDY_TABLE_COLUMNS]
        .drop_duplicates("study_id")
        .sort_values("study_id")
        .reset_index(drop=True)
    )


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
    fd, tmp_name = tempfile.mkstemp(prefix=f".{db_path.name}.", suffix=".tmp", dir=db_path.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        connection = sqlite3.connect(tmp_path)
        try:
            with connection:
                gene_scores.to_sql("genes", connection, index=False)
                evidence.to_sql("gene_evidence", connection, index=False)
                studies.to_sql("studies", connection, index=False)
                meta_rows = [
                    {
                        "key": key,
                        "value": _strict_json_dumps(value, sort_keys=True)
                        if not isinstance(value, str)
                        else value,
                    }
                    for key, value in metadata.items()
                ]
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
        apply_default_file_mode(tmp_path)
        os.replace(tmp_path, db_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


NOMINAL_SIGNIFICANCE_THRESHOLD = 0.05


def _corpus_significance_warnings(gene_scores: pd.DataFrame) -> list[str]:
    """Warn when the evidence tiers are ordering an entirely non-significant corpus.

    DEGORA scores are deliberately relative, so the strongest genes in a corpus
    with no signal still land in tier A. Without a note, "Evidence tier A" over a
    best adjusted p-value of 0.26 reads as a finding.
    """

    if gene_scores.empty or "stouffer_padj" not in gene_scores.columns:
        return []
    padj = pd.to_numeric(gene_scores["stouffer_padj"], errors="coerce").dropna()
    if padj.empty:
        return []
    best = float(padj.min())
    if best <= NOMINAL_SIGNIFICANCE_THRESHOLD:
        return []
    top_tier = ""
    if "evidence_tier" in gene_scores.columns:
        tiers = sorted(str(value) for value in gene_scores["evidence_tier"].dropna().unique())
        top_tier = tiers[0] if tiers else ""
    tier_note = (
        f" The strongest genes are still reported in tier {top_tier}, because tiers rank genes"
        " within this corpus rather than against a significance threshold."
        if top_tier
        else ""
    )
    return [
        f"No gene reached adjusted significance: the smallest stouffer_padj is {best:.3g}, above "
        f"{NOMINAL_SIGNIFICANCE_THRESHOLD}.{tier_note} Treat this ranking as a relative ordering of weak "
        "evidence, not as a set of findings."
    ]


def direction_conflict_warnings(diagnostics: pd.DataFrame | list[dict[str, Any]]) -> list[str]:
    """Name the source units whose direction runs against the rest of the corpus."""

    frame = pd.DataFrame.from_records(diagnostics) if isinstance(diagnostics, list) else diagnostics
    if frame is None or frame.empty or "source_direction_conflict_flag" not in frame.columns:
        return []
    flagged = frame.loc[frame["source_direction_conflict_flag"].astype(bool)]
    if flagged.empty:
        return []
    n_units = int(len(frame))

    def _median_text(record: dict[str, Any]) -> str:
        median = pd.to_numeric(
            pd.Series([record.get("median_pairwise_lfc_spearman")]), errors="coerce"
        ).iloc[0]
        return f"{float(median):.2f}" if pd.notna(median) else "n/a"

    closing = (
        " DEGORA never reverses an effect column, and a reversed source votes against every gene it "
        "shares. Weights and ranks are unchanged."
    )
    records = flagged.to_dict(orient="records")

    # When more than two units are flagged and they are not a minority, the
    # corpus is divided rather than one source being wrong, and DEGORA cannot
    # tell which half carries the intended convention. Repeating "this source
    # disagrees with the other source units" once per unit both contradicts
    # itself and sends the reader to check every one of them; one warning about
    # the corpus is the honest shape.
    if n_units > 2 and len(records) * 2 >= n_units:
        named = ", ".join(
            f"{record.get('source_unit_id')!r} ({_median_text(record)})" for record in records
        )
        return [
            f"{len(records)} of {n_units} source units disagree in direction with the rest, so this "
            f"corpus is split rather than one source being reversed: {named} (median pairwise log2FC "
            "Spearman in brackets). DEGORA cannot tell which half carries the intended convention, so "
            "check how the contrast is written in each of them." + closing
        ]

    warnings: list[str] = []
    for record in records:
        if n_units == 2:
            detail = (
                "with only two source units DEGORA cannot tell which one is reversed; check the contrast "
                "direction of both"
            )
        else:
            detail = "check whether this source's contrast is written control-minus-treatment"
        warnings.append(
            f"source unit {record.get('source_unit_id')!r} disagrees in direction with the other source units "
            f"(median pairwise log2FC Spearman {_median_text(record)}); {detail}." + closing
        )
    return warnings


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

    resolved_output = Path(output_dir).resolve()
    resolved_db = (Path(db_path) if db_path is not None else resolved_output / "degora_scores.db").resolve()
    # The database may intentionally live outside output_dir. Every writer must
    # claim its target identity even for the default location, so a second run
    # that names the same file through --db cannot publish a different database
    # beside this run's CSV/metadata generation. Artifact-output locking always
    # precedes the database target lock, matching the outer CLI and discovery
    # pipeline order; the contexts remain same-thread re-entrant.
    with artifact_output_lock(resolved_output):
        with output_directory_lock(publication_target_lock_path(resolved_db)):
            return _write_score_database_locked(
                harmonized_path,
                resolved_output,
                catalog_path=catalog_path,
                db_path=resolved_db,
                min_studies=min_studies,
                command=command,
                extra_metadata=extra_metadata,
            )


def _write_score_database_locked(
    harmonized_path: Path,
    output_dir: Path,
    *,
    catalog_path: Path | None = None,
    db_path: Path | None = None,
    min_studies: int = 2,
    command: str | None = None,
    extra_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    harmonized_path = harmonized_path.resolve()
    catalog_path = catalog_path.resolve() if catalog_path is not None else None
    output_dir = output_dir.resolve()
    db_path = (db_path or output_dir / "degora_scores.db").resolve()
    path_base = output_dir
    source_path_base = catalog_path.parent if catalog_path is not None else harmonized_path.parent

    if harmonized_path.suffix.lower() in {".parquet", ".pq"}:
        harmonized = pd.read_parquet(harmonized_path)
    else:
        harmonized = pd.read_csv(harmonized_path, low_memory=False)
        harmonized = restore_formula_text_if_marked(harmonized, harmonized_path)
    gene_scores, evidence, metadata = degora_score_table(harmonized, min_studies=min_studies)
    gene_scores = primary_ranked_scores(gene_scores)
    evidence = _portable_source_path_columns(evidence, source_path_base, path_base)
    version_info = runtime_version_info()
    metadata.update(
        {
            **version_info,
            "path_base": "output_directory",
            "source_path_base": "output_directory",
            "external_path_marker": (
                "external-redacted://<sha256-prefix>/<filename> identifies a non-replayable external path; "
                "replace it with the corresponding local file before replay."
            ),
            "harmonized_path": _portable_cli_path(harmonized_path, path_base),
            "catalog_path": _portable_cli_path(catalog_path, path_base) if catalog_path else "",
            "db_path": _portable_cli_path(db_path, path_base),
        }
    )
    if extra_metadata:
        metadata.update(extra_metadata)
    metadata = sanitize_metadata(metadata, path_base)
    studies = _portable_source_path_columns(
        _active_study_table(catalog_path, harmonized),
        source_path_base,
        path_base,
    )

    score_csv = output_dir / "degora_gene_scores.csv"
    metadata_json = output_dir / "degora_score_metadata.json"
    diagnostics_tsv = output_dir / "degora_source_quality_diagnostics.tsv"
    diagnostics_json = output_dir / "degora_source_quality_diagnostics.json"
    summary_path = output_dir / "degora_score_db_summary.json"
    command_args: list[str | Path | int] = [
        "python3",
        "-m",
        "degora.score_db",
        "--harmonized",
        harmonized_path,
        "--output-dir",
        output_dir,
        "--db",
        db_path,
        "--min-studies",
        min_studies,
    ]
    if catalog_path is not None:
        command_args.extend(["--catalog", catalog_path])
    command = command or shell_command(command_args)
    inputs: list[Path] = [harmonized_path]
    if catalog_path is not None:
        inputs.append(catalog_path)
    sidecar_metadata = {"generator": "degora-score-db", **version_info}

    stored_summary = {
        **version_info,
        "path_base": "output_directory",
        "source_path_base": "output_directory",
        "score_csv": _portable_cli_path(score_csv, path_base),
        "metadata_json": _portable_cli_path(metadata_json, path_base),
        "source_quality_diagnostics_tsv": _portable_cli_path(diagnostics_tsv, path_base),
        "source_quality_diagnostics_json": _portable_cli_path(diagnostics_json, path_base),
        "db_path": _portable_cli_path(db_path, path_base),
        "n_gene_scores": int(len(gene_scores)),
        "n_evidence_rows": int(len(evidence)),
        "n_contrasts": int(len(studies)),
        "n_source_units": int(metadata["n_source_units_total"]),
        "n_source_quality_outliers": int(metadata.get("n_source_quality_outliers", 0)),
        "primary_rank_column": metadata.get("primary_rank_column", PRIMARY_RANK_COLUMN),
        "primary_score_column": metadata.get("primary_score_column", PRIMARY_SCORE_COLUMN),
        "top_genes": gene_scores.head(20)["gene_symbol"].tolist(),
        "near_duplicate_source_unit_warnings": list(
            metadata.get("near_duplicate_source_unit_warnings", [])
        ),
        "within_source_direction_warnings": list(
            metadata.get("within_source_direction_warnings", [])
        ),
        "significance_warnings": _corpus_significance_warnings(gene_scores),
        "direction_conflict_warnings": direction_conflict_warnings(metadata.get("source_quality_diagnostics", [])),
    }

    staging_parent = output_dir.parent
    staging_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output_dir.name}.score-db.", dir=staging_parent) as staging_name:
        staging = Path(staging_name)
        staged_score_csv = staging / score_csv.name
        staged_metadata_json = staging / metadata_json.name
        staged_diagnostics_tsv = staging / diagnostics_tsv.name
        staged_diagnostics_json = staging / diagnostics_json.name
        staged_db_path = staging / "database" / db_path.name
        staged_db_path.parent.mkdir()
        staged_summary_path = staging / summary_path.name

        neutralize_formula_text(gene_scores).to_csv(staged_score_csv, index=False)
        staged_metadata_json.write_text(
            _strict_json_dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        diagnostics = pd.DataFrame.from_records(
            metadata.get("source_quality_diagnostics", []),
            columns=SOURCE_QUALITY_DIAGNOSTIC_COLUMNS,
        )
        neutralize_formula_text(diagnostics).to_csv(staged_diagnostics_tsv, sep="\t", index=False)
        staged_diagnostics_json.write_text(
            _strict_json_dumps(diagnostics.to_dict(orient="records"), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        _write_sqlite(staged_db_path, gene_scores, evidence, studies, metadata)
        staged_summary_path.write_text(
            _strict_json_dumps(stored_summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        artifact_pairs = {
            staged_score_csv: score_csv,
            staged_metadata_json: metadata_json,
            staged_diagnostics_tsv: diagnostics_tsv,
            staged_diagnostics_json: diagnostics_json,
            staged_db_path: db_path,
            staged_summary_path: summary_path,
        }
        publication_pairs = dict(artifact_pairs)
        for staged_artifact, final_artifact in artifact_pairs.items():
            staged_source = artifact_source_path(staged_artifact)
            staged_provenance = artifact_provenance_path(staged_artifact)
            artifact_metadata = dict(sidecar_metadata)
            if final_artifact in {score_csv, diagnostics_tsv}:
                artifact_metadata.update(formula_guard_metadata())
            source_text, json_text = source_sidecar_payloads(
                final_artifact,
                command,
                artifact_content_path=staged_artifact,
                inputs=inputs,
                metadata=artifact_metadata,
            )
            staged_source.write_text(source_text, encoding="utf-8")
            if json_text is not None:
                staged_provenance.write_text(json_text, encoding="utf-8")
            publication_pairs[staged_source] = artifact_source_path(final_artifact)
            publication_pairs[staged_provenance] = artifact_provenance_path(final_artifact)

        publish_staged_artifacts(publication_pairs)
    return {
        **stored_summary,
        "score_csv": str(score_csv.resolve()),
        "metadata_json": str(metadata_json.resolve()),
        "source_quality_diagnostics_tsv": str(diagnostics_tsv.resolve()),
        "source_quality_diagnostics_json": str(diagnostics_json.resolve()),
        "db_path": str(db_path.resolve()),
    }


def main(argv: list[str] | None = None) -> int:
    """Build score artifacts directly from an existing harmonized table."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--db", type=Path)
    parser.add_argument("--min-studies", type=int, default=2)
    args = parser.parse_args(argv)
    summary = write_score_database(
        args.harmonized,
        args.output_dir,
        catalog_path=args.catalog,
        db_path=args.db,
        min_studies=args.min_studies,
    )
    print(_strict_json_dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
