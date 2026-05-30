"""Benchmark uncertainty and enrichment utilities for locked gold panels."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import binomtest, hypergeom


@dataclass(frozen=True)
class RecallStats:
    n_gold: int
    n_recovered: int
    recall: float
    ci_low: float
    ci_high: float
    precision: float
    expected_random_recovered: float
    enrichment: float
    hypergeom_pvalue: float


def recall_stats(
    ranked_symbols: list[str],
    positives: set[str],
    k: int,
    *,
    universe_size: int | None = None,
    confidence_level: float = 0.95,
) -> RecallStats:
    """Return recall, exact binomial CI, precision, and top-k enrichment.

    The confidence interval is for the locked positive panel recall. The
    hypergeometric p-value is a list-level enrichment null: among the ranked
    universe, how surprising is at least this many locked positives in top-k.
    """

    n_gold = len(positives)
    universe = int(universe_size or len(ranked_symbols))
    universe = max(universe, len(ranked_symbols), n_gold, 1)
    k_eff = max(0, min(int(k), universe))
    recovered = len(set(ranked_symbols[:k_eff]).intersection(positives))
    if n_gold:
        recall = recovered / n_gold
        ci = binomtest(recovered, n_gold).proportion_ci(
            confidence_level=confidence_level,
            method="exact",
        )
        ci_low = float(ci.low)
        ci_high = float(ci.high)
    else:
        recall = 0.0
        ci_low = 0.0
        ci_high = 0.0
    precision = recovered / k_eff if k_eff else 0.0
    expected = k_eff * n_gold / universe if universe else 0.0
    enrichment = recovered / expected if expected > 0 else np.nan
    if n_gold and k_eff:
        pvalue = float(hypergeom.sf(recovered - 1, universe, n_gold, k_eff))
    else:
        pvalue = 1.0
    return RecallStats(
        n_gold=n_gold,
        n_recovered=recovered,
        recall=float(recall),
        ci_low=ci_low,
        ci_high=ci_high,
        precision=float(precision),
        expected_random_recovered=float(expected),
        enrichment=float(enrichment) if np.isfinite(enrichment) else np.nan,
        hypergeom_pvalue=pvalue,
    )


def background_roc_pr(ranked_symbols: list[str], positives: set[str]) -> dict[str, float]:
    """Use all non-positive ranked genes as background negatives.

    Lower rank is better. AUROC is computed as the probability that a positive
    appears above a background gene. AUPRC is average precision over the ranked
    universe, with prevalence reported so the enrichment over random ranking is
    explicit.
    """

    if not ranked_symbols or not positives:
        return {
            "background_auroc": np.nan,
            "background_auprc": np.nan,
            "background_prevalence": np.nan,
            "background_auprc_enrichment": np.nan,
        }
    seen: set[str] = set()
    unique_symbols: list[str] = []
    for symbol in ranked_symbols:
        if symbol in seen:
            continue
        seen.add(symbol)
        unique_symbols.append(symbol)

    positive_positions = [
        index
        for index, symbol in enumerate(unique_symbols, start=1)
        if symbol in positives
    ]
    n_pos = len(positive_positions)
    n_total = len(unique_symbols)
    n_neg = n_total - n_pos
    auroc_se = np.nan
    auroc_ci_low = np.nan
    auroc_ci_high = np.nan
    if n_pos == 0 or n_neg <= 0:
        auroc = np.nan
    else:
        positives_after = {
            position: n_pos - rank_index
            for rank_index, position in enumerate(positive_positions, start=1)
        }
        # Per-positive DeLong component V10 = fraction of negatives ranked below
        # this positive (lower rank = better). Mean over positives is the AUROC.
        v10 = np.array(
            [(n_total - position - positives_after[position]) / n_neg for position in positive_positions],
            dtype=float,
        )
        auroc = float(v10.mean())
        # Per-negative component V01 = fraction of positives ranked above each
        # negative; one pass tracking cumulative positives by position.
        positive_position_set = set(positive_positions)
        v01_values: list[float] = []
        positives_seen = 0
        for position in range(1, n_total + 1):
            if position in positive_position_set:
                positives_seen += 1
            else:
                v01_values.append(positives_seen / n_pos)
        v01 = np.array(v01_values, dtype=float)
        if n_pos >= 2 and n_neg >= 2:
            var_auroc = v10.var(ddof=1) / n_pos + v01.var(ddof=1) / n_neg
            if np.isfinite(var_auroc) and var_auroc > 0:
                auroc_se = float(np.sqrt(var_auroc))
                auroc_ci_low = float(max(0.0, auroc - 1.96 * auroc_se))
                auroc_ci_high = float(min(1.0, auroc + 1.96 * auroc_se))

    hits = 0
    precisions: list[float] = []
    positive_set = set(positives)
    for index, symbol in enumerate(unique_symbols, start=1):
        if symbol not in positive_set:
            continue
        hits += 1
        precisions.append(hits / index)
    auprc = float(np.sum(precisions) / len(positive_set)) if positive_set else np.nan
    prevalence = len(positive_set.intersection(unique_symbols)) / n_total if n_total else np.nan
    enrichment = auprc / prevalence if np.isfinite(prevalence) and prevalence > 0 else np.nan
    return {
        "background_auroc": float(auroc) if np.isfinite(auroc) else np.nan,
        "background_auroc_se": float(auroc_se) if np.isfinite(auroc_se) else np.nan,
        "background_auroc_ci_low": float(auroc_ci_low) if np.isfinite(auroc_ci_low) else np.nan,
        "background_auroc_ci_high": float(auroc_ci_high) if np.isfinite(auroc_ci_high) else np.nan,
        "background_auprc": float(auprc) if np.isfinite(auprc) else np.nan,
        "background_prevalence": float(prevalence) if np.isfinite(prevalence) else np.nan,
        "background_auprc_enrichment": float(enrichment) if np.isfinite(enrichment) else np.nan,
    }
