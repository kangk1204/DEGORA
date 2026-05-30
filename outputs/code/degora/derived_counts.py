"""Helpers for DEG-like tables derived from public count matrices."""

from __future__ import annotations

from typing import Any

import pandas as pd


LOW_COUNT_FILTER_MIN_COUNT = 10.0
LOW_COUNT_FILTER_MIN_SAMPLES = 2


def low_count_filter_mask(
    counts: pd.DataFrame,
    *,
    min_count: float = LOW_COUNT_FILTER_MIN_COUNT,
    min_samples: int = LOW_COUNT_FILTER_MIN_SAMPLES,
) -> pd.Series:
    """Return an expression-presence mask independent of the treatment label."""

    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return numeric.ge(float(min_count)).sum(axis=1).ge(int(min_samples))


def low_count_filter_summary(
    counts: pd.DataFrame,
    mask: pd.Series,
    *,
    min_count: float = LOW_COUNT_FILTER_MIN_COUNT,
    min_samples: int = LOW_COUNT_FILTER_MIN_SAMPLES,
) -> dict[str, Any]:
    """Return JSON-safe metadata for the count-derived expression filter."""

    before = int(len(counts))
    after = int(mask.sum())
    return {
        "low_count_filter": "raw count >= min_count in at least min_samples selected samples",
        "low_count_filter_min_count": float(min_count),
        "low_count_filter_min_samples": int(min_samples),
        "n_genes_before_low_count_filter": before,
        "n_genes_after_low_count_filter": after,
        "n_genes_removed_low_count_filter": int(before - after),
    }


def attach_low_count_filter_metadata(
    frame: pd.DataFrame,
    summary: dict[str, Any],
) -> pd.DataFrame:
    """Attach filter metadata columns to a derived DEG-like table."""

    out = frame.copy()
    for key, value in summary.items():
        out[key] = value
    return out
