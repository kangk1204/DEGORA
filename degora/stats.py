"""Small statistical utilities used by the iteration-1 slice."""

from __future__ import annotations

import numpy as np


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment with NaN-safe output."""

    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    finite = np.isfinite(pvalues)
    if not finite.any():
        return adjusted

    p = np.clip(pvalues[finite], 0.0, 1.0)
    order = np.argsort(p)
    ranked = p[order]
    n = ranked.size
    factors = n / np.arange(1, n + 1)
    ranked_adjusted = np.minimum.accumulate((ranked * factors)[::-1])[::-1]
    ranked_adjusted = np.clip(ranked_adjusted, 0.0, 1.0)

    finite_indices = np.flatnonzero(finite)
    adjusted[finite_indices[order]] = ranked_adjusted
    return adjusted
