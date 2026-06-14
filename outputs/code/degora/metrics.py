"""Validation metrics for the iteration-1 slice."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd


def recall_at_k(consensus: pd.DataFrame, positives: Iterable[str], k: int) -> dict[str, object]:
    if "gene_symbol" not in consensus.columns:
        raise ValueError("consensus must have a 'gene_symbol' column")
    k_eff = max(0, int(k))
    # Normalize both sides identically (strip + upper, drop NaN/blank, coerce non-strings)
    # so a whitespace/case-dirty or non-string gold panel matches instead of silently
    # under-counting recall or crashing on .upper().
    positives_set = {
        str(gene).strip().upper()
        for gene in positives
        if not pd.isna(gene) and str(gene).strip()
    }
    ranked = consensus["gene_symbol"].astype(str).str.strip().str.upper().tolist()
    top: list[str] = []
    seen: set[str] = set()
    if k_eff:
        for gene in ranked:
            if not gene or gene in seen:
                continue
            seen.add(gene)
            top.append(gene)
            if len(top) >= k_eff:
                break
    recovered = sorted(positives_set.intersection(top))
    return {
        "k": k,
        "n_positives": len(positives_set),
        "n_recovered": len(recovered),
        "recall": len(recovered) / len(positives_set) if positives_set else 0.0,
        "recovered": recovered,
        "missing": sorted(positives_set.difference(top)),
    }
