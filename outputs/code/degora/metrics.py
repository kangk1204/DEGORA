"""Validation metrics for the iteration-1 slice."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd


def recall_at_k(consensus: pd.DataFrame, positives: Iterable[str], k: int) -> dict[str, object]:
    positives_set = {gene.upper() for gene in positives}
    top = consensus.head(k)["gene_symbol"].astype(str).str.upper().tolist()
    recovered = sorted(positives_set.intersection(top))
    return {
        "k": k,
        "n_positives": len(positives_set),
        "n_recovered": len(recovered),
        "recall": len(recovered) / len(positives_set) if positives_set else 0.0,
        "recovered": recovered,
        "missing": sorted(positives_set.difference(top)),
    }
