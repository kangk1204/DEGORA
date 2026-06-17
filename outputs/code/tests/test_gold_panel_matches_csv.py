"""Guard against silent drift between the hardcoded HIF1A panel and its CSV.

degora.gold.HIF1A_UP_TARGETS is a hardcoded frozenset, while figures and reproduction
scripts read data/studies/gold/hypoxia_hif1_gold_panel.csv. Both feed published
recall@k for the SAME locked-blind panel, but nothing else enforces that they agree.
A one-sided edit to either source would silently diverge the numbers while still
running. This test fails loudly if they stop matching.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from degora.gold import HIF1A_UP_TARGETS

# tests/ -> code/ -> outputs/ -> repo root
_REPO_ROOT = Path(__file__).resolve().parents[3]
_GOLD_CSV = _REPO_ROOT / "data" / "studies" / "gold" / "hypoxia_hif1_gold_panel.csv"


def test_hif1a_frozenset_matches_locked_gold_csv() -> None:
    if not _GOLD_CSV.exists():
        pytest.skip(f"gold panel CSV not present in this checkout: {_GOLD_CSV}")
    gold = pd.read_csv(_GOLD_CSV)
    assert "gene_symbol" in gold.columns, "gold panel CSV is missing a gene_symbol column"
    csv_symbols = {
        str(symbol).strip().upper()
        for symbol in gold["gene_symbol"].tolist()
        if not pd.isna(symbol) and str(symbol).strip()
    }
    frozen = set(HIF1A_UP_TARGETS)
    assert csv_symbols == frozen, (
        "HIF1A_UP_TARGETS and hypoxia_hif1_gold_panel.csv have diverged. "
        f"only in frozenset: {sorted(frozen - csv_symbols)}; "
        f"only in CSV: {sorted(csv_symbols - frozen)}"
    )
