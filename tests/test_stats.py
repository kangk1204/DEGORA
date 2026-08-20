from __future__ import annotations

import numpy as np

from degora.stats import bh_adjust


def test_bh_adjust_preserves_order_and_nan_slots() -> None:
    adjusted = bh_adjust(np.array([0.01, 0.04, 0.03, np.nan]))

    assert adjusted[0] == 0.03
    assert adjusted[1] == 0.04
    assert adjusted[2] == 0.04
    assert np.isnan(adjusted[3])
