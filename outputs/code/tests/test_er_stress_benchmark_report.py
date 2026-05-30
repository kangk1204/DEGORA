from __future__ import annotations

import pandas as pd

from scripts.write_er_stress_benchmark_report import build_report


def _summary(path, rows) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def test_er_stress_report_records_primary_and_sensitivity_deltas(tmp_path) -> None:
    primary = tmp_path / "primary.csv"
    full = tmp_path / "full.csv"
    base = {
        "setting_id": "default",
        "run_status": "ok",
        "n_rows": 10,
        "recall_at_10": 0.1,
        "recall_at_20": 0.2,
        "recall_at_50": 0.3,
        "recall_at_100": 0.4,
        "direction_recall_at_50": 0.3,
        "direction_recall_at_100": 0.4,
        "recovered_at_100": "A",
        "direction_recovered_at_100": "A",
        "direction_mismatched_at_100": "",
        "missing_at_100": "B",
        "top10": "A;B",
        "failure_mode": "",
    }
    _summary(
        primary,
        [
            {**base, "method_id": "degora_deg_score", "recall_at_10": 0.4, "recall_at_20": 0.6, "recall_at_50": 0.7},
            {**base, "method_id": "degora_quality_weighted_score", "recall_at_10": 0.5, "recall_at_20": 0.6, "recall_at_50": 0.8},
            {**base, "method_id": "weighted_stouffer", "recall_at_10": 0.3, "recall_at_20": 0.5, "recall_at_50": 0.6},
            {**base, "method_id": "fisher", "recall_at_10": 0.2, "recall_at_20": 0.4, "recall_at_50": 0.5},
        ],
    )
    _summary(
        full,
        [
            {**base, "method_id": "degora_deg_score", "recall_at_10": 0.1},
            {**base, "method_id": "degora_quality_weighted_score", "recall_at_10": 0.3, "recall_at_50": 0.5},
            {**base, "method_id": "fisher", "recall_at_10": 0.5},
        ],
    )

    report = build_report(primary, full)

    assert report["default_degora_method_id"] == "degora_quality_weighted_score"
    assert report["primary_degora_vs_weighted_stouffer_delta"]["recall_at_10"] == 0.2
    assert report["primary_degora_vs_fisher_delta"]["recall_at_50"] == 0.3
    assert report["full_sensitivity_degora_vs_fisher_delta"]["recall_at_10"] == -0.2
    assert report["source_quality_tier_effect"]["degora_recall_primary_minus_full_at_10"] == 0.2
    assert report["quality_weighted_sensitivity_effect"]["quality_weighted_minus_unweighted_score_full_at_10"] == 0.2
    assert "not a basis for choosing the most favorable subset" in report["interpretation"][1]
