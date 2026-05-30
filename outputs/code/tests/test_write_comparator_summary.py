from __future__ import annotations

import pandas as pd

from scripts.write_comparator_summary import SUMMARY_COLUMNS, build_summary


def test_comparator_summary_reports_runnable_and_blocked_methods(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    result = pd.DataFrame(
        {
            "method_id": ["degora_slice"] * 3,
            "setting_id": ["locked"] * 3,
            "gene_id": ["VEGFA", "HK2", "RPL13A"],
            "symbol": ["VEGFA", "HK2", "RPL13A"],
            "rank": [1, 2, 3],
            "score": [10.0, 9.0, 1.0],
            "pvalue": [1e-9, 1e-8, 0.5],
            "padj": [1e-8, 1e-7, 0.8],
            "effect": [2.0, 1.5, 0.0],
            "direction": ["up", "up", "zero"],
            "n_studies": [2, 2, 2],
            "missingness": [0.0, 0.0, 0.0],
            "runtime_s": [0.1, 0.1, 0.1],
            "version": ["test", "test", "test"],
            "status": ["ok", "ok", "ok"],
        }
    )
    result.to_csv(baseline_dir / "hypoxia_degora_slice_locked.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "corpus": ["hypoxia"],
            "method_id": ["awmeta"],
            "setting_id": ["default"],
            "tier": ["tier2"],
            "status": ["open_s1_blocker"],
            "blocker_id": ["awmeta_variance_inputs_missing"],
            "message": ["missing variance"],
            "resolution": ["do not impute"],
        }
    ).to_csv(baseline_dir / "baseline_failure_ledger.csv", index=False)

    summary = build_summary(baseline_dir)

    assert summary.columns.tolist() == SUMMARY_COLUMNS
    assert summary.loc[summary["method_id"].eq("degora_slice"), "recall_at_50"].iloc[0] == 0.1
    awmeta = summary.loc[summary["method_id"].eq("awmeta")].iloc[0]
    assert awmeta["run_status"] == "blocked"
    assert awmeta["failure_mode"] == "awmeta_variance_inputs_missing"
