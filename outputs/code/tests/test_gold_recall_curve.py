from __future__ import annotations

import pandas as pd

from figures.make_gold_recall_curve import _baseline_rankings, _k_summary


def test_k_summary_reports_aurc_at_selected_cutoffs() -> None:
    curve = pd.DataFrame(
        {
            "method_label": ["A", "A", "A", "B", "B", "B"],
            "method_id": ["a", "a", "a", "b", "b", "b"],
            "setting_id": ["default"] * 6,
            "k": [1, 2, 3, 1, 2, 3],
            "recall": [0.0, 0.5, 1.0, 0.5, 0.5, 0.5],
            "direction_recall": [0.0, 0.5, 1.0, 0.5, 0.5, 0.5],
            "n_recovered": [0, 1, 2, 1, 1, 1],
            "n_gold": [2] * 6,
        }
    )

    summary = _k_summary(curve, [3])

    method_a = summary.loc[summary["method_label"].eq("A")].iloc[0]
    method_b = summary.loc[summary["method_label"].eq("B")].iloc[0]
    assert method_a["aurc_at_k"] == 0.5
    assert method_a["direction_aurc_at_k"] == 0.5
    assert method_b["aurc_at_k"] == 0.5


def test_baseline_rankings_uses_manifest_to_ignore_stale_baselines(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    current = baseline_dir / "current_weighted_stouffer_default.tsv"
    stale = baseline_dir / "stale_weighted_stouffer_default.tsv"
    current_frame = pd.DataFrame(
        {
            "method_id": ["weighted_stouffer", "weighted_stouffer"],
            "setting_id": ["default", "default"],
            "gene_id": ["HSPA1A", "ACTB"],
            "symbol": ["HSPA1A", "ACTB"],
            "rank": [1, 2],
            "score": [9.0, 1.0],
            "pvalue": [1e-6, 0.7],
            "padj": [1e-4, 0.9],
            "effect": [2.0, 0.0],
            "direction": ["up", "zero"],
            "n_studies": [2, 2],
            "missingness": [0.0, 0.0],
            "runtime_s": [0.1, 0.1],
            "version": ["test", "test"],
            "status": ["ok", "ok"],
        }
    )
    current_frame.to_csv(current, sep="\t", index=False)
    stale_frame = current_frame.copy()
    stale_frame["symbol"] = ["ACTB", "HSPA1A"]
    stale_frame["gene_id"] = stale_frame["symbol"]
    stale_frame.to_csv(stale, sep="\t", index=False)
    pd.DataFrame(
        {
            "artifact": [current.name],
            "method_id": ["weighted_stouffer"],
            "setting_id": ["default"],
            "artifact_type": ["baseline_result"],
            "status": ["ok"],
        }
    ).to_csv(baseline_dir / "baseline_manifest.csv", index=False)

    rankings = _baseline_rankings(baseline_dir)

    assert rankings[("weighted_stouffer", "default")][0] == ("HSPA1A", "up")
