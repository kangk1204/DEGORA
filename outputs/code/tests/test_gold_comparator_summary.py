from __future__ import annotations

import pandas as pd
import pytest

from scripts.write_gold_comparator_summary import SUMMARY_COLUMNS, build_summary


def test_gold_comparator_summary_includes_degora_score_and_blockers(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    pd.DataFrame(
        {
            "method_id": ["weighted_stouffer"] * 3,
            "setting_id": ["default"] * 3,
            "gene_id": ["ISG15", "RPL13A", "IFIT1"],
            "symbol": ["ISG15", "RPL13A", "IFIT1"],
            "rank": [1, 2, 3],
            "score": [9.0, 1.0, 8.0],
            "pvalue": [1e-6, 0.7, 1e-5],
            "padj": [1e-4, 0.9, 1e-3],
            "effect": [2.0, 0.0, 1.5],
            "direction": ["up", "zero", "up"],
            "n_studies": [2, 2, 2],
            "missingness": [0.0, 0.0, 0.0],
            "runtime_s": [0.1, 0.1, 0.1],
            "version": ["test", "test", "test"],
            "status": ["ok", "ok", "ok"],
        }
    ).to_csv(baseline_dir / "ifn_weighted_stouffer_default.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "corpus": ["ifn"],
            "method_id": ["hstouffer"],
            "setting_id": ["default"],
            "tier": ["tier2"],
            "status": ["open_s1_blocker"],
            "blocker_id": ["hstouffer_deg_table_materializer_blocked"],
            "message": ["blocked"],
            "resolution": ["obtain input"],
        }
    ).to_csv(baseline_dir / "baseline_failure_ledger.csv", index=False)
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["ISG15", "IFIT1"]}).to_csv(gold_path, index=False)
    score_path = tmp_path / "degora_gene_scores.csv"
    pd.DataFrame({"degora_rank": [1, 2, 3], "gene_symbol": ["IFIT1", "ISG15", "RPL13A"]}).to_csv(score_path, index=False)

    summary = build_summary(baseline_dir, gold_path, degora_score_csv=score_path)

    assert summary.columns.tolist() == SUMMARY_COLUMNS
    degora = summary.loc[summary["method_id"].eq("degora_deg_score")].iloc[0]
    assert degora["recall_at_10"] == 1.0
    assert degora["n_recovered_at_10"] == 2
    assert degora["precision_at_10"] == 2 / 3
    assert 0 <= degora["hypergeom_fdr_at_10"] <= 1
    blocked = summary.loc[summary["method_id"].eq("hstouffer")].iloc[0]
    assert blocked["run_status"] == "blocked"


def test_gold_comparator_summary_counts_expected_direction_when_available(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    pd.DataFrame(
        {
            "method_id": ["weighted_stouffer"] * 3,
            "setting_id": ["default"] * 3,
            "gene_id": ["HSPA5", "DDIT3", "ACTB"],
            "symbol": ["HSPA5", "DDIT3", "ACTB"],
            "rank": [1, 2, 3],
            "score": [9.0, 8.0, 1.0],
            "pvalue": [1e-6, 1e-5, 0.7],
            "padj": [1e-4, 1e-3, 0.9],
            "effect": [2.0, -1.5, 0.0],
            "direction": ["up", "down", "zero"],
            "n_studies": [2, 2, 2],
            "missingness": [0.0, 0.0, 0.0],
            "runtime_s": [0.1, 0.1, 0.1],
            "version": ["test", "test", "test"],
            "status": ["ok", "ok", "ok"],
        }
    ).to_csv(baseline_dir / "er_weighted_stouffer_default.tsv", sep="\t", index=False)
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["HSPA5", "DDIT3"], "expected_direction": ["up", "up"]}).to_csv(
        gold_path,
        index=False,
    )
    score_path = tmp_path / "degora_gene_scores.csv"
    pd.DataFrame(
        {
            "degora_rank": [1, 2, 3],
            "gene_symbol": ["DDIT3", "HSPA5", "ACTB"],
            "consensus_direction": ["up", "up", "zero"],
        }
    ).to_csv(score_path, index=False)

    summary = build_summary(baseline_dir, gold_path, degora_score_csv=score_path)

    degora = summary.loc[summary["method_id"].eq("degora_deg_score")].iloc[0]
    assert degora["direction_recall_at_10"] == 1.0
    weighted = summary.loc[summary["method_id"].eq("weighted_stouffer")].iloc[0]
    assert weighted["recall_at_10"] == 1.0
    assert weighted["direction_recall_at_10"] == 0.5
    assert weighted["direction_mismatched_at_100"] == "DDIT3:down!=up"


def test_gold_comparator_summary_adds_quality_weighted_degora_row(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["HSPA5", "DDIT3"], "expected_direction": ["up", "up"]}).to_csv(
        gold_path,
        index=False,
    )
    score_path = tmp_path / "degora_gene_scores.csv"
    pd.DataFrame(
        {
            "degora_rank": [1, 2, 3],
            "quality_weighted_degora_rank": [2, 1, 3],
            "gene_symbol": ["ACTB", "DDIT3", "HSPA5"],
            "consensus_direction": ["up", "down", "up"],
            "quality_weighted_consensus_direction": ["up", "up", "up"],
        }
    ).to_csv(score_path, index=False)

    summary = build_summary(baseline_dir, gold_path, degora_score_csv=score_path)

    assert summary["method_id"].tolist() == ["degora_deg_score", "degora_quality_weighted_score"]
    primary = summary.loc[summary["method_id"].eq("degora_deg_score")].iloc[0]
    secondary = summary.loc[summary["method_id"].eq("degora_quality_weighted_score")].iloc[0]
    assert primary["direction_recall_at_10"] == 0.5
    assert secondary["direction_recall_at_10"] == 1.0


def test_gold_comparator_summary_uses_manifest_to_ignore_stale_baselines(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    current = baseline_dir / "current_weighted_stouffer_default.tsv"
    stale = baseline_dir / "stale_weighted_stouffer_default.tsv"
    columns = {
        "method_id": ["weighted_stouffer"] * 3,
        "setting_id": ["default"] * 3,
        "gene_id": ["ISG15", "IFIT1", "ACTB"],
        "symbol": ["ISG15", "IFIT1", "ACTB"],
        "rank": [1, 2, 3],
        "score": [9.0, 8.0, 1.0],
        "pvalue": [1e-6, 1e-5, 0.7],
        "padj": [1e-4, 1e-3, 0.9],
        "effect": [2.0, 1.5, 0.0],
        "direction": ["up", "up", "zero"],
        "n_studies": [2, 2, 2],
        "missingness": [0.0, 0.0, 0.0],
        "runtime_s": [0.1, 0.1, 0.1],
        "version": ["test", "test", "test"],
        "status": ["ok", "ok", "ok"],
    }
    pd.DataFrame(columns).to_csv(current, sep="\t", index=False)
    stale_columns = {**columns, "gene_id": ["ACTB"], "symbol": ["ACTB"], "rank": [1]}
    for key in ("method_id", "setting_id", "score", "pvalue", "padj", "effect", "direction", "n_studies", "missingness", "runtime_s", "version", "status"):
        stale_columns[key] = stale_columns[key][:1] if isinstance(stale_columns[key], list) else stale_columns[key]
    pd.DataFrame(stale_columns).to_csv(stale, sep="\t", index=False)
    pd.DataFrame(
        {
            "artifact": [current.name],
            "method_id": ["weighted_stouffer"],
            "setting_id": ["default"],
            "artifact_type": ["baseline_result"],
            "status": ["ok"],
        }
    ).to_csv(baseline_dir / "baseline_manifest.csv", index=False)
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["ISG15", "IFIT1"]}).to_csv(gold_path, index=False)

    summary = build_summary(baseline_dir, gold_path)

    weighted = summary.loc[summary["method_id"].eq("weighted_stouffer")].iloc[0]
    assert weighted["n_rows"] == 3
    assert weighted["recall_at_10"] == 1.0
    assert weighted["top10"] == "ISG15;IFIT1;ACTB"


def test_gold_comparator_summary_rejects_empty_gold_panel(tmp_path) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    gold_path = tmp_path / "gold.csv"
    pd.DataFrame({"gene_symbol": ["", None]}).to_csv(gold_path, index=False)

    with pytest.raises(ValueError, match="no positive genes"):
        build_summary(baseline_dir, gold_path)
