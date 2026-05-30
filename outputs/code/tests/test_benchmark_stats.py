from __future__ import annotations

from degora.benchmark_stats import background_roc_pr, recall_stats


def test_recall_stats_reports_exact_ci_precision_and_hypergeom() -> None:
    ranked = ["A", "B", "C", "D", "E"]
    positives = {"A", "C", "X", "Y"}

    stats = recall_stats(ranked, positives, 3, universe_size=10)

    assert stats.n_gold == 4
    assert stats.n_recovered == 2
    assert stats.recall == 0.5
    assert 0 <= stats.ci_low < stats.recall < stats.ci_high <= 1
    assert stats.precision == 2 / 3
    assert stats.expected_random_recovered == 1.2
    assert 0 <= stats.hypergeom_pvalue <= 1


def test_background_roc_pr_uses_ranked_nonpositives_as_background() -> None:
    ranked = ["POS1", "NEG1", "POS2", "NEG2", "NEG3"]
    positives = {"POS1", "POS2"}

    metrics = background_roc_pr(ranked, positives)

    assert metrics["background_auroc"] > 0.5
    assert metrics["background_auprc"] > metrics["background_prevalence"]
    assert metrics["background_auprc_enrichment"] > 1.0
