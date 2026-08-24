from __future__ import annotations

from decimal import Decimal
import inspect
import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora import SCORE_VERSION, __version__
import degora.provenance as provenance
import degora.score_db as score_db
from degora.harmonize import TableMapping, harmonize_frame
from degora.score_db import (
    STUDY_TABLE_COLUMNS,
    _format_top_percent,
    _quality_label,
    _quality_label_frame,
    _replicate_quality_multiplier,
    _source_input_type_weight,
    _source_quality_weight_frame,
    _table_scope_multiplier,
    _write_sqlite,
    degora_score_table,
    write_score_database,
)
from degora.provenance import publish_staged_artifacts


def test_portable_cli_path_is_independent_of_package_layout(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    inside = repo_root / "outputs" / "results" / "scores.db"
    outside = tmp_path / "external" / "scores.db"

    assert score_db._portable_cli_path(inside, repo_root) == "outputs/results/scores.db"
    assert score_db._portable_cli_path(outside, repo_root) == "../external/scores.db"


def test_provenance_paths_commands_and_secrets_are_shareable(tmp_path: Path) -> None:
    artifact_dir = tmp_path / "results"
    artifact_dir.mkdir()
    artifact = artifact_dir / "scores.csv"
    artifact.write_text("gene,score\nA,1\n", encoding="utf-8")
    source = tmp_path / "inputs" / "source=one.csv"
    source.parent.mkdir()
    source.write_text("gene\nA\n", encoding="utf-8")
    command = provenance.shell_command(["tool", "--input", source, f"OUTPUT={artifact}"])

    source_text, json_text = provenance.source_sidecar_payloads(
        artifact,
        command,
        inputs=[source],
        metadata={"ncbi_api_key": "do-not-store", "password": "also-secret"},
    )

    assert json_text is not None
    assert str(tmp_path) not in source_text
    assert str(tmp_path) not in json_text
    record = json.loads(json_text)
    assert record["artifact_path"] == "scores.csv"
    assert record["path_base"] == "artifact_directory"
    assert record["inputs"][0]["path"] == "../inputs/source=one.csv"
    assert record["inputs"][0]["path_replayable"] is True
    assert record["command_replayable"] is True
    assert record["metadata"] == {
        "ncbi_api_key": "[redacted]",
        "password": "[redacted]",
    }

    external = tmp_path.parent.parent / "Users" / "researcher" / "source.csv"
    external_command = provenance.shell_command(["tool", "--input", external])
    external_source, external_json = provenance.source_sidecar_payloads(
        artifact,
        external_command,
        inputs=[external],
    )
    assert external_json is not None
    external_record = json.loads(external_json)
    assert "external-redacted://" in external_source
    assert external_record["command_replayable"] is False
    assert external_record["inputs"][0]["path_replayable"] is False
    assert "replace each external-redacted://" in external_record["replay_warning"]
    assert "researcher" not in external_source
    assert "researcher" not in external_json


def test_vectorized_quality_helpers_match_scalar_contract() -> None:
    frame = pd.DataFrame(
        {
            "source_input_type": [
                "author_deg_table",
                "derived_count_table",
                "normalized_expression_matrix",
                "author_deg_table;derived_count_table",
                None,
                "unlisted_type",
            ],
            "table_scope": ["full_results", "full_results", "significant_only", "", None, "other"],
            "n_ctrl": [3, 2, 1, 0, np.nan, 4],
            "n_treat": [3, 2, 1, 0, 2, 4],
        }
    )
    expected = pd.Series(
        [
            _source_input_type_weight(source)
            * _table_scope_multiplier(scope)
            * _replicate_quality_multiplier(ctrl, treat)
            for source, scope, ctrl, treat in frame.itertuples(index=False, name=None)
        ],
        dtype=float,
    ).clip(0.05, 1.0)

    pd.testing.assert_series_equal(_source_quality_weight_frame(frame).reset_index(drop=True), expected)

    labels = pd.Series([np.nan, -np.inf, 0.59, 0.60, 0.849, 0.85, np.inf])
    expected_labels = pd.Series([_quality_label(value) for value in labels], dtype="string")
    pd.testing.assert_series_equal(_quality_label_frame(labels), expected_labels)


def test_write_sqlite_preserves_existing_db_on_failed_rebuild(tmp_path) -> None:
    db = tmp_path / "degora_scores.db"
    genes = pd.DataFrame({"gene_symbol": ["A", "B"], "degora_rank": [1, 2], "degora_score": [0.9, 0.8]})
    evidence = pd.DataFrame({"gene_symbol": ["A", "B"], "study_id": ["S1", "S1"]})
    studies = pd.DataFrame({"source_unit_id": ["S1"]})

    _write_sqlite(db, genes, evidence, studies, {"corpus": "test"})
    assert db.exists()
    original = db.read_bytes()

    # Duplicate gene_symbol violates the unique index, so the rebuild fails mid-write.
    bad_genes = pd.DataFrame({"gene_symbol": ["A", "A"], "degora_rank": [1, 2], "degora_score": [0.9, 0.8]})
    with pytest.raises(sqlite3.IntegrityError):
        _write_sqlite(db, bad_genes, evidence, studies, {"corpus": "test"})

    # The previous good DB must survive, with no lingering temp file.
    assert db.read_bytes() == original
    assert not (tmp_path / "degora_scores.db.tmp").exists()
    with sqlite3.connect(db) as connection:
        assert connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0] == 2


def test_write_sqlite_does_not_follow_predictable_tmp_symlink(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    db = tmp_path / "degora_scores.db"
    predictable_tmp = tmp_path / "degora_scores.db.tmp"
    victim = tmp_path / "victim.db"
    victim.write_text("do not overwrite")
    predictable_tmp.symlink_to(victim)

    real_connect = sqlite3.connect

    def race_predictable_tmp(path: Path, *args, **kwargs):
        candidate = Path(path)
        if candidate == predictable_tmp:
            candidate.unlink(missing_ok=True)
            candidate.symlink_to(victim)
        return real_connect(path, *args, **kwargs)

    monkeypatch.setattr(sqlite3, "connect", race_predictable_tmp)
    genes = pd.DataFrame({"gene_symbol": ["A"], "degora_rank": [1], "degora_score": [0.9]})
    evidence = pd.DataFrame({"gene_symbol": ["A"], "study_id": ["S1"]})
    studies = pd.DataFrame({"source_unit_id": ["S1"]})

    _write_sqlite(db, genes, evidence, studies, {"corpus": "test"})

    assert db.exists()
    assert victim.read_text() == "do not overwrite"
    assert predictable_tmp.is_symlink()


def _harmonized() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3", "S1", "S2", "S3", "S1", "S3"],
            "paper_id": ["P1", "P1", "P2", "P1", "P1", "P2", "P1", "P2"],
            "gene_symbol": ["VEGFA", "VEGFA", "VEGFA", "RPL13A", "RPL13A", "RPL13A", "HK2", "HK2"],
            "lfc": [2.4, 2.0, 1.8, 0.2, -0.2, 0.1, 2.1, 1.9],
            "signed_z": [6.0, 5.0, 4.5, 0.2, -0.2, 0.1, 4.0, 3.8],
            "pvalue": [1e-8, 1e-7, 1e-6, 0.8, 0.7, 0.9, 1e-5, 2e-5],
            "padj": [1e-6, 1e-5, 1e-4, 0.9, 0.8, 0.9, 1e-3, 1e-3],
            "normalized_rank": [0.02, 0.03, 0.04, 0.8, 0.7, 0.9, 0.06, 0.05],
            "n_ctrl": [3, 3, 4, 3, 3, 4, 3, 4],
            "n_treat": [3, 3, 4, 3, 3, 4, 3, 4],
            "n_genes_in_study": [1000] * 8,
            "pipeline": ["DESeq2", "DESeq2", "edgeR", "DESeq2", "DESeq2", "edgeR", "DESeq2", "edgeR"],
            "assay_type": ["RNA-seq", "RNA-seq", "microarray", "RNA-seq", "RNA-seq", "microarray", "RNA-seq", "microarray"],
            "source_input_type": ["author_deg_table"] * 8,
            "platform": [""] * 8,
            "normalization": [""] * 8,
            "probe_collapse": [""] * 8,
            "species": ["Homo sapiens"] * 8,
            "cell_system": ["A", "A", "B", "A", "A", "B", "A", "B"],
            "hypoxia_modality": ["1% O2"] * 8,
            "duration_h": ["24"] * 8,
            "source_path": ["source.csv"] * 8,
            "source_url": ["https://example.test"] * 8,
        }
    )


def _sidecars(path: Path) -> list[Path]:
    return [path.with_suffix(path.suffix + ".source"), path.with_suffix(path.suffix + ".provenance.json")]


def _score_db_artifacts(outdir: Path, db_path: Path) -> list[Path]:
    artifacts = [
        outdir / "degora_gene_scores.csv",
        outdir / "degora_score_metadata.json",
        outdir / "degora_source_quality_diagnostics.tsv",
        outdir / "degora_source_quality_diagnostics.json",
        db_path,
        outdir / "degora_score_db_summary.json",
    ]
    return [path for artifact in artifacts for path in [artifact, *_sidecars(artifact)]]


def _snapshot(paths: list[Path], *, base: Path) -> dict[str, bytes]:
    return {str(path.relative_to(base)): path.read_bytes() for path in paths if path.exists()}


@pytest.mark.parametrize(
    "bad_min_studies",
    [0, -1, "0", "-1", True, False, "two", 2.5, "2.5", Decimal("2.5"), np.nan, Decimal("NaN")],
)
def test_degora_score_table_min_studies_must_be_integer_at_least_one(bad_min_studies) -> None:
    with pytest.raises(ValueError, match="min_studies must be an integer >= 1"):
        degora_score_table(_harmonized(), min_studies=bad_min_studies)


def test_degora_score_table_min_studies_accepts_integer_string_for_public_api_compatibility() -> None:
    scores, _, metadata = degora_score_table(_harmonized(), min_studies="2")

    assert set(scores["gene_symbol"]) == {"VEGFA", "RPL13A", "HK2"}
    assert metadata["min_studies"] == 2


def test_rank_label_and_top_percent_label_track_primary_quality_rank() -> None:
    # C7: rank_label and top_percent_label are the human-readable companions to the
    # primary quality_weighted rank, so they must be keyed to
    # quality_weighted_degora_rank / quality_weighted_top_percent, not the unweighted
    # degora_rank screening lane.
    scores, _, _ = degora_score_table(_harmonized(), min_studies=2)
    fixture_label_rank = (
        scores["rank_label"].astype(str).str.extract(r"#([\d,]+)")[0].str.replace(",", "", regex=False).astype(int)
    )
    assert (fixture_label_rank == scores["quality_weighted_degora_rank"]).all()
    assert (
        scores["top_percent_label"] == scores["quality_weighted_top_percent"].map(_format_top_percent)
    ).all()
    # Lock the wiring explicitly so a revert to the unweighted lane is visible here even
    # though this compact fixture does not force the two rank lanes apart: the label
    # columns must NOT be sourced from degora_rank / top_percent. (A divergence-exercising
    # end-to-end check needs a full real corpus, which is too slow for this unit; the
    # committed-artifact reconciliation covers that case.)
    from degora import score_db as _score_db_module

    source = inspect.getsource(_score_db_module.degora_score_table)
    assert 'scores["rank_label"] = scores["quality_weighted_degora_rank"].map(' in source
    assert 'scores["top_percent_label"] = scores["quality_weighted_top_percent"].map(' in source
    assert 'scores["evidence_tier"] = _evidence_tier(\n        scores["quality_weighted_top_percent"]' in source
    assert 'scores["quality_weighted_sign_concordance"],' in source


def test_degora_score_prioritizes_repeated_directional_source_unit_support() -> None:
    scores, evidence, metadata = degora_score_table(_harmonized(), min_studies=2)

    assert scores.iloc[0]["gene_symbol"] == "VEGFA"
    vegfa = scores.loc[scores["gene_symbol"].eq("VEGFA")].iloc[0]
    hk2 = scores.loc[scores["gene_symbol"].eq("HK2")].iloc[0]
    rpl13a = scores.loc[scores["gene_symbol"].eq("RPL13A")].iloc[0]

    assert vegfa["n_source_units"] == 2
    assert vegfa["n_contrasts_observed"] == 3
    assert vegfa["sign_concordance"] == 1.0
    assert vegfa["rank_label"] == "#1 / 3"
    assert vegfa["top_percent"] == 33.333333
    assert vegfa["percentile"] == 100.0
    assert vegfa["top_percent_label"] == "top 33.33%"
    assert vegfa["support_label"] == "2 / 2 source units"
    assert vegfa["direction_label"] == "100.0% up-concordant"
    assert vegfa["quality_weighted_degora_rank"] >= 1
    assert vegfa["quality_weighted_degora_score"] > 0
    assert vegfa["source_quality_support_score"] > 0
    assert vegfa["priority_score"] > 0
    assert vegfa["priority_rank"] >= 1
    assert vegfa["evidence_reliability_score"] > 0
    assert vegfa["direction_confidence_index"] > 0.5
    assert vegfa["direction_concordant_source_units"] == 2
    assert vegfa["direction_total_source_units"] == 2
    assert vegfa["direction_posterior_mean"] == 0.75
    # With two source units and min_studies=2, every leave-one-source-out fold
    # is ineligible. This is unavailable, not observed zero stability.
    assert vegfa["loo_total_folds"] == 2
    assert vegfa["loo_rank_evaluable_folds"] == 0
    assert vegfa["loo_penalty_folds"] == 2
    assert not bool(vegfa["loo_component_available"])
    assert pd.isna(vegfa["loo_rank_stability_score"])
    assert vegfa["evidence_reliability_components_used"] == 3
    expected_reliability = 100.0 * np.prod(
        [
            vegfa["support_score"],
            vegfa["source_quality_support_score"],
            vegfa["direction_confidence_index"],
        ]
    ) ** (1.0 / 3.0)
    assert np.isclose(vegfa["evidence_reliability_score"], expected_reliability, atol=1e-6)
    assert vegfa["evidence_reliability_score"] > 50
    assert "heterogeneity_i2" in scores.columns
    assert "heterogeneity_flag" in scores.columns
    assert "re_stouffer_z" in scores.columns
    assert "rra_rho" in scores.columns
    assert "effect_meta_log2fc_re" in scores.columns
    assert 0 <= vegfa["heterogeneity_i2"] <= 1
    assert vegfa["heterogeneity_flag"] in {"low_or_unestimated", "moderate_context_review", "high_context_dependent_review"}
    assert vegfa["re_stouffer_shrinkage_factor"] >= 1.0
    assert 0 <= vegfa["rra_rho"] <= 1
    assert vegfa["rra_rank"] >= 1
    assert vegfa["effect_meta_k"] == 2
    assert vegfa["effect_meta_ci_low"] < vegfa["effect_meta_log2fc_re"] < vegfa["effect_meta_ci_high"]
    assert vegfa["degora_score"] > hk2["degora_score"] > rpl13a["degora_score"]
    assert evidence.loc[evidence["gene_symbol"].eq("VEGFA"), "source_unit_id"].nunique() == 2
    assert set(evidence.loc[evidence["gene_symbol"].eq("VEGFA"), "assay_type"]) == {"RNA-seq", "microarray"}
    assert set(evidence["source_quality_label"]) == {"high"}
    assert len(evidence.loc[evidence["gene_symbol"].eq("VEGFA")]) == 2
    p1 = evidence.loc[evidence["gene_symbol"].eq("VEGFA") & evidence["source_unit_id"].eq("P1")].iloc[0]
    assert np.isclose(p1["signed_z"], 5.5)
    assert p1["min_source_padj"] == 1e-6
    assert p1["contributing_study_ids"] == "S1;S2"
    assert set(scores["evidence_tier"]) <= {"A", "B", "C", "D"}
    assert metadata["score_version"] == SCORE_VERSION
    assert metadata["independent_unit_for_consensus"].startswith("source_unit_id")
    assert metadata["n_contrasts_total"] == 3
    assert metadata["n_source_units_total"] == 2
    assert "no max-|z|" in metadata["source_unit_collapse_rule"]
    assert "stouffer_padj" in metadata["high_confidence_rule"]
    assert "evidence_tier_rules" in metadata
    assert "quality_weighted_score_formula" in metadata
    assert "priority_score_weights" in metadata
    assert "evidence_reliability_score_weights" in metadata
    assert "evidence_reliability_score_rule" in metadata
    assert "available diagnostics" in metadata["evidence_reliability_score_rule"]
    assert "direction_confidence_rule" in metadata
    assert "random_effects_stouffer_rule" in metadata
    assert "rra_rule" in metadata
    assert "effect_meta_rule" in metadata
    assert "heterogeneity_rule" in metadata
    assert "heterogeneity_flag_rule" in metadata
    assert "loo_stability_rule" in metadata
    assert "source_quality_diagnostics" in metadata
    rule = metadata["heterogeneity_rule"]
    assert "not a calibrated Higgins' I2" in rule
    # The rule must not claim a bias direction: Q is not chi-square distributed
    # here, raw (Q-df)/Q is frequently negative, and it is clamped to 0.
    assert "no calibrated bias direction" in rule
    assert "positively biased" not in rule
    assert "missing or non-numeric" in metadata["source_quality_weight_rules"]["replicate_multiplier"]


def test_leave_one_out_reports_unavailable_when_every_fold_is_ineligible() -> None:
    scores, _, metadata = degora_score_table(_harmonized(), min_studies=2)

    loo_numeric = [
        "loo_median_rank",
        "loo_rank_iqr",
        "loo_rank_stability_score",
        "loo_top50_fraction",
        "loo_top100_fraction",
    ]
    assert scores[loo_numeric].isna().all().all()
    assert (scores["loo_total_folds"] == 2).all()
    assert (scores["loo_rank_evaluable_folds"] == 0).all()
    assert (scores["loo_penalty_folds"] == 2).all()
    assert (~scores["loo_component_available"]).all()
    assert (scores["evidence_reliability_components_used"] == 3).all()
    assert "unavailable rather than zero" in metadata["loo_stability_rule"]


def test_reliability_rejects_missing_mandatory_component() -> None:
    frame = pd.DataFrame(
        {
            "support_score": [0.9],
            "source_quality_support_score": [np.nan],
            "direction_confidence_index": [0.8],
            "loo_rank_stability_score": [np.nan],
            "loo_component_available": [False],
        }
    )

    with pytest.raises(ValueError, match="source_quality_support_score"):
        score_db._evidence_reliability_score(frame)


def test_loo_fold_scorer_matches_the_full_priority_lane_before_omission() -> None:
    scores, evidence, _ = degora_score_table(_harmonized(), min_studies=2)
    support_denominator = float(np.log1p(evidence["source_unit_id"].nunique()))

    recomputed = score_db._priority_components_from_evidence(
        evidence,
        support_denominator=support_denominator,
        min_studies=2,
    )
    observed = scores[["gene_symbol", "priority_score", "priority_rank", "direction_confidence_index"]].merge(
        recomputed[["gene_symbol", "priority_score", "n_source_units", "direction_confidence_index"]],
        on="gene_symbol",
        suffixes=("_full", "_fold"),
        validate="one_to_one",
    )
    np.testing.assert_allclose(
        observed["priority_score_full"].to_numpy(dtype=float),
        observed["priority_score_fold"].round(6).to_numpy(dtype=float),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        observed["direction_confidence_index_full"].to_numpy(dtype=float),
        observed["direction_confidence_index_fold"].round(6).to_numpy(dtype=float),
        rtol=0.0,
        atol=0.0,
    )
    observed["priority_score_fold"] = observed["priority_score_fold"].round(6)
    observed["direction_confidence_index_fold"] = observed["direction_confidence_index_fold"].round(6)
    fold_ranked = observed.sort_values(
        ["priority_score_fold", "n_source_units", "direction_confidence_index_fold", "gene_symbol"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    fold_ranks = dict(zip(fold_ranked["gene_symbol"], fold_ranked.index + 1, strict=True))
    assert observed["priority_rank"].tolist() == observed["gene_symbol"].map(fold_ranks).tolist()


def test_loo_fold_scorer_uses_neutral_direction_confidence_for_flat_consensus() -> None:
    evidence = pd.DataFrame(
        {
            "gene_symbol": ["GENE", "GENE"],
            "source_unit_id": ["S1", "S2"],
            "signed_z": [1.0, -1.0],
            "lfc": [1.0, -1.0],
            "normalized_rank": [0.1, 0.1],
            "weight": [1.0, 1.0],
        }
    )

    components = score_db._priority_components_from_evidence(
        evidence,
        support_denominator=float(np.log1p(2)),
        min_studies=2,
    )

    assert components.loc[0, "direction_confidence_index"] == 0.5


def test_loo_fold_ranks_use_the_published_six_decimal_tie_contract(monkeypatch) -> None:
    evidence = pd.DataFrame(
        {
            "gene_symbol": ["AAA", "BBB"],
            "source_unit_id": ["S1", "S2"],
        }
    )
    scores = pd.DataFrame(
        {
            "gene_symbol": ["AAA", "BBB"],
            "priority_rank": [1, 2],
        }
    )

    def fake_fold_components(*args, **kwargs):
        return pd.DataFrame(
            {
                "gene_symbol": ["BBB", "AAA"],
                "priority_score": [0.50000049, 0.50000040],
                "n_source_units": [2, 2],
                "direction_confidence_index": [0.5, 0.5],
            }
        )

    monkeypatch.setattr(score_db, "_priority_components_from_evidence", fake_fold_components)
    stability = score_db._leave_one_source_out_stability(
        evidence,
        scores,
        support_denominator=float(np.log1p(2)),
        min_studies=2,
    ).set_index("gene_symbol")

    assert stability.loc["AAA", "loo_median_rank"] == 1.0
    assert stability.loc["BBB", "loo_median_rank"] == 2.0


def test_evaluated_zero_loo_remains_numeric_and_penalizes_reliability(monkeypatch) -> None:
    evidence = pd.DataFrame(
        {
            "gene_symbol": ["AAA", "BBB", "BBB"],
            "source_unit_id": ["S1", "S2", "S3"],
        }
    )
    scores = pd.DataFrame(
        {
            "gene_symbol": ["AAA", "BBB"],
            "priority_rank": [1, 2],
        }
    )

    def fake_fold_components(subset, **kwargs):
        remaining = set(subset["source_unit_id"].astype(str))
        if remaining == {"S2", "S3"}:
            return pd.DataFrame(
                {
                    "gene_symbol": ["AAA", "BBB"],
                    "priority_score": [0.9, 0.8],
                    "n_source_units": [2, 2],
                    "direction_confidence_index": [0.9, 0.8],
                }
            )
        return pd.DataFrame(
            {
                "gene_symbol": ["BBB"],
                "priority_score": [0.8],
                "n_source_units": [2],
                "direction_confidence_index": [0.8],
            }
        )

    monkeypatch.setattr(score_db, "_priority_components_from_evidence", fake_fold_components)
    stability = score_db._leave_one_source_out_stability(
        evidence,
        scores,
        support_denominator=float(np.log1p(3)),
        min_studies=2,
    ).set_index("gene_symbol")

    aaa = stability.loc["AAA"]
    assert aaa["loo_total_folds"] == 3
    assert aaa["loo_rank_evaluable_folds"] == 1
    assert aaa["loo_penalty_folds"] == 2
    assert bool(aaa["loo_component_available"])
    assert aaa["loo_median_rank"] == 1.0
    assert aaa["loo_rank_iqr"] == 0.0
    assert aaa["loo_top50_fraction"] == 1.0
    assert aaa["loo_top100_fraction"] == 1.0
    assert aaa["loo_rank_stability_score"] == 0.0

    reliability_frame = pd.DataFrame(
        {
            "support_score": [0.9],
            "source_quality_support_score": [0.9],
            "direction_confidence_index": [0.9],
            "loo_rank_stability_score": [0.0],
            "loo_component_available": [True],
        }
    )
    reliability, components_used = score_db._evidence_reliability_score(reliability_frame)
    assert components_used.iloc[0] == 4
    assert reliability.iloc[0] < 10.0


def test_rra_exact_null_uses_positive_zero_neglog10() -> None:
    evidence = pd.DataFrame(
        {
            "gene_symbol": ["NULL_GENE", "NULL_GENE"],
            "source_unit_id": ["S1", "S2"],
            "normalized_rank": [1.0, 1.0],
        }
    )

    result = score_db._rra_beta_layer(evidence, total_source_units=2, min_studies=2)
    value = result.loc[0, "rra_neglog10_rho"]

    assert value == 0.0
    assert not np.signbit(value)


def test_deg_only_tables_emit_explicit_noninferential_pvalue_warning() -> None:
    harmonized = _harmonized().assign(table_scope="deg_only")

    _, _, metadata = degora_score_table(harmonized, min_studies=2)

    warning = metadata["stouffer_inference_warning"]
    assert "DEG-only" in warning
    assert "not valid inferential p-values" in warning


def test_primary_direction_label_uses_quality_weighted_lane() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A", "B", "C"],
            "paper_id": ["A", "B", "C"],
            "gene_symbol": ["GENE"] * 3,
            "lfc": [2.0, 2.0, -5.0],
            "signed_z": [4.0, 4.0, -10.0],
            "pvalue": [1e-4, 1e-4, 1e-20],
            "padj": [1e-3, 1e-3, 1e-18],
            "normalized_rank": [0.01, 0.01, 0.001],
            "n_ctrl": [3, 3, 3],
            "n_treat": [3, 3, 3],
            "n_genes_in_study": [1000] * 3,
            "source_input_type": ["author_deg_table", "author_deg_table", "normalized_expression_matrix"],
            "table_scope": ["full_results"] * 3,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    gene = scores.iloc[0]

    assert gene["consensus_direction"] == "down"
    assert gene["quality_weighted_consensus_direction"] == "up"
    assert gene["direction_label"].endswith("up-concordant")


def test_geometric_score_rejects_missing_or_nonfinite_required_components() -> None:
    from degora.score_db import _weighted_geometric_score_with_weights

    frame = pd.DataFrame({"support_score": [1.0], "direction_score": [float("nan")]})
    with pytest.raises(ValueError, match="non-finite score component"):
        _weighted_geometric_score_with_weights(frame, {"support_score": 0.5, "direction_score": 0.5})


def test_degora_score_component_golden_values() -> None:
    scores, _, _ = degora_score_table(_harmonized(), min_studies=2)
    vegfa = scores.loc[scores["gene_symbol"].eq("VEGFA")].iloc[0]

    assert vegfa["support_score"] == pytest.approx(1.0, abs=1e-6)
    assert vegfa["direction_score"] == pytest.approx(1.0, abs=1e-6)
    assert vegfa["evidence_score"] == pytest.approx(0.583257, abs=1e-6)
    assert vegfa["rank_score_component"] == pytest.approx(0.968377, abs=1e-6)
    assert vegfa["effect_score"] == pytest.approx(0.629470, abs=1e-6)
    # Golden value from the deterministic _harmonized() fixture and degora_score_v1_2_source_unit_mean formula.
    assert vegfa["degora_score"] == pytest.approx(85.305316, abs=1e-6)
    assert vegfa["quality_weighted_degora_score"] == pytest.approx(85.298205, abs=1e-6)
    assert vegfa["source_quality_support_score"] == pytest.approx(1.0, abs=1e-6)
    assert vegfa["quality_weighted_sign_concordance"] == pytest.approx(1.0, abs=1e-6)


def test_effect_meta_layer_skips_numerically_unusable_inverse_variance_weights() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A", "B"],
            "paper_id": ["A", "B"],
            "gene_symbol": ["WEAK", "WEAK"],
            "lfc": [0.1, -0.1],
            "signed_z": [1e-154, -1e-154],
            "pvalue": [1.0, 1.0],
            "padj": [1.0, 1.0],
            "normalized_rank": [0.5, 0.5],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    row = scores.loc[scores["gene_symbol"].eq("WEAK")].iloc[0]
    effect_columns = [
        "effect_meta_log2fc_re",
        "effect_meta_se",
        "effect_meta_ci_low",
        "effect_meta_ci_high",
        "effect_meta_tau2",
        "effect_meta_i2",
    ]

    assert row["effect_meta_k"] == 0
    assert not np.isinf(row[effect_columns].astype(float)).any()


def test_degora_score_does_not_treat_same_source_contrasts_as_independent() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["P1_T1", "P1_T2", "P1_T1", "P2_T1"],
            "paper_id": ["P1", "P1", "P1", "P2"],
            "gene_symbol": ["SAME_SOURCE_ONLY", "SAME_SOURCE_ONLY", "CROSS_SOURCE", "CROSS_SOURCE"],
            "lfc": [4.0, 3.5, 1.1, 1.0],
            "signed_z": [8.0, 7.0, 2.0, 2.2],
            "pvalue": [1e-15, 1e-12, 0.01, 0.02],
            "padj": [1e-12, 1e-10, 0.05, 0.05],
            "normalized_rank": [0.001, 0.002, 0.05, 0.06],
            "n_ctrl": [3, 3, 3, 3],
            "n_treat": [3, 3, 3, 3],
            "n_genes_in_study": [1000] * 4,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)

    assert scores["gene_symbol"].tolist() == ["CROSS_SOURCE"]


def test_study_gene_evidence_preserves_mixed_contrast_provenance() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1"],
            "lfc": [1.0, 2.0],
            "signed_z": [2.0, 4.0],
            "pvalue": [0.05, 0.001],
            "padj": [0.1, 0.01],
            "normalized_rank": [0.2, 0.01],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "pipeline": ["DESeq2", "limma_microarray"],
            "assay_type": ["RNA-seq", "microarray"],
            "source_input_type": ["author_deg_table", "normalized_expression_matrix"],
            "platform": ["", "GPLX"],
            "normalization": ["DESeq2", "RMA/log2"],
            "probe_collapse": ["", "min_pvalue_max_abs_lfc"],
            "species": ["human", "human"],
            "cell_system": ["A", "A"],
            "hypoxia_modality": ["drug 1h", "drug 24h"],
            "duration_h": ["1", "24"],
            "source_path": ["a.csv", "b.csv"],
            "source_url": ["u1", "u2"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["study_id"] == "T1"
    assert row["assay_type"] == "RNA-seq;microarray"
    assert row["platform"] == "GPLX"
    assert row["contributing_study_ids"] == "T1;T2"
    assert row["contributing_assay_types"] == "RNA-seq;microarray"
    assert row["contributing_source_paths"] == "a.csv;b.csv"
    assert row["source_quality_weight"] < 1.0
    assert row["source_quality_label"] == "low"
    assert row["min_source_pvalue"] == 0.001
    assert np.isclose(row["aggregate_pvalue"], 2.0 * 0.0013498980316300933)


def test_study_gene_evidence_metadata_follows_time_course_selection() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2", "T3"],
            "paper_id": ["P1", "P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1", "GENE1"],
            "lfc": [1.0, 8.0, 6.0],
            "signed_z": [2.0, 10.0, 9.0],
            "pvalue": [0.04, 1e-12, 1e-9],
            "padj": [0.08, 1e-10, 1e-7],
            "normalized_rank": [0.2, 0.001, 0.002],
            "n_ctrl": [3, 3, 3],
            "n_treat": [3, 3, 3],
            "n_genes_in_study": [1000, 1000, 1000],
            "duration_h": ["1", "8", "24"],
            "time_course_mode": ["early", "early", "early"],
            "source_path": ["early.csv", "mid.csv", "late.csv"],
            "source_url": ["u1", "u2", "u3"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["study_id"] == "T1"
    assert row["contributing_study_ids"] == "T1"
    assert row["contributing_duration_h"] == "1"
    assert row["contributing_source_paths"] == "early.csv"
    assert row["min_source_pvalue"] == 0.04
    assert row["min_source_padj"] == 0.08
    assert row["lfc"] == 1.0
    assert row["signed_z"] == 2.0


def test_source_quality_uses_conservative_replicate_counts_within_source_unit() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1"],
            "lfc": [1.0, 1.5],
            "signed_z": [2.0, 3.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.05, 0.06],
            "normalized_rank": [0.1, 0.2],
            "n_ctrl": [3, 1],
            "n_treat": [3, 1],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["n_ctrl"] == 1.0
    assert row["n_treat"] == 1.0
    assert row["source_quality_weight"] == 0.5
    assert row["source_quality_label"] == "low"


def test_source_quality_uses_conservative_table_scope_within_source_unit() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["T1", "T2"],
            "paper_id": ["P1", "P1"],
            "gene_symbol": ["GENE1", "GENE1"],
            "lfc": [1.0, 1.5],
            "signed_z": [2.0, 3.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.05, 0.06],
            "normalized_rank": [0.1, 0.2],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "deg_only"],
        }
    )

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)
    row = evidence.iloc[0]

    assert row["table_scope"] == "deg_only;full_results"
    assert row["source_quality_weight"] == 0.65
    assert row["source_quality_label"] == "medium"


def test_study_gene_evidence_preserves_verified_sign_convention() -> None:
    harmonized = _harmonized()
    harmonized["sign_convention"] = [
        "inverted_at_source_flipped_at_ingest",
        "treated_vs_control_as_published",
        "treated_vs_control_as_published",
        "treated_vs_control_as_published",
        "treated_vs_control_as_published",
        "treated_vs_control_as_published",
        "treated_vs_control_as_published",
        "treated_vs_control_as_published",
    ]

    _, evidence, _ = degora_score_table(harmonized, min_studies=1)

    assert "sign_convention" in evidence.columns
    assert evidence["sign_convention"].str.contains(
        "inverted_at_source_flipped_at_ingest",
        regex=False,
    ).any()


def test_write_score_database_emits_sqlite_and_sidecars(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        score_db,
        "runtime_version_info",
        lambda: {
            "degora_version": __version__,
            "degora_code_revision": "abc1234-dirty",
            "degora_code_dirty": "true",
        },
    )
    harmonized_path = tmp_path / "harmonized.csv"
    harmonized = _harmonized()
    harmonized["source_path"] = str(tmp_path / "inputs" / "source.csv")
    harmonized.to_csv(harmonized_path, index=False)

    summary = write_score_database(
        harmonized_path,
        tmp_path,
        db_path=tmp_path / "degora_scores.db",
        extra_metadata={
            "derived_count_pilot": True,
            "claim_scope": "unit_test",
            "output_dir": str(tmp_path),
            "ncbi_api_key": "do-not-store",
        },
    )

    db_path = tmp_path / "degora_scores.db"
    assert db_path.exists()
    assert (tmp_path / "degora_gene_scores.csv").exists()
    assert (tmp_path / "degora_gene_scores.csv.source").exists()
    assert (tmp_path / "degora_score_metadata.json").exists()
    assert (tmp_path / "degora_source_quality_diagnostics.tsv").exists()
    assert (tmp_path / "degora_source_quality_diagnostics.tsv.source").exists()
    assert (tmp_path / "degora_source_quality_diagnostics.json").exists()
    assert (tmp_path / "degora_scores.db.source").exists()
    assert summary["n_gene_scores"] == 3
    assert summary["degora_version"] == __version__
    assert summary["degora_code_revision"] == "abc1234-dirty"
    assert summary["degora_code_dirty"] == "true"
    assert summary["score_csv"] == str((tmp_path / "degora_gene_scores.csv").resolve())
    assert summary["db_path"] == str(db_path.resolve())
    metadata_json = json.loads((tmp_path / "degora_score_metadata.json").read_text())
    assert metadata_json["degora_version"] == __version__
    assert metadata_json["degora_code_revision"] == "abc1234-dirty"
    assert metadata_json["degora_code_dirty"] == "true"
    assert metadata_json["primary_rank_column"] == "quality_weighted_degora_rank"
    assert metadata_json["path_base"] == "output_directory"
    assert metadata_json["source_path_base"] == "output_directory"
    assert metadata_json["harmonized_path"] == "harmonized.csv"
    assert metadata_json["db_path"] == "degora_scores.db"
    assert metadata_json["output_dir"] == "."
    assert metadata_json["ncbi_api_key"] == "[redacted]"
    score_provenance = json.loads((tmp_path / "degora_gene_scores.csv.provenance.json").read_text())
    assert score_provenance["artifact_path"] == "degora_gene_scores.csv"
    assert score_provenance["path_base"] == "artifact_directory"
    assert score_provenance["inputs"][0]["path"] == "harmonized.csv"
    assert str(tmp_path) not in json.dumps(score_provenance)
    assert str(tmp_path) not in (tmp_path / "degora_gene_scores.csv.source").read_text()
    stored_summary = json.loads((tmp_path / "degora_score_db_summary.json").read_text())
    assert stored_summary["path_base"] == "output_directory"
    assert stored_summary["source_path_base"] == "output_directory"
    assert stored_summary["score_csv"] == "degora_gene_scores.csv"
    assert stored_summary["db_path"] == "degora_scores.db"
    assert str(tmp_path) not in json.dumps(stored_summary)
    score_csv = pd.read_csv(tmp_path / "degora_gene_scores.csv")
    assert score_csv["quality_weighted_degora_rank"].tolist() == sorted(score_csv["quality_weighted_degora_rank"].tolist())
    assert score_csv["loo_rank_stability_score"].isna().all()
    assert (score_csv["loo_total_folds"] == 2).all()
    assert (score_csv["loo_rank_evaluable_folds"] == 0).all()
    assert (score_csv["loo_penalty_folds"] == 2).all()
    assert (score_csv["evidence_reliability_components_used"] == 3).all()
    assert "NaN" not in (tmp_path / "degora_gene_scores.csv").read_text(encoding="utf-8")
    assert summary["primary_rank_column"] == "quality_weighted_degora_rank"
    assert summary["top_genes"] == score_csv.head(20)["gene_symbol"].tolist()

    with sqlite3.connect(db_path) as connection:
        top_gene = connection.execute("SELECT gene_symbol FROM genes ORDER BY quality_weighted_degora_rank LIMIT 1").fetchone()[0]
        top_label = connection.execute("SELECT top_percent_label FROM genes ORDER BY quality_weighted_degora_rank LIMIT 1").fetchone()[0]
        quality_rank = connection.execute("SELECT quality_weighted_degora_rank FROM genes ORDER BY quality_weighted_degora_rank LIMIT 1").fetchone()[0]
        priority_score = connection.execute("SELECT priority_score FROM genes ORDER BY quality_weighted_degora_rank LIMIT 1").fetchone()[0]
        reliability_score = connection.execute("SELECT evidence_reliability_score FROM genes ORDER BY quality_weighted_degora_rank LIMIT 1").fetchone()[0]
        loo_state = connection.execute(
            "SELECT loo_rank_stability_score, loo_total_folds, loo_rank_evaluable_folds, "
            "loo_penalty_folds, loo_component_available, evidence_reliability_components_used "
            "FROM genes ORDER BY quality_weighted_degora_rank LIMIT 1"
        ).fetchone()
        evidence_rows = connection.execute("SELECT COUNT(*) FROM gene_evidence WHERE gene_symbol = 'VEGFA'").fetchone()[0]
        quality_columns = [row[1] for row in connection.execute("PRAGMA table_info(gene_evidence)").fetchall()]
        metadata = dict(connection.execute("SELECT key, value FROM meta").fetchall())
        stored_source_paths = {
            row[0] for row in connection.execute("SELECT DISTINCT source_path FROM gene_evidence").fetchall()
        }
        study_source_paths = {
            row[0] for row in connection.execute("SELECT DISTINCT source_path FROM studies").fetchall()
        }

    assert top_gene == "VEGFA"
    assert top_label.startswith("top ")
    assert quality_rank >= 1
    assert priority_score > 0
    assert reliability_score > 0
    assert loo_state == (None, 2, 0, 2, 0, 3)
    assert evidence_rows == 2
    assert "source_quality_weight" in quality_columns
    assert "source_recommended_weight" in quality_columns
    assert "source_reliability_weight" in quality_columns
    assert json.loads(metadata["score_weights"])["support_score"] == 0.30
    assert metadata["degora_version"] == __version__
    assert metadata["degora_code_revision"] == "abc1234-dirty"
    assert metadata["degora_code_dirty"] == "true"
    assert stored_source_paths == {"inputs/source.csv"}
    assert study_source_paths == {"inputs/source.csv"}
    assert "priority_score_weights" in metadata
    assert json.loads(metadata["source_quality_weight_rules"])["source_input_type_weights"]["normalized_expression_matrix"] == 0.35
    assert "evidence_tier_rules" in metadata
    assert metadata["derived_count_pilot"] == "true"
    assert metadata["claim_scope"] == "unit_test"
    assert metadata["output_dir"] == "."
    assert metadata["ncbi_api_key"] == "[redacted]"


def test_score_database_rebases_catalog_sources_to_output_directory(tmp_path: Path) -> None:
    config_dir = tmp_path / "config"
    output_dir = tmp_path / "results"
    config_dir.mkdir()
    harmonized_path = config_dir / "harmonized.csv"
    harmonized = _harmonized()
    harmonized["source_path"] = "deg_tables/source.csv"
    harmonized.to_csv(harmonized_path, index=False)
    catalog_path = config_dir / "catalog.csv"
    pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3"],
            "paper_id": ["P1", "P1", "P2"],
            "source_path": ["deg_tables/source.csv"] * 3,
            "gene_column": ["gene_symbol"] * 3,
            "lfc_column": ["lfc"] * 3,
            "p_column": ["pvalue"] * 3,
            "include": ["yes"] * 3,
        }
    ).to_csv(catalog_path, index=False)

    write_score_database(
        harmonized_path,
        output_dir,
        catalog_path=catalog_path,
        db_path=output_dir / "degora_scores.db",
    )

    expected = {"../config/deg_tables/source.csv"}
    with sqlite3.connect(output_dir / "degora_scores.db") as connection:
        evidence_paths = {
            row[0] for row in connection.execute("SELECT DISTINCT source_path FROM gene_evidence").fetchall()
        }
        study_paths = {
            row[0] for row in connection.execute("SELECT DISTINCT source_path FROM studies").fetchall()
        }
        metadata = dict(connection.execute("SELECT key, value FROM meta").fetchall())

    assert evidence_paths == expected
    assert study_paths == expected
    assert metadata["source_path_base"] == "output_directory"


def test_write_score_database_rolls_back_complete_artifact_set_on_late_sqlite_failure(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    db_path = tmp_path / "degora_scores.db"
    write_score_database(harmonized_path, tmp_path, db_path=db_path, extra_metadata={"generation": "old"})
    artifacts = _score_db_artifacts(tmp_path, db_path)
    before = _snapshot(artifacts, base=tmp_path)

    real_write_sqlite = score_db._write_sqlite

    def fail_after_staged_sqlite(*args, **kwargs) -> None:
        real_write_sqlite(*args, **kwargs)
        raise RuntimeError("forced late sqlite failure")

    monkeypatch.setattr(score_db, "_write_sqlite", fail_after_staged_sqlite)

    with pytest.raises(RuntimeError, match="forced late sqlite failure"):
        write_score_database(harmonized_path, tmp_path, db_path=db_path, extra_metadata={"generation": "new"})

    assert _snapshot(artifacts, base=tmp_path) == before


def test_write_score_database_first_build_failure_leaves_no_generated_artifacts(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    outdir = tmp_path / "score_out"
    db_path = outdir / "degora_scores.db"

    def fail_before_sqlite(*_args, **_kwargs) -> None:
        raise RuntimeError("forced sqlite failure")

    monkeypatch.setattr(score_db, "_write_sqlite", fail_before_sqlite)

    with pytest.raises(RuntimeError, match="forced sqlite failure"):
        write_score_database(harmonized_path, outdir, db_path=db_path)

    assert _snapshot(_score_db_artifacts(outdir, db_path), base=tmp_path) == {}


def test_publish_staged_artifacts_rolls_back_known_targets_on_publish_exception(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    final_a = tmp_path / "a.txt"
    final_b = tmp_path / "b.txt"
    final_a.write_text("old a", encoding="utf-8")
    final_b.write_text("old b", encoding="utf-8")
    staged_a = tmp_path / "stage-a.txt"
    staged_b = tmp_path / "stage-b.txt"
    staged_a.write_text("new a", encoding="utf-8")
    staged_b.write_text("new b", encoding="utf-8")
    real_replace = score_db.os.replace
    calls = 0

    def fail_second_replace(src, dst) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("forced publish failure")
        real_replace(src, dst)

    monkeypatch.setattr(score_db.os, "replace", fail_second_replace)

    with pytest.raises(OSError, match="forced publish failure"):
        publish_staged_artifacts({staged_a: final_a, staged_b: final_b})

    assert final_a.read_text(encoding="utf-8") == "old a"
    assert final_b.read_text(encoding="utf-8") == "old b"


def test_publish_staged_artifacts_cleans_pending_files_on_copy_failure(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    final_a = tmp_path / "a.txt"
    final_b = tmp_path / "b.txt"
    final_a.write_text("old a", encoding="utf-8")
    final_b.write_text("old b", encoding="utf-8")
    staged_a = tmp_path / "stage-a.txt"
    staged_b = tmp_path / "stage-b.txt"
    staged_a.write_text("new a", encoding="utf-8")
    staged_b.write_text("new b", encoding="utf-8")
    real_copy2 = provenance.shutil.copy2
    calls = 0

    def fail_second_copy(src, dst, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("forced pending copy failure")
        return real_copy2(src, dst, *args, **kwargs)

    monkeypatch.setattr(provenance.shutil, "copy2", fail_second_copy)

    with pytest.raises(OSError, match="forced pending copy failure"):
        publish_staged_artifacts({staged_a: final_a, staged_b: final_b})

    assert final_a.read_text(encoding="utf-8") == "old a"
    assert final_b.read_text(encoding="utf-8") == "old b"
    assert not list(tmp_path.glob(".*.pending"))
    assert not list(tmp_path.glob(".*.backup"))


def test_degora_score_caps_nonfinite_lfc_for_browsing_layer() -> None:
    harmonized = _harmonized()
    harmonized.loc[harmonized["gene_symbol"].eq("VEGFA"), "lfc"] = np.inf

    scores, evidence, metadata = degora_score_table(harmonized, min_studies=2)

    assert metadata["n_nonfinite_lfc_capped_for_score"] == 3
    assert np.isfinite(scores.loc[scores["gene_symbol"].eq("VEGFA"), "weighted_lfc"].iloc[0])
    assert evidence.loc[evidence["gene_symbol"].eq("VEGFA"), "lfc"].max() == 10.0


def test_quality_weighted_primary_score_downweights_low_quality_discordant_source() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A1", "B1", "C1", "A1", "B1", "C1"],
            "paper_id": ["A", "B", "C", "A", "B", "C"],
            "gene_symbol": ["UPR1", "UPR1", "UPR1", "CTRL", "CTRL", "CTRL"],
            "lfc": [2.0, 2.0, -8.0, 0.2, 0.2, -0.1],
            "signed_z": [4.0, 4.0, -2.0, 0.3, 0.3, -0.2],
            "pvalue": [1e-4, 1e-4, 0.05, 0.8, 0.8, 0.9],
            "padj": [1e-3, 1e-3, 0.2, 0.9, 0.9, 0.9],
            "normalized_rank": [0.01, 0.01, 0.8, 0.8, 0.8, 0.9],
            "n_ctrl": [3, 3, 2, 3, 3, 2],
            "n_treat": [3, 3, 2, 3, 3, 2],
            "n_genes_in_study": [1000] * 6,
            "pipeline": ["author", "derived", "matrix", "author", "derived", "matrix"],
            "assay_type": ["RNA-seq"] * 6,
            "source_input_type": [
                "author_deg_table",
                "derived_count_table",
                "normalized_expression_matrix",
                "author_deg_table",
                "derived_count_table",
                "normalized_expression_matrix",
            ],
            "platform": [""] * 6,
            "normalization": [""] * 6,
            "probe_collapse": [""] * 6,
            "species": ["human"] * 6,
            "cell_system": ["A", "B", "C", "A", "B", "C"],
            "hypoxia_modality": ["ER stress"] * 6,
            "duration_h": ["24"] * 6,
            "table_scope": ["full_results"] * 6,
            "source_path": ["source.csv"] * 6,
            "source_url": ["https://example.test"] * 6,
        }
    )

    scores, evidence, metadata = degora_score_table(harmonized, min_studies=2)
    upr1 = scores.loc[scores["gene_symbol"].eq("UPR1")].iloc[0]

    # UPR1 is the canonical z-vs-LFC discordant gene: two higher-quality sources are
    # strongly up (signed_z +4) while one low-quality matrix source is strongly down in
    # effect size (lfc -8) but only mildly down in z (-2). The weighted combined signed-z
    # is positive, so sign_concordance is measured against "up" (the two up sources are
    # the concordant ones). consensus_direction must report that same combined-z direction
    # ("up"); deriving it from weighted_lfc (negative here) made the direction_label
    # self-contradictory ("X% down-concordant"). Regression guard for the score_db P1
    # direction fix.
    assert upr1["stouffer_z"] > 0
    assert upr1["consensus_direction"] == "up"
    assert upr1["consensus_direction"] == ("up" if upr1["stouffer_z"] > 0 else "down")
    assert upr1["quality_weighted_consensus_direction"] == "up"
    assert upr1["quality_weighted_degora_score"] > 0
    low_quality = evidence.loc[evidence["source_unit_id"].eq("C")].iloc[0]
    assert low_quality["source_quality_label"] == "low"
    assert low_quality["source_quality_weight"] < 0.4
    assert low_quality["source_reliability_weight"] < 0.65
    assert metadata["n_source_quality_outliers"] == 0


def test_direction_confidence_penalizes_discordant_sources() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["A", "B", "A", "B"],
            "paper_id": ["A", "B", "A", "B"],
            "gene_symbol": ["CONSISTENT", "CONSISTENT", "CONFLICT", "CONFLICT"],
            "lfc": [2.0, 2.0, 2.0, -2.0],
            "signed_z": [4.0, 4.0, 4.0, -4.0],
            "pvalue": [1e-4, 1e-4, 1e-4, 1e-4],
            "padj": [1e-3] * 4,
            "normalized_rank": [0.01, 0.01, 0.02, 0.02],
            "n_ctrl": [3] * 4,
            "n_treat": [3] * 4,
            "n_genes_in_study": [1000] * 4,
            "source_input_type": ["author_deg_table"] * 4,
            "table_scope": ["full_results"] * 4,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    consistent = scores.loc[scores["gene_symbol"].eq("CONSISTENT")].iloc[0]
    conflict = scores.loc[scores["gene_symbol"].eq("CONFLICT")].iloc[0]

    assert consistent["direction_confidence_index"] > conflict["direction_confidence_index"]
    assert consistent["evidence_reliability_score"] > conflict["evidence_reliability_score"]


def test_source_coherence_guardrail_flags_low_quality_outlier_source() -> None:
    rows = []
    for source_unit, source_type, lfc_values in [
        ("A", "author_deg_table", [3.0, 2.0, 1.0]),
        ("B", "author_deg_table", [2.5, 1.5, 0.5]),
        ("C", "normalized_expression_matrix", [-3.0, -2.0, -1.0]),
    ]:
        for gene, lfc in zip(["G1", "G2", "G3"], lfc_values, strict=False):
            rows.append(
                {
                    "study_id": source_unit,
                    "paper_id": source_unit,
                    "gene_symbol": gene,
                    "lfc": lfc,
                    "signed_z": lfc,
                    "pvalue": 0.01,
                    "padj": 0.05,
                    "normalized_rank": 0.01,
                    "n_ctrl": 4,
                    "n_treat": 4,
                    "n_genes_in_study": 1000,
                    "source_input_type": source_type,
                    "table_scope": "full_results",
                }
            )
    harmonized = pd.DataFrame(rows)

    _scores, evidence, metadata = degora_score_table(harmonized, min_studies=1)
    source_c = evidence.loc[evidence["source_unit_id"].eq("C")].iloc[0]

    assert metadata["n_source_quality_outliers"] == 1
    assert bool(source_c["source_outlier_flag"]) is True
    assert source_c["source_coherence_weight"] == 0.5


def test_direction_confidence_counts_sources_against_consensus_direction() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["BIG_DOWN", "SMALL_UP1", "SMALL_UP2"],
            "paper_id": ["BIG", "UP1", "UP2"],
            "gene_symbol": ["GENE", "GENE", "GENE"],
            "lfc": [-1.0, 1.0, 1.0],
            "signed_z": [-3.0, 2.0, 2.0],
            "pvalue": [0.001, 0.01, 0.01],
            "padj": [0.01, 0.05, 0.05],
            "normalized_rank": [0.01, 0.02, 0.02],
            "n_ctrl": [60, 1, 1],
            "n_treat": [60, 1, 1],
            "n_genes_in_study": [1000, 1000, 1000],
            "source_input_type": ["author_deg_table"] * 3,
            "table_scope": ["full_results"] * 3,
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    row = scores.iloc[0]

    assert row["stouffer_z"] < 0
    assert row["consensus_direction"] == "down"
    assert row["direction_concordant_source_units"] == 1
    assert row["direction_total_source_units"] == 3
    assert row["direction_confidence_index"] == 0.4


def test_direction_confidence_is_neutral_when_consensus_z_ties() -> None:
    harmonized = pd.DataFrame(
        {
            "study_id": ["UP", "DOWN"],
            "paper_id": ["UP", "DOWN"],
            "gene_symbol": ["GENE", "GENE"],
            "lfc": [1.0, -1.0],
            "signed_z": [2.0, -2.0],
            "pvalue": [0.01, 0.01],
            "padj": [0.05, 0.05],
            "normalized_rank": [0.02, 0.02],
            "n_ctrl": [1, 1],
            "n_treat": [1, 1],
            "n_genes_in_study": [1000, 1000],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
        }
    )

    scores, _, _ = degora_score_table(harmonized, min_studies=2)
    row = scores.iloc[0]

    assert row["stouffer_z"] == pytest.approx(0.0, abs=1e-12)
    assert row["consensus_direction"] == "flat"
    assert row["direction_concordant_source_units"] == 1
    assert row["direction_total_source_units"] == 2
    assert row["direction_confidence_index"] == 0.5
    assert row["direction_posterior_mean"] == 0.5


def test_studies_table_keeps_every_contrast_without_a_catalog(tmp_path) -> None:
    """The fallback studies table must report contrasts, not collapsed evidence rows.

    Evidence carries one row per (gene, source unit) and names the first contributing
    contrast, so a follow-up contrast that covers fewer genes than its sibling in the
    same source unit is never any gene's first. Deriving the table from evidence
    dropped it, and the reported contrast count then disagreed with n_contrasts_total.
    """

    def contrast(study_id: str, unit: str, genes: list[str]) -> pd.DataFrame:
        frame = pd.DataFrame(
            {"gene": genes, "lfc": [2.0] * len(genes), "pvalue": [0.001] * len(genes)}
        )
        return harmonize_frame(
            frame,
            TableMapping("gene", "lfc", "pvalue"),
            {"study_id": study_id, "paper_id": unit, "source_unit_id": unit, "n_ctrl": 3, "n_treat": 3},
        )

    harmonized = pd.concat(
        [
            contrast("u1_a_full", "u1", ["G1", "G2", "G3"]),
            contrast("u1_b_subset", "u1", ["G1"]),
            contrast("u2_c", "u2", ["G1", "G2", "G3"]),
        ],
        ignore_index=True,
    )
    harmonized_path = tmp_path / "harmonized.csv"
    harmonized.to_csv(harmonized_path, index=False)

    summary = write_score_database(harmonized_path, tmp_path / "out", min_studies=2)

    with sqlite3.connect(summary["db_path"]) as connection:
        studies = pd.read_sql_query("SELECT * FROM studies ORDER BY study_id", connection)
        total = connection.execute(
            "SELECT value FROM meta WHERE key = 'n_contrasts_total'"
        ).fetchone()[0]

    assert list(studies["study_id"]) == ["u1_a_full", "u1_b_subset", "u2_c"]
    assert list(studies["source_unit_id"]) == ["u1", "u1", "u2"]
    assert summary["n_contrasts"] == 3
    assert int(total) == 3
    # The catalog and no-catalog branches must stay schema-compatible.
    assert list(studies.columns) == STUDY_TABLE_COLUMNS


def test_the_condition_field_is_published_under_a_topic_neutral_name(tmp_path) -> None:
    """A general-purpose field must not reach the API under one topic's name.

    The catalog's generic `condition` column is stored as `hypoxia_modality`, and
    that name is pinned into the SQLite schema, two API responses, and the shipped
    workbook headers. Both names are emitted so a reader can move to the neutral
    one before the old one is ever removed.
    """

    def contrast(study_id: str, unit: str) -> pd.DataFrame:
        frame = pd.DataFrame({"gene": ["G1", "G2"], "lfc": [2.0, -1.5], "pvalue": [0.001, 0.002]})
        return harmonize_frame(
            frame,
            TableMapping("gene", "lfc", "pvalue"),
            {
                "study_id": study_id,
                "paper_id": unit,
                "source_unit_id": unit,
                "n_ctrl": 3,
                "n_treat": 3,
                "hypoxia_modality": "Alzheimer disease vs control",
            },
        )

    harmonized = pd.concat([contrast("s1", "u1"), contrast("s2", "u2")], ignore_index=True)
    harmonized_path = tmp_path / "harmonized.csv"
    harmonized.to_csv(harmonized_path, index=False)

    summary = write_score_database(harmonized_path, tmp_path / "out", min_studies=2)

    with sqlite3.connect(summary["db_path"]) as connection:
        for table in ("studies", "gene_evidence"):
            columns = {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}
            assert "hypoxia_modality" in columns, table
            assert "condition" in columns, table
        rows = pd.read_sql_query("SELECT hypoxia_modality, condition FROM studies", connection)

    assert (rows["hypoxia_modality"] == rows["condition"]).all()
    assert rows["condition"].iloc[0] == "Alzheimer disease vs control"
