from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora.ablation import run_ablations
from degora.aggregate import (
    collapse_gene_source_units,
    rank_product_consensus,
    source_unit_rows_for_aggregation,
    time_course_selection_report,
)
from degora.harmonize import TableMapping, harmonize_frame
from degora.score_db import (
    ScoreAblation,
    _priority_components_from_evidence,
    _quality_weighted_consensus,
    _rra_beta_layer,
    _weighted_geometric_score_with_weights,
    degora_score_table,
    study_gene_evidence,
    write_score_database,
)


def _one_source_frame(*, normalized_rank: object = 0.1) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1"],
            "source_unit_id": ["U1"],
            "paper_id": ["P1"],
            "gene_symbol": ["G"],
            "lfc": [1.0],
            "signed_z": [2.0],
            "pvalue": [0.0455],
            "padj": [0.05],
            "normalized_rank": [normalized_rank],
            "n_ctrl": [3],
            "n_treat": [3],
            "n_genes_in_study": [100],
            "source_input_type": ["author_deg_table"],
            "table_scope": ["full_results"],
        }
    )


def _strict_load(path: Path) -> object:
    def reject_constant(token: str) -> None:
        raise AssertionError(f"non-standard JSON constant {token!r} in {path}")

    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)


def test_early_selection_is_fixed_on_raw_rows_before_score_eligibility() -> None:
    raw = pd.DataFrame(
        {
            "study_id": ["U1_1h", "U1_24h", "U1_24h"],
            "source_unit_id": ["U1"] * 3,
            "paper_id": ["P1"] * 3,
            "gene_symbol": ["G", "G", "H"],
            "signed_z": [np.nan, 4.0, 3.0],
            "lfc": [0.0, 2.0, 1.0],
            "pvalue": [1.0, 1e-4, 1e-3],
            "padj": [1.0, 1e-3, 1e-2],
            "normalized_rank": [np.nan, 0.01, 0.02],
            "n_genes_in_study": [100] * 3,
            "n_ctrl": [3] * 3,
            "n_treat": [3] * 3,
            "duration_h": ["1", "24", "24"],
            "time_course_mode": ["early"] * 3,
        }
    )

    report = time_course_selection_report(raw)
    selected_for_scoring = source_unit_rows_for_aggregation(raw)
    scores, evidence, _metadata = degora_score_table(
        raw,
        min_studies=1,
        include_loo_stability=False,
    )

    assert report[0]["rows_after"] == 1
    assert report[0]["genes_after"] == 1
    assert selected_for_scoring.empty
    assert collapse_gene_source_units(raw).empty
    assert evidence.empty
    assert scores.empty


def test_peak_mean_is_selected_exactly_once_for_scores_and_audit() -> None:
    raw = pd.DataFrame(
        {
            # Eight rows are intentional: one selection retains four, whereas an
            # accidental second application would retain only two. A four-row
            # fixture cannot detect the regression because the first pass leaves
            # two and the selector deliberately keeps every group of size <= 2.
            "study_id": [f"T{index}" for index in range(1, 9)],
            "source_unit_id": ["U1"] * 8,
            "paper_id": ["P1"] * 8,
            "gene_symbol": ["G"] * 8,
            "signed_z": [float(index) for index in range(1, 9)],
            "lfc": [float(index) for index in range(1, 9)],
            "pvalue": [0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001],
            "padj": [0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001],
            "normalized_rank": [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1],
            "n_genes_in_study": [100] * 8,
            "n_ctrl": [3] * 8,
            "n_treat": [3] * 8,
            "duration_h": [str(index) for index in range(1, 9)],
            "time_course_mode": ["peak_mean"] * 8,
        }
    )

    report = time_course_selection_report(raw)
    scores, evidence, _metadata = degora_score_table(
        raw,
        min_studies=1,
        include_loo_stability=False,
    )

    assert report[0]["rows_after"] == 4
    assert len(scores) == 1
    assert evidence.loc[0, "n_contrast_rows"] == 4
    assert evidence.loc[0, "contributing_study_ids"] == "T5;T6;T7;T8"
    assert evidence.loc[0, "signed_z"] == pytest.approx(6.5)


@pytest.mark.parametrize(
    "bad_rank",
    [0.0, -0.5, np.nan, np.inf, -np.inf, 1.000001, "not-a-number", True],
)
def test_every_public_scoring_path_rejects_invalid_normalized_rank(bad_rank: object) -> None:
    frame = _one_source_frame(normalized_rank=bad_rank)

    with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
        rank_product_consensus(frame, min_studies=1)
    with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
        study_gene_evidence(frame)
    with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
        degora_score_table(frame, min_studies=1, include_loo_stability=False)


@pytest.mark.parametrize(
    ("lfc", "pvalue", "case_id"),
    [
        (2.0, 1.0, "pvalue_one"),
        (0.0, 0.2, "zero_lfc"),
    ],
)
def test_harmonizer_neutral_rows_do_not_block_valid_scores(
    lfc: float,
    pvalue: float,
    case_id: str,
) -> None:
    harmonized = harmonize_frame(
        pd.DataFrame(
            {
                "gene": [f"NEUTRAL_{case_id}", "ACTIVE"],
                "log2FoldChange": [lfc, 1.0],
                "pvalue": [pvalue, 0.01],
            }
        ),
        TableMapping("gene", "log2FoldChange", "pvalue"),
        {"study_id": "S1", "paper_id": "P1"},
    )
    neutral = harmonized.loc[harmonized["gene_symbol"].str.startswith("NEUTRAL_")].iloc[0]
    assert pd.isna(neutral["signed_z"])
    assert pd.isna(neutral["normalized_rank"])

    scores, evidence, _metadata = degora_score_table(
        harmonized,
        min_studies=1,
        include_loo_stability=False,
    )

    assert scores["gene_symbol"].tolist() == ["ACTIVE"]
    assert evidence["gene_symbol"].tolist() == ["ACTIVE"]


def test_finite_signed_z_cannot_hide_a_missing_rank() -> None:
    frame = _one_source_frame(normalized_rank=np.nan)

    with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
        degora_score_table(frame, min_studies=1, include_loo_stability=False)


def test_malformed_non_neutral_signed_z_cannot_hide_a_missing_rank() -> None:
    frame = _one_source_frame(normalized_rank=np.nan)
    frame["signed_z"] = pd.Series("not-a-z", index=frame.index, dtype=object)
    frame.loc[:, "lfc"] = 2.0
    frame.loc[:, "pvalue"] = 0.01

    operations = [
        lambda: rank_product_consensus(frame, min_studies=1),
        lambda: study_gene_evidence(frame),
        lambda: degora_score_table(frame, min_studies=1, include_loo_stability=False),
    ]
    for operation in operations:
        with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
            operation()


@pytest.mark.parametrize("missing_column", ["pvalue", "lfc"])
def test_score_table_rejects_false_replication_outside_evidence_universe(
    missing_column: str,
) -> None:
    first = _one_source_frame()
    second = _one_source_frame()
    second.loc[:, "study_id"] = "S2"
    second.loc[:, "source_unit_id"] = "U2"
    second.loc[:, "paper_id"] = "P2"
    second.loc[:, missing_column] = np.nan
    frame = pd.concat([first, second], ignore_index=True)

    with pytest.raises(ValueError, match="active rows require"):
        degora_score_table(frame, min_studies=2, include_loo_stability=False)


@pytest.mark.parametrize("bad_pvalue", [-0.1, 1.1, True])
def test_score_table_rejects_pvalues_outside_the_harmonized_contract(
    bad_pvalue: object,
) -> None:
    frame = pd.concat([_one_source_frame(), _one_source_frame()], ignore_index=True)
    frame.loc[1, "study_id"] = "S2"
    frame.loc[1, "source_unit_id"] = "U2"
    frame.loc[1, "paper_id"] = "P2"
    frame["pvalue"] = pd.Series([0.01, bad_pvalue], dtype=object)

    with pytest.raises(ValueError, match=r"present pvalue.*\[0, 1\]"):
        degora_score_table(frame, min_studies=2, include_loo_stability=False)


def test_out_of_range_pvalue_cannot_claim_the_neutral_rank_exception() -> None:
    frame = _one_source_frame(normalized_rank=np.nan)
    frame.loc[:, "signed_z"] = np.nan
    frame.loc[:, "pvalue"] = 2.0

    with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
        degora_score_table(frame, min_studies=1, include_loo_stability=False)


@pytest.mark.parametrize(
    ("neutral_field", "missing_field", "neutral_value"),
    [
        ("pvalue", "lfc", 1.0),
        ("lfc", "pvalue", 0.0),
    ],
)
def test_inactive_neutral_audit_row_does_not_require_unscored_evidence_field(
    neutral_field: str,
    missing_field: str,
    neutral_value: float,
) -> None:
    neutral = _one_source_frame(normalized_rank=np.nan)
    neutral.loc[:, "signed_z"] = np.nan
    neutral.loc[:, neutral_field] = neutral_value
    neutral.loc[:, missing_field] = np.nan
    active = _one_source_frame()
    active.loc[:, "gene_symbol"] = "ACTIVE"
    frame = pd.concat([neutral, active], ignore_index=True)

    scores, evidence, _metadata = degora_score_table(
        frame,
        min_studies=1,
        include_loo_stability=False,
    )

    assert scores["gene_symbol"].tolist() == ["ACTIVE"]
    assert evidence["gene_symbol"].tolist() == ["ACTIVE"]


@pytest.mark.parametrize(
    ("neutral_field", "malformed_field", "neutral_value", "malformed_value"),
    [
        ("lfc", "pvalue", 0.0, 2.0),
        ("pvalue", "lfc", 1.0, "not-an-effect"),
    ],
)
def test_inactive_neutral_row_cannot_hide_present_malformed_evidence(
    neutral_field: str,
    malformed_field: str,
    neutral_value: float,
    malformed_value: object,
) -> None:
    frame = _one_source_frame(normalized_rank=np.nan)
    frame.loc[:, "signed_z"] = np.nan
    frame.loc[:, neutral_field] = neutral_value
    frame[malformed_field] = pd.Series([malformed_value], dtype=object)

    with pytest.raises(ValueError, match="every present lfc"):
        degora_score_table(frame, min_studies=1, include_loo_stability=False)


@pytest.mark.parametrize("bad_rank", [0.0, -0.5, np.inf, -np.inf, 1.000001, "not-a-number", True])
def test_score_ineligible_rows_still_reject_present_malformed_ranks(bad_rank: object) -> None:
    frame = _one_source_frame(normalized_rank=bad_rank)
    frame.loc[:, "signed_z"] = np.nan

    operations = [
        lambda: rank_product_consensus(frame, min_studies=1),
        lambda: study_gene_evidence(frame),
        lambda: degora_score_table(frame, min_studies=1, include_loo_stability=False),
    ]
    for operation in operations:
        with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
            operation()


def test_private_rank_lanes_also_fail_closed_on_invalid_rank() -> None:
    evidence = _one_source_frame(normalized_rank=-0.5)

    operations = [
        lambda: _quality_weighted_consensus(evidence, total_source_quality_weight=1.0),
        lambda: _priority_components_from_evidence(
            evidence,
            support_denominator=float(np.log1p(1)),
            min_studies=1,
        ),
        lambda: _rra_beta_layer(evidence, total_source_units=1, min_studies=1),
    ]
    for operation in operations:
        with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
            operation()


def test_score_database_rejects_invalid_rank_before_writing_outputs(tmp_path: Path) -> None:
    harmonized = tmp_path / "invalid_rank.csv"
    _one_source_frame(normalized_rank=-0.5).to_csv(harmonized, index=False)

    with pytest.raises(ValueError, match=r"0 < normalized_rank <= 1"):
        write_score_database(harmonized, tmp_path / "scores", min_studies=1)


def test_score_diagnostics_are_strict_json_in_memory_files_and_database(tmp_path: Path) -> None:
    frame = _one_source_frame()
    _scores, _evidence, metadata = degora_score_table(
        frame,
        min_studies=1,
        include_loo_stability=False,
    )
    diagnostics = metadata["source_quality_diagnostics"][0]
    undefined_fields = {
        "median_pairwise_lfc_spearman",
        "min_pairwise_lfc_spearman",
        "median_pairwise_sign_agreement",
    }

    assert {name for name in undefined_fields if diagnostics[name] is None} == undefined_fields
    json.dumps(metadata, allow_nan=False)

    harmonized = tmp_path / "harmonized.csv"
    frame.to_csv(harmonized, index=False)
    summary = write_score_database(harmonized, tmp_path / "scores", min_studies=1)
    metadata_file = Path(summary["metadata_json"])
    diagnostics_file = Path(summary["source_quality_diagnostics_json"])
    parsed_metadata = _strict_load(metadata_file)
    parsed_diagnostics = _strict_load(diagnostics_file)

    assert "NaN" not in metadata_file.read_text(encoding="utf-8")
    assert "NaN" not in diagnostics_file.read_text(encoding="utf-8")
    assert parsed_metadata["source_quality_diagnostics"][0]["median_pairwise_lfc_spearman"] is None
    assert parsed_diagnostics[0]["median_pairwise_sign_agreement"] is None

    with sqlite3.connect(summary["db_path"]) as connection:
        stored = connection.execute(
            "SELECT value FROM meta WHERE key = 'source_quality_diagnostics'"
        ).fetchone()[0]
    parsed_stored = json.loads(stored, parse_constant=lambda token: pytest.fail(token))
    assert parsed_stored[0]["min_pairwise_lfc_spearman"] is None


def test_geometric_score_treats_zero_as_absorbing_evidence() -> None:
    components = pd.DataFrame({"support_score": [1.0], "direction_score": [0.0]})
    score = _weighted_geometric_score_with_weights(
        components,
        {"support_score": 0.5, "direction_score": 0.5},
    )

    assert score.iloc[0] == 0.0

    for invalid in (-0.1, 1.1):
        with pytest.raises(ValueError, match=r"0 <= component <= 1"):
            _weighted_geometric_score_with_weights(
                pd.DataFrame({"support_score": [invalid]}),
                {"support_score": 1.0},
            )


def test_ablation_names_are_nonblank_unique_and_reserve_canonical_full() -> None:
    with pytest.raises(ValueError, match="non-blank"):
        ScoreAblation(name="  ")
    with pytest.raises(ValueError, match="reserved"):
        ScoreAblation(name="full", component_weights={"support_score": 1.0})

    unnamed_custom = ScoreAblation(component_weights={"support_score": 1.0})
    assert unnamed_custom.name == "custom"

    duplicate = [
        ScoreAblation(name="same"),
        ScoreAblation(name="same", disable_sample_size_weighting=True),
    ]
    with pytest.raises(ValueError, match="must be unique"):
        run_ablations(pd.DataFrame(), ablations=duplicate)


def test_one_source_support_is_one_in_both_scoring_lanes() -> None:
    scores, _evidence, metadata = degora_score_table(
        _one_source_frame(),
        min_studies=1,
        include_loo_stability=False,
    )
    row = scores.iloc[0]

    assert row["support_score"] == 1.0
    assert row["source_quality_support_score"] == 1.0
    assert "one-source corpus has 1.0 in both support lanes" in metadata["support_normalization_rule"]
