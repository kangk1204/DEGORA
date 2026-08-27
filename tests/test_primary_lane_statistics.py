"""The primary ranking lane must be checkable from what DEGORA writes out.

The unweighted lane has always exported the three statistics its score is built
from. The quality-weighted lane - the one the README calls primary, and the one
the browser orders by - exported none of them, so a reader who saw a positive
``stouffer_z`` beside a ``down`` call had nothing in any output that explained it.
"""

from __future__ import annotations

import sqlite3

import numpy as np
import pandas as pd
import pytest

from degora.score_db import (
    GENE_SCORE_COLUMNS,
    QUALITY_COMPONENT_COLUMNS,
    QUALITY_STATISTIC_COLUMNS,
    write_score_database,
)

PRIMARY_LANE_STATISTICS = ("quality_stouffer_z", "quality_weighted_lfc", "quality_rank_product")


def _harmonized() -> pd.DataFrame:
    """Five source units of unequal quality, harmonized by the real pipeline."""

    from degora.harmonize import TableMapping, harmonize_frame

    rng = np.random.default_rng(7)
    genes = [f"G{index:03d}" for index in range(40)]
    units = [
        ("A", "U1", 6, 6, "author_deg_table"),
        ("B", "U2", 6, 6, "author_deg_table"),
        ("C", "U3", 2, 2, "limma_full_table"),
        ("D", "U4", 3, 3, "derived_count_table"),
        ("E", "U5", 1, 1, "normalized_expression_matrix"),
    ]
    frames = []
    for study, unit, n_ctrl, n_treat, source_input_type in units:
        table = pd.DataFrame(
            {
                "gene": genes,
                "log2FoldChange": rng.normal(0, 1.2, size=len(genes)),
                "pvalue": rng.uniform(1e-9, 1.0, size=len(genes)),
            }
        )
        table["padj"] = np.minimum(table["pvalue"] * 2, 1.0)
        meta = {
            "study_id": study,
            "paper_id": unit,
            "source_unit_id": unit,
            "pipeline": "DESeq2",
            "species": "Homo sapiens",
            "assay_type": "RNA-seq",
            "source_input_type": source_input_type,
            "table_scope": "full_results",
            "n_ctrl": n_ctrl,
            "n_treat": n_treat,
            "source_path": f"deg/{study}.csv",
            "sign_convention": "confirmed_treatment_minus_control",
        }
        frames.append(harmonize_frame(table, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), meta))
    return pd.concat(frames, ignore_index=True)


@pytest.fixture(scope="module")
def scored(tmp_path_factory) -> dict:
    tmp_path = tmp_path_factory.mktemp("primary_lane")
    harmonized = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized, index=False)
    write_score_database(harmonized, tmp_path, db_path=tmp_path / "degora_scores.db", min_studies=2)
    csv = pd.read_csv(tmp_path / "degora_gene_scores.csv")
    with sqlite3.connect(tmp_path / "degora_scores.db") as connection:
        db = pd.read_sql_query("SELECT * FROM genes", connection)
    return {"csv": csv, "db": db}


@pytest.mark.parametrize("column", PRIMARY_LANE_STATISTICS)
def test_the_primary_lane_publishes_the_statistics_it_is_built_from(scored, column: str) -> None:
    assert column in scored["csv"].columns, f"{column} is missing from degora_gene_scores.csv"
    assert column in scored["db"].columns, f"{column} is missing from the score database"


def test_the_primary_direction_can_be_checked_against_its_own_statistic(scored) -> None:
    # The whole point: quality_weighted_consensus_direction is sign(quality_stouffer_z),
    # and a reader can now verify that without reading the source.
    frame = scored["csv"]
    z = pd.to_numeric(frame["quality_stouffer_z"])
    direction = frame["quality_weighted_consensus_direction"].astype(str)
    contradictions = frame.loc[(z.gt(0) & direction.eq("down")) | (z.lt(0) & direction.eq("up")), "gene_symbol"]
    assert contradictions.empty, f"direction contradicts its own statistic for {contradictions.tolist()[:5]}"


def test_the_two_lanes_are_allowed_to_disagree_and_both_stay_explainable(scored) -> None:
    frame = scored["csv"]
    unweighted = pd.to_numeric(frame["stouffer_z"])
    weighted = pd.to_numeric(frame["quality_stouffer_z"])
    for z, column in ((unweighted, "consensus_direction"), (weighted, "quality_weighted_consensus_direction")):
        direction = frame[column].astype(str)
        assert not ((z.gt(0) & direction.eq("down")) | (z.lt(0) & direction.eq("up"))).any(), column


def test_the_weighted_lane_exports_a_counterpart_for_every_unweighted_statistic(scored) -> None:
    columns = set(scored["csv"].columns)
    for unweighted, weighted in (
        ("stouffer_z", "quality_stouffer_z"),
        ("weighted_lfc", "quality_weighted_lfc"),
        ("rank_product", "quality_rank_product"),
    ):
        assert unweighted in columns
        assert weighted in columns, f"{unweighted} is exported but {weighted} is not"


def test_the_new_columns_are_appended_so_no_established_column_moves(scored) -> None:
    columns = list(scored["csv"].columns)
    assert columns[: len(GENE_SCORE_COLUMNS)] == list(GENE_SCORE_COLUMNS)
    tail = columns[len(GENE_SCORE_COLUMNS) :]
    assert set(tail) == set(QUALITY_COMPONENT_COLUMNS) | set(QUALITY_STATISTIC_COLUMNS)
    assert tail[-len(QUALITY_STATISTIC_COLUMNS) :] == list(QUALITY_STATISTIC_COLUMNS)


@pytest.mark.parametrize("column", PRIMARY_LANE_STATISTICS)
def test_the_new_columns_are_numeric_and_finite(scored, column: str) -> None:
    values = pd.to_numeric(scored["csv"][column], errors="coerce")
    assert values.notna().all(), f"{column} has non-numeric entries"
    assert np.isfinite(values.to_numpy()).all(), f"{column} has infinite entries"


def test_the_workbook_documents_every_new_column(tmp_path) -> None:
    import openpyxl

    from degora.excel_export import export_run_workbook

    harmonized = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized, index=False)
    write_score_database(harmonized, tmp_path, db_path=tmp_path / "degora_scores.db", min_studies=2)
    workbook_path = tmp_path / "DEGORA_output.xlsx"
    export_run_workbook(
        tmp_path,
        workbook_path,
        db_path=tmp_path / "degora_scores.db",
        command="pytest",
    )
    workbook = openpyxl.load_workbook(workbook_path)
    dictionary = pd.DataFrame(workbook["Column_dictionary"].values)
    dictionary.columns = dictionary.iloc[0]
    documented = set(dictionary[1:]["column"].astype(str))
    assert set(PRIMARY_LANE_STATISTICS) <= documented
