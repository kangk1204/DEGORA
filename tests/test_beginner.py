from __future__ import annotations

import pandas as pd
import pytest

from degora.beginner import (
    build_catalog,
    catalog_row,
    default_study_id,
    find_source_tables,
    infer_source_table,
    run_init,
)
from degora.slice_runner import validate_catalog_inputs


def _clean_table(path):
    pd.DataFrame(
        {
            "gene": [f"G{index}" for index in range(200)],
            "log2FoldChange": [2.0 - index * 0.01 for index in range(200)],
            "pvalue": [1e-5 if index < 50 else 0.5 for index in range(200)],
            "padj": [1e-4 if index < 50 else 0.9 for index in range(200)],
        }
    ).to_csv(path, index=False)


def _ambiguous_table(path):
    """Two plausible effect columns and no log2 in either name."""

    pd.DataFrame(
        {
            "Gene symbol": [f"G{index}" for index in range(60)],
            "FC": [1.5] * 60,
            "ratio": [1.2] * 60,
            "P.Value": [0.01] * 60,
            "adj.P.Val": [0.05] * 60,
        }
    ).to_csv(path, sep="\t", index=False)


def _not_a_deg_table(path):
    pd.DataFrame({"sample": ["A", "B"], "batch": [1, 2]}).to_csv(path, index=False)


def test_a_clean_table_needs_no_questions(tmp_path) -> None:
    """Everything a DESeq2 export states about itself is read off it."""

    path = tmp_path / "clean_deseq2.csv"
    _clean_table(path)

    inference = infer_source_table(path)

    assert inference.readable
    assert inference.looks_like_a_deg_table
    assert inference.needs_a_question == ()
    assert inference.mapping == {
        "gene_column": "gene",
        "lfc_column": "log2FoldChange",
        "p_column": "pvalue",
        "padj_column": "padj",
    }
    # The scope is inferred too, so the reader is not asked to classify their table.
    assert inference.table_scope == "full_results"
    assert inference.n_rows == 200


def test_a_missing_optional_column_is_not_a_question(tmp_path) -> None:
    """DEGORA runs without an adjusted p-value; asking for an absent one is noise."""

    path = tmp_path / "no_padj.csv"
    pd.DataFrame(
        {"Symbol": [f"G{i}" for i in range(40)], "logFC": [2.0] * 40, "PValue": [1e-6] * 40}
    ).to_csv(path, index=False)

    inference = infer_source_table(path)

    # `mapping` reports the columns that were found, so an absent optional one is
    # simply not a key - and it is not a question either.
    assert "padj_column" not in inference.mapping
    assert inference.mapping["p_column"] == "PValue"
    assert [choice.role for choice in inference.needs_a_question] == []


def test_an_ambiguous_effect_column_is_a_question_not_a_guess(tmp_path) -> None:
    """Picking between two effect columns is the reader's call, not a coin toss."""

    path = tmp_path / "ambiguous.tsv"
    _ambiguous_table(path)

    inference = infer_source_table(path)

    assert inference.looks_like_a_deg_table
    assert "lfc_column" in {choice.role for choice in inference.needs_a_question}


def test_a_file_that_is_not_a_deg_table_is_recognised(tmp_path) -> None:
    """A sample sheet in the same folder must not become a question-by-question walk."""

    path = tmp_path / "sample_sheet.csv"
    _not_a_deg_table(path)

    assert not infer_source_table(path).looks_like_a_deg_table


def test_an_unreadable_file_does_not_end_the_walk(tmp_path) -> None:
    path = tmp_path / "broken.xlsx"
    path.write_text("not a workbook", encoding="utf-8")

    inference = infer_source_table(path)

    assert not inference.readable
    assert inference.problem
    assert not inference.looks_like_a_deg_table


def test_find_source_tables_skips_dotfiles(tmp_path) -> None:
    _clean_table(tmp_path / "real.csv")
    (tmp_path / ".hidden.csv").write_text("gene,lfc\n", encoding="utf-8")

    assert [path.name for path in find_source_tables(tmp_path)] == ["real.csv"]


def test_study_ids_are_readable_and_unique(tmp_path) -> None:
    from pathlib import Path

    assert default_study_id(Path("GSE123 hypoxia (4h).csv"), []) == "GSE123_hypoxia_4h"
    assert default_study_id(Path("a.csv"), ["a"]) == "a_2"
    assert default_study_id(Path("a.csv"), ["a", "a_2"]) == "a_3"


def test_a_reversed_direction_is_excluded_rather_than_flipped(tmp_path) -> None:
    """DEGORA never reverses an effect column, so it must not pretend to here.

    A table whose positive values mean "up in the control group" is unusable as it
    stands. Flipping the sign for the reader would put a correction DEGORA cannot
    verify into the middle of the evidence chain, so the row is written excluded
    with the reason on it.
    """

    from degora.beginner import ContrastAnswers

    path = tmp_path / "reversed.csv"
    _clean_table(path)
    inference = infer_source_table(path)

    row = catalog_row(
        inference,
        ContrastAnswers(positive_means_up_in_treated=False, source_unit_id="P1"),
        study_id="S1",
        catalog_dir=tmp_path,
    )

    assert row["include_in_analysis"] == "no"
    assert "REVERSED" in row["sign_convention"]
    assert "correct the table at source" in row["notes"]
    # The mapping is preserved, so fixing the table is the only remaining step.
    assert row["lfc_column"] == "log2FoldChange"


def test_a_confirmed_direction_is_recorded_not_assumed(tmp_path) -> None:
    from degora.beginner import ContrastAnswers

    path = tmp_path / "forward.csv"
    _clean_table(path)

    row = catalog_row(
        infer_source_table(path),
        ContrastAnswers(positive_means_up_in_treated=True, source_unit_id="P1"),
        study_id="S1",
        catalog_dir=tmp_path,
    )

    assert row["include_in_analysis"] == "yes"
    assert row["sign_convention"] == "confirmed_treatment_minus_control"


def test_the_guided_flow_produces_a_config_that_validates(tmp_path) -> None:
    """The whole point: what init writes has to be something run accepts."""

    deg = tmp_path / "deg"
    deg.mkdir()
    _clean_table(deg / "clean_deseq2.csv")
    _ambiguous_table(deg / "ambiguous.tsv")
    _not_a_deg_table(deg / "sample_sheet.csv")

    answers = iter(
        [
            "human",  # species, asked once
            # ambiguous.tsv - the reader picks the effect column, then the direction
            "FC",
            "yes",
            "drug vs vehicle",
            "PAPER_A",
            "3",
            "3",
            # clean_deseq2.csv
            "yes",
            "hypoxia vs normoxia",
            "PAPER_B",
            "4",
            "4",
        ]
    )
    config = tmp_path / "config.csv"

    summary = run_init(config, deg, ask=lambda question, default="": next(answers), echo=lambda _line: None)

    assert summary["n_contrasts"] == 2
    assert summary["n_source_units"] == 2
    assert [entry["path"] for entry in summary["skipped"]] == ["sample_sheet.csv"]

    validation = validate_catalog_inputs(config)
    assert validation["active_contrasts"] == 2
    assert validation["source_units"] == 2


def test_an_unsure_reader_skips_the_table_rather_than_defaulting(tmp_path) -> None:
    """Pressing enter must never stand in for "yes" on the direction question."""

    deg = tmp_path / "deg"
    deg.mkdir()
    _clean_table(deg / "clean_deseq2.csv")
    _clean_table(deg / "second.csv")

    # "" is unsure, then "yes" to the skip prompt; the second table is confirmed.
    answers = iter(["human", "", "yes", "yes", "hypoxia vs normoxia", "P1", "3", "3"])
    config = tmp_path / "config.csv"

    summary = run_init(config, deg, ask=lambda question, default="": next(answers), echo=lambda _line: None)

    assert summary["n_contrasts"] == 1
    assert summary["skipped"] == [
        {"path": "clean_deseq2.csv", "reason": "contrast direction not confirmed"}
    ]


def test_nothing_confirmed_writes_nothing(tmp_path) -> None:
    deg = tmp_path / "deg"
    deg.mkdir()
    _not_a_deg_table(deg / "sample_sheet.csv")
    config = tmp_path / "config.csv"

    with pytest.raises(ValueError, match="nothing to write"):
        run_init(config, deg, ask=lambda question, default="": "human", echo=lambda _line: None)

    assert not config.exists()


def test_an_existing_config_is_not_overwritten(tmp_path) -> None:
    deg = tmp_path / "deg"
    deg.mkdir()
    _clean_table(deg / "clean.csv")
    config = tmp_path / "config.csv"
    config.write_text("mine\n", encoding="utf-8")

    with pytest.raises(FileExistsError):
        run_init(config, deg, ask=lambda question, default="": "yes", echo=lambda _line: None)

    assert config.read_text(encoding="utf-8") == "mine\n"


def test_build_catalog_keeps_a_stable_column_order() -> None:
    frame = build_catalog([{"study_id": "S1"}])

    assert list(frame.columns)[:3] == ["study_id", "source_unit_id", "source_path"]
    assert "sign_convention" in frame.columns
