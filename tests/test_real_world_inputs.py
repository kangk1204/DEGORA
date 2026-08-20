"""Regressions for shapes that real, published DEG tables actually take.

Every case here was reproduced against a table downloaded from a public
repository, not invented.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from degora.harmonize import (
    ROW_LABEL_COLUMN,
    _normalize_separator,
    _restore_unnamed_row_labels,
)


@pytest.mark.parametrize(
    ("written", "expected"),
    [
        ("tab", "\t"),
        ("TAB", "\t"),
        ("Tab", "\t"),
        ("tsv", "\t"),
        ("\\t", "\t"),
        ("\t", "\t"),
        ("comma", ","),
        ("csv", ","),
        (",", ","),
        ("semicolon", ";"),
        (";", ";"),
        ("pipe", "|"),
        ("|", "|"),
        ("whitespace", r"\s+"),
    ],
)
def test_separator_words_are_accepted_not_treated_as_a_regex(written: str, expected: str) -> None:
    """`sep: tab` reached pandas as a three-character regex and collapsed every row."""

    assert _normalize_separator(written) == expected


def test_a_tab_separated_table_loads_when_sep_is_spelled_out(tmp_path: Path) -> None:
    from degora.harmonize import TableMapping, read_deg_table

    path = tmp_path / "topTags.tsv"
    path.write_text("ORF\tlogFC\tP.Value\nYAL001C\t1.5\t0.001\nYAL002W\t-2.0\t0.02\n", encoding="utf-8")
    mapping = TableMapping(gene_column="ORF", lfc_column="logFC", p_column="P.Value", sep="tab")
    frame = read_deg_table(path, mapping)
    assert list(frame.columns) == ["ORF", "logFC", "P.Value"]
    assert len(frame) == 2


def test_a_wrong_separator_now_explains_itself(tmp_path: Path) -> None:
    """The recovery hint used to fire only when the delimiter was auto-detected."""

    from degora.harmonize import TableMapping, read_deg_table

    path = tmp_path / "results.csv"
    path.write_text("gene,logFC,pvalue\nTP53,1.5,0.001\nMYC,-2.0,0.02\n", encoding="utf-8")
    mapping = TableMapping(gene_column="gene", lfc_column="logFC", p_column="pvalue", sep="tab")
    with pytest.raises(ValueError) as excinfo:
        read_deg_table(path, mapping)
    message = str(excinfo.value)
    assert "single column" in message
    assert "comma-delimited" in message
    assert "auto-detect" in message


def test_r_write_csv_row_labels_become_a_named_column() -> None:
    """`write.csv(res, file)` writes one fewer header field than data fields.

    pandas resolves that by consuming the gene identifiers as an unnamed index,
    which put them somewhere no catalog mapping could reference.
    """

    frame = pd.DataFrame(
        {"baseMean": [10.0, 20.0], "log2FoldChange": [1.0, -1.0], "pvalue": [0.01, 0.2]},
        index=["ENSG00000000003", "ENSG00000000005"],
    )
    assert frame.index.name is None

    restored = _restore_unnamed_row_labels(frame)

    assert list(restored.columns)[0] == ROW_LABEL_COLUMN
    assert list(restored[ROW_LABEL_COLUMN]) == ["ENSG00000000003", "ENSG00000000005"]


def test_a_normal_table_is_left_alone() -> None:
    frame = pd.DataFrame({"gene": ["TP53"], "logFC": [1.0]})
    assert _restore_unnamed_row_labels(frame) is frame


def test_a_named_index_is_left_alone() -> None:
    frame = pd.DataFrame({"logFC": [1.0]}, index=pd.Index(["TP53"], name="gene"))
    assert _restore_unnamed_row_labels(frame) is frame


def test_row_label_restoration_does_not_collide_with_an_existing_column() -> None:
    frame = pd.DataFrame({ROW_LABEL_COLUMN: ["x"], "logFC": [1.0]}, index=["ENSG1"])
    restored = _restore_unnamed_row_labels(frame)
    assert f"{ROW_LABEL_COLUMN}_2" in restored.columns
    assert ROW_LABEL_COLUMN in restored.columns


def test_missing_gene_column_points_at_the_restored_row_labels() -> None:
    from degora.harmonize import _series_as_numeric

    frame = pd.DataFrame({ROW_LABEL_COLUMN: ["ENSG1"], "log2FoldChange": [1.0]})
    with pytest.raises(KeyError) as excinfo:
        _series_as_numeric(frame, "gene")
    assert "write.csv" in str(excinfo.value)
    assert ROW_LABEL_COLUMN in str(excinfo.value)


def test_map_unique_applies_the_function_once_per_distinct_value() -> None:
    """These frames repeat a handful of labels across tens of thousands of rows."""

    from degora.score_db import _map_unique

    calls: list[object] = []

    def record(value: object) -> str:
        calls.append(value)
        return f"seen:{value}"

    values = pd.Series(["a", "b", "a", "b", "a"] * 100)
    mapped = _map_unique(values, record)

    assert len(calls) == 2, f"expected one call per distinct value, got {len(calls)}"
    assert list(mapped[:4]) == ["seen:a", "seen:b", "seen:a", "seen:b"]
    assert len(mapped) == len(values)


def test_map_unique_maps_missing_values_too() -> None:
    from degora.score_db import _map_unique

    values = pd.Series(["a", None, "a"])
    mapped = _map_unique(values, lambda value: "NA" if pd.isna(value) else str(value).upper())
    assert list(mapped) == ["A", "NA", "A"]


def test_map_unique_preserves_the_index() -> None:
    from degora.score_db import _map_unique

    values = pd.Series(["a", "b"], index=[7, 9])
    mapped = _map_unique(values, str.upper)
    assert list(mapped.index) == [7, 9]
