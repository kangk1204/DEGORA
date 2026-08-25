"""Regressions for the defects found auditing 0.4.16."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora.aggregate import time_course_selection_report, time_course_selection_warnings
from degora.discovery import classify_header
from degora.formula_safety import (
    formula_guard_metadata,
    neutralize_formula_text,
    restore_formula_text,
    restore_formula_text_if_marked,
)
from degora.provenance import write_source_sidecar
from degora.slice_runner import CATALOG_ALIASES, _normalize_catalog_columns


def _catalog(**columns) -> pd.DataFrame:
    base = {
        "study_id": ["S1", "S2"],
        "source_path": ["a.csv", "b.csv"],
        "gene_column": ["gene"] * 2,
        "lfc_column": ["lfc"] * 2,
        "p_column": ["p"] * 2,
    }
    return pd.DataFrame({**base, **columns})


@pytest.mark.parametrize(("alias", "canonical"), sorted(CATALOG_ALIASES.items()))
def test_a_wholly_blank_canonical_column_accepts_its_alias(alias: str, canonical: str) -> None:
    """A column of empty cells reads as float64, and pandas 3 refuses the upcast.

    Writing the alias's text into an all-blank canonical column raised
    TypeError: Invalid value ... for dtype 'float64' and ended the run in a raw
    traceback. pandas 3 is inside the declared support range and is what the
    upper-bounds CI job installs. The most ordinary shape that hits it is a blank
    paper_id column beside a filled source_unit_id.
    """

    frame = _catalog(**{canonical: pd.Series([np.nan] * 2, dtype="float64"), alias: ["x", "y"]})

    promoted = _normalize_catalog_columns(frame)

    assert list(promoted[canonical]) == ["x", "y"]


def test_an_explicit_canonical_value_still_wins_over_the_alias() -> None:
    frame = _catalog(paper_id=["kept", ""], source_unit_id=["ignored", "promoted"])

    promoted = _normalize_catalog_columns(frame)

    assert list(promoted["paper_id"]) == ["kept", "promoted"]


def test_a_promoted_legacy_time_course_mode_is_announced() -> None:
    """Blank canonical means `mean`, so promoting `early` moves off the default.

    The release note said "configs that explicitly use early or late can change".
    A reader whose time_course_mode column is blank concludes it does not reach
    them, and their run silently drops every gene absent from the earliest
    contrast.
    """

    frame = _catalog(
        source_unit_id=["U1", "U2"],
        time_course_mode=pd.Series([np.nan] * 2, dtype="float64"),
        temporal_mode=["early", "early"],
    )

    promoted = _normalize_catalog_columns(frame)
    warnings = promoted.attrs["promoted_alias_warnings"]

    assert len(warnings) == 1
    assert "temporal_mode" in warnings[0]
    assert "blank" in warnings[0]
    assert "'mean'" in warnings[0]


def test_promoting_a_legacy_mean_is_not_worth_a_warning() -> None:
    frame = _catalog(
        source_unit_id=["U1", "U2"],
        time_course_mode=pd.Series([np.nan] * 2, dtype="float64"),
        temporal_mode=["mean", ""],
    )

    assert _normalize_catalog_columns(frame).attrs["promoted_alias_warnings"] == []


def test_a_stated_time_course_mode_is_not_reported_as_promoted() -> None:
    frame = _catalog(
        source_unit_id=["U1", "U2"],
        time_course_mode=["early", "early"],
        temporal_mode=["early", "early"],
    )

    assert _normalize_catalog_columns(frame).attrs["promoted_alias_warnings"] == []


def test_a_blank_row_is_reported_even_when_another_row_states_the_canonical_mode() -> None:
    frame = _catalog(
        source_unit_id=["U1", "U2"],
        time_course_mode=["mean", ""],
        temporal_mode=["", "early"],
    )

    promoted = _normalize_catalog_columns(frame)

    assert promoted["time_course_mode"].tolist() == ["mean", "early"]
    assert len(promoted.attrs["promoted_alias_warnings"]) == 1
    assert "early" in promoted.attrs["promoted_alias_warnings"][0]


def _timed(unit: str, genes: int, duration: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": [f"G{index}" for index in range(genes)],
            "study_id": [f"{unit}_{duration}"] * genes,
            "source_unit_id": [unit] * genes,
            "signed_z": [3.0] * genes,
            "duration_h": [duration] * genes,
            "time_course_mode": ["early"] * genes,
        }
    )


def test_early_selection_reports_what_it_dropped() -> None:
    """A 200-gene 24h table beside a 2-gene pilot kept two rows and said nothing."""

    harmonized = pd.concat([_timed("U1", 2, 0.5), _timed("U1", 200, 24.0)], ignore_index=True)

    report = time_course_selection_report(harmonized)

    assert len(report) == 1
    entry = report[0]
    assert entry["source_unit_id"] == "U1"
    assert entry["genes_before"] == 200
    assert entry["genes_after"] == 2
    assert time_course_selection_warnings(report), "a 1% retention must be warned about"


def test_a_unit_keeping_most_of_its_genes_is_not_warned_about() -> None:
    harmonized = pd.concat([_timed("U1", 190, 2.0), _timed("U1", 200, 24.0)], ignore_index=True)

    assert time_course_selection_warnings(time_course_selection_report(harmonized)) == []


def test_mean_mode_produces_no_selection_report() -> None:
    frame = _timed("U1", 10, 4.0)
    frame["time_course_mode"] = "mean"

    assert time_course_selection_report(frame) == []


def test_published_text_cannot_be_read_as_a_formula() -> None:
    """Gene names arrive from downloaded tables, so a formula cell is an input."""

    frame = pd.DataFrame(
        {
            "gene_symbol": ['=HYPERLINK("http://evil","CLICK")', "TP53", "-SUM(A1)", "@cmd", "+x"],
            "degora_score": [-1.5, 2.0, 3.25, 0.0, -0.5],
            "n_studies": [1, 2, 3, 4, 5],
        }
    )

    guarded = neutralize_formula_text(frame)

    assert guarded["gene_symbol"].tolist() == [
        '\'=HYPERLINK("http://evil","CLICK")',
        "TP53",
        "'-SUM(A1)",
        "'@cmd",
        "'+x",
    ]
    # Numbers keep their type and value: prefixing a negative score would corrupt it.
    assert guarded["degora_score"].tolist() == [-1.5, 2.0, 3.25, 0.0, -0.5]
    assert str(guarded["degora_score"].dtype) == "float64"
    assert str(guarded["n_studies"].dtype).startswith("int")


def test_formula_guard_is_reversible_even_for_a_source_apostrophe() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": ["=BAD()", "'=LITERAL", "  #REF!", "TP53"],
            "score": [1.0, -2.0, 3.0, 4.0],
        }
    )

    guarded = neutralize_formula_text(frame)

    assert guarded["gene_symbol"].tolist() == ["'=BAD()", "''=LITERAL", "'  #REF!", "TP53"]
    pd.testing.assert_frame_equal(restore_formula_text(guarded), frame)


def test_formula_guard_is_restored_only_with_matching_provenance(tmp_path: Path) -> None:
    raw = pd.DataFrame({"gene_symbol": ["=BAD()", "'=LITERAL", "TP53"]})
    path = tmp_path / "guarded.csv"
    neutralize_formula_text(raw).to_csv(path, index=False)
    write_source_sidecar(path, "degora test", metadata=formula_guard_metadata())

    loaded = pd.read_csv(path)

    pd.testing.assert_frame_equal(restore_formula_text_if_marked(loaded, path), raw)


def test_tampered_formula_guard_provenance_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "guarded.csv"
    neutralize_formula_text(pd.DataFrame({"gene_symbol": ["=BAD()"]})).to_csv(path, index=False)
    write_source_sidecar(path, "degora test", metadata=formula_guard_metadata())
    path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid formula-guard provenance"):
        restore_formula_text_if_marked(pd.read_csv(path), path)


def test_unmarked_apostrophe_guarded_text_fails_as_ambiguous(tmp_path: Path) -> None:
    path = tmp_path / "guarded.csv"
    pd.DataFrame({"gene_symbol": ["'=BAD()"]}).to_csv(path, index=False)

    with pytest.raises(ValueError, match="no matching DEGORA formula-guard provenance"):
        restore_formula_text_if_marked(pd.read_csv(path), path)


def test_unmarked_raw_formula_text_remains_a_raw_identifier(tmp_path: Path) -> None:
    raw = pd.DataFrame({"gene_symbol": ["=BAD()", "TP53"]})
    path = tmp_path / "raw.csv"
    raw.to_csv(path, index=False)

    pd.testing.assert_frame_equal(restore_formula_text_if_marked(pd.read_csv(path), path), raw)


def test_ordinary_gene_symbols_are_left_exactly_alone() -> None:
    """The score contract depends on this: real data must round-trip unchanged."""

    frame = pd.DataFrame({"gene_symbol": ["TP53", "A1BG-AS1", "HLA-DRB1", "MIR210HG"]})

    assert neutralize_formula_text(frame)["gene_symbol"].tolist() == frame["gene_symbol"].tolist()


@pytest.mark.parametrize(
    "column",
    [
        "feature_id",
        "cluster_id",
        "pathway_id",
        "compound_id",
        "taxon_id",
        "publication_id",
        "dataset_id",
        "experiment_id",
        "protein_id",
        "peak_id",
        "variant_id",
    ],
)
def test_generic_or_non_gene_id_headers_are_not_automatic(column: str) -> None:
    assert classify_header([column, "logFC", "P.Value"])["mapping"]["gene_column"] == ""


@pytest.mark.parametrize(
    "column", ["sample_id", "GSM_id", "patient_id", "subject_id", "run_id", "file_id", "series_id"]
)
def test_sample_scoped_identifiers_stay_refused(column: str) -> None:
    """Admitting a trailing _id must not let a sample column through."""

    assert classify_header([column, "logFC", "P.Value"])["mapping"]["gene_column"] == ""


def test_a_hand_edited_bundle_species_does_not_end_in_a_traceback() -> None:
    """`species` as a bare string raised AttributeError on reader-correctable input."""

    from degora.discovery_run import DiscoveryError, run_discovery_analysis

    for payload in ({"species": "human"}, {"species": None}, {"species": []}, {}):
        with pytest.raises(DiscoveryError):
            run_discovery_analysis(
                payload, [], "/tmp/degora-unused", species="human", min_studies=2, force=True
            )


@pytest.mark.parametrize(
    "column",
    [
        "rank", "row_number", "baseMean", "mean_count", "stat", "score",
        "index", "position", "pathway", "metabolite", "cell_line", "compound",
    ],
)
def test_a_non_gene_column_is_never_auto_accepted_as_the_gene_column(column: str, tmp_path) -> None:
    """v0.4.16 as published accepted all twelve of these without asking.

    Two source units keyed on `rank` produced a config, a clean run, and a
    ranking whose top genes were 1, 2, 3 ... 10. Nothing failed and nothing
    warned, which makes it a scientific false success rather than a crash.
    """

    from degora.beginner import infer_source_table

    numeric = column not in {"pathway", "metabolite", "cell_line", "compound"}
    values = list(range(1, 201)) if numeric else [f"{column}{index}" for index in range(200)]
    path = tmp_path / f"{column}.csv"
    pd.DataFrame(
        {column: values, "log2FoldChange": [1.0] * 200, "pvalue": [0.001] * 200}
    ).to_csv(path, index=False)

    inference = infer_source_table(path)
    asked = {choice.role for choice in inference.needs_a_question}

    assert not (inference.looks_like_a_deg_table and "gene_column" not in asked), (
        f"{column!r} was accepted as a gene column without asking"
    )


@pytest.mark.parametrize(
    ("column", "values"),
    [
        ("gene_symbol", ["TP53", "BRCA1", "EGFR"]),
        ("gene_id", ["ENSG00000141510", "ENSG00000012048", "ENSG00000146648"]),
        ("entrez_id", ["7157", "672", "1956"]),
        ("probe_id", ["1007_s_at", "1053_at", "117_at"]),
    ],
)
def test_a_real_gene_column_is_still_accepted_without_a_question(column: str, values: list[str], tmp_path) -> None:
    """Refusing the twelve must not cost the shapes that should pass."""

    from degora.beginner import infer_source_table

    rows = (values * 70)[:200]
    path = tmp_path / f"{column}.csv"
    pd.DataFrame(
        {column: rows, "log2FoldChange": [1.0] * 200, "pvalue": [0.001] * 200}
    ).to_csv(path, index=False)

    inference = infer_source_table(path)

    assert inference.mapping.get("gene_column") == column
    assert "gene_column" not in {choice.role for choice in inference.needs_a_question}


def test_interrupting_the_cli_is_not_reported_as_a_fault(monkeypatch) -> None:
    """Ctrl+C on a slow search showed a traceback: it is a BaseException.

    The reader stopped a search that had printed nothing for 25 seconds and was
    handed a stack trace for the decision.
    """

    import sys
    import types

    from degora.cli import main

    def interrupted(*_args, **_kwargs):
        raise KeyboardInterrupt

    module = types.ModuleType("degora.discovery_federated")
    module.search_publications = interrupted
    module.page_publication_snapshot = lambda **_kwargs: {"records": [], "total": 0}
    module.resolve_publication_records = lambda *_a, **_k: []
    module.filter_publication_records = lambda records, text_filter="": list(records)
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)

    assert main(
        ["discover", "x", "--species", "human", "--limit", "5", "--output-dir", "/tmp/degora-interrupt"]
    ) == 130


def test_the_sdist_ships_the_script_its_readme_documents() -> None:
    """The README tells sdist readers to run scripts/degora_quickstart.sh."""

    manifest = Path("MANIFEST.in")

    assert manifest.exists(), "no MANIFEST.in, so the sdist carries only the package"
    assert "scripts/degora_quickstart.sh" in manifest.read_text(encoding="utf-8")
