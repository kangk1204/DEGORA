"""Regressions for the v0.4.39 review findings.

Each test is named for what a reader actually hit, not for the function it calls:
a gene that scored twice because two papers spelled it differently, a reversed
contrast nothing warned about, a direction stated in the table that nobody was
shown, a GEO block reported as a missing record, a cuffdiff export refused with
advice that did not apply to it, and a preparation that made the reader run the
whole search again.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora import GENE_SYMBOL_RESOLUTION_VERSION
import degora.harmonize as harmonize_module
from degora.discovery import geo_soft_unavailable_reason
from degora.harmonize import (
    HGNC_SYMBOL_TABLE_PATH,
    TableMapping,
    canonical_gene_symbol,
    contrast_direction_evidence,
    gene_symbol_table_metadata,
    looks_like_cuffdiff,
)
from degora.score_db import (
    WITHIN_SOURCE_DIRECTION_MAX_SIGN_AGREEMENT,
    WITHIN_SOURCE_DIRECTION_SPEARMAN,
    within_source_direction_warnings,
)


# --------------------------------------------------------------------------
# A gene scored twice because one paper wrote CTGF and the other wrote CCN2
# --------------------------------------------------------------------------
RETIREMENTS = [
    ("CTGF", "CCN2"),
    ("CYR61", "CCN1"),
    ("NOV", "CCN3"),
    ("IL8", "CXCL8"),
    ("KIAA0101", "PCLAF"),
    ("MLL", "KMT2A"),
    ("HIST1H1C", "H1-2"),
    ("CARS", "CARS1"),
    ("FAM46C", "TENT5C"),
    ("ADRBK1", "GRK2"),
    ("C11ORF30", "EMSY"),
]


@pytest.mark.parametrize(("previous", "current"), RETIREMENTS)
def test_a_retired_symbol_and_its_current_symbol_are_one_gene(previous: str, current: str) -> None:
    assert canonical_gene_symbol(previous, species="Homo sapiens") == current
    assert canonical_gene_symbol(current, species="Homo sapiens") == current
    assert canonical_gene_symbol(previous.lower(), species="human") == current


def test_the_symbol_table_ships_inside_the_package() -> None:
    assert HGNC_SYMBOL_TABLE_PATH.is_file()
    assert HGNC_SYMBOL_TABLE_PATH.parent.name == "data"
    metadata = gene_symbol_table_metadata()
    assert metadata["gene_symbol_table_snapshot_date"]
    assert len(metadata["gene_symbol_table_source_sha256"]) == 64
    assert len(metadata["gene_symbol_table_sha256"]) == 64


def test_a_truncated_symbol_table_aborts_instead_of_silently_splitting_genes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    corrupt = tmp_path / "hgnc_previous_symbols.tsv"
    corrupt.write_text(
        "# DEGORA gene-symbol retirement table\n"
        "# source\thttps://example.invalid/hgnc.tsv\n"
        "# snapshot_date\t2026-09-02\n"
        f"# source_sha256\t{'0' * 64}\n"
        "previous_symbol\tcurrent_symbol\n"
        "CTGF\tCCN2\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(harmonize_module, "HGNC_SYMBOL_TABLE_PATH", corrupt)
    with pytest.raises(RuntimeError, match="only 1 mappings"):
        harmonize_module._load_retired_symbols()


def test_a_symbol_hgnc_also_uses_for_a_live_gene_is_never_rewritten() -> None:
    # HGNC reuses symbols: BRF1, AK3 and CCNL1 are approved genes and also appear
    # as previous symbols of other genes. Rewriting them would destroy a real gene,
    # so the table excludes every previous symbol that is itself approved.
    for symbol in ("BRF1", "AK3", "CCNL1", "ADCY3", "BRMS1", "CDH3"):
        assert canonical_gene_symbol(symbol) == symbol


def test_an_ambiguous_retirement_keeps_the_documented_choice() -> None:
    # DEC1 is a previous symbol of both BHLHE40 and DELEC1, so it is absent from the
    # table and the documented Excel-damage rule stays the only one that decides it.
    assert canonical_gene_symbol("DEC1") == "BHLHE40"


def test_the_earlier_symbol_rules_still_win_where_they_applied() -> None:
    assert canonical_gene_symbol("SEPT9") == "SEPTIN9"
    assert canonical_gene_symbol("9-Sep") == "SEPTIN9"
    assert canonical_gene_symbol("MARCH1") == "MARCHF1"
    assert canonical_gene_symbol("ENSG00000141510.16") == "ENSG00000141510"
    assert canonical_gene_symbol("7157.0") == "7157"
    # A dotted symbol is still a symbol, not a versioned accession.
    assert canonical_gene_symbol("NKX2.5") == "NKX2.5"
    for blank in ("", "NA", "<NA>", "#N/A", None):
        assert canonical_gene_symbol(blank) == ""


def test_two_tables_that_spell_one_gene_differently_join(tmp_path: Path) -> None:
    from degora.harmonize import harmonize_frame

    old = pd.DataFrame({"gene": ["CTGF", "IL8", "TP53"], "lfc": [2.0, 1.5, 0.3], "p": [1e-8, 1e-6, 0.4]})
    new = pd.DataFrame({"gene": ["CCN2", "CXCL8", "TP53"], "lfc": [1.8, 1.4, 0.2], "p": [1e-7, 1e-5, 0.5]})
    mapping = TableMapping(gene_column="gene", lfc_column="lfc", p_column="p")
    frames = []
    for name, frame in (("old_paper", old), ("new_paper", new)):
        harmonized = harmonize_frame(
            frame,
            mapping,
            {
                "study_id": name,
                "paper_id": name,
                "source_unit_id": name,
                "species": "Homo sapiens",
            },
        )
        frames.append(harmonized[0] if isinstance(harmonized, tuple) else harmonized)
    combined = pd.concat(frames, ignore_index=True)
    per_gene = combined.groupby("gene_symbol")["source_unit_id"].nunique()
    assert per_gene.get("CCN2") == 2, "CTGF and CCN2 must count as one gene with two sources"
    assert per_gene.get("CXCL8") == 2
    assert "CTGF" not in per_gene.index and "IL8" not in per_gene.index
    # The label the source actually carried stays recoverable.
    assert set(combined["input_gene_label"]) >= {"CTGF", "CCN2", "IL8", "CXCL8"}


def test_the_resolution_version_is_recorded_beside_the_score_version() -> None:
    assert "hgnc" in GENE_SYMBOL_RESOLUTION_VERSION


# --------------------------------------------------------------------------
# A contrast written control-minus-treatment inside one source unit
# --------------------------------------------------------------------------
def _reversed_pair_frame(rho_target: str) -> pd.DataFrame:
    rng = np.random.default_rng(20260901)
    n = 4_000
    genes = [f"GENE{i:05d}" for i in range(n)]
    base = rng.normal(0, 1.2, n)
    noise = rng.normal(0, 0.9, n)
    forward = base + rng.normal(0, 0.6, n)
    reversed_ = -(base + noise)
    if rho_target == "aligned":
        reversed_ = base + noise
    rows = []
    for study, values in (("unit_forward", forward), ("unit_other", reversed_)):
        rows.append(
            pd.DataFrame(
                {
                    "study_id": study,
                    "source_unit_id": "GSE1",
                    "paper_id": "GSE1",
                    "gene_symbol": genes,
                    "lfc": values,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def test_a_reversed_contrast_inside_one_source_unit_is_named() -> None:
    warnings = within_source_direction_warnings(_reversed_pair_frame("reversed"))
    assert warnings, "a contrast written control-minus-treatment must be named"
    message = warnings[0]
    assert "'unit_forward'" in message and "'unit_other'" in message
    assert "GSE1" in message
    assert "biologically opposing comparisons" in message
    assert "may have a reversed contrast direction" in message
    assert "does not identify which contrast" in message


def test_two_well_formed_contrasts_of_one_unit_stay_silent() -> None:
    assert within_source_direction_warnings(_reversed_pair_frame("aligned")) == []


def test_the_reversal_bounds_are_reachable_by_a_real_reversed_contrast() -> None:
    # Measured on public data, a real reversed contrast pairs at Spearman -0.23 to
    # -0.50 with 33-43% same-sign agreement; the previous -0.80 / 10% pair could only
    # be reached by a table flipped against itself, so it never named a real one.
    assert -0.6 < WITHIN_SOURCE_DIRECTION_SPEARMAN <= -0.10
    assert 0.30 <= WITHIN_SOURCE_DIRECTION_MAX_SIGN_AGREEMENT < 0.50


def test_a_noisy_twenty_gene_fixture_does_not_trip_the_advisory() -> None:
    rng = np.random.default_rng(7)
    genes = [f"G{i}" for i in range(24)]
    frame = pd.concat(
        [
            pd.DataFrame(
                {
                    "study_id": study,
                    "source_unit_id": "GSE2",
                    "paper_id": "GSE2",
                    "gene_symbol": genes,
                    "lfc": rng.normal(0, 1, len(genes)),
                }
            )
            for study in ("a", "b")
        ],
        ignore_index=True,
    )
    # Independent noise can drift past the correlation bound on 24 genes; the
    # one-sided significance test is what keeps a small fixture quiet.
    assert within_source_direction_warnings(frame) == []


# --------------------------------------------------------------------------
# The direction the table itself states
# --------------------------------------------------------------------------
def test_a_deseq2_results_header_names_its_own_numerator() -> None:
    column = "log2 fold change (MLE): group Starv Control vs Starv 3h TGFb"
    frame = pd.DataFrame({"gene_name": ["A", "B"], column: [1.0, -1.0]})
    message = contrast_direction_evidence(
        frame, TableMapping(gene_column="gene_name", lfc_column=column, p_column="")
    )
    assert "'group Starv Control'" in message
    assert "numerator" in message
    assert "Trust this over the file name" in message


@pytest.mark.parametrize("separator", ["VS", "Vs", "vS"])
def test_a_deseq2_results_header_parses_vs_case_insensitively(separator: str) -> None:
    column = f"log2 fold change (MLE): group Treated {separator} group Control"
    frame = pd.DataFrame({"gene_name": ["A", "B"], column: [1.0, -1.0]})
    message = contrast_direction_evidence(
        frame, TableMapping(gene_column="gene_name", lfc_column=column, p_column="")
    )
    assert "'group Treated' is the numerator" in message
    assert "'group Control' the denominator" in message


def test_a_cuffdiff_table_is_read_from_sample_1_and_sample_2_not_its_file_name() -> None:
    rng = np.random.default_rng(3)
    n = 200
    value_1 = rng.uniform(1, 50, n)
    value_2 = value_1 * rng.uniform(0.2, 5.0, n)
    frame = pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(n)],
            "sample_1": "Dex",
            "sample_2": "Untreated",
            "value_1": value_1,
            "value_2": value_2,
            "log2(fold_change)": np.log2(value_2 / value_1),
            "p_value": rng.uniform(0, 1, n),
        }
    )
    mapping = TableMapping(gene_column="gene", lfc_column="log2(fold_change)", p_column="p_value")
    message = contrast_direction_evidence(frame, mapping)
    assert "'Untreated'" in message and "'Dex'" in message
    assert "positive value means up in 'Untreated'" in message

    inverted = frame.assign(**{"log2(fold_change)": -frame["log2(fold_change)"]})
    message = contrast_direction_evidence(inverted, mapping)
    assert "positive value means up in 'Dex'" in message
    assert "inverted since cuffdiff wrote it" in message


def test_a_multi_comparison_cuffdiff_table_does_not_report_only_the_first_pair() -> None:
    rows: list[dict[str, object]] = []
    for sample_1, sample_2, multiplier in (("A", "B", 2.0), ("C", "D", 0.5)):
        for index in range(25):
            value_1 = float(index + 1)
            value_2 = value_1 * multiplier
            rows.append(
                {
                    "gene": f"{sample_1}{index}",
                    "sample_1": sample_1,
                    "sample_2": sample_2,
                    "value_1": value_1,
                    "value_2": value_2,
                    "log2(fold_change)": np.log2(value_2 / value_1),
                    "p_value": 0.01,
                }
            )
    frame = pd.DataFrame(rows)
    mapping = TableMapping(gene_column="gene", lfc_column="log2(fold_change)", p_column="p_value")
    message = contrast_direction_evidence(frame, mapping)
    assert "2 sample_1/sample_2 comparisons" in message
    assert "A" in message and "B" in message and "C" in message and "D" in message
    assert "no single contrast direction is inferred" in message


def test_a_plain_results_table_says_nothing_rather_than_guessing() -> None:
    frame = pd.DataFrame({"gene_id": ["A"] * 50, "log2FoldChange": [0.5] * 50, "pvalue": [0.1] * 50})
    assert (
        contrast_direction_evidence(
            frame, TableMapping(gene_column="gene_id", lfc_column="log2FoldChange", p_column="pvalue")
        )
        == ""
    )


def test_validate_reports_the_direction_the_table_states(tmp_path: Path) -> None:
    from degora.slice_runner import validate_catalog_inputs

    column = "log2 fold change (MLE): group Ctrl vs Treated"
    rows = []
    for index in (1, 2):
        table = tmp_path / f"deg{index}.csv"
        pd.DataFrame(
            {
                "gene": [f"G{i}" for i in range(40)],
                column: np.linspace(-3, 3, 40) * index,
                "pvalue": np.linspace(0.001, 0.9, 40),
            }
        ).to_csv(table, index=False)
        rows.append(
            {
                "study_id": f"S{index}",
                "source_unit_id": f"U{index}",
                "source_path": table.name,
                "gene_column": "gene",
                "lfc_column": column,
                "p_column": "pvalue",
                "include_in_analysis": "yes",
            }
        )
    catalog = tmp_path / "config.csv"
    pd.DataFrame(rows).to_csv(catalog, index=False)
    result = validate_catalog_inputs(catalog)
    warnings = " ".join(result.get("warnings", []))
    assert "'group Ctrl vs Treated'" in warnings
    assert "confirm this is the treatment-minus-control direction" in warnings


# --------------------------------------------------------------------------
# GEO answering a blocked network with a web page
# --------------------------------------------------------------------------
def test_a_geo_captcha_page_is_reported_as_a_block_not_a_missing_record() -> None:
    body = (
        '<!doctype html><html lang="en-US"><head>'
        '<base href="https://www.google.com/recaptcha/challengepage/"></head><body></body></html>'
    )
    reason = geo_soft_unavailable_reason("GSE52778", body)
    assert "CAPTCHA" in reason
    assert "blocked or rate-limited" in reason
    assert "does not establish whether the Series record exists" in reason


def test_a_generic_geo_html_page_is_not_called_rate_limiting_or_a_valid_record() -> None:
    reason = geo_soft_unavailable_reason("GSE123", "<html><head><title>Maintenance</title></head></html>")
    assert "not enough to diagnose rate limiting" in reason
    assert "verify the accession" in reason
    assert "no Series SOFT record" not in reason


def test_an_empty_and_a_genuinely_wrong_soft_response_keep_their_own_words() -> None:
    assert "empty SOFT response" in geo_soft_unavailable_reason("GSE1", "   ")
    assert geo_soft_unavailable_reason("GSE1", "^PLATFORM = GPL1\n") == (
        "GEO returned no Series SOFT record for GSE1"
    )


# --------------------------------------------------------------------------
# A cuffdiff export refused with advice for a different problem
# --------------------------------------------------------------------------
def test_a_cuffdiff_table_is_recognised() -> None:
    frame = pd.DataFrame(
        columns=["test_id", "gene", "sample_1", "sample_2", "status", "value_1", "value_2", "log2(fold_change)"]
    )
    assert looks_like_cuffdiff(frame)
    assert not looks_like_cuffdiff(pd.DataFrame(columns=["gene", "log2FoldChange", "pvalue"]))


def test_cuffdiff_infinities_are_refused_with_the_cuffdiff_fix(tmp_path: Path) -> None:
    from degora.slice_runner import DegoraConfigError, validate_catalog_inputs

    n = 400
    lfc = np.linspace(-4, 4, n)
    lfc[:60] = 1.7976931348623157e308
    frame = pd.DataFrame(
        {
            "test_id": [f"G{i}" for i in range(n)],
            "gene": [f"G{i}" for i in range(n)],
            "sample_1": "PBS",
            "sample_2": "Drug",
            "status": ["OK"] * (n - 60) + ["NOTEST"] * 60,
            "value_1": np.linspace(1, 100, n),
            "value_2": np.linspace(1, 100, n),
            "log2(fold_change)": lfc,
            "p_value": np.linspace(0.001, 0.9, n),
        }
    )
    table = tmp_path / "gene_exp.diff"
    frame.to_csv(table, sep="\t", index=False)
    catalog = tmp_path / "config.csv"
    pd.DataFrame(
        [
            {
                "study_id": "S1",
                "source_unit_id": "U1",
                "source_path": table.name,
                "gene_column": "gene",
                "lfc_column": "log2(fold_change)",
                "p_column": "p_value",
                "sep": "tab",
                "include_in_analysis": "yes",
            }
        ]
    ).to_csv(catalog, index=False)
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(catalog)
    message = str(excinfo.value)
    assert "cuffdiff gene_exp.diff table" in message
    assert "status" in message
    assert "convert 'log2(fold_change)' to log2" not in message


# --------------------------------------------------------------------------
# Preparing without running the whole search again
# --------------------------------------------------------------------------
def _write_snapshot(directory: Path, *, query: str, species: str) -> None:
    from degora.discovery_export import export_publication_search

    snapshot = {
        "query": query,
        "species": {"key": species, "label": species.title(), "scientific_name": "Homo sapiens"},
        "generated_at": "2026-09-01T00:00:00+00:00",
        "records": [
            {
                "canonical_id": "pmid:1",
                "publication_id": "pmid:1",
                "paper_title": "A study",
                "pmid": "1",
                "pubmed_ids": ["1"],
                "source_unit_id": "GSE1",
                "geo_accessions": ["GSE1"],
                "readiness": "candidate",
                "species": species,
            }
        ],
        "total_records": 1,
        "evaluated_records": 1,
    }
    export_publication_search(snapshot, directory)


def test_from_snapshot_prepares_without_searching_again(tmp_path: Path) -> None:
    from degora.cli import _load_publication_snapshot

    snapshot_dir = tmp_path / "search"
    _write_snapshot(snapshot_dir, query="topic here", species="human")
    loaded = _load_publication_snapshot(snapshot_dir, "topic here", "human")
    assert len(loaded["records"]) == 1
    assert loaded["records"][0]["canonical_id"] == "pmid:1"


def test_from_snapshot_refuses_the_other_species_and_another_query(tmp_path: Path) -> None:
    from degora.cli import CliUsageError, _load_publication_snapshot

    snapshot_dir = tmp_path / "search"
    _write_snapshot(snapshot_dir, query="topic here", species="human")
    with pytest.raises(CliUsageError, match="Human and Mouse are kept in separate workspaces"):
        _load_publication_snapshot(snapshot_dir, "topic here", "mouse")
    with pytest.raises(CliUsageError, match="Pass the snapshot's own query"):
        _load_publication_snapshot(snapshot_dir, "another topic", "human")


def test_from_snapshot_refuses_a_tampered_export(tmp_path: Path) -> None:
    from degora.cli import CliUsageError, _load_publication_snapshot
    from degora.discovery_export import SEARCH_CSV_NAME

    snapshot_dir = tmp_path / "search"
    _write_snapshot(snapshot_dir, query="topic here", species="human")
    (snapshot_dir / SEARCH_CSV_NAME).write_text("tampered\n", encoding="utf-8")
    with pytest.raises(CliUsageError, match="not a complete federated search export"):
        _load_publication_snapshot(snapshot_dir, "topic here", "human")


@pytest.mark.parametrize(
    ("missing_field", "message"),
    [("species", "does not record the search species"), ("query", "does not record the search query")],
)
def test_from_snapshot_refuses_missing_scope_metadata(
    tmp_path: Path,
    missing_field: str,
    message: str,
) -> None:
    from degora.cli import CliUsageError, _load_publication_snapshot
    from degora.discovery_export import SEARCH_JSON_NAME, export_publication_search

    snapshot_dir = tmp_path / "search"
    _write_snapshot(snapshot_dir, query="topic here", species="human")
    snapshot = json.loads((snapshot_dir / SEARCH_JSON_NAME).read_text(encoding="utf-8"))
    snapshot.pop(missing_field)
    export_publication_search(snapshot, snapshot_dir, force=True)
    with pytest.raises(CliUsageError, match=message):
        _load_publication_snapshot(snapshot_dir, "topic here", "human")


def test_from_snapshot_needs_a_selection_and_the_federated_backend() -> None:
    from degora.cli import main

    assert main(["discover", "topic", "--species", "human", "--from-snapshot", "."]) != 0
    assert main(
        ["discover", "topic", "--species", "human", "--source", "geo", "--from-snapshot", ".", "--select", "GSE1"]
    ) != 0


def test_a_search_without_output_dir_still_says_so() -> None:
    from degora.cli import main

    assert main(["discover", "topic", "--species", "human"]) != 0
