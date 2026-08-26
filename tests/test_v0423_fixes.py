"""Regressions for the v0.4.23 fixes: fan-out, identifier spaces, counts, CLI guards."""

from __future__ import annotations


import pandas as pd
import pytest


def test_a_publication_linking_many_series_is_set_aside_with_the_count(tmp_path, monkeypatch) -> None:
    """20 selected publications became 69 series; one consortium paper linked 51."""

    import degora.discovery_prepare as prep
    from degora.discovery import normalize_species

    captured: dict = {}

    def fake_prepare_geo_studies(accessions, species, **kwargs):
        captured["accessions"] = list(accessions)
        return {"studies": [], "excluded_studies": []}

    monkeypatch.setattr(prep, "prepare_geo_studies", fake_prepare_geo_studies)
    mega = {"canonical_id": "pmid:22955616", "title": "consortium", "geo_accessions": [f"GSE{100 + i}" for i in range(8)]}
    normal = {"canonical_id": "accession:GSE343715", "title": "one series", "geo_accessions": ["GSE343715"]}
    excluded: list = []

    prep._prepare_into_staging(
        geo_records=[mega, normal], direct_records=[], spec=normalize_species("human"), query="placenta",
        staging=tmp_path, max_files_per_record=3, transport=None, geo_client=None, excluded=excluded,
    )

    assert captured["accessions"] == ["GSE343715"]
    reasons = [str(item.get("reason") or item.get("exclusion_reason") or item) for item in excluded]
    assert any("links 8 GEO series" in reason and "select them directly" in reason for reason in reasons)


def test_the_inspector_prefers_the_column_whose_values_are_symbols_and_names_the_space() -> None:
    """A column called "Gene Symbol" holding Ensembl IDs joins nothing written in symbols."""

    from degora.discovery import _inspect_rows

    rows = [["GENEID", "GENENAME", "SYMBOL", "logFC", "PValue"]] + [
        [8404 + i, f"protein like {i}", f"SPARCL{i}", 1.0 - i / 100, 0.001] for i in range(30)
    ]
    header = _inspect_rows(rows)
    assert header["mapping"]["gene_column"] == "SYMBOL"
    assert header["gene_identifier_space"] == "gene symbol"

    ensembl_only = [["Gene Symbol", "log2FoldChange", "pvalue"]] + [[f"ENSG{i:011d}", 1.0, 0.001] for i in range(30)]
    header = _inspect_rows(ensembl_only)
    assert header["gene_identifier_space"] == "Ensembl ID"


def test_a_source_unit_in_another_identifier_space_is_flagged_even_with_partial_overlap() -> None:
    """869 genes, all ENSG, reported as a success: overlap alone missed it."""

    from degora.slice_runner import _identifier_space_warnings

    symbols = [f"GENE{i}" for i in range(400)] + [f"ENSG{i:011d}" for i in range(20)]
    ensembl = [f"ENSG{i:011d}" for i in range(450)]
    harmonized = pd.DataFrame(
        {"gene_symbol": symbols + ensembl, "source_unit_id": ["U_symbols"] * len(symbols) + ["U_ensembl"] * len(ensembl)}
    )

    warnings = _identifier_space_warnings(harmonized, min_studies=2)

    assert any("written in Ensembl ID" in w and "gene symbol" in w for w in warnings)


def test_two_symbol_units_with_good_overlap_raise_no_space_warning() -> None:
    from degora.slice_runner import _identifier_space_warnings

    genes = [f"GENE{i}" for i in range(300)]
    harmonized = pd.DataFrame({"gene_symbol": genes + genes[:250], "source_unit_id": ["U1"] * 300 + ["U2"] * 250})

    assert _identifier_space_warnings(harmonized, min_studies=2) == []


def test_the_cli_refuses_a_one_character_query_and_page_zero(tmp_path, capsys) -> None:
    from degora.cli import main

    assert main(["discover", "x", "--species", "human", "--output-dir", str(tmp_path / "a")]) == 2
    assert "at least 2 characters" in capsys.readouterr().err
    with pytest.raises(SystemExit) as excinfo:
        main(["discover", "hypoxia", "--species", "human", "--page", "0", "--output-dir", str(tmp_path / "b")])
    assert excinfo.value.code == 2
    assert "must be 1 or more" in capsys.readouterr().err


def test_serve_without_a_database_says_how_to_get_one(tmp_path) -> None:
    from degora.api import ScoreDatabaseError, _require_degora_score_database

    with pytest.raises(ScoreDatabaseError) as excinfo:
        _require_degora_score_database(tmp_path / "missing.db")

    assert "No run yet?" in str(excinfo.value)
    assert "degora demo" in str(excinfo.value)
