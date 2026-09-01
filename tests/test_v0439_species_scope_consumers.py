from __future__ import annotations

import hashlib
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from degora.api import DegoraRequestHandler
from degora.api import _database_gene_species as api_database_gene_species
from degora.discovery import normalize_species
from degora.discovery_run import _fallback_row
from degora.excel_export import _curated_lookup, _read_gold_from_config
from degora.excel_export import _database_gene_species as export_database_gene_species
from degora.reanalysis import derive_welch_deg


def _write_gene_database(
    path: Path,
    *,
    genes: list[str],
    species: list[str | None] | None,
) -> Path:
    gene_rows = pd.DataFrame(
        {
            "gene_symbol": genes,
            "degora_rank": list(range(1, len(genes) + 1)),
            "degora_score": [1.0 - index * 0.1 for index in range(len(genes))],
            "n_source_units": [1] * len(genes),
            "consensus_direction": ["up"] * len(genes),
        }
    )
    study_rows: dict[str, list[object]] = {
        "study_id": [f"S{index}" for index in range(max(len(species or []), 1))],
        "source_unit_id": [f"P{index}" for index in range(max(len(species or []), 1))],
    }
    evidence_rows: dict[str, list[object]] = {
        "gene_symbol": genes,
        "study_id": ["S0"] * len(genes),
        "source_unit_id": ["P0"] * len(genes),
    }
    if species is not None:
        study_rows["species"] = species
        evidence_rows["species"] = [species[0]] * len(genes)
    with sqlite3.connect(path) as connection:
        gene_rows.to_sql("genes", connection, index=False)
        pd.DataFrame(study_rows).to_sql("studies", connection, index=False)
        pd.DataFrame(evidence_rows).to_sql("gene_evidence", connection, index=False)
    return path


def _handler() -> DegoraRequestHandler:
    handler = object.__new__(DegoraRequestHandler)
    handler.server = SimpleNamespace(server_address=("127.0.0.1", 0))
    return handler


@pytest.mark.parametrize(
    ("species", "expected"),
    [
        (["Homo sapiens"], "Homo sapiens"),
        (["human"], "Homo sapiens"),
        (["Mus musculus"], "Mus musculus"),
        (["human", "mouse"], None),
        ([None], None),
        (None, None),
    ],
)
def test_database_species_inference_is_complete_and_fail_closed(
    tmp_path: Path,
    species: list[str | None] | None,
    expected: str | None,
) -> None:
    database = _write_gene_database(
        tmp_path / "scope.db",
        genes=["ARF4"],
        species=species,
    )

    assert export_database_gene_species(database) == expected
    with sqlite3.connect(database) as connection:
        assert api_database_gene_species(connection) == expected


@pytest.mark.parametrize(
    ("species", "genes", "expected_search", "detail_query", "expected_detail"),
    [
        (["Homo sapiens"], ["ARF4"], ["ARF4"], "Arf2", "ARF4"),
        (["Mus musculus"], ["ARF2", "ARF4"], ["ARF2"], "Arf2", "ARF2"),
    ],
)
def test_api_symbol_lookup_uses_only_the_database_species(
    tmp_path: Path,
    species: list[str],
    genes: list[str],
    expected_search: list[str],
    detail_query: str,
    expected_detail: str,
) -> None:
    database = _write_gene_database(
        tmp_path / "api.db",
        genes=genes,
        species=species,
    )
    handler = _handler()

    search = handler._genes({"q": ["ARF2"]}, db_path=database)
    detail = handler._gene_detail(detail_query, db_path=database)

    assert [row["gene_symbol"] for row in search["genes"]] == expected_search
    assert detail["gene"]["gene_symbol"] == expected_detail


@pytest.mark.parametrize("species", [["human", "mouse"], [None], None])
def test_api_mixed_or_missing_species_never_applies_hgnc(
    tmp_path: Path,
    species: list[str | None] | None,
) -> None:
    database = _write_gene_database(
        tmp_path / "generic.db",
        genes=["ARF4"],
        species=species,
    )

    result = _handler()._genes({"q": ["ARF2"]}, db_path=database)

    assert result["count"] == 0


@pytest.mark.parametrize(
    ("species", "expected_resolved", "expected_present"),
    [
        (["Homo sapiens"], "ARF4", True),
        (["Mus musculus"], "ARF2", False),
        (["human", "mouse"], "ARF2", False),
        (None, "ARF2", False),
    ],
)
def test_gold_panel_resolution_uses_the_full_database_species_scope(
    tmp_path: Path,
    species: list[str] | None,
    expected_resolved: str,
    expected_present: bool,
) -> None:
    database = _write_gene_database(
        tmp_path / "workbook.db",
        genes=["ARF4"],
        species=species,
    )
    config = tmp_path / "config.xlsx"
    with pd.ExcelWriter(config, engine="openpyxl") as writer:
        pd.DataFrame({"gene_symbol": ["ARF2"], "locked": ["yes"]}).to_excel(
            writer,
            sheet_name="GoldPanel",
            index=False,
        )
    inferred = export_database_gene_species(database)

    gold, status, reason = _read_gold_from_config(config, species=inferred)
    lookup = _curated_lookup(
        gold,
        pd.DataFrame({"gene_symbol": ["ARF4"], "quality_weighted_degora_rank": [1]}),
        species=inferred,
    )

    assert status == "locked"
    assert reason == ""
    assert gold.loc[0, "resolved_gene_symbol"] == expected_resolved
    assert bool(lookup.loc[0, "present_in_degora_output"]) is expected_present


@pytest.mark.parametrize(
    ("species", "expected_symbols"),
    [
        ("Homo sapiens", {"ARF4"}),
        ("human", {"ARF4"}),
        ("Mus musculus", {"ARF2", "ARF4"}),
        ("mouse", {"ARF2", "ARF4"}),
        ("Rattus norvegicus", {"ARF2", "ARF4"}),
        (None, {"ARF2", "ARF4"}),
    ],
)
def test_welch_fallback_applies_hgnc_only_for_explicit_human(
    tmp_path: Path,
    species: str | None,
    expected_symbols: set[str],
) -> None:
    matrix = tmp_path / "matrix.csv"
    pd.DataFrame(
        {
            "gene": ["Arf2", "Arf4"],
            "c1": [1.0, 2.0],
            "c2": [1.2, 2.2],
            "t1": [2.0, 3.0],
            "t2": [2.2, 3.2],
        }
    ).to_csv(matrix, index=False)
    output = tmp_path / "derived.csv"

    summary = derive_welch_deg(
        matrix,
        output,
        role="normalized_expression_matrix",
        gene_column="gene",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        normalized_scale="log2",
        species=species,
    )

    assert set(pd.read_csv(output)["gene_symbol"]) == expected_symbols
    assert summary["species"] == str(species or "")


@pytest.mark.parametrize(
    ("species", "expected_symbols"),
    [
        ("human", {"ARF4"}),
        ("mouse", {"ARF2", "ARF4"}),
    ],
)
def test_discovery_fallback_passes_its_scientific_species_to_reanalysis(
    tmp_path: Path,
    species: str,
    expected_symbols: set[str],
) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    matrix = bundle / "matrix.csv"
    pd.DataFrame(
        {
            "gene": ["Arf2", "Arf4"],
            "c1": [1.0, 2.0],
            "c2": [1.2, 2.2],
            "t1": [2.0, 3.0],
            "t2": [2.2, 3.2],
        }
    ).to_csv(matrix, index=False)
    candidate = {
        "candidate_id": "matrix1",
        "name": matrix.name,
        "role": "normalized_expression_matrix",
        "source_url": "https://example.test/matrix.csv",
        "inspection": {
            "status": "upstream_matrix_ready_for_contrast",
            "fetch_scope": "full",
            "local_path": str(matrix),
            "full_file_sha256": hashlib.sha256(matrix.read_bytes()).hexdigest(),
            "sample_columns": ["c1", "c2", "t1", "t2"],
            "gene_column": "gene",
            "header_row": 1,
        },
    }
    entry = {
        "candidate_id": "matrix1",
        "mode": "fallback",
        "direction_confirmed": True,
        "biological_replicates_confirmed": True,
        "control_samples": ["c1", "c2"],
        "treatment_samples": ["t1", "t2"],
        "normalized_scale": "log2",
        "gene_column": "gene",
        "contrast_label": "treatment versus control",
    }

    row, summary = _fallback_row(
        study={"accession": "GSE123456", "files": [candidate]},
        candidate=candidate,
        entry=entry,
        spec=normalize_species(species),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora discovery-analyze",
    )

    assert set(pd.read_csv(row["source_path"])["gene_symbol"]) == expected_symbols
    assert summary["species"] == normalize_species(species).scientific_name
