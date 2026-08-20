from __future__ import annotations

import re
import sys
import types

import pytest

from degora.cli import main


def _snapshot() -> dict:
    return {
        "query": "hypoxia",
        "species": "Human",
        "records": [
            {
                "canonical_id": "pmid:111",
                "source_unit_id": "PMID:111",
                "pubmed_ids": ["111"],
                "doi": "10.1000/example",
                "geo_accessions": ["GSE111"],
                "paper_title": "First",
                "data_readiness": {"verification_state": "verified_ready"},
            },
            {
                "canonical_id": "doi:10.2000/example",
                "source_unit_id": "SRC-2",
                "pubmed_ids": ["222"],
                "doi": "10.2000/example",
                "geo_accessions": ["GSE222"],
                "paper_title": "Second",
                "data_readiness": {"verification_state": "likely_ready"},
            },
        ],
    }


def test_federated_discover_is_default_and_exports_full_snapshot(tmp_path, monkeypatch, capsys) -> None:
    import degora.discovery_export as discovery_export
    import degora.discovery_federated as discovery_federated

    captured = {}

    def fake_search_publications(query, species, limit):
        captured["search"] = {"query": query, "species": species, "limit": limit}
        return _snapshot()

    def fake_export(snapshot, output_dir, *, force=False):
        captured["export"] = {"snapshot": snapshot, "output_dir": output_dir, "force": force}
        return {"search_csv": str(output_dir / "publication_search.csv")}

    monkeypatch.setattr(discovery_federated, "search_publications", fake_search_publications)
    monkeypatch.setattr(discovery_export, "export_publication_search", fake_export)

    output = tmp_path / "search"
    assert main(["discover", "hypoxia", "--species", "human", "--page", "1", "--limit", "37", "--output-dir", str(output)]) == 0

    assert captured["search"] == {"query": "hypoxia", "species": "human", "limit": 37}
    assert captured["export"]["snapshot"]["records"][0]["canonical_id"] == "pmid:111"
    assert captured["export"]["output_dir"] == output.resolve()
    out = capsys.readouterr().out
    assert "federated human PubMed+linked-data snapshot" in out
    assert "page size is 20" in out
    assert "pmid:111" in out
    assert "First" in out
    assert "verified ready" in out


def test_federated_discover_page_two_uses_global_row_numbers(tmp_path, monkeypatch, capsys) -> None:
    import degora.discovery_export as discovery_export
    import degora.discovery_federated as discovery_federated

    snapshot = _snapshot()
    snapshot["records"] = [
        {
            "canonical_id": f"pmid:{index}",
            "source_unit_id": f"PMID:{index}",
            "pubmed_ids": [str(index)],
            "paper_title": f"Publication {index}",
            "data_readiness": {"verification_state": "likely_ready"},
        }
        for index in range(1, 26)
    ]
    monkeypatch.setattr(
        discovery_federated,
        "search_publications",
        lambda query, species, limit: snapshot,
    )
    monkeypatch.setattr(
        discovery_export,
        "export_publication_search",
        lambda snapshot, output_dir, *, force=False: {"search_csv": str(output_dir / "publication_search.csv")},
    )

    assert main(
        [
            "discover",
            "hypoxia",
            "--species",
            "human",
            "--page",
            "2",
            "--output-dir",
            str(tmp_path / "search-page2"),
        ]
    ) == 0

    numbered_rows = [line for line in capsys.readouterr().out.splitlines() if re.match(r"\s*\d+\.", line)]
    assert numbered_rows
    assert numbered_rows[0].lstrip().startswith("21.")


def test_federated_network_unavailability_is_reported_without_traceback(tmp_path, monkeypatch, capsys) -> None:
    import degora.discovery_federated as discovery_federated
    from degora.discovery import DiscoveryUnavailableError

    def unavailable(*args, **kwargs):
        raise DiscoveryUnavailableError("NCBI is temporarily unavailable")

    monkeypatch.setattr(discovery_federated, "search_publications", unavailable)

    assert main(["discover", "hypoxia", "--species", "human", "--output-dir", str(tmp_path / "search")]) == 2
    assert "NCBI is temporarily unavailable" in capsys.readouterr().err


def test_federated_select_matches_identifiers_and_calls_prepare_backend(tmp_path, monkeypatch, capsys) -> None:
    import degora.discovery_federated as discovery_federated

    captured = {}

    def fake_prepare_publication_records(records, species, *, query, max_files_per_record, materialize_dir, force):
        captured.update(
            records=records,
            species=species,
            query=query,
            max_files_per_record=max_files_per_record,
            materialize_dir=materialize_dir,
            force=force,
        )
        return {
            "returned_records": len(records),
            "exports": {"draft_catalog_csv": str(materialize_dir / "DEGORA_discovery_draft_catalog.csv")},
        }

    prepare_module = types.ModuleType("degora.discovery_prepare")
    prepare_module.prepare_publication_records = fake_prepare_publication_records
    monkeypatch.setitem(sys.modules, "degora.discovery_prepare", prepare_module)
    monkeypatch.setattr(discovery_federated, "search_publications", lambda query, species, limit: _snapshot())

    output = tmp_path / "prepared"
    assert main(
        [
            "discover",
            "hypoxia",
            "--species",
            "human",
            "--select",
            "PMID:111",
            "--select",
            "GSE222",
            "--inspection-budget",
            "9",
            "--output-dir",
            str(output),
        ]
    ) == 0

    assert [record["canonical_id"] for record in captured["records"]] == ["pmid:111", "pmid:222"]
    assert captured["species"] == "human"
    assert captured["max_files_per_record"] == 4
    assert captured["materialize_dir"] == output.resolve()
    assert "Review required before analysis" in capsys.readouterr().out


def test_federated_select_rejects_missing_and_duplicate_matches(tmp_path, monkeypatch, capsys) -> None:
    import degora.discovery_federated as discovery_federated

    prepare_module = types.ModuleType("degora.discovery_prepare")
    prepare_module.prepare_publication_records = lambda *args, **kwargs: {}
    monkeypatch.setitem(sys.modules, "degora.discovery_prepare", prepare_module)
    monkeypatch.setattr(discovery_federated, "search_publications", lambda query, species, limit: _snapshot())

    assert main(["discover", "hypoxia", "--species", "human", "--select", "PMID:999", "--output-dir", str(tmp_path / "a")]) == 2
    assert "did not match" in capsys.readouterr().err

    assert main(
        [
            "discover",
            "hypoxia",
            "--species",
            "human",
            "--select",
            "PMID:111",
            "--select",
            "pmid:111",
            "--output-dir",
            str(tmp_path / "b"),
        ]
    ) == 2
    assert "duplicate --select" in capsys.readouterr().err


def test_legacy_geo_mode_keeps_geo_search_contract(tmp_path, monkeypatch, capsys) -> None:
    import degora.discovery as discovery

    captured = {}

    def fake_search_geo(query, species, *, page, page_size, assess_files, global_rank, global_limit):
        captured.update(
            query=query,
            species=species,
            page=page,
            page_size=page_size,
            assess_files=assess_files,
            global_rank=global_rank,
            global_limit=global_limit,
        )
        return {"total_hits": 5, "page": page, "evaluated_studies": 5, "studies": []}

    monkeypatch.setattr(discovery, "search_geo", fake_search_geo)
    monkeypatch.setattr(discovery, "export_search_page", lambda result, output, *, force=False: {"search_csv": str(output / "geo.csv")})

    assert main(["discover", "hypoxia", "--source", "geo", "--species", "human", "--limit", "12", "--output-dir", str(tmp_path / "geo")]) == 0

    assert captured["page_size"] == 20
    assert captured["global_rank"] is True
    assert captured["global_limit"] == 12
    assert "legacy GEO globally ranked" in capsys.readouterr().out


def test_discover_limit_is_validated_at_parser_boundary(tmp_path) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["discover", "hypoxia", "--species", "human", "--limit", "1001", "--output-dir", str(tmp_path)])
    assert exc_info.value.code == 2
