from __future__ import annotations

import sys
import types

import pytest

from degora.discovery import DiscoveryUnavailableError
from degora.discovery_federated import (
    canonical_record_id,
    merge_publication_records,
    page_publication_snapshot,
    rank_publication_records,
    resolve_publication_records,
    search_publications,
)


def test_identifier_graph_deduplicates_connected_publication_records() -> None:
    records = [
        {"pmid": "10", "doi": "10.1/A", "paper_title": "A", "species": "Human", "provider": "pubmed", "provider_id": "p10"},
        {"doi": "https://doi.org/10.1/a", "pmcid": "PMC10", "paper_title": "A full", "species": "Human", "provider": "pmc", "provider_id": "pmc10"},
        {"pmcid": "10", "accession": "GSE10", "paper_title": "A", "species": "Human", "provider": "geo", "provider_id": "gse10"},
    ]

    merged = merge_publication_records(records, "human")

    assert len(merged) == 1
    row = merged[0]
    assert row["canonical_id"] == "pmid:10"
    assert row["pmid"] == "10"
    assert row["doi"] == "10.1/a"
    assert row["pmcid"] == "PMC10"
    assert row["geo_accessions"] == ["GSE10"]
    assert row["provider_ids"] == ["geo:gse10", "pmc:pmc10", "pubmed:p10"]
    assert row["sources"] == ["geo", "pmc", "pubmed"]


def test_canonical_record_id_normalizes_doi_before_accession() -> None:
    row = {"doi": " DOI: https://doi.org/10.1093/NAR/GKAA000. ", "accession": "GSE1"}

    assert canonical_record_id(row) == "doi:10.1093/nar/gkaa000"


def test_no_geo_publication_row_is_preserved_as_metadata_only() -> None:
    merged = merge_publication_records(
        [{"pmid": "123", "paper_title": "Metadata only", "species": "Human", "relevance_rank": 2}],
        "Human",
    )

    assert len(merged) == 1
    assert merged[0]["geo_accessions"] == []
    assert merged[0]["accession"] == ""
    assert merged[0]["source_unit_id"] == "PMID:123"
    assert merged[0]["data_readiness"]["verification_state"] == "metadata_only"


def test_mixed_species_records_are_quarantined_unless_target_file_verified() -> None:
    mixed = merge_publication_records(
        [{"pmid": "1", "accession": "GSE1", "species_evidence": [{"species": "Human"}, {"species": "Mouse"}]}],
        "Human",
    )[0]
    rescued = merge_publication_records(
        [
            {
                "pmid": "2",
                "accession": "GSE2",
                "species_evidence": [{"species": "Human"}, {"species": "Mouse"}],
                "target_species_verified": True,
                "target_species_evidence": "This downloadable file contains Homo sapiens samples only.",
                "candidates": [{"name": "human_DEG.csv", "role": "deg_table", "url": "https://zenodo.org/human_DEG.csv"}],
            }
        ],
        "Human",
    )[0]

    assert mixed["species_decision"] == "mixed_quarantined"
    assert mixed["data_readiness"]["verification_state"] == "mixed_quarantined"
    assert rescued["species_decision"] == "mixed_rescued"
    assert rescued["mixed_rescued"] is True
    assert rescued["data_readiness"]["verification_state"] == "verified_ready"


def test_human_plus_primate_is_quarantined_and_unsubstantiated_rescue_is_rejected() -> None:
    records = [
        {
            "pmid": "3",
            "accession": "GSE3",
            "species_evidence": [{"species": "Homo sapiens"}, {"species": "Macaca mulatta"}],
            "target_species_verified": True,
        }
    ]

    row = merge_publication_records(records, "Human")[0]

    assert row["species"] == ["Human", "Other:Macaca mulatta"]
    assert row["species_decision"] == "mixed_quarantined"
    assert row["mixed_rescued"] is False
    assert row["data_readiness"]["verification_state"] == "mixed_quarantined"


def test_accession_without_verified_file_is_likely_not_content_verified() -> None:
    row = merge_publication_records(
        [{"pmid": "4", "accession": "GSE4", "species": "Human", "target_species_verified": True}],
        "Human",
    )[0]

    assert row["data_readiness"]["verification_state"] == "likely_ready"


def test_rank_120_ready_item_enters_first_top_20_after_global_sort() -> None:
    records = [
        {
            "canonical_id": f"pmid:{index}",
            "pmid": str(index),
            "paper_title": f"paper {index}",
            "year": 2020,
            "relevance_rank": index,
            "data_readiness": {"verification_state": "metadata_only", "priority": 4},
        }
        for index in range(1, 121)
    ]
    records[-1]["data_readiness"] = {"verification_state": "verified_ready", "priority": 0}

    page = page_publication_snapshot({"records": records}, page=1, page_size=20)

    assert len(page["records"]) == 20
    assert page["records"][0]["canonical_id"] == "pmid:120"
    assert any(row["canonical_id"] == "pmid:120" for row in page["records"])


def test_page_publication_snapshot_sorts_globally_for_header_sort() -> None:
    snapshot = {
        "records": [
            {"canonical_id": "pmid:1", "year": 2020},
            {"canonical_id": "pmid:2", "year": 2024},
            {"canonical_id": "pmid:3", "year": 2022},
        ]
    }

    page = page_publication_snapshot(snapshot, page=1, page_size=2, sort_by="year", sort_order="desc")

    assert [row["canonical_id"] for row in page["records"]] == ["pmid:2", "pmid:3"]
    assert page["total_records"] == 3
    assert page["page_size"] == 2


def test_deterministic_tie_order_uses_canonical_id() -> None:
    records = [
        {"canonical_id": "pmid:3", "year": 2020, "relevance_rank": 1, "data_readiness": {"priority": 1}},
        {"canonical_id": "pmid:1", "year": 2020, "relevance_rank": 1, "data_readiness": {"priority": 1}},
        {"canonical_id": "pmid:2", "year": 2020, "relevance_rank": 1, "data_readiness": {"priority": 1}},
    ]

    first = rank_publication_records(records)
    second = rank_publication_records(list(reversed(records)))

    assert [row["canonical_id"] for row in first] == ["pmid:1", "pmid:2", "pmid:3"]
    assert [row["canonical_id"] for row in second] == ["pmid:1", "pmid:2", "pmid:3"]


def test_search_publications_caps_evaluated_records_at_1000() -> None:
    class Provider:
        name = "fixture"

        def search(self, query, species, limit):
            return [
                {
                    "pmid": str(index),
                    "paper_title": f"{query} {index}",
                    "species": species.label,
                    "relevance_rank": index,
                    "accession": f"GSE{index}",
                }
                for index in range(1, limit + 1)
            ]

        def resolve(self, record, species):
            return [{**record, "target_species_verified": True}]

    snapshot = search_publications("ifn", "Human", limit=1200, providers=[Provider()])

    assert snapshot["diagnostics"]["evaluated_limit"] == 1000
    assert snapshot["diagnostics"]["evaluated_records"] == 1000
    assert snapshot["total_records"] == 1000
    assert snapshot["records"][0]["data_readiness"]["verification_state"] == "likely_ready"


def test_search_publications_lazily_imports_default_providers(monkeypatch: pytest.MonkeyPatch) -> None:
    class Provider:
        name = "lazy"

        def search(self, query, species, limit):
            return [{"pmid": "7", "species": species.label, "paper_title": query}]

    module = types.ModuleType("degora.discovery_sources")
    module.default_publication_providers = lambda: [Provider()]
    monkeypatch.setitem(sys.modules, "degora.discovery_sources", module)

    snapshot = search_publications("tnfa", "mouse", providers=None)

    assert snapshot["species"]["key"] == "mouse"
    assert snapshot["species"]["label"] == "Mouse"
    assert snapshot["records"][0]["canonical_id"] == "pmid:7"
    assert snapshot["diagnostics"]["provider_count"] == 1


def test_search_publications_searches_all_providers_even_when_first_fills_cap() -> None:
    calls: list[str] = []

    class Provider:
        def __init__(self, name: str, start: int) -> None:
            self.name = name
            self.start = start

        def search(self, query, species, limit):
            calls.append(self.name)
            return [{"pmid": str(index), "species": species.label} for index in range(self.start, self.start + limit)]

    snapshot = search_publications("ifn", "Human", limit=3, providers=[Provider("first", 1), Provider("second", 100)])

    assert calls == ["first", "second"]
    assert snapshot["diagnostics"]["evaluated_records"] == 3
    assert snapshot["diagnostics"]["searched"][1]["truncated"] is False
    assert snapshot["diagnostics"]["searched"][1]["global_cap_competition"] is True
    assert snapshot["total_records"] == 3
    assert snapshot["provider_events"] == [
        {"provider": "first", "event": "search", "status": "ok", "message": "3 records"},
        {"provider": "second", "event": "search", "status": "ok", "message": "3 records"},
    ]


def test_search_publications_fails_when_every_search_provider_is_unavailable() -> None:
    class FailingProvider:
        name = "failed"

        def search(self, query, species, limit):
            raise OSError("request failed?api_key=do-not-leak")

    with pytest.raises(DiscoveryUnavailableError, match="all publication search providers") as exc_info:
        search_publications("ifn", "human", providers=[FailingProvider()])

    assert "do-not-leak" not in str(exc_info.value)


def test_search_publications_keeps_partial_results_and_redacts_provider_secrets() -> None:
    class FailingProvider:
        name = "failed"

        def search(self, query, species, limit):
            raise RuntimeError("https://service.test?q=x&api_key=do-not-leak")

    class WorkingProvider:
        name = "working"

        def search(self, query, species, limit):
            return [{"pmid": "7", "title": query, "species": species.label}]

    snapshot = search_publications("ifn", "human", providers=[FailingProvider(), WorkingProvider()])

    assert snapshot["provider_status"] == "partial"
    assert snapshot["total"] == 1
    assert "do-not-leak" not in str(snapshot["diagnostics"])


def test_selected_records_are_re_resolved_and_keep_failure_diagnostics() -> None:
    calls: list[str] = []

    class Resolver:
        name = "resolver"

        def resolve(self, record, species):
            calls.append(record["pmid"])
            return [
                {
                    "record_type": "diagnostics",
                    "diagnostics": [
                        {"provider": "repository", "status": "candidate_error", "error": "temporary outage"}
                    ],
                }
            ]

    record = {
        "pmid": "11",
        "species": "Human",
        "candidates": [{"url": "https://zenodo.org/existing.csv", "role": "deg_table"}],
    }
    resolved = resolve_publication_records([record], "human", providers=[Resolver()])

    assert calls == ["11"]
    assert resolved[0]["resolution_state"] == "unavailable"
    assert resolved[0]["resolution_events"][0]["status"] == "candidate_error"


def test_candidate_merge_preserves_distinct_evidence_routes_for_one_accession() -> None:
    merged = merge_publication_records(
        [
            {
                "pmid": "12",
                "accession": "GSE12",
                "species": "Human",
                "candidates": [
                    {
                        "accession": "GSE12",
                        "provider": "geo",
                        "url": "https://ftp.ncbi.nlm.nih.gov/author_DEG.csv",
                        "name": "author_DEG.csv",
                        "role": "deg_table",
                    },
                    {
                        "accession": "GSE12",
                        "provider": "geo",
                        "url": "https://ftp.ncbi.nlm.nih.gov/count_matrix.csv",
                        "name": "count_matrix.csv",
                        "role": "count_matrix",
                    },
                ],
            }
        ],
        "human",
    )

    assert {candidate.get("role") for candidate in merged[0]["candidates"]} >= {"deg_table", "count_matrix"}


def test_actual_provider_shapes_merge_and_repository_resolver_attaches_no_geo_candidate() -> None:
    class PubMedShape:
        name = "ncbi_pubmed"

        def search(self, query, species, limit):
            return [
                {
                    "rank": 1,
                    "pmid": "11",
                    "title": "No GEO paper",
                    "species_evidence": {
                        "requested": species.scientific_name,
                        "basis": "PubMed organism-constrained query",
                    },
                },
                {
                    "rank": 2,
                    "pmid": "22",
                    "title": "Linked GEO paper",
                    "species_evidence": {
                        "requested": species.scientific_name,
                        "basis": "PubMed organism-constrained query",
                    },
                },
            ]

    class GeoShape:
        name = "ncbi_geo"

        def search(self, query, species, limit):
            return [
                {
                    "rank": 1,
                    "pmid": "22",
                    "pmids": ["22"],
                    "accession": "GSE22",
                    "species_evidence": {
                        "requested": species.scientific_name,
                        "observed_taxa": [species.scientific_name],
                        "status": "exact",
                        "basis": "GEO SOFT observed taxa",
                    },
                    "supplementary_file_candidates": [
                        {
                            "url": "https://ftp.ncbi.nlm.nih.gov/GSE22_DESeq2_results.csv.gz",
                            "name": "GSE22_DESeq2_results.csv.gz",
                            "role": "deg_table",
                        }
                    ],
                }
            ]

    class ResolverShape:
        name = "public_repository_resolver"

        def resolve(self, record, species):
            if record.get("pmid") != "11":
                return []
            return [
                {
                    "record_type": "file_candidate",
                    "provider": "europe_pmc",
                    "url": "https://www.ebi.ac.uk/files/no_geo_DEG.csv",
                    "name": "no_geo_DEG.csv",
                    "role": "deg_table",
                }
            ]

    snapshot = search_publications(
        "hypoxia",
        "human",
        limit=10,
        providers=[PubMedShape(), GeoShape(), ResolverShape()],
    )

    by_pmid = {record["pmid"]: record for record in snapshot["records"]}
    assert set(by_pmid) == {"11", "22"}
    assert by_pmid["11"]["candidates"][0]["source_url"].endswith("no_geo_DEG.csv")
    assert by_pmid["11"]["data_readiness"]["verification_state"] == "likely_ready"
    assert by_pmid["22"]["geo_accessions"] == ["GSE22"]
    assert any(candidate.get("role") == "deg_table" for candidate in by_pmid["22"]["candidates"])
    assert snapshot["species"]["key"] == "human"


def test_provider_mixed_taxa_dictionary_is_quarantined_after_cross_source_merge() -> None:
    records = merge_publication_records(
        [
            {"pmid": "9", "species_evidence": {"requested": "Homo sapiens", "basis": "query"}},
            {
                "pmid": "9",
                "accession": "GSE9",
                "species_evidence": {
                    "observed_taxa": ["Homo sapiens", "Mus musculus"],
                    "basis": "GEO SOFT observed taxa",
                },
                "quarantined": True,
            },
        ],
        "human",
    )

    assert records[0]["species_decision"] == "mixed_quarantined"
    assert records[0]["data_readiness"]["verification_state"] == "mixed_quarantined"


def test_search_publications_reports_monotonic_progress_stages() -> None:
    """The browser draws a determinate bar from these fractions."""

    class _Provider:
        name = "stub"

        def search(self, query, species, limit):
            return [
                {"canonical_id": f"PMID:{index}", "pmid": str(index), "title": f"paper {index}"}
                for index in range(1, 4)
            ]

        def resolve(self, record, species):
            return []

    seen: list[tuple[float, str]] = []
    snapshot = search_publications(
        query="hypoxia",
        species="human",
        limit=10,
        providers=[_Provider()],
        progress=lambda fraction, message: seen.append((fraction, message)),
    )

    assert snapshot["total"] >= 1
    assert seen, "no progress was reported"
    fractions = [fraction for fraction, _ in seen]
    assert fractions == sorted(fractions), f"progress must never decrease: {fractions}"
    assert all(0.0 <= fraction <= 1.0 for fraction in fractions)
    messages = " | ".join(message for _, message in seen)
    assert "Querying stub" in messages
    assert "Inspecting linked data" in messages


def test_search_publications_without_progress_callback_is_unchanged() -> None:
    class _Provider:
        name = "stub"

        def search(self, query, species, limit):
            return [{"canonical_id": "PMID:1", "pmid": "1", "title": "paper"}]

        def resolve(self, record, species):
            return []

    with_callback = search_publications(
        query="hypoxia", species="human", limit=5, providers=[_Provider()], progress=lambda *_: None
    )
    without_callback = search_publications(query="hypoxia", species="human", limit=5, providers=[_Provider()])
    assert [row["canonical_id"] for row in with_callback["records"]] == [
        row["canonical_id"] for row in without_callback["records"]
    ]


def test_search_publications_survives_a_failing_progress_callback() -> None:
    class _Provider:
        name = "stub"

        def search(self, query, species, limit):
            return [{"canonical_id": "PMID:1", "pmid": "1", "title": "paper"}]

        def resolve(self, record, species):
            return []

    def boom(_fraction, _message):
        raise RuntimeError("progress sink exploded")

    snapshot = search_publications(
        query="hypoxia", species="human", limit=5, providers=[_Provider()], progress=boom
    )
    assert snapshot["total"] == 1


def test_search_marks_resolution_state_so_prepare_can_reuse_it() -> None:
    class _Provider:
        name = "stub"

        def search(self, query, species, limit):
            return [{"canonical_id": "PMID:7", "pmid": "7", "title": "paper"}]

        def resolve(self, record, species):
            return []

    snapshot = search_publications(query="hypoxia", species="human", limit=5, providers=[_Provider()])
    states = {row.get("resolution_state") for row in snapshot["records"]}
    assert "resolved_no_candidate" in states


def test_resolve_publication_records_skips_publications_search_already_settled() -> None:
    """Preparing a page-one selection used to repeat every provider call."""

    calls: list[str] = []

    class _Provider:
        name = "stub"

        def resolve(self, record, species):
            calls.append(str(record.get("canonical_id")))
            return []

    settled = {"canonical_id": "PMID:1", "pmid": "1", "title": "one", "resolution_state": "resolved_candidates"}
    fresh = {"canonical_id": "PMID:2", "pmid": "2", "title": "two"}

    resolved = resolve_publication_records([settled, fresh], "human", providers=[_Provider()])

    assert [call.lower() for call in calls] == ["pmid:2"], (
        f"only the unresolved publication should be re-resolved, got {calls}"
    )
    ids = {str(row["canonical_id"]).lower() for row in resolved if row.get("canonical_id")}
    assert {"pmid:1", "pmid:2"} <= ids, "reusing a settled record must not drop it"


def test_resolve_publication_records_can_force_a_fresh_resolution() -> None:
    calls: list[str] = []

    class _Provider:
        name = "stub"

        def resolve(self, record, species):
            calls.append(str(record.get("canonical_id")))
            return []

    settled = {"canonical_id": "PMID:1", "pmid": "1", "title": "one", "resolution_state": "resolved_candidates"}
    resolve_publication_records([settled], "human", providers=[_Provider()], reuse_settled=False)
    assert [call.lower() for call in calls] == ["pmid:1"]


def test_unavailable_resolution_is_retried_rather_than_reused() -> None:
    calls: list[str] = []

    class _Provider:
        name = "stub"

        def resolve(self, record, species):
            calls.append(str(record.get("canonical_id")))
            return []

    flaky = {"canonical_id": "PMID:9", "pmid": "9", "title": "nine", "resolution_state": "unavailable"}
    resolve_publication_records([flaky], "human", providers=[_Provider()])
    assert [call.lower() for call in calls] == ["pmid:9"]
