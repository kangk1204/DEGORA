from __future__ import annotations

import csv
import gzip
import http.client
import json
import urllib.error
import urllib.parse
from pathlib import Path

import pytest

from degora.discovery import (
    DiscoveryError,
    DiscoveryUnavailableError,
    GLOBAL_SEARCH_CACHE_SIZE,
    GLOBAL_SEARCH_CACHE_TTL_SECONDS,
    NcbiGeoClient,
    NcbiRequestConfig,
    PUBMED_SUMMARY_BATCH_SIZE,
    SafeNcbiTransport,
    build_geo_query,
    classify_filename,
    export_discovery_bundle,
    export_search_page,
    inspect_candidate_bytes,
    inspect_upstream_bytes,
    normalize_ncbi_url,
    normalize_species,
    parse_geo_soft,
    prepare_geo_studies,
    search_geo,
)


HUMAN_RECORD = {
    "accession": "GSE100001",
    "taxon": "Homo sapiens",
    "title": "Dataset title",
    "summary": "Dataset summary",
    "gdstype": "Expression profiling by high throughput sequencing",
    "pdat": "2024/03/02",
    "n_samples": "6",
    "pubmedids": ["12345"],
}


class FakeGeoClient:
    def __init__(self) -> None:
        self.search_calls: list[tuple[str, str, int, int]] = []
        self.fetch_calls: list[tuple[str, bool]] = []
        self.soft_calls: list[str] = []

    def search_summaries(self, query, species, *, page, page_size):
        self.search_calls.append((query, species.key, page, page_size))
        return 42, [dict(HUMAN_RECORD), {**HUMAN_RECORD, "accession": "GSE100002", "taxon": "Mus musculus"}]

    def accession_summaries(self, accessions, species):
        records = [dict(HUMAN_RECORD)]
        if "GSE100002" in accessions:
            records.append({**HUMAN_RECORD, "accession": "GSE100002", "taxon": "Homo sapiens"})
        return records

    def publication_summaries(self, pmids):
        return {
            "12345": {
                "title": "Paper title",
                "authors": [{"name": "Kim K"}, {"name": "Lee J"}],
                "fulljournalname": "Genome Biology",
                "pubdate": "2023 Dec",
                "sortpubdate": "2023/12/01 00:00",
            }
        }

    def fetch_geo_soft(self, accession):
        self.soft_calls.append(accession)
        if accession == "GSE100002":
            return "\n".join(
                [
                    "^SERIES = GSE100002",
                    "!Series_sample_organism_ch1 = Homo sapiens",
                    "!Series_sample_organism_ch1 = Mus musculus",
                ]
            )
        return "\n".join(
            [
                "^SERIES = GSE100001",
                "!Series_title = Dataset title",
                "!Series_overall_design = treatment versus control",
                "!Series_sample_organism_ch1 = Homo sapiens",
                "!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100001/suppl/GSE100001_DESeq2_results.csv.gz",
                "!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100001/suppl/GSE100001_raw_counts.tsv.gz",
            ]
        )

    def fetch_candidate(self, url, *, full):
        self.fetch_calls.append((url, full))
        if "DESeq2" in url:
            payload = gzip.compress(
                b"gene,log2FoldChange,pvalue,padj\nTP53,2.0,0.001,0.01\nCDKN1A,1.2,0.01,0.04\n"
            )
        else:
            payload = gzip.compress(
                b"gene\tctrl_1\tctrl_2\ttreat_1\ttreat_2\nTP53\t10\t12\t40\t45\nCDKN1A\t4\t5\t30\t28\n"
            )
        return payload, "full" if full else "header_prefix"


def test_species_contract_rejects_cross_species_pooling() -> None:
    assert normalize_species("Human").scientific_name == "Homo sapiens"
    assert normalize_species("Mus musculus").key == "mouse"
    with pytest.raises(DiscoveryError, match="exactly human or mouse"):
        normalize_species("both")


def test_geo_query_keeps_mandatory_species_and_neutralizes_entrez_syntax() -> None:
    query = build_geo_query('hypoxia OR mouse[Organism]', "human")
    assert '"Homo sapiens"[Organism]' in query
    assert '"Mus musculus"[Organism]' not in query
    assert '"hypoxia"[All Fields]' in query


def test_ncbi_url_normalization_converts_only_ncbi_ftp() -> None:
    assert normalize_ncbi_url("ftp://ftp.ncbi.nlm.nih.gov/a.tsv") == "https://ftp.ncbi.nlm.nih.gov/a.tsv"
    with pytest.raises(DiscoveryError, match="host is not allowed"):
        normalize_ncbi_url("https://example.org/a.tsv")
    with pytest.raises(DiscoveryError, match="HTTPS"):
        normalize_ncbi_url("http://www.ncbi.nlm.nih.gov/a.tsv")


def test_transport_rejects_declared_response_larger_than_cap() -> None:
    class Response:
        headers = {"Content-Length": "999999"}

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def geturl(self):
            return "https://www.ncbi.nlm.nih.gov/a.tsv"

        def read(self, _size):
            return b""

    class Opener:
        def open(self, _request, timeout):
            assert timeout > 0
            return Response()

    transport = SafeNcbiTransport()
    transport._opener = Opener()
    with pytest.raises(DiscoveryError, match="safety cap"):
        transport("https://www.ncbi.nlm.nih.gov/a.tsv", {}, 5, 1024)


def test_filename_classifier_separates_deg_and_upstream_tables() -> None:
    assert classify_filename("GSE1_DESeq2_results.csv.gz")["role"] == "deg_table"
    assert classify_filename("GSE1_raw_count_matrix.tsv.gz")["role"] == "count_matrix"
    assert classify_filename("GSE1_series_matrix.txt.gz")["role"] == "normalized_expression_matrix"
    assert classify_filename("GSE1_RAW.tar")["tier"] == "reject"


def test_deg_header_and_values_must_be_usable() -> None:
    good = b"gene,log2FoldChange,pvalue,padj\nA,1.2,0.01,0.03\nB,-2.0,0.02,0.04\n"
    result = inspect_candidate_bytes("result.csv", good)
    assert result["status"] == "ready_for_review"
    assert result["mapping"]["gene_column"] == "gene"

    cuffdiff = (
        b"gene_id\tgene\tlog2(fold_change)\tp_value\tq_value\n"
        b"A\t-\t1.2\t0.01\t0.03\nB\t-\t-2.0\t0.02\t0.04\n"
    )
    cuffdiff_result = inspect_candidate_bytes("gene_exp.diff.tsv", cuffdiff)
    assert cuffdiff_result["status"] == "ready_for_review"
    assert cuffdiff_result["mapping"]["gene_column"] == "gene_id"

    ambiguous = b"gene,foldChange,pvalue\nA,2,0.01\nB,3,0.02\n"
    assert inspect_candidate_bytes("result.csv", ambiguous)["status"] == "requires_lfc_confirmation"

    invalid_p = b"gene,log2FoldChange,pvalue\nA,2,12\nB,3,20\n"
    assert inspect_candidate_bytes("result.csv", invalid_p)["status"] != "ready_for_review"

    multiple_contrasts = (
        b"gene,contrast_A_log2FoldChange,contrast_B_log2FoldChange,contrast_A_pvalue,contrast_B_pvalue\n"
        b"A,2,1,0.01,0.02\nB,-1,-2,0.03,0.04\n"
    )
    ambiguous_mapping = inspect_candidate_bytes("multi_contrast.csv", multiple_contrasts)
    assert ambiguous_mapping["status"] == "requires_column_mapping"
    assert len(ambiguous_mapping["lfc_columns"]) == 2
    assert len(ambiguous_mapping["p_columns"]) == 2


def test_deg_header_prefers_hugo_symbol_over_generic_gene_identifier() -> None:
    payload = (
        b"gene,hugo,logFC,P.Value,adj.P.Val\n"
        b"ENSG000001,A,1.2,0.01,0.03\nENSG000002,B,-2.0,0.02,0.04\n"
    )

    result = inspect_candidate_bytes("author_results.csv", payload)

    assert result["status"] == "ready_for_review"
    assert result["mapping"]["gene_column"] == "hugo"


def test_large_gzip_deg_table_uses_bounded_header_prefix_instead_of_false_rejection() -> None:
    rows = b"gene\tlog2FoldChange\tpvalue\tpadj\nA\t1.2\t0.01\t0.03\nB\t-2.0\t0.02\t0.04\n"
    payload = gzip.compress(rows + (b"C\t0.1\t0.5\t0.8\n" * 400_000))
    result = inspect_candidate_bytes("large_result.tsv.gz", payload)
    assert result["status"] == "ready_for_review"


def test_upstream_matrix_requires_four_numeric_samples() -> None:
    payload = (
        b"gene\tGene Start\tGene Stop\tc1\tc2\tt1\tt2\n"
        b"A\t100\t200\t10\t12\t20\t22\nB\t300\t400\t2\t4\t8\t9\n"
    )
    result = inspect_upstream_bytes("raw_counts.tsv", payload, declared_role="count_matrix")
    assert result["status"] == "upstream_matrix_ready_for_contrast"
    assert result["gene_column"] == "gene"
    assert result["sample_columns"] == ["c1", "c2", "t1", "t2"]


def test_upstream_inspection_does_not_promote_numeric_data_row_to_header() -> None:
    payload = (
        b"Gene ID\tGene Start\tGene Stop\tctrl_1\tctrl_2\ttreat_1\ttreat_2\n"
        b"88191914\t88192060\t88192110\t4\t7\t9\t1\n"
        b"88193000\t88193100\t88193200\t2\t3\t8\t6\n"
    )
    result = inspect_upstream_bytes("processed.tsv", payload, declared_role="unknown_matrix")
    assert result["gene_column"] == "Gene ID"
    assert result["sample_columns"] == ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]


def test_upstream_inspection_prefers_gene_symbol_over_identifier_and_excludes_coordinates() -> None:
    payload = (
        b"Gene ID\tGene Symbol\tGene Start\tGene Stop\tctrl_1\tctrl_2\ttreat_1\ttreat_2\n"
        b"ENSG1\tA\t100\t200\t4\t7\t9\t11\nENSG2\tB\t300\t400\t2\t3\t8\t6\n"
    )
    result = inspect_upstream_bytes("processed.tsv", payload, declared_role="count_matrix")
    assert result["gene_column"] == "Gene Symbol"
    assert result["sample_columns"] == ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]


def test_unknown_statistical_table_is_not_promoted_to_expression_matrix() -> None:
    payload = (
        b"gene\tAveExpr\tt\tadj.P.Val\tB\n"
        b"A\t6.1\t4.2\t0.01\t3.0\nB\t5.4\t-3.1\t0.02\t2.1\n"
    )
    result = inspect_upstream_bytes("processed.tsv", payload, declared_role="unknown_matrix")
    assert result["status"] == "not_upstream_matrix"
    assert result["sample_columns"] == []
    assert set(result["statistic_columns"]) == {"AveExpr", "t", "adj.P.Val", "B"}
    assert "statistic columns" in result["reason"]


def test_search_geo_returns_twenty_item_page_contract_and_publication_metadata() -> None:
    client = FakeGeoClient()
    result = search_geo("hypoxia", "human", page=2, page_size=20, client=client)
    assert client.search_calls == [("hypoxia", "human", 2, 20)]
    assert result["total_hits"] == 42
    assert result["degora_version"]
    if "degora_code_revision" in result:
        assert result["degora_code_revision"]
    assert result["page"] == 2
    assert result["page_size"] == 20
    assert result["total_pages"] is None
    assert result["total_pages_upper_bound"] == 3
    assert result["has_next"] is False
    assert result["returned_studies"] == 1  # mismatched mouse result is rejected
    study = result["studies"][0]
    assert study["paper_title"] == "Paper title"
    assert study["authors_display"] == "Kim K, Lee J"
    assert study["journal"] == "Genome Biology"
    assert study["year"] == 2023
    assert study["source_unit_id"] == "PMID:12345"
    assert study["ncbi_relevance_rank"] == 21
    assert "deg_input_assessment" not in study
    assert client.soft_calls == []
    assert "paper_title" in result["sortable_fields"]


def test_search_geo_keeps_version_provenance_when_git_revision_is_unavailable(monkeypatch) -> None:
    import degora.discovery as discovery

    monkeypatch.setattr(discovery, "runtime_version_info", lambda: {"degora_version": "0.4.0"})

    result = search_geo("hypoxia", "human", page=1, page_size=20, client=FakeGeoClient())

    assert result["degora_version"] == "0.4.0"
    assert "degora_code_revision" not in result


def test_search_assessment_ranks_likely_deg_inputs_and_preserves_ncbi_order() -> None:
    class RankedSearchClient(FakeGeoClient):
        def search_summaries(self, query, species, *, page, page_size):
            self.search_calls.append((query, species.key, page, page_size))
            return 5, [
                {**HUMAN_RECORD, "accession": f"GSE10000{index}", "pubmedids": []}
                for index in range(1, 6)
            ]

        def publication_summaries(self, pmids):
            return {}

        def fetch_geo_soft(self, accession):
            self.soft_calls.append(accession)
            if accession == "GSE100005":
                raise DiscoveryUnavailableError("forced SOFT failure")
            filenames = {
                "GSE100001": ["GSE100001_series_matrix.txt.gz"],
                "GSE100002": ["GSE100002_DESeq2_results.csv.gz"],
                "GSE100003": ["GSE100003_processed_table.tsv.gz"],
                "GSE100004": ["GSE100004_RAW.tar"],
            }[accession]
            return "\n".join(
                [f"^SERIES = {accession}"]
                + [f"!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/{name}" for name in filenames]
            )

        def fetch_candidate(self, url, *, full):  # pragma: no cover - must never be called
            raise AssertionError("search assessment downloaded a candidate payload")

    client = RankedSearchClient()
    result = search_geo("hypoxia", "human", page=1, page_size=5, assess_files=True, client=client)
    studies = result["studies"]

    assert [study["accession"] for study in studies] == [
        "GSE100002",
        "GSE100003",
        "GSE100001",
        "GSE100004",
        "GSE100005",
    ]
    assert [study["ncbi_relevance_rank"] for study in studies] == [2, 3, 1, 4, 5]
    assert [study["deg_input_assessment"]["tier"] for study in studies] == [
        "author_deg_likely",
        "tabular_candidate",
        "matrix_fallback",
        "not_detected",
        "unresolved",
    ]
    assert studies[-1]["deg_input_assessment"]["error"] == "SOFT request failed (DiscoveryUnavailableError)"
    serialized = json.dumps(result).lower()
    assert "ready_for_review" not in serialized
    assert "confirmed deg" not in serialized
    assert result["search_assessment"]["default_sort"] == "deg_input_priority_desc"


def test_global_search_ranks_before_pagination_caps_at_1000_and_reuses_cache(monkeypatch) -> None:
    client = NcbiGeoClient(transport=lambda *args: b"", pace_seconds=0, cache_size=0)
    records = []
    for exact_rank in range(1, 1002):
        accession = f"GSE{200000 + exact_rank}"
        if exact_rank in {10, 500}:
            records.append({**HUMAN_RECORD, "accession": f"GSE{800000 + exact_rank}", "taxon": "Homo sapiens; Mus musculus", "pubmedids": []})
        record = {**HUMAN_RECORD, "accession": accession, "pubmedids": []}
        records.append(record)
        if exact_rank == 11:
            records.append(dict(record))

    gds_calls = []
    soft_calls = []

    def fake_gds_search(term, *, retstart, retmax):
        gds_calls.append((retstart, retmax))
        return len(records), records[retstart : retstart + retmax]

    def fake_publication_summaries(pmids):
        assert list(pmids) == []
        return {}

    def fake_fetch_geo_soft(accession):
        soft_calls.append(accession)
        filename = f"{accession}_DESeq2_results.csv.gz" if accession in {"GSE200120", "GSE201001"} else ""
        lines = [f"^SERIES = {accession}"]
        if filename:
            lines.append(f"!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/{filename}")
        return "\n".join(lines)

    monkeypatch.setattr(client, "_gds_search", fake_gds_search)
    monkeypatch.setattr(client, "publication_summaries", fake_publication_summaries)
    monkeypatch.setattr(client, "fetch_geo_soft", fake_fetch_geo_soft)
    monkeypatch.setattr(
        client,
        "fetch_candidate",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("search downloaded a candidate payload")),
    )

    first = search_geo(
        "hypoxia",
        "human",
        page=1,
        page_size=20,
        assess_files=True,
        global_rank=True,
        global_limit=1000,
        client=client,
    )
    first["studies"][0]["accession"] = "GSE0"  # callers must not mutate the cached snapshot
    page_two = search_geo(
        "hypoxia",
        "human",
        page=2,
        page_size=20,
        assess_files=True,
        global_rank=True,
        global_limit=1000,
        client=client,
    )
    first_again = search_geo(
        "hypoxia",
        "human",
        page=1,
        page_size=20,
        assess_files=True,
        global_rank=True,
        global_limit=1000,
        client=client,
    )

    assert first["cache_hit"] is False
    assert page_two["cache_hit"] is True
    assert first_again["cache_hit"] is True
    assert first_again["studies"][0]["accession"] == "GSE200120"
    assert first_again["studies"][0]["ncbi_relevance_rank"] == 120
    assert [study["ncbi_relevance_rank"] for study in page_two["studies"]] == list(range(20, 40))
    assert first_again["evaluated_studies"] == 1000
    assert first_again["ranking_limit"] == 1000
    assert first_again["ranking_truncated"] is True
    assert first_again["ranking_scope"] == "global_first_1000_exact_single_organism_studies"
    assert first_again["total_pages"] == 50
    assert first_again["has_next"] is True
    assert len(soft_calls) == 1000
    assert len(set(soft_calls)) == 1000
    assert "GSE201001" not in soft_calls
    assert gds_calls


def test_global_search_header_sort_applies_before_page_slice_and_keeps_relevance_ties(monkeypatch) -> None:
    client = NcbiGeoClient(transport=lambda *args: b"", pace_seconds=0, cache_size=0)
    records = [
        {**HUMAN_RECORD, "accession": f"GSE{300000 + index}", "pubmedids": [], "pdat": f"{year}/01/01"}
        for index, year in enumerate([2020, 2024, 2024, 2021, 2023], start=1)
    ]

    monkeypatch.setattr(
        client,
        "_gds_search",
        lambda term, *, retstart, retmax: (len(records), records[retstart : retstart + retmax]),
    )
    monkeypatch.setattr(client, "publication_summaries", lambda pmids: {})
    monkeypatch.setattr(client, "fetch_geo_soft", lambda accession: f"^SERIES = {accession}")

    result = search_geo(
        "hypoxia",
        "human",
        page=1,
        page_size=2,
        assess_files=True,
        global_rank=True,
        global_limit=1000,
        sort_by="year",
        sort_order="desc",
        client=client,
    )

    assert [study["year"] for study in result["studies"]] == [2024, 2024]
    assert [study["ncbi_relevance_rank"] for study in result["studies"]] == [2, 3]
    assert result["sort_by"] == "year"
    assert result["sort_order"] == "desc"

    relevance = search_geo(
        "hypoxia",
        "human",
        page=1,
        page_size=2,
        assess_files=True,
        global_rank=True,
        global_limit=1000,
        sort_by="relevance",
        client=client,
    )
    assert [study["ncbi_relevance_rank"] for study in relevance["studies"]] == [1, 2]
    assert relevance["sort_order"] == "asc"
    assert relevance["cache_hit"] is True
    assert relevance["ranking_truncated"] is False


def test_global_snapshot_cache_is_success_only_ttl_bounded_and_lru_bounded() -> None:
    clock = [100.0]
    client = NcbiGeoClient(
        transport=lambda *args: b"",
        pace_seconds=0,
        cache_size=0,
        monotonic=lambda: clock[0],
    )
    failing_key = ("failure", "human", 1000, 1)
    failures = 0

    def failing_builder():
        nonlocal failures
        failures += 1
        raise DiscoveryUnavailableError("forced build failure")

    with pytest.raises(DiscoveryUnavailableError, match="forced build failure"):
        client.get_or_build_global_search_snapshot(failing_key, failing_builder)
    with pytest.raises(DiscoveryUnavailableError, match="forced build failure"):
        client.get_or_build_global_search_snapshot(failing_key, failing_builder)
    assert failures == 2

    builds = 0

    def successful_builder():
        nonlocal builds
        builds += 1
        return {"studies": [{"accession": f"GSE{builds}"}]}

    key = ("hypoxia", "human", 1000, 1)
    first, first_hit = client.get_or_build_global_search_snapshot(key, successful_builder)
    first["studies"][0]["accession"] = "GSE0"
    cached, cached_hit = client.get_or_build_global_search_snapshot(key, successful_builder)
    assert first_hit is False
    assert cached_hit is True
    assert cached["studies"][0]["accession"] == "GSE1"

    clock[0] += GLOBAL_SEARCH_CACHE_TTL_SECONDS + 1
    expired, expired_hit = client.get_or_build_global_search_snapshot(key, successful_builder)
    assert expired_hit is False
    assert expired["studies"][0]["accession"] == "GSE2"

    for index in range(GLOBAL_SEARCH_CACHE_SIZE + 1):
        cache_key = (f"query-{index}", "human", 1000, 1)
        client.get_or_build_global_search_snapshot(cache_key, lambda index=index: {"index": index})
    assert len(client._global_search_cache) == GLOBAL_SEARCH_CACHE_SIZE


def test_pubmed_summaries_are_requested_in_bounded_batches(monkeypatch) -> None:
    client = NcbiGeoClient(transport=lambda *args: b"", pace_seconds=0, cache_size=0)
    batches = []

    def fake_get_json(url, *, max_bytes=8 * 1024 * 1024):
        del max_bytes
        ids = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)["id"][0].split(",")
        batches.append(ids)
        return {"result": {pmid: {"uid": pmid} for pmid in ids}}

    monkeypatch.setattr(client, "get_json", fake_get_json)
    pmids = [str(100000 + index) for index in range(450)]
    summaries = client.publication_summaries([*pmids, pmids[0]])

    assert [len(batch) for batch in batches] == [PUBMED_SUMMARY_BATCH_SIZE, PUBMED_SUMMARY_BATCH_SIZE, 50]
    assert list(summaries) == pmids


def test_search_export_flattens_assessment_and_escapes_formula_like_filename(tmp_path: Path) -> None:
    result = {
        "studies": [
            {
                "species": "human",
                "accession": "GSE1",
                "paper_title": "Paper",
                "authors": [],
                "pubmed_ids": [],
                "ncbi_relevance_rank": 1,
                "deg_input_assessment": {
                    "tier": "tabular_candidate",
                    "label": "Table to inspect",
                    "priority": 3,
                    "basis": "filename likelihood only",
                    "candidate_files": ["=DEG_results.csv"],
                    "counts": {"deg_like": 0, "tabular": 1, "matrix": 0},
                    "error": "",
                },
            }
        ]
    }

    export_search_page(result, tmp_path)
    with (tmp_path / "geo_search_page.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["deg_input_likelihood_tier"] == "tabular_candidate"
    assert row["deg_input_candidate_filenames"] == "'=DEG_results.csv"
    assert "deg_input_assessment" not in row


def test_ncbi_exact_species_pagination_fills_pages_without_mixed_study_duplicates(monkeypatch) -> None:
    client = NcbiGeoClient()
    records = [
        {
            **HUMAN_RECORD,
            "accession": f"GSE{100000 + index}",
            "taxon": "Homo sapiens; Mus musculus" if index == 15 else "Homo sapiens",
        }
        for index in range(1, 46)
    ]

    def fake_gds_search(term, *, retstart, retmax):
        return len(records), records[retstart : retstart + retmax]

    monkeypatch.setattr(client, "_gds_search", fake_gds_search)
    species = normalize_species("human")
    _, first = client.search_summaries("hypoxia", species, page=1, page_size=20)
    _, second = client.search_summaries("hypoxia", species, page=2, page_size=20)

    assert len(first) == 20
    assert len(second) == 20
    assert all(record["taxon"] == "Homo sapiens" for record in first + second)
    assert {record["accession"] for record in first}.isdisjoint(record["accession"] for record in second)


def test_global_search_diagnostics_record_mixed_organism_exclusion(monkeypatch) -> None:
    client = NcbiGeoClient()
    records = [
        {**HUMAN_RECORD, "accession": "GSE108676", "taxon": "Macaca mulatta; Homo sapiens"},
        {**HUMAN_RECORD, "accession": "GSE225253", "taxon": "Homo sapiens"},
    ]

    def fake_gds_search(term, *, retstart, retmax):
        return len(records), records[retstart : retstart + retmax]

    monkeypatch.setattr(client, "_gds_search", fake_gds_search)
    total, retained, truncated, diagnostics = client.search_summaries_global_detailed(
        "hypoxia",
        normalize_species("human"),
        limit=1000,
    )

    assert total == 2
    assert [record["accession"] for record in retained] == ["GSE225253"]
    assert truncated is False
    assert diagnostics["excluded_reason_counts"] == {"mixed_or_mismatched_organism": 1}
    assert diagnostics["excluded_records_sample"] == [
        {
            "accession": "GSE108676",
            "observed_taxon": "Macaca mulatta; Homo sapiens",
            "reason": "mixed_or_mismatched_organism",
        }
    ]


def test_global_search_diagnostics_distinguish_fetched_processed_and_overflow(monkeypatch) -> None:
    client = NcbiGeoClient()
    records = [
        {**HUMAN_RECORD, "accession": f"GSE{index:06d}", "taxon": "Homo sapiens"}
        for index in range(1, 101)
    ]

    def fake_gds_search(term, *, retstart, retmax):
        return len(records), records[retstart : retstart + retmax]

    monkeypatch.setattr(client, "_gds_search", fake_gds_search)
    _, retained, truncated, diagnostics = client.search_summaries_global_detailed(
        "hypoxia",
        normalize_species("human"),
        limit=1,
    )

    assert [record["accession"] for record in retained] == ["GSE000001"]
    assert truncated is True
    assert diagnostics["raw_records_fetched"] == 100
    assert diagnostics["raw_records_scanned"] == 2
    assert diagnostics["raw_records_processed"] == 2
    assert diagnostics["exact_records_retained"] == 1
    assert diagnostics["exact_records_overflow"] == 1
    assert diagnostics["excluded_record_count"] == 0
    assert diagnostics["records_accounted"] == diagnostics["raw_records_scanned"]


def test_prepare_selected_study_materializes_candidates_but_keeps_catalog_inactive(tmp_path: Path) -> None:
    client = FakeGeoClient()
    result = prepare_geo_studies(
        ["GSE100001", "GSE100002"],
        "human",
        query="hypoxia",
        materialize_dir=tmp_path,
        client=client,
    )
    assert result["returned_studies"] == 1
    assert result["excluded_studies"] == [
        {"accession": "GSE100002", "reason": "mixed or mismatched organism metadata: Homo sapiens, Mus musculus"}
    ]
    study = result["studies"][0]
    assert study["preparation_status"] == "author_deg_ready_for_contrast_review"
    assert study["ready_for_review_count"] == 1
    assert study["upstream_matrix_count"] == 1
    assert all(full for _, full in client.fetch_calls)
    assert list(tmp_path.glob("GSE100001_*"))

    catalog = (tmp_path / "DEGORA_discovery_draft_catalog.csv").read_text()
    assert "include_in_analysis" in catalog
    assert ",no," in catalog
    audit = json.loads((tmp_path / "discovery_audit.json").read_text())
    assert audit["analysis_policy"]["cross_species_pooling"] is False


def test_failed_materialized_preparation_publishes_no_partial_bundle(tmp_path: Path) -> None:
    class FailingSecondClient(FakeGeoClient):
        def accession_summaries(self, accessions, species):
            return [
                {**HUMAN_RECORD, "accession": accession, "pubmedids": [str(12000 + index)]}
                for index, accession in enumerate(accessions, start=1)
            ]

        def fetch_geo_soft(self, accession):
            return "\n".join(
                [
                    f"^SERIES = {accession}",
                    "!Series_sample_organism_ch1 = Homo sapiens",
                    f"!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/{accession}_DESeq2_results.csv.gz",
                ]
            )

        def fetch_candidate(self, url, *, full):
            self.fetch_calls.append((url, full))
            if len(self.fetch_calls) == 2:
                raise OSError("forced second download failure")
            return gzip.compress(b"gene,log2FoldChange,pvalue\nTP53,2.0,0.001\n"), "full"

    target = tmp_path / "prepared"
    with pytest.raises(OSError, match="forced second download failure"):
        prepare_geo_studies(
            ["GSE100001", "GSE100002"],
            "human",
            materialize_dir=target,
            client=FailingSecondClient(),
        )

    assert not target.exists()
    assert not list(tmp_path.glob(".prepared.prepare-*"))


def test_failed_forced_preparation_preserves_previous_complete_bundle(tmp_path: Path) -> None:
    target = tmp_path / "prepared"
    prepare_geo_studies(
        ["GSE100001"],
        "human",
        materialize_dir=target,
        client=FakeGeoClient(),
    )
    before = {
        path.relative_to(target): path.read_bytes()
        for path in target.rglob("*")
        if path.is_file()
    }

    class FailingReplacementClient(FakeGeoClient):
        def fetch_candidate(self, url, *, full):
            raise OSError("forced replacement download failure")

    with pytest.raises(OSError, match="forced replacement download failure"):
        prepare_geo_studies(
            ["GSE100001"],
            "human",
            materialize_dir=target,
            client=FailingReplacementClient(),
            force=True,
        )

    after = {
        path.relative_to(target): path.read_bytes()
        for path in target.rglob("*")
        if path.is_file()
    }
    assert after == before
    assert not list(tmp_path.glob(".prepared.prepare-*"))


def test_export_bundle_refuses_overwrite_without_force(tmp_path: Path) -> None:
    result = {"studies": [], "species": {"key": "human"}}
    export_discovery_bundle(result, tmp_path)
    with pytest.raises(FileExistsError, match="--force"):
        export_discovery_bundle(result, tmp_path)


def test_search_csv_neutralizes_spreadsheet_formulas_but_json_preserves_source_text(tmp_path: Path) -> None:
    title = '=HYPERLINK("https://example.invalid","open")'
    result = {
        "studies": [
            {
                "accession": "GSE1",
                "paper_title": title,
                "authors": ["@SUM(1+1)"],
                "pubmed_ids": ["123"],
            }
        ]
    }

    exports = export_search_page(result, tmp_path)
    with Path(exports["search_csv"]).open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    payload = json.loads(Path(exports["search_json"]).read_text(encoding="utf-8"))

    assert row["paper_title"] == "'" + title
    assert row["authors"].startswith("'@")
    assert payload["studies"][0]["paper_title"] == title


def test_ncbi_client_caches_success_but_never_failure() -> None:
    success_calls = 0

    def success(url, headers, timeout, max_bytes):
        nonlocal success_calls
        success_calls += 1
        return b"ok"

    client = NcbiGeoClient(transport=success, pace_seconds=0, retries=1)
    assert client.get_bytes("https://www.ncbi.nlm.nih.gov/a", max_bytes=10) == b"ok"
    assert client.get_bytes("https://www.ncbi.nlm.nih.gov/a", max_bytes=10) == b"ok"
    assert success_calls == 1

    failure_calls = 0

    def failure(url, headers, timeout, max_bytes):
        nonlocal failure_calls
        failure_calls += 1
        raise urllib.error.URLError("offline")

    failing = NcbiGeoClient(transport=failure, pace_seconds=0, retries=1)
    for _ in range(2):
        with pytest.raises(DiscoveryUnavailableError):
            failing.get_bytes("https://www.ncbi.nlm.nih.gov/a", max_bytes=10)
    assert failure_calls == 2


def test_ncbi_client_retries_incomplete_chunked_response_without_caching_partial_bytes() -> None:
    calls = 0
    sleeps: list[float] = []

    def transient_chunk_failure(url, headers, timeout, max_bytes):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise http.client.IncompleteRead(b"partial")
        return b"complete"

    client = NcbiGeoClient(
        transport=transient_chunk_failure,
        pace_seconds=0,
        retries=2,
        cache_size=2,
        sleep=sleeps.append,
    )

    assert client.get_bytes("https://www.ncbi.nlm.nih.gov/a", max_bytes=20) == b"complete"
    assert client.get_bytes("https://www.ncbi.nlm.nih.gov/a", max_bytes=20) == b"complete"
    assert calls == 2
    assert sleeps == [2.0]


def test_ncbi_identity_is_applied_to_every_eutils_request_without_secret_leakage() -> None:
    captured: list[str] = []
    secret = "test_api_key_123"

    def transport(url, headers, timeout, max_bytes):
        captured.append(url)
        parsed = urllib.parse.urlsplit(url)
        params = urllib.parse.parse_qs(parsed.query)
        if parsed.path.endswith("esearch.fcgi"):
            return json.dumps({"esearchresult": {"count": "1", "idlist": ["1"]}}).encode()
        if params.get("db") == ["gds"]:
            return json.dumps({"result": {"1": HUMAN_RECORD}}).encode()
        return json.dumps({"result": {"12345": {"title": "Paper"}}}).encode()

    config = NcbiRequestConfig(tool="degora_test", email="researcher@example.org", api_key=secret)
    client = NcbiGeoClient(request_config=config, transport=transport, retries=1)
    result = search_geo("hypoxia", "human", page=1, page_size=20, client=client)

    assert client.pace_seconds == pytest.approx(0.11)
    assert config.request_ceiling_per_second == 10
    assert secret not in repr(config)
    assert secret not in json.dumps(result)
    assert len(captured) == 3
    for url in captured:
        params = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
        assert params["tool"] == ["degora_test"]
        assert params["email"] == ["researcher@example.org"]
        assert params["api_key"] == [secret]

    def failure(url, headers, timeout, max_bytes):
        raise urllib.error.URLError(f"failed URL: {url}")

    failing = NcbiGeoClient(request_config=config, transport=failure, pace_seconds=0, retries=1)
    with pytest.raises(DiscoveryUnavailableError) as error:
        failing._gds_search("hypoxia", retstart=0, retmax=1)
    assert secret not in str(error.value)


def test_email_identifies_requests_but_does_not_raise_the_rate_ceiling() -> None:
    email_only = NcbiRequestConfig(email="researcher@example.org")
    assert email_only.pace_seconds == pytest.approx(0.35)
    assert email_only.request_ceiling_per_second == 3
    with pytest.raises(DiscoveryError, match="request ceiling"):
        NcbiGeoClient(request_config=email_only, pace_seconds=0)


def test_soft_metadata_deduplicates_files_and_shared_secondary_pmids_define_one_source_unit(tmp_path: Path) -> None:
    class SharedPmidClient(FakeGeoClient):
        def accession_summaries(self, accessions, species):
            return [
                {**HUMAN_RECORD, "accession": "GSE100001", "pubmedids": ["111"]},
                {**HUMAN_RECORD, "accession": "GSE100002", "pubmedids": ["222"]},
            ]

        def publication_summaries(self, pmids):
            return {
                str(pmid): {
                    "title": f"Paper {pmid}",
                    "authors": [{"name": "Kim K"}],
                    "fulljournalname": "Genome Biology",
                    "pubdate": "2024",
                }
                for pmid in pmids
            }

        def fetch_geo_soft(self, accession):
            url = f"ftp://ftp.ncbi.nlm.nih.gov/{accession}_DESeq2_results.csv.gz"
            return "\n".join(
                [
                    f"^SERIES = {accession}",
                    "!Series_sample_organism_ch1 = Homo sapiens",
                    "!Series_pubmed_id = 999",
                    f"!Series_supplementary_file = {url}",
                    f"!Series_supplementary_file = {url}",
                ]
            )

    result = prepare_geo_studies(
        ["GSE100001", "GSE100002"],
        "human",
        inspection_budget=0,
        materialize_dir=tmp_path,
        client=SharedPmidClient(),
    )

    assert result["returned_studies"] == 2
    assert {study["source_unit_id"] for study in result["studies"]} == {"PMID:111"}
    assert all(study["source_unit_pubmed_ids"] == ["111", "222", "999"] for study in result["studies"])
    for study in result["studies"]:
        supplemental = [item for item in study["files"] if "DESeq2_results" in item["source_url"]]
        assert len(supplemental) == 1


def test_parse_geo_soft_deduplicates_exact_repeated_metadata() -> None:
    text = "\n".join(
        [
            "!Series_pubmed_id = 123",
            "!Series_pubmed_id = 123",
            "!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/a.tsv.gz",
            "!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/a.tsv.gz",
        ]
    )
    parsed = parse_geo_soft(text)
    assert parsed["pubmed_ids"] == ["123"]
    assert parsed["supplementary_files"] == ["https://ftp.ncbi.nlm.nih.gov/a.tsv.gz"]


# --- matrix columns are not always GSM accessions --------------------------


def _labels() -> dict[str, dict]:
    return {
        "GSM6072341": {"title": "4641CERM6M24M", "source": "mammary gland", "characteristics": ["transgene: induced"]},
        "GSM6072342": {"title": "4709CERM6m24M", "source": "mammary gland", "characteristics": ["transgene: induced"]},
        "GSM6072343": {"title": "control 12M", "source": "", "characteristics": ["transgene: uninduced"]},
    }


def test_submitter_column_names_are_matched_to_their_geo_sample() -> None:
    """An author matrix is headed by the submitter's own names, not accessions.

    Keying only on GSM left every column of such a matrix unlabelled, which is
    the state a reader was asked to assign control and treatment from.
    """

    from degora.discovery import match_sample_labels

    resolved = match_sample_labels(
        ["4641CERM6M24M_S2", "4709CERM6m24M_S2", "GSM6072343", "Gene"],
        _labels(),
    )

    assert resolved["4641CERM6M24M_S2"]["accession"] == "GSM6072341"
    # Case and a trailing lane suffix must not defeat the match.
    assert resolved["4709CERM6m24M_S2"]["accession"] == "GSM6072342"
    assert resolved["GSM6072343"]["accession"] == "GSM6072343"
    assert "Gene" not in resolved
    assert resolved["4641CERM6M24M_S2"]["characteristics"] == ["transgene: induced"]


def test_an_ambiguous_column_is_left_unmatched_rather_than_guessed() -> None:
    """A wrong sample label flips a contrast, so silence beats a guess."""

    from degora.discovery import match_sample_labels

    labels = {
        "GSM1": {"title": "wild type", "source": "liver", "characteristics": []},
        "GSM2": {"title": "wild type", "source": "liver", "characteristics": []},
    }

    assert match_sample_labels(["wild type", "liver"], labels) == {}


def test_short_or_empty_column_names_never_match() -> None:
    from degora.discovery import match_sample_labels

    labels = {"GSM1": {"title": "A1", "source": "", "characteristics": []}}

    assert match_sample_labels(["A1", "", "  "], labels) == {}


def test_workbook_candidates_are_inspectable_not_unsupported() -> None:
    """DEGORA read .xls and .xlsx all along; discovery refused to look at them.

    A supplementary workbook - the commonest shape after plain text, and often
    gzipped - was classified unsupported and never offered as a candidate, so the
    search declined files the analysis path could have used.
    """

    from degora.discovery import classify_filename

    for name in ("deg.xlsx", "deg.xls", "deg.xlsx.gz", "deg.xls.gz"):
        assessment = classify_filename(name)
        assert assessment["inspectable"] is True, name
        assert assessment["tier"] != "unsupported", name

    assert classify_filename("paper.pdf")["inspectable"] is False


def test_a_workbook_candidate_is_classified_from_its_contents(tmp_path) -> None:
    """Including the legacy and gzipped shapes, and refusing a file that only claims to be one."""

    import gzip
    import shutil

    from openpyxl import Workbook

    from degora.discovery import DiscoveryError, inspect_candidate_bytes

    book = Workbook()
    sheet = book.active
    sheet.append(["gene", "log2FoldChange", "pvalue"])
    for index in range(1, 6):
        sheet.append([f"G{index}", 2.0 - index * 0.1, 0.001 * index])
    plain = tmp_path / "deg.xlsx"
    book.save(plain)
    packed = tmp_path / "deg.xlsx.gz"
    with plain.open("rb") as source, gzip.open(packed, "wb") as target:
        shutil.copyfileobj(source, target)

    for path in (plain, packed):
        result = inspect_candidate_bytes(path.name, path.read_bytes())
        assert result["status"] == "ready_for_review", path.name
        assert result["mapping"]["lfc_column"] == "log2FoldChange"

    # A renamed text file must not be inspected as a workbook.
    with pytest.raises(DiscoveryError, match="OLE2 signature"):
        inspect_candidate_bytes("renamed.xls", b"gene,lfc\nA,1\n")
    with pytest.raises(DiscoveryError, match="could not be expanded"):
        inspect_candidate_bytes("broken.xlsx.gz", b"not gzip at all")
