from __future__ import annotations

import hashlib
import io
import urllib.error
import urllib.parse
import zipfile
from pathlib import Path

import pytest

from degora.discovery import DiscoveryError, DiscoveryUnavailableError, NcbiRequestConfig, normalize_species
from degora.discovery_sources import (
    NcbiGeoProvider,
    NcbiPubmedProvider,
    PublicRepositoryResolver,
    SafePublicTransport,
    describe_unexpected_payload,
    download_public_candidate,
    inspect_public_archive,
    _link_looks_public_repository,
)


PUBLIC_DNS = [(None, None, None, "", ("8.8.8.8", 443))]
PRIVATE_DNS = [(None, None, None, "", ("127.0.0.1", 443))]


class Response:
    def __init__(self, url: str, payload: bytes, *, status: int = 200, headers: dict[str, str] | None = None) -> None:
        self._url = url
        self._payload = io.BytesIO(payload)
        self.status = status
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def geturl(self) -> str:
        return self._url

    def read(self, size: int = -1) -> bytes:
        return self._payload.read(size)


class MappingOpener:
    def __init__(self, responses: dict[str, Response]) -> None:
        self.responses = responses
        self.urls: list[str] = []

    def open(self, request, timeout: int):
        assert timeout > 0
        url = request.full_url
        self.urls.append(url)
        if url not in self.responses:
            raise AssertionError(f"unexpected URL: {url}")
        return self.responses[url]


class SequenceOpener:
    def __init__(self, responses: list[Response | Exception]) -> None:
        self.responses = list(responses)
        self.urls: list[str] = []

    def open(self, request, timeout: int):
        assert timeout > 0
        self.urls.append(request.full_url)
        if not self.responses:
            raise AssertionError("unexpected extra request")
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


def public_transport(responses: dict[str, Response]) -> SafePublicTransport:
    return SafePublicTransport(resolver=lambda *_args, **_kwargs: PUBLIC_DNS, opener=MappingOpener(responses), sleep=lambda _x: None)


def test_safe_public_transport_rejects_ssrf_url_forms_and_dns() -> None:
    transport = SafePublicTransport(resolver=lambda *_args, **_kwargs: PUBLIC_DNS, opener=MappingOpener({}))
    with pytest.raises(DiscoveryError, match="HTTPS"):
        transport.validate_url("http://www.ncbi.nlm.nih.gov/a.csv")
    with pytest.raises(DiscoveryError, match="credentials"):
        transport.validate_url("https://user:pass@www.ncbi.nlm.nih.gov/a.csv")
    with pytest.raises(DiscoveryError, match="fragments"):
        transport.validate_url("https://www.ncbi.nlm.nih.gov/a.csv#frag")
    with pytest.raises(DiscoveryError, match="default HTTPS port"):
        transport.validate_url("https://www.ncbi.nlm.nih.gov:8443/a.csv")
    with pytest.raises(DiscoveryError, match="host is not allowed"):
        transport.validate_url("https://example.org/a.csv")

    dns_blocked = SafePublicTransport(resolver=lambda *_args, **_kwargs: PRIVATE_DNS, opener=MappingOpener({}))
    with pytest.raises(DiscoveryError, match="DNS result is not public"):
        dns_blocked.validate_url("https://www.ncbi.nlm.nih.gov/a.csv")


def test_safe_public_transport_validates_redirect_and_oversized_payload() -> None:
    redirected = public_transport(
        {
            "https://www.ncbi.nlm.nih.gov/start.csv": Response(
                "https://www.ncbi.nlm.nih.gov/start.csv",
                b"",
                status=302,
                headers={"Location": "https://example.org/blocked.csv"},
            )
        }
    )
    with pytest.raises(DiscoveryError, match="host is not allowed"):
        redirected.get_bytes("https://www.ncbi.nlm.nih.gov/start.csv", max_bytes=20)

    oversized = public_transport(
        {
            "https://ftp.ncbi.nlm.nih.gov/a.csv": Response(
                "https://ftp.ncbi.nlm.nih.gov/a.csv",
                b"abc",
                headers={"Content-Length": "999"},
            )
        }
    )
    with pytest.raises(DiscoveryError, match="safety cap"):
        oversized.get_bytes("https://ftp.ncbi.nlm.nih.gov/a.csv", max_bytes=10)


def test_safe_public_transport_normalizes_network_failures_without_secret_urls() -> None:
    class FailingOpener:
        def __init__(self) -> None:
            self.calls = 0

        def open(self, request, timeout: int):
            self.calls += 1
            raise urllib.error.URLError(f"timeout for {request.full_url}")

    opener = FailingOpener()
    transport = SafePublicTransport(
        resolver=lambda *_args, **_kwargs: PUBLIC_DNS,
        opener=opener,
        sleep=lambda _x: None,
    )
    url = "https://eutils.ncbi.nlm.nih.gov/test?api_key=do-not-leak"

    with pytest.raises(DiscoveryUnavailableError, match="public source request failed") as exc_info:
        transport.get_bytes(url)

    assert "do-not-leak" not in str(exc_info.value)
    assert url not in str(exc_info.value)
    assert opener.calls == 3


def test_safe_public_transport_retries_transient_http_failures() -> None:
    url = "https://eutils.ncbi.nlm.nih.gov/retry"
    opener = SequenceOpener(
        [
            urllib.error.HTTPError(
                url,
                429,
                "Too Many Requests",
                {"Retry-After": "0"},
                None,
            ),
            Response(url, b"recovered"),
        ]
    )
    sleeps: list[float] = []
    transport = SafePublicTransport(
        resolver=lambda *_args, **_kwargs: PUBLIC_DNS,
        opener=opener,
        sleep=sleeps.append,
    )

    assert transport.get_bytes(url) == b"recovered"
    assert opener.urls == [url, url]
    assert sleeps == [0.0]


def test_safe_public_transport_does_not_retry_nontransient_http_failures() -> None:
    url = "https://eutils.ncbi.nlm.nih.gov/missing"
    opener = SequenceOpener(
        [urllib.error.HTTPError(url, 404, "Not Found", {}, None)]
    )
    transport = SafePublicTransport(
        resolver=lambda *_args, **_kwargs: PUBLIC_DNS,
        opener=opener,
        sleep=lambda _x: None,
    )

    with pytest.raises(urllib.error.HTTPError) as exc_info:
        transport.get_bytes(url)

    assert exc_info.value.code == 404
    assert opener.urls == [url]


def test_ncbi_pubmed_provider_normalizes_batched_metadata() -> None:
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {
                "db": "pubmed",
                "term": (
                    '("hypoxia"[All Fields]) AND ("Humans"[MeSH Terms] OR "Homo sapiens"[All Fields] OR '
                    '"human"[Title/Abstract] OR "humans"[Title/Abstract])'
                ),
                "retmax": 2,
                "retmode": "json",
                "sort": "relevance",
                "usehistory": "y",
                "tool": "degora-test",
            }
        )
    )
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "id": "111,222", "retmode": "json", "tool": "degora-test"})
    )
    transport = public_transport(
        {
            search_url: Response(search_url, b'{"esearchresult":{"idlist":["111","222"]}}'),
            summary_url: Response(
                summary_url,
                (
                    b'{"result":{"111":{"title":"A","fulljournalname":"J","sortpubdate":"2024/01/01",'
                    b'"authors":[{"name":"Kim K"}],"articleids":[{"idtype":"doi","value":"10.1/a"},'
                    b'{"idtype":"pmcid","value":"PMC1"}]},"222":{"title":"B","articleids":[]}}}'
                ),
            ),
        }
    )
    provider = NcbiPubmedProvider(
        transport=transport,
        request_config=NcbiRequestConfig(tool="degora-test"),
        sleep=lambda _x: None,
    )
    records = provider.search("hypoxia", normalize_species("human"), 2)
    assert [record["pmid"] for record in records] == ["111", "222"]
    assert records[0]["rank"] == 1
    assert records[0]["doi"] == "10.1/a"
    assert records[0]["pmcid"] == "PMC1"
    assert records[0]["species_evidence"]["basis"] == "PubMed organism-constrained query"


def test_ncbi_pubmed_provider_hydrates_linked_geo_pmid_metadata() -> None:
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "id": "37949939", "retmode": "json", "tool": "degora-test"})
    )
    provider = NcbiPubmedProvider(
        transport=public_transport(
            {
                summary_url: Response(
                    summary_url,
                    (
                        b'{"result":{"37949939":{"title":"Hypoxia study",'
                        b'"fulljournalname":"Genome Biology","sortpubdate":"2023/11/13",'
                        b'"authors":[{"name":"Zhang Y"},{"name":"Kang K"}],'
                        b'"articleids":[{"idtype":"doi","value":"10.1038/example"}]}}}'
                    ),
                )
            }
        ),
        request_config=NcbiRequestConfig(tool="degora-test"),
        sleep=lambda _x: None,
    )

    hydrated = provider.resolve(
        {"pmid": "37949939", "title": "GEO title", "authors": [], "journal": ""},
        normalize_species("human"),
    )

    assert hydrated[0]["authors"] == ["Zhang Y", "Kang K"]
    assert hydrated[0]["journal"] == "Genome Biology"
    assert hydrated[0]["doi"] == "10.1038/example"
    assert "species_evidence" not in hydrated[0]


def test_ncbi_pubmed_provider_skips_hydration_for_complete_or_unlinked_record() -> None:
    provider = NcbiPubmedProvider(
        transport=public_transport({}),
        request_config=NcbiRequestConfig(tool="degora-test"),
        sleep=lambda _x: None,
    )
    complete = {
        "pmid": "1",
        "title": "Complete",
        "authors": ["A"],
        "journal": "J",
        "year": 2024,
    }

    assert provider.resolve(complete, normalize_species("human")) == []
    assert provider.resolve({"title": "No PMID"}, normalize_species("human")) == []


def test_ncbi_geo_provider_quarantines_mixed_records_and_returns_filename_candidates() -> None:
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {
                "db": "gds",
                "term": '("hypoxia"[All Fields]) AND gse[ETYP] AND "Homo sapiens"[Organism]',
                "retmax": 1,
                "retmode": "json",
                "sort": "relevance",
                "usehistory": "y",
                "tool": "degora-test",
            }
        )
    )
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "gds", "id": "999", "retmode": "json", "tool": "degora-test"})
    )
    soft_url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        + urllib.parse.urlencode({"acc": "GSE9", "targ": "self", "form": "text", "view": "brief"})
    )
    soft = "\n".join(
        [
            "^SERIES = GSE9",
            "!Series_title = Mixed record",
            "!Series_sample_organism_ch1 = Homo sapiens",
            "!Series_sample_organism_ch1 = Mus musculus",
            "!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/geo/GSE9_DESeq2_results.csv.gz",
        ]
    ).encode()
    provider = NcbiGeoProvider(
        transport=public_transport(
            {
                search_url: Response(search_url, b'{"esearchresult":{"idlist":["999"]}}'),
                summary_url: Response(summary_url, b'{"result":{"999":{"accession":"GSE9","pubmedids":["111"]}}}'),
                soft_url: Response(soft_url, soft),
            }
        ),
        request_config=NcbiRequestConfig(tool="degora-test"),
        sleep=lambda _x: None,
    )
    records = provider.search("hypoxia", normalize_species("human"), 1)
    assert records[0]["quarantined"] is True
    assert records[0]["species_evidence"]["status"] == "quarantined_mixed_record"
    assert records[0]["supplementary_file_candidates"][0]["role"] == "deg_table"
    assert records[0]["supplementary_file_candidates"][0]["downloaded"] is False


def test_ncbi_pubmed_provider_uses_history_server_for_summary_batches() -> None:
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {
                "db": "pubmed",
                "term": (
                    '("ifn"[All Fields]) AND ("Mice"[MeSH Terms] OR "Mus musculus"[All Fields] OR '
                    '"mouse"[Title/Abstract] OR "mice"[Title/Abstract] OR "murine"[Title/Abstract])'
                ),
                "retmax": 1,
                "retmode": "json",
                "sort": "relevance",
                "usehistory": "y",
                "tool": "degora-test",
            }
        )
    )
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode(
            {
                "db": "pubmed",
                "query_key": "7",
                "WebEnv": "history-token",
                "retstart": 0,
                "retmax": 1,
                "retmode": "json",
                "tool": "degora-test",
            }
        )
    )
    provider = NcbiPubmedProvider(
        transport=public_transport(
            {
                search_url: Response(
                    search_url,
                    b'{"esearchresult":{"idlist":["123"],"querykey":"7","webenv":"history-token"}}',
                ),
                summary_url: Response(summary_url, b'{"result":{"123":{"title":"History-backed","articleids":[]}}}'),
            }
        ),
        request_config=NcbiRequestConfig(tool="degora-test"),
        sleep=lambda _x: None,
    )

    records = provider.search("ifn", normalize_species("mouse"), 1)

    assert records[0]["pmid"] == "123"


def test_europe_pmc_resolver_exposes_bounded_supplementary_zip() -> None:
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
        + urllib.parse.urlencode(
            {"query": "EXT_ID:123 SRC:MED", "format": "json", "resultType": "core", "pageSize": 10}
        )
    )
    resolver = PublicRepositoryResolver(
        transport=public_transport(
            {
                url: Response(
                    url,
                    b'{"resultList":{"result":[{"pmcid":"PMC123","hasSuppl":"Y"}]}}',
                )
            }
        )
    )

    records = resolver._resolve_europe_pmc({"pmid": "123"}, normalize_species("human"))

    assert records == [
        {
            "provider": "europe_pmc",
            "record_type": "file_candidate",
            "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC123/supplementaryFiles?includeInlineImage=n",
            "name": "PMC123_supplementary_files.zip",
            "role": "archive",
            "tier": "archive_candidate",
            "downloaded": False,
        }
    ]


def test_biostudies_resolver_hydrates_public_tabular_files_from_exact_hit() -> None:
    search_url = (
        "https://www.ebi.ac.uk/biostudies/api/v1/search?"
        + urllib.parse.urlencode({"query": "PMC123", "pageSize": 3})
    )
    study_url = "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC123"
    info_url = study_url + "/info"
    resolver = PublicRepositoryResolver(
        transport=public_transport(
            {
                search_url: Response(
                    search_url,
                    b'{"hits":[{"accession":"S-EPMC123","isPublic":true}]}',
                ),
                study_url: Response(
                    study_url,
                    b'{"section":{"files":[{"path":"tables/author_DEG.csv"},{"path":"../unsafe.csv"}]}}',
                ),
                info_url: Response(
                    info_url,
                    b'{"httpLink":"https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC123"}',
                ),
            }
        )
    )

    records = resolver._resolve_biostudies({"pmcid": "PMC123"}, normalize_species("human"))

    assert len(records) == 1
    assert records[0]["url"].endswith("/Files/tables/author_DEG.csv")
    assert records[0]["role"] == "deg_table"


def test_biostudies_resolver_does_not_duplicate_files_directory() -> None:
    search_url = (
        "https://www.ebi.ac.uk/biostudies/api/v1/search?"
        + urllib.parse.urlencode({"query": "PMC123", "pageSize": 3})
    )
    study_url = "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC123"
    info_url = study_url + "/info"
    resolver = PublicRepositoryResolver(
        transport=public_transport(
            {
                search_url: Response(
                    search_url,
                    b'{"hits":[{"accession":"S-EPMC123","isPublic":true}]}',
                ),
                study_url: Response(
                    study_url,
                    b'{"section":{"files":[{"path":"Files/author_DEG.csv"}]}}',
                ),
                info_url: Response(
                    info_url,
                    b'{"httpLink":"https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC123/"}',
                ),
            }
        )
    )

    records = resolver._resolve_biostudies({"pmcid": "PMC123"}, normalize_species("human"))

    assert len(records) == 1
    assert records[0]["url"].endswith("/S-EPMC123/Files/author_DEG.csv")


def test_repository_resolver_converts_provider_failure_to_diagnostics() -> None:
    class FailingResolver(PublicRepositoryResolver):
        def _resolve_europe_pmc(self, record, species):
            raise RuntimeError("forced")

        def _resolve_crossref(self, record, species):
            return [{"provider": "crossref", "record_type": "file_candidate", "url": "https://zenodo.org/file.csv"}]

        def _resolve_datacite(self, record, species):
            return []

    results = FailingResolver(transport=public_transport({})).resolve(
        {"pmid": "1", "doi": "10.1/a"},
        normalize_species("human"),
    )
    assert results[0]["provider"] == "crossref"
    assert results[-1]["record_type"] == "diagnostics"
    assert results[-1]["diagnostics"][0]["status"] == "candidate_error"


def test_repository_metadata_404_is_an_empty_route_not_a_provider_failure() -> None:
    class MissingTransport:
        def get_json(self, url):
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    resolver = PublicRepositoryResolver(transport=MissingTransport())

    assert resolver._resolve_crossref({"doi": "10.1/missing"}, normalize_species("human")) == []
    assert resolver._resolve_datacite({"doi": "10.1/missing"}, normalize_species("human")) == []


def test_repository_link_filter_requires_a_supported_file_signal() -> None:
    assert _link_looks_public_repository("https://www.ebi.ac.uk/files/123", "author_DEG.csv") is True
    assert _link_looks_public_repository("https://www.ebi.ac.uk/files/123") is False
    assert _link_looks_public_repository("https://www.ebi.ac.uk/articles/PMC123") is False


def test_download_public_candidate_is_atomic_and_reports_sha256(tmp_path: Path) -> None:
    url = "https://ftp.ncbi.nlm.nih.gov/GSE9_DESeq2_results.csv"
    payload = b"gene,log2FoldChange,pvalue\nA,1.0,0.01\nB,-1.0,0.02\n"
    result = download_public_candidate(
        {"url": url},
        tmp_path / "candidate.csv",
        transport=public_transport({url: Response(url, payload)}),
    )
    assert result["status"] == "downloaded"
    assert result["sha256"] == hashlib.sha256(payload).hexdigest()
    assert (tmp_path / "candidate.csv").read_bytes() == payload


def test_download_public_candidate_normalizes_stale_http_link_without_url_leak(tmp_path: Path) -> None:
    url = "https://zenodo.org/files/missing.csv?token=do-not-leak"

    class MissingTransport:
        def validate_url(self, candidate_url: str) -> str:
            return candidate_url

        def get_bytes(self, candidate_url: str, *, max_bytes: int) -> bytes:
            raise urllib.error.HTTPError(candidate_url, 404, "Not Found", {}, None)

    with pytest.raises(DiscoveryUnavailableError, match="HTTP 404") as exc_info:
        download_public_candidate(
            {"url": url, "name": "missing.csv"},
            tmp_path / "missing.csv",
            transport=MissingTransport(),
        )

    assert "do-not-leak" not in str(exc_info.value)
    assert url not in str(exc_info.value)
    assert not (tmp_path / "missing.csv").exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_inspect_public_archive_finds_nested_tabular_members() -> None:
    nested = io.BytesIO()
    with zipfile.ZipFile(nested, "w") as archive:
        archive.writestr("inner_DESeq2_results.csv", "gene,log2FoldChange,pvalue\nA,1,0.01\n")
    outer = io.BytesIO()
    with zipfile.ZipFile(outer, "w") as archive:
        archive.writestr("nested.zip", nested.getvalue())
    records = inspect_public_archive(outer.getvalue())
    assert records[0]["member_name"] == "nested.zip!/inner_DESeq2_results.csv"
    assert records[0]["role"] == "deg_table"


def test_inspect_public_archive_rejects_zip_slip_symlink_and_bomb() -> None:
    slip = io.BytesIO()
    with zipfile.ZipFile(slip, "w") as archive:
        archive.writestr("../escape.csv", "x\n")
    with pytest.raises(DiscoveryError, match="unsafe member path"):
        inspect_public_archive(slip.getvalue())

    link = io.BytesIO()
    info = zipfile.ZipInfo("link.csv")
    info.external_attr = 0o120777 << 16
    with zipfile.ZipFile(link, "w") as archive:
        archive.writestr(info, "target")
    with pytest.raises(DiscoveryError, match="symbolic link"):
        inspect_public_archive(link.getvalue())

    bomb = io.BytesIO()
    with zipfile.ZipFile(bomb, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("huge.csv", b"0" * (51 * 1024 * 1024))
    with pytest.raises(DiscoveryError, match="oversized member"):
        inspect_public_archive(bomb.getvalue())


# --- naming what a failed download actually was ----------------------------


@pytest.mark.parametrize(
    ("payload", "expected"),
    [
        (b"<!DOCTYPE html>\n<html><body>404</body></html>", "web page"),
        (b"<html><head><title>Sign in</title></head>", "web page"),
        (b"\x1f\x8b\x08\x00\x00\x00\x00\x00", "gzip"),
        (b"%PDF-1.7\n%\xe2\xe3", "PDF"),
        (b'{"error": "not found"}', "JSON"),
        (b"", "empty response"),
        (b"gene\tlog2FC\tpvalue\nTP53\t2.1\t0.001\n", "plain text"),
        (b"\x00" * 257 + b"ustar\x0000", "tar archive"),
    ],
)
def test_a_failed_download_is_described_by_what_it_actually_is(payload: bytes, expected: str) -> None:
    """"Not a valid ZIP" sends the reader hunting for a corrupt file.

    Naming the payload tells them whether the link moved, needs a login, or was
    simply published in another format.
    """

    assert expected in describe_unexpected_payload(payload)


def test_an_invalid_archive_reports_what_arrived_instead(tmp_path: Path) -> None:
    with pytest.raises(DiscoveryError, match="web page"):
        inspect_public_archive(b"<!DOCTYPE html><html><body>Not found</body></html>")
