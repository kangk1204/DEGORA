"""Public-source discovery adapters with a conservative network boundary.

This module intentionally does not change the existing GEO discovery API.  It
offers standalone public-provider adapters and a transport that can be tested
with injected DNS/openers instead of live network calls.
"""

from __future__ import annotations

import hashlib
import http.client
import io
import ipaddress
import json
import mimetypes
import os
import re
import socket
import stat
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable

from .discovery import (
    API_KEY_PACE_SECONDS,
    DEFAULT_PACE_SECONDS,
    DiscoveryError,
    DiscoveryUnavailableError,
    DiscoveryUnsafeArchiveError,
    MAX_CANDIDATE_BYTES,
    MAX_JSON_BYTES,
    MAX_SOFT_BYTES,
    NcbiRequestConfig,
    SpeciesSpec,
    _query_terms,
    classify_filename,
    inspect_candidate_bytes,
    inspect_upstream_bytes,
    normalize_species,
    parse_geo_soft,
)


PUBLIC_SOURCE_USER_AGENT = "DEGORA-public-sources/0.1 (academic public metadata discovery)"
DEFAULT_TIMEOUT_SECONDS = 20
DEFAULT_TRANSPORT_RETRIES = 3
MAX_RETRY_DELAY_SECONDS = 6.0
MAX_REDIRECTS = 5
MAX_ARCHIVE_DEPTH = 2
MAX_ARCHIVE_MEMBERS = 2_000
MAX_ARCHIVE_EXPANDED_BYTES = 250 * 1024 * 1024
MAX_ARCHIVE_MEMBER_BYTES = 50 * 1024 * 1024
PUBMED_SUMMARY_BATCH_SIZE = 200
GEO_SUMMARY_BATCH_SIZE = 200
GEO_DETAILED_FILE_ASSESSMENT_LIMIT = 20
SUPPORTED_DOWNLOAD_EXTENSIONS = {
    ".csv",
    ".tsv",
    ".txt",
    ".xlsx",
    ".zip",
    ".gz",
}

PROVIDER_HOST_ALLOWLIST = frozenset(
    {
        "eutils.ncbi.nlm.nih.gov",
        "www.ncbi.nlm.nih.gov",
        "ftp.ncbi.nlm.nih.gov",
        "pmc.ncbi.nlm.nih.gov",
        "www.ebi.ac.uk",
        "ebi.ac.uk",
        "europepmc.org",
        "www.europepmc.org",
        "api.crossref.org",
        "doi.crossref.org",
        "api.datacite.org",
        "zenodo.org",
        "sandbox.zenodo.org",
        "figshare.com",
        "api.figshare.com",
        "dryad.org",
        "datadryad.org",
        "mendeley.com",
        "data.mendeley.com",
        "static-content.springer.com",
        "link.springer.com",
        "www.mdpi.com",
        "mdpi.com",
        "mdpi-res.com",
        "www.frontiersin.org",
        "frontiersin.org",
        "journals.plos.org",
        "onlinelibrary.wiley.com",
        "ars.els-cdn.com",
        "www.science.org",
    }
)

PROVIDER_HOST_SUFFIX_ALLOWLIST = (
    ".ncbi.nlm.nih.gov",
    ".ebi.ac.uk",
    ".europepmc.org",
    ".crossref.org",
    ".datacite.org",
    ".zenodo.org",
    ".figshare.com",
    ".dryad.org",
    ".datadryad.org",
    ".mendeley.com",
    ".springer.com",
    ".mdpi.com",
    ".mdpi-res.com",
    ".frontiersin.org",
    ".plos.org",
    ".wiley.com",
    ".els-cdn.com",
    ".science.org",
)

_NCBI_PACE_LOCK = threading.Lock()
_NCBI_LAST_REQUEST = 0.0


def _safe_remote_error(exc: Exception) -> str:
    """Return useful bounded diagnostics without retaining URL credentials."""

    if isinstance(exc, urllib.error.HTTPError):
        return f"HTTPError: remote service returned HTTP {exc.code}"
    text = re.sub(
        r"(?i)(api[_-]?key|access[_-]?token|token|password|secret)=([^&\s]+)",
        r"\1=[redacted]",
        str(exc),
    )
    text = re.sub(r"[\r\n\t]+", " ", text).strip()
    return f"{type(exc).__name__}: {(text or 'remote request failed')[:240]}"


class _NoRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        return None


def _has_allowed_host(host: str) -> bool:
    host = host.lower().rstrip(".")
    return host in PROVIDER_HOST_ALLOWLIST or any(host.endswith(suffix) for suffix in PROVIDER_HOST_SUFFIX_ALLOWLIST)


def _is_public_ip(address: str) -> bool:
    ip = ipaddress.ip_address(address)
    return not (
        ip.is_private
        or ip.is_loopback
        or ip.is_link_local
        or ip.is_reserved
        or ip.is_multicast
        or ip.is_unspecified
    )


def _is_transient_http_status(status: int) -> bool:
    return status in {408, 425, 429, 500, 502, 503, 504}


def _retry_delay_seconds(exc: Exception, attempt: int) -> float:
    if isinstance(exc, urllib.error.HTTPError) and exc.headers is not None:
        retry_after = exc.headers.get("Retry-After")
        try:
            parsed = float(retry_after) if retry_after is not None else None
        except (TypeError, ValueError):
            parsed = None
        if parsed is not None and parsed >= 0:
            return min(parsed, MAX_RETRY_DELAY_SECONDS)
    return min(0.5 * (2**attempt), MAX_RETRY_DELAY_SECONDS)


def _validate_public_url(url: str, *, resolver: Callable[..., Any]) -> str:
    parsed = urllib.parse.urlsplit(str(url).strip())
    host = (parsed.hostname or "").lower().rstrip(".")
    if parsed.scheme.lower() != "https":
        raise DiscoveryError("public-source URLs must use HTTPS")
    if not host or not _has_allowed_host(host):
        raise DiscoveryError(f"public-source URL host is not allowed: {host or '(missing)'}")
    if parsed.username or parsed.password:
        raise DiscoveryError("public-source URLs must not contain credentials")
    if parsed.fragment:
        raise DiscoveryError("public-source URLs must not contain fragments")
    if parsed.port not in (None, 443):
        raise DiscoveryError("public-source URLs must use the default HTTPS port")
    try:
        answers = resolver(host, 443, type=socket.SOCK_STREAM)
    except TypeError:
        answers = resolver(host, 443)
    except OSError as exc:
        raise DiscoveryUnavailableError(f"DNS resolution failed for {host}") from exc
    addresses = {str(item[4][0]) for item in answers}
    if not addresses:
        raise DiscoveryUnavailableError(f"DNS resolution returned no addresses for {host}")
    blocked = sorted(address for address in addresses if not _is_public_ip(address))
    if blocked:
        raise DiscoveryError(f"public-source DNS result is not public: {host}")
    return urllib.parse.urlunsplit(("https", host, parsed.path or "/", parsed.query, ""))


class SafePublicTransport:
    """HTTPS-only bounded transport for known public metadata/file providers."""

    def __init__(
        self,
        *,
        timeout: int = DEFAULT_TIMEOUT_SECONDS,
        max_response_bytes: int = MAX_JSON_BYTES,
        resolver: Callable[..., Any] = socket.getaddrinfo,
        opener: Any | None = None,
        sleep: Callable[[float], None] = time.sleep,
        retries: int = DEFAULT_TRANSPORT_RETRIES,
    ) -> None:
        if timeout <= 0 or timeout > 60:
            raise DiscoveryError("transport timeout must be between 1 and 60 seconds")
        if max_response_bytes <= 0:
            raise DiscoveryError("transport byte cap must be positive")
        if isinstance(retries, bool) or not isinstance(retries, int) or not 1 <= retries <= 5:
            raise DiscoveryError("transport retries must be a whole number from 1 to 5")
        self.timeout = int(timeout)
        self.max_response_bytes = int(max_response_bytes)
        self.resolver = resolver
        self.opener = opener or urllib.request.build_opener(_NoRedirectHandler())
        self.sleep = sleep
        self.retries = retries

    def validate_url(self, url: str) -> str:
        return _validate_public_url(url, resolver=self.resolver)

    def get_bytes(
        self,
        url: str,
        *,
        headers: dict[str, str] | None = None,
        timeout: int | None = None,
        max_bytes: int | None = None,
    ) -> bytes:
        safe_url = self.validate_url(url)
        request_headers = {"User-Agent": PUBLIC_SOURCE_USER_AGENT, **dict(headers or {})}
        limit = int(max_bytes if max_bytes is not None else self.max_response_bytes)
        last_error: Exception | None = None
        for attempt in range(self.retries):
            current_url = safe_url
            try:
                for _redirect in range(MAX_REDIRECTS + 1):
                    request = urllib.request.Request(current_url, headers=request_headers)
                    try:
                        with self.opener.open(request, timeout=timeout or self.timeout) as response:
                            status = getattr(response, "status", getattr(response, "code", 200))
                            if status in {301, 302, 303, 307, 308}:
                                location = response.headers.get("Location")
                                if not location:
                                    raise DiscoveryUnavailableError("redirect response omitted Location")
                                current_url = self.validate_url(urllib.parse.urljoin(current_url, location))
                                continue
                            if int(status) >= 400:
                                raise urllib.error.HTTPError(
                                    current_url,
                                    int(status),
                                    "remote service error",
                                    response.headers,
                                    None,
                                )
                            final_url = getattr(response, "geturl", lambda: current_url)()
                            self.validate_url(final_url)
                            declared = response.headers.get("Content-Length")
                            if declared and "Range" not in request_headers:
                                try:
                                    declared_bytes = int(declared)
                                except ValueError:
                                    declared_bytes = None
                                if declared_bytes is not None and declared_bytes > limit:
                                    raise DiscoveryError(f"remote response exceeds the {limit}-byte safety cap")
                            chunks: list[bytes] = []
                            total = 0
                            while True:
                                chunk = response.read(min(128 * 1024, limit + 1 - total))
                                if not chunk:
                                    break
                                total += len(chunk)
                                if total > limit:
                                    raise DiscoveryError(f"remote response exceeds the {limit}-byte safety cap")
                                chunks.append(chunk)
                            return b"".join(chunks)
                    except urllib.error.HTTPError as exc:
                        if exc.code in {301, 302, 303, 307, 308}:
                            location = exc.headers.get("Location") if exc.headers is not None else None
                            if not location:
                                raise DiscoveryUnavailableError("redirect response omitted Location") from exc
                            current_url = self.validate_url(urllib.parse.urljoin(current_url, location))
                            continue
                        raise
                raise DiscoveryError("public-source redirect chain exceeded the safety limit")
            except urllib.error.HTTPError as exc:
                if not _is_transient_http_status(exc.code):
                    raise
                last_error = exc
                if attempt < self.retries - 1:
                    self.sleep(_retry_delay_seconds(exc, attempt))
            except (urllib.error.URLError, TimeoutError, http.client.HTTPException, OSError) as exc:
                last_error = exc
                if attempt < self.retries - 1:
                    self.sleep(_retry_delay_seconds(exc, attempt))
        assert last_error is not None
        detail = f"HTTP {last_error.code}" if isinstance(last_error, urllib.error.HTTPError) else type(last_error).__name__
        raise DiscoveryUnavailableError(
            f"public source request failed after {self.retries} attempt(s) ({detail})"
        ) from last_error

    def get_json(
        self,
        url: str,
        *,
        headers: dict[str, str] | None = None,
        timeout: int | None = None,
        max_bytes: int | None = None,
    ) -> dict[str, Any]:
        payload = self.get_bytes(url, headers=headers, timeout=timeout, max_bytes=max_bytes or MAX_JSON_BYTES)
        try:
            parsed = json.loads(payload.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise DiscoveryUnavailableError("public source returned invalid JSON") from exc
        if not isinstance(parsed, dict):
            raise DiscoveryUnavailableError("public source returned unexpected JSON")
        return parsed


@dataclass
class ProviderResult:
    records: list[dict[str, Any]] = field(default_factory=list)
    diagnostics: list[dict[str, Any]] = field(default_factory=list)


def _topic_query(query: str) -> str:
    return " AND ".join(f'"{term}"[All Fields]' for term in _query_terms(query))


def _species_query(query: str, species: SpeciesSpec) -> str:
    topic = _topic_query(query)
    if species.key == "human":
        organism = (
            '"Humans"[MeSH Terms] OR "Homo sapiens"[All Fields] OR '
            '"human"[Title/Abstract] OR "humans"[Title/Abstract]'
        )
    else:
        organism = (
            '"Mice"[MeSH Terms] OR "Mus musculus"[All Fields] OR '
            '"mouse"[Title/Abstract] OR "mice"[Title/Abstract] OR "murine"[Title/Abstract]'
        )
    return f"({topic}) AND ({organism})"


def _pace_ncbi(interval: float, sleep: Callable[[float], None]) -> None:
    """Serialize request starts across concurrent Human/Mouse jobs."""

    global _NCBI_LAST_REQUEST
    with _NCBI_PACE_LOCK:
        now = time.monotonic()
        delay = max(0.0, _NCBI_LAST_REQUEST + interval - now)
        if delay:
            sleep(delay)
            now = time.monotonic()
        _NCBI_LAST_REQUEST = now


def _clean_list(values: Iterable[Any]) -> list[str]:
    cleaned: list[str] = []
    for value in values:
        if value is None:
            continue
        text = str(value).strip()
        if not text or text.lower() in {"none", "null", "nan", "n/a"}:
            continue
        cleaned.append(text)
    return list(dict.fromkeys(cleaned))


def _extract_article_ids(summary: dict[str, Any]) -> dict[str, str]:
    ids: dict[str, str] = {}
    for article_id in summary.get("articleids", []) or []:
        if not isinstance(article_id, dict):
            continue
        kind = str(article_id.get("idtype") or "").lower()
        value = str(article_id.get("value") or "").strip()
        if kind in {"doi", "pmcid", "pubmed"} and value:
            ids[kind] = value
    elocation = str(summary.get("elocationid") or "")
    match = re.search(r"10\.\d{4,9}/\S+", elocation)
    if match and "doi" not in ids:
        ids["doi"] = match.group(0).rstrip(".")
    return ids


def _eutils_history_params(search_result: dict[str, Any], *, start: int, size: int) -> dict[str, Any] | None:
    query_key = str(search_result.get("querykey") or search_result.get("query_key") or "").strip()
    web_env = str(search_result.get("webenv") or search_result.get("WebEnv") or "").strip()
    if not query_key or not web_env:
        return None
    return {
        "query_key": query_key,
        "WebEnv": web_env,
        "retstart": start,
        "retmax": size,
    }


def _summary_supplement_candidates(summary: dict[str, Any], accession: str) -> list[dict[str, Any]]:
    """Extract lightweight file hints from the database-specific GDS summary."""

    values: list[str] = []

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            for nested in value.values():
                walk(nested)
        elif isinstance(value, (list, tuple)):
            for nested in value:
                walk(nested)
        elif value not in (None, ""):
            values.extend(part.strip() for part in re.split(r"[;,\n]", str(value)) if part.strip())

    walk(summary.get("suppfile"))
    candidates: dict[str, dict[str, Any]] = {}
    for value in values:
        urls = re.findall(r"https://[^\s<>\"]+", value)
        signals = urls or [value]
        for signal in signals:
            cleaned = signal.rstrip(".)]")
            assessment = classify_filename(cleaned)
            if not assessment.get("inspectable"):
                continue
            url = cleaned if cleaned.lower().startswith("https://") else ""
            key = url or assessment["name"]
            candidates.setdefault(
                key,
                {
                    "provider": "ncbi_geo_summary",
                    "record_type": "file_candidate",
                    "accession": accession,
                    "name": assessment["name"],
                    "url": url,
                    "role": assessment["role"],
                    "tier": assessment["tier"],
                    "downloaded": False,
                    "assessment_scope": "filename_only",
                },
            )
    return list(candidates.values())


class NcbiPubmedProvider:
    name = "ncbi_pubmed"

    def __init__(
        self,
        *,
        transport: SafePublicTransport | None = None,
        request_config: NcbiRequestConfig | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        self.transport = transport or SafePublicTransport()
        self.request_config = request_config or NcbiRequestConfig.from_environment()
        self.sleep = sleep
        self.pace_seconds = API_KEY_PACE_SECONDS if self.request_config.api_key else DEFAULT_PACE_SECONDS

    def _eutils_url(self, endpoint: str, params: dict[str, Any]) -> str:
        merged = {**params, **self.request_config.eutils_params()}
        return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint + "?" + urllib.parse.urlencode(merged)

    def _get_json(self, url: str) -> dict[str, Any]:
        _pace_ncbi(self.pace_seconds, self.sleep)
        payload = self.transport.get_json(url)
        error = payload.get("error") or payload.get("esearchresult", {}).get("errorlist")
        if error:
            raise DiscoveryUnavailableError("NCBI returned an error response")
        return payload

    def _summary_record(
        self,
        pmid: str,
        item: dict[str, Any],
        *,
        rank: int | None = None,
        species: SpeciesSpec | None = None,
    ) -> dict[str, Any]:
        article_ids = _extract_article_ids(item)
        record: dict[str, Any] = {
            "provider": self.name,
            "record_type": "publication",
            "pmid": pmid,
            "doi": article_ids.get("doi", ""),
            "pmcid": article_ids.get("pmcid", ""),
            "title": str(item.get("title") or "").strip(),
            "journal": str(item.get("fulljournalname") or item.get("source") or "").strip(),
            "publication_date": str(item.get("sortpubdate") or item.get("pubdate") or "").strip(),
            "authors": [
                str(author.get("name") or "").strip()
                for author in item.get("authors", []) or []
                if isinstance(author, dict) and str(author.get("name") or "").strip()
            ],
            "identifiers": {"pmid": pmid, **article_ids},
        }
        if rank is not None:
            record["rank"] = rank
        if species is not None:
            record["species_evidence"] = {
                "requested": species.scientific_name,
                "basis": "PubMed organism-constrained query",
            }
        return record

    def search(self, query: str, species: SpeciesSpec, limit: int) -> list[dict[str, Any]]:
        spec = species if isinstance(species, SpeciesSpec) else normalize_species(str(species))
        capped = max(1, min(int(limit), 1000))
        search_url = self._eutils_url(
            "esearch.fcgi",
            {
                "db": "pubmed",
                "term": _species_query(query, spec),
                "retmax": capped,
                "retmode": "json",
                "sort": "relevance",
                "usehistory": "y",
            },
        )
        search = self._get_json(search_url)
        search_result = search.get("esearchresult", {})
        ids = [str(value) for value in search_result.get("idlist", [])]
        records: list[dict[str, Any]] = []
        for start in range(0, len(ids), PUBMED_SUMMARY_BATCH_SIZE):
            batch = ids[start : start + PUBMED_SUMMARY_BATCH_SIZE]
            history = _eutils_history_params(search_result, start=start, size=len(batch))
            summary_url = self._eutils_url(
                "esummary.fcgi",
                {"db": "pubmed", **(history or {"id": ",".join(batch)}), "retmode": "json"},
            )
            summary = self._get_json(summary_url).get("result", {})
            for rank, pmid in enumerate(batch, start=start + 1):
                item = summary.get(pmid)
                if not isinstance(item, dict):
                    continue
                records.append(self._summary_record(pmid, item, rank=rank, species=spec))
        return records

    def resolve(self, record: dict[str, Any], _species: SpeciesSpec) -> list[dict[str, Any]]:
        """Hydrate linked PMID metadata when a GEO-first result lacks citation fields."""

        raw_pmids: list[Any] = [record.get("pmid")]
        for value in (record.get("pubmed_ids"), record.get("pmids")):
            if isinstance(value, (list, tuple, set)):
                raw_pmids.extend(value)
            elif value:
                raw_pmids.append(value)
        identifiers = record.get("identifiers")
        if isinstance(identifiers, dict) and identifiers.get("pmid"):
            raw_pmids.append(identifiers["pmid"])
        pmids = _clean_list(raw_pmids)
        pmid = next((value for value in pmids if re.fullmatch(r"\d+", value)), "")
        if not pmid:
            return []
        if (
            str(record.get("paper_title") or record.get("title") or "").strip()
            and (record.get("authors") or record.get("authors_display"))
            and str(record.get("journal") or "").strip()
            and (record.get("year") or record.get("publication_date"))
        ):
            return []
        summary_url = self._eutils_url(
            "esummary.fcgi",
            {"db": "pubmed", "id": pmid, "retmode": "json"},
        )
        item = self._get_json(summary_url).get("result", {}).get(pmid)
        if not isinstance(item, dict):
            return []
        return [self._summary_record(pmid, item)]


class NcbiGeoProvider:
    name = "ncbi_geo"

    def __init__(
        self,
        *,
        transport: SafePublicTransport | None = None,
        request_config: NcbiRequestConfig | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        self.transport = transport or SafePublicTransport(max_response_bytes=MAX_SOFT_BYTES)
        self.request_config = request_config or NcbiRequestConfig.from_environment()
        self.sleep = sleep
        self.pace_seconds = API_KEY_PACE_SECONDS if self.request_config.api_key else DEFAULT_PACE_SECONDS

    def _eutils_url(self, endpoint: str, params: dict[str, Any]) -> str:
        merged = {**params, **self.request_config.eutils_params()}
        return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint + "?" + urllib.parse.urlencode(merged)

    def _get_json(self, url: str) -> dict[str, Any]:
        _pace_ncbi(self.pace_seconds, self.sleep)
        payload = self.transport.get_json(url)
        return payload

    def search(self, query: str, species: SpeciesSpec, limit: int) -> list[dict[str, Any]]:
        spec = species if isinstance(species, SpeciesSpec) else normalize_species(str(species))
        term = f'({_topic_query(query)}) AND gse[ETYP] AND "{spec.scientific_name}"[Organism]'
        capped = max(1, min(int(limit), 1000))
        search_url = self._eutils_url(
            "esearch.fcgi",
            {
                "db": "gds",
                "term": term,
                "retmax": capped,
                "retmode": "json",
                "sort": "relevance",
                "usehistory": "y",
            },
        )
        search = self._get_json(search_url)
        search_result = search.get("esearchresult", {})
        ids = [str(value) for value in search_result.get("idlist", [])]
        if not ids:
            return []
        summary: dict[str, Any] = {}
        for start in range(0, len(ids), GEO_SUMMARY_BATCH_SIZE):
            batch = ids[start : start + GEO_SUMMARY_BATCH_SIZE]
            history = _eutils_history_params(search_result, start=start, size=len(batch))
            summary_url = self._eutils_url(
                "esummary.fcgi",
                {"db": "gds", **(history or {"id": ",".join(batch)}), "retmode": "json"},
            )
            summary.update(self._get_json(summary_url).get("result", {}))
        records: list[dict[str, Any]] = []
        for rank, uid in enumerate(ids, start=1):
            item = summary.get(uid)
            if not isinstance(item, dict):
                continue
            accession = str(item.get("accession") or "").upper()
            if not re.fullmatch(r"GSE\d+", accession):
                continue
            parsed: dict[str, Any] = {}
            detail_error = ""
            if rank <= GEO_DETAILED_FILE_ASSESSMENT_LIMIT:
                try:
                    _pace_ncbi(self.pace_seconds, self.sleep)
                    soft = self.transport.get_bytes(
                        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
                        + urllib.parse.urlencode({"acc": accession, "targ": "self", "form": "text", "view": "brief"}),
                        max_bytes=MAX_SOFT_BYTES,
                    ).decode("utf-8", "replace")
                    parsed = parse_geo_soft(soft)
                except Exception as exc:  # retain the metadata result when one detail route is unavailable
                    detail_error = _safe_remote_error(exc)
            taxa = parsed.get("taxa", [])
            is_exact = bool(taxa) and taxa == [spec.scientific_name]
            is_mixed = bool(taxa) and not is_exact
            supplement_candidates = _summary_supplement_candidates(item, accession)
            for url in parsed.get("supplementary_files", []):
                assessment = classify_filename(url)
                if assessment.get("inspectable"):
                    soft_candidate = {
                        "provider": self.name,
                        "record_type": "file_candidate",
                        "accession": accession,
                        "name": assessment["name"],
                        "url": url,
                        "role": assessment["role"],
                        "tier": assessment["tier"],
                        "downloaded": False,
                        "assessment_scope": "geo_soft_filename",
                    }
                    supplement_candidates = [
                        candidate
                        for candidate in supplement_candidates
                        if candidate.get("url") != url and candidate.get("name") != assessment["name"]
                    ]
                    supplement_candidates.append(soft_candidate)
            records.append(
                {
                    "provider": self.name,
                    "record_type": "publication",
                    "rank": rank,
                    "accession": accession,
                    "pmid": _clean_list(item.get("pubmedids", []))[:1][0] if item.get("pubmedids") else "",
                    "pmids": _clean_list([*item.get("pubmedids", []), *parsed.get("pubmed_ids", [])]),
                    "title": parsed.get("title") or str(item.get("title") or "").strip(),
                    "summary": parsed.get("summary") or str(item.get("summary") or "").strip(),
                    "publication_date": str(item.get("pdat") or "").strip(),
                    "species_evidence": {
                        "requested": spec.scientific_name,
                        "observed_taxa": taxa,
                        "status": "exact" if is_exact else ("quarantined_mixed_record" if is_mixed else "query_constrained_unverified"),
                        "basis": "GEO SOFT observed taxa" if taxa else "GEO organism-constrained query",
                    },
                    "quarantined": is_mixed,
                    "quarantine_reason": "GEO SOFT record contains mixed or mismatched taxa" if is_mixed else "",
                    "supplementary_file_candidates": supplement_candidates,
                    "detail_assessment": "complete" if parsed else "not_evaluated",
                    "detail_error": detail_error,
                }
            )
        return records


def _candidate_url(candidate: dict[str, Any]) -> str:
    url = str(candidate.get("url") or candidate.get("source_url") or "").strip()
    if not url:
        raise DiscoveryError("candidate is missing a public direct file URL")
    return url


def _is_plausible_direct_file(url: str, content_type: str = "", filename: str = "") -> bool:
    path = urllib.parse.urlsplit(url).path.lower()
    candidate_path = str(filename or path).lower()
    suffixes = "".join(Path(candidate_path).suffixes[-2:])
    ext = Path(candidate_path).suffix.lower()
    if ext in SUPPORTED_DOWNLOAD_EXTENSIONS or suffixes.endswith(".csv.gz") or suffixes.endswith(".tsv.gz") or suffixes.endswith(".txt.gz"):
        return True
    lowered_type = content_type.lower()
    return any(token in lowered_type for token in ("text/csv", "tab-separated-values", "spreadsheet", "zip", "gzip"))


def download_public_candidate(
    candidate: dict[str, Any],
    target: str | Path,
    *,
    max_bytes: int = MAX_CANDIDATE_BYTES,
    transport: SafePublicTransport | None = None,
) -> dict[str, Any]:
    """Download one public direct-file candidate with atomic promotion."""

    if max_bytes <= 0 or max_bytes > MAX_CANDIDATE_BYTES:
        raise DiscoveryError(f"max_bytes must be between 1 and {MAX_CANDIDATE_BYTES}")
    url = _candidate_url(candidate)
    active_transport = transport or SafePublicTransport(max_response_bytes=max_bytes)
    safe_url = active_transport.validate_url(url)
    try:
        payload = active_transport.get_bytes(safe_url, max_bytes=max_bytes)
    except urllib.error.HTTPError as exc:
        # Repository metadata can outlive a moved or withdrawn supplementary
        # file. Treat that candidate as unavailable so preparation can inspect
        # the remaining routes without exposing the signed/query URL.
        raise DiscoveryUnavailableError(
            f"public candidate returned HTTP {exc.code}"
        ) from exc
    guessed_type, _ = mimetypes.guess_type(urllib.parse.urlsplit(safe_url).path)
    if not _is_plausible_direct_file(safe_url, guessed_type or "", str(candidate.get("name") or "")):
        raise DiscoveryError("candidate URL does not look like a public direct tabular/archive file")
    target_path = Path(target)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(payload).hexdigest()
    fd, tmp_name = tempfile.mkstemp(prefix=f".{target_path.name}.", suffix=".tmp", dir=target_path.parent)
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, target_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return {
        "status": "downloaded",
        "url": safe_url,
        "target": str(target_path),
        "bytes": len(payload),
        "sha256": digest,
        "mime_type": guessed_type or "",
        "extension": "".join(target_path.suffixes[-2:]) or target_path.suffix,
    }


# A publisher that has retired a supplementary file usually still answers the
# request - with an HTML error page, a login redirect, or the wrong format
# entirely.  "not a valid ZIP file" sends the reader looking for a corrupt
# download; naming what actually arrived sends them to the right place.
_PAYLOAD_SIGNATURES: tuple[tuple[bytes, str], ...] = (
    (b"\x1f\x8b", "a gzip stream (.gz), not a ZIP"),
    (b"BZh", "a bzip2 stream, not a ZIP"),
    (b"\xfd7zXZ\x00", "an xz stream, not a ZIP"),
    (b"7z\xbc\xaf\x27\x1c", "a 7-Zip archive, not a ZIP"),
    (b"Rar!", "a RAR archive, not a ZIP"),
    (b"%PDF", "a PDF, not a data archive"),
    (b"\x89PNG", "a PNG image, not a data archive"),
    (b"\xff\xd8\xff", "a JPEG image, not a data archive"),
    (b"\xd0\xcf\x11\xe0", "a legacy Office document (.xls/.doc), not a ZIP"),
)
_HTML_MARKERS: tuple[bytes, ...] = (b"<!doctype html", b"<html", b"<head", b"<?xml")


def describe_unexpected_payload(payload: bytes) -> str:
    """Name what a download turned out to be when it is not the ZIP we asked for."""

    if not payload:
        return "an empty response"
    for signature, label in _PAYLOAD_SIGNATURES:
        if payload.startswith(signature):
            return label
    head = payload[:1024].lstrip().lower()
    for marker in _HTML_MARKERS:
        if head.startswith(marker):
            return "a web page, so the file has most likely moved or now needs a login"
    if payload[257:262] == b"ustar":
        return "a tar archive, not a ZIP"
    if head[:1] in (b"{", b"[") and b'"' in head:
        return "a JSON document, not a data archive"
    try:
        text = payload[:256].decode("utf-8")
    except UnicodeDecodeError:
        return "an unrecognized binary format"
    if text.strip():
        return "plain text, not a ZIP"
    return "an unrecognized binary format"


def _safe_archive_member(info: zipfile.ZipInfo) -> None:
    name = info.filename
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts:
        raise DiscoveryUnsafeArchiveError("archive contains an unsafe member path")
    mode = info.external_attr >> 16
    if stat.S_ISLNK(mode):
        raise DiscoveryUnsafeArchiveError("archive contains a symbolic link")
    if info.file_size > MAX_ARCHIVE_MEMBER_BYTES:
        raise DiscoveryUnsafeArchiveError("archive contains an oversized member")


def inspect_public_archive(
    payload: bytes,
    *,
    source_name: str = "candidate.zip",
    max_depth: int = MAX_ARCHIVE_DEPTH,
) -> list[dict[str, Any]]:
    """Inspect ZIP and nested ZIP payloads without extracting to the filesystem."""

    records: list[dict[str, Any]] = []
    total_members = 0
    total_expanded = 0

    def visit(data: bytes, prefix: str, depth: int) -> None:
        nonlocal total_members, total_expanded
        if depth > max_depth:
            raise DiscoveryUnsafeArchiveError("nested ZIP depth exceeds the safety limit")
        try:
            with zipfile.ZipFile(io.BytesIO(data)) as archive:
                for info in archive.infolist():
                    if info.is_dir():
                        continue
                    _safe_archive_member(info)
                    total_members += 1
                    if total_members > MAX_ARCHIVE_MEMBERS:
                        raise DiscoveryUnsafeArchiveError("archive contains too many members")
                    total_expanded += info.file_size
                    if total_expanded > MAX_ARCHIVE_EXPANDED_BYTES:
                        raise DiscoveryUnsafeArchiveError("archive expanded-size cap exceeded")
                    member_name = f"{prefix}{info.filename}"
                    lower = info.filename.lower()
                    if lower.endswith(".zip"):
                        visit(archive.read(info), f"{member_name}!/", depth + 1)
                        continue
                    assessment = classify_filename(info.filename)
                    if assessment.get("inspectable"):
                        records.append(
                            {
                                "provider": "zip_inspector",
                                "record_type": "file_candidate",
                                "archive": source_name,
                                "member_name": member_name,
                                "name": assessment["name"],
                                "role": assessment["role"],
                                "tier": assessment["tier"],
                                "bytes": info.file_size,
                            }
                        )
        except zipfile.BadZipFile as exc:
            raise DiscoveryError(
                f"candidate archive is not a valid ZIP file: the download is {describe_unexpected_payload(data)}"
            ) from exc

    visit(payload, "", 0)
    return records


class PublicRepositoryResolver:
    name = "public_repository_resolver"

    def __init__(self, *, transport: SafePublicTransport | None = None) -> None:
        self.transport = transport or SafePublicTransport()

    def resolve(self, record: dict[str, Any], species: SpeciesSpec) -> list[dict[str, Any]]:
        candidates: list[dict[str, Any]] = []
        diagnostics: list[dict[str, Any]] = []
        primary = (
            ("europe_pmc", self._resolve_europe_pmc),
            ("biostudies", self._resolve_biostudies),
        )
        fallback = (
            ("crossref", self._resolve_crossref),
            ("datacite", self._resolve_datacite),
        )
        for routes in (primary, fallback):
            if routes is fallback and candidates:
                break
            for provider_name, builder in routes:
                try:
                    candidates.extend(builder(record, species))
                except Exception as exc:  # provider failures must not fail the whole resolution
                    diagnostics.append(
                        {
                            "provider": provider_name,
                            "status": "candidate_error",
                            "error": _safe_remote_error(exc),
                        }
                    )
        return candidates + ([{"provider": self.name, "record_type": "diagnostics", "diagnostics": diagnostics}] if diagnostics else [])

    def _resolve_europe_pmc(self, record: dict[str, Any], _species: SpeciesSpec) -> list[dict[str, Any]]:
        query_parts = []
        if record.get("pmid"):
            query_parts.append(f'EXT_ID:{record["pmid"]} SRC:MED')
        if record.get("pmcid"):
            query_parts.append(f'PMCID:{record["pmcid"]}')
        if record.get("doi"):
            query_parts.append(f'DOI:"{record["doi"]}"')
        if not query_parts:
            return []
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
            {"query": " OR ".join(query_parts), "format": "json", "resultType": "core", "pageSize": 10}
        )
        data = self.transport.get_json(url)
        candidates = self._records_from_links("europe_pmc", data)
        results = data.get("resultList", {}).get("result", [])
        if isinstance(results, list):
            for result in results:
                if not isinstance(result, dict):
                    continue
                pmcid = str(result.get("pmcid") or record.get("pmcid") or "").strip().upper()
                has_supplement = str(result.get("hasSuppl") or "").strip().lower() in {"1", "true", "y", "yes"}
                if not pmcid or not has_supplement:
                    continue
                candidates.append(
                    {
                        "provider": "europe_pmc",
                        "record_type": "file_candidate",
                        "url": (
                            "https://www.ebi.ac.uk/europepmc/webservices/rest/"
                            f"{urllib.parse.quote(pmcid, safe='')}/supplementaryFiles?includeInlineImage=n"
                        ),
                        "name": f"{pmcid}_supplementary_files.zip",
                        "role": "archive",
                        "tier": "archive_candidate",
                        "downloaded": False,
                    }
                )
        deduped = {str(candidate.get("url") or candidate.get("name")): candidate for candidate in candidates}
        return list(deduped.values())

    def _resolve_biostudies(self, record: dict[str, Any], _species: SpeciesSpec) -> list[dict[str, Any]]:
        # BioStudies is a file-level confirmation layer, not a reliable DOI or
        # PMID index.  Use one exact stable identifier and hydrate only matching
        # hits instead of issuing a broad title/species search.
        identifiers = _clean_list(
            [record.get("accession"), record.get("pmcid"), record.get("doi"), record.get("pmid")]
        )
        if not identifiers:
            return []
        query = identifiers[0]
        search_url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(
            {"query": query, "pageSize": 3}
        )
        search = self.transport.get_json(search_url)
        hits = search.get("hits", [])
        if not isinstance(hits, list):
            return []
        candidates: list[dict[str, Any]] = []
        for hit in hits[:3]:
            if not isinstance(hit, dict) or hit.get("isPublic") is False:
                continue
            accession = str(hit.get("accession") or "").strip()
            if not accession:
                continue
            encoded_accession = urllib.parse.quote(accession, safe="")
            study_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{encoded_accession}"
            info_url = f"{study_url}/info"
            study = self.transport.get_json(study_url)
            info = self.transport.get_json(info_url)
            candidates.extend(self._records_from_links("biostudies", study))
            base_url = str(info.get("httpLink") or "").rstrip("/")
            if not base_url:
                continue
            for path in _biostudies_file_paths(study)[:100]:
                assessment = classify_filename(path)
                if not assessment.get("inspectable"):
                    continue
                relative_path = PurePosixPath(path)
                path_parts = relative_path.parts
                if path_parts and path_parts[0].lower() == "files":
                    relative_path = PurePosixPath(*path_parts[1:])
                candidates.append(
                    {
                        "provider": "biostudies",
                        "record_type": "file_candidate",
                        "provider_accession": accession,
                        # BioStudies exposes study metadata at the accession
                        # root, while every deposited file lives below its
                        # case-sensitive Files/ directory.
                        "url": base_url + "/Files/" + urllib.parse.quote(str(relative_path), safe="/"),
                        "name": assessment["name"],
                        "role": assessment["role"],
                        "tier": assessment["tier"],
                        "downloaded": False,
                    }
                )
        deduped = {str(candidate.get("url") or candidate.get("name")): candidate for candidate in candidates}
        return list(deduped.values())

    def _resolve_crossref(self, record: dict[str, Any], _species: SpeciesSpec) -> list[dict[str, Any]]:
        doi = str(record.get("doi") or "").strip()
        if not doi:
            return []
        url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
        try:
            data = self.transport.get_json(url)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return []
            raise
        return self._records_from_links("crossref", data)

    def _resolve_datacite(self, record: dict[str, Any], _species: SpeciesSpec) -> list[dict[str, Any]]:
        doi = str(record.get("doi") or "").strip()
        if not doi:
            return []
        url = "https://api.datacite.org/dois/" + urllib.parse.quote(doi, safe="")
        try:
            data = self.transport.get_json(url)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return []
            raise
        return self._records_from_links("datacite", data)

    def _records_from_links(self, provider: str, data: Any) -> list[dict[str, Any]]:
        found: list[dict[str, Any]] = []

        def walk(value: Any) -> None:
            if isinstance(value, dict):
                url = value.get("url") or value.get("URL") or value.get("href")
                filename = value.get("name") or value.get("filename") or value.get("fileName") or ""
                if isinstance(url, str) and _link_looks_public_repository(url, str(filename)):
                    assessment = classify_filename(str(filename) or url)
                    found.append(
                        {
                            "provider": provider,
                            "record_type": "file_candidate",
                            "url": url,
                            "name": assessment["name"],
                            "role": assessment["role"],
                            "tier": assessment["tier"],
                            "downloaded": False,
                        }
                    )
                for item in value.values():
                    walk(item)
            elif isinstance(value, list):
                for item in value:
                    walk(item)

        walk(data)
        deduped: dict[str, dict[str, Any]] = {}
        for item in found:
            deduped.setdefault(item["url"], item)
        return list(deduped.values())


def _link_looks_public_repository(url: str, filename: str = "") -> bool:
    parsed = urllib.parse.urlsplit(str(url))
    host = (parsed.hostname or "").lower()
    path = parsed.path.lower()
    if not _has_allowed_host(host):
        return False
    candidate_name = str(filename or "").lower()
    if any(
        path.endswith(ext) or candidate_name.endswith(ext)
        for ext in (".csv", ".tsv", ".txt", ".xlsx", ".zip", ".gz")
    ):
        return True
    # Repository APIs sometimes expose an opaque numeric download URL next to a
    # separate filename.  The filename gate above still requires a supported
    # tabular/archive suffix before these routes are accepted.
    return bool(candidate_name) and any(
        marker in path
        for marker in ("/ndownloader/files/", "/files/", "/supplementaryfiles")
    )


def _biostudies_file_paths(study: dict[str, Any]) -> list[str]:
    paths: list[str] = []

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            raw_path = value.get("path")
            if raw_path:
                path = PurePosixPath(str(raw_path))
                if not path.is_absolute() and ".." not in path.parts:
                    paths.append(str(path))
            for nested in value.values():
                walk(nested)
        elif isinstance(value, list):
            for nested in value:
                walk(nested)

    walk(study.get("section", {}))
    return list(dict.fromkeys(paths))


def default_publication_providers() -> list[Any]:
    return [NcbiPubmedProvider(), NcbiGeoProvider(), PublicRepositoryResolver()]


def inspect_downloaded_candidate(path: str | Path) -> list[dict[str, Any]]:
    """Return candidate records for a downloaded file or inspect ZIP members."""

    file_path = Path(path)
    payload = file_path.read_bytes()
    if file_path.suffix.lower() == ".zip":
        return inspect_public_archive(payload, source_name=file_path.name)
    assessment = classify_filename(file_path.name)
    if not assessment.get("inspectable"):
        return []
    record: dict[str, Any] = {
        "provider": "download_inspector",
        "record_type": "file_candidate",
        "name": assessment["name"],
        "role": assessment["role"],
        "tier": assessment["tier"],
        "path": str(file_path),
    }
    try:
        if assessment["role"] == "deg_table":
            record["inspection"] = inspect_candidate_bytes(file_path.name, payload)
        elif assessment["role"] in {"count_matrix", "normalized_expression_matrix", "unknown_table"}:
            record["inspection"] = inspect_upstream_bytes(file_path.name, payload, declared_role=assessment["role"])
    except DiscoveryError as exc:
        record["inspection_error"] = str(exc)
    return [record]
