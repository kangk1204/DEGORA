"""Source-neutral publication discovery federation.

This module keeps provider-specific search and resolution separate from the
paper-level records used by DEGORA's discovery UI.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import re
from typing import Any, Callable, Iterable

from .discovery import (
    DiscoveryError,
    DiscoveryUnavailableError,
    SpeciesSpec,
    _query_terms,
    normalize_species as _normalize_species,
)


_MISSING = {"", "none", "null", "nan", "na", "n/a"}
_READINESS_PRIORITY = {
    "verified_ready": 0,
    "target_species_verified": 0,
    "mixed_rescued": 0,
    "likely_ready": 1,
    "likely": 1,
    "candidate": 2,
    "mixed_quarantined": 3,
    "metadata_only": 4,
    "unknown": 5,
}
MAX_SEARCH_RECORDS = 1_000
MAX_DETAILED_RESOLUTION_RECORDS = 20


def normalize_species(species: str | SpeciesSpec) -> SpeciesSpec:
    """Normalize the public discovery species contract to Human or Mouse."""

    if isinstance(species, SpeciesSpec):
        if species.key not in {"human", "mouse"}:
            raise ValueError("species must be Human or Mouse")
        return species
    return _normalize_species(str(species))


def canonical_record_id(record: dict[str, Any]) -> str:
    """Return the preferred stable identifier: PMID, normalized DOI, accession."""

    pmids = _as_list(record.get("pmid")) + _as_list(record.get("pubmed_ids"))
    pmid = _first_sorted(_normalize_pmid(value) for value in pmids)
    if pmid:
        return f"pmid:{pmid}"
    doi = _first_sorted(_normalize_doi(value) for value in _as_list(record.get("doi")) + _as_list(record.get("dois")))
    if doi:
        return f"doi:{doi}"
    accessions = _collect_accessions(record)
    if accessions:
        return f"accession:{accessions[0]}"
    provider_ids = _collect_provider_ids(record)
    if provider_ids:
        provider, value = provider_ids[0].split(":", 1)
        return f"provider:{provider}:{value}"
    title = _clean_text(record.get("paper_title") or record.get("title"))
    if title:
        return f"title:{title.lower()}"
    return "record:unknown"


def merge_publication_records(records: Iterable[dict[str, Any]], species: str | SpeciesSpec) -> list[dict[str, Any]]:
    """Merge publication records connected by publication, accession, or provider IDs."""

    target = normalize_species(species)
    prepared = [_prepare_record(record, target) for record in records]
    if not prepared:
        return []

    parent = list(range(len(prepared)))
    identifiers_by_value: dict[str, int] = {}
    for index, record in enumerate(prepared):
        for identifier in _graph_identifiers(record):
            previous = identifiers_by_value.get(identifier)
            if previous is None:
                identifiers_by_value[identifier] = index
            else:
                _union(parent, previous, index)

    groups: dict[int, list[dict[str, Any]]] = {}
    for index, record in enumerate(prepared):
        groups.setdefault(_find(parent, index), []).append(record)

    merged = [_merge_group(group, target) for group in groups.values()]
    return rank_publication_records(merged, sort_by="canonical_id", sort_order="asc")


def rank_publication_records(
    records: Iterable[dict[str, Any]],
    sort_by: str = "data_readiness",
    sort_order: str | None = None,
) -> list[dict[str, Any]]:
    """Return a deterministic globally ranked publication list."""

    rows = [_copy_record(record) for record in records]
    descending = (sort_order or _default_sort_order(sort_by)).lower() == "desc"

    if sort_by in {None, "", "data_readiness", "readiness", "readiness_score"}:
        # Lower readiness priority is better.  The UI calls the best-first view
        # "descending" even though the internal priority is an ascending rank.
        rows.sort(key=_default_rank_key, reverse=not descending)
    elif sort_by in {"relevance", "relevance_rank", "relevance_score", "pubmed_relevance"}:
        rows.sort(key=lambda row: (_number_or_inf(row.get("relevance_rank")), -_year(row), _canonical_for_sort(row)))
        if not descending:
            rows.reverse()
    elif sort_by in {"year"}:
        rows.sort(key=lambda row: (_year(row), _canonical_for_sort(row)), reverse=descending)
    elif sort_by in {"canonical_id", "id"}:
        rows.sort(key=_canonical_for_sort, reverse=descending)
    elif sort_by in {"paper_title", "title", "journal", "authors", "authors_display", "data_sources"}:
        key = "paper_title" if sort_by == "title" else sort_by
        if key == "authors":
            key = "authors_display"
        rows.sort(key=lambda row: (_sortable_text(row.get(key)), _canonical_for_sort(row)), reverse=descending)
    else:
        rows.sort(key=lambda row: (_sortable_text(row.get(sort_by)), _canonical_for_sort(row)), reverse=descending)
    return rows


def page_publication_snapshot(
    snapshot: dict[str, Any],
    page: int = 1,
    page_size: int = 20,
    sort_by: str | None = None,
    sort_order: str | None = None,
) -> dict[str, Any]:
    """Sort a full snapshot before slicing one page, capped at 20 rows."""

    requested_page = max(1, int(page))
    requested_size = max(1, min(20, int(page_size)))
    records = rank_publication_records(snapshot.get("records", []), sort_by or "data_readiness", sort_order)
    start = (requested_page - 1) * requested_size
    page_records = records[start : start + requested_size]
    out = {key: value for key, value in snapshot.items() if key != "records"}
    out.update(
        {
            "records": page_records,
            "page": requested_page,
            "page_size": requested_size,
            "total_records": len(records),
            "total": len(records),
            "total_pages": (len(records) + requested_size - 1) // requested_size if records else 0,
            "has_next": start + requested_size < len(records),
            "sort_by": sort_by or "data_readiness",
            "sort_order": sort_order or _default_sort_order(sort_by or "data_readiness"),
        }
    )
    return out


def search_publications(
    query: str,
    species: str | SpeciesSpec,
    limit: int = 1000,
    providers: Iterable[Any] | None = None,
    progress: Callable[[float, str], None] | None = None,
) -> dict[str, Any]:
    """Search providers, resolve candidates, merge paper records, and return a full snapshot.

    ``progress`` is an optional ``callback(fraction, message)`` where ``fraction``
    runs from 0.0 to 1.0 across this function's own work. Reporting is advisory:
    a failing callback never interrupts a search.
    """

    def report(fraction: float, message: str) -> None:
        if progress is None:
            return
        try:
            progress(min(1.0, max(0.0, float(fraction))), str(message))
        except Exception:  # noqa: BLE001 - progress reporting must never break a search.
            pass

    target = normalize_species(species)
    if isinstance(limit, bool):
        raise DiscoveryError("limit must be a whole number from 1 to 1000")
    try:
        requested_limit = int(limit)
    except (TypeError, ValueError) as exc:
        raise DiscoveryError("limit must be a whole number from 1 to 1000") from exc
    if requested_limit < 1:
        raise DiscoveryError("limit must be a whole number from 1 to 1000")
    evaluated_limit = min(MAX_SEARCH_RECORDS, requested_limit)
    terms = _query_terms(query)
    safe_query = " ".join(terms)
    active_providers = list(providers) if providers is not None else _default_providers()
    diagnostics: dict[str, Any] = {
        "provider_count": len(active_providers),
        "searched": [],
        "resolved": [],
        "errors": [],
        "evaluated_limit": evaluated_limit,
        "evaluated_records": 0,
    }
    raw_records: list[dict[str, Any]] = []
    provider_events: list[dict[str, Any]] = []
    searchable_providers = 0
    successful_searches = 0

    provider_total = max(1, len(active_providers))
    for provider_index, provider in enumerate(active_providers, start=1):
        provider_name = _provider_name(provider)
        report(
            0.02 + 0.43 * ((provider_index - 1) / provider_total),
            f"Querying {provider_name} ({provider_index} of {provider_total} sources)",
        )
        if not hasattr(provider, "search"):
            diagnostics["searched"].append({"provider": provider_name, "status": "skipped", "reason": "missing_search"})
            continue
        searchable_providers += 1
        try:
            found = provider.search(safe_query, target, evaluated_limit)
        except Exception as exc:  # pragma: no cover - defensive provider isolation
            message = _safe_provider_error(exc)
            diagnostics["errors"].append({"provider": provider_name, "stage": "search", "error": message})
            provider_events.append({"provider": provider_name, "event": "search", "status": "error", "message": message})
            continue
        successful_searches += 1
        found_records = _provider_records(found)
        search_records = [_with_provider_source(record, provider_name) for record in found_records[:evaluated_limit]]
        raw_records.extend(search_records)
        diagnostics["searched"].append(
            {
                "provider": provider_name,
                "status": "ok",
                "records": len(search_records),
                "returned_records": len(found_records),
                "truncated": len(found_records) > evaluated_limit,
                "global_cap_competition": len(raw_records) > evaluated_limit,
            }
        )
        provider_events.append(
            {
                "provider": provider_name,
                "event": "search",
                "status": "ok",
                "message": f"{len(search_records)} records",
            }
        )

    diagnostics["searchable_providers"] = searchable_providers
    diagnostics["successful_searches"] = successful_searches
    if searchable_providers == 0:
        raise DiscoveryUnavailableError("no publication search provider is available")
    if successful_searches == 0:
        raise DiscoveryUnavailableError("all publication search providers are unavailable; retry later")

    # Build one source-neutral publication universe before applying the global
    # cap.  This avoids a full PubMed page starving GEO records from the same
    # query while retaining deterministic provider-relevance ordering.
    report(0.46, f"Merging {len(raw_records)} provider records")
    merged_search_records = merge_publication_records(raw_records, target)
    preliminary = rank_publication_records(
        merged_search_records,
        sort_by="relevance",
        sort_order="desc",
    )[:evaluated_limit]
    report(0.50, f"Ranking {len(preliminary)} unique publications")

    resolved_records: list[dict[str, Any]] = list(preliminary)
    resolution_limit = min(len(preliminary), MAX_DETAILED_RESOLUTION_RECORDS)
    resolvers = [provider for provider in active_providers if hasattr(provider, "resolve")]
    for resolved_index, record in enumerate(preliminary[:resolution_limit], start=1):
        report(
            0.52 + 0.42 * ((resolved_index - 1) / max(1, resolution_limit)),
            f"Inspecting linked data {resolved_index} of {resolution_limit}",
        )
        record_events: list[dict[str, Any]] = []
        for provider in resolvers:
            provider_name = _provider_name(provider)
            try:
                resolved = provider.resolve(record, target)
            except Exception as exc:  # pragma: no cover - defensive provider isolation
                message = _safe_provider_error(exc)
                diagnostics["errors"].append({"provider": provider_name, "stage": "resolve", "error": message})
                event = {"provider": provider_name, "event": "resolve", "status": "error", "message": message}
                provider_events.append(event)
                record_events.append(event)
                continue
            bound_records, events = _bind_resolved_records(record, _provider_records(resolved), provider_name)
            resolved_records.extend(bound_records)
            provider_events.extend(events)
            record_events.extend(events)
            for event in events:
                if str(event.get("status") or "").lower() in {"error", "candidate_error", "unavailable"}:
                    diagnostics["errors"].append(
                        {
                            "provider": event.get("provider") or provider_name,
                            "stage": "resolve",
                            "error": _clean_text(event.get("message")) or "provider resolution failed",
                        }
                    )
            diagnostics["resolved"].append(
                {"provider": provider_name, "publication": record.get("canonical_id"), "records": len(bound_records)}
            )
        # Record the outcome so a later prepare can tell that this publication was
        # already resolved during search and skip repeating the provider calls.
        record["resolution_events"] = record_events
        record["resolution_state"] = _resolution_state(record_events, bool(resolvers))

    diagnostics["evaluated_records"] = len(preliminary)
    diagnostics["detailed_resolution_limit"] = resolution_limit
    diagnostics["detailed_resolution_truncated"] = len(preliminary) > resolution_limit
    report(0.95, "Scoring DEG-input readiness")
    merged = merge_publication_records(resolved_records, target)
    ranked = flag_shared_submission_records(rank_publication_records(merged))
    readiness_counts: dict[str, int] = {}
    for record in ranked:
        state = record.get("data_readiness", {}).get("verification_state", "unknown")
        readiness_counts[state] = readiness_counts.get(state, 0) + 1
    diagnostics["readiness_counts"] = readiness_counts
    return {
        "query": safe_query,
        "species": {"key": target.key, "label": target.label, "scientific_name": target.scientific_name},
        "limit": evaluated_limit,
        "records": ranked,
        "total": len(ranked),
        "total_records": len(ranked),
        "evaluated_records": len(preliminary),
        "ranking_limit": evaluated_limit,
        "ranking_truncated": len(merged_search_records) > evaluated_limit,
        "ranking_contract": "Partially resolved readiness, then provider relevance; 10 rows per browser page, with later selections resolved on demand.",
        "provider_status": "partial" if diagnostics["errors"] else "complete",
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "provider_events": provider_events,
        "diagnostics": diagnostics,
    }


# A publication whose resolution already produced a definite outcome does not
# need the provider calls repeated. "partial" and "unavailable" are retried
# because they indicate a provider that failed or answered incompletely.
SETTLED_RESOLUTION_STATES = frozenset({"resolved_candidates", "resolved_no_candidate"})


def _resolution_state(events: Iterable[dict[str, Any]], has_resolvers: bool) -> str:
    """Classify one publication's resolution outcome from its provider events."""

    collected = list(events)
    found = any(event.get("status") == "found" for event in collected)
    failed = any(event.get("status") in {"error", "candidate_error", "unavailable"} for event in collected)
    if found and not failed:
        return "resolved_candidates"
    if found:
        return "partial"
    if failed:
        return "unavailable"
    return "resolved_no_candidate" if has_resolvers else "no_resolver"


def resolve_publication_records(
    records: Iterable[dict[str, Any]],
    species: str | SpeciesSpec,
    providers: Iterable[Any] | None = None,
    reuse_settled: bool = True,
    progress: Callable[[float, str], None] | None = None,
) -> list[dict[str, Any]]:
    """Resolve missing direct-file routes for at most 20 user-selected papers.

    Search performs detailed public-repository resolution for the first visible
    page.  This bounded on-demand path gives later-page selections the same
    treatment without issuing thousands of provider calls during initial search.

    When ``reuse_settled`` is true, publications that search already resolved to a
    definite outcome are passed through untouched. Selecting the first page and
    preparing it used to repeat every provider call made during search.
    """

    target = normalize_species(species)
    selected = merge_publication_records(records, target)
    if len(selected) > 20:
        raise DiscoveryError("at most 20 publication records can be resolved at once")
    active_providers = list(providers) if providers is not None else _default_providers()
    resolvers = [provider for provider in active_providers if hasattr(provider, "resolve")]
    resolved: list[dict[str, Any]] = []

    def report(fraction: float, message: str) -> None:
        if progress is None:
            return
        try:
            progress(min(1.0, max(0.0, float(fraction))), str(message))
        except Exception:  # noqa: BLE001 - progress reporting must never break resolution.
            pass

    for publication_index, publication in enumerate(selected, start=1):
        report(
            (publication_index - 1) / max(1, len(selected)),
            f"Resolving publication {publication_index} of {len(selected)}",
        )
        if reuse_settled and str(publication.get("resolution_state") or "") in SETTLED_RESOLUTION_STATES:
            resolved.append(dict(publication))
            continue
        publication_events: list[dict[str, Any]] = []
        bound_records: list[dict[str, Any]] = []
        for provider in resolvers:
            provider_name = _provider_name(provider)
            try:
                response = provider.resolve(publication, target)
            except Exception as exc:
                publication_events.append(
                    {
                        "provider": provider_name,
                        "event": "resolve",
                        "status": "error",
                        "message": _safe_provider_error(exc),
                        "publication": publication.get("canonical_id", ""),
                    }
                )
                continue
            bound, events = _bind_resolved_records(publication, _provider_records(response), provider_name)
            bound_records.extend(bound)
            publication_events.extend(events)
        base = dict(publication)
        base["resolution_events"] = publication_events
        base["resolution_state"] = _resolution_state(publication_events, bool(resolvers))
        resolved.append(base)
        resolved.extend(bound_records)
    return rank_publication_records(merge_publication_records(resolved, target))


def _prepare_record(record: dict[str, Any], target: SpeciesSpec) -> dict[str, Any]:
    row = _copy_record(record)
    identifiers = row.get("identifiers") if isinstance(row.get("identifiers"), dict) else {}
    for key in ("pmid", "doi", "pmcid"):
        if not row.get(key) and identifiers.get(key):
            row[key] = identifiers[key]
    if not row.get("pubmed_ids") and row.get("pmids"):
        row["pubmed_ids"] = row.get("pmids")
    if not row.get("candidates") and row.get("supplementary_file_candidates"):
        row["candidates"] = row.get("supplementary_file_candidates")
    if row.get("relevance_rank") in (None, "") and row.get("rank") not in (None, ""):
        row["relevance_rank"] = row.get("rank")
    if row.get("quarantined") is True:
        row["mixed_quarantined"] = True
    row["pubmed_ids"] = sorted(set(filter(None, (_normalize_pmid(value) for value in _as_list(row.get("pubmed_ids")) + _as_list(row.get("pmid"))))))
    row["pmid"] = row["pubmed_ids"][0] if row["pubmed_ids"] else ""
    dois = sorted(set(filter(None, (_normalize_doi(value) for value in _as_list(row.get("doi")) + _as_list(row.get("dois"))))))
    row["doi"] = dois[0] if dois else ""
    if dois:
        row["dois"] = dois
    row["pmcid"] = _first_sorted(_normalize_pmcid(value) for value in _as_list(row.get("pmcid")) + _as_list(row.get("pmcids"))) or ""
    row["geo_accessions"] = _collect_accessions(row)
    row["accession"] = row["geo_accessions"][0] if row["geo_accessions"] else ""
    row["provider_ids"] = _collect_provider_ids(row)
    row["sources"] = _sorted_strings(_as_list(row.get("sources")) + _as_list(row.get("source")) + _as_list(row.get("provider")))
    row["candidates"] = _normalize_candidates(row)
    row["species_evidence"] = _species_evidence(row)
    row["species"] = _species_labels(row["species_evidence"])
    row["species_decision"] = _species_decision(row, target)
    row["mixed_rescued"] = row["species_decision"] == "mixed_rescued"
    row["mixed_quarantined"] = row["species_decision"] == "mixed_quarantined"
    row["data_readiness"] = _readiness(row, target)
    row["canonical_id"] = canonical_record_id(row)
    row["record_kind"] = row.get("record_kind") or "publication"
    row["paper_title"] = _clean_text(row.get("paper_title") or row.get("title"))
    row["authors"] = _normalize_authors(row.get("authors"))
    row["authors_display"] = row.get("authors_display") or _authors_display(row["authors"])
    row["journal"] = _clean_text(row.get("journal"))
    row["year"] = _year(row) or None
    row["relevance_rank"] = _number_or_inf(row.get("relevance_rank"))
    if row["relevance_rank"] == float("inf"):
        row["relevance_rank"] = None
    row["source_unit_id"] = _source_unit_id(row)
    row["shared_submission_units"] = _sorted_strings(_as_list(row.get("shared_submission_units")))
    row["shared_submission_warning"] = _clean_text(row.get("shared_submission_warning"))
    _add_ui_aliases(row)
    return row


def _merge_group(group: list[dict[str, Any]], target: SpeciesSpec) -> dict[str, Any]:
    ordered = sorted(group, key=_default_rank_key)
    base = ordered[0]
    merged: dict[str, Any] = {
        "canonical_id": canonical_record_id(base),
        "record_kind": "publication",
        "paper_title": _first_present(ordered, "paper_title"),
        "authors": _first_nonempty_list(ordered, "authors"),
        "journal": _first_present(ordered, "journal"),
        "year": max((_year(row) for row in ordered), default=0) or None,
        "pubmed_ids": _sorted_strings(value for row in ordered for value in _as_list(row.get("pubmed_ids")) + _as_list(row.get("pmid"))),
        "doi": _first_sorted(_normalize_doi(value) for row in ordered for value in _as_list(row.get("doi")) + _as_list(row.get("dois"))) or "",
        "pmcid": _first_sorted(_normalize_pmcid(value) for row in ordered for value in _as_list(row.get("pmcid")) + _as_list(row.get("pmcids"))) or "",
        "geo_accessions": _sorted_strings(value for row in ordered for value in _collect_accessions(row)),
        "provider_ids": _sorted_strings(value for row in ordered for value in _collect_provider_ids(row)),
        "sources": _sorted_strings(value for row in ordered for value in _as_list(row.get("sources"))),
        "candidates": [],
        "species_evidence": _merge_species_evidence(ordered),
        "relevance_rank": min((_number_or_inf(row.get("relevance_rank")) for row in ordered), default=float("inf")),
        "shared_submission_units": _sorted_strings(
            unit for row in ordered for unit in _as_list(row.get("shared_submission_units"))
        ),
        "shared_submission_warning": _first_text(row.get("shared_submission_warning") for row in ordered),
        "target_species_verified": any(_truthy(row.get("target_species_verified")) for row in ordered),
        "target_species_evidence": _first_present(ordered, "target_species_evidence")
        or _first_present(ordered, "mixed_rescue_evidence")
        or _first_present(ordered, "evidence_text"),
        "mixed_blocked": any(_truthy(row.get("mixed_blocked")) for row in ordered),
        "mixed_quarantined": any(
            _truthy(row.get("mixed_quarantined")) or _truthy(row.get("quarantined"))
            for row in ordered
        ),
        "resolution_events": _merge_resolution_events(ordered),
    }
    merged["pmid"] = merged["pubmed_ids"][0] if merged["pubmed_ids"] else ""
    merged["accession"] = merged["geo_accessions"][0] if merged["geo_accessions"] else ""
    merged["authors_display"] = _authors_display(merged["authors"])
    merged["species"] = _species_labels(merged["species_evidence"])
    merged["species_decision"] = _species_decision({**base, **merged}, target)
    merged["mixed_rescued"] = merged["species_decision"] == "mixed_rescued"
    merged["mixed_quarantined"] = merged["species_decision"] == "mixed_quarantined"
    merged["candidates"] = _merge_candidates(ordered)
    merged["data_readiness"] = _merged_readiness(ordered, merged, target)
    merged["source_unit_id"] = _merged_source_unit_id(ordered, merged)
    explicit_source_units = _sorted_strings(
        row.get("source_unit_id") for row in ordered if _clean_text(row.get("source_unit_id"))
    )
    merged["source_unit_conflict"] = explicit_source_units if len(explicit_source_units) > 1 else []
    resolution_states = _sorted_strings(row.get("resolution_state") for row in ordered)
    if resolution_states:
        has_resolution_error = any(
            str(event.get("status") or "").lower() in {"error", "candidate_error", "unavailable"}
            for event in merged["resolution_events"]
        )
        merged["resolution_state"] = (
            "resolved_candidates"
            if "resolved_candidates" in resolution_states and not has_resolution_error
            else "partial"
            if "resolved_candidates" in resolution_states or "partial" in resolution_states
            else "unavailable"
            if "unavailable" in resolution_states
            else resolution_states[0]
        )
    if merged["relevance_rank"] == float("inf"):
        merged["relevance_rank"] = None
    merged["canonical_id"] = canonical_record_id(merged)
    _add_ui_aliases(merged)
    return merged


def _readiness(record: dict[str, Any], target: SpeciesSpec) -> dict[str, Any]:
    existing = record.get("data_readiness")
    if isinstance(existing, dict):
        tier = existing.get("tier") or existing.get("verification_state") or "unknown"
        state = existing.get("verification_state") or tier
        basis = _sorted_strings(_as_list(existing.get("basis")))
    else:
        tier = str(existing or "").strip() or "unknown"
        state = tier
        basis = []
    decision = _species_decision(record, target)
    if decision == "mixed_quarantined":
        state = "mixed_quarantined"
        tier = "mixed_quarantined"
        basis.append("mixed_species_without_target_file_verification")
    elif (
        decision in {"target_species_verified", "mixed_rescued"}
        and _has_file_candidate(record)
    ) or (state in {"verified_ready", "target_species_verified"} and _has_file_candidate(record)):
        state = "verified_ready"
        tier = "verified_ready"
        basis.append("mixed_species_target_file_verified" if decision == "mixed_rescued" else "target_species_verified")
    elif _has_file_candidate(record):
        state = "likely_ready"
        tier = "likely_ready"
        basis.append("public_data_candidate_present_unverified")
    elif _collect_accessions(record):
        # A repository accession is not evidence that a usable table exists.
        # Treating it as `likely_ready` put every GEO row in the top tier, so
        # the primary sort key stopped discriminating and the first result was
        # routinely one whose files turned out to be browser tracks.
        state = "candidate"
        tier = "candidate"
        basis.append(
            "repository_record_not_inspected"
            if str(record.get("detail_assessment") or "") == "not_evaluated"
            else "repository_record_without_tabular_file"
        )
    elif not tier or tier == "unknown":
        state = "metadata_only"
        tier = "metadata_only"
        basis.append("publication_metadata_only")
    return {
        "tier": tier,
        "priority": _READINESS_PRIORITY.get(state, _READINESS_PRIORITY.get(tier, 5)),
        "verification_state": state,
        "basis": _sorted_strings(basis),
    }


def _merged_readiness(ordered: list[dict[str, Any]], merged: dict[str, Any], target: SpeciesSpec) -> dict[str, Any]:
    if merged["species_decision"] == "mixed_quarantined":
        return {
            "tier": "mixed_quarantined",
            "priority": _READINESS_PRIORITY["mixed_quarantined"],
            "verification_state": "mixed_quarantined",
            "basis": _sorted_strings(
                value
                for row in ordered
                for value in _as_list(row.get("data_readiness", {}).get("basis") if isinstance(row.get("data_readiness"), dict) else [])
            )
            + ["mixed_species_without_target_file_verification"],
        }
    readiness = sorted(
        [*(_readiness(row, target) for row in ordered), _readiness(merged, target)],
        key=lambda value: value["priority"],
    )
    best = dict(readiness[0])
    basis = _sorted_strings(value for item in readiness for value in _as_list(item.get("basis")))
    best["basis"] = basis
    return best


def _species_decision(record: dict[str, Any], target: SpeciesSpec) -> str:
    evidence = _species_evidence(record)
    labels = set(_species_labels(evidence))
    if _truthy(record.get("mixed_blocked")):
        return "mixed_quarantined"
    explicitly_quarantined = _truthy(record.get("mixed_quarantined")) or _truthy(record.get("quarantined"))
    is_mixed = len(labels) > 1 or explicitly_quarantined
    rescue_evidence = _clean_text(
        record.get("target_species_evidence")
        or record.get("mixed_rescue_evidence")
        or record.get("evidence_text")
    )
    if _truthy(record.get("target_species_verified")) and is_mixed:
        return "mixed_rescued" if rescue_evidence else "mixed_quarantined"
    if explicitly_quarantined:
        return "mixed_quarantined"
    if _truthy(record.get("target_species_verified")):
        return "target_species_verified"
    if not labels:
        return "unknown"
    if target.label in labels and len(labels) == 1:
        return "target_species_likely"
    if target.label in labels:
        return "mixed_quarantined"
    return "non_target"


def _species_evidence(record: dict[str, Any]) -> list[dict[str, str]]:
    evidence: list[dict[str, str]] = []
    raw = record.get("species_evidence")
    if isinstance(raw, list):
        for item in raw:
            if isinstance(item, dict):
                label = _species_label(item.get("species") or item.get("label") or item.get("name") or item.get("taxon_id"))
                if label:
                    evidence.append({"species": label, "basis": _clean_text(item.get("basis")) or "provider"})
            else:
                label = _species_label(item)
                if label:
                    evidence.append({"species": label, "basis": "provider"})
    elif isinstance(raw, dict):
        basis = _clean_text(raw.get("basis")) or "provider"
        observed = _as_list(raw.get("observed_taxa"))
        requested = _as_list(raw.get("requested")) if not observed else []
        for value in observed + requested:
            label = _species_label(value)
            if label:
                evidence.append({"species": label, "basis": basis})
    for key in ("species", "organism", "taxon_id"):
        for value in _as_list(record.get(key)):
            label = _species_label(value)
            if label:
                evidence.append({"species": label, "basis": key})
    seen: set[tuple[str, str]] = set()
    unique: list[dict[str, str]] = []
    for item in evidence:
        marker = (item["species"], item["basis"])
        if marker not in seen:
            seen.add(marker)
            unique.append(item)
    return sorted(unique, key=lambda item: (item["species"], item["basis"]))


def _species_label(value: Any) -> str:
    if value is None:
        return ""
    raw = re.sub(r"\s+", " ", str(value).strip())
    text = raw.lower().replace("_", " ").replace("-", " ")
    if text in _MISSING:
        return ""
    if text in {"human", "homo sapiens", "hs", "h sapiens", "9606"}:
        return "Human"
    if text in {"mouse", "mus musculus", "mm", "m musculus", "10090"}:
        return "Mouse"
    if text.startswith("other:"):
        other = raw.split(":", 1)[1].strip()
        return f"Other:{other[:1].upper() + other[1:].lower()}" if other else ""
    if text.startswith("taxon:"):
        taxon = raw.split(":", 1)[1].strip()
        return f"Taxon:{taxon}" if taxon else ""
    if text in {"unknown", "not reported", "not available", "mixed", "multiple", "other"}:
        return ""
    if text.isdigit():
        return f"Taxon:{text}"
    # Preserve non-target taxa so Human+primate and Mouse+other-species records
    # cannot collapse into an apparently single-species result.
    return f"Other:{raw[:1].upper() + raw[1:].lower()}"


def _species_labels(evidence: Iterable[dict[str, str]]) -> list[str]:
    return sorted({item["species"] for item in evidence if item.get("species")})


def _graph_identifiers(record: dict[str, Any]) -> list[str]:
    identifiers: list[str] = []
    identifiers.extend(f"pmid:{normalized}" for value in _as_list(record.get("pubmed_ids")) + _as_list(record.get("pmid")) if (normalized := _normalize_pmid(value)))
    identifiers.extend(f"doi:{normalized}" for value in _as_list(record.get("doi")) + _as_list(record.get("dois")) if (normalized := _normalize_doi(value)))
    identifiers.extend(f"pmcid:{normalized}" for value in _as_list(record.get("pmcid")) + _as_list(record.get("pmcids")) if (normalized := _normalize_pmcid(value)))
    identifiers.extend(f"accession:{value}" for value in _collect_accessions(record))
    identifiers.extend(f"provider:{value}" for value in _collect_provider_ids(record))
    return _sorted_strings(identifiers)


def _collect_accessions(record: dict[str, Any]) -> list[str]:
    values = []
    for key in ("geo_accessions", "accessions", "accession", "gse", "geo"):
        values.extend(_as_list(record.get(key)))
    candidates = record.get("candidates")
    if isinstance(candidates, list):
        for candidate in candidates:
            if isinstance(candidate, dict):
                values.extend(_as_list(candidate.get("accession")) + _as_list(candidate.get("geo_accessions")))
    return _sorted_strings(_normalize_accession(value) for value in values if _normalize_accession(value))


def _collect_provider_ids(record: dict[str, Any]) -> list[str]:
    values: list[str] = []
    raw = record.get("provider_ids")
    if isinstance(raw, dict):
        for provider, ids in raw.items():
            for value in _as_list(ids):
                cleaned = _clean_text(value)
                if cleaned:
                    values.append(f"{_clean_text(provider).lower()}:{cleaned}")
    else:
        for value in _as_list(raw):
            cleaned = _clean_text(value)
            if cleaned:
                values.append(cleaned if ":" in cleaned else f"provider:{cleaned}")
    provider = _clean_text(record.get("provider") or record.get("source")).lower()
    for key in ("provider_id", "id", "uid"):
        value = _clean_text(record.get(key))
        if provider and value:
            values.append(f"{provider}:{value}")
    return _sorted_strings(values)


def _normalize_candidates(record: dict[str, Any]) -> list[dict[str, Any]]:
    raw = record.get("candidates")
    if not isinstance(raw, list):
        raw = []
    supplements = record.get("supplementary_file_candidates")
    if isinstance(supplements, list):
        raw = [*raw, *supplements]
    candidates = [dict(item) for item in raw if isinstance(item, dict)]
    for candidate in candidates:
        if not candidate.get("source_url") and candidate.get("url"):
            candidate["source_url"] = candidate["url"]
        candidate.setdefault("provider", record.get("provider") or record.get("source") or "")
        marker = "|".join(
            str(candidate.get(key) or "")
            for key in ("provider", "accession", "source_url", "name", "role")
        )
        candidate.setdefault("candidate_id", hashlib.sha256(marker.encode("utf-8")).hexdigest()[:16])
    for accession in _collect_accessions(record):
        candidates.append({"accession": accession, "source": "record"})
    return _merge_candidate_list(candidates)


def _has_file_candidate(record: dict[str, Any]) -> bool:
    for candidate in _normalize_candidates(record):
        if _clean_text(candidate.get("source_url") or candidate.get("url")):
            return True
        if _clean_text(candidate.get("name") or candidate.get("filename")) and _clean_text(candidate.get("role")):
            return True
    return False


def _merge_candidates(ordered: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return _merge_candidate_list(candidate for row in ordered for candidate in _normalize_candidates(row))


def _merge_candidate_list(candidates: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        accession = _normalize_accession(candidate.get("accession") or candidate.get("id") or "")
        source_url = _clean_text(candidate.get("source_url") or candidate.get("url"))
        key = repr(
            (
                accession,
                _clean_text(candidate.get("provider") or candidate.get("source")).lower(),
                source_url,
                _clean_text(candidate.get("name") or candidate.get("filename")).lower(),
                _clean_text(candidate.get("role")).lower(),
                _clean_text(candidate.get("record_type") or candidate.get("record_kind")).lower(),
            )
        )
        row = merged.setdefault(key, {})
        for field, value in candidate.items():
            if field == "accession" and accession:
                row[field] = accession
            elif field not in row or not row[field]:
                row[field] = value
        if accession:
            row["accession"] = accession
        if not row.get("source_url") and row.get("url"):
            row["source_url"] = row["url"]
        marker = "|".join(
            str(row.get(field) or "")
            for field in ("provider", "accession", "source_url", "name", "role")
        )
        row.setdefault("candidate_id", hashlib.sha256(marker.encode("utf-8")).hexdigest()[:16])
    return [merged[key] for key in sorted(merged)]


def _merge_resolution_events(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for record in records:
        for event in _as_list(record.get("resolution_events")):
            if not isinstance(event, dict):
                continue
            row = dict(event)
            marker = repr(sorted((str(key), repr(value)) for key, value in row.items()))
            merged.setdefault(marker, row)
    return [merged[key] for key in sorted(merged)]


def _first_text(values: Iterable[Any]) -> str:
    for value in values:
        text = _clean_text(value)
        if text:
            return text
    return ""


def _submission_title_key(record: dict[str, Any]) -> str:
    text = _clean_text(record.get("paper_title") or record.get("title"))
    key = re.sub(r"[^a-z0-9]+", "", text.lower())
    # Short titles collide by accident; a submission title does not.
    return key if len(key) >= 24 else ""


def flag_shared_submission_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Warn where separate repository records are almost certainly one submission.

    Source units collapse on a shared PubMed ID. A submission that has not been
    published yet has no PMID to collapse on, so its arms - deposited as
    separate GEO series under one title - each count as an independent source
    unit. DEGORA's whole replication claim rests on that count, so the reader
    has to be told before treating them as two studies.
    """

    groups: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        if _as_list(record.get("pubmed_ids")) or _clean_text(record.get("pmid")):
            continue  # A publication link would already have collapsed these.
        key = _submission_title_key(record)
        if key:
            groups.setdefault(key, []).append(record)

    for members in groups.values():
        units = _sorted_strings(record.get("source_unit_id") for record in members)
        if len(members) < 2 or len(units) < 2:
            continue
        for record in members:
            own = _clean_text(record.get("source_unit_id"))
            record["shared_submission_units"] = [unit for unit in units if unit != own]
            record["shared_submission_warning"] = (
                "shares its title with "
                f"{len(units) - 1} other repository record"
                f"{'' if len(units) == 2 else 's'} and none is linked to a publication, "
                "so these may be one submission rather than independent studies"
            )
    return records


def _source_unit_id(record: dict[str, Any]) -> str:
    explicit = _clean_text(record.get("source_unit_id"))
    if explicit:
        return explicit
    pmid = _first_sorted(_normalize_pmid(value) for value in _as_list(record.get("pubmed_ids")) + _as_list(record.get("pmid")))
    if pmid:
        return f"PMID:{pmid}"
    doi = _first_sorted(_normalize_doi(value) for value in _as_list(record.get("doi")) + _as_list(record.get("dois")))
    if doi:
        return f"DOI:{doi}"
    accessions = _collect_accessions(record)
    if accessions:
        return accessions[0]
    return canonical_record_id(record)


def _merged_source_unit_id(ordered: list[dict[str, Any]], merged: dict[str, Any]) -> str:
    explicit = _sorted_strings(row.get("source_unit_id") for row in ordered if _clean_text(row.get("source_unit_id")))
    if len(explicit) == 1:
        return explicit[0]
    if len(merged.get("geo_accessions", [])) == 1:
        return merged["geo_accessions"][0]
    return merged["canonical_id"]


def _default_rank_key(record: dict[str, Any]) -> tuple[int, float, int, str]:
    readiness = record.get("data_readiness")
    priority = readiness.get("priority", 5) if isinstance(readiness, dict) else _READINESS_PRIORITY.get(str(readiness), 5)
    return (int(priority), _number_or_inf(record.get("relevance_rank")), -_year(record), _canonical_for_sort(record))


def _default_sort_order(sort_by: str | None) -> str:
    return "desc" if sort_by in {None, "", "data_readiness", "readiness", "readiness_score", "relevance", "relevance_rank", "relevance_score", "year"} else "asc"


def _canonical_for_sort(record: dict[str, Any]) -> str:
    return _clean_text(record.get("canonical_id")) or canonical_record_id(record)


def _normalize_doi(value: Any) -> str:
    text = _clean_text(value).lower()
    if not text:
        return ""
    text = re.sub(r"^doi:\s*", "", text)
    text = re.sub(r"^https?://(dx\.)?doi\.org/", "", text)
    return text.strip().rstrip(".")


def _normalize_pmid(value: Any) -> str:
    text = _clean_text(value)
    match = re.search(r"\d+", text)
    return match.group(0) if match else ""


def _normalize_pmcid(value: Any) -> str:
    text = _clean_text(value).upper()
    if not text:
        return ""
    match = re.search(r"(?:PMC)?(\d+)", text)
    return f"PMC{match.group(1)}" if match else text


def _normalize_accession(value: Any) -> str:
    text = _clean_text(value).upper()
    if not text:
        return ""
    return re.sub(r"\s+", "", text)


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set, frozenset)):
        return list(value)
    return [value]


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in _MISSING else text


def _safe_provider_error(exc: Exception) -> str:
    """Return bounded provider diagnostics without exposing request secrets."""

    text = re.sub(
        r"(?i)(api[_-]?key|access[_-]?token|token|password|secret)=([^&\s]+)",
        r"\1=[redacted]",
        str(exc),
    )
    text = re.sub(r"[\r\n\t]+", " ", text).strip()
    return f"{type(exc).__name__}: {(text or 'provider request failed')[:240]}"


def _sorted_strings(values: Iterable[Any]) -> list[str]:
    return sorted({str(value) for value in values if _clean_text(value)})


def _first_sorted(values: Iterable[str]) -> str:
    items = _sorted_strings(values)
    return items[0] if items else ""


def _first_present(records: list[dict[str, Any]], key: str) -> str:
    values = sorted((_clean_text(record.get(key)) for record in records if _clean_text(record.get(key))), key=lambda value: (-len(value), value.lower()))
    return values[0] if values else ""


def _first_nonempty_list(records: list[dict[str, Any]], key: str) -> list[str]:
    values = [record.get(key) for record in records if _as_list(record.get(key))]
    return _normalize_authors(values[0]) if values else []


def _normalize_authors(value: Any) -> list[str]:
    if isinstance(value, str):
        parts = [part.strip() for part in re.split(r";|,", value) if part.strip()]
        return parts
    return [_clean_text(item) for item in _as_list(value) if _clean_text(item)]


def _authors_display(authors: list[str]) -> str:
    if not authors:
        return ""
    if len(authors) == 1:
        return authors[0]
    return f"{authors[0]} et al."


def _year(record: dict[str, Any]) -> int:
    value = record.get("year") or record.get("publication_year") or record.get("publication_date")
    try:
        return int(value)
    except (TypeError, ValueError):
        match = re.search(r"(19|20)\d{2}", str(value or ""))
        return int(match.group(0)) if match else 0


def _number_or_inf(value: Any) -> float:
    try:
        if value is None:
            return float("inf")
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "verified", "target_species_verified"}


def _copy_record(record: dict[str, Any]) -> dict[str, Any]:
    return dict(record or {})


def _sortable_text(value: Any) -> str:
    if isinstance(value, (list, tuple, set, frozenset)):
        return ", ".join(str(item) for item in value).lower()
    if isinstance(value, dict):
        return str(value.get("verification_state") or value.get("tier") or "").lower()
    return _clean_text(value).lower()


def _add_ui_aliases(row: dict[str, Any]) -> None:
    canonical = _clean_text(row.get("canonical_id")) or canonical_record_id(row)
    row["publication_id"] = canonical
    row["title"] = _clean_text(row.get("paper_title") or row.get("title"))
    readiness = row.get("data_readiness") if isinstance(row.get("data_readiness"), dict) else {}
    row["readiness"] = readiness.get("verification_state") or readiness.get("tier") or "unknown"
    row["readiness_score"] = MAX_SEARCH_RECORDS - int(readiness.get("priority", 5))
    row["readiness_verified"] = row["readiness"] == "verified_ready"
    row["relevance"] = row.get("relevance_rank")
    row["data_sources"] = _sorted_strings(
        [*_as_list(row.get("sources")), *_collect_accessions(row)]
    )
    row["linked_datasets"] = list(_collect_accessions(row))
    pmid = _clean_text(row.get("pmid"))
    doi = _clean_text(row.get("doi"))
    row["pubmed_url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
    row["doi_url"] = f"https://doi.org/{doi}" if doi else ""
    row["source_url"] = row["pubmed_url"] or row["doi_url"] or _clean_text(row.get("source_url"))


def _merge_species_evidence(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    markers = {
        (item["species"], item["basis"])
        for record in records
        for item in _species_evidence(record)
    }
    return [{"species": species, "basis": basis} for species, basis in sorted(markers)]


def _find(parent: list[int], index: int) -> int:
    while parent[index] != index:
        parent[index] = parent[parent[index]]
        index = parent[index]
    return index


def _union(parent: list[int], left: int, right: int) -> None:
    left_root = _find(parent, left)
    right_root = _find(parent, right)
    if left_root != right_root:
        parent[max(left_root, right_root)] = min(left_root, right_root)


def _default_providers() -> list[Any]:
    from degora.discovery_sources import default_publication_providers

    return list(default_publication_providers())


def _provider_name(provider: Any) -> str:
    return _clean_text(getattr(provider, "name", "") or provider.__class__.__name__)


def _provider_records(value: Any) -> list[dict[str, Any]]:
    if value is None:
        return []
    if hasattr(value, "records"):
        value = getattr(value, "records")
    elif isinstance(value, dict) and isinstance(value.get("records"), list):
        value = value["records"]
    elif isinstance(value, dict):
        value = [value]
    try:
        items = list(value)
    except TypeError as exc:
        raise DiscoveryError("discovery provider must return records") from exc
    return [dict(item) for item in items if isinstance(item, dict)]


def _bind_resolved_records(
    publication: dict[str, Any],
    resolved: list[dict[str, Any]],
    provider_name: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Attach source-neutral resolver output to the publication it resolved."""

    records: list[dict[str, Any]] = []
    events: list[dict[str, Any]] = []
    for item in resolved:
        kind = _clean_text(item.get("record_type") or item.get("record_kind")).lower()
        if kind == "diagnostics":
            for diagnostic in _as_list(item.get("diagnostics")):
                if isinstance(diagnostic, dict):
                    events.append(
                        {
                            "provider": diagnostic.get("provider") or provider_name,
                            "event": "resolve",
                            "status": diagnostic.get("status") or "candidate_error",
                            "message": diagnostic.get("error") or diagnostic.get("message") or "",
                            "publication": publication.get("canonical_id", ""),
                        }
                    )
            continue
        if kind == "file_candidate":
            candidate = dict(item)
            candidate.setdefault("provider", provider_name)
            if not candidate.get("source_url") and candidate.get("url"):
                candidate["source_url"] = candidate["url"]
            linked = dict(publication)
            linked["provider"] = provider_name
            linked["sources"] = _sorted_strings([*_as_list(publication.get("sources")), provider_name])
            linked["candidates"] = [*_normalize_candidates(publication), candidate]
            records.append(linked)
            events.append(
                {
                    "provider": provider_name,
                    "event": "candidate",
                    "status": "found",
                    "message": candidate.get("name") or candidate.get("source_url") or "public file candidate",
                    "publication": publication.get("canonical_id", ""),
                }
            )
            continue
        linked = dict(publication)
        for key, value in item.items():
            if value not in (None, "", [], {}):
                linked[key] = value
        linked["provider"] = provider_name
        linked["sources"] = _sorted_strings([*_as_list(publication.get("sources")), *_as_list(item.get("sources")), provider_name])
        records.append(linked)
    return records, events


def _with_provider_source(record: dict[str, Any], provider_name: str) -> dict[str, Any]:
    row = dict(record)
    row.setdefault("provider", provider_name)
    row["sources"] = _sorted_strings(_as_list(row.get("sources")) + [provider_name])
    return row


def _provider_for_record(record: dict[str, Any], providers: list[Any]) -> Any | None:
    source_names = set(_as_list(record.get("sources")) + _as_list(record.get("provider")))
    for provider in providers:
        if _provider_name(provider) in source_names:
            return provider
    return None
