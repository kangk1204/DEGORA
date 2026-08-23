"""Prepare canonical publication discovery records for DEGORA analysis."""

from __future__ import annotations

import hashlib
import io
import json
import re
import shutil
import stat
import tempfile
import urllib.parse
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable

from .discovery import (
    DISCOVERY_BUNDLE_ARTIFACT_TYPE,
    DISCOVERY_BUNDLE_FORMAT_VERSION,
    DISCOVERY_BUNDLE_MARKER,
    DiscoveryError,
    DiscoveryUnavailableError,
    DiscoveryUnsafeArchiveError,
    classify_filename,
    export_discovery_bundle,
    inspect_candidate_bytes,
    inspect_upstream_bytes,
    normalize_species,
    prepare_geo_studies,
)
from .discovery_federated import canonical_record_id, flag_shared_submission_records
from .discovery_sources import (
    MAX_ARCHIVE_DEPTH,
    MAX_ARCHIVE_EXPANDED_BYTES,
    MAX_ARCHIVE_MEMBER_BYTES,
    MAX_ARCHIVE_MEMBERS,
    describe_unexpected_payload,
    download_public_candidate,
)


_TABULAR_MEMBER_RE = re.compile(r"\.(csv|tsv|txt)(\.gz)?$|\.xlsx$", re.IGNORECASE)
_GSE_RE = re.compile(r"^GSE[0-9]+$", re.IGNORECASE)


def prepare_publication_records(
    records: Iterable[dict[str, Any]],
    species: str,
    *,
    query: str = "",
    materialize_dir: str | Path,
    max_records: int = 20,
    max_files_per_record: int = 6,
    transport: Any | None = None,
    geo_client: Any | None = None,
    force: bool = False,
    progress: Callable[[float, str], None] | None = None,
) -> dict[str, Any]:
    """Transactionally prepare publication records for ``run_discovery_analysis``.

    The returned bundle is review-only: author DEG tables and upstream matrices
    are inspected and materialized, but no contrast is activated automatically.

    ``progress`` is an optional ``callback(fraction, message)`` over this
    function's own work. Reporting is advisory and never interrupts a run.
    """

    def report(fraction: float, message: str) -> None:
        if progress is None:
            return
        try:
            progress(min(1.0, max(0.0, float(fraction))), str(message))
        except Exception:  # noqa: BLE001 - progress reporting must never break a prepare.
            pass

    spec = normalize_species(species)
    if isinstance(max_records, bool) or not isinstance(max_records, int) or not 1 <= max_records <= 20:
        raise DiscoveryError("max_records must be a whole number from 1 to 20")
    if (
        isinstance(max_files_per_record, bool)
        or not isinstance(max_files_per_record, int)
        or not 1 <= max_files_per_record <= 12
    ):
        raise DiscoveryError("max_files_per_record must be a whole number from 1 to 12")

    selected, initially_excluded = _select_unique_records(records, max_records=max_records)
    # Recomputed rather than carried: the flag is set on the search snapshot,
    # and re-merging on the way here drops fields the merge does not know about.
    flag_shared_submission_records(selected)
    target = Path(materialize_dir).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    _validate_target(target, force=force)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}.prepare-", dir=target.parent))
    try:
        geo_records, direct_records, excluded = _partition_records(selected, spec)
        excluded = [*initially_excluded, *excluded]
        result = _prepare_into_staging(
            geo_records=geo_records,
            direct_records=direct_records,
            spec=spec,
            query=query,
            staging=staging,
            max_files_per_record=max_files_per_record,
            transport=transport,
            geo_client=geo_client,
            excluded=excluded,
            report=report,
        )
        report(0.92, "Writing the prepared bundle")
        _retarget_paths(result, staging, target)
        result["materialize_dir"] = str(target)
        # Record where the bundle will live, not where it is being staged. The
        # audit JSON is written into staging and then published under `target`,
        # and the staging directory is removed moments later -- so a staging path
        # baked into the persisted document is a link that never resolves for the
        # reader who opens it.
        result["exports"] = _export_paths(target)
        export_discovery_bundle(result, staging, force=True)
        _write_marker(staging, spec.key)
        _publish_prepared_bundle(staging, target, force=force)
        report(1.0, "Prepared bundle published")
        return result
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def _prepare_into_staging(
    *,
    geo_records: list[dict[str, Any]],
    direct_records: list[dict[str, Any]],
    spec: Any,
    query: str,
    staging: Path,
    max_files_per_record: int,
    transport: Any | None,
    geo_client: Any | None,
    excluded: list[dict[str, Any]],
    report: Callable[[float, str], None] = lambda _fraction, _message: None,
) -> dict[str, Any]:
    studies: list[dict[str, Any]] = []
    selected_accessions: list[str] = []
    total_units = max(1, len(geo_records) + len(direct_records))
    if geo_records:
        accessions = _unique_gse_accessions(geo_records)
        selected_accessions = accessions
        record_by_accession = _record_by_geo_accession(geo_records)
        # Several selected publications can report the same GEO series. That is
        # one dataset, so it stays one study - but the publications it absorbs
        # used to vanish from the bundle entirely, leaving a selection of 20
        # that reconciled to neither the prepared nor the excluded list.
        for geo_record in geo_records:
            covered = _geo_accessions(geo_record)
            if any(record_by_accession.get(accession) is geo_record for accession in covered):
                continue
            duplicate = covered[0] if covered else ""
            owner = record_by_accession.get(duplicate) or {}
            excluded.append(
                _excluded(
                    geo_record,
                    f"repository series {duplicate} is already prepared from {_record_id(owner)}; "
                    "publications reporting the same series count as one source unit",
                )
            )
        # The reader selected publications; this phase downloads series. Saying
        # only "12 repository record(s)" against a selection of 20 reads like a
        # miscount.
        report(
            0.05,
            f"Downloading {len(accessions)} repository series "
            f"linked by {len(geo_records)} of {len(geo_records) + len(direct_records)} selected publications",
        )
        try:
            geo_result = prepare_geo_studies(
                accessions,
                spec.key,
                query=query,
                inspection_budget=max(1, len(accessions) * max_files_per_record),
                max_files_per_study=max_files_per_record,
                materialize_dir=staging,
                client=geo_client,
                force=True,
            )
        except (DiscoveryError, DiscoveryUnavailableError) as exc:
            # The repository half failing is not a reason to throw away the
            # publication half the reader also selected.
            for geo_record in geo_records:
                excluded.append(_excluded(geo_record, f"repository preparation failed: {exc}"))
        else:
            excluded.extend(geo_result.get("excluded_studies", []))
            studies.extend(_augment_geo_studies(geo_result.get("studies", []), record_by_accession))

    geo_share = len(geo_records) / total_units
    report(
        0.05 + 0.85 * geo_share,
        f"Inspecting {len(direct_records)} publication(s) with a directly linked file",
    )
    for direct_index, record in enumerate(direct_records, start=1):
        report(
            0.05 + 0.85 * ((len(geo_records) + direct_index - 1) / total_units),
            f"Downloading and inspecting {direct_index} of {len(direct_records)}",
        )
        try:
            study = _prepare_direct_record(record, spec, staging, max_files_per_record, transport)
        except DiscoveryUnsafeArchiveError:
            raise
        except (DiscoveryError, DiscoveryUnavailableError) as exc:
            excluded.append(_excluded(record, f"preparation failed for this record: {exc}"))
            continue
        if study is None:
            excluded.append(_excluded(record, "publication record has no public tabular or archive candidate"))
        elif study["files"]:
            studies.append(study)
        else:
            candidate_errors = list(study.get("candidate_errors") or [])
            excluded.append(
                _excluded(record, _exclusion_reason(candidate_errors), candidate_errors=candidate_errors)
            )

    result = {
        "query": str(query or "").strip(),
        "species": {"key": spec.key, "label": spec.label, "scientific_name": spec.scientific_name},
        "selected_accessions": selected_accessions,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "Canonical federated publication records",
        "returned_studies": len(studies),
        "materialize_dir": str(staging),
        "studies": studies,
        "excluded_studies": excluded,
        "safety_and_review": {
            "auto_run": False,
            "cross_species_pooling": False,
            "max_records": 20,
            "max_files_per_record": max_files_per_record,
            "fetch_scope": "full",
            "archive_policy": (
                "ZIP members are path-normalized, symlinks are rejected, nested ZIP depth is capped at 2, "
                "and only tabular members are materialized."
            ),
            "mixed_species_policy": (
                "mixed_quarantined and mixed_blocked records are excluded; mixed_rescued requires explicit "
                "target_species_verified evidence."
            ),
        },
        "analysis_policy": {
            "auto_run": False,
            "draft_rows_included": False,
            "cross_species_pooling": False,
            "minimum_independent_source_units": 2,
            "source_unit_rule": "Publication records sharing a source_unit_id count as one source unit.",
            "author_table_gate": "Confirm contrast direction and table scope before activation.",
            "upstream_gate": "Choose control/treatment samples and verify biological replicates before fallback analysis.",
            "note": "Human and Mouse catalogs and DEGORA runs are always generated separately.",
        },
    }
    return result


def _select_unique_records(records: Iterable[dict[str, Any]], *, max_records: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    selected: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    seen: set[str] = set()
    for record in records:
        if not isinstance(record, dict):
            excluded.append({"canonical_id": "", "reason": "record is not an object"})
            continue
        key = _record_id(record)
        if key in seen:
            excluded.append(_excluded(record, "duplicate canonical publication record"))
            continue
        seen.add(key)
        if len(selected) >= max_records:
            excluded.append(_excluded(record, "record limit exceeded"))
            continue
        selected.append(record)
    return selected, excluded


def _partition_records(records: list[dict[str, Any]], spec: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    geo_records: list[dict[str, Any]] = []
    direct_records: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for record in records:
        reason = _species_exclusion_reason(record, spec)
        if reason:
            excluded.append(_excluded(record, reason))
            continue
        if _geo_accessions(record):
            geo_records.append(record)
        elif _candidate_records(record):
            direct_records.append(record)
        else:
            excluded.append(_excluded(record, "publication record has no GEO accession or direct file candidate"))
    return geo_records, direct_records, excluded


def _prepare_direct_record(
    record: dict[str, Any],
    spec: Any,
    staging: Path,
    max_files_per_record: int,
    transport: Any | None,
) -> dict[str, Any] | None:
    candidates = _candidate_records(record)[:max_files_per_record]
    if not candidates:
        return None
    record_id = _record_id(record)
    study_dir = staging / "public_files" / _safe_name(record_id)
    files: list[dict[str, Any]] = []
    candidate_errors: list[dict[str, str]] = []
    for index, candidate in enumerate(candidates, start=1):
        downloaded_path = study_dir / f"source_{index:02d}_{_candidate_name(candidate)}"
        # One unreadable supplementary file must cost its own study at most.
        # DiscoveryError and DiscoveryUnavailableError are siblings, so both
        # have to be named here or a rejected archive escapes and takes the
        # whole selection down with it.
        try:
            downloaded = download_public_candidate(candidate, downloaded_path, transport=transport)
            source_url = str(downloaded["url"])
            payload = downloaded_path.read_bytes()
            if downloaded_path.suffix.lower() == ".zip":
                notes: list[str] = []
                files.extend(
                    _materialize_archive_tables(
                        record=record,
                        archive_path=downloaded_path,
                        payload=payload,
                        source_url=source_url,
                        output_dir=study_dir / f"source_{index:02d}_archive",
                        max_files=max_files_per_record - len(files),
                        notes=notes,
                    )
                )
                for note in notes:
                    candidate_errors.append(
                        {
                            "candidate": _candidate_name(candidate),
                            "provider": str(candidate.get("provider") or candidate.get("source") or "public source"),
                            "status": "rejected",
                            "error": note,
                        }
                    )
            else:
                file_record = _file_entry(
                    record=record,
                    name=downloaded_path.name,
                    source_url=source_url,
                    local_path=downloaded_path,
                    payload=payload,
                    declared_role=str(candidate.get("role") or ""),
                )
                if file_record is not None:
                    files.append(file_record)
        except DiscoveryUnsafeArchiveError:
            # A hostile archive is not a flaky download; refuse the whole run.
            raise
        except (DiscoveryError, DiscoveryUnavailableError) as exc:
            candidate_errors.append(
                {
                    "candidate": _candidate_name(candidate),
                    "provider": str(candidate.get("provider") or candidate.get("source") or "public source"),
                    # "unavailable" is worth retrying later; "rejected" never will be.
                    "status": "unavailable" if isinstance(exc, DiscoveryUnavailableError) else "rejected",
                    "error": str(exc),
                }
            )
            downloaded_path.unlink(missing_ok=True)
            shutil.rmtree(study_dir / f"source_{index:02d}_archive", ignore_errors=True)
            continue
        if len(files) >= max_files_per_record:
            break

    if not files:
        return {
            **_publication_study_metadata(record, spec),
            "files": [],
            "candidate_file_count": 0,
            "ready_for_review_count": 0,
            "upstream_matrix_count": 0,
            "preparation_status": _preparation_status(candidate_errors),
            "candidate_errors": candidate_errors,
        }
    ready = sum(item.get("inspection", {}).get("status") == "ready_for_review" for item in files)
    upstream = sum(item.get("inspection", {}).get("status") == "upstream_matrix_ready_for_contrast" for item in files)
    status = (
        "author_deg_ready_for_contrast_review"
        if ready
        else "upstream_matrix_requires_contrast"
        if upstream
        else "author_table_requires_mapping_review"
    )
    return {
        **_publication_study_metadata(record, spec),
        "files": files,
        "candidate_file_count": len(files),
        "ready_for_review_count": ready,
        "upstream_matrix_count": upstream,
        "preparation_status": status,
        "candidate_errors": candidate_errors,
    }


def _materialize_archive_tables(
    *,
    record: dict[str, Any],
    archive_path: Path,
    payload: bytes,
    source_url: str,
    output_dir: Path,
    max_files: int,
    notes: list[str] | None = None,
) -> list[dict[str, Any]]:
    files: list[dict[str, Any]] = []
    materialized_names: set[str] = set()
    total_members = 0
    total_expanded = 0
    nested_notes = notes if notes is not None else []

    def visit(data: bytes, prefix: str, depth: int) -> None:
        nonlocal total_members, total_expanded
        if depth > MAX_ARCHIVE_DEPTH:
            raise DiscoveryUnsafeArchiveError("nested ZIP depth exceeds the safety limit")
        try:
            with zipfile.ZipFile(io.BytesIO(data)) as archive:
                for info in archive.infolist():
                    if info.is_dir():
                        continue
                    _validate_archive_member(info)
                    total_members += 1
                    if total_members > MAX_ARCHIVE_MEMBERS:
                        raise DiscoveryUnsafeArchiveError("archive contains too many members")
                    total_expanded += info.file_size
                    if total_expanded > MAX_ARCHIVE_EXPANDED_BYTES:
                        raise DiscoveryUnsafeArchiveError("archive expanded-size cap exceeded")
                    member_name = f"{prefix}{info.filename}"
                    is_nested_archive = info.filename.lower().endswith(".zip")
                    is_tabular = bool(_TABULAR_MEMBER_RE.search(info.filename))
                    if not is_nested_archive and not is_tabular:
                        continue
                    if len(files) >= max_files:
                        continue
                    raw_member = archive.read(info)
                    if is_nested_archive:
                        visit(raw_member, f"{member_name}!/", depth + 1)
                        continue
                    safe_member = _safe_name(member_name)
                    if safe_member in materialized_names:
                        digest = hashlib.sha256(member_name.encode("utf-8")).hexdigest()[:12]
                        collision_index = 1
                        while True:
                            candidate_name = f"{digest}_{collision_index:03d}_{safe_member}"
                            if candidate_name not in materialized_names:
                                safe_member = candidate_name
                                break
                            collision_index += 1
                    local_path = output_dir / safe_member
                    _write_bytes(local_path, raw_member)
                    file_record = _file_entry(
                        record=record,
                        name=member_name,
                        source_url=source_url,
                        local_path=local_path,
                        payload=raw_member,
                        declared_role="",
                    )
                    if file_record is not None:
                        materialized_names.add(safe_member)
                        files.append(file_record)
                    else:
                        local_path.unlink(missing_ok=True)
        except zipfile.BadZipFile as exc:
            if depth == 0:
                raise DiscoveryError(
                    f"candidate archive is not a valid ZIP file: the download is {describe_unexpected_payload(data)}"
                ) from exc
            # A corrupt member deep inside an otherwise readable archive costs
            # that member, not the tables already extracted beside it.
            nested_notes.append(
                f"nested archive {prefix.rstrip('!/')} could not be read: {describe_unexpected_payload(data)}"
            )

    visit(payload, "", 0)
    if not files:
        archive_path.unlink(missing_ok=True)
    return files


def _file_entry(
    *,
    record: dict[str, Any],
    name: str,
    source_url: str,
    local_path: Path,
    payload: bytes,
    declared_role: str,
) -> dict[str, Any] | None:
    classified = classify_filename(name)
    if not classified.get("inspectable"):
        return None
    role = declared_role or str(classified.get("role") or "")
    if role in {"count_matrix", "normalized_expression_matrix", "unknown_matrix"}:
        inspection = inspect_upstream_bytes(name, payload, declared_role=role)
    else:
        inspection = inspect_candidate_bytes(name, payload)
        if inspection.get("status") == "not_deg_table":
            upstream = inspect_upstream_bytes(name, payload, declared_role="unknown_matrix")
            if upstream.get("status") == "upstream_matrix_ready_for_contrast":
                role = "unknown_matrix"
                inspection = upstream
    inspection = dict(inspection)
    inspection.update(
        {
            "local_path": str(local_path),
            "fetch_scope": "full",
            "full_file_sha256": hashlib.sha256(payload).hexdigest(),
        }
    )
    candidate_id = hashlib.sha256(f"{_record_id(record)}|{source_url}|{name}".encode("utf-8")).hexdigest()[:16]
    return {
        "candidate_id": candidate_id,
        "source_url": source_url,
        "name": Path(name).name,
        "role": role or classified["role"],
        "tier": classified["tier"],
        "reason": classified["reason"],
        "inspection": inspection,
    }


def _augment_geo_studies(studies: list[dict[str, Any]], record_by_accession: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    augmented = []
    for study in studies:
        record = record_by_accession.get(str(study.get("accession") or "").upper())
        if not record:
            augmented.append(study)
            continue
        metadata = _publication_metadata_fields(record)
        merged = {**study, **{key: value for key, value in metadata.items() if value}}
        if record.get("source_unit_id"):
            merged["source_unit_id"] = str(record["source_unit_id"])
        augmented.append(merged)
    return augmented


def _publication_study_metadata(record: dict[str, Any], spec: Any) -> dict[str, Any]:
    metadata = _publication_metadata_fields(record)
    return {
        **metadata,
        "species": spec.key,
        "scientific_name": spec.scientific_name,
        "accession": "",
        "study_type": str(record.get("study_type") or record.get("assay_type") or "publication-linked public table"),
        "release_date": str(record.get("publication_date") or record.get("year") or ""),
        "source_url": str(record.get("source_url") or record.get("url") or ""),
    }


def _publication_metadata_fields(record: dict[str, Any]) -> dict[str, Any]:
    pmids = _strings(record.get("pubmed_ids")) or _strings(record.get("pmids")) or _strings(record.get("pmid"))
    doi = _first(record.get("doi"), record.get("dois"), record.get("publication_doi"))
    return {
        "canonical_id": _record_id(record),
        "provider": str(record.get("provider") or record.get("source_provider") or "publication"),
        "provider_accession": _first(record.get("provider_accession"), record.get("canonical_id")),
        "source_unit_id": str(record.get("source_unit_id") or "").strip() or _source_unit_id(record),
        "source_unit_conflict": list(record.get("source_unit_conflict") or []),
        "shared_submission_units": [str(value) for value in (record.get("shared_submission_units") or []) if value],
        "shared_submission_warning": str(record.get("shared_submission_warning") or ""),
        "source_unit_pubmed_ids": pmids,
        "pubmed_ids": pmids,
        "doi": doi,
        "pmcid": _first(record.get("pmcid"), record.get("pmcids")),
        "paper_title": _first(record.get("paper_title"), record.get("title")),
        "title": _first(record.get("paper_title"), record.get("title")),
        "journal": str(record.get("journal") or ""),
        "year": str(record.get("year") or ""),
        "resolution_state": str(record.get("resolution_state") or ""),
        "resolution_events": list(record.get("resolution_events") or []),
    }


def _species_exclusion_reason(record: dict[str, Any], spec: Any) -> str:
    decision = str(
        record.get("species_decision")
        or record.get("mixed_status")
        or record.get("mixed_activation_status")
        or record.get("organism_status")
        or record.get("species_scope_status")
        or ""
    ).strip().lower()
    if record.get("mixed_blocked") is True:
        decision = "mixed_blocked"
    if record.get("mixed_quarantined") is True or decision == "mixed_quarantined":
        return "mixed_quarantined record excluded from species-specific preparation"
    if decision == "mixed_blocked":
        return "mixed_blocked record excluded from species-specific preparation"
    if record.get("mixed_rescued") is True or decision == "mixed_rescued":
        evidence = _first(
            record.get("target_species_evidence"),
            record.get("mixed_rescue_evidence"),
            record.get("evidence_text"),
            record.get("evidence"),
        )
        if record.get("target_species_verified") is not True or not evidence:
            return "mixed_rescued record lacks target_species_verified=true with textual evidence"
    if decision in {"target_species_verified", "target_species_likely", "likely_ready", "verified_ready", ""}:
        return ""
    target = spec.scientific_name.lower()
    evidence = json.dumps(record.get("species_evidence", ""), ensure_ascii=False).lower()
    if evidence and target not in evidence and spec.key not in evidence:
        return "record species evidence does not match requested species"
    return ""


def _candidate_records(record: dict[str, Any]) -> list[dict[str, Any]]:
    raw: list[Any] = []
    for key in (
        "direct_file_candidates",
        "file_candidates",
        "supplementary_file_candidates",
        "public_file_candidates",
        "candidates",
    ):
        raw.extend(_as_list(record.get(key)))
    candidates: list[dict[str, Any]] = []
    seen: set[str] = set()
    for candidate in raw:
        if not isinstance(candidate, dict):
            continue
        url = str(candidate.get("url") or candidate.get("source_url") or "").strip()
        if not url or url in seen:
            continue
        seen.add(url)
        candidates.append(candidate)
    return sorted(candidates, key=_direct_candidate_priority)


def _direct_candidate_priority(candidate: dict[str, Any]) -> tuple[int, int, str, str]:
    """Prefer author DEG tables before generic tables, archives, and matrices."""

    signal = str(
        candidate.get("name")
        or candidate.get("filename")
        or candidate.get("url")
        or candidate.get("source_url")
        or ""
    )
    assessed = classify_filename(signal)
    role = str(candidate.get("role") or assessed.get("role") or "").lower()
    tier = str(candidate.get("tier") or assessed.get("tier") or "").lower()
    suffix = urllib.parse.urlsplit(str(candidate.get("url") or candidate.get("source_url") or signal)).path.lower()
    is_archive = role == "archive" or suffix.endswith(".zip")
    role_priority = {
        "deg_table": 0,
        "unknown_table": 1,
        "archive": 2,
        "count_matrix": 3,
        "normalized_expression_matrix": 4,
        "unknown_matrix": 5,
    }.get("archive" if is_archive else role, 8)
    tier_priority = {"strong": 0, "weak": 1, "archive_candidate": 2, "upstream": 3}.get(tier, 8)
    url = str(candidate.get("url") or candidate.get("source_url") or "")
    return (role_priority, tier_priority, signal.lower(), url)


def _geo_accessions(record: dict[str, Any]) -> list[str]:
    values = _strings(record.get("geo_accessions"))
    values.extend(_strings(record.get("accessions")))
    values.extend(_strings(record.get("accession")))
    return sorted({value.upper() for value in values if _GSE_RE.fullmatch(value)})


def _unique_gse_accessions(records: list[dict[str, Any]]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for record in sorted(records, key=_record_id):
        for accession in _geo_accessions(record):
            if accession not in seen:
                seen.add(accession)
                output.append(accession)
    return output


def _record_by_geo_accession(records: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    mapping: dict[str, dict[str, Any]] = {}
    for record in records:
        for accession in _geo_accessions(record):
            mapping.setdefault(accession, record)
    return mapping


def _record_id(record: dict[str, Any]) -> str:
    explicit = str(record.get("canonical_id") or "").strip()
    return explicit or canonical_record_id(record)


def _source_unit_id(record: dict[str, Any]) -> str:
    for value in (record.get("source_unit_id"), record.get("pmid"), record.get("doi"), record.get("canonical_id")):
        text = _first(value)
        if text:
            return text
    return _record_id(record)


def _preparation_status(candidate_errors: list[dict[str, str]]) -> str:
    if not candidate_errors:
        return "no_usable_table_resolved"
    if any(str(item.get("status") or "") == "rejected" for item in candidate_errors):
        return "source_rejected"
    return "source_unavailable"


def _exclusion_reason(candidate_errors: list[dict[str, str]]) -> str:
    """Say why a study produced nothing, in terms the reader can act on.

    A source that was briefly down is worth re-running; a file the checks
    rejected never will be, and calling both "temporarily unavailable" sends
    people back to repeat a download that cannot succeed.
    """

    if not candidate_errors:
        return "no usable author DEG table or upstream matrix resolved"
    rejected = [item for item in candidate_errors if str(item.get("status") or "") == "rejected"]
    if not rejected:
        return "public candidate source was temporarily unavailable"
    detail = str(rejected[0].get("error") or "").strip()
    name = str(rejected[0].get("candidate") or "").strip()
    prefix = f"public candidate file {name} was rejected" if name else "public candidate file was rejected"
    return f"{prefix}: {detail}" if detail else prefix


def _excluded(record: dict[str, Any], reason: str, **details: Any) -> dict[str, Any]:
    excluded = {
        "canonical_id": _record_id(record),
        "source_unit_id": str(record.get("source_unit_id") or ""),
        "paper_title": _first(record.get("paper_title"), record.get("title")),
        "resolution_state": str(record.get("resolution_state") or ""),
        "resolution_events": list(record.get("resolution_events") or []),
        "reason": reason,
    }
    excluded.update({key: value for key, value in details.items() if value not in (None, "", [], {})})
    return excluded


def _validate_archive_member(info: zipfile.ZipInfo) -> None:
    path = PurePosixPath(info.filename)
    if not info.filename or path.is_absolute() or ".." in path.parts:
        raise DiscoveryUnsafeArchiveError("archive contains an unsafe member path")
    if stat.S_ISLNK(info.external_attr >> 16):
        raise DiscoveryUnsafeArchiveError("archive contains a symbolic link")
    if info.file_size > MAX_ARCHIVE_MEMBER_BYTES:
        raise DiscoveryUnsafeArchiveError("archive contains an oversized member")


def _write_marker(path: Path, species_key: str) -> None:
    payload = {
        "artifact_type": DISCOVERY_BUNDLE_ARTIFACT_TYPE,
        "format_version": DISCOVERY_BUNDLE_FORMAT_VERSION,
        "species": species_key,
    }
    (path / DISCOVERY_BUNDLE_MARKER).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _publish_prepared_bundle(staging: Path, target: Path, *, force: bool) -> None:
    from .discovery import _publish_prepared_bundle as publish

    publish(staging, target, force=force)


def _validate_target(target: Path, *, force: bool) -> None:
    from .discovery import _validate_preparation_target as validate

    validate(target, force=force)


def _retarget_paths(result: dict[str, Any], staging: Path, target: Path) -> None:
    staging_root = staging.resolve()
    for study in result.get("studies", []):
        for candidate in study.get("files", []):
            inspection = candidate.get("inspection") or {}
            local = inspection.get("local_path")
            if not local:
                continue
            relative = Path(str(local)).resolve().relative_to(staging_root)
            inspection["local_path"] = str(target / relative)


def _export_paths(target: Path) -> dict[str, str]:
    return {
        "output_dir": str(target),
        "audit_json": str(target / "discovery_audit.json"),
        "candidates_csv": str(target / "discovery_candidates.csv"),
        "draft_catalog_csv": str(target / "DEGORA_discovery_draft_catalog.csv"),
    }


def _write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _candidate_name(candidate: dict[str, Any]) -> str:
    name = str(candidate.get("name") or "")
    if not name:
        name = Path(urllib.parse.urlsplit(str(candidate.get("url") or candidate.get("source_url") or "")).path).name
    return _safe_name(name or "candidate.dat")


def _safe_name(value: str) -> str:
    # ':' is deliberately excluded: on Windows "output_dir" / "C:member.csv" resolves to
    # the drive-relative "C:member.csv", so an archive member could name a file outside
    # the bundle even though the member path itself passed the traversal checks.
    safe = re.sub(r"[^A-Za-z0-9_.!-]+", "_", str(value)).strip("._")
    return safe[:180] or hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:16]


def _strings(value: Any) -> list[str]:
    return [str(item).strip() for item in _as_list(value) if str(item).strip()]


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def _first(*values: Any) -> str:
    for value in values:
        for item in _as_list(value):
            text = str(item or "").strip()
            if text:
                return text
    return ""
