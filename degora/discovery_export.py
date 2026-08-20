"""Audit-grade exports for source-neutral publication discovery snapshots."""

from __future__ import annotations

import csv
import io
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


FORMULA_PREFIXES = ("=", "+", "-", "@")
SEARCH_JSON_NAME = "publication_search.json"
SEARCH_CSV_NAME = "publication_search.csv"
SEARCH_XLSX_NAME = "DEGORA_search_results.xlsx"


def _safe_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, (list, tuple, set)):
        value = "; ".join(str(item) for item in value)
    elif isinstance(value, dict):
        value = json.dumps(value, sort_keys=True, ensure_ascii=False)
    else:
        value = str(value)
    return "'" + value if value.startswith(FORMULA_PREFIXES) else value


def _atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_text(path: Path, text: str) -> None:
    _atomic_bytes(path, text.encode("utf-8"))


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if value in (None, ""):
        return []
    return [value]


def _readiness(record: dict[str, Any]) -> dict[str, Any]:
    value = record.get("data_readiness") or record.get("deg_input_assessment") or {}
    return value if isinstance(value, dict) else {}


def publication_rows(snapshot: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in snapshot.get("records", snapshot.get("studies", [])):
        readiness = _readiness(record)
        rows.append(
            {
                "paper_title": record.get("paper_title") or record.get("title", ""),
                "authors": record.get("authors_display") or record.get("authors", []),
                "journal": record.get("journal", ""),
                "year": record.get("year", ""),
                "canonical_id": record.get("canonical_id") or record.get("source_unit_id") or record.get("accession", ""),
                "record_kind": record.get("record_kind", "publication"),
                "species": record.get("species", ""),
                "species_decision": record.get("species_decision", ""),
                "pubmed_ids": record.get("pubmed_ids", []),
                "doi": record.get("doi", ""),
                "pmcid": record.get("pmcid", ""),
                "geo_accessions": record.get("geo_accessions") or _as_list(record.get("accession")),
                "source_unit_id": record.get("source_unit_id", ""),
                "source_unit_conflict": record.get("source_unit_conflict", []),
                "resolution_state": record.get("resolution_state", ""),
                "relevance_rank": record.get("relevance_rank", record.get("ncbi_relevance_rank", "")),
                "readiness_tier": readiness.get("tier", ""),
                "readiness_priority": readiness.get("priority", ""),
                "verification_state": readiness.get("verification_state", "likely"),
                "readiness_basis": readiness.get("basis", ""),
                "candidate_count": len(_as_list(record.get("candidates"))),
            }
        )
    return rows


def _identifier_rows(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        canonical = record.get("canonical_id") or record.get("source_unit_id") or record.get("accession", "")
        values: list[tuple[str, str]] = []
        values.extend(("PMID", str(value)) for value in _as_list(record.get("pubmed_ids")))
        values.extend(("DOI", str(value)) for value in _as_list(record.get("doi")))
        values.extend(("PMCID", str(value)) for value in _as_list(record.get("pmcid")))
        values.extend(("GEO", str(value)) for value in _as_list(record.get("geo_accessions") or record.get("accession")))
        provider_ids = record.get("provider_ids") or {}
        if isinstance(provider_ids, dict):
            for provider, provider_values in provider_ids.items():
                values.extend((str(provider), str(value)) for value in _as_list(provider_values))
        elif isinstance(provider_ids, list):
            for value in provider_ids:
                text = str(value).strip()
                provider, separator, identifier = text.partition(":")
                values.append((provider if separator else "provider", identifier if separator else text))
        seen: set[tuple[str, str]] = set()
        for kind, value in values:
            key = (kind, value.strip())
            if not key[1] or key in seen:
                continue
            seen.add(key)
            rows.append({"canonical_id": canonical, "identifier_type": key[0], "identifier": key[1]})
    return rows


def _dataset_rows(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        canonical = record.get("canonical_id") or record.get("source_unit_id") or record.get("accession", "")
        for accession in _as_list(record.get("geo_accessions") or record.get("accession")):
            if str(accession).strip():
                rows.append({"canonical_id": canonical, "provider": "NCBI GEO", "accession": accession})
        for source in _as_list(record.get("sources")):
            if isinstance(source, dict):
                rows.append(
                    {
                        "canonical_id": canonical,
                        "provider": source.get("provider", ""),
                        "accession": source.get("accession", source.get("id", "")),
                        "landing_url": source.get("landing_url", source.get("source_url", "")),
                    }
                )
    return rows


def _candidate_rows(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        canonical = record.get("canonical_id") or record.get("source_unit_id") or record.get("accession", "")
        for candidate in _as_list(record.get("candidates")):
            if not isinstance(candidate, dict):
                continue
            rows.append(
                {
                    "canonical_id": canonical,
                    "candidate_id": candidate.get("candidate_id", ""),
                    "provider": candidate.get("provider", ""),
                    "name": candidate.get("name", ""),
                    "role": candidate.get("role", ""),
                    "source_input_type": candidate.get("source_input_type", ""),
                    "verification_state": candidate.get("verification_state", "likely"),
                    "source_url": candidate.get("source_url", ""),
                    "landing_url": candidate.get("landing_url", ""),
                    "license": candidate.get("license", ""),
                    "reason": candidate.get("reason", ""),
                }
            )
    return rows


def _species_rows(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        rows.append(
            {
                "canonical_id": record.get("canonical_id") or record.get("source_unit_id") or record.get("accession", ""),
                "requested_species": record.get("species", ""),
                "species_decision": record.get("species_decision", ""),
                "target_species_verified": record.get("target_species_verified", False),
                "species_evidence": record.get("species_evidence", ""),
            }
        )
    return rows


def _event_rows(snapshot: dict[str, Any]) -> list[dict[str, Any]]:
    diagnostics = snapshot.get("provider_events") or snapshot.get("provider_diagnostics") or []
    if isinstance(diagnostics, dict):
        diagnostics = [dict(provider=key, detail=value) for key, value in diagnostics.items()]
    return [value if isinstance(value, dict) else {"detail": value} for value in _as_list(diagnostics)]


def _write_sheet(workbook: Workbook, title: str, rows: list[dict[str, Any]], columns: list[str]) -> None:
    worksheet = workbook.create_sheet(title)
    worksheet.append(columns)
    for row in rows:
        worksheet.append([_safe_cell(row.get(column, "")) for column in columns])
    header_fill = PatternFill("solid", fgColor="18364A")
    for cell in worksheet[1]:
        cell.fill = header_fill
        cell.font = Font(color="FFFFFF", bold=True)
        cell.alignment = Alignment(vertical="center")
    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = worksheet.dimensions
    for index, column in enumerate(columns, start=1):
        values = [str(column)] + [str(_safe_cell(row.get(column, ""))) for row in rows[:200]]
        worksheet.column_dimensions[get_column_letter(index)].width = min(max(max(map(len, values)) + 2, 11), 48)


def build_publication_search_workbook(snapshot: dict[str, Any]) -> bytes:
    records = list(snapshot.get("records", snapshot.get("studies", [])))
    publications = publication_rows(snapshot)
    species_value = snapshot.get("species", "")
    if isinstance(species_value, dict):
        species_value = species_value.get("key") or species_value.get("label") or ""
    query_rows = [
        {"field": "query", "value": snapshot.get("query", "")},
        {"field": "species", "value": species_value},
        {"field": "generated_at", "value": snapshot.get("generated_at", "")},
        {"field": "evaluated_records", "value": snapshot.get("evaluated_records", len(records))},
        {"field": "ranking_limit", "value": snapshot.get("ranking_limit", "")},
        {"field": "ranking_truncated", "value": snapshot.get("ranking_truncated", False)},
        {"field": "ranking_contract", "value": snapshot.get("ranking_contract", "")},
        {"field": "cross_species_pooling", "value": False},
    ]
    workbook = Workbook()
    workbook.remove(workbook.active)
    _write_sheet(workbook, "Query", query_rows, ["field", "value"])
    _write_sheet(workbook, "Publications", publications, list(publications[0]) if publications else [
        "canonical_id", "record_kind", "species", "species_decision", "paper_title", "authors", "journal",
        "year", "pubmed_ids", "doi", "pmcid", "geo_accessions", "source_unit_id", "source_unit_conflict",
        "resolution_state", "relevance_rank",
        "readiness_tier", "readiness_priority", "verification_state", "readiness_basis", "candidate_count",
    ])
    identifiers = _identifier_rows(records)
    _write_sheet(workbook, "Identifiers", identifiers, list(identifiers[0]) if identifiers else ["canonical_id", "identifier_type", "identifier"])
    datasets = _dataset_rows(records)
    _write_sheet(workbook, "Linked datasets", datasets, list(datasets[0]) if datasets else ["canonical_id", "provider", "accession", "landing_url"])
    candidates = _candidate_rows(records)
    _write_sheet(workbook, "Candidate routes", candidates, list(candidates[0]) if candidates else [
        "canonical_id", "candidate_id", "provider", "name", "role", "source_input_type",
        "verification_state", "source_url", "landing_url", "license", "reason",
    ])
    species = _species_rows(records)
    _write_sheet(workbook, "Species decisions", species, list(species[0]) if species else [
        "canonical_id", "requested_species", "species_decision", "target_species_verified", "species_evidence",
    ])
    events = _event_rows(snapshot)
    event_keys = {key for row in events for key in row}
    event_columns = [key for key in ("provider", "event", "status", "message") if key in event_keys]
    event_columns.extend(sorted(event_keys.difference(event_columns)))
    if not event_columns:
        event_columns = ["provider", "event", "status", "message"]
    _write_sheet(workbook, "Provider events", events, event_columns)
    buffer = io.BytesIO()
    workbook.save(buffer)
    return buffer.getvalue()


def _publication_csv(rows: list[dict[str, Any]]) -> str:
    columns = list(rows[0]) if rows else ["canonical_id", "paper_title", "species", "species_decision"]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({column: _safe_cell(row.get(column, "")) for column in columns})
    return buffer.getvalue()


def export_publication_search(snapshot: dict[str, Any], output_dir: str | Path, *, force: bool = False) -> dict[str, str]:
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / SEARCH_JSON_NAME
    csv_path = output / SEARCH_CSV_NAME
    xlsx_path = output / SEARCH_XLSX_NAME
    existing = [path for path in (json_path, csv_path, xlsx_path) if path.exists()]
    if existing and not force:
        raise FileExistsError("federated search output already exists; use --force to replace: " + ", ".join(map(str, existing)))
    rows = publication_rows(snapshot)
    _atomic_text(json_path, json.dumps(snapshot, indent=2, sort_keys=True, ensure_ascii=False) + "\n")
    _atomic_text(csv_path, _publication_csv(rows))
    _atomic_bytes(xlsx_path, build_publication_search_workbook(snapshot))
    return {
        "output_dir": str(output),
        "search_json": str(json_path),
        "search_csv": str(csv_path),
        "search_xlsx": str(xlsx_path),
    }


__all__ = [
    "SEARCH_CSV_NAME",
    "SEARCH_JSON_NAME",
    "SEARCH_XLSX_NAME",
    "build_publication_search_workbook",
    "export_publication_search",
    "publication_rows",
]
