"""Canonical public identifiers used to keep source units independent."""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Any


_MISSING = {"", "none", "null", "nan", "na", "n/a"}
_LIST_SEPARATOR_RE = re.compile(r"[,;|]+")


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in _MISSING else text


def _identifier_chunks(values: Iterable[Any]) -> Iterable[str]:
    for value in values:
        if isinstance(value, (list, tuple, set, frozenset)):
            yield from _identifier_chunks(value)
            continue
        text = _clean_text(value)
        if not text:
            continue
        yield from (chunk.strip() for chunk in _LIST_SEPARATOR_RE.split(text) if chunk.strip())


def normalize_pmids(values: Any) -> list[str]:
    """Return every PMID in scalar, collection, or comma/semicolon/pipe input."""

    normalized = {
        match.group(1)
        for chunk in _identifier_chunks([values])
        if (match := re.fullmatch(r"(?:PMID[\s:_-]*)?(\d+)", chunk, re.IGNORECASE))
    }
    return sorted(normalized)


def normalize_pmcids(values: Any) -> list[str]:
    """Return every canonical ``PMC123`` identifier from common PMC prefixes."""

    normalized = {
        f"PMC{match.group(1)}"
        for chunk in _identifier_chunks([values])
        if (
            match := re.fullmatch(
                r"(?:(?:PMCID)[\s:_-]*)?(?:PMC[\s:_-]*)?(\d+)",
                chunk,
                re.IGNORECASE,
            )
        )
    }
    return sorted(normalized)


def normalize_study_accession(value: Any) -> str:
    """Canonicalize GEO/ArrayExpress aliases without rewriting opaque IDs."""

    text = _clean_text(value).upper()
    if not text:
        return ""
    text = re.sub(r"^ACCESSION[\s:_-]*", "", text)
    compact = re.sub(r"[\s:_-]+", "", text)
    if match := re.fullmatch(r"(?:GSE|EGEOD)(\d+)", compact):
        return f"GSE{match.group(1)}"
    if match := re.fullmatch(r"EMTAB(\d+)", compact):
        return f"E-MTAB-{match.group(1)}"
    return re.sub(r"\s+", "", text)


def recognized_source_unit_key(value: Any) -> str:
    """Return a comparison key only for a recognized public identifier."""

    text = _clean_text(value)
    if not text:
        return ""
    if match := re.fullmatch(r"PMID[\s:_-]*(\d+)", text, re.IGNORECASE):
        return f"pmid:{match.group(1)}"
    if match := re.fullmatch(
        r"(?:PMCID[\s:_-]*(?:PMC[\s:_-]*)?|PMC[\s:_-]*)(\d+)",
        text,
        re.IGNORECASE,
    ):
        return f"pmcid:PMC{match.group(1)}"
    accession = normalize_study_accession(text)
    if re.fullmatch(r"GSE\d+|E-MTAB-\d+", accession):
        return f"accession:{accession}"
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", text, flags=re.IGNORECASE)
    doi = re.sub(r"^DOI\s*:\s*", "", doi, flags=re.IGNORECASE).strip().rstrip(" .")
    if re.match(r"^10\.\d{4,9}/\S+$", doi, re.IGNORECASE):
        return f"doi:{doi.lower()}"
    return ""


def canonical_source_unit_id(value: Any) -> str:
    """Return a stable display ID for recognized aliases, else preserve the text.

    An all-digit value is treated as a PMID only here, where the caller has
    already established that the value is an explicit source-unit identifier.
    Bare numbers in generic accession fields remain ambiguous and are not
    rewritten by :func:`recognized_source_unit_key`.
    """

    text = _clean_text(value)
    if re.fullmatch(r"\d+", text):
        return f"PMID:{text}"
    key = recognized_source_unit_key(text)
    if key.startswith("pmid:"):
        return f"PMID:{key.split(':', 1)[1]}"
    if key.startswith("pmcid:"):
        return f"PMCID:{key.split(':', 1)[1]}"
    if key.startswith("accession:"):
        return key.split(":", 1)[1]
    if key.startswith("doi:"):
        return f"DOI:{key.split(':', 1)[1]}"
    return text


__all__ = [
    "canonical_source_unit_id",
    "normalize_pmcids",
    "normalize_pmids",
    "normalize_study_accession",
    "recognized_source_unit_key",
]
