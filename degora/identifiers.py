"""What space a gene column is written in, judged from its values.

DEGORA matches genes across studies on the identifier itself, so an Ensembl
column and a symbol column have nothing in common even when they describe the
same genes. Every path that reads a gene column - the guided setup, the
preparation inspector, the run - asks this one function, so a column named
"Gene Symbol" that holds Ensembl IDs is called what it is everywhere.
"""

from __future__ import annotations

import re
from typing import Any, Iterable

IDENTIFIER_PATTERNS = (
    ("Ensembl ID", re.compile(r"^ENS[A-Z]*[GT]\d{6,}(\.\d+)?$", re.IGNORECASE)),
    ("RefSeq ID", re.compile(r"^[NX][MRP]_\d+(\.\d+)?$", re.IGNORECASE)),
    ("Affymetrix probe ID", re.compile(r"^\d+_[a-z]?_?at$", re.IGNORECASE)),
    ("Entrez ID", re.compile(r"^\d+$")),
    ("gene symbol", re.compile(r"^[A-Za-z][A-Za-z0-9\-.@_]{0,24}$")),
)
IDENTIFIER_SAMPLE_ROWS = 200
UNKNOWN_IDENTIFIER_SPACE = "unrecognised identifiers"
PROBABILITY_COLUMN_NAME_RE = re.compile(
    r"(^|[^A-Za-z0-9])(adj[._ -]?)?p([._ -]?(value|val))($|[^A-Za-z0-9])|"
    r"^p$|^p[._ -](?!vs($|[._ -])).+|"
    r"(^|[^A-Za-z0-9])p[._ -]?adj(ust(ed)?)?($|[^A-Za-z0-9])|"
    r"q[._ -]?value|fdr|false[._ -]?discovery",
    re.IGNORECASE,
)
BARE_IDENTIFIER_RE = re.compile(r"^(?:id|identifier)$", re.IGNORECASE)
SAMPLE_IDENTIFIER_NAME_RE = re.compile(
    r"(?:^|[_ .-])(?:sample|subject|patient|donor|run|file|series)(?:[_ .-]?(?:id|name|accession|identifier))?$|"
    r"^(?:gsm|srr|err|drr)(?:[_ .-]?(?:id|accession|identifier))?$",
    re.IGNORECASE,
)
SAMPLE_ACCESSION_RE = re.compile(r"^(?:GSM|SRR|ERR|DRR)\d+$", re.IGNORECASE)


def identifier_space_profile(values: Iterable[Any]) -> dict[str, float]:
    """The share of each recognised space among the first IDENTIFIER_SAMPLE_ROWS labels.

    identifier_space() names the majority; a column that is 70% symbols and 30%
    Ensembl IDs is "gene symbol" to it, and the 30% then join nothing in a
    symbol corpus. The profile keeps the minority visible.
    """

    counts: dict[str, int] = {}
    seen = 0
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        seen += 1
        for label, pattern in IDENTIFIER_PATTERNS:
            if pattern.match(text):
                counts[label] = counts.get(label, 0) + 1
                break
        if seen >= IDENTIFIER_SAMPLE_ROWS:
            break
    return {label: hits / seen for label, hits in counts.items()} if seen else {}


def identifier_space(values: Iterable[Any]) -> str:
    """Name the identifier space a gene column is written in.

    DEGORA matches genes across studies by the identifier itself, so a table of
    Ensembl IDs and a table of symbols have nothing in common even when they
    describe the same genes. Recognising that while the config is being built is
    the difference between a sentence now and a run that scores zero genes later.
    """

    counts: dict[str, int] = {}
    seen = 0
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        seen += 1
        for label, pattern in IDENTIFIER_PATTERNS:
            if pattern.match(text):
                counts[label] = counts.get(label, 0) + 1
                break
        if seen >= IDENTIFIER_SAMPLE_ROWS:
            break
    if not counts or not seen:
        return UNKNOWN_IDENTIFIER_SPACE
    label, hits = max(counts.items(), key=lambda item: item[1])
    # A clear majority, or the column is too mixed to name honestly.
    return label if hits >= seen * 0.7 else UNKNOWN_IDENTIFIER_SPACE



# DEGORA matches genes across studies on the identifier itself, and most published
# tables carry symbols, so a symbol column is the default that joins the most.
GENE_SPACE_PREFERENCE = ("gene symbol", "Ensembl ID", "RefSeq ID", "Entrez ID", "Affymetrix probe ID")
