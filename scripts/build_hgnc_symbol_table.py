#!/usr/bin/env python3
"""Regenerate ``degora/data/hgnc_previous_symbols.tsv`` from the HGNC complete set.

DEGORA joins source tables on the gene symbol, so two papers that spell one gene
differently do not meet. The Excel-damage rescue in ``harmonize.py`` covers the
SEPT/MARCH/DEC families; this table covers the ordinary HGNC retirements that
separate an older paper from a newer one (``CTGF`` and ``CCN2``, ``IL8`` and
``CXCL8``, ``KIAA0101`` and ``PCLAF``).

Only *unambiguous* retirements are kept, because a wrong merge is worse than a
missed one:

1. the record is ``status == Approved``;
2. the previous symbol is not itself the approved symbol of any gene - HGNC
   reuses symbols, so ``BRF1`` is both an approved gene and a previous symbol of
   ``ZFP36L1``, and mapping it would destroy the approved gene;
3. the previous symbol names exactly one approved symbol - ``DEC1`` is a
   previous symbol of both ``BHLHE40`` and ``DELEC1`` and is therefore dropped
   here, leaving the documented ``DEC1 -> BHLHE40`` choice in ``harmonize.py``
   as the single rule;
4. the previous symbol is a plain symbol token, never an accession or a bare
   number, so nothing that is really an Ensembl/RefSeq/Entrez identifier is
   rewritten;
5. the previous symbol is not one of the Excel date-damage families, which the
   dedicated repair already resolves before this table is consulted;
6. the previous symbol is not a token that reads as something else in a table -
   HGNC retires ``P`` to ``OCA2`` and ``STAT`` to ``SOAT1``, but such a cell is
   far more often a header fragment or a placeholder than the gene.

Usage:

    python scripts/build_hgnc_symbol_table.py [hgnc_complete_set.txt]

Without an argument the current HGNC complete set is downloaded. The header of
the written file records the source URL, retrieval date and SHA-256 of the exact
input snapshot. HGNC's download URL is mutable, so the retrieval date is not
described as an HGNC release.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
OUTPUT = Path(__file__).resolve().parents[1] / "degora" / "data" / "hgnc_previous_symbols.tsv"

SYMBOL_TOKEN_RE = re.compile(r"^[A-Z][A-Z0-9._@-]{0,29}$")
ACCESSION_LIKE_RE = re.compile(r"^(?:ENS[A-Z]*[GTPR]\d+|[NX][MRPCGWZ]_\d+)$", re.IGNORECASE)
EXCEL_DATE_FAMILY_RE = re.compile(r"^(?:SEPT?|MARCH?|DEC)\d+$")
# Historical symbols whose non-gene meaning dominates in a spreadsheet. HGNC really
# does retire "P" (to OCA2) and "STAT" (to SOAT1), but a cell reading P or STAT in a
# results table is far more often a header fragment, a placeholder or a statistic
# than the gene, and rewriting it would invent evidence. The value of keeping them
# is close to zero; the cost of a wrong merge is not.
NON_GENE_TOKENS = frozenset(
    {
        "ALL", "FALSE", "FC", "GENE", "ID", "INF", "LFC", "MAX", "MEAN", "MIN", "N/A", "NA",
        "NAME", "NAN", "NO", "NONE", "NULL", "P", "RANK", "REF", "STAT", "SUM", "TEST",
        "TOTAL", "TRUE", "YES",
    }
)


def read_rows(source: str | None) -> tuple[list[dict[str, str]], str]:
    if source:
        raw = Path(source).read_bytes()
    else:
        with urllib.request.urlopen(HGNC_URL, timeout=300) as response:
            raw = response.read()
    text = raw.decode("utf-8", "replace")
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    return list(reader), hashlib.sha256(raw).hexdigest()


def build(rows: list[dict[str, str]]) -> tuple[list[tuple[str, str]], dict[str, int]]:
    approved = [row for row in rows if (row.get("status") or "").strip() == "Approved"]
    current = {(row.get("symbol") or "").strip().upper() for row in approved}
    current.discard("")

    targets: dict[str, set[str]] = defaultdict(set)
    for row in approved:
        symbol = (row.get("symbol") or "").strip().upper()
        if not symbol:
            continue
        for previous in (row.get("prev_symbol") or "").split("|"):
            previous = previous.strip().upper()
            if previous:
                targets[previous].add(symbol)

    counts = {
        "approved_records": len(approved),
        "previous_symbols_seen": len(targets),
        "dropped_symbol_reused_as_current": 0,
        "dropped_ambiguous_target": 0,
        "dropped_not_a_symbol_token": 0,
        "dropped_excel_date_family": 0,
        "dropped_non_gene_token": 0,
    }
    pairs: list[tuple[str, str]] = []
    for previous, symbols in targets.items():
        if previous in current:
            counts["dropped_symbol_reused_as_current"] += 1
            continue
        if len(symbols) != 1:
            counts["dropped_ambiguous_target"] += 1
            continue
        if EXCEL_DATE_FAMILY_RE.match(previous):
            counts["dropped_excel_date_family"] += 1
            continue
        if previous in NON_GENE_TOKENS:
            counts["dropped_non_gene_token"] += 1
            continue
        if not SYMBOL_TOKEN_RE.match(previous) or ACCESSION_LIKE_RE.match(previous):
            counts["dropped_not_a_symbol_token"] += 1
            continue
        pairs.append((previous, next(iter(symbols))))
    pairs.sort()
    counts["kept"] = len(pairs)
    return pairs, counts


def main() -> int:
    rows, source_sha256 = read_rows(sys.argv[1] if len(sys.argv) > 1 else None)
    pairs, counts = build(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# DEGORA gene-symbol retirement table\n")
        handle.write(f"# source\t{HGNC_URL}\n")
        handle.write(f"# snapshot_date\t{dt.date.today().isoformat()}\n")
        handle.write(f"# source_sha256\t{source_sha256}\n")
        handle.write(
            "# rule\tunambiguous HGNC prev_symbol -> approved symbol; a previous symbol that is "
            "itself an approved symbol, names more than one approved symbol, is accession- or "
            "number-shaped, belongs to an Excel date-damage family, or is a token whose non-gene "
            "meaning dominates in a spreadsheet is not listed\n"
        )
        handle.write(
            "# counts\t" + "; ".join(f"{key}={value}" for key, value in counts.items()) + "\n"
        )
        handle.write("previous_symbol\tcurrent_symbol\n")
        for previous, symbol in pairs:
            handle.write(f"{previous}\t{symbol}\n")
    for key, value in counts.items():
        print(f"{key}: {value:,}")
    print(f"wrote {OUTPUT} ({OUTPUT.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
