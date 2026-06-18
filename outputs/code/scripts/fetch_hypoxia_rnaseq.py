#!/usr/bin/env python
"""Fetch the hypoxia RNA-seq author DEG tables from their recorded sources.

The hypoxia corpus mixes as-published author DEG tables from GEO, Springer, MDPI,
Mendeley, and a Bioconductor package. This fetcher resolves and downloads the ones
that have a constructable direct file URL (every GEO source: a ``ftp.ncbi`` suppl
path built from the accession plus the stored filename, or an already-direct suppl
URL) plus the Springer ESM workbook, into the catalog's source_path. GEO sources
that need a deterministic transform (sign flip / RefSeq->symbol / sheet extraction)
are downloaded to a ``.raw`` staging file for a separate transform step. Sources
behind a landing page with no constructable file URL (MDPI, Mendeley, Bioconductor)
are reported, not fetched.

    PYTHONPATH=outputs/code python outputs/code/scripts/fetch_hypoxia_rnaseq.py
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

RAW = Path("data/deg/raw")
DOWNLOAD_TIMEOUT_SECONDS = 180
MAX_BYTES = 500 * 1024 * 1024


def _geo_suppl(accession: str, filename: str) -> str:
    # GSE76743 -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76743/suppl/<filename>
    prefix = f"GSE{accession[3:-3]}nnn" if len(accession) > 6 else "GSEnnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/{filename}"


# (study_id, resolved_url, target_path_relative_to_RAW, kind)
#   passthrough  -> downloaded file IS the harmonized source_path (no transform)
#   transform    -> downloaded raw needs a deterministic transform to become source_path
MANIFEST = [
    ("HYP001", _geo_suppl("GSE76743", "GSE76743_HUVEC_HYPOXIA_mRNA_DESeq2.txt.gz"), "GSE76743_HUVEC_HYPOXIA_mRNA_DESeq2.txt.gz", "passthrough"),
    ("HYP009", _geo_suppl("GSE132624", "GSE132624_edgeR_501mel_24h.csv.gz"), "iter6/GSE132624_edgeR_501mel_24h.csv.gz", "passthrough"),
    ("HYP010", _geo_suppl("GSE132624", "GSE132624_edgeR_IGR37_24h.csv.gz"), "iter6/GSE132624_edgeR_IGR37_24h.csv.gz", "passthrough"),
    ("HYP011", _geo_suppl("GSE132624", "GSE132624_edgeR_IGR39_24h.csv.gz"), "iter6/GSE132624_edgeR_IGR39_24h.csv.gz", "passthrough"),
    ("HYP012", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70544/suppl/GSE70544_Normoxia_vs_Hypoxia_Gencode_gene_exp.diff.txt.gz", "iter7/GSE70544_Normoxia_vs_Hypoxia_Gencode_gene_exp.diff.txt.gz", "passthrough"),
    ("HYP014", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311800/suppl/GSE311800_DEG_Results_All_Comparisons.xlsx", "iter8/GSE311800_DEG_Results_All_Comparisons.xlsx", "passthrough"),
    ("HYP015", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293238/suppl/GSE293238_RNAseq_processed.xlsx", "iter9/GSE293238_RNAseq_processed.xlsx", "passthrough"),
    ("HYP018", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313740/suppl/GSE313740_processed_dataxGEO_14dic25.xlsx", "iter12/GSE313740_processed_dataxGEO_14dic25.xlsx", "passthrough"),
    ("HYP019", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225253/suppl/GSE225253_DESeq2_HK-2_hypoxia_vs_HK-2_normoxia_all.xlsx", "iter13/GSE225253_DESeq2_HK-2_hypoxia_vs_HK-2_normoxia_all.xlsx", "passthrough"),
    ("HYP020", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283446/suppl/GSE283446_TableS8.xlsx", "iter15/GSE283446_TableS8.xlsx", "passthrough"),
    # transform sources: download the raw GEO/Springer file to a .raw staging path
    ("HYP013", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108676/suppl/GSE108676_Human_gene_exp.diff.gz", "iter7/GSE108676_Human_gene_exp.diff.gz", "transform:refseq_map"),
    ("HYP017", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255352/suppl/GSE255352_OvsH_deg.xls.gz", "iter11/GSE255352_OvsH_deg.xls.gz", "transform:sign_flip"),
    ("HYP021", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160491/suppl/GSE160491_Processed_all_compare.xls.gz", "iter15/GSE160491_Processed_all_compare.xls.gz", "transform:sign_flip"),
    ("HYP023", "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-026-09875-6/MediaObjects/42003_2026_9875_MOESM2_ESM.xlsx", "iter16/GSE296192_CommBio_ESM2.xlsx", "transform:xlsx_extract"),
    ("HYP024", "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-026-09875-6/MediaObjects/42003_2026_9875_MOESM2_ESM.xlsx", "iter16/GSE296192_CommBio_ESM2.xlsx", "transform:xlsx_extract"),
]
# landing-page sources with no constructable direct file URL (reported, not fetched)
LANDING = {
    "HYP003": "https://bioconductor.org/packages/release/bioc/vignettes/TFEA.ChIP/ (TFEA.ChIP package hypoxia example)",
    "HYP006": "https://www.mdpi.com/1422-0067/23/10/5824 (MDPI supplementary)",
    "HYP007": "https://www.mdpi.com/1422-0067/23/10/5824 (MDPI supplementary)",
    "HYP008": "https://data.mendeley.com/datasets/z42wpkbb8k/1 (Mendeley dataset)",
}


def _download(url: str, target: Path) -> tuple[bool, str]:
    tmp = target.with_name(target.name + ".part")
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(url, headers={"User-Agent": "DEGORA-fetch/1.0"})
        total = 0
        with urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response, open(tmp, "wb") as handle:
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_BYTES:
                    tmp.unlink(missing_ok=True)
                    return False, "exceeds cap"
                handle.write(chunk)
        if total == 0:
            tmp.unlink(missing_ok=True)
            return False, "empty"
        tmp.replace(target)
        return True, f"{total:,} bytes"
    except Exception as exc:  # noqa: BLE001
        tmp.unlink(missing_ok=True)
        return False, f"{type(exc).__name__}: {exc}"


def main() -> int:
    counts = {"downloaded": 0, "present": 0, "failed": 0}
    seen: set[str] = set()
    for study, url, rel, kind in MANIFEST:
        # passthrough rel == the catalog source_path; transform rel == a distinct raw
        # staging filename (so it never overwrites the transformed source_path).
        target = RAW / rel
        if str(target) in seen:
            counts["present"] += 1
            print(f"  shared    {study:<8} {target.name} (already fetched)")
            continue
        seen.add(str(target))
        if target.exists():
            counts["present"] += 1
            print(f"  present   {study:<8} {target.name}")
            continue
        ok, note = _download(url, target)
        counts["downloaded" if ok else "failed"] += 1
        print(f"  {'OK  ' if ok else 'FAIL'}      {study:<8} [{kind}] {target.name}: {note}")
    print("\nlanding-page sources (need manual/web resolution, not fetched):")
    for study, where in LANDING.items():
        print(f"  -- {study}: {where}")
    print(f"\nhypoxia raw fetch: downloaded={counts['downloaded']} present={counts['present']} failed={counts['failed']} | landing={len(LANDING)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
