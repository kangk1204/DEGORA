#!/usr/bin/env python
"""Reconstruct the hypoxia RNA-seq author DEG tables from their recorded sources.

The hypoxia corpus is the messiest by design: it mixes as-published DEG tables
from GEO (direct suppl files, some renamed), a Springer Communications Biology
supplement, an MDPI supplement, a Mendeley dataset, and a Bioconductor package.
This script downloads each source from a constructable/recorded direct URL and
applies the exact deterministic transform needed to land the catalog's
``source_path`` file, so the hypoxia corpus reproduces from scratch with no manual
download. It writes ONLY the per-source DEG tables under data/deg/raw/; it never
touches the locked catalog or gold panel.

After this, `make pipeline TOPIC=hypoxia` re-harmonizes the 19 sources from raw.
Needs network; HYP003 (TFEA.ChIP) additionally needs Rscript to read the .rda.

    PYTHONPATH=outputs/code python outputs/code/scripts/build_hypoxia_sources.py
    make -C outputs/code fetch-hypoxia
"""

from __future__ import annotations

import gzip
import hashlib
import io
import subprocess
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

from degora.harmonize import _repair_excel_date_gene_symbol
from degora.provenance import shell_command, write_source_sidecar

RAW = Path("data/deg/raw")
HGNC = RAW / "ifn/hgnc_complete_set.txt"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
TIMEOUT = 180
COMMAND = "PYTHONPATH=outputs/code python3 outputs/code/scripts/build_hypoxia_sources.py"
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
MAX_DOWNLOAD_BYTES = 200 * 1024 * 1024
MAX_ZIP_MEMBER_BYTES = 100 * 1024 * 1024
EXPECTED_DOWNLOAD_SHA256: dict[str, str] = {}


def _geo_suppl(accession: str, filename: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{accession[3:-3]}nnn/{accession}/suppl/{filename}"


def _get(url: str, *, max_bytes: int = MAX_DOWNLOAD_BYTES, expected_sha256: str | None = None) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "DEGORA-fetch/1.0"})
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        buffer = io.BytesIO()
        hasher = hashlib.sha256()
        total = 0
        while True:
            chunk = response.read(DOWNLOAD_CHUNK_SIZE)
            if not chunk:
                break
            total += len(chunk)
            if total > max_bytes:
                raise RuntimeError(f"download exceeded {max_bytes} bytes: {url}")
            hasher.update(chunk)
            buffer.write(chunk)
    digest = hasher.hexdigest()
    expected = expected_sha256 or EXPECTED_DOWNLOAD_SHA256.get(url)
    if expected and digest.lower() != expected.lower():
        raise RuntimeError(f"sha256 mismatch for {url}: expected {expected}, observed {digest}")
    return buffer.getvalue()


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle:
            tmp = Path(handle.name)
            handle.write(_get(url))
        tmp.replace(target)
    finally:
        if tmp is not None and tmp.exists():
            tmp.unlink()


def _zip_member_bytes(bundle: zipfile.ZipFile, member: str, *, max_bytes: int = MAX_ZIP_MEMBER_BYTES) -> bytes:
    info = bundle.getinfo(member)
    if info.file_size > max_bytes:
        raise RuntimeError(f"zip member exceeded {max_bytes} bytes: {member}")
    data = bundle.read(info)
    if len(data) > max_bytes:
        raise RuntimeError(f"zip member exceeded {max_bytes} bytes after read: {member}")
    return data


def _repair_gene_column(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    out = frame.copy()
    out[column] = out[column].map(_repair_excel_date_gene_symbol)
    return out


def _record_generated(path: Path, *, source_url: str, transform: str, inputs: list[Path] | None = None) -> None:
    write_source_sidecar(
        path,
        COMMAND,
        inputs=inputs or [],
        metadata={"generator": "hypoxia-source-builder", "source_url": source_url, "transform": transform},
    )


def _record_download(path: Path, url: str, *, study_id: str, transform: str = "passthrough") -> None:
    write_source_sidecar(
        path,
        shell_command(["curl", "-L", "-o", path, url]),
        metadata={"generator": "hypoxia-source-download", "source_url": url, "study_id": study_id, "transform": transform},
    )


def _save(df: pd.DataFrame, rel: str, *, source_url: str, transform: str, inputs: list[Path] | None = None) -> None:
    out = RAW / rel
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    _record_generated(out, source_url=source_url, transform=transform, inputs=inputs)
    print(f"  wrote {rel} ({len(df)} rows)")


def _repair_existing_csv(rel: str, column: str, *, source_url: str, transform: str) -> bool:
    path = RAW / rel
    if not path.exists():
        return False
    text = path.read_text()
    lines = text.splitlines()
    if lines and lines[0].split(",", 1)[0] == column:
        repaired = [lines[0]]
        for line in lines[1:]:
            gene, sep, rest = line.partition(",")
            repaired_gene = _repair_excel_date_gene_symbol(gene)
            repaired.append(f"{repaired_gene}{sep}{rest}")
        path.write_text("\n".join(repaired) + ("\n" if text.endswith(("\n", "\r\n")) else ""))
        _record_generated(path, source_url=source_url, transform=transform)
        print(f"  repaired {rel} ({len(lines) - 1} rows)")
    else:
        frame = pd.read_csv(path)
        frame = _repair_gene_column(frame, column)
        _save(frame, rel, source_url=source_url, transform=transform)
    return True


# ---- GEO passthrough: download IS the catalog source_path -----------------------------
GEO_PASSTHROUGH = [
    ("GSE76743", "GSE76743_HUVEC_HYPOXIA_mRNA_DESeq2.txt.gz", "GSE76743_HUVEC_HYPOXIA_mRNA_DESeq2.txt.gz"),
    ("GSE70544", "GSE70544_Normoxia_vs_Hypoxia_Gencode_gene_exp.diff.txt.gz", "iter7/GSE70544_Normoxia_vs_Hypoxia_Gencode_gene_exp.diff.txt.gz"),
    ("GSE311800", "GSE311800_DEG_Results_All_Comparisons.xlsx", "iter8/GSE311800_DEG_Results_All_Comparisons.xlsx"),
    ("GSE293238", "GSE293238_RNAseq_processed.xlsx", "iter9/GSE293238_RNAseq_processed.xlsx"),
    ("GSE313740", "GSE313740_processed_dataxGEO_14dic25.xlsx", "iter12/GSE313740_processed_dataxGEO_14dic25.xlsx"),
    ("GSE225253", "GSE225253_DESeq2_HK-2_hypoxia_vs_HK-2_normoxia_all.xlsx", "iter13/GSE225253_DESeq2_HK-2_hypoxia_vs_HK-2_normoxia_all.xlsx"),
    ("GSE283446", "GSE283446_TableS8.xlsx", "iter15/GSE283446_TableS8.xlsx"),
    # GSE132624 edgeR per-cell-line files (24 h) downloaded under the catalog's renamed paths
    ("GSE132624", "GSE132624_edgeR_501mel_DMEM_1_O2_24h_relative_to_0h_QLF-test_filtered_1cpm_3rep.csv.gz", "iter6/GSE132624_edgeR_501mel_24h.csv.gz"),
    ("GSE132624", "GSE132624_edgeR_IGR37_DMEM_1_O2_24h_relative_to_0h_QLF-test_filtered_1cpm_3rep.csv.gz", "iter6/GSE132624_edgeR_IGR37_24h.csv.gz"),
    ("GSE132624", "GSE132624_edgeR_IGR39_DMEM_1_O2_24h_relative_to_0h_QLF-test_filtered_1cpm_3rep.csv.gz", "iter6/GSE132624_edgeR_IGR39_24h.csv.gz"),
]

MENDELEY_URL = "https://data.mendeley.com/public-files/datasets/z42wpkbb8k/files/66d0f53c-1a5c-4ff5-8f96-4d8f1b713c0c/file_downloaded"
EUROPEPMC_BAUER = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9144826/supplementaryFiles"
TFEA_RDA_URL = "https://raw.githubusercontent.com/LauraPS1/TFEA.ChIP/master/data/hypoxia.rda"
EUROPEPMC_KINDRICK = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13181004/supplementaryFiles"
HYP008_MENDELEY_REL = "iter3/mendeley_z42wpkbb8k_hypoxia_vs_normoxia.csv"
HYP008_TRANSFORM = (
    "Mendeley CSV subset; source log2FoldChange sign-flipped from normoxia-vs-hypoxia "
    "to hypoxia-vs-normoxia; repaired Excel-date MARCH/SEPT/DEC gene symbols"
)
HYP008_DIRECTION_SENTINELS = frozenset({"BNIP3", "CA9", "EGLN3", "NDRG1", "SLC2A1", "VEGFA"})


def _refseq_to_symbol() -> dict[str, str]:
    if not HGNC.exists():
        _download(HGNC_URL, HGNC)
    _record_download(HGNC, HGNC_URL, study_id="HGNC", transform="reference_symbol_map")
    hg = pd.read_csv(HGNC, sep="\t", usecols=["symbol", "refseq_accession"], dtype=str).dropna(subset=["refseq_accession"])
    return dict(zip(hg["refseq_accession"].str.split(".").str[0], hg["symbol"]))


def _read_tsv_gz(rel: str) -> pd.DataFrame:
    with gzip.open(RAW / rel, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def _hyp008_marker_median(frame: pd.DataFrame) -> float | None:
    values = pd.to_numeric(
        frame.loc[frame["hgnc_symbol"].astype(str).isin(HYP008_DIRECTION_SENTINELS), "hypoxia_log2FoldChange"],
        errors="coerce",
    ).dropna()
    if values.empty:
        return None
    return float(values.median())


def _normalize_hyp008_mendeley(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "log2FoldChange" in out.columns:
        out["hypoxia_log2FoldChange"] = -pd.to_numeric(out["log2FoldChange"], errors="coerce")
    elif "hypoxia_log2FoldChange" in out.columns:
        out["hypoxia_log2FoldChange"] = pd.to_numeric(out["hypoxia_log2FoldChange"], errors="coerce")
        marker_median = _hyp008_marker_median(out)
        if marker_median is not None and marker_median < 0:
            out["hypoxia_log2FoldChange"] = -out["hypoxia_log2FoldChange"]
    else:
        raise ValueError("HYP008 Mendeley table must contain log2FoldChange or hypoxia_log2FoldChange")
    out = _repair_gene_column(out, "hgnc_symbol")
    return out[["hgnc_symbol", "hypoxia_log2FoldChange", "pvalue", "padj"]]


def _write_hyp008_mendeley() -> None:
    target = RAW / HYP008_MENDELEY_REL
    if target.exists():
        source = pd.read_csv(target)
    else:
        source = pd.read_csv(io.BytesIO(_get(MENDELEY_URL)))
    _save(
        _normalize_hyp008_mendeley(source),
        HYP008_MENDELEY_REL,
        source_url=MENDELEY_URL,
        transform=HYP008_TRANSFORM,
    )


def main() -> int:
    print("[GEO passthrough]")
    for accession, filename, rel in GEO_PASSTHROUGH:
        url = _geo_suppl(accession, filename)
        target = RAW / rel
        if target.exists():
            _record_download(target, url, study_id=accession)
            print(f"  present {rel}")
            continue
        _download(url, target)
        _record_download(target, url, study_id=accession)
        print(f"  fetched {rel}")

    print("[Mendeley HYP008]")
    _write_hyp008_mendeley()

    print("[MDPI Bauer HYP006/HYP007 via Europe PMC]")
    bauer_outputs = [
        ("log2FC_AH_N", "padj_AH_N", "iter3/bauer_2022_thp1_total_mrna_AH_vs_N.csv"),
        ("log2FoldChange_CH_N", "padj_CH_N", "iter3/bauer_2022_thp1_total_mrna_CH_vs_N.csv"),
    ]
    if all(
        _repair_existing_csv(
            rel,
            "gene_symbol",
            source_url=EUROPEPMC_BAUER,
            transform=f"Europe PMC MDPI supplement nested zip; DGE_Counts_log2FC columns {lfc_col}/{padj_col}; repaired Excel-date MARCH/SEPT/DEC gene symbols",
        )
        for lfc_col, padj_col, rel in bauer_outputs
    ):
        pass
    else:
        bundle = zipfile.ZipFile(io.BytesIO(_get(EUROPEPMC_BAUER)))
        inner = zipfile.ZipFile(io.BytesIO(_zip_member_bytes(bundle, "ijms-23-05824-s001.zip")))
        b = pd.read_excel(io.BytesIO(_zip_member_bytes(inner, "Bauer et al - Table S1.xlsx")), sheet_name="DGE_Counts_log2FC")
        for lfc_col, padj_col, rel in bauer_outputs:
            out = pd.DataFrame(
                {
                    "gene_symbol": b["external_gene_name"],
                    "hypoxia_log2FoldChange": pd.to_numeric(b[lfc_col], errors="coerce"),
                    "pvalue": pd.to_numeric(b[padj_col], errors="coerce"),  # only adjusted p is provided
                    "padj": pd.to_numeric(b[padj_col], errors="coerce"),
                }
            ).dropna(subset=["gene_symbol"])
            out = _repair_gene_column(out, "gene_symbol")
            _save(
                out,
                rel,
                source_url=EUROPEPMC_BAUER,
                transform=f"Europe PMC MDPI supplement nested zip; DGE_Counts_log2FC columns {lfc_col}/{padj_col}; repaired Excel-date MARCH/SEPT/DEC gene symbols",
            )

    print("[Bioconductor TFEA.ChIP HYP003 via Rscript]")
    target = RAW / "TFEA_hypoxia_symbol.csv"
    if target.exists():
        _record_generated(target, source_url=TFEA_RDA_URL, transform="Rscript load hypoxia data frame from TFEA.ChIP .rda")
        print(f"  present TFEA_hypoxia_symbol.csv ({sum(1 for _ in open(target)) - 1} rows)")
    else:
        with tempfile.TemporaryDirectory(prefix="degora_tfea_") as tmpdir:
            rda = Path(tmpdir) / "hypoxia.rda"
            rda.write_bytes(_get(TFEA_RDA_URL))
            subprocess.run(
                ["Rscript", "-e", f'load("{rda}"); write.csv(hypoxia, "{target}", row.names=FALSE)'],
                check=True,
            )
        _record_generated(target, source_url=TFEA_RDA_URL, transform="Rscript load hypoxia data frame from TFEA.ChIP .rda")
        print(f"  wrote TFEA_hypoxia_symbol.csv ({sum(1 for _ in open(target)) - 1} rows)")

    print("[transform HYP013 RefSeq->symbol]")
    ref2sym = _refseq_to_symbol()
    gse108676_url = _geo_suppl("GSE108676", "GSE108676_Human_gene_exp.diff.gz")
    gse108676_raw = RAW / "iter7/GSE108676_Human_gene_exp.diff.gz"
    if not gse108676_raw.exists():
        _download(gse108676_url, gse108676_raw)
    _record_download(gse108676_raw, gse108676_url, study_id="GSE108676", transform="raw_refseq_table_for_symbol_mapping")
    c = _read_tsv_gz("iter7/GSE108676_Human_gene_exp.diff.gz")
    c["gene_symbol"] = c["gene"].astype(str).str.split(".").str[0].map(ref2sym)
    c = c.dropna(subset=["gene_symbol"]).rename(
        columns={"log2.fold_change.": "hypoxia_log2FoldChange", "p_value": "pvalue", "q_value": "padj"}
    )
    _save(
        c[["gene_symbol", "hypoxia_log2FoldChange", "pvalue", "padj"]],
        "iter7/GSE108676_human_refseq_mapped_hypoxia_vs_normoxia.csv",
        source_url=gse108676_url,
        transform="RefSeq accession stripped and mapped through HGNC complete set",
        inputs=[gse108676_raw, HGNC],
    )

    print("[transform HYP017/HYP021 sign-flip]")
    gse255352_url = _geo_suppl("GSE255352", "GSE255352_OvsH_deg.xls.gz")
    gse255352_raw = RAW / "iter11/GSE255352_OvsH_deg.xls.gz"
    if not gse255352_raw.exists():
        _download(gse255352_url, gse255352_raw)
    _record_download(gse255352_raw, gse255352_url, study_id="GSE255352", transform="raw_normoxia_vs_hypoxia_table")
    a = _read_tsv_gz("iter11/GSE255352_OvsH_deg.xls.gz")
    a["hypoxia_vs_normoxia_log2FoldChange"] = -pd.to_numeric(a["log2FoldChange"], errors="coerce")
    a = _repair_gene_column(a, "gene_name")
    _save(
        a[["gene_name", "hypoxia_vs_normoxia_log2FoldChange", "pvalue", "padj", "gene_biotype"]],
        "iter11/GSE255352_HvsO_deg_signflipped.csv",
        source_url=gse255352_url,
        transform="sign-flipped O(normoxia)-vs-H(hypoxia) log2FoldChange to H-vs-O; repaired Excel-date MARCH/SEPT/DEC gene symbols",
        inputs=[gse255352_raw],
    )

    gse160491_url = _geo_suppl("GSE160491", "GSE160491_Processed_all_compare.xls.gz")
    gse160491_raw = RAW / "iter15/GSE160491_Processed_all_compare.xls.gz"
    if not gse160491_raw.exists():
        _download(gse160491_url, gse160491_raw)
    _record_download(gse160491_raw, gse160491_url, study_id="GSE160491", transform="raw_normoxia_vs_hypoxia_table")
    g = _read_tsv_gz("iter15/GSE160491_Processed_all_compare.xls.gz")
    g["hypoxia_vs_normoxia_log2FoldChange"] = -pd.to_numeric(g["N_siConvsH_siCon_log2FoldChange"], errors="coerce")
    g = g.rename(columns={"N_siConvsH_siCon_pvalue": "pvalue", "N_siConvsH_siCon_padj": "padj"})
    g = _repair_gene_column(g, "gene_name")
    _save(
        g[["gene_name", "hypoxia_vs_normoxia_log2FoldChange", "pvalue", "padj", "gene_biotype"]],
        "iter15/GSE160491_siControl_hypoxia_vs_normoxia_signflipped.csv",
        source_url=gse160491_url,
        transform="sign-flipped N_siCon-vs-H_siCon log2FoldChange to hypoxia-vs-normoxia; repaired Excel-date MARCH/SEPT/DEC gene symbols",
        inputs=[gse160491_raw],
    )

    print("[Springer Kindrick HYP023/HYP024 -- PC3/HCT116 RNA-seq from Figure 4 source data]")
    kindrick_outputs = [
        ("Figure 4A&C source data ", "iter16/GSE296192_PC3_hypoxia_vs_normoxia.csv"),
        ("Figure 4B&D source data", "iter16/GSE296192_HCT116_hypoxia_vs_normoxia.csv"),
    ]
    if all(
        _repair_existing_csv(
            rel,
            "gene_symbol",
            source_url=EUROPEPMC_KINDRICK,
            transform=f"Europe PMC/Springer ESM2 sheet {sheet!r}; repaired Excel-date MARCH/SEPT/DEC gene symbols",
        )
        for sheet, rel in kindrick_outputs
    ):
        pass
    else:
        bundle = zipfile.ZipFile(io.BytesIO(_get(EUROPEPMC_KINDRICK)))
        esm = pd.ExcelFile(io.BytesIO(_zip_member_bytes(bundle, "42003_2026_9875_MOESM2_ESM.xlsx")))
        for sheet, rel in kindrick_outputs:
            df = pd.read_excel(esm, sheet_name=sheet, header=1).rename(columns=lambda c: str(c))
            df = df.rename(columns={df.columns[0]: "gene_symbol"})
            adj = next(c for c in df.columns if "adjusted" in c.lower() and "p" in c.lower())
            out = pd.DataFrame(
                {
                    "gene_symbol": df["gene_symbol"].astype(str).str.strip(),
                    "hypoxia_log2FoldChange": pd.to_numeric(df["foldchange"], errors="coerce"),
                    "pvalue": pd.to_numeric(df[adj], errors="coerce"),
                    "padj": pd.to_numeric(df[adj], errors="coerce"),
                }
            ).dropna(subset=["gene_symbol", "hypoxia_log2FoldChange"])
            out = _repair_gene_column(out, "gene_symbol")
            out = out[out["gene_symbol"].ne("nan") & out["gene_symbol"].ne("")]
            _save(
                out,
                rel,
                source_url=EUROPEPMC_KINDRICK,
                transform=f"Europe PMC/Springer ESM2 sheet {sheet!r}; repaired Excel-date MARCH/SEPT/DEC gene symbols",
            )

    print("\nhypoxia RNA-seq sources reconstructed. Run `make pipeline TOPIC=hypoxia` to re-harmonize from raw.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
