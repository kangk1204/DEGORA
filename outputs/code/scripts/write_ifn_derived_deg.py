#!/usr/bin/env python
"""Derive a preregistered IFN-response pilot corpus from public count matrices.

This is a derived-count pilot, not an as-published DEG-table corpus. It exists
to test whether DEGORA recovers canonical interferon-stimulated genes when
the topic changes from hypoxia to IFN response.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from degora.derived_counts import attach_low_count_filter_metadata, low_count_filter_mask, low_count_filter_summary
from degora.provenance import write_source_sidecar


IFN_GOLD_GENES = [
    "ISG15",
    "RSAD2",
    "CMPK2",
    "MX1",
    "MX2",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "IFI44",
    "IFI44L",
    "IFI6",
    "IFITM1",
    "USP18",
    "DDX58",
    "IFIH1",
    "HERC5",
    "STAT1",
    "IRF7",
    "IFI27",
    "BST2",
    "ISG20",
    "IRF9",
    "STAT2",
    "EIF2AK2",
]


RAW_INPUTS = {
    "gse147507": {
        "path": "data/deg/raw/ifn/GSE147507_RawReadCounts_Human.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz",
    },
    "gse221804": {
        "path": "data/deg/raw/ifn/GSE221804_raw_counts_WT_IFNa.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221804/suppl/GSE221804_raw_counts_WT_IFNa.csv.gz",
    },
    "hgnc": {
        "path": "data/deg/raw/ifn/hgnc_complete_set.txt",
        "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    },
}


@dataclass(frozen=True)
class ContrastSpec:
    study_id: str
    paper_id: str
    source_url: str
    source_path: Path
    raw_source_key: str
    species: str
    cell_system: str
    condition: str
    duration_h: str
    n_ctrl: int
    n_treat: int
    control_columns: tuple[str, ...]
    treat_columns: tuple[str, ...]
    output_name: str


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).clip(0.0, 1.0)
    n = len(values)
    order = np.argsort(values.to_numpy(dtype=float))
    ranked = values.to_numpy(dtype=float)[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(out, index=pvalues.index)


def _log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    libraries = numeric.sum(axis=0).replace(0.0, np.nan)
    cpm = numeric.divide(libraries, axis=1) * 1_000_000.0
    return np.log2(cpm.fillna(0.0) + 1.0)


def _collapse_by_gene(counts: pd.DataFrame) -> pd.DataFrame:
    frame = counts.copy()
    frame.index = frame.index.astype("string").str.replace(r"\.\d+$", "", regex=True).str.upper()
    frame = frame.loc[frame.index.notna() & (frame.index != "")]
    return frame.groupby(frame.index).sum(numeric_only=True)


def _read_symbol_counts(path: Path, *, sep: str) -> pd.DataFrame:
    counts = pd.read_csv(path, sep=sep, index_col=0)
    return _collapse_by_gene(counts)


def _read_hgnc_reference(path: Path) -> tuple[dict[str, str], set[str]]:
    hgnc = pd.read_csv(path, sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id", "locus_group"])
    hgnc = hgnc.dropna(subset=["symbol", "ensembl_gene_id"])
    hgnc = hgnc.loc[hgnc["locus_group"].astype("string").str.lower().eq("protein-coding gene")].copy()
    hgnc["ensembl_gene_id"] = hgnc["ensembl_gene_id"].str.replace(r"\.\d+$", "", regex=True).str.upper()
    hgnc["symbol"] = hgnc["symbol"].str.upper()
    return dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"], strict=False)), set(hgnc["symbol"])


def _read_hgnc_map(path: Path) -> dict[str, str]:
    mapping, _ = _read_hgnc_reference(path)
    return mapping


def _read_ensembl_counts(path: Path, hgnc_map: dict[str, str]) -> pd.DataFrame:
    raw = pd.read_csv(path, index_col=0)
    ensembl_ids = raw.index.astype("string").str.replace(r"\.\d+$", "", regex=True).str.upper()
    symbols = pd.Series(ensembl_ids.map(hgnc_map), index=raw.index)
    raw = raw.loc[symbols.notna()].copy()
    raw.index = symbols.loc[symbols.notna()].astype(str)
    return _collapse_by_gene(raw)


def _derive_contrast(counts: pd.DataFrame, spec: ContrastSpec) -> pd.DataFrame:
    missing = sorted(set(spec.control_columns + spec.treat_columns).difference(counts.columns))
    if missing:
        raise ValueError(f"{spec.study_id} missing columns: {missing}")

    selected_counts = counts[list(spec.control_columns + spec.treat_columns)]
    expressed = low_count_filter_mask(selected_counts)
    filter_summary = low_count_filter_summary(selected_counts, expressed)
    log_counts = _log_cpm(selected_counts).loc[expressed]
    control = log_counts[list(spec.control_columns)]
    treat = log_counts[list(spec.treat_columns)]
    log2fc = treat.mean(axis=1) - control.mean(axis=1)
    test = stats.ttest_ind(treat.to_numpy(dtype=float), control.to_numpy(dtype=float), axis=1, equal_var=False, nan_policy="omit")
    pvalue = pd.Series(test.pvalue, index=log_counts.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    out = pd.DataFrame(
        {
            "gene_symbol": log_counts.index.astype(str),
            "log2FoldChange": log2fc,
            "pvalue": pvalue,
        }
    )
    out["padj"] = _bh_adjust(out["pvalue"])
    out = attach_low_count_filter_metadata(out, filter_summary)
    out = out.sort_values(["pvalue", "gene_symbol"]).reset_index(drop=True)
    return out


def _write_gold_panel(path: Path, command: str) -> None:
    now = "2026-05-30"
    rows = [
        {
            "gene_symbol": gene,
            "expected_direction": "up",
            "panel_role": "core_canonical",
            "locked_at": now,
            "evidence_basis": (
                "Expanded canonical panel (MSigDB Hallmark + topic-specific target/ChIP "
                "literature); locked blind to DEGORA output."
            ),
            "notes": (
                "Type I/II interferon (ISG) core canonical target; expanded gold panel locked "
                "blind to DEGORA output; used only for recall/prioritization diagnostics."
            ),
        }
        for gene in IFN_GOLD_GENES
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, lineterminator="\n")
    write_source_sidecar(path, command, metadata={"generator": "ifn-derived-gold-panel"})


def _ensure_gold_panel(path: Path, command: str) -> list[str]:
    """Return locked IFN gold symbols, creating a panel only when it is absent."""

    if path.exists():
        frame = pd.read_csv(path)
        if "gene_symbol" not in frame.columns:
            raise ValueError(f"IFN gold panel {path} is missing gene_symbol")
        return frame["gene_symbol"].dropna().astype(str).str.strip().str.upper().loc[lambda s: s.ne("")].tolist()
    _write_gold_panel(path, command)
    return list(IFN_GOLD_GENES)


def _contrast_specs(raw_dir: Path) -> list[ContrastSpec]:
    gse147 = Path(RAW_INPUTS["gse147507"]["path"])
    gse221 = Path(RAW_INPUTS["gse221804"]["path"])
    mock = (
        "Series9_NHBE_Mock_1",
        "Series9_NHBE_Mock_2",
        "Series9_NHBE_Mock_3",
        "Series9_NHBE_Mock_4",
    )
    specs = [
        ContrastSpec(
            study_id="IFN_GSE147507_NHBE_IFNB_4h",
            paper_id="GSE147507_NHBE_IFNB",
            source_url=RAW_INPUTS["gse147507"]["url"],
            source_path=raw_dir / "GSE147507_IFNB_4h_vs_mock_derived_logcpm_welch.csv",
            raw_source_key="gse147507",
            species="Homo sapiens",
            cell_system="NHBE",
            condition="IFNB vs mock",
            duration_h="4",
            n_ctrl=4,
            n_treat=2,
            control_columns=mock,
            treat_columns=("Series9_NHBE_IFNB_4h_1", "Series9_NHBE_IFNB_4h_2"),
            output_name="GSE147507_IFNB_4h_vs_mock_derived_logcpm_welch.csv",
        ),
        ContrastSpec(
            study_id="IFN_GSE147507_NHBE_IFNB_6h",
            paper_id="GSE147507_NHBE_IFNB",
            source_url=RAW_INPUTS["gse147507"]["url"],
            source_path=raw_dir / "GSE147507_IFNB_6h_vs_mock_derived_logcpm_welch.csv",
            raw_source_key="gse147507",
            species="Homo sapiens",
            cell_system="NHBE",
            condition="IFNB vs mock",
            duration_h="6",
            n_ctrl=4,
            n_treat=2,
            control_columns=mock,
            treat_columns=("Series9_NHBE_IFNB_6h_1", "Series9_NHBE_IFNB_6h_2"),
            output_name="GSE147507_IFNB_6h_vs_mock_derived_logcpm_welch.csv",
        ),
        ContrastSpec(
            study_id="IFN_GSE147507_NHBE_IFNB_12h",
            paper_id="GSE147507_NHBE_IFNB",
            source_url=RAW_INPUTS["gse147507"]["url"],
            source_path=raw_dir / "GSE147507_IFNB_12h_vs_mock_derived_logcpm_welch.csv",
            raw_source_key="gse147507",
            species="Homo sapiens",
            cell_system="NHBE",
            condition="IFNB vs mock",
            duration_h="12",
            n_ctrl=4,
            n_treat=2,
            control_columns=mock,
            treat_columns=("Series9_NHBE_IFNB_12h_1", "Series9_NHBE_IFNB_12h_2"),
            output_name="GSE147507_IFNB_12h_vs_mock_derived_logcpm_welch.csv",
        ),
    ]
    un = ("hs_HuhWT_UN_rep1", "hs_HuhWT_UN_rep2")
    for duration in ("2", "4", "8", "24"):
        specs.append(
            ContrastSpec(
                study_id=f"IFN_GSE221804_HUHWT_IFNA_{duration}h",
                paper_id="GSE221804_HUHWT_IFNA",
                source_url=RAW_INPUTS["gse221804"]["url"],
                source_path=raw_dir / f"GSE221804_IFNA_{duration}h_vs_un_derived_logcpm_welch.csv",
                raw_source_key="gse221804",
                species="Homo sapiens",
                cell_system="Huh7.5 WT",
                condition="IFNa vs untreated",
                duration_h=duration,
                n_ctrl=2,
                n_treat=2,
                control_columns=un,
                treat_columns=(f"hs_HuhWT_IFNa_{duration}h_rep1", f"hs_HuhWT_IFNa_{duration}h_rep2"),
                output_name=f"GSE221804_IFNA_{duration}h_vs_un_derived_logcpm_welch.csv",
            )
        )
    return specs


def _catalog_rows(specs: list[ContrastSpec]) -> list[dict[str, Any]]:
    rows = []
    for spec in specs:
        rows.append(
            {
                "study_id": spec.study_id,
                "paper_id": spec.paper_id,
                "source_url": spec.source_url,
                "source_path": str(spec.source_path),
                "pipeline": "logCPM_Welch_derived_from_public_counts",
                "species": spec.species,
                "cell_system": spec.cell_system,
                "hypoxia_modality": spec.condition,
                "duration_h": spec.duration_h,
                "n_ctrl": spec.n_ctrl,
                "n_treat": spec.n_treat,
                "gene_column": "gene_symbol",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
                "sep": "",
                "sheet_name": "",
                "gene_type_column": "",
                "gene_type_keep": "",
                "include_in_analysis": True,
                "notes": (
                    "Derived-count IFN pilot contrast; log2FC and p-values computed from public count matrix "
                    "using logCPM and Welch t-test after an independent low-count expression filter. "
                    "Not an as-published DEG table."
                ),
                "source_input_type": "derived_count_table",
                "assay_type": "RNA-seq",
                "table_scope": "full_results",
            }
        )
    return rows


def write_ifn_derived_corpus(raw_dir: Path, catalog_path: Path, gold_path: Path, summary_path: Path, *, command: str) -> dict[str, Any]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    for source in RAW_INPUTS.values():
        path = Path(source["path"])
        if not path.exists():
            raise FileNotFoundError(f"Missing required raw input: {path}")
        write_source_sidecar(
            path,
            f"curl -L -o {path} {source['url']}",
            inputs=[],
            metadata={"generator": "manual-download-source", "source_url": source["url"]},
        )

    gold_genes = _ensure_gold_panel(gold_path, command)
    gold_gene_set = set(gold_genes)
    hgnc_map, protein_coding_symbols = _read_hgnc_reference(Path(RAW_INPUTS["hgnc"]["path"]))
    gse147_raw = _read_symbol_counts(Path(RAW_INPUTS["gse147507"]["path"]), sep="\t")
    gse147 = gse147_raw.loc[gse147_raw.index.isin(protein_coding_symbols)].copy()
    gse221 = _read_ensembl_counts(Path(RAW_INPUTS["gse221804"]["path"]), hgnc_map)
    counts_by_source = {"gse147507": gse147, "gse221804": gse221}
    specs = _contrast_specs(raw_dir)
    summaries = []
    for spec in specs:
        deg = _derive_contrast(counts_by_source[spec.raw_source_key], spec)
        deg_path = raw_dir / spec.output_name
        deg.to_csv(deg_path, index=False)
        raw_input = Path(RAW_INPUTS[spec.raw_source_key]["path"])
        write_source_sidecar(
            deg_path,
            command,
            inputs=[raw_input, RAW_INPUTS["hgnc"]["path"]] if spec.raw_source_key == "gse221804" else [raw_input],
            metadata={
                "generator": "ifn-derived-logcpm-welch",
                "study_id": spec.study_id,
                "raw_source_key": spec.raw_source_key,
                "derived_count_pilot": True,
            },
        )
        top20 = set(deg.head(20)["gene_symbol"])
        summaries.append(
            {
                "study_id": spec.study_id,
                "paper_id": spec.paper_id,
                "rows": int(len(deg)),
                "top20_gold_hits": sorted(top20.intersection(gold_gene_set)),
                "top_gene": str(deg.iloc[0]["gene_symbol"]),
                "output": str(deg_path),
            }
        )

    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    catalog = pd.DataFrame(_catalog_rows(specs))
    catalog.to_csv(catalog_path, index=False)
    write_source_sidecar(catalog_path, command, inputs=[Path(RAW_INPUTS["gse147507"]["path"]), Path(RAW_INPUTS["gse221804"]["path"]), Path(RAW_INPUTS["hgnc"]["path"])], metadata={"generator": "ifn-derived-catalog"})

    summary = {
        "corpus": "IFN response derived-count pilot",
        "claim_guardrail": "Derived from public count matrices; not an as-published DEG-table validation corpus.",
        "n_contrasts": len(specs),
        "n_source_units": int(catalog["paper_id"].nunique()),
        "gold_genes": gold_genes,
        "contrast_summaries": summaries,
        "catalog": str(catalog_path),
        "gold_panel": str(gold_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(summary_path, command, inputs=[catalog_path, gold_path], metadata={"generator": "ifn-derived-summary"})
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/deg/raw/ifn"))
    parser.add_argument("--catalog", type=Path, default=Path("data/studies/ifn_derived_catalog.csv"))
    parser.add_argument("--gold", type=Path, default=Path("data/studies/gold/ifn_gold_panel.csv"))
    parser.add_argument("--summary", type=Path, default=Path("outputs/results/ifn-pilot/ifn_derived_deg_summary.json"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_ifn_derived_deg.py"
        f" --raw-dir {args.raw_dir}"
        f" --catalog {args.catalog}"
        f" --gold {args.gold}"
        f" --summary {args.summary}"
    )
    summary = write_ifn_derived_corpus(args.raw_dir, args.catalog, args.gold, args.summary, command=command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
