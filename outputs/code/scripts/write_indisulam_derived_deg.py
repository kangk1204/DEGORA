#!/usr/bin/env python
"""Derive a guarded Indisulam DEGORA pilot from public count matrices."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from degora.provenance import write_source_sidecar


RAW_INPUTS = {
    "gse223011": {
        "path": "data/deg/raw/indisulam/GSE223011_Nijhuis_KURAMOCHI_VC_INDISULAM_humanCounts.tsv.xlsx",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223011/suppl/GSE223011_Nijhuis_KURAMOCHI_VC_INDISULAM_humanCounts.tsv.xlsx",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223011",
    },
    "gse268568": {
        "path": "data/deg/raw/indisulam/GSE268568_RawCount_Matrix_ann.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268568/suppl/GSE268568_RawCount_Matrix_ann.txt.gz",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268568",
    },
    "hgnc": {
        "path": "data/deg/raw/ifn/hgnc_complete_set.txt",
        "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        "geo": "",
    },
}

GUARDRAIL_FLAGS = {
    "derived_count_pilot": True,
    "as_published_author_deg_validation": False,
    "primary_min_studies": 2,
    "claim_scope": "exploratory_cross_source_pilot",
    "microarray_admissibility": (
        "Microarray contrasts are admissible as a separate secondary evidence tier only when a direct "
        "Indisulam/E7070 treatment contrast has raw array data, a normalized expression matrix suitable "
        "for limma, or a full author DEG table with gene-level logFC and p-value."
    ),
}

INDISULAM_ANCHORS = [
    ("RBM39", "primary_degradation_substrate", "Han2017_Science_PMID28302793; Faust2020_NatChemBiol_PMID31844297"),
    ("DCAF15", "core_glue_machinery", "Han2017_Science_PMID28302793; Faust2020_NatChemBiol_PMID31844297"),
    ("DDB1", "core_glue_machinery", "Faust2020_NatChemBiol_PMID31844297"),
    ("DDA1", "core_glue_machinery", "Faust2020_NatChemBiol_PMID31844297"),
    ("CUL4A", "core_glue_machinery", "Faust2020_NatChemBiol_PMID31844297"),
    ("CUL4B", "core_glue_machinery", "Faust2020_NatChemBiol_PMID31844297"),
    ("RBM23", "related_substrate_selectivity", "Faust2020_NatChemBiol_PMID31844297"),
    ("CDK4", "downstream_cell_cycle_or_metabolic_response", "Sikka2022_NatCommun_PMID35296644"),
    ("TYMS", "downstream_cell_cycle_or_metabolic_response", "Sikka2022_NatCommun_PMID35296644"),
    ("THOC1", "downstream_splicing_or_transcription_context", "Indisulam_RBM39_splicing_literature"),
    ("EZH2", "downstream_cancer_context", "Indisulam_RBM39_splicing_literature"),
    ("ZMYND8", "downstream_cancer_context", "Indisulam_RBM39_splicing_literature"),
    ("MYC", "PARP_DDR_or_MYC_context", "Sikka2022_NatCommun_PMID35296644"),
    ("ATM", "PARP_DDR_context", "Xu2023_CellReports_GSE223011"),
    ("BRCA1", "PARP_DDR_context", "Xu2023_CellReports_GSE223011"),
    ("RAD51", "PARP_DDR_context", "Xu2023_CellReports_GSE223011"),
]


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _normalize_symbol(symbols: pd.Series) -> pd.Series:
    return symbols.astype("string").str.strip().str.upper()


def _collapse_by_symbol(counts: pd.DataFrame) -> pd.DataFrame:
    frame = counts.copy()
    frame.index = _normalize_symbol(pd.Series(frame.index, index=frame.index))
    frame = frame.loc[frame.index.notna() & (frame.index != "")]
    return frame.groupby(frame.index).sum(numeric_only=True)


def _read_hgnc_map(path: Path) -> dict[str, str]:
    hgnc = pd.read_csv(path, sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id"])
    hgnc = hgnc.dropna(subset=["symbol", "ensembl_gene_id"])
    hgnc["ensembl_gene_id"] = hgnc["ensembl_gene_id"].str.replace(r"\.\d+$", "", regex=True).str.upper()
    hgnc["symbol"] = hgnc["symbol"].str.upper()
    return dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"], strict=False))


def _map_ensembl_to_symbol(frame: pd.DataFrame, gene_column: str, hgnc_map: dict[str, str]) -> tuple[pd.Series, int]:
    ids = frame[gene_column].astype("string").str.replace(r"\.\d+$", "", regex=True).str.upper()
    symbols = ids.map(hgnc_map)
    n_unmapped = int(symbols.isna().sum())
    return symbols, n_unmapped


def _read_gse223011_counts(path: Path, hgnc_map: dict[str, str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = pd.read_excel(path, sheet_name="humanCounts")
    count_columns = [
        "../Human/3. Merged BAM/KURA_INDY-1_merged.bam",
        "../Human/3. Merged BAM/KURA_INDY-2_merged.bam",
        "../Human/3. Merged BAM/KURA_INDY-3_merged.bam",
        "../Human/3. Merged BAM/KURA_VC-1_merged.bam",
        "../Human/3. Merged BAM/KURA_VC-2_merged.bam",
        "../Human/3. Merged BAM/KURA_VC-3_merged.bam",
    ]
    missing = sorted(set(["Geneid", *count_columns]).difference(frame.columns))
    if missing:
        raise ValueError(f"GSE223011 missing columns: {missing}")
    symbols, n_unmapped = _map_ensembl_to_symbol(frame, "Geneid", hgnc_map)
    counts = frame.loc[symbols.notna(), count_columns].copy()
    counts.index = symbols.loc[symbols.notna()].astype(str)
    counts = _collapse_by_symbol(counts)
    return counts, {
        "raw_rows": int(len(frame)),
        "mapped_rows": int(symbols.notna().sum()),
        "n_unmapped_genes": n_unmapped,
        "n_gene_symbols_after_collapse": int(len(counts)),
    }


def _read_gse268568_counts(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = pd.read_csv(path, sep="\t")
    count_columns = ["WT1", "WT2", "WT3", "WTindi1", "WTindi2", "WTindi3", "KO1", "KO2", "KO3", "KOindi1", "KOindi2", "KOindi3"]
    missing = sorted(set(["gene_sym", *count_columns]).difference(frame.columns))
    if missing:
        raise ValueError(f"GSE268568 missing columns: {missing}")
    symbols = _normalize_symbol(frame["gene_sym"])
    counts = frame.loc[symbols.notna() & symbols.ne(""), count_columns].copy()
    counts.index = symbols.loc[symbols.notna() & symbols.ne("")].astype(str)
    counts = _collapse_by_symbol(counts)
    return counts, {
        "raw_rows": int(len(frame)),
        "mapped_rows": int(symbols.notna().sum()),
        "n_unmapped_genes": int(symbols.isna().sum() + symbols.eq("").sum()),
        "n_gene_symbols_after_collapse": int(len(counts)),
    }


def _derive_contrast(counts: pd.DataFrame, spec: ContrastSpec) -> pd.DataFrame:
    columns = list(spec.control_columns + spec.treat_columns)
    missing = sorted(set(columns).difference(counts.columns))
    if missing:
        raise ValueError(f"{spec.study_id} missing columns: {missing}")
    log_counts = _log_cpm(counts[columns])
    control = log_counts[list(spec.control_columns)]
    treat = log_counts[list(spec.treat_columns)]
    log2fc = treat.mean(axis=1) - control.mean(axis=1)
    test = stats.ttest_ind(treat.to_numpy(dtype=float), control.to_numpy(dtype=float), axis=1, equal_var=False, nan_policy="omit")
    pvalue = pd.Series(test.pvalue, index=log_counts.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    out = pd.DataFrame({"gene_symbol": log_counts.index.astype(str), "log2FoldChange": log2fc, "pvalue": pvalue})
    out["padj"] = _bh_adjust(out["pvalue"])
    return out.sort_values(["pvalue", "gene_symbol"]).reset_index(drop=True)


def _contrast_specs(raw_dir: Path) -> list[ContrastSpec]:
    return [
        ContrastSpec(
            study_id="INDISULAM_GSE223011_KURAMOCHI_24h",
            paper_id="GSE223011_Xu2023",
            source_url=RAW_INPUTS["gse223011"]["geo"],
            source_path=raw_dir / "GSE223011_KURAMOCHI_indisulam_24h_vs_vehicle_derived_logcpm_welch.csv",
            raw_source_key="gse223011",
            species="Homo sapiens",
            cell_system="KURAMOCHI high-grade serous ovarian carcinoma",
            condition="Indisulam 10uM vs vehicle",
            duration_h="24",
            n_ctrl=3,
            n_treat=3,
            control_columns=(
                "../Human/3. Merged BAM/KURA_VC-1_merged.bam",
                "../Human/3. Merged BAM/KURA_VC-2_merged.bam",
                "../Human/3. Merged BAM/KURA_VC-3_merged.bam",
            ),
            treat_columns=(
                "../Human/3. Merged BAM/KURA_INDY-1_merged.bam",
                "../Human/3. Merged BAM/KURA_INDY-2_merged.bam",
                "../Human/3. Merged BAM/KURA_INDY-3_merged.bam",
            ),
            output_name="GSE223011_KURAMOCHI_indisulam_24h_vs_vehicle_derived_logcpm_welch.csv",
        ),
        ContrastSpec(
            study_id="INDISULAM_GSE268568_CAL27_WT_24h",
            paper_id="GSE268568_Ando2024",
            source_url=RAW_INPUTS["gse268568"]["geo"],
            source_path=raw_dir / "GSE268568_CAL27_WT_indisulam_24h_vs_DMSO_derived_logcpm_welch.csv",
            raw_source_key="gse268568",
            species="Homo sapiens",
            cell_system="CAL27 WT",
            condition="Indisulam 1uM vs DMSO",
            duration_h="24",
            n_ctrl=3,
            n_treat=3,
            control_columns=("WT1", "WT2", "WT3"),
            treat_columns=("WTindi1", "WTindi2", "WTindi3"),
            output_name="GSE268568_CAL27_WT_indisulam_24h_vs_DMSO_derived_logcpm_welch.csv",
        ),
        ContrastSpec(
            study_id="INDISULAM_GSE268568_CAL27_LATS12KO_24h",
            paper_id="GSE268568_Ando2024",
            source_url=RAW_INPUTS["gse268568"]["geo"],
            source_path=raw_dir / "GSE268568_CAL27_LATS12KO_indisulam_24h_vs_DMSO_derived_logcpm_welch.csv",
            raw_source_key="gse268568",
            species="Homo sapiens",
            cell_system="CAL27 LATS1/2 KO",
            condition="Indisulam 1uM vs DMSO",
            duration_h="24",
            n_ctrl=3,
            n_treat=3,
            control_columns=("KO1", "KO2", "KO3"),
            treat_columns=("KOindi1", "KOindi2", "KOindi3"),
            output_name="GSE268568_CAL27_LATS12KO_indisulam_24h_vs_DMSO_derived_logcpm_welch.csv",
        ),
    ]


def _catalog_rows(specs: list[ContrastSpec], row_counts: dict[str, int]) -> list[dict[str, Any]]:
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
                "table_scope": "full_results",
                "rank_universe_size": row_counts.get(spec.study_id, ""),
                "include_in_analysis": True,
                "notes": (
                    "Derived-count Indisulam pilot contrast; log2FC and p-values computed from public count matrix "
                    "using logCPM and Welch t-test. Not an as-published author DEG table."
                ),
            }
        )
    return rows


def _write_gold_panel(path: Path, command: str) -> None:
    rows = [
        {
            "gene_symbol": symbol,
            "anchor_role": role,
            "claim_scope": "mechanism_or_pathway_anchor_not_directional_expression_claim",
            "expected_direction": "not_asserted",
            "evidence_basis": evidence,
            "source_reference": evidence,
        }
        for symbol, role, evidence in INDISULAM_ANCHORS
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    write_source_sidecar(path, command, metadata={"generator": "indisulam-anchor-panel", **GUARDRAIL_FLAGS})


def _availability_records(raw_paths: dict[str, Path]) -> list[dict[str, Any]]:
    records = [
        {
            "source_id": "GSE223011",
            "assay_type": "RNA-seq",
            "status": "active_primary_count_source",
            "direct_treatment_design": "KURAMOCHI vehicle vs 10uM Indisulam 24h, 3 vs 3",
            "author_deg_table_found": False,
            "local_path": str(raw_paths["gse223011"]),
            "local_size_bytes": raw_paths["gse223011"].stat().st_size if raw_paths["gse223011"].exists() else 0,
            "sha256": _sha256(raw_paths["gse223011"]) if raw_paths["gse223011"].exists() else "",
            "source_url": RAW_INPUTS["gse223011"]["geo"],
            "use_decision": "scored_as_derived_count_pilot",
            "reason": "public count matrix has clear treated/control samples; no author full DEG table used here",
        },
        {
            "source_id": "GSE268568",
            "assay_type": "RNA-seq",
            "status": "active_primary_count_source",
            "direct_treatment_design": "CAL27 WT and LATS1/2 KO DMSO vs 1uM Indisulam 24h, 3 vs 3 each",
            "author_deg_table_found": False,
            "local_path": str(raw_paths["gse268568"]),
            "local_size_bytes": raw_paths["gse268568"].stat().st_size if raw_paths["gse268568"].exists() else 0,
            "sha256": _sha256(raw_paths["gse268568"]) if raw_paths["gse268568"].exists() else "",
            "source_url": RAW_INPUTS["gse268568"]["geo"],
            "use_decision": "scored_as_derived_count_pilot",
            "reason": "public annotated raw-count matrix has clear treated/control samples",
        },
        {
            "source_id": "GSE164505",
            "assay_type": "RNA-seq",
            "status": "deferred_secondary_source",
            "direct_treatment_design": "neuroblastoma Indisulam/RBM39 perturbation RNA-seq, multiple samples",
            "author_deg_table_found": "not_evaluated_for_this_run",
            "local_path": "",
            "local_size_bytes": "",
            "sha256": "",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164505",
            "use_decision": "deferred",
            "reason": "larger multi-condition design requires separate contrast parsing",
        },
        {
            "source_id": "GSE254842",
            "assay_type": "RNA-seq_or_processed_expression",
            "status": "deferred_secondary_source",
            "direct_treatment_design": "MDA-MB-231 DMSO/Indisulam/Pladienolide B candidate",
            "author_deg_table_found": "not_evaluated_for_this_run",
            "local_path": "",
            "local_size_bytes": "",
            "sha256": "",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254842",
            "use_decision": "deferred",
            "reason": "processed-matrix derivation and treatment labeling should be audited separately",
        },
        {
            "source_id": "E7070_microarray_literature_PMID12467223",
            "assay_type": "microarray",
            "status": "candidate_secondary_evidence_tier",
            "direct_treatment_design": "E7070/Indisulam sulfonamide expression profiling reported in literature",
            "author_deg_table_found": "not_recovered_in_this_run",
            "local_path": "",
            "local_size_bytes": "",
            "sha256": "",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/12467223/",
            "use_decision": "eligible_only_after_raw_or_full_gene_level_table_recovery",
            "reason": "microarray is usable through limma/rank-level harmonization, but not mixed into this primary count-derived pilot without source table recovery",
        },
    ]
    return records


def write_indisulam_derived_corpus(raw_dir: Path, catalog_path: Path, gold_path: Path, output_dir: Path, *, command: str) -> dict[str, Any]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_paths = {key: Path(info["path"]) for key, info in RAW_INPUTS.items()}
    for key, info in RAW_INPUTS.items():
        path = Path(info["path"])
        if not path.exists():
            raise FileNotFoundError(f"Missing required raw input: {path}")
        write_source_sidecar(
            path,
            f"curl -L -o {path} {info['url']}",
            inputs=[],
            metadata={"generator": "manual-download-source", "source_url": info["url"]},
        )

    hgnc_map = _read_hgnc_map(raw_paths["hgnc"])
    gse223, gse223_map = _read_gse223011_counts(raw_paths["gse223011"], hgnc_map)
    gse268, gse268_map = _read_gse268568_counts(raw_paths["gse268568"])
    counts_by_source = {"gse223011": gse223, "gse268568": gse268}
    mapping_summaries = {"gse223011": gse223_map, "gse268568": gse268_map}

    specs = _contrast_specs(raw_dir)
    contrast_summaries = []
    row_counts = {}
    for spec in specs:
        deg = _derive_contrast(counts_by_source[spec.raw_source_key], spec)
        deg_path = raw_dir / spec.output_name
        deg.to_csv(deg_path, index=False)
        row_counts[spec.study_id] = int(len(deg))
        raw_input = raw_paths[spec.raw_source_key]
        inputs = [raw_input, raw_paths["hgnc"]] if spec.raw_source_key == "gse223011" else [raw_input]
        write_source_sidecar(
            deg_path,
            command,
            inputs=inputs,
            metadata={"generator": "indisulam-derived-logcpm-welch", "study_id": spec.study_id, **GUARDRAIL_FLAGS},
        )
        contrast_summaries.append(
            {
                "study_id": spec.study_id,
                "paper_id": spec.paper_id,
                "rows": int(len(deg)),
                "top_gene": str(deg.iloc[0]["gene_symbol"]) if not deg.empty else "",
                "anchor_hits_top100": sorted(set(deg.head(100)["gene_symbol"]).intersection({row[0] for row in INDISULAM_ANCHORS})),
                "output": str(deg_path),
            }
        )

    catalog = pd.DataFrame(_catalog_rows(specs, row_counts))
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(catalog_path, index=False)
    write_source_sidecar(catalog_path, command, inputs=[raw_paths["gse223011"], raw_paths["gse268568"], raw_paths["hgnc"]], metadata={"generator": "indisulam-derived-catalog", **GUARDRAIL_FLAGS})
    _write_gold_panel(gold_path, command)

    availability = pd.DataFrame(_availability_records(raw_paths))
    availability_tsv = output_dir / "indisulam_source_availability.tsv"
    availability_json = output_dir / "indisulam_source_availability.json"
    availability.to_csv(availability_tsv, sep="\t", index=False)
    availability_json.write_text(json.dumps(availability.to_dict(orient="records"), indent=2, sort_keys=True) + "\n")

    extra_metadata = {
        **GUARDRAIL_FLAGS,
        "topic": "Indisulam / E7070 RBM39 degrader response",
        "primary_assay_type": "RNA-seq count-derived DEG-like tables",
        "secondary_assay_types_allowed": ["microarray_limma_full_table", "author_full_deg_table"],
        "anchor_expected_direction_policy": "not_asserted",
    }
    extra_metadata_path = output_dir / "indisulam_score_metadata_extra.json"
    extra_metadata_path.write_text(json.dumps(extra_metadata, indent=2, sort_keys=True) + "\n")

    summary = {
        "corpus": "Indisulam derived-count pilot",
        **GUARDRAIL_FLAGS,
        "n_contrasts": len(specs),
        "n_source_units": int(catalog["paper_id"].nunique()),
        "mapping_summaries": mapping_summaries,
        "contrast_summaries": contrast_summaries,
        "catalog": str(catalog_path),
        "gold_panel": str(gold_path),
        "source_availability_tsv": str(availability_tsv),
        "score_metadata_extra": str(extra_metadata_path),
    }
    summary_path = output_dir / "indisulam_derived_deg_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    for artifact in (availability_tsv, availability_json, extra_metadata_path, summary_path):
        write_source_sidecar(
            artifact,
            command,
            inputs=[catalog_path, gold_path, raw_paths["gse223011"], raw_paths["gse268568"], raw_paths["hgnc"]],
            metadata={"generator": "indisulam-derived-summary", **GUARDRAIL_FLAGS},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/deg/raw/indisulam"))
    parser.add_argument("--catalog", type=Path, default=Path("data/studies/indisulam_derived_catalog.csv"))
    parser.add_argument("--gold", type=Path, default=Path("data/studies/gold/indisulam_anchor_panel.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/results/indisulam-pilot"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_indisulam_derived_deg.py"
        f" --raw-dir {args.raw_dir}"
        f" --catalog {args.catalog}"
        f" --gold {args.gold}"
        f" --output-dir {args.output_dir}"
    )
    summary = write_indisulam_derived_corpus(args.raw_dir, args.catalog, args.gold, args.output_dir, command=command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
