#!/usr/bin/env python
"""Prepare a locked heat-shock/HSF1 benchmark corpus for DEGORA.

The benchmark uses public RNA-seq count matrices with clear heat-shock versus
control labels. The locked gold panel is written before any DEGORA or baseline
score is produced.
"""

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

from degora.derived_counts import attach_low_count_filter_metadata, low_count_filter_mask, low_count_filter_summary
from degora.provenance import shell_command, write_source_sidecar


RAW_INPUTS = {
    "gse164834": {
        "path": "data/deg/raw/heat_shock/GSE164834_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
            "&acc=GSE164834&file=GSE164834_raw_counts_GRCh38.p13_NCBI.tsv.gz"
        ),
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164834",
    },
    "gse123980": {
        "path": "data/deg/raw/heat_shock/GSE123980_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
            "&acc=GSE123980&file=GSE123980_raw_counts_GRCh38.p13_NCBI.tsv.gz"
        ),
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123980",
    },
    "gse124609": {
        "path": "data/deg/raw/heat_shock/GSE124609_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
            "&acc=GSE124609&file=GSE124609_raw_counts_GRCh38.p13_NCBI.tsv.gz"
        ),
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124609",
    },
    "gse132447": {
        "path": "data/deg/raw/heat_shock/GSE132447_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
            "&acc=GSE132447&file=GSE132447_raw_counts_GRCh38.p13_NCBI.tsv.gz"
        ),
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132447",
    },
    "gse73471": {
        "path": "data/deg/raw/heat_shock/GSE73471_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
            "&acc=GSE73471&file=GSE73471_raw_counts_GRCh38.p13_NCBI.tsv.gz"
        ),
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73471",
    },
    "gse130493": {
        "path": "data/deg/raw/heat_shock/GSE130493_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130493/suppl/GSE130493_RAW.tar",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130493",
    },
    "gse57397": {
        "path": "data/deg/raw/heat_shock/GSE57397_raw_counts_GRCh38.p13_NCBI.tsv.gz",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
            "&acc=GSE57397&file=GSE57397_raw_counts_GRCh38.p13_NCBI.tsv.gz"
        ),
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57397",
    },
    "human_annot": {
        "path": "data/deg/raw/heat_shock/Human.GRCh38.p13.annot.tsv.gz",
        "url": "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz",
        "geo": "",
    },
}


HEAT_SHOCK_GOLD_GENES = [
    ("HSPA1A", "inducible_HSP70", "up"),
    ("HSPA1B", "inducible_HSP70", "up"),
    ("HSPA6", "inducible_HSP70", "up"),
    ("HSPH1", "HSP110_disaggregase", "up"),
    ("DNAJB1", "HSP40_cochaperone", "up"),
    ("HSPB1", "small_HSP", "up"),
    ("BAG3", "HSP70_cochaperone", "up"),
    ("HSPB8", "small_HSP", "up"),
    ("SERPINH1", "collagen_chaperone", "up"),
    ("HSPA4", "HSP70_family", "up"),
    ("HSPA4L", "HSP70_family", "up"),
    ("HSP90AA1", "HSP90_chaperone", "up"),
    ("HSP90AB1", "HSP90_chaperone", "up"),
    ("DNAJB6", "HSP40_cochaperone", "up"),
    ("DNAJC7", "HSP40_cochaperone", "up"),
    ("HSPA5", "stress_chaperone", "up"),
]


GUARDRAIL_FLAGS = {
    "benchmark_topic": "heat shock / HSF1 transcriptional response",
    "gold_panel_locked_before_scoring": True,
    "primary_min_studies": 2,
    "claim_scope": "gold_panel_prioritization_in_a_prelocked_heat_shock_hsf1_benchmark",
    "score_interpretation": "prioritization_index_not_calibrated_probability",
    "time_course_policy": (
        "Related time points or related phenotypes from one accession share a paper_id/source-unit and are "
        "mean-aggregated by the score; no max-signal source-unit selection is used."
    ),
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
    pipeline: str = "logCPM_Welch_derived_from_NCBI_GEO_counts"
    assay_type: str = "RNA-seq"
    source_input_type: str = "derived_count_table"
    platform: str = "NCBI_GEO_RNAseq_counts_GRCh38.p13"
    normalization: str = "logCPM_from_raw_counts"
    notes: str = ""


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


def _clean_symbol(values: pd.Series) -> pd.Series:
    return values.astype("string").str.strip().str.upper().replace({"": pd.NA, "NAN": pd.NA, "NONE": pd.NA})


def _collapse_by_symbol(counts: pd.DataFrame) -> pd.DataFrame:
    frame = counts.copy()
    symbols = _clean_symbol(pd.Series(frame.index, index=frame.index)).to_numpy()
    mask = pd.notna(symbols) & (symbols != "")
    frame = frame.loc[mask]
    frame.index = symbols[mask]
    return frame.groupby(frame.index).sum(numeric_only=True)


def _log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    libraries = numeric.sum(axis=0).replace(0.0, np.nan)
    cpm = numeric.divide(libraries, axis=1) * 1_000_000.0
    return np.log2(cpm.fillna(0.0) + 1.0)


def _read_ncbi_geneid_symbol_map(path: Path) -> dict[str, str]:
    annot = pd.read_csv(path, sep="\t", dtype=str)
    if "GeneType" in annot.columns:
        annot = annot.loc[annot["GeneType"].astype("string").str.lower().eq("protein-coding")].copy()
    annot = annot.dropna(subset=["GeneID", "Symbol"])
    annot["GeneID"] = annot["GeneID"].astype(str)
    annot["Symbol"] = _clean_symbol(annot["Symbol"])
    annot = annot.dropna(subset=["Symbol"]).drop_duplicates("GeneID", keep="first")
    return dict(zip(annot["GeneID"], annot["Symbol"], strict=False))


def _read_ncbi_geneid_counts(path: Path, geneid_to_symbol: dict[str, str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    raw = pd.read_csv(path, sep="\t", index_col=0)
    gene_ids = raw.index.astype("string").str.strip()
    symbols = pd.Series(gene_ids.map(geneid_to_symbol), index=raw.index)
    counts = raw.loc[symbols.notna()].copy()
    counts.index = symbols.loc[symbols.notna()].astype(str)
    collapsed = _collapse_by_symbol(counts)
    return collapsed, {
        "raw_rows": int(len(raw)),
        "mapped_rows": int(symbols.notna().sum()),
        "n_unmapped_genes": int(symbols.isna().sum()),
        "n_gene_symbols_after_collapse": int(len(collapsed)),
        "sample_columns": list(raw.columns),
    }


def _derive_logcpm_welch_contrast(counts: pd.DataFrame, spec: ContrastSpec) -> pd.DataFrame:
    columns = list(spec.control_columns + spec.treat_columns)
    missing = sorted(set(columns).difference(counts.columns))
    if missing:
        raise ValueError(f"{spec.study_id} missing columns: {missing}")
    selected_counts = counts[columns]
    expressed = low_count_filter_mask(selected_counts)
    filter_summary = low_count_filter_summary(selected_counts, expressed)
    log_counts = _log_cpm(selected_counts).loc[expressed]
    control = log_counts[list(spec.control_columns)]
    treat = log_counts[list(spec.treat_columns)]
    log2fc = treat.mean(axis=1) - control.mean(axis=1)
    test = stats.ttest_ind(treat.to_numpy(dtype=float), control.to_numpy(dtype=float), axis=1, equal_var=False, nan_policy="omit")
    pvalue = pd.Series(test.pvalue, index=log_counts.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    out = pd.DataFrame({"gene_symbol": log_counts.index.astype(str), "log2FoldChange": log2fc, "pvalue": pvalue})
    out["padj"] = _bh_adjust(out["pvalue"])
    out["source_input_type"] = spec.source_input_type
    out = attach_low_count_filter_metadata(out, filter_summary)
    return out.sort_values(["pvalue", "gene_symbol"]).reset_index(drop=True)


def _contrast_specs(raw_dir: Path) -> list[ContrastSpec]:
    common_note = (
        "Derived-count heat-shock benchmark contrast; log2FC and p-values computed from NCBI GEO raw "
        "gene counts using logCPM and Welch t-test. Source-unit collapse is by paper_id."
    )
    return [
        ContrastSpec(
            study_id="HS_GSE164834_SW620_HS_1h_vs_Ctrl",
            paper_id="GSE164834_SW620_HS",
            source_url=RAW_INPUTS["gse164834"]["geo"],
            source_path=raw_dir / "GSE164834_SW620_HS_1h_vs_Ctrl_logcpm_welch.csv",
            raw_source_key="gse164834",
            species="Homo sapiens",
            cell_system="SW620 colon cancer cells, scrambled siRNA",
            condition="42C heat shock 1h vs 37C control",
            duration_h="1",
            n_ctrl=3,
            n_treat=3,
            control_columns=("GSM5020463", "GSM5020464", "GSM5020465"),
            treat_columns=("GSM5020469", "GSM5020470", "GSM5020471"),
            output_name="GSE164834_SW620_HS_1h_vs_Ctrl_logcpm_welch.csv",
            notes=common_note,
        ),
        ContrastSpec(
            study_id="HS_GSE164834_SW620_HS_1h_REC_1h_vs_Ctrl",
            paper_id="GSE164834_SW620_HS",
            source_url=RAW_INPUTS["gse164834"]["geo"],
            source_path=raw_dir / "GSE164834_SW620_HS_1h_REC_1h_vs_Ctrl_logcpm_welch.csv",
            raw_source_key="gse164834",
            species="Homo sapiens",
            cell_system="SW620 colon cancer cells, scrambled siRNA",
            condition="42C heat shock 1h plus 37C recovery 1h vs 37C control",
            duration_h="2",
            n_ctrl=3,
            n_treat=3,
            control_columns=("GSM5020463", "GSM5020464", "GSM5020465"),
            treat_columns=("GSM5020475", "GSM5020476", "GSM5020477"),
            output_name="GSE164834_SW620_HS_1h_REC_1h_vs_Ctrl_logcpm_welch.csv",
            notes=common_note,
        ),
        ContrastSpec(
            study_id="HS_GSE124609_WI38_Young_HS_vs_Ctrl",
            paper_id="GSE124609_WI38_HS",
            source_url=RAW_INPUTS["gse124609"]["geo"],
            source_path=raw_dir / "GSE124609_WI38_Young_HS_vs_Ctrl_logcpm_welch.csv",
            raw_source_key="gse124609",
            species="Homo sapiens",
            cell_system="WI38 young human lung fibroblasts",
            condition="heat shock vs control",
            duration_h="not_reported_in_matrix",
            n_ctrl=2,
            n_treat=2,
            control_columns=("GSM3537582", "GSM3537583"),
            treat_columns=("GSM3537584", "GSM3537585"),
            output_name="GSE124609_WI38_Young_HS_vs_Ctrl_logcpm_welch.csv",
            notes=common_note,
        ),
        ContrastSpec(
            study_id="HS_GSE124609_WI38_Senescent_HS_vs_Ctrl",
            paper_id="GSE124609_WI38_HS",
            source_url=RAW_INPUTS["gse124609"]["geo"],
            source_path=raw_dir / "GSE124609_WI38_Senescent_HS_vs_Ctrl_logcpm_welch.csv",
            raw_source_key="gse124609",
            species="Homo sapiens",
            cell_system="WI38 senescent human lung fibroblasts",
            condition="heat shock vs control",
            duration_h="not_reported_in_matrix",
            n_ctrl=2,
            n_treat=2,
            control_columns=("GSM3537586", "GSM3537587"),
            treat_columns=("GSM3537588", "GSM3537589"),
            output_name="GSE124609_WI38_Senescent_HS_vs_Ctrl_logcpm_welch.csv",
            notes=common_note,
        ),
        ContrastSpec(
            study_id="HS_GSE132447_PrimaryNeurons_HS_3h_REC_vs_Untreated",
            paper_id="GSE132447_PrimaryNeurons_HS",
            source_url=RAW_INPUTS["gse132447"]["geo"],
            source_path=raw_dir / "GSE132447_PrimaryNeurons_HS_3h_REC_vs_Untreated_logcpm_welch.csv",
            raw_source_key="gse132447",
            species="Homo sapiens",
            cell_system="primary neurons",
            condition="heat shock with 3h recovery vs untreated",
            duration_h="3h_recovery",
            n_ctrl=3,
            n_treat=3,
            control_columns=("GSM3864912", "GSM3864913", "GSM3864914"),
            treat_columns=("GSM3864918", "GSM3864919", "GSM3864920"),
            output_name="GSE132447_PrimaryNeurons_HS_3h_REC_vs_Untreated_logcpm_welch.csv",
            notes=(
                common_note
                + " GSE132447 1h-recovery replicates are incomplete in the NCBI generated count file, so the primary benchmark uses the complete 3h-recovery arm."
            ),
        ),
    ]


def _catalog_rows(specs: list[ContrastSpec], row_counts: dict[str, int]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for spec in specs:
        rows.append(
            {
                "study_id": spec.study_id,
                "paper_id": spec.paper_id,
                "source_url": spec.source_url,
                "source_path": str(spec.source_path),
                "pipeline": spec.pipeline,
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
                "assay_type": spec.assay_type,
                "source_input_type": spec.source_input_type,
                "platform": spec.platform,
                "normalization": spec.normalization,
                "probe_id_column": "",
                "probe_collapse": "",
                "table_scope": "full_results",
                "rank_universe_size": row_counts.get(spec.study_id, ""),
                "include_in_analysis": True,
                "notes": spec.notes,
            }
        )
    return rows


def _write_gold_panel(path: Path, command: str) -> None:
    rows = [
        {
            "gene_symbol": symbol,
            "expected_direction": direction,
            "evidence_class": evidence_class,
            "locked": "yes",
            "evidence_basis": (
                "Pre-output compact heat-shock/HSF1 sentinel panel: canonical heat-inducible HSP70/HSP90, "
                "HSP40/co-chaperone, and small-HSP genes supported by Reactome cellular-response-to-heat-stress "
                "membership and HSF1 target literature."
            ),
            "notes": "Locked before heat-shock DEGORA or baseline scoring; used only for recall/prioritization diagnostics.",
        }
        for symbol, evidence_class, direction in HEAT_SHOCK_GOLD_GENES
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    write_source_sidecar(path, command, metadata={"generator": "heat-shock-hsf1-gold-panel", **GUARDRAIL_FLAGS})


def _availability_records(raw_paths: dict[str, Path], summaries: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    def file_record(key: str, status: str, use_decision: str, reason: str) -> dict[str, Any]:
        path = raw_paths[key]
        return {
            "source_id": key.upper(),
            "status": status,
            "use_decision": use_decision,
            "local_path": str(path),
            "local_file_exists": bool(path.exists()),
            "local_size_bytes": path.stat().st_size if path.exists() else 0,
            "sha256": _sha256(path) if path.exists() and path.is_file() else "",
            "source_url": RAW_INPUTS[key]["geo"] or RAW_INPUTS[key]["url"],
            "reason": reason,
            "summary": json.dumps(summaries.get(key, {}), sort_keys=True),
        }

    return [
        file_record("gse164834", "active_primary_source", "included", "SW620 42C heat-shock and recovery arms have 3 control and 3 treatment replicates"),
        file_record("gse124609", "active_primary_source", "included", "WI38 young and senescent RNA-seq heat-shock contrasts have two replicates per arm"),
        file_record("gse132447", "active_primary_source", "included", "primary neuron 3h-recovery RNA-seq arm has complete 3 vs 3 count columns"),
        file_record(
            "gse123980",
            "deferred_candidate",
            "excluded_metadata_ambiguous_for_primary",
            "K562 total RNA-seq sample titles indicate HS15/HS30 but treatment characteristics for total RNA rows do not explicitly encode heat exposure; nascent TT-seq rows are not a primary expression table",
        ),
        file_record("gse73471", "deferred_candidate", "excluded_low_control_replication", "WI38 heat-shock source has only one untreated control in the available matrix"),
        file_record("gse130493", "deferred_candidate", "excluded_no_matched_control", "downloaded RAW archive contains heat-shock samples but no matched control samples"),
        file_record("gse57397", "deferred_candidate", "excluded_no_untreated_control", "HCT116 source contains heat-shock and recovery samples but no matched untreated/non-heat control arm"),
    ]


def _write_markdown_summary(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# Heat Shock / HSF1 Benchmark Preparation",
        "",
        "Gold panel was locked before scoring.",
        "",
        f"- Active contrasts: {summary['n_contrasts']}",
        f"- Independent source units: {summary['n_source_units']}",
        f"- Gold genes: {', '.join(summary['gold_genes'])}",
        "",
        "## Inclusion Decisions",
        "",
        "| source | status | decision | reason |",
        "| --- | --- | --- | --- |",
    ]
    for row in summary["source_availability"]:
        lines.append(f"| {row['source_id']} | {row['status']} | {row['use_decision']} | {row['reason']} |")
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- All active inputs are full NCBI GEO RNA-seq count matrices, not DEG-only lists.",
            "- Same-accession time-point or phenotype rows share a source unit.",
            "- The locked panel asserts expected up-regulation under heat shock.",
            "- The score is a prioritization index, not a calibrated probability.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def write_heat_shock_benchmark(
    raw_dir: Path,
    catalog_path: Path,
    gold_path: Path,
    output_dir: Path,
    *,
    command: str,
) -> dict[str, Any]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_paths = {key: Path(info["path"]) for key, info in RAW_INPUTS.items()}
    required_keys = ["gse164834", "gse124609", "gse132447", "human_annot"]
    optional_keys = ["gse123980", "gse73471", "gse130493", "gse57397"]
    for key in required_keys:
        info = RAW_INPUTS[key]
        path = Path(info["path"])
        if not path.exists():
            raise FileNotFoundError(f"Missing required raw input: {path}")
        write_source_sidecar(
            path,
            shell_command(["curl", "-L", "-o", path, info["url"]]),
            inputs=[],
            metadata={"generator": "manual-download-source", "source_url": info["url"]},
        )
    for key in optional_keys:
        info = RAW_INPUTS[key]
        path = Path(info["path"])
        if path.exists():
            write_source_sidecar(
                path,
                shell_command(["curl", "-L", "-o", path, info["url"]]),
                inputs=[],
                metadata={"generator": "manual-download-source", "source_url": info["url"]},
            )

    _write_gold_panel(gold_path, command)
    geneid_to_symbol = _read_ncbi_geneid_symbol_map(raw_paths["human_annot"])
    counts_by_source: dict[str, pd.DataFrame] = {}
    source_summaries: dict[str, dict[str, Any]] = {}
    for key in ["gse164834", "gse124609", "gse132447"]:
        counts, summary = _read_ncbi_geneid_counts(raw_paths[key], geneid_to_symbol)
        counts_by_source[key] = counts
        source_summaries[key] = summary

    specs = _contrast_specs(raw_dir)
    contrast_summaries: list[dict[str, Any]] = []
    row_counts: dict[str, int] = {}
    gold_symbols = {gene for gene, _, _ in HEAT_SHOCK_GOLD_GENES}
    for spec in specs:
        deg = _derive_logcpm_welch_contrast(counts_by_source[spec.raw_source_key], spec)
        deg_path = raw_dir / spec.output_name
        deg.to_csv(deg_path, index=False)
        row_counts[spec.study_id] = int(len(deg))
        write_source_sidecar(
            deg_path,
            command,
            inputs=[raw_paths[spec.raw_source_key], raw_paths["human_annot"]],
            metadata={"generator": "heat-shock-logcpm-welch", "study_id": spec.study_id, **GUARDRAIL_FLAGS},
        )
        top100_gold_hits = sorted(set(deg.head(100)["gene_symbol"]).intersection(gold_symbols))
        contrast_summaries.append(
            {
                "study_id": spec.study_id,
                "paper_id": spec.paper_id,
                "rows": int(len(deg)),
                "top_gene": str(deg.iloc[0]["gene_symbol"]) if not deg.empty else "",
                "gold_hits_top100": top100_gold_hits,
                "output": str(deg_path),
                "source_input_type": spec.source_input_type,
            }
        )

    catalog = pd.DataFrame(_catalog_rows(specs, row_counts))
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(catalog_path, index=False)
    catalog_inputs = [raw_paths[key] for key in required_keys] + [gold_path]
    write_source_sidecar(
        catalog_path,
        command,
        inputs=catalog_inputs,
        metadata={"generator": "heat-shock-hsf1-catalog", **GUARDRAIL_FLAGS},
    )

    availability = _availability_records(raw_paths, source_summaries)
    availability_tsv = output_dir / "heat_shock_source_availability.tsv"
    availability_json = output_dir / "heat_shock_source_availability.json"
    pd.DataFrame(availability).to_csv(availability_tsv, sep="\t", index=False)
    availability_json.write_text(json.dumps(availability, indent=2, sort_keys=True) + "\n")

    extra_metadata = {
        **GUARDRAIL_FLAGS,
        "gold_panel": str(gold_path),
        "active_corpus": "GSE164834 + GSE124609 + GSE132447 NCBI GEO RNA-seq count-derived heat-shock contrasts",
        "deferred_sources": {
            "GSE123980": "total RNA heat-shock metadata ambiguous for primary benchmark; nascent TT-seq not used as primary expression input",
            "GSE73471": "one untreated control only",
            "GSE130493": "heat-shock-only archive without matched controls",
            "GSE57397": "heat-shock and recovery only; no untreated/non-heat control arm",
        },
        "gold_external_basis": [
            "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS.html",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC6888953/",
        ],
    }
    extra_metadata_path = output_dir / "heat_shock_score_metadata_extra.json"
    extra_metadata_path.write_text(json.dumps(extra_metadata, indent=2, sort_keys=True) + "\n")

    summary: dict[str, Any] = {
        "corpus": "heat shock / HSF1 response benchmark",
        **GUARDRAIL_FLAGS,
        "n_contrasts": len(specs),
        "n_source_units": int(catalog["paper_id"].nunique()),
        "gold_genes": [gene for gene, _, _ in HEAT_SHOCK_GOLD_GENES],
        "catalog": str(catalog_path),
        "gold_panel": str(gold_path),
        "contrast_summaries": contrast_summaries,
        "source_summaries": source_summaries,
        "source_availability": availability,
        "source_availability_tsv": str(availability_tsv),
        "score_metadata_extra": str(extra_metadata_path),
    }
    summary_json = output_dir / "heat_shock_benchmark_preparation_summary.json"
    summary_md = output_dir / "heat_shock_benchmark_preparation_summary.md"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    _write_markdown_summary(summary, summary_md)

    for artifact in (availability_tsv, availability_json, extra_metadata_path, summary_json, summary_md):
        write_source_sidecar(
            artifact,
            command,
            inputs=[catalog_path, gold_path, *catalog_inputs],
            metadata={"generator": "heat-shock-hsf1-preparation-summary", **GUARDRAIL_FLAGS},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/deg/raw/heat_shock"))
    parser.add_argument("--catalog", type=Path, default=Path("data/studies/heat_shock_hsf1_catalog.csv"))
    parser.add_argument("--gold", type=Path, default=Path("data/studies/gold/heat_shock_hsf1_gold_panel.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/results/heat-shock-benchmark"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = shell_command(
        [
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_heat_shock_benchmark.py",
            "--raw-dir",
            args.raw_dir,
            "--catalog",
            args.catalog,
            "--gold",
            args.gold,
            "--output-dir",
            args.output_dir,
        ]
    )
    summary = write_heat_shock_benchmark(args.raw_dir, args.catalog, args.gold, args.output_dir, command=command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
