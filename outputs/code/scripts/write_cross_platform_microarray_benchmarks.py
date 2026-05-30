#!/usr/bin/env python
"""Build ER-stress and IFN RNA-seq+microarray DEGORA benchmarks."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from io import StringIO
from pathlib import Path
from typing import Any
from urllib.request import urlretrieve

import pandas as pd

from derive_microarray_deg import derive_microarray_deg
from degora.baselines import write_baseline_outputs
from degora.provenance import shell_command, write_source_sidecar
from degora.score_db import write_score_database
from degora.slice_runner import CATALOG_COLUMNS, run_slice
from write_gold_comparator_summary import build_summary, write_markdown


RAW_DIR = Path("data/deg/raw/microarray_cross_platform")
RESULT_ROOT = Path("outputs/results/cross-platform-microarray")
HARMONIZED_DIR = Path("data/deg/harmonized")

GSE71634_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71634/matrix/GSE71634_series_matrix.txt.gz"
GSE19519_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19519/matrix/GSE19519_series_matrix.txt.gz"
GSE3045_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3045/matrix/GSE3045_series_matrix.txt.gz"
GSE22282_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22282/matrix/GSE22282_series_matrix.txt.gz"
GPL10558_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
GPL570_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

GEO_URLS = {
    "GSE71634": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71634",
    "GSE19519": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19519",
    "GSE3045": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3045",
    "GSE22282": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22282",
}

TOPICS = {
    "IFN": {
        "rna_only_result": Path("outputs/results/ifn-pilot"),
        "rna_only_summary": Path("outputs/results/ifn-pilot/ifn_gold_comparator_summary.csv"),
        "mixed_result": Path("outputs/results/ifn-cross-platform"),
        "mixed_summary_name": "ifn_cross_platform_gold_comparator_summary",
        "gold": Path("data/studies/gold/ifn_gold_panel.csv"),
        "base_catalog": Path("data/studies/ifn_derived_catalog.csv"),
        "mixed_catalog": Path("data/studies/ifn_cross_platform_catalog.csv"),
    },
    "ER stress": {
        "rna_only_result": Path("outputs/results/er-stress-benchmark-primary"),
        "rna_only_summary": Path("outputs/results/er-stress-benchmark-primary/er_stress_primary_gold_comparator_summary.csv"),
        "mixed_result": Path("outputs/results/er-stress-cross-platform"),
        "mixed_summary_name": "er_stress_cross_platform_gold_comparator_summary",
        "gold": Path("data/studies/gold/er_stress_upr_gold_panel.csv"),
        "base_catalog": Path("data/studies/er_stress_upr_primary_catalog.csv"),
        "mixed_catalog": Path("data/studies/er_stress_upr_cross_platform_catalog.csv"),
    },
    "Hypoxia": {
        "rna_only_result": Path("outputs/results/hypoxia-hif1-benchmark"),
        "rna_only_summary": Path("outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv"),
        "mixed_result": Path("outputs/results/hypoxia-cross-platform"),
        "mixed_summary_name": "hypoxia_cross_platform_gold_comparator_summary",
        "gold": Path("data/studies/gold/hypoxia_hif1_gold_panel.csv"),
        "base_catalog": Path("data/studies/hypoxia_catalog.csv"),
        "mixed_catalog": Path("data/studies/hypoxia_cross_platform_catalog.csv"),
    },
}


def _download(path: Path, url: str, *, command: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.stat().st_size == 0:
        urlretrieve(url, path)
    write_source_sidecar(path, f"curl -L -o {path} {url}", metadata={"generator": "geo-download", "source_url": url})


def _read_text_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read().splitlines()


def _strip_quotes(value: str) -> str:
    return value.strip().strip('"')


def _series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    lines = _read_text_lines(path)
    metadata: dict[str, list[str]] = {}
    table_start = lines.index("!series_matrix_table_begin")
    table_end = lines.index("!series_matrix_table_end")
    for line in lines[:table_start]:
        if not line.startswith("!Sample_"):
            continue
        fields = line.split("\t")
        metadata[fields[0]] = [_strip_quotes(value) for value in fields[1:]]
    expression = pd.read_csv(StringIO("\n".join(lines[table_start + 1 : table_end])), sep="\t", dtype=str)
    expression = expression.rename(columns={"ID_REF": "probe_id"})
    accessions = metadata.get("!Sample_geo_accession", list(expression.columns[1:]))
    titles = metadata.get("!Sample_title", [""] * len(accessions))
    samples = pd.DataFrame({"sample": accessions, "title": titles})
    return expression, samples


def _annotation(path: Path) -> pd.DataFrame:
    lines = _read_text_lines(path)
    data_lines = [line for line in lines if line and not line.startswith(("#", "^", "!"))]
    table = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", dtype=str, low_memory=False)
    symbol_column = None
    for candidate in ["Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "ILMN_Gene"]:
        if candidate in table.columns:
            symbol_column = candidate
            break
    if symbol_column is None:
        raise ValueError(f"Could not find gene-symbol column in {path}; columns={list(table.columns)[:20]}")
    out = table.loc[:, ["ID", symbol_column]].rename(columns={"ID": "probe_id", symbol_column: "gene_symbol"})
    out["gene_symbol"] = out["gene_symbol"].map(_first_symbol)
    out = out.loc[out["gene_symbol"].ne("")].drop_duplicates("probe_id")
    return out


def _first_symbol(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "null", "---"}:
        return ""
    for part in re.split(r"///|//|;|,", text):
        symbol = part.strip().upper()
        if symbol and symbol not in {"---", "NA", "N/A", "NULL"}:
            return symbol
    return ""


def _write_gene_matrix(
    expression: pd.DataFrame,
    samples: pd.DataFrame,
    annotation: pd.DataFrame,
    output: Path,
    *,
    command: str,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    merged = expression.merge(annotation, on="probe_id", how="inner")
    ordered = ["gene_symbol", "probe_id", *samples["sample"].tolist()]
    merged = merged.loc[:, ordered].copy()
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    write_source_sidecar(output, command, inputs=[], metadata={"generator": "geo-series-matrix-symbol-merge", **metadata})
    return {
        "output": str(output),
        "n_probe_rows_with_symbol": int(len(merged)),
        "n_gene_symbols": int(merged["gene_symbol"].nunique()),
    }


def _sample_names(samples: pd.DataFrame, predicate) -> list[str]:
    subset = samples.loc[samples["title"].map(lambda value: predicate(str(value).lower()))]
    return subset["sample"].astype(str).tolist()


def _derive_microarray_inputs(command: str) -> dict[str, Any]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    gse71634_matrix = RAW_DIR / "GSE71634_series_matrix.txt.gz"
    gse19519_matrix = RAW_DIR / "GSE19519_series_matrix.txt.gz"
    gse3045_matrix = RAW_DIR / "GSE3045_series_matrix.txt.gz"
    gse22282_matrix = RAW_DIR / "GSE22282_series_matrix.txt.gz"
    gpl10558 = RAW_DIR / "GPL10558.annot.gz"
    gpl570 = RAW_DIR / "GPL570.annot.gz"
    for path, url in [
        (gse71634_matrix, GSE71634_MATRIX_URL),
        (gse19519_matrix, GSE19519_MATRIX_URL),
        (gse3045_matrix, GSE3045_MATRIX_URL),
        (gse22282_matrix, GSE22282_MATRIX_URL),
        (gpl10558, GPL10558_ANNOT_URL),
        (gpl570, GPL570_ANNOT_URL),
    ]:
        _download(path, url, command=command)

    ifn_expression, ifn_samples = _series_matrix(gse71634_matrix)
    er_expression, er_samples = _series_matrix(gse19519_matrix)
    hyp_astro_expression, hyp_astro_samples = _series_matrix(gse3045_matrix)
    hyp_dc_expression, hyp_dc_samples = _series_matrix(gse22282_matrix)
    ifn_matrix = RAW_DIR / "GSE71634_IFNB_PBMC_probe_gene_matrix.csv"
    er_matrix = RAW_DIR / "GSE19519_ER_STRESS_LCL_probe_gene_matrix.csv"
    hyp_astro_matrix = RAW_DIR / "GSE3045_hypoxia_astrocyte_probe_gene_matrix.csv"
    hyp_dc_matrix = RAW_DIR / "GSE22282_hypoxia_mDC_probe_gene_matrix.csv"
    ifn_matrix_summary = _write_gene_matrix(
        ifn_expression,
        ifn_samples,
        _annotation(gpl10558),
        ifn_matrix,
        command=command,
        metadata={"accession": "GSE71634", "platform": "GPL10558", "source_url": GEO_URLS["GSE71634"]},
    )
    er_matrix_summary = _write_gene_matrix(
        er_expression,
        er_samples,
        _annotation(gpl570),
        er_matrix,
        command=command,
        metadata={"accession": "GSE19519", "platform": "GPL570", "source_url": GEO_URLS["GSE19519"]},
    )
    hyp_astro_matrix_summary = _write_gene_matrix(
        hyp_astro_expression,
        hyp_astro_samples,
        _annotation(gpl570),
        hyp_astro_matrix,
        command=command,
        metadata={"accession": "GSE3045", "platform": "GPL570", "source_url": GEO_URLS["GSE3045"]},
    )
    hyp_dc_matrix_summary = _write_gene_matrix(
        hyp_dc_expression,
        hyp_dc_samples,
        _annotation(gpl570),
        hyp_dc_matrix,
        command=command,
        metadata={"accession": "GSE22282", "platform": "GPL570", "source_url": GEO_URLS["GSE22282"]},
    )

    ifn_ctrl = _sample_names(ifn_samples, lambda title: "unstimulated" in title)
    ifn_treat = _sample_names(ifn_samples, lambda title: "stimulated" in title and "unstimulated" not in title)
    er_tun_ctrl = _sample_names(er_samples, lambda title: "60grandparents" in title and "dmso" in title and "8_hours" in title)
    er_tun_treat = _sample_names(er_samples, lambda title: "60grandparents" in title and "tunicamycin" in title and "8_hours" in title)
    er_tg_ctrl = _sample_names(er_samples, lambda title: "twins" in title and "dmso" in title and "4h" in title)
    er_tg_treat = _sample_names(er_samples, lambda title: "twins" in title and "thapsigargin" in title and "4h" in title)
    hyp_astro_ctrl = _sample_names(hyp_astro_samples, lambda title: "normoxia" in title)
    hyp_astro_treat = _sample_names(hyp_astro_samples, lambda title: "hypoxia" in title)
    hyp_dc_ctrl = _sample_names(hyp_dc_samples, lambda title: "normoxia" in title)
    hyp_dc_treat = _sample_names(hyp_dc_samples, lambda title: "hypoxia" in title)

    derived: dict[str, Any] = {
        "downloaded_inputs": [str(gse71634_matrix), str(gse19519_matrix), str(gse3045_matrix), str(gse22282_matrix), str(gpl10558), str(gpl570)],
        "ifn_matrix": ifn_matrix_summary,
        "er_matrix": er_matrix_summary,
        "hypoxia_astrocyte_matrix": hyp_astro_matrix_summary,
        "hypoxia_mdc_matrix": hyp_dc_matrix_summary,
        "sample_groups": {
            "ifn_ctrl": len(ifn_ctrl),
            "ifn_treat": len(ifn_treat),
            "er_tunicamycin_ctrl": len(er_tun_ctrl),
            "er_tunicamycin_treat": len(er_tun_treat),
            "er_thapsigargin_ctrl": len(er_tg_ctrl),
            "er_thapsigargin_treat": len(er_tg_treat),
            "hypoxia_astrocyte_ctrl": len(hyp_astro_ctrl),
            "hypoxia_astrocyte_treat": len(hyp_astro_treat),
            "hypoxia_mdc_ctrl": len(hyp_dc_ctrl),
            "hypoxia_mdc_treat": len(hyp_dc_treat),
        },
    }

    specs = [
        (
            ifn_matrix,
            RAW_DIR / "GSE71634_PBMC_IFNB_18h_vs_unstimulated_microarray_welch.csv",
            "gene_symbol",
            "probe_id",
            ifn_ctrl,
            ifn_treat,
            True,
            RESULT_ROOT / "gse71634_ifnb_microarray_summary.json",
            {
                "study_id": "IFN_GSE71634_PBMC_IFNB_18h",
                "paper_id": "GSE71634_PBMC_IFNB",
                "source_url": GEO_URLS["GSE71634"],
                "platform": "GPL10558",
                "normalization": "GEO_series_matrix_log2_transform",
            },
        ),
        (
            er_matrix,
            RAW_DIR / "GSE19519_LCL_tunicamycin_8h_vs_DMSO_microarray_welch.csv",
            "gene_symbol",
            "probe_id",
            er_tun_ctrl,
            er_tun_treat,
            False,
            RESULT_ROOT / "gse19519_tunicamycin_microarray_summary.json",
            {
                "study_id": "ER_GSE19519_LCL_TUNICAMYCIN_8h",
                "paper_id": "GSE19519_ER_STRESS",
                "source_url": GEO_URLS["GSE19519"],
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_log2_processed",
            },
        ),
        (
            er_matrix,
            RAW_DIR / "GSE19519_LCL_thapsigargin_4h_vs_DMSO_microarray_welch.csv",
            "gene_symbol",
            "probe_id",
            er_tg_ctrl,
            er_tg_treat,
            False,
            RESULT_ROOT / "gse19519_thapsigargin_microarray_summary.json",
            {
                "study_id": "ER_GSE19519_LCL_THAPSIGARGIN_4h",
                "paper_id": "GSE19519_ER_STRESS",
                "source_url": GEO_URLS["GSE19519"],
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_log2_processed",
            },
        ),
        (
            hyp_astro_matrix,
            RAW_DIR / "GSE3045_primary_astrocyte_hypoxia_24h_vs_normoxia_microarray_welch.csv",
            "gene_symbol",
            "probe_id",
            hyp_astro_ctrl,
            hyp_astro_treat,
            True,
            RESULT_ROOT / "gse3045_hypoxia_astrocyte_microarray_summary.json",
            {
                "study_id": "HYP_GSE3045_ASTROCYTE_HYPOXIA_24h",
                "paper_id": "GSE3045_ASTROCYTE_HYPOXIA",
                "source_url": GEO_URLS["GSE3045"],
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_MAS5_log2_transform",
            },
        ),
        (
            hyp_dc_matrix,
            RAW_DIR / "GSE22282_mDC_hypoxia_vs_normoxia_microarray_welch.csv",
            "gene_symbol",
            "probe_id",
            hyp_dc_ctrl,
            hyp_dc_treat,
            True,
            RESULT_ROOT / "gse22282_hypoxia_mdc_microarray_summary.json",
            {
                "study_id": "HYP_GSE22282_MDC_HYPOXIA",
                "paper_id": "GSE22282_MDC_HYPOXIA",
                "source_url": GEO_URLS["GSE22282"],
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_MAS5_log2_transform",
            },
        ),
    ]
    for matrix, output, gene_col, probe_col, ctrl, treat, log2_transform, summary_path, metadata in specs:
        summary = derive_microarray_deg(
            matrix,
            output,
            gene_column=gene_col,
            probe_column=probe_col,
            control_samples=ctrl,
            treatment_samples=treat,
            log2_transform=log2_transform,
            command=command,
            summary_path=summary_path,
            metadata=metadata,
        )
        derived[metadata["study_id"]] = summary
    return derived


def _append_catalog(base_catalog: Path, rows: list[dict[str, Any]], output_catalog: Path, command: str) -> None:
    base = pd.read_csv(base_catalog)
    catalog = pd.concat([base, pd.DataFrame(rows)], ignore_index=True, sort=False)
    for column in CATALOG_COLUMNS:
        if column not in catalog.columns:
            catalog[column] = ""
    catalog = catalog[CATALOG_COLUMNS]
    output_catalog.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(output_catalog, index=False)
    write_source_sidecar(output_catalog, command, inputs=[base_catalog, *[Path(row["source_path"]) for row in rows]], metadata={"generator": "cross-platform-catalog"})


def _catalog_rows() -> dict[str, list[dict[str, Any]]]:
    common = {
        "pipeline": "welch_microarray_normalized_matrix",
        "species": "Homo sapiens",
        "gene_column": "gene_symbol",
        "lfc_column": "log2FoldChange",
        "p_column": "pvalue",
        "padj_column": "padj",
        "sep": "",
        "sheet_name": "",
        "gene_type_column": "",
        "gene_type_keep": "",
        "assay_type": "microarray",
        "source_input_type": "normalized_expression_matrix",
        "probe_id_column": "probe_id",
        "probe_collapse": "min_pvalue_max_abs_lfc",
        "time_course_mode": "mean",
        "temporal_mode": "",
        "table_scope": "full_results",
        "rank_universe_size": "",
        "include_in_analysis": True,
    }
    return {
        "IFN": [
            {
                **common,
                "study_id": "IFN_GSE71634_PBMC_IFNB_18h",
                "paper_id": "GSE71634_PBMC_IFNB",
                "source_url": GEO_URLS["GSE71634"],
                "source_path": str(RAW_DIR / "GSE71634_PBMC_IFNB_18h_vs_unstimulated_microarray_welch.csv"),
                "cell_system": "PBMC healthy donors",
                "hypoxia_modality": "IFN-beta vs unstimulated",
                "duration_h": "18",
                "n_ctrl": 20,
                "n_treat": 20,
                "platform": "GPL10558",
                "normalization": "GEO_series_matrix_log2_transform",
                "notes": "Cross-platform IFN microarray source added after RNA-seq-only benchmark lock.",
            }
        ],
        "ER stress": [
            {
                **common,
                "study_id": "ER_GSE19519_LCL_TUNICAMYCIN_8h",
                "paper_id": "GSE19519_ER_STRESS",
                "source_url": GEO_URLS["GSE19519"],
                "source_path": str(RAW_DIR / "GSE19519_LCL_tunicamycin_8h_vs_DMSO_microarray_welch.csv"),
                "cell_system": "immortalized B cells / LCL",
                "hypoxia_modality": "tunicamycin vs DMSO",
                "duration_h": "8",
                "n_ctrl": 60,
                "n_treat": 60,
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_log2_processed",
                "notes": "GSE19519 ER-stress microarray source; shares source unit with thapsigargin row. Derived from a public normalized matrix with a Welch fallback that treats samples as independent; family/pairing structure is not modeled and is downweighted through source_input_type=normalized_expression_matrix.",
            },
            {
                **common,
                "study_id": "ER_GSE19519_LCL_THAPSIGARGIN_4h",
                "paper_id": "GSE19519_ER_STRESS",
                "source_url": GEO_URLS["GSE19519"],
                "source_path": str(RAW_DIR / "GSE19519_LCL_thapsigargin_4h_vs_DMSO_microarray_welch.csv"),
                "cell_system": "immortalized B cells / LCL",
                "hypoxia_modality": "thapsigargin vs DMSO",
                "duration_h": "4",
                "n_ctrl": 52,
                "n_treat": 52,
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_log2_processed",
                "notes": "GSE19519 ER-stress microarray source; shares source unit with tunicamycin row. Derived from a public normalized matrix with a Welch fallback that treats samples as independent; family/pairing structure is not modeled and is downweighted through source_input_type=normalized_expression_matrix.",
            },
        ],
        "Hypoxia": [
            {
                **common,
                "study_id": "HYP_GSE3045_ASTROCYTE_HYPOXIA_24h",
                "paper_id": "GSE3045_ASTROCYTE_HYPOXIA",
                "source_url": GEO_URLS["GSE3045"],
                "source_path": str(RAW_DIR / "GSE3045_primary_astrocyte_hypoxia_24h_vs_normoxia_microarray_welch.csv"),
                "cell_system": "primary human astrocytes",
                "hypoxia_modality": "1% O2 hypoxia 24h vs normoxia",
                "duration_h": "24",
                "n_ctrl": 3,
                "n_treat": 3,
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_MAS5_log2_transform",
                "notes": "Human hypoxia microarray source added after RNA-seq-only benchmark lock; GEO design reports three independent normoxia and three hypoxia astrocyte batches.",
            },
            {
                **common,
                "study_id": "HYP_GSE22282_MDC_HYPOXIA",
                "paper_id": "GSE22282_MDC_HYPOXIA",
                "source_url": GEO_URLS["GSE22282"],
                "source_path": str(RAW_DIR / "GSE22282_mDC_hypoxia_vs_normoxia_microarray_welch.csv"),
                "cell_system": "human mature dendritic cells",
                "hypoxia_modality": "hypoxia vs normoxia",
                "duration_h": "not_reported_in_matrix",
                "n_ctrl": 3,
                "n_treat": 3,
                "platform": "GPL570",
                "normalization": "GEO_series_matrix_MAS5_log2_transform",
                "notes": "Human hypoxia microarray source added after RNA-seq-only benchmark lock; GEO design reports three healthy-donor biological replicates under normoxia and hypoxia.",
            },
        ],
    }


def _run_pipeline(topic: str, command: str) -> dict[str, Any]:
    info = TOPICS[topic]
    outdir = info["mixed_result"]
    outdir.mkdir(parents=True, exist_ok=True)
    catalog = info["mixed_catalog"]
    metrics = run_slice(catalog, outdir, HARMONIZED_DIR, min_studies=2)
    harmonized = HARMONIZED_DIR / f"{outdir.name}_harmonized.csv"
    score_summary = write_score_database(
        harmonized,
        outdir,
        catalog_path=catalog,
        db_path=outdir / "degora_scores.db",
        min_studies=2,
        command=command,
        extra_metadata={
            "benchmark_mode": "rna_seq_plus_microarray",
            "indisulam_excluded": True,
            "cross_platform_claim": "microarray added to a locked RNA-seq benchmark topic to assess marker-rank sharpening",
        },
    )
    baseline_summary = write_baseline_outputs(harmonized, outdir / "baselines", corpus=outdir.name, min_studies=2)
    comparator = build_summary(
        outdir / "baselines",
        info["gold"],
        degora_score_csv=outdir / "degora_gene_scores.csv",
        degora_score_label="quality_weighted_default",
    )
    summary_stem = str(info["mixed_summary_name"])
    summary_csv = outdir / f"{summary_stem}.csv"
    summary_md = outdir / f"{summary_stem}.md"
    comparator.to_csv(summary_csv, index=False)
    write_markdown(comparator, summary_md, title=f"{topic} RNA-seq+microarray comparator summary", gold_path=info["gold"])
    for artifact in (summary_csv, summary_md):
        write_source_sidecar(
            artifact,
            command,
            inputs=[outdir / "baselines", info["gold"], outdir / "degora_gene_scores.csv"],
            metadata={"generator": "cross-platform-gold-comparator-summary", "topic": topic},
        )
    return {
        "topic": topic,
        "catalog": str(catalog),
        "harmonized": str(harmonized),
        "outdir": str(outdir),
        "slice_metrics": metrics,
        "score_summary": score_summary,
        "baseline_summary": baseline_summary,
        "comparator_summary": str(summary_csv),
    }


def _read_degora_quality_weighted(summary_path: Path) -> pd.Series:
    frame = pd.read_csv(summary_path)
    subset = frame.loc[frame["method_id"].eq("degora_quality_weighted_score")].copy()
    if subset.empty:
        subset = frame.loc[frame["method_id"].eq("degora_deg_score")].copy()
    for column in ["recall_at_10", "recall_at_20", "recall_at_50", "recall_at_100", "direction_recall_at_100"]:
        subset[column] = pd.to_numeric(subset[column], errors="coerce")
    if subset.empty:
        raise ValueError(f"{summary_path} does not contain a DEGORA score row")
    return subset.sort_values(["method_id", "setting_id"]).iloc[0]


def _benchmark_summary(command: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for topic, info in TOPICS.items():
        for mode, result_dir, summary_path in [
            ("RNA-seq only", info["rna_only_result"], info["rna_only_summary"]),
            ("RNA-seq + microarray", info["mixed_result"], info["mixed_result"] / f"{info['mixed_summary_name']}.csv"),
        ]:
            best = _read_degora_quality_weighted(summary_path)
            metadata = json.loads((result_dir / "degora_score_metadata.json").read_text())
            rows.append(
                {
                    "topic": topic,
                    "mode": mode,
                    "method_id": best["method_id"],
                    "setting_id": best["setting_id"],
                    "n_gene_scores": metadata.get("n_gene_scores"),
                    "n_source_units_total": metadata.get("n_source_units_total"),
                    "n_contrasts_total": metadata.get("n_contrasts_total"),
                    "recall_at_10": best["recall_at_10"],
                    "recall_at_20": best["recall_at_20"],
                    "recall_at_50": best["recall_at_50"],
                    "recall_at_100": best["recall_at_100"],
                    "direction_recall_at_100": best["direction_recall_at_100"],
                    "top10": best["top10"],
                }
            )
    frame = pd.DataFrame(rows)
    output = RESULT_ROOT / "cross_platform_benchmark_summary.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    write_source_sidecar(output, command, inputs=[TOPICS["IFN"]["rna_only_summary"], TOPICS["ER stress"]["rna_only_summary"]], metadata={"generator": "cross-platform-benchmark-summary"})
    return frame


def _rank_column(frame: pd.DataFrame) -> str:
    if "quality_weighted_degora_rank" in frame.columns:
        return "quality_weighted_degora_rank"
    return "degora_rank"


def _direction_column(frame: pd.DataFrame) -> str:
    if "quality_weighted_consensus_direction" in frame.columns:
        return "quality_weighted_consensus_direction"
    return "consensus_direction"


def _marker_rank_shifts(command: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for topic, info in TOPICS.items():
        gold = pd.read_csv(info["gold"])
        genes = gold["gene_symbol"].dropna().astype(str).str.upper().tolist()
        rna = pd.read_csv(info["rna_only_result"] / "degora_gene_scores.csv")
        mixed = pd.read_csv(info["mixed_result"] / "degora_gene_scores.csv")
        rna_rank_col = _rank_column(rna)
        mixed_rank_col = _rank_column(mixed)
        rna_dir_col = _direction_column(rna)
        mixed_dir_col = _direction_column(mixed)
        rna_map = rna.set_index(rna["gene_symbol"].astype(str).str.upper())
        mixed_map = mixed.set_index(mixed["gene_symbol"].astype(str).str.upper())
        for gene in genes:
            if gene not in rna_map.index or gene not in mixed_map.index:
                continue
            rna_row = rna_map.loc[gene]
            mixed_row = mixed_map.loc[gene]
            if isinstance(rna_row, pd.DataFrame):
                rna_row = rna_row.iloc[0]
            if isinstance(mixed_row, pd.DataFrame):
                mixed_row = mixed_row.iloc[0]
            rna_rank = float(rna_row[rna_rank_col])
            mixed_rank = float(mixed_row[mixed_rank_col])
            rows.append(
                {
                    "topic": topic,
                    "gene_symbol": gene,
                    "expected_direction": gold.loc[gold["gene_symbol"].astype(str).str.upper().eq(gene), "expected_direction"].iloc[0]
                    if "expected_direction" in gold.columns
                    else "",
                    "rna_only_rank": rna_rank,
                    "mixed_rank": mixed_rank,
                    "rank_delta": rna_rank - mixed_rank,
                    "rna_only_top_percent": float(rna_row.get("quality_weighted_top_percent", rna_row.get("top_percent"))),
                    "mixed_top_percent": float(mixed_row.get("quality_weighted_top_percent", mixed_row.get("top_percent"))),
                    "top_percent_delta": float(rna_row.get("quality_weighted_top_percent", rna_row.get("top_percent")))
                    - float(mixed_row.get("quality_weighted_top_percent", mixed_row.get("top_percent"))),
                    "rna_only_direction": rna_row.get(rna_dir_col, ""),
                    "mixed_direction": mixed_row.get(mixed_dir_col, ""),
                    "rna_only_source_units": int(rna_row.get("n_source_units", 0)),
                    "mixed_source_units": int(mixed_row.get("n_source_units", 0)),
                    "rna_only_reliability": float(rna_row.get("evidence_reliability_score", 0)),
                    "mixed_reliability": float(mixed_row.get("evidence_reliability_score", 0)),
                }
            )
    frame = pd.DataFrame(rows)
    output = RESULT_ROOT / "cross_platform_marker_rank_shifts.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    write_source_sidecar(output, command, metadata={"generator": "cross-platform-marker-rank-shifts"})
    return frame


def _dataset_sources(command: str, preparation: dict[str, Any]) -> pd.DataFrame:
    rows = [
        {
            "topic": "IFN",
            "accession": "GSE71634",
            "assay_type": "microarray",
            "platform": "GPL10558",
            "samples": preparation["sample_groups"]["ifn_ctrl"] + preparation["sample_groups"]["ifn_treat"],
            "contrast": "IFN-beta 18h vs unstimulated PBMC",
            "source_url": GEO_URLS["GSE71634"],
            "use": "added to RNA-seq-only IFN benchmark",
        },
        {
            "topic": "ER stress",
            "accession": "GSE19519",
            "assay_type": "microarray",
            "platform": "GPL570",
            "samples": preparation["sample_groups"]["er_tunicamycin_ctrl"] + preparation["sample_groups"]["er_tunicamycin_treat"],
            "contrast": "tunicamycin 8h vs DMSO LCL",
            "source_url": GEO_URLS["GSE19519"],
            "use": "added to RNA-seq-only ER primary benchmark",
        },
        {
            "topic": "ER stress",
            "accession": "GSE19519",
            "assay_type": "microarray",
            "platform": "GPL570",
            "samples": preparation["sample_groups"]["er_thapsigargin_ctrl"] + preparation["sample_groups"]["er_thapsigargin_treat"],
            "contrast": "thapsigargin 4h vs DMSO LCL",
            "source_url": GEO_URLS["GSE19519"],
            "use": "same source unit as GSE19519 tunicamycin row",
        },
        {
            "topic": "Hypoxia",
            "accession": "GSE3045",
            "assay_type": "microarray",
            "platform": "GPL570",
            "samples": preparation["sample_groups"]["hypoxia_astrocyte_ctrl"] + preparation["sample_groups"]["hypoxia_astrocyte_treat"],
            "contrast": "primary human astrocyte 1% O2 24h vs normoxia",
            "source_url": GEO_URLS["GSE3045"],
            "use": "added to human-only hypoxia benchmark",
        },
        {
            "topic": "Hypoxia",
            "accession": "GSE22282",
            "assay_type": "microarray",
            "platform": "GPL570",
            "samples": preparation["sample_groups"]["hypoxia_mdc_ctrl"] + preparation["sample_groups"]["hypoxia_mdc_treat"],
            "contrast": "mature dendritic-cell hypoxia vs normoxia",
            "source_url": GEO_URLS["GSE22282"],
            "use": "added to human-only hypoxia benchmark",
        },
    ]
    frame = pd.DataFrame(rows)
    output = RESULT_ROOT / "cross_platform_dataset_sources.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    write_source_sidecar(output, command, metadata={"generator": "cross-platform-dataset-sources"})
    return frame


def write_cross_platform_benchmarks(output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    command = shell_command(["PYTHONPATH=outputs/code", "python", "outputs/code/scripts/write_cross_platform_microarray_benchmarks.py", "--output-dir", output_dir])
    preparation = _derive_microarray_inputs(command)
    for topic, rows in _catalog_rows().items():
        _append_catalog(TOPICS[topic]["base_catalog"], rows, TOPICS[topic]["mixed_catalog"], command)
    pipeline = [_run_pipeline(topic, command) for topic in TOPICS]
    benchmark = _benchmark_summary(command)
    rank_shifts = _marker_rank_shifts(command)
    datasets = _dataset_sources(command, preparation)
    summary = {
        "indisulam_excluded": True,
        "output_dir": str(output_dir),
        "preparation": preparation,
        "pipeline": pipeline,
        "tables": {
            "benchmark_summary": str(output_dir / "cross_platform_benchmark_summary.csv"),
            "marker_rank_shifts": str(output_dir / "cross_platform_marker_rank_shifts.csv"),
            "dataset_sources": str(output_dir / "cross_platform_dataset_sources.csv"),
        },
        "n_benchmark_rows": int(len(benchmark)),
        "n_rank_shift_rows": int(len(rank_shifts)),
        "n_dataset_rows": int(len(datasets)),
    }
    summary_path = output_dir / "cross_platform_microarray_benchmark_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(summary_path, command, metadata={"generator": "cross-platform-microarray-benchmark-summary"})
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = write_cross_platform_benchmarks(args.output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
