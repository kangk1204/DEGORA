#!/usr/bin/env python
"""Prepare GSE93829 microarray sensitivity inputs for the Indisulam/RBM39 axis."""

from __future__ import annotations

import argparse
import gzip
import json
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd

from derive_microarray_deg import derive_microarray_deg
from degora.provenance import write_source_sidecar


GSE93829_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93829/matrix/GSE93829_series_matrix.txt.gz"
GPL17077_TEXT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17077&targ=self&form=text&view=full"
GSE93829_GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93829"

CONTROL_DMSO = ["GSM2463166", "GSM2463167", "GSM2463168"]
TREAT_E7820 = ["GSM2463169", "GSM2463170", "GSM2463171"]
CONTROL_SIRNA = ["GSM2463172", "GSM2463173", "GSM2463174"]
TREAT_SIRBM39 = ["GSM2463175", "GSM2463176", "GSM2463177"]

GUARDRAILS = {
    "claim_scope": "cross_platform_rbm39_sulfonamide_mechanism_sensitivity",
    "as_published_author_deg_validation": False,
    "primary_indisulam_result_replaced": False,
    "includes_microarray": True,
    "microarray_direct_indisulam": False,
    "microarray_compound": "E7820",
    "microarray_interpretation": (
        "GSE93829 is an E7820/RBM39-axis microarray sensitivity layer, not a direct Indisulam treatment dataset."
    ),
}


def _read_lines(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read().splitlines()
    return path.read_text(errors="replace").splitlines()


def _between_markers(lines: list[str], start_marker: str, end_marker: str) -> pd.DataFrame:
    start = lines.index(start_marker) + 1
    end = lines.index(end_marker)
    return pd.read_csv(StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)


def _series_matrix(series_matrix_path: Path) -> pd.DataFrame:
    matrix = _between_markers(_read_lines(series_matrix_path), "!series_matrix_table_begin", "!series_matrix_table_end")
    matrix = matrix.rename(columns={"ID_REF": "probe_id"})
    return matrix


def _platform_annotation(platform_path: Path) -> pd.DataFrame:
    annotation = _between_markers(_read_lines(platform_path), "!platform_table_begin", "!platform_table_end")
    return annotation.loc[:, ["ID", "GENE_SYMBOL"]].rename(columns={"ID": "probe_id", "GENE_SYMBOL": "gene_symbol"})


def _write_probe_matrix(series_matrix_path: Path, platform_path: Path, output_path: Path, *, command: str) -> dict[str, Any]:
    expression = _series_matrix(series_matrix_path)
    annotation = _platform_annotation(platform_path)
    merged = expression.merge(annotation, on="probe_id", how="left")
    merged["gene_symbol"] = merged["gene_symbol"].astype("string").str.strip().str.upper()
    merged = merged.loc[merged["gene_symbol"].notna() & merged["gene_symbol"].ne("")].copy()
    ordered = ["gene_symbol", "probe_id", *CONTROL_DMSO, *TREAT_E7820, *CONTROL_SIRNA, *TREAT_SIRBM39]
    merged = merged.loc[:, ordered]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False)
    write_source_sidecar(
        output_path,
        command,
        inputs=[series_matrix_path, platform_path],
        metadata={"generator": "gse93829-probe-gene-matrix", **GUARDRAILS},
    )
    return {
        "probe_matrix": str(output_path),
        "n_probe_rows_with_symbol": int(len(merged)),
        "n_gene_symbols": int(merged["gene_symbol"].nunique()),
    }


def _microarray_catalog_rows(raw_dir: Path) -> list[dict[str, Any]]:
    common = {
        "paper_id": "GSE93829_Uehara2017",
        "source_url": GSE93829_GEO_URL,
        "pipeline": "welch_microarray_normalized_matrix",
        "assay_type": "microarray",
        "source_input_type": "normalized_expression_matrix",
        "platform": "GPL17077",
        "normalization": "GeneSpring_12.5_processed_series_matrix_log_scale",
        "probe_id_column": "probe_id",
        "probe_collapse": "min_pvalue_max_abs_lfc",
        "species": "Homo sapiens",
        "cell_system": "HCT116",
        "duration_h": "",
        "n_ctrl": 3,
        "n_treat": 3,
        "gene_column": "gene_symbol",
        "lfc_column": "log2FoldChange",
        "p_column": "pvalue",
        "padj_column": "padj",
        "sep": "",
        "sheet_name": "",
        "gene_type_column": "",
        "gene_type_keep": "",
        "table_scope": "full_results",
        "rank_universe_size": "",
        "include_in_analysis": True,
    }
    return [
        {
            **common,
            "study_id": "SULFONAMIDE_GSE93829_HCT116_E7820_24h",
            "source_path": str(raw_dir / "GSE93829_HCT116_E7820_24h_vs_DMSO_microarray_welch.csv"),
            "hypoxia_modality": "E7820 1uM 24h vs DMSO",
            "duration_h": "24",
            "notes": "Microarray sensitivity layer for aryl-sulfonamide/RBM39 axis; not direct Indisulam.",
        },
        {
            **common,
            "study_id": "RBM39KD_GSE93829_HCT116_siRBM39_48h",
            "source_path": str(raw_dir / "GSE93829_HCT116_siRBM39_48h_vs_nontarget_microarray_welch.csv"),
            "hypoxia_modality": "siRBM39 48h vs non-target siRNA",
            "duration_h": "48",
            "notes": "Genetic RBM39 perturbation microarray sensitivity layer; mechanism anchor, not drug treatment.",
        },
    ]


def write_microarray_sensitivity(
    series_matrix: Path,
    platform: Path,
    base_catalog: Path,
    output_catalog: Path,
    raw_dir: Path,
    output_dir: Path,
    *,
    command: str,
) -> dict[str, Any]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    for path, url, key in [
        (series_matrix, GSE93829_MATRIX_URL, "gse93829_series_matrix"),
        (platform, GPL17077_TEXT_URL, "gpl17077_platform"),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")
        write_source_sidecar(path, f"curl -L -o {path} {url}", metadata={"generator": "manual-download-source", "source_key": key, "source_url": url})

    probe_matrix = raw_dir / "GSE93829_normalized_probe_gene_matrix.csv"
    probe_summary = _write_probe_matrix(series_matrix, platform, probe_matrix, command=command)

    e7820_path = raw_dir / "GSE93829_HCT116_E7820_24h_vs_DMSO_microarray_welch.csv"
    sirbm39_path = raw_dir / "GSE93829_HCT116_siRBM39_48h_vs_nontarget_microarray_welch.csv"
    e7820_summary = derive_microarray_deg(
        probe_matrix,
        e7820_path,
        gene_column="gene_symbol",
        probe_column="probe_id",
        control_samples=CONTROL_DMSO,
        treatment_samples=TREAT_E7820,
        command=command,
        summary_path=output_dir / "gse93829_e7820_microarray_summary.json",
        metadata={"study_id": "SULFONAMIDE_GSE93829_HCT116_E7820_24h", "paper_id": "GSE93829_Uehara2017", "platform": "GPL17077", "normalization": "GeneSpring_12.5_processed_series_matrix_log_scale", **GUARDRAILS},
    )
    sirbm39_summary = derive_microarray_deg(
        probe_matrix,
        sirbm39_path,
        gene_column="gene_symbol",
        probe_column="probe_id",
        control_samples=CONTROL_SIRNA,
        treatment_samples=TREAT_SIRBM39,
        command=command,
        summary_path=output_dir / "gse93829_sirbm39_microarray_summary.json",
        metadata={"study_id": "RBM39KD_GSE93829_HCT116_siRBM39_48h", "paper_id": "GSE93829_Uehara2017", "platform": "GPL17077", "normalization": "GeneSpring_12.5_processed_series_matrix_log_scale", **GUARDRAILS},
    )

    base = pd.read_csv(base_catalog)
    catalog = pd.concat([base, pd.DataFrame(_microarray_catalog_rows(raw_dir))], ignore_index=True, sort=False)
    output_catalog.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(output_catalog, index=False)
    write_source_sidecar(output_catalog, command, inputs=[base_catalog, e7820_path, sirbm39_path], metadata={"generator": "indisulam-cross-platform-catalog", **GUARDRAILS})

    availability = [
        {
            "source_id": "GSE93829",
            "assay_type": "microarray",
            "platform": "GPL17077",
            "status": "active_sensitivity_source",
            "direct_indisulam": False,
            "compound_or_perturbation": "E7820 and siRBM39",
            "source_url": GSE93829_GEO_URL,
            "use_decision": "included_as_cross_platform_mechanism_sensitivity",
            "reason": "E7820/RBM39-axis microarray; related aryl-sulfonamide mechanism, not direct Indisulam treatment.",
        }
    ]
    availability_tsv = output_dir / "indisulam_microarray_source_availability.tsv"
    availability_json = output_dir / "indisulam_microarray_source_availability.json"
    pd.DataFrame(availability).to_csv(availability_tsv, sep="\t", index=False)
    availability_json.write_text(json.dumps(availability, indent=2, sort_keys=True) + "\n")

    extra_metadata = {
        **GUARDRAILS,
        "topic": "Indisulam/E7070 with E7820/RBM39-axis microarray sensitivity",
        "primary_score_relation": "sensitivity_layer_not_replacement",
        "microarray_source_unit": "GSE93829_Uehara2017",
    }
    extra_metadata_path = output_dir / "indisulam_microarray_score_metadata_extra.json"
    extra_metadata_path.write_text(json.dumps(extra_metadata, indent=2, sort_keys=True) + "\n")

    summary = {
        **extra_metadata,
        "probe_summary": probe_summary,
        "e7820_summary": e7820_summary,
        "sirbm39_summary": sirbm39_summary,
        "catalog": str(output_catalog),
        "source_availability_tsv": str(availability_tsv),
        "score_metadata_extra": str(extra_metadata_path),
    }
    summary_path = output_dir / "indisulam_microarray_preparation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    for artifact in (availability_tsv, availability_json, extra_metadata_path, summary_path):
        write_source_sidecar(artifact, command, inputs=[series_matrix, platform, output_catalog], metadata={"generator": "indisulam-microarray-sensitivity", **GUARDRAILS})
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series-matrix", type=Path, default=Path("data/deg/raw/microarray/GSE93829_series_matrix.txt.gz"))
    parser.add_argument("--platform", type=Path, default=Path("data/deg/raw/microarray/GPL17077_platform.txt"))
    parser.add_argument("--base-catalog", type=Path, default=Path("data/studies/indisulam_derived_catalog.csv"))
    parser.add_argument("--output-catalog", type=Path, default=Path("data/studies/indisulam_cross_platform_catalog.csv"))
    parser.add_argument("--raw-dir", type=Path, default=Path("data/deg/raw/microarray"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/results/indisulam-microarray"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_indisulam_microarray_sensitivity.py"
        f" --series-matrix {args.series_matrix}"
        f" --platform {args.platform}"
        f" --base-catalog {args.base_catalog}"
        f" --output-catalog {args.output_catalog}"
        f" --raw-dir {args.raw_dir}"
        f" --output-dir {args.output_dir}"
    )
    summary = write_microarray_sensitivity(
        args.series_matrix,
        args.platform,
        args.base_catalog,
        args.output_catalog,
        args.raw_dir,
        args.output_dir,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
