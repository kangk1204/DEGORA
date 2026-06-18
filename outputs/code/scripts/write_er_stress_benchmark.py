#!/usr/bin/env python
"""Prepare a locked ER-stress/UPR benchmark corpus for DEGORA.

The gold panel is written before any score output is produced. Active sources
mix one as-published Cuffdiff table with two public-matrix derived sources; the
derived sources are explicitly labeled so they cannot be mistaken for a primary
as-published DEG-table validation corpus.
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
    "gse84450": {
        "path": "data/deg/raw/er_stress/GSE84450_WT_UTvsWT_DTT_gene_exp.diff.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84450/suppl/GSE84450_WT_UTvsWT_DTT_gene_exp.diff.gz",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84450",
    },
    "gse102505": {
        "path": "data/deg/raw/er_stress/GSE102505_glioma_nh125_exp_count_matrix.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102505/suppl/GSE102505_glioma_nh125_exp_count_matrix.csv.gz",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102505",
    },
    "gse103667": {
        "path": "data/deg/raw/er_stress/GSE103667_TE.norm.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103667/suppl/GSE103667_TE.norm.txt.gz",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103667",
    },
    "gse84989": {
        "path": "data/deg/raw/er_stress/GSE84989_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84989/suppl/GSE84989_RAW.tar",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84989",
    },
    "gse296996": {
        "path": "data/deg/raw/er_stress/GSE296996_count_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296996/suppl/GSE296996_count_matrix.txt.gz",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296996",
    },
    "gse245918": {
        "path": "data/deg/raw/er_stress/GSE245918_rsem_counts.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245918/suppl/GSE245918_rsem_counts.txt.gz",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245918",
    },
    "hgnc": {
        "path": "data/deg/raw/ifn/hgnc_complete_set.txt",
        "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        "geo": "",
    },
}

ER_STRESS_GOLD_GENES = [
    ("HSPA5", "core_chaperone", "up"),
    ("HSP90B1", "core_chaperone", "up"),
    ("DNAJB9", "ERAD_chaperone", "up"),
    ("HERPUD1", "ERAD", "up"),
    ("DDIT3", "PERK_ATF4_CHOP", "up"),
    ("ATF4", "PERK_ATF4_CHOP", "up"),
    ("ATF3", "integrated_stress_response", "up"),
    ("PPP1R15A", "integrated_stress_response", "up"),
    ("CHAC1", "integrated_stress_response", "up"),
    ("ASNS", "amino_acid_stress_response", "up"),
    ("PDIA4", "protein_folding", "up"),
    ("EDEM1", "ERAD", "up"),
    ("SEL1L", "ERAD", "up"),
    ("DERL3", "ERAD", "up"),
    ("MANF", "ER_stress_secreted_factor", "up"),
    ("TRIB3", "integrated_stress_response", "up"),
    ("DDIT4", "stress_response", "up"),
    ("XBP1", "IRE1_XBP1_axis", "up"),
]

GUARDRAIL_FLAGS = {
    "benchmark_topic": "ER stress / unfolded protein response",
    "gold_panel_locked_before_scoring": True,
    "primary_min_studies": 2,
    "claim_scope": "gold_panel_prioritization_in_a_prelocked_er_stress_upr_benchmark",
    "score_interpretation": "prioritization_index_not_calibrated_probability",
    "time_course_policy": (
        "Related time points or cell-line contrasts from one accession share a paper_id/source-unit and are "
        "mean-aggregated by the primary score; no max-signal source-unit selection is used."
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
    pipeline: str
    assay_type: str
    source_input_type: str
    platform: str = ""
    normalization: str = ""
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


def _collapse_deg_rows(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["gene_symbol"] = _clean_symbol(out["gene_symbol"])
    out = out.loc[out["gene_symbol"].notna()].copy()
    out["abs_lfc"] = pd.to_numeric(out["log2FoldChange"], errors="coerce").abs()
    out["n_source_rows_for_gene"] = out.groupby("gene_symbol")["gene_symbol"].transform("size")
    out = (
        out.sort_values(["gene_symbol", "pvalue", "abs_lfc", "padj"], ascending=[True, True, False, True])
        .drop_duplicates("gene_symbol", keep="first")
        .drop(columns=["abs_lfc"])
        .sort_values(["pvalue", "gene_symbol"])
        .reset_index(drop=True)
    )
    return out


def _log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    libraries = numeric.sum(axis=0).replace(0.0, np.nan)
    cpm = numeric.divide(libraries, axis=1) * 1_000_000.0
    return np.log2(cpm.fillna(0.0) + 1.0)


def _read_hgnc_map(path: Path) -> dict[str, str]:
    hgnc = pd.read_csv(path, sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id", "locus_group"])
    hgnc = hgnc.dropna(subset=["symbol", "ensembl_gene_id"])
    hgnc = hgnc.loc[hgnc["locus_group"].astype("string").str.lower().eq("protein-coding gene")].copy()
    hgnc["ensembl_gene_id"] = hgnc["ensembl_gene_id"].str.replace(r"\.\d+$", "", regex=True).str.upper()
    hgnc["symbol"] = hgnc["symbol"].str.upper()
    return dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"], strict=False))


def _read_human_ensembl_counts(path: Path, hgnc_map: dict[str, str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    raw = pd.read_csv(path, index_col=0)
    ensembl_ids = raw.index.astype("string").str.replace(r"\.\d+$", "", regex=True).str.upper()
    symbols = pd.Series(ensembl_ids.map(hgnc_map), index=raw.index)
    counts = raw.loc[symbols.notna()].copy()
    counts.index = symbols.loc[symbols.notna()].astype(str)
    collapsed = _collapse_by_symbol(counts)
    return collapsed, {
        "raw_rows": int(len(raw)),
        "mapped_rows": int(symbols.notna().sum()),
        "n_unmapped_genes": int(symbols.isna().sum()),
        "n_gene_symbols_after_collapse": int(len(collapsed)),
    }


def _read_ensembl_count_table(
    path: Path,
    hgnc_map: dict[str, str],
    *,
    sep: str = ",",
    ensembl_column: str | None = None,
    drop_columns: tuple[str, ...] = (),
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read a tab/CSV ENSG count table by the same ENSG->HGNC recipe as GSE102505.

    ``ensembl_column`` selects the gene-id column when the table is not indexed by
    the first column; ``drop_columns`` removes annotation columns (such as a
    secondary ``gene_name`` field) before the count matrix is assembled.
    """

    raw = pd.read_csv(path, sep=sep)
    if ensembl_column is None:
        ensembl_column = str(raw.columns[0])
    ensembl_series = raw[ensembl_column].astype("string").str.replace(r"\.\d+$", "", regex=True).str.upper()
    count_columns = [
        column
        for column in raw.columns
        if column != ensembl_column and column not in set(drop_columns)
    ]
    counts = raw[count_columns].copy()
    counts.index = ensembl_series
    symbols = pd.Series(ensembl_series.map(hgnc_map).to_numpy(), index=counts.index)
    counts = counts.loc[symbols.notna()].copy()
    counts.index = symbols.loc[symbols.notna()].astype(str)
    collapsed = _collapse_by_symbol(counts)
    return collapsed, {
        "raw_rows": int(len(raw)),
        "mapped_rows": int(symbols.notna().sum()),
        "n_unmapped_genes": int(symbols.isna().sum()),
        "n_gene_symbols_after_collapse": int(len(collapsed)),
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


def _clean_gse84450_cuffdiff(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    raw = pd.read_csv(path, sep="\t")
    ok = raw.loc[raw["status"].astype("string").str.upper().eq("OK")].copy()
    lfc = pd.to_numeric(ok["log2(fold_change)"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    clean = pd.DataFrame(
        {
            "gene_symbol": ok["gene"],
            "log2FoldChange": lfc.clip(lower=-10.0, upper=10.0),
            "pvalue": pd.to_numeric(ok["p_value"], errors="coerce").clip(lower=0.0, upper=1.0),
            "padj": pd.to_numeric(ok["q_value"], errors="coerce").clip(lower=0.0, upper=1.0),
            "cuffdiff_status": ok["status"],
            "cuffdiff_value_1": ok["value_1"],
            "cuffdiff_value_2": ok["value_2"],
        }
    )
    clean = clean.dropna(subset=["log2FoldChange", "pvalue"]).copy()
    clean = _collapse_deg_rows(clean)
    summary = {
        "raw_rows": int(len(raw)),
        "status_ok_rows": int(len(ok)),
        "n_gene_symbols_after_collapse": int(len(clean)),
        "n_dropped_non_ok_or_invalid": int(len(raw) - len(clean)),
    }
    return clean, summary


def _derive_normalized_expression_contrast(path: Path, spec: ContrastSpec) -> tuple[pd.DataFrame, dict[str, Any]]:
    matrix = pd.read_csv(path, sep="\t")
    required = ["name", *spec.control_columns, *spec.treat_columns]
    missing = sorted(set(required).difference(matrix.columns))
    if missing:
        raise ValueError(f"{spec.study_id} missing columns: {missing}")
    expression = matrix[list(spec.control_columns + spec.treat_columns)].apply(pd.to_numeric, errors="coerce")
    genes = _clean_symbol(matrix["name"])
    valid = genes.notna() & expression.notna().all(axis=1)
    expression = expression.loc[valid].copy()
    genes = genes.loc[valid]
    control = expression[list(spec.control_columns)]
    treat = expression[list(spec.treat_columns)]
    log2fc = treat.mean(axis=1) - control.mean(axis=1)
    test = stats.ttest_ind(treat.to_numpy(dtype=float), control.to_numpy(dtype=float), axis=1, equal_var=False, nan_policy="omit")
    pvalue = pd.Series(test.pvalue, index=expression.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    out = pd.DataFrame({"gene_symbol": genes.astype(str), "log2FoldChange": log2fc, "pvalue": pvalue})
    out["padj"] = _bh_adjust(out["pvalue"])
    out = _collapse_deg_rows(out)
    out["source_input_type"] = spec.source_input_type
    summary = {
        "raw_rows": int(len(matrix)),
        "valid_rows": int(valid.sum()),
        "n_gene_symbols_after_collapse": int(len(out)),
        "n_collapsed_duplicate_symbol_rows": int(valid.sum() - len(out)),
    }
    return out, summary


def _gse102505_specs(raw_dir: Path) -> list[ContrastSpec]:
    common = {
        "paper_id": "GSE102505_glioma_TUN",
        "source_url": RAW_INPUTS["gse102505"]["geo"],
        "raw_source_key": "gse102505",
        "species": "Homo sapiens",
        "condition": "tunicamycin vs DMSO",
        "duration_h": "not_reported_in_matrix",
        "n_ctrl": 3,
        "n_treat": 3,
        "pipeline": "logCPM_Welch_derived_from_public_counts",
        "assay_type": "RNA-seq",
        "source_input_type": "derived_count_table",
        "normalization": "logCPM_from_raw_counts",
        "notes": "Count-derived ER-stress benchmark contrast; same GSE102505 paper_id collapses cell-line contrasts to one source unit.",
    }
    specs = []
    for prefix, label in (("t4213", "T4213 glioma"), ("u251", "U251 glioma"), ("nha", "normal human astrocytes")):
        specs.append(
            ContrastSpec(
                study_id=f"ER_GSE102505_{prefix.upper()}_TUN_vs_DMSO",
                source_path=raw_dir / f"GSE102505_{prefix.upper()}_tunicamycin_vs_DMSO_derived_logcpm_welch.csv",
                cell_system=label,
                control_columns=(f"{prefix}_dmso_1", f"{prefix}_dmso_2", f"{prefix}_dmso_3"),
                treat_columns=(f"{prefix}_tun_1", f"{prefix}_tun_2", f"{prefix}_tun_3"),
                output_name=f"GSE102505_{prefix.upper()}_tunicamycin_vs_DMSO_derived_logcpm_welch.csv",
                **common,
            )
        )
    return specs


def _gse296996_specs(raw_dir: Path) -> list[ContrastSpec]:
    return [
        ContrastSpec(
            study_id="ER_GSE296996_HEPG2_TUN_vs_DMSO",
            paper_id="GSE296996_HepG2_TUN",
            source_url=RAW_INPUTS["gse296996"]["geo"],
            source_path=raw_dir / "GSE296996_HepG2_tunicamycin_vs_DMSO_derived_logcpm_welch.csv",
            raw_source_key="gse296996",
            species="Homo sapiens",
            cell_system="HepG2 hepatoblastoma",
            condition="tunicamycin vs DMSO",
            duration_h="not_reported_in_matrix",
            n_ctrl=5,
            n_treat=5,
            control_columns=("HepG2-Ctrl_rep1", "HepG2-Ctrl_rep2", "HepG2-Ctrl_rep3", "HepG2-Ctrl_rep4", "HepG2-Ctrl_rep5"),
            treat_columns=(
                "HepG2-Tunicamycin_rep1",
                "HepG2-Tunicamycin_rep2",
                "HepG2-Tunicamycin_rep3",
                "HepG2-Tunicamycin_rep4",
                "HepG2-Tunicamycin_rep5",
            ),
            output_name="GSE296996_HepG2_tunicamycin_vs_DMSO_derived_logcpm_welch.csv",
            pipeline="logCPM_Welch_derived_from_public_counts",
            assay_type="RNA-seq",
            source_input_type="derived_count_table",
            normalization="logCPM_from_raw_counts",
            notes="Count-derived human ER-stress benchmark contrast (HepG2 tunicamycin vs DMSO, GSM8981325-29 vs GSM8981330-34) added to reach >=2 independent human RNA-seq source units.",
        )
    ]


def _gse245918_specs(raw_dir: Path) -> list[ContrastSpec]:
    common = {
        "paper_id": "GSE245918_TUN",
        "source_url": RAW_INPUTS["gse245918"]["geo"],
        "raw_source_key": "gse245918",
        "species": "Homo sapiens",
        "condition": "tunicamycin vs DMSO",
        "duration_h": "24",
        "n_ctrl": 3,
        "n_treat": 3,
        "pipeline": "logCPM_Welch_derived_from_public_counts",
        "assay_type": "RNA-seq",
        "source_input_type": "derived_count_table",
        "normalization": "logCPM_from_raw_counts",
        "notes": "Count-derived human ER-stress benchmark contrast; same GSE245918 paper_id collapses OCI-AML3 and HEK293T cell-line contrasts to one source unit.",
    }
    specs = []
    cell_specs = (
        (
            "OCI_AML3",
            "OCI-AML3 acute myeloid leukemia",
            ("A10", "A11", "A12"),  # OCI-AML3 24H DMSO (GSM7851003-05)
            ("A13", "A14", "A15"),  # OCI-AML3 24H tunicamycin (GSM7851006-08)
        ),
        (
            "HEK293T",
            "HEK293T embryonic kidney",
            ("A16", "A17", "A18"),  # HEK293T DMSO 24h (GSM7851009-11)
            ("A19", "A20", "A21"),  # HEK293T tunicamycin 24h (GSM7851012-14)
        ),
    )
    for prefix, label, control_columns, treat_columns in cell_specs:
        specs.append(
            ContrastSpec(
                study_id=f"ER_GSE245918_{prefix}_TUN_vs_DMSO",
                source_path=raw_dir / f"GSE245918_{prefix}_tunicamycin_vs_DMSO_derived_logcpm_welch.csv",
                cell_system=label,
                control_columns=control_columns,
                treat_columns=treat_columns,
                output_name=f"GSE245918_{prefix}_tunicamycin_vs_DMSO_derived_logcpm_welch.csv",
                **common,
            )
        )
    return specs


def _all_specs(raw_dir: Path) -> list[ContrastSpec]:
    return [
        ContrastSpec(
            study_id="ER_GSE84450_WT_DTT_vs_UT",
            paper_id="GSE84450_DTT",
            source_url=RAW_INPUTS["gse84450"]["geo"],
            source_path=raw_dir / "GSE84450_WT_DTT_vs_UT_cuffdiff_clean.csv",
            raw_source_key="gse84450",
            species="Mus musculus",
            cell_system="mouse embryonic fibroblast WT",
            condition="DTT vs untreated",
            duration_h="not_reported_in_table",
            n_ctrl=2,
            n_treat=2,
            control_columns=(),
            treat_columns=(),
            output_name="GSE84450_WT_DTT_vs_UT_cuffdiff_clean.csv",
            pipeline="Cuffdiff2_author_table",
            assay_type="RNA-seq",
            source_input_type="author_deg_table",
            normalization="Cuffdiff_FPKM",
            notes="As-published Cuffdiff gene_exp.diff full-result table; status != OK rows excluded before scoring.",
        ),
        *_gse102505_specs(raw_dir),
        *_gse296996_specs(raw_dir),
        *_gse245918_specs(raw_dir),
        ContrastSpec(
            study_id="ER_GSE103667_THAP_vs_DMSO",
            paper_id="GSE103667_THAP",
            source_url=RAW_INPUTS["gse103667"]["geo"],
            source_path=raw_dir / "GSE103667_THAP_vs_DMSO_normexpr_welch.csv",
            raw_source_key="gse103667",
            species="Mus musculus",
            cell_system="TE matrix source",
            condition="thapsigargin vs DMSO",
            duration_h="not_reported_in_matrix",
            n_ctrl=2,
            n_treat=2,
            control_columns=("TE.DMSO1", "TE.DMSO2"),
            treat_columns=("TE.THAP1", "TE.THAP2"),
            output_name="GSE103667_THAP_vs_DMSO_normexpr_welch.csv",
            pipeline="welch_normalized_expression_matrix",
            assay_type="RNA-seq",
            source_input_type="normalized_expression_matrix",
            normalization="TE.norm",
            notes=(
                "Secondary ER-stress source from normalized expression matrix with two replicates per group; "
                "included as lower-confidence support, not a full author DEG table."
            ),
        ),
    ]


def _catalog_rows(specs: list[ContrastSpec], row_counts: dict[str, int]) -> list[dict[str, Any]]:
    rows = []
    for spec in specs:
        is_human = spec.species == "Homo sapiens"
        notes = spec.notes
        if not is_human:
            notes = f"{notes} Excluded from the manuscript benchmark after the human-only species gate."
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
                "include_in_analysis": is_human,
                "notes": notes,
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
                "Pre-output ER stress/UPR sentinel panel centered on Hallmark UPR and canonical "
                "PERK/ATF4, IRE1/XBP1, chaperone, and ERAD transcriptional response genes."
            ),
            "notes": "Locked before ER stress DEGORA or baseline scoring; used only for recall/prioritization diagnostics.",
        }
        for symbol, evidence_class, direction in ER_STRESS_GOLD_GENES
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    write_source_sidecar(path, command, metadata={"generator": "er-stress-upr-gold-panel", **GUARDRAIL_FLAGS})


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
        file_record("gse84450", "species_sensitivity_source", "excluded_from_manuscript_primary", "as-published Cuffdiff full-result DTT ER-stress table; mouse source retained only for non-primary sensitivity after the human-only species gate"),
        file_record("gse102505", "active_primary_source", "included", "public count matrix with clear tunicamycin vs DMSO columns in three cell systems"),
        file_record("gse296996", "active_primary_source", "included", "public HepG2 count matrix with five tunicamycin vs five DMSO replicate columns"),
        file_record("gse245918", "active_primary_source", "included", "public RSEM count matrix with OCI-AML3 and HEK293T tunicamycin vs DMSO 24h arms mapped from the GEO series matrix"),
        file_record("gse103667", "species_sensitivity_source", "excluded_from_manuscript_primary", "normalized thapsigargin vs DMSO matrix with two replicates per arm; mouse source retained only for non-primary sensitivity after the human-only species gate"),
        file_record("gse84989", "deferred_candidate", "excluded_pending_metadata_lock", "sample-label metadata was not locked for this benchmark run"),
    ]


def _write_markdown_summary(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# ER Stress / UPR Benchmark Preparation",
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
            "- hStouffer/AWmeta failures remain blocker evidence, not wins.",
            "- Same-accession cell-line/time-course rows share a source unit.",
            "- Derived count or normalized-matrix rows are labeled as derived inputs.",
            "- The score is a prioritization index, not a calibrated probability.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def _write_primary_only_catalog(catalog: pd.DataFrame, path: Path, command: str, inputs: list[Path]) -> pd.DataFrame:
    primary = catalog.copy()
    secondary_mask = primary["source_input_type"].eq("normalized_expression_matrix")
    primary.loc[secondary_mask, "include_in_analysis"] = False
    primary.loc[secondary_mask, "notes"] = (
        primary.loc[secondary_mask, "notes"].astype(str)
        + " Excluded from primary-only ER benchmark because it is a lower-confidence normalized-matrix source; retained in full sensitivity."
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    primary.to_csv(path, index=False)
    write_source_sidecar(
        path,
        command,
        inputs=inputs,
        metadata={"generator": "er-stress-upr-primary-only-catalog", **GUARDRAIL_FLAGS},
    )
    return primary


def write_er_stress_benchmark(
    raw_dir: Path,
    catalog_path: Path,
    gold_path: Path,
    output_dir: Path,
    *,
    primary_catalog_path: Path,
    command: str,
) -> dict[str, Any]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_paths = {key: Path(info["path"]) for key, info in RAW_INPUTS.items()}
    for key, info in RAW_INPUTS.items():
        path = Path(info["path"])
        if not path.exists():
            raise FileNotFoundError(f"Missing required raw input: {path}")
        write_source_sidecar(
            path,
            shell_command(["curl", "-L", "-o", path, info["url"]]),
            inputs=[],
            metadata={"generator": "manual-download-source", "source_url": info["url"]},
        )

    _write_gold_panel(gold_path, command)
    hgnc_map = _read_hgnc_map(raw_paths["hgnc"])
    gse102_counts, gse102_mapping = _read_human_ensembl_counts(raw_paths["gse102505"], hgnc_map)
    gse296_counts, gse296_mapping = _read_ensembl_count_table(
        raw_paths["gse296996"], hgnc_map, sep="\t", ensembl_column="Geneid", drop_columns=("gene_name",)
    )
    gse245_counts, gse245_mapping = _read_ensembl_count_table(
        raw_paths["gse245918"], hgnc_map, sep="\t", ensembl_column="Ensembl.id"
    )
    specs = _all_specs(raw_dir)

    contrast_summaries: list[dict[str, Any]] = []
    row_counts: dict[str, int] = {}
    source_summaries: dict[str, dict[str, Any]] = {
        "gse102505": gse102_mapping,
        "gse296996": gse296_mapping,
        "gse245918": gse245_mapping,
    }
    for spec in specs:
        if spec.raw_source_key == "gse84450":
            deg, source_summary = _clean_gse84450_cuffdiff(raw_paths["gse84450"])
            source_summaries["gse84450"] = source_summary
            inputs = [raw_paths["gse84450"]]
            generator = "er-stress-gse84450-cuffdiff-clean"
        elif spec.raw_source_key == "gse102505":
            deg = _derive_logcpm_welch_contrast(gse102_counts, spec)
            inputs = [raw_paths["gse102505"], raw_paths["hgnc"]]
            generator = "er-stress-gse102505-logcpm-welch"
        elif spec.raw_source_key == "gse296996":
            deg = _derive_logcpm_welch_contrast(gse296_counts, spec)
            inputs = [raw_paths["gse296996"], raw_paths["hgnc"]]
            generator = "er-stress-gse296996-logcpm-welch"
        elif spec.raw_source_key == "gse245918":
            deg = _derive_logcpm_welch_contrast(gse245_counts, spec)
            inputs = [raw_paths["gse245918"], raw_paths["hgnc"]]
            generator = "er-stress-gse245918-logcpm-welch"
        elif spec.raw_source_key == "gse103667":
            deg, source_summary = _derive_normalized_expression_contrast(raw_paths["gse103667"], spec)
            source_summaries["gse103667"] = source_summary
            inputs = [raw_paths["gse103667"]]
            generator = "er-stress-gse103667-normalized-expression-welch"
        else:  # pragma: no cover - guarded by _all_specs.
            raise AssertionError(f"Unhandled raw source: {spec.raw_source_key}")

        deg_path = raw_dir / spec.output_name
        deg.to_csv(deg_path, index=False)
        row_counts[spec.study_id] = int(len(deg))
        write_source_sidecar(
            deg_path,
            command,
            inputs=inputs,
            metadata={"generator": generator, "study_id": spec.study_id, **GUARDRAIL_FLAGS},
        )
        top100_gold_hits = sorted(set(deg.head(100)["gene_symbol"]).intersection({gene for gene, _, _ in ER_STRESS_GOLD_GENES}))
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
    catalog_inputs = [
        raw_paths["gse84450"],
        raw_paths["gse102505"],
        raw_paths["gse296996"],
        raw_paths["gse245918"],
        raw_paths["gse103667"],
        raw_paths["hgnc"],
        gold_path,
    ]
    write_source_sidecar(
        catalog_path,
        command,
        inputs=catalog_inputs,
        metadata={"generator": "er-stress-upr-catalog", **GUARDRAIL_FLAGS},
    )
    primary_catalog = _write_primary_only_catalog(catalog, primary_catalog_path, command, catalog_inputs)

    availability = _availability_records(raw_paths, source_summaries)
    availability_tsv = output_dir / "er_stress_source_availability.tsv"
    availability_json = output_dir / "er_stress_source_availability.json"
    pd.DataFrame(availability).to_csv(availability_tsv, sep="\t", index=False)
    availability_json.write_text(json.dumps(availability, indent=2, sort_keys=True) + "\n")

    extra_metadata = {
        **GUARDRAIL_FLAGS,
        "gold_panel": str(gold_path),
        "active_corpus": "GSE84450 DTT author table + GSE102505/GSE296996/GSE245918 human tunicamycin count-derived contrasts + GSE103667 thapsigargin normalized matrix",
        "deferred_sources": ["GSE84989 sample-label metadata not locked"],
    }
    extra_metadata_path = output_dir / "er_stress_score_metadata_extra.json"
    extra_metadata_path.write_text(json.dumps(extra_metadata, indent=2, sort_keys=True) + "\n")

    summary: dict[str, Any] = {
        "corpus": "ER stress / unfolded protein response benchmark",
        **GUARDRAIL_FLAGS,
        "n_contrasts": len(specs),
        "n_source_units": int(catalog["paper_id"].nunique()),
        "n_primary_active_contrasts": int(pd.Series(primary_catalog["include_in_analysis"]).astype(str).str.lower().isin({"true", "1", "yes"}).sum()),
        "n_primary_source_units": int(primary_catalog.loc[primary_catalog["include_in_analysis"].astype(str).str.lower().isin({"true", "1", "yes"}), "paper_id"].nunique()),
        "gold_genes": [gene for gene, _, _ in ER_STRESS_GOLD_GENES],
        "catalog": str(catalog_path),
        "primary_catalog": str(primary_catalog_path),
        "gold_panel": str(gold_path),
        "contrast_summaries": contrast_summaries,
        "source_summaries": source_summaries,
        "source_availability": availability,
        "source_availability_tsv": str(availability_tsv),
        "score_metadata_extra": str(extra_metadata_path),
    }
    summary_json = output_dir / "er_stress_benchmark_preparation_summary.json"
    summary_md = output_dir / "er_stress_benchmark_preparation_summary.md"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    _write_markdown_summary(summary, summary_md)

    for artifact in (availability_tsv, availability_json, extra_metadata_path, summary_json, summary_md):
        write_source_sidecar(
            artifact,
            command,
            inputs=[catalog_path, gold_path, raw_paths["gse84450"], raw_paths["gse102505"], raw_paths["gse296996"], raw_paths["gse245918"], raw_paths["gse103667"], raw_paths["hgnc"]],
            metadata={"generator": "er-stress-upr-preparation-summary", **GUARDRAIL_FLAGS},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/deg/raw/er_stress"))
    parser.add_argument("--catalog", type=Path, default=Path("data/studies/er_stress_upr_catalog.csv"))
    parser.add_argument("--primary-catalog", type=Path, default=Path("data/studies/er_stress_upr_primary_catalog.csv"))
    parser.add_argument("--gold", type=Path, default=Path("data/studies/gold/er_stress_upr_gold_panel.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/results/er-stress-benchmark"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = shell_command(
        [
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_er_stress_benchmark.py",
            "--raw-dir",
            args.raw_dir,
            "--catalog",
            args.catalog,
            "--primary-catalog",
            args.primary_catalog,
            "--gold",
            args.gold,
            "--output-dir",
            args.output_dir,
        ]
    )
    summary = write_er_stress_benchmark(
        args.raw_dir,
        args.catalog,
        args.gold,
        args.output_dir,
        primary_catalog_path=args.primary_catalog,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
