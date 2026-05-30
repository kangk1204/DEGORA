"""Uniform baseline adapters and R preflight/blocker reporting.

Executable Tier 0 methods share the same long schema. Direct-prior-art methods
are never silently skipped: runnable adapters emit the same schema, while
missing packages, wrappers, or faithful DEG-table adaptations are written to the
failure ledger as ``open_s1_blocker`` rows.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues, norm

from .aggregate import _eligible_study_gene_stats, rank_product_consensus, slice_consensus
from .stats import bh_adjust

BASELINE_RESULT_COLUMNS = [
    "method_id",
    "setting_id",
    "gene_id",
    "symbol",
    "rank",
    "score",
    "pvalue",
    "padj",
    "effect",
    "direction",
    "n_studies",
    "missingness",
    "runtime_s",
    "version",
    "status",
]

FAILURE_LEDGER_COLUMNS = [
    "corpus",
    "method_id",
    "setting_id",
    "tier",
    "status",
    "blocker_id",
    "message",
    "resolution",
]

PARITY_MATRIX_COLUMNS = [
    "method_id",
    "method_family",
    "implementation_source",
    "version_or_commit",
    "setting_id",
    "setting_type",
    "parameters",
    "input_requirements",
    "output_schema_status",
    "run_status",
    "failure_mode",
    "equal_tuning_status",
    "claim_allowed",
    "notes",
]

PUBLIC_SUMMARY_TOOL_INPUT_COLUMNS = [
    "method_id",
    "method_family",
    "public_summary_deg_status",
    "can_run_from_current_public_files",
    "current_pipeline_status",
    "required_inputs",
    "current_public_inputs_available",
    "missing_or_nonfaithful_inputs",
    "faithful_adapter_decision",
    "manuscript_use",
    "source_or_package",
    "notes",
]

R_PVALUE_FLOOR = 1e-300
BASELINE_LFC_CAP = 10.0


def _resolve_manifest_artifact_path(artifact: object, baseline_dir: Path) -> Path:
    raw = str(artifact).strip()
    path = Path(raw)
    if path.is_absolute():
        return path

    candidates = [
        baseline_dir / path,
        Path.cwd() / path,
        _repo_root() / path,
        baseline_dir / path.name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def baseline_result_paths(baseline_dir: Path) -> list[Path]:
    """Return the current successful baseline result TSVs.

    Generated baseline directories can retain stale TSVs from previous corpus
    names. When the manifest exists, it is the authority for which baseline
    result files belong to the current run.
    """

    baseline_dir = Path(baseline_dir)
    manifest_path = baseline_dir / "baseline_manifest.csv"
    if not manifest_path.exists():
        return sorted(baseline_dir.glob("*.tsv"))

    manifest = pd.read_csv(manifest_path)
    required = {"artifact", "artifact_type", "status"}
    missing_columns = required.difference(manifest.columns)
    if missing_columns:
        raise ValueError(
            f"baseline manifest {manifest_path} missing required columns: {sorted(missing_columns)}"
        )

    rows = manifest.loc[
        manifest["artifact_type"].astype(str).eq("baseline_result")
        & manifest["status"].astype(str).eq("ok")
    ]
    paths: list[Path] = []
    seen: set[Path] = set()
    missing_artifacts: list[str] = []
    for artifact in rows["artifact"].dropna().tolist():
        path = _resolve_manifest_artifact_path(artifact, baseline_dir)
        if not path.exists():
            missing_artifacts.append(str(artifact))
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        paths.append(path)
    if missing_artifacts:
        raise FileNotFoundError(
            f"baseline manifest {manifest_path} references missing artifacts: {missing_artifacts}"
        )
    return paths



AWMETA_VARIANCE_COLUMNS = {
    "se",
    "stderr",
    "standard_error",
    "standarderror",
    "lfc_se",
    "lfcse",
    "logfc_se",
    "logfcse",
    "beta_se",
    "betase",
    "var",
    "variance",
    "vi",
    "weight",
    "weights",
}

HSTOUFFER_REQUIRED_ORIGINAL_COLUMNS = ["id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

HSTOUFFER_SOURCE_LEVEL_EVIDENCE = [
    {
        "source": "CR4CKID/hStouffer hStouffer.py",
        "commit": "306e38c26919f19e7c3dfd6cd646005c502b3310",
        "evidence": "The pinned script reads each input with pandas read_csv(..., sep='\\t'), uses the first column as gene id, column index 2 as log2FC, column index 5 as p-value, and the named lfcSE column for REM filtering.",
        "implication": "A faithful adapter needs original DESeq2-like tabular results retaining lfcSE and the canonical column order/fields, not a harmonized p-value-derived surrogate table.",
    },
    {
        "source": "CR4CKID/hStouffer hStouffer.py",
        "commit": "306e38c26919f19e7c3dfd6cd646005c502b3310",
        "evidence": "The REM stage squares lfcSE to form per-study variances and drops rows without positive SE/log2FC before deciding DEG direction.",
        "implication": "Imputing lfcSE from DEGORA signed_z or p-values would change the published method's variance model and is blocked.",
    },
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _resolve_source_path(source_path: Any) -> Path | None:
    if pd.isna(source_path):
        return None
    raw = str(source_path)
    if not raw:
        return None
    path = Path(raw)
    if path.exists():
        return path
    data_marker = "/data/"
    if data_marker in raw:
        candidate = _repo_root() / "data" / raw.split(data_marker, 1)[1]
        if candidate.exists():
            return candidate
    candidate = _repo_root() / raw
    if candidate.exists():
        return candidate
    return path


def _read_source_columns(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {"exists": False, "error": "missing source_path", "sheets": []}
    if not path.exists():
        return {"exists": False, "error": f"source file not found: {path}", "sheets": []}
    try:
        if path.suffix.lower() in {".xlsx", ".xls"}:
            workbook = pd.ExcelFile(path)
            sheets = []
            for sheet_name in workbook.sheet_names[:10]:
                header = pd.read_excel(workbook, sheet_name=sheet_name, nrows=0)
                sheets.append({"sheet": str(sheet_name), "columns": [str(column) for column in header.columns]})
            return {"exists": True, "error": "", "sheets": sheets}
        header = pd.read_csv(path, sep=None, engine="python", nrows=0, compression="infer")
        return {"exists": True, "error": "", "sheets": [{"sheet": "file", "columns": [str(column) for column in header.columns]}]}
    except Exception as exc:  # pragma: no cover - exact parser failures depend on optional engines and input files.
        return {"exists": True, "error": f"{type(exc).__name__}: {exc}", "sheets": []}


def _hstouffer_source_input_audit(harmonized: pd.DataFrame) -> list[dict[str, Any]]:
    required_lower = {column.lower(): column for column in HSTOUFFER_REQUIRED_ORIGINAL_COLUMNS}
    fields = ["study_id", "paper_id", "pipeline", "source_path", "source_url"]
    if not set(fields).issubset(harmonized.columns):
        return []
    records: list[dict[str, Any]] = []
    for row in harmonized[fields].drop_duplicates().sort_values("study_id").itertuples(index=False):
        row_dict = row._asdict()
        path = _resolve_source_path(row_dict["source_path"])
        header = _read_source_columns(path)
        all_columns: list[str] = []
        for sheet in header["sheets"]:
            all_columns.extend(sheet["columns"])
        lower_columns = {column.lower() for column in all_columns}
        present = [original for lowered, original in required_lower.items() if lowered in lower_columns]
        missing = [original for lowered, original in required_lower.items() if lowered not in lower_columns]
        pipeline = str(row_dict["pipeline"])
        status = "compatible_header_detected" if not missing and pipeline.lower() == "deseq2" else "blocked"
        if pipeline.lower() != "deseq2":
            status = "blocked_non_deseq2_pipeline"
        elif missing:
            status = "blocked_missing_required_columns"
        records.append(
            {
                "study_id": row_dict["study_id"],
                "paper_id": row_dict["paper_id"],
                "pipeline": pipeline,
                "source_path": str(row_dict["source_path"]),
                "source_url": row_dict["source_url"],
                "local_file_exists": bool(header["exists"]),
                "header_error": header["error"],
                "detected_sheets": header["sheets"],
                "required_columns_present": present,
                "required_columns_missing": missing,
                "audit_status": status,
            }
        )
    return records


AWMETA_EFFECT_COLUMNS = {
    "lfc",
    "logfc",
    "log2fc",
    "log2foldchange",
    "log2_fold_change",
    "effect",
    "estimate",
    "beta",
}


def _normalized_column_name(column: object) -> str:
    return str(column).strip().lower().replace(" ", "_").replace("-", "_").replace(".", "_").replace("(", "").replace(")", "")


def _resolve_awmeta_source_path(path_value: object) -> Path | None:
    if path_value is None or (isinstance(path_value, float) and math.isnan(path_value)):
        return None
    path = Path(str(path_value))
    if path.exists():
        return path
    marker = "/data/"
    text = str(path)
    if marker in text:
        relative = Path("data") / Path(text.split(marker, 1)[1])
        candidate = Path.cwd() / relative
        if candidate.exists():
            return candidate
    return path


def _source_header(path: Path) -> dict[str, Any]:
    """Read only source-table headers for an AWmeta input audit."""

    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith((".xlsx", ".xls")):
        workbook = pd.ExcelFile(path)
        sheet = workbook.sheet_names[0]
        columns = pd.read_excel(path, sheet_name=sheet, nrows=0).columns.tolist()
        return {"read_status": "ok", "sheet": sheet, "columns": [str(column) for column in columns]}
    if suffixes.endswith((".csv", ".csv.gz")):
        columns = pd.read_csv(path, nrows=0).columns.tolist()
        return {"read_status": "ok", "sheet": "", "columns": [str(column) for column in columns]}
    columns = pd.read_csv(path, nrows=0, sep=None, engine="python").columns.tolist()
    return {"read_status": "ok", "sheet": "", "columns": [str(column) for column in columns]}


def awmeta_source_input_audit(harmonized: pd.DataFrame) -> list[dict[str, Any]]:
    """Audit original source tables for faithful AWmeta variance/weight inputs.

    The audit is intentionally conservative: it only credits explicit original
    columns named like standard errors, variances, or weights. It does not infer
    standard errors from p-values, z-scores, test statistics, baseMean/logCPM, or
    replicate expression columns.
    """

    required = {"study_id", "source_path", "pipeline"}
    if not required.issubset(harmonized.columns):
        return []
    records: list[dict[str, Any]] = []
    study_sources = harmonized.loc[:, ["study_id", "pipeline", "source_path"]].drop_duplicates(["study_id", "source_path"])
    for row in study_sources.itertuples(index=False):
        source = _resolve_awmeta_source_path(row.source_path)
        record: dict[str, Any] = {
            "study_id": str(row.study_id),
            "pipeline": str(row.pipeline),
            "source_path": str(row.source_path),
            "resolved_source_exists": bool(source and source.exists()),
            "source_read_status": "not_found",
            "sheet": "",
            "columns_inspected": [],
            "effect_columns": [],
            "variance_or_weight_columns": [],
            "acquisition_requirement": "obtain original per-gene variance/SE or documented AWmeta-equivalent weights for this study",
        }
        if source and source.exists():
            try:
                header = _source_header(source)
                columns = list(header["columns"])
                normalized = [_normalized_column_name(column) for column in columns]
                exact = dict(zip(normalized, columns, strict=False))
                effect_columns = [str(exact[name]) for name in normalized if name in AWMETA_EFFECT_COLUMNS]
                variance_columns = [str(exact[name]) for name in normalized if name in AWMETA_VARIANCE_COLUMNS]
                record.update(
                    {
                        "source_read_status": str(header["read_status"]),
                        "sheet": str(header["sheet"]),
                        "columns_inspected": [str(column) for column in columns],
                        "effect_columns": effect_columns,
                        "variance_or_weight_columns": variance_columns,
                    }
                )
                if variance_columns:
                    record["acquisition_requirement"] = "variance/SE candidate present; verify it is the original per-gene within-study effect variance before AWmeta use"
            except Exception as exc:  # pragma: no cover - depends on optional engine support and source-table format
                record["source_read_status"] = f"read_failed:{type(exc).__name__}: {exc}"
        records.append(record)
    return records

def awmeta_deg_table_feasibility(harmonized: pd.DataFrame) -> dict[str, Any]:
    """Assess whether harmonized DEG tables can faithfully feed AWmeta/AW-REM.

    AWmeta's AW-REM branch is an effect-size meta-analysis. A faithful adapter
    needs per-study effect sizes plus within-study variances/standard errors
    (or an equivalent documented weight input). The DEGORA harmonized as-published
    DEG schema deliberately keeps p-values and log-fold changes from heterogeneous
    pipelines, but it does not currently preserve a statistically comparable
    variance column. Inferring SE as abs(lfc) / abs(z_from_p) would bake in
    pipeline-specific test-statistic assumptions and is therefore not a faithful
    prior-art comparator.
    """

    columns = {_normalized_column_name(column) for column in harmonized.columns}
    variance_columns = sorted(columns.intersection(AWMETA_VARIANCE_COLUMNS))
    required_basic = {"study_id", "gene_symbol", "lfc", "pvalue"}
    missing_basic = sorted(required_basic.difference(columns))
    source_audit = awmeta_source_input_audit(harmonized)
    studies_with_source_variance = sorted(
        {
            str(row["study_id"])
            for row in source_audit
            if row.get("variance_or_weight_columns")
        }
    )
    active_studies = sorted(harmonized["study_id"].dropna().astype(str).unique()) if "study_id" in harmonized.columns else []
    studies_with_harmonized_variance: list[str] = []
    if variance_columns and "study_id" in harmonized.columns:
        variance_frame = harmonized.loc[:, ["study_id", *[column for column in harmonized.columns if _normalized_column_name(column) in variance_columns]]].copy()
        variance_frame["_has_variance"] = variance_frame.drop(columns=["study_id"]).notna().any(axis=1)
        studies_with_harmonized_variance = sorted(
            variance_frame.loc[variance_frame["_has_variance"], "study_id"].dropna().astype(str).unique()
        )
    studies_with_any_variance = sorted(set(studies_with_source_variance).union(studies_with_harmonized_variance))
    studies_missing_source_variance = sorted(set(active_studies).difference(studies_with_any_variance))
    common = {
        "required_inputs": {
            "per_study_gene": "gene identifier/symbol",
            "effect": "original per-study log fold-change/effect estimate",
            "variance": "original per-study effect variance, standard error, or documented equivalent weight",
        },
        "prohibited_approximations": [
            "derive SE from p-values",
            "derive SE from signed_z",
            "mix raw logFC with heterogeneous pipeline test statistics as equivalent variances",
            "silently filter studies lacking variance/SE inputs",
        ],
        "active_studies": active_studies,
        "source_input_audit": source_audit,
        "studies_with_harmonized_variance_candidates": studies_with_harmonized_variance,
        "studies_with_source_variance_candidates": studies_with_source_variance,
        "studies_missing_source_variance_candidates": studies_missing_source_variance,
        "data_acquisition_requirement": (
            "For every active study, acquire the original per-gene effect-size variance/standard-error column "
            "or a published AWmeta-equivalent weight table tied to the same effect estimates; otherwise keep AWmeta/AW-REM blocked."
        ),
    }
    if missing_basic:
        return {
            "faithful": False,
            "blocker_id": "awmeta_required_columns_missing",
            "message": f"AWmeta DEG-table adapter requires {sorted(required_basic)} plus variance/SE fields; missing basic columns: {missing_basic}.",
            "variance_columns": variance_columns,
            "missing_basic_columns": missing_basic,
            **common,
        }
    if not variance_columns:
        return {
            "faithful": False,
            "blocker_id": "awmeta_variance_inputs_missing",
            "message": (
                "AWmeta/AW-REM requires per-study effect-size variances, standard errors, or documented equivalent weights. "
                "The harmonized as-published DEG tables provide study_id, gene_symbol, lfc, and pvalue but no uniform variance/SE field; "
                "deriving SE from p-values/signed z-scores would be a non-faithful approximation across heterogeneous DEG pipelines."
            ),
            "variance_columns": variance_columns,
            "missing_basic_columns": [],
            **common,
        }
    if studies_missing_source_variance:
        return {
            "faithful": False,
            "blocker_id": "awmeta_incomplete_original_variance_inputs",
            "message": (
                "Variance-like columns are present for only a subset of source tables. A faithful AWmeta/AW-REM run would "
                "silently filter the corpus unless original variance/SE or documented equivalent weights are acquired for "
                f"all active studies; missing studies: {studies_missing_source_variance}."
            ),
            "variance_columns": variance_columns,
            "missing_basic_columns": [],
            **common,
        }
    return {
        "faithful": False,
        "blocker_id": "awmeta_weight_optimization_unimplemented",
        "message": (
            "Variance-like columns are present, but a documented AW-Fisher/AW-REM adaptive-weight implementation "
            "and default weight initialization are still not implemented in this thin-slice codebase."
        ),
        "variance_columns": variance_columns,
        "missing_basic_columns": [],
        **common,
    }

BASELINE_MANIFEST_COLUMNS = [
    "artifact",
    "method_id",
    "setting_id",
    "artifact_type",
    "status",
    "source_command",
    "input_harmonized",
    "rows",
    "notes",
]


PREFLIGHT_PACKAGES = [
    "metafor",
    "MetaVolcanoR",
    "RobustRankAggreg",
    "AWFisher",
    "metaRNASeq",
    "RankProd",
    "MetaDE",
    "DExMA",
    "MetaIntegrator",
]

PREFLIGHT_PACKAGE_SOURCES = {
    "metafor": "CRAN/metafor",
    "MetaVolcanoR": "Bioconductor/MetaVolcanoR",
    "RobustRankAggreg": "CRAN/RobustRankAggreg",
    "AWFisher": "Bioconductor/AWFisher",
    "metaRNASeq": "R-universe/CRAN metaRNASeq",
    "RankProd": "Bioconductor/RankProd",
    "MetaDE": "CRAN archive/MetaDE",
    "DExMA": "Bioconductor/DExMA",
    "MetaIntegrator": "CRAN/MetaIntegrator",
}

HSTOUFFER_SOURCE_PIN = {
    "repo_url": "https://github.com/CR4CKID/hStouffer",
    "commit": "306e38c26919f19e7c3dfd6cd646005c502b3310",
    "script_path": "hStouffer/hStouffer.py",
    "license": "MIT",
    "required_input_shape": "directory of DESeq2 result .txt files with id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj",
    "published_default_command": "hStouffer.py -d <deseq2-result-dir> -o <prefix>",
}


@dataclass(frozen=True)
class MethodSpec:
    method_id: str
    tier: str
    method_family: str
    implementation_source: str
    input_requirements: str
    parameters: str
    setting_id: str = "default"
    setting_type: str = "default"
    r_package: str | None = None
    resolution: str = ""


R_BACKED_METHODS = [
    MethodSpec(
        "robustrankaggreg",
        "tier1",
        "direct_prior_art",
        "CRAN RobustRankAggreg via Rscript",
        "within-study ranked gene lists",
        "method=aggregateRanks; min_studies=2",
        r_package="RobustRankAggreg",
        resolution="Install R and RobustRankAggreg.",
    ),
    MethodSpec(
        "metavolcanor",
        "tier1",
        "direct_prior_art",
        "Bioconductor MetaVolcanoR via R/rpy2",
        "per-study gene, logFC, p-value tables with comparable effect definitions",
        "mode=p_comb default; min_studies=2",
        r_package="MetaVolcanoR",
        resolution="Install R, BiocManager, rpy2, and MetaVolcanoR.",
    ),
    MethodSpec(
        "awfisher",
        "tier1",
        "direct_prior_art",
        "Bioconductor AWFisher::AWFisher_pvalue via Rscript",
        "per-study p-value matrix; study contribution weights are inferred by AW-Fisher",
        "AWFisher_pvalue on p_from_signed_z; min_studies=2",
        r_package="AWFisher",
        resolution="Install R/Bioconductor AWFisher and run the uniform-schema adapter.",
    ),
    MethodSpec(
        "metarnaseq_fisher",
        "tier1",
        "direct_prior_art",
        "metaRNASeq::fishercomb via Rscript",
        "per-study p-value vectors for each gene",
        "fishercomb on p_from_signed_z; min_studies=2",
        r_package="metaRNASeq",
        resolution="Install metaRNASeq and run the Fisher combination adapter.",
    ),
    MethodSpec(
        "metarnaseq_invnorm",
        "tier1",
        "direct_prior_art",
        "metaRNASeq::invnorm via Rscript",
        "per-study p-value vectors plus study replicate-count weights",
        "invnorm on p_from_signed_z with source-unit replicate weights and neutral p=0.5 for missing study-gene entries; min_studies=2",
        r_package="metaRNASeq",
        resolution="Install metaRNASeq and run the inverse-normal combination adapter.",
    ),
    MethodSpec(
        "rankprod_exact",
        "tier2",
        "direct_prior_art",
        "Bioconductor RankProd official package",
        "replicate-level expression/rank-product inputs with class labels and origin metadata",
        "official RankProd defaults; min_studies=2",
        r_package="RankProd",
        resolution="Acquire expression matrices/class labels or keep the current rank_product_approx explicitly labeled as a summary-rank approximation.",
    ),
    MethodSpec(
        "metade",
        "tier2",
        "direct_prior_art",
        "MetaDE R package",
        "p-value matrix for p-value modes or effect-size plus variance/raw-expression inputs for ES/FEM/REM modes",
        "MetaDE default methods where inputs are faithful; min_studies=2",
        r_package="MetaDE",
        resolution="Install MetaDE; run only p-value/rank modes from public DEG summaries and require variance/raw expression for effect-size modes.",
    ),
    MethodSpec(
        "dexma",
        "tier2",
        "direct_prior_art",
        "Bioconductor DExMA workflow",
        "expression matrices or GEO accessions with sample phenotype metadata for differential-expression meta-analysis",
        "DExMA workflow defaults; min_studies=2",
        r_package="DExMA",
        resolution="Acquire raw/processed expression matrices and phenotype metadata, or document DExMA as out-of-scope for summary-only public DEG files.",
    ),
    MethodSpec(
        "metaintegrator",
        "tier2",
        "direct_prior_art",
        "CRAN MetaIntegrator",
        "study-level expression matrices, phenotype labels, effect sizes and variances for disease signature meta-analysis",
        "MetaIntegrator defaults; min_studies=2",
        r_package="MetaIntegrator",
        resolution="Acquire expression matrices and phenotype labels; do not run from already-filtered public DEG summaries.",
    ),
    MethodSpec(
        "hstouffer",
        "tier2",
        "direct_prior_art",
        "published CR4CKID/hStouffer Python script",
        "uniform-pipeline study p-values/effects; dynamic cutoff parameters",
        "published defaults; min_studies=2",
        resolution="Pin the CR4CKID/hStouffer Python implementation and materialize compatible DESeq2-like inputs.",
    ),
    MethodSpec(
        "awmeta",
        "tier2",
        "direct_prior_art",
        "AWmeta/AW-Fisher + AW-REM faithful implementation from preprint",
        "per-study p-values, effects, variance/weight inputs with documented initialization",
        "published defaults; min_studies=2",
        r_package="metafor",
        resolution="Implement AW-Fisher/AW-REM from the preprint with documented choices and install metafor dependency.",
    ),
]

TIER0_METHOD_METADATA = {
    "degora_slice": {
        "method_family": "degora_locked_reference",
        "implementation_source": "repo degora.aggregate.slice_consensus",
        "parameters": "locked project defaults; min_studies={min_studies}",
        "input_requirements": "harmonized DEGORA DEG table",
        "version_or_commit": "degora-thin-slice",
        "setting_id": "locked",
    },
    "weighted_stouffer": {
        "method_family": "classical",
        "implementation_source": "scipy.stats.norm",
        "parameters": "sample-size weighted signed-z Stouffer; min_studies={min_studies}",
        "input_requirements": "per-study signed z-scores and sample sizes",
        "version_or_commit": "scipy.stats.norm",
        "setting_id": "default",
    },
    "unweighted_stouffer": {
        "method_family": "classical",
        "implementation_source": "scipy.stats.norm",
        "parameters": "unweighted signed-z Stouffer; min_studies={min_studies}",
        "input_requirements": "per-study signed z-scores",
        "version_or_commit": "scipy.stats.norm",
        "setting_id": "default",
    },
    "fisher": {
        "method_family": "classical",
        "implementation_source": "scipy.stats.combine_pvalues",
        "parameters": "Fisher combined p-values from signed-z p-values; min_studies={min_studies}",
        "input_requirements": "per-study p-values",
        "version_or_commit": "scipy.stats.combine_pvalues",
        "setting_id": "default",
    },
    "rank_product_approx": {
        "method_family": "classical_rank_approximation",
        "implementation_source": "repo degora.aggregate.rank_product_consensus",
        "parameters": "rank-product approximation; min_studies={min_studies}",
        "input_requirements": "within-study normalized ranks",
        "version_or_commit": "degora rank-product approximation",
        "setting_id": "default",
    },
    "sign_vote": {
        "method_family": "classical",
        "implementation_source": "repo sign-test approximation",
        "parameters": "two-sided sign vote; min_studies={min_studies}",
        "input_requirements": "per-study signed z-score sign calls",
        "version_or_commit": "two-sided sign-test approximation",
        "setting_id": "default",
    },
}


DIRECT_METHOD_METADATA = {
    "metavolcanor": {
        "method_family": "direct_prior_art",
        "implementation_source": "Bioconductor MetaVolcanoR::combining_mv via Rscript",
        "parameters": "pcriteria=pvalue; foldchangecol=Log2FC; metafc=Mean; metathr=0.01; min_studies={min_studies}",
        "input_requirements": "per-study Symbol, Log2FC, pvalue tables from harmonized DEG rows",
        "version_or_commit": "{version}",
        "setting_id": "default",
    },
    "robustrankaggreg": {
        "method_family": "direct_prior_art",
        "implementation_source": "CRAN RobustRankAggreg::aggregateRanks via Rscript",
        "parameters": "method=RRA; ranked by within-study normalized rank; N=union genes; min_studies={min_studies}",
        "input_requirements": "within-study ranked gene lists from harmonized DEG tables",
        "version_or_commit": "{version}",
        "setting_id": "default",
    },
    "awfisher": {
        "method_family": "direct_prior_art",
        "implementation_source": "Bioconductor AWFisher::AWFisher_pvalue via Rscript",
        "parameters": "AWFisher_pvalue on two-sided p_from_signed_z matrix; missing study-gene entries set to p=1; min_studies={min_studies}",
        "input_requirements": "per-source-unit gene p-values from harmonized DEG tables",
        "version_or_commit": "{version}",
        "setting_id": "default",
    },
    "metarnaseq_fisher": {
        "method_family": "direct_prior_art",
        "implementation_source": "metaRNASeq::fishercomb via Rscript",
        "parameters": "fishercomb on two-sided p_from_signed_z vectors; missing study-gene entries set to p=1; min_studies={min_studies}",
        "input_requirements": "per-source-unit gene p-values from harmonized DEG tables",
        "version_or_commit": "{version}",
        "setting_id": "default",
    },
    "metarnaseq_invnorm": {
        "method_family": "direct_prior_art",
        "implementation_source": "metaRNASeq::invnorm via Rscript",
        "parameters": "invnorm on two-sided p_from_signed_z vectors using source-unit sample-size weights; missing study-gene entries set to neutral p=0.5; min_studies={min_studies}",
        "input_requirements": "per-source-unit gene p-values plus source-unit sample-size weights from harmonized DEG tables",
        "version_or_commit": "{version}",
        "setting_id": "default",
    },
}

INVNORM_UNINFORMATIVE_BLOCKER_ID = "metarnaseq_invnorm_uninformative_sparse_public_summary"
INVNORM_UNINFORMATIVE_PREFIX = f"{INVNORM_UNINFORMATIVE_BLOCKER_ID}:"

SUMMARY_ONLY_BLOCKERS = {
    "rankprod_exact": {
        "blocker_id": "rankprod_exact_requires_expression_or_origin_labels",
        "message": (
            "The official RankProd workflow is an expression/rank-product analysis with replicate/origin labels. "
            "DEGORA can run a summary-rank approximation from public DEG tables, but treating it as official RankProd would be non-faithful."
        ),
    },
    "metade": {
        "blocker_id": "metade_effect_size_modes_require_variance_or_raw_expression",
        "message": (
            "MetaDE p-value/rank ideas are represented by runnable Fisher/Stouffer/AWFisher/rank baselines, but MetaDE effect-size/FEM/REM modes "
            "need per-study variances or raw expression-derived effect-size inputs that are absent from heterogeneous public DEG summaries."
        ),
    },
    "dexma": {
        "blocker_id": "dexma_requires_expression_matrices_and_phenotype_metadata",
        "message": (
            "DExMA is a gene-expression meta-analysis workflow over expression matrices/GEO-derived objects and sample phenotype metadata; "
            "already-filtered public supplementary DEG tables are not sufficient for a faithful DExMA run."
        ),
    },
    "metaintegrator": {
        "blocker_id": "metaintegrator_requires_expression_and_phenotype_objects",
        "message": (
            "MetaIntegrator is designed around study expression data, phenotype labels, and effect-size/variance objects for signature meta-analysis; "
            "public DEG-only tables do not contain the required cohort-level inputs."
        ),
    },
}

PUBLIC_WORKFLOW_PRIOR_ART_ROWS = [
    {
        "method_id": "omicc",
        "method_family": "public expression reuse / meta-analysis platform",
        "public_status": "workflow_level_public_expression_platform",
        "can_run": "not_from_current_deg_summary_files",
        "required": "public expression studies, sample-group annotations, expression matrices, and OMiCC platform curation",
        "missing": "sample-level expression matrices and phenotype/group labels are not present in already-filtered public DEG summary tables",
        "decision": "cite and compare qualitatively as a public-expression reuse platform; do not treat as a same-input gene-rank comparator",
        "manuscript_use": "workflow/resource prior art",
        "source": "OMiCC / MetaIntegrator / RankProd public expression platform",
        "notes": "OMiCC supports public-data reuse and meta-analysis from expression studies; DEGORA instead starts from heterogeneous as-published DEG evidence rows.",
    },
    {
        "method_id": "imageo",
        "method_family": "GEO expression meta-analysis web workflow",
        "public_status": "workflow_level_geo_reanalysis_platform",
        "can_run": "not_from_current_deg_summary_files",
        "required": "GEO dataset identifiers, platform expression matrices, sample labels, and GEO-derived preprocessing",
        "missing": "already-filtered supplementary DEG tables do not contain the expression objects and sample metadata required by the ImaGEO workflow",
        "decision": "cite and compare qualitatively; use as a contrast to DEGORA's as-published DEG-table operating regime",
        "manuscript_use": "workflow/resource prior art",
        "source": "ImaGEO web tool",
        "notes": "ImaGEO is a GEO-starting workflow, while DEGORA is a summary-evidence ingestion and audit workflow.",
    },
    {
        "method_id": "networkanalyst",
        "method_family": "statistical/visual/network expression meta-analysis platform",
        "public_status": "workflow_level_expression_or_gene_table_platform",
        "can_run": "not_as_same_input_gene_rank_comparator",
        "required": "normalized expression tables or compatible multi-table expression inputs, sample labels, and downstream network-analysis settings",
        "missing": "DEGORA's benchmark contract supplies DEG evidence rows, not complete expression matrices for NetworkAnalyst's statistical workflow",
        "decision": "cite and compare as a web workflow; not a replacement for the current public DEG-summary comparator table",
        "manuscript_use": "workflow/resource prior art",
        "source": "NetworkAnalyst",
        "notes": "Relevant for user-facing expression meta-analysis and downstream networks; orthogonal to DEGORA's source-unit evidence-card output.",
    },
    {
        "method_id": "crossmeta",
        "method_family": "GEO retrieval / cross-platform microarray meta-analysis package",
        "public_status": "workflow_level_raw_or_processed_expression_reanalysis",
        "can_run": "not_from_current_deg_summary_files",
        "required": "GEO accessions, expression matrices, platform annotation, normalization, and phenotype labels",
        "missing": "public DEG-only files omit the expression-level inputs crossmeta downloads and normalizes before meta-analysis",
        "decision": "cite as an alternative expression-reanalysis workflow; not a same-input comparator",
        "manuscript_use": "workflow/resource prior art",
        "source": "Bioconductor crossmeta",
        "notes": "crossmeta avoids as-published DEG heterogeneity by returning to GEO expression data; DEGORA accepts the heterogeneity as input.",
    },
    {
        "method_id": "deet",
        "method_family": "uniformly processed differential-expression atlas",
        "public_status": "resource_prior_art_not_same_input",
        "can_run": "not_as_method_comparator",
        "required": "user gene list compared against a uniformly processed DEG atlas",
        "missing": "DEET does not ingest arbitrary heterogeneous supplementary DEG tables to build a topic-specific source-auditable DB",
        "decision": "cite as resource prior art; differentiate uniform precomputed atlas from DEGORA's user-built evidence database",
        "manuscript_use": "database/resource prior art",
        "source": "DEET differential expression atlas",
        "notes": "Important nearest resource prior art; not a per-corpus meta-analysis baseline under the current input contract.",
    },
    {
        "method_id": "creeds",
        "method_family": "crowd-extracted expression signature database",
        "public_status": "resource_prior_art_not_same_input",
        "can_run": "not_as_method_comparator",
        "required": "curated or crowd-extracted up/down signatures from GEO studies",
        "missing": "CREEDS is a signature resource rather than a same-input source-unit meta-analysis method for newly collected DEG tables",
        "decision": "cite as expression-signature database prior art; compare output philosophy and provenance limits qualitatively",
        "manuscript_use": "database/resource prior art",
        "source": "CREEDS / Ma'ayan Lab expression signatures",
        "notes": "Closest resource-style threat to the database framing; DEGORA's differentiator is local, topic-specific, source-unit-auditable DEG evidence construction.",
    },
    {
        "method_id": "generic_pvalue_combiners",
        "method_family": "generic p-value combination packages",
        "public_status": "method_family_covered_by_runnable_baselines",
        "can_run": "yes_in_principle",
        "required": "per-study p-values, optional weights, and method-specific independence/dependence assumptions",
        "missing": "no missing input for basic Fisher/Stouffer-style families; package-specific variants are redundant with the covered family-level baselines",
        "decision": "cover by Fisher, weighted/unweighted Stouffer, AWFisher, and metaRNASeq rows rather than adding many redundant p-value packages",
        "manuscript_use": "method-family coverage row",
        "source": "metapod / metap / TFisher class methods",
        "notes": "Use this row to preempt reviewer requests for generic p-value packages; add package-specific runs only if a reviewer requests a distinct test family.",
    },
]


def _method_metadata(method_id: str) -> dict[str, str]:
    if method_id in TIER0_METHOD_METADATA:
        return TIER0_METHOD_METADATA[method_id]
    return DIRECT_METHOD_METADATA[method_id]


def validate_baseline_result(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and order a result frame according to the PRD baseline schema."""

    missing = [column for column in BASELINE_RESULT_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"baseline result missing required columns: {missing}")
    out = frame.copy()
    for column in ["method_id", "setting_id", "gene_id", "symbol", "direction", "version", "status"]:
        out[column] = out[column].astype("string")
    for column in ["rank", "n_studies"]:
        out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
    for column in ["score", "pvalue", "padj", "effect", "missingness", "runtime_s"]:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    bad_status = set(out["status"].dropna().astype(str)).difference({"ok", "failed", "open_s1_blocker"})
    if bad_status:
        raise ValueError(f"unsupported baseline status values: {sorted(bad_status)}")
    return out[BASELINE_RESULT_COLUMNS]


def _study_level(harmonized: pd.DataFrame, min_studies: int) -> pd.DataFrame:
    by_study = _eligible_study_gene_stats(harmonized, min_studies)
    if by_study.empty:
        return by_study
    by_study = by_study.copy()
    lfc = pd.to_numeric(by_study["lfc"], errors="coerce") if "lfc" in by_study.columns else pd.Series(np.nan, index=by_study.index)
    signed_z = pd.to_numeric(by_study["signed_z"], errors="coerce") if "signed_z" in by_study.columns else pd.Series(np.nan, index=by_study.index)
    by_study["lfc"] = np.where(
        np.isfinite(lfc),
        lfc.clip(lower=-BASELINE_LFC_CAP, upper=BASELINE_LFC_CAP),
        np.sign(signed_z).fillna(0.0) * BASELINE_LFC_CAP,
    )
    by_study = by_study.copy()
    by_study["p_from_signed_z"] = 2.0 * norm.sf(np.abs(by_study["signed_z"].to_numpy(dtype=float)))
    return by_study


def _study_matrix_inputs(
    harmonized: pd.DataFrame,
    min_studies: int,
    *,
    missing_pvalue_fill: float = 1.0,
) -> dict[str, Any]:
    """Build gene-by-source-unit matrices for R comparators that combine p-values."""

    by_study = _study_level(harmonized, min_studies)
    if by_study.empty:
        return {
            "by_study": by_study,
            "p_matrix": pd.DataFrame(),
            "lfc_matrix": pd.DataFrame(),
            "stats": pd.DataFrame(columns=["gene_symbol", "n_studies", "effect"]),
            "studies": [],
            "replicate_counts": [],
            "total_studies": 0,
        }
    required = {"study_id", "gene_symbol", "lfc", "p_from_signed_z", "weight"}
    missing = required.difference(by_study.columns)
    if missing:
        raise ValueError(f"matrix comparator input missing required columns: {sorted(missing)}")

    matrix_input = by_study.loc[:, ["study_id", "gene_symbol", "lfc", "signed_z", "p_from_signed_z", "weight"]].copy()
    matrix_input["study_id"] = matrix_input["study_id"].astype(str)
    matrix_input["gene_symbol"] = matrix_input["gene_symbol"].astype(str)
    for column in ["lfc", "signed_z", "p_from_signed_z", "weight"]:
        matrix_input[column] = pd.to_numeric(matrix_input[column], errors="coerce")
    matrix_input = matrix_input.dropna(subset=["study_id", "gene_symbol", "p_from_signed_z"])
    matrix_input = matrix_input.loc[np.isfinite(matrix_input["p_from_signed_z"].to_numpy(dtype=float))].copy()
    matrix_input["p_from_signed_z"] = matrix_input["p_from_signed_z"].clip(lower=R_PVALUE_FLOOR, upper=1.0)
    matrix_input = matrix_input.sort_values(["gene_symbol", "study_id", "p_from_signed_z"])
    matrix_input = matrix_input.drop_duplicates(["gene_symbol", "study_id"], keep="first")
    if matrix_input.empty:
        return {
            "by_study": matrix_input,
            "p_matrix": pd.DataFrame(),
            "lfc_matrix": pd.DataFrame(),
            "stats": pd.DataFrame(columns=["gene_symbol", "n_studies", "effect"]),
            "studies": [],
            "replicate_counts": [],
            "total_studies": 0,
        }

    studies = sorted(matrix_input["study_id"].unique())
    genes = sorted(matrix_input["gene_symbol"].unique())
    fill_value = float(missing_pvalue_fill)
    if not np.isfinite(fill_value) or fill_value <= 0.0 or fill_value > 1.0:
        raise ValueError(f"missing_pvalue_fill must be in (0, 1], got {missing_pvalue_fill!r}")
    p_matrix = (
        matrix_input.pivot(index="gene_symbol", columns="study_id", values="p_from_signed_z")
        .reindex(index=genes, columns=studies)
        .astype(float)
        .fillna(fill_value)
        .clip(lower=R_PVALUE_FLOOR, upper=1.0)
    )
    lfc_matrix = (
        matrix_input.pivot(index="gene_symbol", columns="study_id", values="lfc")
        .reindex(index=genes, columns=studies)
        .astype(float)
    )
    stats = matrix_input.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        effect=("lfc", "mean"),
        mean_signed_z=("signed_z", "mean"),
        mean_weight=("weight", "mean"),
    )
    stats["effect"] = stats["effect"].fillna(stats["mean_signed_z"]).fillna(0.0)
    replicate_counts = (
        matrix_input.groupby("study_id")["weight"]
        .median()
        .reindex(studies)
        .pow(2)
        .replace([np.inf, -np.inf], np.nan)
        .fillna(1.0)
        .round()
        .clip(lower=1)
        .astype(int)
        .tolist()
    )
    return {
        "by_study": matrix_input,
        "p_matrix": p_matrix,
        "lfc_matrix": lfc_matrix,
        "stats": stats,
        "studies": studies,
        "replicate_counts": replicate_counts,
        "total_studies": len(studies),
    }


def _write_matrix_with_gene_column(frame: pd.DataFrame, path: Path) -> None:
    out = frame.reset_index().rename(columns={"gene_symbol": "gene_symbol"})
    out.to_csv(path, index=False)


def _direction(effect: pd.Series) -> pd.Series:
    return pd.Series(
        np.select([effect.gt(0), effect.lt(0)], ["up", "down"], default="zero"),
        index=effect.index,
        dtype="string",
    )


def _finalize(
    frame: pd.DataFrame,
    *,
    method_id: str,
    setting_id: str,
    effect: str | pd.Series,
    pvalue: str | pd.Series,
    n_studies: str | pd.Series,
    runtime_s: float,
    version: str,
    tie_breaker: str | pd.Series | None = None,
    tie_breaker_desc: bool = False,
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    out = pd.DataFrame(
        {
            "method_id": method_id,
            "setting_id": setting_id,
            "gene_id": frame["gene_symbol"].astype("string"),
            "symbol": frame["gene_symbol"].astype("string"),
            "effect": frame[effect] if isinstance(effect, str) else effect,
            "pvalue": frame[pvalue] if isinstance(pvalue, str) else pvalue,
            "n_studies": frame[n_studies] if isinstance(n_studies, str) else n_studies,
            "missingness": 0.0,
            "runtime_s": float(runtime_s),
            "version": version,
            "status": "ok",
        }
    )
    out["padj"] = bh_adjust(out["pvalue"].to_numpy(dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        out["score"] = -np.log10(out["padj"].clip(lower=np.nextafter(0.0, 1.0)))
    out["direction"] = _direction(out["effect"])
    sort_columns = ["padj", "pvalue", "symbol"]
    if tie_breaker is not None:
        tie_values = frame[tie_breaker] if isinstance(tie_breaker, str) else tie_breaker
        out["_tie_breaker"] = pd.to_numeric(tie_values, errors="coerce")
        if tie_breaker_desc:
            out["_tie_breaker"] = -out["_tie_breaker"]
        sort_columns = ["padj", "pvalue", "_tie_breaker", "symbol"]
    out = out.sort_values(sort_columns, na_position="last").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1, dtype=int)
    if "_tie_breaker" in out.columns:
        out = out.drop(columns=["_tie_breaker"])
    return validate_baseline_result(out)


def _metarnaseq_invnorm_uninformative_result(result: pd.DataFrame) -> bool:
    """Detect the sparse-summary failure mode where invnorm emits no rank signal."""

    if result.empty or not {"pvalue", "test_statistic"}.issubset(result.columns):
        return False
    pvalues = pd.to_numeric(result["pvalue"], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    test_stats = (
        pd.to_numeric(result["test_statistic"], errors="coerce")
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .abs()
    )
    if pvalues.empty or test_stats.empty:
        return False
    all_neutral_p = bool(pvalues.ge(1.0 - 1e-12).all())
    no_statistic_signal = bool(test_stats.le(1e-12).all() or test_stats.nunique(dropna=True) <= 1)
    if all_neutral_p and no_statistic_signal:
        return True

    raw = result.copy()
    raw["_pvalue"] = pd.to_numeric(raw["pvalue"], errors="coerce")
    raw["_test_statistic"] = pd.to_numeric(raw["test_statistic"], errors="coerce")
    top_n = min(100, len(raw))
    if top_n < 10:
        return False
    top = raw.sort_values(["_pvalue", "gene_symbol"], na_position="last").head(top_n)
    top_p = top["_pvalue"].replace([np.inf, -np.inf], np.nan)
    top_stat = top["_test_statistic"].replace([np.inf, -np.inf], np.nan)
    saturated_p = top_p.isna() | top_p.le(0.0) | top_p.le(R_PVALUE_FLOOR)
    if top_n >= 100 and saturated_p.all():
        return True
    finite_stat = top_stat.dropna()
    no_top_tie_breaker = finite_stat.empty or finite_stat.abs().le(1e-12).all() or finite_stat.nunique(dropna=True) <= 1
    return bool(saturated_p.all() and no_top_tie_breaker)


def degora_slice_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    consensus = slice_consensus(harmonized, min_studies=min_studies)
    runtime_s = time.perf_counter() - start
    return _finalize(
        consensus,
        method_id="degora_slice",
        setting_id="locked",
        effect="stouffer_z",
        pvalue="stouffer_p",
        n_studies="n_studies",
        runtime_s=runtime_s,
        version="degora-thin-slice",
    )


def weighted_stouffer_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    by_study["_wz"] = by_study["weight"] * by_study["signed_z"]
    by_study["_w2"] = by_study["weight"] ** 2
    out = by_study.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        sum_wz=("_wz", "sum"),
        sum_w2=("_w2", "sum"),
    )
    out["effect"] = out["sum_wz"] / np.sqrt(out["sum_w2"])
    out["pvalue"] = 2.0 * norm.sf(np.abs(out["effect"]))
    return _finalize(out, method_id="weighted_stouffer", setting_id="default", effect="effect", pvalue="pvalue", n_studies="n_studies", runtime_s=time.perf_counter() - start, version="scipy.stats.norm")


def unweighted_stouffer_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    out = by_study.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        sum_z=("signed_z", "sum"),
    )
    out["effect"] = out["sum_z"] / np.sqrt(out["n_studies"])
    out["pvalue"] = 2.0 * norm.sf(np.abs(out["effect"]))
    return _finalize(out, method_id="unweighted_stouffer", setting_id="default", effect="effect", pvalue="pvalue", n_studies="n_studies", runtime_s=time.perf_counter() - start, version="scipy.stats.norm")


def fisher_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    records: list[dict[str, Any]] = []
    for gene, gene_rows in by_study.groupby("gene_symbol", sort=False):
        pvalues = gene_rows["p_from_signed_z"].clip(lower=np.nextafter(0.0, 1.0), upper=1.0).to_numpy(dtype=float)
        _, pvalue = combine_pvalues(pvalues, method="fisher")
        records.append({"gene_symbol": gene, "n_studies": int(gene_rows["study_id"].nunique()), "effect": float(gene_rows["lfc"].mean()), "pvalue": float(pvalue)})
    return _finalize(pd.DataFrame.from_records(records), method_id="fisher", setting_id="default", effect="effect", pvalue="pvalue", n_studies="n_studies", runtime_s=time.perf_counter() - start, version="scipy.stats.combine_pvalues")


def rank_product_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    ranks = rank_product_consensus(harmonized, min_studies=min_studies)
    if ranks.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    out = ranks.rename(columns={"n_studies_rank": "n_studies", "rank_score": "effect"}).copy()
    out["pvalue"] = np.minimum(1.0, out["rank_product"].to_numpy(dtype=float))
    return _finalize(out, method_id="rank_product_approx", setting_id="default", effect="effect", pvalue="pvalue", n_studies="n_studies", runtime_s=time.perf_counter() - start, version="degora rank-product approximation")


def sign_vote_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    records: list[dict[str, Any]] = []
    for gene, gene_rows in by_study.groupby("gene_symbol", sort=False):
        signs = np.sign(gene_rows["signed_z"].to_numpy(dtype=float))
        positive = int(np.sum(signs > 0))
        negative = int(np.sum(signs < 0))
        n = int(len(signs))
        effect = (positive - negative) / n if n else 0.0
        tail = min(sum(math.comb(n, k) for k in range(positive, n + 1)), sum(math.comb(n, k) for k in range(negative, n + 1)))
        pvalue = min(1.0, 2.0 * tail / (2**n))
        records.append({"gene_symbol": gene, "n_studies": n, "effect": effect, "pvalue": pvalue})
    return _finalize(pd.DataFrame.from_records(records), method_id="sign_vote", setting_id="default", effect="effect", pvalue="pvalue", n_studies="n_studies", runtime_s=time.perf_counter() - start, version="two-sided sign-test approximation")


def _r_subprocess_env() -> dict[str, str]:
    """Return an R subprocess env with adjacent conda shared libraries visible."""

    env = os.environ.copy()
    rscript = shutil.which("Rscript")
    lib_dirs: list[str] = []
    if rscript:
        prefix = Path(rscript).resolve().parents[1]
        for candidate in [prefix / "lib", *prefix.glob("pkgs/*/lib")]:
            if candidate.is_dir() and any(candidate.glob("libfftw3*.so*")):
                lib_dirs.append(str(candidate))
    existing = env.get("LD_LIBRARY_PATH", "")
    if lib_dirs:
        parts = [*lib_dirs, *([existing] if existing else [])]
        env["LD_LIBRARY_PATH"] = ":".join(dict.fromkeys(parts))
    return env


def _run_rscript(args: list[str], *, timeout: int) -> subprocess.CompletedProcess[str]:
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript executable not found on PATH")
    return subprocess.run([rscript, *args], check=False, capture_output=True, text=True, timeout=timeout, env=_r_subprocess_env())


def _r_package_version(package: str) -> str:
    if not shutil.which("Rscript"):
        raise RuntimeError("Rscript executable not found on PATH")
    r_code = "pkg <- commandArgs(trailingOnly=TRUE)[1]; cat(as.character(packageVersion(pkg)))"
    completed = _run_rscript(["-e", r_code, package], timeout=15)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(detail or f"could not query R package version for {package}")
    return completed.stdout.strip()


def robustrankaggreg_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    """Run RobustRankAggreg on per-study ranked gene lists via Rscript."""

    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    required = {"study_id", "gene_symbol", "normalized_rank", "signed_z", "lfc"}
    missing = required.difference(by_study.columns)
    if missing:
        raise ValueError(f"RobustRankAggreg input missing required columns: {sorted(missing)}")

    ranked = by_study.loc[:, ["study_id", "gene_symbol", "normalized_rank", "signed_z", "lfc"]].copy()
    ranked = ranked.dropna(subset=["study_id", "gene_symbol", "normalized_rank"])
    ranked = ranked.sort_values(["study_id", "normalized_rank", "gene_symbol"])
    ranked = ranked.drop_duplicates(["study_id", "gene_symbol"], keep="first")
    if ranked.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    ranked["rank_order"] = ranked.groupby("study_id").cumcount() + 1
    n_universe = int(ranked["gene_symbol"].nunique())

    r_code = """
    args <- commandArgs(trailingOnly=TRUE)
    input <- args[1]
    output <- args[2]
    n_universe <- as.integer(args[3])
    ranks <- read.csv(input, stringsAsFactors=FALSE)
    ranks <- ranks[order(ranks$study_id, ranks$rank_order, ranks$gene_symbol), ]
    glist <- split(ranks$gene_symbol, ranks$study_id)
    result <- RobustRankAggreg::aggregateRanks(glist=glist, N=n_universe, method="RRA")
    write.csv(result, output, row.names=FALSE)
    """
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript executable not found on PATH")
    with tempfile.TemporaryDirectory(prefix="degora_rra_") as tmp:
        tmp_path = Path(tmp)
        rank_path = tmp_path / "ranks.csv"
        result_path = tmp_path / "rra.csv"
        ranked.loc[:, ["study_id", "gene_symbol", "rank_order"]].to_csv(rank_path, index=False)
        completed = _run_rscript(["-e", r_code, str(rank_path), str(result_path), str(n_universe)], timeout=120)
        if completed.returncode != 0:
            detail = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)
            raise RuntimeError(detail or "RobustRankAggreg aggregateRanks failed")
        result = pd.read_csv(result_path)

    name_column = "Name" if "Name" in result.columns else result.columns[0]
    score_column = "Score" if "Score" in result.columns else result.columns[-1]
    stats = ranked.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        effect=("signed_z", "mean"),
        mean_lfc=("lfc", "mean"),
    )
    out = result.loc[:, [name_column, score_column]].rename(columns={name_column: "gene_symbol", score_column: "pvalue"})
    out = out.merge(stats, on="gene_symbol", how="left")
    out["effect"] = out["effect"].fillna(out["mean_lfc"])
    version = _r_package_version("RobustRankAggreg")
    return _finalize(
        out,
        method_id="robustrankaggreg",
        setting_id="default",
        effect="effect",
        pvalue="pvalue",
        n_studies="n_studies",
        runtime_s=time.perf_counter() - start,
        version=f"RobustRankAggreg {version}",
    )


def metavolcanor_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    """Run MetaVolcanoR's Fisher p-value combining method on harmonized DEG tables."""

    start = time.perf_counter()
    by_study = _study_level(harmonized, min_studies)
    if by_study.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)
    required = {"study_id", "gene_symbol", "lfc", "p_from_signed_z"}
    missing = required.difference(by_study.columns)
    if missing:
        raise ValueError(f"MetaVolcanoR input missing required columns: {sorted(missing)}")

    meta_input = by_study.loc[:, ["study_id", "gene_symbol", "lfc", "p_from_signed_z"]].copy()
    meta_input = meta_input.dropna(subset=["study_id", "gene_symbol", "lfc", "p_from_signed_z"])
    meta_input["p_from_signed_z"] = meta_input["p_from_signed_z"].clip(lower=np.nextafter(0.0, 1.0), upper=1.0)
    meta_input = meta_input.loc[
        np.isfinite(meta_input["lfc"].to_numpy(dtype=float))
        & np.isfinite(meta_input["p_from_signed_z"].to_numpy(dtype=float))
    ].copy()
    meta_input = meta_input.sort_values(["study_id", "gene_symbol", "p_from_signed_z"])
    meta_input = meta_input.drop_duplicates(["study_id", "gene_symbol"], keep="first")
    if meta_input.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)

    total_studies = int(meta_input["study_id"].nunique())
    stats = meta_input.groupby("gene_symbol", as_index=False).agg(n_studies=("study_id", "nunique"))

    r_code = r"""
    args <- commandArgs(trailingOnly=TRUE)
    input <- args[1]
    output <- args[2]
    ranks <- read.csv(input, stringsAsFactors=FALSE)
    diffexp <- split(ranks[, c("Symbol", "Log2FC", "pvalue")], ranks$study_id)
    result <- MetaVolcanoR::combining_mv(
      diffexp=diffexp,
      pcriteria="pvalue",
      foldchangecol="Log2FC",
      genenamecol="Symbol",
      geneidcol="Symbol",
      metafc="Mean",
      metathr=0.01,
      collaps=FALSE,
      jobname="degora_metavolcanor",
      outputfolder=tempdir(),
      draw="PDF"
    )
    meta <- result@metaresult
    write.csv(meta, output, row.names=FALSE)
    """
    with tempfile.TemporaryDirectory(prefix="degora_metavolcanor_") as tmp:
        tmp_path = Path(tmp)
        input_path = tmp_path / "metavolcanor_input.csv"
        result_path = tmp_path / "metavolcanor_result.csv"
        meta_input.rename(columns={"gene_symbol": "Symbol", "lfc": "Log2FC", "p_from_signed_z": "pvalue"}).to_csv(input_path, index=False)
        completed = _run_rscript(["-e", r_code, str(input_path), str(result_path)], timeout=180)
        if completed.returncode != 0:
            detail = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)
            raise RuntimeError(detail or "MetaVolcanoR combining_mv failed")
        result = pd.read_csv(result_path)

    if result.empty:
        raise RuntimeError("MetaVolcanoR combining_mv returned no rows; check runtime dependencies and input compatibility")
    required_result = {"Symbol", "metap", "metafc"}
    missing_result = required_result.difference(result.columns)
    if missing_result:
        raise RuntimeError(f"MetaVolcanoR result missing required columns: {sorted(missing_result)}")

    out = result.loc[:, ["Symbol", "metap", "metafc"]].rename(columns={"Symbol": "gene_symbol", "metap": "pvalue", "metafc": "effect"})
    out = out.merge(stats, on="gene_symbol", how="left")
    finalized = _finalize(
        out,
        method_id="metavolcanor",
        setting_id="default",
        effect="effect",
        pvalue="pvalue",
        n_studies="n_studies",
        runtime_s=time.perf_counter() - start,
        version=f"MetaVolcanoR {_r_package_version('MetaVolcanoR')}",
    )
    finalized["missingness"] = 1.0 - (finalized["n_studies"].astype(float) / float(total_studies))
    return validate_baseline_result(finalized)


def awfisher_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    """Run AWFisher on the source-unit p-value matrix via Rscript."""

    start = time.perf_counter()
    inputs = _study_matrix_inputs(harmonized, min_studies)
    p_matrix = inputs["p_matrix"]
    lfc_matrix = inputs["lfc_matrix"]
    stats = inputs["stats"]
    total_studies = int(inputs["total_studies"])
    if p_matrix.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)

    r_code = r"""
    args <- commandArgs(trailingOnly=TRUE)
    p_input <- args[1]
    lfc_input <- args[2]
    output <- args[3]
    pdat <- read.csv(p_input, check.names=FALSE, stringsAsFactors=FALSE)
    ldat <- read.csv(lfc_input, check.names=FALSE, stringsAsFactors=FALSE)
    genes <- pdat[[1]]
    pmat <- as.matrix(pdat[, -1, drop=FALSE])
    lmat <- as.matrix(ldat[, -1, drop=FALSE])
    storage.mode(pmat) <- "double"
    storage.mode(lmat) <- "double"
    pmat[!is.finite(pmat) | pmat <= 0] <- .Machine$double.xmin
    pmat[pmat > 1] <- 1
    result <- AWFisher::AWFisher_pvalue(pmat)
    weights <- result$weights
    selected_effect <- numeric(nrow(lmat))
    selected_n <- integer(nrow(lmat))
    for (i in seq_len(nrow(lmat))) {
      idx <- weights[i, ] == 1 & is.finite(lmat[i, ])
      selected_n[i] <- sum(idx)
      if (selected_n[i] > 0) {
        selected_effect[i] <- mean(lmat[i, idx])
      } else {
        selected_effect[i] <- mean(lmat[i, is.finite(lmat[i, ])])
      }
    }
    out <- data.frame(
      gene_symbol=genes,
      pvalue=result$pvalues,
      aw_selected_studies=selected_n,
      effect=selected_effect,
      stringsAsFactors=FALSE
    )
    write.csv(out, output, row.names=FALSE)
    """
    with tempfile.TemporaryDirectory(prefix="degora_awfisher_") as tmp:
        tmp_path = Path(tmp)
        p_path = tmp_path / "p_matrix.csv"
        lfc_path = tmp_path / "lfc_matrix.csv"
        result_path = tmp_path / "awfisher_result.csv"
        _write_matrix_with_gene_column(p_matrix, p_path)
        _write_matrix_with_gene_column(lfc_matrix, lfc_path)
        completed = _run_rscript(["-e", r_code, str(p_path), str(lfc_path), str(result_path)], timeout=180)
        if completed.returncode != 0:
            detail = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)
            raise RuntimeError(detail or "AWFisher AWFisher_pvalue failed")
        result = pd.read_csv(result_path)

    out = result.merge(stats.loc[:, ["gene_symbol", "n_studies", "effect"]], on="gene_symbol", how="left", suffixes=("_aw", "_mean"))
    out["effect"] = out["effect_aw"].fillna(out["effect_mean"]).fillna(0.0)
    finalized = _finalize(
        out,
        method_id="awfisher",
        setting_id="default",
        effect="effect",
        pvalue="pvalue",
        n_studies="n_studies",
        runtime_s=time.perf_counter() - start,
        version=f"AWFisher {_r_package_version('AWFisher')}",
    )
    finalized["missingness"] = 1.0 - (finalized["n_studies"].astype(float) / float(total_studies))
    return validate_baseline_result(finalized)


def _metarnaseq_adapter(harmonized: pd.DataFrame, *, method_id: str, min_studies: int = 2) -> pd.DataFrame:
    start = time.perf_counter()
    missing_fill = 0.5 if method_id == "metarnaseq_invnorm" else 1.0
    inputs = _study_matrix_inputs(harmonized, min_studies, missing_pvalue_fill=missing_fill)
    p_matrix = inputs["p_matrix"]
    stats = inputs["stats"]
    total_studies = int(inputs["total_studies"])
    replicate_counts = inputs["replicate_counts"]
    if p_matrix.empty:
        return pd.DataFrame(columns=BASELINE_RESULT_COLUMNS)

    if method_id == "metarnaseq_fisher":
        r_function = "fishercomb"
        extra_args = ""
    elif method_id == "metarnaseq_invnorm":
        r_function = "invnorm"
        extra_args = ", nrep=nrep"
    else:
        raise ValueError(f"unsupported metaRNASeq adapter method: {method_id}")
    r_code = rf"""
    args <- commandArgs(trailingOnly=TRUE)
    input <- args[1]
    output <- args[2]
    nrep <- as.integer(strsplit(args[3], ",")[[1]])
    pdat <- read.csv(input, check.names=FALSE, stringsAsFactors=FALSE)
    genes <- pdat[[1]]
    mat <- as.matrix(pdat[, -1, drop=FALSE])
    storage.mode(mat) <- "double"
    mat[!is.finite(mat) | mat <= 0] <- .Machine$double.xmin
    mat[mat >= 1] <- 1 - .Machine$double.eps
    mat[mat > 1] <- 1 - .Machine$double.eps
    indpval <- as.list(as.data.frame(mat, check.names=FALSE))
    names(indpval) <- colnames(mat)
    for (i in seq_along(indpval)) names(indpval[[i]]) <- genes
    result <- metaRNASeq::{r_function}(indpval{extra_args})
    out <- data.frame(
      gene_symbol=names(result$rawpval),
      pvalue=as.numeric(result$rawpval),
      test_statistic=as.numeric(result$TestStatistic),
      stringsAsFactors=FALSE
    )
    write.csv(out, output, row.names=FALSE)
    """
    with tempfile.TemporaryDirectory(prefix=f"degora_{method_id}_") as tmp:
        tmp_path = Path(tmp)
        p_path = tmp_path / "p_matrix.csv"
        result_path = tmp_path / f"{method_id}_result.csv"
        _write_matrix_with_gene_column(p_matrix, p_path)
        completed = _run_rscript(
            ["-e", r_code, str(p_path), str(result_path), ",".join(str(value) for value in replicate_counts)],
            timeout=180,
        )
        if completed.returncode != 0:
            detail = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)
            raise RuntimeError(detail or f"metaRNASeq {r_function} failed")
        result = pd.read_csv(result_path)

    if method_id == "metarnaseq_invnorm" and _metarnaseq_invnorm_uninformative_result(result):
        raise RuntimeError(
            f"{INVNORM_UNINFORMATIVE_PREFIX} metaRNASeq::invnorm returned an uninformative "
            "or saturated p-value rank after sparse public-summary missing entries were set to "
            "neutral p=0.5; reporting this as a blocked comparator is more faithful than an "
            "all-tie or floor-saturated rank."
        )

    out = result.merge(stats.loc[:, ["gene_symbol", "n_studies", "effect"]], on="gene_symbol", how="left")
    out["pvalue"] = pd.to_numeric(out["pvalue"], errors="coerce").clip(lower=R_PVALUE_FLOOR, upper=1.0)
    test_statistic = pd.to_numeric(out["test_statistic"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    out["test_statistic"] = test_statistic.fillna(pd.to_numeric(out["effect"], errors="coerce").abs()).fillna(0.0)
    finalized = _finalize(
        out,
        method_id=method_id,
        setting_id="default",
        effect="effect",
        pvalue="pvalue",
        n_studies="n_studies",
        runtime_s=time.perf_counter() - start,
        version=f"metaRNASeq {_r_package_version('metaRNASeq')}",
        tie_breaker="test_statistic",
        tie_breaker_desc=True,
    )
    finalized["missingness"] = 1.0 - (finalized["n_studies"].astype(float) / float(total_studies))
    return validate_baseline_result(finalized)


def metarnaseq_fisher_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    return _metarnaseq_adapter(harmonized, method_id="metarnaseq_fisher", min_studies=min_studies)


def metarnaseq_invnorm_adapter(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    return _metarnaseq_adapter(harmonized, method_id="metarnaseq_invnorm", min_studies=min_studies)


TIER0_ADAPTERS: dict[str, Callable[[pd.DataFrame, int], pd.DataFrame]] = {
    "degora_slice": degora_slice_adapter,
    "weighted_stouffer": weighted_stouffer_adapter,
    "unweighted_stouffer": unweighted_stouffer_adapter,
    "fisher": fisher_adapter,
    "rank_product_approx": rank_product_adapter,
    "sign_vote": sign_vote_adapter,
}

DIRECT_ADAPTERS: dict[str, Callable[[pd.DataFrame, int], pd.DataFrame]] = {
    "metavolcanor": metavolcanor_adapter,
    "robustrankaggreg": robustrankaggreg_adapter,
    "awfisher": awfisher_adapter,
    "metarnaseq_fisher": metarnaseq_fisher_adapter,
    "metarnaseq_invnorm": metarnaseq_invnorm_adapter,
}


def run_tier0_baselines(harmonized: pd.DataFrame, min_studies: int = 2) -> dict[str, pd.DataFrame]:
    return {method_id: adapter(harmonized, min_studies) for method_id, adapter in TIER0_ADAPTERS.items()}


def run_available_direct_baselines(
    harmonized: pd.DataFrame,
    preflight: dict[str, Any],
    min_studies: int = 2,
    runtime_failures: list[dict[str, str]] | None = None,
) -> dict[str, pd.DataFrame]:
    """Run direct-prior-art adapters whose dependencies are available."""

    outputs: dict[str, pd.DataFrame] = {}
    for method_id, adapter in DIRECT_ADAPTERS.items():
        spec = next(spec for spec in R_BACKED_METHODS if spec.method_id == method_id)
        package = spec.r_package or ""
        package_status = preflight.get("packages", {}).get(package, {})
        if package and not package_status.get("available"):
            continue
        try:
            outputs[method_id] = adapter(harmonized, min_studies)
        except Exception as exc:  # pragma: no cover - depends on local R runtime state
            if runtime_failures is None:
                raise
            runtime_failures.append({"method_id": method_id, "message": str(exc)})
    return outputs


def _r_version() -> tuple[str | None, str | None]:
    rscript = shutil.which("Rscript")
    if not rscript:
        return None, "Rscript executable not found on PATH"
    completed = _run_rscript(["--version"], timeout=10)
    text = (completed.stdout or completed.stderr or "").strip()
    return (text if completed.returncode == 0 else None), (None if completed.returncode == 0 else text)


def _r_package_available(package: str) -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if not rscript:
        return False, "Rscript executable not found on PATH"
    command = [rscript, "-e", f"quit(status = ifelse(requireNamespace('{package}', quietly=TRUE), 0, 1))"]
    completed = _run_rscript(command[1:], timeout=15)
    if completed.returncode == 0:
        return True, "available"
    detail = (completed.stderr or completed.stdout or "not available").strip().splitlines()
    return False, detail[-1] if detail else "not available"


def _r_package_details(package: str) -> dict[str, Any]:
    """Return version/path/source and the exact R probe log for a package."""

    source = PREFLIGHT_PACKAGE_SOURCES.get(package, "R package")
    rscript = shutil.which("Rscript")
    if not rscript:
        return {
            "available": False,
            "source": source,
            "version": None,
            "install_path": None,
            "failure_log": "Rscript executable not found on PATH",
            "message": "Rscript executable not found on PATH",
        }
    r_code = (
        "pkg <- commandArgs(trailingOnly=TRUE)[1]; "
        "if (!requireNamespace(pkg, quietly=TRUE)) { "
        "cat('requireNamespace returned FALSE\\n'); quit(status=1) }; "
        "desc <- packageDescription(pkg); "
        "cat(paste0('version=', as.character(desc$Version), '\\n')); "
        "cat(paste0('path=', find.package(pkg), '\\n'))"
    )
    completed = _run_rscript(["-e", r_code, package], timeout=15)
    stdout = (completed.stdout or "").strip()
    stderr = (completed.stderr or "").strip()
    combined = "\n".join(part for part in [stdout, stderr] if part)
    if completed.returncode != 0:
        lines = combined.splitlines()
        return {
            "available": False,
            "source": source,
            "version": None,
            "install_path": None,
            "failure_log": combined or "not available",
            "message": lines[-1] if lines else "not available",
        }
    details: dict[str, str] = {}
    for line in stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            details[key] = value
    return {
        "available": True,
        "source": source,
        "version": details.get("version"),
        "install_path": details.get("path"),
        "failure_log": "",
        "message": "available",
    }


def r_preflight_report(packages: Iterable[str] = PREFLIGHT_PACKAGES) -> dict[str, Any]:
    version, r_error = _r_version()
    rpy2_available = importlib.util.find_spec("rpy2") is not None
    package_status: dict[str, dict[str, Any]] = {}
    for package in packages:
        package_status[package] = _r_package_details(package)
    return {
        "r_version": version,
        "r_error": r_error,
        "rscript_path": shutil.which("Rscript"),
        "rpy2_available": rpy2_available,
        "interop_status": "ok" if version and rpy2_available else "blocked",
        "packages": package_status,
        "external_sources": {"hstouffer": HSTOUFFER_SOURCE_PIN},
    }


def hstouffer_materializer_feasibility(harmonized: pd.DataFrame) -> dict[str, Any]:
    """Assess whether the pinned hStouffer script can be run faithfully.

    hStouffer consumes original DESeq2 result tables and its final REM filter
    depends on ``lfcSE``. DEGORA's harmonized table carries enough fields for
    classical Stouffer p-value combination, but not enough to reconstruct the
    published hStouffer inputs without imputation or silent study filtering.
    """

    required_original = HSTOUFFER_REQUIRED_ORIGINAL_COLUMNS
    harmonized_columns = set(harmonized.columns)
    pipeline_counts = harmonized.get("pipeline", pd.Series(dtype="string")).fillna("unknown").astype(str).value_counts().to_dict()
    active_studies = int(harmonized["study_id"].nunique()) if "study_id" in harmonized.columns else 0
    missing_original_fields = [field for field in ["baseMean", "lfcSE", "stat"] if field not in harmonized_columns]
    has_surrogate_stat = "signed_z" in harmonized_columns
    non_deseq2_pipelines = sorted(pipeline for pipeline in pipeline_counts if pipeline.lower() != "deseq2")
    source_input_audit = _hstouffer_source_input_audit(harmonized)
    compatible_source_studies = sorted(
        item["study_id"] for item in source_input_audit if item["audit_status"] == "compatible_header_detected"
    )
    blocked_source_studies = sorted(item["study_id"] for item in source_input_audit if item["audit_status"] != "compatible_header_detected")
    blockers: list[str] = []
    if missing_original_fields:
        blockers.append("harmonized table lacks original DESeq2 fields required by hStouffer REM/filtering: " + ", ".join(missing_original_fields))
    if has_surrogate_stat:
        blockers.append("signed_z is a harmonized p-derived surrogate and is not a pinned original DESeq2 stat column")
    if non_deseq2_pipelines:
        blockers.append("included corpus contains non-DESeq2/unspecified pipeline rows: " + _summarize_blocked_values(non_deseq2_pipelines, 8))
    if blocked_source_studies:
        blockers.append(
            "source-level audit found studies without a compatible original DESeq2-like header or pipeline: "
            + _summarize_blocked_values(blocked_source_studies, 12)
        )
    return {
        "status": "blocked",
        "source_pin": HSTOUFFER_SOURCE_PIN,
        "source_level_method_evidence": HSTOUFFER_SOURCE_LEVEL_EVIDENCE,
        "active_studies": active_studies,
        "pipeline_counts": pipeline_counts,
        "required_original_columns": required_original,
        "harmonized_columns_present": sorted(harmonized_columns),
        "missing_original_fields": missing_original_fields,
        "compatible_source_studies": compatible_source_studies,
        "blocked_source_studies": blocked_source_studies,
        "source_input_audit": source_input_audit,
        "can_materialize_without_imputation_or_filtering": False,
        "blockers": blockers,
        "decision": "Do not emit a runnable hStouffer baseline from DEGORA harmonized/as-published DEG tables until original per-study DESeq2 result files with lfcSE/baseMean/stat are available for every included study or the comparator is explicitly reframed as an approximation.",
    }


def _summarize_blocked_values(values: list[str], limit: int) -> str:
    if len(values) <= limit:
        return ", ".join(values)
    shown = ", ".join(values[:limit])
    remaining = len(values) - limit
    return f"{shown}, ... ({remaining} more; total {len(values)})"


def failure_ledger(
    corpus: str,
    preflight: dict[str, Any],
    hstouffer_feasibility: dict[str, Any] | None = None,
    harmonized: pd.DataFrame | None = None,
    runtime_failures: list[dict[str, str]] | None = None,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    runtime_failure_by_method = {item["method_id"]: item["message"] for item in (runtime_failures or [])}
    for spec in R_BACKED_METHODS:
        blocker_id = ""
        message = ""
        if not preflight.get("r_version") and spec.method_id != "hstouffer":
            blocker_id = "rscript_missing"
            message = str(preflight.get("r_error") or "Rscript unavailable")
        elif spec.method_id == "hstouffer":
            blocker_id = "hstouffer_deg_table_materializer_blocked" if hstouffer_feasibility else "hstouffer_input_materializer_unverified"
            if hstouffer_feasibility:
                source_pin = hstouffer_feasibility.get("source_pin", HSTOUFFER_SOURCE_PIN)
                blockers = "; ".join(hstouffer_feasibility.get("blockers", []))
                message = (
                    f"hStouffer source pinned to CR4CKID/hStouffer@{source_pin.get('commit')}; "
                    "faithful materialization from the harmonized/as-published DEG table is blocked. "
                    f"{blockers}"
                )
            else:
                message = (
                    f"hStouffer source pinned to CR4CKID/hStouffer@{HSTOUFFER_SOURCE_PIN['commit']}, "
                    "but no harmonized input feasibility report was provided for this preflight-only run."
                )
        elif spec.method_id == "awmeta":
            if harmonized is not None:
                feasibility = awmeta_deg_table_feasibility(harmonized)
                blocker_id = str(feasibility["blocker_id"])
                message = str(feasibility["message"])
            else:
                blocker_id = "awmeta_deg_table_adapter_missing"
                message = "AWmeta source is available, but a faithful DEG-table adapter for AW-Fisher/AW-REM is not present in this thin-slice codebase."
        elif spec.method_id in runtime_failure_by_method:
            message = runtime_failure_by_method[spec.method_id]
            if message.startswith(INVNORM_UNINFORMATIVE_PREFIX):
                blocker_id = INVNORM_UNINFORMATIVE_BLOCKER_ID
                message = message.removeprefix(INVNORM_UNINFORMATIVE_PREFIX).strip()
            else:
                blocker_id = f"{spec.method_id}_runtime_failed"
        elif spec.method_id in SUMMARY_ONLY_BLOCKERS:
            blocker = SUMMARY_ONLY_BLOCKERS[spec.method_id]
            pkg = preflight.get("packages", {}).get(spec.r_package or "", {})
            package_note = ""
            if spec.r_package:
                package_note = f" Package status for {spec.r_package}: {pkg.get('message', 'unresolved')}."
            blocker_id = str(blocker["blocker_id"])
            message = str(blocker["message"]) + package_note
        elif spec.r_package:
            pkg = preflight.get("packages", {}).get(spec.r_package, {})
            if not pkg.get("available"):
                blocker_id = f"r_package_missing:{spec.r_package}"
                message = str(pkg.get("failure_log") or pkg.get("message") or "R package unavailable")
            elif spec.method_id not in DIRECT_ADAPTERS:
                blocker_id = f"{spec.method_id}_wrapper_missing"
                message = f"{spec.method_id} package is available, but no uniform-schema wrapper has been implemented and validated."
        if blocker_id:
            records.append(
                {
                    "corpus": corpus,
                    "method_id": spec.method_id,
                    "setting_id": spec.setting_id,
                    "tier": spec.tier,
                    "status": "open_s1_blocker",
                    "blocker_id": blocker_id,
                    "message": message,
                    "resolution": spec.resolution,
                }
            )
    return pd.DataFrame.from_records(records, columns=FAILURE_LEDGER_COLUMNS)


def _claim_allowed(*, run_status: str, failure_mode: str = "") -> str:
    """Return whether this row can support a same-input benchmark claim."""

    if run_status == "ok" and not failure_mode:
        return "yes"
    return "no"


def baseline_parity_matrix(
    *,
    corpus: str,
    outputs: dict[str, pd.DataFrame],
    preflight: dict[str, Any],
    ledger: pd.DataFrame,
    min_studies: int,
) -> pd.DataFrame:
    """Build the PRD parity gate matrix for resolved and blocked baselines."""

    records: list[dict[str, Any]] = []
    emitted_settings: set[tuple[str, str]] = set()

    for method_id, result in outputs.items():
        meta = _method_metadata(method_id)
        setting_id = str(meta["setting_id"])
        emitted_settings.add((method_id, setting_id))
        version_or_commit = str(meta["version_or_commit"]).format(version=result["version"].iloc[0] if not result.empty else "unknown")
        if method_id == "degora_slice":
            equal_tuning_status = "locked_reference_no_extra_tuning"
        elif method_id in DIRECT_METHOD_METADATA:
            equal_tuning_status = "default_direct_prior_art_wrapper_no_extra_tuning"
        else:
            equal_tuning_status = "default_only_classical_baseline"
        records.append(
            {
                "method_id": method_id,
                "method_family": meta["method_family"],
                "implementation_source": meta["implementation_source"],
                "version_or_commit": version_or_commit,
                "setting_id": setting_id,
                "setting_type": "locked" if method_id == "degora_slice" else "default",
                "parameters": str(meta["parameters"]).format(min_studies=min_studies),
                "input_requirements": meta["input_requirements"],
                "output_schema_status": "uniform_schema_validated",
                "run_status": "ok",
                "failure_mode": "",
                "equal_tuning_status": equal_tuning_status,
                "claim_allowed": _claim_allowed(run_status="ok"),
                "notes": f"{len(result)} rows emitted for {corpus}.",
            }
        )

    ledger_by_method = {str(row.method_id): row for row in ledger.itertuples(index=False)}
    for spec in R_BACKED_METHODS:
        pkg = preflight.get("packages", {}).get(spec.r_package or "", {})
        row = ledger_by_method.get(spec.method_id)
        blocked = row is not None
        failure_mode = str(row.blocker_id) if blocked else ""
        notes = str(row.message) if blocked else "Package available; wrapper execution still requires implementation validation."
        if spec.method_id == "hstouffer":
            hstouffer_pin = preflight.get("external_sources", {}).get("hstouffer", HSTOUFFER_SOURCE_PIN)
            version_or_commit = f"{hstouffer_pin.get('repo_url')}@{hstouffer_pin.get('commit')}"
        else:
            version_or_commit = str(pkg.get("version") or "unresolved")
        if pkg.get("install_path"):
            version_or_commit = f"{version_or_commit}; {pkg['install_path']}"
        for setting_id, setting_type, parameters in [
            (spec.setting_id, spec.setting_type, spec.parameters),
            ("tuned_min_studies_3", "one_tuned_setting", spec.parameters.replace(f"min_studies={min_studies}", "min_studies=3") if f"min_studies={min_studies}" in spec.parameters else f"{spec.parameters}; min_studies=3"),
        ]:
            if (spec.method_id, setting_id) in emitted_settings:
                continue
            records.append(
                {
                    "method_id": spec.method_id,
                    "method_family": spec.method_family,
                    "implementation_source": spec.implementation_source,
                    "version_or_commit": version_or_commit,
                    "setting_id": setting_id,
                    "setting_type": setting_type,
                    "parameters": parameters,
                    "input_requirements": spec.input_requirements,
                    "output_schema_status": "not_emitted_blocked" if blocked else "pending_wrapper_execution",
                    "run_status": "blocked" if blocked else "not_run",
                    "failure_mode": failure_mode,
                    "equal_tuning_status": "blocked_before_tuning_no_result_dependent_changes" if blocked else "eligible_for_default_plus_one_tuned_setting",
                    "claim_allowed": _claim_allowed(run_status="blocked" if blocked else "not_run", failure_mode=failure_mode),
                    "notes": notes,
                }
            )
    return pd.DataFrame.from_records(records, columns=PARITY_MATRIX_COLUMNS)


def _ledger_status(method_id: str, outputs: dict[str, pd.DataFrame], ledger: pd.DataFrame) -> str:
    if ledger.empty:
        return "run_ok" if method_id in outputs else "not_run"
    rows = ledger.loc[ledger["method_id"].astype(str).eq(method_id)]
    if rows.empty:
        return "run_ok" if method_id in outputs else "not_run"
    row = rows.iloc[0]
    return f"blocked:{row['blocker_id']}"


def public_summary_tool_input_requirements(
    *,
    corpus: str,
    harmonized: pd.DataFrame,
    outputs: dict[str, pd.DataFrame],
    preflight: dict[str, Any],
    ledger: pd.DataFrame,
    hstouffer_report: dict[str, Any],
    awmeta_report: dict[str, Any],
) -> pd.DataFrame:
    """Build a reviewer-facing table for public-file comparator feasibility."""

    available_columns = sorted(str(column) for column in harmonized.columns)
    compact_available = ", ".join(column for column in ["study_id", "source_unit_id", "gene_symbol", "lfc", "pvalue", "padj", "signed_z", "normalized_rank", "n_ctrl", "n_treat", "pipeline", "source_path"] if column in harmonized.columns)
    packages = preflight.get("packages", {})

    def package_source(package: str | None) -> str:
        if not package:
            return "DEGORA repository implementation"
        status = packages.get(package, {})
        version = status.get("version") or "not installed"
        return f"{PREFLIGHT_PACKAGE_SOURCES.get(package, package)}; version={version}"

    def record(
        method_id: str,
        method_family: str,
        public_status: str,
        can_run: str,
        required: str,
        missing: str,
        decision: str,
        manuscript_use: str,
        source: str,
        notes: str,
        current_status: str | None = None,
    ) -> dict[str, Any]:
        return {
            "method_id": method_id,
            "method_family": method_family,
            "public_summary_deg_status": public_status,
            "can_run_from_current_public_files": can_run,
            "current_pipeline_status": current_status or _ledger_status(method_id, outputs, ledger),
            "required_inputs": required,
            "current_public_inputs_available": compact_available,
            "missing_or_nonfaithful_inputs": missing,
            "faithful_adapter_decision": decision,
            "manuscript_use": manuscript_use,
            "source_or_package": source,
            "notes": notes,
        }

    invnorm_status = _ledger_status("metarnaseq_invnorm", outputs, ledger)
    invnorm_blocked = invnorm_status.startswith("blocked:")
    invnorm_missing = (
        "official inverse-normal adapter produced an uninformative or saturated rank on this public-summary corpus"
        if invnorm_blocked
        else "true per-study one-sided p-values are absent; source-unit sample-size weights are used where available"
    )
    invnorm_decision = (
        "blocked for this corpus after the official adapter returned an uninformative sparse public-summary result"
        if invnorm_blocked
        else "official metaRNASeq invnorm adapter emitted with documented weight proxy"
    )
    invnorm_use = "blocked comparator row" if invnorm_blocked else "direct prior-art comparator with caveat"
    invnorm_note = (
        "DEGORA marks this comparator blocked instead of reporting all-tie or floor-saturated ranks."
        if invnorm_blocked
        else "Missing study-gene entries are set to neutral p=0.5 to avoid treating absent public-summary evidence as infinite inverse-normal evidence; if the official function returns an all-p=1/no-statistic result, DEGORA marks the comparator blocked instead of emitting alphabetical ranks."
    )

    records = [
        record(
            "degora_slice",
            "DEGORA reference",
            "native_public_summary_method",
            "yes",
            "study/source-unit, gene symbol, p-value or signed z, log2FC, optional sample sizes, optional quality metadata",
            "",
            "native reference output; no external comparator adaptation needed",
            "primary method",
            package_source(None),
            f"{corpus}: public harmonized columns inspected={len(available_columns)}.",
            current_status="run_ok" if "degora_slice" in outputs else "not_run",
        ),
        record(
            "weighted_stouffer",
            "classical p-value/z method",
            "summary_compatible",
            "yes",
            "per-study signed z and sample-size weights",
            "",
            "runnable as transparent Tier 0 baseline",
            "classical comparator",
            "SciPy norm.sf",
            "Uses DEGORA's conservative source-unit collapsed signed z values.",
        ),
        record(
            "fisher",
            "classical p-value method",
            "summary_compatible",
            "yes",
            "per-study p-values",
            "",
            "runnable as transparent Tier 0 baseline",
            "classical comparator",
            "SciPy combine_pvalues",
            "Uses two-sided p-values recovered from signed z to avoid one-sided inflation.",
        ),
        record(
            "robustrankaggreg",
            "rank aggregation",
            "summary_compatible_official",
            "yes",
            "within-study ranked gene lists",
            "",
            "official R adapter emitted when package is installed",
            "direct prior-art comparator",
            package_source("RobustRankAggreg"),
            "Does not need raw expression; ranks are reconstructed from public DEG tables.",
        ),
        record(
            "metavolcanor",
            "p-value/effect meta-analysis",
            "summary_compatible_official",
            "yes",
            "per-study gene, logFC, p-value tables",
            "",
            "official R adapter emitted when package is installed",
            "direct prior-art comparator",
            package_source("MetaVolcanoR"),
            "Uses MetaVolcanoR p-combination mode only; effect-size variance modes are not claimed.",
        ),
        record(
            "awfisher",
            "adaptive weighted p-value method",
            "summary_compatible_official",
            "yes",
            "gene-by-study p-value matrix",
            "",
            "official AWFisher adapter emitted when package is installed",
            "direct prior-art comparator",
            package_source("AWFisher"),
            "Missing study-gene entries are set to p=1; AWFisher returns adaptive study-contribution weights.",
        ),
        record(
            "metarnaseq_fisher",
            "RNA-seq p-value combination",
            "summary_compatible_official",
            "yes",
            "per-study p-value vectors",
            "",
            "official metaRNASeq fishercomb adapter emitted when package is installed",
            "direct prior-art comparator",
            package_source("metaRNASeq"),
            "Uses two-sided p_from_signed_z because public summary tables do not provide uniform one-sided p-values.",
        ),
        record(
            "metarnaseq_invnorm",
            "RNA-seq inverse-normal p-value combination",
            "summary_compatible_official",
            "blocked_on_current_public_files" if invnorm_blocked else "yes_with_weight_proxy",
            "per-study p-value vectors and replicate-count weights",
            invnorm_missing,
            invnorm_decision,
            invnorm_use,
            package_source("metaRNASeq"),
            invnorm_note,
            current_status=invnorm_status,
        ),
        record(
            "rank_product_approx",
            "rank aggregation approximation",
            "summary_compatible_approximation",
            "yes_as_approximation",
            "within-study normalized ranks",
            "official RankProd replicate/origin inputs are not present",
            "keep as explicit approximation, not official RankProd",
            "sensitivity comparator",
            "DEGORA repository implementation",
            "Useful because public DEG files often retain ranks but not replicate-level expression.",
        ),
        record(
            "rankprod_exact",
            "official rank product",
            "blocked_for_public_summary_only",
            "no",
            "replicate-level expression or rank-product inputs with class/origin labels",
            "raw expression/class labels/origin metadata; current R package also unavailable if RankProd preflight is missing",
            "do not run from public DEG-only files",
            "blocked comparator row",
            package_source("RankProd"),
            "Prevents overclaiming that summary-rank approximation is the official package.",
        ),
        record(
            "metade",
            "microarray meta-analysis suite",
            "partial_summary_coverage",
            "partial",
            "p-value matrices for p-value modes; effect sizes plus variances or raw expression for FEM/REM/effect-size modes",
            "per-study effect variances and raw expression objects for non-p-value modes",
            "p-value/rank concepts are covered by Fisher/Stouffer/AWFisher/rank baselines; effect-size modes blocked",
            "blocked/covered comparator row",
            package_source("MetaDE"),
            "Avoids double-counting MetaDE methods already represented by primitive p-value/rank comparators.",
        ),
        record(
            "dexma",
            "gene-expression meta-analysis workflow",
            "blocked_for_public_summary_only",
            "no",
            "expression matrices or GEO-derived expression objects plus phenotype/sample metadata",
            "sample-level expression matrices, phenotype metadata, normalization/QC objects",
            "do not run from already-filtered supplementary DEG tables",
            "blocked comparator row",
            package_source("DExMA"),
            "DExMA is relevant for a future raw-expression benchmark, not the current public summary-file benchmark.",
        ),
        record(
            "metaintegrator",
            "signature meta-analysis",
            "blocked_for_public_summary_only",
            "no",
            "study expression data, phenotype labels, effect sizes and variances",
            "cohort-level expression and phenotype objects; comparable variance inputs",
            "do not run from public DEG-only files",
            "blocked comparator row",
            package_source("MetaIntegrator"),
            "Useful reviewer-facing no-go row because it is a common expression-signature meta-analysis package.",
        ),
        record(
            "hstouffer",
            "RNA-seq Stouffer/REM workflow",
            "blocked_for_public_summary_only",
            "no",
            "DESeq2-like full result tables with id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj",
            "; ".join(hstouffer_report.get("blockers", [])) or "original DESeq2 lfcSE/baseMean/stat fields",
            "blocked until original compatible DESeq2-like files exist for every active source",
            "blocked comparator row",
            f"{HSTOUFFER_SOURCE_PIN['repo_url']}@{HSTOUFFER_SOURCE_PIN['commit']}",
            "The feasibility JSON gives per-source header evidence.",
        ),
        record(
            "awmeta",
            "adaptive weighted effect-size meta-analysis",
            "blocked_for_public_summary_only",
            "no",
            "per-study effect sizes plus original within-study variance/SE or documented equivalent weights",
            str(awmeta_report.get("message", "variance/SE inputs unavailable")),
            "blocked; never infer SE from p-values or signed z",
            "blocked comparator row",
            package_source("metafor"),
            "The feasibility JSON gives per-source variance/SE audit evidence.",
        ),
        record(
            "gene_set_pathway_tools",
            "pathway enrichment / gene-set methods",
            "orthogonal_not_gene_rank_comparator",
            "not_as_gene_level_comparator",
            "ranked full gene universe or expression matrix plus gene sets",
            "these tools answer pathway-level enrichment, not gene-level cross-study ranking",
            "report as orthogonal downstream analysis, not a replacement comparator for DEGORA gene ranking",
            "scope guardrail",
            "fgsea/GSEA/GSVA class methods",
            "Include only if the manuscript adds a pathway endpoint.",
            current_status="out_of_scope",
        ),
    ]
    for prior in PUBLIC_WORKFLOW_PRIOR_ART_ROWS:
        records.append(
            record(
                str(prior["method_id"]),
                str(prior["method_family"]),
                str(prior["public_status"]),
                str(prior["can_run"]),
                str(prior["required"]),
                str(prior["missing"]),
                str(prior["decision"]),
                str(prior["manuscript_use"]),
                str(prior["source"]),
                str(prior["notes"]),
                current_status="workflow_or_resource_prior_art",
            )
        )
    return pd.DataFrame.from_records(records, columns=PUBLIC_SUMMARY_TOOL_INPUT_COLUMNS)


def _write_markdown_table(frame: pd.DataFrame, path: Path) -> None:
    lines = [
        "| " + " | ".join(frame.columns) + " |",
        "| " + " | ".join(["---"] * len(frame.columns)) + " |",
    ]
    for row in frame.itertuples(index=False):
        values = [str(value).replace("\n", " ").replace("|", "\\|") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    path.write_text("\n".join(lines) + "\n")


def baseline_manifest(
    *,
    outputs: dict[str, pd.DataFrame],
    output_paths: dict[str, Path],
    command: str,
    harmonized_path: Path,
    preflight_path: Path,
    ledger_path: Path,
    parity_path: Path,
    direct_blockers: int,
    extra_artifacts: list[tuple[Path, str, int, str]] | None = None,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for method_id, result in outputs.items():
        meta = _method_metadata(method_id)
        path = output_paths[method_id]
        records.append(
            {
                "artifact": str(path),
                "method_id": method_id,
                "setting_id": str(meta["setting_id"]),
                "artifact_type": "baseline_result",
                "status": "ok",
                "source_command": command,
                "input_harmonized": str(harmonized_path.resolve()),
                "rows": int(len(result)),
                "notes": "uniform schema result with sibling .source",
            }
        )
    support_artifacts = [
        (preflight_path, "r_preflight_report", 1, "R/R-package probe with version/source/path/failure logs"),
        (ledger_path, "failure_ledger", direct_blockers, "direct-prior-art blockers; empty means no open blockers"),
        (parity_path, "baseline_parity_matrix", -1, "binary manuscript-claim gate for default plus one tuned setting where runnable"),
    ]
    support_artifacts.extend(extra_artifacts or [])
    for path, artifact_type, rows, notes in support_artifacts:
        records.append(
            {
                "artifact": str(path),
                "method_id": "all",
                "setting_id": "all",
                "artifact_type": artifact_type,
                "status": "ok",
                "source_command": command,
                "input_harmonized": str(harmonized_path.resolve()),
                "rows": rows,
                "notes": notes,
            }
        )
    return pd.DataFrame.from_records(records, columns=BASELINE_MANIFEST_COLUMNS)


def write_baseline_outputs(
    harmonized_path: Path,
    output_dir: Path,
    *,
    corpus: str,
    min_studies: int = 2,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    outputs = run_tier0_baselines(harmonized, min_studies=min_studies)
    command = f"make -C outputs/code baseline CORPUS={corpus} HARMONIZED={harmonized_path.resolve()} BASELINE_OUTDIR={output_dir.resolve()}"
    preflight = r_preflight_report()
    runtime_failures: list[dict[str, str]] = []
    outputs.update(run_available_direct_baselines(harmonized, preflight, min_studies=min_studies, runtime_failures=runtime_failures))

    written: list[str] = []
    output_paths: dict[str, Path] = {}
    for method_id, result in outputs.items():
        result = validate_baseline_result(result)
        path = output_dir / f"{corpus}_{method_id}_default.tsv"
        if method_id == "degora_slice":
            path = output_dir / f"{corpus}_{method_id}_locked.tsv"
        result.to_csv(path, sep="\t", index=False)
        path.with_suffix(path.suffix + ".source").write_text(command + "\n")
        written.append(str(path))
        output_paths[method_id] = path

    preflight_path = output_dir / "r_preflight_report.json"
    preflight_path.write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    preflight_path.with_suffix(preflight_path.suffix + ".source").write_text(command + "\n")

    hstouffer_report = hstouffer_materializer_feasibility(harmonized)
    hstouffer_path = output_dir / "hstouffer_feasibility_report.json"
    hstouffer_path.write_text(json.dumps(hstouffer_report, indent=2, sort_keys=True) + "\n")
    hstouffer_path.with_suffix(hstouffer_path.suffix + ".source").write_text(command + "\n")

    awmeta_report = awmeta_deg_table_feasibility(harmonized)
    awmeta_path = output_dir / "awmeta_feasibility_report.json"
    awmeta_path.write_text(json.dumps(awmeta_report, indent=2, sort_keys=True) + "\n")
    awmeta_path.with_suffix(awmeta_path.suffix + ".source").write_text(command + "\n")

    ledger = failure_ledger(
        corpus,
        preflight,
        hstouffer_feasibility=hstouffer_report,
        harmonized=harmonized,
        runtime_failures=runtime_failures,
    )
    ledger_path = output_dir / "baseline_failure_ledger.csv"
    ledger.to_csv(ledger_path, index=False)
    ledger_path.with_suffix(ledger_path.suffix + ".source").write_text(command + "\n")

    parity = baseline_parity_matrix(corpus=corpus, outputs=outputs, preflight=preflight, ledger=ledger, min_studies=min_studies)
    parity_path = output_dir / "baseline_parity_matrix.csv"
    parity.to_csv(parity_path, index=False)
    parity_path.with_suffix(parity_path.suffix + ".source").write_text(command + "\n")

    input_requirements = public_summary_tool_input_requirements(
        corpus=corpus,
        harmonized=harmonized,
        outputs=outputs,
        preflight=preflight,
        ledger=ledger,
        hstouffer_report=hstouffer_report,
        awmeta_report=awmeta_report,
    )
    input_requirements_path = output_dir / "public_summary_tool_input_requirements.csv"
    input_requirements.to_csv(input_requirements_path, index=False)
    input_requirements_path.with_suffix(input_requirements_path.suffix + ".source").write_text(command + "\n")
    input_requirements_md_path = output_dir / "public_summary_tool_input_requirements.md"
    _write_markdown_table(input_requirements, input_requirements_md_path)
    input_requirements_md_path.with_suffix(input_requirements_md_path.suffix + ".source").write_text(command + "\n")

    direct_blockers = int(len(ledger))
    manifest_path = output_dir / "baseline_manifest.csv"
    manifest = baseline_manifest(
        outputs=outputs,
        output_paths=output_paths,
        command=command,
        harmonized_path=harmonized_path,
        preflight_path=preflight_path,
        ledger_path=ledger_path,
        parity_path=parity_path,
        direct_blockers=direct_blockers,
        extra_artifacts=[
            (
                input_requirements_path,
                "public_summary_tool_input_requirements",
                int(len(input_requirements)),
                "supplementary reviewer table: why each comparator can or cannot run from public supplementary DEG files",
            ),
            (
                input_requirements_md_path,
                "public_summary_tool_input_requirements_markdown",
                int(len(input_requirements)),
                "markdown rendering of the public-file comparator feasibility table",
            ),
        ],
    )
    manifest.loc[len(manifest)] = {
        "artifact": str(hstouffer_path),
        "method_id": "hstouffer",
        "setting_id": "default",
        "artifact_type": "feasibility_report",
        "status": "blocked",
        "source_command": command,
        "input_harmonized": str(harmonized_path.resolve()),
        "rows": 1,
        "notes": "pinned source plus faithful input-materializer blocker dossier",
    }
    manifest.loc[len(manifest)] = {
        "artifact": str(awmeta_path),
        "method_id": "awmeta",
        "setting_id": "default",
        "artifact_type": "feasibility_report",
        "status": "blocked",
        "source_command": command,
        "input_harmonized": str(harmonized_path.resolve()),
        "rows": 1,
        "notes": "faithful AWmeta/AW-REM blocker dossier for missing variance/SE inputs",
    }
    # Now that the manifest path is known, record its own row and rewrite.
    manifest.loc[len(manifest)] = {
        "artifact": str(manifest_path),
        "method_id": "all",
        "setting_id": "all",
        "artifact_type": "baseline_manifest",
        "status": "ok",
        "source_command": command,
        "input_harmonized": str(harmonized_path.resolve()),
        "rows": int(len(manifest) + 1),
        "notes": "manifest of generated baseline artifacts",
    }
    manifest.to_csv(manifest_path, index=False)
    manifest_path.with_suffix(manifest_path.suffix + ".source").write_text(command + "\n")

    return {
        "corpus": corpus,
        "harmonized_path": str(harmonized_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "tier0_outputs": written,
        "n_open_s1_blockers": int(len(ledger)),
        "parity_matrix": str(parity_path),
        "public_summary_tool_input_requirements": str(input_requirements_path),
        "baseline_manifest": str(manifest_path),
    }
