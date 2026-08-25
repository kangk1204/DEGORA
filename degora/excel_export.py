"""Excel audit workbook export for DEGORA run outputs."""

from __future__ import annotations

import json
import os
import sqlite3
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.comments import Comment

from . import runtime_version_info
from .excel_io import locked_panel_mask, read_config_sheet
from .harmonize import canonical_gene_symbol
from .formula_safety import EXCEL_ERROR_LITERALS, restore_formula_text_if_marked
from .provenance import (
    normalize_ooxml_zip,
    portable_command,
    portable_path,
    set_reproducible_workbook_properties,
    write_source_sidecar,
)
from .score_db import HETEROGENEITY_RULE, RANDOM_EFFECTS_STOUFFER_RULE


EXCEL_MAX_ROWS = 1_048_576
DEFAULT_WORKBOOK_NAME = "DEGORA_output.xlsx"
COMMENT_AUTHOR = "DEGORA"
FORMULA_PREFIX_WHITESPACE = " \t\r\n"

# Large corpora (e.g. the cross-platform hypoxia benchmark) can produce hundreds of
# thousands of per-source evidence rows. Writing them all into a single Excel sheet is slow
# and can exhaust memory or fail in constrained environments, so the Gene_evidence sheet is
# capped to the top-ranked genes' evidence. Complete source-level evidence remains in SQLite,
# while the gene-scores CSV remains the complete ranked-gene table. Override with
# DEGORA_EXCEL_EVIDENCE_ROW_CAP (0 disables it).
EVIDENCE_SHEET_ROW_CAP = 100_000


def _evidence_row_cap() -> int:
    raw = os.environ.get("DEGORA_EXCEL_EVIDENCE_ROW_CAP", "")
    if raw == "":
        return EVIDENCE_SHEET_ROW_CAP
    try:
        return max(int(raw), 0)
    except (TypeError, ValueError):
        return EVIDENCE_SHEET_ROW_CAP


def _cap_evidence_for_sheet(evidence: pd.DataFrame, genes: pd.DataFrame, cap: int) -> tuple[pd.DataFrame, bool]:
    """Return (evidence_for_sheet, was_capped). Full source-level evidence stays in SQLite."""
    if cap <= 0 or len(evidence) <= cap:
        return evidence, False
    ordered = evidence
    rank_col = "quality_weighted_degora_rank"
    if "gene_symbol" in evidence.columns and {"gene_symbol", rank_col} <= set(genes.columns):
        ranks = genes.drop_duplicates("gene_symbol").set_index("gene_symbol")[rank_col]
        ordered = evidence.assign(_gene_rank=evidence["gene_symbol"].map(ranks))
        ordered = ordered.sort_values("_gene_rank", kind="mergesort", na_position="last").drop(columns=["_gene_rank"])
    return ordered.head(cap), True


SHEET_GUIDE: dict[str, dict[str, str]] = {
    "Workbook_guide": {
        "row_grain": "one row per workbook sheet",
        "what_it_contains": "Short navigation guide for the DEGORA Excel output.",
        "how_to_use": "Start here to identify the ranked output and source-level audit tabs.",
    },
    "Column_dictionary": {
        "row_grain": "one row per exported column",
        "what_it_contains": "Column meanings, value scales, and blank-value interpretation.",
        "how_to_use": "Use this tab when a header is not self-explanatory.",
    },
    "Run_summary": {
        "row_grain": "one row per run-level metric",
        "what_it_contains": "Result folder, row counts, source-unit counts, curated-gene recall, and top genes.",
        "how_to_use": "Quickly confirm run size and benchmark recovery before opening detailed tables.",
    },
    "Gene_scores": {
        "row_grain": "one row per ranked gene",
        "what_it_contains": "The primary DEGORA ranked gene table with score components and audit metrics.",
        "how_to_use": "Use this as the main prioritized gene list; higher score and lower rank indicate stronger priority.",
    },
    "Gene_evidence": {
        "row_grain": "one row per gene per collapsed source unit",
        "what_it_contains": "Traceable source-unit evidence behind each gene score.",
        "how_to_use": "Audit why a gene ranked highly and which studies, assays, and contrasts contributed evidence.",
    },
    "Source_units": {
        "row_grain": "one row per source/contrast entry",
        "what_it_contains": "Input source metadata, sample counts, platforms, and source links.",
        "how_to_use": "Check which independent source units and contrasts were included.",
    },
    "Curated_lookup": {
        "row_grain": "one row per curated GoldPanel gene",
        "what_it_contains": "Curated benchmark genes matched back to DEGORA ranks and top-k hit flags.",
        "how_to_use": "Inspect benchmark recall gene by gene.",
    },
    "Source_quality": {
        "row_grain": "one row per source-quality diagnostic",
        "what_it_contains": "Diagnostics used to review source reliability when available.",
        "how_to_use": "Review source-level quality warnings without changing the primary rank by hand.",
    },
    "Metadata": {
        "row_grain": "one row per metadata field",
        "what_it_contains": "Flattened DEGORA score metadata JSON.",
        "how_to_use": "Confirm settings, paths, versions, and provenance recorded for the run.",
    },
    "SQLite_meta": {
        "row_grain": "one row per SQLite metadata key",
        "what_it_contains": "Raw key/value metadata stored in the DEGORA SQLite database.",
        "how_to_use": "Use for low-level audit or reproducibility checks.",
    },
}


COLUMN_DEFINITIONS: dict[str, tuple[str, str, str]] = {
    "sheet": ("Workbook sheet name.", "text", "not expected"),
    "row_grain": ("What one row represents in the sheet.", "text", "not expected"),
    "what_it_contains": ("Short description of the sheet contents.", "text", "not expected"),
    "how_to_use": ("Suggested interpretation or first-use action.", "text", "not expected"),
    "column": ("Column header exactly as exported.", "text", "not expected"),
    "meaning": ("Plain-language meaning of the column.", "text", "not expected"),
    "scale_or_values": ("Expected value scale or examples.", "text", "not expected"),
    "missing_or_blank": ("How to read blank or missing values.", "text", "not expected"),
    "field": ("Run or metadata field name.", "text", "not expected"),
    "value": ("Value associated with the field/key.", "text or number", "blank means unavailable or not applicable"),
    "key": ("SQLite metadata key.", "text", "not expected"),
    "degora_rank": ("Rank by the original DEGORA score.", "1 is highest priority", "blank means not ranked"),
    "rank_label": ("Human-readable rank label.", "text such as '#1 / 20000'", "blank means not ranked"),
    "gene_symbol": ("HGNC-style gene symbol after harmonization.", "uppercase gene symbol", "not expected for scored rows"),
    "evidence_tier": ("Qualitative evidence tier assigned from DEGORA support and score metrics.", "text label", "blank means not classified"),
    "degora_score": ("Original DEGORA prioritization score; relative index, not a probability.", "higher is stronger", "blank means not scored"),
    "top_percent": ("Percentile position from the top of the ranked gene list.", "0-100, lower is better", "blank means not ranked"),
    "percentile": ("Rank percentile in the scored gene list.", "0-100", "blank means not ranked"),
    "top_percent_label": ("Human-readable top-percent label.", "text", "blank means not ranked"),
    "consensus_direction": ("Consensus up/down direction from contributing source units.", "up, down, mixed, or similar text", "blank means unavailable"),
    "n_source_units": ("Number of independent collapsed source units supporting the gene.", "non-negative integer", "0/blank means no retained support"),
    "n_contrasts_observed": ("Number of input DEG contrast rows in which the gene was observed.", "non-negative integer", "blank means unavailable"),
    "support_label": ("Human-readable source-unit support summary.", "text", "blank means unavailable"),
    "source_units": ("Semicolon-delimited source-unit IDs contributing to the gene.", "text list", "blank means unavailable"),
    "sign_concordance": ("Fraction of source-unit evidence agreeing with the consensus direction.", "0-1, higher is more concordant", "blank means unavailable"),
    "direction_label": ("Human-readable direction concordance summary.", "text", "blank means unavailable"),
    "support_score": ("Score component from log-scaled independent source-unit support.", "0-1, higher is stronger", "blank means unavailable"),
    "direction_score": ("Score component from agreement with consensus direction.", "0-1, higher is stronger", "blank means unavailable"),
    "evidence_score": ("Score component from source-unit signed-z evidence strength.", "0-1, higher is stronger", "blank means unavailable"),
    "rank_score_component": ("Score component from high within-source ranks after source-unit collapse.", "0-1, higher is stronger", "blank means unavailable"),
    "effect_score": ("Score component from absolute weighted log2 fold-change when available.", "0-1, higher is stronger", "blank means unavailable"),
    "priority_rank": ("Rank by the priority score variant.", "1 is highest priority", "blank means unavailable"),
    "priority_score": ("Alternative priority score retained for audit.", "higher is stronger", "blank means unavailable"),
    "priority_top_percent": ("Top-percent position for the priority score variant.", "0-100, lower is better", "blank means unavailable"),
    "evidence_reliability_score": ("Conditional reliability summary over available diagnostics; it does not determine the primary rank.", "0-100, higher is stronger", "compare together with evidence_reliability_components_used"),
    "evidence_reliability_components_used": ("Number of reliability diagnostics included after availability filtering.", "3 or 4", "3 means LOO was unavailable; 4 means LOO was evaluated"),
    "direction_confidence_index": ("Direction-confidence index from concordance and support.", "0-1, higher is stronger", "blank means unavailable"),
    "quality_weighted_direction_confidence_index": ("Quality-weighted version of direction confidence.", "0-1, higher is stronger", "blank means unavailable"),
    "direction_concordant_source_units": ("Number of source units agreeing with the consensus direction.", "integer", "blank means unavailable"),
    "direction_total_source_units": ("Number of source units used for direction assessment.", "integer", "blank means unavailable"),
    "direction_posterior_mean": (
        "Legacy field name for the Beta(1,1)-shrunk direction-agreement index; not a calibrated posterior probability.",
        "0-1",
        "blank means unavailable",
    ),
    "loo_total_folds": ("Number of global source-unit omission folds in the LOO diagnostic.", "non-negative integer", "equals evaluable plus penalty folds"),
    "loo_rank_evaluable_folds": ("Number of LOO folds in which the gene remained eligible and received a priority rank.", "non-negative integer", "0 means numeric LOO diagnostics are unavailable"),
    "loo_penalty_folds": ("Number of LOO folds in which the gene fell below min_studies and received the penalty-rank state.", "non-negative integer", "disjoint from evaluable folds"),
    "loo_component_available": ("Whether at least one LOO fold was rank-evaluable for this gene.", "TRUE/FALSE", "FALSE means numeric LOO diagnostics are blank"),
    "loo_median_rank": ("Median rank across rank-evaluable leave-one-source-unit-out folds.", "rank, lower is better", "penalty folds are excluded; blank means zero folds were rank-evaluable"),
    "loo_rank_iqr": ("Interquartile range across rank-evaluable leave-one-source-unit-out folds.", "rank units, lower is more stable", "penalty folds are excluded; blank means zero folds were rank-evaluable"),
    "loo_rank_stability_score": ("Rank-stability score under leave-one-source-unit-out analysis; ineligible folds remain negative evidence through the penalty-rank state.", "0-1, higher is more stable", "blank means zero folds were rank-evaluable"),
    "loo_top50_fraction": ("Fraction of rank-evaluable leave-one-out folds where the gene remains in the top 50.", "0-1", "penalty folds are excluded; blank means zero folds were rank-evaluable"),
    "loo_top100_fraction": ("Fraction of rank-evaluable leave-one-out folds where the gene remains in the top 100.", "0-1", "penalty folds are excluded; blank means zero folds were rank-evaluable"),
    "quality_weighted_degora_rank": ("Primary DEGORA rank used for prioritization in current outputs.", "1 is highest priority", "blank means not ranked"),
    "quality_weighted_degora_score": ("Primary fixed quality-weighted DEGORA score; relative index, not a probability.", "higher is stronger", "blank means not scored"),
    "quality_weighted_top_percent": ("Top-percent position for the primary quality-weighted rank.", "0-100, lower is better", "blank means not ranked"),
    "quality_weighted_consensus_direction": ("Consensus direction under quality-weighted source evidence.", "up, down, mixed, or similar text", "blank means unavailable"),
    "quality_weighted_sign_concordance": ("Quality-weighted fraction of evidence agreeing with consensus direction.", "0-1, higher is more concordant", "blank means unavailable"),
    "source_quality_support_score": ("Support contribution after source-quality weighting.", "0-1 or score units", "blank means unavailable"),
    "source_quality_weight_sum": ("Total source-quality weight contributing to the gene.", "non-negative number", "blank means unavailable"),
    "stouffer_z": ("Fixed-effect Stouffer combined z statistic.", "signed z score", "blank means unavailable"),
    "stouffer_p": ("P value from fixed-effect Stouffer combination.", "0-1", "blank means unavailable"),
    "stouffer_neglog10_p": ("Underflow-safe negative log10 of the fixed-effect Stouffer p value, used to order p-value ties.", "higher is stronger", "blank means unavailable"),
    "stouffer_padj": ("BH-adjusted Stouffer p value; a ranking aid, not a calibrated FDR.", "0-1", "blank means unavailable"),
    "heterogeneity_q": ("Cochran-like heterogeneity Q across source units.", "non-negative number", "blank means unavailable"),
    "heterogeneity_df": ("Degrees of freedom for heterogeneity summary.", "integer", "blank means unavailable"),
    "heterogeneity_i2": (HETEROGENEITY_RULE, "0-1, higher is more heterogeneous", "blank means unavailable"),
    "heterogeneity_flag": ("Text flag summarizing source-unit heterogeneity.", "text label", "blank means not flagged"),
    "re_stouffer_z": (f"Z statistic from the {RANDOM_EFFECTS_STOUFFER_RULE}.", "signed z score", "blank means unavailable"),
    "re_stouffer_p": (f"P value from the {RANDOM_EFFECTS_STOUFFER_RULE}.", "0-1", "blank means unavailable"),
    "re_stouffer_padj": (f"BH-adjusted p from the {RANDOM_EFFECTS_STOUFFER_RULE}.", "0-1", "blank means unavailable"),
    "re_stouffer_shrinkage_factor": (
        f"Divisor applied in the {RANDOM_EFFECTS_STOUFFER_RULE}.",
        ">=1, higher means more descriptive shrinkage",
        "blank means unavailable",
    ),
    "rra_rho": ("RobustRankAggreg-style rho statistic; a ranking aid, not a calibrated FDR or p-value. Use rra_neglog10_rho for top genes when rho underflows.", "0-1, lower is stronger", "blank means unavailable"),
    "rra_neglog10_rho": ("Negative log10 transform of RRA rho.", "higher is stronger", "blank means unavailable"),
    "rra_rank": ("Rank implied by RobustRankAggreg-style evidence.", "1 is strongest", "blank means unavailable"),
    "effect_meta_log2fc_re": ("Random-effects summary log2 fold change.", "log2 fold-change", "blank means insufficient effect data"),
    "effect_meta_se": ("Standard error of random-effects log2 fold-change summary.", "non-negative number", "blank means unavailable"),
    "effect_meta_ci_low": ("Lower bound of random-effects log2 fold-change confidence interval.", "log2 fold-change", "blank means unavailable"),
    "effect_meta_ci_high": ("Upper bound of random-effects log2 fold-change confidence interval.", "log2 fold-change", "blank means unavailable"),
    "effect_meta_tau2": ("Between-source variance estimate for effect summary.", "non-negative number", "blank means unavailable"),
    "effect_meta_i2": ("Inverse-variance Cochran-Q Higgins' I2 of the random-effects log2FC meta-summary.", "0-1", "blank means unavailable"),
    "effect_meta_k": ("Number of source units contributing to effect-size summary.", "integer", "blank means unavailable"),
    "effect_meta_se_source": ("How the effect-size standard error was derived.", "text label", "blank means unavailable"),
    "weighted_lfc": ("Weighted log2 fold-change summary used as an audit field.", "log2 fold-change", "blank means unavailable"),
    "rank_product": ("Rank-product style aggregate rank statistic.", "positive number, lower is stronger", "blank means unavailable"),
    "rank_score": ("Rank-based DEGORA reference score.", "higher is stronger", "blank means unavailable"),
    "slice_rank": ("Rank from the unweighted DEGORA slice/reference output.", "1 is highest priority", "blank means unavailable"),
    "high_confidence": ("Whether the gene passes the high-confidence rule used by DEGORA.", "1 or 0 (stored as integer in SQLite)", "blank means not flagged"),
    "study_id": ("Input contrast or study identifier.", "text", "not expected"),
    "source_unit_id": ("Independent collapsed source unit; related contrasts share this ID to avoid overcounting.", "text", "not expected"),
    "paper_id": ("Publication or dataset identifier used as source-unit fallback when needed.", "text", "blank means unavailable"),
    "pipeline": ("Analysis pipeline or source method for the DEG table.", "text", "blank means unavailable"),
    "assay_type": ("Assay platform class.", "RNA-seq, microarray, or text", "blank means unavailable"),
    "source_input_type": ("Type of input table used by DEGORA.", "text label", "blank means unavailable"),
    "table_scope": ("Whether the source table is genome-wide, DEG-only, or another declared scope.", "text label", "blank means unavailable"),
    "platform": ("Assay platform or accession, when reported.", "text", "blank means unavailable"),
    "normalization": ("Normalization method reported for the source.", "text", "blank means unavailable"),
    "probe_collapse": ("Probe-to-gene collapse method for microarray sources.", "text", "blank means not applicable or unavailable"),
    "species": ("Species represented by the source data.", "text", "blank means unavailable"),
    "cell_system": ("Cell type, tissue, or experimental system.", "text", "blank means unavailable"),
    "condition": (
        "Condition or perturbation descriptor retained from the source catalog.",
        "text",
        "Blank when the catalog did not record one.",
    ),
    "hypoxia_modality": ("Condition or perturbation descriptor retained from the source catalog.", "text", "blank means unavailable"),
    "duration_h": ("Treatment or exposure duration in hours when available.", "number or text", "blank means unavailable"),
    "time_course_mode": (
        "How related time-course rows were preselected before collapse: mean keeps all, early/late keep "
        "the globally smallest/largest duration_h, peak_mean keeps each gene's strongest half by "
        "p-value-derived |signed_z| rather than by fold change.",
        "text label",
        "blank means not a time-course source",
    ),
    "temporal_mode": ("Temporal aggregation mode used for source-unit collapse.", "text label", "blank means not applicable"),
    "n_ctrl": ("Number of control samples represented by the source row.", "integer", "blank means unavailable"),
    "n_treat": ("Number of treatment/test samples represented by the source row.", "integer", "blank means unavailable"),
    "lfc": ("Source-level log2 fold change after harmonization.", "log2 fold-change", "blank means unavailable"),
    "signed_z": ("Signed evidence statistic derived from source p value and direction.", "signed z score", "blank means unavailable"),
    "aggregate_pvalue": ("Two-sided p of the source-unit weighted-mean signed z; descriptive only, not a Stouffer/combined p.", "0-1", "blank means unavailable"),
    "min_source_pvalue": ("Smallest source p value among rows collapsed into the source unit.", "0-1", "blank means unavailable"),
    "min_source_padj": ("Smallest adjusted p value among collapsed source rows.", "0-1", "blank means unavailable"),
    "normalized_rank": ("Within-source normalized rank of the gene.", "0-1, lower is stronger", "blank means unavailable"),
    "n_genes_in_study": ("Declared or inferred rank universe size for the source table.", "integer", "blank means unavailable"),
    "weight": ("Source-unit analysis weight used by DEGORA.", "non-negative number", "blank means unavailable"),
    "source_quality_weight": ("Quality weight assigned to a source unit.", "non-negative number", "blank means unavailable"),
    "source_quality_label": ("Text label explaining source quality status.", "text", "blank means unavailable"),
    "source_coherence_weight": ("Weight component from source coherence diagnostics.", "non-negative number", "blank means unavailable"),
    "source_recommended_weight": ("Recommended source weight after diagnostics.", "non-negative number", "blank means unavailable"),
    "source_reliability_weight": ("Weight component from source reliability diagnostics.", "non-negative number", "blank means unavailable"),
    "source_reliability_label": ("Text label for source reliability.", "text", "blank means unavailable"),
    "source_outlier_flag": ("Whether source diagnostics flagged the row/source as an outlier.", "1 or 0 (stored as integer in SQLite)", "blank means not flagged"),
    "direction": ("Direction of gene change in the source unit.", "up or down", "blank means unavailable"),
    "source_path": ("Local input file path used for provenance.", "path text", "blank means unavailable"),
    "source_url": ("Public source URL or accession link when available.", "URL/text", "blank means unavailable"),
    "contributing_study_ids": ("Study IDs collapsed into this source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_pipelines": ("Pipelines represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_assay_types": ("Assay types represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_source_input_types": ("Input table types represented in the collapsed source-unit row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_platforms": ("Platforms represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_normalizations": ("Normalizations represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_probe_collapse": ("Probe-collapse methods represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_duration_h": ("Durations represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_time_course_modes": ("Time-course modes represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_source_paths": ("Input paths represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "contributing_source_urls": ("Source URLs represented in the collapsed source-unit evidence row.", "semicolon-delimited text", "blank means unavailable"),
    "n_contrast_rows": ("Number of contrast rows collapsed into the source-unit evidence row.", "integer", "blank means unavailable"),
    "n_studies_in_source_unit": ("Number of study IDs collapsed into the source unit.", "integer", "blank means unavailable"),
    "notes": ("Free-text notes retained from the source catalog.", "text", "blank means no note"),
    "resolved_gene_symbol": (
        "The symbol DEGORA scored this panel gene as, after legacy and Excel-date symbol repair.",
        "text",
        "differs from gene_symbol when the panel used a legacy symbol such as SEPT9",
    ),
    "present_in_degora_output": ("Whether a curated gene was present in the DEGORA ranked output.", "boolean", "false means absent from scored output"),
    "n_genes": ("Number of distinct genes contributed by this source unit.", "integer", "blank means unavailable"),
    "n_pairwise_comparisons": ("Number of other source units this one was compared against for coherence.", "integer", "blank means unavailable"),
    "median_pairwise_lfc_spearman": ("Median Spearman correlation of log2FC between this source unit and other source units.", "-1 to 1, higher is more coherent", "blank means too few overlapping genes"),
    "min_pairwise_lfc_spearman": ("Minimum pairwise log2FC Spearman correlation against other source units.", "-1 to 1", "blank means too few overlapping genes"),
    "median_pairwise_sign_agreement": ("Median fraction of overlapping genes with the same log2FC sign across source-unit pairs.", "0-1, higher is more concordant", "blank means unavailable"),
    "recommended_role": ("Suggested role for the source unit in analysis (primary or sensitivity).", "text label", "blank means unavailable"),
}


SHEET_COLUMN_OVERRIDES: dict[tuple[str, str], tuple[str, str, str]] = {
    ("Run_summary", "field"): ("Run-level metric name.", "text", "not expected"),
    ("Run_summary", "value"): ("Run-level metric value.", "text or number", "blank means unavailable"),
    ("Metadata", "field"): ("Metadata field from degora_score_metadata.json.", "text", "not expected"),
    ("Metadata", "value"): ("Serialized metadata value.", "text or number", "blank means unavailable"),
    ("SQLite_meta", "key"): ("Raw metadata key stored in the SQLite database.", "text", "not expected"),
    ("SQLite_meta", "value"): ("Raw metadata value stored in the SQLite database.", "text", "blank means unavailable"),
}


def _table_exists(db_path: Path, table_name: str) -> bool:
    if not db_path.exists() or db_path.stat().st_size == 0:
        return False
    with sqlite3.connect(db_path) as connection:
        row = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (table_name,),
        ).fetchone()
    return row is not None


def _read_table(db_path: Path, table_name: str) -> pd.DataFrame:
    if not _table_exists(db_path, table_name):
        return pd.DataFrame()
    with sqlite3.connect(db_path) as connection:
        return pd.read_sql_query(f"SELECT * FROM {table_name}", connection)


def _read_evidence_for_workbook(db_path: Path, cap: int) -> tuple[pd.DataFrame, int, bool]:
    """Read only the evidence rows that can be written to the workbook.

    The complete evidence table can be much larger than the Excel export cap.
    Apply the cap in SQLite so export memory is bounded by the workbook payload,
    while returning the uncapped row count for the audit summary.
    """

    if not _table_exists(db_path, "gene_evidence"):
        return pd.DataFrame(), 0, False
    with sqlite3.connect(db_path) as connection:
        total = int(connection.execute("SELECT COUNT(*) FROM gene_evidence").fetchone()[0])
        if cap <= 0 or total <= cap:
            return pd.read_sql_query("SELECT * FROM gene_evidence", connection), total, False

        gene_columns = {
            str(row[1]) for row in connection.execute("PRAGMA table_info(genes)").fetchall()
        }
        evidence_columns = {
            str(row[1]) for row in connection.execute("PRAGMA table_info(gene_evidence)").fetchall()
        }
        rank_column = next(
            (
                candidate
                for candidate in ("quality_weighted_degora_rank", "degora_rank")
                if candidate in gene_columns
            ),
            None,
        )
        if "gene_symbol" in gene_columns and "gene_symbol" in evidence_columns and rank_column:
            query = (
                "SELECT e.* FROM gene_evidence AS e "
                "LEFT JOIN genes AS g ON g.gene_symbol = e.gene_symbol "
                f"ORDER BY g.{rank_column} IS NULL, g.{rank_column}, e.rowid LIMIT ?"
            )
        else:
            query = "SELECT * FROM gene_evidence ORDER BY rowid LIMIT ?"
        evidence = pd.read_sql_query(query, connection, params=(cap,))
    return evidence, total, True


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return value if isinstance(value, dict) else {}


def _read_gold_from_config(config_path: Path | None) -> tuple[pd.DataFrame, str, str]:
    """Read optional curated genes without disguising a broken panel as an absent one."""

    if config_path is None:
        return pd.DataFrame(), "not_provided", "no run config was supplied for GoldPanel lookup"
    if not config_path.exists():
        return pd.DataFrame(), "not_provided", "the run config is unavailable for GoldPanel lookup"
    if config_path.suffix.lower() not in {".xlsx", ".xls"}:
        return pd.DataFrame(), "not_provided", "the run config is not an Excel workbook"
    try:
        with pd.ExcelFile(config_path) as workbook:
            if "GoldPanel" not in workbook.sheet_names:
                return pd.DataFrame(), "not_provided", "the Excel config has no GoldPanel sheet"
            gold = read_config_sheet(workbook, "GoldPanel")
    except (ImportError, KeyError, OSError, TypeError, ValueError, zipfile.BadZipFile) as exc:
        # Do not include the exception message: parser failures can echo a local path
        # or workbook cell text into an otherwise portable public manifest.
        return (
            pd.DataFrame(),
            "read_error",
            f"GoldPanel could not be read ({type(exc).__name__}); curated recall was not calculated",
        )
    if "gene_symbol" not in gold.columns:
        return (
            pd.DataFrame(),
            "invalid",
            "GoldPanel is missing the required gene_symbol column; curated recall was not calculated",
        )
    gold = gold.copy()
    named_rows = int(gold["gene_symbol"].astype("string").fillna("").str.strip().ne("").sum())
    if "locked" in gold.columns:
        gold = gold.loc[locked_panel_mask(gold["locked"])].copy()
    # gene_symbol keeps what the user typed so their panel stays recognisable;
    # resolved_gene_symbol is what DEGORA actually scores and is what the ranked
    # output is matched on. Without the second column a panel written with legacy
    # symbols reported every gene as absent.
    gold["gene_symbol"] = gold["gene_symbol"].astype("string").fillna("").str.upper().str.strip()
    gold["resolved_gene_symbol"] = [canonical_gene_symbol(value) for value in gold["gene_symbol"].tolist()]
    gold = (
        gold.loc[gold["resolved_gene_symbol"].astype("string").fillna("").ne("")]
        .drop_duplicates("resolved_gene_symbol")
        .reset_index(drop=True)
    )
    if gold.empty:
        if named_rows:
            return (
                gold,
                "not_provided",
                (
                    f"GoldPanel lists {named_rows:,} gene(s) but none are locked; set locked=yes "
                    "on the rows that should be used for curated recall"
                ),
            )
        return gold, "not_provided", "GoldPanel contains no locked gene symbols"
    return gold, "locked", ""


def _curated_lookup(gold: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    if gold.empty or genes.empty or "gene_symbol" not in genes.columns:
        return pd.DataFrame()
    rank_columns = [
        "gene_symbol",
        "quality_weighted_degora_rank",
        "quality_weighted_degora_score",
        "quality_weighted_consensus_direction",
        "quality_weighted_sign_concordance",
        "n_source_units",
        "n_contrasts_observed",
        "source_units",
    ]
    present = [column for column in rank_columns if column in genes.columns]
    if "resolved_gene_symbol" not in gold.columns:
        gold = gold.assign(
            resolved_gene_symbol=[canonical_gene_symbol(value) for value in gold["gene_symbol"].tolist()]
        )
    lookup = gold.merge(
        genes[present].rename(columns={"gene_symbol": "resolved_gene_symbol"}),
        on="resolved_gene_symbol",
        how="left",
    )
    # Coerce to a real Series first: when the rank column is absent (e.g. a legacy
    # score CSV via the fallback path), lookup.get(...) returns None and
    # pd.to_numeric(None) yields a scalar, whose .notna()/.le() would raise.
    raw_rank = (
        lookup["quality_weighted_degora_rank"]
        if "quality_weighted_degora_rank" in lookup.columns
        else pd.Series(pd.NA, index=lookup.index)
    )
    rank = pd.to_numeric(raw_rank, errors="coerce")
    lookup["present_in_degora_output"] = rank.notna()
    for cutoff in [10, 20, 50, 100]:
        lookup[f"top{cutoff}_hit"] = rank.le(cutoff).fillna(False)
    return lookup


def _summary_rows(
    result_dir: Path,
    genes: pd.DataFrame,
    evidence: pd.DataFrame,
    studies: pd.DataFrame,
    gold: pd.DataFrame,
    version_info: dict[str, str] | None = None,
    *,
    evidence_row_count: int | None = None,
    path_base: Path | None = None,
    gold_panel_status: str = "not_provided",
    gold_panel_reason: str = "",
) -> pd.DataFrame:
    version_info = version_info or runtime_version_info()
    path_base = path_base or result_dir
    rank_col = "quality_weighted_degora_rank" if "quality_weighted_degora_rank" in genes.columns else "degora_rank"
    top = (
        genes.sort_values(rank_col).head(20)["gene_symbol"].tolist()
        if rank_col in genes.columns and "gene_symbol" in genes.columns
        else []
    )
    rows: list[dict[str, Any]] = [
        {"field": "path_base", "value": "workbook_directory"},
        {"field": "result_dir", "value": portable_path(result_dir, path_base)},
        {"field": "degora_version", "value": version_info.get("degora_version", "")},
        {"field": "degora_code_revision", "value": version_info.get("degora_code_revision", "")},
        {"field": "n_scored_genes", "value": int(len(genes))},
        {
            "field": "n_gene_evidence_rows",
            "value": int(len(evidence) if evidence_row_count is None else evidence_row_count),
        },
        {
            "field": "n_source_units",
            "value": int(studies["source_unit_id"].nunique()) if "source_unit_id" in studies.columns else 0,
        },
        {"field": "n_studies", "value": int(len(studies))},
        {"field": "n_curated_genes", "value": int(len(gold))},
        {"field": "gold_panel_status", "value": gold_panel_status},
        {"field": "gold_panel_reason", "value": gold_panel_reason},
        {"field": "top_genes", "value": ";".join(map(str, top))},
    ]
    if not gold.empty and "gene_symbol" in genes.columns:
        ranked = genes.sort_values(rank_col) if rank_col in genes.columns else genes
        curated_column = "resolved_gene_symbol" if "resolved_gene_symbol" in gold.columns else "gene_symbol"
        curated = set(gold[curated_column].astype(str))
        for cutoff in [10, 20, 50, 100]:
            hits = len(set(ranked.head(cutoff)["gene_symbol"].astype(str)) & curated)
            rows.append({"field": f"curated_hits_top{cutoff}", "value": hits})
            rows.append({"field": f"curated_recall_top{cutoff}", "value": hits / len(curated) if curated else ""})
    return pd.DataFrame(rows)


def _metadata_table(metadata: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for key, value in sorted(metadata.items()):
        text = json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else "" if value is None else str(value)
        rows.append({"field": key, "value": text})
    return pd.DataFrame(rows)


def _base_sheet_name(sheet_name: str) -> str:
    if sheet_name in SHEET_GUIDE:
        return sheet_name
    for base_name in SHEET_GUIDE:
        if sheet_name.startswith(f"{base_name}_"):
            return base_name
    return sheet_name


def _fallback_definition(sheet_name: str, column: str) -> tuple[str, str, str]:
    if column.startswith("top") and column.endswith("_hit"):
        cutoff = column.removeprefix("top").removesuffix("_hit")
        return (
            f"Whether this curated gene appears within the top {cutoff} primary DEGORA ranks.",
            "boolean",
            "false means not recovered at this cutoff or absent from output",
        )
    if column.startswith("curated_hits_top"):
        cutoff = column.removeprefix("curated_hits_top")
        return (f"Number of curated genes recovered within the top {cutoff} ranks.", "integer", "blank means no GoldPanel")
    if column.startswith("curated_recall_top"):
        cutoff = column.removeprefix("curated_recall_top")
        return (f"Fraction of curated genes recovered within the top {cutoff} ranks.", "0-1", "blank means no GoldPanel")
    if column.endswith("_path"):
        return ("Path-like provenance field retained from DEGORA inputs or outputs.", "path text", "blank means unavailable")
    if column.endswith("_url"):
        return ("URL-like provenance field retained from DEGORA inputs or outputs.", "URL/text", "blank means unavailable")
    if column.startswith("n_"):
        return ("Count field exported from DEGORA for audit.", "integer", "blank means unavailable")
    if column.endswith("_score"):
        return ("Score-like DEGORA audit field.", "numeric; direction depends on field", "blank means unavailable")
    if column.endswith("_rank"):
        return ("Rank-like DEGORA audit field.", "rank; usually lower is stronger", "blank means unavailable")
    return (
        f"Exported {sheet_name} field preserved from DEGORA outputs for audit.",
        "text or numeric",
        "blank means unavailable or not applicable",
    )


def _column_definition(sheet_name: str, column: str) -> tuple[str, str, str]:
    base_name = _base_sheet_name(sheet_name)
    return (
        SHEET_COLUMN_OVERRIDES.get((base_name, column))
        or COLUMN_DEFINITIONS.get(column)
        or _fallback_definition(base_name, column)
    )


def _workbook_guide() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for sheet, info in SHEET_GUIDE.items():
        rows.append(
            {
                "sheet": sheet,
                "row_grain": info["row_grain"],
                "what_it_contains": info["what_it_contains"],
                "how_to_use": info["how_to_use"],
            }
        )
    return pd.DataFrame(rows)


def _column_dictionary(sheet_frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    dictionary_frames = {"Workbook_guide": _workbook_guide()}
    dictionary_frames.update(sheet_frames)
    for sheet, frame in dictionary_frames.items():
        row_grain = SHEET_GUIDE.get(sheet, {}).get("row_grain", "")
        for column in frame.columns:
            meaning, scale, missing = _column_definition(sheet, str(column))
            rows.append(
                {
                    "sheet": sheet,
                    "row_grain": row_grain,
                    "column": str(column),
                    "meaning": meaning,
                    "scale_or_values": scale,
                    "missing_or_blank": missing,
                }
            )
    return pd.DataFrame(rows)


def _comment_text(sheet_name: str, column: str) -> str:
    meaning, scale, missing = _column_definition(sheet_name, column)
    return f"{meaning}\nValues: {scale}\nBlank: {missing}"


def _write_sheet_chunks(writer: pd.ExcelWriter, frame: pd.DataFrame, base_name: str) -> None:
    if frame.empty:
        pd.DataFrame().to_excel(writer, sheet_name=base_name[:31], index=False)
        return
    max_data_rows = EXCEL_MAX_ROWS - 1
    if len(frame) <= max_data_rows:
        frame.to_excel(writer, sheet_name=base_name[:31], index=False)
        return
    for idx, start in enumerate(range(0, len(frame), max_data_rows), start=1):
        frame.iloc[start : start + max_data_rows].to_excel(writer, sheet_name=f"{base_name[:27]}_{idx}"[:31], index=False)


def _annotate_headers(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        if worksheet.max_row < 1:
            continue
        for cell in worksheet[1]:
            column = str(cell.value or "").strip()
            if not column:
                continue
            cell.comment = Comment(_comment_text(worksheet.title, column), COMMENT_AUTHOR)


def _force_formula_like_text(writer: pd.ExcelWriter) -> None:
    """Store formula-like strings as literal XLSX text without changing their value.

    Only text cells can be mistaken for a formula, so the numeric grid - which is
    the overwhelming majority of a DEGORA workbook - is skipped without touching
    the individual cells.
    """

    for worksheet in writer.book.worksheets:
        for row in worksheet.iter_rows():
            for cell in row:
                if cell.data_type not in {"s", "f", "str", "e"}:
                    continue
                value = cell.value
                if not isinstance(value, str):
                    continue
                stripped = value.lstrip(FORMULA_PREFIX_WHITESPACE)
                if stripped.startswith(("=", "+", "-", "@")) or stripped.upper() in EXCEL_ERROR_LITERALS:
                    cell.data_type = "s"


def _autosize(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        # Only the header cell is read, but `worksheet.columns` materialises every
        # cell of every column to get it - a second full pass over the workbook.
        for cell in worksheet[1]:
            header = str(cell.value or "")
            worksheet.column_dimensions[cell.column_letter].width = min(max(len(header) + 2, 12), 48)


def export_run_workbook(
    result_dir: Path,
    output: Path | None = None,
    *,
    config_path: Path | None = None,
    db_path: Path | None = None,
    command: str,
) -> dict[str, Any]:
    """Export a DEGORA run folder to an Excel workbook.

    The workbook is an audit convenience built from canonical DEGORA outputs:
    the ranked score CSV, SQLite evidence database, metadata JSON, and optional
    GoldPanel sheet from the run config.
    """

    result_dir = Path(result_dir)
    output = output or (result_dir / DEFAULT_WORKBOOK_NAME)
    db_path = db_path or (result_dir / "degora_scores.db")
    score_csv = result_dir / "degora_gene_scores.csv"
    metadata_json = result_dir / "degora_score_metadata.json"
    diagnostics_tsv = result_dir / "degora_source_quality_diagnostics.tsv"

    genes = _read_table(db_path, "genes")
    if genes.empty and score_csv.exists():
        genes = pd.read_csv(score_csv)
        genes = restore_formula_text_if_marked(genes, score_csv)
    evidence_sheet, evidence_total_rows, evidence_capped = _read_evidence_for_workbook(
        db_path,
        _evidence_row_cap(),
    )
    studies = _read_table(db_path, "studies")
    meta = _read_table(db_path, "meta")
    metadata = _read_json(metadata_json)
    version_info = {
        **runtime_version_info(),
        **{key: str(metadata.get(key, "")) for key in ("degora_version", "degora_code_revision") if metadata.get(key)},
    }
    diagnostics = pd.DataFrame()
    if diagnostics_tsv.exists():
        try:
            diagnostics = pd.read_csv(diagnostics_tsv, sep="\t")
            diagnostics = restore_formula_text_if_marked(diagnostics, diagnostics_tsv)
        except (pd.errors.EmptyDataError, pd.errors.ParserError):
            # Mirror the JSON/gold readers: a truncated or hand-edited TSV must not
            # abort the whole workbook export; fall back to an empty Source_quality sheet.
            diagnostics = pd.DataFrame()
    gold, gold_panel_status, gold_panel_reason = _read_gold_from_config(config_path)
    lookup = _curated_lookup(gold, genes)
    summary = _summary_rows(
        result_dir,
        genes,
        evidence_sheet,
        studies,
        gold,
        version_info,
        evidence_row_count=evidence_total_rows,
        path_base=output.parent,
        gold_panel_status=gold_panel_status,
        gold_panel_reason=gold_panel_reason,
    )
    if evidence_capped:
        summary = pd.concat(
            [
                summary,
                pd.DataFrame(
                    [
                        {"field": "gene_evidence_rows_total", "value": evidence_total_rows},
                        {"field": "gene_evidence_rows_in_workbook", "value": int(len(evidence_sheet))},
                        {
                            "field": "gene_evidence_note",
                            "value": (
                                f"Gene_evidence is limited to the top {int(len(evidence_sheet))} rows "
                                "(by gene rank) for this large corpus; the complete evidence for all genes "
                                "is in degora_scores.db (gene_evidence). degora_gene_scores.csv contains "
                                "the complete ranked-gene table, not source-level evidence."
                            ),
                        },
                    ]
                ),
            ],
            ignore_index=True,
        )
    metadata_frame = _metadata_table(metadata)
    sheet_frames = {
        "Run_summary": summary,
        "Gene_scores": genes,
        "Gene_evidence": evidence_sheet,
        "Source_units": studies,
        "Curated_lookup": lookup,
        "Source_quality": diagnostics,
        "Metadata": metadata_frame,
        "SQLite_meta": meta,
    }
    guide = _workbook_guide()
    dictionary = _column_dictionary(sheet_frames)

    output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        guide.to_excel(writer, sheet_name="Workbook_guide", index=False)
        dictionary.to_excel(writer, sheet_name="Column_dictionary", index=False)
        for sheet_name, frame in sheet_frames.items():
            _write_sheet_chunks(writer, frame, sheet_name)
        _force_formula_like_text(writer)
        _annotate_headers(writer)
        _autosize(writer)
        # Pin the workbook's own created/modified stamps before openpyxl saves them.
        set_reproducible_workbook_properties(writer.book)
    # The saved archive still carries per-member write timestamps, so rewrite it with
    # fixed timestamps and a stable member order. Identical inputs then produce a
    # byte-identical workbook, matching the CSV and SQLite outputs of the same run.
    normalize_ooxml_zip(output)

    manifest = output.with_suffix(".manifest.json")
    validation = output.with_suffix(".validation.txt")
    manifest_base = manifest.parent
    inputs = [
        path
        for path in [db_path, score_csv, metadata_json, diagnostics_tsv, config_path]
        if path is not None and path.exists()
    ]
    manifest_data = {
        **version_info,
        "generated_at": "deterministic",
        "script": "degora.excel_export.export_run_workbook",
        "path_base": "manifest_directory",
        "command": portable_command(command, manifest_base),
        "gold_panel": {
            "status": gold_panel_status,
            "reason": gold_panel_reason,
            "gene_count": int(len(gold)),
        },
        "inputs": [portable_path(path, manifest_base) for path in inputs],
        "outputs": [portable_path(path, manifest_base) for path in [output, manifest, validation]],
        "sheets": {
            "Workbook_guide": "Reader-facing guide to workbook tabs and row grain.",
            "Column_dictionary": "Definitions, scales, and missing-value notes for exported columns.",
            "Run_summary": "Run-level counts and top genes.",
            "Gene_scores": "Full DEGORA ranked gene table.",
            "Gene_evidence": "Source-unit evidence rows from the SQLite database (capped to the top-ranked genes for very large corpora; see Run_summary; full evidence is in the SQLite database and CSV).",
            "Source_units": "Source/contrast metadata.",
            "Curated_lookup": "Optional GoldPanel markers with DEGORA ranks.",
            "Source_quality": "Source-quality diagnostics, when present.",
            "Metadata": "Score metadata JSON flattened to field/value rows.",
            "SQLite_meta": "Raw meta table from SQLite, when present.",
        },
    }
    manifest.write_text(json.dumps(manifest_data, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    validation.write_text(
        "\n".join(
            [
                f"n_scored_genes={len(genes)}",
                f"degora_version={version_info.get('degora_version', '')}",
                f"degora_code_revision={version_info.get('degora_code_revision', '')}",
                f"n_gene_evidence_rows={evidence_total_rows}",
                f"n_source_units={studies['source_unit_id'].nunique() if 'source_unit_id' in studies.columns else 0}",
                f"n_curated_genes={len(gold)}",
                f"gold_panel_status={gold_panel_status}",
                f"gold_panel_reason={gold_panel_reason}",
            ]
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    sidecar_metadata = {"generator": "degora-run-workbook", **version_info}
    for artifact in [output, manifest, validation]:
        write_source_sidecar(artifact, command, inputs=inputs, metadata=sidecar_metadata)
    return {
        "output": output.as_posix(),
        "manifest": manifest.as_posix(),
        "validation": validation.as_posix(),
        "rows_gene_scores": int(len(genes)),
        "rows_gene_evidence": evidence_total_rows,
        "rows_curated_lookup": int(len(lookup)),
        "gold_panel_status": gold_panel_status,
        "gold_panel_reason": gold_panel_reason,
    }
