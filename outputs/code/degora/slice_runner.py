"""Command implementation for the iteration-1 vertical slice."""

from __future__ import annotations

import argparse
import difflib
import json
from pathlib import Path
from typing import Any

import pandas as pd

from .aggregate import slice_consensus
from .harmonize import TableMapping, harmonize_frame, normalize_table_scope, read_deg_table
from .metrics import recall_at_k
from .provenance import shell_command, write_source_sidecar


CATALOG_COLUMNS = [
    "study_id",
    "paper_id",
    "source_url",
    "source_path",
    "pipeline",
    "species",
    "cell_system",
    "hypoxia_modality",
    "duration_h",
    "n_ctrl",
    "n_treat",
    "gene_column",
    "lfc_column",
    "p_column",
    "padj_column",
    "sep",
    "sheet_name",
    "gene_type_column",
    "gene_type_keep",
    "assay_type",
    "source_input_type",
    "platform",
    "normalization",
    "probe_id_column",
    "probe_collapse",
    "time_course_mode",
    "temporal_mode",
    "table_scope",
    "rank_universe_size",
    "include_in_analysis",
    "notes",
]

CATALOG_ALIASES = {
    "source_unit_id": "paper_id",
    "condition": "hypoxia_modality",
    "contrast_label": "hypoxia_modality",
    "time_h": "duration_h",
    "result_scope": "table_scope",
    "input_scope": "table_scope",
    "tested_gene_count": "rank_universe_size",
    "rank_universe": "rank_universe_size",
    "technology": "assay_type",
    "assay": "assay_type",
    "input_type": "source_input_type",
    "platform_id": "platform",
    "array_platform": "platform",
    "probe_column": "probe_id_column",
    "probe_collapse_rule": "probe_collapse",
    "time_mode": "time_course_mode",
    "temporal_aggregation": "time_course_mode",
    "include": "include_in_analysis",
}

ESSENTIAL_CATALOG_COLUMNS = [
    "study_id",
    "paper_id",
    "source_path",
    "gene_column",
    "lfc_column",
    "p_column",
]

CATALOG_COLUMN_HELP = {
    "study_id": "unique row ID for one DEG contrast, such as IFN_GSE001_4h",
    "paper_id": "independent source-unit ID; reuse it for time-course rows from the same paper/dataset",
    "source_path": "local path to the DEG table file",
    "gene_column": "column in the DEG table that contains gene symbols or IDs",
    "lfc_column": "column in the DEG table that contains log2 fold change",
    "p_column": "column in the DEG table that contains nominal p-value; use padj only when no p-value exists and note it",
    "padj_column": "optional adjusted p-value/FDR column",
    "table_scope": "auto, full_results, or deg_only; use deg_only when the table only lists significant genes",
    "rank_universe_size": "optional number of genes originally tested; important for DEG-only lists",
    "assay_type": "RNA-seq, microarray, proteomics, or other source assay; blank means unknown",
    "source_input_type": "author_deg_table, derived_count_table, normalized_expression_matrix, or similar",
    "platform": "microarray platform such as GPL570, sequencing platform, or blank if not needed",
    "normalization": "normalization used by the source, e.g. RMA/log2, quantile/log2, DESeq2, or edgeR",
    "probe_collapse": "for microarray sources, how probes were collapsed to gene symbols",
    "time_course_mode": "mean, early, late, or peak_mean for same-source time-course contrasts",
    "include_in_analysis": "yes/no flag; blank means yes",
}

OPTIONAL_CATALOG_DEFAULTS = {
    "source_url": "",
    "pipeline": "unknown_pipeline",
    "species": "",
    "cell_system": "",
    "hypoxia_modality": "",
    "duration_h": "",
    "n_ctrl": "",
    "n_treat": "",
    "padj_column": "",
    "sep": "",
    "sheet_name": "",
    "gene_type_column": "",
    "gene_type_keep": "",
    "assay_type": "",
    "source_input_type": "",
    "platform": "",
    "normalization": "",
    "probe_id_column": "",
    "probe_collapse": "",
    "time_course_mode": "mean",
    "temporal_mode": "",
    "table_scope": "auto",
    "rank_universe_size": "",
    "include_in_analysis": True,
    "notes": "",
}

TIME_COURSE_MODE_ALIASES = {
    "": "mean",
    "auto": "mean",
    "all": "mean",
    "source_mean": "mean",
    "average": "mean",
    "mean": "mean",
    "first": "early",
    "earliest": "early",
    "early": "early",
    "last": "late",
    "latest": "late",
    "late": "late",
    "peak": "peak_mean",
    "peak_mean": "peak_mean",
    "strongest_window": "peak_mean",
}


class DegoraConfigError(ValueError):
    """Beginner-readable configuration error with concrete repair hints."""

    def __init__(
        self,
        title: str,
        *,
        problems: list[str],
        fixes: list[str] | None = None,
        context: str | None = None,
    ) -> None:
        self.title = title
        self.problems = problems
        self.fixes = fixes or []
        self.context = context
        super().__init__(self._format())

    def _format(self) -> str:
        lines = [f"DEGORA config error: {self.title}"]
        if self.context:
            lines.extend(["", f"Context: {self.context}"])
        if self.problems:
            lines.extend(["", "Problems:"])
            lines.extend(f"- {problem}" for problem in self.problems)
        if self.fixes:
            lines.extend(["", "How to fix:"])
            lines.extend(f"- {fix}" for fix in self.fixes)
        return "\n".join(lines)


def _format_columns(columns: list[Any]) -> str:
    return ", ".join(map(str, columns)) if columns else "(no columns found)"


def _user_row_number(index: Any) -> str:
    if isinstance(index, int):
        return str(index + 2)
    return str(index)


def _nonempty(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value)
    if text in {r"\t", "\t"}:
        return "\t"
    text = text.strip()
    return text if text != "" else None


def _normalize_time_course_setting(value: Any) -> str | None:
    text = _nonempty(value)
    label = "" if text is None else text.strip().lower().replace("-", "_").replace(" ", "_")
    return TIME_COURSE_MODE_ALIASES.get(label)


def _count_labels(series: pd.Series, *, unknown_label: str = "unknown") -> dict[str, int]:
    """Return JSON-safe value counts with missing or blank labels collapsed."""

    labels = series.map(lambda value: _nonempty(value) or unknown_label)
    return {str(label): int(count) for label, count in labels.value_counts(dropna=False).items()}


def _read_catalog_frame(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        with pd.ExcelFile(path) as workbook:
            sheet_names = workbook.sheet_names
            if "Contrasts" in sheet_names:
                return pd.read_excel(workbook, sheet_name="Contrasts")
            if len(sheet_names) == 1:
                return pd.read_excel(workbook, sheet_name=sheet_names[0])
        raise DegoraConfigError(
            "Excel workbook is missing a Contrasts sheet",
            context=f"config file: {path}",
            problems=[
                f"Found sheets: {_format_columns(sheet_names)}",
                "DEGORA reads DEG input rows from a sheet named 'Contrasts'.",
            ],
            fixes=[
                "Rename the sheet that lists DEG tables to 'Contrasts'.",
                "Keep Project and AdvancedSettings sheets if you want them; GoldPanel is used only for optional locked recall metrics.",
            ],
        )
    return pd.read_csv(path)


def _read_locked_gold_panel(path: Path) -> dict[str, Any]:
    """Read an optional locked gold panel from beginner Excel configs."""

    if path.suffix.lower() not in {".xlsx", ".xls"}:
        return {
            "status": "not_provided",
            "source": "",
            "gene_column": "",
            "genes": [],
            "reason": "catalog is not an Excel workbook with a GoldPanel sheet",
        }
    with pd.ExcelFile(path) as workbook:
        if "GoldPanel" not in workbook.sheet_names:
            return {
                "status": "not_provided",
                "source": "GoldPanel",
                "gene_column": "",
                "genes": [],
                "reason": "Excel workbook has no GoldPanel sheet",
            }
        gold = pd.read_excel(workbook, sheet_name="GoldPanel")
    if "gene_symbol" not in gold.columns:
        return {
            "status": "invalid",
            "source": "GoldPanel",
            "gene_column": "gene_symbol",
            "genes": [],
            "reason": "GoldPanel sheet is missing the gene_symbol column",
        }
    if "locked" in gold.columns:
        locked = gold["locked"].astype("string").fillna("").str.strip().str.lower()
        gold = gold.loc[locked.isin({"1", "true", "t", "yes", "y", "locked"}) | locked.eq("")]
    genes = (
        gold["gene_symbol"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .loc[lambda series: series.ne("")]
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    return {
        "status": "locked" if genes else "not_provided",
        "source": "GoldPanel",
        "gene_column": "gene_symbol",
        "genes": genes,
        "reason": "" if genes else "GoldPanel contains no locked gene symbols",
    }


def _normalize_catalog_columns(catalog: pd.DataFrame) -> pd.DataFrame:
    frame = catalog.copy()
    for alias, canonical in CATALOG_ALIASES.items():
        if canonical not in frame.columns and alias in frame.columns:
            frame[canonical] = frame[alias]
    for column, default in OPTIONAL_CATALOG_DEFAULTS.items():
        if column not in frame.columns:
            frame[column] = default
    return frame


def _validate_catalog_columns(catalog: pd.DataFrame, path: Path) -> None:
    missing = [column for column in ESSENTIAL_CATALOG_COLUMNS if column not in catalog.columns]
    if not missing:
        return
    problems = [
        f"Missing required column '{column}': {CATALOG_COLUMN_HELP.get(column, 'required catalog field')}."
        for column in missing
    ]
    raise DegoraConfigError(
        "catalog is missing required column(s)",
        context=f"config file: {path}; available columns: {_format_columns(list(catalog.columns))}",
        problems=problems,
        fixes=[
            "Add the missing column(s) to the CSV file or to the Excel Contrasts sheet.",
            "For time-course data, put one row per time point and reuse paper_id/source_unit_id for related rows.",
            "Supported beginner aliases include source_unit_id->paper_id, time_h->duration_h, condition->hypoxia_modality, and include->include_in_analysis.",
        ],
    )


def _validate_catalog_required_values(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    problems: list[str] = []
    active = catalog.loc[include_mask]
    for index, row in active.iterrows():
        for column in ESSENTIAL_CATALOG_COLUMNS:
            if _nonempty(row.get(column)) is None:
                problems.append(
                    f"Row {_user_row_number(index)} has an empty '{column}' value "
                    f"({CATALOG_COLUMN_HELP.get(column, 'required catalog field')})."
                )
    if not problems:
        return
    raise DegoraConfigError(
        "active contrast row(s) have empty required values",
        context=f"config file: {path}",
        problems=problems,
        fixes=[
            "Fill the highlighted cells in the Contrasts sheet.",
            "If a row is only a note or placeholder, set include_in_analysis/include to 'no'.",
        ],
    )


def _validate_optional_scope_values(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    problems: list[str] = []
    fixes: list[str] = []
    active = catalog.loc[include_mask]
    for index, row in active.iterrows():
        try:
            normalize_table_scope(row.get("table_scope", "auto"))
        except ValueError as exc:
            problems.append(
                f"Row {_user_row_number(index)} has unsupported table_scope={row.get('table_scope')!r}: {exc}"
            )
            fixes.append("Use table_scope=auto, full_results, deg_only, or ambiguous.")

        time_course_mode = row.get("time_course_mode", "mean")
        if _normalize_time_course_setting(time_course_mode) is None:
            problems.append(
                f"Row {_user_row_number(index)} has time_course_mode={time_course_mode!r}; "
                "it must be mean, early, late, peak_mean, or blank."
            )
            fixes.append(
                "Use time_course_mode=mean unless you predeclare that the source unit should use early, late, or peak_mean time points."
            )

        raw_universe = row.get("rank_universe_size", "")
        if _nonempty(raw_universe) is not None:
            universe = pd.to_numeric(pd.Series([raw_universe]), errors="coerce").iloc[0]
            if pd.isna(universe) or float(universe) <= 0:
                problems.append(
                    f"Row {_user_row_number(index)} has rank_universe_size={raw_universe!r}; "
                    "it must be a positive number or blank."
                )
                fixes.append(
                    "For DEG-only lists, enter the number of genes originally tested if the paper reports it; otherwise leave blank."
                )

    if problems:
        raise DegoraConfigError(
            "table-scope settings are not valid",
            context=f"config file: {path}",
            problems=problems,
            fixes=fixes,
        )


def _microarray_warnings(catalog: pd.DataFrame) -> list[str]:
    """Return non-fatal warnings for active microarray metadata."""

    if "assay_type" not in catalog.columns:
        return []
    warnings: list[str] = []
    for _, row in catalog.iterrows():
        assay = str(row.get("assay_type", "")).strip().lower()
        if assay != "microarray":
            continue
        study_id = str(row.get("study_id", ""))
        pipeline = str(row.get("pipeline", "")).strip().lower()
        if "limma" not in pipeline and "welch_microarray" not in pipeline:
            warnings.append(
                f"{study_id}: assay_type=microarray is best paired with pipeline=limma_microarray "
                "or welch_microarray_normalized_matrix."
            )
        for column in ["platform", "normalization", "probe_collapse"]:
            if _nonempty(row.get(column)) is None:
                warnings.append(f"{study_id}: microarray row is missing {column}; keep this source exploratory until documented.")
        if _nonempty(row.get("source_input_type")) is None:
            warnings.append(
                f"{study_id}: microarray row is missing source_input_type; specify author_deg_table, "
                "limma_full_table, or normalized_expression_matrix."
            )
    return warnings


def read_catalog(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    catalog = _normalize_catalog_columns(_read_catalog_frame(path))
    _validate_catalog_columns(catalog, path)
    return catalog[CATALOG_COLUMNS]


def catalog_include_mask(catalog: pd.DataFrame) -> pd.Series:
    """Return the active-analysis mask from the pre-registered catalog flag."""

    truthy = {"1", "true", "t", "yes", "y", "include", "included"}
    falsy = {"0", "false", "f", "no", "n", "exclude", "excluded"}
    values = catalog["include_in_analysis"]
    mask = []
    for index, value in values.items():
        if value is None or pd.isna(value) or str(value).strip() == "":
            mask.append(True)
            continue
        text = str(value).strip().lower()
        if text in truthy:
            mask.append(True)
        elif text in falsy:
            mask.append(False)
        else:
            raise DegoraConfigError(
                "include flag has an unsupported value",
                problems=[
                    f"Row {_user_row_number(index)} has include_in_analysis={value!r}.",
                    "DEGORA only accepts yes/no-style include flags.",
                ],
                fixes=[
                    "Use yes/true/include/1 to keep a row.",
                    "Use no/false/exclude/0 to exclude a row.",
                    "Leave the cell blank if the row should be included.",
                ],
            )
    return pd.Series(mask, index=catalog.index, dtype=bool)


def apply_gene_type_filter(frame: pd.DataFrame, column: str | None, keep: str | None) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Filter a source table by gene type when the catalog explicitly requests it."""

    summary = {
        "applied": False,
        "column": column,
        "keep": keep,
        "rows_before": int(len(frame)),
        "rows_after": int(len(frame)),
    }
    if not column or not keep:
        return frame, summary
    if column not in frame.columns:
        suggestion = difflib.get_close_matches(column, list(map(str, frame.columns)), n=1)
        fix = f"Did you mean '{suggestion[0]}'?" if suggestion else "Use one of the available column names exactly."
        raise DegoraConfigError(
            "gene type filter column is missing from a source table",
            problems=[
                f"gene_type_column={column!r} was requested, but the source table does not contain that column.",
                f"Available columns: {_format_columns(list(frame.columns))}",
            ],
            fixes=[
                fix,
                "If you do not need a gene-type filter, leave gene_type_column and gene_type_keep blank.",
            ],
        )

    keep_values = {value.strip().lower() for value in keep.split("|") if value.strip()}
    mask = frame[column].astype("string").str.strip().str.lower().isin(keep_values)
    filtered = frame.loc[mask].copy()
    summary.update({"applied": True, "rows_after": int(len(filtered))})
    return filtered, summary


def _validate_source_columns(frame: pd.DataFrame, mapping: TableMapping, row: dict[str, Any]) -> None:
    requested = [
        ("gene_column", mapping.gene_column, "gene symbols or gene IDs"),
        ("lfc_column", mapping.lfc_column, "log2 fold change"),
        ("p_column", mapping.p_column, "p-value"),
    ]
    if mapping.padj_column:
        requested.append(("padj_column", mapping.padj_column, "adjusted p-value/FDR"))

    available = list(map(str, frame.columns))
    problems: list[str] = []
    fixes: list[str] = []
    for catalog_column, source_column, meaning in requested:
        if source_column in frame.columns:
            continue
        suggestion = difflib.get_close_matches(str(source_column), available, n=1)
        problems.append(
            f"{row['study_id']}: {catalog_column}={source_column!r} should identify the source-table column for {meaning}, "
            "but that column was not found."
        )
        if suggestion:
            fixes.append(f"{row['study_id']}: change {catalog_column} to '{suggestion[0]}'.")
    if problems:
        raise DegoraConfigError(
            "source table column mapping is wrong",
            context=f"source file: {row.get('source_path', '')}; available columns: {_format_columns(available)}",
            problems=problems,
            fixes=fixes
            + [
                "Open the DEG table and copy the column name exactly into the Contrasts sheet.",
                "Column names are case-sensitive.",
            ],
        )


def _resolve_source_path(raw_source_path: Any, catalog_path: Path, catalog_repo_root: Path) -> Path:
    source_path = Path(str(raw_source_path))
    if source_path.is_absolute():
        return source_path

    candidates = [
        Path.cwd() / source_path,
        catalog_path.parent / source_path,
        catalog_repo_root / source_path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return catalog_repo_root / source_path


def _infer_project_root(catalog_path: Path) -> Path:
    """Find the repository root for resolving beginner config paths."""

    for parent in [catalog_path.parent, *catalog_path.parents]:
        if (parent / "outputs" / "code").exists() and (parent / "data").exists():
            return parent
    if len(catalog_path.parents) > 2:
        return catalog_path.parents[2]
    return Path.cwd()


def validate_catalog_inputs(catalog_path: Path) -> dict[str, Any]:
    """Validate catalog/config rows and source-table column mappings without writing outputs."""

    catalog_path = catalog_path.resolve()
    catalog_repo_root = _infer_project_root(catalog_path)
    full_catalog = read_catalog(catalog_path)
    include_mask = catalog_include_mask(full_catalog)
    _validate_catalog_required_values(full_catalog, include_mask, catalog_path)
    _validate_optional_scope_values(full_catalog, include_mask, catalog_path)
    catalog = full_catalog.loc[include_mask].copy()
    if catalog.empty:
        raise DegoraConfigError(
            "catalog contains no active contrasts",
            problems=["Every row is excluded by include_in_analysis/include, or the catalog is empty."],
            fixes=["Set include_in_analysis/include to yes for at least one contrast row."],
        )

    source_units = set()
    checked_sources: list[str] = []
    for row in catalog.to_dict(orient="records"):
        source_path = _resolve_source_path(row["source_path"], catalog_path, catalog_repo_root)
        if not source_path.exists():
            raise DegoraConfigError(
                "source DEG table file was not found",
                problems=[
                    f"{row['study_id']}: source_path points to {source_path}, but that file does not exist.",
                ],
                fixes=[
                    "Check the file path in the Contrasts sheet.",
                    "Relative paths can point from your current folder, the Excel file folder, or the project root.",
                    "If the file was moved, update source_path rather than editing analysis outputs by hand.",
                ],
            )

        mapping = TableMapping(
            gene_column=row["gene_column"],
            lfc_column=row["lfc_column"],
            p_column=row["p_column"],
            padj_column=_nonempty(row.get("padj_column")),
            sep=_nonempty(row.get("sep")),
            sheet_name=_nonempty(row.get("sheet_name")),
        )
        try:
            raw_frame = read_deg_table(source_path, mapping)
        except Exception as exc:
            raise DegoraConfigError(
                "source DEG table could not be read",
                context=f"{row['study_id']}: {source_path}",
                problems=[str(exc)],
                fixes=[
                    "Check that the file is a supported CSV/TSV/TXT/XLS/XLSX table.",
                    "For Excel sources, set sheet_name to the exact sheet containing the DEG table.",
                    "For TSV files, set sep to \\t if DEGORA cannot infer it.",
                ],
            ) from exc
        _validate_source_columns(raw_frame, mapping, row)
        apply_gene_type_filter(
            raw_frame,
            _nonempty(row.get("gene_type_column")),
            _nonempty(row.get("gene_type_keep")),
        )
        paper_id = _nonempty(row.get("paper_id"))
        study_id = _nonempty(row.get("study_id"))
        source_units.add(paper_id or study_id)
        checked_sources.append(str(source_path))

    return {
        "config_path": str(catalog_path),
        "active_contrasts": int(len(catalog)),
        "excluded_contrasts": int((~include_mask).sum()),
        "source_units": int(len(source_units)),
        "checked_sources": checked_sources,
    }


def run_slice(catalog_path: Path, output_dir: Path, harmonized_dir: Path, min_studies: int) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    harmonized_dir.mkdir(parents=True, exist_ok=True)

    catalog_path = catalog_path.resolve()
    catalog_repo_root = _infer_project_root(catalog_path)
    full_catalog = read_catalog(catalog_path)
    include_mask = catalog_include_mask(full_catalog)
    _validate_catalog_required_values(full_catalog, include_mask, catalog_path)
    _validate_optional_scope_values(full_catalog, include_mask, catalog_path)
    catalog = full_catalog.loc[include_mask].copy()
    excluded_catalog = full_catalog.loc[~include_mask].copy()
    if catalog.empty:
        raise ValueError("Catalog contains no active studies")

    harmonized_tables = []
    source_inputs: list[Path] = []
    input_warnings: list[str] = _microarray_warnings(catalog)
    filter_summaries: dict[str, dict[str, Any]] = {}

    for row in catalog.to_dict(orient="records"):
        source_path = _resolve_source_path(row["source_path"], catalog_path, catalog_repo_root)
        if not source_path.exists():
            raise DegoraConfigError(
                "source DEG table file was not found",
                problems=[
                    f"{row['study_id']}: source_path points to {source_path}, but that file does not exist.",
                ],
                fixes=[
                    "Check the file path in the Contrasts sheet.",
                    "Relative paths can point from your current folder, the Excel file folder, or the project root.",
                    "If the file was moved, update source_path rather than editing analysis outputs by hand.",
                ],
            )
        source_inputs.append(source_path)

        mapping = TableMapping(
            gene_column=row["gene_column"],
            lfc_column=row["lfc_column"],
            p_column=row["p_column"],
            padj_column=_nonempty(row.get("padj_column")),
            sep=_nonempty(row.get("sep")),
            sheet_name=_nonempty(row.get("sheet_name")),
        )
        try:
            raw_frame = read_deg_table(source_path, mapping)
        except Exception as exc:
            raise DegoraConfigError(
                "source DEG table could not be read",
                context=f"{row['study_id']}: {source_path}",
                problems=[str(exc)],
                fixes=[
                    "Check that the file is a supported CSV/TSV/TXT/XLS/XLSX table.",
                    "For Excel sources, set sheet_name to the exact sheet containing the DEG table.",
                    "For TSV files, set sep to \\t if DEGORA cannot infer it.",
                ],
            ) from exc
        _validate_source_columns(raw_frame, mapping, row)
        filtered_frame, filter_summary = apply_gene_type_filter(
            raw_frame,
            _nonempty(row.get("gene_type_column")),
            _nonempty(row.get("gene_type_keep")),
        )
        filter_summaries[row["study_id"]] = filter_summary
        row["source_path"] = str(source_path)
        harmonized = harmonize_frame(filtered_frame, mapping, row)
        if harmonized.empty:
            input_warnings.append(f"{row['study_id']} produced zero harmonized rows")
        harmonized_tables.append(harmonized)

    if not harmonized_tables:
        raise ValueError("Catalog contains no usable studies")

    all_harmonized = pd.concat(harmonized_tables, ignore_index=True)
    harmonized_stem = f"{output_dir.name}_harmonized"
    harmonized_csv = harmonized_dir / f"{harmonized_stem}.csv"
    harmonized_parquet = harmonized_dir / f"{harmonized_stem}.parquet"
    all_harmonized.to_csv(harmonized_csv, index=False)
    all_harmonized.to_parquet(harmonized_parquet, index=False)

    consensus = slice_consensus(all_harmonized, min_studies=min_studies)
    consensus_path = output_dir / "slice_consensus.csv"
    consensus.to_csv(consensus_path, index=False)

    result_harmonized_path = output_dir / "slice_harmonized.csv"
    all_harmonized.to_csv(result_harmonized_path, index=False)
    rank_universe_warnings = sorted(
        {
            str(value)
            for value in all_harmonized.get("rank_universe_warning", pd.Series(dtype=str)).dropna().unique()
            if str(value).strip()
        }
    )
    gene_symbol_collapse_warnings = sorted(
        {
            str(value)
            for value in all_harmonized.get("gene_symbol_collapse_warning", pd.Series(dtype=str)).dropna().unique()
            if str(value).strip()
        }
    )
    input_warnings.extend(gene_symbol_collapse_warnings)

    gold_panel = _read_locked_gold_panel(catalog_path)
    if gold_panel["status"] == "locked":
        recall50 = recall_at_k(consensus, gold_panel["genes"], 50)
        recall100 = recall_at_k(consensus, gold_panel["genes"], 100)
    else:
        recall50 = {
            "status": "not_applicable",
            "reason": "no locked topic-specific GoldPanel was provided; do not interpret hypoxia-specific recall for this run",
        }
        recall100 = {
            "status": "not_applicable",
            "reason": "no locked topic-specific GoldPanel was provided; do not interpret hypoxia-specific recall for this run",
        }

    source_inputs_sorted = [Path(value) for value in sorted({str(path) for path in source_inputs})]
    metrics = {
        "catalog_path": str(catalog_path),
        "source_input_files": [str(path) for path in source_inputs_sorted],
        "n_catalog_rows": int(len(full_catalog)),
        "n_active_catalog_rows": int(len(catalog)),
        "excluded_catalog_rows": excluded_catalog[["study_id", "notes"]].to_dict(orient="records"),
        "n_harmonized_rows": int(len(all_harmonized)),
        "n_consensus_genes": int(len(consensus)),
        "min_studies": min_studies,
        "pvalue_clipped_rows": int(all_harmonized["pvalue_was_clipped"].sum()),
        "study_row_counts": all_harmonized.groupby("study_id").size().astype(int).to_dict(),
        "source_filter_summary": filter_summaries,
        "table_scope_counts": all_harmonized[["study_id", "table_scope"]]
        .drop_duplicates()
        .groupby("table_scope")
        .size()
        .astype(int)
        .to_dict()
        if "table_scope" in all_harmonized.columns
        else {},
        "rank_universe_warnings": rank_universe_warnings,
        "gene_symbol_collapse_warnings": gene_symbol_collapse_warnings,
        "pipeline_counts": _count_labels(catalog["pipeline"]),
        "assay_type_counts": _count_labels(catalog["assay_type"]) if "assay_type" in catalog.columns else {},
        "source_input_type_counts": _count_labels(catalog["source_input_type"])
        if "source_input_type" in catalog.columns
        else {},
        "warnings": input_warnings,
        "gold_panel_status": gold_panel["status"],
        "gold_panel_source": gold_panel["source"],
        "gold_panel_gene_count": len(gold_panel["genes"]),
        "gold_panel_reason": gold_panel["reason"],
        "recall_at_50": recall50,
        "recall_at_100": recall100,
    }
    metrics_path = output_dir / "slice_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    command = shell_command(
        [
            "make",
            "-C",
            "outputs/code",
            "slice",
            f"CATALOG={catalog_path}",
            f"OUTDIR={output_dir}",
            f"HARMONIZED_DIR={harmonized_dir}",
            f"SLICE_MIN_STUDIES={min_studies}",
        ]
    )
    for artifact in (harmonized_csv, harmonized_parquet, consensus_path, result_harmonized_path, metrics_path):
        write_source_sidecar(
            artifact,
            command,
            inputs=[catalog_path, *source_inputs_sorted],
            metadata={"generator": "slice", "min_studies": min_studies},
        )

    return metrics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--harmonized-dir", type=Path, required=True)
    parser.add_argument("--min-studies", type=int, default=2)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    metrics = run_slice(args.catalog, args.output_dir, args.harmonized_dir, args.min_studies)
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
