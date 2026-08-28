"""Command implementation for the iteration-1 vertical slice."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import os
import re
import stat
import zipfile
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from .aggregate import (
    TIME_COURSE_RETENTION_WARNING_FRACTION,
    _source_unit_series,
    duration_hours,
    slice_consensus,
    time_course_selection_report,
    time_course_selection_warnings,
    validate_min_studies,
)
from .excel_io import locked_panel_mask, read_config_sheet
from .identifiers import UNKNOWN_IDENTIFIER_SPACE, identifier_space, identifier_space_profile
from .formula_safety import (
    formula_guard_metadata,
    neutralize_formula_text,
    restore_formula_text_if_marked,
)
from .harmonize import (
    is_workbook_path,
    TableMapping,
    bounded_value_examples,
    canonical_gene_symbol,
    duplicate_source_headers,
    harmonize_frame,
    normalize_header_row,
    normalize_lfc_scale,
    normalize_table_scope,
    read_deg_table,
    resolve_column_name,
    validate_table_mapping_roles,
)
from .metrics import recall_at_k
from .provenance import artifact_output_lock, portable_path, shell_command, write_source_sidecar


CATALOG_COLUMNS = [
    "study_id",
    "paper_id",
    "source_unit_id",
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
    "header_row",
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
    "sign_convention",
    "lfc_scale",
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
    "temporal_mode": "time_course_mode",
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

BEGINNER_REQUIRED_CONTRAST_COLUMNS = [
    "study_id",
    "source_unit_id (or paper_id)",
    "source_path",
    "gene_column",
    "lfc_column",
    "p_column",
]

REQUIRED_SOURCE_TABLE_MAPPINGS = [
    ("gene_column", "gene symbols or gene IDs"),
    ("lfc_column", "numeric log2 fold change"),
    ("p_column", "numeric p-value in [0, 1]"),
]

OPTIONAL_SOURCE_TABLE_MAPPINGS = [
    ("padj_column", "adjusted p-value/FDR in [0, 1]"),
]

CATALOG_COLUMN_HELP = {
    "study_id": "unique row ID for one DEG contrast, such as IFN_GSE001_4h",
    "paper_id": "independent source-unit ID; use source_unit_id in the beginner Excel sheet",
    "source_path": "local path to the DEG table file",
    "gene_column": "column in the DEG table that contains gene symbols or IDs",
    "lfc_column": "column in the DEG table that contains log2 fold change",
    "p_column": "column in the DEG table that contains nominal p-value in [0, 1]",
    "padj_column": "optional adjusted p-value/FDR column; leave blank if unavailable",
    "table_scope": "auto, full_results, or deg_only; use deg_only when the table only lists significant genes",
    "rank_universe_size": "optional number of genes originally tested; important for DEG-only lists",
    "sign_convention": "optional verified effect direction, including any explicit sign inversion applied before ingest",
    "lfc_scale": "optional; write log2 to state that lfc_column is already a log2 fold change when its header does not say so",
    "assay_type": "RNA-seq, microarray, proteomics, or other source assay; blank means unknown",
    "source_input_type": "author_deg_table, derived_count_table, normalized_expression_matrix, or similar",
    "platform": "microarray platform such as GPL570, sequencing platform, or blank if not needed",
    "normalization": "normalization used by the source, e.g. RMA/log2, quantile/log2, DESeq2, or edgeR",
    "probe_collapse": "for microarray sources, how probes were collapsed to gene symbols",
    "time_course_mode": (
        "how same-source time-course contrasts are preselected: mean keeps all, early/late keep the "
        "globally smallest/largest duration_h, peak_mean keeps each gene's strongest half by "
        "p-value-derived |signed_z| (not by fold change)"
    ),
    "include_in_analysis": "yes/no flag; blank means yes",
    "header_row": "optional 1-based row number of the line holding the column names; blank means the first row",
}

OPTIONAL_CATALOG_DEFAULTS = {
    "source_unit_id": "",
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
    "header_row": "",
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
    "sign_convention": "",
    "lfc_scale": "",
    "include_in_analysis": True,
    "notes": "",
}

# A declared biological group size above this is a typing slip, not a design.
# `degora/discovery_run.py` applies the same bound to the browser's review panel;
# the catalog path had no equivalent check until 0.4.33.
MAX_DECLARED_GROUP_SIZE = 10_000

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


def _is_external_absolute_path(path: Path, catalog_path: Path) -> bool:
    if not path.is_absolute():
        return False
    try:
        path.resolve().relative_to(catalog_path.parent.resolve())
    except ValueError:
        return True
    return False


def _source_context(row: dict[str, Any], source_path: Path, catalog_path: Path, available: list[Any] | None = None) -> str:
    if _is_external_absolute_path(source_path, catalog_path):
        # Both the path and the column names stay hidden. A column name can itself
        # be private - an assay name, a cohort label - which is why the test for
        # this uses a sentinel header rather than an ordinary one.
        return (
            f"{row.get('study_id', 'source table')}: source_path is outside the config folder; "
            "path and available columns hidden"
        )
    if available is None:
        return f"source file: {row.get('source_path', '')}"
    return f"source file: {row.get('source_path', '')}; available columns: {_format_columns(list(available))}"


def _source_path_context(row: dict[str, Any], source_path: Path, catalog_path: Path) -> str:
    if _is_external_absolute_path(source_path, catalog_path):
        name = source_path.name or "external source file"
        return f"{row.get('study_id', 'source table')}: external source file {name}"
    return f"{row.get('study_id', 'source table')}: {source_path}"


def _source_path_problem(row: dict[str, Any], source_path: Path, catalog_path: Path, message: str) -> str:
    if _is_external_absolute_path(source_path, catalog_path):
        name = source_path.name or "external source file"
        return f"{row['study_id']}: source_path points to external source file {name}, {message}"
    return f"{row['study_id']}: source_path points to {source_path}, {message}"


def _readable_source_read_failure(source_path: Path, catalog_path: Path, exc: Exception) -> str:
    message = _readable_read_failure(source_path, exc)
    if not _is_external_absolute_path(source_path, catalog_path):
        return message
    full = str(source_path)
    return message.replace(full, source_path.name or "external source file")


def _display_catalog_column(column: str) -> str:
    return "source_unit_id (or paper_id)" if column == "paper_id" else column


def _format_source_mapping_contract(mappings: list[tuple[str, str]]) -> list[str]:
    return [f"{column} -> {meaning}" for column, meaning in mappings]


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


# Two genome-wide tables in one identifier space overlap far above this; two
# DEG-only lists can legitimately overlap little, so this stays a warning.
IDENTIFIER_OVERLAP_WARNING_SHARE = 0.10
# A minority space at or above this share of a unit's identifiers is reported.
MIXED_IDENTIFIER_WARNING_SHARE = 0.10


def _spread_sample(identifiers: set[str], size: int = 200) -> list[str]:
    """Up to `size` identifiers taken evenly across the sorted set."""

    ordered = sorted(identifiers)
    if len(ordered) <= size:
        return ordered
    step = len(ordered) / size
    return [ordered[int(index * step)] for index in range(size)]


def run_warning_messages(metrics: dict[str, Any]) -> list[str]:
    messages: list[str] = []
    seen: set[str] = set()

    warning_values: list[Any] = []
    for key in ("warnings", "identifier_space_warnings", "rank_universe_warnings"):
        warning_values.extend(metrics.get(key, []) or [])
    if metrics.get("gold_panel_status") in {"invalid", "read_error"}:
        warning_values.append(
            metrics.get("gold_panel_reason") or "GoldPanel could not support curated recall metrics"
        )
    elif metrics.get("gold_panel_status") == "not_provided" and int(metrics.get("n_panel_rows", 0) or 0):
        # A filled GoldPanel whose rows all say locked=no is dropped in full. That
        # was visible only inside DEGORA_output.validation.txt, so a user who set
        # up a panel had no way to learn from the run that it was never used.
        warning_values.append(
            metrics.get("gold_panel_reason")
            or "GoldPanel rows were found but none are locked; curated recall was not calculated"
        )

    # Rows dropped for a missing or unparsable gene, effect or p-value used to be
    # reported only past a 10% share. Eight '<1E-16' rows in a 300-row table are
    # under that share and are the eight most significant genes in the table, so
    # every non-zero count is named here, compactly, whatever its share.
    dropped_counts = metrics.get("unusable_row_counts") or {}
    input_counts = metrics.get("input_row_counts") or {}
    already_detailed = " ".join(str(value) for value in warning_values)
    if isinstance(dropped_counts, dict):
        compact: list[str] = []
        for study_id, count in dropped_counts.items():
            try:
                dropped = int(count)
            except (TypeError, ValueError):
                continue
            if not dropped or f"{study_id}: {dropped:,} of" in already_detailed:
                continue
            total = input_counts.get(study_id) if isinstance(input_counts, dict) else None
            of_total = f" of {int(total):,}" if isinstance(total, (int, float)) and total else ""
            compact.append(f"{study_id}: {dropped:,}{of_total} row(s)")
        if compact:
            warning_values.append(
                "Rows without a usable gene identifier, numeric log2 fold change or numeric p-value were "
                "dropped before ranking - " + "; ".join(compact) + ". See unusable_row_counts in slice_metrics.json."
            )

    try:
        clipped_rows = int(metrics.get("pvalue_clipped_rows", 0) or 0)
    except (TypeError, ValueError):
        clipped_rows = 0
    if clipped_rows:
        text = (
            f"{clipped_rows} row(s) reported pvalue < 1e-300; values were floored "
            "to 1e-300 before signed-z scoring."
        )
        warning_values.append(text)

    for value in warning_values:
        text = str(value).strip()
        if text and text not in seen:
            messages.append(text)
            seen.add(text)

    return messages


def _identifier_space_warnings(harmonized: pd.DataFrame, *, min_studies: int = 2) -> list[str]:
    """Flag source units whose identifiers do not meet the rest of the corpus.

    Mixing Ensembl IDs with gene symbols produces a run that succeeds and reports
    every source unit as independent, while the mismatched unit supports no gene
    at all. Nothing else in the pipeline notices, because a unit that never joins
    simply never appears in a consensus row.
    """

    if min_studies <= 1 or harmonized.empty or "source_unit_id" not in harmonized.columns:
        return []
    by_unit = {
        str(unit): set(group["gene_symbol"].dropna().astype(str))
        for unit, group in harmonized.groupby("source_unit_id", sort=True)
    }
    if len(by_unit) < 2:
        return []

    warnings: list[str] = []
    # The space each unit is written in, judged from its values. Overlap alone
    # missed a real case: a column named "Gene Symbol" that held Ensembl IDs
    # overlapped a symbol table on the 4% of rows where both fell back to
    # Ensembl, and the run reported success over 869 genes that were all ENSG.
    # Sampled across the sorted range, not from its head: the first two hundred
    # sorted identifiers of a symbol table with a few ENSG fallbacks are all ENSG.
    spaces = {unit: identifier_space(_spread_sample(identifiers)) for unit, identifiers in by_unit.items() if identifiers}
    # A unit can be written in two spaces at once - 70% symbols, 30% Ensembl -
    # and the majority rule above calls it symbols. The minority then joins
    # nothing in a symbol corpus, silently. Say the split; do not convert or refuse.
    for unit, identifiers in by_unit.items():
        profile = identifier_space_profile(_spread_sample(identifiers))
        minority = sorted(
            ((share, label) for label, share in profile.items() if label != spaces.get(unit) and share >= MIXED_IDENTIFIER_WARNING_SHARE),
            reverse=True,
        )
        if minority:
            share, label = minority[0]
            warnings.append(
                f"source_unit_id={unit!r} mixes identifier spaces: about {share:.0%} of its gene identifiers "
                f"are {label} while the rest are {spaces.get(unit, 'unrecognised')}. The minority joins only "
                "sources written the same way; map the unit onto one space before trusting its contribution."
            )
    named = [space for space in spaces.values() if space != UNKNOWN_IDENTIFIER_SPACE]
    dominant = max(set(named), key=named.count) if named else UNKNOWN_IDENTIFIER_SPACE
    for unit, identifiers in by_unit.items():
        if not identifiers:
            continue
        space = spaces.get(unit, UNKNOWN_IDENTIFIER_SPACE)
        if space != UNKNOWN_IDENTIFIER_SPACE and dominant != UNKNOWN_IDENTIFIER_SPACE and space != dominant:
            warnings.append(
                f"source_unit_id={unit!r} is written in {space} while the rest of the corpus is written in "
                f"{dominant}; the genes it shares with the others are only those both sides happened to "
                "write the same way. Map every source onto one identifier space before trusting the ranking."
            )
            continue
        best_unit, best_overlap = "", 0
        for other_unit, other in by_unit.items():
            if other_unit == unit:
                continue
            overlap = len(identifiers & other)
            if overlap > best_overlap:
                best_unit, best_overlap = other_unit, overlap
        share = best_overlap / len(identifiers)
        if best_overlap == 0:
            example = sorted(identifiers)[0]
            warnings.append(
                f"source_unit_id={unit!r} shares no gene identifier with any other source unit "
                f"({len(identifiers):,} identifiers, e.g. {example!r}); it cannot contribute to any "
                "score. Map every source onto one identifier space (all symbols, or all Ensembl IDs)."
            )
        elif share < IDENTIFIER_OVERLAP_WARNING_SHARE:
            warnings.append(
                f"source_unit_id={unit!r} shares only {best_overlap:,} of its {len(identifiers):,} "
                f"gene identifiers with any other source unit (best match {best_unit!r}); check that "
                "every source uses the same identifier space."
            )
    return warnings


def _count_labels(series: pd.Series, *, unknown_label: str = "unknown") -> dict[str, int]:
    """Return JSON-safe value counts with missing or blank labels collapsed."""

    labels = series.map(lambda value: _nonempty(value) or unknown_label)
    return {str(label): int(count) for label, count in labels.value_counts(dropna=False).items()}


def _readable_read_failure(path: Path, exc: Exception) -> str:
    """Describe why a config file could not be opened, in the reader's terms.

    A file that is a valid ZIP but not a workbook reaches pandas as an unknown
    spreadsheet engine, and the only line describing the problem was an internal
    option key. Name the shape of the file instead; the original message is still
    available in the traceback the error chains from.
    """

    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        # A missing reader and a damaged file are different problems with different
        # fixes, and reporting the first as the second sent readers looking for a
        # corrupt download when a valid workbook simply had no engine to open it.
        message = str(exc).lower()
        if "xlrd" in message or "openpyxl" in message or "install" in message:
            engine = "xlrd" if suffix == ".xls" else "openpyxl"
            return (
                f"{path.name} needs the {engine} reader, which is not installed in this environment. "
                f"Reinstall DEGORA (pip install -e . or pip install degora), or install {engine} directly. "
                "The file itself was not read, so this says nothing about whether it is valid."
            )
        if suffix == ".xls":
            # Legacy .xls is an OLE2 compound file, not a ZIP.
            header = b""
            try:
                with path.open("rb") as handle:
                    header = handle.read(8)
            except OSError:
                pass
            if header.startswith(b"\xd0\xcf\x11\xe0"):
                return (
                    f"{path.name} is a legacy Excel workbook, but its contents could not be read: {exc}. "
                    "Open it in a spreadsheet tool and save it as .xlsx, or export the sheet as CSV."
                )
            return (
                f"{path.name} is not a legacy Excel workbook: it does not start with the OLE2 signature "
                "such a file has. It may be a CSV or TSV that was renamed, or a partial download."
            )
        try:
            with zipfile.ZipFile(path) as archive:
                names = set(archive.namelist())
        except (OSError, zipfile.BadZipFile):
            return (
                f"{path.name} is not a readable Excel workbook: it is not a ZIP archive. "
                "It may be a CSV or TSV that was renamed, or a partial download."
            )
        if "[Content_Types].xml" not in names:
            sample = ", ".join(sorted(names)[:3]) or "(empty archive)"
            return (
                f"{path.name} is a ZIP archive but not an Excel workbook: it has no [Content_Types].xml "
                f"(it contains {sample}). Save the config from Excel or a spreadsheet tool as .xlsx, "
                "or use a CSV config instead."
            )
    encoding_hint = _text_encoding_hint(path, exc)
    if encoding_hint:
        return encoding_hint
    return str(exc)


def _text_encoding_hint(path: Path, exc: Exception) -> str:
    """Name the encoding problem behind a raw codec error, when it can be seen.

    A UTF-16 export (Excel's "Unicode Text") starts with a byte-order mark and
    used to surface as "'utf-8' codec can't decode byte 0xff", which says nothing
    a reader can act on.
    """

    if not isinstance(exc, UnicodeDecodeError) and "codec can't decode" not in str(exc):
        return ""
    header = b""
    try:
        with path.open("rb") as handle:
            header = handle.read(4)
    except OSError:
        return ""
    if header.startswith((b"\xff\xfe", b"\xfe\xff")):
        return (
            f"{path.name} is UTF-16 encoded (it starts with a UTF-16 byte-order mark), and DEGORA reads "
            "text tables as UTF-8. Open it in a spreadsheet tool and save it as CSV UTF-8, or convert it "
            "with `iconv -f UTF-16 -t UTF-8`."
        )
    return (
        f"{path.name} is not UTF-8 encoded text ({exc}). Save it as CSV UTF-8 from a spreadsheet tool, "
        "or convert it with iconv, and try again."
    )


def _read_catalog_frame(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise DegoraConfigError(
            "config file was not found",
            context=f"config file: {path}",
            problems=[f"No file exists at {path}."],
            fixes=[
                "Check the path you passed to `degora validate` / `degora run`.",
                "Create a starter config with `degora template <name>.xlsx`, or run `degora demo` for a worked example.",
            ],
        )
    suffix = path.suffix.lower()
    try:
        if suffix in {".xlsx", ".xls"}:
            with pd.ExcelFile(path) as workbook:
                sheet_names = workbook.sheet_names
                if "Contrasts" in sheet_names:
                    return read_config_sheet(workbook, "Contrasts")
                if len(sheet_names) == 1:
                    return read_config_sheet(workbook, sheet_names[0])
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
        frame = pd.read_csv(path)
        try:
            return restore_formula_text_if_marked(frame, path)
        except ValueError as exc:
            raise DegoraConfigError(
                "CSV config formula-guard provenance is invalid or incomplete",
                context=f"config file: {path}",
                problems=[str(exc)],
                fixes=[
                    "Keep the generated .provenance.json sidecar beside the CSV config.",
                    "If the config was edited, regenerate it or use the original raw identifiers rather than guarded spreadsheet text.",
                ],
            ) from exc
    except DegoraConfigError:
        raise
    except Exception as exc:
        raise DegoraConfigError(
            "config file could not be read",
            context=f"config file: {path}",
            problems=[_readable_read_failure(path, exc)],
            fixes=[
                "Make sure the file is a valid CSV or Excel (.xlsx) workbook.",
                "For Excel configs, keep the DEG rows on a sheet named 'Contrasts'.",
                "For TSV-style files, check the delimiter and text encoding.",
            ],
        ) from exc


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
    try:
        with pd.ExcelFile(path) as workbook:
            if "GoldPanel" not in workbook.sheet_names:
                return {
                    "status": "not_provided",
                    "source": "GoldPanel",
                    "gene_column": "",
                    "genes": [],
                    "reason": "Excel workbook has no GoldPanel sheet",
                }
            gold = read_config_sheet(workbook, "GoldPanel")
    except (ImportError, KeyError, OSError, TypeError, ValueError, zipfile.BadZipFile) as exc:
        # Parser messages can contain local paths or workbook cell text. Record only
        # the exception type so persisted metrics remain safe to share.
        return {
            "status": "read_error",
            "source": "GoldPanel",
            "gene_column": "gene_symbol",
            "genes": [],
            "reason": (
                f"GoldPanel could not be read ({type(exc).__name__}); "
                "curated recall was not calculated"
            ),
        }
    if "gene_symbol" not in gold.columns:
        return {
            "status": "invalid",
            "source": "GoldPanel",
            "gene_column": "gene_symbol",
            "genes": [],
            "reason": (
                "GoldPanel is missing the required gene_symbol column; "
                "curated recall was not calculated"
            ),
        }
    gold = gold.copy()
    named_rows = int(gold["gene_symbol"].astype("string").fillna("").str.strip().ne("").sum())
    if "locked" in gold.columns:
        gold = gold.loc[locked_panel_mask(gold["locked"])]
    # Resolve panel symbols the same way source tables are resolved. A panel written
    # with legacy symbols (SEPT9, MARCH1, DEC1) previously reported every one of its
    # genes as absent, because the scored output holds the current symbol instead.
    genes = sorted(
        {
            symbol
            for symbol in (
                canonical_gene_symbol(value) for value in gold["gene_symbol"].tolist()
            )
            if symbol
        }
    )
    if genes:
        reason = ""
    elif named_rows:
        reason = (
            f"GoldPanel lists {named_rows:,} gene(s) but none are locked; set locked=yes "
            "on the rows that should be used for curated recall"
        )
    else:
        reason = "GoldPanel contains no locked gene symbols"
    return {
        "status": "locked" if genes else "not_provided",
        "source": "GoldPanel",
        "gene_column": "gene_symbol",
        "genes": genes,
        "reason": reason,
        "n_panel_rows": named_rows,
    }


def _normalize_header_label(label: Any) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(label).strip().lower()).strip("_")


def _canonicalize_catalog_headers(catalog: pd.DataFrame) -> pd.DataFrame:
    """Map catalog column headers to their canonical names ignoring case/spaces/hyphens.

    Only headers that resolve to a known catalog column or alias are renamed (e.g.
    ``Gene Column``/``GENE_COLUMN``/``P-Column`` -> ``gene_column``/``p_column``). This is
    a header convenience only; the *values* that map to source-table columns
    (gene_column, lfc_column, ...) are still matched against the DEG file exactly. A
    header is left untouched if its canonical form is already present, so an exact
    column always wins over a differently-cased duplicate.
    """

    known = set(CATALOG_COLUMNS) | set(CATALOG_ALIASES.keys())
    claimed = set(catalog.columns)
    rename: dict[str, str] = {}
    for column in catalog.columns:
        canonical = _normalize_header_label(column)
        if column == canonical or canonical not in known or canonical in claimed:
            continue
        rename[column] = canonical
        claimed.add(canonical)
    return catalog.rename(columns=rename) if rename else catalog


def _promoted_alias_warnings(catalog: pd.DataFrame) -> list[str]:
    """Report a legacy header that decided a setting the canonical column left blank.

    Promotion is silent by design for headers that only rename a value. It is not
    safe to be silent for time_course_mode: a blank canonical column means `mean`,
    so promoting a legacy `temporal_mode` of `early` moves a run off the default
    and can drop every gene a source unit did not measure at its earliest point.
    A reader told only that "configs that explicitly use early or late can change"
    reads their blank column and concludes the change does not reach them.
    """

    warnings: list[str] = []
    for alias, canonical in CATALOG_ALIASES.items():
        if canonical != "time_course_mode" or alias not in catalog.columns:
            continue
        promoted_rows = (
            pd.Series(True, index=catalog.index)
            if canonical not in catalog.columns
            else catalog[canonical].map(_nonempty).isna()
        )
        promoted = sorted(
            {
                str(value).strip()
                for value in catalog.loc[promoted_rows, alias].tolist()
                if _nonempty(value) and _normalize_time_course_setting(value) not in {None, "mean"}
            }
        )
        if promoted:
            warnings.append(
                f"The legacy '{alias}' column set time_course_mode to {', '.join(promoted)} because "
                "the 'time_course_mode' column is absent or blank. A blank column means 'mean', so "
                "this run is not using the default: early and late keep only the source unit's "
                "earliest or latest timed contrast, and genes measured at no other time drop out. "
                f"Put the value in 'time_course_mode' to state it, or clear '{alias}' to run as mean."
            )
    return warnings


def _normalize_catalog_columns(catalog: pd.DataFrame) -> pd.DataFrame:
    frame = _canonicalize_catalog_headers(catalog.copy())
    # Computed before promotion fills the canonical column, and carried on the
    # frame so validate and run can both report it without re-reading the file.
    promoted_warnings = _promoted_alias_warnings(frame)
    for alias, canonical in CATALOG_ALIASES.items():
        if alias not in frame.columns:
            continue
        if canonical not in frame.columns:
            frame[canonical] = frame[alias]
            continue
        # An explicit canonical header wins only where it has a value. This keeps
        # the documented source_unit_id alias usable when a beginner template also
        # contains a blank legacy paper_id column.
        blank_canonical = frame[canonical].map(_nonempty).isna()
        # A column that is blank in every row arrives typed as float64, because
        # that is what a CSV or workbook of empty cells reads as. Writing the
        # alias's text into it is an upcast, which pandas 3 refuses outright
        # instead of widening - so the column is widened here, deliberately,
        # before the values land. Without this a blank paper_id column beside a
        # filled source_unit_id ends a run in a raw TypeError.
        if blank_canonical.any():
            frame[canonical] = frame[canonical].astype("object")
        frame.loc[blank_canonical, canonical] = frame.loc[blank_canonical, alias]
    for column, default in OPTIONAL_CATALOG_DEFAULTS.items():
        if column not in frame.columns:
            frame[column] = default
    frame.attrs["promoted_alias_warnings"] = promoted_warnings
    return frame


def _validate_catalog_columns(
    catalog: pd.DataFrame,
    path: Path,
    *,
    file_columns: list[Any] | None = None,
) -> None:
    missing = [column for column in ESSENTIAL_CATALOG_COLUMNS if column not in catalog.columns]
    if not missing:
        return
    problems = [
        f"Missing required Contrasts column '{_display_catalog_column(column)}': "
        f"{CATALOG_COLUMN_HELP.get(column, 'required catalog field')}."
        for column in missing
    ]
    # Report the headers the file actually carried. Optional catalog defaults are
    # injected before this check, so listing catalog.columns told the reader that
    # source_unit_id was "available" in the same message that said it was missing.
    available = list(catalog.columns) if file_columns is None else list(file_columns)
    fixes = [
        f"Required Contrasts columns are: {_format_columns(BEGINNER_REQUIRED_CONTRAST_COLUMNS)}.",
        "Add the missing column(s) to the CSV file or to the Excel Contrasts sheet.",
        "For time-course data, put one row per time point and reuse paper_id/source_unit_id for related rows.",
        "Supported beginner aliases include source_unit_id->paper_id, time_h->duration_h, condition->hypoxia_modality, and include->include_in_analysis.",
    ]
    if _looks_like_deg_table(available):
        # The commonest way to reach this error is passing a DEG table where the
        # config belongs, and the required-column list alone does not say so.
        fixes.insert(
            0,
            "This file looks like a DEG results table, not a DEGORA config. A DEG table belongs in "
            "source_path inside the config; create the config with `degora template <name>.xlsx` or "
            "`degora init <name>.csv --deg-dir <folder>`.",
        )
    if len(available) == 1 and any(delimiter in str(available[0]) for delimiter in (";", "\t", "|")):
        # A European Excel "Save as CSV" writes semicolons; the whole header then
        # arrived as one column name and the message said study_id was missing.
        delimiter = next(d for d in (";", "\t", "|") if d in str(available[0]))
        label = {";": "semicolon", "\t": "tab", "|": "pipe"}[delimiter]
        fixes.insert(
            0,
            f"The file has a single column whose name contains the {label} character, so it is "
            f"{label}-delimited rather than comma-delimited. Save the config as a comma-separated CSV "
            "(or as .xlsx) and try again.",
        )
    raise DegoraConfigError(
        "catalog is missing required column(s)",
        context=f"config file: {path}; available columns: {_format_columns(available)}",
        problems=problems,
        fixes=fixes,
    )


DEG_TABLE_HEADER_HINTS = frozenset(
    {
        "adj.p.val",
        "adjpval",
        "avgexpr",
        "basemean",
        "fdr",
        "gene",
        "gene_id",
        "gene_symbol",
        "log2foldchange",
        "logcpm",
        "logfc",
        "p.value",
        "padj",
        "pvalue",
    }
)


def _looks_like_deg_table(columns: list[Any]) -> bool:
    """Report whether these headers look like a DEG results table."""

    seen = {str(column).strip().lower() for column in columns}
    return len(seen & DEG_TABLE_HEADER_HINTS) >= 2


def _validate_catalog_required_values(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    problems: list[str] = []
    active = catalog.loc[include_mask]
    for index, row in active.iterrows():
        for column in ESSENTIAL_CATALOG_COLUMNS:
            if _nonempty(row.get(column)) is None:
                problems.append(
                    f"Row {_user_row_number(index)} has an empty '{_display_catalog_column(column)}' value "
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
            f"Required Contrasts columns are: {_format_columns(BEGINNER_REQUIRED_CONTRAST_COLUMNS)}.",
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

        try:
            normalize_lfc_scale(row.get("lfc_scale", ""))
        except ValueError as exc:
            problems.append(
                f"Row {_user_row_number(index)} has unsupported lfc_scale={row.get('lfc_scale')!r}: {exc}"
            )
            fixes.append(
                "Write lfc_scale=log2 to state that lfc_column already holds log2 fold changes, "
                "or leave the cell blank. DEGORA never converts a scale."
            )

        time_course_mode = row.get("time_course_mode", "mean")
        if _normalize_time_course_setting(time_course_mode) is None:
            problems.append(
                f"Row {_user_row_number(index)} has time_course_mode={time_course_mode!r}; "
                "it must be mean, early, late, peak_mean, or blank."
            )
            fixes.append(
                "Use time_course_mode=mean unless you predeclare that the source unit should use early, late, or peak_mean time points."
            )

        try:
            normalize_header_row(row.get("header_row"))
        except ValueError as exc:
            problems.append(f"Row {_user_row_number(index)}: {exc}")
            fixes.append("Use header_row only when the column names are not on the first line of the file.")

        raw_universe = row.get("rank_universe_size", "")
        if _nonempty(raw_universe) is not None:
            universe = pd.to_numeric(pd.Series([raw_universe]), errors="coerce").iloc[0]
            invalid = pd.isna(universe)
            if not invalid:
                value = float(universe)
                invalid = value <= 0 or value in (float("inf"), float("-inf")) or value != int(value)
            if invalid:
                problems.append(
                    f"Row {_user_row_number(index)} has rank_universe_size={raw_universe!r}; "
                    "it must be a positive whole number of genes, or blank."
                )
                fixes.append(
                    "For DEG-only lists, enter the number of genes originally tested if the paper reports it; otherwise leave blank."
                )

    if problems:
        raise DegoraConfigError(
            "optional contrast settings are not valid (table_scope, time_course_mode, rank_universe_size)",
            context=f"config file: {path}",
            problems=problems,
            fixes=fixes,
        )


def _validate_source_unit_time_course_modes(
    catalog: pd.DataFrame,
    include_mask: pd.Series,
    path: Path,
) -> None:
    """Require one temporal preselection policy per independent source unit."""

    active = catalog.loc[include_mask].copy()
    if active.empty:
        return
    active["_source_unit_id"] = _source_unit_series(active)
    active["_normalized_time_course_mode"] = active["time_course_mode"].map(
        _normalize_time_course_setting
    )
    problems: list[str] = []
    for source_unit_id, group in active.groupby("_source_unit_id", sort=True):
        modes = sorted(set(group["_normalized_time_course_mode"].dropna().astype(str)))
        if len(modes) <= 1:
            continue
        rows = ", ".join(str(_user_row_number(index)) for index in group.index)
        blank_note = (
            " A blank time_course_mode cell means 'mean', so leaving one row blank beside "
            "another row's value is a conflict, not an omission."
            if "mean" in modes and group["time_course_mode"].map(_nonempty).isna().any()
            else ""
        )
        problems.append(
            f"source_unit_id={source_unit_id!r} uses conflicting normalized modes {modes} "
            f"across config rows {rows}.{blank_note}"
        )
    if problems:
        raise DegoraConfigError(
            "source unit has conflicting time_course_mode values",
            context=f"config file: {path}",
            problems=problems,
            fixes=[
                "Use one time_course_mode for every active row sharing a source_unit_id.",
                "Use mean to keep every related contrast, early/late to keep the globally earliest/latest "
                "duration_h, or peak_mean for the documented gene-specific strongest-window summary.",
            ],
        )


def _validate_time_course_durations(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    """Require a plain numeric duration_h on every row an early/late unit selects from.

    `early` and `late` keep the contrast at the unit's smallest or largest
    duration_h. A duration written as "30min" used to parse as 30, "4h" as 4, and
    a blank one was silently discarded once any sibling had a number - so the
    unit's earliest point could be the wrong contrast, with every direction that
    changed between the two time points inverted and nothing said.
    """

    active = catalog.loc[include_mask].copy()
    if active.empty:
        return
    active["_source_unit_id"] = _source_unit_series(active)
    active["_mode"] = active["time_course_mode"].map(_normalize_time_course_setting)
    problems: list[str] = []
    for source_unit_id, group in active.groupby("_source_unit_id", sort=True):
        modes = set(group["_mode"].dropna().astype(str))
        if not modes & {"early", "late"}:
            continue
        for index, row in group.iterrows():
            raw = row.get("duration_h")
            if not np.isfinite(duration_hours(raw)):
                shown = "" if _nonempty(raw) is None else str(raw).strip()
                problems.append(
                    f"Row {_user_row_number(index)} (source_unit_id={source_unit_id!r}, time_course_mode="
                    f"{sorted(modes)[0]}) has duration_h={shown!r}; early/late need a plain number of hours "
                    "on every row of the unit."
                )
    if problems:
        raise DegoraConfigError(
            "time-course durations are not numeric",
            context=f"config file: {path}",
            problems=problems,
            fixes=[
                "Write duration_h as a number of hours (0.5, 4, 24), without units or text.",
                "Fill duration_h on every active row that shares the source_unit_id, or switch that unit to "
                "time_course_mode=mean, which does not select by duration.",
            ],
        )


def _ignored_workbook_settings(row: dict[str, Any], source_path: Path) -> list[str]:
    """Say when a filled sheet_name applies to nothing.

    A sheet_name on a CSV row was silently ignored. Silence hides the two
    mistakes it usually means - the reader intended a workbook, or pointed the
    row at the wrong file - so it is reported, and the run goes on.
    """

    sheet = _nonempty(row.get("sheet_name"))
    if sheet is None or is_workbook_path(source_path):
        return []
    return [
        f"{row['study_id']}: sheet_name={str(sheet).strip()!r} is set but {source_path.name} is not a workbook, so "
        "it is ignored. If the DEG table lives in an Excel sheet, point source_path at the .xlsx/.xls file; "
        "otherwise clear sheet_name."
    ]


def _source_table_identity(row: dict[str, Any], catalog_path: Path) -> tuple[str, str, str]:
    """(resolved path, sheet, header row) - what a source table is, before column mapping."""

    source_path = _resolve_source_path(row["source_path"], catalog_path)
    try:
        resolved = str(source_path.resolve())
    except OSError:
        resolved = str(source_path)
    # A sheet name means nothing for a CSV; carrying it into the identity let
    # the same CSV declared twice, once with a stray sheet_name, pass the
    # duplicate-table check as two different tables.
    sheet = (_nonempty(row.get("sheet_name")) or "") if is_workbook_path(source_path) else ""
    # header_row blank and header_row=1 read the same file the same way; as raw
    # text they made two identities of one table.
    header_row = normalize_header_row(row.get("header_row"))
    header = "" if header_row in (None, 1) else str(header_row)
    return resolved, sheet, header


def _source_mapping_identity(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return tuple(  # type: ignore[return-value]
        (_nonempty(row.get(column)) or "") for column in ("gene_column", "lfc_column", "p_column", "padj_column")
    )


def _file_digest(path: Path) -> str:
    # Only a regular file is hashed: opening a FIFO here would wait for a writer
    # that never comes, which is the hang _require_readable_source_file exists to
    # refuse a few lines later.
    try:
        if not stat.S_ISREG(os.stat(path).st_mode):
            return ""
    except OSError:
        return ""
    digest = hashlib.sha256()
    try:
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError:
        return ""
    return digest.hexdigest()


def _reject_duplicate_source_tables(catalog: pd.DataFrame, catalog_path: Path) -> list[str]:
    """Refuse one result table counted as two independent source units.

    DEGORA's replication rule counts independent source units, and nothing in the
    scoring can tell that two units were the same file: the same table declared
    under U1 and U2 validated, ran, and gave every gene "2 / 2 source units" with
    perfect concordance. The same table mapped twice inside one unit (two
    contrasts read from two column sets of one workbook) is ordinary and is left
    alone. Returns advisory warnings for the borderline case - one file, two
    units, but different column mappings - which can be a workbook that really
    does carry two studies side by side.
    """

    rows = catalog.to_dict(orient="records")
    units = _source_unit_series(catalog).tolist()
    by_table: dict[tuple[str, str, str], list[tuple[Any, ...]]] = {}
    digests: dict[str, str] = {}
    by_digest: dict[tuple[str, str, str], list[tuple[Any, ...]]] = {}
    for index, (row, unit) in enumerate(zip(rows, units)):
        identity = _source_table_identity(row, catalog_path)
        mapping = _source_mapping_identity(row)
        study_id = str(row.get("study_id", ""))
        by_table.setdefault(identity, []).append((study_id, unit, mapping, catalog.index[index]))
        path = identity[0]
        if path not in digests:
            digests[path] = _file_digest(Path(path))
        if digests[path]:
            # Reported as the catalog wrote it: an absolute path outside the config
            # folder is private, and the existing readers hide it for that reason.
            by_digest.setdefault((digests[path], identity[1], identity[2]), []).append(
                (study_id, unit, mapping, catalog.index[index], str(row.get("source_path", "")), path)
            )

    problems: list[str] = []
    warnings: list[str] = []
    seen_pairs: set[tuple[str, str]] = set()

    def _report(entries: list[tuple], *, same_bytes_only: bool) -> None:
        for left_index, left in enumerate(entries):
            for right in entries[left_index + 1 :]:
                left_study, left_unit, left_mapping, left_row = left[:4]
                right_study, right_unit, right_mapping, right_row = right[:4]
                if left_unit == right_unit:
                    continue
                pair = tuple(sorted((left_study, right_study)))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                if same_bytes_only:
                    what = (
                        f"{left_study} ({left[4]}) and {right_study} ({right[4]}) are byte-identical files "
                        f"declared under different source units ({left_unit!r} and {right_unit!r})"
                    )
                else:
                    what = (
                        f"{left_study} and {right_study} read the same table (rows {_user_row_number(left_row)} and "
                        f"{_user_row_number(right_row)}) under different source units "
                        f"({left_unit!r} and {right_unit!r})"
                    )
                if left_mapping == right_mapping:
                    problems.append(what + " with identical column mappings, so one result table would count twice.")
                else:
                    warnings.append(
                        what + " with different column mappings. That is fine only if the file really carries two "
                        "independent studies side by side; if both column sets describe one experiment, give the "
                        "rows one source_unit_id."
                    )

    for entries in by_table.values():
        if len(entries) > 1:
            _report(entries, same_bytes_only=False)
    for entries in by_digest.values():
        if len(entries) > 1 and len({entry[5] for entry in entries}) > 1:
            _report(entries, same_bytes_only=True)
    if problems:
        raise DegoraConfigError(
            "one result table is declared as two independent source units",
            context=f"config file: {catalog_path}",
            problems=problems,
            fixes=[
                "Give the rows that come from one experiment the same source_unit_id; DEGORA counts independent "
                "source units, and one table cannot replicate itself.",
                "If the second row was meant to be a different file, fix its source_path.",
            ],
        )
    return warnings


def _validate_optional_replicate_counts(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    """Reject malformed declared replicate counts without inventing missing metadata.

    Legacy summary tables may omit both counts and remain analyzable with the
    documented missing-count quality weight. Once either count is supplied,
    however, both must be positive whole numbers for the exact contrast.
    """

    problems: list[str] = []
    active = catalog.loc[include_mask]
    for index, row in active.iterrows():
        raw_control = _nonempty(row.get("n_ctrl"))
        raw_treatment = _nonempty(row.get("n_treat"))
        if raw_control is None and raw_treatment is None:
            continue
        if raw_control is None or raw_treatment is None:
            missing = "n_ctrl" if raw_control is None else "n_treat"
            problems.append(
                f"Row {_user_row_number(index)} supplies only one replicate count; {missing} is blank."
            )
            continue
        for label, raw_value in (("n_ctrl", raw_control), ("n_treat", raw_treatment)):
            numeric = pd.to_numeric(pd.Series([raw_value]), errors="coerce").iloc[0]
            invalid = pd.isna(numeric)
            if not invalid:
                value = float(numeric)
                invalid = value <= 0 or value in (float("inf"), float("-inf")) or value != int(value)
            if invalid:
                problems.append(
                    f"Row {_user_row_number(index)} has {label}={raw_value!r}; it must be a positive whole number."
                )
                continue
            # The same implausible-size bound the discovery review panel applies.
            # Without it here, a typing slip in a workbook (3 becomes 999999)
            # passed validation and saturated the contrast's sample-size weight
            # at its min(sqrt(n_ctrl + n_treat), 4) ceiling, moving almost every
            # rank with no warning anywhere in the run.
            if float(numeric) > MAX_DECLARED_GROUP_SIZE:
                problems.append(
                    f"Row {_user_row_number(index)} has {label}={raw_value!r}; that is not a plausible "
                    f"number of biological replicates for one contrast (maximum "
                    f"{MAX_DECLARED_GROUP_SIZE:,}). Enter the samples in this group, not cells or reads."
                )
    if problems:
        raise DegoraConfigError(
            "replicate counts are not valid",
            context=f"config file: {path}",
            problems=problems,
            fixes=[
                "Enter positive whole-number biological replicate counts for both n_ctrl and n_treat.",
                "Use counts for the exact selected contrast, not the total number of samples in the study.",
                "If counts are genuinely unavailable for a legacy published table, leave both cells blank.",
            ],
        )


def _reject_duplicate_active_study_ids(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    """Each active row must carry a unique study_id.

    Source-unit collapse groups by (gene_symbol, source_unit_id), so two active
    rows that share a study_id would silently double-count one contrast inside
    the within-source-unit weighted mean instead of being rejected as a
    copy-paste error. study_id is documented as a unique per-contrast identifier,
    so enforce it for active rows.
    """

    if "study_id" not in catalog.columns:
        return
    ids = catalog.loc[include_mask, "study_id"].astype("string").str.strip()
    collide = ids[ids.duplicated(keep=False) & ids.ne("")]
    if not collide.empty:
        duplicated = sorted(collide.dropna().unique())
        raise DegoraConfigError(
            "catalog has duplicate study_id values among active rows",
            context=f"config file: {path}",
            problems=[f"study_id {value!r} appears on more than one active row" for value in duplicated],
            fixes=[
                "Give each active contrast row a unique study_id; group related time points "
                "or contrasts from one source with a shared paper_id or source_unit_id (not study_id)."
            ],
        )


# Exported provenance columns - contributing_study_ids, source_units,
# contributing_source_paths - are semicolon-joined lists. An identifier holding the
# delimiter makes those lists unparseable and inflates every count derived by
# splitting them.
IDENTIFIER_LIST_DELIMITER = ";"


def _reject_delimiter_in_identifiers(catalog: pd.DataFrame, include_mask: pd.Series, path: Path) -> None:
    """Keep the list delimiter out of the identifiers those lists are built from."""

    problems: list[str] = []
    active = catalog.loc[include_mask]
    for column in ("study_id", "paper_id", "source_unit_id"):
        if column not in active.columns:
            continue
        for index, value in active[column].items():
            text = _nonempty(value)
            if text and IDENTIFIER_LIST_DELIMITER in text:
                problems.append(
                    f"Row {_user_row_number(index)} has {_display_catalog_column(column)}={text!r}, "
                    f"which contains {IDENTIFIER_LIST_DELIMITER!r}."
                )
    if not problems:
        return
    raise DegoraConfigError(
        "identifier contains the list delimiter",
        context=f"config file: {path}",
        problems=problems,
        fixes=[
            f"Remove {IDENTIFIER_LIST_DELIMITER!r} from study_id, paper_id and source_unit_id.",
            "DEGORA joins contributing study IDs and source units into semicolon-separated lists in "
            "the score table, the evidence table and the workbook; an identifier holding that "
            "character makes those lists ambiguous and the counts derived from them wrong.",
        ],
    )


def _mixed_species_warnings(catalog: pd.DataFrame) -> list[str]:
    """Say when one run mixes species, because scoring will not notice.

    The Search workflow keeps Human and Mouse in separate workspaces, but the
    scoring path is not species-specific: a hand-written config naming one human
    and one mouse source satisfies min_studies=2 and produces a pooled ranking.
    The species is recorded on every evidence row, so the mixing is visible after
    the fact - it was only never announced.
    """

    if "species" not in catalog.columns:
        return []
    labels = sorted({text for text in (_nonempty(value) for value in catalog["species"]) if text})
    if len(labels) < 2:
        return []
    return [
        "this run mixes "
        + ", ".join(repr(label) for label in labels)
        + " in one ranking. DEGORA scoring matches on gene symbol and is not species-specific, so "
        "these sources are pooled and can satisfy the min_studies replication rule between them. "
        "Run one species at a time unless cross-species pooling is what you intend."
    ]


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


def _reject_ambiguous_headers(catalog: pd.DataFrame, path: Path) -> None:
    """Reject catalogs where several headers resolve to the same catalog field.

    Header matching ignores case/spaces/separators, so two literal 'gene_column' columns
    -- or 'Gene Column' plus 'gene_column' -- both map to one field; silently keeping one
    and dropping the other would discard a user's data without any warning.
    """

    known = set(CATALOG_COLUMNS) | set(CATALOG_ALIASES.keys())
    seen: dict[str, list[str]] = {}
    for column in catalog.columns:
        canonical = _normalize_header_label(column)
        if canonical in known:
            seen.setdefault(canonical, []).append(str(column))
    ambiguous = {canon: cols for canon, cols in seen.items() if len(cols) > 1}
    if ambiguous:
        raise DegoraConfigError(
            "catalog has duplicate or ambiguous column headers",
            context=f"config file: {path}",
            problems=[
                f"Multiple headers map to '{canon}': {_format_columns(cols)}."
                for canon, cols in sorted(ambiguous.items())
            ],
            fixes=[
                "Keep exactly one column per field; remove or rename the duplicate header(s).",
                "Header matching ignores case, spaces, and separators (so 'Gene Column' and 'gene_column' collide).",
            ],
        )


def read_catalog(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    frame = _read_catalog_frame(path)
    _reject_ambiguous_headers(frame, path)
    catalog = _normalize_catalog_columns(frame)
    _validate_catalog_columns(catalog, path, file_columns=list(frame.columns))
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
        # Excel/CSV float-promotes an integer include column that has any blank cell, so a
        # user's 1/0 flags arrive as "1.0"/"0.0"; normalize integer-valued floats back.
        if text not in truthy and text not in falsy:
            try:
                number = float(text)
            except ValueError:
                number = None
            if number is not None and number.is_integer():
                text = str(int(number))
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


def _validate_source_columns(
    frame: pd.DataFrame,
    mapping: TableMapping,
    row: dict[str, Any],
    source_path: Path,
    catalog_path: Path,
) -> None:
    try:
        validate_table_mapping_roles(frame, mapping, study_id=row.get("study_id"))
    except ValueError as exc:
        raise DegoraConfigError(
            "source table column mapping reuses incompatible roles",
            context=_source_context(row, source_path, catalog_path, list(frame.columns)),
            problems=[str(exc)],
            fixes=[
                "Use different source-table columns for gene_column, lfc_column, and p_column.",
                "Do not point lfc_column at a p-value/FDR column or reuse gene identifiers as statistical values.",
                "padj_column may equal p_column only when the same unit-interval column is intentionally used as both p-value and adjusted p-value/FDR.",
            ],
        ) from exc

    requested = [
        ("required", column, getattr(mapping, column), meaning)
        for column, meaning in REQUIRED_SOURCE_TABLE_MAPPINGS
    ]
    if mapping.padj_column:
        requested.extend(
            ("optional_if_filled", column, getattr(mapping, column), meaning)
            for column, meaning in OPTIONAL_SOURCE_TABLE_MAPPINGS
            if getattr(mapping, column)
        )

    available = list(map(str, frame.columns))
    problems: list[str] = []
    fixes: list[str] = []
    for requirement, catalog_column, source_column, meaning in requested:
        # Match harmonize_frame: a catalog written before restored row labels got
        # one name may still spell that column the way pandas did.
        if resolve_column_name(frame, source_column) in frame.columns:
            continue
        suggestion = difflib.get_close_matches(str(source_column), available, n=1)
        if requirement == "required":
            problems.append(
                f"{row['study_id']}: required source-table mapping {catalog_column}={source_column!r} "
                f"should identify the column for {meaning}, but that column was not found."
            )
        else:
            problems.append(
                f"{row['study_id']}: optional {catalog_column} is filled with {source_column!r}, "
                f"so DEGORA expects a source-table column for {meaning}; that column was not found. "
                f"Leave {catalog_column} blank if the source table does not provide it."
            )
        if suggestion:
            fixes.append(f"{row['study_id']}: change {catalog_column} to '{suggestion[0]}'.")
    if problems:
        raise DegoraConfigError(
            "source table column mapping is wrong",
            context=_source_context(row, source_path, catalog_path, available),
            problems=problems,
            fixes=fixes
            + [
                "Required source-table mappings are: "
                + _format_columns(_format_source_mapping_contract(REQUIRED_SOURCE_TABLE_MAPPINGS))
                + ".",
                "Optional source-table mappings such as padj_column are checked only when you fill them.",
                "Open the DEG table and copy the column name exactly into the Contrasts sheet.",
                "Column names are case-sensitive.",
            ],
        )
    _reject_ambiguous_source_headers(frame, mapping, row, source_path, catalog_path)


def _reject_ambiguous_source_headers(
    frame: pd.DataFrame,
    mapping: TableMapping,
    row: dict[str, Any],
    source_path: Path,
    catalog_path: Path,
) -> None:
    """Refuse a mapping onto a header the source table carries more than once.

    A supplementary table with two contrasts side by side is written
    `gene, logFC, P.Value, logFC, P.Value`. pandas renames the second block to
    `logFC.1` before anyone sees it, so a catalog that says `logFC` bound to the
    first block in silence - and when the intended contrast was the second one,
    every direction in the results was inverted with nothing to notice.
    """

    duplicates = duplicate_source_headers(frame)
    if not duplicates:
        return
    problems: list[str] = []
    for catalog_column in ("gene_column", "lfc_column", "p_column", "padj_column"):
        source_column = getattr(mapping, catalog_column)
        if not source_column:
            continue
        count = duplicates.get(str(source_column))
        if count:
            mangled = ", ".join(f"'{source_column}.{index}'" for index in range(1, count))
            problems.append(
                f"{row['study_id']}: {catalog_column}={source_column!r} names a header that appears {count} times "
                f"in the table; pandas exposes the later copies as {mangled}, and DEGORA cannot tell which block "
                "you mean."
            )
    if not problems:
        return
    raise DegoraConfigError(
        "source table has repeated column headers",
        context=_source_context(row, source_path, catalog_path, list(map(str, frame.columns))),
        problems=problems,
        fixes=[
            "Rename the repeated headers in the source table so each contrast has its own column names, "
            "or split the table into one file per contrast.",
            "To use a later block deliberately, set the mapping to the pandas name of that copy "
            "(for example lfc_column=logFC.1), which selects it explicitly.",
        ],
    )


# A handful of unparsable cells is ordinary in a published table; a column that is
# largely text is a mapping mistake, and validate is where it should surface. Both
# conditions are needed: the share alone rejects a ten-row table over one odd cell,
# and the count alone rejects a twenty-thousand-row table over three.
NON_NUMERIC_REJECT_SHARE = 0.10
NON_NUMERIC_REJECT_MINIMUM = 3




def _count_normalized_gene_symbols(harmonized: pd.DataFrame) -> int:
    """Count genes whose stored symbol differs from the label the source carried.

    Repairing "9-Sep" to SEPTIN9 changes the name the reader will search for, so
    the run has to be able to say how often it happened.
    """

    if not {"input_gene_label", "gene_symbol"}.issubset(harmonized.columns):
        return 0
    renamed: set[str] = set()
    labels = harmonized["input_gene_label"].tolist()
    symbols = harmonized["gene_symbol"].tolist()
    for label, symbol in zip(labels, symbols):
        text = "" if label is None or (not isinstance(label, str) and pd.isna(label)) else str(label)
        target = "" if symbol is None or (not isinstance(symbol, str) and pd.isna(symbol)) else str(symbol)
        if not text.strip() or not target.strip():
            continue
        parts = [part.strip().upper() for part in text.split(";") if part.strip()]
        if any(part != target for part in parts):
            renamed.add(target)
    return len(renamed)


def _report_all_missing_source_files(catalog: pd.DataFrame, catalog_path: Path) -> None:
    """Report every active row whose source table is missing, not just the first.

    Leaving the Excel template's example rows in place is the commonest beginner
    mistake, and it leaves more than one broken path behind. Reporting them one at
    a time turned a single edit into a fix-and-re-run loop, so collect the whole
    set before raising. Rows that exist but are not regular files are left to the
    per-row check, which explains what the path is instead.
    """

    problems: list[str] = []
    for row in catalog.to_dict(orient="records"):
        source_path = _resolve_source_path(row["source_path"], catalog_path)
        if source_path.exists():
            continue
        problems.append(_source_path_problem(row, source_path, catalog_path, "but that file does not exist."))
    if not problems:
        return
    raise DegoraConfigError(
        "source DEG table file was not found",
        problems=problems,
        fixes=[
            "Check the file path in the Contrasts sheet.",
            "Relative paths are resolved from the Excel/CSV config folder, then the current folder.",
            "If the file was moved, update source_path rather than editing analysis outputs by hand.",
            "Delete the template's example rows, or set their include column to no, before your own rows are run.",
        ],
    )


def _require_readable_source_file(source_path: Path, row: dict[str, Any], catalog_path: Path) -> None:
    """Refuse anything that is not a regular file before a reader blocks on it.

    ``exists()`` is true for a FIFO, a device and a socket, and pandas then waits
    for a writer that never comes: the command produced no output, used no CPU and
    never returned, which is indistinguishable from a hang.
    """

    if source_path.is_file():
        return
    if not source_path.exists():
        raise DegoraConfigError(
            "source DEG table file was not found",
            problems=[
                _source_path_problem(row, source_path, catalog_path, "but that file does not exist."),
            ],
            fixes=[
                "Check the file path in the Contrasts sheet.",
                "Relative paths are resolved from the Excel/CSV config folder, then the current folder.",
                "If the file was moved, update source_path rather than editing analysis outputs by hand.",
            ],
        )
    kind = "directory" if source_path.is_dir() else "not a regular file (for example a pipe, socket or device)"
    raise DegoraConfigError(
        "source DEG table path is not a readable file",
        problems=[_source_path_problem(row, source_path, catalog_path, f"which is {kind}.")],
        fixes=[
            "Point source_path at a CSV, TSV, TXT, XLS or XLSX file on disk.",
            "DEGORA reads each source table more than once, so a stream or pipe cannot be used.",
        ],
    )


# Effect-scale sanity checks. Direction is the sign of the effect, so a linear
# fold change (2.5 = up, 0.4 = down) makes every gene "up" and nothing downstream
# can notice. A column with no negative values and values on both sides of 1 -
# below 0.5 as well as above 1 - is that signature; a log2 table would carry
# negatives. |log2FC| beyond 30 is a 10^9-fold change and is not a log2 value.
LINEAR_FOLD_CHANGE_MIN_ROWS = 10
LOG2_FOLD_CHANGE_ABSURD = 30.0
# Past this share of impossible values - and this many of them, so one outlier
# in a thirty-row table stays a warning - the column is refused, not warned about.
ABSURD_SHARE_REFUSED = 0.05
ABSURD_COUNT_REFUSED = 20


def _effect_scale_problems(
    raw: pd.Series,
    numeric: pd.Series,
    *,
    study_id: str,
    source_column: str,
    declared_log2: bool = False,
) -> tuple[list[str], list[str]]:
    """Return (problems, warnings) about whether an effect column is log2-scaled.

    ``declared_log2`` is the catalog's ``lfc_scale=log2``: the reader stating the
    scale where the header does not. It carries the same weight as a header that
    says log2, because it is the same claim made by a person instead of a file.
    """

    values = numeric[np.isfinite(numeric.to_numpy(dtype=float))]
    if values.empty:
        return [], []
    header_says_log = bool(re.search(r"log", str(source_column), re.IGNORECASE)) or declared_log2
    problems: list[str] = []
    warnings: list[str] = []
    if len(values) < LINEAR_FOLD_CHANGE_MIN_ROWS and not declared_log2:
        # The shape tests below need ten values; a smaller table is neither
        # refused nor checked, and that silence read as approval.
        warnings.append(
            f"{study_id}: lfc_column={source_column!r} has only {len(values):,} numeric value(s), fewer than the "
            f"{LINEAR_FOLD_CHANGE_MIN_ROWS} DEGORA needs to tell a log2 fold change from a linear one, so the "
            "scale is taken as declared or named. If these are linear fold changes, convert them to log2; if they "
            "are already log2 and the header does not say so, set lfc_scale=log2 on this row."
        )
    n_negative = int((values < 0).sum())
    n_below_half = int(((values > 0) & (values < 0.5)).sum())
    n_above_one = int((values > 1.0).sum())
    if (
        len(values) >= LINEAR_FOLD_CHANGE_MIN_ROWS
        and n_negative == 0
        and n_below_half > 0
        and n_above_one > 0
    ):
        message = (
            f"{study_id}: lfc_column={source_column!r} has no negative values but {n_below_half:,} value(s) "
            f"between 0 and 0.5 and {n_above_one:,} value(s) above 1 (range {float(values.min()):g} to "
            f"{float(values.max()):g}). That is the shape of a linear fold change (0.4 = 2.5-fold down), not a "
            "log2 fold change. DEGORA reads direction from the sign, so every gene in this table would be "
            "called up."
        )
        # A header that says log2 is an explicit claim; a table of up-regulated
        # genes only can look like this. Without that claim the shape decides.
        if declared_log2:
            warnings.append(message + " The catalog declares lfc_scale=log2, so this is reported, not refused.")
        elif header_says_log:
            warnings.append(message + " The column name says log2, so this is reported, not refused; check it.")
        else:
            problems.append(
                message + " If these values are already log2, say so: set lfc_scale to log2 on this row."
            )
    elif n_negative == 0 and len(values) >= LINEAR_FOLD_CHANGE_MIN_ROWS and not header_says_log:
        warnings.append(
            f"{study_id}: lfc_column={source_column!r} has no negative values and its name does not say log2. "
            "If it holds a linear fold change, convert it to log2 before running; if this is a table of "
            "up-regulated genes only, no action is needed."
        )
    elif n_negative and not header_says_log and float(values.abs().min()) >= 1.0 and len(values) >= LINEAR_FOLD_CHANGE_MIN_ROWS:
        warnings.append(
            f"{study_id}: lfc_column={source_column!r} is not named as a log2 value and every |value| is >= 1. "
            "A signed linear fold change (-2.5 = 2.5-fold down) has that shape; so does a log2 table filtered at "
            "|log2FC| >= 1. Confirm the scale - a linear value inflates every effect size DEGORA reports."
        )
    absurd = values.abs().max()
    if float(absurd) > LOG2_FOLD_CHANGE_ABSURD:
        n_absurd = int((values.abs() > LOG2_FOLD_CHANGE_ABSURD).sum())
        message = (
            f"{study_id}: lfc_column={source_column!r} has {n_absurd:,} value(s) with |log2FC| > "
            f"{LOG2_FOLD_CHANGE_ABSURD:g} (largest {float(absurd):g}). A log2 fold change that large is a "
            "10^9-fold change or more; check that the column is on a log2 scale and not a ratio, a count or a statistic."
        )
        # A few values past the line are outliers worth a look. Thousands of them
        # are a count or an intensity mapped as the effect size - a column of
        # read counts reached scoring this way with a warning and nothing else -
        # and no scale declaration can make a count a fold change.
        if n_absurd >= ABSURD_COUNT_REFUSED and n_absurd / len(values) > ABSURD_SHARE_REFUSED:
            problems.append(
                message
                + f" {n_absurd / len(values):.0%} of the values are past that line, so this column is not a log2 "
                "fold change; map the log2 fold-change column instead."
            )
        else:
            warnings.append(message)
    return problems, warnings


def _looks_like_adjusted_p(p_column: str, mapping: TableMapping) -> bool:
    """True when the p_column is the adjusted column - by name, or because it is also padj_column."""

    from .discovery import PADJ_RE

    if mapping.padj_column and str(mapping.padj_column).strip().lower() == p_column.strip().lower():
        return True
    # PADJ_RE alone, as the classifier init and the browser use: p_val_adj and
    # pvals_adj also match the nominal pattern, and requiring them not to hid
    # exactly the Seurat and scanpy headers this is for. No nominal header
    # (pvalue, P.Value, pval, p) matches PADJ_RE, so nothing is misclassified.
    return bool(PADJ_RE.search(p_column))


def _validate_numeric_source_columns(
    frame: pd.DataFrame,
    mapping: TableMapping,
    row: dict[str, Any],
    source_path: Path,
    catalog_path: Path,
) -> list[str]:
    """Reject an effect or p-value column that does not hold usable numbers.

    The p-value column was already checked for being inside [0, 1], which is a
    range check on values that parsed. A column of 'UP'/'n/a'/'#DIV/0!' parses to
    nothing at all, passed validate, and then lost its rows silently during the
    run - so the mistake that costs the most rows was the one validate could not
    see. A p-value written as a bound ('<1E-16') is rejected outright whatever its
    share: those rows are the most significant genes in the table, and dropping
    them as unparsable removed exactly the genes the table was published for.
    Returns non-fatal warnings about the effect column's scale.
    """

    problems: list[str] = []
    scale_problem_seen = False
    fixes: list[str] = []
    warnings: list[str] = []
    for catalog_column, meaning in (("lfc_column", "log2 fold change"), ("p_column", "p-value"), ("padj_column", "adjusted p-value")):
        source_column = getattr(mapping, catalog_column)
        if not source_column:
            continue
        resolved = resolve_column_name(frame, source_column)
        if resolved not in frame.columns:
            continue
        raw = frame[resolved]
        numeric = pd.to_numeric(raw, errors="coerce")
        if catalog_column == "lfc_column":
            # "inf" parses as a number; validate said "config OK" and the run then
            # set the row aside. Say it here, where the reader is looking.
            nonfinite = numeric.notna() & ~np.isfinite(numeric.to_numpy(dtype=float))
            n_nonfinite = int(nonfinite.sum())
            if n_nonfinite:
                shown = ", ".join(str(value) for value in raw.loc[nonfinite].head(3))
                warnings.append(
                    f"{row['study_id']}: lfc_column={source_column!r} holds {n_nonfinite:,} infinite value(s) "
                    f"(examples: {shown}); an infinite effect has no usable magnitude, so the run sets those rows "
                    "aside as unusable. Check how the table was exported (a division by zero or a filtered-out "
                    "group writes Inf)."
                )
        if catalog_column == "p_column" and _looks_like_adjusted_p(str(source_column), mapping):
            # degora init refuses to pick this by default and explains it; the
            # browser requires adjusted_p_as_pvalue_confirmed; a hand-written
            # catalog said nothing. Same fact, same words, on the third path.
            warnings.append(
                f"{row['study_id']}: p_column={source_column!r} is an adjusted p-value/FDR column, so every "
                "p-value DEGORA reads from this table is already adjusted. That is a usable answer when the "
                "table has no unadjusted p-value; if it has one, map p_column to it and keep the adjusted "
                "column in padj_column."
            )
        if catalog_column == "lfc_column":
            scale_problems, scale_warnings = _effect_scale_problems(
                raw,
                numeric,
                study_id=str(row["study_id"]),
                source_column=str(source_column),
                declared_log2=normalize_lfc_scale(row.get("lfc_scale")) == "log2",
            )
            problems.extend(scale_problems)
            warnings.extend(scale_warnings)
            if scale_problems:
                scale_problem_seen = True
                fixes.append(
                    f"{row['study_id']}: convert {source_column!r} to log2 (log2 of the ratio, or sign x log2|value| "
                    "for a signed linear fold change) and save the table, map a log2 column instead, or - if "
                    "the values are already log2 - set lfc_scale to log2 on this row."
                )
        unparsed = numeric.isna() & raw.notna()
        n_unparsed = int(unparsed.sum())
        n_present = int(raw.notna().sum())
        if not n_unparsed or not n_present:
            continue
        bounded = bounded_value_examples(raw.loc[unparsed]) if catalog_column != "lfc_column" else []
        if bounded:
            n_bounded = int(
                sum(1 for value in raw.loc[unparsed] if str(value).strip() and bounded_value_examples(pd.Series([value])))
            )
            problems.append(
                f"{row['study_id']}: {catalog_column}={source_column!r} has {n_bounded:,} {meaning}(s) written as a "
                f"bound (examples: {', '.join(repr(text) for text in bounded[:3])}). Those are the most "
                "significant rows in the table, and DEGORA would drop them as unreadable."
            )
            fixes.append(
                f"{row['study_id']}: replace each bound with its numeric value (for example '<1E-16' -> 1E-16) so "
                "the row keeps a conservative p-value, or clear the cell to drop the row deliberately."
            )
            continue
        share = n_unparsed / n_present
        examples = ", ".join(
            repr(text)
            for text in dict.fromkeys(str(value).strip() for value in raw.loc[unparsed] if str(value).strip())
            if text
        )
        if n_unparsed >= NON_NUMERIC_REJECT_MINIMUM and share >= NON_NUMERIC_REJECT_SHARE:
            problems.append(
                f"{row['study_id']}: {catalog_column}={source_column!r} should hold numeric {meaning}, but "
                f"{n_unparsed:,} of {n_present:,} non-empty values ({share:.1%}) are not numbers "
                f"(examples: {examples[:160]})."
            )
    if not problems:
        return warnings
    raise DegoraConfigError(
        "source table effect column is not on a log2 scale"
        if scale_problem_seen
        else "source table effect or p-value column is not numeric",
        context=_source_context(row, source_path, catalog_path),
        problems=problems,
        fixes=fixes
        + [
            "Point lfc_column and p_column at the numeric columns of the DEG table.",
            "Replace placeholder text such as 'UP', 'NA' or a spreadsheet error value with an empty cell.",
            "Rows without a numeric effect and p-value cannot be ranked and are dropped during the run.",
        ],
    )


def _resolve_source_path(raw_source_path: Any, catalog_path: Path) -> Path:
    source_path = Path(str(raw_source_path))
    if source_path.is_absolute():
        return source_path

    # Resolve against the config directory first. That makes a config and its
    # data bundle self-contained and prevents an unrelated same-named file in
    # the process working directory from being selected accidentally.
    candidates = [
        catalog_path.parent / source_path,
        Path.cwd() / source_path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return catalog_path.parent / source_path


def _validate_generated_replicate_provenance(
    catalog: pd.DataFrame,
    catalog_path: Path,
) -> None:
    """Reconcile catalog counts with count-bearing generated-table provenance."""

    problems: list[str] = []
    for row in catalog.to_dict(orient="records"):
        source_path = _resolve_source_path(row["source_path"], catalog_path)
        provenance_path = source_path.with_suffix(source_path.suffix + ".provenance.json")
        if not provenance_path.is_file():
            continue
        try:
            payload = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            problems.append(f"{row['study_id']}: cannot read generated-source provenance {provenance_path}: {exc}")
            continue
        metadata = payload.get("metadata", {})
        if not isinstance(metadata, dict):
            continue
        recorded_control = metadata.get("n_ctrl")
        recorded_treatment = metadata.get("n_treat")
        if recorded_control is None and recorded_treatment is None:
            continue
        if recorded_control is None or recorded_treatment is None:
            problems.append(f"{row['study_id']}: generated-source provenance supplies only one of n_ctrl/n_treat.")
            continue
        declared_control = _nonempty(row.get("n_ctrl"))
        declared_treatment = _nonempty(row.get("n_treat"))
        if declared_control is None or declared_treatment is None:
            problems.append(
                f"{row['study_id']}: catalog replicate counts are blank, but generated-source provenance records "
                f"n_ctrl={recorded_control}, n_treat={recorded_treatment}."
            )
            continue
        try:
            raw_counts = (recorded_control, recorded_treatment, declared_control, declared_treatment)
            numeric_counts = tuple(float(value) for value in raw_counts)
            if any(
                value <= 0 or value in (float("inf"), float("-inf")) or value != int(value)
                for value in numeric_counts
            ):
                raise ValueError("replicate counts must be positive whole numbers")
            recorded_pair = (int(numeric_counts[0]), int(numeric_counts[1]))
            declared_pair = (int(numeric_counts[2]), int(numeric_counts[3]))
        except (TypeError, ValueError, OverflowError):
            problems.append(f"{row['study_id']}: generated-source provenance contains invalid replicate counts.")
            continue
        if declared_pair != recorded_pair:
            problems.append(
                f"{row['study_id']}: catalog n_ctrl/n_treat={declared_pair[0]}/{declared_pair[1]} disagree with "
                f"generated-source provenance {recorded_pair[0]}/{recorded_pair[1]}."
            )
        effect_direction = str(metadata.get("effect_direction", "")).strip()
        if effect_direction and effect_direction != "treatment_minus_control":
            problems.append(
                f"{row['study_id']}: generated-source provenance declares effect_direction={effect_direction!r}; "
                "DEGORA matrix fallbacks require treatment_minus_control."
            )
    if problems:
        raise DegoraConfigError(
            "catalog metadata disagree with generated-source provenance",
            context=f"config file: {catalog_path}",
            problems=problems,
            fixes=[
                "Copy n_ctrl and n_treat from the derivation summary/provenance into the catalog.",
                "Do not edit replicate counts to change source weights after a DEG table has been generated.",
                "Regenerate a matrix fallback if the selected control or treatment columns were wrong.",
            ],
        )


def _try_write_parquet(frame: pd.DataFrame, path: Path) -> str | None:
    """Write an optional parquet mirror, returning a warning on engine failure."""

    try:
        frame.to_parquet(path, index=False)
    except (ImportError, ValueError) as exc:
        return (
            "Optional harmonized parquet was not written because no usable parquet "
            f"engine is available: {exc}"
        )
    return None


def _check_identifier_space_overlap(
    identifiers_by_unit: dict[str, set[str]], path: Path
) -> list[str]:
    """Warn when source units do not share a gene identifier space.

    `degora run` already diagnoses this precisely, but only after harmonizing and
    scoring: a config mixing Ensembl IDs in one study with symbols in another
    passed `validate` with "config OK" and then failed at the end of the run.
    The preflight has the same identifier sets in hand, so it says so first.

    This warns rather than refuses on purpose. The run-time message is computed
    from the harmonized table and is the better of the two - it names the largest
    pairwise overlap and a sample identifier from each unit - so replacing it
    with an earlier refusal would cost more than it gained. A run at
    `min_studies=1`, or one whose overlap is partial, is also legitimate. The
    preflight's job here is to tell the reader before the work, not to decide
    for them.
    """

    units = {unit: genes for unit, genes in identifiers_by_unit.items() if genes}
    if len(units) < 2:
        return []
    counts = ", ".join(f"{unit} ({len(genes):,})" for unit, genes in list(units.items())[:6])
    shared = set.intersection(*units.values())
    if shared:
        return []

    items = list(units.items())
    best_pair = ""
    largest = 0
    for i, (left_name, left) in enumerate(items):
        for right_name, right in items[i + 1 :]:
            overlap = len(left & right)
            if overlap > largest:
                largest = overlap
                best_pair = f"{left_name} vs {right_name}"

    samples = "; ".join(
        f"{unit}: {sorted(genes)[0]!r}" for unit, genes in items[:3] if genes
    )
    if largest == 0:
        return [
            f"no two of the {len(units)} source units share a single gene identifier. "
            f"Identifiers per source unit: {counts}. First identifier in each: {samples}. "
            "No gene can reach a replication floor of 2, so a run at min_studies=2 or above will "
            "score nothing. Map every source onto one identifier space - all gene symbols, or all "
            "Ensembl IDs - or check that gene_column points at the identifier column and not at a "
            "probe ID, row number or rank."
        ]
    return [
        f"no identifier is shared by all {len(units)} source units; the largest overlap between "
        f"any two is {largest:,} ({best_pair}). Identifiers per source unit: {counts}. Genes are "
        "matched by identifier, so sources written in different identifier spaces contribute no "
        "shared evidence."
    ]


def validate_catalog_inputs(
    catalog_path: Path,
    *,
    progress: Callable[[int, int, str], None] | None = None,
) -> dict[str, Any]:
    """Validate catalog/config rows and source-table column mappings without writing outputs."""

    catalog_path = catalog_path.resolve()
    full_catalog = read_catalog(catalog_path)
    include_mask = catalog_include_mask(full_catalog)
    _validate_catalog_required_values(full_catalog, include_mask, catalog_path)
    _validate_optional_scope_values(full_catalog, include_mask, catalog_path)
    _validate_source_unit_time_course_modes(full_catalog, include_mask, catalog_path)
    _validate_time_course_durations(full_catalog, include_mask, catalog_path)
    _validate_optional_replicate_counts(full_catalog, include_mask, catalog_path)
    _reject_duplicate_active_study_ids(full_catalog, include_mask, catalog_path)
    _reject_delimiter_in_identifiers(full_catalog, include_mask, catalog_path)
    catalog = full_catalog.loc[include_mask].copy()
    if catalog.empty:
        raise DegoraConfigError(
            "catalog contains no active contrasts",
            problems=["Every row is excluded by include_in_analysis/include, or the catalog is empty."],
            fixes=["Set include_in_analysis/include to yes for at least one contrast row."],
        )
    _validate_generated_replicate_provenance(catalog, catalog_path)
    _report_all_missing_source_files(catalog, catalog_path)
    source_warnings: list[str] = _reject_duplicate_source_tables(catalog, catalog_path)

    checked_sources: list[str] = []
    rows = catalog.to_dict(orient="records")
    # Validation already opens and harmonizes every source table, so the gene
    # identifier each one actually carries is in hand here. Collecting them per
    # source unit lets the preflight say that two sources share no identifier
    # space - the commonest way a config that validates cleanly goes on to score
    # zero genes at the end of a full run.
    identifiers_by_unit: dict[str, set[str]] = {}
    unit_ids = _source_unit_series(catalog).tolist()
    total_rows = len(rows)
    if progress is not None:
        progress(0, total_rows, "")
    for index, row in enumerate(rows, start=1):
        source_path = _resolve_source_path(row["source_path"], catalog_path)
        _require_readable_source_file(source_path, row, catalog_path)

        mapping = TableMapping(
            gene_column=row["gene_column"],
            lfc_column=row["lfc_column"],
            p_column=row["p_column"],
            padj_column=_nonempty(row.get("padj_column")),
            sep=_nonempty(row.get("sep")),
            sheet_name=_nonempty(row.get("sheet_name")),
            header_row=normalize_header_row(row.get("header_row")),
        )
        try:
            raw_frame = read_deg_table(source_path, mapping)
        except Exception as exc:
            raise DegoraConfigError(
                "source DEG table could not be read",
                context=_source_path_context(row, source_path, catalog_path),
                problems=[_readable_source_read_failure(source_path, catalog_path, exc)],
                fixes=[
                    "Check that the file is a supported CSV/TSV/TXT/XLS/XLSX table.",
                    "For Excel sources, set sheet_name to the exact sheet containing the DEG table.",
                    "For TSV files, set sep to \\t if DEGORA cannot infer it.",
                ],
            ) from exc
        _validate_source_columns(raw_frame, mapping, row, source_path, catalog_path)
        source_warnings.extend(_validate_numeric_source_columns(raw_frame, mapping, row, source_path, catalog_path))
        source_warnings.extend(_ignored_workbook_settings(row, source_path))
        filtered_frame, _ = apply_gene_type_filter(
            raw_frame,
            _nonempty(row.get("gene_type_column")),
            _nonempty(row.get("gene_type_keep")),
        )
        row["source_path"] = portable_path(source_path, catalog_path.parent)
        try:
            harmonized = harmonize_frame(filtered_frame, mapping, row)
        except (KeyError, TypeError, ValueError) as exc:
            raise DegoraConfigError(
                "source DEG table statistical values are invalid",
                context=_source_path_context(row, source_path, catalog_path),
                problems=[str(exc)],
                fixes=[
                    "Check that gene_column, lfc_column, p_column, and padj_column (if filled) map to the intended source columns.",
                    "p_column and padj_column must contain values in [0, 1]; do not map -log10 values, test statistics, or percentages.",
                    "lfc_column must contain numeric log2 fold changes and gene_column must contain usable identifiers.",
                ],
            ) from exc
        if harmonized.empty:
            raise DegoraConfigError(
                "source DEG table produced no usable rows",
                context=_source_path_context(row, source_path, catalog_path),
                problems=[
                    "No row retained a nonblank gene identifier together with numeric log2 fold change and p-value values."
                ],
                fixes=[
                    "Verify the selected source-table columns and any gene-type filter.",
                    "Use the full results sheet rather than a notes or metadata sheet.",
                ],
            )
        unit_id = str(unit_ids[index - 1]) if index - 1 < len(unit_ids) else ""
        if unit_id and "gene_symbol" in harmonized.columns:
            identifiers_by_unit.setdefault(unit_id, set()).update(
                harmonized["gene_symbol"].dropna().astype(str)
            )
        checked_sources.append(str(source_path))
        if progress is not None:
            progress(index, total_rows, str(row.get("study_id") or source_path.name))

    # Count independent source units with the same precedence the scoring layer uses
    # (explicit source_unit_id > paper_id > study_id), so the preflight count cannot
    # disagree with aggregate.py / score_db.py.
    unit_series = _source_unit_series(catalog)
    source_units = set(unit_series[unit_series.ne("")].tolist())
    identifier_overlap_warnings = _check_identifier_space_overlap(identifiers_by_unit, catalog_path)
    gold_panel = _read_locked_gold_panel(catalog_path)
    warnings = [
        *catalog.attrs.get("promoted_alias_warnings", []),
        *_microarray_warnings(catalog),
        *_mixed_species_warnings(catalog),
        *identifier_overlap_warnings,
        *source_warnings,
    ]
    if gold_panel["status"] in {"invalid", "read_error"}:
        warnings.append(gold_panel["reason"])

    return {
        "config_path": str(catalog_path),
        "active_contrasts": int(len(catalog)),
        "excluded_contrasts": int((~include_mask).sum()),
        "source_units": int(len(source_units)),
        "checked_sources": checked_sources,
        "statistical_values_checked": True,
        "required_contrasts_columns": BEGINNER_REQUIRED_CONTRAST_COLUMNS,
        "required_source_table_mappings": _format_source_mapping_contract(REQUIRED_SOURCE_TABLE_MAPPINGS),
        "optional_source_table_mappings": _format_source_mapping_contract(OPTIONAL_SOURCE_TABLE_MAPPINGS),
        "gold_panel_status": gold_panel["status"],
        "n_panel_rows": int(gold_panel.get("n_panel_rows", 0) or 0),
        "gold_panel_source": gold_panel["source"],
        "gold_panel_gene_count": len(gold_panel["genes"]),
        "gold_panel_reason": gold_panel["reason"],
        # Surface the same non-fatal microarray advisories that run_slice emits, so the
        # `degora validate` preflight flags them before a full run rather than after.
        "warnings": list(dict.fromkeys(str(message) for message in warnings if str(message).strip())),
    }


def _harmonized_copy_stem(output_dir: Path, harmonized_dir: Path) -> str:
    """Name the harmonized copy after the output directory without clobbering another run's.

    The copy is `<output_dir.name>_harmonized`, so two runs whose output folders
    are both called `results` shared one copy and the second overwrote the first.
    The sidecar beside an existing copy records the `--output-dir` that wrote it;
    when that is a different directory, this run's copy gets a short suffix
    derived from its own output path instead.
    """

    stem = f"{output_dir.name}_harmonized"
    existing = harmonized_dir / f"{stem}.csv"
    sidecar = existing.with_suffix(existing.suffix + ".source")
    if not existing.exists() or not sidecar.exists():
        return stem
    try:
        command = sidecar.read_text(encoding="utf-8").splitlines()[0]
    except (OSError, IndexError, UnicodeDecodeError):
        return stem
    match = re.search(r"--output-dir\s+(\S+)", command)
    if not match:
        return stem
    recorded = match.group(1)
    recorded_path = Path(recorded)
    if not recorded_path.is_absolute():
        recorded_path = harmonized_dir.resolve() / recorded_path
    try:
        same = recorded_path.resolve() == output_dir.resolve()
    except OSError:
        same = False
    if same:
        return stem
    digest = hashlib.sha256(str(output_dir.resolve()).encode("utf-8")).hexdigest()[:8]
    return f"{output_dir.name}_{digest}_harmonized"


def run_slice(catalog_path: Path, output_dir: Path, harmonized_dir: Path, min_studies: int) -> dict[str, Any]:
    """Harmonize the catalog and write the slice outputs, holding the output directory."""

    try:
        min_studies = validate_min_studies(min_studies)
    except ValueError as exc:
        # min_studies <= 0 disables the replication floor entirely (every gene passes
        # n_studies.ge(min_studies)); reject it so a meaningless value is not recorded
        # as legitimate. 1 is intentionally supported (single-source scoring).
        raise DegoraConfigError(
            "min_studies must be at least 1",
            problems=[f"Got min_studies={min_studies!r}."],
            fixes=["Use 1 to score single-source genes, or 2+ to require independent replication."],
        ) from exc
    # Preserve the public, actionable configuration error before the generic
    # lock helper tries to create the output directory.  The check is repeated
    # under the lock below so a path change between preflight and acquisition
    # still fails closed.
    for label, directory in (("output", output_dir), ("harmonized", harmonized_dir)):
        if directory.exists() and not directory.is_dir():
            raise DegoraConfigError(
                f"the {label} path is a file, not a folder",
                context=f"{label} directory: {directory}",
                problems=["A file already exists at this path, so no folder can be created there."],
                fixes=["Point the output at a folder (for example outputs/results/degora-run), or move the file out of the way."],
            )
    with artifact_output_lock(output_dir):
        return _run_slice_transaction(catalog_path, output_dir, harmonized_dir, min_studies)


def _run_slice_transaction(
    catalog_path: Path,
    output_dir: Path,
    harmonized_dir: Path,
    min_studies: int,
) -> dict[str, Any]:
    # min_studies validation remains disk-free. Create both directories only after
    # the complete publication scope is locked; later path or catalog failures may
    # leave an empty directory because atomicity applies to published artifacts, not
    # to the directory's existence.
    for label, directory in (("output", output_dir), ("harmonized", harmonized_dir)):
        if directory.exists() and not directory.is_dir():
            raise DegoraConfigError(
                f"the {label} path is a file, not a folder",
                context=f"{label} directory: {directory}",
                problems=["A file already exists at this path, so no folder can be created there."],
                fixes=["Point the output at a folder (for example outputs/results/degora-run), or move the file out of the way."],
            )
        try:
            directory.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise DegoraConfigError(
                f"could not create the {label} directory",
                context=f"{label} directory: {directory}",
                problems=[
                    "A file stands where a folder is needed on this path."
                    if isinstance(exc, (FileExistsError, NotADirectoryError))
                    else str(exc)
                ],
                fixes=[
                    "Choose a path that is writable and is not an existing file.",
                    f"Check the permissions on the parent folder of {directory}.",
                ],
            ) from exc

    return _run_slice_locked(catalog_path, output_dir, harmonized_dir, min_studies)


def _run_slice_locked(catalog_path: Path, output_dir: Path, harmonized_dir: Path, min_studies: int) -> dict[str, Any]:
    catalog_path = catalog_path.resolve()
    full_catalog = read_catalog(catalog_path)
    include_mask = catalog_include_mask(full_catalog)
    _validate_catalog_required_values(full_catalog, include_mask, catalog_path)
    _validate_optional_scope_values(full_catalog, include_mask, catalog_path)
    _validate_source_unit_time_course_modes(full_catalog, include_mask, catalog_path)
    _validate_time_course_durations(full_catalog, include_mask, catalog_path)
    _validate_optional_replicate_counts(full_catalog, include_mask, catalog_path)
    _reject_duplicate_active_study_ids(full_catalog, include_mask, catalog_path)
    _reject_delimiter_in_identifiers(full_catalog, include_mask, catalog_path)
    catalog = full_catalog.loc[include_mask].copy()
    excluded_catalog = full_catalog.loc[~include_mask].copy()
    if catalog.empty:
        raise DegoraConfigError(
            "catalog contains no active contrasts",
            context=f"config file: {catalog_path}",
            problems=["Every row is excluded by include_in_analysis/include, or the catalog is empty."],
            fixes=["Set include_in_analysis/include to yes for at least one contrast row."],
        )
    _validate_generated_replicate_provenance(catalog, catalog_path)
    _report_all_missing_source_files(catalog, catalog_path)

    harmonized_tables = []
    source_inputs: list[Path] = []
    input_warnings: list[str] = [
        *catalog.attrs.get("promoted_alias_warnings", []),
        *_microarray_warnings(catalog),
        *_mixed_species_warnings(catalog),
        *_reject_duplicate_source_tables(catalog, catalog_path),
    ]
    filter_summaries: dict[str, dict[str, Any]] = {}

    for row in catalog.to_dict(orient="records"):
        source_path = _resolve_source_path(row["source_path"], catalog_path)
        _require_readable_source_file(source_path, row, catalog_path)
        source_inputs.append(source_path)

        mapping = TableMapping(
            gene_column=row["gene_column"],
            lfc_column=row["lfc_column"],
            p_column=row["p_column"],
            padj_column=_nonempty(row.get("padj_column")),
            sep=_nonempty(row.get("sep")),
            sheet_name=_nonempty(row.get("sheet_name")),
            header_row=normalize_header_row(row.get("header_row")),
        )
        try:
            raw_frame = read_deg_table(source_path, mapping)
        except Exception as exc:
            raise DegoraConfigError(
                "source DEG table could not be read",
                context=_source_path_context(row, source_path, catalog_path),
                problems=[_readable_source_read_failure(source_path, catalog_path, exc)],
                fixes=[
                    "Check that the file is a supported CSV/TSV/TXT/XLS/XLSX table.",
                    "For Excel sources, set sheet_name to the exact sheet containing the DEG table.",
                    "For TSV files, set sep to \\t if DEGORA cannot infer it.",
                ],
            ) from exc
        _validate_source_columns(raw_frame, mapping, row, source_path, catalog_path)
        input_warnings.extend(_validate_numeric_source_columns(raw_frame, mapping, row, source_path, catalog_path))
        input_warnings.extend(_ignored_workbook_settings(row, source_path))
        filtered_frame, filter_summary = apply_gene_type_filter(
            raw_frame,
            _nonempty(row.get("gene_type_column")),
            _nonempty(row.get("gene_type_keep")),
        )
        filter_summaries[row["study_id"]] = filter_summary
        row["source_path"] = portable_path(source_path, catalog_path.parent)
        try:
            harmonized = harmonize_frame(filtered_frame, mapping, row)
        except ValueError as exc:
            raise DegoraConfigError(
                "source DEG table failed harmonization",
                context=_source_path_context(row, source_path, catalog_path),
                problems=[str(exc)],
                fixes=[
                    "Check that gene_column, lfc_column, p_column, and padj_column (if filled) map to the intended source columns.",
                    "p_column must contain raw p-values in [0, 1]; do not map -log10(p), test statistics, or percentages.",
                    "Exact p=0 is accepted and handled automatically, but negative p-values and values above 1 are rejected.",
                ],
            ) from exc
        if harmonized.empty:
            input_warnings.append(f"{row['study_id']} produced zero harmonized rows")
        harmonized_tables.append(harmonized)

    if not harmonized_tables:
        raise ValueError("Catalog contains no usable studies")

    all_harmonized = pd.concat(harmonized_tables, ignore_index=True)
    harmonized_stem = _harmonized_copy_stem(output_dir, harmonized_dir)
    harmonized_csv = harmonized_dir / f"{harmonized_stem}.csv"
    harmonized_parquet = harmonized_dir / f"{harmonized_stem}.parquet"
    neutralize_formula_text(all_harmonized).to_csv(harmonized_csv, index=False)
    optional_output_warnings: list[str] = []
    parquet_warning = _try_write_parquet(all_harmonized, harmonized_parquet)
    if parquet_warning:
        optional_output_warnings.append(parquet_warning)

    consensus = slice_consensus(all_harmonized, min_studies=min_studies)
    consensus_path = output_dir / "slice_consensus.csv"
    neutralize_formula_text(consensus).to_csv(consensus_path, index=False)

    result_harmonized_path = output_dir / "slice_harmonized.csv"
    neutralize_formula_text(all_harmonized).to_csv(result_harmonized_path, index=False)
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
    unusable_row_warnings = sorted(
        {
            str(value)
            for value in all_harmonized.get("unusable_row_warning", pd.Series(dtype=str)).dropna().unique()
            if str(value).strip()
        }
    )
    input_warnings.extend(gene_symbol_collapse_warnings)
    input_warnings.extend(unusable_row_warnings)
    if int(min_studies) == 1:
        # The default floor of two is what makes a score replicated evidence;
        # a run below it used to say nothing about what it had given up.
        input_warnings.append(
            "min_studies=1: a gene is scored from a single source unit, so this ranking shows no replication; "
            "treat it as exploratory prioritisation, not as replicated evidence."
        )

    # What early/late preselection kept, and where it left a source unit
    # contributing a small minority of the genes it started with.
    time_course_report = time_course_selection_report(all_harmonized)
    time_course_report_warnings = time_course_selection_warnings(time_course_report)
    input_warnings.extend(time_course_report_warnings)

    gold_panel = _read_locked_gold_panel(catalog_path)
    if gold_panel["status"] == "locked":
        recall50 = recall_at_k(consensus, gold_panel["genes"], 50)
        recall100 = recall_at_k(consensus, gold_panel["genes"], 100)
    else:
        recall50 = {
            "status": "not_applicable",
            "reason": "no locked topic-specific GoldPanel was provided; this run has no benchmark recall to interpret",
        }
        recall100 = {
            "status": "not_applicable",
            "reason": "no locked topic-specific GoldPanel was provided; this run has no benchmark recall to interpret",
        }

    source_inputs_sorted = [Path(value) for value in sorted({str(path) for path in source_inputs})]
    metrics = {
        # Every recorded path is relative to the config directory so metrics
        # remain shareable without leaking a username or workstation layout.
        "path_base": "catalog_directory",
        "catalog_path": portable_path(catalog_path, catalog_path.parent),
        "source_input_files": [
            portable_path(path, catalog_path.parent) for path in source_inputs_sorted
        ],
        "n_catalog_rows": int(len(full_catalog)),
        "n_active_catalog_rows": int(len(catalog)),
        "excluded_catalog_rows": excluded_catalog[["study_id", "notes"]].fillna("").to_dict(orient="records"),
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
        "identifier_space_warnings": _identifier_space_warnings(all_harmonized, min_studies=min_studies),
        "time_course_selection": time_course_report,
        "time_course_selection_warnings": time_course_report_warnings,
        "time_course_selection_warning_fraction": TIME_COURSE_RETENTION_WARNING_FRACTION,
        "rank_universe_warnings": rank_universe_warnings,
        "gene_symbol_collapse_warnings": gene_symbol_collapse_warnings,
        "n_gene_symbols_normalized": _count_normalized_gene_symbols(all_harmonized),
        "unusable_row_warnings": unusable_row_warnings,
        # rows_before minus the harmonized count used to be the only way to see
        # this, and it could not say which column was responsible.
        "unusable_row_counts": {
            str(study_id): int(count)
            for study_id, count in all_harmonized.groupby("study_id")["n_rows_dropped_unusable"]
            .max()
            .items()
            if int(count)
        }
        if "n_rows_dropped_unusable" in all_harmonized.columns
        else {},
        "input_row_counts": {
            str(study_id): int(count)
            for study_id, count in all_harmonized.groupby("study_id")["n_input_rows"].max().items()
        }
        if "n_input_rows" in all_harmonized.columns
        else {},
        "pipeline_counts": _count_labels(catalog["pipeline"]),
        "assay_type_counts": _count_labels(catalog["assay_type"]) if "assay_type" in catalog.columns else {},
        "source_input_type_counts": _count_labels(catalog["source_input_type"])
        if "source_input_type" in catalog.columns
        else {},
        "warnings": input_warnings,
        "optional_output_warnings": optional_output_warnings,
        "gold_panel_status": gold_panel["status"],
        "n_panel_rows": int(gold_panel.get("n_panel_rows", 0) or 0),
        "gold_panel_source": gold_panel["source"],
        "gold_panel_gene_count": len(gold_panel["genes"]),
        "gold_panel_reason": gold_panel["reason"],
        "recall_rank_source": "slice_consensus_order",
        "recall_rank_note": "slice_metrics recall is computed before score_db and is not the primary quality-weighted DEGORA rank.",
        "recall_at_50": recall50,
        "recall_at_100": recall100,
    }
    metrics_path = output_dir / "slice_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    command = shell_command(
        [
            "python3",
            "-m",
            "degora.slice_runner",
            "--catalog",
            catalog_path,
            "--output-dir",
            output_dir.resolve(),
            "--harmonized-dir",
            harmonized_dir.resolve(),
            "--min-studies",
            min_studies,
        ]
    )
    artifacts = [harmonized_csv, consensus_path, result_harmonized_path, metrics_path]
    if harmonized_parquet.exists():
        artifacts.insert(1, harmonized_parquet)
    for artifact in artifacts:
        metadata = {"generator": "slice", "min_studies": min_studies}
        if artifact.suffix.lower() == ".csv":
            metadata.update(formula_guard_metadata())
        write_source_sidecar(
            artifact,
            command,
            inputs=[catalog_path, *source_inputs_sorted],
            metadata=metadata,
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
