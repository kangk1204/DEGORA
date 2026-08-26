"""Activate reviewed discovery candidates and run one species-specific DEGORA analysis."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .discovery import DiscoveryError, normalize_species
from .excel_export import DEFAULT_WORKBOOK_NAME, export_run_workbook
from .formula_safety import formula_guard_metadata, neutralize_formula_text
from .harmonize import _read_excel_any, _restore_unnamed_row_labels
from .provenance import shell_command, write_source_sidecar
from .reanalysis import derive_welch_deg
from .score_db import write_score_database
from .slice_runner import CATALOG_COLUMNS, run_slice, validate_catalog_inputs

MAX_ACTIVE_CANDIDATES = 40
MAX_CONTRAST_LABEL = 180


def _candidate_index(prepared: dict[str, Any]) -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    index: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
    for study in prepared.get("studies", []):
        for candidate in study.get("files", []):
            candidate_id = str(candidate.get("candidate_id") or "")
            if not candidate_id or candidate_id in index:
                raise DiscoveryError("prepared bundle contains a missing or duplicate candidate_id")
            index[candidate_id] = (study, candidate)
    return index


def _text(value: Any, *, field: str, required: bool = False, maximum: int = MAX_CONTRAST_LABEL) -> str:
    text = str(value or "").strip()
    if required and not text:
        raise DiscoveryError(f"{field} is required")
    if len(text) > maximum:
        raise DiscoveryError(f"{field} is too long; maximum length is {maximum} characters")
    if any(ord(char) < 32 for char in text):
        raise DiscoveryError(f"{field} contains control characters")
    return text


def _optional_count(value: Any, *, field: str, required: bool = False) -> int | str:
    message = f"{field} must be a positive whole number" + ("" if required else " or blank")
    if value is None:
        if required:
            raise DiscoveryError(message)
        return ""
    if isinstance(value, bool):
        raise DiscoveryError(message)
    text = str(value).strip()
    if not text:
        if required:
            raise DiscoveryError(message)
        return ""
    if not re.fullmatch(r"[0-9]+", text):
        raise DiscoveryError(message)
    number = int(text)
    if number < 1:
        raise DiscoveryError(message)
    return number


def _assay_type(study_type: str) -> str:
    value = study_type.lower()
    if "array" in value:
        return "microarray"
    if "sequencing" in value or "high throughput" in value:
        return "RNA-seq"
    return "expression profiling"


def _author_pipeline(filename: str) -> str:
    value = filename.lower()
    if "deseq" in value:
        return "DESeq2"
    if "edger" in value:
        return "edgeR"
    if "limma" in value or "toptable" in value:
        return "limma"
    if "cuffdiff" in value or "gene_exp.diff" in value:
        return "Cuffdiff"
    return "author_reported_pipeline_unspecified"


def _first_text(*values: Any) -> str:
    for value in values:
        if isinstance(value, (list, tuple)):
            nested = _first_text(*value)
            if nested:
                return nested
            continue
        text = str(value or "").strip()
        if text:
            return text
    return ""


def _doi_unit(value: Any) -> str:
    doi = _first_text(value)
    if not doi:
        return ""
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi, flags=re.IGNORECASE).strip()
    doi = re.sub(r"^doi:\s*", "", doi, flags=re.IGNORECASE).strip()
    doi = doi.rstrip(" .")
    if not doi:
        return ""
    return f"DOI:{doi.lower()}"


def _provider_accession_unit(study: dict[str, Any]) -> str:
    canonical_id = _first_text(study.get("canonical_id"), study.get("provider_accession"), study.get("accession"))
    if not canonical_id:
        return ""
    provider = _first_text(study.get("provider"), study.get("source_provider"))
    if provider and not str(canonical_id).upper().startswith(("GSE", "PMID:", "DOI:", "PMC")):
        return f"{provider.upper()}:{canonical_id}"
    return str(canonical_id).upper()


def _contained_local_path(candidate: dict[str, Any], bundle_root: Path) -> Path:
    inspection = candidate.get("inspection", {})
    raw = inspection.get("local_path", "")
    if not raw:
        raise DiscoveryError("selected candidate was inspected but not fully downloaded")
    path = Path(str(raw)).resolve()
    if not path.is_relative_to(bundle_root):
        raise DiscoveryError("selected candidate path falls outside its discovery bundle")
    if not path.is_file():
        raise FileNotFoundError(f"selected candidate file is missing: {path.name}")
    if inspection.get("fetch_scope") != "full":
        raise DiscoveryError("selected candidate is only a header preview; prepare it with full materialization first")
    expected_sha256 = str(inspection.get("full_file_sha256") or "").lower()
    if not re.fullmatch(r"[a-f0-9]{64}", expected_sha256):
        raise DiscoveryError("selected candidate is missing its preparation SHA-256 integrity record")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != expected_sha256:
        raise DiscoveryError("selected candidate changed after preparation; download and inspect it again")
    return path


def _without_list_delimiter(value: str) -> str:
    """Strip the semicolon DEGORA joins identifier lists with.

    A DOI may legitimately contain one, and it would then reach the catalog as a
    source_unit_id that makes every semicolon-joined provenance list ambiguous.
    """

    return value.replace(";", "_")


def _paper_source_unit(study: dict[str, Any]) -> str:
    explicit = str(study.get("source_unit_id") or "").strip()
    if explicit:
        return _without_list_delimiter(
            _text(explicit, field="source_unit_id", required=True, maximum=160)
        )
    pmids = [
        str(value).strip()
        for value in [*study.get("source_unit_pubmed_ids", []), *study.get("pubmed_ids", [])]
        if str(value).strip()
    ]
    if pmids:
        return _without_list_delimiter(f"PMID:{pmids[0]}")
    doi = _doi_unit(_first_text(study.get("doi"), study.get("dois"), study.get("publication_doi")))
    if doi:
        return _without_list_delimiter(doi)
    provider_accession = _provider_accession_unit(study)
    if provider_accession:
        return _without_list_delimiter(provider_accession)
    raise DiscoveryError("prepared study is missing a stable publication or data identifier")


def _study_accession_key(study: dict[str, Any]) -> str:
    key = _first_text(study.get("accession"), study.get("canonical_id"), study.get("provider_accession"))
    if not key:
        key = _paper_source_unit(study)
    key = re.sub(r"[^A-Za-z0-9_.:-]+", "_", key).strip("_")
    if not key:
        raise DiscoveryError("prepared study is missing a deterministic catalog accession")
    return key[:80]


def _publication_metadata_note(study: dict[str, Any]) -> str:
    parts = []
    doi = _doi_unit(_first_text(study.get("doi"), study.get("dois"), study.get("publication_doi")))
    pmcid = _first_text(study.get("pmcid"), study.get("pmcids"))
    if doi:
        parts.append(doi)
    if pmcid:
        pmcid = pmcid.removeprefix("PMCID:").removeprefix("PMC")
        parts.append(f"PMCID:PMC{pmcid}")
    return "; ".join(parts)


def _source_url(study: dict[str, Any], candidate: dict[str, Any]) -> str:
    url = _first_text(candidate.get("source_url"), study.get("source_url"))
    if url:
        return url
    doi = _doi_unit(_first_text(study.get("doi"), study.get("dois"), study.get("publication_doi")))
    if doi:
        return "https://doi.org/" + doi.removeprefix("DOI:")
    return ""


def _mixed_status(study: dict[str, Any]) -> str:
    status = _first_text(
        study.get("mixed_status"),
        study.get("mixed_activation_status"),
        study.get("organism_status"),
        study.get("species_scope_status"),
        study.get("target_species_status"),
        study.get("degora_mixed_status"),
    )
    if status:
        return status.lower()
    if study.get("mixed_blocked") is True:
        return "mixed_blocked"
    if study.get("mixed_quarantined") is True:
        return "mixed_quarantined"
    if study.get("mixed_rescued") is True:
        return "mixed_rescued"
    return ""


def _validate_mixed_activation(study: dict[str, Any]) -> None:
    status = _mixed_status(study)
    if status in {"mixed_blocked", "mixed_quarantined"}:
        raise DiscoveryError(f"{status} study activation is not allowed")
    if status == "mixed_rescued":
        evidence = _first_text(
            study.get("target_species_evidence"),
            study.get("mixed_rescue_evidence"),
            study.get("species_evidence"),
            study.get("evidence_text"),
            study.get("evidence"),
        )
        if study.get("target_species_verified") is not True or not evidence:
            raise DiscoveryError(
                "mixed_rescued study activation requires target_species_verified=true and nonempty evidence text"
            )


def _validate_prepared_source_units(studies: Iterable[dict[str, Any]]) -> None:
    pmid_units: dict[str, str] = {}
    for study in studies:
        source_unit = _paper_source_unit(study)
        pmids = [
            str(value).strip()
            for value in [*study.get("source_unit_pubmed_ids", []), *study.get("pubmed_ids", [])]
            if str(value).strip()
        ]
        for pmid in pmids:
            previous = pmid_units.setdefault(pmid, source_unit)
            if previous != source_unit:
                raise DiscoveryError(
                    "prepared studies sharing a PubMed ID must share one source_unit_id; prepare the bundle again"
                )


def _base_catalog_row(
    *,
    study: dict[str, Any],
    candidate: dict[str, Any],
    entry: dict[str, Any],
    spec,
    source_path: Path,
    sequence: int,
) -> dict[str, Any]:
    accession = _study_accession_key(study)
    source_unit_id = _paper_source_unit(study)
    contrast = _text(entry.get("contrast_label"), field="contrast_label", required=True)
    row = {column: "" for column in CATALOG_COLUMNS}
    row.update(
        {
            "study_id": f"{spec.key}_{accession}_{sequence:03d}_{str(candidate['candidate_id'])[:8]}",
            "paper_id": source_unit_id,
            "source_unit_id": source_unit_id,
            "source_url": _source_url(study, candidate),
            "source_path": str(source_path),
            "species": spec.scientific_name,
            "cell_system": _text(entry.get("cell_system"), field="cell_system"),
            "hypoxia_modality": contrast,
            "duration_h": _text(entry.get("duration_h"), field="duration_h", maximum=32),
            "n_ctrl": _optional_count(entry.get("n_ctrl"), field="n_ctrl"),
            "n_treat": _optional_count(entry.get("n_treat"), field="n_treat"),
            "assay_type": _text(entry.get("assay_type"), field="assay_type", maximum=80)
            or _assay_type(str(study.get("study_type") or "")),
            "platform": _text(entry.get("platform"), field="platform", maximum=80),
            "time_course_mode": "mean",
            "table_scope": _text(entry.get("table_scope") or "auto", field="table_scope", maximum=32),
            "include_in_analysis": "yes",
        }
    )
    return row


AUTHOR_REVIEWABLE_STATUSES = frozenset(
    {
        "ready_for_review",
        "requires_column_mapping",
        "requires_lfc_confirmation",
        "requires_pvalue_mapping",
    }
)
AUTHOR_MAPPING_FIELDS = ("gene_column", "lfc_column", "p_column", "padj_column")
AUTHOR_DUPLICATE_GENE_POLICIES = frozenset({"harmonizer", "keep_first"})


def _author_mapping(inspection: dict[str, Any], entry: dict[str, Any]) -> dict[str, str]:
    detected = inspection.get("mapping", {})
    mapping = {
        field: _text(
            entry.get(field) if field in entry else detected.get(field),
            field=field,
            maximum=160,
        )
        for field in AUTHOR_MAPPING_FIELDS
    }
    for required in ("gene_column", "lfc_column", "p_column"):
        if not mapping[required]:
            raise DiscoveryError(f"{required} is required for an author DEG table")

    status = str(inspection.get("status") or "")
    explicit_changes = {
        field
        for field in AUTHOR_MAPPING_FIELDS
        if field in entry
        and str(entry.get(field) or "").strip() != str(detected.get(field) or "").strip()
    }
    sheet_changed = bool(
        str(entry.get("sheet_name") or "").strip()
        and str(entry.get("sheet_name")).strip() != str(inspection.get("sheet_name") or "").strip()
    )
    if (status != "ready_for_review" or explicit_changes or sheet_changed) and entry.get("column_mapping_confirmed") is not True:
        raise DiscoveryError(
            "column_mapping_confirmed=true is required when activating a table that needs explicit column or sheet review"
        )
    if status == "requires_pvalue_mapping" and not str(entry.get("p_column") or "").strip():
        raise DiscoveryError("requires_pvalue_mapping candidates need an explicit p_column selection")
    if mapping["padj_column"] and mapping["p_column"] == mapping["padj_column"]:
        if entry.get("adjusted_p_as_pvalue_confirmed") is not True:
            raise DiscoveryError(
                "adjusted_p_as_pvalue_confirmed=true is required when the adjusted-p/FDR column is also used as p_column"
            )
    if status == "requires_lfc_confirmation" and entry.get("lfc_scale_confirmed_log2") is not True:
        raise DiscoveryError("lfc_scale_confirmed_log2=true is required for an ambiguous effect-size column")
    return mapping


def _read_author_frame(path: Path, *, sheet_name: str, header_row: int) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    header = max(header_row - 1, 0)
    if suffixes.endswith((".xlsx", ".xls", ".xlsx.gz", ".xls.gz")):
        selected_sheet: str | int = sheet_name or 0
        return _read_excel_any(path, sheet_name=selected_sheet, header_row=header_row)
    separator = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return pd.read_csv(path, sep=separator, header=header)


def _author_filter(entry: dict[str, Any]) -> tuple[str, str]:
    column = _text(entry.get("row_filter_column"), field="row_filter_column", maximum=160)
    value = _text(entry.get("row_filter_value"), field="row_filter_value", maximum=240)
    if bool(column) != bool(value):
        raise DiscoveryError("row_filter_column and row_filter_value must be supplied together")
    if column and entry.get("row_filter_confirmed") is not True:
        raise DiscoveryError("row_filter_confirmed=true is required for an author-table row subset")
    return column, value


def _author_duplicate_gene_policy(entry: dict[str, Any]) -> str:
    policy = _text(
        entry.get("duplicate_gene_policy") or "harmonizer",
        field="duplicate_gene_policy",
        maximum=32,
    )
    if policy not in AUTHOR_DUPLICATE_GENE_POLICIES:
        raise DiscoveryError("duplicate_gene_policy must be harmonizer or keep_first")
    if policy == "keep_first" and entry.get("duplicate_gene_policy_confirmed") is not True:
        raise DiscoveryError(
            "duplicate_gene_policy_confirmed=true is required when preserving the first usable source row"
        )
    return policy


def _author_activation_key(candidate_id: str, inspection: dict[str, Any], entry: dict[str, Any]) -> tuple[str, ...]:
    if inspection.get("status") not in AUTHOR_REVIEWABLE_STATUSES:
        raise DiscoveryError("author DEG candidate must pass header review before activation")
    mapping = _author_mapping(inspection, entry)
    filter_column, filter_value = _author_filter(entry)
    _author_duplicate_gene_policy(entry)
    sheet_name = _text(
        entry.get("sheet_name") or inspection.get("sheet_name"),
        field="sheet_name",
        maximum=160,
    )
    return (
        candidate_id,
        sheet_name,
        *(mapping[field] for field in AUTHOR_MAPPING_FIELDS),
        filter_column,
        filter_value,
    )


def _materialize_author_table(
    *,
    study: dict[str, Any],
    candidate: dict[str, Any],
    entry: dict[str, Any],
    source_path: Path,
    derived_dir: Path,
    sequence: int,
    replay_command: str,
) -> tuple[Path, dict[str, Any], dict[str, str]]:
    inspection = candidate.get("inspection", {})
    mapping = _author_mapping(inspection, entry)
    sheet_name = _text(
        entry.get("sheet_name") or inspection.get("sheet_name"),
        field="sheet_name",
        maximum=160,
    )
    try:
        header_row = int(inspection.get("header_row") or 1)
    except (TypeError, ValueError) as exc:
        raise DiscoveryError("prepared author candidate has an invalid header_row") from exc
    if not 1 <= header_row <= 100:
        raise DiscoveryError("prepared author candidate header_row must be between 1 and 100")

    try:
        frame = _read_author_frame(source_path, sheet_name=sheet_name, header_row=header_row)
    except (KeyError, ValueError) as exc:
        raise DiscoveryError(f"selected author table sheet could not be read: {exc}") from exc
    frame.columns = [str(value).strip() for value in frame.columns]
    # The inspector and read_deg_table both call an R export's unnamed label
    # column `row_name`; the materialised frame has to carry the same name, or
    # the mapping the reader just confirmed is "not found" in the file it came from.
    frame = _restore_unnamed_row_labels(frame)
    required_columns = [mapping["gene_column"], mapping["lfc_column"], mapping["p_column"]]
    if mapping["padj_column"]:
        required_columns.append(mapping["padj_column"])
    missing = sorted({column for column in required_columns if column not in frame.columns})
    if missing:
        raise DiscoveryError("selected author table column(s) were not found: " + ", ".join(missing))

    filter_column, filter_value = _author_filter(entry)
    n_input_rows = int(len(frame))
    if filter_column:
        if filter_column not in frame.columns:
            raise DiscoveryError(f"row_filter_column was not found in the selected author table: {filter_column}")
        normalized = frame[filter_column].astype("string").fillna("").str.strip()
        frame = frame.loc[normalized.eq(filter_value)].copy()
        if frame.empty:
            raise DiscoveryError(
                f"row filter matched no author-table rows: {filter_column}={filter_value!r}"
            )

    genes = frame[mapping["gene_column"]].astype("string").str.strip()
    lfc = pd.to_numeric(frame[mapping["lfc_column"]], errors="coerce")
    pvalue = pd.to_numeric(frame[mapping["p_column"]], errors="coerce")
    valid = genes.notna() & genes.ne("") & lfc.notna() & pvalue.notna() & pvalue.between(0.0, 1.0)
    if int(valid.sum()) < 2:
        raise DiscoveryError("selected author-table mapping produced fewer than two usable gene/effect/p-value rows")
    selected = pd.DataFrame(
        {
            "gene_symbol": genes,
            "log2FoldChange": lfc,
            "pvalue": pvalue,
        }
    )
    if mapping["padj_column"]:
        selected["padj"] = pd.to_numeric(frame[mapping["padj_column"]], errors="coerce")
    selected = selected.loc[valid].reset_index(drop=True)
    duplicate_gene_policy = _author_duplicate_gene_policy(entry)
    n_usable_rows_before_duplicate_policy = int(len(selected))
    duplicate_gene_rows = selected["gene_symbol"].duplicated(keep=False)
    n_duplicate_gene_rows = int(duplicate_gene_rows.sum())
    n_duplicate_genes = int(selected.loc[duplicate_gene_rows, "gene_symbol"].nunique())
    if duplicate_gene_policy == "keep_first":
        selected = selected.drop_duplicates("gene_symbol", keep="first").reset_index(drop=True)

    derived_dir.mkdir(parents=True, exist_ok=True)
    accession = _study_accession_key(study)
    output = derived_dir / f"{sequence:02d}_{accession}_{str(candidate['candidate_id'])[:10]}_author.csv"
    neutralize_formula_text(selected).to_csv(output, index=False, lineterminator="\n")
    provenance = {
        "generator": "degora.discovery_run._materialize_author_table",
        "operation": "author_table_column_projection_and_optional_exact_row_subset",
        "statistical_reanalysis": False,
        "candidate_id": str(candidate.get("candidate_id") or ""),
        "sheet_name": sheet_name,
        "header_row": header_row,
        "source_mapping": mapping,
        "output_mapping": {
            "gene_column": "gene_symbol",
            "lfc_column": "log2FoldChange",
            "p_column": "pvalue",
            "padj_column": "padj" if mapping["padj_column"] else "",
        },
        "row_filter_column": filter_column,
        "row_filter_value": filter_value,
        "n_input_rows": n_input_rows,
        "n_rows_after_filter": int(len(frame)),
        "duplicate_gene_policy": duplicate_gene_policy,
        "duplicate_gene_policy_interpretation": (
            "preserve rows for downstream min_pvalue_max_abs_lfc harmonization"
            if duplicate_gene_policy == "harmonizer"
            else "keep the first usable source row for each gene in original row order"
        ),
        "n_usable_rows_before_duplicate_policy": n_usable_rows_before_duplicate_policy,
        "n_duplicate_gene_rows": n_duplicate_gene_rows,
        "n_duplicate_genes": n_duplicate_genes,
        "n_usable_output_rows": int(len(selected)),
        **formula_guard_metadata(),
    }
    write_source_sidecar(output, replay_command, inputs=[source_path], metadata=provenance)
    return output, provenance, provenance["output_mapping"]


def _author_row(
    *,
    study: dict[str, Any],
    candidate: dict[str, Any],
    entry: dict[str, Any],
    spec,
    bundle_root: Path,
    derived_dir: Path,
    sequence: int,
    replay_command: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    inspection = candidate.get("inspection", {})
    if inspection.get("status") not in AUTHOR_REVIEWABLE_STATUSES:
        raise DiscoveryError("author DEG candidate must pass header review before activation")
    if entry.get("direction_confirmed") is not True:
        raise DiscoveryError("direction_confirmed=true is required; positive log2FC must mean treatment/case minus control")
    if entry.get("invert_lfc"):
        raise DiscoveryError("automatic log2FC inversion is not supported; use a correctly oriented author table")
    author_entry = {
        **entry,
        "n_ctrl": _optional_count(entry.get("n_ctrl"), field="n_ctrl", required=True),
        "n_treat": _optional_count(entry.get("n_treat"), field="n_treat", required=True),
    }
    original_source_path = _contained_local_path(candidate, bundle_root)
    source_path, derivation, mapping = _materialize_author_table(
        study=study,
        candidate=candidate,
        entry=entry,
        source_path=original_source_path,
        derived_dir=derived_dir,
        sequence=sequence,
        replay_command=replay_command,
    )
    row = _base_catalog_row(
        study=study,
        candidate=candidate,
        entry=author_entry,
        spec=spec,
        source_path=source_path,
        sequence=sequence,
    )
    row.update(
        {
            "pipeline": _text(entry.get("pipeline"), field="pipeline", maximum=80)
            or _author_pipeline(str(candidate.get("name") or "")),
            "gene_column": str(mapping.get("gene_column") or ""),
            "lfc_column": str(mapping.get("lfc_column") or ""),
            "p_column": str(mapping.get("p_column") or ""),
            "padj_column": str(mapping.get("padj_column") or ""),
            "sheet_name": "",
            "source_input_type": "author_deg_table",
            "normalization": "author_reported",
            "notes": (
                "Selected from the DEGORA discovery browser. The user explicitly confirmed that positive log2FC "
                "represents treatment/case minus control. DEGORA projected the reviewed author columns and, when "
                "declared, applied one exact row filter; no differential expression was recomputed."
            ),
        }
    )
    metadata_note = _publication_metadata_note(study)
    if metadata_note:
        row["notes"] = f"{row['notes']} Publication metadata: {metadata_note}."
    source_mapping = derivation.get("source_mapping", {})
    if source_mapping.get("padj_column") and source_mapping.get("p_column") == source_mapping.get("padj_column"):
        row["notes"] = (
            f"{row['notes']} The reviewed adjusted-p/FDR column was explicitly used as the p-value input because "
            "the author table did not report a nominal p-value; this conservative semantic substitution is recorded."
        )
    if derivation.get("row_filter_column"):
        row["notes"] = (
            f"{row['notes']} Exact source-row subset: {derivation['row_filter_column']}="
            f"{derivation['row_filter_value']!r}."
        )
    if derivation.get("duplicate_gene_policy") == "keep_first":
        row["notes"] = (
            f"{row['notes']} Duplicate gene symbols retained the first usable source row in original row order "
            "because that legacy/manual extraction rule was explicitly selected."
        )
    return row, derivation


COUNT_SAMPLE_ROWS = 2000
COUNT_WHOLE_NUMBER_SHARE = 0.95


def _require_whole_number_counts(source_path: Path, sample_columns: list[str]) -> None:
    """Refuse matrix_type=count_matrix when the selected columns are not whole numbers."""

    suffixes = "".join(Path(source_path).suffixes).lower()
    try:
        if suffixes.endswith((".xlsx", ".xls", ".xlsx.gz", ".xls.gz")):
            frame = _read_excel_any(Path(source_path), 0, header_row=1).head(COUNT_SAMPLE_ROWS)
        else:
            separator = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
            frame = pd.read_csv(source_path, sep=separator, nrows=COUNT_SAMPLE_ROWS)
    except Exception:  # noqa: BLE001 - the derivation reports an unreadable file itself
        return
    columns = [name for name in sample_columns if name in frame.columns]
    if not columns:
        return
    values = pd.to_numeric(frame[columns].stack(), errors="coerce").dropna()
    if values.empty:
        return
    whole = float(((values - values.round()).abs() < 1e-9).mean())
    if whole < COUNT_WHOLE_NUMBER_SHARE:
        raise DiscoveryError(
            f"matrix_type=count_matrix, but only {whole:.0%} of the selected columns' values in the first "
            f"{COUNT_SAMPLE_ROWS:,} rows are whole numbers (preflight; the derivation checks every row); "
            "this looks like a normalized matrix (FPKM, TPM, CPM or a log scale). Select it as "
            "normalized_expression_matrix with normalized_scale log2 or linear, or choose the raw count file."
        )


def _fallback_row(
    *,
    study: dict[str, Any],
    candidate: dict[str, Any],
    entry: dict[str, Any],
    spec,
    bundle_root: Path,
    derived_dir: Path,
    sequence: int,
    replay_command: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    inspection = candidate.get("inspection", {})
    if inspection.get("status") != "upstream_matrix_ready_for_contrast":
        raise DiscoveryError("fallback candidate must pass matrix inspection before activation")
    if entry.get("direction_confirmed") is not True:
        raise DiscoveryError("direction_confirmed=true is required for treatment-minus-control fallback analysis")
    if entry.get("biological_replicates_confirmed") is not True:
        raise DiscoveryError(
            "biological_replicates_confirmed=true is required; technical replicates, paired samples, and cells "
            "must not be treated as independent Welch replicates"
        )
    source_path = _contained_local_path(candidate, bundle_root)
    allowed_samples = set(map(str, inspection.get("sample_columns", [])))
    control_samples = [str(value) for value in entry.get("control_samples", [])]
    treatment_samples = [str(value) for value in entry.get("treatment_samples", [])]
    requested = set(control_samples + treatment_samples)
    unknown = sorted(requested.difference(allowed_samples))
    if unknown:
        raise DiscoveryError("selected sample column(s) were not found in the inspected matrix: " + ", ".join(unknown))
    role = str(candidate.get("role") or inspection.get("declared_role") or "")
    if role == "unknown_matrix":
        role = _text(entry.get("matrix_type"), field="matrix_type", required=True, maximum=40)
        if role not in {"count_matrix", "normalized_expression_matrix"}:
            # Reader-correctable input; the derivation raised a bare ValueError for it.
            raise DiscoveryError(
                f"matrix_type={role!r} is not recognised; use count_matrix for raw counts or "
                "normalized_expression_matrix for a normalized matrix (with normalized_scale log2 or linear)"
            )
    if role == "count_matrix":
        # A fractional matrix (FPKM, TPM, a log scale) selected as raw counts
        # would be handed to a count model as if it were counts. The values say
        # which it is before any test is run.
        _require_whole_number_counts(source_path, control_samples + treatment_samples)
    normalized_scale = ""
    if role == "normalized_expression_matrix":
        normalized_scale = _text(entry.get("normalized_scale"), field="normalized_scale", required=True, maximum=16)
        if normalized_scale not in {"log2", "linear"}:
            raise DiscoveryError("normalized_scale must be log2 or linear for a normalized expression matrix")
    gene_column = _text(entry.get("gene_column") or inspection.get("gene_column"), field="gene_column", required=True)
    accession = _study_accession_key(study)
    derived_path = derived_dir / f"{spec.key}_{accession}_{str(candidate['candidate_id'])[:10]}_welch.csv"
    summary = derive_welch_deg(
        source_path,
        derived_path,
        role=role,
        gene_column=gene_column,
        control_samples=control_samples,
        treatment_samples=treatment_samples,
        normalized_scale=normalized_scale or None,
        sheet_name=inspection.get("sheet_name") or None,
        command=replay_command,
        metadata={
            "accession": study.get("accession", ""),
            "species": spec.scientific_name,
            "source_url": _source_url(study, candidate),
            "biological_replicates_confirmed": True,
            "inference_scope": "exploratory_screening_not_confirmatory",
        },
    )
    row = _base_catalog_row(
        study=study,
        candidate=candidate,
        entry={**entry, "n_ctrl": summary["n_ctrl"], "n_treat": summary["n_treat"]},
        spec=spec,
        source_path=derived_path,
        sequence=sequence,
    )
    row.update(
        {
            "pipeline": summary["pipeline"],
            "gene_column": summary["gene_column"],
            "lfc_column": summary["lfc_column"],
            "p_column": summary["p_column"],
            "padj_column": summary["padj_column"],
            "source_input_type": summary["source_input_type"],
            "normalization": summary["normalization"],
            "probe_collapse": "median_expression" if role != "count_matrix" else "",
            "table_scope": "full_results",
            "notes": (
                "Labeled fallback derived from a public matrix by the documented Welch workflow. "
                "The user attested that selected columns are independent biological replicates; effect direction "
                "is treatment minus control. This fallback is exploratory screening, not confirmatory inference."
            ),
        }
    )
    metadata_note = _publication_metadata_note(study)
    if metadata_note:
        row["notes"] = f"{row['notes']} Publication metadata: {metadata_note}."
    return row, summary


DISCOVERY_RUN_MARKER = ".degora-discovery-run.json"
DISCOVERY_RUN_ARTIFACT_TYPE = "degora_discovery_analysis"
DISCOVERY_RUN_FORMAT_VERSION = 1


def _recognized_discovery_output(output: Path) -> bool:
    marker = output / DISCOVERY_RUN_MARKER
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    return bool(
        isinstance(payload, dict)
        and payload.get("artifact_type") == DISCOVERY_RUN_ARTIFACT_TYPE
        and payload.get("format_version") == DISCOVERY_RUN_FORMAT_VERSION
    )


def _begin_output_transaction(output: Path, *, force: bool) -> tuple[bool, Path | None]:
    """Reserve the final directory while keeping a forced prior run recoverable."""

    if output.parent == output or output == Path.home().resolve() or output == Path.cwd().resolve():
        raise DiscoveryError("analysis output must be a dedicated subdirectory")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not output.is_dir():
        raise FileExistsError(f"analysis output exists and is not a directory: {output}")
    existed_empty = output.exists() and not any(output.iterdir())
    backup: Path | None = None
    if output.exists() and not existed_empty:
        if not force:
            raise FileExistsError(f"analysis output already exists and is not empty: {output}")
        if not _recognized_discovery_output(output):
            raise DiscoveryError("refusing --force because the output is not a recognized DEGORA discovery run")
        backup = Path(tempfile.mkdtemp(prefix=f".{output.name}.backup-", dir=output.parent))
        backup.rmdir()
        output.replace(backup)
    try:
        output.mkdir(parents=True, exist_ok=True)
    except BaseException:
        if backup is not None and backup.exists():
            backup.replace(output)
        raise
    return existed_empty, backup


def _rollback_output_transaction(output: Path, *, existed_empty: bool, backup: Path | None) -> None:
    if output.exists():
        shutil.rmtree(output)
    if backup is not None and backup.exists():
        backup.replace(output)
    elif existed_empty:
        output.mkdir(parents=True, exist_ok=True)


def _commit_output_transaction(backup: Path | None) -> None:
    if backup is not None and backup.exists():
        shutil.rmtree(backup)


def _execute_discovery_analysis(
    prepared: dict[str, Any],
    selections: Iterable[dict[str, Any]],
    output_dir: str | Path,
    *,
    species: str,
    min_studies: int = 2,
    force: bool = False,
    extra_metadata: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Build an active catalog and run DEGORA for exactly one species."""

    spec = normalize_species(species)
    # A hand-edited bundle can carry species as a bare string or null. Reading
    # .get on that raised AttributeError and ended the command in a traceback,
    # for input the reader can correct.
    prepared_species_field = prepared.get("species")
    prepared_species = (
        str(prepared_species_field.get("key") or "")
        if isinstance(prepared_species_field, dict)
        else ""
    )
    if prepared_species != spec.key:
        raise DiscoveryError("prepared bundle species does not match the requested analysis species")
    entries = list(selections)
    if not entries:
        raise DiscoveryError("select at least one prepared candidate for analysis")
    if len(entries) > MAX_ACTIVE_CANDIDATES:
        raise DiscoveryError(f"at most {MAX_ACTIVE_CANDIDATES} candidates can be activated in one run")
    candidate_index = _candidate_index(prepared)
    _validate_prepared_source_units(prepared.get("studies", []))
    bundle_root_text = str(prepared.get("materialize_dir") or "")
    if not bundle_root_text:
        raise DiscoveryError("prepared bundle has no materialized file directory")
    bundle_root = Path(bundle_root_text).resolve()
    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()) and not force:
        raise FileExistsError(f"analysis output already exists and is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    derived_dir = output / "derived_tables"
    prepared_path = output / "prepared_bundle.json"
    request_path = output / "analysis_request.json"
    replay_args: list[str | Path] = [
        "degora",
        "discovery-analyze",
        prepared_path,
        request_path,
        "--species",
        spec.key,
        "--output-dir",
        output,
        "--min-studies",
        str(min_studies),
    ]
    if force:
        replay_args.append("--force")
    replay_command = shell_command(replay_args)
    prepared_path.write_text(json.dumps(prepared, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    rows: list[dict[str, Any]] = []
    author_derivations: list[dict[str, Any]] = []
    fallback_summaries: list[dict[str, Any]] = []
    seen_activations: set[tuple[str, ...]] = set()
    for sequence, entry in enumerate(entries, start=1):
        candidate_id = _text(entry.get("candidate_id"), field="candidate_id", required=True, maximum=64)
        pair = candidate_index.get(candidate_id)
        if pair is None:
            raise DiscoveryError(f"candidate does not belong to this prepared bundle: {candidate_id}")
        study, candidate = pair
        if str(study.get("species") or "") != spec.key:
            raise DiscoveryError("cross-species candidate detected in a species-specific bundle")
        _validate_mixed_activation(study)
        mode = _text(entry.get("mode"), field="mode", required=True, maximum=24)
        if mode == "author":
            activation_key = _author_activation_key(candidate_id, candidate.get("inspection", {}), entry)
            if activation_key in seen_activations:
                raise DiscoveryError(
                    "the same author candidate extraction was selected more than once; use a distinct exact row filter"
                )
            seen_activations.add(activation_key)
            row, summary = _author_row(
                study=study,
                candidate=candidate,
                entry=entry,
                spec=spec,
                bundle_root=bundle_root,
                derived_dir=derived_dir,
                sequence=sequence,
                replay_command=replay_command,
            )
            author_derivations.append(summary)
        elif mode == "fallback":
            activation_key = (candidate_id, "fallback")
            if activation_key in seen_activations:
                raise DiscoveryError(f"candidate selected more than once: {candidate_id}")
            seen_activations.add(activation_key)
            row, summary = _fallback_row(
                study=study,
                candidate=candidate,
                entry=entry,
                spec=spec,
                bundle_root=bundle_root,
                derived_dir=derived_dir,
                sequence=sequence,
                replay_command=replay_command,
            )
            fallback_summaries.append(summary)
        else:
            raise DiscoveryError("mode must be author or fallback")
        rows.append(row)

    source_units = sorted({str(row["source_unit_id"]) for row in rows})
    if len(source_units) < min_studies:
        raise DiscoveryError(
            f"DEGORA requires at least {min_studies} independent selected source units; found {len(source_units)}"
        )
    request_path.write_text(
        json.dumps(
            {
                "species": spec.key,
                "cross_species_pooling": False,
                "min_studies": min_studies,
                "source_units": source_units,
                "selections": entries,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    catalog_path = output / f"DEGORA_{spec.key}_selected_catalog.csv"
    catalog = pd.DataFrame(rows, columns=CATALOG_COLUMNS)
    neutralize_formula_text(catalog).to_csv(catalog_path, index=False, lineterminator="\n")
    write_source_sidecar(
        catalog_path,
        replay_command,
        inputs=[request_path],
        metadata={"generator": "discovery_selected_catalog", **formula_guard_metadata()},
    )
    validation = validate_catalog_inputs(catalog_path)
    results_dir = output / "results"
    harmonized_dir = output / "harmonized"
    metrics = run_slice(catalog_path, results_dir, harmonized_dir, min_studies=min_studies)
    harmonized_path = results_dir / "slice_harmonized.csv"
    db_path = results_dir / f"degora_{spec.key}_scores.db"
    score_summary = write_score_database(
        harmonized_path,
        results_dir,
        catalog_path=catalog_path,
        db_path=db_path,
        min_studies=min_studies,
        command=replay_command,
        extra_metadata={
            "discovery_species": spec.key,
            "discovery_cross_species_pooling": "false",
            "discovery_source_units": ",".join(source_units),
            **{str(key): str(value) for key, value in (extra_metadata or {}).items()},
        },
    )
    # A run that scored nothing is a failure, not a run with an empty table. The
    # standard CLI already refuses it; without the same refusal here a discovery
    # run registered status "complete" with top_genes [] and a workbook nobody
    # could read anything out of, and the enclosing transaction kept that partial
    # output as a finished run.
    if int(score_summary.get("n_gene_scores", 0) or 0) == 0:
        raise DiscoveryError(
            f"DEGORA scored zero genes at min_studies={min_studies}. No gene had usable, "
            "directional evidence from enough independent source units. Check that the "
            "selected sources share a gene identifier space and that contrast direction "
            "and table scope are correct, then prepare and analyze again."
        )
    excel_workbook = export_run_workbook(
        results_dir,
        results_dir / DEFAULT_WORKBOOK_NAME,
        config_path=catalog_path,
        db_path=db_path,
        command=replay_command,
    )
    (output / DISCOVERY_RUN_MARKER).write_text(
        json.dumps(
            {
                "artifact_type": DISCOVERY_RUN_ARTIFACT_TYPE,
                "format_version": DISCOVERY_RUN_FORMAT_VERSION,
                "species": spec.key,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "status": "complete",
        "species": {"key": spec.key, "label": spec.label, "scientific_name": spec.scientific_name},
        "output_dir": str(output),
        "catalog_path": str(catalog_path),
        "analysis_request": str(request_path),
        "results_dir": str(results_dir),
        "db_path": str(db_path),
        "score_csv": str(score_summary["score_csv"]),
        "excel_workbook": excel_workbook,
        "top_genes": score_summary.get("top_genes", []),
        "source_units": source_units,
        "n_source_units": len(source_units),
        "n_active_contrasts": len(rows),
        "validation": validation,
        "metrics": metrics,
        "author_derivations": author_derivations,
        "fallback_derivations": fallback_summaries,
        "cross_species_pooling": False,
    }


def run_discovery_analysis(
    prepared: dict[str, Any],
    selections: Iterable[dict[str, Any]],
    output_dir: str | Path,
    *,
    species: str,
    min_studies: int = 2,
    force: bool = False,
    extra_metadata: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Run a species-specific activation with rollback on every failed attempt."""

    output = Path(output_dir).resolve()
    existed_empty, backup = _begin_output_transaction(output, force=force)
    try:
        result = _execute_discovery_analysis(
            prepared,
            selections,
            output,
            species=species,
            min_studies=min_studies,
            force=force,
            extra_metadata=extra_metadata,
        )
    except BaseException:
        _rollback_output_transaction(output, existed_empty=existed_empty, backup=backup)
        raise
    _commit_output_transaction(backup)
    return result


__all__ = ["run_discovery_analysis"]
