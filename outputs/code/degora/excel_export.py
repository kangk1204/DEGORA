"""Excel audit workbook export for DEGORA run outputs."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from .excel_io import read_config_sheet
from .provenance import write_source_sidecar


EXCEL_MAX_ROWS = 1_048_576
DEFAULT_WORKBOOK_NAME = "DEGORA_output.xlsx"


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


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return value if isinstance(value, dict) else {}


def _read_gold_from_config(config_path: Path | None) -> pd.DataFrame:
    if config_path is None or not config_path.exists() or config_path.suffix.lower() not in {".xlsx", ".xls"}:
        return pd.DataFrame()
    try:
        with pd.ExcelFile(config_path) as workbook:
            if "GoldPanel" not in workbook.sheet_names:
                return pd.DataFrame()
            gold = read_config_sheet(workbook, "GoldPanel")
    except Exception:
        return pd.DataFrame()
    if "gene_symbol" not in gold.columns:
        return pd.DataFrame()
    gold = gold.copy()
    if "locked" in gold.columns:
        locked = gold["locked"].astype("string").fillna("").str.strip().str.lower()
        keep = locked.isin({"", "1", "true", "t", "yes", "y", "locked"})
        gold = gold.loc[keep].copy()
    gold["gene_symbol"] = gold["gene_symbol"].astype(str).str.upper().str.strip()
    return gold.loc[gold["gene_symbol"].ne("")].drop_duplicates("gene_symbol").reset_index(drop=True)


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
    lookup = gold.merge(genes[present], on="gene_symbol", how="left")
    rank = pd.to_numeric(lookup.get("quality_weighted_degora_rank"), errors="coerce")
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
) -> pd.DataFrame:
    rank_col = "quality_weighted_degora_rank" if "quality_weighted_degora_rank" in genes.columns else "degora_rank"
    top = (
        genes.sort_values(rank_col).head(20)["gene_symbol"].tolist()
        if rank_col in genes.columns and "gene_symbol" in genes.columns
        else []
    )
    rows: list[dict[str, Any]] = [
        {"field": "result_dir", "value": str(result_dir.resolve())},
        {"field": "n_scored_genes", "value": int(len(genes))},
        {"field": "n_gene_evidence_rows", "value": int(len(evidence))},
        {
            "field": "n_source_units",
            "value": int(studies["source_unit_id"].nunique()) if "source_unit_id" in studies.columns else "",
        },
        {"field": "n_studies", "value": int(len(studies))},
        {"field": "n_curated_genes", "value": int(len(gold))},
        {"field": "top_genes", "value": ";".join(map(str, top))},
    ]
    if not gold.empty and "gene_symbol" in genes.columns:
        ranked = genes.sort_values(rank_col) if rank_col in genes.columns else genes
        curated = set(gold["gene_symbol"].astype(str))
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


def _autosize(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for cells in worksheet.columns:
            header = str(cells[0].value or "")
            worksheet.column_dimensions[cells[0].column_letter].width = min(max(len(header) + 2, 12), 48)


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
    evidence = _read_table(db_path, "gene_evidence")
    studies = _read_table(db_path, "studies")
    meta = _read_table(db_path, "meta")
    metadata = _read_json(metadata_json)
    diagnostics = pd.read_csv(diagnostics_tsv, sep="\t") if diagnostics_tsv.exists() else pd.DataFrame()
    gold = _read_gold_from_config(config_path)
    lookup = _curated_lookup(gold, genes)
    summary = _summary_rows(result_dir, genes, evidence, studies, gold)
    metadata_frame = _metadata_table(metadata)

    output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Run_summary", index=False)
        _write_sheet_chunks(writer, genes, "Gene_scores")
        _write_sheet_chunks(writer, evidence, "Gene_evidence")
        _write_sheet_chunks(writer, studies, "Source_units")
        _write_sheet_chunks(writer, lookup, "Curated_lookup")
        _write_sheet_chunks(writer, diagnostics, "Source_quality")
        _write_sheet_chunks(writer, metadata_frame, "Metadata")
        _write_sheet_chunks(writer, meta, "SQLite_meta")
        _autosize(writer)

    manifest = output.with_suffix(".manifest.json")
    validation = output.with_suffix(".validation.txt")
    inputs = [
        path
        for path in [db_path, score_csv, metadata_json, diagnostics_tsv, config_path]
        if path is not None and path.exists()
    ]
    manifest_data = {
        "generated_at": "deterministic",
        "script": "degora.excel_export.export_run_workbook",
        "command": command,
        "inputs": [path.as_posix() for path in inputs],
        "outputs": [output.as_posix(), manifest.as_posix(), validation.as_posix()],
        "sheets": {
            "Run_summary": "Run-level counts and top genes.",
            "Gene_scores": "Full DEGORA ranked gene table.",
            "Gene_evidence": "Source-unit evidence rows from the SQLite database.",
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
                f"n_gene_evidence_rows={len(evidence)}",
                f"n_source_units={studies['source_unit_id'].nunique() if 'source_unit_id' in studies.columns else 0}",
                f"n_curated_genes={len(gold)}",
            ]
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    for artifact in [output, manifest, validation]:
        write_source_sidecar(artifact, command, inputs=inputs, metadata={"generator": "degora-run-workbook"})
    return {
        "output": output.as_posix(),
        "manifest": manifest.as_posix(),
        "validation": validation.as_posix(),
        "rows_gene_scores": int(len(genes)),
        "rows_gene_evidence": int(len(evidence)),
        "rows_curated_lookup": int(len(lookup)),
    }
