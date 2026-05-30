"""Export manuscript-facing DEGORA score resource artifacts."""

from __future__ import annotations

import json
import sqlite3
import urllib.request
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pandas as pd

from .provenance import shell_command, write_source_sidecar


TOP_GENE_COLUMNS = [
    "rank_label",
    "gene_symbol",
    "evidence_tier",
    "degora_score",
    "priority_score",
    "priority_top_percent",
    "evidence_reliability_score",
    "direction_confidence_index",
    "loo_rank_stability_score",
    "top_percent_label",
    "support_label",
    "direction_label",
    "weighted_lfc",
    "high_confidence",
]

SOURCE_FAMILY_COLUMNS = [
    "gene_symbol",
    "n_studies",
    "n_source_units",
    "support_label",
    "source_units",
    "raw_evidence_row_count",
    "direction_label",
]


def validate_columns(frame: pd.DataFrame, required: list[str], *, label: str) -> None:
    """Fail fast when a manuscript-facing artifact would miss required fields."""

    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def top_gene_resource_table(score_csv: str | Path, *, top_n: int = 20) -> pd.DataFrame:
    """Return the human-readable top-gene DEGORA score table."""

    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    scores = pd.read_csv(score_csv)
    validate_columns(scores, TOP_GENE_COLUMNS, label="score CSV")
    return scores.loc[:, TOP_GENE_COLUMNS].head(top_n).copy()


def source_family_collapse_table(
    score_csv: str | Path,
    db_path: str | Path,
    *,
    top_n: int = 20,
) -> pd.DataFrame:
    """Return independent source-unit support with raw evidence rows separated."""

    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    scores = pd.read_csv(score_csv)
    validate_columns(scores, ["gene_symbol", "support_label", "direction_label"], label="score CSV")
    top_scores = scores.head(top_n).copy()
    genes = top_scores["gene_symbol"].astype(str).str.upper().tolist()
    if not genes:
        return pd.DataFrame(columns=SOURCE_FAMILY_COLUMNS)
    placeholders = ",".join("?" for _ in genes)

    with sqlite3.connect(db_path) as connection:
        available = {
            row[1]
            for row in connection.execute("PRAGMA table_info(gene_evidence)").fetchall()
        }
        n_contrast_expr = "n_contrast_rows" if "n_contrast_rows" in available else "1 AS n_contrast_rows"
        n_studies_expr = (
            "n_studies_in_source_unit"
            if "n_studies_in_source_unit" in available
            else "1 AS n_studies_in_source_unit"
        )
        evidence = pd.read_sql_query(
            f"""
            SELECT gene_symbol, study_id, source_unit_id, {n_contrast_expr}, {n_studies_expr}
            FROM gene_evidence
            WHERE gene_symbol IN ({placeholders})
            """,
            connection,
            params=genes,
        )

    if evidence.empty:
        grouped = pd.DataFrame(columns=["gene_symbol", "n_studies", "n_source_units", "source_units", "raw_evidence_row_count"])
    else:
        evidence["gene_symbol"] = evidence["gene_symbol"].astype(str).str.upper()
        grouped = evidence.groupby("gene_symbol", as_index=False).agg(
            n_studies=("n_studies_in_source_unit", "sum"),
            n_source_units=("source_unit_id", "nunique"),
            source_units=("source_unit_id", lambda values: ";".join(sorted(set(map(str, values))))),
            raw_evidence_row_count=("n_contrast_rows", "sum"),
        )

    collapse = top_scores[["gene_symbol", "support_label", "direction_label"]].merge(grouped, on="gene_symbol", how="left")
    for column in ("n_studies", "n_source_units", "raw_evidence_row_count"):
        collapse[column] = collapse[column].fillna(0).astype(int)
    collapse["source_units"] = collapse["source_units"].fillna("")
    collapse = collapse[SOURCE_FAMILY_COLUMNS]
    validate_columns(collapse, SOURCE_FAMILY_COLUMNS, label="source-family collapse table")
    return collapse


def capture_gene_api_json(
    api_url: str,
    output_path: str | Path,
    *,
    db_path: str | Path,
    command: str,
) -> dict[str, Any]:
    """Capture one local gene API response as a provenance-tracked JSON artifact."""

    parsed = urlparse(api_url)
    host = (parsed.hostname or "").lower()
    if parsed.scheme not in {"http", "https"} or host not in {"127.0.0.1", "localhost", "::1"}:
        raise ValueError("api_url must point to the local DEGORA API on localhost or 127.0.0.1")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(api_url, timeout=10) as response:
        payload = json.loads(response.read().decode("utf-8"))
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(
        output,
        command,
        inputs=[db_path],
        metadata={"generator": "degora-resource-package", "api_url": api_url},
    )
    return payload


def _api_gene_slug(api_url: str) -> str:
    parsed = urlparse(api_url)
    slug = Path(parsed.path).name.strip() or "response"
    safe = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in slug)
    return safe or "response"


def export_score_resource_package(
    score_csv: str | Path,
    db_path: str | Path,
    output_dir: str | Path,
    *,
    top_n: int = 20,
    api_url: str | None = None,
    command: str | None = None,
) -> dict[str, Any]:
    """Write the DEGORA manuscript/resource tables and optional API evidence."""

    score_csv = Path(score_csv)
    db_path = Path(db_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    command = command or shell_command(
        [
            "make",
            "-C",
            "outputs/code",
            "resource-package",
            f"RESOURCE_SCORE_CSV={score_csv.resolve()}",
            f"RESOURCE_SCORE_DB={db_path.resolve()}",
            f"RESOURCE_OUTDIR={output_dir.resolve()}",
            f"RESOURCE_TOP_N={top_n}",
            f"RESOURCE_API_URL={api_url or ''}",
        ]
    )

    top_table = top_gene_resource_table(score_csv, top_n=top_n)
    collapse_table = source_family_collapse_table(score_csv, db_path, top_n=top_n)

    top_path = output_dir / "degora_score_resource_top_genes.tsv"
    collapse_path = output_dir / "degora_score_source_family_collapse_top_genes.tsv"
    top_table.to_csv(top_path, sep="\t", index=False)
    collapse_table.to_csv(collapse_path, sep="\t", index=False)

    inputs = [score_csv, db_path]
    write_source_sidecar(top_path, command, inputs=inputs, metadata={"generator": "degora-resource-package", "top_n": top_n})
    write_source_sidecar(
        collapse_path,
        command,
        inputs=inputs,
        metadata={"generator": "degora-resource-package", "top_n": top_n, "raw_evidence_row_count": "gene_evidence rows before source-unit collapse"},
    )

    summary: dict[str, Any] = {
        "top_gene_table": str(top_path.resolve()),
        "source_family_collapse_table": str(collapse_path.resolve()),
        "top_n": int(top_n),
        "top_genes": top_table["gene_symbol"].tolist(),
    }

    if api_url:
        api_path = output_dir / f"api_genes_{_api_gene_slug(api_url)}.json"
        payload = capture_gene_api_json(api_url, api_path, db_path=db_path, command=command)
        summary["api_gene_json"] = str(api_path.resolve())
        summary["api_gene_symbol"] = payload.get("gene", {}).get("gene_symbol")

    summary_path = output_dir / "degora_score_resource_package_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(summary_path, command, inputs=inputs, metadata={"generator": "degora-resource-package"})
    return summary
