#!/usr/bin/env python
"""Write source-support-aware DEGORA top-gene reports."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from degora.provenance import shell_command, write_source_sidecar


SUMMARY_COLUMNS = [
    "metric",
    "cutoff",
    "denominator",
    "n_recovered",
    "recall",
    "recovered_genes",
    "missing_genes",
    "criteria",
]

TOP_GENE_COLUMNS = [
    "rank",
    "gene_symbol",
    "locked_gold",
    "interpretive_marker",
    "known_or_locked_marker",
    "expected_direction",
    "consensus_direction",
    "direction_match",
    "n_source_units",
    "sign_concordance",
    "source_support_pass",
    "degora_score",
    "priority_score",
    "evidence_reliability_score",
    "direction_confidence_index",
    "quality_weighted_degora_score",
    "source_units",
    "marker_role",
    "validation_role",
    "evidence_basis",
    "source_url",
    "notes",
]

SCORE_OPTIONAL_COLUMNS = [
    "degora_score",
    "priority_score",
    "evidence_reliability_score",
    "direction_confidence_index",
    "quality_weighted_degora_score",
]


def _normalize_symbol(value: object) -> str:
    return str(value).strip().upper()


def _normalize_direction(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"up", "+", "1", "increase", "increased", "upregulated", "up-regulated"}:
        return "up"
    if text in {"down", "-", "-1", "decrease", "decreased", "downregulated", "down-regulated"}:
        return "down"
    return ""


def _read_gold(path: Path, gene_column: str) -> tuple[set[str], dict[str, str]]:
    frame = pd.read_csv(path)
    if gene_column not in frame.columns:
        raise ValueError(f"gold file {path} missing gene column {gene_column!r}")
    genes = frame[gene_column].dropna().map(_normalize_symbol)
    positives = {gene for gene in genes if gene}
    directions: dict[str, str] = {}
    if "expected_direction" in frame.columns:
        for row in frame.to_dict(orient="records"):
            gene = _normalize_symbol(row.get(gene_column, ""))
            direction = _normalize_direction(row.get("expected_direction", ""))
            if gene in positives and direction:
                directions[gene] = direction
    return positives, directions


def _read_marker_panel(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    frame = pd.read_csv(path)
    if "gene_symbol" not in frame.columns:
        raise ValueError(f"marker panel {path} missing gene_symbol column")
    records: dict[str, dict[str, str]] = {}
    for row in frame.fillna("").to_dict(orient="records"):
        gene = _normalize_symbol(row.get("gene_symbol", ""))
        if not gene or gene in records:
            continue
        records[gene] = {str(key): str(value).strip() for key, value in row.items()}
    return records


def _ranked_scores(score_csv: Path, rank_column: str) -> pd.DataFrame:
    frame = pd.read_csv(score_csv)
    if "gene_symbol" not in frame.columns:
        raise ValueError(f"score file {score_csv} missing gene_symbol column")
    if rank_column not in frame.columns:
        raise ValueError(f"score file {score_csv} missing rank column {rank_column!r}")
    ranked = frame.copy()
    ranked["_symbol"] = ranked["gene_symbol"].map(_normalize_symbol)
    ranked[rank_column] = pd.to_numeric(ranked[rank_column], errors="coerce")
    ranked = ranked.dropna(subset=["_symbol", rank_column])
    ranked = ranked.loc[ranked["_symbol"].ne("")]
    return ranked.sort_values([rank_column, "_symbol"]).reset_index(drop=True)


def _num(row: pd.Series, column: str) -> float | None:
    if column not in row.index:
        return None
    value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
    if pd.isna(value):
        return None
    return float(value)


def _string(row: pd.Series, column: str) -> str:
    if column not in row.index or pd.isna(row[column]):
        return ""
    return str(row[column]).strip()


def _fmt_float(value: object) -> str:
    if value == "" or value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return f"{float(value):.3f}"


def _join(values: Iterable[str]) -> str:
    return ";".join(values)


def _source_supported(
    row: pd.Series,
    expected_direction: str,
    consensus_direction: str,
    *,
    min_source_units: int,
    min_sign_concordance: float,
) -> bool:
    if expected_direction and consensus_direction != expected_direction:
        return False
    n_source_units = _num(row, "n_source_units")
    sign_concordance = _num(row, "sign_concordance")
    if n_source_units is None or n_source_units < min_source_units:
        return False
    if sign_concordance is None or sign_concordance < min_sign_concordance:
        return False
    return True


def _top_gene_rows(
    ranked: pd.DataFrame,
    positives: set[str],
    expected_directions: dict[str, str],
    markers: dict[str, dict[str, str]],
    *,
    rank_column: str,
    direction_column: str,
    top_n: int,
    min_source_units: int,
    min_sign_concordance: float,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, series in ranked.head(top_n).iterrows():
        symbol = str(series["_symbol"])
        marker = markers.get(symbol, {})
        expected = expected_directions.get(symbol) or _normalize_direction(marker.get("expected_direction", ""))
        consensus = _normalize_direction(series.get(direction_column, "")) if direction_column in series.index else ""
        direction_match: bool | str
        direction_match = "" if not expected else consensus == expected
        source_support = (
            _source_supported(
                series,
                expected,
                consensus,
                min_source_units=min_source_units,
                min_sign_concordance=min_sign_concordance,
            )
            if expected
            else ""
        )
        record: dict[str, object] = {
            "rank": int(float(series[rank_column])),
            "gene_symbol": symbol,
            "locked_gold": symbol in positives,
            "interpretive_marker": symbol in markers,
            "known_or_locked_marker": symbol in positives or symbol in markers,
            "expected_direction": expected,
            "consensus_direction": consensus,
            "direction_match": direction_match,
            "n_source_units": _num(series, "n_source_units") if "n_source_units" in series.index else "",
            "sign_concordance": _num(series, "sign_concordance") if "sign_concordance" in series.index else "",
            "source_support_pass": source_support,
            "source_units": _string(series, "source_units"),
            "marker_role": marker.get("marker_role", ""),
            "validation_role": marker.get("validation_role", ""),
            "evidence_basis": marker.get("evidence_basis", ""),
            "source_url": marker.get("source_url", ""),
            "notes": marker.get("notes", ""),
        }
        for column in SCORE_OPTIONAL_COLUMNS:
            record[column] = _num(series, column) if column in series.index else ""
        rows.append(record)
    return pd.DataFrame(rows, columns=TOP_GENE_COLUMNS)


def _summary_rows(
    ranked: pd.DataFrame,
    positives: set[str],
    expected_directions: dict[str, str],
    markers: dict[str, dict[str, str]],
    *,
    rank_column: str,
    direction_column: str,
    cutoffs: list[int],
    min_source_units: int,
    min_sign_concordance: float,
) -> pd.DataFrame:
    all_known = positives.union(markers)
    rows: list[dict[str, object]] = []
    symbol_to_row = {str(row["_symbol"]): pd.Series(row) for row in ranked.to_dict(orient="records")}

    def add_row(metric: str, cutoff: int, denominator: set[str], recovered: set[str], criteria: str) -> None:
        missing = sorted(denominator.difference(recovered))
        rows.append(
            {
                "metric": metric,
                "cutoff": cutoff,
                "denominator": len(denominator),
                "n_recovered": len(recovered),
                "recall": len(recovered) / len(denominator) if denominator else 0.0,
                "recovered_genes": _join(sorted(recovered)),
                "missing_genes": _join(missing),
                "criteria": criteria,
            }
        )

    for cutoff in cutoffs:
        top_symbols = ranked.head(cutoff)["_symbol"].astype(str).tolist()
        top_set = set(top_symbols)
        locked_membership = positives.intersection(top_set)
        locked_direction = {
            gene
            for gene in locked_membership
            if not expected_directions.get(gene)
            or _normalize_direction(symbol_to_row[gene].get(direction_column, "")) == expected_directions[gene]
        }
        locked_source_supported = {
            gene
            for gene in locked_direction
            if _source_supported(
                symbol_to_row[gene],
                expected_directions.get(gene, ""),
                _normalize_direction(symbol_to_row[gene].get(direction_column, "")),
                min_source_units=min_source_units,
                min_sign_concordance=min_sign_concordance,
            )
        }
        add_row("locked_membership_recall", cutoff, positives, locked_membership, "gene is in locked gold panel")
        add_row(
            "locked_direction_recall",
            cutoff,
            positives,
            locked_direction,
            "locked gold gene is recovered with expected direction when declared",
        )
        add_row(
            "locked_source_supported_recall",
            cutoff,
            positives,
            locked_source_supported,
            f"locked direction match plus n_source_units>={min_source_units} and sign_concordance>={min_sign_concordance}",
        )
        if markers:
            known_membership = all_known.intersection(top_set)
            known_source_supported = set(locked_source_supported)
            for gene in markers:
                if gene not in top_set or gene not in symbol_to_row:
                    continue
                marker_expected = _normalize_direction(markers[gene].get("expected_direction", ""))
                consensus = _normalize_direction(symbol_to_row[gene].get(direction_column, ""))
                if marker_expected and consensus != marker_expected:
                    continue
                if _source_supported(
                    symbol_to_row[gene],
                    marker_expected,
                    consensus,
                    min_source_units=min_source_units,
                    min_sign_concordance=min_sign_concordance,
                ):
                    known_source_supported.add(gene)
            add_row(
                "known_or_locked_marker_coverage",
                cutoff,
                all_known,
                known_membership,
                "locked gold plus interpretive marker panel; not a locked validation metric",
            )
            add_row(
                "known_or_locked_marker_source_supported",
                cutoff,
                all_known,
                known_source_supported,
                "known/locked marker recovered with direction/source support; interpretive only",
            )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def build_report(
    score_csv: Path,
    gold_path: Path,
    *,
    marker_panel: Path | None = None,
    gold_gene_column: str = "gene_symbol",
    rank_column: str = "degora_rank",
    direction_column: str = "consensus_direction",
    top_n: int = 10,
    cutoffs: list[int] | None = None,
    min_source_units: int = 2,
    min_sign_concordance: float = 1.0,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    if min_source_units < 1:
        raise ValueError("min_source_units must be at least 1")
    if not 0.0 <= min_sign_concordance <= 1.0:
        raise ValueError("min_sign_concordance must be between 0 and 1")
    cutoffs = sorted(set(cutoffs or [10, 20, 50, 100]))
    positives, expected_directions = _read_gold(gold_path, gold_gene_column)
    if not positives:
        raise ValueError(f"gold file {gold_path} contains no positive genes")
    markers = _read_marker_panel(marker_panel)
    ranked = _ranked_scores(score_csv, rank_column)
    top_genes = _top_gene_rows(
        ranked,
        positives,
        expected_directions,
        markers,
        rank_column=rank_column,
        direction_column=direction_column,
        top_n=top_n,
        min_source_units=min_source_units,
        min_sign_concordance=min_sign_concordance,
    )
    summary = _summary_rows(
        ranked,
        positives,
        expected_directions,
        markers,
        rank_column=rank_column,
        direction_column=direction_column,
        cutoffs=cutoffs,
        min_source_units=min_source_units,
        min_sign_concordance=min_sign_concordance,
    )
    metadata = {
        "score_csv": str(score_csv),
        "gold_path": str(gold_path),
        "marker_panel": str(marker_panel) if marker_panel else "",
        "rank_column": rank_column,
        "direction_column": direction_column,
        "top_n": top_n,
        "cutoffs": cutoffs,
        "min_source_units": min_source_units,
        "min_sign_concordance": min_sign_concordance,
        "n_scored_genes": int(len(ranked)),
        "n_locked_gold": int(len(positives)),
        "n_interpretive_markers": int(len(markers)),
    }
    return summary, top_genes, metadata


def write_markdown(summary: pd.DataFrame, top_genes: pd.DataFrame, output: Path, *, title: str, metadata: dict[str, Any]) -> None:
    lines = [
        f"# {title}",
        "",
        f"Score table: `{metadata['score_csv']}`.",
        f"Locked gold panel: `{metadata['gold_path']}`.",
    ]
    if metadata.get("marker_panel"):
        lines.append(f"Interpretive marker panel: `{metadata['marker_panel']}`.")
    lines.extend(
        [
            "",
            "Guardrail: locked recall rows are validation diagnostics. Interpretive marker rows are biology annotation only and must not be reported as pre-locked benchmark performance.",
            "",
            f"Source-support criterion: `n_source_units >= {metadata['min_source_units']}` and `sign_concordance >= {metadata['min_sign_concordance']}` with expected direction when declared.",
            "",
            "## Summary",
            "",
            "| metric | K | recovered / denominator | recall | recovered genes |",
            "| --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.metric} | {row.cutoff} | {row.n_recovered} / {row.denominator} | "
            f"{float(row.recall):.2f} | {row.recovered_genes} |"
        )
    lines.extend(
        [
            "",
            "## Top Genes",
            "",
            "| rank | gene | locked | marker | direction | sources | concordance | reliability | role |",
            "| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in top_genes.itertuples(index=False):
        locked = "yes" if row.locked_gold else "no"
        marker = "yes" if row.interpretive_marker else "no"
        direction = f"{row.consensus_direction or 'NA'}"
        if row.expected_direction:
            direction += f" / expected {row.expected_direction}"
        reliability = _fmt_float(row.evidence_reliability_score)
        lines.append(
            f"| {row.rank} | {row.gene_symbol} | {locked} | {marker} | {direction} | "
            f"{_fmt_float(row.n_source_units)} | {_fmt_float(row.sign_concordance)} | {reliability} | {row.marker_role} |"
        )
    output.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-csv", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--marker-panel", type=Path)
    parser.add_argument("--gold-gene-column", default="gene_symbol")
    parser.add_argument("--rank-column", default="degora_rank")
    parser.add_argument("--direction-column", default="consensus_direction")
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--cutoffs", type=int, nargs="+", default=[10, 20, 50, 100])
    parser.add_argument("--min-source-units", type=int, default=2)
    parser.add_argument("--min-sign-concordance", type=float, default=1.0)
    parser.add_argument("--title", default="DEGORA Source-Support Report")
    parser.add_argument("--output-summary-csv", type=Path, required=True)
    parser.add_argument("--output-top-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args(argv)

    summary, top_genes, metadata = build_report(
        args.score_csv,
        args.gold,
        marker_panel=args.marker_panel,
        gold_gene_column=args.gold_gene_column,
        rank_column=args.rank_column,
        direction_column=args.direction_column,
        top_n=args.top_n,
        cutoffs=args.cutoffs,
        min_source_units=args.min_source_units,
        min_sign_concordance=args.min_sign_concordance,
    )
    for path in (args.output_summary_csv, args.output_top_tsv, args.output_json, args.output_md):
        path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_summary_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    top_genes.to_csv(args.output_top_tsv, sep="\t", index=False)
    args.output_json.write_text(
        json.dumps(
            {
                "metadata": metadata,
                "summary": summary.to_dict(orient="records"),
                "top_genes": top_genes.to_dict(orient="records"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    write_markdown(summary, top_genes, args.output_md, title=args.title, metadata=metadata)
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_source_support_report.py",
            "--score-csv",
            args.score_csv,
            "--gold",
            args.gold,
            *(
                [
                    "--marker-panel",
                    args.marker_panel,
                ]
                if args.marker_panel
                else []
            ),
            "--gold-gene-column",
            args.gold_gene_column,
            "--rank-column",
            args.rank_column,
            "--direction-column",
            args.direction_column,
            "--top-n",
            args.top_n,
            "--cutoffs",
            *args.cutoffs,
            "--min-source-units",
            args.min_source_units,
            "--min-sign-concordance",
            args.min_sign_concordance,
            "--title",
            args.title,
            "--output-summary-csv",
            args.output_summary_csv,
            "--output-top-tsv",
            args.output_top_tsv,
            "--output-json",
            args.output_json,
            "--output-md",
            args.output_md,
        ]
    )
    inputs: list[Path] = [args.score_csv, args.gold]
    if args.marker_panel:
        inputs.append(args.marker_panel)
    for artifact in (args.output_summary_csv, args.output_top_tsv, args.output_json, args.output_md):
        write_source_sidecar(
            artifact,
            command,
            inputs=inputs,
            metadata={"generator": "source-support-report", **metadata},
        )
    print(json.dumps({"summary_rows": len(summary), "top_gene_rows": len(top_genes)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
