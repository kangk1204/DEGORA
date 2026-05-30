#!/usr/bin/env python
"""Write a concise ER-stress benchmark comparison report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.provenance import shell_command, write_source_sidecar


METHOD_LABELS = {
    "degora_deg_score": "DEGORA score",
    "degora_quality_weighted_score": "DEGORA quality-weighted score",
    "weighted_stouffer": "weighted Stouffer",
    "fisher": "Fisher",
    "metavolcanor": "MetaVolcanoR",
    "robustrankaggreg": "RobustRankAggreg",
}
DEFAULT_DEGORA_METHOD_ID = "degora_quality_weighted_score"


def _method_rows(path: Path) -> dict[str, dict[str, Any]]:
    frame = pd.read_csv(path)
    rows: dict[str, dict[str, Any]] = {}
    for row in frame.to_dict(orient="records"):
        method = str(row["method_id"])
        if method in METHOD_LABELS:
            rows[method] = row
    return rows


def _recall(row: dict[str, Any], k: int) -> float:
    value = row.get(f"recall_at_{k}", 0.0)
    if value == "":
        return 0.0
    return float(value)


def _direction_recall(row: dict[str, Any], k: int) -> float | str:
    value = row.get(f"direction_recall_at_{k}", "")
    if value == "":
        return ""
    return float(value)


def _text(row: dict[str, Any], key: str) -> str:
    value = row.get(key, "")
    if value == "" or pd.isna(value):
        return ""
    return str(value)


def _compact_rows(rows: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    compact = []
    for method, label in METHOD_LABELS.items():
        if method not in rows:
            continue
        row = rows[method]
        compact.append(
            {
                "method_id": method,
                "label": label,
                "status": row.get("run_status", ""),
                "recall_at_10": _recall(row, 10),
                "recall_at_20": _recall(row, 20),
                "recall_at_50": _recall(row, 50),
                "recall_at_100": _recall(row, 100),
                "direction_recall_at_50": _direction_recall(row, 50),
                "direction_recall_at_100": _direction_recall(row, 100),
                "top10": _text(row, "top10"),
                "recovered_at_100": _text(row, "recovered_at_100"),
                "direction_recovered_at_100": _text(row, "direction_recovered_at_100"),
                "direction_mismatched_at_100": _text(row, "direction_mismatched_at_100"),
            }
        )
    return compact


def _default_degora_row(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    if DEFAULT_DEGORA_METHOD_ID in rows:
        return rows[DEFAULT_DEGORA_METHOD_ID]
    if "degora_deg_score" in rows:
        return rows["degora_deg_score"]
    raise ValueError("summary is missing a DEGORA score row")


def _delta(degora: dict[str, Any], other: dict[str, Any], k: int) -> float:
    return round(_recall(degora, k) - _recall(other, k), 10)


def _tier_delta(primary_row: dict[str, Any], full_row: dict[str, Any], k: int) -> float:
    return round(_recall(primary_row, k) - _recall(full_row, k), 10)


def build_report(primary_summary: Path, full_summary: Path) -> dict[str, Any]:
    primary = _method_rows(primary_summary)
    full = _method_rows(full_summary)
    degora_primary = _default_degora_row(primary)
    degora_full = _default_degora_row(full)

    unweighted_full = full.get("degora_deg_score")
    quality_full = full.get("degora_quality_weighted_score")
    return {
        "benchmark": "ER stress / UPR",
        "default_degora_method_id": str(degora_primary.get("method_id", DEFAULT_DEGORA_METHOD_ID)),
        "claim": (
            "ER stress / UPR is reported with the predeclared quality-weighted DEGORA ranking, alongside "
            "primary and full sensitivity source sets. The claim is directional evidence integration and "
            "auditable cross-platform support, not cherry-picked superiority from a favorable subset."
        ),
        "primary_set": _compact_rows(primary),
        "full_sensitivity_set": _compact_rows(full),
        "primary_degora_vs_weighted_stouffer_delta": {
            f"recall_at_{k}": _delta(degora_primary, primary["weighted_stouffer"], k)
            for k in (10, 20, 50, 100)
            if "weighted_stouffer" in primary
        },
        "primary_degora_vs_fisher_delta": {
            f"recall_at_{k}": _delta(degora_primary, primary["fisher"], k)
            for k in (10, 20, 50, 100)
            if "fisher" in primary
        },
        "full_sensitivity_degora_vs_fisher_delta": {
            f"recall_at_{k}": _delta(degora_full, full["fisher"], k)
            for k in (10, 20, 50, 100)
            if "fisher" in full
        },
        "source_quality_tier_effect": {
            f"degora_recall_primary_minus_full_at_{k}": _tier_delta(degora_primary, degora_full, k)
            for k in (10, 20, 50, 100)
        },
        "quality_weighted_sensitivity_effect": {
            f"quality_weighted_minus_unweighted_score_full_at_{k}": _delta(quality_full, unweighted_full, k)
            for k in (10, 20, 50, 100)
        }
        if quality_full and unweighted_full
        else {},
        "interpretation": [
            "Use the quality-weighted DEGORA row as the fixed manuscript-facing ranking for both primary and sensitivity source sets.",
            "Report primary and full source sets together; differences are source-quality sensitivity evidence, not a basis for choosing the most favorable subset.",
            "Keep the unweighted DEGORA score as a reference index so the effect of source-quality weighting remains visible.",
            "The comparison claim is early-rank prioritization and evidence-DB usability, not universal statistical superiority.",
            "hStouffer and AWmeta remain faithful-input blockers, not defeated baselines.",
        ],
        "comparison_strategy": [
            "Pre-lock gold panel and source-quality tiers before scoring.",
            "Report recall@10, @20, @50, and @100; do not only report the most favorable cutoff.",
            "Use direction-aware recall when a gold panel declares expected directions; this exposes direction mismatches that simple membership recall hides.",
            "Use the source-support-aware summaries because DEGORA's intended DB use is directional consistency across independent source units.",
            "Report source-quality diagnostics because normalized-matrix sources can conflict with full-table sources and should not silently flip the primary evidence story.",
            "Keep primary-quality and full-sensitivity datasets separate instead of hiding lower-confidence sources.",
            "Preserve per-gene evidence rows so top hits can be audited source by source.",
        ],
    }


def write_markdown(report: dict[str, Any], output: Path) -> None:
    lines = [
        "# ER Stress / UPR Benchmark Report",
        "",
        report["claim"],
        "",
        "## Primary Set",
        "",
        "| method | status | recall@10 | recall@20 | recall@50 | recall@100 | dir@50 | dir@100 | top10 |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in report["primary_set"]:
        dir50 = "" if row["direction_recall_at_50"] == "" else f"{row['direction_recall_at_50']:.2f}"
        dir100 = "" if row["direction_recall_at_100"] == "" else f"{row['direction_recall_at_100']:.2f}"
        lines.append(
            f"| {row['label']} | {row['status']} | {row['recall_at_10']:.2f} | {row['recall_at_20']:.2f} | "
            f"{row['recall_at_50']:.2f} | {row['recall_at_100']:.2f} | {dir50} | {dir100} | {row['top10']} |"
        )
    lines.extend(
        [
            "",
            "## Full Sensitivity",
            "",
            "| method | status | recall@10 | recall@20 | recall@50 | recall@100 | dir@50 | dir@100 | top10 |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report["full_sensitivity_set"]:
        dir50 = "" if row["direction_recall_at_50"] == "" else f"{row['direction_recall_at_50']:.2f}"
        dir100 = "" if row["direction_recall_at_100"] == "" else f"{row['direction_recall_at_100']:.2f}"
        lines.append(
            f"| {row['label']} | {row['status']} | {row['recall_at_10']:.2f} | {row['recall_at_20']:.2f} | "
            f"{row['recall_at_50']:.2f} | {row['recall_at_100']:.2f} | {dir50} | {dir100} | {row['top10']} |"
        )
    lines.extend(["", "## Source-Quality Tier Effect", ""])
    for key, value in report["source_quality_tier_effect"].items():
        lines.append(f"- `{key}`: {value:.2f}")
    if report.get("quality_weighted_sensitivity_effect"):
        lines.extend(["", "## Quality-Weighted Secondary Score Effect", ""])
        for key, value in report["quality_weighted_sensitivity_effect"].items():
            lines.append(f"- `{key}`: {value:.2f}")
    lines.extend(["", "## Interpretation", ""])
    lines.extend(f"- {item}" for item in report["interpretation"])
    lines.extend(["", "## Stronger Comparison Strategy", ""])
    lines.extend(f"- {item}" for item in report["comparison_strategy"])
    output.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-summary", type=Path, required=True)
    parser.add_argument("--full-summary", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args(argv)

    report = build_report(args.primary_summary, args.full_summary)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    write_markdown(report, args.output_md)
    command = shell_command(
        [
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/write_er_stress_benchmark_report.py",
            "--primary-summary",
            args.primary_summary,
            "--full-summary",
            args.full_summary,
            "--output-json",
            args.output_json,
            "--output-md",
            args.output_md,
        ]
    )
    for artifact in (args.output_json, args.output_md):
        write_source_sidecar(
            artifact,
            command,
            inputs=[args.primary_summary, args.full_summary],
            metadata={"generator": "er-stress-upr-benchmark-report"},
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
