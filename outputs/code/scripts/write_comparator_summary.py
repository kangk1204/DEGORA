#!/usr/bin/env python
"""Summarize baseline comparator recovery on the locked HIF target panel."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

from degora.gold import HIF1A_UP_TARGETS


SUMMARY_COLUMNS = [
    "method_id",
    "setting_id",
    "run_status",
    "n_rows",
    "recall_at_50",
    "recall_at_100",
    "recovered_at_50",
    "recovered_at_100",
    "missing_at_100",
    "top10",
    "failure_mode",
]


def _join(values: list[str]) -> str:
    return ";".join(values)


def _recall(frame: pd.DataFrame, k: int) -> tuple[float, list[str], list[str]]:
    positives = {gene.upper() for gene in HIF1A_UP_TARGETS}
    top = frame.sort_values("rank").head(k)["symbol"].astype(str).str.upper().tolist()
    recovered = sorted(positives.intersection(top))
    missing = sorted(positives.difference(top))
    return len(recovered) / len(positives), recovered, missing


def _baseline_rows(baseline_dir: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for path in sorted(baseline_dir.glob("hypoxia_*_*.tsv")):
        frame = pd.read_csv(path, sep="\t")
        if frame.empty:
            continue
        method_id = str(frame["method_id"].iloc[0])
        setting_id = str(frame["setting_id"].iloc[0])
        recall50, recovered50, _ = _recall(frame, 50)
        recall100, recovered100, missing100 = _recall(frame, 100)
        top10 = frame.sort_values("rank").head(10)["symbol"].astype(str).tolist()
        records.append(
            {
                "method_id": method_id,
                "setting_id": setting_id,
                "run_status": "ok",
                "n_rows": int(len(frame)),
                "recall_at_50": recall50,
                "recall_at_100": recall100,
                "recovered_at_50": _join(recovered50),
                "recovered_at_100": _join(recovered100),
                "missing_at_100": _join(missing100),
                "top10": _join(top10),
                "failure_mode": "",
            }
        )
    return records


def _blocked_rows(baseline_dir: Path, emitted: set[tuple[str, str]]) -> list[dict[str, object]]:
    failure_path = baseline_dir / "baseline_failure_ledger.csv"
    if not failure_path.exists():
        return []
    failure = pd.read_csv(failure_path)
    records: list[dict[str, object]] = []
    for row in failure.itertuples(index=False):
        key = (str(row.method_id), str(row.setting_id))
        if key in emitted:
            continue
        records.append(
            {
                "method_id": str(row.method_id),
                "setting_id": str(row.setting_id),
                "run_status": "blocked",
                "n_rows": 0,
                "recall_at_50": "",
                "recall_at_100": "",
                "recovered_at_50": "",
                "recovered_at_100": "",
                "missing_at_100": "",
                "top10": "",
                "failure_mode": str(row.blocker_id),
            }
        )
    return records


def build_summary(baseline_dir: Path) -> pd.DataFrame:
    records = _baseline_rows(baseline_dir)
    emitted = {(str(row["method_id"]), str(row["setting_id"])) for row in records}
    records.extend(_blocked_rows(baseline_dir, emitted))
    frame = pd.DataFrame.from_records(records, columns=SUMMARY_COLUMNS)
    if frame.empty:
        return frame
    order = {
        "degora_slice": 0,
        "weighted_stouffer": 1,
        "unweighted_stouffer": 2,
        "fisher": 3,
        "rank_product_approx": 4,
        "sign_vote": 5,
        "metavolcanor": 6,
        "robustrankaggreg": 7,
        "hstouffer": 8,
        "awmeta": 9,
    }
    frame["_order"] = frame["method_id"].map(order).fillna(999)
    return frame.sort_values(["_order", "setting_id"]).drop(columns=["_order"]).reset_index(drop=True)


def _summary_title(baseline_dir: Path) -> str:
    for part in reversed(baseline_dir.parts):
        match = re.fullmatch(r"iter-(\d+)", part)
        if match:
            return f"Iteration {match.group(1)} Comparator Recall Summary"
    return "Comparator Recall Summary"


def write_markdown(summary: pd.DataFrame, output: Path, title: str) -> None:
    lines = [
        f"# {title}",
        "",
        "Locked positives: HIF1A_UP_TARGETS from `outputs/code/degora/gold.py`.",
        "",
        "| method | status | rows | recall@50 | recall@100 | failure | top10 |",
        "| --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in summary.itertuples(index=False):
        r50 = "" if row.recall_at_50 == "" else f"{float(row.recall_at_50):.2f}"
        r100 = "" if row.recall_at_100 == "" else f"{float(row.recall_at_100):.2f}"
        lines.append(
            f"| {row.method_id}/{row.setting_id} | {row.run_status} | {row.n_rows} | {r50} | {r100} | {row.failure_mode} | {row.top10} |"
        )
    lines.extend(
        [
            "",
            "Interpretation guardrail: comparative superiority claims remain disallowed while hStouffer and AWmeta are blocked.",
            "",
        ]
    )
    output.write_text("\n".join(lines))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=Path("../../outputs/results/iter-13/baselines"))
    parser.add_argument("--output-csv", type=Path, default=Path("../../outputs/results/iter-13/comparator_recall_summary.csv"))
    parser.add_argument("--output-md", type=Path, default=Path("../../outputs/results/iter-13/comparator_recall_summary.md"))
    args = parser.parse_args(argv)

    baseline_dir = args.baseline_dir.resolve()
    output_csv = args.output_csv.resolve()
    output_md = args.output_md.resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary = build_summary(baseline_dir)
    summary.to_csv(output_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    write_markdown(summary, output_md, _summary_title(baseline_dir))
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/write_comparator_summary.py "
        f"--baseline-dir {baseline_dir} --output-csv {output_csv} --output-md {output_md}"
    )
    output_csv.with_suffix(output_csv.suffix + ".source").write_text(command + "\n")
    output_md.with_suffix(output_md.suffix + ".source").write_text(command + "\n")
    print(summary.to_json(orient="records"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
