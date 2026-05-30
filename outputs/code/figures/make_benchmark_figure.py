#!/usr/bin/env python
"""Create Figure 2: benchmark comparison across locked mechanism panels."""

from __future__ import annotations

import argparse
import json
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

from degora.provenance import shell_command, write_source_sidecar


FIGURE_ID = "FIGURE_2_BENCHMARK"
DEFAULT_OUTPUT = Path("outputs/figures/manuscript/FIGURE_2_BENCHMARK")

BENCHMARKS = [
    ("IFN", "ifn-pilot", Path("outputs/results/ifn-pilot/ifn_gold_comparator_summary.csv"), Path("outputs/results/ifn-pilot/deep-metrics/ifn-pilot_degora_advantage_metrics.csv")),
    ("ER stress", "er-stress-cross-platform", Path("outputs/results/er-stress-cross-platform/er_stress_cross_platform_gold_comparator_summary.csv"), Path("outputs/results/er-stress-cross-platform/deep-metrics/er-stress-cross-platform_degora_advantage_metrics.csv")),
    ("Heat shock", "heat-shock-benchmark", Path("outputs/results/heat-shock-benchmark/heat_shock_gold_comparator_summary.csv"), Path("outputs/results/heat-shock-benchmark/deep-metrics/heat-shock-benchmark_degora_advantage_metrics.csv")),
    ("Hypoxia", "hypoxia-hif1-benchmark", Path("outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv"), Path("outputs/results/hypoxia-hif1-benchmark/deep-metrics/hypoxia-hif1-benchmark_degora_advantage_metrics.csv")),
]

METHOD_LABELS = {
    "degora_deg_score": "DEGORA",
    "degora_quality_weighted_score": "DEGORA quality",
    "weighted_stouffer": "Weighted Stouffer",
    "fisher": "Fisher",
    "awfisher": "AWFisher",
    "metarnaseq_fisher": "metaRNASeq Fisher",
    "metavolcanor": "MetaVolcanoR",
    "robustrankaggreg": "RobustRankAggreg",
    "rank_product_approx": "Rank product approx.",
}

METHOD_ORDER = [
    "DEGORA",
    "DEGORA quality",
    "Weighted Stouffer",
    "Fisher",
    "AWFisher",
    "metaRNASeq Fisher",
    "MetaVolcanoR",
    "RobustRankAggreg",
    "Rank product approx.",
]


def _style() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10, "svg.fonttype": "none"})


def _read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames = []
    evidence_frames = []
    for label, corpus, summary_path, evidence_path in BENCHMARKS:
        frame = pd.read_csv(summary_path)
        frame["benchmark"] = label
        frame["corpus"] = corpus
        frame["method_label"] = frame["method_id"].map(METHOD_LABELS).fillna(frame["method_id"].astype(str))
        for col in ["recall_at_10", "recall_at_20", "recall_at_50", "recall_at_100", "direction_recall_at_100"]:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
        frames.append(frame)
        evidence = pd.read_csv(evidence_path)
        evidence["benchmark"] = label
        evidence["corpus"] = corpus
        evidence_frames.append(evidence)
    recall = pd.concat(frames, ignore_index=True)
    evidence = pd.concat(evidence_frames, ignore_index=True)
    fixed_rows = []
    for benchmark, group in recall.loc[recall["run_status"].eq("ok")].groupby("benchmark", sort=False):
        degora = group.loc[group["method_id"].eq("degora_quality_weighted_score")]
        if degora.empty:
            degora = group.loc[group["method_id"].eq("degora_deg_score")]
        if not degora.empty:
            fixed_rows.append(degora.iloc[0])
    degora_best = pd.DataFrame(fixed_rows).reset_index(drop=True)
    return recall, evidence, degora_best


def _recall_long(degora_best: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in degora_best.itertuples(index=False):
        for k in [10, 20, 50, 100]:
            rows.append({"benchmark": row.benchmark, "k": k, "recall": float(getattr(row, f"recall_at_{k}")), "top10": row.top10})
    return pd.DataFrame(rows)


def _evidence_top100(evidence: pd.DataFrame) -> pd.DataFrame:
    subset = evidence.loc[evidence["subset"].eq("top100_quality_weighted")].copy()
    cols = [
        "benchmark",
        "median_n_source_units",
        "median_sign_concordance",
        "fraction_full_concordance",
        "median_loo_rank_stability_score",
        "median_source_quality_weight_sum",
    ]
    return subset[cols].reset_index(drop=True)


def _heatmap_data(recall: pd.DataFrame) -> pd.DataFrame:
    frame = recall.loc[recall["run_status"].eq("ok") & recall["method_label"].isin(METHOD_ORDER)].copy()
    frame["method_label"] = pd.Categorical(frame["method_label"], METHOD_ORDER, ordered=True)
    return frame.sort_values(["method_label", "benchmark"]).reset_index(drop=True)


def _save(fig: plt.Figure, stem: Path) -> list[Path]:
    outputs = []
    for suffix in [".png", ".pdf", ".svg"]:
        out = stem.with_suffix(suffix)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        outputs.append(out)
    plt.close(fig)
    return outputs


def _plot_panel_a(recall_long: pd.DataFrame, stem: Path) -> list[Path]:
    _style()
    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    for benchmark, group in recall_long.groupby("benchmark", sort=False):
        ax.plot(group["k"], group["recall"], marker="o", linewidth=2.2, label=benchmark)
    ax.set_title("A. DEGORA recall by cutoff", loc="left", fontweight="bold")
    ax.set_xlabel("Top K genes")
    ax.set_ylabel("Locked-panel recall")
    ax.set_xticks([10, 20, 50, 100])
    ax.set_ylim(0, 1.02)
    ax.grid(True, color="#e5e7eb")
    ax.legend(frameon=False, loc="lower right")
    return _save(fig, stem)


def _plot_panel_b(heatmap: pd.DataFrame, stem: Path) -> list[Path]:
    _style()
    benchmarks = [item[0] for item in BENCHMARKS]
    methods = [m for m in METHOD_ORDER if m in set(heatmap["method_label"].astype(str))]
    matrix = np.full((len(methods), len(benchmarks)), np.nan)
    for i, method in enumerate(methods):
        for j, benchmark in enumerate(benchmarks):
            subset = heatmap.loc[heatmap["method_label"].astype(str).eq(method) & heatmap["benchmark"].eq(benchmark)]
            if not subset.empty:
                matrix[i, j] = float(subset["recall_at_100"].iloc[0])
    fig, ax = plt.subplots(figsize=(7.8, 3.9))
    cmap = LinearSegmentedColormap.from_list("recall", ["#f8fafc", "#a7d8c9", "#0f766e"])
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(benchmarks)))
    ax.set_xticklabels(benchmarks, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(benchmarks)):
            value = matrix[i, j]
            ax.text(j, i, "" if np.isnan(value) else f"{value:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("B. Recall@100 versus summary-level comparators", loc="left", fontweight="bold")
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02, label="Recall@100")
    return _save(fig, stem)


def _plot_panel_c(evidence: pd.DataFrame, stem: Path) -> list[Path]:
    _style()
    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    x = np.arange(len(evidence))
    ax.bar(x - 0.18, evidence["median_n_source_units"], width=0.36, color="#265d97", label="Median source units")
    ax2 = ax.twinx()
    ax2.bar(x + 0.18, evidence["median_sign_concordance"], width=0.36, color="#0f766e", label="Median sign concordance")
    ax.set_xticks(x)
    ax.set_xticklabels(evidence["benchmark"], rotation=25, ha="right")
    ax.set_ylabel("Source units")
    ax2.set_ylabel("Sign concordance")
    ax2.set_ylim(0, 1.05)
    ax.set_title("C. Top-100 DEGORA evidence support", loc="left", fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles + handles2, labels + labels2, frameon=False, loc="upper left")
    return _save(fig, stem)


def _plot_panel_d(degora_best: pd.DataFrame, stem: Path) -> list[Path]:
    _style()
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    labels = degora_best["benchmark"].tolist()
    y = np.arange(len(labels))
    ax.barh(y, degora_best["recall_at_10"], color="#dbeafe", label="Recall@10")
    ax.barh(y, degora_best["recall_at_50"] - degora_best["recall_at_10"], left=degora_best["recall_at_10"], color="#93c5fd", label="Gain to @50")
    ax.barh(y, degora_best["recall_at_100"] - degora_best["recall_at_50"], left=degora_best["recall_at_50"], color="#265d97", label="Gain to @100")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Recall")
    ax.set_title("D. Early versus broader prioritization", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb")
    ax.legend(frameon=False, loc="lower right")
    return _save(fig, stem)


def _plot_combined(recall_long: pd.DataFrame, heatmap: pd.DataFrame, evidence: pd.DataFrame, degora_best: pd.DataFrame, stem: Path) -> list[Path]:
    _style()
    fig = plt.figure(figsize=(11.2, 8.2))
    gs = fig.add_gridspec(2, 2, wspace=0.35, hspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    for benchmark, group in recall_long.groupby("benchmark", sort=False):
        ax_a.plot(group["k"], group["recall"], marker="o", linewidth=2.0, label=benchmark)
    ax_a.set_title("A. DEGORA recall", loc="left", fontweight="bold")
    ax_a.set_xlabel("Top K genes")
    ax_a.set_ylabel("Recall")
    ax_a.set_xticks([10, 20, 50, 100])
    ax_a.set_ylim(0, 1.02)
    ax_a.grid(True, color="#e5e7eb")
    ax_a.legend(frameon=False, fontsize=7)

    ax_b = fig.add_subplot(gs[0, 1])
    benchmarks = [item[0] for item in BENCHMARKS]
    methods = [m for m in METHOD_ORDER if m in set(heatmap["method_label"].astype(str))]
    matrix = np.full((len(methods), len(benchmarks)), np.nan)
    for i, method in enumerate(methods):
        for j, benchmark in enumerate(benchmarks):
            subset = heatmap.loc[heatmap["method_label"].astype(str).eq(method) & heatmap["benchmark"].eq(benchmark)]
            if not subset.empty:
                matrix[i, j] = float(subset["recall_at_100"].iloc[0])
    cmap = LinearSegmentedColormap.from_list("recall", ["#f8fafc", "#a7d8c9", "#0f766e"])
    ax_b.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax_b.set_xticks(np.arange(len(benchmarks)))
    ax_b.set_xticklabels(benchmarks, rotation=25, ha="right")
    ax_b.set_yticks(np.arange(len(methods)))
    ax_b.set_yticklabels(methods, fontsize=7)
    for i in range(len(methods)):
        for j in range(len(benchmarks)):
            value = matrix[i, j]
            ax_b.text(j, i, "" if np.isnan(value) else f"{value:.2f}", ha="center", va="center", fontsize=6.5)
    ax_b.set_title("B. Recall@100", loc="left", fontweight="bold")

    ax_c = fig.add_subplot(gs[1, 0])
    x = np.arange(len(evidence))
    ax_c.bar(x, evidence["median_n_source_units"], color="#265d97")
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(evidence["benchmark"], rotation=25, ha="right")
    ax_c.set_ylabel("Median source units in top100")
    ax_c.set_title("C. Evidence support", loc="left", fontweight="bold")
    ax_c.grid(axis="y", color="#e5e7eb")

    ax_d = fig.add_subplot(gs[1, 1])
    y = np.arange(len(degora_best))
    ax_d.barh(y, degora_best["recall_at_10"], color="#dbeafe", label="@10")
    ax_d.barh(y, degora_best["recall_at_50"] - degora_best["recall_at_10"], left=degora_best["recall_at_10"], color="#93c5fd", label="@50 gain")
    ax_d.barh(y, degora_best["recall_at_100"] - degora_best["recall_at_50"], left=degora_best["recall_at_50"], color="#265d97", label="@100 gain")
    ax_d.set_yticks(y)
    ax_d.set_yticklabels(degora_best["benchmark"])
    ax_d.set_xlim(0, 1)
    ax_d.set_title("D. Recall gains by cutoff", loc="left", fontweight="bold")
    ax_d.legend(frameon=False, fontsize=7)
    ax_d.grid(axis="x", color="#e5e7eb")
    fig.suptitle("DEGORA benchmark against summary-level comparators", x=0.01, y=0.995, ha="left", fontsize=12, fontweight="bold")
    return _save(fig, stem)


def _write_xlsx(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    wb = Workbook()
    first = True
    for name, frame in sheets.items():
        ws = wb.active if first else wb.create_sheet()
        first = False
        ws.title = name[:31]
        for row in dataframe_to_rows(frame, index=False, header=True):
            ws.append(row)
    wb.save(path)


def _docx_escape(text: str) -> str:
    return escape(text).replace("\n", "</w:t><w:br/><w:t>")


def _write_docx(path: Path, paragraphs: list[str]) -> None:
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>
"""
    body = "".join(f"<w:p><w:r><w:t>{_docx_escape(p)}</w:t></w:r></w:p>" for p in paragraphs)
    document = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:body>{body}<w:sectPr/></w:body>
</w:document>
"""
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", rels)
        archive.writestr("word/document.xml", document)


def build(output_dir: Path, command: str) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    recall, evidence, degora_best = _read_inputs()
    recall_long = _recall_long(degora_best)
    evidence_top100 = _evidence_top100(evidence)
    heatmap = _heatmap_data(recall)
    source_xlsx = output_dir / "figure2_benchmark_source_data.xlsx"
    recall_csv = output_dir / "figure2_benchmark_degora_recall.csv"
    heatmap_csv = output_dir / "figure2_benchmark_method_recall100.csv"
    legend = output_dir / "figure2_benchmark_legend.docx"
    manifest = output_dir / "figure2_benchmark_manifest.json"
    validation = output_dir / "figure2_benchmark_validation.txt"
    recall_long.to_csv(recall_csv, index=False)
    heatmap.to_csv(heatmap_csv, index=False)
    _write_xlsx(
        source_xlsx,
        {"all_method_recall": recall, "degora_best": degora_best, "degora_recall_long": recall_long, "top100_evidence": evidence_top100, "method_heatmap": heatmap},
    )
    figure_outputs = []
    figure_outputs.extend(_plot_combined(recall_long, heatmap, evidence_top100, degora_best, output_dir / "figure2_benchmark"))
    figure_outputs.extend(_plot_panel_a(recall_long, output_dir / "figure2_benchmark_A_recall"))
    figure_outputs.extend(_plot_panel_b(heatmap, output_dir / "figure2_benchmark_B_methods"))
    figure_outputs.extend(_plot_panel_c(evidence_top100, output_dir / "figure2_benchmark_C_support"))
    figure_outputs.extend(_plot_panel_d(degora_best, output_dir / "figure2_benchmark_D_cutoffs"))
    _write_docx(
        legend,
        [
            "Figure 2. DEGORA benchmark against summary-level ranking methods.",
            "A, DEGORA recall across top-K cutoffs for IFN, ER stress, heat shock, and hypoxia/HIF1 locked panels. B, Recall@100 for runnable public-summary comparators. C, Median source-unit support among top-100 quality-weighted DEGORA genes. D, Recall gains from top10 to top50 and top100.",
        ],
    )
    outputs = [*figure_outputs, source_xlsx, recall_csv, heatmap_csv, legend, manifest, validation]
    validation.write_text("\n".join([f"{p}: exists={p.exists()} size_bytes={p.stat().st_size if p.exists() else 0}" for p in outputs if p != manifest]) + "\n")
    manifest.write_text(
        json.dumps(
            {
                "figure_id": FIGURE_ID,
                "generated_at": datetime.now(UTC).isoformat(),
                "script": "outputs/code/figures/make_benchmark_figure.py",
                "command": command,
                "inputs": [str(item[2]) for item in BENCHMARKS] + [str(item[3]) for item in BENCHMARKS],
                "outputs": [str(p) for p in outputs],
                "source_data": [str(source_xlsx), str(recall_csv), str(heatmap_csv)],
                "panel_claims": ["Recall and support values are read from generated comparator and deep-metric tables."],
                "transformations": {"degora_default": "fixed quality-weighted DEGORA row per benchmark when present, otherwise the unweighted DEGORA reference row; no recall-based row selection"},
                "validation": {"benchmarks": degora_best[["benchmark", "recall_at_10", "recall_at_50", "recall_at_100"]].to_dict(orient="records")},
                "known_limitations": ["Recall panels use locked positive gene panels and do not estimate false-positive rates."],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    inputs = [item[2] for item in BENCHMARKS] + [item[3] for item in BENCHMARKS]
    for artifact in outputs:
        write_source_sidecar(artifact, command, inputs=inputs, metadata={"generator": "benchmark-figure", "figure_id": FIGURE_ID})
    return {"figure_id": FIGURE_ID, "outputs": [str(p) for p in outputs]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    command = shell_command(["PYTHONPATH=outputs/code", "python", "outputs/code/figures/make_benchmark_figure.py", "--output-dir", args.output_dir])
    print(json.dumps(build(args.output_dir, command), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
