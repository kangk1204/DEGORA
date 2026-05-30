#!/usr/bin/env python
"""Create the cross-platform RNA-seq plus microarray benchmark figure package."""

from __future__ import annotations

import argparse
import json
import math
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

from degora.provenance import shell_command, write_source_sidecar


FIGURE_ID = "DEGORA_CROSS_PLATFORM"
DEFAULT_RESULTS = Path("outputs/results/cross-platform-microarray")
DEFAULT_OUTPUT = Path("outputs/figures/manuscript/DEGORA_CROSS_PLATFORM")
TOPIC_ORDER = ["IFN", "ER stress", "Hypoxia"]
MODE_ORDER = ["RNA-seq only", "RNA-seq + microarray"]
COLORS = {"RNA-seq only": "#265d97", "RNA-seq + microarray": "#0f766e"}
TOPIC_MARKERS = {"IFN": "o", "ER stress": "s", "Hypoxia": "^"}
TOPIC_LINESTYLES = {"IFN": "-", "ER stress": "--", "Hypoxia": ":"}


def _set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 8.5,
            "legend.fontsize": 8,
            "svg.fonttype": "none",
        }
    )


def _save(fig: plt.Figure, stem: Path) -> list[Path]:
    outputs = []
    for suffix in (".png", ".pdf", ".svg"):
        out = stem.with_suffix(suffix)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        outputs.append(out)
    plt.close(fig)
    return outputs


def _recall_long(benchmark: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in benchmark.itertuples(index=False):
        for k in [10, 20, 50, 100]:
            rows.append(
                {
                    "topic": row.topic,
                    "mode": row.mode,
                    "k": k,
                    "recall": float(getattr(row, f"recall_at_{k}")),
                    "direction_recall_at_100": float(row.direction_recall_at_100),
                    "n_source_units_total": int(row.n_source_units_total),
                    "n_contrasts_total": int(row.n_contrasts_total),
                    "top10": row.top10,
                }
            )
    frame = pd.DataFrame(rows)
    frame["topic"] = pd.Categorical(frame["topic"], TOPIC_ORDER, ordered=True)
    frame["mode"] = pd.Categorical(frame["mode"], MODE_ORDER, ordered=True)
    return frame.sort_values(["topic", "mode", "k"]).reset_index(drop=True)


def _support_summary(shifts: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for topic, group in shifts.groupby("topic", sort=False):
        rows.extend(
            [
                {
                    "topic": topic,
                    "mode": "RNA-seq only",
                    "median_source_units": float(group["rna_only_source_units"].median()),
                    "median_reliability": float(group["rna_only_reliability"].median()),
                    "median_top_percent": float(group["rna_only_top_percent"].median()),
                },
                {
                    "topic": topic,
                    "mode": "RNA-seq + microarray",
                    "median_source_units": float(group["mixed_source_units"].median()),
                    "median_reliability": float(group["mixed_reliability"].median()),
                    "median_top_percent": float(group["mixed_top_percent"].median()),
                },
            ]
        )
    frame = pd.DataFrame(rows)
    frame["topic"] = pd.Categorical(frame["topic"], TOPIC_ORDER, ordered=True)
    frame["mode"] = pd.Categorical(frame["mode"], MODE_ORDER, ordered=True)
    return frame.sort_values(["topic", "mode"]).reset_index(drop=True)


def _source_summary(datasets: pd.DataFrame, benchmark: pd.DataFrame) -> pd.DataFrame:
    micro = datasets.groupby("topic", as_index=False).agg(
        microarray_contrasts=("contrast", "count"),
        microarray_samples=("samples", "sum"),
        microarray_accessions=("accession", lambda values: ";".join(sorted(set(map(str, values))))),
    )
    mixed = benchmark.loc[benchmark["mode"].eq("RNA-seq + microarray"), ["topic", "n_source_units_total", "n_contrasts_total"]]
    out = mixed.merge(micro, on="topic", how="left")
    return out.sort_values("topic", key=lambda s: s.map({topic: i for i, topic in enumerate(TOPIC_ORDER)})).reset_index(drop=True)


def _plot_panel_a(recall: pd.DataFrame, stem: Path) -> list[Path]:
    _set_style()
    fig, axes = plt.subplots(1, len(TOPIC_ORDER), figsize=(3.25 * len(TOPIC_ORDER), 3.2), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, topic in zip(axes, TOPIC_ORDER, strict=True):
        subset = recall.loc[recall["topic"].astype(str).eq(topic)]
        for mode in MODE_ORDER:
            group = subset.loc[subset["mode"].astype(str).eq(mode)]
            ax.plot(group["k"], group["recall"], marker="o", linewidth=2.2, color=COLORS[mode], label=mode)
        ax.set_title(topic, loc="left", fontweight="bold")
        ax.set_xlabel("Top K genes")
        ax.set_xticks([10, 20, 50, 100])
        ax.set_ylim(0, 1.02)
        ax.grid(True, color="#e5e7eb")
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Gold-panel recall")
    axes[-1].legend(frameon=False, loc="lower right")
    fig.suptitle("A. RNA-seq-only versus RNA-seq+microarray recall", x=0.01, y=1.02, ha="left", fontsize=11, fontweight="bold")
    return _save(fig, stem)


def _plot_panel_b(shifts: pd.DataFrame, stem: Path) -> list[Path]:
    _set_style()
    fig, axes = plt.subplots(1, len(TOPIC_ORDER), figsize=(3.25 * len(TOPIC_ORDER), 3.3), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, topic in zip(axes, TOPIC_ORDER, strict=True):
        subset = shifts.loc[shifts["topic"].eq(topic)].copy()
        ax.scatter(
            subset["rna_only_top_percent"],
            subset["mixed_top_percent"],
            s=34,
            color="#0f766e",
            alpha=0.82,
            edgecolor="#1f2937",
            linewidth=0.35,
        )
        limit = max(float(subset["rna_only_top_percent"].max()), float(subset["mixed_top_percent"].max()), 0.1)
        ax.plot([0, limit], [0, limit], color="#9ca3af", linestyle="--", linewidth=1.0)
        labels = subset.sort_values("top_percent_delta", ascending=False).head(5)
        for row in labels.itertuples(index=False):
            ax.text(row.rna_only_top_percent, row.mixed_top_percent, row.gene_symbol, fontsize=7)
        ax.set_title(topic, loc="left", fontweight="bold")
        ax.set_xlabel("RNA-seq-only top percent")
        ax.grid(True, color="#e5e7eb")
    axes[0].set_ylabel("RNA-seq+microarray top percent")
    fig.suptitle("B. Locked marker rank shifts", x=0.01, y=1.02, ha="left", fontsize=11, fontweight="bold")
    return _save(fig, stem)


def _plot_panel_c(support: pd.DataFrame, stem: Path) -> list[Path]:
    _set_style()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2), sharey=False)
    x = np.arange(len(TOPIC_ORDER))
    width = 0.34
    for i, mode in enumerate(MODE_ORDER):
        subset = support.loc[support["mode"].astype(str).eq(mode)].set_index("topic").reindex(TOPIC_ORDER)
        axes[0].bar(x + (i - 0.5) * width, subset["median_source_units"], width=width, color=COLORS[mode], label=mode)
        axes[1].bar(x + (i - 0.5) * width, subset["median_reliability"], width=width, color=COLORS[mode], label=mode)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(TOPIC_ORDER)
        ax.grid(axis="y", color="#e5e7eb")
        ax.set_axisbelow(True)
    axes[0].set_title("Median source units", loc="left", fontweight="bold")
    axes[1].set_title("Median reliability score", loc="left", fontweight="bold")
    axes[0].set_ylabel("Gold-panel markers")
    axes[1].set_ylim(0, 105)
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("C. Marker-level evidence support", x=0.01, y=1.02, ha="left", fontsize=11, fontweight="bold")
    return _save(fig, stem)


def _plot_panel_d(sources: pd.DataFrame, stem: Path) -> list[Path]:
    _set_style()
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    x = np.arange(len(sources))
    ax.bar(x - 0.17, sources["microarray_samples"], width=0.34, color="#9a6a12", label="Microarray samples")
    ax2 = ax.twinx()
    ax2.bar(x + 0.17, sources["microarray_contrasts"], width=0.34, color="#6b7280", label="Microarray contrasts")
    ax.set_xticks(x)
    ax.set_xticklabels(sources["topic"])
    ax.set_ylabel("Samples")
    ax2.set_ylabel("Contrasts")
    ax.set_title("D. Added microarray evidence", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#e5e7eb")
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles + handles2, labels + labels2, frameon=False, loc="upper left")
    return _save(fig, stem)


def _plot_combined(
    recall: pd.DataFrame,
    shifts: pd.DataFrame,
    support: pd.DataFrame,
    sources: pd.DataFrame,
    stem: Path,
) -> list[Path]:
    _set_style()
    fig = plt.figure(figsize=(11.4, 8.4))
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.43)

    ax_a = fig.add_subplot(gs[0, 0])
    for topic in TOPIC_ORDER:
        for mode in MODE_ORDER:
            group = recall.loc[recall["topic"].astype(str).eq(topic) & recall["mode"].astype(str).eq(mode)]
            ax_a.plot(
                group["k"],
                group["recall"],
                marker=TOPIC_MARKERS.get(topic, "o"),
                linestyle=TOPIC_LINESTYLES.get(topic, "-"),
                color=COLORS[mode],
                label=f"{topic} {mode}",
            )
    ax_a.set_title("A. Gold-panel recall", loc="left", fontweight="bold")
    ax_a.set_xlabel("Top K genes")
    ax_a.set_ylabel("Recall")
    ax_a.set_ylim(0, 1.02)
    ax_a.set_xticks([10, 20, 50, 100])
    ax_a.grid(True, color="#e5e7eb")
    ax_a.legend(frameon=False, fontsize=6.8, loc="lower right")

    ax_b = fig.add_subplot(gs[0, 1])
    for topic in TOPIC_ORDER:
        subset = shifts.loc[shifts["topic"].eq(topic)]
        ax_b.scatter(
            subset["rna_only_top_percent"],
            subset["mixed_top_percent"],
            s=32,
            marker=TOPIC_MARKERS.get(topic, "o"),
            alpha=0.78,
            label=topic,
        )
    limit = max(float(shifts["rna_only_top_percent"].max()), float(shifts["mixed_top_percent"].max()), 0.1)
    ax_b.plot([0, limit], [0, limit], color="#9ca3af", linestyle="--", linewidth=1.0)
    for row in shifts.sort_values("top_percent_delta", ascending=False).head(8).itertuples(index=False):
        ax_b.text(row.rna_only_top_percent, row.mixed_top_percent, row.gene_symbol, fontsize=6.7)
    ax_b.set_title("B. Marker top-percent shift", loc="left", fontweight="bold")
    ax_b.set_xlabel("RNA-seq only")
    ax_b.set_ylabel("RNA-seq + microarray")
    ax_b.grid(True, color="#e5e7eb")
    ax_b.legend(frameon=False, loc="upper left")

    ax_c = fig.add_subplot(gs[1, 0])
    x = np.arange(len(TOPIC_ORDER))
    width = 0.34
    for i, mode in enumerate(MODE_ORDER):
        subset = support.loc[support["mode"].astype(str).eq(mode)].set_index("topic").reindex(TOPIC_ORDER)
        ax_c.bar(x + (i - 0.5) * width, subset["median_source_units"], width=width, color=COLORS[mode], label=mode)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(TOPIC_ORDER)
    ax_c.set_ylabel("Median source units")
    ax_c.set_title("C. Marker support", loc="left", fontweight="bold")
    ax_c.grid(axis="y", color="#e5e7eb")
    ax_c.legend(frameon=False)

    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.bar(np.arange(len(sources)), sources["microarray_samples"], color="#9a6a12")
    ax_d.set_xticks(np.arange(len(sources)))
    ax_d.set_xticklabels(sources["topic"])
    ax_d.set_ylabel("Added microarray samples")
    ax_d.set_title("D. Added microarray data", loc="left", fontweight="bold")
    ax_d.grid(axis="y", color="#e5e7eb")

    fig.suptitle("DEGORA cross-platform benchmark", x=0.01, y=0.995, ha="left", fontsize=12, fontweight="bold")
    return _save(fig, stem)


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
    body = "".join(f"<w:p><w:r><w:t>{_docx_escape(paragraph)}</w:t></w:r></w:p>" for paragraph in paragraphs)
    document = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:body>{body}<w:sectPr/></w:body>
</w:document>
"""
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", rels)
        archive.writestr("word/document.xml", document)


def _write_xlsx(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    workbook = Workbook()
    first = True
    for name, frame in sheets.items():
        ws = workbook.active if first else workbook.create_sheet()
        first = False
        ws.title = name[:31]
        for row in dataframe_to_rows(frame, index=False, header=True):
            ws.append(row)
    workbook.save(path)


def _write_validation(path: Path, outputs: list[Path], sheets: dict[str, pd.DataFrame]) -> None:
    lines = ["Cross-platform microarray figure validation", ""]
    for out in outputs:
        lines.append(f"{out}: exists={out.exists()} size_bytes={out.stat().st_size if out.exists() else 0}")
    lines.append("")
    for name, frame in sheets.items():
        lines.append(f"{name}: rows={len(frame)} columns={len(frame.columns)}")
    path.write_text("\n".join(lines) + "\n")


def _write_manifest(
    path: Path,
    *,
    command: str,
    inputs: list[Path],
    outputs: list[Path],
    source_data: list[Path],
    benchmark: pd.DataFrame,
    shifts: pd.DataFrame,
) -> None:
    manifest = {
        "figure_id": FIGURE_ID,
        "generated_at": datetime.now(UTC).isoformat(),
        "script": "outputs/code/figures/make_cross_platform_microarray_figure.py",
        "command": command,
        "inputs": [str(path) for path in inputs],
        "outputs": [str(path) for path in outputs],
        "source_data": [str(path) for path in source_data],
        "panel_claims": [
            "A, recall values are read from cross_platform_benchmark_summary.csv.",
            "B, marker top-percent shifts are computed from quality-weighted DEGORA ranks in RNA-seq-only and mixed corpora.",
            "C, source support is summarized across locked gold-panel marker rows.",
            "D, added microarray sample counts are read from cross_platform_dataset_sources.csv.",
        ],
        "transformations": {
            "recall": "long-form recall@10/20/50/100 from the selected DEGORA ranking per benchmark mode",
            "rank_shift": "positive rank_delta or top_percent_delta means the marker moved upward after microarray addition",
            "support": "median source units and evidence reliability among locked marker genes",
        },
        "validation": {
            "benchmark_rows": benchmark.to_dict(orient="records"),
            "median_rank_delta": shifts.groupby("topic")["rank_delta"].median().to_dict(),
            "median_top_percent_delta": shifts.groupby("topic")["top_percent_delta"].median().to_dict(),
        },
        "known_limitations": [
            "Microarray DEG values are derived from GEO processed series matrices with a Welch fallback; limma full tables remain preferred when available.",
            "GSE19519 tunicamycin and thapsigargin rows share one source unit to avoid source inflation.",
            "Recall is against compact locked panels and does not estimate false-positive rates.",
        ],
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def build_figure_package(results_dir: Path, output_dir: Path, command: str) -> dict[str, object]:
    benchmark = pd.read_csv(results_dir / "cross_platform_benchmark_summary.csv")
    shifts = pd.read_csv(results_dir / "cross_platform_marker_rank_shifts.csv")
    datasets = pd.read_csv(results_dir / "cross_platform_dataset_sources.csv")
    recall = _recall_long(benchmark)
    support = _support_summary(shifts)
    sources = _source_summary(datasets, benchmark)

    output_dir.mkdir(parents=True, exist_ok=True)
    source_xlsx = output_dir / "degora_cross_platform_source_data.xlsx"
    legend = output_dir / "degora_cross_platform_legend.docx"
    manifest = output_dir / "degora_cross_platform_manifest.json"
    validation = output_dir / "degora_cross_platform_validation.txt"
    recall_csv = output_dir / "degora_cross_platform_recall_long.csv"
    shifts_csv = output_dir / "degora_cross_platform_marker_rank_shifts.csv"

    recall.to_csv(recall_csv, index=False)
    shifts.to_csv(shifts_csv, index=False)
    sheets = {
        "benchmark_summary": benchmark,
        "recall_long": recall,
        "marker_rank_shifts": shifts,
        "support_summary": support,
        "dataset_sources": datasets,
        "source_summary": sources,
    }
    _write_xlsx(source_xlsx, sheets)

    figure_outputs: list[Path] = []
    figure_outputs.extend(_plot_combined(recall, shifts, support, sources, output_dir / "degora_cross_platform"))
    figure_outputs.extend(_plot_panel_a(recall, output_dir / "degora_cross_platform_A_recall"))
    figure_outputs.extend(_plot_panel_b(shifts, output_dir / "degora_cross_platform_B_rank_shift"))
    figure_outputs.extend(_plot_panel_c(support, output_dir / "degora_cross_platform_C_support"))
    figure_outputs.extend(_plot_panel_d(sources, output_dir / "degora_cross_platform_D_sources"))

    _write_docx(
        legend,
        [
            "DEGORA_CROSS_PLATFORM. RNA-seq and microarray integration benchmark.",
            "A, Gold-panel recall for RNA-seq-only and RNA-seq+microarray DEGORA rankings in IFN, ER-stress, and hypoxia topics. B, Top-percent ranks for locked marker genes before and after microarray addition; points below the diagonal moved upward in the mixed analysis. C, Median source-unit support and evidence reliability for locked markers. D, Microarray samples added from GSE71634, GSE19519, GSE3045, and GSE22282.",
        ],
    )

    outputs = [*figure_outputs, source_xlsx, recall_csv, shifts_csv, legend, manifest, validation]
    inputs = [
        results_dir / "cross_platform_benchmark_summary.csv",
        results_dir / "cross_platform_marker_rank_shifts.csv",
        results_dir / "cross_platform_dataset_sources.csv",
    ]
    _write_validation(validation, [*figure_outputs, source_xlsx, recall_csv, shifts_csv, legend], sheets)
    _write_manifest(
        manifest,
        command=command,
        inputs=inputs,
        outputs=outputs,
        source_data=[source_xlsx, recall_csv, shifts_csv],
        benchmark=benchmark,
        shifts=shifts,
    )
    for artifact in outputs:
        write_source_sidecar(
            artifact,
            command,
            inputs=inputs,
            metadata={"generator": "cross-platform-microarray-figure", "figure_id": FIGURE_ID},
        )

    return {
        "figure_id": FIGURE_ID,
        "output_dir": str(output_dir),
        "n_benchmark_rows": int(len(benchmark)),
        "n_marker_rows": int(len(shifts)),
        "outputs": [str(path) for path in outputs],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = shell_command(
        [
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/figures/make_cross_platform_microarray_figure.py",
            "--results-dir",
            args.results_dir,
            "--output-dir",
            args.output_dir,
        ]
    )
    summary = build_figure_package(args.results_dir, args.output_dir, command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
