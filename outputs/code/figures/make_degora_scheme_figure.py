#!/usr/bin/env python
"""Create Figure 1: DEGORA workflow scheme from structured source data."""

from __future__ import annotations

import argparse
import json
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

from degora.provenance import shell_command, write_source_sidecar


FIGURE_ID = "FIGURE_1_SCHEME"
DEFAULT_OUTPUT = Path("outputs/figures/manuscript/FIGURE_1_SCHEME")


WORKFLOW_ROWS = [
    {"panel": "A", "order": 1, "module": "Public evidence", "label": "RNA-seq DEG tables", "detail": "gene, log2FC, p-value, source metadata"},
    {"panel": "A", "order": 2, "module": "Public evidence", "label": "Microarray DEG tables", "detail": "probe-collapsed gene summaries"},
    {"panel": "A", "order": 3, "module": "Harmonize", "label": "Common gene-level schema", "detail": "symbol, signed z, normalized rank, direction"},
    {"panel": "A", "order": 4, "module": "Source units", "label": "Collapse related contrasts", "detail": "paper/source-unit mean aggregation"},
    {"panel": "A", "order": 5, "module": "Score", "label": "DEGORA ranking", "detail": "support, direction, evidence, rank, effect"},
    {"panel": "A", "order": 6, "module": "Outputs", "label": "Evidence database", "detail": "CSV, SQLite, API, HTML atlas"},
]

SCORE_ROWS = [
    {"panel": "B", "component": "support_score", "weight": 0.30, "meaning": "independent source-unit support"},
    {"panel": "B", "component": "direction_score", "weight": 0.25, "meaning": "same-direction sign concordance"},
    {"panel": "B", "component": "evidence_score", "weight": 0.20, "meaning": "combined signed evidence strength"},
    {"panel": "B", "component": "rank_score_component", "weight": 0.15, "meaning": "within-source rank signal"},
    {"panel": "B", "component": "effect_score", "weight": 0.10, "meaning": "absolute weighted log2FC"},
]

OUTPUT_ROWS = [
    {"panel": "C", "artifact": "degora_gene_scores.csv", "purpose": "ranked gene table"},
    {"panel": "C", "artifact": "degora_scores.db", "purpose": "local SQLite evidence database"},
    {"panel": "C", "artifact": "degora_atlas.html", "purpose": "searchable browser output"},
    {"panel": "C", "artifact": "figure/source-data packages", "purpose": "manuscript and supplementary outputs"},
]


def _style() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none"})


def _draw_box(ax, xy, width, height, title, detail, color):
    rect = Rectangle(xy, width, height, facecolor=color, edgecolor="#1f2937", linewidth=0.8)
    ax.add_patch(rect)
    ax.text(xy[0] + width / 2, xy[1] + height * 0.62, title, ha="center", va="center", fontweight="bold", fontsize=8.5)
    ax.text(xy[0] + width / 2, xy[1] + height * 0.32, detail, ha="center", va="center", fontsize=7.2, color="#374151", wrap=True)


def _draw_arrow(ax, start, end):
    arrow = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, linewidth=1.0, color="#4b5563")
    ax.add_patch(arrow)


def _plot_panel_a(path: Path) -> list[Path]:
    _style()
    fig, ax = plt.subplots(figsize=(8.4, 2.9))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 3)
    ax.axis("off")
    rows = pd.DataFrame(WORKFLOW_ROWS)
    xs = [0.2, 0.2, 2.5, 4.8, 7.1, 9.4]
    ys = [1.65, 0.35, 1.0, 1.0, 1.0, 1.0]
    colors = ["#dbeafe", "#fce7f3", "#ecfdf5", "#fef3c7", "#ede9fe", "#e0f2fe"]
    for row, x, y, color in zip(rows.itertuples(index=False), xs, ys, colors, strict=True):
        _draw_box(ax, (x, y), 1.85, 0.9, row.label, row.detail, color)
    for start_y in [2.1, 0.8]:
        _draw_arrow(ax, (2.05, start_y), (2.45, 1.45))
    for x1, x2 in [(4.35, 4.75), (6.65, 7.05), (8.95, 9.35)]:
        _draw_arrow(ax, (x1, 1.45), (x2, 1.45))
    ax.text(0.02, 2.82, "A. Workflow", fontweight="bold", fontsize=10)
    return _save(fig, path)


def _plot_panel_b(path: Path) -> list[Path]:
    _style()
    rows = pd.DataFrame(SCORE_ROWS)
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.barh(rows["component"], rows["weight"], color=["#0f766e", "#265d97", "#9a6a12", "#7c3aed", "#b4234b"])
    ax.set_xlim(0, 0.34)
    ax.set_xlabel("Primary score weight")
    ax.set_title("B. Score components", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb")
    for i, row in enumerate(rows.itertuples(index=False)):
        ax.text(row.weight + 0.008, i, f"{row.weight:.2f}", va="center", fontsize=8)
    fig.tight_layout()
    return _save(fig, path)


def _plot_panel_c(path: Path) -> list[Path]:
    _style()
    fig, ax = plt.subplots(figsize=(4.9, 3.2))
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 4.8)
    ax.axis("off")
    rows = pd.DataFrame(OUTPUT_ROWS)
    for i, row in enumerate(rows.itertuples(index=False)):
        y = 3.7 - i * 1.05
        _draw_box(ax, (0.25, y), 4.5, 0.72, row.artifact, row.purpose, "#f8fafc")
    ax.text(0.04, 4.55, "C. User-facing outputs", fontweight="bold", fontsize=10)
    return _save(fig, path)


def _plot_combined(path: Path) -> list[Path]:
    _style()
    fig = plt.figure(figsize=(11.2, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], wspace=0.25, hspace=0.35)
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.set_xlim(0, 12)
    ax_a.set_ylim(0, 3)
    ax_a.axis("off")
    rows = pd.DataFrame(WORKFLOW_ROWS)
    xs = [0.2, 0.2, 2.5, 4.8, 7.1, 9.4]
    ys = [1.65, 0.35, 1.0, 1.0, 1.0, 1.0]
    colors = ["#dbeafe", "#fce7f3", "#ecfdf5", "#fef3c7", "#ede9fe", "#e0f2fe"]
    for row, x, y, color in zip(rows.itertuples(index=False), xs, ys, colors, strict=True):
        _draw_box(ax_a, (x, y), 1.85, 0.9, row.label, row.detail, color)
    for start_y in [2.1, 0.8]:
        _draw_arrow(ax_a, (2.05, start_y), (2.45, 1.45))
    for x1, x2 in [(4.35, 4.75), (6.65, 7.05), (8.95, 9.35)]:
        _draw_arrow(ax_a, (x1, 1.45), (x2, 1.45))
    ax_a.text(0.02, 2.82, "A. Workflow", fontweight="bold", fontsize=10)

    ax_b = fig.add_subplot(gs[1, 0])
    score = pd.DataFrame(SCORE_ROWS)
    ax_b.barh(score["component"], score["weight"], color=["#0f766e", "#265d97", "#9a6a12", "#7c3aed", "#b4234b"])
    ax_b.set_xlim(0, 0.34)
    ax_b.set_xlabel("Primary score weight")
    ax_b.set_title("B. Score components", loc="left", fontweight="bold")
    ax_b.grid(axis="x", color="#e5e7eb")

    ax_c = fig.add_subplot(gs[1, 1])
    ax_c.set_xlim(0, 5)
    ax_c.set_ylim(0, 4.8)
    ax_c.axis("off")
    for i, row in enumerate(pd.DataFrame(OUTPUT_ROWS).itertuples(index=False)):
        _draw_box(ax_c, (0.25, 3.7 - i * 1.05), 4.5, 0.72, row.artifact, row.purpose, "#f8fafc")
    ax_c.text(0.04, 4.55, "C. User-facing outputs", fontweight="bold", fontsize=10)
    fig.suptitle("DEGORA public DEG evidence integration", x=0.01, y=0.995, ha="left", fontsize=12, fontweight="bold")
    return _save(fig, path)


def _save(fig: plt.Figure, stem: Path) -> list[Path]:
    outputs = []
    for suffix in [".png", ".pdf", ".svg"]:
        out = stem.with_suffix(suffix)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        outputs.append(out)
    plt.close(fig)
    return outputs


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
    workflow = pd.DataFrame(WORKFLOW_ROWS)
    score = pd.DataFrame(SCORE_ROWS)
    outputs_table = pd.DataFrame(OUTPUT_ROWS)
    source_xlsx = output_dir / "figure1_scheme_source_data.xlsx"
    legend = output_dir / "figure1_scheme_legend.docx"
    manifest = output_dir / "figure1_scheme_manifest.json"
    validation = output_dir / "figure1_scheme_validation.txt"
    _write_xlsx(source_xlsx, {"workflow": workflow, "score_components": score, "outputs": outputs_table})
    figure_outputs = []
    figure_outputs.extend(_plot_combined(output_dir / "figure1_scheme"))
    figure_outputs.extend(_plot_panel_a(output_dir / "figure1_scheme_A_workflow"))
    figure_outputs.extend(_plot_panel_b(output_dir / "figure1_scheme_B_score"))
    figure_outputs.extend(_plot_panel_c(output_dir / "figure1_scheme_C_outputs"))
    _write_docx(
        legend,
        [
            "Figure 1. DEGORA workflow.",
            "A, Public RNA-seq and microarray DEG tables are harmonized into a common gene-level schema and collapsed by independent source units before scoring. B, Primary score weights for source support, direction concordance, evidence strength, rank signal, and effect size. C, Output files used for analysis, local browsing, and manuscript reporting.",
        ],
    )
    all_outputs = [*figure_outputs, source_xlsx, legend, manifest, validation]
    validation.write_text("\n".join([f"{p}: exists={p.exists()} size_bytes={p.stat().st_size if p.exists() else 0}" for p in all_outputs if p != manifest]) + "\n")
    manifest.write_text(
        json.dumps(
            {
                "figure_id": FIGURE_ID,
                "generated_at": datetime.now(UTC).isoformat(),
                "script": "outputs/code/figures/make_degora_scheme_figure.py",
                "command": command,
                "inputs": [],
                "outputs": [str(p) for p in all_outputs],
                "source_data": [str(source_xlsx)],
                "panel_claims": ["Workflow steps, score weights, and outputs are declared in structured source-data sheets."],
                "transformations": {"layout": "schematic generated from workflow, score, and output tables"},
                "validation": {"workflow_rows": len(workflow), "score_rows": len(score), "output_rows": len(outputs_table)},
                "known_limitations": ["Schematic panel shows analysis structure, not quantitative benchmark performance."],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    for artifact in all_outputs:
        write_source_sidecar(artifact, command, metadata={"generator": "degora-scheme-figure", "figure_id": FIGURE_ID})
    return {"figure_id": FIGURE_ID, "outputs": [str(p) for p in all_outputs]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    command = shell_command(["PYTHONPATH=outputs/code", "python", "outputs/code/figures/make_degora_scheme_figure.py", "--output-dir", args.output_dir])
    print(json.dumps(build(args.output_dir, command), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
