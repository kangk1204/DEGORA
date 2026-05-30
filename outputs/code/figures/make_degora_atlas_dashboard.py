#!/usr/bin/env python
"""Build a static searchable DEGORA atlas and manuscript figure package."""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

from degora.provenance import shell_command, write_source_sidecar


FIGURE_ID = "DEGORA_ATLAS"
OUTPUT_STEM = "degora_atlas"


@dataclass(frozen=True)
class CorpusInput:
    corpus_id: str
    label: str
    topic: str
    result_dir: Path
    comparator_summary: Path | None = None
    metadata: Path | None = None

    @property
    def score_csv(self) -> Path:
        return self.result_dir / "degora_gene_scores.csv"

    @property
    def score_db(self) -> Path:
        return self.result_dir / "degora_scores.db"

    @property
    def metadata_json(self) -> Path:
        return self.metadata or self.result_dir / "degora_score_metadata.json"


DEFAULT_CORPORA = [
    CorpusInput(
        "ifn-pilot",
        "IFN",
        "interferon response",
        Path("outputs/results/ifn-pilot"),
        Path("outputs/results/ifn-pilot/ifn_gold_comparator_summary.csv"),
    ),
    CorpusInput(
        "er-stress-cross-platform",
        "ER stress",
        "ER stress / UPR RNA-seq + microarray",
        Path("outputs/results/er-stress-cross-platform"),
        Path("outputs/results/er-stress-cross-platform/er_stress_cross_platform_gold_comparator_summary.csv"),
    ),
    CorpusInput(
        "heat-shock-benchmark",
        "Heat shock",
        "heat shock / HSF1",
        Path("outputs/results/heat-shock-benchmark"),
        Path("outputs/results/heat-shock-benchmark/heat_shock_gold_comparator_summary.csv"),
    ),
    CorpusInput(
        "hypoxia-hif1-benchmark",
        "Hypoxia",
        "hypoxia / HIF1",
        Path("outputs/results/hypoxia-hif1-benchmark"),
        Path("outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv"),
    ),
    CorpusInput(
        "ifn-cross-platform",
        "IFN mixed",
        "interferon response RNA-seq + microarray",
        Path("outputs/results/ifn-cross-platform"),
        Path("outputs/results/ifn-cross-platform/ifn_cross_platform_gold_comparator_summary.csv"),
    ),
    CorpusInput(
        "hypoxia-cross-platform",
        "Hypoxia mixed",
        "hypoxia RNA-seq + microarray",
        Path("outputs/results/hypoxia-cross-platform"),
        Path("outputs/results/hypoxia-cross-platform/hypoxia_cross_platform_gold_comparator_summary.csv"),
    ),
]


GENE_INDEX_COLUMNS = [
    "corpus_id",
    "corpus_label",
    "topic",
    "gene_symbol",
    "display_rank",
    "display_score",
    "display_top_percent",
    "display_direction",
    "display_sign_concordance",
    "degora_rank",
    "degora_score",
    "priority_score",
    "priority_top_percent",
    "quality_weighted_degora_rank",
    "quality_weighted_degora_score",
    "quality_weighted_top_percent",
    "n_source_units",
    "n_contrasts_observed",
    "weighted_lfc",
    "evidence_tier",
    "evidence_reliability_score",
    "direction_confidence_index",
    "quality_weighted_direction_confidence_index",
    "loo_rank_stability_score",
    "loo_top50_fraction",
    "loo_top100_fraction",
    "high_confidence",
]


METHOD_LABELS = {
    "degora_deg_score": "DEGORA",
    "degora_quality_weighted_score": "DEGORA quality",
    "weighted_stouffer": "Weighted Stouffer",
    "fisher": "Fisher",
    "awfisher": "AWFisher",
    "metarnaseq_fisher": "metaRNASeq Fisher",
    "metarnaseq_invnorm": "metaRNASeq invnorm",
    "metavolcanor": "MetaVolcanoR",
    "robustrankaggreg": "RobustRankAggreg",
    "rank_product_approx": "Rank product approx.",
    "degora_slice": "DEGORA slice",
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


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _num(value: object, digits: int = 4) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return round(number, digits)


def _read_gene_scores(corpus: CorpusInput) -> pd.DataFrame:
    frame = pd.read_csv(corpus.score_csv)
    frame = frame.copy()
    frame["corpus_id"] = corpus.corpus_id
    frame["corpus_label"] = corpus.label
    frame["topic"] = corpus.topic

    if "quality_weighted_degora_rank" in frame.columns:
        frame["display_rank"] = pd.to_numeric(frame["quality_weighted_degora_rank"], errors="coerce")
        frame["display_score"] = pd.to_numeric(frame.get("quality_weighted_degora_score"), errors="coerce")
        frame["display_top_percent"] = pd.to_numeric(frame.get("quality_weighted_top_percent"), errors="coerce")
        frame["display_direction"] = frame.get("quality_weighted_consensus_direction", frame.get("consensus_direction", ""))
        frame["display_sign_concordance"] = pd.to_numeric(
            frame.get("quality_weighted_sign_concordance", frame.get("sign_concordance")), errors="coerce"
        )
    else:
        frame["display_rank"] = pd.to_numeric(frame["degora_rank"], errors="coerce")
        frame["display_score"] = pd.to_numeric(frame["degora_score"], errors="coerce")
        frame["display_top_percent"] = pd.to_numeric(frame["top_percent"], errors="coerce")
        frame["display_direction"] = frame.get("consensus_direction", "")
        frame["display_sign_concordance"] = pd.to_numeric(frame.get("sign_concordance"), errors="coerce")

    for column in GENE_INDEX_COLUMNS:
        if column not in frame.columns:
            frame[column] = np.nan
    frame = frame[GENE_INDEX_COLUMNS].copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.upper()
    frame = frame.dropna(subset=["gene_symbol", "display_rank"])
    frame = frame.sort_values(["corpus_id", "display_rank", "gene_symbol"])
    return frame


def _read_method_recall(corpus: CorpusInput) -> pd.DataFrame:
    if corpus.comparator_summary is None or not corpus.comparator_summary.exists():
        return pd.DataFrame()
    frame = pd.read_csv(corpus.comparator_summary)
    frame = frame.copy()
    frame["corpus_id"] = corpus.corpus_id
    frame["corpus_label"] = corpus.label
    frame["topic"] = corpus.topic
    frame["method_label"] = frame["method_id"].map(METHOD_LABELS).fillna(frame["method_id"].astype(str))
    for column in ["recall_at_10", "recall_at_20", "recall_at_50", "recall_at_100", "direction_recall_at_100"]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def _corpus_summary(corpora: list[CorpusInput], gene_index: pd.DataFrame, method_recall: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for corpus in corpora:
        subset = gene_index.loc[gene_index["corpus_id"].eq(corpus.corpus_id)]
        metadata = _read_json(corpus.metadata_json)
        recall_subset = method_recall.loc[
            method_recall["corpus_id"].eq(corpus.corpus_id)
            & method_recall["method_id"].eq("degora_deg_score")
            & method_recall["setting_id"].astype(str).str.contains("quality", case=False, na=False)
        ]
        if recall_subset.empty:
            recall_subset = method_recall.loc[
                method_recall["corpus_id"].eq(corpus.corpus_id)
                & method_recall["method_id"].isin(["degora_deg_score", "degora_quality_weighted_score"])
            ]
        recall_at_100 = _num(recall_subset["recall_at_100"].iloc[0]) if not recall_subset.empty else None
        recall_at_50 = _num(recall_subset["recall_at_50"].iloc[0]) if not recall_subset.empty else None
        top = subset.sort_values(["display_rank", "gene_symbol"]).head(1)
        rows.append(
            {
                "corpus_id": corpus.corpus_id,
                "corpus_label": corpus.label,
                "topic": corpus.topic,
                "n_genes": int(len(subset)),
                "n_source_units_total": int(metadata.get("n_source_units_total", 0) or 0),
                "n_contrasts_total": int(metadata.get("n_contrasts_total", 0) or 0),
                "n_evidence_rows": int(metadata.get("n_source_unit_gene_evidence_rows", 0) or 0),
                "top_gene": "" if top.empty else str(top["gene_symbol"].iloc[0]),
                "top_gene_rank": None if top.empty else int(float(top["display_rank"].iloc[0])),
                "degora_recall_at_50": recall_at_50,
                "degora_recall_at_100": recall_at_100,
                "score_version": str(metadata.get("score_version", "")),
                "score_warning": str(metadata.get("score_warning", "")),
            }
        )
    return pd.DataFrame(rows)


def _top100_evidence(gene_index: pd.DataFrame) -> pd.DataFrame:
    subset = gene_index.loc[pd.to_numeric(gene_index["display_rank"], errors="coerce").le(100)].copy()
    rows: list[dict[str, Any]] = []
    for (corpus_id, corpus_label), group in subset.groupby(["corpus_id", "corpus_label"], sort=False):
        rows.append(
            {
                "corpus_id": corpus_id,
                "corpus_label": corpus_label,
                "n_top100": int(len(group)),
                "median_source_units": _num(pd.to_numeric(group["n_source_units"], errors="coerce").median(), 3),
                "median_sign_concordance": _num(pd.to_numeric(group["display_sign_concordance"], errors="coerce").median(), 3),
                "median_evidence_reliability": _num(pd.to_numeric(group["evidence_reliability_score"], errors="coerce").median(), 3),
                "median_loo_stability": _num(pd.to_numeric(group["loo_rank_stability_score"], errors="coerce").median(), 3),
            }
        )
    return pd.DataFrame(rows)


def _prior_art_counts(prior_art: pd.DataFrame) -> pd.DataFrame:
    if prior_art.empty:
        return pd.DataFrame()
    frame = prior_art.copy()
    status_col = "public_summary_deg_status" if "public_summary_deg_status" in frame.columns else "current_pipeline_status"
    manuscript_col = "manuscript_use" if "manuscript_use" in frame.columns else status_col
    rows = (
        frame.groupby([status_col, manuscript_col], dropna=False)
        .size()
        .reset_index(name="n_methods_or_resources")
        .rename(columns={status_col: "public_summary_status", manuscript_col: "manuscript_use"})
    )
    return rows.sort_values(["public_summary_status", "manuscript_use"]).reset_index(drop=True)


def _html_payload(
    *,
    corpus_summary: pd.DataFrame,
    method_recall: pd.DataFrame,
    prior_art: pd.DataFrame,
    gene_index: pd.DataFrame,
) -> dict[str, Any]:
    search_columns = [
        "corpus_id",
        "corpus_label",
        "topic",
        "gene_symbol",
        "display_rank",
        "display_score",
        "display_top_percent",
        "display_direction",
        "display_sign_concordance",
        "n_source_units",
        "weighted_lfc",
        "evidence_tier",
        "evidence_reliability_score",
        "direction_confidence_index",
        "loo_rank_stability_score",
        "high_confidence",
    ]
    compact_genes = gene_index[search_columns].copy()
    for col in compact_genes.columns:
        if col not in {"corpus_id", "corpus_label", "topic", "gene_symbol", "display_direction", "evidence_tier"}:
            compact_genes[col] = compact_genes[col].map(lambda value: _num(value, 4))
    recall_cols = [
        "corpus_id",
        "corpus_label",
        "method_id",
        "method_label",
        "run_status",
        "recall_at_10",
        "recall_at_50",
        "recall_at_100",
        "direction_recall_at_100",
        "top10",
    ]
    recall_for_html = method_recall[[col for col in recall_cols if col in method_recall.columns]].copy()
    return {
        "generated_at": datetime.now(UTC).isoformat(),
        "corpora": corpus_summary.to_dict(orient="records"),
        "methods": recall_for_html.fillna("").to_dict(orient="records"),
        "prior_art": prior_art.fillna("").to_dict(orient="records"),
        "gene_columns": search_columns,
        "genes": compact_genes.fillna("").values.tolist(),
    }


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>DEGORA Atlas</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f7f8f6;
      --ink: #20252b;
      --muted: #59656f;
      --panel: #ffffff;
      --line: #dfe4dc;
      --green: #0f766e;
      --red: #b4234b;
      --blue: #265d97;
      --gold: #9a6a12;
      --gray: #6b7280;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 14px;
    }
    header {
      padding: 18px 24px 14px;
      background: var(--panel);
      border-bottom: 1px solid var(--line);
      position: sticky;
      top: 0;
      z-index: 3;
    }
    .topline {
      display: flex;
      align-items: flex-end;
      justify-content: space-between;
      gap: 16px;
      flex-wrap: wrap;
    }
    h1 { margin: 0; font-size: 22px; letter-spacing: 0; }
    .subtitle { color: var(--muted); margin-top: 4px; }
    .meta { color: var(--muted); font-variant-numeric: tabular-nums; }
    .summary {
      display: grid;
      grid-template-columns: repeat(4, minmax(140px, 1fr));
      gap: 10px;
      padding: 14px 24px;
      border-bottom: 1px solid var(--line);
      background: #fbfcfa;
    }
    .metric {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 10px 12px;
      min-height: 70px;
    }
    .metric span { display: block; color: var(--muted); font-size: 12px; }
    .metric strong { display: block; margin-top: 4px; font-size: 20px; font-variant-numeric: tabular-nums; }
    main {
      display: grid;
      grid-template-columns: minmax(0, 1.05fr) minmax(420px, .95fr);
      gap: 16px;
      padding: 16px 24px 24px;
    }
    section {
      min-width: 0;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
    }
    .section-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
    }
    h2 { margin: 0; font-size: 15px; letter-spacing: 0; }
    .controls {
      display: grid;
      grid-template-columns: minmax(170px, 1fr) 150px 120px 118px;
      gap: 8px;
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
      background: #fbfcfa;
    }
    input, select, button {
      width: 100%;
      min-width: 0;
      height: 34px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      color: var(--ink);
      font: inherit;
      padding: 0 10px;
    }
    button { background: var(--green); border-color: var(--green); color: #fff; font-weight: 700; cursor: pointer; }
    table { width: 100%; border-collapse: collapse; table-layout: fixed; }
    th, td {
      padding: 8px 9px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
    th { color: var(--muted); font-size: 12px; background: #fbfcfa; }
    tbody tr { cursor: pointer; }
    tbody tr:hover { background: #edf7f4; }
    .num { text-align: right; font-variant-numeric: tabular-nums; }
    .gene { font-weight: 800; }
    .badge {
      display: inline-flex;
      min-width: 46px;
      height: 22px;
      padding: 0 7px;
      align-items: center;
      justify-content: center;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 800;
      background: #e7f5f2;
      color: var(--green);
    }
    .badge.down { background: #fde7ed; color: var(--red); }
    .badge.flat, .badge.zero { background: #f3efe2; color: var(--gold); }
    .viz-grid { display: grid; grid-template-columns: 1fr; gap: 14px; padding: 14px; }
    .chart { min-height: 210px; border-bottom: 1px solid var(--line); padding-bottom: 14px; }
    .chart:last-child { border-bottom: 0; padding-bottom: 0; }
    .chart-title { display: flex; justify-content: space-between; gap: 10px; margin-bottom: 8px; color: var(--muted); }
    svg { width: 100%; height: auto; display: block; }
    .detail { padding: 14px; }
    .detail-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 10px; margin-bottom: 12px; }
    .note { color: var(--muted); line-height: 1.45; margin: 10px 0 0; }
    @media (max-width: 1080px) {
      main { grid-template-columns: 1fr; padding: 12px; }
      .summary { grid-template-columns: repeat(2, minmax(140px, 1fr)); padding: 12px; }
    }
    @media (max-width: 640px) {
      .controls { grid-template-columns: 1fr; }
      .detail-grid { grid-template-columns: 1fr; }
      .summary { grid-template-columns: 1fr; }
      th:nth-child(4), td:nth-child(4), th:nth-child(7), td:nth-child(7) { display: none; }
    }
  </style>
</head>
<body>
  <header>
    <div class="topline">
      <div>
        <h1>DEGORA Atlas</h1>
        <div class="subtitle">Searchable direction-aware evidence database built from public DEG summary files.</div>
      </div>
      <div class="meta" id="generated"></div>
    </div>
  </header>
  <div class="summary" id="summary"></div>
  <main>
    <section>
      <div class="section-head">
        <h2>Gene Search</h2>
        <div class="meta" id="status"></div>
      </div>
      <div class="controls">
        <input id="query" placeholder="Gene symbol, e.g. CMPK2, VEGFA, HSPA1A" autocomplete="off">
        <select id="corpus"></select>
        <select id="direction">
          <option value="">All directions</option>
          <option value="up">Up</option>
          <option value="down">Down</option>
          <option value="flat">Flat/zero</option>
        </select>
        <button id="search">Search</button>
      </div>
      <table>
        <thead>
          <tr>
            <th style="width:118px">Corpus</th>
            <th style="width:82px" class="num">Rank</th>
            <th>Gene</th>
            <th style="width:74px">Dir</th>
            <th style="width:82px" class="num">Top %</th>
            <th style="width:80px" class="num">Units</th>
            <th style="width:90px" class="num">Rel.</th>
          </tr>
        </thead>
        <tbody id="genes"></tbody>
      </table>
    </section>
    <section>
      <div class="section-head">
        <h2>Atlas Views</h2>
        <div class="meta">quality-weighted ranks</div>
      </div>
      <div class="viz-grid">
        <div class="chart">
          <div class="chart-title"><strong>Recall@100 by method</strong><span>locked/anchor panels</span></div>
          <div id="heatmap"></div>
        </div>
        <div class="chart">
          <div class="chart-title"><strong>Top-100 evidence support</strong><span>median source units and reliability</span></div>
          <div id="supportPlot"></div>
        </div>
        <div id="detail" class="detail">
          <p class="note">Select a gene row to inspect corpus-level evidence metrics. The full SQLite atlas database is generated next to this HTML file.</p>
        </div>
      </div>
    </section>
  </main>
  <script>
    const ATLAS = __ATLAS_JSON__;
    const cols = Object.fromEntries(ATLAS.gene_columns.map((name, index) => [name, index]));
    const rows = ATLAS.genes;
    const fmt = (value, digits=2) => value === "" || value === null || Number.isNaN(Number(value)) ? "" : Number(value).toFixed(digits);
    const pct100 = (value, digits=1) => {
      const numeric = Number(value);
      if (!Number.isFinite(numeric)) return "";
      return `${numeric.toFixed(digits)}%`;
    };
    const fraction100 = (value) => {
      const numeric = Number(value);
      if (!Number.isFinite(numeric)) return 0;
      return Math.max(0, Math.min(1, numeric / 100));
    };
    const esc = (value) => String(value ?? "").replace(/[&<>"']/g, ch => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[ch]));
    const cls = (value) => String(value ?? "").replace(/[^A-Za-z0-9_-]/g, "");
    const badge = (value) => `<span class="badge ${cls(value)}">${esc(value || "flat")}</span>`;
    const geneVal = (row, name) => row[cols[name]];

    function init() {
      document.getElementById("generated").textContent = `Generated ${ATLAS.generated_at.slice(0, 10)} · ${rows.length.toLocaleString()} scored gene rows`;
      document.getElementById("corpus").innerHTML = `<option value="">All corpora</option>` + ATLAS.corpora.map(c => `<option value="${esc(c.corpus_id)}">${esc(c.corpus_label)}</option>`).join("");
      document.getElementById("summary").innerHTML = [
        ["Corpora", ATLAS.corpora.length],
        ["Scored rows", rows.length.toLocaleString()],
        ["Methods/resources", ATLAS.prior_art.length],
        ["Best IFN top gene", (ATLAS.corpora.find(c => c.corpus_id === "ifn-pilot") || {}).top_gene || ""]
      ].map(([k, v]) => `<div class="metric"><span>${esc(k)}</span><strong>${esc(v)}</strong></div>`).join("");
      drawHeatmap();
      drawSupport();
      search();
    }

    function search() {
      const q = document.getElementById("query").value.trim().toUpperCase();
      const corpus = document.getElementById("corpus").value;
      const direction = document.getElementById("direction").value;
      let filtered = rows.filter(row => {
        const symbol = String(geneVal(row, "gene_symbol"));
        const dir = String(geneVal(row, "display_direction") || "").toLowerCase();
        if (q && !symbol.includes(q)) return false;
        if (corpus && geneVal(row, "corpus_id") !== corpus) return false;
        if (direction) {
          if (direction === "flat" && !["flat", "zero", ""].includes(dir)) return false;
          if (direction !== "flat" && dir !== direction) return false;
        }
        return true;
      });
      filtered.sort((a, b) => Number(geneVal(a, "display_rank")) - Number(geneVal(b, "display_rank")));
      const shown = filtered.slice(0, 120);
      document.getElementById("status").textContent = `${filtered.length.toLocaleString()} matches · ${shown.length} shown`;
      document.getElementById("genes").innerHTML = shown.map((row, i) => `
        <tr data-i="${i}">
          <td>${esc(geneVal(row, "corpus_label"))}</td>
          <td class="num">${fmt(geneVal(row, "display_rank"), 0)}</td>
          <td class="gene">${esc(geneVal(row, "gene_symbol"))}</td>
          <td>${badge(geneVal(row, "display_direction"))}</td>
          <td class="num">${fmt(geneVal(row, "display_top_percent"), 3)}</td>
          <td class="num">${fmt(geneVal(row, "n_source_units"), 0)}</td>
          <td class="num">${pct100(geneVal(row, "evidence_reliability_score"), 1)}</td>
        </tr>
      `).join("");
      Array.from(document.querySelectorAll("#genes tr")).forEach((tr, i) => tr.addEventListener("click", () => showDetail(shown[i])));
      if (shown.length) showDetail(shown[0]);
    }

    function showDetail(row) {
      const symbol = geneVal(row, "gene_symbol");
      const matches = rows.filter(r => geneVal(r, "gene_symbol") === symbol)
        .sort((a, b) => Number(geneVal(a, "display_rank")) - Number(geneVal(b, "display_rank")));
      document.getElementById("detail").innerHTML = `
        <h2>${esc(symbol)}</h2>
        <div class="detail-grid">
          <div class="metric"><span>Best corpus</span><strong>${esc(geneVal(row, "corpus_label"))}</strong></div>
          <div class="metric"><span>Best rank</span><strong>${fmt(geneVal(row, "display_rank"), 0)}</strong></div>
          <div class="metric"><span>Weighted LFC</span><strong>${fmt(geneVal(row, "weighted_lfc"), 2)}</strong></div>
          <div class="metric"><span>Reliability</span><strong>${pct100(geneVal(row, "evidence_reliability_score"), 1)}</strong></div>
          <div class="metric"><span>Direction confidence</span><strong>${fmt(Number(geneVal(row, "direction_confidence_index")) * 100, 1)}%</strong></div>
          <div class="metric"><span>LOO stability</span><strong>${fmt(Number(geneVal(row, "loo_rank_stability_score")) * 100, 1)}%</strong></div>
        </div>
        <table>
          <thead><tr><th>Corpus</th><th class="num">Rank</th><th>Direction</th><th class="num">Top %</th><th class="num">Units</th><th class="num">Score</th></tr></thead>
          <tbody>${matches.map(r => `<tr><td>${esc(geneVal(r, "corpus_label"))}</td><td class="num">${fmt(geneVal(r, "display_rank"), 0)}</td><td>${badge(geneVal(r, "display_direction"))}</td><td class="num">${fmt(geneVal(r, "display_top_percent"), 3)}</td><td class="num">${fmt(geneVal(r, "n_source_units"), 0)}</td><td class="num">${fmt(geneVal(r, "display_score"), 2)}</td></tr>`).join("")}</tbody>
        </table>
        <p class="note">Ranks and scores are relative prioritization indices, not calibrated probabilities. Detailed per-source evidence is stored in degora_atlas.db and each corpus-level degora_scores.db.</p>
      `;
    }

    function color(value) {
      const v = Math.max(0, Math.min(1, Number(value) || 0));
      const r = Math.round(245 - 210 * v);
      const g = Math.round(246 - 84 * v);
      const b = Math.round(248 - 118 * v);
      return `rgb(${r},${g},${b})`;
    }

    function drawHeatmap() {
      const corpora = ATLAS.corpora.map(c => c.corpus_label);
      const methods = [...new Set(ATLAS.methods.map(m => m.method_label))].filter(m => ["DEGORA","DEGORA quality","Weighted Stouffer","Fisher","AWFisher","metaRNASeq Fisher","MetaVolcanoR","RobustRankAggreg","Rank product approx."].includes(m));
      const cellW = 112, cellH = 28, left = 150, top = 42;
      const width = left + corpora.length * cellW + 20, height = top + methods.length * cellH + 36;
      const body = [];
      corpora.forEach((c, j) => body.push(`<text x="${left + j * cellW + cellW/2}" y="22" text-anchor="middle" font-size="11">${esc(c)}</text>`));
      methods.forEach((m, i) => {
        body.push(`<text x="8" y="${top + i * cellH + 18}" font-size="11">${esc(m)}</text>`);
        corpora.forEach((c, j) => {
          const rec = ATLAS.methods.find(x => x.corpus_label === c && x.method_label === m && x.run_status === "ok");
          const v = rec ? Number(rec.recall_at_100 || 0) : 0;
          body.push(`<rect x="${left + j * cellW}" y="${top + i * cellH}" width="${cellW-4}" height="${cellH-4}" rx="4" fill="${color(v)}"></rect>`);
          body.push(`<text x="${left + j * cellW + cellW/2}" y="${top + i * cellH + 17}" text-anchor="middle" font-size="11" fill="#17202a">${rec ? fmt(v, 2) : ""}</text>`);
        });
      });
      document.getElementById("heatmap").innerHTML = `<svg viewBox="0 0 ${width} ${height}" role="img">${body.join("")}</svg>`;
    }

    function drawSupport() {
      const data = ATLAS.corpora.map(c => {
        const topRows = rows.filter(r => geneVal(r, "corpus_id") === c.corpus_id && Number(geneVal(r, "display_rank")) <= 100);
        const med = (arr) => {
          const v = arr.map(Number).filter(Number.isFinite).sort((a,b)=>a-b);
          return v.length ? v[Math.floor(v.length/2)] : 0;
        };
        return {label: c.corpus_label, units: med(topRows.map(r => geneVal(r, "n_source_units"))), rel: fraction100(med(topRows.map(r => geneVal(r, "evidence_reliability_score"))))};
      });
      const width = 760, height = 230, left = 92, bottom = 38, top = 16;
      const maxUnits = Math.max(1, ...data.map(d => d.units));
      const barW = 36, gap = 58;
      const body = [`<line x1="${left}" y1="${height-bottom}" x2="${width-10}" y2="${height-bottom}" stroke="#cfd6cc"></line>`];
      data.forEach((d, i) => {
        const x = left + i * (barW + gap);
        const h = (height - bottom - top) * d.units / maxUnits;
        const y = height - bottom - h;
        body.push(`<rect x="${x}" y="${y}" width="${barW}" height="${h}" rx="4" fill="#265d97"></rect>`);
        body.push(`<circle cx="${x + barW/2}" cy="${height - bottom - (height-bottom-top)*d.rel}" r="5" fill="#b4234b"></circle>`);
        body.push(`<text x="${x + barW/2}" y="${height-14}" text-anchor="middle" font-size="11">${esc(d.label)}</text>`);
      });
      body.push(`<text x="8" y="22" font-size="11" fill="#59656f">bars: median units</text>`);
      body.push(`<text x="8" y="38" font-size="11" fill="#59656f">red: reliability (0-100%)</text>`);
      document.getElementById("supportPlot").innerHTML = `<svg viewBox="0 0 ${width} ${height}" role="img">${body.join("")}</svg>`;
    }

    document.getElementById("search").addEventListener("click", search);
    document.getElementById("query").addEventListener("keydown", e => { if (e.key === "Enter") search(); });
    document.getElementById("corpus").addEventListener("change", search);
    document.getElementById("direction").addEventListener("change", search);
    init();
  </script>
</body>
</html>
"""


def _write_html(path: Path, payload: dict[str, Any]) -> None:
    encoded = (
        json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )
    path.write_text(HTML_TEMPLATE.replace("__ATLAS_JSON__", encoded), encoding="utf-8")


def _write_sqlite(path: Path, tables: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    with sqlite3.connect(path) as connection:
        for table, frame in tables.items():
            frame.to_sql(table, connection, index=False)
        connection.execute("CREATE INDEX idx_gene_scores_symbol ON gene_scores(gene_symbol)")
        connection.execute("CREATE INDEX idx_gene_scores_corpus_rank ON gene_scores(corpus_id, display_rank)")
        connection.execute("CREATE INDEX idx_method_recall_corpus_method ON method_recall(corpus_id, method_id)")


def _ordered_recall_heatmap(method_recall: pd.DataFrame) -> pd.DataFrame:
    frame = method_recall.loc[method_recall["run_status"].astype(str).eq("ok")].copy()
    frame = frame.loc[frame["method_label"].isin(METHOD_ORDER)]
    frame["method_label"] = pd.Categorical(frame["method_label"], categories=METHOD_ORDER, ordered=True)
    frame["recall_at_100"] = pd.to_numeric(frame["recall_at_100"], errors="coerce")
    return frame.sort_values(["method_label", "corpus_label"]).reset_index(drop=True)


def _set_plot_style() -> None:
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


def _save_figure(fig: plt.Figure, stem: Path) -> list[Path]:
    outputs = []
    for suffix in (".png", ".pdf", ".svg"):
        path = stem.with_suffix(suffix)
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def _plot_panel_a(corpus_summary: pd.DataFrame, stem: Path) -> list[Path]:
    _set_plot_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    data = corpus_summary.copy()
    x = np.arange(len(data))
    ax.bar(x - 0.18, data["n_genes"] / 1000, width=0.36, color="#265d97", label="Scored genes (thousand)")
    ax2 = ax.twinx()
    ax2.bar(x + 0.18, data["n_source_units_total"], width=0.36, color="#0f766e", label="Source units")
    ax.set_xticks(x)
    ax.set_xticklabels(data["corpus_label"], rotation=35, ha="right")
    ax.set_ylabel("Scored genes (thousand)")
    ax2.set_ylabel("Independent source units")
    ax.set_title("A. Atlas scale", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#e5e7eb")
    ax.set_axisbelow(True)
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles + handles2, labels + labels2, loc="upper left", frameon=False)
    return _save_figure(fig, stem)


def _plot_panel_b(recall_heatmap: pd.DataFrame, corpus_summary: pd.DataFrame, stem: Path) -> list[Path]:
    _set_plot_style()
    corpora = corpus_summary["corpus_label"].tolist()
    methods = [method for method in METHOD_ORDER if method in set(recall_heatmap["method_label"].astype(str))]
    matrix = np.full((len(methods), len(corpora)), np.nan)
    for i, method in enumerate(methods):
        for j, corpus in enumerate(corpora):
            sub = recall_heatmap.loc[
                recall_heatmap["method_label"].astype(str).eq(method) & recall_heatmap["corpus_label"].eq(corpus)
            ]
            if not sub.empty:
                matrix[i, j] = float(sub["recall_at_100"].iloc[0])
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    cmap = LinearSegmentedColormap.from_list("degora_recall", ["#f8fafc", "#a7d8c9", "#0f766e"])
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(corpora)))
    ax.set_xticklabels(corpora, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(corpora)):
            value = matrix[i, j]
            label = "" if np.isnan(value) else f"{value:.2f}"
            ax.text(j, i, label, ha="center", va="center", color="#17202a", fontsize=7.5)
    ax.set_title("B. Recall@100 by method", loc="left", fontweight="bold")
    fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="Recall@100")
    return _save_figure(fig, stem)


def _plot_panel_c(top100: pd.DataFrame, stem: Path) -> list[Path]:
    _set_plot_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.scatter(
        top100["median_source_units"],
        top100["median_evidence_reliability"],
        s=80,
        c=top100["median_sign_concordance"],
        cmap="viridis",
        vmin=0,
        vmax=1,
        edgecolor="#1f2937",
        linewidth=0.4,
    )
    for row in top100.itertuples(index=False):
        ax.text(row.median_source_units + 0.04, row.median_evidence_reliability, row.corpus_label, fontsize=7.5)
    ax.set_xlabel("Median source units among top 100")
    ax.set_ylabel("Median evidence reliability (%)")
    ax.set_ylim(0, 105)
    ax.grid(True, color="#e5e7eb")
    ax.set_axisbelow(True)
    ax.set_title("C. Top-100 evidence quality", loc="left", fontweight="bold")
    return _save_figure(fig, stem)


def _plot_panel_d(prior_counts: pd.DataFrame, stem: Path) -> list[Path]:
    _set_plot_style()
    counts = prior_counts.groupby("public_summary_status")["n_methods_or_resources"].sum().sort_values()
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    colors = ["#6b7280", "#9a6a12", "#265d97", "#0f766e", "#b4234b", "#7c3aed"][: len(counts)]
    ax.barh(np.arange(len(counts)), counts.values, color=colors)
    ax.set_yticks(np.arange(len(counts)))
    ax.set_yticklabels(counts.index)
    ax.set_xlabel("Methods/resources")
    ax.set_title("D. Prior-art input status", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb")
    ax.set_axisbelow(True)
    return _save_figure(fig, stem)


def _plot_combined(
    corpus_summary: pd.DataFrame,
    recall_heatmap: pd.DataFrame,
    top100: pd.DataFrame,
    prior_counts: pd.DataFrame,
    stem: Path,
) -> list[Path]:
    _set_plot_style()
    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.08], height_ratios=[1, 1], wspace=0.36, hspace=0.45)

    ax_a = fig.add_subplot(gs[0, 0])
    data = corpus_summary.copy()
    x = np.arange(len(data))
    ax_a.bar(x, data["n_genes"] / 1000, color="#265d97")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(data["corpus_label"], rotation=35, ha="right")
    ax_a.set_ylabel("Scored genes (thousand)")
    ax_a.set_title("A. Atlas scale", loc="left", fontweight="bold")
    ax_a.grid(axis="y", color="#e5e7eb")

    ax_b = fig.add_subplot(gs[0, 1])
    corpora = corpus_summary["corpus_label"].tolist()
    methods = [method for method in METHOD_ORDER if method in set(recall_heatmap["method_label"].astype(str))]
    matrix = np.full((len(methods), len(corpora)), np.nan)
    for i, method in enumerate(methods):
        for j, corpus in enumerate(corpora):
            sub = recall_heatmap.loc[
                recall_heatmap["method_label"].astype(str).eq(method) & recall_heatmap["corpus_label"].eq(corpus)
            ]
            if not sub.empty:
                matrix[i, j] = float(sub["recall_at_100"].iloc[0])
    cmap = LinearSegmentedColormap.from_list("degora_recall", ["#f8fafc", "#a7d8c9", "#0f766e"])
    ax_b.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax_b.set_xticks(np.arange(len(corpora)))
    ax_b.set_xticklabels(corpora, rotation=35, ha="right")
    ax_b.set_yticks(np.arange(len(methods)))
    ax_b.set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(corpora)):
            value = matrix[i, j]
            ax_b.text(j, i, "" if np.isnan(value) else f"{value:.2f}", ha="center", va="center", fontsize=6.8)
    ax_b.set_title("B. Recall@100", loc="left", fontweight="bold")

    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.scatter(
        top100["median_source_units"],
        top100["median_evidence_reliability"],
        s=70,
        c=top100["median_sign_concordance"],
        cmap="viridis",
        vmin=0,
        vmax=1,
        edgecolor="#1f2937",
        linewidth=0.4,
    )
    for row in top100.itertuples(index=False):
        ax_c.text(row.median_source_units + 0.04, row.median_evidence_reliability, row.corpus_label, fontsize=7)
    ax_c.set_xlabel("Median source units in top 100")
    ax_c.set_ylabel("Median reliability (%)")
    ax_c.set_ylim(0, 105)
    ax_c.set_title("C. Evidence quality", loc="left", fontweight="bold")
    ax_c.grid(True, color="#e5e7eb")

    ax_d = fig.add_subplot(gs[1, 1])
    counts = prior_counts.groupby("public_summary_status")["n_methods_or_resources"].sum().sort_values()
    ax_d.barh(np.arange(len(counts)), counts.values, color="#9a6a12")
    ax_d.set_yticks(np.arange(len(counts)))
    ax_d.set_yticklabels(counts.index)
    ax_d.set_xlabel("Methods/resources")
    ax_d.set_title("D. Prior-art input status", loc="left", fontweight="bold")
    ax_d.grid(axis="x", color="#e5e7eb")

    fig.suptitle("DEGORA searchable evidence atlas", x=0.01, y=0.995, ha="left", fontsize=12, fontweight="bold")
    return _save_figure(fig, stem)


def _docx_escape(text: str) -> str:
    return escape(text).replace("\n", "</w:t><w:br/><w:t>")


def _write_minimal_docx(path: Path, paragraphs: list[str]) -> None:
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
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", rels)
        archive.writestr("word/document.xml", document)


def _write_xlsx(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    workbook = Workbook()
    first = True
    for name, frame in sheets.items():
        worksheet = workbook.active if first else workbook.create_sheet()
        first = False
        worksheet.title = name[:31]
        for row in dataframe_to_rows(frame, index=False, header=True):
            worksheet.append(row)
    workbook.save(path)


def _write_manifest(
    path: Path,
    *,
    command: str,
    inputs: list[Path],
    outputs: list[Path],
    source_data: list[Path],
    corpus_summary: pd.DataFrame,
    top100: pd.DataFrame,
) -> None:
    manifest = {
        "figure_id": FIGURE_ID,
        "generated_at": datetime.now(UTC).isoformat(),
        "script": "outputs/code/figures/make_degora_atlas_dashboard.py",
        "command": command,
        "inputs": [str(path) for path in inputs],
        "outputs": [str(path) for path in outputs],
        "source_data": [str(path) for path in source_data],
        "panel_claims": [
            "A, corpus scale is computed from degora_gene_scores.csv and degora_score_metadata.json.",
            "B, recall@100 is read from corpus-level gold or anchor comparator summaries.",
            "C, top-100 evidence quality is computed from quality-weighted DEGORA ranks.",
            "D, prior-art input status is summarized from prior_art_coverage_summary.csv.",
        ],
        "transformations": {
            "gene_rank_used_for_dashboard": "quality_weighted_degora_rank when present, otherwise degora_rank",
            "top100_evidence": "median metrics among genes with dashboard rank <= 100 per corpus",
            "static_html": "self-contained HTML with embedded compact score index; detailed evidence remains in SQLite",
        },
        "validation": {
            "n_corpora": int(len(corpus_summary)),
            "n_total_scored_gene_rows": int(corpus_summary["n_genes"].sum()),
            "top100_rows": top100.to_dict(orient="records"),
        },
        "known_limitations": [
            "The static HTML embeds score-level rows, not every source-level evidence row.",
            "Recall summaries use locked or anchor panels and do not estimate false-positive rates.",
            "DEGORA scores are relative prioritization indices, not calibrated probabilities.",
        ],
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def _write_validation(path: Path, outputs: list[Path], tables: dict[str, pd.DataFrame]) -> None:
    lines = ["DEGORA atlas validation", ""]
    for artifact in outputs:
        lines.append(f"{artifact}: exists={artifact.exists()} size_bytes={artifact.stat().st_size if artifact.exists() else 0}")
    lines.append("")
    for name, frame in tables.items():
        lines.append(f"{name}: rows={len(frame)} columns={len(frame.columns)}")
    path.write_text("\n".join(lines) + "\n")


def build_atlas_package(
    *,
    corpora: list[CorpusInput],
    prior_art_path: Path,
    output_dir: Path,
    command: str,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    missing = [str(corpus.score_csv) for corpus in corpora if not corpus.score_csv.exists()]
    if missing:
        raise FileNotFoundError(f"missing score CSV files: {missing}")
    if not prior_art_path.exists():
        raise FileNotFoundError(f"missing prior-art summary: {prior_art_path}")

    gene_index = pd.concat([_read_gene_scores(corpus) for corpus in corpora], ignore_index=True)
    method_recall = pd.concat([_read_method_recall(corpus) for corpus in corpora], ignore_index=True)
    prior_art = pd.read_csv(prior_art_path)
    corpus_summary = _corpus_summary(corpora, gene_index, method_recall)
    top100 = _top100_evidence(gene_index)
    prior_counts = _prior_art_counts(prior_art)
    recall_heatmap = _ordered_recall_heatmap(method_recall)

    html_path = output_dir / "degora_atlas.html"
    sqlite_path = output_dir / "degora_atlas.db"
    source_xlsx = output_dir / "degora_atlas_source_data.xlsx"
    gene_index_csv = output_dir / "degora_atlas_gene_index.csv"
    legend_docx = output_dir / "degora_atlas_legend.docx"
    manifest_json = output_dir / "degora_atlas_manifest.json"
    validation_txt = output_dir / "degora_atlas_validation.txt"

    _write_html(
        html_path,
        _html_payload(
            corpus_summary=corpus_summary,
            method_recall=method_recall,
            prior_art=prior_art,
            gene_index=gene_index,
        ),
    )
    _write_sqlite(
        sqlite_path,
        {
            "gene_scores": gene_index,
            "corpus_summary": corpus_summary,
            "method_recall": method_recall,
            "prior_art_coverage": prior_art,
            "prior_art_counts": prior_counts,
        },
    )
    gene_index.to_csv(gene_index_csv, index=False)
    _write_xlsx(
        source_xlsx,
        {
            "corpus_summary": corpus_summary,
            "recall_heatmap": recall_heatmap,
            "top100_evidence": top100,
            "prior_art_counts": prior_counts,
            "dashboard_gene_index": gene_index,
        },
    )

    figure_outputs: list[Path] = []
    figure_outputs.extend(_plot_combined(corpus_summary, recall_heatmap, top100, prior_counts, output_dir / OUTPUT_STEM))
    figure_outputs.extend(_plot_panel_a(corpus_summary, output_dir / "degora_atlas_A_corpus_scale"))
    figure_outputs.extend(_plot_panel_b(recall_heatmap, corpus_summary, output_dir / "degora_atlas_B_recall_heatmap"))
    figure_outputs.extend(_plot_panel_c(top100, output_dir / "degora_atlas_C_evidence_quality"))
    figure_outputs.extend(_plot_panel_d(prior_counts, output_dir / "degora_atlas_D_prior_art_status"))

    _write_minimal_docx(
        legend_docx,
        [
            "DEGORA_ATLAS. Searchable DEGORA evidence atlas.",
            "A, scored genes and independent source units across the included DEGORA corpora. B, recall@100 for runnable ranking methods against locked or anchor panels. C, median source support and evidence reliability among the top 100 quality-weighted DEGORA genes; point color encodes median sign concordance. D, method/resource counts by public-summary input status.",
        ],
    )

    outputs = [html_path, sqlite_path, source_xlsx, gene_index_csv, legend_docx, manifest_json, validation_txt, *figure_outputs]
    source_data = [source_xlsx, gene_index_csv]
    inputs = [
        prior_art_path,
        *[corpus.score_csv for corpus in corpora],
        *[corpus.comparator_summary for corpus in corpora if corpus.comparator_summary is not None],
        *[corpus.metadata_json for corpus in corpora if corpus.metadata_json.exists()],
    ]
    _write_validation(
        validation_txt,
        outputs=[html_path, sqlite_path, source_xlsx, gene_index_csv, legend_docx, *figure_outputs],
        tables={
            "gene_index": gene_index,
            "corpus_summary": corpus_summary,
            "method_recall": method_recall,
            "prior_art": prior_art,
            "top100": top100,
            "prior_counts": prior_counts,
        },
    )
    _write_manifest(
        manifest_json,
        command=command,
        inputs=inputs,
        outputs=outputs,
        source_data=source_data,
        corpus_summary=corpus_summary,
        top100=top100,
    )
    for artifact in outputs:
        write_source_sidecar(
            artifact,
            command,
            inputs=inputs,
            metadata={"generator": "degora-atlas-dashboard", "figure_id": FIGURE_ID},
        )
    return {
        "figure_id": FIGURE_ID,
        "n_corpora": int(len(corpus_summary)),
        "n_gene_rows": int(len(gene_index)),
        "n_method_rows": int(len(method_recall)),
        "n_prior_art_rows": int(len(prior_art)),
        "outputs": [str(path) for path in outputs],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/figures/manuscript/DEGORA_ATLAS"))
    parser.add_argument("--prior-art", type=Path, default=Path("outputs/results/prior_art_coverage_summary.csv"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = shell_command(
        [
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/figures/make_degora_atlas_dashboard.py",
            "--output-dir",
            args.output_dir,
            "--prior-art",
            args.prior_art,
        ]
    )
    summary = build_atlas_package(
        corpora=DEFAULT_CORPORA,
        prior_art_path=args.prior_art,
        output_dir=args.output_dir,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
