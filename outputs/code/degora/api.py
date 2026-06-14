"""Local HTTP API and browser UI for a DEGORA SQLite score database."""

from __future__ import annotations

import errno
import json
import math
import sqlite3
import sys
from contextlib import closing
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, quote, unquote, urlparse


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>DEGORA</title>
  <link rel="icon" href="data:,">
  <style>
    :root {
      color-scheme: light;
      --bg: #f6f7f5;
      --panel: #ffffff;
      --line: #d9ded8;
      --ink: #1f2933;
      --muted: #64706b;
      --accent: #0f766e;
      --warn: #a16207;
      --danger: #9f1239;
    }
    * { box-sizing: border-box; }
    html { height: 100%; }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 14px;
      height: 100vh;
      overflow: hidden;
      display: grid;
      grid-template-rows: auto minmax(0, 1fr) auto;
    }
    header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      padding: 18px 24px 12px;
      border-bottom: 1px solid var(--line);
      background: var(--panel);
      position: sticky;
      top: 0;
      z-index: 2;
    }
    h1 {
      margin: 0;
      font-size: 20px;
      font-weight: 700;
      letter-spacing: 0;
    }
    .meta {
      display: flex;
      gap: 16px;
      flex-wrap: wrap;
      color: var(--muted);
      justify-content: flex-end;
    }
    main {
      display: grid;
      grid-template-columns: minmax(320px, 60%) 18px minmax(320px, 1fr);
      gap: 6px;
      padding: 16px 24px 24px;
      flex: 1;
      min-height: 0;
      height: 100%;
      overflow: hidden;
    }
    section {
      min-width: 0;
      min-height: 0;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
    }
    .genes-panel, .evidence-panel {
      display: flex;
      flex-direction: column;
      height: 100%;
    }
    .splitter {
      align-self: stretch;
      width: 18px;
      min-width: 18px;
      border: 0;
      border-radius: 8px;
      background: transparent;
      cursor: col-resize;
      position: relative;
      touch-action: none;
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .splitter::before {
      content: "";
      position: absolute;
      inset: 8px 7px;
      border-radius: 999px;
      background: #b9c8c0;
      transition: background 0.12s ease, inset 0.12s ease;
    }
    .splitter::after {
      content: "";
      width: 4px;
      height: 28px;
      border-radius: 999px;
      background:
        radial-gradient(circle, #5a6f64 1.6px, transparent 1.8px) 0 0 / 4px 7px repeat-y;
      position: relative;
      z-index: 1;
      opacity: 0.9;
      transition: opacity 0.12s ease, transform 0.12s ease;
    }
    .splitter:hover::before,
    .splitter:focus-visible::before,
    body.is-resizing .splitter::before {
      inset: 8px 5px;
      background: var(--accent);
    }
    .splitter:hover::after,
    .splitter:focus-visible::after,
    body.is-resizing .splitter::after {
      opacity: 1;
      transform: scale(1.08);
    }
    .splitter:focus-visible {
      outline: 2px solid var(--accent);
      outline-offset: 2px;
    }
    body.is-resizing {
      cursor: col-resize;
      user-select: none;
    }
    .section-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
    }
    h2 {
      margin: 0;
      font-size: 15px;
      font-weight: 700;
      letter-spacing: 0;
    }
    .controls {
      display: grid;
      grid-template-columns: minmax(150px, 1fr) 120px 120px 92px;
      gap: 8px;
      width: 100%;
      padding: 12px 14px;
      border-bottom: 1px solid var(--line);
      background: #fbfcfb;
    }
    input, select, button {
      width: 100%;
      height: 34px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      color: var(--ink);
      font: inherit;
      padding: 0 10px;
    }
    button {
      background: var(--accent);
      border-color: var(--accent);
      color: #fff;
      font-weight: 650;
      cursor: pointer;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      table-layout: fixed;
    }
    th, td {
      padding: 9px 10px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: middle;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
    th {
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
      background: #fbfcfb;
    }
    .gene-table-shell {
      flex: 1 1 0;
      min-height: 0;
      position: relative;
    }
    .gene-table-scroll {
      height: 100%;
      overflow: auto;
      overscroll-behavior: contain;
    }
    .gene-table-scroll table { min-width: 760px; }
    .gene-rank-col { width: 58px; }
    .gene-tier-col { width: 54px; }
    .gene-symbol-col { width: 160px; }
    .gene-top-col { width: 72px; }
    .gene-score-col { width: 78px; }
    .gene-units-col { width: 62px; }
    .gene-sign-col { width: 82px; }
    .gene-lfc-col { width: 72px; }
    .gene-table-scroll thead th {
      position: sticky;
      top: 0;
      z-index: 1;
    }
    th.sortable { padding: 0; }
    .sort-head {
      display: flex;
      align-items: center;
      gap: 5px;
      width: 100%;
      min-height: 36px;
      height: auto;
      padding: 9px 10px;
      border: 0;
      border-radius: 0;
      background: transparent;
      color: inherit;
      font: inherit;
      font-weight: 700;
      text-align: left;
      cursor: pointer;
    }
    th.num .sort-head { justify-content: flex-end; text-align: right; }
    .sort-head:hover, .sort-head:focus-visible { background: #eef7f5; outline: none; }
    .sort-head:disabled {
      cursor: wait;
      color: #8a9b98;
    }
    .sort-indicator {
      display: inline-block;
      width: 10px;
      color: var(--accent);
      font-weight: 900;
    }
    .table-loading {
      position: absolute;
      inset: 0;
      z-index: 3;
      display: flex;
      align-items: center;
      justify-content: center;
      background: rgba(248, 250, 249, 0.74);
      backdrop-filter: blur(1px);
    }
    .table-loading[hidden] { display: none; }
    .loading-card {
      width: min(320px, calc(100% - 32px));
      padding: 14px 16px;
      border: 1px solid #c9ded9;
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 16px 36px rgba(21, 50, 45, 0.14);
    }
    .loading-title {
      display: block;
      margin-bottom: 8px;
      color: var(--ink);
      font-weight: 800;
    }
    .loading-note {
      display: block;
      margin-top: 8px;
      color: var(--muted);
      font-size: 12px;
    }
    .loading-bar {
      height: 4px;
      overflow: hidden;
      border-radius: 999px;
      background: #dfe9e6;
    }
    .loading-bar::before {
      content: "";
      display: block;
      width: 46%;
      height: 100%;
      border-radius: inherit;
      background: var(--accent);
      animation: loading-slide 1s ease-in-out infinite;
    }
    .genes-panel.is-loading .gene-table-scroll table { opacity: 0.45; }
    .genes-panel.is-loading .gene-table-scroll { cursor: wait; }
    @keyframes loading-slide {
      0% { transform: translateX(-115%); }
      50% { transform: translateX(85%); }
      100% { transform: translateX(250%); }
    }
    tbody tr { cursor: pointer; }
    tbody tr:hover { background: #eef7f5; }
    #genes tr.selected { background: #d4ece7; box-shadow: inset 3px 0 0 var(--accent); }
    #genes tr:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
    footer.legend {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 12px 22px;
      padding: 0 24px 22px;
      color: var(--muted);
      font-size: 12px;
    }
    footer.legend b { color: var(--ink); font-weight: 700; margin-right: 4px; }
    footer.legend .tier, footer.legend .badge { margin: 0 2px; }
    .more { padding: 10px 14px; border-top: 1px solid var(--line); text-align: center; }
    .more button { width: auto; min-width: 160px; }
    .more button[hidden] { display: none; }
    .num { text-align: right; font-variant-numeric: tabular-nums; }
    .gene { font-weight: 750; }
    .badge {
      display: inline-flex;
      align-items: center;
      min-width: 44px;
      justify-content: center;
      height: 22px;
      padding: 0 8px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
      background: #e7f5f2;
      color: #0f5d56;
    }
    .badge.down { background: #fde8ee; color: var(--danger); }
    .badge.flat { background: #f4f2e8; color: var(--warn); }
    .tier {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      width: 26px;
      height: 24px;
      border-radius: 6px;
      font-weight: 800;
      background: #eef2ff;
      color: #3730a3;
    }
    .tier.B { background: #e7f5f2; color: #0f5d56; }
    .tier.C { background: #f4f2e8; color: var(--warn); }
    .tier.D { background: #f1f5f9; color: var(--muted); }
    .detail-body { padding: 14px; }
    .evidence-panel .detail-body {
      flex: 1 1 0;
      min-height: 0;
      overflow: auto;
    }
    .kv {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 10px;
      margin-bottom: 14px;
    }
    .metric {
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 9px 10px;
      min-height: 58px;
    }
    .metric span {
      display: block;
      color: var(--muted);
      font-size: 12px;
      margin-bottom: 4px;
    }
    .metric strong {
      display: block;
      font-size: 18px;
      font-variant-numeric: tabular-nums;
      overflow-wrap: anywhere;
    }
    .sources {
      color: var(--muted);
      line-height: 1.45;
      overflow-wrap: anywhere;
      margin: 0 0 14px;
    }
    .empty {
      padding: 22px 14px;
      color: var(--muted);
    }
    .status {
      color: var(--muted);
      font-variant-numeric: tabular-nums;
    }
    .hint { color: var(--muted); font-size: 12px; }
    .metric[data-tip], th[data-tip], td[data-tip] { cursor: help; }
    .ev-scroll { max-height: 45vh; overflow: auto; border: 1px solid var(--line); border-radius: 6px; }
    .ev-scroll table { border-collapse: separate; border-spacing: 0; }
    .ev-scroll thead th { position: sticky; top: 0; z-index: 1; }
    th[data-tip] { text-decoration: underline dotted; text-underline-offset: 3px; }
    .metric[data-tip] span { text-decoration: underline dotted; text-underline-offset: 2px; }
    #tip {
      position: fixed;
      z-index: 50;
      max-width: 320px;
      padding: 8px 10px;
      border-radius: 6px;
      background: #1f2933;
      color: #fff;
      font-size: 12.5px;
      line-height: 1.45;
      box-shadow: 0 6px 24px rgba(0, 0, 0, .18);
      pointer-events: none;
      display: none;
    }
    @media (max-width: 980px) {
      header { align-items: flex-start; flex-direction: column; }
      .meta { justify-content: flex-start; }
      body { height: auto; min-height: 100vh; overflow: auto; }
      main { grid-template-columns: 1fr; padding: 12px; overflow: visible; min-height: auto; }
      .splitter { display: none; }
      section { min-height: auto; height: auto; }
      .gene-table-shell { flex: none; height: 60vh; max-height: 60vh; }
      .evidence-panel .detail-body { overflow: visible; }
      .controls { grid-template-columns: 1fr 1fr; }
      .kv { grid-template-columns: 1fr 1fr; }
    }
    @media (max-width: 580px) {
      .controls { grid-template-columns: 1fr; }
      .kv { grid-template-columns: 1fr; }
      th:nth-child(4), td:nth-child(4),
      th:nth-child(6), td:nth-child(6) { display: none; }
    }
  </style>
</head>
<body>
  <header>
    <h1>DEGORA</h1>
    <div class="meta" id="meta"></div>
  </header>
  <main id="layoutMain">
    <section class="genes-panel">
      <div class="section-head">
        <h2>Genes</h2>
        <div class="status" id="status" role="status" aria-live="polite"></div>
      </div>
      <div class="controls">
        <input id="query" placeholder="Gene symbol" autocomplete="off" aria-label="Gene symbol search" data-tip="Type part of a gene symbol to filter the list; leave blank to show all genes.">
        <input id="minUnits" type="number" min="1" value="1" aria-label="Min source units" data-tip="Show only genes supported by at least this many independent source units (one paper = one unit).">
        <select id="direction" aria-label="Direction" data-tip="Filter by consensus regulation direction: up, down, flat, or all.">
          <option value="">All directions</option>
          <option value="up">Up</option>
          <option value="down">Down</option>
          <option value="flat">Flat</option>
        </select>
        <button id="load" data-tip="Apply the filters above and reload the gene list.">Search</button>
      </div>
      <div class="gene-table-shell">
        <div class="gene-table-scroll" id="geneTableScroll">
          <table>
            <colgroup>
              <col class="gene-rank-col">
              <col class="gene-tier-col">
              <col class="gene-symbol-col">
              <col class="gene-top-col">
              <col class="gene-score-col">
              <col class="gene-units-col">
              <col class="gene-sign-col">
              <col class="gene-lfc-col">
            </colgroup>
            <thead>
              <tr>
                <th class="sortable" data-tip="DEGORA rank among all scored genes (1 = strongest)." aria-sort="ascending"><button class="sort-head" type="button" data-sort="rank">Rank <span class="sort-indicator" aria-hidden="true">^</span></button></th>
                <th class="sortable" data-tip="Confidence tier from rank, support, and direction (A strongest to D weakest)." aria-sort="none"><button class="sort-head" type="button" data-sort="tier">Tier <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="sortable" data-tip="Gene symbol." aria-sort="none"><button class="sort-head" type="button" data-sort="gene">Gene <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Position as a percent of all scored genes (e.g. top 1%)." aria-sort="none"><button class="sort-head" type="button" data-sort="top">Top <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="DEGORA quality-weighted prioritization score: a relative index, not a probability." aria-sort="none"><button class="sort-head" type="button" data-sort="score">Score <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Number of independent source units supporting the gene (one paper = one unit)." aria-sort="none"><button class="sort-head" type="button" data-sort="units">Units <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Direction concordance: percent of supporting evidence agreeing on the consensus direction." aria-sort="none"><button class="sort-head" type="button" data-sort="sign">Sign <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Sample-size-weighted mean log2 fold-change across supporting source units." aria-sort="none"><button class="sort-head" type="button" data-sort="lfc">LFC <span class="sort-indicator" aria-hidden="true"></span></button></th>
              </tr>
            </thead>
            <tbody id="genes"></tbody>
          </table>
        </div>
        <div class="table-loading" id="tableLoading" role="status" aria-live="polite" hidden>
          <div class="loading-card">
            <strong class="loading-title" id="loadingTitle">Loading genes...</strong>
            <div class="loading-bar" aria-hidden="true"></div>
            <span class="loading-note" id="loadingNote">Refreshing the table</span>
          </div>
        </div>
      </div>
      <div class="more"><button id="loadMore" type="button" hidden>Load more</button></div>
    </section>
    <div
      class="splitter"
      id="layoutSplitter"
      role="separator"
      aria-orientation="vertical"
      aria-label="Resize gene and evidence panels"
      aria-valuemin="25"
      aria-valuemax="75"
      aria-valuenow="61"
      tabindex="0"
      data-tip="Drag left or right to resize the gene and evidence panels. Double-click to reset."
    ></div>
    <section class="evidence-panel">
      <div class="section-head">
        <h2>Evidence</h2>
        <span class="hint">Hover any label or column for its meaning</span>
      </div>
      <div class="detail-body" id="detail">
        <div class="empty">Select a gene.</div>
      </div>
    </section>
  </main>
  <footer class="legend">
    <span><b>Tier</b> <span class="tier">A</span><span class="tier B">B</span><span class="tier C">C</span><span class="tier D">D</span> A strongest to D weakest</span>
    <span><b>Direction</b> <span class="badge up">up</span><span class="badge down">down</span><span class="badge flat">flat</span></span>
    <span>Hover any label, column header, or truncated cell for its meaning.</span>
  </footer>
  <div id="tip" role="tooltip"></div>
  <script>
    const $ = (id) => document.getElementById(id);
    const fmt = (value, digits = 3) => {
      if (value === null || value === undefined || Number.isNaN(Number(value))) return "";
      const n = Number(value);
      if (Math.abs(n) >= 100) return n.toFixed(1);
      if (Math.abs(n) < 0.001 && n !== 0) return n.toExponential(2);
      return n.toFixed(digits);
    };
    const esc = (value) => String(value ?? "").replace(/[&<>"']/g, (ch) => ({
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#39;"
    }[ch]));
    const safeClass = (value) => String(value ?? "").replace(/[^A-Za-z0-9_-]/g, "");
    const badge = (direction) => `<span class="badge ${safeClass(direction)}">${esc(direction)}</span>`;
    const tier = (value) => `<span class="tier ${safeClass(value)}">${esc(value)}</span>`;
    const SPLIT_STORAGE_KEY = "degoraPanelSplitPercentV2";
    const DEFAULT_SPLIT_PERCENT = 60;

    function splitBounds() {
      const main = $("layoutMain");
      const width = main ? main.getBoundingClientRect().width : 0;
      const minFromPixels = width > 0 ? Math.min(42, Math.max(25, (320 / width) * 100)) : 25;
      const maxFromPixels = width > 0 ? Math.max(58, Math.min(75, 100 - (320 / width) * 100)) : 75;
      if (minFromPixels >= maxFromPixels) return { min: 35, max: 65 };
      return { min: minFromPixels, max: maxFromPixels };
    }

    function storedSplitPercent() {
      try {
        const stored = window.localStorage.getItem(SPLIT_STORAGE_KEY);
        if (stored === null || stored === "") return DEFAULT_SPLIT_PERCENT;
        const value = Number(stored);
        return Number.isFinite(value) ? value : DEFAULT_SPLIT_PERCENT;
      } catch (_) {
        return DEFAULT_SPLIT_PERCENT;
      }
    }

    function setPanelSplit(percent, persist = true) {
      const main = $("layoutMain");
      const splitter = $("layoutSplitter");
      if (!main || !splitter) return;
      if (window.matchMedia("(max-width: 980px)").matches) {
        main.style.gridTemplateColumns = "";
        return;
      }
      const bounds = splitBounds();
      const clamped = Math.min(bounds.max, Math.max(bounds.min, percent));
      main.style.gridTemplateColumns = `minmax(320px, ${clamped}%) 18px minmax(320px, 1fr)`;
      splitter.setAttribute("aria-valuenow", String(Math.round(clamped)));
      if (persist) {
        try { window.localStorage.setItem(SPLIT_STORAGE_KEY, String(clamped)); } catch (_) {}
      }
    }

    function initPanelResize() {
      const main = $("layoutMain");
      const splitter = $("layoutSplitter");
      if (!main || !splitter) return;
      let activePointerId = null;
      let activeMouse = false;

      const applyClientX = (clientX) => {
        const rect = main.getBoundingClientRect();
        if (!rect.width) return;
        setPanelSplit(((clientX - rect.left) / rect.width) * 100);
      };
      const stopResize = () => {
        activePointerId = null;
        activeMouse = false;
        document.body.classList.remove("is-resizing");
      };

      setPanelSplit(storedSplitPercent(), false);

      splitter.addEventListener("pointerdown", (event) => {
        if (window.matchMedia("(max-width: 980px)").matches) return;
        event.preventDefault();
        activePointerId = event.pointerId;
        try { splitter.setPointerCapture(event.pointerId); } catch (_) {}
        document.body.classList.add("is-resizing");
        applyClientX(event.clientX);
      });
      const trackPointer = (event) => {
        if (event.pointerId !== activePointerId) return;
        applyClientX(event.clientX);
      };
      const stopPointer = (event) => {
        if (event.pointerId === activePointerId) stopResize();
      };
      splitter.addEventListener("pointermove", trackPointer);
      document.addEventListener("pointermove", trackPointer);
      splitter.addEventListener("pointerup", stopPointer);
      splitter.addEventListener("pointercancel", stopPointer);
      document.addEventListener("pointerup", stopPointer);
      document.addEventListener("pointercancel", stopPointer);
      splitter.addEventListener("mousedown", (event) => {
        if (window.matchMedia("(max-width: 980px)").matches || event.button !== 0 || activePointerId !== null) return;
        event.preventDefault();
        activeMouse = true;
        document.body.classList.add("is-resizing");
        applyClientX(event.clientX);
      });
      document.addEventListener("mousemove", (event) => {
        if (!activeMouse) return;
        applyClientX(event.clientX);
      });
      document.addEventListener("mouseup", () => {
        if (activeMouse) stopResize();
      });
      splitter.addEventListener("dblclick", () => setPanelSplit(DEFAULT_SPLIT_PERCENT));
      splitter.addEventListener("keydown", (event) => {
        const current = Number(splitter.getAttribute("aria-valuenow")) || DEFAULT_SPLIT_PERCENT;
        const bounds = splitBounds();
        if (event.key === "ArrowLeft") { event.preventDefault(); setPanelSplit(current - 2); }
        else if (event.key === "ArrowRight") { event.preventDefault(); setPanelSplit(current + 2); }
        else if (event.key === "Home") { event.preventDefault(); setPanelSplit(bounds.min); }
        else if (event.key === "End") { event.preventDefault(); setPanelSplit(bounds.max); }
        else if (event.key === "Enter" || event.key === " ") { event.preventDefault(); setPanelSplit(DEFAULT_SPLIT_PERCENT); }
      });
      window.addEventListener("resize", () => {
        setPanelSplit(storedSplitPercent(), false);
      });
    }

    async function getJson(path) {
      const response = await fetch(path);
      if (!response.ok) {
        let message = await response.text();
        try {
          const payload = JSON.parse(message);
          if (payload && payload.error) message = payload.error;
        } catch (_) {}
        throw new Error(message);
      }
      return response.json();
    }

    async function loadMeta() {
      const health = await getJson("/api/health");
      $("meta").innerHTML = [
        health.db_name,
        `${health.gene_count.toLocaleString()} genes`,
        `${health.study_count.toLocaleString()} studies`,
        `${health.source_unit_count.toLocaleString()} source units`
      ].map((text) => `<span>${esc(text)}</span>`).join("");
    }

    function markSelected(symbol) {
      document.querySelectorAll("#genes tr").forEach((row) => {
        row.classList.toggle("selected", row.dataset.gene === symbol);
      });
    }

    const PAGE_SIZE = 100;
    const SORT_DEFAULTS = {
      rank: "asc",
      tier: "asc",
      gene: "asc",
      top: "asc",
      score: "desc",
      units: "desc",
      sign: "desc",
      lfc: "desc"
    };
    let sortState = { sort: "rank", order: "asc" };
    let page = { query: null, loaded: 0, total: 0 };
    let genesLoading = false;

    function geneRowHtml(gene) {
      return `
        <tr data-gene="${esc(gene.gene_symbol)}" tabindex="0" aria-label="Show evidence for ${esc(gene.gene_symbol)}">
          <td class="num">${esc(gene.degora_rank)}</td>
          <td>${tier(gene.evidence_tier)}</td>
          <td class="gene">${esc(gene.gene_symbol)}</td>
          <td class="num">${esc(gene.top_percent_label)}</td>
          <td class="num">${fmt(gene.degora_score, 2)}</td>
          <td class="num">${esc(gene.n_source_units)}</td>
          <td class="num">${fmt(gene.sign_concordance * 100, 1)}%</td>
          <td class="num">${fmt(gene.weighted_lfc, 2)}</td>
        </tr>`;
    }

    function currentQuery() {
      const params = new URLSearchParams();
      const q = $("query").value.trim();
      const minUnits = $("minUnits").value.trim();
      const direction = $("direction").value;
      if (q) params.set("q", q);
      if (minUnits) params.set("min_units", minUnits);
      if (direction) params.set("direction", direction);
      params.set("sort", sortState.sort);
      params.set("order", sortState.order);
      return params;
    }

    function updateSortHeaders() {
      document.querySelectorAll("[data-sort]").forEach((button) => {
        const active = button.dataset.sort === sortState.sort;
        const th = button.closest("th");
        const indicator = button.querySelector(".sort-indicator");
        if (th) th.setAttribute("aria-sort", active ? (sortState.order === "asc" ? "ascending" : "descending") : "none");
        if (indicator) indicator.textContent = active ? (sortState.order === "asc" ? "^" : "v") : "";
      });
    }

    function setGeneLoading(isLoading, title = "Loading genes...", note = "Refreshing the table") {
      genesLoading = isLoading;
      document.querySelector(".genes-panel").classList.toggle("is-loading", isLoading);
      document.querySelector(".genes-panel").setAttribute("aria-busy", isLoading ? "true" : "false");
      $("tableLoading").hidden = !isLoading;
      $("loadingTitle").textContent = title;
      $("loadingNote").textContent = note;
      $("load").disabled = isLoading;
      $("loadMore").disabled = isLoading;
      document.querySelectorAll("[data-sort]").forEach((button) => {
        button.disabled = isLoading;
      });
    }

    async function setSort(sort) {
      if (genesLoading) return;
      if (sortState.sort === sort) {
        sortState = { sort, order: sortState.order === "asc" ? "desc" : "asc" };
      } else {
        sortState = { sort, order: SORT_DEFAULTS[sort] || "asc" };
      }
      updateSortHeaders();
      $("geneTableScroll").scrollTop = 0;
      await loadGenes("Sorting genes...", "Applying the selected table order");
    }

    function updateGeneStatus() {
      $("status").textContent = page.loaded < page.total
        ? `Showing ${page.loaded.toLocaleString()} of ${page.total.toLocaleString()}`
        : `${page.total.toLocaleString()} gene${page.total === 1 ? "" : "s"}`;
      const more = $("loadMore");
      more.hidden = page.loaded >= page.total;
      if (!more.hidden) more.textContent = `Load ${Math.min(PAGE_SIZE, page.total - page.loaded).toLocaleString()} more`;
    }

    async function fetchGenePage(append, loadingTitle = "Loading genes...", loadingNote = "Refreshing the table") {
      if (genesLoading) return;
      if (append && page.loaded >= page.total) return;
      const title = append ? "Loading more genes..." : loadingTitle;
      const note = append ? "Appending the next page" : loadingNote;
      setGeneLoading(true, title, note);
      $("status").textContent = title;
      const params = new URLSearchParams(page.query);
      params.set("limit", String(PAGE_SIZE));
      params.set("offset", String(append ? page.loaded : 0));
      try {
        const data = await getJson(`/api/genes?${params.toString()}`);
        page.total = data.count;
        const html = data.genes.map(geneRowHtml).join("");
        if (append) {
          $("genes").insertAdjacentHTML("beforeend", html);
          page.loaded += data.genes.length;
        } else {
          $("genes").innerHTML = html;
          page.loaded = data.genes.length;
        }
        updateGeneStatus();
        if (!append) {
          if (page.loaded) loadGene(data.genes[0].gene_symbol);
          else $("detail").innerHTML = `<div class="empty">No genes match these filters. Try lowering <b>Min units</b> or clearing the gene search.</div>`;
        }
      } catch (error) {
        $("status").textContent = "error";
        if (!append) {
          $("genes").innerHTML = "";
          $("loadMore").hidden = true;
          $("detail").innerHTML = `<div class="empty">Could not load genes: ${esc(error.message)}</div>`;
        }
        return;
      } finally {
        setGeneLoading(false);
      }
    }

    async function loadGenes(loadingTitle = "Loading genes...", loadingNote = "Refreshing the table") {
      if (genesLoading) return;
      page = { query: currentQuery(), loaded: 0, total: 0 };
      await fetchGenePage(false, loadingTitle, loadingNote);
    }

    async function loadGene(symbol) {
      markSelected(symbol);
      let data;
      try {
        data = await getJson(`/api/genes/${encodeURIComponent(symbol)}`);
      } catch (error) {
        $("detail").innerHTML = `<div class="empty">Could not load ${esc(symbol)}: ${esc(error.message)}</div>`;
        return;
      }
      const gene = data.gene;
      const evidenceRows = data.evidence.map((row) => `
        <tr>
          <td>${esc(row.contributing_study_ids || row.study_id)}</td>
          <td>${esc(row.source_unit_id)}</td>
          <td>${esc(row.contributing_pipelines || row.pipeline || "")}</td>
          <td>${esc(row.contributing_assay_types || row.assay_type || "")}</td>
          <td>${esc(row.contributing_time_course_modes || row.time_course_mode || "")}</td>
          <td class="num">${fmt(row.source_reliability_weight, 2)}</td>
          <td class="num">${fmt(row.lfc, 2)}</td>
          <td class="num">${fmt(row.signed_z, 2)}</td>
          <td class="num">${fmt(row.normalized_rank, 4)}</td>
        </tr>
      `).join("");
      $("detail").innerHTML = `
        <div class="kv">
          <div class="metric" data-tip="DEGORA rank of this gene among all scored genes (1 = strongest)."><span>Rank</span><strong>${esc(gene.rank_label)}</strong></div>
          <div class="metric" data-tip="Where this gene sits as a percent of all scored genes (e.g. top 1%)."><span>Top fraction</span><strong>${esc(gene.top_percent_label)}</strong></div>
          <div class="metric" data-tip="Coarse confidence tier (A strongest to D weakest) from rank, support, and direction."><span>Evidence tier</span><strong>${tier(gene.evidence_tier)}</strong></div>
          <div class="metric" data-tip="The DEGORA quality-weighted prioritization score: a relative index, not a probability."><span>Score</span><strong>${fmt(gene.degora_score, 2)}</strong></div>
          <div class="metric" data-tip="How many independent source units support this gene (one paper = one unit)."><span>Source support</span><strong>${esc(gene.support_label)}</strong></div>
          <div class="metric" data-tip="Consensus regulation direction (up/down/flat) across the supporting sources."><span>Direction</span><strong>${esc(gene.direction_label)}</strong></div>
          <div class="metric" data-tip="Evidence-strength component combining repeated support and signal magnitude."><span>Evidence</span><strong>${fmt(gene.evidence_score, 2)}</strong></div>
          <div class="metric" data-tip="Contribution from how highly this gene ranked within each source's DEG list."><span>Rank signal</span><strong>${fmt(gene.rank_score_component, 2)}</strong></div>
          <div class="metric" data-tip="Sample-size-weighted mean log2 fold-change across supporting source units."><span>Weighted LFC</span><strong>${fmt(gene.weighted_lfc, 2)}</strong></div>
          <div class="metric" data-tip="Effect/rank/direction-focused prioritization score."><span>Priority</span><strong>${fmt(gene.priority_score, 2)}</strong></div>
          <div class="metric" data-tip="Combines support, source quality, direction confidence, and leave-one-source-out rank stability."><span>Reliability</span><strong>${fmt(gene.evidence_reliability_score, 2)}</strong></div>
          <div class="metric" data-tip="Direction-consistency index, shrunk toward 50% when evidence is weak or discordant."><span>Direction confidence</span><strong>${fmt(gene.direction_confidence_index * 100, 1)}%</strong></div>
          <div class="metric" data-tip="How stable the rank stays when each source unit is left out one at a time (higher = more robust)."><span>LOO stability</span><strong>${fmt(gene.loo_rank_stability_score * 100, 1)}%</strong></div>
        </div>
        <p class="sources">${badge(gene.consensus_direction)} DEGORA score is a relative prioritization score, not a probability.</p>
        <p class="sources">${esc(gene.source_units || "")}</p>
        <div class="ev-scroll">
        <table>
          <thead>
            <tr>
              <th data-tip="Study/contrast ID(s) contributing this row; one source unit may bundle several contrasts.">Study</th>
              <th data-tip="Independent source unit (paper/dataset); contrasts from one paper collapse into one unit.">Unit</th>
              <th data-tip="Analysis pipeline of the source DEG table (e.g. DESeq2, limma).">Pipeline</th>
              <th data-tip="Assay type of the source (e.g. RNA-seq, microarray).">Assay</th>
              <th data-tip="How same-source time-course contrasts were summarized (mean/early/late/peak_mean).">Time mode</th>
              <th class="num" data-tip="Source reliability weight: predeclared, sample-size-aware source-quality weight.">Rel</th>
              <th class="num" data-tip="log2 fold-change contributed by this source unit.">LFC</th>
              <th class="num" data-tip="Signed z-score (direction x significance) from the source p-value and sign.">z</th>
              <th class="num" data-tip="Within-study normalized rank of the gene (near 0 = top of that study's list).">Rank</th>
            </tr>
          </thead>
          <tbody>${evidenceRows}</tbody>
        </table>
        </div>
      `;
    }

    $("load").addEventListener("click", loadGenes);
    $("loadMore").addEventListener("click", () => fetchGenePage(true));
    document.querySelectorAll("[data-sort]").forEach((button) => {
      button.addEventListener("click", () => setSort(button.dataset.sort));
    });
    $("query").addEventListener("keydown", (event) => {
      if (event.key === "Enter") loadGenes();
    });
    $("genes").addEventListener("click", (event) => {
      const row = event.target.closest("tr");
      if (row && row.dataset.gene) loadGene(row.dataset.gene);
    });
    $("genes").addEventListener("keydown", (event) => {
      if (event.key !== "Enter" && event.key !== " ") return;
      const row = event.target.closest("tr");
      if (row && row.dataset.gene) { event.preventDefault(); loadGene(row.dataset.gene); }
    });

    const tip = $("tip");
    const placeTip = (text, x, y) => {
      tip.textContent = text;
      tip.style.display = "block";
      const r = tip.getBoundingClientRect();
      let left = x + 14;
      let top = y + 14;
      if (left + r.width > window.innerWidth - 6) left = x - r.width - 14;
      if (top + r.height > window.innerHeight - 6) top = y - r.height - 14;
      tip.style.left = Math.max(6, left) + "px";
      tip.style.top = Math.max(6, top) + "px";
    };
    const hideTip = () => { tip.style.display = "none"; };
    document.addEventListener("mousemove", (event) => {
      // Curated explanations (metric cards, column headers) take priority.
      const tipped = event.target.closest("[data-tip]");
      if (tipped) { placeTip(tipped.getAttribute("data-tip"), event.clientX, event.clientY); return; }
      // Otherwise, if a table cell is truncated, reveal its full content on hover.
      const cell = event.target.closest("td, th");
      if (cell && cell.scrollWidth > cell.clientWidth + 1) {
        const text = cell.textContent.trim();
        if (text) { placeTip(text, event.clientX, event.clientY); return; }
      }
      hideTip();
    });
    document.addEventListener("mouseleave", hideTip);
    window.addEventListener("scroll", hideTip, true);

    updateSortHeaders();
    initPanelResize();
    loadMeta().then(loadGenes).catch((error) => {
      $("status").textContent = "error";
      $("detail").innerHTML = `<div class="empty">${esc(error.message)}</div>`;
    });
  </script>
</body>
</html>
"""


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def _connect(db_path: Path) -> sqlite3.Connection:
    if not db_path.exists():
        # Expose only the filename over the wire: the absolute path can leak local
        # directory names to a client if the server is bound to a non-loopback host.
        raise FileNotFoundError(f"DEGORA database is not available: {Path(db_path).name}")
    uri = f"file:{quote(str(db_path.resolve()), safe='/')}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def _looks_like_local_path(value: Any) -> bool:
    text = str(value).strip()
    if not text:
        return False
    if Path(text).is_absolute():
        return True
    if len(text) >= 3 and text[1] == ":" and text[2] in {"/", "\\"}:
        return True
    return text.startswith(("/", "\\\\"))


LOCAL_PATH_REDACTION = "[redacted: local path]"


def _contains_local_path(value: Any) -> bool:
    text = str(value).strip()
    if not text:
        return False
    if _looks_like_local_path(text):
        return True
    return any(_looks_like_local_path(part) for part in text.replace("\n", ";").split(";"))


def _redact_meta_for_network(meta: dict[str, str]) -> dict[str, str]:
    redacted = dict(meta)
    for key, value in meta.items():
        if key.endswith(("_path", "_paths", "_dir", "_dirs")) and _contains_local_path(value):
            redacted[key] = LOCAL_PATH_REDACTION
    return redacted


def _redact_record_paths_for_network(record: dict[str, Any]) -> dict[str, Any]:
    redacted = dict(record)
    for key, value in record.items():
        if key.endswith(("_path", "_paths", "_dir", "_dirs")) and _contains_local_path(value):
            redacted[key] = LOCAL_PATH_REDACTION
    return redacted


def _redact_records_paths_for_network(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [_redact_record_paths_for_network(record) for record in records]


def _row_dicts(cursor: sqlite3.Cursor) -> list[dict[str, Any]]:
    return [dict(row) for row in cursor.fetchall()]


def _one_row(cursor: sqlite3.Cursor) -> dict[str, Any] | None:
    row = cursor.fetchone()
    return dict(row) if row is not None else None


def _int_param(params: dict[str, list[str]], name: str, default: int, *, minimum: int = 0, maximum: int = 500) -> int:
    raw = params.get(name, [str(default)])[0]
    try:
        value = int(raw)
    except ValueError:
        raise ValueError(f"{name} must be an integer") from None
    if value < minimum or value > maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return value


def _float_param(params: dict[str, list[str]], name: str, default: float, *, minimum: float = 0.0) -> float:
    raw = params.get(name, [str(default)])[0]
    try:
        value = float(raw)
    except ValueError:
        raise ValueError(f"{name} must be a number") from None
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum:g}")
    return value


def _text_param(params: dict[str, list[str]], name: str, default: str = "", *, maximum: int = 128) -> str:
    value = params.get(name, [default])[0].strip()
    if len(value) > maximum:
        raise ValueError(f"{name} is too long; maximum length is {maximum} characters")
    return value


def _escape_like_literal(value: str) -> str:
    return value.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


GENE_SORT_COLUMNS = {
    "rank": ("degora_rank", "ASC"),
    "tier": ("evidence_tier", "ASC"),
    "gene": ("gene_symbol", "ASC"),
    "top": ("top_percent", "ASC"),
    "score": ("degora_score", "DESC"),
    "units": ("n_source_units", "DESC"),
    "sign": ("sign_concordance", "DESC"),
    "lfc": ("weighted_lfc", "DESC"),
}


def _gene_order_clause(params: dict[str, list[str]]) -> tuple[str, str, str]:
    sort_key = _text_param(params, "sort", "rank", maximum=32).lower()
    if sort_key not in GENE_SORT_COLUMNS:
        allowed = ", ".join(sorted(GENE_SORT_COLUMNS))
        raise ValueError(f"sort must be one of: {allowed}")
    column, default_order = GENE_SORT_COLUMNS[sort_key]
    raw_order = _text_param(params, "order", default_order.lower(), maximum=8).lower()
    if raw_order in {"asc", "ascending"}:
        order = "ASC"
    elif raw_order in {"desc", "descending"}:
        order = "DESC"
    else:
        raise ValueError("order must be asc or desc")
    return f"{column} IS NULL ASC, {column} {order}, degora_rank ASC, gene_symbol ASC", sort_key, order.lower()


class DegoraRequestHandler(BaseHTTPRequestHandler):
    """Serve the static browser UI and JSON endpoints."""

    server: "DegoraHttpServer"

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        if self.server.quiet:
            return
        super().log_message(format, *args)

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        try:
            if parsed.path in {"/", "/index.html"}:
                self._send_html(INDEX_HTML)
            elif parsed.path == "/api/health":
                self._send_json(self._health())
            elif parsed.path == "/api/meta":
                self._send_json({"meta": self._meta()})
            elif parsed.path == "/api/studies":
                self._send_json({"studies": self._studies()})
            elif parsed.path == "/api/genes":
                self._send_json(self._genes(parse_qs(parsed.query)))
            elif parsed.path.startswith("/api/genes/"):
                symbol = unquote(parsed.path.removeprefix("/api/genes/")).upper()
                self._send_json(self._gene_detail(symbol))
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except FileNotFoundError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.NOT_FOUND)
        except ValueError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
        except sqlite3.Error as exc:
            self._send_json({"error": f"database error: {exc}"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def _send_html(self, html: str) -> None:
        encoded = html.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)

    def _send_json(self, payload: dict[str, Any], *, status: HTTPStatus = HTTPStatus.OK) -> None:
        encoded = json.dumps(_jsonable(payload), sort_keys=True).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)

    def _health(self) -> dict[str, Any]:
        with closing(_connect(self.server.db_path)) as connection:
            gene_count = connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
            study_count = connection.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
            source_unit_count = connection.execute("SELECT COUNT(DISTINCT source_unit_id) FROM studies").fetchone()[0]
            top_gene = connection.execute("SELECT gene_symbol FROM genes ORDER BY degora_rank LIMIT 1").fetchone()
        return {
            "status": "ok",
            "db_name": self.server.db_path.name,
            "gene_count": gene_count,
            "study_count": study_count,
            "source_unit_count": source_unit_count,
            "top_gene": top_gene[0] if top_gene else None,
        }

    def _meta(self) -> dict[str, str]:
        with closing(_connect(self.server.db_path)) as connection:
            rows = _row_dicts(connection.execute("SELECT key, value FROM meta ORDER BY key"))
        meta = {row["key"]: row["value"] for row in rows}
        if self.server.server_address[0] not in LOOPBACK_HOSTS:
            return _redact_meta_for_network(meta)
        return meta

    def _studies(self) -> list[dict[str, Any]]:
        with closing(_connect(self.server.db_path)) as connection:
            rows = _row_dicts(connection.execute("SELECT * FROM studies ORDER BY source_unit_id, study_id"))
        if self.server.server_address[0] not in LOOPBACK_HOSTS:
            return _redact_records_paths_for_network(rows)
        return rows

    def _genes(self, params: dict[str, list[str]]) -> dict[str, Any]:
        limit = _int_param(params, "limit", 50, minimum=1, maximum=500)
        offset = _int_param(params, "offset", 0, minimum=0, maximum=1_000_000)
        min_units = _int_param(params, "min_units", 1, minimum=1, maximum=10_000)
        min_score = _float_param(params, "min_score", 0.0, minimum=0.0)
        direction = _text_param(params, "direction", maximum=16).lower()
        query = _text_param(params, "q", maximum=128).upper()
        order_clause, sort_key, sort_order = _gene_order_clause(params)

        where = ["n_source_units >= ?", "degora_score >= ?"]
        values: list[Any] = [min_units, min_score]
        if direction:
            if direction not in {"up", "down", "flat"}:
                raise ValueError("direction must be up, down, or flat")
            where.append("consensus_direction = ?")
            values.append(direction)
        if query:
            where.append("gene_symbol LIKE ? ESCAPE '\\'")
            values.append(f"%{_escape_like_literal(query)}%")

        where_clause = " AND ".join(where)
        with closing(_connect(self.server.db_path)) as connection:
            count = connection.execute(f"SELECT COUNT(*) FROM genes WHERE {where_clause}", values).fetchone()[0]
            rows = _row_dicts(
                connection.execute(
                    f"SELECT * FROM genes WHERE {where_clause} ORDER BY {order_clause} LIMIT ? OFFSET ?",
                    values + [limit, offset],
                )
            )
        return {"count": count, "limit": limit, "offset": offset, "sort": sort_key, "order": sort_order, "genes": rows}

    def _gene_detail(self, symbol: str) -> dict[str, Any]:
        if not symbol:
            raise ValueError("gene symbol is required")
        if len(symbol) > 128:
            raise ValueError("gene symbol is too long; maximum length is 128 characters")
        with closing(_connect(self.server.db_path)) as connection:
            gene = _one_row(connection.execute("SELECT * FROM genes WHERE gene_symbol = ?", [symbol]))
            if gene is None:
                raise FileNotFoundError(f"gene not found: {symbol[:64]}")
            evidence = _row_dicts(
                connection.execute(
                    "SELECT * FROM gene_evidence WHERE gene_symbol = ? ORDER BY source_unit_id, study_id",
                    [symbol],
                )
            )
        if self.server.server_address[0] not in LOOPBACK_HOSTS:
            evidence = _redact_records_paths_for_network(evidence)
        return {"gene": gene, "evidence": evidence}


class DegoraHttpServer(ThreadingHTTPServer):
    """HTTP server carrying the database path for request handlers."""

    # Keep auto-port fallback deterministic on Windows too: SO_REUSEADDR there
    # can allow two listeners on the same port, so a busy-port probe would lie.
    allow_reuse_address = False

    def __init__(self, server_address: tuple[str, int], db_path: str | Path, *, quiet: bool = False) -> None:
        super().__init__(server_address, DegoraRequestHandler)
        self.db_path = Path(db_path).resolve()
        self.quiet = quiet


MAX_PORT_ATTEMPTS = 20


def create_server(
    db_path: str | Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    *,
    quiet: bool = False,
    auto_port: bool = True,
) -> DegoraHttpServer:
    """Bind the local server, auto-avoiding a port already held by another server.

    If the requested port is busy (e.g. a previous `degora serve` is still running),
    try the next few ports and finally an OS-assigned free port, so a leftover server
    no longer crashes a new run with `OSError: Address already in use`. Pass
    ``port=0`` or ``auto_port=False`` to bind exactly one port.
    """

    if port == 0 or not auto_port:
        return DegoraHttpServer((host, port), db_path, quiet=quiet)

    candidates = [port + offset for offset in range(MAX_PORT_ATTEMPTS)] + [0]
    last_error: OSError | None = None
    for index, candidate in enumerate(candidates):
        try:
            server = DegoraHttpServer((host, candidate), db_path, quiet=quiet)
        except OSError as exc:
            if exc.errno not in (errno.EADDRINUSE, errno.EACCES):
                raise
            last_error = exc
            continue
        if index > 0 and not quiet:
            print(
                f"Port {port} is already in use (another local DEGORA server may be running); "
                f"using port {server.server_address[1]} instead.",
                flush=True,
            )
        return server
    raise OSError(
        errno.EADDRINUSE,
        f"Could not bind {host} on port {port} or the next {MAX_PORT_ATTEMPTS - 1} ports. "
        "Stop the other server, or pass a free --port.",
    ) from last_error


LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "::1", "::ffff:127.0.0.1"})


def serve(db_path: str | Path, host: str = "127.0.0.1", port: int = 8765, *, quiet: bool = False) -> None:
    db_path = Path(db_path)
    if not db_path.exists():
        # Fail before binding so a missing DB never serves (and never leaks its path over HTTP).
        raise FileNotFoundError(f"DEGORA database does not exist: {db_path}")
    if host not in LOOPBACK_HOSTS:
        print(
            f"WARNING: serving on {host} exposes this no-authentication, read-only DEGORA "
            "database browser to your network. Anyone who can reach this port can read every "
            "gene, evidence row, study, and source URL. Use the default 127.0.0.1 unless you "
            "intend public access.",
            file=sys.stderr,
        )
    server = create_server(db_path, host=host, port=port, quiet=quiet)
    address, bound_port = server.server_address
    print(f"DEGORA browser/API: http://{address}:{bound_port}", flush=True)
    print(f"Database: {server.db_path}", flush=True)
    try:
        server.serve_forever()
    finally:
        server.server_close()
