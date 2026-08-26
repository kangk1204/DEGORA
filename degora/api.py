"""Local HTTP API and browser UI for a DEGORA SQLite score database."""

from __future__ import annotations

import errno
import inspect
import ipaddress
import json
import math
import os
import re
import secrets
import shutil
import socket
import sqlite3
import sys
import threading
import time
from contextlib import closing
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable
from urllib.parse import parse_qs, parse_qsl, quote, unquote, urlencode, urlparse, urlsplit, urlunsplit

from . import format_version_info, runtime_version_info
from .excel_export import EXCEL_ERROR_LITERALS
from .harmonize import canonical_gene_symbol
from .score_db import (
    PRIMARY_CONCORDANCE_COLUMN,
    PRIMARY_DIRECTION_COLUMN,
    PRIMARY_RANK_COLUMN,
    PRIMARY_SCORE_COLUMN,
    PRIMARY_TOP_PERCENT_COLUMN,
)
TOKEN_REDACTION = "[redacted]"
_TOKEN_ARG_RE = re.compile(r"(?i)(--api-token(?:=|\s+))('[^']*'|\"[^\"]*\"|[^\s]+)")
_TOKEN_QUERY_KEY_RE = r"(?:t|%74)(?:o|%6f)(?:k|%6b)(?:e|%65)(?:n|%6e)"
_TOKEN_QUERY_RE = re.compile(rf"(?i)(^|[?&#;\s])({_TOKEN_QUERY_KEY_RE})=[^&#;\s]*")
_TOKEN_HEADER_RE = re.compile(r"(?i)(X-DEGORA-Token\s*[:=]\s*)('[^']*'|\"[^\"]*\"|[^\s]+)")


def redact_token_from_url_text(text: str) -> str:
    """Redact token values from request lines, URLs, and command snippets."""

    def _replace_query_token(match: re.Match[str]) -> str:
        return f"{match.group(1)}{match.group(2)}={TOKEN_REDACTION}"

    text = _TOKEN_QUERY_RE.sub(_replace_query_token, text)
    text = _TOKEN_ARG_RE.sub(rf"\1{TOKEN_REDACTION}", text)
    return _TOKEN_HEADER_RE.sub(rf"\1{TOKEN_REDACTION}", text)


def strip_token_query_param(url: str) -> str:
    """Remove token query parameters while preserving all non-token parameters."""

    parts = urlsplit(url)
    query = urlencode(
        [(key, value) for key, value in parse_qsl(parts.query, keep_blank_values=True) if key.lower() != "token"],
        doseq=True,
    )
    return urlunsplit((parts.scheme, parts.netloc, parts.path, query, parts.fragment))


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="referrer" content="no-referrer">
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
      min-width: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 14px;
      height: 100vh;
      overflow: hidden;
      display: grid;
      grid-template-columns: minmax(0, 1fr);
      grid-template-rows: auto minmax(0, 1fr) auto;
    }
    header {
      display: flex;
      min-width: 0;
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
    [hidden] { display: none !important; }
    .brand { display: flex; align-items: center; gap: 12px; min-width: 0; }
    .brand-mark {
      display: inline-grid;
      place-items: center;
      width: 34px;
      height: 34px;
      border-radius: 10px;
      color: #fff;
      background: linear-gradient(145deg, #0f766e, #115e59);
      box-shadow: 0 8px 20px rgba(15, 118, 110, .2);
      font-weight: 850;
    }
    .brand-copy { min-width: 0; }
    .brand-copy small { display: none; }
    .top-nav { display: flex; gap: 6px; margin-left: auto; }
    .top-nav button {
      width: auto;
      height: 32px;
      padding: 0 13px;
      border-color: transparent;
      color: var(--muted);
      background: transparent;
    }
    .top-nav button.active { color: var(--accent); background: #e7f5f2; border-color: #cce7e1; }
    .discovery-view {
      grid-row: 2 / span 2;
      min-height: 0;
      min-width: 0;
      overflow: auto;
      padding: 18px 24px 42px;
      background: var(--bg);
    }
    .discovery-shell { width: min(1380px, 100%); min-width: 0; margin: 0 auto; }
    .search-hero { max-width: 940px; margin: 0 auto 22px; padding: 76px 0 28px; text-align: center; }
    .discovery-view.has-results { padding-top: 12px; }
    .discovery-view.has-results .search-hero { margin-bottom: 12px; padding: 0; }
    .discovery-view.has-results .search-hero > .eyebrow,
    .discovery-view.has-results .search-hero > h2,
    .discovery-view.has-results .search-hero > p,
    .discovery-view.has-results .policy-row,
    .discovery-view.has-results .cross-species-action { display: none; }
    .eyebrow {
      display: none;
      color: var(--accent);
      font-size: 11px;
      font-weight: 850;
      letter-spacing: .16em;
      text-transform: uppercase;
    }
    .search-hero h2 { margin: 0 0 8px; font-size: clamp(30px, 4vw, 44px); letter-spacing: -.04em; }
    .search-hero p { max-width: 700px; margin: 0 auto 18px; color: var(--muted); font-size: 14px; line-height: 1.45; }
    .species-tabs {
      display: flex;
      align-self: center;
      gap: 2px;
      padding: 4px;
      margin: 0;
      border: 0;
      border-radius: 10px;
      background: #eef0ed;
    }
    .species-tabs button {
      width: auto;
      min-width: 78px;
      height: 38px;
      border: 0;
      border-radius: 8px;
      color: var(--muted);
      background: transparent;
    }
    .species-tabs button.active {
      color: #fff;
      background: var(--accent);
      box-shadow: 0 5px 14px rgba(15, 118, 110, .2);
    }
    .discovery-search {
      display: grid;
      grid-template-columns: auto minmax(0, 1fr) 116px;
      gap: 8px;
      max-width: 900px;
      margin: 0 auto;
      padding: 7px;
      border: 1px solid #cbd8d3;
      border-radius: 16px;
      background: #fff;
      box-shadow: 0 16px 42px rgba(31, 41, 51, .09);
    }
    .discovery-search input { height: 46px; border: 0; font-size: 16px; padding: 0 14px; outline: none; }
    .discovery-search input:focus-visible {
      outline: 2px solid var(--accent);
      outline-offset: -2px;
      border-radius: 10px;
    }
    .discovery-search button { height: 46px; border-radius: 11px; }
    .cross-species-action { display: flex; justify-content: center; margin-top: 9px; }
    .cross-species-action button { width: auto; height: 30px; padding: 0 10px; border-color: transparent; background: transparent; color: var(--muted); font-size: 11px; }
    .policy-row { display: flex; justify-content: center; flex-wrap: wrap; gap: 7px; margin-top: 9px; }
    .policy-chip {
      display: inline-flex;
      align-items: center;
      gap: 6px;
      padding: 5px 9px;
      border: 0;
      border-radius: 999px;
      background: transparent;
      color: var(--muted);
      font-size: 11px;
    }
    .policy-chip::before { content: none; }
    .discovery-card {
      border: 1px solid var(--line);
      border-radius: 14px;
      background: rgba(255,255,255,.96);
      box-shadow: 0 10px 30px rgba(31, 41, 51, .045);
      overflow: hidden;
    }
    #discoveryResultsCard { overflow: visible; }
    #discoveryResultsCard .discovery-card-head {
      border-radius: 14px 14px 0 0;
      background: rgba(255,255,255,.98);
    }
    .discovery-card + .discovery-card { margin-top: 18px; }
    .discovery-card-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 16px 18px;
      border-bottom: 1px solid var(--line);
    }
    .discovery-card-head > div { min-width: 0; }
    .discovery-card-head h3 { margin: 0; font-size: 15px; }
    .discovery-card-head p { margin: 3px 0 0; color: var(--muted); font-size: 12px; }
    .study-action-bar {
      position: sticky;
      top: 0;
      z-index: 5;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      min-height: 52px;
      padding: 8px 16px;
      border-bottom: 1px solid var(--line);
      background: rgba(247, 250, 249, .97);
      box-shadow: 0 8px 18px rgba(31, 41, 51, .06);
      backdrop-filter: blur(8px);
    }
    .study-order-control { display: flex; align-items: center; flex-wrap: wrap; gap: 5px 8px; }
    .study-order-status { color: var(--muted); font-size: 12px; }
    .readiness-note { flex-basis: 100%; color: var(--muted); font-size: 11px; line-height: 1.35; }
    .study-order-reset { width: auto; height: 28px; padding: 0 8px; border-color: var(--line); background: #fff; color: var(--accent); font-size: 11px; }
    .study-action-bar .selection-action { margin-left: auto; }
    .study-action-bar .selection-action button { min-width: 0; }
    .study-action-bar .selection-action .action-secondary { min-width: 72px; }
    .results-scroll { overflow: auto; }
    .result-filter { display: flex; align-items: center; gap: 8px; margin: 0 0 8px; flex-wrap: wrap; }
    .result-filter label { color: var(--muted); font-size: 12px; }
    .result-filter input { flex: 1; min-width: 200px; max-width: 420px; }
    .result-filter .filter-count { color: var(--muted); font-size: 12px; }
    .study-table { min-width: 1120px; table-layout: fixed; }
    .study-table th:first-child, .study-table td:first-child { width: 46px; text-align: center; }
    .study-table th:nth-child(2) { width: 29%; }
    .study-table th:nth-child(3) { width: 18%; }
    .study-table th:nth-child(4) { width: 12%; }
    .study-table th:nth-child(5) { width: 74px; }
    .study-table th:nth-child(6) { width: 14%; }
    .study-table th:nth-child(7) { width: 17%; }
    .study-table th:nth-child(8) { width: 92px; }
    .study-title { font-weight: 720; white-space: normal; line-height: 1.35; }
    .dataset-title { display: block; margin-top: 4px; color: var(--muted); font-size: 11px; white-space: normal; }
    .study-publication-meta { display: none; margin-top: 5px; color: var(--muted); font-size: 11px; line-height: 1.35; white-space: normal; }
    .study-inspect { width: auto; min-width: 70px; height: 32px; padding: 0 10px; }
    /* Every cell clips with an ellipsis; a cell holding only a button must not,
       or the column shows "Inspect …" with the dots meaning nothing. */
    td.inspect-cell { overflow: visible; text-overflow: clip; }
    .readiness-basis { display: block; margin-top: 2px; }
    .mobile-field-label { display: none; }
    .study-table input[type="checkbox"], .candidate-row input[type="checkbox"] { width: 17px; height: 17px; accent-color: var(--accent); }
    /* Browsers grey a disabled 17px checkbox almost imperceptibly, so a row the
       cap has locked looked identical to one waiting to be ticked. */
    .study-table input[type="checkbox"]:disabled { opacity: 0.3; cursor: not-allowed; }
    .selection-limit {
      margin: 0 0 8px;
      padding: 8px 11px;
      border: 1px solid #f5d68a;
      border-radius: 8px;
      background: #fdf6e3;
      color: #7a5a12;
      font-size: 12px;
      line-height: 1.5;
    }
    .table-footer {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 12px 16px;
      border-top: 1px solid var(--line);
      background: #fbfcfb;
    }
    .pager { display: flex; align-items: center; gap: 8px; }
    .pager button, .action-secondary { width: auto; background: #fff; color: var(--ink); border-color: var(--line); }
    .pager span { min-width: 90px; text-align: center; color: var(--muted); font-variant-numeric: tabular-nums; }
    .selection-action { display: flex; align-items: center; gap: 10px; }
    .selection-action button { width: auto; min-width: 170px; }
    .discovery-empty { padding: 42px 20px; text-align: center; color: var(--muted); line-height: 1.6; }
    .discovery-loading { padding: 34px 20px; text-align: center; color: var(--muted); }
    .candidate-study { padding: 16px 18px; border-bottom: 1px solid var(--line); }
    .candidate-study:last-child { border-bottom: 0; }
    .candidate-study h4 { margin: 0 0 4px; font-size: 14px; }
    .candidate-study > p { margin: 0 0 12px; color: var(--muted); font-size: 12px; }
    .candidate-row {
      display: grid;
      grid-template-columns: 28px minmax(180px, 1.4fr) minmax(220px, 1fr) minmax(220px, 1fr);
      gap: 10px;
      align-items: start;
      padding: 12px;
      margin-top: 8px;
      border: 1px solid var(--line);
      border-radius: 10px;
      background: #fbfcfb;
    }
    .candidate-name { font-weight: 720; overflow-wrap: anywhere; }
    .candidate-note { display: block; margin-top: 4px; color: var(--muted); font-size: 11px; line-height: 1.4; }
    .candidate-fields { display: grid; gap: 7px; }
    /* The author card carries six children but the track list above defines
       four columns, so these last two groups wrapped onto an implicit second
       row - the metadata fields landing inside the 28px checkbox column and
       colliding with the confirmation lines. Give them the full width under
       the header instead. */
    .candidate-row > .candidate-fields-wide,
    .candidate-row > .candidate-confirms { grid-column: 2 / -1; }
    .candidate-fields-wide { grid-template-columns: repeat(auto-fit, minmax(168px, 1fr)); column-gap: 10px; }
    .candidate-confirms { margin-top: 2px; padding-top: 9px; border-top: 1px solid var(--line); }
    .candidate-fields label { display: grid; gap: 4px; color: var(--muted); font-size: 11px; }
    .candidate-fields input, .candidate-fields select { height: 32px; font-size: 12px; }
    .candidate-fields [aria-invalid="true"], .sample-groups[aria-invalid="true"] { outline: 2px solid #fda4af; outline-offset: 1px; }
    .confirm-line { display: flex !important; grid-template-columns: none !important; align-items: center; gap: 7px !important; color: var(--ink) !important; }
    .confirm-line input { width: 16px; height: 16px; }
    /* A confirmation that does not gate this row is not a question for this row. */
    .confirm-line.not-required { display: none !important; }
    .candidate-tools { display: flex; flex-wrap: wrap; align-items: center; gap: 7px; }
    .candidate-tools button { width: auto; min-width: 0; height: 30px; padding: 0 10px; border-radius: 6px; font-size: 11px; }
    .sample-groups { display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 5px; max-height: 300px; overflow: auto; }
    .sample-counts { grid-column: 1 / -1; display: flex; gap: 8px; align-items: center; color: var(--muted); font-size: 11px; }
    .sample-counts strong { color: var(--ink); font-variant-numeric: tabular-nums; }
    .sample-item {
      display: grid;
      grid-template-columns: minmax(0, 1fr) 92px;
      grid-template-areas: "id select" "label select" "traits select";
      align-items: center;
      gap: 1px 6px;
      padding: 4px 5px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
    }
    .sample-item span { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-size: 11px; }
    .sample-id { grid-area: id; font-variant-numeric: tabular-nums; color: var(--muted); font-size: 10px; }
    .sample-label { grid-area: label; color: var(--ink); font-weight: 620; }
    .sample-traits { grid-area: traits; color: var(--muted); font-size: 10px; }
    .sample-label-missing { color: var(--muted); font-weight: 500; font-style: italic; }
    .prepared-blocked {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 8px;
      margin-bottom: 10px;
      padding: 10px 12px;
      border: 1px solid #f0b8b8;
      border-radius: 9px;
      background: #fff3f3;
      color: #8b3333;
      font-size: 12px;
      line-height: 1.55;
    }
    .prepared-blocked strong { color: #6f2626; }
    .shared-submission {
      display: block;
      margin-top: 3px;
      color: #8a6206;
      font-size: 10.5px;
      line-height: 1.4;
    }
    .independence-warning {
      grid-column: 1 / -1;
      margin: 0 0 10px;
      padding: 9px 11px;
      border: 1px solid #f5d68a;
      border-radius: 8px;
      background: #fdf6e3;
      color: #7a5a12;
      font-size: 12px;
      line-height: 1.5;
    }
    .prepared-blocked button { width: auto; min-width: 0; height: 28px; padding: 0 10px; border-radius: 7px; font-size: 11px; }
    /* Everything stays readable and selectable; only the controls go inert. */
    /* The two collapsible panels and the confirmation lines are children of the
       four-column candidate grid. Without a span they fill whatever cell comes
       next - the 28px checkbox column, on the second row - and every word of
       "Columns DEGORA read from the file" lands on its own line with the inputs
       clipped to four characters. */
    .candidate-row > .candidate-advanced,
    .candidate-row > .candidate-confirms,
    .candidate-row > .series-sample-note { grid-column: 1 / -1; }
    .candidate-advanced > summary { cursor: pointer; list-style: none; font-weight: 700; color: var(--muted); font-size: 12px; }
    .candidate-advanced > summary::-webkit-details-marker { display: none; }
    .candidate-advanced > summary::before { content: "\25B8 "; }
    .candidate-advanced[open] > summary::before { content: "\25BE "; }
    .unanalyzable-group { margin-top: 18px; border-top: 1px dashed var(--line); padding-top: 10px; }
    .preferred-note { margin: 4px 0 2px; color: var(--muted); font-size: 12px; }
    .alternative-candidates { margin-top: 8px; }
    .alternative-candidates > summary { cursor: pointer; color: var(--muted); font-size: 12px; font-weight: 700; }
    .unanalyzable-group > summary { cursor: pointer; font-weight: 700; color: var(--muted); }
    .is-unanalyzable .candidate-row { background: #f7f8f8; }
    .is-unanalyzable .candidate-row input:disabled,
    .is-unanalyzable .candidate-row select:disabled,
    .is-unanalyzable .candidate-row button:disabled { opacity: 0.5; cursor: not-allowed; }
    .unusable-files ul { margin: 5px 0 6px; padding-left: 18px; }
    .unusable-files li { margin: 2px 0; }
    .unusable-files code { font-size: 11px; overflow-wrap: anywhere; }
    .sample-labels-missing { grid-column: 1 / -1; color: #7a5a12; }
    .sample-bulk { grid-column: 1 / -1; display: flex; flex-wrap: wrap; align-items: center; gap: 5px; padding-bottom: 2px; }
    .sample-bulk .sample-filter { flex: 1 1 130px; min-width: 110px; height: 26px; font-size: 11px; }
    .sample-bulk button { width: auto; min-width: 0; height: 26px; padding: 0 8px; border-radius: 6px; font-size: 11px; }
    .sample-bulk-count { color: var(--muted); font-size: 10px; }
    .sample-item.is-filtered-out { display: none; }
    .sample-item select { grid-area: select; }
    .sample-item select { height: 28px; padding: 0 5px; font-size: 11px; }
    .status-pill { display: inline-flex; padding: 4px 8px; border-radius: 999px; background: #e7f5f2; color: #0f5d56; font-size: 11px; font-weight: 720; }
    .status-pill.blocked { background: #f4f2e8; color: var(--warn); }
    .analysis-banner { padding: 18px; background: linear-gradient(120deg, #ecfdf8, #eff6ff); }
    .analysis-banner h3 { margin: 0 0 6px; }
    .analysis-banner p { margin: 0 0 12px; color: var(--muted); }
    .analysis-banner button { width: auto; }
    .analysis-actions { display: flex; flex-wrap: wrap; gap: 10px; }
    .error-box, .notice-box { margin: 12px 18px 18px; padding: 11px 12px; border-radius: 8px; }
    .error-box { border: 1px solid #fecdd3; background: #fff1f2; color: var(--danger); }
    .notice-box { border: 1px solid #bfe4dc; background: #effaf7; color: #0f5d56; }
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
    .section-title-group { display: flex; align-items: center; gap: 9px; min-width: 0; }
    .export-action {
      width: auto;
      height: 28px;
      padding: 0 9px;
      border-color: var(--line);
      background: #fff;
      color: var(--accent);
      font-size: 11px;
      white-space: nowrap;
    }
    .export-action:hover, .export-action:focus-visible { border-color: #8bc8bc; background: #eef8f6; }
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
    button:disabled { cursor: not-allowed; opacity: .55; }
    a { color: #0f6f68; text-decoration: none; }
    a:hover, a:focus-visible { text-decoration: underline; text-underline-offset: 2px; }
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
    .sort-head:hover { background: #eef7f5; }
    /* The hover tint alone made keyboard focus indistinguishable from hover. */
    .sort-head:focus-visible {
      background: #eef7f5;
      outline: 2px solid var(--accent);
      outline-offset: -2px;
    }
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
    /* Determinate variant: the span carries an inline width percentage. */
    .loading-bar.is-determinate::before { content: none; }
    .loading-bar.is-determinate > span {
      display: block;
      height: 100%;
      border-radius: inherit;
      background: var(--accent);
      transition: width 320ms ease-out;
    }
    /* While a search, prepare or analysis is running, every control that would
       start competing work is dimmed and locked. The progress panel itself is
       deliberately excluded so it stays legible. */
    .discovery-view.is-busy .discovery-search,
    .discovery-view.is-busy .cross-species-action,
    .discovery-view.is-busy .study-action-bar,
    .discovery-view.is-busy .table-footer,
    .discovery-view.is-busy .study-table {
      opacity: 0.45;
      filter: saturate(0.65);
      pointer-events: none;
      transition: opacity 180ms ease-out, filter 180ms ease-out;
    }
    .discovery-view.is-busy .discovery-search { cursor: progress; }
    .discovery-view .discovery-search,
    .discovery-view .cross-species-action,
    .discovery-view .study-action-bar,
    .discovery-view .table-footer,
    .discovery-view .study-table {
      transition: opacity 180ms ease-out, filter 180ms ease-out;
    }
    @media (prefers-reduced-motion: reduce) {
      .discovery-view .discovery-search,
      .discovery-view .cross-species-action,
      .discovery-view .study-action-bar,
      .discovery-view .table-footer,
      .discovery-view .study-table { transition: none; }
    }
    .search-progress {
      width: min(420px, 100%);
      margin: 0 auto;
      text-align: left;
    }
    .search-progress .loading-note + .loading-note { margin-top: 4px; }
    .job-cancel {
      width: auto;
      min-width: 0;
      height: 30px;
      margin-top: 14px;
      padding: 0 12px;
      border: 1px solid var(--line);
      border-radius: 7px;
      background: transparent;
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
      cursor: pointer;
    }
    .job-cancel:hover:not(:disabled) { color: var(--ink); border-color: var(--ink); }
    .job-cancel:disabled { cursor: default; opacity: 0.6; }
    @media (prefers-reduced-motion: reduce) {
      .loading-bar::before { animation: none; width: 100%; opacity: 0.55; }
      .loading-bar.is-determinate > span { transition: none; }
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
    .mobile-study-tools, .study-mobile-meta { display: none; }
    .evidence-kind {
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 3px 7px;
      background: #e7f5f2;
      color: #0f5d56;
      font-size: 11px;
      font-weight: 750;
      white-space: nowrap;
    }
    .evidence-kind.standardized { background: #eef2ff; color: #3730a3; }
    .evidence-kind.exploratory { background: #fff4d6; color: #8a4b08; }
    .deg-input {
      display: inline-flex;
      align-items: center;
      max-width: 144px;
      border: 1px solid transparent;
      border-radius: 999px;
      padding: 4px 8px;
      font-size: 10.5px;
      font-weight: 760;
      line-height: 1.2;
      white-space: normal;
    }
    .deg-input.author_deg_likely { border-color: #a7d8c8; background: #e7f5ef; color: #0d624f; }
    .deg-input.tabular_candidate { border-color: #b9cee8; background: #eef5fc; color: #245985; }
    .deg-input.matrix_fallback { border-color: #e7cd96; background: #fff6df; color: #855512; }
    .deg-input.not_detected { border-color: #d6dce2; background: #f4f6f8; color: #5c6773; }
    .deg-input.unresolved { border-color: #e3c4c4; background: #fff0f0; color: #8b3333; }
    .ev-scroll { max-height: 45vh; overflow: auto; border: 1px solid var(--line); border-radius: 6px; }
    .ev-scroll table {
      min-width: 1180px;
      table-layout: auto;
      border-collapse: separate;
      border-spacing: 0;
    }
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
    @media (max-width: 1120px) {
      header { align-items: center; flex-wrap: wrap; padding: 12px 18px 10px; }
      .meta { flex-basis: 100%; justify-content: flex-start; margin-left: 46px; gap: 12px; font-size: 12px; }
    }
    @media (max-width: 980px) {
      body { height: 100vh; min-height: 0; overflow: hidden; }
      main { grid-template-columns: 1fr; align-content: start; padding: 12px; overflow: auto; min-height: 0; height: auto; }
      .splitter { display: none; }
      section { min-height: auto; height: auto; }
      /* `.genes-panel, .evidence-panel { height: 100% }` out-specifies the
         `section` rule above, and `section { overflow: hidden }` then clipped
         the stacked evidence panel down to a ~28px sliver. Match the class
         specificity so the panel grows with its content instead. */
      .genes-panel, .evidence-panel { height: auto; overflow: visible; }
      .gene-table-shell { flex: none; height: 60vh; max-height: 60vh; }
      /* `flex: 1 1 0` collapsed the stacked detail body to its minimum height
         once the panel stopped filling the viewport. Size it to its content. */
      .evidence-panel .detail-body { overflow: visible; flex: none; height: auto; }
      .controls { grid-template-columns: 1fr 1fr; }
      .kv { grid-template-columns: 1fr 1fr; }
      .discovery-view { grid-row: 2; padding: 14px 12px 30px; overflow: auto; }
      .candidate-row { grid-template-columns: 28px 1fr; }
      .candidate-fields, .sample-groups { grid-column: 2; }
    }
    @media (max-width: 700px) {
      header { gap: 8px 10px; padding: 10px 12px 9px; }
      .brand { flex: 1; }
      .brand-copy small { font-size: 10px; }
      .top-nav { order: 3; display: grid; grid-template-columns: repeat(2, 1fr); width: 100%; margin: 0; }
      .meta {
        order: 4;
        flex-basis: auto;
        width: 100%;
        margin: 0;
        gap: 12px;
        flex-wrap: nowrap;
        overflow-x: auto;
        font-size: 11px;
        scrollbar-width: thin;
      }
      .meta span { white-space: nowrap; }
      .search-hero h2 { font-size: 26px; }
      .search-hero p { font-size: 13px; }
      .controls { grid-template-columns: 1fr; }
      .kv { grid-template-columns: 1fr; }
      .discovery-search { grid-template-columns: 1fr; }
      .discovery-search input { width: 100%; min-width: 0; }
      .species-tabs { display: grid; grid-template-columns: repeat(2, 1fr); }
      .species-tabs button { min-width: 0; }
      .policy-row { justify-content: flex-start; flex-wrap: nowrap; overflow-x: auto; padding: 1px 1px 4px; scrollbar-width: thin; }
      .policy-chip { flex: 0 0 auto; }
      .discovery-card-head { align-items: flex-start; flex-wrap: wrap; }
      .study-action-bar { align-items: flex-start; flex-direction: column; gap: 6px; padding: 8px 12px; }
      .study-action-bar .selection-action { width: 100%; margin-left: 0; display: grid; grid-template-columns: minmax(0, 1fr) auto; gap: 7px; }
      .study-action-bar .selection-action .status { grid-column: 1 / -1; }
      .study-action-bar .selection-action button { min-width: 0; }
      .table-footer { align-items: flex-start; flex-direction: column; }
      .table-footer .selection-action { width: 100%; flex-wrap: wrap; }
      .table-footer .selection-action button { flex: 1 1 150px; min-width: 0; }
      .pager { width: 100%; justify-content: space-between; }
      .pager span { min-width: 0; }
      .section-head { align-items: flex-start; }
      .section-title-group { align-items: flex-start; flex-direction: column; gap: 5px; }
      .sample-groups { grid-template-columns: 1fr; }
      .mobile-study-tools {
        display: flex;
        align-items: end;
        gap: 8px;
        padding: 10px 0;
      }
      .mobile-study-tools label { flex: 1; min-width: 0; color: var(--muted); font-size: 12px; }
      .mobile-study-tools select { width: 100%; margin-top: 4px; }
      .mobile-study-tools button { width: auto; min-width: 118px; }
      .results-scroll { overflow: visible; }
      .study-table { display: block; min-width: 0; }
      .study-table thead { display: none; }
      .study-table tbody { display: grid; gap: 10px; }
      .study-table tr {
        display: grid;
        grid-template-columns: 32px minmax(0, 1fr) auto;
        gap: 5px 9px;
        padding: 12px;
        border: 1px solid var(--line);
        border-radius: 8px;
        background: #fff;
      }
      .study-table td { display: block; padding: 0; border: 0; white-space: normal; }
      .study-table td:nth-child(1) { grid-column: 1; grid-row: 1 / span 3; padding-top: 3px; }
      .study-table td:nth-child(2) { grid-column: 2 / 4; grid-row: 1; }
      .study-table td:nth-child(3),
      .study-table td:nth-child(4),
      .study-table td:nth-child(5) { display: none; }
      .study-table td:nth-child(6) { grid-column: 2; grid-row: 2; }
      .study-table td:nth-child(7) { grid-column: 2; grid-row: 3; color: var(--muted); font-size: 12px; }
      .study-table td:nth-child(8) { grid-column: 3; grid-row: 2 / span 2; align-self: center; }
      .study-mobile-meta, .study-publication-meta { display: block; margin-top: 6px; color: var(--muted); font-size: 12px; line-height: 1.4; }
      .mobile-field-label { display: block; margin-bottom: 3px; color: var(--muted); font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: .04em; }
    }
    @media (max-width: 420px) {
      .policy-row { flex-wrap: wrap; overflow: visible; }
      .policy-chip { white-space: normal; }
    }
  </style>
</head>
<body>
  <header>
    <div class="brand">
      <span class="brand-mark" aria-hidden="true">D</span>
      <div class="brand-copy"><h1>DEGORA</h1><small>Source-resolved differential-expression evidence</small></div>
    </div>
    <nav class="top-nav" aria-label="Application views">
      <button id="discoverNav" class="active" type="button" aria-pressed="true">Discover</button>
      <button id="atlasNav" type="button" aria-pressed="false">Evidence atlas</button>
    </nav>
    <div class="meta" id="meta"></div>
  </header>
  <div class="discovery-view" id="discoveryView">
    <div class="discovery-shell">
      <div class="search-hero">
        <div class="eyebrow">Publication + linked public data discovery</div>
        <h2>Find Human or Mouse studies</h2>
        <p id="searchGuidance">Search papers and linked repositories, rank likely DEG inputs, inspect selected files, then run species-specific DEGORA.</p>
        <div class="discovery-search">
          <div class="species-tabs" role="tablist" aria-label="Species workspaces">
            <button id="humanSpeciesTab" role="tab" type="button" class="active" tabindex="0" aria-selected="true" aria-pressed="true" aria-controls="discoveryResultsCard" data-species="human">Human</button>
            <button id="mouseSpeciesTab" role="tab" type="button" tabindex="-1" aria-selected="false" aria-pressed="false" aria-controls="discoveryResultsCard" data-species="mouse">Mouse</button>
          </div>
          <input id="discoveryQuery" autocomplete="off" required minlength="2" maxlength="200" placeholder="e.g. hypoxia, Alzheimer’s, irradiation" aria-label="Search papers and data" aria-describedby="searchGuidance">
          <button id="discoverySearch" type="button">Search</button>
        </div>
        <div class="cross-species-action"><button id="discoverySearchBoth" class="action-secondary" type="button">Run separate Human + Mouse searches</button></div>
        <div class="policy-row" aria-label="Analysis policy">
          <span class="policy-chip">Human and Mouse never pooled</span>
          <span class="policy-chip">10 per page · up to 1,000 ranked</span>
        </div>
      </div>

      <section class="discovery-card" id="discoveryResultsCard" role="tabpanel" aria-labelledby="humanSpeciesTab">
        <div class="discovery-card-head">
          <div><h3 id="resultsTitle">Study results</h3><p id="resultsSubtitle" role="status" aria-live="polite" aria-atomic="true">Choose Human or Mouse, enter a keyword, and search.</p></div>
          <span class="status-pill" id="speciesState">Human workspace</span>
        </div>
        <div class="study-action-bar" id="discoveryActions" hidden>
          <div class="study-order-control"><span class="study-order-status" id="studyOrderStatus">Sort: DEG readiness</span><button class="study-order-reset" id="resetStudySort" type="button" hidden>Restore readiness ranking</button><span class="readiness-note">Readiness is an estimate until Prepare selection inspects file content.</span></div>
          <div class="selection-action">
            <span class="status" id="selectedStatus" role="status" aria-live="polite">0 / 20 selected · Human</span>
            <button class="action-secondary" id="downloadSearchExcel" type="button" disabled>Download search Excel</button>
            <button class="action-secondary" id="clearSelected" type="button">Clear</button>
            <button id="prepareSelected" type="button">Prepare selection</button>
          </div>
        </div>
        <div id="discoveryResults"><div class="discovery-empty">Search papers and linked data to see provisional publication matches.</div></div>
        <div class="table-footer" id="discoveryFooter" hidden>
          <div class="pager">
            <button id="discoveryPrev" type="button">Previous</button>
            <span id="discoveryPage" role="status" aria-live="polite" aria-atomic="true">Page 1 / 1</span>
            <button id="discoveryNext" type="button">Next</button>
          </div>
        </div>
        <div id="discoveryNotice" class="notice-box" role="status" tabindex="-1" hidden></div>
      </section>

      <section class="discovery-card" id="preparedCard" hidden>
        <div class="discovery-card-head">
          <div><h3>Prepared evidence</h3><p>Review contrast direction for author tables, or assign control and treatment samples for a labeled fallback.</p></div>
          <span class="status-pill" id="preparedStatus">Review required</span>
        </div>
        <div id="preparedCandidates"></div>
        <div class="table-footer">
          <label class="confirm-line" id="speciesConfirmLine" hidden><input id="speciesConfirmed" type="checkbox"> <span id="speciesConfirmText"></span></label>
          <span class="status" id="analysisEligibility" role="status" aria-live="polite">Select candidates from at least two independent studies.</span>
          <div class="selection-action">
            <button class="action-secondary" id="backToResults" type="button">Back to studies</button>
            <button id="runDiscoveryAnalysis" type="button">Run species-specific DEGORA</button>
          </div>
        </div>
        <div id="discoveryError" class="error-box" role="alert" hidden></div>
      </section>

      <section class="discovery-card" id="analysisCompleteCard" hidden>
        <div class="analysis-banner">
          <h3 id="analysisCompleteTitle">Analysis complete</h3>
          <p id="analysisCompleteText"></p>
          <div class="analysis-actions">
            <button id="openAnalysis" type="button">Open evidence atlas</button>
            <button class="action-secondary" id="downloadAnalysisExcel" type="button" disabled>Download Excel</button>
          </div>
        </div>
      </section>
    </div>
  </div>
  <main id="layoutMain" hidden>
    <section class="genes-panel">
      <div class="section-head">
        <div class="section-title-group"><h2>Genes</h2><button class="export-action" id="exportGenes" type="button" disabled>Download ranking CSV</button></div>
        <div class="status" id="status" role="status" aria-live="polite"></div>
      </div>
      <div class="controls">
        <input id="query" placeholder="Gene symbol" autocomplete="off" maxlength="128" aria-label="Gene symbol search" data-tip="Type part of a gene symbol to filter the list; leave blank to show all genes.">
        <input id="minUnits" type="number" min="1" max="10000" step="1" value="1" aria-label="Min source units" data-tip="Show only genes supported by at least this many independent source units (one paper = one unit).">
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
                <th class="sortable" data-tip="Primary quality-weighted DEGORA rank among all scored genes (1 = strongest)." aria-sort="ascending"><button class="sort-head" type="button" data-sort="rank">Rank <span class="sort-indicator" aria-hidden="true">^</span></button></th>
                <th class="sortable" data-tip="Confidence tier from rank, support, and direction (A strongest to D weakest)." aria-sort="none"><button class="sort-head" type="button" data-sort="tier">Tier <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="sortable" data-tip="Gene symbol." aria-sort="none"><button class="sort-head" type="button" data-sort="gene">Gene <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Primary quality-weighted rank position as a percent of all scored genes (e.g. top 1%)." aria-sort="none"><button class="sort-head" type="button" data-sort="top">Top <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="DEGORA quality-weighted prioritization score: a relative index, not a probability." aria-sort="none"><button class="sort-head" type="button" data-sort="score">Score <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Number of independent source units supporting the gene (one paper = one unit)." aria-sort="none"><button class="sort-head" type="button" data-sort="units">Units <span class="sort-indicator" aria-hidden="true"></span></button></th>
                <th class="num sortable" data-tip="Quality-weighted direction concordance: percent of supporting evidence agreeing on the consensus direction." aria-sort="none"><button class="sort-head" type="button" data-sort="sign">Sign <span class="sort-indicator" aria-hidden="true"></span></button></th>
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
      aria-valuenow="60"
      tabindex="0"
      data-tip="Drag left or right to resize the gene and evidence panels. Double-click to reset."
    ></div>
    <section class="evidence-panel">
      <div class="section-head">
        <div class="section-title-group"><h2>Evidence</h2><button class="export-action" id="exportEvidence" type="button" disabled>Download evidence CSV</button></div>
      </div>
      <div class="detail-body" id="detail">
        <div class="empty">Select a gene.</div>
      </div>
    </section>
  </main>
  <footer class="legend" id="atlasLegend" hidden>
    <span><b>Tier</b> <span class="tier">A</span><span class="tier B">B</span><span class="tier C">C</span><span class="tier D">D</span> A strongest to D weakest</span>
    <span><b>Direction</b> <span class="badge up">up</span><span class="badge down">down</span><span class="badge flat">flat</span></span>
    <span>Hover or focus any label or column header for its meaning.</span>
  </footer>
  <div id="tip" role="tooltip"></div>
  <script>
    const $ = (id) => document.getElementById(id);
    function readApiToken() {
      const url = new URL(window.location.href);
      const hashText = url.hash.startsWith("#") ? url.hash.slice(1) : url.hash;
      const hashParams = new URLSearchParams(hashText);
      let storedToken = "";
      try { storedToken = window.sessionStorage.getItem("degoraApiToken") || ""; } catch (_) {}
      const token = hashParams.get("token") || url.searchParams.get("token") || storedToken;
      if (hashParams.has("token") || url.searchParams.has("token")) {
        try { window.sessionStorage.setItem("degoraApiToken", token); } catch (_) {}
      }
      if (hashParams.has("token") || url.searchParams.has("token")) {
        url.searchParams.delete("token");
        hashParams.delete("token");
        const cleanHash = hashParams.toString();
        const cleanUrl = `${url.pathname}${url.search}${cleanHash ? `#${cleanHash}` : ""}`;
        window.history.replaceState(null, "", cleanUrl);
      }
      return token;
    }
    const API_TOKEN = readApiToken();
    const fmt = (value, digits = 3) => {
      if (value === null || value === undefined || Number.isNaN(Number(value))) return "";
      const n = Number(value);
      if (Math.abs(n) >= 100) return n.toFixed(1);
      if (Math.abs(n) < 0.001 && n !== 0) return n.toExponential(2);
      return n.toFixed(digits);
    };
    const fmtNullablePercent = (value, digits = 1) => {
      if (value === null || value === undefined || !Number.isFinite(Number(value))) return "N/A";
      return `${fmt(Number(value) * 100, digits)}%`;
    };
    const reliabilityBasisLabel = (gene) => {
      const used = Number(gene?.evidence_reliability_components_used);
      return Number.isInteger(used) && used >= 1 && used <= 4 ? `${used}/4 diagnostics` : "available diagnostics";
    };
    const looStabilityLabel = (gene) => {
      const value = fmtNullablePercent(gene?.loo_rank_stability_score, 1);
      const evaluable = Number(gene?.loo_rank_evaluable_folds);
      const total = Number(gene?.loo_total_folds);
      const foldText = Number.isInteger(evaluable) && Number.isInteger(total)
        ? `${evaluable}/${total} evaluable folds`
        : "fold count unavailable";
      return `${value} · ${foldText}`;
    };
    const esc = (value) => String(value ?? "").replace(/[&<>"']/g, (ch) => ({
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#39;"
    }[ch]));
    const safeHttpUrl = (value) => {
      try {
        const text = String(value || "").trim();
        if (!text) return "";
        const url = new URL(text);
        return url.protocol === "http:" || url.protocol === "https:" ? url.href : "";
      } catch (_) {
        return "";
      }
    };
    const externalLink = (url, label, className = "") => {
      const href = safeHttpUrl(url);
      return href
        ? `<a${className ? ` class="${esc(className)}"` : ""} href="${esc(href)}" target="_blank" rel="noopener noreferrer">${esc(label)}</a>`
        : esc(label);
    };
    const formatElapsed = (seconds) => {
      const total = Math.max(0, Math.round(seconds));
      if (total < 60) return `${total}s`;
      const minutes = Math.floor(total / 60);
      return `${minutes}m ${String(total % 60).padStart(2, "0")}s`;
    };
    const csvCell = (value) => {
      if (value === null || value === undefined) return "";
      let text = String(value);
      if (typeof value === "string" && /^\\s*[=+\\-@]/.test(text)) text = `'${text}`;
      return `"${text.replace(/"/g, '""')}"`;
    };
    function downloadCsv(filename, headers, rows) {
      const lines = [headers.map(csvCell).join(","), ...rows.map((row) => row.map(csvCell).join(","))];
      const blob = new Blob(["\ufeff", lines.join("\\r\\n"), "\\r\\n"], { type: "text/csv;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = filename.replace(/[^A-Za-z0-9._-]+/g, "_");
      link.hidden = true;
      document.body.appendChild(link);
      link.click();
      link.remove();
      window.setTimeout(() => URL.revokeObjectURL(url), 1000);
    }
    const safeClass = (value) => String(value ?? "").replace(/[^A-Za-z0-9_-]/g, "");
    const badge = (direction) => `<span class="badge ${safeClass(direction)}">${esc(direction)}</span>`;
    const tier = (value) => `<span class="tier ${safeClass(value)}">${esc(value)}</span>`;
    function evidenceKind(row) {
      const types = String(row.contributing_source_input_types || row.source_input_type || "")
        .split(";").map((value) => value.trim()).filter(Boolean);
      if (types.some((value) => ["derived_count_table", "normalized_expression_matrix"].includes(value))) {
        return { label: types.includes("author_deg_table") ? "Mixed · exploratory" : "Exploratory fallback", className: "exploratory" };
      }
      if (types.some((value) => value === "limma_full_table")) {
        return { label: types.includes("author_deg_table") ? "Mixed standardized" : "Standardized reanalysis", className: "standardized" };
      }
      if (types.length && types.every((value) => value === "author_deg_table")) {
        return { label: "Author DEG", className: "author" };
      }
      return { label: types.join("; ") || "Input type unavailable", className: "standardized" };
    }
    const firstPresent = (...values) => values.find((value) => value !== null && value !== undefined && value !== "");
    const primaryRank = (gene) => firstPresent(gene.quality_weighted_degora_rank, gene.degora_rank);
    const primaryScore = (gene) => firstPresent(gene.quality_weighted_degora_score, gene.degora_score);
    const primaryTopPercent = (gene) => firstPresent(gene.quality_weighted_top_percent, gene.top_percent);
    const primaryConcordance = (gene) => firstPresent(gene.quality_weighted_sign_concordance, gene.sign_concordance);
    const primaryDirection = (gene) => firstPresent(gene.quality_weighted_consensus_direction, gene.consensus_direction, "flat");
    const primaryDirectionConfidence = (gene) => firstPresent(gene.quality_weighted_direction_confidence_index, gene.direction_confidence_index);
    function primaryRankLabel(gene) {
      const rank = primaryRank(gene);
      if (rank === null || rank === undefined || rank === "") return gene.rank_label || "";
      return `#${Number(rank).toLocaleString()}`;
    }
    function topPercentLabel(gene) {
      const value = primaryTopPercent(gene);
      if (value === null || value === undefined || value === "" || Number.isNaN(Number(value))) {
        return gene.top_percent_label || "";
      }
      const n = Number(value);
      if (n < 0.01) return `top ${n.toFixed(4)}%`;
      if (n < 1.0) return `top ${n.toFixed(3)}%`;
      return `top ${n.toFixed(2)}%`;
    }
    function topPercentTableLabel(gene) {
      const value = primaryTopPercent(gene);
      if (value === null || value === undefined || value === "" || Number.isNaN(Number(value))) return "—";
      const n = Number(value);
      if (n < 0.01) return `${n.toFixed(4)}%`;
      if (n < 1.0) return `${n.toFixed(3)}%`;
      return `${n.toFixed(2)}%`;
    }
    function primaryDirectionLabel(gene) {
      const concordance = primaryConcordance(gene);
      const direction = primaryDirection(gene);
      if (concordance === null || concordance === undefined || concordance === "" || Number.isNaN(Number(concordance))) {
        return gene.direction_label || direction;
      }
      return `${fmt(Number(concordance) * 100, 1)}% ${direction}-concordant`;
    }
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

      // The grid track percentage is of main's CONTENT box, but the split used to
      // be computed from its border box, so the handle landed padding-width away
      // from the cursor.
      const contentGeometry = () => {
        const rect = main.getBoundingClientRect();
        const styles = window.getComputedStyle(main);
        const left = rect.left
          + (parseFloat(styles.paddingLeft) || 0)
          + (parseFloat(styles.borderLeftWidth) || 0);
        const right = rect.right
          - (parseFloat(styles.paddingRight) || 0)
          - (parseFloat(styles.borderRightWidth) || 0);
        return { left, width: right - left };
      };
      const applyClientX = (clientX) => {
        const geometry = contentGeometry();
        if (geometry.width <= 0) return;
        setPanelSplit(((clientX - geometry.left) / geometry.width) * 100);
      };
      // Pressing the handle used to snap the split to the cursor immediately, so
      // the first click of a double-click moved the handle out from under the
      // second one and "double-click to reset" could never fire.
      const DRAG_THRESHOLD_PX = 3;
      let pressClientX = null;
      let dragging = false;
      const beginPress = (clientX) => {
        pressClientX = clientX;
        dragging = false;
      };
      const trackClientX = (clientX) => {
        if (!dragging) {
          if (pressClientX === null || Math.abs(clientX - pressClientX) < DRAG_THRESHOLD_PX) return;
          dragging = true;
        }
        applyClientX(clientX);
      };
      const stopResize = () => {
        activePointerId = null;
        activeMouse = false;
        pressClientX = null;
        dragging = false;
        document.body.classList.remove("is-resizing");
      };

      setPanelSplit(storedSplitPercent(), false);

      splitter.addEventListener("pointerdown", (event) => {
        if (window.matchMedia("(max-width: 980px)").matches) return;
        event.preventDefault();
        activePointerId = event.pointerId;
        try { splitter.setPointerCapture(event.pointerId); } catch (_) {}
        document.body.classList.add("is-resizing");
        beginPress(event.clientX);
      });
      const trackPointer = (event) => {
        if (event.pointerId !== activePointerId) return;
        trackClientX(event.clientX);
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
        beginPress(event.clientX);
      });
      document.addEventListener("mousemove", (event) => {
        if (!activeMouse) return;
        trackClientX(event.clientX);
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

    async function getJson(path, { action = false } = {}) {
      const headers = {};
      if (API_TOKEN) headers["X-DEGORA-Token"] = API_TOKEN;
      if (action) headers["X-DEGORA-Action"] = "1";
      const options = Object.keys(headers).length ? { headers } : {};
      const response = await fetch(path, options);
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

    async function postJson(path, payload) {
      const headers = { "Content-Type": "application/json", "X-DEGORA-Action": "1" };
      if (API_TOKEN) headers["X-DEGORA-Token"] = API_TOKEN;
      const response = await fetch(path, { method: "POST", headers, body: JSON.stringify(payload) });
      if (!response.ok) {
        let message = await response.text();
        try {
          const parsed = JSON.parse(message);
          if (parsed && parsed.error) message = parsed.error;
        } catch (_) {}
        throw new Error(message);
      }
      return response.json();
    }

    const DISCOVERY_PAGE_SIZE = 10;
    const DISCOVERY_GLOBAL_RANK_LIMIT = 1000;
    const MAX_SELECTED_STUDIES = 20;
    const newSpeciesState = () => ({
      query: "",
      suggestedQuery: "",
      page: 1,
      totalPages: 0,
      totalHits: 0,
      evaluatedStudies: 0,
      rankingLimit: DISCOVERY_GLOBAL_RANK_LIMIT,
      rankingTruncated: false,
      cacheHit: false,
      hasNext: false,
      studies: [],
      searchId: "",
      jobId: "",
      jobStatus: "",
      jobProgress: null,
      jobCancelled: false,
      cancelling: false,
      jobMessage: "",
      jobStartedAt: 0,
      prepareProgress: null,
      prepareCancelled: false,
      prepareJobId: "",
      prepareMessage: "",
      prepareJobStartedAt: 0,
      providerStatus: "",
      providerErrors: [],
      noticeLevel: "info",
      verified: false,
      selected: new Set(),
      sort: { key: "readiness", order: "desc" },
      textFilter: "",
      totalUnfiltered: 0,
      prepared: null,
      bundleId: "",
      run: null,
      draft: {},
      cloneCounter: 0,
      loading: false,
      error: "",
      notice: "",
      analysisError: "",
      preparing: false,
      analyzing: false,
      searchRequest: 0,
      prepareRequest: 0,
      analysisRequest: 0
    });
    const discoveryStates = { human: newSpeciesState(), mouse: newSpeciesState() };
    let activeSpecies = "human";
    let preferredDiscoverySpecies = "human";
    let discoveryOpened = false;
    let activeRunId = "";
    const speciesLabel = (key) => key === "mouse" ? "Mouse" : "Human";
    const activeDiscoveryState = () => discoveryStates[activeSpecies];
    const currentAtlasContext = () => {
      if (activeRunId) return { key: `${activeSpecies}:run:${activeRunId}`, kind: "run", species: activeSpecies, runId: activeRunId };
      if (activeSpecies === "human") return { key: "human:base", kind: "base", species: "human", runId: "" };
      return { key: `${activeSpecies}:none`, kind: "none", species: activeSpecies, runId: "" };
    };
    const atlasApi = (suffix, context = currentAtlasContext()) => context.kind === "run"
      ? `/api/discovery/runs/${context.runId}${suffix}`
      : `/api${suffix}`;

    function renderDiscoveryHeaderMeta(state = activeDiscoveryState()) {
      let activity = "ready to search";
      if (state.analyzing) activity = "analysis running";
      else if (state.run) activity = "analysis ready";
      else if (state.preparing) activity = "preparing studies";
      else if (state.prepared) activity = "studies prepared";
      else if (state.loading) activity = "searching NCBI";
      else if (state.query) activity = `${state.evaluatedStudies.toLocaleString()} assessed studies`;
      $("meta").innerHTML = [
        `${speciesLabel(activeSpecies)} workspace`,
        activity,
        "species kept separate"
      ].map((text) => `<span>${esc(text)}</span>`).join("");
    }

    function isDiscoveryView() {
      return !$("discoveryView").hidden;
    }

    function showView(name) {
      const discovery = name === "discover";
      $("discoveryView").hidden = !discovery;
      $("layoutMain").hidden = discovery;
      $("atlasLegend").hidden = discovery;
      $("discoverNav").classList.toggle("active", discovery);
      $("atlasNav").classList.toggle("active", !discovery);
      $("discoverNav").setAttribute("aria-pressed", discovery ? "true" : "false");
      $("atlasNav").setAttribute("aria-pressed", discovery ? "false" : "true");
      if (discovery) {
        renderDiscoveryHeaderMeta();
      } else {
        setPanelSplit(storedSplitPercent(), false);
        // The header kept a Discover statistic ("48 assessed studies") while the
        // atlas was showing a different corpus. ensureAtlasContext() short-circuits
        // when the context is unchanged, so refresh the header unconditionally.
        void loadMeta(currentAtlasContext(), atlasContextGeneration);
        void ensureAtlasContext();
      }
    }

    function degAssessment(study) {
      const value = study && study.deg_input_assessment;
      return value && typeof value === "object" ? value : {};
    }

    function degInputBadge(study) {
      const assessment = degAssessment(study);
      const allowed = new Set(["author_deg_likely", "tabular_candidate", "matrix_fallback", "not_detected", "unresolved"]);
      const tier = allowed.has(assessment.tier) ? assessment.tier : "unresolved";
      const label = assessment.label || "File availability unresolved";
      const files = Array.isArray(assessment.candidate_files) ? assessment.candidate_files : [];
      const detail = [assessment.basis, files.length ? `Candidates: ${files.join(", ")}` : ""].filter(Boolean).join(" ");
      return `<span class="deg-input ${tier}" title="${esc(detail)}">${esc(label)}</span>`;
    }

    function publicationKey(study) {
      return String(study.canonical_id || study.publication_id || study.source_unit_id || study.accession || study.pmid || study.doi || study.title || study.paper_title || "");
    }

    function asListText(value) {
      if (Array.isArray(value)) return value.filter(Boolean).join(", ");
      return value == null ? "" : String(value);
    }

    // The search badge is an estimate from metadata. Once preparation has
    // actually opened the files, the row should report what was found instead
    // of the guess that sent the reader there - otherwise the same study gets
    // picked again next week.
    // DEGORA needs two independent source units. When a preparation cannot
    // reach that, every field below is work that will be discarded the moment
    // a new selection is prepared - state.draft is reset on each run - so the
    // review is not merely blocked at the end, it is pointless from the start.
    function usableSourceUnits(prepared) {
      const units = new Set();
      (prepared && prepared.studies || []).forEach((study) => {
        if (!(study.files || []).some(eligibleCandidate)) return;
        units.add(String(study.source_unit_id || study.accession || study.canonical_id || ""));
      });
      units.delete("");
      return units.size;
    }

    function preparedOutcomes(state) {
      const prepared = state.prepared;
      if (!prepared) return {};
      const outcomes = {};
      const put = (keys, outcome) => {
        keys.filter(Boolean).forEach((key) => { outcomes[String(key)] = outcome; });
      };
      (prepared.studies || []).forEach((study) => {
        const usable = (study.files || []).filter(eligibleCandidate);
        const outcome = !usable.length
          ? { label: "no usable table", ok: false, reason: "Preparation opened every linked file and none was a DEG table or expression matrix." }
          : usable.some(isAuthorReviewable)
          ? { label: "author DEG ready", ok: true, reason: "An author DEG table was resolved; confirm the contrast direction." }
          : { label: "needs group assignment", ok: true, reason: "Only an expression matrix was resolved; assign control and treatment samples." };
        put([publicationKey(study), study.accession, study.canonical_id, study.source_unit_id], outcome);
      });
      (prepared.excluded_studies || []).forEach((item) => {
        put(
          [item.canonical_id, item.accession, item.source_unit_id, item.paper_title],
          { label: "excluded", ok: false, reason: item.reason || "" },
        );
      });
      return outcomes;
    }

    function readinessBadge(study, outcomes) {
      const outcome = outcomes
        ? outcomes[publicationKey(study)] || outcomes[String(study.accession || "")] || null
        : null;
      if (outcome) {
        return `<span class="deg-input ${outcome.ok ? "author_deg_likely" : "unresolved"}"`
          + ` title="${esc(outcome.reason)}">${esc(`prepared · ${outcome.label}`)}</span>`;
      }
      const detail = study.data_readiness && typeof study.data_readiness === "object" ? study.data_readiness : {};
      const readiness = detail.verification_state || detail.tier || study.readiness || study.readiness_label || study.status || "provisional";
      const verified = Boolean(study.verified || study.final || study.readiness_verified || readiness === "verified_ready");
      // "likely_ready" reads as a promise. It is an estimate from the record's
      // metadata, and the top-ranked one can hold no usable file at all, so say
      // what the estimate rests on: how many candidate files were actually seen.
      const files = Array.isArray(degAssessment(study).candidate_files) ? degAssessment(study).candidate_files : [];
      // What the record is likely to yield, judged from its file names - the
      // same order preparation ranks what it has opened. An estimate, said as one.
      const likely = detail.likely_input && detail.likely_input !== "no tabular file seen"
        ? `likely ${detail.likely_input} · `
        : "";
      const seen = files.length
        ? `${likely}${files.length} candidate file${files.length === 1 ? "" : "s"}, nothing inspected yet`
        : "nothing inspected yet";
      const headline = verified
        ? "data confirmed"
        : readiness === "likely_ready"
        ? "may have usable data"
        : "not inspected yet";
      const basis = [degAssessment(study).basis, detail.basis, `state: ${readiness}`]
        .filter(Boolean)
        .join(" · ");
      // One phrase on the badge. "data confirmed · nothing inspected yet" put a
      // claim and its caveat in the same pill, which read as a contradiction;
      // the caveat belongs on the line under it, with the relevance figure.
      return `<span class="deg-input ${verified ? "author_deg_likely" : "matrix_fallback"}"`
        + ` title="${esc(basis)}">${esc(headline)}</span>`
        + `<span class="dataset-title readiness-basis">${esc(seen)}</span>`;
    }

    // Source units collapse on a shared PubMed ID, so an unpublished submission
    // deposited as several series counts as several independent studies. That
    // count is what the replication claim rests on, so say it on the row.
    function sharedSubmissionNote(study) {
      const peers = Array.isArray(study.shared_submission_units) ? study.shared_submission_units : [];
      if (!peers.length) return "";
      const warning = study.shared_submission_warning
        || "shares its title with other repository records that are not linked to a publication";
      return `<span class="shared-submission" title="${esc(warning)}">`
        + `May be one submission with ${esc(peers.slice(0, 3).join(", "))}`
        + `${peers.length > 3 ? ` and ${peers.length - 3} more` : ""}</span>`;
    }

    function publicationProvenance(study) {
      const identifiers = [];
      const pmids = Array.isArray(study.pubmed_ids) ? study.pubmed_ids : (study.pmid ? [study.pmid] : []);
      if (pmids.length) identifiers.push(`PMID ${pmids.slice(0, 2).join(", ")}`);
      if (study.doi) identifiers.push(`DOI ${study.doi}`);
      const accessions = Array.isArray(study.geo_accessions) ? study.geo_accessions : (study.accession ? [study.accession] : []);
      if (accessions.length) identifiers.push(accessions.slice(0, 3).join(", "));
      if (study.source_unit_id) identifiers.push(`source unit ${study.source_unit_id}`);
      const conflict = Array.isArray(study.source_unit_conflict) ? study.source_unit_conflict : [];
      if (conflict.length) identifiers.push(`source-unit conflict: ${conflict.join(" / ")}`);
      // The raw decision key reads as "species target species verified"; say it once.
      const speciesWords = {
        target_species_verified: "species verified",
        query_constrained: "species from query",
        mixed_rescued: "mixed series · target species file",
        mixed_species: "mixed species",
      };
      if (study.species_decision) {
        const key = String(study.species_decision);
        identifiers.push(speciesWords[key] || `species ${key.replaceAll("_", " ")}`);
      }
      return identifiers.join(" · ");
    }

    function relevanceText(study) {
      const value = study.relevance ?? study.relevance_score ?? study.score ?? "";
      if (typeof value === "number") return value.toFixed(value > 10 ? 0 : 3);
      return String(value || "—");
    }

    function studySortHead(label, key) {
      const state = activeDiscoveryState();
      const active = state.sort.key === key;
      const direction = active ? (state.sort.order === "asc" ? "ascending" : "descending") : "not sorted";
      return `<button class="sort-head" type="button" data-study-sort="${esc(key)}" aria-label="Sort all assessed studies by ${esc(label)}; ${direction}">${esc(label)} <span class="sort-indicator" aria-hidden="true">${active ? (state.sort.order === "asc" ? "\u25B4" : "\u25BE") : ""}</span></button>`;
    }

    function studyAriaSort(key) {
      const state = activeDiscoveryState();
      if (state.sort.key !== key) return "none";
      return state.sort.order === "asc" ? "ascending" : "descending";
    }

    const STUDY_SORT_LABELS = {
      relevance: "Relevance",
      readiness: "Readiness",
      title: "Publication",
      authors: "Authors",
      journal: "Journal",
      year: "Year",
      data_sources: "Data sources"
    };

    function studyOrderLabel(state) {
      if (state.sort.key === "readiness") return "Sort: DEG readiness · relevance tie-break";
      if (state.sort.key === "relevance") return "Global publication relevance";
      const arrow = state.sort.order === "asc" ? "ascending" : "descending";
      return `${STUDY_SORT_LABELS[state.sort.key] || state.sort.key} ${arrow} · persisted publication snapshot`;
    }

    function renderDiscoveryNotice(state) {
      // A real failure outranks a stale success notice, but a *newer* failure
      // (a rejected prepare) must not be masked by an older search error.
      const notice = state.noticeLevel === "error" && state.notice
        ? state.notice
        : state.error
        ? `Search failed: ${state.error}`
        : (state.notice || "");
      const box = $("discoveryNotice");
      const failed = Boolean(state.error) || /(?:failed|error|unavailable)/i.test(notice);
      box.textContent = notice;
      box.hidden = !notice;
      box.className = failed ? "error-box" : "notice-box";
      box.setAttribute("role", failed ? "alert" : "status");
    }

    // A publication with no PMID/DOI/accession/title collapses to a shared
    // fallback id, so one checkbox would tick every row that shares it and the
    // server rejects the prepare as ambiguous. Mark those rows unselectable.
    function pageSelectability(state) {
      const counts = new Map();
      state.studies.forEach((study) => {
        const key = publicationKey(study);
        counts.set(key, (counts.get(key) || 0) + 1);
      });
      return (key) => {
        if (!key.trim()) return "no-id";
        if ((counts.get(key) || 0) > 1) return "ambiguous-id";
        return "";
      };
    }

    function selectableKeys(state) {
      const blocked = pageSelectability(state);
      return state.studies.map(publicationKey).filter((key) => key && !blocked(key));
    }

    // renderDiscoveryResults() replaces the whole table, which moves focus to
    // <body> and leaves a keyboard user unable to tick a second row.
    function renderDiscoveryResultsKeepingFocus() {
      const active = document.activeElement;
      const accession = active && active.classList && active.classList.contains("study-select")
        ? String(active.dataset.accession || "")
        : "";
      const wasSelectAll = Boolean(active) && active.id === "selectPageStudies";
      renderDiscoveryResults();
      const restored = wasSelectAll
        ? $("selectPageStudies")
        : accession
        ? document.querySelector(`.study-select[data-accession="${CSS.escape(accession)}"]`)
        : null;
      if (restored && !restored.disabled) restored.focus();
    }

    function updatePageSelectionCheckbox(state) {
      const selectAll = $("selectPageStudies");
      if (!selectAll) return;
      const pageAccessions = selectableKeys(state);
      selectAll.checked = pageAccessions.length > 0 && pageAccessions.every((value) => state.selected.has(value));
      selectAll.indeterminate = !selectAll.checked && pageAccessions.some((value) => state.selected.has(value));
      const capBlocksEveryRow = !selectAll.checked
        && state.selected.size >= MAX_SELECTED_STUDIES
        && pageAccessions.every((value) => !state.selected.has(value));
      selectAll.disabled = pageAccessions.length === 0 || capBlocksEveryRow;
      selectAll.title = capBlocksEveryRow
        ? `Selection limit of ${MAX_SELECTED_STUDIES} reached. Clear a selection to choose another.`
        : "";
    }

    // Every control that could start competing work is locked while a search,
    // prepare or analysis is in flight. Disabling matters as much as the dimming:
    // pointer-events alone still leaves the controls reachable by keyboard.
    // Controls whose enabled state nothing else computes: they must be switched
    // back on when the work finishes, or the search box would stay dead forever.
    const BUSY_OWNED_CONTROLS = [
      "discoverySearch", "discoveryQuery", "discoverySearchBoth",
      "humanSpeciesTab", "mouseSpeciesTab", "mobileStudySort", "mobileStudyOrder", "resultFilter",
    ];
    // Controls the render pass already computes: only ever force these off.
    const BUSY_FORCED_CONTROLS = [
      "discoveryPrev", "discoveryNext", "downloadSearchExcel",
      "clearSelected", "prepareSelected", "selectPageStudies", "resetStudySort",
    ];

    function applyDiscoveryBusyState(state) {
      const busy = Boolean(state.loading || state.preparing || state.analyzing);
      const view = $("discoveryView");
      view.classList.toggle("is-busy", busy);
      view.setAttribute("aria-busy", busy ? "true" : "false");
      BUSY_OWNED_CONTROLS.forEach((id) => {
        const element = document.getElementById(id);
        if (element) element.disabled = busy;
      });
      if (!busy) return;
      BUSY_FORCED_CONTROLS.forEach((id) => {
        const element = document.getElementById(id);
        if (element) element.disabled = true;
      });
      document
        .querySelectorAll("#discoveryResults .study-select, #discoveryResults .study-inspect, #discoveryResults .sort-head")
        .forEach((element) => { element.disabled = true; });
    }

    function renderDiscoveryResults() {
      renderDiscoveryResultsView();
      applyDiscoveryBusyState(activeDiscoveryState());
    }

    function renderDiscoveryResultsView() {
      const state = activeDiscoveryState();
      $("discoveryView").classList.toggle("has-results", Boolean(state.query));
      renderDiscoveryHeaderMeta(state);
      renderDiscoveryNotice(state);
      $("speciesState").textContent = `${speciesLabel(activeSpecies)} workspace`;
      $("resultsTitle").textContent = state.query && !state.loading && !state.error
        ? `${state.totalHits.toLocaleString()} studies`
        : "Study results";
      if (state.loading) {
        const percent = typeof state.jobProgress === "number"
          ? Math.max(0, Math.min(100, Math.round(state.jobProgress * 100)))
          : null;
        const stage = state.jobMessage
          || (state.jobStatus === "running" ? "Search running" : "Search queued");
        const elapsed = state.jobStartedAt
          ? Math.max(0, Math.round((Date.now() - state.jobStartedAt) / 1000))
          : 0;
        const elapsedText = elapsed >= 1 ? ` · ${formatElapsed(elapsed)} elapsed` : "";
        const percentText = percent === null ? "" : ` · ${percent}%`;
        // Only the stage and percentage go to the live region: the elapsed
        // second-counter would otherwise be announced once per second.
        $("resultsSubtitle").textContent = `${stage}${percentText}`;
        const bar = percent === null
          ? `<div class="loading-bar" aria-hidden="true"></div>`
          : `<div class="loading-bar is-determinate" aria-hidden="true"><span style="width: ${percent}%"></span></div>`;
        // The card is built once and then updated in place: rebuilding it on every
        // 650 ms poll replaced the Stop button under the reader's pointer, so a
        // press whose mousedown and mouseup landed on different instances was lost.
        const existingCard = $("discoveryResults").querySelector(".search-progress");
        if (existingCard && existingCard.dataset.species === activeSpecies) {
          existingCard.querySelector(".loading-stage").textContent = `${stage}${percentText}${elapsedText}`;
          const fill = existingCard.querySelector(".loading-bar span");
          if (fill && percent !== null) fill.style.width = `${percent}%`;
          if (percent !== null) existingCard.setAttribute("aria-valuenow", String(percent));
          const stopButton = existingCard.querySelector("#cancelSearchJob");
          if (stopButton) {
            stopButton.disabled = Boolean(state.cancelling);
            stopButton.textContent = state.cancelling ? "Stopping..." : "Stop this search";
          }
        } else {
          $("discoveryResults").innerHTML = `<div class="discovery-loading">`
            + `<div class="loading-card search-progress" data-species="${esc(activeSpecies)}" aria-busy="true"`
            + `${percent === null ? "" : ` role="progressbar" aria-valuenow="${percent}" aria-valuemin="0" aria-valuemax="100"`}>`
            + `<strong class="loading-title">Building the ${esc(speciesLabel(activeSpecies))} publication snapshot</strong>`
            + bar
            + `<span class="loading-note loading-stage">${esc(stage)}${esc(percentText)}${esc(elapsedText)}</span>`
            + `<span class="loading-note">Human and Mouse searches run as independent jobs and are never pooled.</span>`
            + `<button type="button" class="job-cancel" id="cancelSearchJob"${state.cancelling ? " disabled" : ""}>`
            + `${state.cancelling ? "Stopping..." : "Stop this search"}</button>`
            + `</div></div>`;
        }
        $("discoveryActions").hidden = true;
        $("discoveryFooter").hidden = true;
        updateSelectedStatus();
        return;
      }
      if (!state.query) {
        $("resultsSubtitle").textContent = "Choose Human or Mouse, enter a keyword, and search.";
        $("discoveryResults").innerHTML = `<div class="discovery-empty">Search papers and linked data to see provisional publication matches.</div>`;
        $("discoveryActions").hidden = true;
        $("discoveryFooter").hidden = true;
        updateSelectedStatus();
        return;
      }
      if (state.error) {
        $("resultsSubtitle").textContent = `Search failed for “${state.query}”: ${state.error}`;
        $("discoveryResults").innerHTML = `<div class="discovery-empty"><strong>Publication discovery could not complete.</strong><br>Keep the query and retry when network services are available.<br><button class="action-secondary study-inspect" type="button" data-retry-search>Retry search</button></div>`;
        $("discoveryActions").hidden = state.selected.size === 0;
        if (!$("discoveryActions").hidden) updateSelectedStatus();
        $("discoveryFooter").hidden = true;
        return;
      }
      // A stopped search and a search that found nothing look identical once the
      // progress bar is gone, and only one of them means the query has no matches.
      if (state.jobCancelled && !state.studies.length) {
        $("resultsSubtitle").textContent = `Search stopped for \u201c${state.query}\u201d`;
        $("discoveryResults").innerHTML = `<div class="discovery-empty">`
          + `<strong>You stopped this search, so it has no results.</strong><br>`
          + `This is not a finding about the query. Anything already downloaded is kept `
          + `and a later run reuses it.<br>`
          + `<button class="action-secondary study-inspect" type="button" data-retry-search>Search again</button></div>`;
        $("discoveryActions").hidden = true;
        $("discoveryFooter").hidden = true;
        updateSelectedStatus();
        return;
      }
      const verification = state.providerStatus === "partial"
        ? "partial snapshot · some sources unavailable"
        : state.verified
        ? "complete snapshot"
        : "search in progress";
      $("resultsSubtitle").textContent = `${speciesLabel(activeSpecies)} only · author DEG tables preferred · ≥2 source units · ${verification}`;
      $("discoveryActions").hidden = false;
      $("studyOrderStatus").textContent = studyOrderLabel(state);
      $("resetStudySort").hidden = state.sort.key === "readiness";
      $("downloadSearchExcel").disabled = !state.searchId || !state.verified;
      if (!state.studies.length) {
        const degraded = state.providerStatus === "partial" && state.providerErrors.length;
        $("discoveryResults").innerHTML = degraded
          ? `<div class="discovery-empty"><strong>Some data sources did not answer, so this result set is incomplete.</strong><br>`
            + `Unavailable: ${esc(state.providerErrors.join(", "))}.<br>`
            + `Retry when those services are reachable before concluding that the query has no matches.<br>`
            + `<button class="action-secondary study-inspect" type="button" data-retry-search>Retry search</button></div>`
          : `<div class="discovery-empty">No ${esc(speciesLabel(activeSpecies))} publication records were returned on this page.</div>`;
      } else {
        const blockedReason = pageSelectability(state);
        const outcomes = preparedOutcomes(state);
        const capReached = state.selected.size >= MAX_SELECTED_STUDIES;
        const rows = state.studies.map((study) => {
          const key = publicationKey(study);
          const isSelected = state.selected.has(key);
          const checked = isSelected ? "checked" : "";
          const blocked = blockedReason(key);
          // Show unavailable rows as disabled instead of silently rejecting the
          // click: at the cap, and for rows without a usable unique identifier.
          const disabledReason = blocked === "no-id"
            ? "This publication has no usable identifier, so it cannot be prepared."
            : blocked === "ambiguous-id"
            ? "Several results share this identifier, so a single one cannot be prepared."
            : !isSelected && capReached
            ? `Selection limit of ${MAX_SELECTED_STUDIES} reached. Clear a selection to choose another.`
            : "";
          const selectAttrs = disabledReason
            ? ` disabled title="${esc(disabledReason)}"`
            : "";
          const paperHref = study.pubmed_url || study.source_url || study.url || study.doi_url;
          const title = study.title || study.paper_title || study.dataset_title || key || "Untitled publication";
          const authors = asListText(study.authors_display || study.authors);
          const dataSources = asListText(study.data_sources || study.linked_datasets || study.datasets);
          const provenance = publicationProvenance(study);
          const publicationMeta = [authors, study.journal, study.year].filter(Boolean).map(esc).join(" · ");
          return `<tr>
            <td><input class="study-select" type="checkbox" data-accession="${esc(key)}" aria-label="Select ${esc(title)}" ${checked}${selectAttrs}></td>
            <td>${externalLink(paperHref, title, "study-title")}<span class="study-publication-meta">${publicationMeta || "Publication metadata unavailable"}</span><span class="dataset-title">${esc(provenance)}</span>${sharedSubmissionNote(study)}</td>
            <td>${esc(authors || "—")}</td>
            <td>${esc(study.journal || "—")}</td>
            <td>${esc(study.year || "—")}</td>
            <td><span class="mobile-field-label">Linked data</span>${esc(dataSources || "—")}</td>
            <td><span class="mobile-field-label">DEG readiness</span>${readinessBadge(study, outcomes)}<span class="dataset-title">relevance ${esc(relevanceText(study))}</span></td>
            <td class="inspect-cell"><button class="action-secondary study-inspect" type="button" data-study-inspect="${esc(key)}" aria-label="Inspect DEG inputs for ${esc(title)}">Inspect</button></td>
          </tr>`;
        }).join("");
        const sortOptions = [
          ["readiness", "Readiness"], ["relevance", "Relevance"],
          ["title", "Publication"], ["authors", "Authors"], ["journal", "Journal"],
          ["year", "Year"], ["data_sources", "Data sources"]
        ].map(([key, label]) => `<option value="${key}" ${state.sort.key === key ? "selected" : ""}>${label}</option>`).join("");
        // A disabled checkbox fires no click, so the notice the click handler
        // sets can never appear - the reader presses "select all" on page 2 and
        // gets silence. Say it before they press anything.
        const selectableHere = selectableKeys(state).filter((value) => !state.selected.has(value)).length;
        const limitBanner = capReached && selectableHere > 0
          ? `<div class="selection-limit">Selection limit reached: ${MAX_SELECTED_STUDIES} publications across all pages`
            + ` (${state.selected.size - state.studies.reduce((total, study) => total + (state.selected.has(publicationKey(study)) ? 1 : 0), 0)} of them on other pages).`
            + ` Untick a row or press Clear to choose different ones here.</div>`
          : "";
        // Narrowing a thousand-record snapshot used to mean running the whole
        // search again against live providers. This filters what is already here.
        const filterCount = state.textFilter
          ? `<span class="filter-count">${state.totalHits.toLocaleString()} of ${state.totalUnfiltered.toLocaleString()} match</span>`
          : "";
        const filterBar = `<div class="result-filter"><label for="resultFilter">Narrow these results</label>`
          + `<input id="resultFilter" type="search" placeholder="title, author, journal or year" `
          + `maxlength="100" value="${esc(state.textFilter)}" `
          + `data-tip="Filters the records already found. It does not run a new search.">`
          + `${filterCount}</div>`;
        $("discoveryResults").innerHTML = limitBanner + filterBar + `<div class="mobile-study-tools"><label>Sort persisted records<select id="mobileStudySort">${sortOptions}</select></label><button id="mobileStudyOrder" type="button" aria-label="Reverse global sort order">${state.sort.order === "asc" ? "Ascending" : "Descending"}</button></div><div class="results-scroll"><table class="study-table">
          <thead><tr>
            <th><input id="selectPageStudies" type="checkbox" aria-label="Select all studies on this page"></th>
            <th class="sortable" aria-sort="${studyAriaSort("title")}">${studySortHead("Publication", "title")}</th>
            <th class="sortable" aria-sort="${studyAriaSort("authors")}">${studySortHead("Authors", "authors")}</th>
            <th class="sortable" aria-sort="${studyAriaSort("journal")}">${studySortHead("Journal", "journal")}</th>
            <th class="sortable" aria-sort="${studyAriaSort("year")}">${studySortHead("Year", "year")}</th>
            <th class="sortable" aria-sort="${studyAriaSort("data_sources")}">${studySortHead("Linked data", "data_sources")}</th>
            <th class="sortable" aria-sort="${studyAriaSort("readiness")}">${studySortHead("DEG readiness", "readiness")}</th>
            <th>Inspect</th>
          </tr></thead><tbody>${rows}</tbody></table></div>`;
        updatePageSelectionCheckbox(state);
      }
      $("discoveryFooter").hidden = false;
      const totalPages = Math.max(state.totalPages, 1);
      const firstShown = state.studies.length ? ((state.page - 1) * DISCOVERY_PAGE_SIZE) + 1 : 0;
      const lastShown = state.studies.length ? firstShown + state.studies.length - 1 : 0;
      $("discoveryPage").textContent = `Showing ${firstShown.toLocaleString()}–${lastShown.toLocaleString()} of ${state.totalHits.toLocaleString()} · Page ${state.page.toLocaleString()} / ${totalPages.toLocaleString()}`;
      $("discoveryPrev").disabled = state.page <= 1;
      $("discoveryNext").disabled = !state.hasNext;
      updateSelectedStatus();
    }

    function updateSelectedStatus() {
      const state = activeDiscoveryState();
      const count = state.selected.size;
      // The cap counts across every page, and the page size equals the cap, so
      // "20 / 20" used to read as "all 20 rows here" while the selection in fact
      // lived on another page. Say where the selection actually is.
      const onPage = state.studies.reduce(
        (total, study) => total + (state.selected.has(publicationKey(study)) ? 1 : 0),
        0,
      );
      const scope = count > 0 && state.studies.length ? ` · ${onPage} on this page` : "";
      $("selectedStatus").textContent =
        `${count.toLocaleString()} selected of max ${MAX_SELECTED_STUDIES}${scope} · ${speciesLabel(activeSpecies)}`;
      $("selectedStatus").title = count > 0 && onPage < count
        ? `${count - onPage} selected publication(s) are on other pages and will be included in Prepare selection.`
        : "";
      const busy = Boolean(state.loading || state.preparing || state.analyzing);
      $("prepareSelected").textContent = state.preparing ? "Downloading and inspecting..." : "Prepare selection";
      $("prepareSelected").disabled = count === 0 || busy;
      $("clearSelected").disabled = count === 0 || busy;
    }

    async function refreshSearchPage(species, requestId) {
      const state = discoveryStates[species];
      if (!state.searchId) return;
      const params = new URLSearchParams({
        page: String(state.page),
        page_size: String(DISCOVERY_PAGE_SIZE),
        sort_by: state.sort.key,
        sort_order: state.sort.order
      });
      if (state.textFilter) params.set("filter", state.textFilter);
      const data = await getJson(`/api/discovery/searches/${state.searchId}/records?${params.toString()}`);
      if (requestId !== state.searchRequest) return;
      state.studies = data.records || [];
      state.totalHits = Number(data.total || data.search?.total || 0);
      state.totalPages = Number(data.total_pages || 0);
      state.totalUnfiltered = Number(data.total_unfiltered || data.total || 0);
      state.evaluatedStudies = state.totalUnfiltered;
      state.rankingLimit = Number(data.search?.limit || DISCOVERY_GLOBAL_RANK_LIMIT);
      state.rankingTruncated = state.totalHits > state.rankingLimit;
      state.cacheHit = true;
      state.hasNext = Boolean(data.has_next);
      state.page = Number(data.page || state.page);
      state.verified = data.search?.status === "complete";
      // The server already reports a degraded provider set, but nothing read it,
      // so an outage was indistinguishable from "this query has no results".
      const snapshot = data.search?.snapshot || data;
      state.providerStatus = String(snapshot.provider_status || "");
      const providerErrors = ((snapshot.diagnostics || {}).errors || [])
        .map((entry) => `${entry.provider || "a source"} (${entry.stage || "search"})`);
      state.providerErrors = [...new Set(providerErrors)];
      state.error = "";
    }

    // Stopping a job is a request to the server, not just a client-side give-up:
    // the worker is mid-download and would otherwise keep fetching from public
    // repositories after the reader has walked away from the answer.
    async function cancelDiscoveryJob(kind) {
      const species = activeSpecies;
      const state = discoveryStates[species];
      const jobId = kind === "prepare" ? state.prepareJobId : state.jobId;
      if (!jobId || state.cancelling) return;
      state.cancelling = true;
      if (activeSpecies === species) {
        kind === "prepare" ? renderPreparedState() : renderDiscoveryResults();
      }
      let outcome = null;
      try {
        outcome = await postJson(`/api/discovery/jobs/${jobId}/cancel`, {});
      } catch (error) {
        state.cancelling = false;
        state.notice = `Could not stop the job: ${error.message}`;
        state.noticeLevel = "error";
        if (activeSpecies === species) {
          renderDiscoveryResults();
          renderPreparedState();
        }
        return;
      }
      state.cancelling = false;
      if (outcome && outcome.cancelled === false) {
        const job = outcome.job || {};
        const notice = outcome.reason || "The job finished before it could be stopped.";
        if (kind === "search" && job.status === "complete") {
          state.searchRequest += 1;
          const adoptRequest = state.searchRequest;
          if (job.search_id) state.searchId = job.search_id;
          state.loading = false;
          state.jobCancelled = false;
          state.notice = notice;
          state.noticeLevel = "info";
          try {
            await refreshSearchPage(species, adoptRequest);
          } catch (error) {
            state.error = error.message;
            state.notice = `Search completed, but the saved result could not be loaded: ${error.message}`;
            state.noticeLevel = "error";
          }
          if (activeSpecies === species) renderDiscoveryResults();
          return;
        }
        if (kind === "prepare" && job.status === "complete") {
          state.prepareRequest += 1;
          state.preparing = false;
          state.prepareCancelled = false;
          state.notice = notice;
          state.noticeLevel = "info";
          if (job.result) {
            state.prepared = job.result;
            state.bundleId = job.result.bundle_id || state.bundleId;
          }
          if (activeSpecies === species) {
            renderDiscoveryResults();
            renderPreparedState();
          }
          return;
        }
        state.notice = notice;
        state.noticeLevel = "info";
        if (activeSpecies === species) {
          renderDiscoveryResults();
          renderPreparedState();
        }
        return;
      }
      // Retire the in-flight work client-side too. The poll loop reads the
      // cancelled status on its next tick, but the reader should not watch a
      // progress bar keep moving for another 650ms after pressing stop.
      if (kind === "prepare") {
        state.prepareRequest += 1;
        state.preparing = false;
        state.prepareCancelled = true;
      } else {
        state.searchRequest += 1;
        state.loading = false;
        state.jobCancelled = true;
      }
      state.notice = outcome && outcome.cancelled === false
        ? "The job finished before it could be stopped."
        : (outcome && outcome.reason) || "Stopped.";
      state.noticeLevel = outcome && outcome.cancelled === false ? "info" : "info";
      if (activeSpecies === species) {
        renderDiscoveryResults();
        renderPreparedState();
      }
    }

    async function pollSearchJob(species, requestId) {
      const state = discoveryStates[species];
      while (requestId === state.searchRequest && state.jobId) {
        const payload = await getJson(`/api/discovery/jobs/${state.jobId}`);
        if (requestId !== state.searchRequest) return;
        const job = payload.job || payload;
        state.jobStatus = job.status || "";
        state.jobProgress = typeof job.progress === "number" ? job.progress : null;
        state.jobMessage = typeof job.message === "string" ? job.message : "";
        if (!state.jobStartedAt) {
          const started = Date.parse(job.created_at || "");
          state.jobStartedAt = Number.isFinite(started) ? started : Date.now();
        }
        if (job.status === "complete") {
          await refreshSearchPage(species, requestId);
          return;
        }
        // Cancelled is not a failure: someone asked for this. It can also arrive
        // from another tab, so the loop recognises the state rather than relying
        // on the button that usually causes it.
        if (job.status === "cancelled") {
          state.jobCancelled = true;
          return;
        }
        if (job.status === "failed" || job.status === "interrupted") {
          throw new Error(job.error || "publication search failed");
        }
        if (activeSpecies === species) renderDiscoveryResults();
        await new Promise((resolve) => setTimeout(resolve, 650));
      }
    }

    async function searchStudies({ resetPage = false, speciesOverride = "" } = {}) {
      const requestSpecies = speciesOverride || activeSpecies;
      const state = discoveryStates[requestSpecies];
      const input = $("discoveryQuery");
      const query = input.value.trim();
      if (query.length < 2) {
        input.setCustomValidity(query ? "Enter at least 2 characters." : "Enter a condition, perturbation, disease, or pathway.");
        input.reportValidity();
        input.focus();
        return;
      }
      input.setCustomValidity("");
      const queryChanged = query !== state.query;
      if (resetPage || queryChanged) state.page = 1;
      // Any new search re-mints the snapshot the bundle was built from, so the
      // bundle, its review draft and its run stop applying - not only when the
      // query text changed. Bumping the request ids retires work that is still
      // in flight; without it a late response reinstated a bundle built from
      // the previous snapshot, and Run then launched against a stale selection.
      state.prepareRequest += 1;
      state.analysisRequest += 1;
      state.preparing = false;
      state.analyzing = false;
      state.jobCancelled = false;
      state.prepareCancelled = false;
      state.selected.clear();
      state.prepared = null;
      state.bundleId = "";
      state.run = null;
      state.draft = {};
      state.cloneCounter = 0;
      state.analysisError = "";
      if (queryChanged) {
        state.sort = { key: "readiness", order: "desc" };
      }
      if (resetPage || queryChanged) {
        state.textFilter = "";
        state.totalUnfiltered = 0;
      }
      if (activeSpecies === requestSpecies) {
        activeRunId = "";
        invalidateAtlasContext();
        renderPreparedState();
      }
      state.query = query;
      if (resetPage || queryChanged) {
        // The rows of the previous snapshot must not outlive it: when a re-run
        // search was stopped, they came back under the new search id, looked
        // like an answer to the query in the box, and could not be prepared.
        state.studies = [];
        state.totalHits = 0;
        state.totalPages = 0;
        state.hasNext = false;
        state.evaluatedStudies = 0;
        state.rankingTruncated = false;
        state.cacheHit = false;
      }
      state.loading = true;
      state.error = "";
      state.notice = "";
      state.noticeLevel = "info";
      state.verified = false;
      state.jobStatus = "queued";
      state.jobProgress = null;
      state.jobMessage = "";
      state.jobStartedAt = Date.now();
      const requestId = ++state.searchRequest;
      const requestPage = state.page;
      if (activeSpecies === requestSpecies) renderDiscoveryResults();
      try {
        if (!resetPage && state.searchId && !queryChanged) {
          await refreshSearchPage(requestSpecies, requestId);
        } else {
          const data = await postJson("/api/discovery/searches", {
            query,
            species: requestSpecies,
            limit: DISCOVERY_GLOBAL_RANK_LIMIT
          });
          if (requestId !== state.searchRequest) return;
          // A new snapshot re-mints record ids, so a selection carried over from
          // the previous one points at rows that exist on no page and fails the
          // prepare server-side while still consuming cap slots.
          state.selected.clear();
          state.searchId = data.search_id;
          state.jobId = data.job_id;
          state.jobStatus = data.status || "queued";
          state.page = requestPage;
          await pollSearchJob(requestSpecies, requestId);
        }
      } catch (error) {
        if (requestId !== state.searchRequest) return;
        state.studies = [];
        state.totalHits = 0;
        state.totalPages = 0;
        state.evaluatedStudies = 0;
        state.rankingTruncated = false;
        state.cacheHit = false;
        state.hasNext = false;
        state.verified = false;
        state.error = error.message;
        if (activeSpecies === requestSpecies) renderDiscoveryResults();
      } finally {
        if (requestId !== state.searchRequest) return;
        state.loading = false;
        if (activeSpecies === requestSpecies) renderDiscoveryResults();
      }
    }

    let resultFilterTimer = 0;

    // The panel re-renders by replacing innerHTML, which destroys the element the
    // reader is typing into. Put the caret back where it was, or the filter loses
    // focus on the first keystroke and every one after it.
    function restoreFocus(elementId) {
      if (!elementId) return;
      const control = $(elementId);
      if (!control) return;
      const end = control.value.length;
      control.focus();
      try {
        control.setSelectionRange(end, end);
      } catch (_) {
        // Some input types refuse a selection range; focus alone is enough.
      }
    }

    async function reloadSearchRecords({ resetPage = false, keepFocus = "" } = {}) {
      const state = activeDiscoveryState();
      if (resetPage) state.page = 1;
      if (!state.searchId) {
        await searchStudies({ resetPage: true });
        return;
      }
      state.loading = true;
      state.error = "";
      // Without this the panel replayed the previous job's "Job completed. 100%"
      // caption, with an elapsed counter that grew for the tab's lifetime.
      state.jobStatus = "";
      state.jobProgress = null;
      state.jobMessage = "";
      state.jobStartedAt = Date.now();
      const requestId = ++state.searchRequest;
      renderDiscoveryResults();
      restoreFocus(keepFocus);
      try {
        await refreshSearchPage(activeSpecies, requestId);
      } catch (error) {
        if (requestId !== state.searchRequest) return;
        state.error = error.message;
      } finally {
        if (requestId !== state.searchRequest) return;
        state.loading = false;
        renderDiscoveryResults();
        restoreFocus(keepFocus);
      }
    }

    async function searchBothSpecies() {
      await Promise.all([
        searchStudies({ resetPage: true, speciesOverride: "human" }),
        searchStudies({ resetPage: true, speciesOverride: "mouse" })
      ]);
    }

    function setSpecies(species) {
      capturePreparedDraft();
      activeSpecies = species;
      activeRunId = activeDiscoveryState().run?.run_id || "";
      invalidateAtlasContext();
      document.querySelectorAll("[data-species]").forEach((button) => {
        const active = button.dataset.species === species;
        button.classList.toggle("active", active);
        button.setAttribute("aria-pressed", active ? "true" : "false");
        button.setAttribute("aria-selected", active ? "true" : "false");
        button.setAttribute("tabindex", active ? "0" : "-1");
      });
      $("discoveryResultsCard").setAttribute("aria-labelledby", `${species}SpeciesTab`);
      const state = activeDiscoveryState();
      $("discoveryQuery").value = state.query || state.suggestedQuery;
      renderDiscoveryResults();
      renderPreparedState();
    }

    async function loadDiscoveryDefaults() {
      try {
        const payload = await getJson("/api/meta");
        const meta = payload && payload.meta ? payload.meta : {};
        const keyword = String(meta.demo_search_keyword || "").trim();
        const species = ["human", "mouse"].includes(String(meta.demo_search_species || "").toLowerCase())
          ? String(meta.demo_search_species).toLowerCase() : "human";
        preferredDiscoverySpecies = species;
        if (keyword) discoveryStates[species].suggestedQuery = keyword;
        const state = activeDiscoveryState();
        if (!state.query && !$("discoveryQuery").value.trim()) {
          $("discoveryQuery").value = state.suggestedQuery;
        }
        return meta;
      } catch (_) {
        // Demo search defaults are optional; the atlas and manual search remain usable.
        return {};
      }
    }

    function capturePreparedDraft() {
      const state = activeDiscoveryState();
      if (!state.prepared || $("preparedCard").hidden) return;
      const draft = {};
      $("preparedCandidates").querySelectorAll(".candidate-row").forEach((row) => {
        const key = row.dataset.activationKey || row.dataset.candidateId;
        const samples = {};
        row.querySelectorAll("[data-sample]").forEach((select) => { samples[select.dataset.sample] = select.value; });
        draft[key] = {
          enabled: Boolean(row.querySelector(".candidate-enable")?.checked),
          candidateId: row.dataset.candidateId,
          clone: row.dataset.clone === "true",
          contrast: row.querySelector(".contrast-label")?.value || "",
          direction: Boolean(row.querySelector(".direction-confirmed")?.checked),
          biologicalReplicates: Boolean(row.querySelector(".biological-replicates-confirmed")?.checked),
          tableScope: row.querySelector(".table-scope")?.value || "auto",
          cellSystem: row.querySelector(".cell-system")?.value || "",
          durationH: row.querySelector(".duration-h")?.value || "",
          nCtrl: row.querySelector(".n-ctrl")?.value || "",
          nTreat: row.querySelector(".n-treat")?.value || "",
          platform: row.querySelector(".platform")?.value || "",
          assayType: row.querySelector(".assay-type")?.value || "",
          pipeline: row.querySelector(".pipeline")?.value || "",
          sheetName: row.querySelector(".sheet-name")?.value || "",
          geneColumn: row.querySelector(".gene-column")?.value || "",
          lfcColumn: row.querySelector(".lfc-column")?.value || "",
          pColumn: row.querySelector(".p-column")?.value || "",
          padjColumn: row.querySelector(".padj-column")?.value || "",
          columnMappingConfirmed: Boolean(row.querySelector(".column-mapping-confirmed")?.checked),
          adjustedPAsPvalueConfirmed: Boolean(row.querySelector(".adjusted-p-as-pvalue-confirmed")?.checked),
          lfcScaleConfirmedLog2: Boolean(row.querySelector(".lfc-scale-confirmed-log2")?.checked),
          rowFilterColumn: row.querySelector(".row-filter-column")?.value || "",
          rowFilterValue: row.querySelector(".row-filter-value")?.value || "",
          rowFilterConfirmed: Boolean(row.querySelector(".row-filter-confirmed")?.checked),
          duplicateGenePolicy: row.querySelector(".duplicate-gene-policy")?.value || "harmonizer",
          duplicateGenePolicyConfirmed: Boolean(row.querySelector(".duplicate-gene-policy-confirmed")?.checked),
          matrixType: row.querySelector(".matrix-type")?.value || "",
          normalizedScale: row.querySelector(".normalized-scale")?.value || "",
          samples
        };
      });
      state.draft = draft;
    }

    function isPositiveWholeNumber(value) {
      return /^[1-9][0-9]*$/.test(String(value || "").trim());
    }

    const AUTHOR_REVIEWABLE_STATUSES = new Set([
      "ready_for_review",
      "requires_column_mapping",
      "requires_lfc_confirmation",
      "requires_pvalue_mapping"
    ]);

    function isAuthorReviewable(candidate) {
      const status = candidate.inspection?.status || "";
      return AUTHOR_REVIEWABLE_STATUSES.has(status);
    }

    // "No usable table was resolved within the safety limits" names no file and
    // no limit, so a reader cannot tell a study that published only browser
    // tracks from one whose table was a megabyte over the cap. List what was
    // actually there.
    function candidateRejection(candidate) {
      const inspection = candidate.inspection || {};
      const status = inspection.status || "";
      // A rejected file was never inspected, so its inspection note only says
      // it was skipped. What the reader needs is why it was never a candidate.
      if (candidate.tier === "reject" || candidate.role === "unsupported") {
        return candidate.reason || "not a DEG table or expression matrix";
      }
      if (status && status !== "not_inspected") return inspection.reason || status;
      return inspection.reason || candidate.reason || "not inspected within the safety limits";
    }

    function candidateLabel(candidate) {
      if (candidate.name) return candidate.name;
      if (candidate.member_name) return candidate.member_name;
      const url = String(candidate.source_url || "");
      const tail = url.split("?")[0].split("/").filter(Boolean).pop();
      return tail || "unnamed file";
    }

    function unusableStudyHtml(study) {
      const files = study.files || [];
      if (!files.length) {
        return `<div class="candidate-note">The repository listed no supplementary file for this study.</div>`;
      }
      const shown = files.slice(0, 8);
      const rest = files.length - shown.length;
      const items = shown
        .map((candidate) => `<li><code>${esc(candidateLabel(candidate))}</code> — ${esc(candidateRejection(candidate))}</li>`)
        .join("");
      return `<div class="candidate-note unusable-files">`
        + `<strong>${files.length} file${files.length === 1 ? "" : "s"} were found and none could be used:</strong>`
        + `<ul>${items}${rest > 0 ? `<li>and ${rest} more</li>` : ""}</ul>`
        + `Pick another study, or open the repository page to check whether a DEG table exists elsewhere.</div>`;
    }

    function eligibleCandidate(candidate) {
      const status = candidate.inspection?.status || "";
      return isAuthorReviewable(candidate) || status === "upstream_matrix_ready_for_contrast";
    }

    function detectedAuthorValue(candidate, field) {
      const inspection = candidate.inspection || {};
      const mapping = inspection.mapping || {};
      return field === "sheet_name" ? (inspection.sheet_name || "") : (mapping[field] || "");
    }

    function authorDraftKeys(candidateId) {
      return Object.keys(activeDiscoveryState().draft).filter((key) => {
        const item = activeDiscoveryState().draft[key] || {};
        return item.candidateId === candidateId || key === candidateId || key.startsWith(`${candidateId}::clone::`);
      });
    }

    // A control that already carries a value must not be hidden behind a closed
    // panel: the reader would not see what their run is actually going to do.
    function anyValueSet(...values) {
      return values.some((value) => String(value ?? "").trim() !== "");
    }

    function authorCandidateHtml(study, candidate, defaultChecked, activationKey, clone) {
      activationKey = activationKey || candidate.candidate_id;
      const mapping = candidate.inspection?.mapping || {};
      const inspection = candidate.inspection || {};
      const status = inspection.status || "";
      const draft = activeDiscoveryState().draft[activationKey] || {};
      const checked = draft.enabled === undefined ? defaultChecked : draft.enabled;
      const contrast = draft.contrast === undefined ? "" : draft.contrast;
      const cellSystem = draft.cellSystem === undefined ? "" : draft.cellSystem;
      const durationH = draft.durationH === undefined ? "" : draft.durationH;
      const nCtrl = draft.nCtrl === undefined ? "" : draft.nCtrl;
      const nTreat = draft.nTreat === undefined ? "" : draft.nTreat;
      const platform = draft.platform === undefined ? "" : draft.platform;
      const assayType = draft.assayType === undefined ? "" : draft.assayType;
      const pipeline = draft.pipeline === undefined ? "" : draft.pipeline;
      const sheetName = draft.sheetName === undefined ? (inspection.sheet_name || "") : draft.sheetName;
      const geneColumn = draft.geneColumn === undefined ? (mapping.gene_column || "") : draft.geneColumn;
      const lfcColumn = draft.lfcColumn === undefined ? (mapping.lfc_column || "") : draft.lfcColumn;
      const pColumn = draft.pColumn === undefined ? (mapping.p_column || "") : draft.pColumn;
      const padjColumn = draft.padjColumn === undefined ? (mapping.padj_column || "") : draft.padjColumn;
      const rowFilterColumn = draft.rowFilterColumn === undefined ? "" : draft.rowFilterColumn;
      const rowFilterValue = draft.rowFilterValue === undefined ? "" : draft.rowFilterValue;
      const duplicateGenePolicy = draft.duplicateGenePolicy || "harmonizer";
      // Open the column panel when the reader has to look at it: the file was not
      // read cleanly, or they have already changed something in it.
      // The split between groups is not derivable from a results table - it has one
      // row per gene - and it feeds the source weight directly, so it is never
      // guessed. What the linked series holds in total is knowable, and is the
      // number a reader is otherwise squinting at the paper to find.
      const seriesSamples = Number(study.n_samples || 0);
      const seriesSampleNote = seriesSamples > 0
        ? `<p class="candidate-note">The linked series lists ${seriesSamples} sample${seriesSamples === 1 ? "" : "s"} in total. `
          + `Enter the numbers for this contrast only; together they cannot exceed ${seriesSamples}.</p>`
        : "";
      const columnsOpen = status !== "ready_for_review"
        || anyValueSet(draft.sheetName, draft.geneColumn, draft.lfcColumn, draft.pColumn, draft.padjColumn);
      // And never hide a setting that is already carrying a value, or the reader
      // cannot see what their run is actually going to do.
      const advancedOpen = anyValueSet(
        rowFilterColumn,
        rowFilterValue,
        assayType,
        pipeline,
        cellSystem,
        durationH,
        platform
      ) || duplicateGenePolicy === "keep_first";
      return `<div class="candidate-row" data-candidate-id="${esc(candidate.candidate_id)}" data-activation-key="${esc(activationKey)}" data-clone="${clone ? "true" : "false"}" data-accession="${esc(study.accession)}" data-source-unit="${esc(study.source_unit_id || study.accession)}" data-mode="author" data-author-status="${esc(status)}" data-series-samples="${esc(String(seriesSamples || ""))}" data-detected-sheet-name="${esc(detectedAuthorValue(candidate, "sheet_name"))}" data-detected-gene-column="${esc(detectedAuthorValue(candidate, "gene_column"))}" data-detected-lfc-column="${esc(detectedAuthorValue(candidate, "lfc_column"))}" data-detected-p-column="${esc(detectedAuthorValue(candidate, "p_column"))}" data-detected-padj-column="${esc(detectedAuthorValue(candidate, "padj_column"))}">
        <input class="candidate-enable" type="checkbox" aria-label="Use ${esc(candidate.name)}" ${checked ? "checked" : ""}>
        <div><span class="candidate-name">${esc(candidate.name)}${clone ? " · additional cohort/contrast" : ""}</span><span class="candidate-note">Author DEG candidate · ${esc(status)} · ${esc(geneColumn)} / ${esc(lfcColumn)} / ${esc(pColumn)}</span><div class="candidate-tools"><span class="status-pill">${status === "ready_for_review" ? "One-click mapping accepted" : "Mapping review required"}</span><button class="action-secondary clone-author-candidate" type="button">Add cohort/contrast</button></div></div>
        <div class="candidate-fields">
          <label>Contrast label<input class="contrast-label" value="${esc(contrast)}" placeholder="Enter exact paper/GEO contrast" autocomplete="off"></label>
          <label>Table scope<select class="table-scope"><option value="auto" ${draft.tableScope === "auto" || draft.tableScope === undefined ? "selected" : ""}>Auto review</option><option value="full_results" ${draft.tableScope === "full_results" ? "selected" : ""}>Full results</option><option value="deg_only" ${draft.tableScope === "deg_only" ? "selected" : ""}>Significant genes only</option></select></label>
          <label>Control biological n<input class="n-ctrl" type="number" min="1" step="1" inputmode="numeric" value="${esc(nCtrl)}" required title="Independent biological control samples in this exact contrast, not total GEO samples"></label>
          <label>Treatment biological n<input class="n-treat" type="number" min="1" step="1" inputmode="numeric" value="${esc(nTreat)}" required title="Independent biological treatment/case samples in this exact contrast, not total GEO samples"></label>
        </div>
        ${seriesSampleNote}
        <details class="candidate-advanced" ${columnsOpen ? "open" : ""}>
          <summary>Columns DEGORA read from the file${columnsOpen ? "" : " - detected, open to change"}</summary>
        <div class="candidate-fields">
          <label>Sheet name<input class="sheet-name" value="${esc(sheetName)}" maxlength="160" autocomplete="off"></label>
          <label>Gene column<input class="gene-column" value="${esc(geneColumn)}" maxlength="160" autocomplete="off"></label>
          <label>log2FC column<input class="lfc-column" value="${esc(lfcColumn)}" maxlength="160" autocomplete="off"></label>
          <label>P-value column<input class="p-column" value="${esc(pColumn)}" maxlength="160" autocomplete="off"></label>
          <label>Adjusted-p column<input class="padj-column" value="${esc(padjColumn)}" maxlength="160" autocomplete="off"></label>
        </div>
        </details>
        <details class="candidate-advanced" ${advancedOpen ? "open" : ""}>
          <summary>Advanced settings${advancedOpen ? "" : " - none set"}</summary>
        <div class="candidate-fields candidate-fields-wide">
          <label>Row filter column<input class="row-filter-column" value="${esc(rowFilterColumn)}" maxlength="160" autocomplete="off"></label>
          <label>Row filter value<input class="row-filter-value" value="${esc(rowFilterValue)}" maxlength="240" autocomplete="off"></label>
          <label>Duplicate genes<select class="duplicate-gene-policy"><option value="harmonizer" ${duplicateGenePolicy === "harmonizer" ? "selected" : ""}>DEGORA default: lowest p, then largest |log2FC|</option><option value="keep_first" ${duplicateGenePolicy === "keep_first" ? "selected" : ""}>Legacy/manual: keep first source row</option></select></label>
          <label>Assay type<input class="assay-type" value="${esc(assayType)}" maxlength="80" autocomplete="off"></label>
          <label>Pipeline<input class="pipeline" value="${esc(pipeline)}" maxlength="80" autocomplete="off"></label>
          <label>Cell or tissue system (optional)<input class="cell-system" value="${esc(cellSystem)}" maxlength="160" autocomplete="off"></label>
          <label>Duration h (optional)<input class="duration-h" value="${esc(durationH)}" maxlength="32" placeholder="e.g. 24 or 24-48" autocomplete="off"></label>
          <label>Platform (optional)<input class="platform" value="${esc(platform)}" maxlength="80" autocomplete="off"></label>
        </div>
        </details>
        <div class="candidate-fields candidate-confirms">
          <label class="confirm-line"><input class="direction-confirmed" type="checkbox" ${draft.direction ? "checked" : ""}> A positive value here means the gene went UP in the treated group, not the control group</label>
          <label class="confirm-line"><input class="column-mapping-confirmed" type="checkbox" ${draft.columnMappingConfirmed ? "checked" : ""}> The gene, effect and p-value columns chosen above are the right ones</label>
          <label class="confirm-line"><input class="adjusted-p-as-pvalue-confirmed" type="checkbox" ${draft.adjustedPAsPvalueConfirmed ? "checked" : ""}> This table has no separate raw p-value, so its adjusted p-value is being used as one</label>
          <label class="confirm-line"><input class="lfc-scale-confirmed-log2" type="checkbox" ${draft.lfcScaleConfirmedLog2 ? "checked" : ""}> The effect column is already a log2 fold change, not a plain ratio</label>
          <label class="confirm-line"><input class="row-filter-confirmed" type="checkbox" ${draft.rowFilterConfirmed ? "checked" : ""}> The filter above picks one comparison and does not mix several together</label>
          <label class="confirm-line"><input class="duplicate-gene-policy-confirmed" type="checkbox" ${draft.duplicateGenePolicyConfirmed ? "checked" : ""}> Keeping the first row for each repeated gene reproduces how this table was originally read</label>
        </div>
      </div>`;
    }

    // Twenty-odd samples assigned one dropdown at a time is where attention
    // runs out. Filtering by the label text keeps the choice explicit and
    // reviewable - the reader still sees which rows moved - while removing the
    // clicking. Any bulk move clears the direction attestations, because they
    // were made about a different assignment.
    function sampleBulkHtml() {
      return `<div class="sample-bulk">
        <input class="sample-filter" type="search" placeholder="Filter by label, e.g. uninduced" aria-label="Filter samples by label" autocomplete="off">
        <button class="action-secondary sample-bulk-apply" type="button" data-group="control">Set Control</button>
        <button class="action-secondary sample-bulk-apply" type="button" data-group="treatment">Set Treatment</button>
        <button class="action-secondary sample-bulk-apply" type="button" data-group="">Set Ignore</button>
        <span class="sample-bulk-count" aria-live="polite"></span>
      </div>`;
    }

    function sampleFilterText(item) {
      return String(item.textContent || "").toLowerCase();
    }

    function matchingSampleItems(row) {
      const field = row.querySelector(".sample-filter");
      const needle = String(field ? field.value : "").trim().toLowerCase();
      const items = [...row.querySelectorAll(".sample-item")];
      return needle ? items.filter((item) => sampleFilterText(item).includes(needle)) : items;
    }

    function refreshSampleFilter(row) {
      const field = row.querySelector(".sample-filter");
      if (!field) return;
      const needle = field.value.trim().toLowerCase();
      const items = [...row.querySelectorAll(".sample-item")];
      let matched = 0;
      items.forEach((item) => {
        const hit = !needle || sampleFilterText(item).includes(needle);
        item.classList.toggle("is-filtered-out", Boolean(needle) && !hit);
        if (hit) matched += 1;
      });
      const counter = row.querySelector(".sample-bulk-count");
      if (counter) {
        counter.textContent = needle
          ? `${matched} of ${items.length} match`
          : `${items.length} sample${items.length === 1 ? "" : "s"}`;
      }
      // Matching is plain substring, so "induced" also matches "uninduced".
      // The visible rows are the preview; the count on the button makes the
      // press a decision about a stated number rather than a hopeful one.
      const names = { control: "Control", treatment: "Treatment", "": "Ignore" };
      row.querySelectorAll(".sample-bulk-apply").forEach((button) => {
        button.textContent = `Set ${names[button.dataset.group] || "Ignore"} (${matched})`;
        button.disabled = matched === 0;
      });
    }

    function applySampleBulk(row, group) {
      const targets = matchingSampleItems(row);
      if (!targets.length) return;
      targets.forEach((item) => {
        const select = item.querySelector("[data-sample]");
        if (select) select.value = group;
      });
      // The attestations below were ticked about the previous assignment.
      ["direction-confirmed", "biological-replicates-confirmed"].forEach((name) => {
        const box = row.querySelector(`.${name}`);
        if (box && box.checked) box.checked = false;
      });
      updateSampleCounts(row);
      capturePreparedDraft();
      updateAnalysisEligibility();
    }

    function fallbackCandidateHtml(study, candidate) {
      const inspection = candidate.inspection || {};
      const draft = activeDiscoveryState().draft[candidate.candidate_id] || {};
      const contrast = draft.contrast === undefined ? "" : draft.contrast;
      const role = candidate.role || inspection.declared_role || "";
      const matrixType = draft.matrixType || "";
      const matrixTypeControl = role === "unknown_matrix" ? `<label>Matrix type<select class="matrix-type"><option value="" ${matrixType === "" ? "selected" : ""}>Choose matrix type</option><option value="normalized_expression_matrix" ${matrixType === "normalized_expression_matrix" ? "selected" : ""}>Normalized expression</option><option value="count_matrix" ${matrixType === "count_matrix" ? "selected" : ""}>Raw counts</option></select></label>` : "";
      const normalizedScale = draft.normalizedScale || "";
      const normalizedScaleControl = role === "count_matrix" ? "" : `<label>Normalized value scale<select class="normalized-scale"><option value="" ${normalizedScale === "" ? "selected" : ""}>Choose confirmed scale</option><option value="log2" ${normalizedScale === "log2" ? "selected" : ""}>Already log2 scale</option><option value="linear" ${normalizedScale === "linear" ? "selected" : ""}>Linear, apply log2(x + 1)</option></select></label>`;
      // A GSM accession on its own tells a reader nothing about which arm it
      // belongs to, and this is the one control they must not get wrong.
      // Author matrices use the submitter's own column names, so preparation
      // resolves those to samples and stores the result on the candidate.
      const columnLabels = inspection.sample_labels || {};
      const studyLabels = study.sample_labels || {};
      const columns = inspection.sample_columns || [];
      const labelledCount = columns.filter(
        (name) => columnLabels[name] || studyLabels[String(name).toUpperCase()],
      ).length;
      const samples = columns.map((sample) => {
        const group = draft.samples?.[sample] || "";
        const meta = columnLabels[sample] || studyLabels[String(sample).toUpperCase()] || {};
        const descriptor = meta.title || meta.source || "";
        const traits = (meta.characteristics || []).join(" · ");
        const tip = [sample, descriptor, traits].filter(Boolean).join(" — ");
        // When the submitter's title is just the column name again, printing
        // both wastes the line that should carry the characteristics.
        const flat = (value) => String(value).toLowerCase().replace(/[^a-z0-9]+/g, "");
        const echoesColumn = descriptor
          && (flat(sample).startsWith(flat(descriptor)) || flat(descriptor).startsWith(flat(sample)));
        const matched = Boolean(meta.accession || descriptor || traits);
        const labelLine = descriptor && !echoesColumn
          ? `<span class="sample-label">${esc(descriptor)}</span>`
          : !matched && labelledCount
          ? `<span class="sample-label sample-label-missing">not matched to a GEO sample</span>`
          : "";
        return `<label class="sample-item" title="${esc(tip)}"><span class="sample-id">${esc(sample)}</span>${labelLine}${traits ? `<span class="sample-traits">${esc(traits)}</span>` : ""}<select data-sample="${esc(sample)}"><option value="" ${group === "" ? "selected" : ""}>Ignore</option><option value="control" ${group === "control" ? "selected" : ""}>Control</option><option value="treatment" ${group === "treatment" ? "selected" : ""}>Treatment</option></select></label>`;
      }).join("");
      return `<div class="candidate-row" data-candidate-id="${esc(candidate.candidate_id)}" data-accession="${esc(study.accession)}" data-source-unit="${esc(study.source_unit_id || study.accession)}" data-mode="fallback" data-role="${esc(role)}">
        <input class="candidate-enable" type="checkbox" aria-label="Use ${esc(candidate.name)}" ${draft.enabled ? "checked" : ""}>
        <div><span class="candidate-name">${esc(candidate.name)}</span><span class="candidate-note">Labeled fallback · ${esc(role)} · choose biological groups explicitly</span></div>
        <div class="candidate-fields">
          <label>Contrast label<input class="contrast-label" value="${esc(contrast)}" placeholder="Enter exact paper/GEO contrast" autocomplete="off"></label>
          ${matrixTypeControl}
          ${normalizedScaleControl}
          <label>Gene column<input class="gene-column" value="${esc(draft.geneColumn === undefined ? (inspection.gene_column || "") : draft.geneColumn)}"></label>
          <label class="confirm-line"><input class="direction-confirmed" type="checkbox" ${draft.direction ? "checked" : ""}> The control and treated groups above are correct, and the comparison runs treated minus control</label>
          <label class="confirm-line"><input class="biological-replicates-confirmed" type="checkbox" ${draft.biologicalReplicates ? "checked" : ""}> Each selected column is a separate biological sample, not a repeat measurement of the same one</label>
        </div>
        <div class="sample-groups" role="group" aria-label="Assign independent biological control and treatment samples"><div class="sample-counts" aria-live="polite"><span>Control <strong data-control-count>0</strong></span><span>Treatment <strong data-treatment-count>0</strong></span><span>Required: 2 + 2</span>${columns.length && !labelledCount ? `<span class="sample-labels-missing">GEO returned no matching sample labels for these columns — check the series page before assigning.</span>` : ""}</div>${columns.length > 4 ? sampleBulkHtml() : ""}${samples || `<span class="candidate-note">No numeric sample columns detected.</span>`}</div>
      </div>`;
    }

    // Mirrors discovery.candidate_preference for bundles prepared before the
    // server stamped preference_rank on each file. Lower is shown first.
    function candidatePreferenceRank(candidate) {
      if (Number.isFinite(candidate.preference_rank)) return candidate.preference_rank;
      const status = String((candidate.inspection || {}).status || "");
      const name = String(candidate.name || "");
      if (status === "ready_for_review") return 0;
      if (["requires_pvalue_mapping", "requires_lfc_confirmation", "candidate_header", "requires_column_mapping"].includes(status)) return 1;
      if (status.startsWith("upstream_matrix")) {
        if (/log2|log_2|tmm|vst|rlog|voom/i.test(name)) return 3;
        if (/fpkm|tpm|rpkm|cpm|normali[sz]ed|expression/i.test(name)) return 4;
        if (/raw|counts?|htseq|featurecounts|salmon|kallisto|rsem/i.test(name)) return 2;
        return 4;
      }
      return 9;
    }

    function renderPreparedState() {
      const state = activeDiscoveryState();
      $("downloadAnalysisExcel").disabled = !state.run?.excel_workbook;
      $("discoveryError").textContent = state.analysisError || "";
      $("discoveryError").hidden = !state.analysisError;
      if (!state.prepared) {
        if (state.preparing) {
          // Preparation issues dozens of paced network round trips. Show the same
          // determinate progress the search does instead of a dead button label.
          const percent = typeof state.prepareProgress === "number"
            ? Math.max(0, Math.min(100, Math.round(state.prepareProgress * 100)))
            : null;
          const stage = state.prepareMessage || "Preparation queued";
          const elapsed = state.prepareJobStartedAt
            ? Math.max(0, Math.round((Date.now() - state.prepareJobStartedAt) / 1000))
            : 0;
          const elapsedText = elapsed >= 1 ? ` · ${formatElapsed(elapsed)} elapsed` : "";
          const percentText = percent === null ? "" : ` · ${percent}%`;
          const bar = percent === null
            ? `<div class="loading-bar" aria-hidden="true"></div>`
            : `<div class="loading-bar is-determinate" aria-hidden="true"><span style="width: ${percent}%"></span></div>`;
          $("preparedCard").hidden = false;
          $("preparedStatus").textContent = "Preparing";
          $("preparedCandidates").innerHTML = `<div class="discovery-loading">`
            + `<div class="loading-card search-progress" aria-busy="true"`
            + `${percent === null ? "" : ` role="progressbar" aria-valuenow="${percent}" aria-valuemin="0" aria-valuemax="100"`}>`
            + `<strong class="loading-title">Downloading and inspecting the selected publications</strong>`
            + bar
            + `<span class="loading-note">${esc(stage)}${esc(percentText)}${esc(elapsedText)}</span>`
            + `<span class="loading-note">Public repositories are queried at a paced rate, so this takes longer for larger selections.</span>`
            + `<button type="button" class="job-cancel" id="cancelPrepareJob"${state.cancelling ? " disabled" : ""}>`
            + `${state.cancelling ? "Stopping..." : "Stop preparing"}</button>`
            + `</div></div>`;
          $("analysisEligibility").textContent = "Preparation in progress.";
          $("runDiscoveryAnalysis").disabled = true;
          $("analysisCompleteCard").hidden = !state.run;
          renderDiscoveryHeaderMeta();
          return;
        }
        $("preparedCard").hidden = true;
        $("analysisCompleteCard").hidden = !state.run;
        renderDiscoveryHeaderMeta();
        return;
      }
      $("preparedCard").hidden = false;
      const allStudies = state.prepared.studies || [];
      // A reader scanning for something to activate should not have to read
      // past three "no usable table" cards to find the fourth that has one.
      // Studies with at least one candidate come first, in their order; the
      // rest are grouped under one heading, collapsed, with their reasons intact.
      const analyzable = allStudies.filter((study) => (study.files || []).some(eligibleCandidate));
      const unanalyzable = allStudies.filter((study) => !(study.files || []).some(eligibleCandidate));
      const renderStudy = (study) => {
        // One file in front, the rest behind it. A series that ships the same
        // samples as raw counts, log2 TMM and FPKM is one experiment three
        // times; showing all three as peers invited selecting all three, which
        // the server now refuses. The preferred one is the least processed
        // evidence (server-ranked; ranked here too for bundles prepared before).
        const candidates = (study.files || []).filter(eligibleCandidate)
          .map((candidate, index) => ({ candidate, index, rank: candidatePreferenceRank(candidate) }))
          .sort((a, b) => a.rank - b.rank || a.index - b.index)
          .map((entry) => entry.candidate);
        const hasAuthor = candidates.some(isAuthorReviewable);
        const firstAuthor = candidates.findIndex(isAuthorReviewable);
        const renderCandidate = (candidate, index) => {
          if (!isAuthorReviewable(candidate)) return fallbackCandidateHtml(study, candidate);
          const keys = authorDraftKeys(candidate.candidate_id);
          const base = authorCandidateHtml(study, candidate, hasAuthor && index === firstAuthor, candidate.candidate_id, false);
          const clones = keys.filter((key) => key !== candidate.candidate_id).map((key) => authorCandidateHtml(study, candidate, false, key, true)).join("");
          return base + clones;
        };
        const primary = candidates.length ? renderCandidate(candidates[0], 0) : "";
        const preferredNote = candidates.length && candidates[0].preference_reason
          ? `<p class="preferred-note">Shown first: ${esc(candidates[0].preference_reason)}.</p>`
          : "";
        const others = candidates.slice(1);
        const alternatives = others.length
          ? `<details class="alternative-candidates"><summary>${others.length} other file${others.length === 1 ? "" : "s"} from this series`
            + ` - usually the same samples in another normalization; open only if the file above is not the one to use</summary>`
            + others.map((candidate, offset) => renderCandidate(candidate, offset + 1)).join("") + `</details>`
          : "";
        const rows = preferredNote + primary + alternatives;
        return `<div class="candidate-study"><h4>${esc([study.accession, study.paper_title || study.title || "Untitled study"].filter(Boolean).join(" · "))}</h4><p>${esc(study.preparation_status || "review required")}</p>${rows || unusableStudyHtml(study)}</div>`;
      };
      const html = analyzable.map(renderStudy).join("")
        + (unanalyzable.length
          ? `<details class="unanalyzable-group"><summary>${unanalyzable.length} stud${unanalyzable.length === 1 ? "y" : "ies"} with no usable table`
            + ` - open to see why each was set aside</summary>${unanalyzable.map(renderStudy).join("")}</details>`
          : "");
      const excluded = (state.prepared.excluded_studies || []).map((item) => {
        // The server sends canonical_id/paper_title/source_unit_id, never `accession`,
        // so every excluded card used to render an empty heading.
        const label = item.paper_title || item.canonical_id || item.source_unit_id || item.accession || "Unidentified publication";
        return `<div class="candidate-study"><h4>${esc(label)}</h4><p>${esc(item.reason)}</p></div>`;
      }).join("");
      const units = usableSourceUnits(state.prepared);
      const excludedCount = (state.prepared.excluded_studies || []).length;
      const blocked = units < 2;
      const blockedNotice = blocked
        ? `<div class="prepared-blocked" role="status">`
          + `<strong>This preparation cannot be analysed.</strong> `
          + `${units} of ${studies.length + excludedCount} prepared stud${studies.length + excludedCount === 1 ? "y" : "ies"} produced a usable candidate, `
          + `and DEGORA needs two independent source units. The review fields below are switched off because preparing a new selection clears them — `
          + `go back to the results, add another study, and prepare again.`
          + `<button class="action-secondary" type="button" data-back-to-results>Back to studies</button></div>`
        : "";
      $("preparedCandidates").innerHTML = blockedNotice + (html + excluded || `<div class="discovery-empty">No candidates were prepared.</div>`);
      $("preparedCandidates").classList.toggle("is-unanalyzable", blocked);
      $("preparedCandidates").querySelectorAll(".candidate-row").forEach(refreshSampleFilter);
      $("preparedCandidates")
        .querySelectorAll(".candidate-row input, .candidate-row select, .candidate-row button")
        .forEach((control) => { control.disabled = blocked; });
      $("preparedStatus").textContent = blocked
        ? `${studies.length} prepared · ${units} usable`
        : `${studies.length} studies prepared`;
      $("analysisCompleteCard").hidden = !state.run;
      if (state.run) {
        $("analysisCompleteTitle").textContent = `${speciesLabel(activeSpecies)} DEGORA analysis complete`;
        $("analysisCompleteText").textContent = `${state.run.n_source_units || state.run.source_units?.length || 0} independent source units were analyzed separately from the other species. Top genes: ${(state.run.top_genes || []).slice(0, 8).join(", ")}.`;
      }
      updateAnalysisEligibility();
      renderDiscoveryHeaderMeta();
    }

    async function pollPrepareJob(species, requestId, jobId) {
      const state = discoveryStates[species];
      state.prepareJobStartedAt = Date.now();
      while (requestId === state.prepareRequest) {
        const payload = await getJson(`/api/discovery/jobs/${jobId}`);
        if (requestId !== state.prepareRequest) return null;
        const job = payload.job || payload;
        state.prepareProgress = typeof job.progress === "number" ? job.progress : null;
        state.prepareMessage = typeof job.message === "string" ? job.message : "";
        if (job.status === "complete") return job.result || {};
        if (job.status === "cancelled") {
          state.prepareCancelled = true;
          return null;
        }
        if (job.status === "failed" || job.status === "interrupted") {
          throw new Error(job.error || "preparation failed");
        }
        if (activeSpecies === species) renderPreparedState();
        await new Promise((resolve) => setTimeout(resolve, 650));
      }
      return null;
    }

    // Prepare takes tens of seconds and its progress card sits below a full
    // page of results, so pressing the button looked like nothing happened.
    // Bring the card into view when the work starts, not only when it ends.
    function revealPreparedCard() {
      const card = $("preparedCard");
      if (!card || card.hidden) return;
      const box = card.getBoundingClientRect();
      const viewport = window.innerHeight || document.documentElement.clientHeight || 0;
      // Already comfortably on screen: scrolling again would only jitter the
      // view under a reader who is mid-sentence.
      if (box.top >= 0 && box.top < viewport * 0.75) return;
      const reduced = window.matchMedia
        && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
      card.scrollIntoView({ behavior: reduced ? "auto" : "smooth", block: "start" });
    }

    async function prepareSelectedStudies({ recordIds: explicitIds = null } = {}) {
      const requestSpecies = activeSpecies;
      const state = discoveryStates[requestSpecies];
      // Inspect prepares a single publication without touching the selection:
      // it used to overwrite it, silently discarding a full 20-item selection
      // along with any bundle and review work already done.
      const recordIds = explicitIds ? [...explicitIds] : [...state.selected];
      if (!recordIds.length) return;
      const requestId = ++state.prepareRequest;
      const query = state.query;
      state.preparing = true;
      state.prepareProgress = null;
      state.prepareMessage = "";
      state.prepareJobStartedAt = Date.now();
      state.prepared = null;
      state.bundleId = "";
      state.run = null;
      state.draft = {};
      state.cloneCounter = 0;
      state.analysisError = "";
      state.notice = "";
      state.noticeLevel = "info";
      let preparationFailed = false;
      if (activeSpecies === requestSpecies) {
        activeRunId = "";
        invalidateAtlasContext();
        renderPreparedState();
        updateSelectedStatus();
        revealPreparedCard();
      }
      try {
        const started = await postJson("/api/discovery/prepare-jobs", {
          species: requestSpecies,
          query,
          search_id: state.searchId,
          record_ids: recordIds
        });
        if (requestId !== state.prepareRequest) return;
        state.prepareJobId = started.job_id;
        const data = await pollPrepareJob(requestSpecies, requestId, started.job_id);
        if (requestId !== state.prepareRequest || data === null) return;
        state.prepared = data;
        state.bundleId = data.bundle_id;
        state.run = null;
        state.draft = {};
        state.cloneCounter = 0;
        if (activeSpecies === requestSpecies) {
          renderPreparedState();
          // Usually a no-op now: the card was brought into view when the work
          // started, so a reader who scrolled away on purpose is left alone.
          revealPreparedCard();
        }
      } catch (error) {
        if (requestId !== state.prepareRequest) return;
        state.notice = `Preparation failed: ${error.message}`;
        state.noticeLevel = "error";
        preparationFailed = true;
        if (activeSpecies === requestSpecies) renderPreparedState();
      } finally {
        if (requestId !== state.prepareRequest) return;
        state.preparing = false;
        if (activeSpecies === requestSpecies) {
          renderDiscoveryResults();
          renderPreparedState();
          updateSelectedStatus();
          if (preparationFailed) $("discoveryNotice").focus();
        }
      }
    }

    function selectedCandidateRows() {
      return [...$("preparedCandidates").querySelectorAll(".candidate-row")].filter((row) => row.querySelector(".candidate-enable")?.checked);
    }

    function updateSampleCounts(row) {
      if (row.dataset.mode !== "fallback") return;
      const values = [...row.querySelectorAll("[data-sample]")].map((select) => select.value);
      const control = values.filter((value) => value === "control").length;
      const treatment = values.filter((value) => value === "treatment").length;
      const controlNode = row.querySelector("[data-control-count]");
      const treatmentNode = row.querySelector("[data-treatment-count]");
      if (controlNode) controlNode.textContent = String(control);
      if (treatmentNode) treatmentNode.textContent = String(treatment);
    }

    // Two ticked units that belong to one submission are not two studies. The
    // gate counts units, so without this the reader clears the ">= 2 independent
    // source units" bar with correlated arms of a single experiment.
    function sharedSubmissionConflicts(unitIds) {
      const prepared = activeDiscoveryState().prepared;
      if (!prepared) return [];
      const selected = new Set([...unitIds].filter(Boolean).map(String));
      const conflicts = [];
      (prepared.studies || []).forEach((study) => {
        const unit = String(study.source_unit_id || study.accession || "");
        if (!selected.has(unit)) return;
        const peers = (study.shared_submission_units || []).filter((peer) => selected.has(String(peer)));
        if (peers.length) conflicts.push({ unit, peers: peers.map(String) });
      });
      return conflicts;
    }

    function renderIndependenceWarning(unitIds) {
      const holder = $("preparedCandidates");
      const existing = holder.querySelector(".independence-warning");
      if (existing) existing.remove();
      const conflicts = sharedSubmissionConflicts(unitIds);
      if (!conflicts.length) return;
      const named = conflicts
        .map((item) => [item.unit, ...item.peers].sort().join(" + "))
        .filter((value, index, all) => all.indexOf(value) === index);
      const node = document.createElement("div");
      node.className = "independence-warning";
      node.setAttribute("role", "status");
      node.innerHTML = `<strong>These may not be independent studies.</strong> `
        + `${esc(named.join("; "))} share a title and none is linked to a publication, `
        + `so they are probably one submission deposited as several series. DEGORA counts them as `
        + `separate source units, which would treat correlated data as replication. `
        + `Untick all but one unless you have confirmed they are separate experiments.`;
      holder.insertBefore(node, holder.firstChild);
    }

    // A record found only through the literature search carries the organism filter
    // that produced the search and no per-record evidence, so DEGORA labels it
    // query_constrained. The README asks the reader to confirm the species of such
    // a record; nothing in the browser ever collected that answer or recorded it.
    function unconfirmedSpeciesStudies(state) {
      const studies = (state.prepared && state.prepared.studies) || [];
      return studies.filter((study) => {
        const decision = String(study.species_decision || "").toLowerCase();
        return decision === "query_constrained" || decision === "unknown" || decision === "";
      });
    }

    function updateSpeciesConfirmation(state) {
      const pending = unconfirmedSpeciesStudies(state);
      const line = $("speciesConfirmLine");
      const control = $("speciesConfirmed");
      if (!pending.length) {
        line.hidden = true;
        control.checked = false;
        return true;
      }
      line.hidden = false;
      $("speciesConfirmText").textContent =
        `${pending.length} of these record${pending.length === 1 ? " was" : "s were"} matched by the `
        + `${speciesLabel(activeSpecies)} search filter, not by a per-record organism check. `
        + `I have confirmed ${pending.length === 1 ? "it is" : "they are"} ${speciesLabel(activeSpecies)} data.`;
      return Boolean(control.checked);
    }

    function updateAnalysisEligibility() {
      const rows = selectedCandidateRows();
      const speciesConfirmed = updateSpeciesConfirmation(activeDiscoveryState());
      $("preparedCandidates").querySelectorAll(".candidate-row").forEach(updateSampleCounts);
      // A switched-off form must not also be painted red for being incomplete.
      const blockedState = activeDiscoveryState();
      if (blockedState.prepared && usableSourceUnits(blockedState.prepared) < 2) {
        $("preparedCandidates").querySelectorAll("[aria-invalid]").forEach((node) => {
          node.removeAttribute("aria-invalid");
        });
        $("runDiscoveryAnalysis").disabled = true;
        $("runDiscoveryAnalysis").textContent = "Run species-specific DEGORA";
        $("analysisEligibility").textContent =
          `Not analysable: ${usableSourceUnits(blockedState.prepared)} usable stud`
          + `${usableSourceUnits(blockedState.prepared) === 1 ? "y" : "ies"} in this preparation, and DEGORA needs two `
          + `independent ${speciesLabel(activeSpecies)} source units. Add another study and prepare again.`;
        return;
      }
      const units = new Set(rows.map((row) => row.dataset.sourceUnit || row.dataset.accession));
      const hasFallback = rows.some((row) => row.dataset.mode === "fallback");
      const textValue = (row, selector) => (row.querySelector(selector)?.value || "").trim();
      const boolValue = (row, selector) => Boolean(row.querySelector(selector)?.checked);
      const authorMappingEdited = (row) => (
        textValue(row, ".sheet-name") !== (row.dataset.detectedSheetName || "")
        || textValue(row, ".gene-column") !== (row.dataset.detectedGeneColumn || "")
        || textValue(row, ".lfc-column") !== (row.dataset.detectedLfcColumn || "")
        || textValue(row, ".p-column") !== (row.dataset.detectedPColumn || "")
        || textValue(row, ".padj-column") !== (row.dataset.detectedPadjColumn || "")
      );
      const authorRowFilterComplete = (row) => {
        const hasColumn = Boolean(textValue(row, ".row-filter-column"));
        const hasValue = Boolean(textValue(row, ".row-filter-value"));
        return hasColumn === hasValue && (!hasColumn || boolValue(row, ".row-filter-confirmed"));
      };
      const authorDuplicatePolicyComplete = (row) => (
        textValue(row, ".duplicate-gene-policy") !== "keep_first"
        || boolValue(row, ".duplicate-gene-policy-confirmed")
      );
      const authorReviewComplete = (row) => {
        const status = row.dataset.authorStatus || "";
        const pColumn = textValue(row, ".p-column");
        const padjColumn = textValue(row, ".padj-column");
        const needsMapping = status !== "ready_for_review" || authorMappingEdited(row);
        return Boolean(textValue(row, ".gene-column"))
          && Boolean(textValue(row, ".lfc-column"))
          && Boolean(pColumn)
          && (!needsMapping || boolValue(row, ".column-mapping-confirmed"))
          && (!(pColumn && padjColumn && pColumn === padjColumn) || boolValue(row, ".adjusted-p-as-pvalue-confirmed"))
          && (status !== "requires_lfc_confirmation" || boolValue(row, ".lfc-scale-confirmed-log2"))
          && authorRowFilterComplete(row)
          && authorDuplicatePolicyComplete(row);
      };
      const reviewComplete = rows.length > 0 && speciesConfirmed && rows.every((row) => {
        if (!(row.querySelector(".contrast-label")?.value || "").trim()) return false;
        if (!row.querySelector(".direction-confirmed")?.checked) return false;
        if (row.querySelector(".lfc-scale-confirmed") && !row.querySelector(".lfc-scale-confirmed").checked) return false;
        if (row.dataset.mode !== "fallback") {
          return isPositiveWholeNumber(row.querySelector(".n-ctrl")?.value)
            && isPositiveWholeNumber(row.querySelector(".n-treat")?.value)
            && authorReviewComplete(row);
        }
        if (!(row.querySelector(".gene-column")?.value || "").trim()) return false;
        if (!row.querySelector(".biological-replicates-confirmed")?.checked) return false;
        if (row.dataset.role === "unknown_matrix" && !row.querySelector(".matrix-type")?.value) return false;
        const resolvedRole = row.dataset.role === "unknown_matrix" ? row.querySelector(".matrix-type")?.value : row.dataset.role;
        if (resolvedRole === "normalized_expression_matrix" && !row.querySelector(".normalized-scale")?.value) return false;
        const groups = [...row.querySelectorAll("[data-sample]")].map((select) => select.value);
        return groups.filter((value) => value === "control").length >= 2
          && groups.filter((value) => value === "treatment").length >= 2;
      });
      const markInvalid = (control, invalid) => {
        if (!control) return;
        if (invalid) control.setAttribute("aria-invalid", "true");
        else control.removeAttribute("aria-invalid");
      };
      // The activation gate already asks for each confirmation only where it
      // applies. The panel asked for all of them on every row, so a table needing
      // one assertion presented six, and the one that mattered read like the rest.
      const showWhenRequired = (control, required) => {
        const line = control?.closest(".confirm-line");
        if (!line) return;
        line.classList.toggle("not-required", !required);
        if (!required && control.checked) control.checked = false;
      };
      $("preparedCandidates").querySelectorAll(".candidate-row").forEach((row) => {
        const enabled = Boolean(row.querySelector(".candidate-enable")?.checked);
        const contrast = row.querySelector(".contrast-label");
        const direction = row.querySelector(".direction-confirmed");
        const contrastValid = Boolean((contrast?.value || "").trim());
        const directionValid = Boolean(direction?.checked);
        markInvalid(contrast, enabled && !contrastValid);
        markInvalid(direction, enabled && !directionValid);
        let rowValid = contrastValid && directionValid;
        if (row.dataset.mode === "author") {
          const nCtrl = row.querySelector(".n-ctrl");
          const nTreat = row.querySelector(".n-treat");
          const seriesTotal = Number(row.dataset.seriesSamples || 0);
          const enteredTotal = Number(nCtrl?.value || 0) + Number(nTreat?.value || 0);
          // Entering the series total in both boxes is the common mistake, and it
          // doubles the source's weight against every other study in the corpus.
          const fitsSeries = !seriesTotal || enteredTotal <= seriesTotal;
          const ctrlValid = isPositiveWholeNumber(nCtrl?.value) && fitsSeries;
          const treatValid = isPositiveWholeNumber(nTreat?.value) && fitsSeries;
          const gene = row.querySelector(".gene-column");
          const lfc = row.querySelector(".lfc-column");
          const pColumn = row.querySelector(".p-column");
          const padjColumn = row.querySelector(".padj-column");
          const mappingConfirmed = row.querySelector(".column-mapping-confirmed");
          const adjustedPConfirmed = row.querySelector(".adjusted-p-as-pvalue-confirmed");
          const lfcScaleConfirmed = row.querySelector(".lfc-scale-confirmed-log2");
          const rowFilterColumn = row.querySelector(".row-filter-column");
          const rowFilterValue = row.querySelector(".row-filter-value");
          const rowFilterConfirmed = row.querySelector(".row-filter-confirmed");
          const duplicateGenePolicy = row.querySelector(".duplicate-gene-policy");
          const duplicateGenePolicyConfirmed = row.querySelector(".duplicate-gene-policy-confirmed");
          const needsMapping = (row.dataset.authorStatus || "") !== "ready_for_review" || authorMappingEdited(row);
          const pValue = textValue(row, ".p-column");
          const padjValue = textValue(row, ".padj-column");
          const adjustedPValid = !(pValue && padjValue && pValue === padjValue) || Boolean(adjustedPConfirmed?.checked);
          const lfcScaleValid = row.dataset.authorStatus !== "requires_lfc_confirmation" || Boolean(lfcScaleConfirmed?.checked);
          const rowFilterPairValid = Boolean(textValue(row, ".row-filter-column")) === Boolean(textValue(row, ".row-filter-value"));
          const rowFilterValid = rowFilterPairValid && (!textValue(row, ".row-filter-column") || Boolean(rowFilterConfirmed?.checked));
          const duplicatePolicyValid = duplicateGenePolicy?.value !== "keep_first" || Boolean(duplicateGenePolicyConfirmed?.checked);
          markInvalid(nCtrl, enabled && !ctrlValid);
          markInvalid(nTreat, enabled && !treatValid);
          markInvalid(gene, enabled && !textValue(row, ".gene-column"));
          markInvalid(lfc, enabled && !textValue(row, ".lfc-column"));
          markInvalid(pColumn, enabled && !pValue);
          showWhenRequired(mappingConfirmed, needsMapping);
          showWhenRequired(adjustedPConfirmed, Boolean(pValue && padjValue && pValue === padjValue));
          showWhenRequired(lfcScaleConfirmed, row.dataset.authorStatus === "requires_lfc_confirmation");
          showWhenRequired(rowFilterConfirmed, Boolean(textValue(row, ".row-filter-column")));
          showWhenRequired(duplicateGenePolicyConfirmed, duplicateGenePolicy?.value === "keep_first");
          markInvalid(mappingConfirmed, enabled && needsMapping && !mappingConfirmed?.checked);
          markInvalid(adjustedPConfirmed, enabled && !adjustedPValid);
          markInvalid(lfcScaleConfirmed, enabled && !lfcScaleValid);
          markInvalid(rowFilterColumn, enabled && !rowFilterPairValid);
          markInvalid(rowFilterValue, enabled && !rowFilterPairValid);
          markInvalid(rowFilterConfirmed, enabled && !rowFilterValid);
          markInvalid(duplicateGenePolicyConfirmed, enabled && !duplicatePolicyValid);
          rowValid = rowValid && ctrlValid && treatValid && authorReviewComplete(row);
        } else {
          const gene = row.querySelector(".gene-column");
          const biological = row.querySelector(".biological-replicates-confirmed");
          const matrixType = row.querySelector(".matrix-type");
          const normalizedScale = row.querySelector(".normalized-scale");
          const geneValid = Boolean((gene?.value || "").trim());
          const biologicalValid = Boolean(biological?.checked);
          const matrixTypeValid = row.dataset.role !== "unknown_matrix" || Boolean(matrixType?.value);
          const resolvedRole = row.dataset.role === "unknown_matrix" ? matrixType?.value : row.dataset.role;
          const scaleValid = resolvedRole !== "normalized_expression_matrix" || Boolean(normalizedScale?.value);
          const groups = [...row.querySelectorAll("[data-sample]")].map((select) => select.value);
          const groupsValid = groups.filter((value) => value === "control").length >= 2
            && groups.filter((value) => value === "treatment").length >= 2;
          markInvalid(gene, enabled && !geneValid);
          markInvalid(biological, enabled && !biologicalValid);
          markInvalid(matrixType, enabled && !matrixTypeValid);
          markInvalid(normalizedScale, enabled && !scaleValid);
          markInvalid(row.querySelector(".sample-groups"), enabled && !groupsValid);
          rowValid = rowValid && geneValid && biologicalValid && matrixTypeValid && scaleValid && groupsValid;
        }
        markInvalid(row, enabled && !rowValid);
      });
      renderIndependenceWarning(units);
      const eligible = units.size >= 2 && reviewComplete;
      const state = activeDiscoveryState();
      $("runDiscoveryAnalysis").disabled = !eligible || state.analyzing;
      $("runDiscoveryAnalysis").textContent = state.analyzing ? "Running DEGORA..." : "Run species-specific DEGORA";
      if (units.size < 2) {
        // "Select candidates from at least two independent studies" is true and
        // useless when only one of the prepared studies has a candidate to
        // select. Say how much supply there is.
        const prepared = state.prepared || {};
        const preparedStudies = prepared.studies || [];
        const usable = preparedStudies.filter((study) => (study.files || []).some(eligibleCandidate)).length;
        const total = preparedStudies.length + (prepared.excluded_studies || []).length;
        const supply = total
          ? ` ${usable} of ${total} prepared stud${total === 1 ? "y" : "ies"} produced a usable candidate.`
          : "";
        $("analysisEligibility").textContent = units.size === 1
          ? `One independent ${speciesLabel(activeSpecies)} study is selected; DEGORA needs two.${supply}`
          : `Select candidates from at least two independent ${speciesLabel(activeSpecies)} studies.${supply}`;
      } else if (!reviewComplete) {
        $("analysisEligibility").textContent = hasFallback
          ? "Complete each exact contrast and direction confirmation; fallback matrices also require scale, biological-replicate attestation, and 2 + 2 sample assignment."
          : "Enter exact contrasts, mappings, positive whole-number biological group sizes, and all required author-table confirmations.";
      } else {
        $("analysisEligibility").textContent = `${rows.length} candidate${rows.length === 1 ? "" : "s"} from ${units.size} independent ${speciesLabel(activeSpecies)} studies; review complete.`;
      }
    }

    function collectAnalysisSelections() {
      return selectedCandidateRows().map((row) => {
        const common = {
          candidate_id: row.dataset.candidateId,
          mode: row.dataset.mode,
          contrast_label: row.querySelector(".contrast-label")?.value.trim() || "",
          direction_confirmed: Boolean(row.querySelector(".direction-confirmed")?.checked)
        };
        if (row.dataset.mode === "author") {
          common.table_scope = row.querySelector(".table-scope")?.value || "auto";
          common.cell_system = row.querySelector(".cell-system")?.value.trim() || "";
          common.duration_h = row.querySelector(".duration-h")?.value.trim() || "";
          common.n_ctrl = Number(row.querySelector(".n-ctrl")?.value || 0);
          common.n_treat = Number(row.querySelector(".n-treat")?.value || 0);
          common.platform = row.querySelector(".platform")?.value.trim() || "";
          common.assay_type = row.querySelector(".assay-type")?.value.trim() || "";
          common.pipeline = row.querySelector(".pipeline")?.value.trim() || "";
          common.sheet_name = row.querySelector(".sheet-name")?.value.trim() || "";
          common.gene_column = row.querySelector(".gene-column")?.value.trim() || "";
          common.lfc_column = row.querySelector(".lfc-column")?.value.trim() || "";
          common.p_column = row.querySelector(".p-column")?.value.trim() || "";
          common.padj_column = row.querySelector(".padj-column")?.value.trim() || "";
          common.column_mapping_confirmed = Boolean(row.querySelector(".column-mapping-confirmed")?.checked);
          common.adjusted_p_as_pvalue_confirmed = Boolean(row.querySelector(".adjusted-p-as-pvalue-confirmed")?.checked);
          common.lfc_scale_confirmed_log2 = Boolean(row.querySelector(".lfc-scale-confirmed-log2")?.checked);
          common.row_filter_column = row.querySelector(".row-filter-column")?.value.trim() || "";
          common.row_filter_value = row.querySelector(".row-filter-value")?.value.trim() || "";
          common.row_filter_confirmed = Boolean(row.querySelector(".row-filter-confirmed")?.checked);
          common.duplicate_gene_policy = row.querySelector(".duplicate-gene-policy")?.value || "harmonizer";
          common.duplicate_gene_policy_confirmed = Boolean(row.querySelector(".duplicate-gene-policy-confirmed")?.checked);
        } else {
          common.gene_column = row.querySelector(".gene-column")?.value.trim() || "";
          common.normalized_scale = row.querySelector(".normalized-scale")?.value || "";
          common.biological_replicates_confirmed = Boolean(row.querySelector(".biological-replicates-confirmed")?.checked);
          common.control_samples = [];
          common.treatment_samples = [];
          row.querySelectorAll("[data-sample]").forEach((select) => {
            if (select.value === "control") common.control_samples.push(select.dataset.sample);
            if (select.value === "treatment") common.treatment_samples.push(select.dataset.sample);
          });
          if (row.dataset.role === "unknown_matrix") common.matrix_type = row.querySelector(".matrix-type")?.value || "";
        }
        return common;
      });
    }

    async function runSelectedAnalysis() {
      const requestSpecies = activeSpecies;
      const state = discoveryStates[requestSpecies];
      const selections = collectAnalysisSelections();
      const requestId = ++state.analysisRequest;
      state.analyzing = true;
      state.analysisError = "";
      state.run = null;
      if (activeSpecies === requestSpecies) {
        activeRunId = "";
        invalidateAtlasContext();
        renderPreparedState();
      }
      try {
        const pendingSpecies = unconfirmedSpeciesStudies(state);
        const result = await postJson("/api/discovery/analyze", {
          bundle_id: state.bundleId,
          species: requestSpecies,
          selections,
          species_confirmed: pendingSpecies.length === 0 || Boolean($("speciesConfirmed").checked),
          species_confirmation_required_for: pendingSpecies.length
        });
        if (requestId !== state.analysisRequest) return;
        state.run = result;
        if (activeSpecies === requestSpecies) {
          activeRunId = result.run_id;
          invalidateAtlasContext();
          // invalidateAtlasContext() empties the atlas and waits for the next
          // visit to reload it. If the atlas is the visible view, nothing else
          // triggers that reload and it stays blank until the user bounces tabs.
          if (!isDiscoveryView()) void ensureAtlasContext();
          renderDiscoveryHeaderMeta();
          renderPreparedState();
          $("analysisCompleteCard").scrollIntoView({ behavior: "smooth", block: "start" });
        }
      } catch (error) {
        if (requestId !== state.analysisRequest) return;
        state.analysisError = error.message;
      } finally {
        if (requestId !== state.analysisRequest) return;
        state.analyzing = false;
        if (activeSpecies === requestSpecies) renderPreparedState();
      }
    }

    function openDiscoveryAnalysis() {
      const state = activeDiscoveryState();
      if (!state.run) return;
      activeRunId = state.run.run_id;
      invalidateAtlasContext();
      showView("atlas");
    }

    async function downloadAnalysisExcel() {
      const state = activeDiscoveryState();
      if (!state.run?.run_id || !state.run?.excel_workbook) return;
      const button = $("downloadAnalysisExcel");
      button.disabled = true;
      button.textContent = "Preparing Excel...";
      try {
        const headers = {};
        if (API_TOKEN) headers["X-DEGORA-Token"] = API_TOKEN;
        const options = Object.keys(headers).length ? { headers } : {};
        const response = await fetch(`/api/discovery/runs/${state.run.run_id}/export.xlsx`, options);
        if (!response.ok) {
          let message = await response.text();
          try {
            const payload = JSON.parse(message);
            if (payload && payload.error) message = payload.error;
          } catch (_) {}
          throw new Error(message);
        }
        const blob = await response.blob();
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = `DEGORA_${state.run.species?.key || activeSpecies}_output.xlsx`;
        link.hidden = true;
        document.body.appendChild(link);
        link.click();
        link.remove();
        window.setTimeout(() => URL.revokeObjectURL(url), 1000);
      } catch (error) {
        state.analysisError = `Excel download failed: ${error.message}`;
        renderPreparedState();
      } finally {
        button.textContent = "Download Excel";
        button.disabled = !state.run?.excel_workbook;
      }
    }

    async function downloadSearchExcel() {
      const state = activeDiscoveryState();
      if (!state.searchId || !state.verified) return;
      const button = $("downloadSearchExcel");
      button.disabled = true;
      button.textContent = "Preparing search Excel...";
      try {
        const headers = {};
        if (API_TOKEN) headers["X-DEGORA-Token"] = API_TOKEN;
        const options = Object.keys(headers).length ? { headers } : {};
        const response = await fetch(`/api/discovery/searches/${state.searchId}/export.xlsx`, options);
        if (!response.ok) {
          let message = await response.text();
          try {
            const payload = JSON.parse(message);
            if (payload && payload.error) message = payload.error;
          } catch (_) {}
          throw new Error(message);
        }
        const blob = await response.blob();
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = `DEGORA_${activeSpecies}_${state.searchId}_search.xlsx`;
        link.hidden = true;
        document.body.appendChild(link);
        link.click();
        link.remove();
        window.setTimeout(() => URL.revokeObjectURL(url), 1000);
        state.notice = "Search workbook downloaded with publications, identifiers, linked datasets, candidate routes, species decisions, and provider events.";
      } catch (error) {
        state.notice = `Search Excel download failed: ${error.message}`;
      } finally {
        button.textContent = "Download search Excel";
        button.disabled = !state.searchId || !state.verified;
        renderDiscoveryNotice(state);
      }
    }

    async function loadMeta(context, generation) {
      if (context.kind === "run") {
        const run = await getJson(`/api/discovery/runs/${context.runId}/summary`);
        if (!atlasContextIsCurrent(context, generation)) return false;
        $("meta").innerHTML = [
          `${speciesLabel(context.species)} discovery run`,
          `${run.n_source_units || run.source_units?.length || 0} source units`,
          "species-specific"
        ].map((text) => `<span>${esc(text)}</span>`).join("");
        return true;
      }
      const health = await getJson("/api/health");
      if (!atlasContextIsCurrent(context, generation)) return false;
      $("meta").innerHTML = [
        health.db_name,
        `DEGORA ${health.database_degora_version || health.degora_version}`,
        `${health.gene_count.toLocaleString()} genes`,
        `${health.study_count.toLocaleString()} studies`,
        `${health.source_unit_count.toLocaleString()} source units`
      ].map((text) => `<span>${esc(text)}</span>`).join("");
      return true;
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
    let atlasContextGeneration = 0;
    let atlasLoadedContextKey = "";
    let geneListRequestSerial = 0;
    let geneDetailRequestSerial = 0;
    let currentGeneDetail = null;

    function atlasContextIsCurrent(context, generation) {
      return generation === atlasContextGeneration && context.key === currentAtlasContext().key;
    }

    function invalidateAtlasContext() {
      atlasContextGeneration += 1;
      atlasLoadedContextKey = "";
      geneListRequestSerial += 1;
      geneDetailRequestSerial += 1;
      page = { query: null, loaded: 0, total: 0 };
      currentGeneDetail = null;
      if (genesLoading) setGeneLoading(false);
      $("genes").innerHTML = "";
      $("loadMore").hidden = true;
      $("exportGenes").disabled = true;
      $("exportEvidence").disabled = true;
      $("status").textContent = "";
      $("detail").innerHTML = `<div class="empty">Atlas context changed. Loading the matching evidence...</div>`;
    }

    async function ensureAtlasContext() {
      const context = currentAtlasContext();
      if (atlasLoadedContextKey === context.key) return;
      const generation = ++atlasContextGeneration;
      geneListRequestSerial += 1;
      geneDetailRequestSerial += 1;
      page = { query: null, loaded: 0, total: 0 };
      currentGeneDetail = null;
      // Clear the filters with the context, not only the rows. A gene search left
      // over from the previous run silently applied to the new one, so an analysis
      // of 11,886 genes opened showing nine and read as an analysis that found nine.
      $("query").value = "";
      $("minUnits").value = "1";
      $("direction").value = "";
      $("genes").innerHTML = "";
      $("loadMore").hidden = true;
      $("exportGenes").disabled = true;
      $("exportEvidence").disabled = true;
      if (context.kind === "none") {
        $("meta").innerHTML = [
          `${speciesLabel(context.species)} workspace`,
          "No analysis run",
          "species kept separate"
        ].map((text) => `<span>${esc(text)}</span>`).join("");
        $("status").textContent = `No ${speciesLabel(context.species)} discovery run`;
        $("detail").innerHTML = `<div class="empty">Run a ${esc(speciesLabel(context.species))} discovery analysis before opening its evidence atlas. Human results are not shown in the Mouse workspace.</div>`;
        atlasLoadedContextKey = context.key;
        return;
      }
      $("status").textContent = "Loading matching evidence...";
      $("detail").innerHTML = `<div class="empty">Loading ${esc(speciesLabel(context.species))} evidence...</div>`;
      try {
        const metaReady = await loadMeta(context, generation);
        if (!metaReady || !atlasContextIsCurrent(context, generation)) return;
        const loaded = await loadGenes(
          context.kind === "run" ? "Loading discovery ranking..." : "Loading genes...",
          context.kind === "run" ? `${speciesLabel(context.species)} sources only` : "Refreshing the submitted evidence atlas",
          context,
          generation,
        );
        if (loaded && atlasContextIsCurrent(context, generation)) atlasLoadedContextKey = context.key;
      } catch (error) {
        if (!atlasContextIsCurrent(context, generation)) return;
        $("status").textContent = "error";
        $("detail").innerHTML = `<div class="empty">Could not load this atlas context: ${esc(error.message)}</div>`;
      }
    }

    function geneRowHtml(gene) {
      return `
        <tr data-gene="${esc(gene.gene_symbol)}" tabindex="0" aria-label="Show evidence for ${esc(gene.gene_symbol)}">
          <td class="num">${esc(primaryRank(gene))}</td>
          <td>${tier(gene.evidence_tier)}</td>
          <td class="gene">${esc(gene.gene_symbol)}</td>
          <td class="num">${esc(topPercentTableLabel(gene))}</td>
          <td class="num">${fmt(primaryScore(gene), 2)}</td>
          <td class="num">${esc(gene.n_source_units)}</td>
          <td class="num">${fmt(Number(primaryConcordance(gene)) * 100, 1)}%</td>
          <td class="num">${fmt(gene.weighted_lfc, 2)}</td>
        </tr>`;
    }

    function currentQuery() {
      const params = new URLSearchParams();
      const q = $("query").value.trim().slice(0, 128);
      const rawMinUnits = $("minUnits").value.trim();
      const direction = $("direction").value;
      if (q) params.set("q", q);
      // The server accepts an integer in [1, 10000]; sending anything else
      // returned a raw 400 that wiped the table. Clamp instead.
      if (rawMinUnits) {
        const parsed = Math.trunc(Number(rawMinUnits));
        if (Number.isFinite(parsed)) {
          const clamped = String(Math.min(10000, Math.max(1, parsed)));
          params.set("min_units", clamped);
          // Reflect the clamp so the field agrees with the results on screen.
          if (clamped !== rawMinUnits) $("minUnits").value = clamped;
        } else {
          $("minUnits").value = "1";
        }
      }
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
        if (indicator) indicator.textContent = active ? (sortState.order === "asc" ? "\u25B4" : "\u25BE") : "";
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
      $("exportGenes").disabled = isLoading || page.total === 0 || currentAtlasContext().kind === "none";
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
      $("exportGenes").disabled = page.total === 0 || genesLoading || currentAtlasContext().kind === "none";
    }

    async function fetchGenePage(
      append,
      loadingTitle = "Loading genes...",
      loadingNote = "Refreshing the table",
      context = currentAtlasContext(),
      generation = atlasContextGeneration,
    ) {
      if (context.kind === "none") return false;
      if (append && (genesLoading || page.loaded >= page.total)) return false;
      const requestId = append ? geneListRequestSerial : ++geneListRequestSerial;
      const title = append ? "Loading more genes..." : loadingTitle;
      const note = append ? "Appending the next page" : loadingNote;
      setGeneLoading(true, title, note);
      $("status").textContent = title;
      const params = new URLSearchParams(page.query);
      params.set("limit", String(PAGE_SIZE));
      params.set("offset", String(append ? page.loaded : 0));
      try {
        const data = await getJson(`${atlasApi("/genes", context)}?${params.toString()}`);
        if (requestId !== geneListRequestSerial || !atlasContextIsCurrent(context, generation)) return false;
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
          if (page.loaded) void loadGene(data.genes[0].gene_symbol, context, generation);
          else $("detail").innerHTML = `<div class="empty">No genes match these filters. Try lowering <b>Min units</b> or clearing the gene search.</div>`;
        }
        return true;
      } catch (error) {
        if (requestId !== geneListRequestSerial || !atlasContextIsCurrent(context, generation)) return false;
        $("status").textContent = "error";
        if (!append) {
          $("genes").innerHTML = "";
          $("loadMore").hidden = true;
          $("detail").innerHTML = `<div class="empty">Could not load genes: ${esc(error.message)}</div>`;
        }
        return false;
      } finally {
        if (requestId === geneListRequestSerial && atlasContextIsCurrent(context, generation)) setGeneLoading(false);
      }
    }

    async function loadGenes(
      loadingTitle = "Loading genes...",
      loadingNote = "Refreshing the table",
      context = currentAtlasContext(),
      generation = atlasContextGeneration,
    ) {
      if (context.kind === "none") return false;
      currentGeneDetail = null;
      $("exportEvidence").disabled = true;
      page = { query: currentQuery(), loaded: 0, total: 0 };
      return fetchGenePage(false, loadingTitle, loadingNote, context, generation);
    }

    async function loadGene(symbol, context = currentAtlasContext(), generation = atlasContextGeneration) {
      if (context.kind === "none") return;
      const requestId = ++geneDetailRequestSerial;
      markSelected(symbol);
      currentGeneDetail = null;
      $("exportEvidence").disabled = true;
      let data;
      try {
        data = await getJson(`${atlasApi("/genes", context)}/${encodeURIComponent(symbol)}`);
      } catch (error) {
        if (requestId !== geneDetailRequestSerial || !atlasContextIsCurrent(context, generation)) return;
        $("detail").innerHTML = `<div class="empty">Could not load ${esc(symbol)}: ${esc(error.message)}</div>`;
        return;
      }
      if (requestId !== geneDetailRequestSerial || !atlasContextIsCurrent(context, generation)) return;
      currentGeneDetail = data;
      $("exportEvidence").disabled = false;
      const gene = data.gene;
      const evidenceRows = data.evidence.map((row) => {
        const kind = evidenceKind(row);
        const sampleSize = `${firstPresent(row.n_ctrl, "?")} / ${firstPresent(row.n_treat, "?")}`;
        const sourceLabel = row.source_url ? "Open source" : "Unavailable";
        return `<tr>
          <td>${esc(row.contributing_study_ids || row.study_id)}</td>
          <td>${esc(row.source_unit_id)}</td>
          <td><span class="evidence-kind ${kind.className}">${esc(kind.label)}</span></td>
          <td>${esc(row.table_scope || "unavailable")}</td>
          <td>${esc(row.contributing_pipelines || row.pipeline || "")}</td>
          <td>${esc(row.contributing_assay_types || row.assay_type || "")}</td>
          <td class="num">${esc(sampleSize)}</td>
          <td class="num">${fmt(row.source_reliability_weight, 2)}</td>
          <td class="num">${fmt(row.lfc, 2)}</td>
          <td class="num">${fmt(row.min_source_pvalue, 3)}</td>
          <td class="num">${fmt(row.min_source_padj, 3)}</td>
          <td class="num">${fmt(row.normalized_rank, 4)}</td>
          <td>${externalLink(row.source_url, sourceLabel)}</td>
        </tr>`;
      }).join("");
      $("detail").innerHTML = `
        <div class="kv">
          <div class="metric" tabindex="0" data-tip="Primary quality-weighted DEGORA rank of this gene among all scored genes (1 = strongest)."><span>Rank</span><strong>${esc(primaryRankLabel(gene))}</strong></div>
          <div class="metric" tabindex="0" data-tip="Where the primary quality-weighted rank sits as a percent of all scored genes (e.g. top 1%)."><span>Top fraction</span><strong>${esc(topPercentLabel(gene))}</strong></div>
          <div class="metric" tabindex="0" data-tip="Coarse confidence tier (A strongest to D weakest) from rank, support, and direction."><span>Evidence tier</span><strong>${tier(gene.evidence_tier)}</strong></div>
          <div class="metric" tabindex="0" data-tip="The DEGORA quality-weighted prioritization score: a relative index, not a probability."><span>Score</span><strong>${fmt(primaryScore(gene), 2)}</strong></div>
          <div class="metric" tabindex="0" data-tip="How many independent source units support this gene (one paper = one unit)."><span>Source support</span><strong>${esc(gene.support_label)}</strong></div>
          <div class="metric" tabindex="0" data-tip="Primary quality-weighted consensus regulation direction (up/down/flat) across the supporting sources."><span>Direction</span><strong>${esc(primaryDirectionLabel(gene))}</strong></div>
          <div class="metric" tabindex="0" data-tip="Evidence-strength component combining repeated support and signal magnitude."><span>Evidence</span><strong>${fmt(gene.evidence_score, 2)}</strong></div>
          <div class="metric" tabindex="0" data-tip="Contribution from how highly this gene ranked within each source's DEG list."><span>Rank signal</span><strong>${fmt(gene.rank_score_component, 2)}</strong></div>
          <div class="metric" tabindex="0" data-tip="Sample-size-weighted mean log2 fold-change across supporting source units."><span>Weighted LFC</span><strong>${fmt(gene.weighted_lfc, 2)}</strong></div>
          <div class="metric" tabindex="0" data-tip="Effect/rank/direction-focused prioritization score."><span>Priority</span><strong>${fmt(gene.priority_score, 2)}</strong></div>
          <div class="metric" tabindex="0" data-tip="Conditional 0-100 summary over available diagnostics. Support, source quality, and direction confidence are mandatory; LOO participates only when evaluable, with available weights renormalized. This field does not determine the primary rank."><span>Reliability · ${esc(reliabilityBasisLabel(gene))}</span><strong>${fmt(gene.evidence_reliability_score, 2)}</strong></div>
          <div class="metric" tabindex="0" data-tip="Quality-weighted direction-consistency index, shrunk toward 50% when evidence is weak or discordant."><span>Direction confidence</span><strong>${fmt(Number(primaryDirectionConfidence(gene)) * 100, 1)}%</strong></div>
          <div class="metric" tabindex="0" data-tip="Priority-rank stability across global leave-one-source-unit-out folds. N/A means no fold kept this gene eligible under min_studies; it is not zero stability."><span>LOO stability</span><strong>${esc(looStabilityLabel(gene))}</strong></div>
        </div>
        <p class="sources">${badge(primaryDirection(gene))} DEGORA quality-weighted score is a relative prioritization score, not a probability.</p>
        <p class="sources"><strong>Source units:</strong> ${esc(gene.source_units || "unavailable")}</p>
        <div class="ev-scroll">
        <table>
          <thead>
            <tr>
              <th data-tip="Study/contrast ID(s) contributing this row; one source unit may bundle several contrasts.">Study</th>
              <th data-tip="Independent source unit (paper/dataset); contrasts from one paper collapse into one unit.">Unit</th>
              <th data-tip="Evidence provenance class. Derived matrices are explicitly exploratory; author DEG tables remain distinct.">Evidence type</th>
              <th data-tip="Whether the source table contains full tested results, significant genes only, or an unresolved scope.">Scope</th>
              <th data-tip="Analysis pipeline of the source DEG table (e.g. DESeq2, limma).">Pipeline</th>
              <th data-tip="Assay type of the source (e.g. RNA-seq, microarray).">Assay</th>
              <th class="num" data-tip="Control / treatment biological replicate counts reported for the source contrast.">n C/T</th>
              <th class="num" data-tip="Source reliability weight: fixed heuristic, sample-size-aware source-quality weight.">Rel</th>
              <th class="num" data-tip="log2 fold-change contributed by this source unit.">LFC</th>
              <th class="num" data-tip="Smallest source p-value among rows collapsed into this source unit.">p</th>
              <th class="num" data-tip="Smallest adjusted p-value among rows collapsed into this source unit.">FDR</th>
              <th class="num" data-tip="Within-study normalized rank of the gene (near 0 = top of that study's list).">Rank</th>
              <th data-tip="Original public source URL when available.">Source</th>
            </tr>
          </thead>
          <tbody>${evidenceRows}</tbody>
        </table>
        </div>
      `;
    }

    async function exportGeneRanking() {
      const context = currentAtlasContext();
      const generation = atlasContextGeneration;
      if (context.kind === "none" || genesLoading) return;
      const button = $("exportGenes");
      button.disabled = true;
      button.textContent = "Preparing CSV...";
      $("status").textContent = "Preparing complete filtered ranking...";
      const exportQuery = currentQuery();
      const records = [];
      let offset = 0;
      try {
        while (true) {
          const params = new URLSearchParams(exportQuery);
          params.set("limit", "500");
          params.set("offset", String(offset));
          const data = await getJson(`${atlasApi("/genes", context)}?${params.toString()}`);
          if (!atlasContextIsCurrent(context, generation)) return;
          const batch = data.genes || [];
          records.push(...batch);
          offset += batch.length;
          if (!batch.length || offset >= Number(data.count || 0)) break;
        }
        const headers = [
          "quality_weighted_rank", "evidence_tier", "gene_symbol", "quality_weighted_top_percent",
          "quality_weighted_score", "n_source_units", "quality_weighted_sign_concordance",
          "weighted_lfc", "quality_weighted_consensus_direction", "support_label"
        ];
        const rows = records.map((gene) => [
          primaryRank(gene), gene.evidence_tier, gene.gene_symbol, primaryTopPercent(gene),
          primaryScore(gene), gene.n_source_units, primaryConcordance(gene), gene.weighted_lfc,
          primaryDirection(gene), gene.support_label
        ]);
        const runLabel = context.kind === "run" ? context.runId : "submitted";
        downloadCsv(`degora_${context.species}_${runLabel}_ranking.csv`, headers, rows);
        $("status").textContent = `Downloaded ${records.length.toLocaleString()} filtered genes`;
      } catch (error) {
        if (atlasContextIsCurrent(context, generation)) $("status").textContent = `CSV export failed: ${error.message}`;
      } finally {
        button.textContent = "Download ranking CSV";
        if (atlasContextIsCurrent(context, generation)) {
          window.setTimeout(() => {
            if (atlasContextIsCurrent(context, generation)) updateGeneStatus();
          }, 1200);
        }
      }
    }

    function exportCurrentEvidence() {
      if (!currentGeneDetail) return;
      const gene = currentGeneDetail.gene;
      const headers = [
        "gene_symbol", "quality_weighted_rank", "evidence_tier", "quality_weighted_score",
        "consensus_direction", "study_ids", "source_unit_id", "source_input_type", "table_scope",
        "pipeline", "assay_type", "n_ctrl", "n_treat", "source_reliability_weight", "lfc",
        "min_source_pvalue", "min_source_padj", "normalized_rank", "source_url"
      ];
      const rows = (currentGeneDetail.evidence || []).map((row) => [
        gene.gene_symbol, primaryRank(gene), gene.evidence_tier, primaryScore(gene), primaryDirection(gene),
        row.contributing_study_ids || row.study_id, row.source_unit_id,
        row.contributing_source_input_types || row.source_input_type, row.table_scope,
        row.contributing_pipelines || row.pipeline, row.contributing_assay_types || row.assay_type,
        row.n_ctrl, row.n_treat, row.source_reliability_weight, row.lfc, row.min_source_pvalue,
        row.min_source_padj, row.normalized_rank, row.source_url
      ]);
      downloadCsv(`degora_${currentAtlasContext().species}_${gene.gene_symbol}_evidence.csv`, headers, rows);
    }

    $("discoverNav").addEventListener("click", () => {
      if (!discoveryOpened) {
        discoveryOpened = true;
        setSpecies(preferredDiscoverySpecies);
      }
      showView("discover");
    });
    $("atlasNav").addEventListener("click", () => showView("atlas"));
    document.querySelectorAll("[data-species]").forEach((button) => {
      button.addEventListener("click", () => setSpecies(button.dataset.species));
    });
    document.querySelector(".species-tabs").addEventListener("keydown", (event) => {
      if (!["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
      event.preventDefault();
      const nextSpecies = event.key === "Home" || event.key === "ArrowLeft" ? "human" : "mouse";
      setSpecies(nextSpecies);
      $(`${nextSpecies}SpeciesTab`).focus();
    });
    $("discoverySearch").addEventListener("click", () => {
      if (activeDiscoveryState().loading) return;
      void searchStudies({ resetPage: true });
    });
    $("discoverySearchBoth").addEventListener("click", searchBothSpecies);
    $("discoveryQuery").addEventListener("keydown", (event) => {
      // Same guard as the button: Enter used to start a second server-side job.
      if (event.key !== "Enter" || activeDiscoveryState().loading) return;
      void searchStudies({ resetPage: true });
    });
    $("discoveryQuery").addEventListener("input", (event) => event.target.setCustomValidity(""));
    $("discoveryPrev").addEventListener("click", () => {
      const state = activeDiscoveryState();
      if (state.page > 1) { state.page -= 1; void reloadSearchRecords(); }
    });
    $("discoveryNext").addEventListener("click", () => {
      const state = activeDiscoveryState();
      if (state.hasNext) { state.page += 1; void reloadSearchRecords(); }
    });
    $("discoveryResults").addEventListener("click", (event) => {
      if (event.target.closest("#cancelSearchJob")) {
        void cancelDiscoveryJob("search");
        return;
      }
      const retry = event.target.closest("[data-retry-search]");
      if (retry) {
        void searchStudies({ resetPage: true });
        return;
      }
      const inspect = event.target.closest("[data-study-inspect]");
      if (inspect) {
        const state = activeDiscoveryState();
        const key = String(inspect.dataset.studyInspect || "");
        const blocked = pageSelectability(state)(key);
        if (blocked) {
          // Inspect used to bypass the selection rules, clearing an existing
          // selection and then failing server-side on the ambiguous identifier.
          state.notice = blocked === "no-id"
            ? "This publication has no usable identifier, so it cannot be prepared."
            : "Several results share this identifier, so a single one cannot be prepared.";
          renderDiscoveryNotice(state);
          return;
        }
        void prepareSelectedStudies({ recordIds: [key] });
        return;
      }
      const mobileOrder = event.target.closest("#mobileStudyOrder");
      if (mobileOrder) {
        const state = activeDiscoveryState();
        state.sort.order = state.sort.order === "asc" ? "desc" : "asc";
        void reloadSearchRecords({ resetPage: true });
        return;
      }
      const sort = event.target.closest("[data-study-sort]");
      if (!sort) return;
      const state = activeDiscoveryState();
      const key = sort.dataset.studySort;
      state.sort = state.sort.key === key
        ? { key, order: state.sort.order === "asc" ? "desc" : "asc" }
        : { key, order: ["relevance", "readiness", "year"].includes(key) ? "desc" : "asc" };
      void reloadSearchRecords({ resetPage: true });
    });
    // `input` fires per keystroke, which is what the debounce below was written
    // for; `change` alone waited for Enter or blur and left the box looking inert.
    const applyResultFilter = (event) => {
      const state = activeDiscoveryState();
      if (event.target.id === "resultFilter") {
        const wanted = event.target.value.trim().slice(0, 100);
        if (wanted === state.textFilter) return;
        state.textFilter = wanted;
        // Debounced: the snapshot is already on the server, but each keystroke
        // would otherwise be a request and the results would flicker under a
        // reader who is still typing.
        window.clearTimeout(resultFilterTimer);
        resultFilterTimer = window.setTimeout(() => {
          void reloadSearchRecords({ resetPage: true, keepFocus: "resultFilter" });
        }, 250);
        return;
      }
      if (event.target.id === "mobileStudySort") {
        const key = event.target.value;
        state.sort = {
          key,
          order: ["relevance", "readiness", "year"].includes(key) ? "desc" : "asc"
        };
        void reloadSearchRecords({ resetPage: true });
        return;
      }
      if (event.target.id === "selectPageStudies") {
        let capped = false;
        selectableKeys(state).forEach((key) => {
          if (!event.target.checked) state.selected.delete(key);
          else if (state.selected.has(key)) return;
          else if (state.selected.size < MAX_SELECTED_STUDIES) state.selected.add(key);
          else capped = true;
        });
        if (capped) {
          state.notice = `Selection limit reached: at most ${MAX_SELECTED_STUDIES} publications across all pages. Clear some to choose others.`;
        } else if (state.noticeLevel !== "error") {
          state.notice = "";
        }
        renderDiscoveryResultsKeepingFocus();
        return;
      }
      if (event.target.classList.contains("study-select")) {
        const accession = String(event.target.dataset.accession || "");
        if (!accession.trim()) {
          // Guard the empty fallback id: it used to enter the set, inflate the
          // counter, and make the whole prepare fail server-side.
          event.target.checked = false;
          state.notice = "This publication has no usable identifier, so it cannot be prepared.";
          renderDiscoveryNotice(state);
          return;
        }
        if (!event.target.checked) state.selected.delete(accession);
        else if (state.selected.size < MAX_SELECTED_STUDIES) state.selected.add(accession);
        else {
          event.target.checked = false;
          state.notice = `Selection limit reached: at most ${MAX_SELECTED_STUDIES} publications across all pages. Clear some to choose others.`;
        }
        if (state.selected.size < MAX_SELECTED_STUDIES && !event.target.checked && state.noticeLevel !== "error") {
          state.notice = "";
        }
        // Re-render so newly capped rows become visibly disabled instead of
        // silently rejecting the next click, without losing keyboard focus.
        renderDiscoveryResultsKeepingFocus();
      }
    };
    $("discoveryResults").addEventListener("change", applyResultFilter);
    $("discoveryResults").addEventListener("input", (event) => {
      if (event.target.id === "resultFilter") applyResultFilter(event);
    });
    $("clearSelected").addEventListener("click", () => {
      const state = activeDiscoveryState();
      state.selected.clear();
      state.notice = "";
      state.noticeLevel = "info";
      renderDiscoveryResults();
    });
    $("resetStudySort").addEventListener("click", () => {
      const state = activeDiscoveryState();
      state.sort = { key: "readiness", order: "desc" };
      void reloadSearchRecords({ resetPage: true });
    });
    $("downloadSearchExcel").addEventListener("click", downloadSearchExcel);
    $("prepareSelected").addEventListener("click", prepareSelectedStudies);
    $("preparedCandidates").addEventListener("change", () => { capturePreparedDraft(); updateAnalysisEligibility(); });
    // The species attestation lives in the card footer, outside #preparedCandidates,
    // so the listener above never saw it: ticking the box left Run disabled until
    // some other field was touched, and unticking it left Run enabled.
    $("speciesConfirmed").addEventListener("change", updateAnalysisEligibility);
    $("preparedCandidates").addEventListener("input", (event) => {
      const filter = event.target.closest(".sample-filter");
      if (filter) {
        const filterRow = filter.closest(".candidate-row");
        if (filterRow) refreshSampleFilter(filterRow);
        return;
      }
      capturePreparedDraft();
      updateAnalysisEligibility();
    });
    $("preparedCandidates").addEventListener("click", (event) => {
      if (event.target.closest("#cancelPrepareJob")) {
        void cancelDiscoveryJob("prepare");
        return;
      }
      const bulk = event.target.closest(".sample-bulk-apply");
      if (bulk) {
        const bulkRow = bulk.closest(".candidate-row");
        if (bulkRow) applySampleBulk(bulkRow, bulk.dataset.group || "");
        return;
      }
      const button = event.target.closest(".clone-author-candidate");
      if (!button) return;
      const row = button.closest(".candidate-row");
      if (!row) return;
      const state = activeDiscoveryState();
      capturePreparedDraft();
      const key = `${row.dataset.candidateId}::clone::${++state.cloneCounter}`;
      state.draft[key] = {
        ...(state.draft[row.dataset.activationKey || row.dataset.candidateId] || {}),
        enabled: true,
        candidateId: row.dataset.candidateId,
        clone: true,
        contrast: "",
        rowFilterValue: ""
      };
      renderPreparedState();
      const clone = $(`preparedCandidates`).querySelector(`[data-activation-key="${CSS.escape(key)}"]`);
      if (clone) clone.scrollIntoView({ behavior: "smooth", block: "nearest" });
    });
    $("backToResults").addEventListener("click", () => $("discoveryResultsCard").scrollIntoView({ behavior: "smooth" }));
    $("preparedCandidates").addEventListener("click", (event) => {
      if (event.target.closest("[data-back-to-results]")) {
        $("discoveryResultsCard").scrollIntoView({ behavior: "smooth", block: "start" });
      }
    });
    $("runDiscoveryAnalysis").addEventListener("click", runSelectedAnalysis);
    $("openAnalysis").addEventListener("click", openDiscoveryAnalysis);
    $("downloadAnalysisExcel").addEventListener("click", downloadAnalysisExcel);

    $("load").addEventListener("click", () => loadGenes());
    $("loadMore").addEventListener("click", () => fetchGenePage(true));
    $("exportGenes").addEventListener("click", exportGeneRanking);
    $("exportEvidence").addEventListener("click", exportCurrentEvidence);
    document.querySelectorAll("[data-sort]").forEach((button) => {
      button.addEventListener("click", () => setSort(button.dataset.sort));
    });
    // Enter worked in the gene box but not in the source-unit box next to it.
    ["query", "minUnits"].forEach((id) => {
      $(id).addEventListener("keydown", (event) => {
        if (event.key === "Enter") loadGenes();
      });
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
    document.addEventListener("focusin", (event) => {
      const tipped = event.target.closest("[data-tip]");
      if (!tipped) return;
      const rect = tipped.getBoundingClientRect();
      placeTip(tipped.getAttribute("data-tip"), rect.left + Math.min(rect.width, 24), rect.bottom);
    });
    document.addEventListener("focusout", hideTip);
    window.addEventListener("scroll", hideTip, true);

    updateSortHeaders();
    initPanelResize();
    setSpecies("human");
    discoveryOpened = true;
    // `degora serve <db>` follows a run, so a database that already holds scored
    // genes is what the user asked to look at: open the atlas on it. Discover
    // stays the landing view for a database with nothing scored yet.
    void (async () => {
      let scored = 0;
      try {
        const meta = await loadDiscoveryDefaults();
        scored = Number((meta && meta.n_gene_scores) || 0);
      } catch (_) {
        scored = 0;
      }
      // The demo's pre-filled search lives in the workspace of the species the
      // demo was built for; `degora demo --species mouse` used to open on Human
      // with an empty box and the keyword parked in an invisible state object.
      if (preferredDiscoverySpecies !== activeSpecies) setSpecies(preferredDiscoverySpecies);
      showView(Number.isFinite(scored) && scored > 0 ? "atlas" : "discover");
    })();
  </script>
</body>
</html>
"""


def _iso_timestamp(value: Any) -> Any:
    """Render an epoch-seconds float as ISO-8601 UTC; pass strings and None through."""

    if value is None or isinstance(value, str):
        return value
    try:
        seconds = float(value)
    except (TypeError, ValueError):
        return value
    if not math.isfinite(seconds):
        return None
    from datetime import datetime, timezone

    return datetime.fromtimestamp(seconds, tz=timezone.utc).isoformat().replace("+00:00", "Z")


def _jsonable(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        # json.dumps serialises a tuple as an array but a NaN inside it as the bare
        # token NaN, which JSON.parse in the browser refuses.
        return [_jsonable(item) for item in value]
    return value


def _connect(db_path: Path) -> sqlite3.Connection:
    if not db_path.exists():
        # Expose only the filename over the wire: the absolute path can leak local
        # directory names to a client if the server is bound to a non-loopback host.
        raise FileNotFoundError(f"DEGORA database is not available: {Path(db_path).name}")
    uri = f"{db_path.resolve().as_uri()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def _looks_like_local_path(value: Any) -> bool:
    text = str(value).strip()
    if not text:
        return False
    # Strip a file:/file:// scheme so file:///absolute/path reaches the
    # is_absolute() check instead of passing through as an opaque URL string.
    lowered = text.lower()
    if lowered.startswith("file://"):
        text = text[len("file://"):]
    elif lowered.startswith("file:"):
        text = text[len("file:"):]
    if not text:
        return False
    if Path(text).is_absolute():
        return True
    if len(text) >= 3 and text[1] == ":" and text[2] in {"/", "\\"}:
        return True
    return text.startswith(("/", "\\\\"))


LOCAL_PATH_REDACTION = "[redacted: local path]"
# A POSIX path worth redacting descends ("/etc/passwd"), names a file
# ("/catalog.csv") or starts at a known root. The looser "anything after a
# slash" rule this used to be turned "1%/21% O2", "ratio (A)/(B)" and
# "log2FC>1 &/or padj<0.05" - ordinary experimental metadata - into
# "[redacted: local path]" on a network-shared browser. discovery_store learned
# the same lesson first; its rule is reused here rather than copied loosely.
_LOCAL_PATH_TOKEN_RE = re.compile(
    r"(?i)(?:"
    r"file:///[^\s,;\"')\]]+|"
    r"file:/[^\s,;\"')\]]+|"
    r"(?<![A-Za-z0-9])[A-Z]:[\\/][^\s,;\"')\]]+|"
    r"\\\\[^\s,;\"')\]]+\\[^\s,;\"')\]]+"
    r")"
)


def _contains_local_path(value: Any) -> bool:
    from .discovery_store import _POSIX_ABSOLUTE_PATH_RE

    text = str(value).strip()
    if not text:
        return False
    if _looks_like_local_path(text):
        return True
    if any(_looks_like_local_path(part) for part in text.replace("\n", ";").split(";")):
        return True
    if _LOCAL_PATH_TOKEN_RE.search(text) is not None:
        return True
    return _POSIX_ABSOLUTE_PATH_RE.search(text) is not None


# Redact by VALUE, not by key suffix: free-text fields such as source_url,
# contributing_source_urls, and notes can hold an absolute local path or a
# file:// URL, which a key-suffix filter (_path/_dir only) would leak to a
# network client. _looks_like_local_path requires absolute-path/drive-letter/UNC
# shapes, so legitimate https:// URLs and plain text are preserved.
def _redact_meta_for_network(meta: dict[str, str]) -> dict[str, str]:
    redacted = dict(meta)
    for key, value in meta.items():
        if _contains_local_path(value):
            redacted[key] = LOCAL_PATH_REDACTION
    return redacted


def _redact_record_paths_for_network(record: dict[str, Any]) -> dict[str, Any]:
    redacted = dict(record)
    for key, value in record.items():
        if _contains_local_path(value):
            redacted[key] = LOCAL_PATH_REDACTION
    return redacted


def _redact_records_paths_for_network(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [_redact_record_paths_for_network(record) for record in records]


def _redact_local_paths_for_network(value: Any) -> Any:
    """Recursively redact local paths from persisted discovery artifacts."""

    if isinstance(value, dict):
        return {key: _redact_local_paths_for_network(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_redact_local_paths_for_network(item) for item in value]
    if isinstance(value, tuple):
        return [_redact_local_paths_for_network(item) for item in value]
    return LOCAL_PATH_REDACTION if _contains_local_path(value) else value


def _row_dicts(cursor: sqlite3.Cursor) -> list[dict[str, Any]]:
    return [dict(row) for row in cursor.fetchall()]


def _one_row(cursor: sqlite3.Cursor) -> dict[str, Any] | None:
    row = cursor.fetchone()
    return dict(row) if row is not None else None


def _normalize_gene_api_record(record: dict[str, Any]) -> dict[str, Any]:
    """Expose SQLite-backed logical fields with JSON boolean semantics."""

    if "loo_component_available" in record and record["loo_component_available"] is not None:
        record["loo_component_available"] = bool(record["loo_component_available"])
    return record


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
    "rank": (PRIMARY_RANK_COLUMN, "ASC", "degora_rank"),
    "tier": ("evidence_tier", "ASC", "evidence_tier"),
    "gene": ("gene_symbol", "ASC", "gene_symbol"),
    "top": (PRIMARY_TOP_PERCENT_COLUMN, "ASC", "top_percent"),
    "score": (PRIMARY_SCORE_COLUMN, "DESC", "degora_score"),
    "units": ("n_source_units", "DESC", "n_source_units"),
    "sign": (PRIMARY_CONCORDANCE_COLUMN, "DESC", "sign_concordance"),
    "lfc": ("weighted_lfc", "DESC", "weighted_lfc"),
}


def _table_columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})").fetchall()}


def _column_or_fallback(preferred: str, fallback: str, available: set[str]) -> str:
    if preferred in available:
        return preferred
    if fallback in available:
        return fallback
    raise ValueError(f"database is missing required gene column: {preferred} or {fallback}")


def _gene_order_clause(params: dict[str, list[str]], available_columns: set[str]) -> tuple[str, str, str]:
    sort_key = _text_param(params, "sort", "rank", maximum=32).lower()
    if sort_key not in GENE_SORT_COLUMNS:
        allowed = ", ".join(sorted(GENE_SORT_COLUMNS))
        raise ValueError(f"sort must be one of: {allowed}")
    preferred_column, default_order, fallback_column = GENE_SORT_COLUMNS[sort_key]
    column = _column_or_fallback(preferred_column, fallback_column, available_columns)
    raw_order = _text_param(params, "order", default_order.lower(), maximum=8).lower()
    if raw_order in {"asc", "ascending"}:
        order = "ASC"
    elif raw_order in {"desc", "descending"}:
        order = "DESC"
    else:
        raise ValueError("order must be asc or desc")
    rank_column = _column_or_fallback(PRIMARY_RANK_COLUMN, "degora_rank", available_columns)
    return f"{column} IS NULL ASC, {column} {order}, {rank_column} ASC, gene_symbol ASC", sort_key, order.lower()


DISCOVERY_SEARCH_SORT_COLUMNS = {
    "relevance",
    "relevance_score",
    "readiness",
    "readiness_score",
    "title",
    "paper_title",
    "authors",
    "authors_display",
    "journal",
    "year",
    "data_sources",
}
DISCOVERY_SEARCH_EXPORT_SHEETS = [
    "Query",
    "Publications",
    "Identifiers",
    "Linked datasets",
    "Candidate routes",
    "Species decisions",
    "Provider events",
]
FORMULA_PREFIXES = ("=", "+", "-", "@")


def _canonical_species(value: str) -> str:
    species = value.strip().lower()
    if species not in {"human", "mouse"}:
        raise ValueError("species must be human or mouse")
    return species


def _search_result_records(snapshot: dict[str, Any]) -> list[dict[str, Any]]:
    for key in ("records", "publications", "studies", "items", "results"):
        value = snapshot.get(key)
        if isinstance(value, list):
            return [dict(item) for item in value if isinstance(item, dict)]
    return []


def _record_sort_value(record: dict[str, Any], sort_by: str) -> Any:
    aliases = {
        "title": ("title", "paper_title"),
        "paper_title": ("paper_title", "title"),
        "authors": ("authors_display", "authors"),
        "authors_display": ("authors_display", "authors"),
        "readiness": ("readiness_score", "readiness"),
        "relevance": ("relevance_score", "relevance"),
        "data_sources": ("data_sources", "linked_datasets", "datasets"),
    }
    for key in aliases.get(sort_by, (sort_by,)):
        value = record.get(key)
        if value is not None:
            if isinstance(value, list):
                return ", ".join(str(item) for item in value)
            return value
    return ""


def _page_snapshot_fallback(
    snapshot: dict[str, Any],
    *,
    page: int,
    page_size: int,
    sort_by: str,
    sort_order: str,
) -> dict[str, Any]:
    records = _search_result_records(snapshot)
    reverse = sort_order == "desc"
    if sort_by != "relevance":
        records = sorted(
            records,
            key=lambda record: (
                _record_sort_value(record, sort_by) in {"", None},
                str(_record_sort_value(record, sort_by)).lower()
                if not isinstance(_record_sort_value(record, sort_by), (int, float))
                else _record_sort_value(record, sort_by),
            ),
            reverse=reverse,
        )
    total = len(records)
    offset = (page - 1) * page_size
    return {
        "records": records[offset : offset + page_size],
        "total": total,
        "page": page,
        "page_size": page_size,
        "total_pages": max(1, math.ceil(total / page_size)) if total else 0,
        "sort_by": sort_by,
        "sort_order": sort_order,
        "has_next": offset + page_size < total,
    }


def _call_page_publication_snapshot(
    snapshot: dict[str, Any],
    *,
    page: int,
    page_size: int,
    sort_by: str,
    sort_order: str,
    text_filter: str = "",
) -> dict[str, Any]:
    try:
        from .discovery_federated import page_publication_snapshot
    except ImportError:
        return _page_snapshot_fallback(
            snapshot,
            page=page,
            page_size=page_size,
            sort_by=sort_by,
            sort_order=sort_order,
        )
    arguments: dict[str, Any] = {
        "page": page,
        "page_size": page_size,
        "sort_by": sort_by,
        "sort_order": sort_order,
    }
    # An implementation without the filter parameter still pages correctly; passing
    # it unconditionally would turn "cannot narrow the results" into a 500.
    if _accepts_keyword(page_publication_snapshot, "text_filter"):
        arguments["text_filter"] = text_filter
    result = page_publication_snapshot(snapshot, **arguments)
    if isinstance(result, dict):
        return result
    raise ValueError("page_publication_snapshot must return a dictionary")


def _accepts_keyword(func: Any, name: str) -> bool:
    """Whether ``func`` will accept ``name`` as a keyword argument."""

    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):
        return False
    parameters = signature.parameters
    if any(parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in parameters.values()):
        return True
    return name in parameters


def _accepts_progress_callback(func: Any) -> bool:
    """Report whether ``func`` takes a ``progress`` keyword.

    Test doubles replace ``search_publications`` with strict keyword-only stubs,
    so an unconditional ``progress=`` would raise ``TypeError`` against them.
    """

    try:
        parameters = inspect.signature(func).parameters
    except (TypeError, ValueError):  # pragma: no cover - builtins and C callables.
        return False
    if "progress" in parameters:
        return True
    return any(parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in parameters.values())


def _call_with_optional_progress(func: Any, *args: Any, progress: Any = None, **kwargs: Any) -> Any:
    """Call ``func`` forwarding ``progress`` only when it accepts the keyword.

    Test doubles replace these functions with strict signatures, so an
    unconditional ``progress=`` would raise ``TypeError`` against them.
    """

    if progress is not None and _accepts_progress_callback(func):
        kwargs["progress"] = progress
    return func(*args, **kwargs)


def _call_with_optional_keywords(func: Any, *args: Any, optional_keywords: dict[str, Any], **kwargs: Any) -> Any:
    for key, value in optional_keywords.items():
        if key == "before_publish" and _accepts_keyword(func, key):
            kwargs[key] = value
        elif value is not None and _accepts_keyword(func, key):
            kwargs[key] = value
    return func(*args, **kwargs)


def _commit_discovery_job(manager: Any, job_id: str) -> None:
    commit = getattr(manager, "commit", None)
    if not callable(commit):
        raise RuntimeError("the discovery job manager does not provide a commit barrier")
    commit(job_id)


def _call_search_publications(
    query: str,
    species: str,
    *,
    limit: int,
    progress: Callable[[float, str], None] | None = None,
) -> dict[str, Any]:
    from .discovery_federated import search_publications

    call_kwargs: dict[str, Any] = {"query": query, "species": species, "limit": limit}
    if progress is not None and _accepts_progress_callback(search_publications):
        call_kwargs["progress"] = progress
    result = search_publications(**call_kwargs)
    if isinstance(result, dict):
        snapshot = dict(result)
    elif isinstance(result, list):
        snapshot = {"records": result}
    else:
        raise ValueError("search_publications must return a dictionary or list")
    snapshot.setdefault("query", query)
    snapshot.setdefault("species", {"key": species, "label": species.title()})
    snapshot.setdefault("limit", limit)
    snapshot.setdefault("status", "complete")
    snapshot.setdefault("records", _search_result_records(snapshot))
    snapshot.setdefault("total", len(_search_result_records(snapshot)))
    return snapshot


def _formula_neutral(value: Any) -> Any:
    if isinstance(value, (dict, list)):
        value = json.dumps(_jsonable(value), sort_keys=True)
    if not isinstance(value, str):
        return value
    stripped = value.lstrip(" \t\r\n")
    if stripped.startswith(FORMULA_PREFIXES) or stripped.upper() in EXCEL_ERROR_LITERALS:
        return "'" + value
    return value


def _first_value(record: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        value = record.get(key)
        if value not in (None, ""):
            return value
    return ""


def _normalize_publication_identifier(value: Any) -> str:
    text = str(value or "").strip()
    lowered = text.lower()
    if lowered.startswith("pmid:"):
        digits = "".join(char for char in text.split(":", 1)[1] if char.isdigit())
        return f"pmid:{digits}" if digits else ""
    if lowered.startswith("doi:"):
        doi = text.split(":", 1)[1].strip().lower()
        doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi)
        return f"doi:{doi}" if doi else ""
    if re.fullmatch(r"GSE[0-9]+", text, flags=re.IGNORECASE):
        return text.lower()
    return lowered


def _value_list(value: Any) -> list[Any]:
    if value in (None, ""):
        return []
    if isinstance(value, (list, tuple, set, frozenset)):
        return list(value)
    return [value]


def _publication_record_keys(record: dict[str, Any]) -> set[str]:
    values: list[Any] = [
        record.get("canonical_id"),
        record.get("publication_id"),
        record.get("source_unit_id"),
        record.get("accession"),
    ]
    values.extend(_value_list(record.get("geo_accessions")))
    values.extend(f"PMID:{value}" for value in _value_list(record.get("pubmed_ids")) + _value_list(record.get("pmid")))
    values.extend(f"DOI:{value}" for value in _value_list(record.get("dois")) + _value_list(record.get("doi")))
    return {
        normalized
        for value in values
        if (normalized := _normalize_publication_identifier(value))
    }


def _select_publication_snapshot_records(
    snapshot: dict[str, Any],
    identifiers: list[Any],
) -> list[dict[str, Any]]:
    if not 1 <= len(identifiers) <= 20:
        raise ValueError("record_ids must contain between 1 and 20 selections")
    normalized = [_normalize_publication_identifier(value) for value in identifiers]
    if any(not value for value in normalized):
        raise ValueError("record_ids contains an invalid publication identifier")
    if len(set(normalized)) != len(normalized):
        raise ValueError("record_ids contains duplicate selections")
    records = _search_result_records(snapshot)
    selected: list[dict[str, Any]] = []
    for original, key in zip(identifiers, normalized, strict=True):
        matches = [record for record in records if key in _publication_record_keys(record)]
        if not matches:
            raise ValueError(f"selected publication was not found in the persisted search: {original}")
        if len(matches) > 1:
            raise ValueError(f"selected identifier is ambiguous in the persisted search: {original}")
        selected.append(matches[0])
    return selected


def _load_discovery_store_classes() -> tuple[type[Any], type[Any]]:
    from .discovery_store import DiscoveryJobManager, DiscoveryStateStore

    return DiscoveryStateStore, DiscoveryJobManager


def _api_job_status(status: str) -> str:
    return "complete" if status == "completed" else status


def _api_job_progress(value: Any) -> float | None:
    """Clamp a stored job progress fraction to [0, 1] for the browser."""

    if value is None or isinstance(value, bool):
        return None
    try:
        fraction = float(value)
    except (TypeError, ValueError):
        return None
    if fraction != fraction or fraction in {float("inf"), float("-inf")}:
        return None
    return round(min(1.0, max(0.0, fraction)), 4)


def _clean_job_message(value: Any) -> str:
    """Collapse a stored job message to a single short display line."""

    if not isinstance(value, str):
        return ""
    collapsed = " ".join(value.split())
    return collapsed[:160]


class DegoraRequestHandler(BaseHTTPRequestHandler):
    """Serve the static browser UI and JSON endpoints."""

    server: "DegoraHttpServer"
    # Reap half-open / slow-loris connections so a peer cannot pin a worker thread
    # indefinitely (each connection gets its own thread; see DegoraHttpServer).
    timeout = 30

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        if self.server.quiet:
            return
        redacted_args = tuple(redact_token_from_url_text(arg) if isinstance(arg, str) else arg for arg in args)
        super().log_message(format, *redacted_args)

    def _token_authorized(self, parsed: Any) -> bool:
        expected = self.server.access_token
        if not expected:
            return True
        supplied = self.headers.get("X-DEGORA-Token", "")
        if not supplied:
            supplied = parse_qs(parsed.query).get("token", [""])[0]
        # compare_digest refuses non-ASCII str; a mistyped or mangled token used to
        # raise TypeError out of the handler and drop the connection instead of 401.
        return secrets.compare_digest(
            str(supplied).encode("utf-8", "surrogateescape"),
            str(expected).encode("utf-8", "surrogateescape"),
        )

    def _host_header_authorized(self) -> bool:
        """Reject DNS-rebinding Host headers when the server is loopback-only."""

        bound_host = str(self.server.server_address[0])
        if not _is_loopback_host(bound_host):
            return True
        raw_host = str(self.headers.get("Host") or "").strip()
        if not raw_host:
            return False
        try:
            parsed_host = urlsplit(f"//{raw_host}")
            hostname = parsed_host.hostname or ""
            header_port = parsed_host.port
        except ValueError:
            return False
        if not _is_loopback_host(hostname):
            return False
        return header_port in {None, int(self.server.server_address[1])}

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        try:
            if not self._host_header_authorized():
                self._send_json(
                    {"error": "invalid Host header for loopback DEGORA server"},
                    status=HTTPStatus.MISDIRECTED_REQUEST,
                )
                return
            is_index = parsed.path in {"/", "/index.html"}
            if not is_index and not self._token_authorized(parsed):
                self._send_json({"error": "missing or invalid DEGORA access token"}, status=HTTPStatus.UNAUTHORIZED)
                return
            if is_index:
                self._send_html(INDEX_HTML)
            elif parsed.path == "/api/health":
                self._send_json(self._health())
            elif parsed.path == "/api/meta":
                self._send_json({"meta": self._meta()})
            elif parsed.path == "/api/studies":
                self._send_json({"studies": self._studies()})
            elif parsed.path == "/api/discovery/search":
                self._require_action_header()
                self._send_json(self._discovery_search(parse_qs(parsed.query)))
            elif re.fullmatch(r"/api/discovery/jobs/[a-f0-9]{16}", parsed.path):
                self._require_loopback_discovery()
                job_id = parsed.path.split("/")[4]
                self._send_json({"job": self._discovery_job(job_id)})
            elif re.fullmatch(r"/api/discovery/searches/[a-f0-9]{16}", parsed.path):
                self._require_loopback_discovery()
                search_id = parsed.path.split("/")[4]
                self._send_json({"search": self._discovery_publication_search(search_id)})
            elif re.fullmatch(r"/api/discovery/searches/[a-f0-9]{16}/records", parsed.path):
                self._require_loopback_discovery()
                search_id = parsed.path.split("/")[4]
                self._send_json(self._discovery_publication_records(search_id, parse_qs(parsed.query)))
            elif re.fullmatch(r"/api/discovery/searches/[a-f0-9]{16}/export\.xlsx", parsed.path):
                self._require_loopback_discovery()
                search_id = parsed.path.split("/")[4]
                self._send_xlsx_bytes(
                    self._discovery_publication_export(search_id),
                    filename=f"DEGORA_discovery_{search_id}.xlsx",
                )
            elif re.fullmatch(r"/api/discovery/runs/[a-f0-9]{16}/summary", parsed.path):
                run_id = parsed.path.split("/")[4]
                self._send_json(self._discovery_run_summary(run_id))
            elif re.fullmatch(r"/api/discovery/runs/[a-f0-9]{16}/genes", parsed.path):
                run_id = parsed.path.split("/")[4]
                self._send_json(self._genes(parse_qs(parsed.query), db_path=self._discovery_run_db(run_id)))
            elif re.fullmatch(r"/api/discovery/runs/[a-f0-9]{16}/genes/[^/]+", parsed.path):
                parts = parsed.path.split("/")
                run_id = parts[4]
                symbol = unquote(parts[6]).upper()
                self._send_json(self._gene_detail(symbol, db_path=self._discovery_run_db(run_id)))
            elif re.fullmatch(r"/api/discovery/runs/[a-f0-9]{16}/export\.xlsx", parsed.path):
                if not _is_loopback_host(self.server.server_address[0]):
                    self._send_json(
                        {"error": "discovery artifact downloads are available only on a loopback server"},
                        status=HTTPStatus.FORBIDDEN,
                    )
                    return
                run_id = parsed.path.split("/")[4]
                self._send_file_download(
                    self._discovery_run_excel(run_id),
                    content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    filename="DEGORA_output.xlsx",
                )
            elif parsed.path == "/api/genes":
                self._send_json(self._genes(parse_qs(parsed.query)))
            elif parsed.path.startswith("/api/genes/"):
                symbol = unquote(parsed.path.removeprefix("/api/genes/")).upper()
                self._send_json(self._gene_detail(symbol))
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except FileNotFoundError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.NOT_FOUND)
        except PermissionError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.FORBIDDEN)
        except ValueError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
        except RuntimeError as exc:
            if exc.__class__.__name__ == "DiscoveryUnavailableError":
                self._send_json({"error": str(exc)}, status=HTTPStatus.SERVICE_UNAVAILABLE)
                return
            raise
        except sqlite3.Error as exc:
            self._send_json({"error": f"database error: {exc}"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)
        except (ConnectionError, socket.timeout, TimeoutError):
            # The peer went away (a closed tab, an aborted download) or a body
            # stopped arriving. There is nobody to answer and nothing on disk to
            # blame; writing a second response here used to raise inside the
            # handler and print a nested traceback for every closed tab.
            return
        except OSError as exc:
            detail = exc.strerror or exc.__class__.__name__
            self._send_json({"error": f"discovery filesystem error: {detail}"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def do_POST(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        try:
            if not self._host_header_authorized():
                self._send_json(
                    {"error": "invalid Host header for loopback DEGORA server"},
                    status=HTTPStatus.MISDIRECTED_REQUEST,
                )
                return
            if not self._token_authorized(parsed):
                self._send_json({"error": "missing or invalid DEGORA access token"}, status=HTTPStatus.UNAUTHORIZED)
                return
            if not _is_loopback_host(self.server.server_address[0]):
                self._send_json(
                    {"error": "discovery download and analysis actions are available only on a loopback server"},
                    status=HTTPStatus.FORBIDDEN,
                )
                return
            payload = self._read_json_action()
            if parsed.path == "/api/discovery/searches":
                self._send_json(self._discovery_create_publication_search(payload), status=HTTPStatus.ACCEPTED)
            elif parsed.path == "/api/discovery/prepare":
                self._send_json(self._discovery_prepare(payload), status=HTTPStatus.CREATED)
            elif parsed.path == "/api/discovery/prepare-jobs":
                # Same work as /api/discovery/prepare, run as a job so the browser
                # can show real stage progress instead of a blocking request.
                self._send_json(self._discovery_prepare_job(payload), status=HTTPStatus.ACCEPTED)
            elif parsed.path == "/api/discovery/analyze":
                self._send_json(self._discovery_analyze(payload), status=HTTPStatus.CREATED)
            elif re.fullmatch(r"/api/discovery/jobs/[a-f0-9]{16}/cancel", parsed.path):
                self._send_json(self._discovery_cancel_job(parsed.path.split("/")[4]))
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except FileNotFoundError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.NOT_FOUND)
        except FileExistsError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.CONFLICT)
        except PermissionError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.FORBIDDEN)
        except ValueError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
        except RuntimeError as exc:
            if exc.__class__.__name__ == "DiscoveryUnavailableError":
                self._send_json({"error": str(exc)}, status=HTTPStatus.SERVICE_UNAVAILABLE)
                return
            raise
        except sqlite3.Error as exc:
            self._send_json({"error": f"database error: {exc}"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)
        except (ConnectionError, socket.timeout, TimeoutError):
            # The peer went away (a closed tab, an aborted download) or a body
            # stopped arriving. There is nobody to answer and nothing on disk to
            # blame; writing a second response here used to raise inside the
            # handler and print a nested traceback for every closed tab.
            return
        except OSError as exc:
            detail = exc.strerror or exc.__class__.__name__
            self._send_json({"error": f"discovery filesystem error: {detail}"}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def _read_json_action(self) -> dict[str, Any]:
        content_type = self.headers.get("Content-Type", "").split(";", 1)[0].strip().lower()
        if content_type != "application/json":
            raise ValueError("discovery actions require Content-Type: application/json")
        self._require_action_header()
        raw_length = self.headers.get("Content-Length", "")
        try:
            length = int(raw_length)
        except ValueError as exc:
            raise ValueError("a valid Content-Length header is required") from exc
        if not 1 <= length <= 64 * 1024:
            raise ValueError("JSON action body must be between 1 and 65536 bytes")
        try:
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError("request body is not valid UTF-8 JSON") from exc
        if not isinstance(payload, dict):
            raise ValueError("request body must be a JSON object")
        return payload

    def _require_action_header(self) -> None:
        if self.headers.get("X-DEGORA-Action") != "1":
            raise ValueError("discovery actions require the X-DEGORA-Action: 1 header")

    def _require_loopback_discovery(self) -> None:
        if not _is_loopback_host(self.server.server_address[0]):
            raise PermissionError("persistent discovery actions are available only on a loopback server")

    def _send_html(self, html: str) -> None:
        encoded = html.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Referrer-Policy", "no-referrer")
        # The page is fully self-contained (one inline style + script, data: favicon),
        # so a strict CSP hardens --allow-network mode against MIME sniffing,
        # clickjacking, and any future escaping regression without breaking the UI.
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("X-Frame-Options", "DENY")
        self.send_header(
            "Content-Security-Policy",
            "default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; "
            "connect-src 'self'; img-src data:; base-uri 'none'; form-action 'none'",
        )
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)

    def _send_json(self, payload: dict[str, Any], *, status: HTTPStatus = HTTPStatus.OK) -> None:
        # allow_nan=False: a non-finite float that escaped _jsonable would be
        # written as the bare token NaN, which no browser JSON parser accepts.
        # Failing loudly here beats shipping a response the page cannot read.
        encoded = json.dumps(_jsonable(payload), sort_keys=True, allow_nan=False).encode("utf-8")
        try:
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(encoded)))
            self.end_headers()
            self.wfile.write(encoded)
        except ConnectionError:
            # The client closed the socket first; there is nobody left to tell.
            return

    def _send_file_download(self, path: Path, *, content_type: str, filename: str) -> None:
        size = path.stat().st_size
        with path.open("rb") as stream:
            try:
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", content_type)
                self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
                self.send_header("X-Content-Type-Options", "nosniff")
                self.send_header("Cache-Control", "no-store")
                self.send_header("Content-Length", str(size))
                self.end_headers()
                shutil.copyfileobj(stream, self.wfile, length=1024 * 1024)
            except ConnectionError:
                # Closing a browser tab or cancelling a large download is a normal
                # client-side event after the response may already have started.
                return

    def _send_xlsx_bytes(self, payload: bytes, *, filename: str) -> None:
        try:
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
        except ConnectionError:
            return

    def _health(self) -> dict[str, Any]:
        version_info = runtime_version_info()
        with closing(_connect(self.server.db_path)) as connection:
            columns = _table_columns(connection, "genes")
            rank_column = _column_or_fallback(PRIMARY_RANK_COLUMN, "degora_rank", columns)
            gene_count = connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
            study_count = connection.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
            source_unit_count = connection.execute("SELECT COUNT(DISTINCT source_unit_id) FROM studies").fetchone()[0]
            top_gene = connection.execute(
                f"SELECT gene_symbol FROM genes ORDER BY {rank_column} IS NULL ASC, {rank_column} ASC, gene_symbol ASC LIMIT 1"
            ).fetchone()
            meta_rows = connection.execute(
                "SELECT key, value FROM meta WHERE key IN ('degora_version', 'degora_code_revision')"
            ).fetchall()
        db_version_info = {str(key): str(value) for key, value in meta_rows}
        return {
            "status": "ok",
            "db_name": self.server.db_path.name,
            **version_info,
            "database_degora_version": db_version_info.get("degora_version", ""),
            "database_degora_code_revision": db_version_info.get("degora_code_revision", ""),
            "gene_count": gene_count,
            "study_count": study_count,
            "source_unit_count": source_unit_count,
            "top_gene": top_gene[0] if top_gene else None,
        }

    def _meta(self) -> dict[str, str]:
        with closing(_connect(self.server.db_path)) as connection:
            rows = _row_dicts(connection.execute("SELECT key, value FROM meta ORDER BY key"))
        meta = {row["key"]: row["value"] for row in rows}
        if not _is_loopback_host(self.server.server_address[0]):
            return _redact_meta_for_network(meta)
        return meta

    def _studies(self) -> list[dict[str, Any]]:
        with closing(_connect(self.server.db_path)) as connection:
            rows = _row_dicts(connection.execute("SELECT * FROM studies ORDER BY source_unit_id, study_id"))
        if not _is_loopback_host(self.server.server_address[0]):
            return _redact_records_paths_for_network(rows)
        return rows

    def _discovery_search(self, params: dict[str, list[str]]) -> dict[str, Any]:
        from .discovery import DEFAULT_GLOBAL_RANK_LIMIT, DEFAULT_PAGE_SIZE, search_geo

        query = _text_param(params, "q", maximum=200)
        species = _text_param(params, "species", maximum=32)
        page = _int_param(params, "page", 1, minimum=1, maximum=500)
        sort_by = _text_param(params, "sort", "deg_input_priority", maximum=32)
        sort_order = _text_param(params, "order", maximum=4) or None
        return search_geo(
            query,
            species,
            page=page,
            page_size=DEFAULT_PAGE_SIZE,
            assess_files=True,
            global_rank=True,
            global_limit=DEFAULT_GLOBAL_RANK_LIMIT,
            sort_by=sort_by,
            sort_order=sort_order,
        )

    def _discovery_create_publication_search(self, payload: dict[str, Any]) -> dict[str, Any]:
        query = str(payload.get("query") or "").strip()
        if len(query) < 2:
            raise ValueError("query must be at least 2 characters")
        if len(query) > 200:
            raise ValueError("query is too long; maximum length is 200 characters")
        species = _canonical_species(str(payload.get("species") or ""))
        raw_limit = payload.get("limit", 1000)
        if isinstance(raw_limit, bool):
            raise ValueError("limit must be an integer")
        try:
            limit = int(raw_limit)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("limit must be an integer") from exc
        if not 1 <= limit <= 1000:
            raise ValueError("limit must be between 1 and 1000")
        store = self.server.discovery_search_store
        manager = self.server.discovery_job_manager
        search_id = secrets.token_hex(8)
        now = time.time()
        search_payload = {
            "id": search_id,
            "query": query,
            "species": species,
            "limit": limit,
            "status": "queued",
            "error": "",
            "snapshot": {},
            "total": 0,
            "created_at": now,
            "updated_at": now,
        }
        store.save_search(search_id, search_payload)

        def worker(_job_id: str, payload: dict[str, Any], progress: Any) -> dict[str, Any]:
            from .discovery_store import DiscoveryJobCancelled, _sanitize_text

            def mark_search_failed(exc: BaseException) -> None:
                try:
                    current = store.get_search(search_id)
                except BaseException:  # noqa: BLE001 - the job failure is still recorded by the manager.
                    current = None
                failed = dict(current) if isinstance(current, dict) else dict(search_payload)
                failed.update(
                    {
                        "status": "failed",
                        "error": _sanitize_text(str(exc)),
                        "updated_at": time.time(),
                    }
                )
                try:
                    store.save_search(search_id, failed)
                except BaseException:  # noqa: BLE001 - best-effort projection must not mask the real failure.
                    return

            progress(0.02, "Starting publication search.")

            def stage(fraction: float, message: str) -> None:
                # search_publications reports 0..1 over its own work; map that into
                # the 0.03..0.95 band so the job's own bookends stay monotonic.
                progress(0.03 + 0.92 * min(1.0, max(0.0, float(fraction))), message)

            try:
                snapshot = _call_search_publications(
                    payload["query"],
                    payload["species"],
                    limit=int(payload["limit"]),
                    progress=stage,
                )
            except DiscoveryJobCancelled:
                raise
            except BaseException as exc:  # noqa: BLE001 - keep search projection terminal for provider hard failures.
                mark_search_failed(exc)
                raise
            complete = dict(search_payload)
            complete.update(
                {
                    "status": "complete",
                    "error": "",
                    "snapshot": snapshot,
                    "total": snapshot.get("total", len(_search_result_records(snapshot))),
                    "updated_at": time.time(),
                }
            )
            progress(0.97, "Saving the publication snapshot.")
            _commit_discovery_job(manager, _job_id)
            try:
                store.save_search(search_id, complete)
            except DiscoveryJobCancelled:
                raise
            except BaseException as exc:  # noqa: BLE001 - persist hard storage failures too.
                mark_search_failed(exc)
                raise
            progress(1.0, "Publication snapshot persisted.")
            return {"search_id": search_id}

        job = manager.submit(
            "publication_search",
            {"search_id": search_id, "query": query, "species": species, "limit": limit},
            worker,
        )
        return {"job_id": job["job_id"], "search_id": search_id, "status": "queued"}

    @staticmethod
    def _mark_search_cancelled(store: Any, payload: dict[str, Any]) -> None:
        search_id = str(payload.get("search_id") or "")
        if not search_id:
            return
        try:
            search = store.get_search(search_id)
        except Exception:  # noqa: BLE001 - the job is cancelled either way.
            return
        if not isinstance(search, dict) or search.get("status") not in {"queued", "running"}:
            return
        stopped = dict(search)
        stopped.update({"status": "cancelled", "error": "", "updated_at": time.time()})
        try:
            store.save_search(search_id, stopped)
        except Exception:  # noqa: BLE001 - reporting the cancellation matters more.
            return

    def _discovery_cancel_job(self, job_id: str) -> dict[str, Any]:
        """Stop a running search or preparation at the reader's request.

        A job that has already finished is not an error to report - the reader
        pressed a button a moment too late, and the honest answer is that the
        work is done, along with the state it actually reached.
        """

        store = self.server.discovery_search_store
        if store is None:
            raise RuntimeError("discovery is not available on this server")
        if store.get_job(job_id) is None:
            raise FileNotFoundError("discovery search job was not found")

        manager = self.server.discovery_job_manager
        cancel = getattr(manager, "cancel", None) if manager is not None else None
        if cancel is None:
            # A manager without cancellation cannot stop the worker thread, and
            # marking the job cancelled while it kept running and writing would be
            # a worse answer than saying plainly that it cannot be stopped.
            raise RuntimeError("this server cannot cancel a running job")

        cancelled = cancel(job_id)
        if cancelled is not None:
            # A search job persists its snapshot record separately, and that
            # record is what /searches/{id} answers from. Left alone it would
            # report "queued" for the rest of the server's life.
            self._mark_search_cancelled(store, cancelled.get("payload") or {})
        job = self._discovery_job(job_id)
        if cancelled is not None and job["status"] == "complete":
            return {
                "cancelled": False,
                "reason": "The job completed while the stop request was being processed; the completed result was kept.",
                "job": job,
            }
        if cancelled is None:
            if job["status"] == "cancelled":
                return {
                    "cancelled": False,
                    "reason": "The job was already cancelled.",
                    "job": job,
                }
            if job["status"] not in {"complete", "failed", "interrupted", "cancelled"}:
                return {
                    "cancelled": False,
                    "reason": "The job was already saving its result and could not be stopped.",
                    "job": job,
                }
            return {
                "cancelled": False,
                "reason": "The job had already finished before it could be cancelled.",
                "job": job,
            }
        return {
            "cancelled": True,
            # Said here rather than left for the reader to wonder about: stopping
            # does not undo downloads, and it does guarantee no partial result.
            "reason": (
                "Stopped. Files already downloaded are kept and a later run reuses "
                "them; no partial result was recorded for this job."
            ),
            "job": job,
        }

    def _discovery_job(self, job_id: str) -> dict[str, Any]:
        job = self.server.discovery_search_store.get_job(job_id)
        if job is None:
            raise FileNotFoundError("discovery search job was not found")
        raw_error = job.get("error", "")
        if isinstance(raw_error, dict):
            raw_error = raw_error.get("message", "") or json.dumps(raw_error, sort_keys=True)
        status = _api_job_status(str(job["status"]))
        result = job.get("result") if str(job.get("status")) == "completed" else None
        if result is not None and not _is_loopback_host(str(self.server.server_address[0])):
            result = _redact_local_paths_for_network(result)
        search_id = job.get("search_id") or job.get("payload", {}).get("search_id")
        return {
            "id": job.get("id", job.get("job_id")),
            "search_id": search_id,
            "status": status,
            "error": raw_error,
            # progress/message are already persisted by the job manager on every
            # update; surfacing them lets the browser show real stage feedback
            # instead of a constant "search running" placeholder.
            "progress": 1.0 if status == "complete" else _api_job_progress(job.get("progress")),
            "message": "Job completed." if status == "complete" else _clean_job_message(job.get("message")),
            "created_at": job.get("created_at"),
            "updated_at": job.get("updated_at"),
            # Search jobs persist their payload separately and only return an id
            # here; prepare jobs return the whole bundle, which the browser
            # collects from this field when the job completes.
            "result": result,
        }

    def _discovery_publication_search(self, search_id: str) -> dict[str, Any]:
        search = self.server.discovery_search_store.get_search(search_id)
        if search is None:
            raise FileNotFoundError("discovery search was not found")
        snapshot = search.get("snapshot") if isinstance(search.get("snapshot"), dict) else {}
        return {
            "id": search["id"],
            "query": search["query"],
            "species": search["species"],
            "limit": search["limit"],
            "status": search["status"],
            "error": search.get("error", ""),
            "total": search.get("total", snapshot.get("total", len(_search_result_records(snapshot)))),
            # Job records carry ISO-8601 UTC timestamps; the search record stores
            # epoch floats internally. One field name, one type, on the wire.
            "created_at": _iso_timestamp(search.get("created_at")),
            "updated_at": _iso_timestamp(search.get("updated_at")),
        }

    def _discovery_publication_records(self, search_id: str, params: dict[str, list[str]]) -> dict[str, Any]:
        from .discovery import DEFAULT_PAGE_SIZE

        search = self.server.discovery_search_store.get_search(search_id)
        if search is None:
            raise FileNotFoundError("discovery search was not found")
        page = _int_param(params, "page", 1, minimum=1, maximum=500)
        page_size = _int_param(params, "page_size", DEFAULT_PAGE_SIZE, minimum=1, maximum=20)
        sort_by = _text_param(params, "sort_by", "readiness", maximum=40).lower()
        if sort_by not in DISCOVERY_SEARCH_SORT_COLUMNS:
            allowed = ", ".join(sorted(DISCOVERY_SEARCH_SORT_COLUMNS))
            raise ValueError(f"sort_by must be one of: {allowed}")
        sort_order = _text_param(params, "sort_order", "desc", maximum=8).lower()
        if sort_order not in {"asc", "desc"}:
            raise ValueError("sort_order must be asc or desc")
        text_filter = _text_param(params, "filter", "", maximum=100)
        snapshot = search.get("snapshot") if isinstance(search.get("snapshot"), dict) else {}
        page_data = _call_page_publication_snapshot(
            snapshot,
            page=page,
            page_size=page_size,
            sort_by=sort_by,
            sort_order=sort_order,
            text_filter=text_filter,
        )
        page_data.setdefault("records", [])
        page_data.setdefault("total", len(_search_result_records(snapshot)))
        page_data.setdefault("page", page)
        page_data.setdefault("page_size", page_size)
        page_data.setdefault("total_pages", max(1, math.ceil(int(page_data["total"]) / page_size)) if page_data["total"] else 0)
        page_data.setdefault("sort_by", sort_by)
        page_data.setdefault("sort_order", sort_order)
        page_data.setdefault("has_next", page < int(page_data.get("total_pages") or 0))
        return {
            "search": self._discovery_publication_search(search_id),
            **page_data,
        }

    def _append_xlsx_rows(self, sheet: Any, rows: list[dict[str, Any]], headers: list[str]) -> None:
        sheet.append(headers)
        for row in rows:
            sheet.append([_formula_neutral(row.get(header, "")) for header in headers])

    def _discovery_publication_export(self, search_id: str) -> bytes:
        from .discovery_export import build_publication_search_workbook

        search = self.server.discovery_search_store.get_search(search_id)
        if search is None:
            raise FileNotFoundError("discovery search was not found")
        snapshot = search.get("snapshot") if isinstance(search.get("snapshot"), dict) else {}
        enriched = dict(snapshot)
        enriched.setdefault("search_id", search_id)
        enriched.setdefault("query", search.get("query", ""))
        enriched.setdefault("species", search.get("species", ""))
        enriched.setdefault("limit", search.get("limit", 1000))
        return build_publication_search_workbook(enriched)

    def _discovery_prepare_job(self, payload: dict[str, Any]) -> dict[str, Any]:
        """Run a preparation as a job so its stages can be polled.

        Preparation issues dozens of paced network round trips; as a blocking
        request it held a connection open for tens of seconds with no output.
        """

        manager = self.server.discovery_job_manager
        request = dict(payload)

        def worker(_job_id: str, _payload: dict[str, Any], progress: Any) -> dict[str, Any]:
            progress(0.02, "Starting preparation.")
            commit_calls = 0

            def stage(fraction: float, message: str) -> None:
                progress(0.05 + 0.90 * min(1.0, max(0.0, float(fraction))), message)

            def before_publish() -> None:
                nonlocal commit_calls
                if commit_calls:
                    raise RuntimeError("the preparation attempted to publish more than once")
                _commit_discovery_job(manager, _job_id)
                commit_calls += 1

            result = _call_with_optional_keywords(
                self._discovery_prepare,
                request,
                progress=stage,
                before_publish=before_publish,
                optional_keywords={"remember_result": False},
            )
            if commit_calls != 1:
                raise RuntimeError("the preparation returned without reaching its commit barrier")
            self._remember_discovery_bundle(result)
            return result

        job = manager.submit("publication_prepare", {"species": request.get("species")}, worker)
        return {"job_id": job["job_id"], "status": "queued"}

    def _discovery_prepare(
        self,
        payload: dict[str, Any],
        progress: Callable[[float, str], None] | None = None,
        before_publish: Callable[[], None] | None = None,
        remember_result: bool = True,
    ) -> dict[str, Any]:
        from .discovery import normalize_species, prepare_geo_studies

        species = str(payload.get("species") or "")
        spec = normalize_species(species)
        query = str(payload.get("query") or "")
        bundle_id = secrets.token_hex(8)
        target = self.server.discovery_root / spec.key / "bundles" / bundle_id
        publish_reached = before_publish is None

        def guarded_before_publish() -> None:
            nonlocal publish_reached
            if before_publish is None:
                return
            before_publish()
            publish_reached = True

        publish_callback = guarded_before_publish if before_publish is not None else None
        search_id = str(payload.get("search_id") or "").strip()
        record_ids = payload.get("record_ids")
        try:
            if search_id or record_ids is not None:
                if not re.fullmatch(r"[a-f0-9]{16}", search_id):
                    raise ValueError("search_id is invalid")
                if not isinstance(record_ids, list):
                    raise ValueError("record_ids must be a JSON list")
                search = self.server.discovery_search_store.get_search(search_id)
                if search is None:
                    raise FileNotFoundError("persisted discovery search was not found")
                if str(search.get("species") or "").lower() != spec.key:
                    raise ValueError("search species does not match the preparation workspace")
                persisted_query = str(search.get("query") or "").strip()
                if query and query != persisted_query:
                    raise ValueError("preparation query does not match the persisted search")
                from .discovery_federated import resolve_publication_records
                from .discovery_prepare import prepare_publication_records

                snapshot = search.get("snapshot") if isinstance(search.get("snapshot"), dict) else {}
                records = _select_publication_snapshot_records(snapshot, record_ids)
                # Resolution is roughly the first third of the work, preparation the rest.
                records = _call_with_optional_progress(
                    resolve_publication_records,
                    records,
                    spec,
                    progress=(lambda fraction, message: progress(0.30 * fraction, message)) if progress else None,
                )
                result = _call_with_optional_keywords(
                    prepare_publication_records,
                    records,
                    spec.key,
                    query=persisted_query,
                    materialize_dir=target,
                    optional_keywords={
                        "progress": (
                            lambda fraction, message: progress(0.30 + 0.70 * fraction, message)
                        ) if progress else None,
                        "before_publish": publish_callback,
                    },
                )
                result["search_id"] = search_id
            else:
                # Backward-compatible GSE-only action for older clients.
                accessions = payload.get("accessions")
                if not isinstance(accessions, list):
                    raise ValueError("record_ids or legacy accessions must be a JSON list")
                result = _call_with_optional_keywords(
                    prepare_geo_studies,
                    accessions,
                    spec.key,
                    query=query,
                    materialize_dir=target,
                    optional_keywords={"before_publish": publish_callback},
                )
            if not publish_reached:
                raise RuntimeError("the preparation implementation skipped its required commit barrier")
        except BaseException:
            if before_publish is not None:
                shutil.rmtree(target, ignore_errors=True)
            raise
        result["bundle_id"] = bundle_id
        if remember_result:
            self._remember_discovery_bundle(result)
        return result

    def _remember_discovery_bundle(self, result: dict[str, Any]) -> None:
        bundle_id = str(result.get("bundle_id") or "")
        if not re.fullmatch(r"[a-f0-9]{16}", bundle_id):
            raise ValueError("prepared discovery result has an invalid bundle_id")
        with self.server.discovery_lock:
            self.server.remember_discovery(self.server.discovery_bundles, bundle_id, result)
        save_artifact = getattr(self.server.discovery_search_store, "save_artifact", None)
        if callable(save_artifact):
            save_artifact("bundle", bundle_id, result)

    def _discovery_analyze(self, payload: dict[str, Any]) -> dict[str, Any]:
        from .discovery import normalize_species
        from .discovery_run import run_discovery_analysis

        bundle_id = str(payload.get("bundle_id") or "")
        if not re.fullmatch(r"[a-f0-9]{16}", bundle_id):
            raise ValueError("bundle_id is invalid")
        with self.server.discovery_lock:
            prepared = self.server.discovery_bundles.get(bundle_id)
        if prepared is None:
            get_artifact = getattr(self.server.discovery_search_store, "get_artifact", None)
            if callable(get_artifact):
                prepared = get_artifact("bundle", bundle_id)
        if prepared is None:
            raise FileNotFoundError("prepared discovery bundle was not found")
        species = str(payload.get("species") or "")
        spec = normalize_species(species)
        selections = payload.get("selections")
        if not isinstance(selections, list):
            raise ValueError("selections must be a JSON list")
        # A record matched only by the search's organism filter carries no
        # per-record species evidence. The reader is asked to confirm it; the answer
        # is recorded with the run, so an audit can see it was confirmed rather than
        # assumed - which is what the documentation already promised.
        raw_pending = payload.get("species_confirmation_required_for") or 0
        if isinstance(raw_pending, bool) or not isinstance(raw_pending, (int, float)) or not math.isfinite(raw_pending):
            raise ValueError("species_confirmation_required_for must be a whole number")
        pending_species = int(raw_pending)
        if pending_species < 0 or pending_species > 10_000:
            raise ValueError("species_confirmation_required_for is out of range")
        species_confirmed = bool(payload.get("species_confirmed"))
        if pending_species and not species_confirmed:
            raise ValueError(
                f"{pending_species} selected record(s) were matched by the {spec.label} search filter "
                "rather than a per-record organism check; confirm their species before analysing"
            )
        run_id = secrets.token_hex(8)
        output = self.server.discovery_root / spec.key / "runs" / run_id
        result = run_discovery_analysis(
            prepared,
            selections,
            output,
            species=spec.key,
            min_studies=2,
            extra_metadata={
                "discovery_species_confirmed_by_reviewer": "true" if species_confirmed else "false",
                "discovery_species_records_needing_confirmation": str(pending_species),
            },
        )
        record = {"run_id": run_id, "bundle_id": bundle_id, **result}
        with self.server.discovery_lock:
            self.server.remember_discovery(self.server.discovery_runs, run_id, record)
        save_artifact = getattr(self.server.discovery_search_store, "save_artifact", None)
        if callable(save_artifact):
            save_artifact("run", run_id, record)
        return record

    def _load_discovery_run_summary(self, run_id: str) -> dict[str, Any]:
        with self.server.discovery_lock:
            result = self.server.discovery_runs.get(run_id)
        if result is None:
            get_artifact = getattr(self.server.discovery_search_store, "get_artifact", None)
            if callable(get_artifact):
                result = get_artifact("run", run_id)
        if result is None:
            raise FileNotFoundError("discovery analysis run was not found")
        return result

    def _discovery_run_summary(self, run_id: str) -> dict[str, Any]:
        result = self._load_discovery_run_summary(run_id)
        if not _is_loopback_host(self.server.server_address[0]):
            return _redact_local_paths_for_network(result)
        return result

    def _discovery_run_db(self, run_id: str) -> Path:
        result = self._load_discovery_run_summary(run_id)
        path = Path(str(result.get("db_path") or "")).resolve()
        expected_root = (self.server.discovery_root / str(result["species"]["key"]) / "runs" / run_id).resolve()
        if not path.is_relative_to(expected_root) or not path.is_file():
            raise FileNotFoundError("discovery analysis database is not available")
        return path

    def _discovery_run_excel(self, run_id: str) -> Path:
        result = self._load_discovery_run_summary(run_id)
        species = str(result.get("species", {}).get("key") or "")
        if species not in {"human", "mouse"}:
            raise FileNotFoundError("discovery analysis Excel workbook is not available")
        expected_root = (self.server.discovery_root / species / "runs" / run_id).resolve()
        output_dir = Path(str(result.get("output_dir") or "")).resolve()
        candidate = expected_root / "results" / "DEGORA_output.xlsx"
        resolved = candidate.resolve()
        if (
            output_dir != expected_root
            or resolved != candidate
            or not resolved.is_relative_to(expected_root)
            or not resolved.is_file()
        ):
            raise FileNotFoundError("discovery analysis Excel workbook is not available")
        return resolved

    def _genes(self, params: dict[str, list[str]], *, db_path: Path | None = None) -> dict[str, Any]:
        limit = _int_param(params, "limit", 50, minimum=1, maximum=500)
        offset = _int_param(params, "offset", 0, minimum=0, maximum=1_000_000)
        min_units = _int_param(params, "min_units", 1, minimum=1, maximum=10_000)
        min_score = _float_param(params, "min_score", 0.0, minimum=0.0)
        direction = _text_param(params, "direction", maximum=16).lower()
        query = _text_param(params, "q", maximum=128).upper()

        with closing(_connect(db_path or self.server.db_path)) as connection:
            columns = _table_columns(connection, "genes")
            score_column = _column_or_fallback(PRIMARY_SCORE_COLUMN, "degora_score", columns)
            direction_column = _column_or_fallback(PRIMARY_DIRECTION_COLUMN, "consensus_direction", columns)
            order_clause, sort_key, sort_order = _gene_order_clause(params, columns)

            where = ["n_source_units >= ?", f"{score_column} >= ?"]
            values: list[Any] = [min_units, min_score]
            if direction:
                if direction not in {"up", "down", "flat"}:
                    raise ValueError("direction must be up, down, or flat")
                where.append(f"{direction_column} = ?")
                values.append(direction)
            if query:
                # A search for SEPT9 has to find the gene DEGORA actually scored it
                # as (SEPTIN9). Both forms are matched so partial queries such as
                # "SEPT" keep working unchanged.
                terms = [query]
                resolved = canonical_gene_symbol(query)
                if resolved and resolved != query:
                    terms.append(resolved)
                where.append("(" + " OR ".join("gene_symbol LIKE ? ESCAPE '\\'" for _ in terms) + ")")
                values.extend(f"%{_escape_like_literal(term)}%" for term in terms)

            where_clause = " AND ".join(where)
            count = connection.execute(f"SELECT COUNT(*) FROM genes WHERE {where_clause}", values).fetchone()[0]
            rows = _row_dicts(
                connection.execute(
                    f"SELECT * FROM genes WHERE {where_clause} ORDER BY {order_clause} LIMIT ? OFFSET ?",
                    values + [limit, offset],
                )
            )
            rows = [_normalize_gene_api_record(row) for row in rows]
        return {"count": count, "limit": limit, "offset": offset, "sort": sort_key, "order": sort_order, "genes": rows}

    def _gene_detail(self, symbol: str, *, db_path: Path | None = None) -> dict[str, Any]:
        if not symbol:
            raise ValueError("gene symbol is required")
        if len(symbol) > 128:
            raise ValueError("gene symbol is too long; maximum length is 128 characters")
        with closing(_connect(db_path or self.server.db_path)) as connection:
            matched = symbol
            gene = _one_row(connection.execute("SELECT * FROM genes WHERE gene_symbol = ?", [symbol]))
            if gene is None:
                # A legacy or Excel-date symbol names the same gene as the current
                # symbol the run scored; resolve it rather than reporting a 404.
                resolved = canonical_gene_symbol(symbol)
                if resolved and resolved != symbol:
                    gene = _one_row(
                        connection.execute("SELECT * FROM genes WHERE gene_symbol = ?", [resolved])
                    )
                    if gene is not None:
                        matched = resolved
            if gene is None:
                raise FileNotFoundError(f"gene not found: {symbol[:64]}")
            gene = _normalize_gene_api_record(gene)
            evidence = _row_dicts(
                connection.execute(
                    "SELECT * FROM gene_evidence WHERE gene_symbol = ? ORDER BY source_unit_id, study_id",
                    [matched],
                )
            )
        if not _is_loopback_host(self.server.server_address[0]):
            evidence = _redact_records_paths_for_network(evidence)
        payload: dict[str, Any] = {"gene": gene, "evidence": evidence}
        if matched != symbol:
            payload["requested_gene_symbol"] = symbol
            payload["resolved_gene_symbol"] = matched
        return payload


_DISCOVERY_PROCESS_LOCK_GUARD = threading.Lock()
_DISCOVERY_PROCESS_LOCKS: set[str] = set()


class DiscoveryWorkspaceInUseError(RuntimeError):
    """Raised when another DEGORA server owns the same discovery workspace."""


class _DiscoveryRootProcessLock:
    """Hold one cross-process owner lock for a discovery workspace."""

    def __init__(self, path: Path) -> None:
        self.path = path.resolve()
        self._key = str(self.path)
        self._handle: Any | None = None
        self._claimed = False

    def acquire(self) -> None:
        with _DISCOVERY_PROCESS_LOCK_GUARD:
            if self._key in _DISCOVERY_PROCESS_LOCKS:
                raise DiscoveryWorkspaceInUseError("another DEGORA server is already using this discovery workspace")
            _DISCOVERY_PROCESS_LOCKS.add(self._key)
            self._claimed = True
        handle = None
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            handle = self.path.open("a+b")
            if os.name == "nt":  # pragma: no cover - exercised on native Windows only.
                import msvcrt

                handle.seek(0, os.SEEK_END)
                if handle.tell() == 0:
                    handle.write(b"\0")
                    handle.flush()
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except (ImportError, OSError) as exc:
            if handle is not None:
                handle.close()
            self._release_claim()
            raise DiscoveryWorkspaceInUseError("another DEGORA server is already using this discovery workspace") from exc
        self._handle = handle

    def _release_claim(self) -> None:
        # Only drop the shared key this instance actually claimed. A caller that
        # releases after a failed acquire would otherwise erase the entry belonging
        # to the server that already owns the workspace.
        with _DISCOVERY_PROCESS_LOCK_GUARD:
            if self._claimed:
                self._claimed = False
                _DISCOVERY_PROCESS_LOCKS.discard(self._key)

    def release(self) -> None:
        handle = self._handle
        self._handle = None
        if handle is not None:
            try:
                if os.name == "nt":  # pragma: no cover - exercised on native Windows only.
                    import msvcrt

                    handle.seek(0)
                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl

                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except (ImportError, OSError):
                pass
            finally:
                handle.close()
        self._release_claim()


class DegoraHttpServer(ThreadingHTTPServer):
    """HTTP server carrying the database path for request handlers."""

    # On Windows SO_REUSEADDR lets two listeners share a port, so a busy-port
    # probe would lie there and the flag stays off. On POSIX it only lets a new
    # listener bind while the previous run's connections sit in TIME_WAIT;
    # without it, Ctrl-C and restart moved the server to the next port for a
    # minute and blamed a second DEGORA that did not exist.
    allow_reuse_address = os.name != "nt"
    # Per-connection worker threads must not outlive the process, so a stuck
    # handler cannot keep the interpreter alive after Ctrl-C / shutdown.
    daemon_threads = True
    discovery_registry_limit = 100

    def __init__(
        self,
        server_address: tuple[str, int],
        db_path: str | Path,
        *,
        quiet: bool = False,
        access_token: str | None = None,
        discovery_root: str | Path | None = None,
    ) -> None:
        # ::1 is in LOOPBACK_HOSTS, so it has to bind as IPv6 rather than fail
        # inside socket.getaddrinfo with an address-family traceback.
        try:
            if ipaddress.ip_address(str(server_address[0]).strip()).version == 6:
                self.address_family = socket.AF_INET6
        except ValueError:
            pass
        # The discovery workspace is created before the socket is bound: a
        # permissions problem here is a filesystem error, and used to be swallowed
        # by the port-retry loop as "port is already in use".
        root = (
            Path(discovery_root).resolve()
            if discovery_root is not None
            else (Path(db_path).resolve().parent / "degora_discovery").resolve()
        )
        try:
            root.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise DiscoveryWorkspaceError(
                f"DEGORA cannot create its discovery workspace {root} ({exc.strerror or exc}). "
                "The results folder has to be writable by the user running `degora serve`; "
                "move or copy the results, or fix the folder permissions."
            ) from exc
        super().__init__(server_address, DegoraRequestHandler)
        self._degora_closed = False
        self._discovery_process_lock: _DiscoveryRootProcessLock | None = None
        self.discovery_search_store: Any | None = None
        self.discovery_job_manager: Any | None = None
        self.db_path = Path(db_path).resolve()
        self.quiet = quiet
        self.access_token = access_token
        self.discovery_root = (
            Path(discovery_root).resolve()
            if discovery_root is not None
            else (self.db_path.parent / "degora_discovery").resolve()
        )
        self.discovery_lock = threading.RLock()
        self.discovery_bundles: dict[str, dict[str, Any]] = {}
        self.discovery_runs: dict[str, dict[str, Any]] = {}
        self.discovery_root.mkdir(parents=True, exist_ok=True)
        self._discovery_process_lock = _DiscoveryRootProcessLock(self.discovery_root / ".degora-server.lock")
        try:
            self._discovery_process_lock.acquire()
            store_class, manager_class = _load_discovery_store_classes()
            self.discovery_search_store = store_class(self.discovery_root / "discovery.sqlite3")
            recover_jobs = getattr(self.discovery_search_store, "recover_interrupted_jobs", None)
            if callable(recover_jobs):
                recover_jobs()
            try:
                self.discovery_job_manager = manager_class(self.discovery_search_store, max_workers=2)
            except TypeError:
                self.discovery_job_manager = manager_class(self.discovery_search_store)
        except Exception:
            close_manager = getattr(self.discovery_job_manager, "close", None)
            if callable(close_manager):
                close_manager()
            shutdown_manager = getattr(self.discovery_job_manager, "shutdown", None)
            if callable(shutdown_manager):
                shutdown_manager()
            close_store = getattr(self.discovery_search_store, "close", None)
            if callable(close_store):
                close_store()
            self._discovery_process_lock.release()
            super().server_close()
            raise

    def remember_discovery(self, registry: dict[str, dict[str, Any]], key: str, value: dict[str, Any]) -> None:
        registry[key] = value
        while len(registry) > self.discovery_registry_limit:
            registry.pop(next(iter(registry)))

    def server_close(self) -> None:
        if getattr(self, "_degora_closed", False):
            return
        self._degora_closed = True
        try:
            manager = getattr(self, "discovery_job_manager", None)
            shutdown_manager = getattr(manager, "shutdown", None)
            if callable(shutdown_manager):
                try:
                    shutdown_manager(wait=False, cancel_futures=True, interrupt=True)
                except TypeError:
                    shutdown_manager()
            else:
                close_manager = getattr(manager, "close", None)
                if callable(close_manager):
                    try:
                        close_manager(wait=False, cancel_futures=True, interrupt=True)
                    except TypeError:
                        close_manager()
            close_store = getattr(getattr(self, "discovery_search_store", None), "close", None)
            if callable(close_store):
                close_store()
        finally:
            try:
                process_lock = getattr(self, "_discovery_process_lock", None)
                if process_lock is not None:
                    process_lock.release()
            finally:
                super().server_close()


MAX_PORT_ATTEMPTS = 20


def create_server(
    db_path: str | Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    *,
    quiet: bool = False,
    auto_port: bool = True,
    access_token: str | None = None,
    discovery_root: str | Path | None = None,
) -> DegoraHttpServer:
    """Bind the local server, auto-avoiding a port held by an unrelated service.

    If the requested port is busy, try the next few ports and finally an
    OS-assigned free port. A live DEGORA server that owns the same discovery
    workspace is rejected by the workspace lock even when another port is free.
    Pass ``port=0`` or ``auto_port=False`` to bind exactly one port.
    """

    if port == 0 or not auto_port:
        return DegoraHttpServer(
            (host, port),
            db_path,
            quiet=quiet,
            access_token=access_token,
            discovery_root=discovery_root,
        )

    candidates = [port + offset for offset in range(MAX_PORT_ATTEMPTS)] + [0]
    last_error: OSError | None = None
    for index, candidate in enumerate(candidates):
        try:
            server = DegoraHttpServer(
                (host, candidate),
                db_path,
                quiet=quiet,
                access_token=access_token,
                discovery_root=discovery_root,
            )
        except DiscoveryWorkspaceError:
            raise
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


def _is_loopback_host(value: str) -> bool:
    host = str(value or "").strip().lower().rstrip(".")
    if host in LOOPBACK_HOSTS:
        return True
    try:
        address = ipaddress.ip_address(host)
    except ValueError:
        return False
    mapped = getattr(address, "ipv4_mapped", None)
    return bool(address.is_loopback or (mapped is not None and mapped.is_loopback))


class ScoreDatabaseError(ValueError):
    """The path given to `degora serve` is not a DEGORA score database."""


class DiscoveryWorkspaceError(ValueError):
    """The discovery workspace beside the database cannot be created."""


REQUIRED_SCORE_DB_TABLES = ("genes", "gene_evidence", "studies", "meta")
# The columns the dashboard's first requests read, per table. /api/health counts
# studies and /api/genes orders by the primary rank, so a database that has the
# tables but not these columns still fails behind a page that loaded normally.
REQUIRED_SCORE_DB_COLUMNS = {
    "genes": ("gene_symbol", "degora_rank", "degora_score"),
    "gene_evidence": ("gene_symbol", "source_unit_id", "study_id"),
    "studies": ("study_id", "source_unit_id"),
    "meta": ("key", "value"),
}


def _require_degora_score_database(db_path: Path) -> None:
    """Refuse to serve a file that is not a DEGORA score database.

    Existence alone was checked before, so a mistyped path to a DEG table started a
    server that returned the dashboard with HTTP 200 and then failed every single
    API call with a 500. That is the hardest failure to diagnose, because the page
    loads. Fail here instead, before binding, with the command that makes the file.
    """

    fix = (
        "Pass the degora_scores.db written by `degora run` "
        "(for example: degora serve outputs/results/degora-run/degora_scores.db). "
        "No run yet? `degora demo my_demo` then `degora run my_demo/degora_demo_config.xlsx` "
        "writes one in a minute, and `degora launch <config>` runs and opens the browser in one step."
    )
    if db_path.is_dir():
        raise ScoreDatabaseError(f"DEGORA database path is a directory, not a database file: {db_path}\n{fix}")
    if not db_path.exists():
        # The first thing a new reader types is `degora serve` with no run behind
        # it; the answer has to say where a database comes from, not only that
        # this path has none.
        raise ScoreDatabaseError(f"DEGORA database does not exist: {db_path}\n{fix}")
    # The same read-only connection every request uses. Building the URI by hand
    # from the raw path let SQLite read `#`, `?` and `%XX` as URI syntax: it opened
    # a *different* file, read-write, created it when it did not exist, and then
    # reported the reader's perfectly good database as not a DEGORA database.
    try:
        with closing(_connect(db_path)) as connection:
            names = {
                str(row[0])
                for row in connection.execute("SELECT name FROM sqlite_master WHERE type IN ('table', 'view')")
            }
            missing = [table for table in REQUIRED_SCORE_DB_TABLES if table not in names]
            if missing:
                raise ScoreDatabaseError(
                    f"{db_path} is a SQLite file but not a DEGORA score database; it is missing "
                    f"{', '.join(missing)}.\n{fix}"
                )
            missing_columns: list[str] = []
            for table, columns in REQUIRED_SCORE_DB_COLUMNS.items():
                present = {str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})")}
                missing_columns.extend(f"{table}.{column}" for column in columns if column not in present)
            if missing_columns:
                raise ScoreDatabaseError(
                    f"{db_path} is a DEGORA score database with missing column(s): "
                    f"{', '.join(missing_columns)}. It was probably written by a different DEGORA version "
                    f"or edited by hand; rerun `degora run` to rebuild it.\n{fix}"
                )
    except sqlite3.DatabaseError as exc:
        raise ScoreDatabaseError(
            f"{db_path} is not a DEGORA score database (SQLite could not read it: {exc}).\n{fix}"
        ) from exc


def serve(
    db_path: str | Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    *,
    quiet: bool = False,
    allow_network: bool = False,
    access_token: str | None = None,
) -> None:
    db_path = Path(db_path)
    if not db_path.exists():
        # Fail before binding so a missing DB never serves (and never leaks its path over HTTP).
        # The preflight's message says where a database comes from; a bare "does
        # not exist" here reached the reader first and said nothing about that.
        _require_degora_score_database(db_path)
        raise FileNotFoundError(f"DEGORA database does not exist: {db_path}")
    _require_degora_score_database(db_path)
    token = access_token
    if not _is_loopback_host(host):
        if not allow_network:
            raise PermissionError(
                f"Refusing to serve DEGORA on non-loopback host {host!r} without --allow-network. "
                "Use the default 127.0.0.1 for local use, or pass --allow-network to expose the "
                "read-only browser/API with token protection."
            )
        token = token or secrets.token_urlsafe(24)
        print(
            f"WARNING: serving on {host} exposes this read-only DEGORA database browser to your "
            "network. A per-run access token is required; keep the printed URL private.",
            file=sys.stderr,
        )
    server = create_server(db_path, host=host, port=port, quiet=quiet, access_token=token)
    address, bound_port = server.server_address[:2]
    host_text = f"[{address}]" if ":" in str(address) else str(address)
    url = f"http://{host_text}:{bound_port}"
    if token:
        url = f"{url}#token={quote(token, safe='')}"
    print(f"DEGORA browser/API: {url}", flush=True)
    print(f"DEGORA version: {format_version_info()}", flush=True)
    print(f"Database: {server.db_path}", flush=True)
    try:
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nStopped DEGORA browser/API.", flush=True)
    finally:
        server.server_close()
