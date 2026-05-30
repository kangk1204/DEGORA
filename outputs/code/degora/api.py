"""Local HTTP API and browser UI for a DEGORA SQLite score database."""

from __future__ import annotations

import json
import math
import sqlite3
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
    body {
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 14px;
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
      grid-template-columns: minmax(0, 1.35fr) minmax(360px, 0.85fr);
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
    tbody tr { cursor: pointer; }
    tbody tr:hover { background: #eef7f5; }
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
    @media (max-width: 980px) {
      header { align-items: flex-start; flex-direction: column; }
      .meta { justify-content: flex-start; }
      main { grid-template-columns: 1fr; padding: 12px; }
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
  <main>
    <section>
      <div class="section-head">
        <h2>Genes</h2>
        <div class="status" id="status"></div>
      </div>
      <div class="controls">
        <input id="query" placeholder="Gene symbol" autocomplete="off">
        <input id="minUnits" type="number" min="1" value="1" aria-label="Min source units">
        <select id="direction" aria-label="Direction">
          <option value="">All directions</option>
          <option value="up">Up</option>
          <option value="down">Down</option>
          <option value="flat">Flat</option>
        </select>
        <button id="load">Load</button>
      </div>
      <table>
        <thead>
          <tr>
            <th style="width:58px">Rank</th>
            <th style="width:54px">Tier</th>
            <th>Gene</th>
            <th class="num">Top</th>
            <th class="num">Score</th>
            <th class="num">Units</th>
            <th class="num">Sign</th>
            <th class="num">LFC</th>
          </tr>
        </thead>
        <tbody id="genes"></tbody>
      </table>
    </section>
    <section>
      <div class="section-head">
        <h2>Evidence</h2>
      </div>
      <div class="detail-body" id="detail">
        <div class="empty">Select a gene.</div>
      </div>
    </section>
  </main>
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

    async function getJson(path) {
      const response = await fetch(path);
      if (!response.ok) throw new Error(await response.text());
      return response.json();
    }

    async function loadMeta() {
      const health = await getJson("/api/health");
      $("meta").innerHTML = [
        `${health.gene_count.toLocaleString()} genes`,
        `${health.study_count.toLocaleString()} studies`,
        `${health.source_unit_count.toLocaleString()} source units`
      ].map((text) => `<span>${text}</span>`).join("");
    }

    async function loadGenes() {
      const params = new URLSearchParams();
      params.set("limit", "100");
      const q = $("query").value.trim();
      const minUnits = $("minUnits").value.trim();
      const direction = $("direction").value;
      if (q) params.set("q", q);
      if (minUnits) params.set("min_units", minUnits);
      if (direction) params.set("direction", direction);
      const data = await getJson(`/api/genes?${params.toString()}`);
      $("status").textContent = `${data.count.toLocaleString()} shown`;
      $("genes").innerHTML = data.genes.map((gene) => `
        <tr data-gene="${esc(gene.gene_symbol)}">
          <td class="num">${gene.degora_rank}</td>
          <td>${tier(gene.evidence_tier)}</td>
          <td class="gene">${esc(gene.gene_symbol)}</td>
          <td class="num">${esc(gene.top_percent_label)}</td>
          <td class="num">${fmt(gene.degora_score, 2)}</td>
          <td class="num">${gene.n_source_units}</td>
          <td class="num">${fmt(gene.sign_concordance * 100, 1)}%</td>
          <td class="num">${fmt(gene.weighted_lfc, 2)}</td>
        </tr>
      `).join("");
      document.querySelectorAll("#genes tr").forEach((row) => {
        row.addEventListener("click", () => loadGene(row.dataset.gene));
      });
      if (data.genes.length) loadGene(data.genes[0].gene_symbol);
      else $("detail").innerHTML = `<div class="empty">No genes.</div>`;
    }

    async function loadGene(symbol) {
      const data = await getJson(`/api/genes/${encodeURIComponent(symbol)}`);
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
          <div class="metric"><span>Rank</span><strong>${esc(gene.rank_label)}</strong></div>
          <div class="metric"><span>Top fraction</span><strong>${esc(gene.top_percent_label)}</strong></div>
          <div class="metric"><span>Evidence tier</span><strong>${tier(gene.evidence_tier)}</strong></div>
          <div class="metric"><span>Score</span><strong>${fmt(gene.degora_score, 2)}</strong></div>
          <div class="metric"><span>Source support</span><strong>${esc(gene.support_label)}</strong></div>
          <div class="metric"><span>Direction</span><strong>${esc(gene.direction_label)}</strong></div>
          <div class="metric"><span>Evidence</span><strong>${fmt(gene.evidence_score, 2)}</strong></div>
          <div class="metric"><span>Rank signal</span><strong>${fmt(gene.rank_score_component, 2)}</strong></div>
          <div class="metric"><span>Weighted LFC</span><strong>${fmt(gene.weighted_lfc, 2)}</strong></div>
          <div class="metric"><span>Priority</span><strong>${fmt(gene.priority_score, 2)}</strong></div>
          <div class="metric"><span>Reliability</span><strong>${fmt(gene.evidence_reliability_score, 2)}</strong></div>
          <div class="metric"><span>Direction confidence</span><strong>${fmt(gene.direction_confidence_index * 100, 1)}%</strong></div>
          <div class="metric"><span>LOO stability</span><strong>${fmt(gene.loo_rank_stability_score * 100, 1)}%</strong></div>
        </div>
        <p class="sources">${badge(gene.consensus_direction)} DEGORA score is a relative prioritization score, not a probability.</p>
        <p class="sources">${esc(gene.source_units || "")}</p>
        <table>
          <thead>
            <tr>
              <th>Study</th>
              <th>Unit</th>
              <th>Pipeline</th>
              <th>Assay</th>
              <th>Time mode</th>
              <th class="num">Rel</th>
              <th class="num">LFC</th>
              <th class="num">z</th>
              <th class="num">Rank</th>
            </tr>
          </thead>
          <tbody>${evidenceRows}</tbody>
        </table>
      `;
    }

    $("load").addEventListener("click", loadGenes);
    $("query").addEventListener("keydown", (event) => {
      if (event.key === "Enter") loadGenes();
    });
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
        raise FileNotFoundError(f"DEGORA database does not exist: {db_path}")
    uri = f"file:{quote(str(db_path.resolve()), safe='/')}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    return connection


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
        value = default
    return min(max(value, minimum), maximum)


def _float_param(params: dict[str, list[str]], name: str, default: float, *, minimum: float = 0.0) -> float:
    raw = params.get(name, [str(default)])[0]
    try:
        value = float(raw)
    except ValueError:
        value = default
    if not math.isfinite(value):
        return default
    return max(value, minimum)


def _text_param(params: dict[str, list[str]], name: str, default: str = "", *, maximum: int = 128) -> str:
    value = params.get(name, [default])[0].strip()
    if len(value) > maximum:
        raise ValueError(f"{name} is too long; maximum length is {maximum} characters")
    return value


def _escape_like_literal(value: str) -> str:
    return value.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


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
        with _connect(self.server.db_path) as connection:
            gene_count = connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
            study_count = connection.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
            source_unit_count = connection.execute("SELECT COUNT(DISTINCT source_unit_id) FROM studies").fetchone()[0]
            top_gene = connection.execute("SELECT gene_symbol FROM genes ORDER BY degora_rank LIMIT 1").fetchone()
        return {
            "status": "ok",
            "db_path": str(self.server.db_path),
            "gene_count": gene_count,
            "study_count": study_count,
            "source_unit_count": source_unit_count,
            "top_gene": top_gene[0] if top_gene else None,
        }

    def _meta(self) -> dict[str, str]:
        with _connect(self.server.db_path) as connection:
            rows = _row_dicts(connection.execute("SELECT key, value FROM meta ORDER BY key"))
        return {row["key"]: row["value"] for row in rows}

    def _studies(self) -> list[dict[str, Any]]:
        with _connect(self.server.db_path) as connection:
            return _row_dicts(connection.execute("SELECT * FROM studies ORDER BY source_unit_id, study_id"))

    def _genes(self, params: dict[str, list[str]]) -> dict[str, Any]:
        limit = _int_param(params, "limit", 50, minimum=1, maximum=500)
        offset = _int_param(params, "offset", 0, minimum=0, maximum=1_000_000)
        min_units = _int_param(params, "min_units", 1, minimum=1, maximum=10_000)
        min_score = _float_param(params, "min_score", 0.0, minimum=0.0)
        direction = _text_param(params, "direction", maximum=16).lower()
        query = _text_param(params, "q", maximum=128).upper()

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
        with _connect(self.server.db_path) as connection:
            count = connection.execute(f"SELECT COUNT(*) FROM genes WHERE {where_clause}", values).fetchone()[0]
            rows = _row_dicts(
                connection.execute(
                    f"SELECT * FROM genes WHERE {where_clause} ORDER BY degora_rank LIMIT ? OFFSET ?",
                    values + [limit, offset],
                )
            )
        return {"count": count, "limit": limit, "offset": offset, "genes": rows}

    def _gene_detail(self, symbol: str) -> dict[str, Any]:
        if not symbol:
            raise ValueError("gene symbol is required")
        with _connect(self.server.db_path) as connection:
            gene = _one_row(connection.execute("SELECT * FROM genes WHERE gene_symbol = ?", [symbol]))
            if gene is None:
                raise FileNotFoundError(f"gene not found: {symbol}")
            evidence = _row_dicts(
                connection.execute(
                    "SELECT * FROM gene_evidence WHERE gene_symbol = ? ORDER BY source_unit_id, study_id",
                    [symbol],
                )
            )
        return {"gene": gene, "evidence": evidence}


class DegoraHttpServer(ThreadingHTTPServer):
    """HTTP server carrying the database path for request handlers."""

    def __init__(self, server_address: tuple[str, int], db_path: str | Path, *, quiet: bool = False) -> None:
        super().__init__(server_address, DegoraRequestHandler)
        self.db_path = Path(db_path).resolve()
        self.quiet = quiet


def create_server(db_path: str | Path, host: str = "127.0.0.1", port: int = 8765, *, quiet: bool = False) -> DegoraHttpServer:
    return DegoraHttpServer((host, port), db_path, quiet=quiet)


def serve(db_path: str | Path, host: str = "127.0.0.1", port: int = 8765, *, quiet: bool = False) -> None:
    server = create_server(db_path, host=host, port=port, quiet=quiet)
    address, bound_port = server.server_address
    print(f"DEGORA browser/API: http://{address}:{bound_port}")
    print(f"Database: {server.db_path}")
    try:
        server.serve_forever()
    finally:
        server.server_close()
