#!/usr/bin/env python
"""Build the DEGORA Evidence Explorer: a fully static, offline, single-file
interactive evidence database.

The explorer is one self-contained HTML file. It embeds:
  * the sql.js (SQLite-WASM) engine (vendored offline, no CDN at runtime),
  * a gzip-compressed lean SQLite database (all scored genes plus per-source
    provenance evidence for the top-ranked genes of every corpus),
and runs real SQL in the browser. It works from ``file://`` with no server and
no network, so the archived Zenodo copy stays functional. This is the
provenance-first, auditable resource the manuscript describes: every gene can
be drilled down to the individual source units, pipelines, platforms, cell
systems, directions, and source URLs that support it.

Design choices (grounded in the SOTA-resource strategy review):
  * Self-contained single file -> opens from file:// / Zenodo unzip / Pages.
  * Real client-side SQL over the gene_scores x gene_evidence join -> the
    auditable drill-down that ranking-only methods lack.
  * Gene permalinks via the URL hash, CSV export per gene and per query, and a
    pointer to the full per-corpus SQLite databases for exhaustive evidence.
  * Reproducible: generated from the committed degora_scores.db files, with a
    .source/.provenance sidecar like every other artifact.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import sqlite3
import tempfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from degora.provenance import shell_command, write_source_sidecar


FIGURE_ID = "DEGORA_EVIDENCE_EXPLORER"
OUTPUT_STEM = "degora_evidence_explorer"
VENDOR_DIR = Path(__file__).resolve().parent / "vendor"
TOP_EVIDENCE_RANK = 1000  # per corpus: genes whose full provenance is embedded


@dataclass(frozen=True)
class CorpusInput:
    corpus_id: str
    label: str
    topic: str
    input_mode: str
    result_dir: Path
    comparator_summary: Path

    @property
    def score_db(self) -> Path:
        return self.result_dir / "degora_scores.db"


DEFAULT_CORPORA = [
    CorpusInput("ifn-pilot", "IFN", "interferon response", "RNA-seq",
                Path("outputs/results/ifn-pilot"),
                Path("outputs/results/ifn-pilot/ifn_gold_comparator_summary.csv")),
    CorpusInput("ifn-cross-platform", "IFN", "interferon response", "RNA-seq + microarray",
                Path("outputs/results/ifn-cross-platform"),
                Path("outputs/results/ifn-cross-platform/ifn_cross_platform_gold_comparator_summary.csv")),
    CorpusInput("er-stress-cross-platform", "ER stress", "ER stress / UPR", "RNA-seq + microarray",
                Path("outputs/results/er-stress-cross-platform"),
                Path("outputs/results/er-stress-cross-platform/er_stress_cross_platform_gold_comparator_summary.csv")),
    CorpusInput("heat-shock-benchmark", "Heat shock", "heat shock / HSF1", "RNA-seq",
                Path("outputs/results/heat-shock-benchmark"),
                Path("outputs/results/heat-shock-benchmark/heat_shock_gold_comparator_summary.csv")),
    CorpusInput("hypoxia-hif1-benchmark", "Hypoxia", "hypoxia / HIF1", "RNA-seq",
                Path("outputs/results/hypoxia-hif1-benchmark"),
                Path("outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv")),
    CorpusInput("hypoxia-cross-platform", "Hypoxia", "hypoxia / HIF1", "RNA-seq + microarray",
                Path("outputs/results/hypoxia-cross-platform"),
                Path("outputs/results/hypoxia-cross-platform/hypoxia_cross_platform_gold_comparator_summary.csv")),
]

SCORE_COLUMNS = [
    "gene_symbol", "quality_weighted_degora_rank", "quality_weighted_degora_score",
    "quality_weighted_top_percent", "consensus_direction", "quality_weighted_consensus_direction",
    "n_source_units", "n_contrasts_observed", "sign_concordance", "quality_weighted_sign_concordance",
    "evidence_tier", "evidence_reliability_score", "direction_posterior_mean",
    "heterogeneity_i2", "heterogeneity_flag", "re_stouffer_z", "rra_rho",
    "effect_meta_log2fc_re", "effect_meta_ci_low", "effect_meta_ci_high", "effect_meta_i2",
    "effect_meta_k", "loo_rank_stability_score", "weighted_lfc", "high_confidence",
]
EVIDENCE_COLUMNS = [
    "gene_symbol", "source_unit_id", "paper_id", "pipeline", "assay_type", "platform",
    "cell_system", "species", "n_ctrl", "n_treat", "lfc", "signed_z", "normalized_rank",
    "direction", "source_url", "source_quality_weight", "source_quality_label",
]
METHOD_LABELS = {
    "degora_quality_weighted_score": "DEGORA (quality-weighted)", "degora_deg_score": "DEGORA",
    "degora_slice": "DEGORA slice", "weighted_stouffer": "Weighted Stouffer",
    "unweighted_stouffer": "Unweighted Stouffer", "fisher": "Fisher", "awfisher": "AWFisher",
    "metarnaseq_fisher": "metaRNASeq Fisher", "metarnaseq_invnorm": "metaRNASeq invnorm",
    "metavolcanor": "MetaVolcanoR", "robustrankaggreg": "RobustRankAggreg",
    "rank_product_approx": "Rank product approx.", "sign_vote": "Sign vote",
}


def _existing(table: sqlite3.Connection, name: str, wanted: list[str]) -> list[str]:
    have = {row[1] for row in table.execute(f"pragma table_info({name})")}
    return [c for c in wanted if c in have]


def _build_lean_db(corpora: list[CorpusInput]) -> tuple[bytes, dict[str, Any]]:
    """Assemble a compact SQLite DB: all gene_scores plus top-rank evidence."""
    out = sqlite3.connect(":memory:")
    out.execute(
        "create table gene_scores(corpus text, corpus_label text, input_mode text, "
        + ", ".join(f'"{c}" ' for c in SCORE_COLUMNS) + ")"
    )
    out.execute(
        "create table gene_evidence(corpus text, " + ", ".join(f'"{c}" ' for c in EVIDENCE_COLUMNS) + ")"
    )
    out.execute("create table corpora(corpus text, corpus_label text, topic text, input_mode text, n_genes int, n_source_units int)")
    out.execute("create table method_recall(corpus text, corpus_label text, method_id text, method_label text, recall_at_10 real, recall_at_50 real, recall_at_100 real)")

    n_scores = n_evidence = 0
    corpora_rows: list[dict[str, Any]] = []
    import csv as _csv

    for corpus in corpora:
        src = sqlite3.connect(f"file:{corpus.score_db}?mode=ro", uri=True)
        sc = _existing(src, "genes", SCORE_COLUMNS)
        rank_col = "quality_weighted_degora_rank" if "quality_weighted_degora_rank" in sc else "degora_rank"
        for row in src.execute(f"select {','.join(sc)} from genes"):
            d = dict(zip(sc, row))
            out.execute(
                "insert into gene_scores values (" + ",".join(["?"] * (3 + len(SCORE_COLUMNS))) + ")",
                [corpus.corpus_id, corpus.label, corpus.input_mode] + [d.get(c) for c in SCORE_COLUMNS],
            )
            n_scores += 1
        keep = {g.upper() for (g,) in src.execute(f"select gene_symbol from genes where {rank_col} <= ?", (TOP_EVIDENCE_RANK,))}
        ec = _existing(src, "gene_evidence", EVIDENCE_COLUMNS)
        for row in src.execute(f"select {','.join(ec)} from gene_evidence"):
            d = dict(zip(ec, row))
            if str(d.get("gene_symbol", "")).upper() in keep:
                out.execute(
                    "insert into gene_evidence values (" + ",".join(["?"] * (1 + len(EVIDENCE_COLUMNS))) + ")",
                    [corpus.corpus_id] + [d.get(c) for c in EVIDENCE_COLUMNS],
                )
                n_evidence += 1
        n_units = src.execute("select count(distinct source_unit_id) from studies").fetchone()[0]
        out.execute(
            "insert into corpora values (?,?,?,?,?,?)",
            (corpus.corpus_id, corpus.label, corpus.topic, corpus.input_mode,
             len(keep | {g.upper() for (g,) in src.execute('select gene_symbol from genes')}), n_units),
        )
        corpora_rows.append({"corpus": corpus.corpus_id, "label": corpus.label, "input_mode": corpus.input_mode})
        src.close()

        if corpus.comparator_summary.exists():
            with corpus.comparator_summary.open(encoding="utf-8") as handle:
                for mrow in _csv.DictReader(handle):
                    if mrow.get("run_status") != "ok":
                        continue
                    mid = mrow.get("method_id", "")
                    def _f(key: str) -> float | None:
                        try:
                            return round(float(mrow.get(key, "") or "nan"), 4)
                        except ValueError:
                            return None
                    out.execute(
                        "insert into method_recall values (?,?,?,?,?,?,?)",
                        (corpus.corpus_id, corpus.label, mid, METHOD_LABELS.get(mid, mid),
                         _f("recall_at_10"), _f("recall_at_50"), _f("recall_at_100")),
                    )

    out.execute("create index ix_scores_sym on gene_scores(gene_symbol)")
    out.execute("create index ix_scores_corpus on gene_scores(corpus, quality_weighted_degora_rank)")
    out.execute("create index ix_ev_sym on gene_evidence(gene_symbol)")
    out.commit()

    tmp = Path(tempfile.mkstemp(suffix=".db")[1])
    disk = sqlite3.connect(tmp)
    out.backup(disk)
    disk.close()
    out.close()
    raw = tmp.read_bytes()
    tmp.unlink()
    stats = {
        "n_scores": n_scores, "n_evidence": n_evidence, "n_corpora": len(corpora),
        "top_evidence_rank": TOP_EVIDENCE_RANK, "raw_db_bytes": len(raw),
    }
    return raw, stats


def _render_html(db_raw: bytes, stats: dict[str, Any]) -> str:
    wasm_b64 = base64.b64encode((VENDOR_DIR / "sql-wasm.wasm").read_bytes()).decode("ascii")
    sqljs_loader = (VENDOR_DIR / "sql-wasm.js").read_text(encoding="utf-8")
    db_gz_b64 = base64.b64encode(gzip.compress(db_raw, 9)).decode("ascii")
    generated = datetime.now(UTC).isoformat()
    meta = {
        "generated_at": generated, "n_scores": stats["n_scores"], "n_evidence": stats["n_evidence"],
        "n_corpora": stats["n_corpora"], "top_evidence_rank": stats["top_evidence_rank"],
    }
    return _HTML_TEMPLATE.format(
        meta_json=json.dumps(meta),
        sqljs_loader=sqljs_loader,
        wasm_b64=wasm_b64,
        db_gz_b64=db_gz_b64,
        app_js=_APP_JS,
        css=_CSS,
    )


def build_explorer(*, corpora: list[CorpusInput], output_dir: Path, command: str) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    missing = [str(c.score_db) for c in corpora if not c.score_db.exists()]
    if missing:
        raise FileNotFoundError(f"missing per-corpus score DB: {missing}")
    if not (VENDOR_DIR / "sql-wasm.wasm").exists() or not (VENDOR_DIR / "sql-wasm.js").exists():
        raise FileNotFoundError(f"vendored sql.js missing under {VENDOR_DIR}")

    db_raw, stats = _build_lean_db(corpora)
    html = _render_html(db_raw, stats)

    html_path = output_dir / f"{OUTPUT_STEM}.html"
    db_path = output_dir / f"{OUTPUT_STEM}.db"
    manifest_path = output_dir / f"{OUTPUT_STEM}_manifest.json"
    html_path.write_text(html, encoding="utf-8")
    db_path.write_bytes(db_raw)  # also ship the lean DB uncompressed for power users

    manifest = {
        "figure_id": FIGURE_ID,
        "generated_at": datetime.now(UTC).isoformat(),
        "script": "outputs/code/figures/make_degora_evidence_explorer.py",
        "command": command,
        "outputs": [str(html_path), str(db_path), str(manifest_path)],
        "stats": stats,
        "design": {
            "engine": "sql.js (SQLite compiled to WebAssembly), vendored offline",
            "deployment": "single self-contained HTML; embeds gzip-compressed SQLite + wasm; "
                          "opens from file:// with no server or network",
            "tables": ["gene_scores (all scored genes)", "gene_evidence (per-source provenance for top-ranked genes)",
                       "corpora", "method_recall"],
            "interactions": ["client-side SQL gene search", "corpus/direction/tier filters",
                             "per-gene per-source provenance drill-down with source URLs",
                             "cross-corpus comparison", "recall heatmap", "gene permalinks via URL hash",
                             "CSV export per gene and per query"],
        },
        "known_limitations": [
            f"Embedded provenance covers the top {TOP_EVIDENCE_RANK} ranked genes per corpus; "
            "lower-ranked genes show score-level metrics and point to the full per-corpus degora_scores.db.",
            "DEGORA scores are relative prioritization indices, not calibrated probabilities.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    inputs = [c.score_db for c in corpora] + [c.comparator_summary for c in corpora if c.comparator_summary.exists()]
    inputs += [VENDOR_DIR / "sql-wasm.wasm", VENDOR_DIR / "sql-wasm.js"]
    for artifact in (html_path, db_path, manifest_path):
        write_source_sidecar(artifact, command, inputs=inputs,
                             metadata={"generator": "degora-evidence-explorer", "figure_id": FIGURE_ID})

    return {"figure_id": FIGURE_ID, "outputs": [str(html_path), str(db_path), str(manifest_path)], **stats,
            "html_bytes": html_path.stat().st_size}


_CSS = """
:root{color-scheme:light;--bg:#f6f8f6;--ink:#1b2227;--muted:#5b6770;--panel:#fff;--line:#dde3dc;
--green:#0f766e;--red:#b4234b;--blue:#265d97;--gold:#9a6a12;--chip:#e7f5f2}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;font-size:14px}
header{padding:16px 22px;background:var(--panel);border-bottom:1px solid var(--line);position:sticky;top:0;z-index:5}
.topline{display:flex;justify-content:space-between;align-items:flex-end;gap:16px;flex-wrap:wrap}
h1{margin:0;font-size:21px}.sub{color:var(--muted);margin-top:3px}
.meta{color:var(--muted);font-variant-numeric:tabular-nums;text-align:right;font-size:12px}
.summary{display:grid;grid-template-columns:repeat(5,minmax(120px,1fr));gap:10px;padding:12px 22px;border-bottom:1px solid var(--line);background:#fbfcfa}
.metric{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:9px 11px}
.metric span{display:block;color:var(--muted);font-size:11px}.metric strong{display:block;margin-top:3px;font-size:18px;font-variant-numeric:tabular-nums}
main{display:grid;grid-template-columns:minmax(0,1fr) minmax(430px,.92fr);gap:16px;padding:16px 22px 28px}
section{min-width:0;background:var(--panel);border:1px solid var(--line);border-radius:8px;overflow:hidden}
.head{display:flex;justify-content:space-between;align-items:center;gap:10px;padding:11px 14px;border-bottom:1px solid var(--line)}
h2{margin:0;font-size:15px}.controls{display:grid;grid-template-columns:minmax(150px,1fr) 140px 120px 110px;gap:8px;padding:12px 14px;border-bottom:1px solid var(--line);background:#fbfcfa}
input,select,button{height:34px;border:1px solid var(--line);border-radius:6px;background:#fff;color:var(--ink);font:inherit;padding:0 10px;width:100%;min-width:0}
button{background:var(--green);border-color:var(--green);color:#fff;font-weight:700;cursor:pointer}
button.ghost{background:#fff;color:var(--green)}
table{width:100%;border-collapse:collapse;table-layout:fixed}
th,td{padding:8px 9px;border-bottom:1px solid var(--line);text-align:left;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
th{color:var(--muted);font-size:12px;background:#fbfcfa;position:sticky;top:0}
.tbl-wrap{max-height:60vh;overflow:auto}
tbody tr{cursor:pointer}tbody tr:hover{background:#eef7f4}tbody tr.sel{background:#e0f0eb}
.num{text-align:right;font-variant-numeric:tabular-nums}.gene{font-weight:800}
.badge{display:inline-flex;min-width:42px;height:21px;padding:0 7px;align-items:center;justify-content:center;border-radius:999px;font-size:12px;font-weight:800;background:var(--chip);color:var(--green)}
.badge.down{background:#fde7ed;color:var(--red)}.badge.flat,.badge.zero{background:#f3efe2;color:var(--gold)}
.tier{display:inline-flex;width:24px;height:22px;align-items:center;justify-content:center;border-radius:6px;font-weight:800;background:#eef2ff;color:#3730a3}
.tier.B{background:#e7f5f2;color:#0f5d56}.tier.C{background:#f4f2e8;color:var(--gold)}.tier.D{background:#f1f5f9;color:var(--muted)}
.detail{padding:14px;max-height:78vh;overflow:auto}
.kv{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:9px;margin-bottom:12px}
.evcard{border:1px solid var(--line);border-radius:8px;padding:10px 12px;margin-bottom:9px;background:#fcfdfb}
.evcard .row1{display:flex;justify-content:space-between;gap:8px;align-items:baseline}
.evcard .unit{font-weight:700}.evcard .tags{color:var(--muted);font-size:12px;margin-top:3px;line-height:1.5}
.evcard a{color:var(--blue)}
.note{color:var(--muted);line-height:1.5;margin:8px 0 0}
.bar{height:8px;border-radius:4px;background:#e9efe9;overflow:hidden;margin-top:5px}
.bar>i{display:block;height:100%;background:var(--green)}
svg{width:100%;height:auto;display:block}
.toolbar{display:flex;gap:8px;flex-wrap:wrap;padding:0 14px 12px}
.toolbar button{width:auto;padding:0 12px}
.loading{padding:40px 22px;color:var(--muted);text-align:center}
@media(max-width:1040px){main{grid-template-columns:1fr}.summary{grid-template-columns:repeat(2,1fr)}}
@media(max-width:620px){.controls{grid-template-columns:1fr}.kv{grid-template-columns:1fr}}
"""

_APP_JS = r"""
const META = __META__;
let SQL, db;
const $ = (id)=>document.getElementById(id);
const esc=(v)=>String(v??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const cls=(v)=>String(v??"").replace(/[^A-Za-z0-9_-]/g,"");
const fmt=(v,d=2)=>{const n=Number(v);return (v===null||v===""||Number.isNaN(n))?"":n.toFixed(d);};
const pct=(v,d=1)=>{const n=Number(v);return Number.isFinite(n)?(n*100).toFixed(d)+"%":"";};
const badge=(v)=>`<span class="badge ${cls((v||"flat").toLowerCase())}">${esc(v||"flat")}</span>`;

function b64ToBytes(b64){const bin=atob(b64);const a=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)a[i]=bin.charCodeAt(i);return a;}
async function gunzip(bytes){
  if(typeof DecompressionStream!=="undefined"){
    const ds=new DecompressionStream("gzip");
    const out=new Response(new Blob([bytes]).stream().pipeThrough(ds));
    return new Uint8Array(await out.arrayBuffer());
  }
  throw new Error("This browser lacks DecompressionStream; please use a current browser or the bundled .db file.");
}
function q(sql,params=[]){const r=db.exec(sql,params);if(!r.length)return [];const {columns,values}=r[0];return values.map(row=>Object.fromEntries(row.map((v,i)=>[columns[i],v])));}

async function boot(){
  try{
    SQL = await initSqlJs({ wasmBinary: b64ToBytes(window.__WASM_B64__) });
    const dbBytes = await gunzip(b64ToBytes(window.__DB_GZ_B64__));
    db = new SQL.Database(dbBytes);
    delete window.__DB_GZ_B64__; delete window.__WASM_B64__;
    initUI();
  }catch(e){ document.body.innerHTML = `<div class="loading">Failed to load evidence database: ${esc(e.message)}</div>`; }
}

function initUI(){
  const corpora=q("select corpus,corpus_label,input_mode,n_source_units from corpora");
  $("app").style.display="";
  $("loading").style.display="none";
  $("generated").textContent=`Generated ${META.generated_at.slice(0,10)} · ${META.n_scores.toLocaleString()} scored genes · ${META.n_evidence.toLocaleString()} provenance rows · ${META.n_corpora} corpora`;
  $("corpus").innerHTML=`<option value="">All corpora</option>`+corpora.map(c=>`<option value="${esc(c.corpus)}">${esc(c.corpus_label)} (${esc(c.input_mode)})</option>`).join("");
  const totUnits=corpora.reduce((s,c)=>s+(c.n_source_units||0),0);
  $("summary").innerHTML=[
    ["Scored genes",META.n_scores.toLocaleString()],
    ["Provenance rows",META.n_evidence.toLocaleString()],
    ["Corpus configs",META.n_corpora],
    ["Source units",totUnits],
    ["Engine","SQLite-WASM"]
  ].map(([k,v])=>`<div class="metric"><span>${esc(k)}</span><strong>${esc(v)}</strong></div>`).join("");
  drawHeatmap();
  ["query","corpus","direction","minUnits"].forEach(id=>{
    $(id).addEventListener(id==="query"||id==="minUnits"?"input":"change",()=>{runSearch();syncHash();});
  });
  $("dlQuery").addEventListener("click",downloadQuery);
  window.addEventListener("hashchange",applyHash);
  if(!applyHash()) runSearch();
}

function filters(){
  const wh=[],pa=[];
  const qv=$("query").value.trim().toUpperCase();
  const cv=$("corpus").value, dv=$("direction").value, mu=parseInt($("minUnits").value||"1",10);
  if(qv){wh.push("upper(gene_symbol) like ?");pa.push("%"+qv+"%");}
  if(cv){wh.push("corpus = ?");pa.push(cv);}
  if(dv){wh.push("lower(coalesce(quality_weighted_consensus_direction,consensus_direction,'')) = ?");pa.push(dv);}
  if(mu>1){wh.push("n_source_units >= ?");pa.push(mu);}
  return [wh.length?("where "+wh.join(" and ")):"",pa];
}
function runSearch(){
  const [where,pa]=filters();
  const rows=q(`select corpus,corpus_label,input_mode,gene_symbol,quality_weighted_degora_rank as rank,
    coalesce(quality_weighted_consensus_direction,consensus_direction) as dir,quality_weighted_top_percent as top,
    n_source_units as units,evidence_tier as tier,evidence_reliability_score as rel,quality_weighted_degora_score as score
    from gene_scores ${where} order by (rank is null), rank limit 300`,pa);
  const total=q(`select count(*) n from gene_scores ${where}`,pa)[0].n;
  $("status").textContent=`${total.toLocaleString()} matches · ${rows.length} shown`;
  $("rows").innerHTML=rows.map((r,i)=>`<tr data-g="${esc(r.gene_symbol)}" data-c="${esc(r.corpus)}">
    <td>${esc(r.corpus_label)}</td><td class="num">${r.rank??""}</td><td class="gene">${esc(r.gene_symbol)}</td>
    <td>${badge(r.dir)}</td><td class="num">${fmt(r.top,3)}</td><td class="num">${r.units??""}</td>
    <td><span class="tier ${cls(r.tier)}">${esc(r.tier)}</span></td></tr>`).join("");
  [...document.querySelectorAll("#rows tr")].forEach(tr=>tr.addEventListener("click",()=>{
    document.querySelectorAll("#rows tr.sel").forEach(x=>x.classList.remove("sel"));tr.classList.add("sel");
    showGene(tr.dataset.g,tr.dataset.c);syncHash(tr.dataset.g,tr.dataset.c);
  }));
  if(rows.length){const f=document.querySelector("#rows tr");f.classList.add("sel");showGene(rows[0].gene_symbol,rows[0].corpus);}
  else $("detail").innerHTML=`<p class="note">No genes match the current filters.</p>`;
}

function showGene(sym,corpus){
  const all=q(`select * from gene_scores where gene_symbol = ? order by (quality_weighted_degora_rank is null), quality_weighted_degora_rank`,[sym]);
  if(!all.length){$("detail").innerHTML=`<p class="note">No record for ${esc(sym)}.</p>`;return;}
  const row=all.find(r=>r.corpus===corpus)||all[0];
  const m=(k,v,d=2)=>`<div class="metric"><span>${k}</span><strong>${v}</strong></div>`;
  const i2=Number(row.heterogeneity_i2);
  const head=`<div class="head" style="padding-left:0;border:0"><h2>${esc(sym)} · ${esc(row.corpus_label)} <span class="badge ${cls((row.quality_weighted_consensus_direction||row.consensus_direction||'flat').toLowerCase())}">${esc(row.quality_weighted_consensus_direction||row.consensus_direction||'flat')}</span></h2></div>`;
  const kv=`<div class="kv">
    ${m("DEGORA rank",row.quality_weighted_degora_rank??"")}
    ${m("Top %",fmt(row.quality_weighted_top_percent,3))}
    ${m("Evidence tier",`<span class="tier ${cls(row.evidence_tier)}">${esc(row.evidence_tier)}</span>`)}
    ${m("Source units",row.n_source_units??"")}
    ${m("Sign concordance",pct(row.quality_weighted_sign_concordance??row.sign_concordance))}
    ${m("Direction posterior",pct(row.direction_posterior_mean))}
    ${m("Heterogeneity I²",Number.isFinite(i2)?i2.toFixed(2):"")}
    ${m("RE log2FC (95% CI)",`${fmt(row.effect_meta_log2fc_re)} <span style="font-size:12px;color:var(--muted)">[${fmt(row.effect_meta_ci_low)}, ${fmt(row.effect_meta_ci_high)}]</span>`)}
    ${m("LOO stability",pct(row.loo_rank_stability_score))}
  </div>`;
  const het=(row.heterogeneity_flag||"").includes("high")?`<p class="note">⚠ Flagged <b>high_context_dependent_review</b> (I²≥0.75): inspect the per-source evidence before treating this as a context-stable marker.</p>`:"";
  const ev=q(`select * from gene_evidence where gene_symbol = ? and corpus = ? order by source_quality_weight desc, source_unit_id`,[sym,corpus]);
  let evHtml;
  if(ev.length){
    evHtml=ev.map(e=>{
      const url=e.source_url?`<a href="${esc(e.source_url)}" target="_blank" rel="noopener">source ↗</a>`:"";
      const reps=(e.n_ctrl!=null&&e.n_treat!=null)?`${e.n_ctrl}v${e.n_treat}`:"";
      return `<div class="evcard"><div class="row1"><span class="unit">${esc(e.paper_id||e.source_unit_id)}</span>${badge(e.direction)}</div>
        <div class="tags">${esc(e.assay_type||"")}${e.platform?" · "+esc(e.platform):""}${e.cell_system?" · "+esc(e.cell_system):""}${reps?" · n="+reps:""}<br>
        pipeline: ${esc(e.pipeline||"n/a")} · log2FC ${fmt(e.lfc)} · signed-z ${fmt(e.signed_z)} · rank ${fmt(e.normalized_rank,4)} · ${esc(e.source_quality_label||"")} ${url}</div></div>`;
    }).join("");
  } else {
    evHtml=`<p class="note">Per-source provenance for ${esc(sym)} is beyond the top ${META.top_evidence_rank} embedded rows for this corpus. The complete evidence is in the corpus-level <code>degora_scores.db</code> (gene_evidence table) shipped with this resource.</p>`;
  }
  const cross=all.length>1?`<table><thead><tr><th>Corpus</th><th class="num">Rank</th><th>Dir</th><th class="num">Top%</th><th class="num">Units</th></tr></thead><tbody>${all.map(r=>`<tr><td>${esc(r.corpus_label)} <span style="color:var(--muted);font-size:11px">${esc(r.input_mode)}</span></td><td class="num">${r.quality_weighted_degora_rank??""}</td><td>${badge(r.quality_weighted_consensus_direction||r.consensus_direction)}</td><td class="num">${fmt(r.quality_weighted_top_percent,3)}</td><td class="num">${r.n_source_units??""}</td></tr>`).join("")}</tbody></table>`:"";
  $("detail").innerHTML=`${head}${kv}${het}
    <div class="toolbar"><button class="ghost" onclick="dlGene('${esc(sym)}','${esc(corpus)}')">Download evidence (CSV)</button>
    <button class="ghost" onclick="copyLink()">Copy permalink</button></div>
    ${cross?`<h2 style="font-size:13px;margin:6px 0">Cross-corpus ranks</h2>${cross}`:""}
    <h2 style="font-size:13px;margin:14px 0 6px">Per-source evidence (${ev.length})</h2>${evHtml}
    <p class="note">Scores are relative prioritization indices, not calibrated probabilities.</p>`;
}

function toCSV(rows){if(!rows.length)return "";const cols=Object.keys(rows[0]);
  return [cols.join(",")].concat(rows.map(r=>cols.map(c=>{const v=r[c]==null?"":String(r[c]);return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v;}).join(","))).join("\n");}
function dl(name,text){const b=new Blob([text],{type:"text/csv"});const a=document.createElement("a");a.href=URL.createObjectURL(b);a.download=name;a.click();URL.revokeObjectURL(a.href);}
window.dlGene=(sym,corpus)=>{const ev=q("select * from gene_evidence where gene_symbol=? and corpus=?",[sym,corpus]);dl(`DEGORA_${sym}_${corpus}_evidence.csv`,toCSV(ev.length?ev:q("select * from gene_scores where gene_symbol=? and corpus=?",[sym,corpus])));};
function downloadQuery(){const [where,pa]=filters();dl("DEGORA_query.csv",toCSV(q(`select * from gene_scores ${where} order by (quality_weighted_degora_rank is null),quality_weighted_degora_rank limit 5000`,pa)));}
window.copyLink=()=>{navigator.clipboard&&navigator.clipboard.writeText(location.href);};

function syncHash(sym,corpus){
  const p=new URLSearchParams();
  if($("query").value)p.set("q",$("query").value);
  if($("corpus").value)p.set("corpus",$("corpus").value);
  if($("direction").value)p.set("dir",$("direction").value);
  if(sym)p.set("gene",sym);if(corpus)p.set("gc",corpus);
  history.replaceState(null,"","#"+p.toString());
}
function applyHash(){
  if(!location.hash||location.hash.length<2)return false;
  const p=new URLSearchParams(location.hash.slice(1));
  $("query").value=p.get("q")||"";$("corpus").value=p.get("corpus")||"";$("direction").value=p.get("dir")||"";
  runSearch();
  const g=p.get("gene");if(g)showGene(g,p.get("gc")||$("corpus").value||"ifn-pilot");
  return true;
}

function drawHeatmap(){
  const rows=q("select corpus_label,method_label,recall_at_100 from method_recall where method_label is not null");
  if(!rows.length){$("heatmap").innerHTML="";return;}
  const corpora=[...new Set(rows.map(r=>r.corpus_label))];
  const order=["DEGORA (quality-weighted)","DEGORA","Weighted Stouffer","Fisher","AWFisher","metaRNASeq Fisher","MetaVolcanoR","RobustRankAggreg","Rank product approx."];
  const methods=order.filter(m=>rows.some(r=>r.method_label===m));
  const cw=104,ch=26,left=178,top=44,W=left+corpora.length*cw+12,H=top+methods.length*ch+30;
  const col=(v)=>{v=Math.max(0,Math.min(1,Number(v)||0));return `rgb(${Math.round(248-210*v)},${Math.round(250-92*v)},${Math.round(248-126*v)})`;};
  let b="";
  corpora.forEach((c,j)=>b+=`<text x="${left+j*cw+cw/2}" y="24" text-anchor="middle" font-size="11">${esc(c)}</text>`);
  methods.forEach((m,i)=>{b+=`<text x="8" y="${top+i*ch+17}" font-size="11">${esc(m)}</text>`;
    corpora.forEach((c,j)=>{const r=rows.find(x=>x.corpus_label===c&&x.method_label===m);const v=r?Number(r.recall_at_100):null;
      b+=`<rect x="${left+j*cw}" y="${top+i*ch}" width="${cw-4}" height="${ch-4}" rx="4" fill="${v==null?'#f3f4f3':col(v)}"></rect>`;
      b+=`<text x="${left+j*cw+(cw-4)/2}" y="${top+i*ch+17}" text-anchor="middle" font-size="11" fill="#16202a">${v==null?'':v.toFixed(2)}</text>`;});});
  $("heatmap").innerHTML=`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Recall at 100 by method and corpus">${b}</svg>`;
}
boot();
"""

_HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DEGORA Evidence Explorer</title>
<style>{css}</style>
</head>
<body>
<header>
  <div class="topline">
    <div><h1>DEGORA Evidence Explorer</h1>
    <div class="sub">Auditable gene-evidence database over public DEG summaries — runs entirely in your browser, offline.</div></div>
    <div class="meta" id="generated"></div>
  </div>
</header>
<div id="loading" class="loading">Loading the embedded SQLite evidence database…</div>
<div id="app" style="display:none">
  <div class="summary" id="summary"></div>
  <main>
    <section>
      <div class="head"><h2>Gene search</h2><div class="meta" id="status"></div></div>
      <div class="controls">
        <input id="query" placeholder="Gene symbol, e.g. VEGFA, DDIT3, CMPK2" autocomplete="off" aria-label="Gene symbol">
        <select id="corpus" aria-label="Corpus"></select>
        <select id="direction" aria-label="Direction"><option value="">All directions</option><option value="up">Up</option><option value="down">Down</option><option value="flat">Flat</option></select>
        <input id="minUnits" type="number" min="1" value="1" aria-label="Minimum source units">
      </div>
      <div class="toolbar"><button id="dlQuery" class="ghost">Download current query (CSV)</button></div>
      <div class="tbl-wrap"><table><thead><tr>
        <th style="width:120px">Corpus</th><th style="width:64px" class="num">Rank</th><th>Gene</th>
        <th style="width:70px">Dir</th><th style="width:76px" class="num">Top %</th><th style="width:64px" class="num">Units</th><th style="width:52px">Tier</th>
      </tr></thead><tbody id="rows"></tbody></table></div>
    </section>
    <section>
      <div class="head"><h2>Evidence &amp; provenance</h2><div class="meta">click a gene</div></div>
      <div id="heatmap" style="padding:12px 14px;border-bottom:1px solid var(--line)"></div>
      <div class="detail" id="detail"></div>
    </section>
  </main>
</div>
<script>{sqljs_loader}</script>
<script>window.__WASM_B64__="{wasm_b64}";window.__DB_GZ_B64__="{db_gz_b64}";</script>
<script>const __META__={meta_json};{app_js}</script>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/figures/manuscript/DEGORA_EVIDENCE_EXPLORER"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = shell_command(
        ["PYTHONPATH=outputs/code", "python", "outputs/code/figures/make_degora_evidence_explorer.py",
         "--output-dir", args.output_dir]
    )
    summary = build_explorer(corpora=DEFAULT_CORPORA, output_dir=args.output_dir, command=command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
