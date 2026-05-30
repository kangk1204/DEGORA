"""Tests for the static DEGORA Evidence Explorer generator.

These exercise the lean-DB builder and the self-contained HTML round-trip on a
small synthetic corpus, so they run fast and need no committed artifacts.
"""

from __future__ import annotations

import base64
import gzip
import re
import sqlite3
from pathlib import Path

import pytest

from figures import make_degora_evidence_explorer as explorer


def _make_corpus_db(path: Path, corpus_id: str) -> None:
    conn = sqlite3.connect(path)
    conn.execute(
        "create table genes (gene_symbol text, quality_weighted_degora_rank int, "
        "quality_weighted_degora_score real, consensus_direction text, "
        "quality_weighted_consensus_direction text, n_source_units int, "
        "sign_concordance real, evidence_tier text, heterogeneity_i2 real, "
        "heterogeneity_flag text, effect_meta_log2fc_re real, effect_meta_ci_low real, "
        "effect_meta_ci_high real, weighted_lfc real)"
    )
    conn.execute(
        "create table gene_evidence (gene_symbol text, source_unit_id text, paper_id text, "
        "pipeline text, assay_type text, platform text, cell_system text, species text, "
        "n_ctrl int, n_treat int, lfc real, signed_z real, normalized_rank real, "
        "direction text, source_url text, source_quality_weight real, source_quality_label text)"
    )
    conn.execute("create table studies (study_id text, source_unit_id text)")
    conn.executemany(
        "insert into genes (gene_symbol, quality_weighted_degora_rank, quality_weighted_degora_score, "
        "consensus_direction, quality_weighted_consensus_direction, n_source_units, sign_concordance, "
        "evidence_tier, heterogeneity_i2, heterogeneity_flag, effect_meta_log2fc_re, effect_meta_ci_low, "
        "effect_meta_ci_high, weighted_lfc) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [
            ("VEGFA", 1, 98.0, "up", "up", 3, 1.0, "A", 0.2, "low", 2.1, 1.5, 2.7, 2.0),
            ("EGLN3", 2, 95.0, "up", "up", 2, 1.0, "B", 0.1, "low", 1.8, 1.0, 2.6, 1.7),
        ],
    )
    conn.executemany(
        "insert into gene_evidence (gene_symbol, source_unit_id, paper_id, pipeline, assay_type, "
        "platform, cell_system, species, n_ctrl, n_treat, lfc, signed_z, normalized_rank, direction, "
        "source_url, source_quality_weight, source_quality_label) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [
            ("VEGFA", "GSE1", "GSE1", "DESeq2", "RNA-seq", "", "HeLa", "Homo sapiens", 3, 3,
             2.0, 5.0, 0.01, "up", "https://example.org/GSE1", 1.0, "high"),
            ("VEGFA", "GSE2", "GSE2", "edgeR", "RNA-seq", "", "A549", "Homo sapiens", 3, 3,
             2.2, 4.5, 0.02, "up", "https://example.org/GSE2", 0.9, "high"),
        ],
    )
    conn.execute("insert into studies values ('S1','GSE1')")
    conn.execute("insert into studies values ('S2','GSE2')")
    conn.commit()
    conn.close()


def test_build_lean_db_and_html_roundtrip(tmp_path: Path) -> None:
    result_dir = tmp_path / "results" / "demo"
    result_dir.mkdir(parents=True)
    _make_corpus_db(result_dir / "degora_scores.db", "demo")
    corpus = explorer.CorpusInput(
        "demo", "Demo", "demo topic", "RNA-seq", result_dir, result_dir / "missing_summary.csv"
    )

    raw, stats = explorer._build_lean_db([corpus])
    assert stats["n_scores"] == 2
    assert stats["n_evidence"] == 2

    # The lean DB must be a valid SQLite carrying both tables and the provenance join.
    db_path = tmp_path / "lean.db"
    db_path.write_bytes(raw)
    conn = sqlite3.connect(db_path)
    assert conn.execute("select count(*) from gene_scores").fetchone()[0] == 2
    rows = conn.execute(
        "select e.source_url from gene_evidence e join gene_scores s "
        "on e.gene_symbol = s.gene_symbol where s.gene_symbol = 'VEGFA'"
    ).fetchall()
    assert any(r[0] and r[0].startswith("http") for r in rows)
    conn.close()


def test_render_html_is_self_contained(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    if not (explorer.VENDOR_DIR / "sql-wasm.wasm").exists():
        pytest.skip("vendored sql.js not present")
    result_dir = tmp_path / "results" / "demo"
    result_dir.mkdir(parents=True)
    _make_corpus_db(result_dir / "degora_scores.db", "demo")
    corpus = explorer.CorpusInput(
        "demo", "Demo", "demo topic", "RNA-seq", result_dir, result_dir / "missing_summary.csv"
    )
    raw, stats = explorer._build_lean_db([corpus])
    html = explorer._render_html(raw, stats)

    # No runtime network dependency, embeds the engine and the gzip-compressed DB.
    assert "initSqlJs" in html
    assert "__DB_GZ_B64__" in html and "__WASM_B64__" in html
    assert "https://cdn" not in html and "http://cdn" not in html

    # The embedded DB must decode and query (the offline guarantee).
    db_b64 = re.search(r'__DB_GZ_B64__="([A-Za-z0-9+/=]+)"', html).group(1)
    decoded = gzip.decompress(base64.b64decode(db_b64))
    db_path = tmp_path / "embedded.db"
    db_path.write_bytes(decoded)
    conn = sqlite3.connect(db_path)
    assert conn.execute("select count(*) from gene_evidence").fetchone()[0] == 2
    conn.close()
