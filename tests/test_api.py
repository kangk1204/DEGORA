from __future__ import annotations

import json
import http.client
import sqlite3
from pathlib import Path
import threading
import time
import urllib.request
from urllib.error import HTTPError

import pandas as pd
import pytest

from degora import __version__
import degora.api as api
from degora.api import (
    ScoreDatabaseError,
    LOOPBACK_HOSTS,
    DiscoveryWorkspaceInUseError,
    create_server,
    redact_token_from_url_text,
    serve,
    strip_token_query_param,
)
from degora.provenance import EXTERNAL_PATH_PREFIX
from degora.score_db import write_score_database


def test_discovery_job_manager_interrupts_running_search(tmp_path) -> None:
    """A stopped server must not block on an in-flight search, and must persist why.

    This drives the same store/manager classes the running server uses, so the
    interruption contract is exercised on the live path rather than on a stand-in.
    """

    store_class, manager_class = api._load_discovery_store_classes()
    store = store_class(tmp_path / "discovery.sqlite3")
    manager = manager_class(store, max_workers=1)

    started = threading.Event()
    release = threading.Event()
    search_id = "a1b2c3d4e5f60718"
    store.save_search(search_id, {"id": search_id, "query": "hypoxia", "status": "queued"})

    def worker(_job_id, _payload, _progress):
        started.set()
        release.wait(timeout=5)
        return {"search_id": search_id}

    job = manager.submit("publication_search", {"search_id": search_id}, worker)
    assert started.wait(timeout=2)

    before = time.monotonic()
    manager.shutdown(wait=False, cancel_futures=True, interrupt=True)
    assert time.monotonic() - before < 0.5
    assert store.get_job(job["job_id"])["status"] == "interrupted"

    release.set()
    manager.shutdown(wait=True)
    assert store.get_job(job["job_id"])["status"] == "interrupted"


def _harmonized() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S1", "S2"],
            "paper_id": ["P1", "P2", "P1", "P2"],
            "gene_symbol": ["VEGFA", "VEGFA", "RPL13A", "RPL13A"],
            "lfc": [2.0, 1.8, 0.1, -0.1],
            "signed_z": [5.0, 4.5, 0.1, -0.1],
            "pvalue": [1e-7, 1e-6, 0.9, 0.8],
            "padj": [1e-5, 1e-4, 0.9, 0.9],
            "normalized_rank": [0.02, 0.03, 0.9, 0.8],
            "n_ctrl": [3, 4, 3, 4],
            "n_treat": [3, 4, 3, 4],
            "n_genes_in_study": [1000] * 4,
            "pipeline": ["DESeq2", "edgeR", "DESeq2", "edgeR"],
            "assay_type": ["RNA-seq", "microarray", "RNA-seq", "microarray"],
            "source_input_type": ["author_deg_table", "limma_full_table", "author_deg_table", "limma_full_table"],
            "platform": ["", "GPL570", "", "GPL570"],
            "normalization": ["", "RMA/log2", "", "RMA/log2"],
            "probe_collapse": ["", "author_gene_level", "", "author_gene_level"],
            "species": ["Homo sapiens"] * 4,
            "cell_system": ["A", "B", "A", "B"],
            "hypoxia_modality": ["1% O2"] * 4,
            "duration_h": ["24"] * 4,
            "source_path": ["source.csv"] * 4,
            "source_url": ["https://example.test"] * 4,
        }
    )


def _get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=5) as response:
        return json.loads(response.read().decode("utf-8"))


def test_local_api_serves_health_gene_list_and_detail(tmp_path, monkeypatch) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")
    monkeypatch.setattr(
        api,
        "runtime_version_info",
        lambda: {
            "degora_version": __version__,
            "degora_code_revision": "abc1234-dirty",
            "degora_code_dirty": "true",
        },
    )

    server = create_server(tmp_path / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"

    try:
        with urllib.request.urlopen(base_url, timeout=5) as response:
            html = response.read().decode("utf-8")
        with urllib.request.urlopen(f"{base_url}/api/health", timeout=5) as response:
            assert response.headers.get("Access-Control-Allow-Origin") is None
        health = _get_json(f"{base_url}/api/health")
        genes = _get_json(f"{base_url}/api/genes?q=VEGF&limit=5")
        gene_asc = _get_json(f"{base_url}/api/genes?sort=gene&order=asc&limit=5")
        score_desc = _get_json(f"{base_url}/api/genes?sort=score&order=desc&limit=5")
        wildcard = _get_json(f"{base_url}/api/genes?q=_&limit=5")
        detail = _get_json(f"{base_url}/api/genes/VEGFA")
        try:
            _get_json(f"{base_url}/api/genes?q={'A' * 129}")
        except HTTPError as exc:
            assert exc.code == 400
        else:  # pragma: no cover - assertion branch
            raise AssertionError("long gene query should be rejected")
        for query in ("sort=source_url", "order=sideways"):
            try:
                _get_json(f"{base_url}/api/genes?{query}")
            except HTTPError as exc:
                assert exc.code == 400
            else:  # pragma: no cover - assertion branch
                raise AssertionError(f"{query} should be rejected")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert "<title>DEGORA</title>" in html
    assert '<meta name="viewport" content="width=device-width, initial-scale=1">' in html
    assert "const esc =" in html
    assert "const url = new URL(text);" in html
    assert "new URL(String(value || \"\"), window.location.href)" not in html
    assert 'id="layoutSplitter"' in html
    assert 'role="separator"' in html
    assert "radial-gradient" in html
    assert "18px minmax" in html
    assert "initPanelResize();" in html
    assert "JSON.parse(message)" in html
    assert 'data-species="human"' in html
    assert 'data-species="mouse"' in html
    assert 'aria-pressed="true"' in html
    assert "10 per page · up to 1,000 ranked" in html
    assert "MAX_SELECTED_STUDIES = 20" in html
    assert "Find Human or Mouse studies" in html
    assert ".discovery-search input { width: 100%; min-width: 0; }" in html
    assert "grid-template-columns: minmax(0, 1fr);" in html
    assert 'required minlength="2" maxlength="200"' in html
    assert "data-study-sort" in html
    assert 'aria-sort="${studyAriaSort' in html
    assert '${state.evaluatedStudies.toLocaleString()} assessed' in html
    assert "Likely author DEG table" not in html  # labels arrive from assessed search JSON
    assert 'study.deg_input_assessment' in html
    assert "Human and Mouse never pooled" in html
    assert 'id="discoveryActions" hidden' in html
    assert 'id="clearSelected"' in html
    assert 'id="resetStudySort"' in html
    assert 'state.sort = { key: "readiness", order: "desc" };' in html
    assert "Prepare selection" in html
    assert "Enter exact paper/GEO contrast" in html
    assert 'class="n-ctrl" type="number" min="1" step="1"' in html
    assert 'class="n-treat" type="number" min="1" step="1"' in html
    assert 'class="cell-system"' in html
    assert 'class="duration-h"' in html
    assert 'class="platform"' in html
    assert "function isPositiveWholeNumber(value)" in html
    assert 'common.n_ctrl = Number(row.querySelector(".n-ctrl")?.value || 0);' in html
    assert 'common.n_treat = Number(row.querySelector(".n-treat")?.value || 0);' in html
    assert "positive whole-number biological group sizes" in html
    assert "Choose confirmed scale" in html
    assert "biological-replicates-confirmed" in html
    # The assertion still has to be made; it is phrased for a reader who does not
    # already know what "independent biological replicate" rules out.
    assert "each column is a separate biological sample" in html
    assert "Enter exact contrasts, mappings, positive whole-number biological group sizes" in html
    assert "fallback matrices also require scale, biological-replicate attestation, and 2 + 2 sample assignment" in html
    assert "data-source-unit" in html
    assert "row.dataset.sourceUnit || row.dataset.accession" in html
    assert ".study-table th:nth-child(4)" in html
    assert '"X-DEGORA-Action": "1"' in html
    assert '<button id="discoverNav" class="active"' in html
    assert '<button id="atlasNav" type="button" aria-pressed="false"' in html
    assert '<div class="discovery-view" id="discoveryView">' in html
    assert '<main id="layoutMain" hidden>' in html
    assert 'showView("discover");' in html
    assert '<button id="discoverySearch" type="button">Search</button>' in html
    assert "Run separate Human + Mouse searches" in html
    assert 'data-study-inspect="${esc(key)}"' in html
    assert 'data-retry-search' in html
    assert "const queryChanged = query !== state.query;" in html
    assert "state.selected.clear();" in html
    assert 'postJson("/api/discovery/searches"' in html
    assert 'getJson(`/api/discovery/jobs/${state.jobId}`)' in html
    assert 'getJson(`/api/discovery/searches/${state.searchId}/records?${params.toString()}`)' in html
    assert "const requestSpecies = activeSpecies;" in html
    assert "requestId !== state.searchRequest" in html
    assert "atlasContextGeneration" in html
    assert "function renderDiscoveryHeaderMeta" in html
    assert 'else if (state.query) activity = `${state.evaluatedStudies.toLocaleString()} assessed studies`' in html
    assert "Human results are not shown in the Mouse workspace" in html
    assert 'id="discoveryNotice" class="notice-box" role="status"' in html
    assert 'box.className = failed ? "error-box" : "notice-box"' in html
    assert "Sort all assessed studies by" in html
    assert "Exploratory fallback" in html
    assert "row.min_source_padj" in html
    assert "row.source_url" in html
    assert 'class="metric" tabindex="0"' in html
    assert "const fmtNullablePercent" in html
    assert 'return "N/A"' in html
    assert "Conditional 0-100 summary over available diagnostics" in html
    assert "evidence_reliability_components_used" in html
    assert "loo_rank_evaluable_folds" in html
    assert "N/A means no fold kept this gene eligible" in html
    assert 'document.addEventListener("focusin"' in html
    assert ".study-table td:nth-child(6) { grid-column: 2" in html
    assert ".study-table td:nth-child(8) { grid-column: 3" in html
    assert "topPercentTableLabel(gene)" in html
    assert 'const payload = await getJson("/api/meta");' in html
    assert "demo_search_keyword" in html
    assert "demo_search_species" in html
    # The dashboard still primes Discover defaults at boot; it now also uses the
    # scored-gene count from that same payload to choose the landing view.
    assert "loadDiscoveryDefaults()" in html
    assert 'showView(Number.isFinite(scored) && scored > 0 ? "atlas" : "discover")' in html
    assert ".mobile-study-tools button { width: auto; min-width: 118px; }" in html
    assert 'id="exportGenes"' in html
    assert 'id="exportEvidence"' in html
    assert 'id="downloadAnalysisExcel"' in html
    assert "function downloadAnalysisExcel()" in html
    assert '`/api/discovery/runs/${state.run.run_id}/export.xlsx`' in html
    assert "function exportGeneRanking()" in html
    assert "function exportCurrentEvidence()" in html
    assert "const exportQuery = currentQuery();" in html
    assert "const params = new URLSearchParams(exportQuery);" in html
    assert "Source units:</strong>" in html
    assert "downloadCsv(`degora_${context.species}_${runLabel}_ranking.csv`" in html
    assert health["status"] == "ok"
    assert health["degora_version"] == __version__
    assert health["degora_code_revision"] == "abc1234-dirty"
    assert health["degora_code_dirty"] == "true"
    assert health["database_degora_version"] == __version__
    assert health["gene_count"] == 2
    # Health must not leak the absolute on-disk db path; only the filename is exposed.
    assert "db_path" not in health
    assert health["db_name"] == "degora_scores.db"
    assert genes["genes"][0]["gene_symbol"] == "VEGFA"
    assert wildcard["count"] == 0
    assert genes["genes"][0]["rank_label"] == "#1 / 2"
    assert "top_percent_label" in genes["genes"][0]
    assert "evidence_tier" in genes["genes"][0]
    assert gene_asc["sort"] == "gene"
    assert gene_asc["order"] == "asc"
    assert [gene["gene_symbol"] for gene in gene_asc["genes"]] == ["RPL13A", "VEGFA"]
    assert score_desc["genes"][0]["gene_symbol"] == "VEGFA"
    assert detail["gene"]["gene_symbol"] == "VEGFA"
    assert detail["gene"]["support_label"] == "2 / 2 source units"
    assert detail["gene"]["direction_label"] == "100.0% up-concordant"
    assert detail["gene"]["loo_rank_stability_score"] is None
    assert detail["gene"]["loo_total_folds"] == 2
    assert detail["gene"]["loo_rank_evaluable_folds"] == 0
    assert detail["gene"]["loo_penalty_folds"] == 2
    assert detail["gene"]["loo_component_available"] is False
    assert detail["gene"]["evidence_reliability_components_used"] == 3
    assert len(detail["evidence"]) == 2
    assert {row["assay_type"] for row in detail["evidence"]} == {"RNA-seq", "microarray"}
    assert all("contributing_study_ids" in row for row in detail["evidence"])


def test_api_uses_quality_weighted_primary_rank_and_score(tmp_path) -> None:
    db = tmp_path / "degora_scores.db"
    with sqlite3.connect(db) as connection:
        pd.DataFrame(
            [
                {
                    "gene_symbol": "OLD_TOP",
                    "degora_rank": 1,
                    "degora_score": 100.0,
                    "quality_weighted_degora_rank": 2,
                    "quality_weighted_degora_score": 10.0,
                    "quality_weighted_top_percent": 100.0,
                    "quality_weighted_consensus_direction": "up",
                    "quality_weighted_sign_concordance": 0.5,
                    "n_source_units": 1,
                    "consensus_direction": "up",
                    "sign_concordance": 0.5,
                },
                {
                    "gene_symbol": "PRIMARY_TOP",
                    "degora_rank": 2,
                    "degora_score": 1.0,
                    "quality_weighted_degora_rank": 1,
                    "quality_weighted_degora_score": 99.0,
                    "quality_weighted_top_percent": 50.0,
                    "quality_weighted_consensus_direction": "down",
                    "quality_weighted_sign_concordance": 1.0,
                    "n_source_units": 1,
                    "consensus_direction": "down",
                    "sign_concordance": 1.0,
                },
            ]
        ).to_sql("genes", connection, index=False)
        pd.DataFrame({"source_unit_id": ["P1"], "study_id": ["S1"]}).to_sql("studies", connection, index=False)
        pd.DataFrame({"key": ["degora_version"], "value": [__version__]}).to_sql("meta", connection, index=False)

    server = create_server(db, port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"

    try:
        health = _get_json(f"{base_url}/api/health")
        by_rank = _get_json(f"{base_url}/api/genes?sort=rank&order=asc&limit=2")
        by_score = _get_json(f"{base_url}/api/genes?sort=score&order=desc&limit=2")
        min_score = _get_json(f"{base_url}/api/genes?min_score=90&limit=2")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert health["top_gene"] == "PRIMARY_TOP"
    assert [gene["gene_symbol"] for gene in by_rank["genes"]] == ["PRIMARY_TOP", "OLD_TOP"]
    assert [gene["gene_symbol"] for gene in by_score["genes"]] == ["PRIMARY_TOP", "OLD_TOP"]
    assert [gene["gene_symbol"] for gene in min_score["genes"]] == ["PRIMARY_TOP"]


def test_meta_redacts_local_paths_when_bound_non_loopback(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(
        harmonized_path,
        tmp_path,
        db_path=tmp_path / "degora_scores.db",
        extra_metadata={"output_dir": str(tmp_path)},
    )

    server = create_server(tmp_path / "degora_scores.db", host="0.0.0.0", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    _, port = server.server_address
    base_url = f"http://127.0.0.1:{port}"

    try:
        meta = _get_json(f"{base_url}/api/meta")["meta"]
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert meta["db_path"] == "degora_scores.db"
    assert meta["harmonized_path"] == "harmonized.csv"
    assert meta["output_dir"] == "."
    assert str(tmp_path) not in json.dumps(meta)


def test_serve_requires_explicit_network_allow_for_non_loopback(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    db = tmp_path / "degora_scores.db"
    write_score_database(harmonized_path, tmp_path, db_path=db)

    with pytest.raises(PermissionError, match="--allow-network"):
        serve(db, host="0.0.0.0", port=0, quiet=True)



def _write_minimal_score_db(db_path) -> None:
    """Write the smallest file `degora serve` accepts as a score database."""

    with sqlite3.connect(db_path) as connection:
        connection.execute("CREATE TABLE genes (gene_symbol TEXT, degora_rank INTEGER, degora_score REAL)")
        connection.execute("CREATE TABLE gene_evidence (gene_symbol TEXT, source_unit_id TEXT, study_id TEXT)")
        connection.execute("CREATE TABLE studies (study_id TEXT, source_unit_id TEXT)")
        connection.execute("CREATE TABLE meta (key TEXT, value TEXT)")


def test_serve_handles_keyboard_interrupt_cleanly(tmp_path, monkeypatch, capsys) -> None:
    db = tmp_path / "degora_scores.db"
    # serve() refuses a path that is not a DEGORA score database before it binds,
    # so this fixture needs a real one rather than a stub byte string.
    _write_minimal_score_db(db)
    state = {"closed": False}

    class InterruptingServer:
        db_path = db
        server_address = ("127.0.0.1", 8765)

        def serve_forever(self) -> None:
            raise KeyboardInterrupt

        def server_close(self) -> None:
            state["closed"] = True

    monkeypatch.setattr(api, "create_server", lambda *args, **kwargs: InterruptingServer())

    serve(db, quiet=True)

    captured = capsys.readouterr()
    assert "DEGORA browser/API: http://127.0.0.1:8765" in captured.out
    assert "Stopped DEGORA browser/API." in captured.out
    assert state["closed"]


def test_access_token_protects_api_when_configured(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    server = create_server(tmp_path / "degora_scores.db", port=0, quiet=True, access_token="secret-token")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"

    try:
        with urllib.request.urlopen(base_url, timeout=5) as response:
            html = response.read().decode("utf-8")
        with pytest.raises(HTTPError) as exc_info:
            _get_json(f"{base_url}/api/health")
        assert exc_info.value.code == 401

        request = urllib.request.Request(f"{base_url}/api/health", headers={"X-DEGORA-Token": "secret-token"})
        with urllib.request.urlopen(request, timeout=5) as response:
            health = json.loads(response.read().decode("utf-8"))
        query_health = _get_json(f"{base_url}/api/health?token=secret-token")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert "readApiToken()" in html
    assert health["status"] == "ok"
    assert query_health["status"] == "ok"


def test_loopback_server_rejects_nonlocal_host_header(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    db = tmp_path / "degora_scores.db"
    write_score_database(harmonized_path, tmp_path, db_path=db)
    server = create_server(db, port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address

    try:
        connection = http.client.HTTPConnection(host, port, timeout=5)
        connection.putrequest("GET", "/api/health", skip_host=True)
        connection.putheader("Host", "attacker.example")
        connection.endheaders()
        response = connection.getresponse()
        payload = json.loads(response.read().decode("utf-8"))
        connection.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert response.status == 421
    assert "Host header" in payload["error"]


def test_query_token_is_redacted_from_access_log(monkeypatch) -> None:
    handler = object.__new__(api.DegoraRequestHandler)
    handler.server = type("Server", (), {"quiet": False})()
    captured: dict[str, object] = {}

    def fake_log_message(self, fmt: str, *args: object) -> None:
        captured["format"] = fmt
        captured["args"] = args

    monkeypatch.setattr(api.BaseHTTPRequestHandler, "log_message", fake_log_message)

    handler.log_message('"%s" %s', "GET /api/health?token=secret-token HTTP/1.1", "200")

    assert captured["args"] == ("GET /api/health?token=[redacted] HTTP/1.1", "200")


def test_percent_encoded_token_query_is_redacted() -> None:
    text = "GET /api/health?to%6Ben=secret%2Dtoken&limit=1 HTTP/1.1"

    redacted = redact_token_from_url_text(text)

    assert "secret%2Dtoken" not in redacted
    assert "to%6Ben=[redacted]" in redacted
    assert strip_token_query_param("http://127.0.0.1/api/health?to%6Ben=secret&limit=1").endswith(
        "/api/health?limit=1"
    )


def test_sqlite_readonly_uri_uses_path_as_uri(tmp_path, monkeypatch) -> None:
    db = tmp_path / "space dir" / "scores.db"
    db.parent.mkdir()
    db.write_bytes(b"")
    captured: dict[str, object] = {}

    class FakeConnection:
        row_factory = None

    def fake_connect(database: str, *, uri: bool = False):
        captured["uri"] = database
        captured["uri_flag"] = uri
        return FakeConnection()

    monkeypatch.setattr(api.sqlite3, "connect", fake_connect)

    connection = api._connect(db)

    assert isinstance(connection, FakeConnection)
    assert captured["uri"] == f"{db.resolve().as_uri()}?mode=ro"
    assert captured["uri_flag"] is True


def test_serve_prints_token_as_fragment_not_query(tmp_path, monkeypatch, capsys) -> None:
    db = tmp_path / "degora_scores.db"
    # serve() refuses a path that is not a DEGORA score database before it binds,
    # so this fixture needs a real one rather than a stub byte string.
    _write_minimal_score_db(db)
    state = {"closed": False}

    class InterruptingServer:
        db_path = db
        server_address = ("127.0.0.1", 8765)

        def serve_forever(self) -> None:
            raise KeyboardInterrupt

        def server_close(self) -> None:
            state["closed"] = True

    monkeypatch.setattr(api, "create_server", lambda *args, **kwargs: InterruptingServer())

    serve(db, quiet=True, access_token="secret token")

    captured = capsys.readouterr()
    assert "DEGORA browser/API: http://127.0.0.1:8765#token=secret%20token" in captured.out
    assert "?token=" not in captured.out
    assert state["closed"]


def test_network_api_redacts_source_paths_in_studies_and_gene_evidence(tmp_path) -> None:
    secret_source = tmp_path / "Users" / "researcher" / "project" / "source.csv"
    harmonized = _harmonized()
    harmonized["source_path"] = str(secret_source)
    harmonized_path = tmp_path / "harmonized.csv"
    harmonized.to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    server = create_server(tmp_path / "degora_scores.db", host="0.0.0.0", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    _, port = server.server_address
    base_url = f"http://127.0.0.1:{port}"

    try:
        studies = _get_json(f"{base_url}/api/studies")["studies"]
        detail = _get_json(f"{base_url}/api/genes/VEGFA")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    payload = json.dumps({"studies": studies, "detail": detail})
    assert str(tmp_path) not in payload
    stored_paths = [
        *[study["source_path"] for study in studies],
        *[row["source_path"] for row in detail["evidence"]],
        *[row["contributing_source_paths"] for row in detail["evidence"]],
    ]
    assert all(value.startswith(EXTERNAL_PATH_PREFIX) for value in stored_paths)
    assert all(value.endswith("/source.csv") for value in stored_paths)


def test_meta_redaction_catches_posix_paths_on_windows() -> None:
    from degora.api import _formula_neutral, _redact_meta_for_network, _redact_record_paths_for_network

    meta = _redact_meta_for_network(
        {
            "db_path": "/mnt/c/Projects/DEGORA/outputs/degora_scores.db",
            "output_dir": "C:\\Projects\\DEGORA\\outputs",
            "source_url": "https://example.test/data",
            "notes": "created from /Users/researcher/private/source.tsv before upload",
            "file_url": "see file:///Users/researcher/private/source.tsv",
            "windows_note": "source C:\\Users\\researcher\\private\\source.tsv",
            "unc_note": "source \\\\server\\share\\source.tsv",
        }
    )

    assert meta["db_path"] == "[redacted: local path]"
    assert meta["output_dir"] == "[redacted: local path]"
    assert meta["source_url"] == "https://example.test/data"
    assert meta["notes"] == "[redacted: local path]"
    assert meta["file_url"] == "[redacted: local path]"
    assert meta["windows_note"] == "[redacted: local path]"
    assert meta["unc_note"] == "[redacted: local path]"
    record = _redact_record_paths_for_network(
        {
            "source_path": "/Users/researcher/private/cohort_2019_patients.csv",
            "output_dir": "C:\\Projects\\DEGORA\\patient-results",
            "source_url": "https://example.test/data",
            "notes": "created from /Users/researcher/private/source.tsv before upload",
            "windows_note": "source C:\\Users\\researcher\\private\\source.tsv",
            "unc_note": "source \\\\server\\share\\source.tsv",
        }
    )
    assert record["source_url"] == "https://example.test/data"
    assert record["source_path"] == "[redacted: local path]"
    assert record["output_dir"] == "[redacted: local path]"
    assert record["notes"] == "[redacted: local path]"
    assert record["windows_note"] == "[redacted: local path]"
    assert record["unc_note"] == "[redacted: local path]"
    assert _formula_neutral("  =SUM(A1:A2)") == "'  =SUM(A1:A2)"
    assert _formula_neutral("\t+cmd") == "'\t+cmd"
    assert _formula_neutral("#VALUE!") == "'#VALUE!"
    assert _formula_neutral(" \t#REF!") == "' \t#REF!"


def test_genes_pagination_offset_count_and_no_overlap(tmp_path) -> None:
    # 150 scored genes (each supported by 3 source units) so the list spans >1 page of 100.
    rows = []
    for i in range(150):
        for unit in ("P1", "P2", "P3"):
            rows.append(
                {
                    "study_id": f"{unit}_S",
                    "paper_id": unit,
                    "gene_symbol": f"GENE{i:03d}",
                    "lfc": 2.0,
                    "signed_z": 6.0 - i * 0.01,
                    "pvalue": 1e-6,
                    "padj": 1e-5,
                    "normalized_rank": 0.001 + i * 0.0001,
                    "n_ctrl": 3,
                    "n_treat": 3,
                    "n_genes_in_study": 20000,
                    "pipeline": "DESeq2",
                    "assay_type": "RNA-seq",
                    "source_input_type": "author_deg_table",
                    "platform": "",
                    "normalization": "DESeq2",
                    "probe_collapse": "",
                    "species": "Homo sapiens",
                    "cell_system": "A",
                    "hypoxia_modality": "x",
                    "duration_h": "24",
                    "source_path": "s.csv",
                    "source_url": "u",
                }
            )
    harmonized_path = tmp_path / "harmonized.csv"
    pd.DataFrame(rows).to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    server = create_server(tmp_path / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"
    try:
        page1 = _get_json(f"{base_url}/api/genes?limit=100&offset=0")
        page2 = _get_json(f"{base_url}/api/genes?limit=100&offset=100")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    # count is the full match total; pages are capped at limit and addressable by offset.
    assert page1["count"] == 150
    assert len(page1["genes"]) == 100
    assert len(page2["genes"]) == 50
    first = [gene["gene_symbol"] for gene in page1["genes"]]
    second = [gene["gene_symbol"] for gene in page2["genes"]]
    # 'Load more' appends a fresh, non-overlapping page...
    assert set(first).isdisjoint(second)
    # ...and ranks stay globally ordered and distinct across both pages.
    ranks = [gene["degora_rank"] for gene in page1["genes"]] + [gene["degora_rank"] for gene in page2["genes"]]
    assert ranks == sorted(ranks)
    assert len(set(ranks)) == 150


def test_genes_api_rejects_invalid_numeric_query_params(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    server = create_server(tmp_path / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"
    try:
        for query, expected in [
            ("limit=abc", "limit must be an integer"),
            ("limit=0", "limit must be between 1 and 500"),
            ("offset=-1", "offset must be between 0 and 1000000"),
            ("min_units=1.5", "min_units must be an integer"),
            ("min_score=nan", "min_score must be finite"),
            ("min_score=-1", "min_score must be >= 0"),
        ]:
            with pytest.raises(HTTPError) as exc_info:
                _get_json(f"{base_url}/api/genes?{query}")
            assert exc_info.value.code == 400
            body = exc_info.value.read().decode("utf-8")
            assert expected in body
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_health_does_not_leak_absolute_db_path_when_db_missing(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    db = tmp_path / "secret_dir" / "degora_scores.db"
    db.parent.mkdir()
    write_score_database(harmonized_path, db.parent, db_path=db)

    server = create_server(db, port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"
    try:
        db.unlink()  # database removed while the server is running
        try:
            _get_json(f"{base_url}/api/health")
            raise AssertionError("expected an error after the database was removed")
        except HTTPError as exc:
            body = exc.read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    # The HTTP error must not leak the absolute path / directory names, only the filename.
    assert "secret_dir" not in body
    assert str(tmp_path) not in body
    assert "not available" in body


def test_gene_detail_rejects_overlong_symbol_with_400(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    server = create_server(tmp_path / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"
    try:
        try:
            _get_json(f"{base_url}/api/genes/{'A' * 5000}")
            raise AssertionError("expected 400 for an overlong gene symbol")
        except HTTPError as exc:
            code = exc.code
            body = exc.read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert code == 400
    assert "too long" in body
    assert len(body) < 500  # the 5000-char symbol is not echoed back


def test_serve_rejects_missing_database(tmp_path) -> None:
    with pytest.raises((FileNotFoundError, ScoreDatabaseError)):
        serve(tmp_path / "does_not_exist.db")


def test_loopback_hosts_membership() -> None:
    assert "127.0.0.1" in LOOPBACK_HOSTS
    assert "localhost" in LOOPBACK_HOSTS
    assert "0.0.0.0" not in LOOPBACK_HOSTS


def test_discovery_readiness_sort_prefers_numeric_score() -> None:
    record = {"readiness": "Likely author DEG", "readiness_score": 0.93}
    assert api._record_sort_value(record, "readiness") == 0.93


def test_create_server_falls_back_when_port_in_use(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")
    db = tmp_path / "degora_scores.db"

    first = create_server(db, port=0, quiet=True, discovery_root=tmp_path / "discovery-first")
    try:
        busy_port = first.server_address[1]
        # Requesting the already-bound port must not raise; it auto-falls-back.
        second = create_server(db, port=busy_port, quiet=True, discovery_root=tmp_path / "discovery-second")
        try:
            assert second.server_address[1] != busy_port
            assert second.server_address[1] != 0
        finally:
            second.server_close()
    finally:
        first.server_close()


def test_create_server_rejects_duplicate_discovery_workspace(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")
    db = tmp_path / "degora_scores.db"
    discovery_root = tmp_path / "shared-discovery"

    first = create_server(db, port=0, quiet=True, discovery_root=discovery_root)
    try:
        with pytest.raises(DiscoveryWorkspaceInUseError, match="already using this discovery workspace"):
            create_server(db, port=0, quiet=True, discovery_root=discovery_root)
    finally:
        first.server_close()


def test_dashboard_shows_real_search_progress_instead_of_a_static_panel() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # Determinate bar driven by the job's progress fraction.
    assert ".loading-bar.is-determinate" in html
    assert 'class="loading-bar is-determinate"' in html
    assert "state.jobProgress" in html
    assert "state.jobMessage" in html
    # Elapsed time gives movement even while a stage is long.
    assert "formatElapsed" in html
    assert "state.jobStartedAt" in html
    # The poller must read the new fields off the job payload.
    assert "typeof job.progress" in html
    assert "typeof job.message" in html
    # Animation is suppressed for users who ask for reduced motion.
    assert "prefers-reduced-motion" in html


def test_dashboard_reports_where_the_selection_actually_is() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # "20 / 20" used to read as "all 20 rows on this page" while the selection
    # lived on another page and the cap was global.
    assert "selected of max" in html
    assert "on this page" in html
    assert "20 / ${MAX_SELECTED_STUDIES} selected" not in html


def test_dashboard_disables_unselectable_publications_instead_of_ignoring_clicks() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "pageSelectability" in html
    assert "selectableKeys" in html
    assert "ambiguous-id" in html and "no-id" in html
    assert "Several results share this identifier" in html
    assert "has no usable identifier" in html
    assert "Selection limit of ${MAX_SELECTED_STUDIES} reached" in html


def test_dashboard_guards_against_double_submitting_a_search() -> None:
    from degora.api import INDEX_HTML

    assert "if (activeDiscoveryState().loading) return;" in INDEX_HTML


def test_dashboard_ranks_notices_by_severity_and_recency() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # A stale success notice must not occupy the alert region while a real
    # failure is hidden, and a newer prepare failure must not be masked by an
    # older search error.
    assert 'state.noticeLevel === "error" && state.notice' in html
    assert 'noticeLevel: "info"' in html
    assert 'state.noticeLevel = "error";' in html
    # Clearing a notice must not silently discard an error-level one.
    assert 'state.noticeLevel !== "error"' in html


def test_dashboard_clamps_atlas_filters_before_calling_the_api() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert 'maxlength="128"' in html
    assert 'max="10000"' in html
    assert "Math.min(10000, Math.max(1, parsed))" in html


def test_stacked_evidence_panel_is_not_clipped_on_narrow_viewports() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert ".genes-panel, .evidence-panel { height: auto; overflow: visible; }" in html
    assert ".evidence-panel .detail-body { overflow: visible; flex: none; height: auto; }" in html


def test_keyboard_focus_is_visible_on_sort_headers_and_the_search_input() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert ".sort-head:focus-visible {" in html
    assert ".discovery-search input:focus-visible {" in html
    # The hover tint alone must no longer stand in for a focus ring.
    assert ".sort-head:hover, .sort-head:focus-visible { background: #eef7f5; outline: none; }" not in html


def test_dashboard_guards_enter_as_well_as_the_search_button() -> None:
    from degora.api import INDEX_HTML

    # Enter used to bypass the double-submit guard and start extra server jobs.
    assert 'if (event.key !== "Enter" || activeDiscoveryState().loading) return;' in INDEX_HTML


def test_dashboard_keeps_keyboard_focus_when_toggling_a_selection() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "renderDiscoveryResultsKeepingFocus" in html
    assert "CSS.escape(accession)" in html


def test_inspect_respects_the_unselectable_rules() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # Inspect used to clear the selection and then fail server-side on a
    # publication whose identifier is shared by several rows.
    assert "const blocked = pageSelectability(state)(key);" in html


def test_repaging_resets_the_previous_jobs_progress() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert html.count('state.jobStartedAt = Date.now();') >= 2
    assert 'state.jobMessage = "";' in html


def test_select_all_is_disabled_once_the_cap_blocks_every_row() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "capBlocksEveryRow" in html
    assert "selectAll.disabled = pageAccessions.length === 0 || capBlocksEveryRow;" in html


def test_a_new_snapshot_drops_identifiers_from_the_previous_one() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "// A new snapshot re-mints record ids" in html


def test_excluded_study_cards_are_labelled() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "item.paper_title || item.canonical_id || item.source_unit_id" in html
    assert "<h4>${esc(item.accession)}</h4>" not in html


def test_progress_uses_one_live_region_and_does_not_announce_each_second() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # The card is a progressbar, not a second polite live region duplicating
    # the subtitle once per poll.
    assert 'class="loading-card search-progress" aria-busy="true"' in html
    assert 'role="progressbar"' in html
    assert '$("resultsSubtitle").textContent = `${stage}${percentText}`;' in html


def test_atlas_filter_clamp_is_reflected_in_the_field() -> None:
    from degora.api import INDEX_HTML

    assert '$("minUnits").value = clamped;' in INDEX_HTML


def test_enter_applies_atlas_filters_from_either_input() -> None:
    from degora.api import INDEX_HTML

    assert '["query", "minUnits"].forEach((id) => {' in INDEX_HTML


def test_dashboard_locks_and_fades_controls_while_work_is_in_flight() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "function applyDiscoveryBusyState(state)" in html
    assert 'view.classList.toggle("is-busy", busy);' in html
    assert 'view.setAttribute("aria-busy", busy ? "true" : "false");' in html
    # Dimming alone would still leave the controls reachable by keyboard.
    assert "BUSY_OWNED_CONTROLS" in html and "BUSY_FORCED_CONTROLS" in html
    assert "if (element) element.disabled = busy;" in html
    assert ".discovery-view.is-busy .discovery-search" in html
    assert "pointer-events: none;" in html
    # The progress panel must not be dimmed with everything else.
    assert ".discovery-view.is-busy .search-progress" not in html


def test_a_new_search_retires_in_flight_prepare_and_analysis_work() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # A late prepare/analysis response used to reinstate a bundle built from the
    # snapshot the new search had just replaced.
    assert "state.prepareRequest += 1;" in html
    assert "state.analysisRequest += 1;" in html


def test_finishing_an_analysis_reloads_a_visible_atlas() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "function isDiscoveryView()" in html
    assert "if (!isDiscoveryView()) void ensureAtlasContext();" in html


def test_inspect_does_not_overwrite_the_selection() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "async function prepareSelectedStudies({ recordIds: explicitIds = null } = {})" in html
    assert "void prepareSelectedStudies({ recordIds: [key] });" in html
    # The old body cleared the whole selection to prepare one publication.
    assert "state.selected.clear();\n        state.selected.add(key);" not in html


def test_atlas_view_refreshes_the_header_metadata() -> None:
    from degora.api import INDEX_HTML

    assert "void loadMeta(currentAtlasContext(), atlasContextGeneration);" in INDEX_HTML


def test_header_activity_follows_the_prepare_and_analysis_pipeline() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert html.count("renderDiscoveryHeaderMeta();") >= 3


def test_prepared_study_heading_has_no_leading_separator() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert '${esc(study.accession)} · ${esc(study.paper_title' not in html
    assert 'filter(Boolean).join(" · ")' in html


def test_hidden_action_bar_never_keeps_the_other_workspaces_numbers() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # Both early-return branches must refresh the selection status so the bar
    # cannot be revealed later carrying the previous species' counts.
    assert html.count('$("discoveryFooter").hidden = true;\n        updateSelectedStatus();') >= 2


def test_prepare_shows_the_same_determinate_progress_as_search() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "async function pollPrepareJob(species, requestId, jobId)" in html
    assert '"/api/discovery/prepare-jobs"' in html
    assert "state.prepareProgress" in html and "state.prepareMessage" in html
    assert "Downloading and inspecting the selected publications" in html
    # The bundle now arrives through the job record rather than the POST response.
    assert "job.result" in html


def test_a_degraded_provider_set_is_reported_rather_than_shown_as_no_results() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "state.providerStatus" in html and "state.providerErrors" in html
    assert "partial snapshot · some sources unavailable" in html
    assert "Some data sources did not answer" in html
    assert "before concluding that the query has no matches" in html


def test_splitter_does_not_snap_on_press_and_uses_the_content_box() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    # A press used to jump the handle to the cursor, so the first click of a
    # double-click moved it out from under the second.
    assert "const DRAG_THRESHOLD_PX = 3;" in html
    assert "const beginPress = (clientX) => {" in html
    assert "const trackClientX = (clientX) => {" in html
    # The grid track percentage is of the content box, not the border box.
    assert "const contentGeometry = () => {" in html
    assert "parseFloat(styles.paddingLeft)" in html
    assert "setPanelSplit(((clientX - rect.left) / rect.width) * 100);" not in html


def test_dashboard_explains_the_selection_cap_before_a_click_that_cannot_work() -> None:
    """A disabled checkbox fires no event, so the click handler's notice is unreachable.

    Filling page 1 uses the whole cross-page budget, and the reader's next move
    is to press "select all" on page 2 and get silence. The page has to say why
    before anything is pressed.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "selection-limit" in html
    assert "Selection limit reached: ${MAX_SELECTED_STUDIES} publications across all pages" in html
    assert "Untick a row or press Clear to choose different ones here." in html
    # And a locked row has to look locked; browsers barely grey a 17px checkbox.
    assert '.study-table input[type="checkbox"]:disabled' in html
    assert "cursor: not-allowed" in html


def test_dashboard_orders_the_first_page_by_deg_readiness() -> None:
    """Readiness first, provider relevance only as a tie-break."""

    from degora.api import INDEX_HTML

    assert 'sort: { key: "readiness", order: "desc" }' in INDEX_HTML
    assert "Sort: DEG readiness" in INDEX_HTML


def test_dashboard_names_the_files_a_study_could_not_use() -> None:
    """"No usable table within the safety limits" names no file and no limit.

    A study that published only browser tracks and one whose table was a
    megabyte over the cap produced the same sentence, so a reader could not
    tell which of the two they were looking at.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "unusableStudyHtml" in html
    assert "candidateRejection" in html
    assert "were found and none could be used" in html
    assert "The repository listed no supplementary file for this study." in html
    # A rejected file was never inspected, so its inspection note is useless here.
    assert 'candidate.tier === "reject" || candidate.role === "unsupported"' in html


def test_dashboard_labels_matrix_columns_from_the_candidate_first() -> None:
    """Author matrices are headed by submitter names, resolved during preparation."""

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "inspection.sample_labels" in html
    assert "GEO returned no matching sample labels for these columns" in html
    # Printing the submitter title when it merely repeats the column wastes the
    # line that should carry the characteristics.
    assert "echoesColumn" in html


def test_dashboard_replaces_the_search_estimate_with_the_prepared_outcome() -> None:
    """A row that promised likely_ready and delivered nothing must stop promising it.

    The search badge is scored from metadata. Once preparation has opened the
    files, the row should report what was found, or the same study gets picked
    again next week.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "preparedOutcomes" in html
    assert "prepared · ${outcome.label}" in html
    for label in ("author DEG ready", "needs group assignment", "no usable table", "excluded"):
        assert label in html


def test_dashboard_says_how_much_usable_evidence_a_preparation_produced() -> None:
    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "prepared stud${total === 1 ? \"y\" : \"ies\"} produced a usable candidate" in html
    assert "One independent ${speciesLabel(activeSpecies)} study is selected; DEGORA needs two." in html


def test_dashboard_offers_filtered_bulk_sample_assignment_that_clears_attestations() -> None:
    """Twenty dropdowns one at a time is where attention runs out.

    Bulk assignment is only safe if it stays explicit: the reader sees the rows
    it will touch and the count on the button, and the direction attestations -
    made about the previous assignment - are cleared by the move.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "sampleBulkHtml" in html and "applySampleBulk" in html
    assert "matchingSampleItems" in html and "refreshSampleFilter" in html
    assert "Filter by label" in html
    assert '["direction-confirmed", "biological-replicates-confirmed"].forEach' in html
    # The button states its target count, so pressing is a decision about a number.
    assert "`Set ${names[button.dataset.group] || \"Ignore\"} (${matched})`" in html
    # Only offered where it earns its space.
    assert "columns.length > 4 ? sampleBulkHtml() : \"\"" in html


def test_dashboard_shows_the_prepare_progress_where_the_reader_is_looking() -> None:
    """Prepare runs for tens of seconds behind a full page of results.

    Pressing the button and seeing nothing move is indistinguishable from a
    button that did not work.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "revealPreparedCard" in html
    # Brought into view when the work starts, not only when it finishes.
    assert "updateSelectedStatus();\n        revealPreparedCard();" in html
    # And left alone when it is already on screen, or when motion is unwanted.
    assert "viewport * 0.75" in html
    assert "prefers-reduced-motion: reduce" in html


def test_dashboard_stops_a_review_that_can_never_be_analysed() -> None:
    """One usable study cannot reach two independent source units.

    Every field would be filled for nothing, and worse, state.draft is cleared
    on the next preparation, so the work is provably discarded.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "usableSourceUnits" in html
    assert "This preparation cannot be analysed." in html
    assert "preparing a new selection clears them" in html
    assert "is-unanalyzable" in html
    assert "data-back-to-results" in html
    # The controls go inert; the file names and reasons stay readable.
    assert '.querySelectorAll(".candidate-row input, .candidate-row select, .candidate-row button")' in html
    assert "Not analysable: " in html


def test_a_browser_page_holds_ten_publications() -> None:
    """The page size is stated in four places and they have to agree.

    The number appears in the dashboard constant, the caption a reader sees,
    the CLI's own paging, and the ranking contract written into every audit
    bundle - so a change in one place that misses the others publishes a claim
    the tool does not honour.
    """

    import inspect

    from degora import discovery_federated
    from degora.api import INDEX_HTML
    from degora.cli import DISCOVERY_PAGE_SIZE
    from degora.discovery import DEFAULT_PAGE_SIZE

    assert DEFAULT_PAGE_SIZE == 10
    assert DISCOVERY_PAGE_SIZE == DEFAULT_PAGE_SIZE
    assert f"const DISCOVERY_PAGE_SIZE = {DEFAULT_PAGE_SIZE};" in INDEX_HTML
    assert f"{DEFAULT_PAGE_SIZE} per page · up to 1,000 ranked" in INDEX_HTML
    assert f"{DEFAULT_PAGE_SIZE} rows per browser page" in inspect.getsource(discovery_federated)


def test_no_prose_states_a_page_size_the_tool_does_not_honour() -> None:
    """Catch the stale page-size claim the four-constant check above cannot see.

    That check compares constants and two exact dashboard strings. It passed
    while the README still promised a "20-row page" and the CLI docstring still
    described a "20-row publication page", because both used a phrasing no
    assertion named. Scan the shipped prose for any stated page size instead, so
    the wording can change without the number drifting loose again.
    """

    import re

    from degora.discovery import DEFAULT_PAGE_SIZE

    root = Path(__file__).resolve().parents[1]
    sources = [root / "README.md", root / "CHANGELOG.md", *sorted((root / "degora").glob("*.py"))]
    # "N-row page", "N rows per page", "N per page" -- every way the docs have
    # spelled a page size so far.
    claim = re.compile(r"(\d+)(?:-row page|\s+rows?\s+per\s+(?:page|browser page)|\s+per page)")

    offenders = []
    for path in sources:
        for number, line_no, line in (
            (match.group(1), index, text)
            for index, text in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
            for match in claim.finditer(text)
        ):
            if int(number) != DEFAULT_PAGE_SIZE:
                offenders.append(f"{path.name}:{line_no}: {line.strip()[:120]}")

    assert not offenders, "prose states a page size other than {}:\n{}".format(
        DEFAULT_PAGE_SIZE, "\n".join(offenders)
    )


def test_a_page_no_longer_exhausts_the_selection_budget_on_its_own() -> None:
    """Ten per page against a cap of twenty is why the cap now spans pages.

    With both at twenty, ticking select-all on page one spent the whole budget
    and page two was locked before the reader ever reached it.
    """

    from degora.api import INDEX_HTML
    from degora.discovery import DEFAULT_PAGE_SIZE, MAX_SELECTED_STUDIES

    assert DEFAULT_PAGE_SIZE < MAX_SELECTED_STUDIES
    assert f"const MAX_SELECTED_STUDIES = {MAX_SELECTED_STUDIES};" in INDEX_HTML


def test_dashboard_warns_before_counting_one_submission_as_replication() -> None:
    """The gate counts source units, and an unpublished submission split across
    several series counts as several.

    Nothing in the interface said so, and the whole replication claim is that
    count.
    """

    from degora.api import INDEX_HTML

    html = INDEX_HTML
    assert "sharedSubmissionNote" in html
    assert "sharedSubmissionConflicts" in html
    assert "renderIndependenceWarning" in html
    assert "These may not be independent studies." in html
    assert "would treat correlated data as replication" in html
    assert "May be one submission with" in html


def test_changing_the_analysis_context_clears_the_gene_filters() -> None:
    """A filter left over from one run must not silently apply to the next.

    Switching context reset the rows, the page and the detail pane but not the
    gene search box, so an analysis of 11,886 genes opened showing the nine that
    matched a `TP53` filter typed against the previous run - which reads as an
    analysis that found nine genes.
    """

    from degora.api import INDEX_HTML

    start = INDEX_HTML.index("async function ensureAtlasContext()")
    end = INDEX_HTML.index("$(\"genes\").innerHTML = \"\";", start)
    reset = INDEX_HTML[start:end]

    assert '$("query").value = "";' in reset
    assert '$("minUnits").value = "1";' in reset
    assert '$("direction").value = "";' in reset
    # The reset values have to be the control's own defaults.
    assert '<option value="">All directions</option>' in INDEX_HTML
    assert 'id="minUnits" type="number" min="1" max="10000" step="1" value="1"' in INDEX_HTML


def test_a_confirmation_is_only_shown_where_it_gates() -> None:
    """Six assertions on every row made the one that mattered read like the rest.

    The activation gate already asks for each confirmation only where it applies -
    a mapping confirmation when the mapping was edited, a log2 confirmation only
    for a table whose effect column does not say it is log2, and so on. The review
    panel rendered all six on every row regardless, so a table needing one
    assertion presented six bioinformatics judgements.

    Direction is deliberately not in this list: it gates every row, always.
    """

    from degora.api import INDEX_HTML

    assert ".confirm-line.not-required { display: none !important; }" in INDEX_HTML
    assert "const showWhenRequired" in INDEX_HTML

    # Each conditional confirmation is toggled by the same condition the gate uses.
    for control, condition in (
        ("mappingConfirmed", "needsMapping"),
        ("adjustedPConfirmed", 'Boolean(pValue && padjValue && pValue === padjValue)'),
        ("lfcScaleConfirmed", 'row.dataset.authorStatus === "requires_lfc_confirmation"'),
        ("rowFilterConfirmed", 'Boolean(textValue(row, ".row-filter-column"))'),
        ("duplicateGenePolicyConfirmed", 'duplicateGenePolicy?.value === "keep_first"'),
    ):
        assert f"showWhenRequired({control}, {condition})" in INDEX_HTML, control

    # A hidden confirmation must not carry a stale tick into the gate.
    assert "if (!required && control.checked) control.checked = false;" in INDEX_HTML

    # Direction is never hidden, and says what it means without naming log2FC.
    assert "showWhenRequired(direction" not in INDEX_HTML
    assert "A positive value here means the gene went UP in the treated group" in INDEX_HTML


def test_the_review_confirmations_are_written_for_a_non_specialist() -> None:
    """The wording is the interface here: an unread assertion is an unmade one."""

    from degora.api import INDEX_HTML

    for phrase in (
        "The gene, effect and p-value columns chosen above are the right ones",
        "This table has no separate raw p-value, so its adjusted p-value is being used as one",
        "The effect column is already a log2 fold change, not a plain ratio",
        "The filter above picks one comparison and does not mix several together",
        "Keeping the first row for each repeated gene reproduces how this table was originally read",
    ):
        assert phrase in INDEX_HTML, phrase

    # The jargon these replaced must be gone, not merely joined.
    for jargon in (
        "Sheet and gene/log2FC/p mappings are intentionally confirmed",
        "Adjusted-p/FDR column may be used as p-value for this activation",
        "documented source-order workflow",
    ):
        assert jargon not in INDEX_HTML, jargon


def test_every_hidden_confirmation_sits_in_something_that_can_be_hidden() -> None:
    """The toggle finds its row with closest(".confirm-line"), so each must have one.

    A control rendered outside that wrapper would be toggled by a call that
    silently does nothing, and the row would keep asking a question that does not
    gate it - the exact failure this change exists to remove, reintroduced quietly.
    """

    import re

    from degora.api import INDEX_HTML

    wrapped = set(
        re.findall(r'<label class="confirm-line[^"]*"[^>]*><input class="([a-z0-9-]+)"', INDEX_HTML)
    )
    toggled = re.findall(r"showWhenRequired\((\w+),", INDEX_HTML)
    control_classes = {
        "mappingConfirmed": "column-mapping-confirmed",
        "adjustedPConfirmed": "adjusted-p-as-pvalue-confirmed",
        "lfcScaleConfirmed": "lfc-scale-confirmed-log2",
        "rowFilterConfirmed": "row-filter-confirmed",
        "duplicateGenePolicyConfirmed": "duplicate-gene-policy-confirmed",
    }

    assert set(toggled) == set(control_classes), toggled
    for name in toggled:
        assert control_classes[name] in wrapped, name

    # The helper still looks for that wrapper, so the two halves cannot drift apart.
    helper = re.search(r"const showWhenRequired = \(control, required\) => \{(.*?)\};", INDEX_HTML, re.S)
    assert helper is not None
    assert 'closest(".confirm-line")' in helper.group(1)


def test_the_readiness_badge_says_what_the_estimate_rests_on() -> None:
    """"likely_ready" reads as a promise; it is an estimate from metadata.

    An audit inspected the top-ranked likely_ready candidate and found no usable
    file at all, and five publications across eight repository records yielded two
    analysable datasets. The label already said "search estimate", but what the
    estimate was built from - how many candidate files had actually been seen -
    appeared nowhere on the row.
    """

    from degora.api import INDEX_HTML

    assert "may have usable data" in INDEX_HTML
    assert "not inspected yet" in INDEX_HTML
    assert "data confirmed" in INDEX_HTML
    # The count is in the visible label, not only a tooltip nobody hovers.
    assert 'candidate file${files.length === 1 ? "" : "s"}' in INDEX_HTML
    assert "nothing inspected yet" in INDEX_HTML
    # The machine state stays reachable for an audit, in the tooltip.
    assert "`state: ${readiness}`" in INDEX_HTML
    # The raw tier no longer stands alone as the reader-facing label.
    assert "`search estimate · ${readiness}`" not in INDEX_HTML


def test_every_interpolated_name_in_a_candidate_row_exists() -> None:
    """A template can reference a variable nobody declared and still match a string test.

    While collapsing the advanced settings, the markup referencing `columnsOpen`
    landed without the declaration: every string assertion still passed, and the
    row would have thrown a ReferenceError in the browser. This checks the two
    functions that build a prepared-candidate row for bare interpolations whose
    name is not declared in the function, its parameters, or module scope.
    """

    import re

    from degora.api import INDEX_HTML

    builtins = {"esc", "fmt", "badge", "tier", "Boolean", "String", "Number", "Math", "JSON", "Object", "Array"}
    module_level = set(re.findall(r"\n    (?:function|const|let)\s+([A-Za-z_$][\w$]*)", INDEX_HTML))

    for name in ("authorCandidateHtml", "fallbackCandidateHtml"):
        start = INDEX_HTML.index(f"function {name}(")
        end = INDEX_HTML.index("\n    function ", start + 10)
        body = INDEX_HTML[start:end]

        declared = set(re.findall(r"\b(?:const|let|var)\s+([A-Za-z_$][\w$]*)", body))
        parameters = re.search(rf"function {name}\(([^)]*)\)", body).group(1)
        declared |= {part.strip() for part in parameters.split(",") if part.strip()}

        # ${identifier} and ${identifier ? ... } - the shapes that name a variable
        # outright rather than reading a property off one.
        used = set(re.findall(r"\$\{([A-Za-z_$][\w$]*)\s*(?:\?|\})", body))
        assert used, name
        undefined = sorted(used - declared - module_level - builtins)
        assert not undefined, f"{name} interpolates undeclared: {undefined}"


def test_advanced_settings_start_collapsed_but_never_hide_a_set_value() -> None:
    """Eight optional controls competed with the two or three that need deciding.

    Collapsing them is only safe if a setting already carrying a value opens the
    panel: otherwise the reader cannot see what their run is going to do.
    """

    from degora.api import INDEX_HTML

    assert "<details class=\"candidate-advanced\"" in INDEX_HTML
    assert "function anyValueSet(" in INDEX_HTML
    # Advanced opens when anything in it is set, including the non-default policy.
    assert "const advancedOpen = anyValueSet(" in INDEX_HTML
    assert 'duplicateGenePolicy === "keep_first"' in INDEX_HTML
    # Columns open when the file was not read cleanly or the reader edited them -
    # edited, not merely captured: the draft holds the prefilled detected mapping
    # after any re-render, and that must not count as "set" (0.4.29).
    assert 'const columnsOpen = status !== "ready_for_review"' in INDEX_HTML
    assert "const mappingEdited = " in INDEX_HTML
    assert "|| mappingEdited" in INDEX_HTML
    assert "anyValueSet(draft.sheetName, draft.geneColumn" not in INDEX_HTML
    # Collapsed, not removed: the activation gate still queries these controls.
    assert 'class="row-filter-column"' in INDEX_HTML
    assert 'class="duplicate-gene-policy"' in INDEX_HTML


def test_group_sizes_are_bounded_by_the_series_rather_than_guessed() -> None:
    """The split between groups is not derivable, and it feeds the source weight.

    A results table has one row per gene, so nothing in it says how many samples
    were in each arm; inferring it would put an unverifiable number straight into
    min(sqrt(n_ctrl + n_treat), 4). What the linked series holds in total is
    knowable, and is what a reader is otherwise squinting at the paper to find.
    """

    from degora.api import INDEX_HTML

    assert "The linked series lists" in INDEX_HTML
    assert "together they cannot exceed" in INDEX_HTML
    assert 'data-series-samples="' in INDEX_HTML

    # Entering the series total in both boxes doubles the source's weight against
    # every other study, so the pair is checked against the series.
    assert "const fitsSeries = !seriesTotal || enteredTotal <= seriesTotal;" in INDEX_HTML
    assert "isPositiveWholeNumber(nCtrl?.value) && fitsSeries" in INDEX_HTML
    assert "isPositiveWholeNumber(nTreat?.value) && fitsSeries" in INDEX_HTML

    # And nothing prefills them: a guessed group size is not an improvement.
    assert 'class="n-ctrl" type="number" min="1" step="1" inputmode="numeric" value="${esc(nCtrl)}"' in INDEX_HTML


def test_the_stop_button_is_wired_to_something_that_exists() -> None:
    """A button rendered into innerHTML has no handler until delegation finds it.

    Both progress cards are rebuilt on every poll tick, so the buttons cannot
    carry their own listener; they are reached through the container delegation.
    A mismatch between the rendered id and the delegated selector would produce a
    button that looks live and does nothing, which is worse than no button.
    """

    from degora.api import INDEX_HTML

    for button_id, container in (
        ("cancelSearchJob", '$("discoveryResults").addEventListener("click"'),
        ("cancelPrepareJob", '$("preparedCandidates").addEventListener("click"'),
    ):
        assert f'id="{button_id}"' in INDEX_HTML, f"{button_id} is never rendered"
        assert f'closest("#{button_id}")' in INDEX_HTML, f"{button_id} has no delegated handler"
        # The handler must be inside the container that actually holds the button.
        handler_at = INDEX_HTML.index(container)
        delegated_at = INDEX_HTML.index(f'closest("#{button_id}")')
        assert 0 < delegated_at - handler_at < 400, f"{button_id} is delegated from the wrong container"

    # Both call the one cancel function, and it is defined.
    assert "async function cancelDiscoveryJob(kind)" in INDEX_HTML
    assert INDEX_HTML.count('cancelDiscoveryJob("search")') == 1
    assert INDEX_HTML.count('cancelDiscoveryJob("prepare")') == 1


def test_a_cancelled_job_ends_both_poll_loops() -> None:
    """Neither loop may treat "cancelled" as a failure or keep polling forever."""

    from degora.api import INDEX_HTML

    # Search, prepare and (since 0.4.31) analysis each poll a job.
    assert INDEX_HTML.count('job.status === "cancelled"') == 3
    # The cancelled check precedes the failure check in every loop, so a
    # cancellation is never reported to the reader as an error.
    for failure in ('throw new Error(job.error || "publication search failed")',
                    'throw new Error(job.error || "preparation failed")',
                    'throw new Error(job.error || "analysis failed")'):
        failure_at = INDEX_HTML.index(failure)
        cancelled_at = INDEX_HTML.rindex('job.status === "cancelled"', 0, failure_at)
        assert failure_at - cancelled_at < 400


def test_a_too_late_cancel_adopts_the_completed_browser_result() -> None:
    from degora.api import INDEX_HTML

    false_branch = INDEX_HTML.index("outcome && outcome.cancelled === false")
    retire_branch = INDEX_HTML.index("Retire the in-flight work client-side too")
    false_body = INDEX_HTML[false_branch:retire_branch]
    assert false_branch < retire_branch
    assert "await refreshSearchPage(species, adoptRequest);" in INDEX_HTML
    assert "Search completed, but the saved result could not be loaded" in INDEX_HTML
    assert "state.prepared = job.result;" in INDEX_HTML
    assert "state.notice = notice;" in false_body
    assert "return;" in false_body
    assert "state.jobCancelled = true;" not in false_body
    assert "state.prepareCancelled = true;" not in false_body


def test_a_stopped_search_is_not_reported_as_an_empty_one() -> None:
    """"No records were returned" is a claim about the query, not about stopping."""

    from degora.api import INDEX_HTML

    assert "You stopped this search, so it has no results." in INDEX_HTML
    assert "This is not a finding about the query." in INDEX_HTML
