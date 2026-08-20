from __future__ import annotations

import json
import http.client
import sqlite3
import threading
import time
import urllib.request
from urllib.error import HTTPError

import pandas as pd
import pytest

from degora import __version__
import degora.api as api
from degora.api import (
    LOOPBACK_HOSTS,
    DiscoveryWorkspaceInUseError,
    create_server,
    redact_token_from_url_text,
    serve,
    strip_token_query_param,
)
from degora.provenance import EXTERNAL_PATH_PREFIX
from degora.score_db import write_score_database


def test_local_discovery_fallback_interrupts_running_search(tmp_path, monkeypatch) -> None:
    started = threading.Event()
    release = threading.Event()

    def blocked_search(query: str, species: str, *, limit: int) -> dict:
        started.set()
        release.wait(timeout=5)
        return {"query": query, "species": species, "limit": limit, "records": []}

    monkeypatch.setattr(api, "_call_search_publications", blocked_search)
    store = api._LocalDiscoveryStateStore(tmp_path / "discovery.sqlite3")
    manager = api._LocalDiscoveryJobManager(store, max_workers=1)
    search = store.create_search(query="hypoxia", species="human", limit=20)
    job = store.create_job(search_id=search["id"])
    manager.submit_search(
        job_id=job["id"],
        search_id=search["id"],
        query="hypoxia",
        species="human",
        limit=20,
    )
    assert started.wait(timeout=2)

    before = time.monotonic()
    manager.close(wait=False, cancel_futures=True, interrupt=True)
    assert time.monotonic() - before < 0.5
    assert store.get_job(job["id"])["status"] == "interrupted"
    assert store.get_search(search["id"])["status"] == "interrupted"

    release.set()
    manager.close(wait=True)
    assert store.get_job(job["id"])["status"] == "interrupted"
    assert store.get_search(search["id"])["status"] == "interrupted"


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


def test_local_api_serves_health_gene_list_and_detail(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

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
    assert "20 per page · up to 1,000 ranked" in html
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
    assert "independent biological replicates" in html
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
    assert "void loadDiscoveryDefaults();" in html
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


def test_serve_handles_keyboard_interrupt_cleanly(tmp_path, monkeypatch, capsys) -> None:
    db = tmp_path / "degora_scores.db"
    db.write_bytes(b"stub")
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
    db.write_bytes(b"stub")
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
    from degora.api import _redact_meta_for_network

    meta = _redact_meta_for_network(
        {
            "db_path": "/mnt/c/Projects/DEGORA/outputs/degora_scores.db",
            "output_dir": "C:\\Projects\\DEGORA\\outputs",
            "source_url": "https://example.test/data",
        }
    )

    assert meta["db_path"] == "[redacted: local path]"
    assert meta["output_dir"] == "[redacted: local path]"
    assert meta["source_url"] == "https://example.test/data"


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
    with pytest.raises(FileNotFoundError):
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
