from __future__ import annotations

import json
import sqlite3
import threading
import urllib.request
from urllib.error import HTTPError

import pandas as pd
import pytest

from degora import __version__
import degora.api as api
from degora.api import LOCAL_PATH_REDACTION, LOOPBACK_HOSTS, create_server, serve
from degora.score_db import write_score_database


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
    assert "const esc =" in html
    assert 'id="layoutSplitter"' in html
    assert 'role="separator"' in html
    assert "radial-gradient" in html
    assert "18px minmax" in html
    assert "initPanelResize();" in html
    assert "JSON.parse(message)" in html
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

    assert meta["db_path"] == "[redacted: local path]"
    assert meta["harmonized_path"] == "[redacted: local path]"
    assert meta["output_dir"] == "[redacted: local path]"
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
    secret_source = tmp_path / "Users" / "keunsoo" / "Projects_main" / "source.csv"
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
    assert all(study["source_path"] == LOCAL_PATH_REDACTION for study in studies)
    assert all(row["source_path"] == LOCAL_PATH_REDACTION for row in detail["evidence"])
    assert all(row["contributing_source_paths"] == LOCAL_PATH_REDACTION for row in detail["evidence"])


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


def test_create_server_falls_back_when_port_in_use(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")
    db = tmp_path / "degora_scores.db"

    first = create_server(db, port=0, quiet=True)
    try:
        busy_port = first.server_address[1]
        # Requesting the already-bound port must not raise; it auto-falls-back.
        second = create_server(db, port=busy_port, quiet=True)
        try:
            assert second.server_address[1] != busy_port
            assert second.server_address[1] != 0
        finally:
            second.server_close()
    finally:
        first.server_close()
