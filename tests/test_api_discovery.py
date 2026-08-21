from __future__ import annotations

import json
import shutil
import threading
import urllib.parse
import urllib.request
from pathlib import Path
from urllib.error import HTTPError

import pandas as pd
import pytest

import degora.discovery as discovery
import degora.discovery_run as discovery_run
from degora.api import DegoraRequestHandler, INDEX_HTML, create_server
from degora.score_db import write_score_database


def _db(tmp_path: Path) -> Path:
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "gene_symbol": ["TP53", "TP53"],
            "lfc": [2.0, 1.5],
            "signed_z": [4.0, 3.0],
            "pvalue": [0.001, 0.01],
            "padj": [0.01, 0.02],
            "normalized_rank": [0.01, 0.02],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "pipeline": ["DESeq2", "edgeR"],
            "species": ["Homo sapiens", "Homo sapiens"],
            "source_path": ["a.csv", "b.csv"],
            "source_url": ["https://example.test/a", "https://example.test/b"],
        }
    )
    path = tmp_path / "harmonized.csv"
    harmonized.to_csv(path, index=False)
    db = tmp_path / "degora_scores.db"
    write_score_database(path, tmp_path, db_path=db)
    return db


def _request_json(
    url: str,
    *,
    payload: dict | None = None,
    action: bool = False,
    token: str | None = None,
) -> tuple[int, dict]:
    data = None
    headers = {}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    if action:
        headers["X-DEGORA-Action"] = "1"
    if token:
        headers["X-DEGORA-Token"] = token
    request = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(request, timeout=10) as response:
        return response.status, json.loads(response.read().decode())


def _request_bytes(url: str, *, token: str | None = None) -> tuple[int, dict[str, str], bytes]:
    headers = {"X-DEGORA-Token": token} if token else {}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=10) as response:
        return response.status, dict(response.headers.items()), response.read()


def test_file_download_treats_client_disconnect_as_normal(tmp_path: Path) -> None:
    workbook = tmp_path / "DEGORA_output.xlsx"
    workbook.write_bytes(b"PK\x03\x04synthetic-xlsx")

    class DisconnectingWriter:
        def write(self, data: bytes) -> int:
            raise BrokenPipeError("client closed the download")

    handler = object.__new__(DegoraRequestHandler)
    handler.wfile = DisconnectingWriter()
    handler.send_response = lambda *args, **kwargs: None
    handler.send_header = lambda *args, **kwargs: None
    handler.end_headers = lambda: None

    handler._send_file_download(
        workbook,
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        filename="DEGORA_output.xlsx",
    )


def test_discovery_api_search_prepare_analyze_and_run_query_are_species_scoped(tmp_path: Path, monkeypatch) -> None:
    db = _db(tmp_path)

    search_calls = []

    def fake_search(
        query,
        species,
        *,
        page,
        page_size,
        assess_files,
        global_rank,
        global_limit,
        sort_by,
        sort_order,
    ):
        search_calls.append(
            {
                "assess_files": assess_files,
                "global_rank": global_rank,
                "global_limit": global_limit,
                "sort_by": sort_by,
                "sort_order": sort_order,
            }
        )
        if species not in {"human", "mouse"}:
            raise discovery.DiscoveryError("species must be exactly human or mouse")
        return {
            "query": query,
            "species": {"key": species},
            "page": page,
            "page_size": page_size,
            "total_hits": 1,
            "total_pages": 1,
            "evaluated_studies": 1,
            "ranking_limit": global_limit,
            "ranking_truncated": False,
            "studies": [{"accession": "GSE1", "paper_title": "Paper"}],
        }

    def fake_prepare(accessions, species, *, query, materialize_dir):
        Path(materialize_dir).mkdir(parents=True)
        return {
            "query": query,
            "species": {"key": species},
            "materialize_dir": str(materialize_dir),
            "selected_accessions": accessions,
            "studies": [],
        }

    def fake_run(prepared, selections, output, *, species, min_studies):
        output = Path(output)
        target = output / "results" / f"degora_{species}_scores.db"
        target.parent.mkdir(parents=True)
        shutil.copy2(db, target)
        workbook = target.parent / "DEGORA_output.xlsx"
        workbook.write_bytes(b"PK\x03\x04synthetic-xlsx")
        return {
            "status": "complete",
            "species": {"key": species, "label": species.title(), "scientific_name": "Homo sapiens"},
            "db_path": str(target),
            "output_dir": str(output),
            "top_genes": ["TP53"],
            "excel_workbook": {"output": str(workbook)},
        }

    monkeypatch.setattr(discovery, "search_geo", fake_search)
    monkeypatch.setattr(discovery, "prepare_geo_studies", fake_prepare)
    monkeypatch.setattr(discovery_run, "run_discovery_analysis", fake_run)

    server = create_server(db, port=0, quiet=True, discovery_root=tmp_path / "discovery")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base = f"http://{host}:{port}"
    try:
        _, search = _request_json(
            f"{base}/api/discovery/search?q=hypoxia&species=human&page=2&sort=year&order=desc",
            action=True,
        )
        _, relevance = _request_json(
            f"{base}/api/discovery/search?q=hypoxia&species=human&page=1&sort=relevance",
            action=True,
        )
        status, prepared = _request_json(
            f"{base}/api/discovery/prepare",
            payload={"query": "hypoxia", "species": "human", "accessions": ["GSE1"]},
            action=True,
        )
        analyze_status, run = _request_json(
            f"{base}/api/discovery/analyze",
            payload={
                "bundle_id": prepared["bundle_id"],
                "species": "human",
                "selections": [{"candidate_id": "x"}],
            },
            action=True,
        )
        _, genes = _request_json(f"{base}/api/discovery/runs/{run['run_id']}/genes?limit=10")
        excel_status, excel_headers, excel_bytes = _request_bytes(
            f"{base}/api/discovery/runs/{run['run_id']}/export.xlsx"
        )
        with pytest.raises(HTTPError) as exc_info:
            _request_json(f"{base}/api/discovery/search?q=x&species=human&page=1")
        assert exc_info.value.code == 400
        with pytest.raises(HTTPError) as exc_info:
            _request_json(f"{base}/api/discovery/search?q=x&species=both&page=1", action=True)
        assert exc_info.value.code == 400
        with pytest.raises(HTTPError) as exc_info:
            _request_json(
                f"{base}/api/discovery/prepare",
                payload={"query": "x", "species": "human", "accessions": ["GSE1"]},
                action=False,
            )
        assert exc_info.value.code == 400
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert search["species"]["key"] == "human"
    assert search["page"] == 2
    assert search["page_size"] == 10
    assert search_calls[0] == {
        "assess_files": True,
        "global_rank": True,
        "global_limit": 1000,
        "sort_by": "year",
        "sort_order": "desc",
    }
    assert search_calls[1]["sort_by"] == "relevance"
    assert search_calls[1]["sort_order"] is None
    assert relevance["page"] == 1
    assert all(call["global_rank"] for call in search_calls)
    assert status == 201
    assert analyze_status == 201
    assert run["species"]["key"] == "human"
    assert genes["genes"][0]["gene_symbol"] == "TP53"
    assert excel_status == 200
    assert excel_headers["Content-Type"] == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    assert excel_headers["Content-Disposition"] == 'attachment; filename="DEGORA_output.xlsx"'
    assert excel_headers["Cache-Control"] == "no-store"
    assert excel_headers["X-Content-Type-Options"] == "nosniff"
    assert excel_bytes == b"PK\x03\x04synthetic-xlsx"


def test_discovery_excel_download_requires_token_and_rejects_symlink_escape(tmp_path: Path) -> None:
    db = _db(tmp_path)
    discovery_root = tmp_path / "discovery"
    run_id = "a" * 16
    run_root = discovery_root / "human" / "runs" / run_id
    results = run_root / "results"
    results.mkdir(parents=True)
    workbook = results / "DEGORA_output.xlsx"
    workbook.write_bytes(b"PK\x03\x04authorized")

    server = create_server(
        db,
        port=0,
        quiet=True,
        access_token="secret",
        discovery_root=discovery_root,
    )
    server.discovery_runs[run_id] = {
        "run_id": run_id,
        "species": {"key": "human"},
        "output_dir": str(run_root),
        "excel_workbook": {"output": str(workbook)},
    }
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    url = f"http://{host}:{port}/api/discovery/runs/{run_id}/export.xlsx"
    try:
        with pytest.raises(HTTPError) as exc_info:
            _request_bytes(url)
        assert exc_info.value.code == 401

        status, _, payload = _request_bytes(url, token="secret")
        assert status == 200
        assert payload == b"PK\x03\x04authorized"

        outside = tmp_path / "outside.xlsx"
        outside.write_bytes(b"PK\x03\x04outside")
        workbook.unlink()
        workbook.symlink_to(outside)
        with pytest.raises(HTTPError) as exc_info:
            _request_bytes(url, token="secret")
        assert exc_info.value.code == 404
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_discovery_post_is_forbidden_on_non_loopback_server(tmp_path: Path) -> None:
    db = _db(tmp_path)
    discovery_root = tmp_path / "discovery"
    server = create_server(
        db,
        host="0.0.0.0",
        port=0,
        quiet=True,
        access_token="secret",
        discovery_root=discovery_root,
    )
    run_id = "b" * 16
    run_root = discovery_root / "human" / "runs" / run_id
    results = run_root / "results"
    results.mkdir(parents=True)
    workbook = results / "DEGORA_output.xlsx"
    workbook.write_bytes(b"PK\x03\x04network-blocked")
    server.discovery_runs[run_id] = {
        "run_id": run_id,
        "species": {"key": "human"},
        "output_dir": str(run_root),
        "excel_workbook": {"output": str(workbook)},
        "nested": {"candidate_paths": [str(run_root / "private-input.csv")]},
    }
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    _, port = server.server_address
    try:
        requests = [
            (
                "/api/discovery/prepare",
                {"query": "x", "species": "human", "accessions": ["GSE1"]},
            ),
            (
                "/api/discovery/analyze",
                {"bundle_id": "0" * 16, "species": "human", "selections": []},
            ),
        ]
        for endpoint, payload in requests:
            with pytest.raises(HTTPError) as exc_info:
                _request_json(
                    f"http://127.0.0.1:{port}{endpoint}",
                    payload=payload,
                    action=True,
                    token="secret",
                )
            assert exc_info.value.code == 403
        with pytest.raises(HTTPError) as exc_info:
            _request_bytes(
                f"http://127.0.0.1:{port}/api/discovery/runs/{run_id}/export.xlsx",
                token="secret",
            )
        assert exc_info.value.code == 403
        summary_url = f"http://127.0.0.1:{port}/api/discovery/runs/{run_id}/summary"
        with pytest.raises(HTTPError) as exc_info:
            _request_json(summary_url)
        assert exc_info.value.code == 401
        status, summary = _request_json(summary_url, token="secret")
        assert status == 200
        assert str(tmp_path) not in json.dumps(summary)
        assert summary["output_dir"] == "[redacted: local path]"
        assert summary["nested"]["candidate_paths"] == ["[redacted: local path]"]
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_discovery_dashboard_author_activation_contract_is_exposed_and_serialized() -> None:
    required_author_fields = [
        "sheet_name",
        "gene_column",
        "lfc_column",
        "p_column",
        "padj_column",
        "column_mapping_confirmed",
        "adjusted_p_as_pvalue_confirmed",
        "lfc_scale_confirmed_log2",
        "row_filter_column",
        "row_filter_value",
        "row_filter_confirmed",
        "duplicate_gene_policy",
        "duplicate_gene_policy_confirmed",
        "assay_type",
        "pipeline",
    ]
    for field in required_author_fields:
        assert f"common.{field}" in INDEX_HTML

    for css_class in [
        "sheet-name",
        "gene-column",
        "lfc-column",
        "p-column",
        "padj-column",
        "column-mapping-confirmed",
        "adjusted-p-as-pvalue-confirmed",
        "lfc-scale-confirmed-log2",
        "row-filter-column",
        "row-filter-value",
        "row-filter-confirmed",
        "duplicate-gene-policy",
        "duplicate-gene-policy-confirmed",
        "assay-type",
        "pipeline",
    ]:
        assert css_class in INDEX_HTML

    for status in [
        "ready_for_review",
        "requires_column_mapping",
        "requires_lfc_confirmation",
        "requires_pvalue_mapping",
    ]:
        assert status in INDEX_HTML

    assert "data-activation-key" in INDEX_HTML
    assert "clone-author-candidate" in INDEX_HTML
    assert "state.cloneCounter" in INDEX_HTML
    assert "AUTHOR_REVIEWABLE_STATUSES.has(status)" in INDEX_HTML
    assert "status !== \"ready_for_review\" || authorMappingEdited(row)" in INDEX_HTML
    assert "pColumn && padjColumn && pColumn === padjColumn" in INDEX_HTML
    assert "rowFilterPairValid" in INDEX_HTML
    assert 'duplicateGenePolicy?.value !== "keep_first"' in INDEX_HTML
