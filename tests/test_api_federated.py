from __future__ import annotations

import io
import json
import re
import sys
import threading
import time
import types
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from openpyxl import load_workbook
import pytest

import degora.discovery as discovery
from degora.api import INDEX_HTML, create_server


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


def _start_server(tmp_path: Path):
    db = tmp_path / "degora.db"
    db.touch()
    server = create_server(db, port=0, quiet=True, discovery_root=tmp_path / "discovery")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, thread, f"http://{host}:{port}"


def _stop_server(server, thread) -> None:
    server.shutdown()
    server.server_close()
    thread.join(timeout=5)


def _install_federated_module(monkeypatch, *, calls: list[dict]) -> None:
    module = types.ModuleType("degora.discovery_federated")

    def search_publications(*, query, species, limit):
        records = []
        for index in range(120):
            records.append(
                {
                    "publication_id": f"{species}-{index:03d}",
                    "title": "=formula risk" if index == 0 else f"{species} paper {index:03d}",
                    "authors": f"Author {index}",
                    "journal": "Journal",
                    "year": 2000 + (index % 20),
                    "readiness": "ready" if index % 2 == 0 else "candidate",
                    "relevance": 120 - index,
                    "data_sources": ["GEO", f"GSE{index:04d}"],
                    "candidate_routes": ["author_deg_table"],
                    "species_decision": species,
                }
            )
        return {
            "query": query,
            "species": {"key": species},
            "limit": limit,
            "records": records,
            "total": len(records),
            "provider_events": [{"provider": "fake", "event": "search", "status": "ok", "message": "+neutralize"}],
        }

    def page_publication_snapshot(snapshot, *, page, page_size, sort_by, sort_order):
        calls.append({"page": page, "page_size": page_size, "sort_by": sort_by, "sort_order": sort_order})
        records = list(snapshot["records"])
        reverse = sort_order == "desc"
        if sort_by == "year":
            records.sort(key=lambda record: record["year"], reverse=reverse)
        elif sort_by == "title":
            records.sort(key=lambda record: record["title"], reverse=reverse)
        elif sort_by == "relevance":
            records.sort(key=lambda record: record["relevance"], reverse=reverse)
        offset = (page - 1) * page_size
        return {
            "records": records[offset : offset + page_size],
            "total": len(records),
            "page": page,
            "page_size": page_size,
            "total_pages": 6,
            "sort_by": sort_by,
            "sort_order": sort_order,
            "has_next": page < 6,
        }

    module.search_publications = search_publications
    module.page_publication_snapshot = page_publication_snapshot
    module.resolve_publication_records = lambda records, species: list(records)
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)


def _create_search(base: str, species: str) -> dict:
    status, created = _request_json(
        f"{base}/api/discovery/searches",
        payload={"query": "hypoxia", "species": species, "limit": 120},
        action=True,
    )
    assert status == 202
    deadline = time.time() + 5
    while time.time() < deadline:
        _, job_payload = _request_json(f"{base}/api/discovery/jobs/{created['job_id']}")
        if job_payload["job"]["status"] == "complete":
            return created
        time.sleep(0.05)
    raise AssertionError("search job did not complete")


def test_federated_search_async_pagination_sort_restart_and_export(tmp_path: Path, monkeypatch) -> None:
    page_calls: list[dict] = []
    _install_federated_module(monkeypatch, calls=page_calls)
    server, thread, base = _start_server(tmp_path)
    try:
        created = _create_search(base, "human")
        _, search_payload = _request_json(f"{base}/api/discovery/searches/{created['search_id']}")
        _, first_page = _request_json(f"{base}/api/discovery/searches/{created['search_id']}/records?page=1&page_size=20")
        _, sorted_page = _request_json(
            f"{base}/api/discovery/searches/{created['search_id']}/records?"
            f"{urllib.parse.urlencode({'page': 1, 'page_size': 20, 'sort_by': 'year', 'sort_order': 'asc'})}"
        )
        status, headers, xlsx = _request_bytes(f"{base}/api/discovery/searches/{created['search_id']}/export.xlsx")
    finally:
        _stop_server(server, thread)

    assert search_payload["search"]["status"] == "complete"
    assert search_payload["search"]["total"] == 120
    assert len(first_page["records"]) == 20
    assert first_page["records"][0]["publication_id"] == "human-000"
    assert {"page": 1, "page_size": 20, "sort_by": "readiness", "sort_order": "desc"} in page_calls
    assert sorted_page["records"][0]["year"] == 2000
    assert {"page": 1, "page_size": 20, "sort_by": "year", "sort_order": "asc"} in page_calls
    assert status == 200
    assert headers["Content-Type"] == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

    workbook = load_workbook(io.BytesIO(xlsx), data_only=False)
    assert workbook.sheetnames == [
        "Query",
        "Publications",
        "Identifiers",
        "Linked datasets",
        "Candidate routes",
        "Species decisions",
        "Provider events",
    ]
    assert workbook["Publications"]["A2"].value == "'=formula risk"
    assert workbook["Provider events"]["D2"].value == "'+neutralize"

    restarted, restarted_thread, restarted_base = _start_server(tmp_path)
    try:
        _, persisted = _request_json(f"{restarted_base}/api/discovery/searches/{created['search_id']}")
        _, persisted_page = _request_json(
            f"{restarted_base}/api/discovery/searches/{created['search_id']}/records?page=6&page_size=20"
        )
    finally:
        _stop_server(restarted, restarted_thread)
    assert persisted["search"]["status"] == "complete"
    assert len(persisted_page["records"]) == 20


def test_federated_search_keeps_human_and_mouse_isolated(tmp_path: Path, monkeypatch) -> None:
    _install_federated_module(monkeypatch, calls=[])
    server, thread, base = _start_server(tmp_path)
    try:
        human = _create_search(base, "human")
        mouse = _create_search(base, "mouse")
        _, human_page = _request_json(f"{base}/api/discovery/searches/{human['search_id']}/records?page=1&page_size=20")
        _, mouse_page = _request_json(f"{base}/api/discovery/searches/{mouse['search_id']}/records?page=1&page_size=20")
    finally:
        _stop_server(server, thread)

    assert human["search_id"] != mouse["search_id"]
    assert human_page["records"][0]["publication_id"].startswith("human-")
    assert mouse_page["records"][0]["publication_id"].startswith("mouse-")
    assert "Search papers and linked repositories" in INDEX_HTML
    assert "Run separate Human + Mouse searches" in INDEX_HTML
    assert "Human and Mouse never pooled" in INDEX_HTML
    assert "Publication + linked public data discovery" in INDEX_HTML
    assert 'sort: { key: "readiness", order: "desc" }' in INDEX_HTML
    assert "Sort: DEG readiness · relevance tie-break" in INDEX_HTML
    assert 'id="downloadSearchExcel"' in INDEX_HTML
    assert '`/api/discovery/searches/${state.searchId}/export.xlsx`' in INDEX_HTML
    assert "complete snapshot" in INDEX_HTML
    assert "source-unit conflict" in INDEX_HTML
    assert 'class="mobile-field-label">DEG readiness' in INDEX_HTML
    assert 'role="tabpanel"' in INDEX_HTML


def test_federated_prepare_uses_persisted_canonical_selection_and_persists_bundle(tmp_path: Path, monkeypatch) -> None:
    _install_federated_module(monkeypatch, calls=[])
    captured: dict = {}
    prepare_module = types.ModuleType("degora.discovery_prepare")

    def fake_prepare(records, species, *, query, materialize_dir):
        captured.update(records=records, species=species, query=query, materialize_dir=str(materialize_dir))
        return {
            "species": {"key": species},
            "query": query,
            "studies": [],
            "excluded_studies": [],
            "returned_studies": 0,
            "materialize_dir": str(materialize_dir),
            "exports": {},
        }

    prepare_module.prepare_publication_records = fake_prepare
    monkeypatch.setitem(sys.modules, "degora.discovery_prepare", prepare_module)
    server, thread, base = _start_server(tmp_path)
    try:
        created = _create_search(base, "human")
        status, prepared = _request_json(
            f"{base}/api/discovery/prepare",
            payload={
                "species": "human",
                "query": "hypoxia",
                "search_id": created["search_id"],
                "record_ids": ["human-000"],
            },
            action=True,
        )
    finally:
        _stop_server(server, thread)

    assert status == 201
    assert captured["records"][0]["publication_id"] == "human-000"
    assert captured["species"] == "human"
    assert prepared["search_id"] == created["search_id"]

    restarted, restarted_thread, _base = _start_server(tmp_path)
    try:
        persisted = restarted.discovery_search_store.get_artifact("bundle", prepared["bundle_id"])
    finally:
        _stop_server(restarted, restarted_thread)
    assert persisted["bundle_id"] == prepared["bundle_id"]
    assert persisted["search_id"] == created["search_id"]


def test_federated_prepare_rejects_query_relabeling(tmp_path: Path, monkeypatch) -> None:
    _install_federated_module(monkeypatch, calls=[])
    server, thread, base = _start_server(tmp_path)
    try:
        created = _create_search(base, "human")
        request = urllib.request.Request(
            f"{base}/api/discovery/prepare",
            data=json.dumps(
                {
                    "species": "human",
                    "query": "unrelated label",
                    "search_id": created["search_id"],
                    "record_ids": ["human-000"],
                }
            ).encode(),
            headers={"Content-Type": "application/json", "X-DEGORA-Action": "1"},
        )
        try:
            urllib.request.urlopen(request, timeout=10)
        except urllib.error.HTTPError as exc:
            assert exc.code == 400
            payload = json.loads(exc.read().decode())
        else:  # pragma: no cover - the request must be rejected.
            raise AssertionError("query relabeling was accepted")
    finally:
        _stop_server(server, thread)

    assert "does not match" in payload["error"]


def test_legacy_discovery_search_endpoint_remains_synchronous_geo_wrapper(tmp_path: Path, monkeypatch) -> None:
    def fake_search_geo(query, species, **kwargs):
        return {
            "query": query,
            "species": {"key": species},
            "page": kwargs["page"],
            "page_size": kwargs["page_size"],
            "studies": [{"accession": "GSE1", "paper_title": "legacy"}],
        }

    monkeypatch.setattr(discovery, "search_geo", fake_search_geo)
    server, thread, base = _start_server(tmp_path)
    try:
        status, payload = _request_json(f"{base}/api/discovery/search?q=hypoxia&species=human&page=1", action=True)
    finally:
        _stop_server(server, thread)

    assert status == 200
    assert payload["studies"] == [{"accession": "GSE1", "paper_title": "legacy"}]


def test_discovery_workspace_has_a_single_live_server_owner(tmp_path: Path) -> None:
    db = tmp_path / "degora.db"
    db.touch()
    root = tmp_path / "discovery"
    first = create_server(db, port=0, quiet=True, discovery_root=root)
    try:
        with pytest.raises(RuntimeError, match="already using this discovery workspace"):
            create_server(db, port=0, quiet=True, discovery_root=root)
    finally:
        first.server_close()

    restarted = create_server(db, port=0, quiet=True, discovery_root=root)
    restarted.server_close()


def test_internal_type_error_is_not_retried_as_a_second_provider_call(tmp_path: Path, monkeypatch) -> None:
    calls: list[tuple[str, str, int]] = []
    module = types.ModuleType("degora.discovery_federated")

    def search_publications(*, query, species, limit):
        calls.append((query, species, limit))
        raise TypeError("provider implementation bug")

    module.search_publications = search_publications
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)
    server, thread, base = _start_server(tmp_path)
    try:
        status, created = _request_json(
            f"{base}/api/discovery/searches",
            payload={"query": "hypoxia", "species": "human", "limit": 10},
            action=True,
        )
        assert status == 202
        deadline = time.time() + 5
        while time.time() < deadline:
            _, payload = _request_json(f"{base}/api/discovery/jobs/{created['job_id']}")
            if payload["job"]["status"] == "failed":
                break
            time.sleep(0.05)
        else:
            raise AssertionError("failing search job did not reach terminal state")
    finally:
        _stop_server(server, thread)

    assert calls == [("hypoxia", "human", 10)]
    assert payload["job"]["error"] == "provider implementation bug"


def test_persistent_discovery_endpoints_require_token_and_action_header(tmp_path: Path, monkeypatch) -> None:
    _install_federated_module(monkeypatch, calls=[])
    db = tmp_path / "degora.db"
    db.touch()
    server = create_server(
        db,
        port=0,
        quiet=True,
        access_token="secret-token",
        discovery_root=tmp_path / "discovery",
    )
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base = f"http://{host}:{port}"
    endpoint = f"{base}/api/discovery/searches"
    try:
        with pytest.raises(urllib.error.HTTPError) as exc_info:
            _request_json(endpoint, payload={"query": "hypoxia", "species": "human", "limit": 20}, action=True)
        assert exc_info.value.code == 401

        with pytest.raises(urllib.error.HTTPError) as exc_info:
            _request_json(
                endpoint,
                payload={"query": "hypoxia", "species": "human", "limit": 20},
                token="secret-token",
            )
        assert exc_info.value.code == 400

        status, created = _request_json(
            endpoint,
            payload={"query": "hypoxia", "species": "human", "limit": 20},
            action=True,
            token="secret-token",
        )
        assert status == 202
        job_url = f"{base}/api/discovery/jobs/{created['job_id']}"
        with pytest.raises(urllib.error.HTTPError) as exc_info:
            _request_json(job_url)
        assert exc_info.value.code == 401
        deadline = time.time() + 5
        while time.time() < deadline:
            _, job = _request_json(job_url, token="secret-token")
            if job["job"]["status"] == "complete":
                break
            time.sleep(0.05)
        else:
            raise AssertionError("authorized discovery search did not complete")
        export_url = f"{base}/api/discovery/searches/{created['search_id']}/export.xlsx"
        with pytest.raises(urllib.error.HTTPError) as exc_info:
            _request_bytes(export_url)
        assert exc_info.value.code == 401
        assert _request_bytes(export_url, token="secret-token")[0] == 200
    finally:
        _stop_server(server, thread)


def _install_progress_reporting_module(monkeypatch, *, seen: list[tuple[float, str]]) -> None:
    """A stub that accepts the optional progress callback, unlike the strict one."""

    module = types.ModuleType("degora.discovery_federated")

    def search_publications(*, query, species, limit, progress=None):
        for fraction, message in ((0.25, "Querying stub (1 of 2 sources)"), (0.75, "Inspecting linked data 5 of 20")):
            if progress is not None:
                progress(fraction, message)
                seen.append((fraction, message))
            time.sleep(0.05)
        return {
            "query": query,
            "species": {"key": species},
            "limit": limit,
            "records": [{"publication_id": "p1", "title": "paper", "year": 2020}],
            "total": 1,
        }

    module.search_publications = search_publications
    module.page_publication_snapshot = lambda snapshot, **kwargs: {
        "records": list(snapshot["records"]),
        "total": 1,
        "page": 1,
        "page_size": 20,
        "total_pages": 1,
        "sort_by": "data_readiness",
        "sort_order": "desc",
        "has_next": False,
    }
    module.resolve_publication_records = lambda records, species: list(records)
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)


def test_discovery_job_exposes_progress_and_message(tmp_path: Path, monkeypatch) -> None:
    """The browser draws a determinate bar from these two fields."""

    reported: list[tuple[float, str]] = []
    _install_progress_reporting_module(monkeypatch, seen=reported)
    server, thread, base = _start_server(tmp_path)
    try:
        status, created = _request_json(
            f"{base}/api/discovery/searches",
            payload={"query": "hypoxia", "species": "human", "limit": 20},
            action=True,
        )
        assert status == 202
        observed: list[tuple[float, str]] = []
        deadline = time.time() + 10
        while time.time() < deadline:
            _, payload = _request_json(f"{base}/api/discovery/jobs/{created['job_id']}")
            job = payload["job"]
            assert "progress" in job and "message" in job
            if job["progress"] is not None:
                observed.append((job["progress"], job["message"]))
            if job["status"] == "complete":
                break
            time.sleep(0.02)
        else:  # pragma: no cover - only on a stalled job.
            raise AssertionError("search job did not complete")

        assert reported, "the optional progress callback was not forwarded to search_publications"
        assert observed, "no progress fractions reached the API"
        fractions = [fraction for fraction, _ in observed]
        assert fractions == sorted(fractions), f"job progress must never decrease: {fractions}"
        assert observed[-1][0] == 1.0
        assert observed[-1][1]
    finally:
        _stop_server(server, thread)


def test_discovery_job_progress_is_clamped_and_message_is_bounded(tmp_path: Path, monkeypatch) -> None:
    from degora.api import _api_job_progress, _clean_job_message

    assert _api_job_progress(1.4) == 1.0
    assert _api_job_progress(-2) == 0.0
    assert _api_job_progress(None) is None
    assert _api_job_progress("nope") is None
    assert _api_job_progress(True) is None
    assert _api_job_progress(float("nan")) is None
    assert _clean_job_message("  many\n\nspaces  ") == "many spaces"
    assert _clean_job_message(None) == ""
    assert len(_clean_job_message("x" * 500)) == 160


def test_strict_search_stub_without_progress_still_runs(tmp_path: Path, monkeypatch) -> None:
    """The forwarding guard must not break callers that take no progress kwarg."""

    calls: list[dict] = []
    _install_federated_module(monkeypatch, calls=calls)
    server, thread, base = _start_server(tmp_path)
    try:
        created = _create_search(base, "human")
        _, payload = _request_json(f"{base}/api/discovery/jobs/{created['job_id']}")
        assert payload["job"]["status"] == "complete"
        assert payload["job"]["progress"] == 1.0
    finally:
        _stop_server(server, thread)


def test_prepare_job_route_reports_progress_and_returns_the_bundle(tmp_path: Path, monkeypatch) -> None:
    """Preparation used to be a blocking request with no way to show stages."""

    calls: list[dict] = []
    _install_federated_module(monkeypatch, calls=calls)
    stages: list[tuple[float, str]] = []

    prepare_module = types.ModuleType("degora.discovery_prepare")

    def prepare_publication_records(records, species, *, query, materialize_dir, progress=None, **kwargs):
        for fraction, message in ((0.2, "Downloading 1 of 3"), (0.8, "Downloading 3 of 3")):
            if progress is not None:
                progress(fraction, message)
                stages.append((fraction, message))
            time.sleep(0.05)
        Path(materialize_dir).mkdir(parents=True, exist_ok=True)
        return {"bundle_id": "abc", "studies": [{"accession": "GSE1"}], "excluded_studies": []}

    prepare_module.prepare_publication_records = prepare_publication_records
    monkeypatch.setitem(sys.modules, "degora.discovery_prepare", prepare_module)

    server, thread, base = _start_server(tmp_path)
    try:
        created = _create_search(base, "human")
        _, records_page = _request_json(f"{base}/api/discovery/searches/{created['search_id']}/records?page=1&page_size=20")
        record_ids = [row["publication_id"] for row in records_page["records"][:2]]

        status, started = _request_json(
            f"{base}/api/discovery/prepare-jobs",
            payload={
                "species": "human",
                "query": "hypoxia",
                "search_id": created["search_id"],
                "record_ids": record_ids,
            },
            action=True,
        )
        assert status == 202
        assert started["job_id"]

        observed: list[float] = []
        deadline = time.time() + 15
        while time.time() < deadline:
            _, payload = _request_json(f"{base}/api/discovery/jobs/{started['job_id']}")
            job = payload["job"]
            if job["progress"] is not None:
                observed.append(job["progress"])
            if job["status"] == "complete":
                assert job["result"], "the prepared bundle must reach the browser through the job"
                # The server assigns the bundle id; the worker's value is replaced.
                assert re.fullmatch(r"[a-f0-9]{16}", job["result"]["bundle_id"])
                assert job["result"]["studies"] == [{"accession": "GSE1"}]
                break
            if job["status"] == "failed":
                raise AssertionError(job["error"])
            time.sleep(0.02)
        else:  # pragma: no cover - only on a stalled job.
            raise AssertionError("prepare job did not complete")

        assert stages, "the progress callback was not forwarded to prepare_publication_records"
        assert observed == sorted(observed), f"prepare progress must never decrease: {observed}"
        assert observed[-1] == 1.0
    finally:
        _stop_server(server, thread)


def test_job_result_is_withheld_until_the_job_completes(tmp_path: Path, monkeypatch) -> None:
    calls: list[dict] = []
    _install_federated_module(monkeypatch, calls=calls)
    server, thread, base = _start_server(tmp_path)
    try:
        status, created = _request_json(
            f"{base}/api/discovery/searches",
            payload={"query": "hypoxia", "species": "human", "limit": 20},
            action=True,
        )
        assert status == 202
        _, payload = _request_json(f"{base}/api/discovery/jobs/{created['job_id']}")
        job = payload["job"]
        if job["status"] != "complete":
            assert job["result"] is None
    finally:
        _stop_server(server, thread)
