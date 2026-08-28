"""Focused branch coverage for the v0.4.34 durable discovery-state hardening."""

from __future__ import annotations

import concurrent.futures
import json
import sqlite3
import threading
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from degora.discovery_store import (
    DEFAULT_DB_NAME,
    DiscoveryJobCancelled,
    DiscoveryJobManager,
    DiscoveryQueueFullError,
    DiscoveryStateStore,
    DiscoveryStoreError,
    _db_path,
    _json_loads,
    _sanitize_url,
    sanitize_discovery_payload,
)


def test_url_and_json_sanitizer_covers_all_credential_locations() -> None:
    payload = {
        "userinfo": "https://alice:PASS@example.test/data?id=GSE1",
        "query": "https://example.test/data?api-key=QUERYSECRET&id=GSE2",
        "oauth_fragment": "https://example.test/callback#access_token=OAUTHSECRET&state=SAFE",
        "router_fragment": (
            "https://example.test/app#/callback?X-Amz-Security-Token=AWSSECRET&record_id=PMID1"
        ),
        "service authorization": "AUTHSECRET",
        "build_token": "BUILDSECRET",
        "api key": "APISECRET",
        "next_token": "SAFE_NEXT",
        "page_token": "SAFE_PAGE",
        "token_count": 3,
        "nested": ("https://example.test/data?signature=SIGSECRET&accession=GSE3", 7),
        "https://example.test/key?token=KEYSECRET": "safe value",
    }

    cleaned = sanitize_discovery_payload(payload)
    encoded = json.dumps(cleaned, sort_keys=True)

    for secret in (
        "PASS",
        "QUERYSECRET",
        "OAUTHSECRET",
        "AWSSECRET",
        "AUTHSECRET",
        "BUILDSECRET",
        "APISECRET",
        "SIGSECRET",
        "KEYSECRET",
    ):
        assert secret not in encoded
    assert "[redacted]@example.test" in cleaned["userinfo"]
    assert parse_qs(urlsplit(cleaned["query"]).query) == {
        "api-key": ["[redacted]"],
        "id": ["GSE2"],
    }
    assert parse_qs(urlsplit(cleaned["oauth_fragment"]).fragment) == {
        "access_token": ["[redacted]"],
        "state": ["SAFE"],
    }
    route, route_query = urlsplit(cleaned["router_fragment"]).fragment.split("?", 1)
    assert route == "/callback"
    assert parse_qs(route_query) == {
        "X-Amz-Security-Token": ["[redacted]"],
        "record_id": ["PMID1"],
    }
    assert cleaned["service authorization"] == "[redacted]"
    assert cleaned["build_token"] == "[redacted]"
    assert cleaned["api key"] == "[redacted]"
    assert cleaned["next_token"] == "SAFE_NEXT"
    assert cleaned["page_token"] == "SAFE_PAGE"
    assert cleaned["token_count"] == 3
    assert isinstance(cleaned["nested"], list)
    assert cleaned["nested"][1] == 7
    assert any("%5Bredacted%5D" in key for key in cleaned if key.startswith("https://"))


def test_url_sanitizer_fails_closed_without_damaging_non_urls() -> None:
    malformed = "http://[broken"
    malformed_secret = "http://[broken?token=SECRET"
    no_authority = "https:///local/path?token=SECRET"
    unsupported = "mailto:reader@example.test?token=ordinary-text"
    safe = "https://example.test/data?id=GSE1&page_token=NEXT"

    assert _sanitize_url(malformed) == malformed
    assert _sanitize_url(malformed_secret) == "[redacted: credential URL]"
    assert _sanitize_url(no_authority) == "[redacted: credential URL]"
    assert _sanitize_url(unsupported) == unsupported
    assert _sanitize_url(safe) == safe


def test_strict_json_failures_roll_back_and_leave_the_store_reusable(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state")

    with pytest.raises(DiscoveryStoreError, match="JSON-serializable"):
        store.create_job("search", {"score": float("nan")})
    assert store.list_jobs() == []

    valid = store.create_job("search", {"score": 1.0})
    with pytest.raises(DiscoveryStoreError, match="JSON-serializable"):
        store.update_job(valid["job_id"], status="running", result={"bad": object()})
    assert store.get_job(valid["job_id"])["status"] == "queued"

    with pytest.raises(DiscoveryStoreError, match="JSON-serializable"):
        store.save_search("1" * 16, {"score": float("inf")})
    assert store.get_search("1" * 16) is None
    store.save_search("1" * 16, {"score": 1.0})

    with pytest.raises(DiscoveryStoreError, match="JSON-serializable"):
        store.save_artifact("table", "2" * 16, {"payload": object()})
    assert store.get_artifact("table", "2" * 16) is None
    store.save_artifact("table", "2" * 16, {"rows": []})

    assert store.get_search("1" * 16) == {"score": 1.0}
    assert store.get_artifact("table", "2" * 16) == {"rows": []}
    assert _json_loads(None) is None


def test_database_path_resolution_and_schema_version_rejection(tmp_path: Path) -> None:
    directory = tmp_path / "existing"
    directory.mkdir()
    explicit = tmp_path / "explicit.sqlite3"
    implicit = tmp_path / "implicit"

    assert _db_path(directory) == directory / DEFAULT_DB_NAME
    assert _db_path(explicit) == explicit
    assert _db_path(implicit) == implicit / DEFAULT_DB_NAME

    store = DiscoveryStateStore(explicit)
    with sqlite3.connect(store.db_path) as connection:
        connection.execute("UPDATE metadata SET value='999' WHERE key='schema_version'")

    with pytest.raises(DiscoveryStoreError, match="unsupported discovery store schema version: 999"):
        DiscoveryStateStore(explicit)


@pytest.mark.parametrize("limit", [True, 0, -1, 1.5])
def test_create_job_rejects_invalid_active_limits(tmp_path: Path, limit: object) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")

    with pytest.raises(DiscoveryStoreError, match="active_job_limit must be a positive integer"):
        store.create_job("search", {}, active_job_limit=limit)  # type: ignore[arg-type]


def test_store_validation_and_lifecycle_boundaries(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")

    with pytest.raises(DiscoveryStoreError, match="kind is required"):
        store.create_job(" ", {})
    with pytest.raises(DiscoveryStoreError, match="16-character lowercase hex"):
        store.create_job("search", {}, job_id="ABC")
    with pytest.raises(DiscoveryStoreError, match="kind is required"):
        store.save_artifact(" ", "1" * 16, {})
    with pytest.raises(DiscoveryStoreError, match="limit must be a positive integer"):
        store.list_jobs(limit=False)

    job = store.create_job("search", {}, job_id="a" * 16)
    assert store.get_job("b" * 16) is None
    assert store.get_search("c" * 16) is None
    assert store.get_artifact("table", "d" * 16) is None

    with pytest.raises(KeyError, match="unknown discovery job"):
        store.update_job("b" * 16, progress=0.1)
    with pytest.raises(DiscoveryStoreError, match="invalid job status"):
        store.update_job(job["job_id"], status="unknown")
    with pytest.raises(DiscoveryStoreError, match="queued -> completed"):
        store.update_job(job["job_id"], status="completed")
    with pytest.raises(DiscoveryStoreError, match="between 0.0 and 1.0"):
        store.update_job(job["job_id"], progress=1.1)

    running = store.update_job(job["job_id"], status="running", progress=0.4)
    assert running["started_at"] is not None
    with pytest.raises(DiscoveryStoreError, match="cannot decrease"):
        store.update_job(job["job_id"], progress=0.3)
    completed = store.update_job(job["job_id"], status="completed", result={"ok": True})
    assert completed["progress"] == 1.0
    assert completed["completed_at"] is not None


def test_concurrent_store_admission_never_exceeds_the_durable_limit(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    contender_count = 8
    limit = 3
    gate = threading.Barrier(contender_count)

    def submit(index: int) -> tuple[str, str]:
        gate.wait(timeout=5)
        job_id = f"{index:016x}"
        try:
            store.create_job("search", {"index": index}, job_id=job_id, active_job_limit=limit)
        except DiscoveryQueueFullError:
            return "full", job_id
        return "admitted", job_id

    with concurrent.futures.ThreadPoolExecutor(max_workers=contender_count) as executor:
        outcomes = list(executor.map(submit, range(contender_count)))

    admitted = [job_id for outcome, job_id in outcomes if outcome == "admitted"]
    refused = [job_id for outcome, job_id in outcomes if outcome == "full"]
    assert len(admitted) == limit
    assert len(refused) == contender_count - limit
    assert len(store.list_jobs()) == limit

    for job_id in admitted:
        store.update_job(job_id, status="running")
        store.update_job(job_id, status="completed", result={})
    replacement = store.create_job("search", {}, active_job_limit=limit)
    assert replacement["status"] == "queued"


def test_interrupt_recovery_handles_malformed_search_links_and_exclusions(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    active_search_id = "1" * 16
    terminal_search_id = "2" * 16
    list_search_id = "3" * 16
    missing_search_id = "4" * 16
    excluded_search_id = "5" * 16
    store.save_search(active_search_id, {"status": "running", "records": []})
    store.save_search(terminal_search_id, {"status": "complete", "records": []})
    store.save_search(list_search_id, ["not", "a", "mapping"])
    store.save_search(excluded_search_id, {"status": "running", "records": []})

    jobs = {
        "active": store.create_job("publication_search", {"search_id": active_search_id}),
        "terminal": store.create_job("publication_search", {"search_id": terminal_search_id}),
        "list_search": store.create_job("publication_search", {"search_id": list_search_id}),
        "missing_search": store.create_job("publication_search", {"search_id": missing_search_id}),
        "bad_id": store.create_job("publication_search", {"search_id": "invalid"}),
        "missing_id": store.create_job("publication_search", {"query": "hypoxia"}),
        "list_payload": store.create_job("publication_search", ["invalid", "payload"]),
        "other_kind": store.create_job("prepare", {}),
        "excluded": store.create_job("publication_search", {"search_id": excluded_search_id}),
    }

    interrupted = store.interrupt_active_jobs(
        "server stopped",
        exclude_job_ids=[jobs["excluded"]["job_id"]],
    )

    assert {job["job_id"] for job in interrupted} == {
        job["job_id"] for name, job in jobs.items() if name != "excluded"
    }
    assert store.get_job(jobs["excluded"]["job_id"])["status"] == "queued"
    assert store.get_search(active_search_id)["status"] == "interrupted"
    assert store.get_search(active_search_id)["error"] == "server stopped"
    assert store.get_search(terminal_search_id)["status"] == "complete"
    assert store.get_search(list_search_id) == ["not", "a", "mapping"]
    assert store.get_search(missing_search_id) is None
    assert store.get_search(excluded_search_id)["status"] == "running"


def test_running_cancellation_retains_queue_capacity_until_worker_exit(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    manager = DiscoveryJobManager(store, max_workers=1, max_pending_jobs=1)
    first_started = threading.Event()
    release_first = threading.Event()

    def blocked_worker(_job_id, _payload, _progress):
        first_started.set()
        assert release_first.wait(timeout=5)
        return {"slot": 1}

    first = manager.submit("search", {"slot": 1}, blocked_worker)
    assert first_started.wait(timeout=5)
    cancelled = manager.cancel(first["job_id"])
    assert cancelled is not None and cancelled["status"] == "cancelled"

    with pytest.raises(DiscoveryQueueFullError, match="queue is full"):
        manager.submit(
            "search",
            {"slot": 2},
            lambda _job_id, payload, _progress: {"slot": payload["slot"]},
        )
    release_first.set()
    manager.shutdown(wait=True)

    assert store.get_job(first["job_id"])["result"] is None
    assert store.get_job(first["job_id"])["status"] == "cancelled"


@pytest.mark.parametrize("limit", [True, 0, -1, 1.5])
def test_manager_rejects_invalid_pending_limits(tmp_path: Path, limit: object) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")

    with pytest.raises(DiscoveryStoreError, match="max_pending_jobs must be a positive integer"):
        DiscoveryJobManager(store, max_pending_jobs=limit)  # type: ignore[arg-type]


def test_shutdown_racing_with_submission_interrupts_the_created_job(tmp_path: Path, monkeypatch) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    manager = DiscoveryJobManager(store, max_workers=1)
    original_create_job = store.create_job
    created_job_id: list[str] = []

    def create_then_close(*args, **kwargs):
        job = original_create_job(*args, **kwargs)
        created_job_id.append(job["job_id"])
        manager._closing.set()
        return job

    monkeypatch.setattr(store, "create_job", create_then_close)

    with pytest.raises(DiscoveryStoreError, match="manager is shutting down"):
        manager.submit("search", {}, lambda _job_id, _payload, _progress: {})

    assert len(created_job_id) == 1
    final = store.get_job(created_job_id[0])
    assert final is not None and final["status"] == "interrupted"
    manager.shutdown(wait=True)


def test_commit_fence_distinguishes_committing_cancelled_and_closing_states(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    manager = DiscoveryJobManager(store, max_workers=1)
    committing_id = "1" * 16
    cancelled_id = "2" * 16
    closing_id = "3" * 16

    manager.commit(committing_id)
    manager.commit(committing_id)
    assert manager._job_stop_reason(committing_id) == ""

    with manager._condition:
        manager._state[cancelled_id] = "cancelled"
    assert "reader" in manager._job_stop_reason(cancelled_id)
    with pytest.raises(DiscoveryJobCancelled, match="before its result was published"):
        manager.commit(cancelled_id)

    manager._closing.set()
    assert "stopping" in manager._job_stop_reason(closing_id)
    with pytest.raises(DiscoveryJobCancelled, match="before its result was published"):
        manager.commit(closing_id)
    manager.shutdown(wait=True)


def test_worker_failure_persistence_redacts_credential_bearing_url(tmp_path: Path) -> None:
    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(_job_id, _payload, _progress):
        raise RuntimeError(
            "provider rejected https://alice:PASS@example.test/data?access_token=TOKENSECRET"
            "&id=GSE1#session=SESSIONSECRET"
        )

    job = manager.submit("search", {}, worker)
    manager.shutdown(wait=True)

    failed = store.get_job(job["job_id"])
    assert failed is not None and failed["status"] == "failed"
    message = failed["error"]["message"]
    assert "PASS" not in message
    assert "TOKENSECRET" not in message
    assert "SESSIONSECRET" not in message
    assert "id=GSE1" in message
    assert "[redacted]@example.test" in message
