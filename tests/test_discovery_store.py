from __future__ import annotations

import concurrent.futures
import sqlite3
import threading
import time

import pytest

from degora.discovery_store import DiscoveryJobManager, DiscoveryStateStore, DiscoveryStoreError


def test_store_persists_jobs_across_reopen(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    job = store.create_job("search", {"query": "hypoxia", "api_key": "do-not-store"})
    store.update_job(job["job_id"], status="running", progress=0.25, message="started")
    store.update_job(job["job_id"], status="completed", result={"hits": [1, 2]})

    reopened = DiscoveryStateStore(tmp_path)
    recovered = reopened.get_job(job["job_id"])

    assert recovered is not None
    assert recovered["status"] == "completed"
    assert recovered["progress"] == 1.0
    assert recovered["payload"]["api_key"] == "[redacted]"
    assert recovered["result"] == {"hits": [1, 2]}


def test_store_uses_schema_metadata_and_wal_when_available(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)

    with sqlite3.connect(store.db_path) as connection:
        schema_version = connection.execute("SELECT value FROM metadata WHERE key='schema_version'").fetchone()[0]
        journal_mode = connection.execute("PRAGMA journal_mode").fetchone()[0]

    assert schema_version == "1"
    assert journal_mode in {"wal", "delete", "truncate", "persist", "memory", "off"}


def test_concurrent_manager_submissions_complete(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=4)

    def worker(job_id, payload, progress):
        progress(0.5, f"half {job_id}")
        return {"value": payload["value"]}

    jobs = [manager.submit("search", {"value": index}, worker) for index in range(12)]
    manager.shutdown()

    completed = [store.get_job(job["job_id"]) for job in jobs]
    assert {job["status"] for job in completed if job is not None} == {"completed"}
    assert sorted(job["result"]["value"] for job in completed if job is not None) == list(range(12))


def test_direct_concurrent_create_job_ids_are_unique(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)

    def create(index: int) -> str:
        return store.create_job("search", {"index": index})["job_id"]

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        job_ids = list(executor.map(create, range(40)))

    assert len(job_ids) == len(set(job_ids))
    assert len(store.list_jobs(limit=100)) == 40


def test_manager_records_worker_failure(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(job_id, payload, progress):
        progress(0.2, "before failure")
        raise RuntimeError("boom")

    job = manager.submit("search", {"query": "x"}, worker)
    manager.shutdown()

    failed = store.get_job(job["job_id"])
    assert failed is not None
    assert failed["status"] == "failed"
    assert failed["progress"] == 0.2
    assert failed["error"] == {"message": "boom", "type": "RuntimeError"}


def test_progress_must_be_monotonic(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    job = store.create_job("search", {})
    store.update_job(job["job_id"], status="running", progress=0.4)

    with pytest.raises(DiscoveryStoreError, match="progress cannot decrease"):
        store.update_job(job["job_id"], progress=0.3)


def test_invalid_transition_is_rejected(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    job = store.create_job("search", {})

    with pytest.raises(DiscoveryStoreError, match="queued -> completed"):
        store.update_job(job["job_id"], status="completed")

    assert store.get_job(job["job_id"])["status"] == "queued"


def test_recover_interrupted_jobs_preserves_completed_and_artifacts(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    queued = store.create_job("search", {"state": "queued"})
    running = store.create_job("search", {"state": "running"})
    completed = store.create_job("search", {"state": "completed"})
    store.update_job(running["job_id"], status="running", progress=0.6, message="working")
    store.update_job(completed["job_id"], status="running", progress=0.5)
    store.update_job(completed["job_id"], status="completed", result={"ok": True})
    store.save_artifact("bundle", "aaaaaaaaaaaaaaaa", {"path": "bundle.json"})

    recovered = store.recover_interrupted_jobs()

    assert {job["job_id"] for job in recovered} == {queued["job_id"], running["job_id"]}
    assert store.get_job(queued["job_id"])["status"] == "interrupted"
    running_after = store.get_job(running["job_id"])
    assert running_after["status"] == "interrupted"
    assert "interrupted before completion" in running_after["message"]
    assert store.get_job(completed["job_id"])["status"] == "completed"
    assert store.get_artifact("bundle", "aaaaaaaaaaaaaaaa") == {"path": "bundle.json"}


def test_search_and_artifact_roundtrip(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    store.save_search("1111111111111111", {"query": "hypoxia", "items": [1, 2]})
    store.save_artifact("table", "2222222222222222", {"rows": [{"gene": "VEGFA"}]})

    reopened = DiscoveryStateStore(tmp_path)

    assert reopened.get_search("1111111111111111") == {"items": [1, 2], "query": "hypoxia"}
    assert reopened.get_artifact("table", "2222222222222222") == {"rows": [{"gene": "VEGFA"}]}


def test_shutdown_waits_for_running_job(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(job_id, payload, progress):
        time.sleep(0.02)
        progress(0.7, "almost")
        return {"done": True}

    job = manager.submit("search", {}, worker)
    manager.shutdown()

    assert store.get_job(job["job_id"])["status"] == "completed"


def test_interrupting_shutdown_returns_promptly_and_preserves_terminal_state(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    release = threading.Event()

    def worker(job_id, payload, progress):
        started.set()
        release.wait(timeout=5)
        progress(0.8, "late progress")
        return {"done": True}

    job = manager.submit("search", {}, worker)
    assert started.wait(timeout=2)

    before = time.monotonic()
    manager.shutdown(wait=False, cancel_futures=True, interrupt=True)
    elapsed = time.monotonic() - before

    assert elapsed < 0.5
    interrupted = store.get_job(job["job_id"])
    assert interrupted is not None
    assert interrupted["status"] == "interrupted"
    assert "server was stopped" in interrupted["message"]

    release.set()
    manager.shutdown(wait=True)
    assert store.get_job(job["job_id"])["status"] == "interrupted"
