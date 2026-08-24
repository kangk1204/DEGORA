from __future__ import annotations

import concurrent.futures
import sqlite3
import threading
import time

import pytest

from degora.discovery_store import (
    DiscoveryJobCancelled,
    DiscoveryJobManager,
    DiscoveryStateStore,
    DiscoveryStoreError,
)


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


def test_worker_failure_message_redacts_secret_and_absolute_path(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(job_id, payload, progress):
        progress(0.2, "before failure")
        raise RuntimeError("provider failed with api_key=SECRET123 at /private/tmp/catalog.csv")

    job = manager.submit("search", {"query": "x"}, worker)
    manager.shutdown()

    failed = store.get_job(job["job_id"])
    assert failed is not None
    assert failed["status"] == "failed"
    message = failed["error"]["message"]
    assert "SECRET123" not in message
    assert "/private/tmp/catalog.csv" not in message
    assert "api_key=[redacted]" in message
    assert "[redacted: local path]" in message


@pytest.mark.parametrize(
    "secret_text,path_text",
    [
        ("client_secret=SECRET123", "/mnt/private/catalog.csv"),
        ("Authorization: Bearer SECRET123", "file:///Users/reviewer/catalog.csv"),
        ("token=SECRET123", "/catalog.csv"),
        ("password=SECRET123", "file:///catalog.csv"),
        ("ncbi_api_key:SECRET123", r"C:\\Users\\reviewer\\catalog.csv"),
        ("access_token=SECRET123", r"\\\\server\\share\\catalog.csv"),
    ],
)
def test_worker_failure_redaction_handles_common_secret_and_path_shapes(
    tmp_path, secret_text: str, path_text: str
) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(job_id, payload, progress):
        raise RuntimeError(f"provider failed with {secret_text} at {path_text}")

    job = manager.submit("search", {"query": "x"}, worker)
    manager.shutdown()

    message = store.get_job(job["job_id"])["error"]["message"]
    assert "SECRET123" not in message
    assert path_text not in message
    assert "[redacted]" in message
    assert "[redacted: local path]" in message


def test_artifact_paths_are_preserved_for_restart_and_local_reuse(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    materialize_dir = tmp_path / "prepared" / "bundle"
    source_path = materialize_dir / "study.csv"
    payload = {
        "materialize_dir": str(materialize_dir),
        "studies": [{"source_path": str(source_path)}],
    }

    store.save_artifact("bundle", "a" * 16, payload)
    reopened = DiscoveryStateStore(tmp_path / "state")

    assert reopened.get_artifact("bundle", "a" * 16) == payload


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


def test_shutdown_drains_committing_job_before_interrupting_others(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=2)
    committing_ready = threading.Event()
    running_ready = threading.Event()
    release_committing = threading.Event()

    def committing_worker(job_id, payload, progress):
        progress(0.6, "Ready to publish.")
        manager.commit(job_id)
        committing_ready.set()
        release_committing.wait(timeout=5)
        progress(0.9, "Publishing.")
        return {"done": "committed"}

    def running_worker(job_id, payload, progress):
        running_ready.set()
        threading.Event().wait(timeout=5)
        return {"done": "running"}

    committed_job = manager.submit("search", {"kind": "committing"}, committing_worker)
    running_job = manager.submit("search", {"kind": "running"}, running_worker)
    assert committing_ready.wait(timeout=5)
    assert running_ready.wait(timeout=5)

    shutdown_done = threading.Event()
    shutdown_thread = threading.Thread(
        target=lambda: (manager.shutdown(wait=False, cancel_futures=True, interrupt=True), shutdown_done.set())
    )
    shutdown_thread.start()

    assert not shutdown_done.wait(timeout=0.1)
    release_committing.set()
    assert shutdown_done.wait(timeout=5)
    manager.shutdown(wait=True)

    assert store.get_job(committed_job["job_id"])["status"] == "completed"
    assert store.get_job(running_job["job_id"])["status"] == "interrupted"


def test_shutdown_preserves_committing_job_after_drain_timeout(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)
    manager._commit_drain_timeout = 0.01
    committed = threading.Event()
    release = threading.Event()

    def worker(job_id, payload, progress):
        progress(0.6, "Ready to publish.")
        manager.commit(job_id)
        committed.set()
        release.wait(timeout=5)
        return {"done": True}

    job = manager.submit("search", {}, worker)
    assert committed.wait(timeout=5)

    manager.shutdown(wait=False, cancel_futures=True, interrupt=True)

    assert store.get_job(job["job_id"])["status"] == "running"
    release.set()
    manager.shutdown(wait=True)
    final = store.get_job(job["job_id"])
    assert final["status"] == "completed"
    assert final["result"] == {"done": True}


def test_submit_after_shutdown_does_not_create_queued_orphan(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)
    manager.shutdown(wait=False, cancel_futures=True, interrupt=True)

    with pytest.raises(DiscoveryStoreError, match="shutting down"):
        manager.submit("search", {}, lambda job_id, payload, progress: {"ok": True})

    assert store.list_jobs() == []


def test_submit_scheduling_failure_interrupts_created_job(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)
    manager.shutdown(wait=True)

    with pytest.raises(RuntimeError):
        manager.submit(
            "search",
            {},
            lambda job_id, payload, progress: {"ok": True},
            job_id="1234567890abcdef",
        )

    final = store.get_job("1234567890abcdef")
    assert final is not None
    assert final["status"] == "interrupted"


def test_cancelled_queued_future_cleanup_when_executor_drops_it(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    release = threading.Event()

    def blocking_worker(job_id, payload, progress):
        started.set()
        release.wait(timeout=5)
        return {"ok": True}

    first = manager.submit("search", {"slot": 1}, blocking_worker)
    second = manager.submit("search", {"slot": 2}, lambda job_id, payload, progress: {"ok": True})
    assert started.wait(timeout=5)

    assert manager.cancel(second["job_id"]) is not None
    manager.shutdown(wait=False, cancel_futures=True)
    release.set()
    manager.shutdown(wait=True)

    assert store.get_job(first["job_id"])["status"] == "completed"
    assert store.get_job(second["job_id"])["status"] == "cancelled"
    assert manager._state == {}


def test_custom_reused_job_id_completes_and_cleans_manager_state(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path)
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(job_id, payload, progress):
        progress(0.3, "Working.")
        manager.commit(job_id)
        return {"ok": True}

    job = manager.submit("search", {}, worker, job_id="abcdef1234567890")
    manager.shutdown(wait=True)

    assert job["job_id"] == "abcdef1234567890"
    assert store.get_job(job["job_id"])["status"] == "completed"
    assert manager._state == {}


def test_a_cancelled_worker_stops_instead_of_writing_after_shutdown(tmp_path) -> None:
    """Shutdown has to reach a worker that is already running, not just the queue.

    cancel_futures only drops work that has not started. A prepare already
    downloading kept going after the server was stopped and the job was recorded
    as interrupted, so files appeared in the workspace after the run that was
    supposed to own them had ended. Every stage already reports progress, so that
    is where the worker now notices.
    """

    store = DiscoveryStateStore(tmp_path / "discovery.sqlite3")
    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    release = threading.Event()
    side_effect = tmp_path / "written-after-shutdown.txt"
    outcome: list[str] = []

    def worker(_job_id, _payload, progress):
        progress(0.1, "stage one")
        started.set()
        release.wait(timeout=5)
        try:
            progress(0.5, "stage two")
        except BaseException as exc:  # noqa: BLE001 - record what reached the worker
            outcome.append(type(exc).__name__)
            raise
        side_effect.write_text("worker kept going", encoding="utf-8")
        outcome.append("completed")
        return {"ok": True}

    job = manager.submit("publication_prepare", {}, worker)
    assert started.wait(timeout=2)

    manager.shutdown(wait=False, cancel_futures=True, interrupt=True)
    assert store.get_job(job["job_id"])["status"] == "interrupted"

    release.set()
    for _ in range(50):
        if outcome:
            break
        time.sleep(0.05)

    assert outcome == [DiscoveryJobCancelled.__name__], outcome
    assert not side_effect.exists(), "a cancelled worker still wrote to the workspace"
    assert store.get_job(job["job_id"])["status"] == "interrupted"
