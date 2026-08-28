"""Regressions for bounded cancellation, URL redaction and commit draining."""

from __future__ import annotations

import json
import sqlite3
import threading

import pytest

from degora.discovery_store import (
    DiscoveryJobManager,
    DiscoveryQueueFullError,
    DiscoveryStateStore,
    DiscoveryStoreError,
    _sanitize_text,
    sanitize_discovery_payload,
)


def test_cancel_churn_never_exceeds_physical_or_tracked_admission_bound(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    manager = DiscoveryJobManager(store, max_workers=1, max_pending_jobs=1)
    executor_blocked = threading.Event()
    release_executor = threading.Event()

    def occupy_only_executor_worker() -> None:
        executor_blocked.set()
        assert release_executor.wait(timeout=10)

    blocker = manager._executor.submit(occupy_only_executor_worker)
    assert executor_blocked.wait(timeout=5)

    job = manager.submit(
        "search",
        {"index": 0},
        lambda _job_id, _payload, _progress: {"unexpected": True},
    )
    cancelled = manager.cancel(job["job_id"])
    assert cancelled is not None and cancelled["status"] == "cancelled"

    # Future.cancel() would report the queued Future done but leave its WorkItem
    # in ThreadPoolExecutor's queue. The retained admission slot makes hundreds
    # of cancel/re-submit attempts fail closed instead of accumulating tombstones.
    for index in range(1, 251):
        with pytest.raises(DiscoveryQueueFullError, match="queue is full"):
            manager.submit(
                "search",
                {"index": index},
                lambda _job_id, _payload, _progress: {"unexpected": True},
            )
        assert len(manager._futures) == 1
        assert len(manager._future_job_ids) == 1
        assert manager._executor._work_queue.qsize() <= 1

    assert [item["status"] for item in store.list_jobs(limit=300)] == ["cancelled"]
    release_executor.set()
    blocker.result(timeout=5)
    manager.shutdown(wait=True)
    assert len(manager._futures) == 0
    assert len(manager._future_job_ids) == 0


def test_cancelling_running_job_retains_admission_until_worker_exits(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    manager = DiscoveryJobManager(store, max_workers=1, max_pending_jobs=1)
    started = threading.Event()
    release = threading.Event()

    def worker(_job_id, _payload, _progress):
        started.set()
        assert release.wait(timeout=5)
        return {"done": True}

    first = manager.submit("search", {}, worker)
    assert started.wait(timeout=5)
    assert manager.cancel(first["job_id"])["status"] == "cancelled"

    with pytest.raises(DiscoveryQueueFullError, match="queue is full"):
        manager.submit("search", {}, worker)

    release.set()
    for _ in range(100):
        if not manager._future_job_ids:
            break
        threading.Event().wait(0.01)
    second = manager.submit("search", {}, lambda *_args: {"done": True})
    manager.shutdown(wait=True)

    assert store.get_job(first["job_id"])["status"] == "cancelled"
    assert store.get_job(second["job_id"])["status"] == "completed"


def test_url_redaction_fails_closed_and_preserves_surrounding_punctuation() -> None:
    text = (
        "broken http://[invalid?token=FIRSTSECRET)., then "
        "valid https://example.test/data?id=GSE1&access_token=SECONDSECRET), done"
    )

    cleaned = _sanitize_text(text)

    assert "FIRSTSECRET" not in cleaned and "SECONDSECRET" not in cleaned
    assert "[redacted: credential URL]).," in cleaned
    assert "id=GSE1&access_token=%5Bredacted%5D)," in cleaned
    assert cleaned.endswith(" done")


def test_semicolon_url_credentials_are_redacted_before_persistence(tmp_path) -> None:
    query_secret = "QUERY_SEMICOLON_SECRET"
    fragment_secret = "FRAGMENT_SEMICOLON_SECRET"
    payload = {
        "query": (
            "https://example.test/data?id=GSE1;"
            f"access_token={query_secret}&safe=1"
        ),
        "fragment": (
            "https://example.test/callback#state=ok;"
            f"session_id={fragment_secret}"
        ),
    }

    cleaned = sanitize_discovery_payload(payload)

    assert cleaned == {
        "query": (
            "https://example.test/data?id=GSE1;"
            "access_token=%5Bredacted%5D&safe=1"
        ),
        "fragment": (
            "https://example.test/callback#state=ok;"
            "session_id=%5Bredacted%5D"
        ),
    }
    assert sanitize_discovery_payload(cleaned) == cleaned

    store = DiscoveryStateStore(tmp_path / "state")
    store.save_search("1" * 16, payload)
    with sqlite3.connect(store.db_path) as connection:
        persisted = connection.execute(
            "SELECT payload_json FROM searches WHERE search_id=?", ("1" * 16,)
        ).fetchone()[0]

    assert query_secret not in persisted
    assert fragment_secret not in persisted
    assert "id=GSE1;access_token=%5Bredacted%5D&safe=1" in persisted
    assert "state=ok;session_id=%5Bredacted%5D" in persisted


@pytest.mark.parametrize(
    "source_url",
    [
        "https://x.test/cb?id=1%2526access_token=DOUBLEENCODED-DELIMITER",
        "https://x.test/cb?%2561ccess_token=DOUBLEENCODED-KEY",
        "https://x.test/cb#state=ok%253Bsession_id=DOUBLEENCODED-FRAGMENT",
    ],
)
def test_double_encoded_credential_structure_fails_closed_before_store(
    tmp_path, source_url: str
) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    store.save_artifact("bundle", "a" * 16, {"source_url": source_url})

    with sqlite3.connect(store.db_path) as connection:
        persisted = connection.execute(
            "SELECT payload_json FROM artifacts WHERE artifact_id=?", ("a" * 16,)
        ).fetchone()[0]

    assert "DOUBLEENCODED" not in persisted
    assert "[redacted: credential URL]" in persisted


@pytest.mark.parametrize(
    "source_url",
    [
        "https://x.test/cb?id=GSE1&token=FIRST&%2574oken=MIXED-SECOND",
        "https://x.test/cb?access_token=FIRST&id=GSE1%2526token=MIXED-SECOND",
        "https://x.test/cb?token=FIRST#state=ok%2526token=MIXED-SECOND",
    ],
)
def test_mixed_plain_and_encoded_credentials_fail_closed_after_partial_redaction(
    tmp_path, source_url: str
) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    store.save_artifact("bundle", "b" * 16, {"source_url": source_url})

    with sqlite3.connect(store.db_path) as connection:
        persisted = connection.execute(
            "SELECT payload_json FROM artifacts WHERE artifact_id=?", ("b" * 16,)
        ).fetchone()[0]

    assert "FIRST" not in persisted
    assert "MIXED-SECOND" not in persisted
    assert "[redacted: credential URL]" in persisted


@pytest.mark.parametrize("delimiter", ["&", ";"])
def test_sanitized_mapping_key_collision_is_rejected_without_leaking_secrets(
    tmp_path, delimiter: str
) -> None:
    payload = {
        f"https://example.test/data?id=GSE1{delimiter}token=FIRSTSECRET": {"record": 1},
        f"https://example.test/data?id=GSE1{delimiter}token=SECONDSECRET": {"record": 2},
    }

    with pytest.raises(DiscoveryStoreError, match="mapping keys collide") as exc_info:
        sanitize_discovery_payload(payload)
    assert "FIRSTSECRET" not in str(exc_info.value)
    assert "SECONDSECRET" not in str(exc_info.value)

    store = DiscoveryStateStore(tmp_path / "state")
    with pytest.raises(DiscoveryStoreError, match="mapping keys collide"):
        store.save_search("1" * 16, payload)
    assert store.get_search("1" * 16) is None

    manager = DiscoveryJobManager(store, max_workers=1, max_pending_jobs=1)
    with pytest.raises(DiscoveryStoreError, match="mapping keys collide"):
        manager.submit("search", payload, lambda *_args: {})
    accepted = manager.submit("search", {"safe": True}, lambda *_args: {"ok": True})
    manager.shutdown(wait=True)
    assert store.get_job(accepted["job_id"])["status"] == "completed"


def test_stringified_mapping_key_collision_is_rejected() -> None:
    with pytest.raises(DiscoveryStoreError, match="mapping keys collide"):
        sanitize_discovery_payload({1: "numeric", "1": "text"})


def test_common_pat_private_key_and_session_cookie_fields_are_redacted_before_store(
    tmp_path,
) -> None:
    secrets_by_key = {
        "github_pat": "GITHUB-PAT-SECRET",
        "gitPat": "GIT-PAT-SECRET",
        "private_key": "PRIVATE-KEY-SECRET",
        "secretKey": "SECRET-KEY-VALUE",
        "session_cookie": "SESSION-COOKIE-SECRET",
        "authCookie": "AUTH-COOKIE-SECRET",
    }
    store = DiscoveryStateStore(tmp_path / "state")
    store.save_search("2" * 16, {"provider_state": secrets_by_key, "id": "GSE123"})

    persisted = store.get_search("2" * 16)
    assert persisted["provider_state"] == {
        key: "[redacted]" for key in secrets_by_key
    }
    assert persisted["id"] == "GSE123"


def test_shutdown_cannot_return_while_commit_owner_is_still_publishing(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    manager = DiscoveryJobManager(store, max_workers=1)
    manager._commit_drain_timeout = 0.01
    commit_acquired = threading.Event()
    release_publication = threading.Event()
    published = tmp_path / "published.json"

    def worker(job_id, _payload, _progress):
        manager.commit(job_id)
        commit_acquired.set()
        assert release_publication.wait(timeout=5)
        published.write_text(json.dumps({"complete": True}), encoding="utf-8")
        return {"path": str(published)}

    job = manager.submit("analysis", {}, worker)
    assert commit_acquired.wait(timeout=5)
    shutdown_done = threading.Event()
    shutdown_thread = threading.Thread(
        target=lambda: (
            manager.shutdown(wait=False, cancel_futures=True, interrupt=True),
            shutdown_done.set(),
        )
    )
    shutdown_thread.start()

    assert not shutdown_done.wait(timeout=0.1)
    assert not published.exists()
    assert store.get_job(job["job_id"])["status"] == "running"

    release_publication.set()
    assert shutdown_done.wait(timeout=5)
    shutdown_thread.join(timeout=5)
    assert published.exists()
    final = store.get_job(job["job_id"])
    assert final["status"] == "completed"
    assert final["result"] == {"path": str(published)}
