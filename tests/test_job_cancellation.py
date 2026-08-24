"""Cancelling a discovery job, at the store, the manager and the HTTP surface."""

from __future__ import annotations

import threading
import time

import pytest

from degora.discovery_store import (
    DiscoveryJobCancelled,
    DiscoveryJobManager,
    DiscoveryStateStore,
    DiscoveryStoreError,
)


@pytest.fixture()
def store(tmp_path) -> DiscoveryStateStore:
    return DiscoveryStateStore(tmp_path / "discovery")


def test_a_queued_job_can_be_cancelled_before_it_starts(store) -> None:
    job = store.create_job("search", {"query": "hypoxia"})

    cancelled = store.cancel_job(job["job_id"], "Job was cancelled by the reader.")

    assert cancelled is not None
    assert cancelled["status"] == "cancelled"
    assert store.get_job(job["job_id"])["completed_at"]


def test_cancelling_a_finished_job_reports_that_it_was_too_late(store) -> None:
    """The reader pressed the button a moment late; that is not an error."""

    job = store.create_job("search", {"query": "hypoxia"})
    store.update_job(job["job_id"], status="running")
    store.update_job(job["job_id"], status="completed", result={"hits": 3})

    assert store.cancel_job(job["job_id"], "Job was cancelled by the reader.") is None
    # And the result the reader can already see is still there.
    assert store.get_job(job["job_id"])["result"] == {"hits": 3}


def test_a_cancelled_job_is_final(store) -> None:
    job = store.create_job("search", {"query": "hypoxia"})
    store.cancel_job(job["job_id"], "Job was cancelled by the reader.")

    with pytest.raises(DiscoveryStoreError, match="invalid job transition"):
        store.update_job(job["job_id"], status="completed", result={"hits": 3})


def test_cancelling_an_unknown_job_is_an_error(store) -> None:
    with pytest.raises(DiscoveryStoreError, match="unknown job"):
        store.cancel_job("0" * 16, "Job was cancelled by the reader.")


def test_a_running_worker_stops_at_its_next_progress_report(store) -> None:
    """The worker is mid-download; the progress callback is where it notices."""

    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    may_continue = threading.Event()
    stages_run = []

    def worker(job_id, payload, progress):
        progress(0.1, "Stage one.")
        stages_run.append("one")
        started.set()
        may_continue.wait(timeout=5)
        progress(0.5, "Stage two.")  # must raise
        stages_run.append("two")
        return {"hits": 3}

    job = manager.submit("search", {}, worker)
    assert started.wait(timeout=5)

    cancelled = manager.cancel(job["job_id"])
    may_continue.set()
    manager.shutdown(wait=True)

    assert cancelled is not None and cancelled["status"] == "cancelled"
    assert stages_run == ["one"], "the worker kept going past the cancellation"
    final = store.get_job(job["job_id"])
    assert final["status"] == "cancelled"
    assert final["result"] is None, "a cancelled job must never record a result"


def test_commit_then_cancel_keeps_result_publishable(store) -> None:
    manager = DiscoveryJobManager(store, max_workers=1)
    committed = threading.Event()
    release = threading.Event()

    def worker(job_id, payload, progress):
        progress(0.5, "Ready to publish.")
        manager.commit(job_id)
        committed.set()
        release.wait(timeout=5)
        return {"hits": 3}

    job = manager.submit("search", {}, worker)
    assert committed.wait(timeout=5)

    assert manager.cancel(job["job_id"]) is None
    assert store.get_job(job["job_id"])["status"] == "running"

    release.set()
    manager.shutdown(wait=True)

    final = store.get_job(job["job_id"])
    assert final["status"] == "completed"
    assert final["result"] == {"hits": 3}
    assert manager._state == {}


def test_cancel_then_commit_stops_result_publish(store) -> None:
    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    may_commit = threading.Event()
    outcome: list[str] = []

    def worker(job_id, payload, progress):
        progress(0.2, "Working.")
        started.set()
        may_commit.wait(timeout=5)
        try:
            manager.commit(job_id)
        except BaseException as exc:  # noqa: BLE001 - assert the control-flow boundary
            outcome.append(type(exc).__name__)
            raise
        return {"hits": 3}

    job = manager.submit("search", {}, worker)
    assert started.wait(timeout=5)

    cancelled = manager.cancel(job["job_id"])
    may_commit.set()
    manager.shutdown(wait=True)

    assert cancelled is not None and cancelled["status"] == "cancelled"
    assert outcome == [DiscoveryJobCancelled.__name__]
    final = store.get_job(job["job_id"])
    assert final["status"] == "cancelled"
    assert final["result"] is None
    assert manager._state == {}


def test_duplicate_cancel_cannot_clear_cancelled_state_before_worker_commit(store) -> None:
    manager = DiscoveryJobManager(store, max_workers=1)
    worker_started = threading.Event()
    may_commit = threading.Event()
    first_cancel_in_store = threading.Event()
    release_first_cancel = threading.Event()
    cancel_calls = 0
    cancel_calls_lock = threading.Lock()
    commit_outcome: list[str] = []
    original_cancel_job = store.cancel_job

    def blocking_cancel_job(job_id, message):
        nonlocal cancel_calls
        with cancel_calls_lock:
            cancel_calls += 1
            call_number = cancel_calls
        if call_number == 1:
            first_cancel_in_store.set()
            assert release_first_cancel.wait(timeout=5)
        return original_cancel_job(job_id, message)

    store.cancel_job = blocking_cancel_job

    def worker(job_id, payload, progress):
        progress(0.5, "Ready to publish.")
        worker_started.set()
        may_commit.wait(timeout=5)
        try:
            manager.commit(job_id)
        except DiscoveryJobCancelled:
            commit_outcome.append("cancelled")
            raise
        return {"unsafe": True}

    job = manager.submit("search", {}, worker)
    assert worker_started.wait(timeout=5)

    outcomes: list[dict | None] = []
    first = threading.Thread(target=lambda: outcomes.append(manager.cancel(job["job_id"])))
    second = threading.Thread(target=lambda: outcomes.append(manager.cancel(job["job_id"])))
    first.start()
    assert first_cancel_in_store.wait(timeout=5)
    may_commit.set()
    second.start()
    time.sleep(0.05)
    release_first_cancel.set()
    first.join(timeout=5)
    second.join(timeout=5)
    manager.shutdown(wait=True)

    assert not first.is_alive() and not second.is_alive()
    # If the worker wins the wake-up race it removes its terminal in-memory state
    # before the duplicate request rechecks it; that second store call is then a
    # harmless terminal read. In either ordering it must not reopen publication.
    assert cancel_calls in {1, 2}
    assert sum(outcome is not None for outcome in outcomes) == 1
    assert commit_outcome == ["cancelled"]
    final = store.get_job(job["job_id"])
    assert final["status"] == "cancelled"
    assert final["result"] is None
    assert manager._state == {}


def test_a_worker_that_finishes_during_the_cancel_does_not_write_a_result(store) -> None:
    """The race that matters: cancel arrives between the last report and the return.

    Without the guard after the worker returns, the manager tried to write
    'completed' over a job the store had already moved to 'cancelled'. The store
    refuses that transition, so the job survived - but as an unexplained failure
    to a reader who had just been told it stopped cleanly.
    """

    manager = DiscoveryJobManager(store, max_workers=1)
    reported = threading.Event()
    may_return = threading.Event()

    def worker(job_id, payload, progress):
        progress(0.9, "Almost done.")
        reported.set()
        may_return.wait(timeout=5)
        return {"hits": 3}  # finishes with no further progress report

    job = manager.submit("search", {}, worker)
    assert reported.wait(timeout=5)

    manager.cancel(job["job_id"])
    may_return.set()
    manager.shutdown(wait=True)

    final = store.get_job(job["job_id"])
    assert final["status"] == "cancelled"
    assert final["result"] is None


def test_a_cancelled_worker_that_raises_is_not_reported_as_a_failure(store) -> None:
    """A cancelled download usually dies on the way out. That is the cancellation."""

    manager = DiscoveryJobManager(store, max_workers=1)
    reported = threading.Event()
    may_raise = threading.Event()

    def worker(job_id, payload, progress):
        progress(0.4, "Downloading.")
        reported.set()
        may_raise.wait(timeout=5)
        raise ConnectionResetError("connection closed mid-download")

    job = manager.submit("search", {}, worker)
    assert reported.wait(timeout=5)

    manager.cancel(job["job_id"])
    may_raise.set()
    manager.shutdown(wait=True)

    final = store.get_job(job["job_id"])
    assert final["status"] == "cancelled"
    assert final["error"] in (None, {}, "")


def test_cancelling_one_job_leaves_the_others_running(store) -> None:
    """Human and Mouse searches are separate jobs; stopping one is not stopping both."""

    manager = DiscoveryJobManager(store, max_workers=2)
    release = threading.Event()

    def worker(job_id, payload, progress):
        progress(0.2, "Working.")
        release.wait(timeout=5)
        progress(0.8, "Still working.")
        return {"species": payload["species"]}

    human = manager.submit("search", {"species": "human"}, worker)
    mouse = manager.submit("search", {"species": "mouse"}, worker)
    time.sleep(0.2)

    manager.cancel(human["job_id"])
    release.set()
    manager.shutdown(wait=True)

    assert store.get_job(human["job_id"])["status"] == "cancelled"
    assert store.get_job(mouse["job_id"])["status"] == "completed"
    assert store.get_job(mouse["job_id"])["result"] == {"species": "mouse"}


def test_terminal_or_unknown_cancel_leaves_no_manager_state(store) -> None:
    manager = DiscoveryJobManager(store, max_workers=1)
    job = store.create_job("search", {})
    store.update_job(job["job_id"], status="running")
    store.update_job(job["job_id"], status="completed", result={"ok": True})

    assert manager.cancel(job["job_id"]) is None
    assert manager._state == {}

    with pytest.raises(DiscoveryStoreError, match="unknown job"):
        manager.cancel("0" * 16)
    assert manager._state == {}


def test_cancel_store_baseexception_cleans_pending_state_and_wakes_commit(tmp_path) -> None:
    class RaisingCancelStore(DiscoveryStateStore):
        def __init__(self, root):
            super().__init__(root)
            self.entered_cancel = threading.Event()
            self.release_cancel = threading.Event()

        def cancel_job(self, job_id: str, message: str):
            self.entered_cancel.set()
            self.release_cancel.wait(timeout=5)
            raise KeyboardInterrupt("operator interrupted cancellation")

    store = RaisingCancelStore(tmp_path / "discovery")
    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    may_commit = threading.Event()
    committed = threading.Event()
    cancel_errors: list[str] = []

    def worker(job_id, payload, progress):
        progress(0.2, "Working.")
        started.set()
        may_commit.wait(timeout=5)
        manager.commit(job_id)
        committed.set()
        return {"hits": 3}

    job = manager.submit("search", {}, worker)
    assert started.wait(timeout=5)

    def cancel_job():
        try:
            manager.cancel(job["job_id"])
        except BaseException as exc:  # noqa: BLE001 - assert manager cleanup for BaseException.
            cancel_errors.append(type(exc).__name__)
        else:
            cancel_errors.append("no-error")

    cancel_thread = threading.Thread(target=cancel_job)
    cancel_thread.start()
    assert store.entered_cancel.wait(timeout=5)

    may_commit.set()
    time.sleep(0.05)
    assert not committed.is_set()

    store.release_cancel.set()
    cancel_thread.join(timeout=5)
    assert not cancel_thread.is_alive()
    assert committed.wait(timeout=5)
    manager.shutdown(wait=True)

    assert cancel_errors == ["KeyboardInterrupt"]
    final = store.get_job(job["job_id"])
    assert final["status"] == "completed"
    assert final["result"] == {"hits": 3}
    assert manager._state == {}


def test_repeated_cancelled_jobs_do_not_retain_manager_state(store) -> None:
    manager = DiscoveryJobManager(store, max_workers=2)
    releases = [threading.Event() for _ in range(5)]
    started = [threading.Event() for _ in range(5)]

    def worker(job_id, payload, progress):
        progress(0.1, "Started.")
        started[payload["index"]].set()
        releases[payload["index"]].wait(timeout=5)
        progress(0.9, "Should stop if cancelled.")
        return {"index": payload["index"]}

    jobs = [manager.submit("search", {"index": index}, worker) for index in range(5)]
    for event in started[:2]:
        assert event.wait(timeout=5)
    for job in jobs:
        manager.cancel(job["job_id"])
    for release in releases:
        release.set()
    manager.shutdown(wait=True)

    assert {store.get_job(job["job_id"])["status"] for job in jobs} == {"cancelled"}
    assert manager._state == {}
    assert manager._future_job_ids == {}
    assert manager._futures == set()


def test_cancel_pending_cleanup_when_worker_finishes_before_store_cancel(tmp_path) -> None:
    class SlowCancelStore(DiscoveryStateStore):
        def __init__(self, root):
            super().__init__(root)
            self.entered_cancel = threading.Event()
            self.release_cancel = threading.Event()

        def cancel_job(self, job_id: str, message: str):
            self.entered_cancel.set()
            self.release_cancel.wait(timeout=5)
            return super().cancel_job(job_id, message)

    store = SlowCancelStore(tmp_path / "discovery")
    manager = DiscoveryJobManager(store, max_workers=1)
    started = threading.Event()
    release_worker = threading.Event()
    cancel_result: list[dict | None] = []

    def worker(job_id, payload, progress):
        progress(0.2, "Working.")
        started.set()
        release_worker.wait(timeout=5)
        return {"hits": 3}

    job = manager.submit("search", {}, worker)
    assert started.wait(timeout=5)

    cancel_thread = threading.Thread(target=lambda: cancel_result.append(manager.cancel(job["job_id"])))
    cancel_thread.start()
    assert store.entered_cancel.wait(timeout=5)

    release_worker.set()
    manager.shutdown(wait=True)
    store.release_cancel.set()
    cancel_thread.join(timeout=5)

    assert cancel_result and cancel_result[0] is not None
    assert store.get_job(job["job_id"])["status"] == "cancelled"
    assert manager._state == {}
