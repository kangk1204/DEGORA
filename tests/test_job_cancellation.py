"""Cancelling a discovery job, at the store, the manager and the HTTP surface."""

from __future__ import annotations

import threading
import time

import pytest

from degora.discovery_store import (
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
