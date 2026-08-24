"""Durable SQLite state for asynchronous discovery work."""

from __future__ import annotations

import concurrent.futures
import json
import re
import secrets
import sqlite3
import threading
from collections.abc import Callable, Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1"
DEFAULT_DB_NAME = "discovery_state.sqlite3"
# "cancelled" is kept apart from "interrupted" because they answer different
# questions about a job that did not finish. "interrupted" means the server
# stopped underneath it and the work is worth resuming; "cancelled" means a
# person decided it should not finish, and re-running it is their call to make.
TERMINAL_STATUSES = frozenset({"completed", "failed", "interrupted", "cancelled"})
VALID_STATUSES = frozenset({"queued", "running", *TERMINAL_STATUSES})
SEARCH_TERMINAL_STATUSES = frozenset({"complete", "failed", "interrupted", "cancelled"})
ALLOWED_TRANSITIONS = {
    "queued": frozenset({"running", "interrupted", "cancelled"}),
    "running": frozenset({"completed", "failed", "interrupted", "cancelled"}),
    "completed": frozenset(),
    "failed": frozenset(),
    "interrupted": frozenset(),
    "cancelled": frozenset(),
}
SENSITIVE_KEY_PARTS = (
    "api_key",
    "apikey",
    "authorization",
    "bearer",
    "client_secret",
    "credential",
    "ncbi_api_key",
    "password",
    "refresh_token",
    "secret",
    "token",
)
_MISSING = object()
_SENSITIVE_ASSIGNMENT_RE = re.compile(
    r"(?i)\b(api[_-]?key|apikey|ncbi[_-]?api[_-]?key|client[_-]?secret|access[_-]?token|"
    r"refresh[_-]?token|token|secret|password|credentials?)"
    r"(\s*[=:]\s*)([^,\s;]+)"
)
_BEARER_TOKEN_RE = re.compile(r"(?i)\b(Bearer\s+)([^,\s;]+)")
_FILE_URI_RE = re.compile(r'''(?i)file:/+[^\s,;)\]}"']+''')
_POSIX_ABSOLUTE_PATH_RE = re.compile(
    r'''(?<![\w:/.-])/(?!/)[^/\s,;)\]}"']+(?:/[^/\s,;)\]}"']+)*'''
)
_WINDOWS_ABSOLUTE_PATH_RE = re.compile(r"(?i)\b[A-Z]:[\\/](?:[^\\/\s,;)]+[\\/]?)+")
_UNC_PATH_RE = re.compile(r"\\\\[^\\\s,;)]+\\[^,\s;)]+")


class DiscoveryStoreError(ValueError):
    """Invalid discovery state operation."""


def _sanitize_text(value: str) -> str:
    text = _BEARER_TOKEN_RE.sub(r"\1[redacted]", value)
    text = _SENSITIVE_ASSIGNMENT_RE.sub(r"\1\2[redacted]", text)
    text = _FILE_URI_RE.sub("[redacted: local path]", text)
    text = _UNC_PATH_RE.sub("[redacted: local path]", text)
    text = _WINDOWS_ABSOLUTE_PATH_RE.sub("[redacted: local path]", text)
    text = _POSIX_ABSOLUTE_PATH_RE.sub("[redacted: local path]", text)
    return text


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def _new_id() -> str:
    return secrets.token_hex(8)


def _validate_id(value: str, *, field: str) -> str:
    text = str(value)
    if len(text) != 16 or any(char not in "0123456789abcdef" for char in text):
        raise DiscoveryStoreError(f"{field} must be a 16-character lowercase hex identifier")
    return text


def _sanitize_json(value: Any) -> Any:
    if isinstance(value, dict):
        cleaned = {}
        for key, item in value.items():
            text_key = str(key)
            lowered = text_key.lower()
            if any(part in lowered for part in SENSITIVE_KEY_PARTS):
                cleaned[text_key] = "[redacted]"
            else:
                cleaned[text_key] = _sanitize_json(item)
        return cleaned
    if isinstance(value, list):
        return [_sanitize_json(item) for item in value]
    if isinstance(value, tuple):
        return [_sanitize_json(item) for item in value]
    return value


def _json_dumps(value: Any) -> str:
    try:
        return json.dumps(
            _sanitize_json(value),
            allow_nan=False,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )
    except (TypeError, ValueError) as exc:
        raise DiscoveryStoreError("payload must be JSON-serializable") from exc


def _json_loads(value: str | None) -> Any:
    if value is None:
        return None
    return json.loads(value)


def _db_path(root_or_db: str | Path) -> Path:
    path = Path(root_or_db)
    if path.exists() and path.is_dir():
        return path / DEFAULT_DB_NAME
    if path.suffix:
        return path
    return path / DEFAULT_DB_NAME


class DiscoveryStateStore:
    """Small per-operation SQLite store for discovery jobs and outputs."""

    def __init__(self, root_or_db: str | Path):
        self.db_path = _db_path(root_or_db)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize()

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.db_path, timeout=30, isolation_level=None)
        connection.row_factory = sqlite3.Row
        return connection

    def _initialize(self) -> None:
        with self._connect() as connection:
            try:
                connection.execute("PRAGMA journal_mode=WAL")
            except sqlite3.DatabaseError:
                pass
            connection.execute("PRAGMA foreign_keys=ON")
            connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS jobs (
                    job_id TEXT PRIMARY KEY,
                    kind TEXT NOT NULL,
                    status TEXT NOT NULL,
                    progress REAL NOT NULL,
                    message TEXT NOT NULL,
                    payload_json TEXT NOT NULL,
                    result_json TEXT,
                    error_json TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    started_at TEXT,
                    completed_at TEXT
                );
                CREATE TABLE IF NOT EXISTS searches (
                    search_id TEXT PRIMARY KEY,
                    payload_json TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS artifacts (
                    kind TEXT NOT NULL,
                    artifact_id TEXT NOT NULL,
                    payload_json TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    PRIMARY KEY (kind, artifact_id)
                );
                CREATE INDEX IF NOT EXISTS idx_jobs_updated_at ON jobs(updated_at);
                """
            )
            try:
                connection.execute("BEGIN IMMEDIATE")
                version = connection.execute("SELECT value FROM metadata WHERE key='schema_version'").fetchone()
                if version is None:
                    connection.execute(
                        "INSERT INTO metadata(key, value) VALUES('schema_version', ?)",
                        (SCHEMA_VERSION,),
                    )
                elif version["value"] != SCHEMA_VERSION:
                    raise DiscoveryStoreError(
                        f"unsupported discovery store schema version: {version['value']}"
                    )
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise

    def _read_job_in_transaction(self, connection: sqlite3.Connection, job_id: str) -> sqlite3.Row:
        row = connection.execute(
            "SELECT * FROM jobs WHERE job_id=?",
            (_validate_id(job_id, field="job_id"),),
        ).fetchone()
        if row is None:
            raise KeyError(f"unknown discovery job: {job_id}")
        return row

    @staticmethod
    def _job_dict(row: sqlite3.Row) -> dict[str, Any]:
        return {
            "job_id": row["job_id"],
            "kind": row["kind"],
            "status": row["status"],
            "progress": row["progress"],
            "message": row["message"],
            "payload": _json_loads(row["payload_json"]),
            "result": _json_loads(row["result_json"]),
            "error": _json_loads(row["error_json"]),
            "created_at": row["created_at"],
            "updated_at": row["updated_at"],
            "started_at": row["started_at"],
            "completed_at": row["completed_at"],
        }

    def create_job(self, kind: str, payload: Any, job_id: str | None = None) -> dict[str, Any]:
        text_kind = str(kind).strip()
        if not text_kind:
            raise DiscoveryStoreError("kind is required")
        resolved_job_id = _validate_id(job_id, field="job_id") if job_id is not None else _new_id()
        now = _utc_now()
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """
                    INSERT INTO jobs(
                        job_id, kind, status, progress, message, payload_json,
                        created_at, updated_at
                    ) VALUES(?, ?, 'queued', 0.0, '', ?, ?, ?)
                    """,
                    (resolved_job_id, text_kind, _json_dumps(payload), now, now),
                )
                row = self._read_job_in_transaction(connection, resolved_job_id)
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        return self._job_dict(row)

    def update_job(
        self,
        job_id: str,
        *,
        status: str | None = None,
        progress: float | None = None,
        message: str | None = None,
        result: Any = _MISSING,
        error: Any = _MISSING,
    ) -> dict[str, Any]:
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                row = self._read_job_in_transaction(connection, job_id)
                current_status = row["status"]
                next_status = current_status if status is None else str(status)
                if next_status not in VALID_STATUSES:
                    raise DiscoveryStoreError(f"invalid job status: {next_status}")
                if next_status != current_status and next_status not in ALLOWED_TRANSITIONS[current_status]:
                    raise DiscoveryStoreError(f"invalid job transition: {current_status} -> {next_status}")

                current_progress = float(row["progress"])
                next_progress = current_progress if progress is None else float(progress)
                if next_progress < current_progress:
                    raise DiscoveryStoreError("job progress cannot decrease")
                if not 0.0 <= next_progress <= 1.0:
                    raise DiscoveryStoreError("job progress must be between 0.0 and 1.0")
                if next_status == "completed":
                    next_progress = 1.0

                now = _utc_now()
                started_at = row["started_at"]
                completed_at = row["completed_at"]
                if current_status == "queued" and next_status == "running":
                    started_at = now
                if next_status in TERMINAL_STATUSES and current_status != next_status:
                    completed_at = now

                connection.execute(
                    """
                    UPDATE jobs
                    SET status=?, progress=?, message=?, result_json=?, error_json=?,
                        updated_at=?, started_at=?, completed_at=?
                    WHERE job_id=?
                    """,
                    (
                        next_status,
                        next_progress,
                        row["message"] if message is None else str(message),
                        row["result_json"] if result is _MISSING else _json_dumps(result),
                        row["error_json"] if error is _MISSING else _json_dumps(error),
                        now,
                        started_at,
                        completed_at,
                        row["job_id"],
                    ),
                )
                updated = self._read_job_in_transaction(connection, row["job_id"])
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        return self._job_dict(updated)

    def get_job(self, job_id: str) -> dict[str, Any] | None:
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM jobs WHERE job_id=?",
                (_validate_id(job_id, field="job_id"),),
            ).fetchone()
        return None if row is None else self._job_dict(row)

    def list_jobs(self, limit: int = 100) -> list[dict[str, Any]]:
        if isinstance(limit, bool) or int(limit) < 1:
            raise DiscoveryStoreError("limit must be a positive integer")
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT * FROM jobs ORDER BY updated_at DESC, created_at DESC LIMIT ?",
                (int(limit),),
            ).fetchall()
        return [self._job_dict(row) for row in rows]

    def save_search(self, search_id: str, payload: Any) -> None:
        now = _utc_now()
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """
                    INSERT INTO searches(search_id, payload_json, created_at, updated_at)
                    VALUES(?, ?, ?, ?)
                    ON CONFLICT(search_id) DO UPDATE SET
                        payload_json=excluded.payload_json,
                        updated_at=excluded.updated_at
                    """,
                    (_validate_id(search_id, field="search_id"), _json_dumps(payload), now, now),
                )
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise

    def get_search(self, search_id: str) -> Any:
        with self._connect() as connection:
            row = connection.execute(
                "SELECT payload_json FROM searches WHERE search_id=?",
                (_validate_id(search_id, field="search_id"),),
            ).fetchone()
        return None if row is None else _json_loads(row["payload_json"])

    def save_artifact(self, kind: str, artifact_id: str, payload: Any) -> None:
        text_kind = str(kind).strip()
        if not text_kind:
            raise DiscoveryStoreError("kind is required")
        now = _utc_now()
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                connection.execute(
                    """
                    INSERT INTO artifacts(kind, artifact_id, payload_json, created_at, updated_at)
                    VALUES(?, ?, ?, ?, ?)
                    ON CONFLICT(kind, artifact_id) DO UPDATE SET
                        payload_json=excluded.payload_json,
                        updated_at=excluded.updated_at
                    """,
                    (text_kind, _validate_id(artifact_id, field="artifact_id"), _json_dumps(payload), now, now),
                )
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise

    def get_artifact(self, kind: str, artifact_id: str) -> Any:
        with self._connect() as connection:
            row = connection.execute(
                "SELECT payload_json FROM artifacts WHERE kind=? AND artifact_id=?",
                (str(kind).strip(), _validate_id(artifact_id, field="artifact_id")),
            ).fetchone()
        return None if row is None else _json_loads(row["payload_json"])

    def interrupt_active_jobs(
        self,
        message: str,
        *,
        exclude_job_ids: Iterable[str] = (),
    ) -> list[dict[str, Any]]:
        """Atomically interrupt active jobs except results already past commit."""

        excluded = {_validate_id(job_id, field="job_id") for job_id in exclude_job_ids}
        now = _utc_now()
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                rows = [
                    row
                    for row in connection.execute(
                        "SELECT * FROM jobs WHERE status IN ('queued', 'running')"
                    ).fetchall()
                    if row["job_id"] not in excluded
                ]
                for row in rows:
                    connection.execute(
                        """
                        UPDATE jobs
                        SET status='interrupted', message=?, updated_at=?, completed_at=COALESCE(completed_at, ?)
                        WHERE job_id=?
                        """,
                        (message, now, now, row["job_id"]),
                    )
                    if row["kind"] == "publication_search":
                        payload = _json_loads(row["payload_json"])
                        if isinstance(payload, dict):
                            search_id = payload.get("search_id")
                            if isinstance(search_id, str) and search_id:
                                try:
                                    search_id = _validate_id(search_id, field="search_id")
                                except DiscoveryStoreError:
                                    continue
                                search_row = connection.execute(
                                    "SELECT payload_json FROM searches WHERE search_id=?",
                                    (search_id,),
                                ).fetchone()
                                if search_row is not None:
                                    search = _json_loads(search_row["payload_json"])
                                    if (
                                        isinstance(search, dict)
                                        and search.get("status") not in SEARCH_TERMINAL_STATUSES
                                    ):
                                        interrupted = dict(search)
                                        interrupted.update(
                                            {"status": "interrupted", "error": message, "updated_at": now}
                                        )
                                        connection.execute(
                                            """
                                            UPDATE searches
                                            SET payload_json=?, updated_at=?
                                            WHERE search_id=?
                                            """,
                                            (_json_dumps(interrupted), now, search_id),
                                        )
                recovered = connection.execute(
                    "SELECT * FROM jobs WHERE job_id IN (%s) ORDER BY created_at"
                    % ",".join("?" for _ in rows),
                    tuple(row["job_id"] for row in rows),
                ).fetchall() if rows else []
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        return [self._job_dict(row) for row in recovered]

    def cancel_job(self, job_id: str, message: str) -> dict[str, Any] | None:
        """Move a queued or running job to cancelled, or report that it is too late.

        Returns the cancelled job, or None when the job had already reached a
        terminal state. The read and the write share one transaction because the
        worker may be completing in the same instant, and a cancel that silently
        overwrote a completed result would discard work the reader can already see.
        """

        _validate_id(job_id, field="job_id")
        now = _utc_now()
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                row = connection.execute(
                    "SELECT * FROM jobs WHERE job_id=?", (job_id,)
                ).fetchone()
                if row is None:
                    # The outer handler rolls back; doing it here too would leave
                    # it rolling back a transaction that is no longer open.
                    raise DiscoveryStoreError(f"unknown job: {job_id}")
                if row["status"] in TERMINAL_STATUSES:
                    connection.execute("ROLLBACK")
                    return None
                connection.execute(
                    """
                    UPDATE jobs
                    SET status='cancelled', message=?, updated_at=?,
                        completed_at=COALESCE(completed_at, ?)
                    WHERE job_id=?
                    """,
                    (message, now, now, job_id),
                )
                cancelled = self._read_job_in_transaction(connection, job_id)
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        return self._job_dict(cancelled)

    def recover_interrupted_jobs(self) -> list[dict[str, Any]]:
        return self.interrupt_active_jobs(
            "Job was interrupted before completion and marked interrupted during store recovery."
        )


class DiscoveryJobCancelled(BaseException):
    """Raised inside a worker when its job must stop: a shutdown, or a person.

    It derives from BaseException, not Exception, for the same reason
    KeyboardInterrupt does: the discovery workers wrap progress reporting in
    ``except Exception`` so that a reporting failure can never break a search,
    and a cancellation that those guards could swallow would not cancel anything.
    """


class DiscoveryJobManager:
    """Submit discovery jobs to a thread pool while persisting state."""

    def __init__(self, store: DiscoveryStateStore, max_workers: int = 4):
        self.store = store
        self._executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)
        self._futures: set[concurrent.futures.Future[Any]] = set()
        self._lock = threading.Lock()
        self._condition = threading.Condition(self._lock)
        self._closing = threading.Event()
        self._state: dict[str, str] = {}
        self._future_job_ids: dict[concurrent.futures.Future[Any], str] = {}
        self._commit_drain_timeout = 2.0

    def _job_stop_reason(self, job_id: str) -> str:
        with self._condition:
            state = self._state.get(job_id)
            if state == "committing":
                return ""
            if self._closing.is_set():
                return "discovery job was cancelled because the local DEGORA server is stopping"
            if state in {"cancel_pending", "cancelled"}:
                return "discovery job was cancelled by the reader"
            return ""

    def _finish_job_state(self, job_id: str) -> None:
        with self._condition:
            self._state.pop(job_id, None)
            self._condition.notify_all()

    def commit(self, job_id: str) -> None:
        """Fence off the point after which a job must publish its terminal state.

        Workers call this immediately before making a completed result visible.
        The decision is linearized with reader cancellation while keeping SQLite
        writes outside the manager lock.
        """

        _validate_id(job_id, field="job_id")
        with self._condition:
            while True:
                state = self._state.get(job_id)
                if state == "committing":
                    return
                if self._closing.is_set() or state == "cancelled":
                    raise DiscoveryJobCancelled("discovery job was cancelled before its result was published")
                if state == "cancel_pending":
                    self._condition.wait()
                    continue
                self._state[job_id] = "committing"
                self._condition.notify_all()
                return

    def submit(
        self,
        kind: str,
        payload: Any,
        worker: Callable[[str, Any, Callable[[float, str | None], dict[str, Any]]], Any],
        job_id: str | None = None,
    ) -> dict[str, Any]:
        if self._closing.is_set():
            raise DiscoveryStoreError("discovery manager is shutting down")
        job = self.store.create_job(kind, payload, job_id=job_id)
        if self._closing.is_set():
            self.store.update_job(
                job["job_id"],
                status="interrupted",
                message="Job was interrupted because the local DEGORA server was stopped.",
            )
            raise DiscoveryStoreError("discovery manager is shutting down")

        def progress_callback(progress: float, message: str | None = None) -> dict[str, Any]:
            # Cancellation point. cancel_futures only drops work that has not
            # started; a worker already downloading or writing needs somewhere to
            # notice, and every stage already reports here.
            stop_reason = self._job_stop_reason(job["job_id"])
            if stop_reason:
                raise DiscoveryJobCancelled(stop_reason)
            return self.store.update_job(job["job_id"], progress=progress, message=message)

        def complete_job(result: Any) -> None:
            try:
                self.store.update_job(job["job_id"], status="completed", result=result, message="Job completed.")
            except DiscoveryStoreError:
                terminal = self.store.get_job(job["job_id"])
                if terminal and terminal["status"] in TERMINAL_STATUSES:
                    return
                raise

        def fail_job(exc: BaseException) -> None:
            try:
                self.store.update_job(
                    job["job_id"],
                    status="failed",
                    error={"type": type(exc).__name__, "message": _sanitize_text(str(exc))},
                    message="Job failed.",
                )
            except DiscoveryStoreError:
                terminal = self.store.get_job(job["job_id"])
                if terminal and terminal["status"] in TERMINAL_STATUSES:
                    return
                raise

        def should_stop_without_terminal_write() -> bool:
            return bool(self._job_stop_reason(job["job_id"]))

        def run() -> None:
            try:
                if self._closing.is_set():
                    return
                self.store.update_job(job["job_id"], status="running", message="Job started.")
                result = worker(job["job_id"], payload, progress_callback)
                if should_stop_without_terminal_write():
                    return
                complete_job(result)
            except DiscoveryJobCancelled:
                # interrupt_active_jobs has already written the terminal state, and
                # anything this worker would have written next is exactly what the
                # cancellation exists to prevent.
                return
            except BaseException as exc:  # noqa: BLE001 - persist arbitrary worker failures.
                if should_stop_without_terminal_write():
                    return
                fail_job(exc)
            finally:
                self._finish_job_state(job["job_id"])

        try:
            # Register the future while holding the same condition used by
            # commit/shutdown. A very fast worker can otherwise enter its commit
            # section before shutdown can discover which future must be drained.
            with self._condition:
                future = self._executor.submit(run)
                self._futures.add(future)
                self._future_job_ids[future] = job["job_id"]
                self._condition.notify_all()
        except RuntimeError:
            self.store.update_job(
                job["job_id"],
                status="interrupted",
                message="Job was interrupted because the local DEGORA server was stopped.",
            )
            raise
        future.add_done_callback(self._discard_future)
        return job

    def cancel(self, job_id: str) -> dict[str, Any] | None:
        """Stop a job at the reader's request.

        Returns the cancelled job, or None if it had already finished - the
        caller needs that difference, because "too late" and "stopped" are
        different things to tell someone who just pressed a button.

        Files already downloaded are left where they are. They are valid inputs
        that cost a download to fetch, and a later run finds them cached. What
        cancelling guarantees is narrower and more important: a cancelled job
        never records a result, so partial work cannot be read as a finished
        search.
        """

        _validate_id(job_id, field="job_id")
        with self._condition:
            while True:
                state = self._state.get(job_id)
                if state == "committing":
                    return None
                if state == "cancelled":
                    return None
                if state == "cancel_pending":
                    self._condition.wait()
                    continue
                self._state[job_id] = "cancel_pending"
                self._condition.notify_all()
                break
        try:
            cancelled = self.store.cancel_job(job_id, "Job was cancelled by the reader.")
        except BaseException:
            with self._condition:
                if self._state.get(job_id) == "cancel_pending":
                    self._state.pop(job_id, None)
                    self._condition.notify_all()
            raise
        with self._condition:
            live_future = any(
                mapped_job_id == job_id and not future.done()
                for future, mapped_job_id in self._future_job_ids.items()
            )
            if cancelled is None or not live_future:
                if self._state.get(job_id) == "cancel_pending":
                    self._state.pop(job_id, None)
            else:
                self._state[job_id] = "cancelled"
            self._condition.notify_all()
        return cancelled

    def _discard_future(self, future: concurrent.futures.Future[Any]) -> None:
        with self._condition:
            self._futures.discard(future)
            job_id = self._future_job_ids.pop(future, None)
            if job_id is not None:
                self._state.pop(job_id, None)
            self._condition.notify_all()

    def _await_committing(self) -> set[str]:
        with self._condition:
            committing = {
                future
                for future, job_id in self._future_job_ids.items()
                if self._state.get(job_id) == "committing" and not future.done()
            }
        if committing:
            concurrent.futures.wait(committing, timeout=self._commit_drain_timeout)
        with self._condition:
            return {
                job_id
                for future, job_id in self._future_job_ids.items()
                if self._state.get(job_id) == "committing" and not future.done()
            }

    def shutdown(
        self,
        *,
        wait: bool = True,
        cancel_futures: bool = False,
        interrupt: bool = False,
    ) -> None:
        """Stop accepting work; optionally return promptly and persist interruption."""

        if interrupt:
            self._closing.set()
            still_committing = self._await_committing()
            self.store.interrupt_active_jobs(
                "Job was interrupted because the local DEGORA server was stopped.",
                exclude_job_ids=still_committing,
            )
        self._executor.shutdown(wait=wait, cancel_futures=cancel_futures)
