"""Regression tests for the v0.4.34 security, concurrency and export audit."""

from __future__ import annotations

import json
import multiprocessing
import os
import sqlite3
import threading
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import pytest

from degora import api, discovery_export, discovery_run, provenance
from degora.discovery_export import (
    SEARCH_CSV_NAME,
    SEARCH_JSON_NAME,
    SEARCH_MANIFEST_NAME,
    SEARCH_XLSX_NAME,
    export_publication_search,
    verify_publication_search_export,
)
from degora.discovery_store import (
    DiscoveryJobManager,
    DiscoveryQueueFullError,
    DiscoveryStateStore,
    _sanitize_json,
)
from degora.excel_export import export_run_workbook, verify_run_workbook_export
from degora.provenance import source_sidecar_payloads


def _minimal_score_db(path: Path) -> None:
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE genes (gene_symbol TEXT, degora_rank INTEGER, degora_score REAL)")
        connection.execute("CREATE TABLE gene_evidence (gene_symbol TEXT, source_unit_id TEXT, study_id TEXT)")
        connection.execute("CREATE TABLE studies (study_id TEXT, source_unit_id TEXT)")
        connection.execute("CREATE TABLE meta (key TEXT, value TEXT)")
        connection.execute("INSERT INTO genes VALUES ('ISG15', 1, 10.0)")
        connection.execute("INSERT INTO gene_evidence VALUES ('ISG15', 'U1', 'S1')")
        connection.execute("INSERT INTO studies VALUES ('S1', 'U1')")
        connection.execute("INSERT INTO meta VALUES ('degora_version', 'test')")


def _request_json(
    url: str,
    *,
    method: str = "GET",
    payload: dict | None = None,
    headers: dict[str, str] | None = None,
) -> tuple[int, dict]:
    body = None if payload is None else json.dumps(payload).encode("utf-8")
    request_headers = dict(headers or {})
    if body is not None:
        request_headers.setdefault("Content-Type", "application/json")
        request_headers.setdefault("X-DEGORA-Action", "1")
    request = urllib.request.Request(url, data=body, headers=request_headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode("utf-8"))


def _snapshot(query: str = "hypoxia") -> dict:
    return {
        "query": query,
        "species": {"key": "human"},
        "records": [
            {
                "canonical_id": "PMID:123",
                "paper_title": query,
                "candidates": [
                    {
                        "candidate_id": "c1",
                        "source_url": (
                            "https://storage.googleapis.test/data.csv?X-Goog-Credential=GOOGSECRET"
                            "&X-Goog-Signature=SIGNED&id=GSE123#session_token=SESSIONSECRET"
                        ),
                    }
                ],
            }
        ],
    }


def _minimal_result_dir(path: Path) -> Path:
    path.mkdir()
    with sqlite3.connect(path / "degora_scores.db") as connection:
        pd.DataFrame(
            {
                "gene_symbol": ["ISG15"],
                "quality_weighted_degora_rank": [1],
                "quality_weighted_degora_score": [1.0],
            }
        ).to_sql("genes", connection, index=False)
        pd.DataFrame({"gene_symbol": ["ISG15"], "study_id": ["S1"]}).to_sql(
            "gene_evidence", connection, index=False
        )
        pd.DataFrame({"study_id": ["S1"], "source_unit_id": ["U1"]}).to_sql(
            "studies", connection, index=False
        )
        pd.DataFrame({"key": ["corpus"], "value": ["test"]}).to_sql(
            "meta", connection, index=False
        )
    (path / "degora_score_metadata.json").write_text("{}\n", encoding="utf-8")
    return path


@pytest.mark.skipif(os.name == "nt", reason="the inode handoff regression is POSIX-specific")
@pytest.mark.filterwarnings("ignore:This process .* multi-threaded:DeprecationWarning")
def test_lock_waiter_and_third_process_share_one_stable_inode(tmp_path) -> None:
    context = multiprocessing.get_context("fork")
    output = tmp_path / "results"
    first_held = context.Event()
    release_first = context.Event()
    first_unlocked = context.Event()
    finish_first = context.Event()
    second_held = context.Event()
    release_second = context.Event()
    third_refused = context.Event()
    errors = context.Queue()

    def first_process() -> None:
        try:
            import fcntl

            real_flock = fcntl.flock

            def pause_after_unlock(fd, operation):
                result = real_flock(fd, operation)
                if operation == fcntl.LOCK_UN:
                    first_unlocked.set()
                    assert finish_first.wait(timeout=5)
                return result

            fcntl.flock = pause_after_unlock
            with provenance.output_directory_lock(output):
                first_held.set()
                assert release_first.wait(timeout=5)
        except BaseException as exc:  # noqa: BLE001 - relayed to parent
            errors.put(repr(exc))

    def second_process() -> None:
        try:
            assert first_unlocked.wait(timeout=5)
            with provenance.output_directory_lock(output):
                second_held.set()
                assert release_second.wait(timeout=5)
        except BaseException as exc:  # noqa: BLE001 - relayed to parent
            errors.put(repr(exc))

    def third_process() -> None:
        try:
            with provenance.output_directory_lock(output):
                errors.put("third process entered while the second held the lock")
        except provenance.OutputDirectoryBusyError:
            third_refused.set()
        except BaseException as exc:  # noqa: BLE001 - relayed to parent
            errors.put(repr(exc))

    first = context.Process(target=first_process)
    second = context.Process(target=second_process)
    first.start()
    second.start()
    assert first_held.wait(timeout=5)
    inode = (output / ".degora-run.lock").stat().st_ino
    release_first.set()
    assert first_unlocked.wait(timeout=5)
    assert second_held.wait(timeout=5)
    finish_first.set()
    first.join(timeout=5)
    assert not first.is_alive()
    assert (output / ".degora-run.lock").stat().st_ino == inode

    third = context.Process(target=third_process)
    third.start()
    third.join(timeout=5)
    assert not third.is_alive() and third_refused.is_set()
    release_second.set()
    second.join(timeout=5)
    assert not second.is_alive()
    assert errors.empty(), [errors.get() for _ in range(errors.qsize())]


def test_cancel_before_analysis_publish_rolls_back_the_complete_staged_run(tmp_path, monkeypatch) -> None:
    staged = threading.Event()
    release = threading.Event()
    output = tmp_path / "run"

    def fake_execute(_prepared, _selections, output_dir, **_kwargs):
        target = Path(output_dir)
        (target / ".degora-discovery-run.json").write_text('{"status":"complete"}\n', encoding="utf-8")
        staged.set()
        assert release.wait(timeout=5)
        return {"status": "complete", "output_dir": str(target)}

    monkeypatch.setattr(discovery_run, "_execute_discovery_analysis", fake_execute)
    store = DiscoveryStateStore(tmp_path / "state")
    manager = DiscoveryJobManager(store, max_workers=1)

    def worker(job_id, _payload, _progress):
        return discovery_run.run_discovery_analysis(
            {},
            [],
            output,
            species="human",
            before_publish=lambda: manager.commit(job_id),
        )

    job = manager.submit("discovery_analyze", {}, worker)
    assert staged.wait(timeout=5)
    cancelled = manager.cancel(job["job_id"])
    release.set()
    manager.shutdown(wait=True)

    assert cancelled is not None and cancelled["status"] == "cancelled"
    assert store.get_job(job["job_id"])["result"] is None
    assert not output.exists(), "a cancelled run published an orphan output directory"


def test_credential_urls_are_redacted_recursively_but_safe_ids_are_unchanged(tmp_path) -> None:
    safe = "https://example.test/data?id=GSE123&accession=PMID%3A456&page_token=NEXT_PAGE"
    payload = {
        "safe": safe,
        "safe_fields": {"page_token": "NEXT_PAGE", "token_count": 12, "source_unit_id": "PMID:123"},
        "nested": [
            {
                "aws": "https://s3.test/x?X-Amz-Credential=AWSSECRET&X-Amz-Signature=SIG&id=GSE1",
                "google": "https://g.test/x?X-Goog-Credential=GOOGSECRET&X-Goog-Signature=SIG",
                "oauth": "https://oauth.test/cb#access_token=JWTSECRET&state=SAFE_STATE",
                "session_url": "see https://app.test/x?session_id=SESSIONSECRET&record_id=PMID123 now",
            }
        ],
    }

    cleaned = _sanitize_json(payload)
    encoded = json.dumps(cleaned, sort_keys=True)
    assert all(secret not in encoded for secret in ("AWSSECRET", "GOOGSECRET", "JWTSECRET", "SESSIONSECRET"))
    assert cleaned["safe"] == safe
    assert cleaned["safe_fields"] == {
        "page_token": "NEXT_PAGE",
        "token_count": 12,
        "source_unit_id": "PMID:123",
    }
    assert "id=GSE1" in cleaned["nested"][0]["aws"]
    assert "state=SAFE_STATE" in cleaned["nested"][0]["oauth"]
    assert "record_id=PMID123" in cleaned["nested"][0]["session_url"]

    store = DiscoveryStateStore(tmp_path / "state.sqlite3")
    store.save_search("a" * 16, payload)
    persisted = json.dumps(store.get_search("a" * 16), sort_keys=True)
    assert all(secret not in persisted for secret in ("AWSSECRET", "GOOGSECRET", "JWTSECRET", "SESSIONSECRET"))


def test_search_export_redacts_credentials_and_manifest_detects_tampering(tmp_path) -> None:
    result = export_publication_search(_snapshot(), tmp_path)
    verify_publication_search_export(tmp_path)
    combined = Path(result["search_json"]).read_text(encoding="utf-8") + Path(result["search_csv"]).read_text(
        encoding="utf-8"
    )
    assert "GOOGSECRET" not in combined and "SESSIONSECRET" not in combined
    assert "GSE123" in combined

    Path(result["search_csv"]).write_text("tampered\n", encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        verify_publication_search_export(tmp_path)


def test_search_export_failure_restores_the_previous_complete_generation(tmp_path, monkeypatch) -> None:
    export_publication_search(_snapshot("old"), tmp_path)
    targets = [
        tmp_path / SEARCH_JSON_NAME,
        tmp_path / SEARCH_CSV_NAME,
        tmp_path / SEARCH_XLSX_NAME,
        tmp_path / SEARCH_MANIFEST_NAME,
    ]
    before = {path.name: path.read_bytes() for path in targets}
    real_replace = provenance.os.replace
    failed = False

    def fail_xlsx_publish(source, destination):
        nonlocal failed
        if (
            not failed
            and Path(destination) == tmp_path / SEARCH_XLSX_NAME
            and Path(source).suffix == ".pending"
        ):
            failed = True
            raise OSError("synthetic XLSX publication failure")
        return real_replace(source, destination)

    monkeypatch.setattr(provenance.os, "replace", fail_xlsx_publish)
    with pytest.raises(OSError, match="synthetic XLSX"):
        export_publication_search(_snapshot("new"), tmp_path, force=True)

    assert failed
    assert {path.name: path.read_bytes() for path in targets} == before
    verify_publication_search_export(tmp_path)


def test_workbook_manifest_detects_a_mixed_generation(tmp_path) -> None:
    result = export_run_workbook(_minimal_result_dir(tmp_path / "results"), command="pytest v0434")
    verify_run_workbook_export(result["output"])
    Path(result["validation"]).write_text("tampered=true\n", encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        verify_run_workbook_export(result["output"])


def test_workbook_set_publication_failure_restores_manifest_and_sidecars(tmp_path, monkeypatch) -> None:
    result_dir = _minimal_result_dir(tmp_path / "results")
    result = export_run_workbook(result_dir, command="pytest old generation")
    primary = [Path(result[key]) for key in ("output", "manifest", "validation")]
    targets = [
        candidate
        for artifact in primary
        for candidate in (
            artifact,
            provenance.artifact_source_path(artifact),
            provenance.artifact_provenance_path(artifact),
        )
    ]
    before = {path.name: path.read_bytes() for path in targets}
    (result_dir / "degora_score_metadata.json").write_text(
        json.dumps({"generation": "new"}) + "\n",
        encoding="utf-8",
    )
    real_replace = provenance.os.replace
    failed = False

    def fail_validation_publish(source, destination):
        nonlocal failed
        if (
            not failed
            and Path(destination) == Path(result["validation"])
            and Path(source).suffix == ".pending"
        ):
            failed = True
            raise OSError("synthetic validation publication failure")
        return real_replace(source, destination)

    monkeypatch.setattr(provenance.os, "replace", fail_validation_publish)
    with pytest.raises(OSError, match="synthetic validation"):
        export_run_workbook(result_dir, command="pytest new generation")

    assert failed
    assert {path.name: path.read_bytes() for path in targets} == before
    verify_run_workbook_export(result["output"])


def test_provenance_json_refuses_nonfinite_metadata(tmp_path) -> None:
    artifact = tmp_path / "artifact.csv"
    artifact.write_text("x\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Out of range float values"):
        source_sidecar_payloads(artifact, "pytest", metadata={"metric": float("nan")})


def test_durable_queue_admission_is_bounded_atomically(tmp_path) -> None:
    store = DiscoveryStateStore(tmp_path / "state")
    manager = DiscoveryJobManager(store, max_workers=1, max_pending_jobs=2)
    started = threading.Event()
    release = threading.Event()

    def worker(_job_id, _payload, _progress):
        started.set()
        assert release.wait(timeout=5)
        return {}

    manager.submit("search", {"index": 1}, worker)
    assert started.wait(timeout=5)
    manager.submit("search", {"index": 2}, worker)
    with pytest.raises(DiscoveryQueueFullError, match="queue is full"):
        manager.submit("search", {"index": 3}, worker)
    with sqlite3.connect(store.db_path) as connection:
        assert connection.execute(
            "SELECT COUNT(*) FROM jobs WHERE status IN ('queued', 'running')"
        ).fetchone()[0] == 2
    release.set()
    manager.shutdown(wait=True)


def test_http_queue_overload_is_429(tmp_path, monkeypatch) -> None:
    db = tmp_path / "scores.db"
    _minimal_score_db(db)
    started = threading.Event()
    release = threading.Event()

    def blocked_search(*_args, **_kwargs):
        started.set()
        assert release.wait(timeout=5)
        return {"records": [], "total": 0}

    monkeypatch.setattr(api, "_call_search_publications", blocked_search)
    server = api.create_server(db, port=0, quiet=True, max_pending_jobs=1)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address[:2]
    url = f"http://{host}:{port}/api/discovery/searches"
    try:
        first_status, _ = _request_json(
            url,
            method="POST",
            payload={"query": "hypoxia", "species": "human", "limit": 5},
        )
        assert first_status == 202 and started.wait(timeout=5)
        second_status, second = _request_json(
            url,
            method="POST",
            payload={"query": "fibrosis", "species": "human", "limit": 5},
        )
        assert second_status == 429
        assert "queue is full" in second["error"]
    finally:
        release.set()
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_serve_defaults_to_fragment_token_and_protects_api_end_to_end(tmp_path, monkeypatch) -> None:
    db = tmp_path / "scores.db"
    _minimal_score_db(db)
    original_create_server = api.create_server
    ready = threading.Event()
    captured: dict[str, object] = {}

    def capturing_create_server(*args, **kwargs):
        server = original_create_server(*args, **kwargs)
        captured.update({"server": server, "token": kwargs.get("access_token")})
        ready.set()
        return server

    monkeypatch.setattr(api, "create_server", capturing_create_server)
    thread = threading.Thread(target=lambda: api.serve(db, port=0, quiet=True), daemon=True)
    thread.start()
    assert ready.wait(timeout=5)
    server = captured["server"]
    token = str(captured["token"])
    assert token and token != "None"
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    fragment_url = f"{base}#token={urllib.parse.quote(token, safe='')}"
    supplied = urllib.parse.parse_qs(urllib.parse.urlsplit(fragment_url).fragment)["token"][0]
    try:
        assert _request_json(f"{base}/api/health")[0] == 401
        status, health = _request_json(
            f"{base}/api/health",
            headers={"X-DEGORA-Token": supplied},
        )
        assert status == 200 and health["status"] == "ok"
    finally:
        server.shutdown()
        thread.join(timeout=5)


def test_open_browser_receives_only_a_private_bootstrap_file_url(tmp_path, monkeypatch) -> None:
    db = tmp_path / "scores.db"
    _minimal_score_db(db)
    captured: dict[str, object] = {}

    class InterruptingServer:
        db_path = db
        server_address = ("127.0.0.1", 8765)

        def serve_forever(self):
            raise KeyboardInterrupt

        def server_close(self):
            return None

    def fake_create_server(*_args, **kwargs):
        captured["token"] = kwargs["access_token"]
        return InterruptingServer()

    def capture_bootstrap(url: str) -> bool:
        captured["url"] = url
        parsed = urllib.parse.urlsplit(url)
        bootstrap = Path(urllib.request.url2pathname(parsed.path))
        captured["bootstrap"] = bootstrap
        captured["document"] = bootstrap.read_text(encoding="utf-8")
        captured["file_mode"] = bootstrap.stat().st_mode & 0o777
        captured["directory_mode"] = bootstrap.parent.stat().st_mode & 0o777
        return True

    monkeypatch.setattr(api, "create_server", fake_create_server)
    monkeypatch.setattr(api.webbrowser, "open", capture_bootstrap)

    api.serve(db, quiet=True, open_browser=True)

    opened_url = str(captured["url"])
    parsed = urllib.parse.urlsplit(opened_url)
    token = str(captured["token"])
    assert parsed.scheme == "file" and not parsed.query and not parsed.fragment
    assert token not in opened_url and urllib.parse.quote(token, safe="") not in opened_url
    assert captured["file_mode"] == 0o600
    assert captured["directory_mode"] == 0o700
    assert f"#token={urllib.parse.quote(token, safe='')}" in str(captured["document"])
    assert not Path(captured["bootstrap"]).exists()


def test_search_manifest_verifier_rejects_every_incomplete_commit_shape(tmp_path) -> None:
    result = export_publication_search(_snapshot(), tmp_path)
    manifest_path = Path(result["manifest"])
    original_text = manifest_path.read_text(encoding="utf-8")
    original = json.loads(original_text)

    manifest_path.unlink()
    with pytest.raises(ValueError, match="missing or invalid"):
        verify_publication_search_export(tmp_path)

    for invalid_text in ("{not-json", '{"value": NaN}\n'):
        manifest_path.write_text(invalid_text, encoding="utf-8")
        with pytest.raises(ValueError, match="missing or invalid"):
            verify_publication_search_export(tmp_path)

    unsupported_values = (
        [],
        {**original, "artifact_type": "other"},
        {**original, "format_version": 999},
    )
    for value in unsupported_values:
        manifest_path.write_text(json.dumps(value), encoding="utf-8")
        with pytest.raises(ValueError, match="unsupported"):
            verify_publication_search_export(tmp_path)

    for files in ([], {SEARCH_JSON_NAME: original["files"][SEARCH_JSON_NAME]}):
        manifest_path.write_text(json.dumps({**original, "files": files}), encoding="utf-8")
        with pytest.raises(ValueError, match="complete artifact set"):
            verify_publication_search_export(tmp_path)

    invalid_entry = json.loads(original_text)
    invalid_entry["files"][SEARCH_JSON_NAME] = None
    manifest_path.write_text(json.dumps(invalid_entry), encoding="utf-8")
    with pytest.raises(ValueError, match="artifact is missing"):
        verify_publication_search_export(tmp_path)

    missing_artifact = Path(result["search_json"])
    missing_bytes = missing_artifact.read_bytes()
    manifest_path.write_text(original_text, encoding="utf-8")
    missing_artifact.unlink()
    with pytest.raises(ValueError, match="artifact is missing"):
        verify_publication_search_export(tmp_path)
    missing_artifact.write_bytes(missing_bytes)

    wrong_size = json.loads(original_text)
    wrong_size["files"][SEARCH_JSON_NAME]["size_bytes"] += 1
    manifest_path.write_text(json.dumps(wrong_size), encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        verify_publication_search_export(tmp_path)

    wrong_generation = json.loads(original_text)
    wrong_generation["generation_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(wrong_generation), encoding="utf-8")
    with pytest.raises(ValueError, match="generation digest"):
        verify_publication_search_export(tmp_path)


def test_workbook_manifest_verifier_rejects_every_incomplete_commit_shape(tmp_path) -> None:
    result = export_run_workbook(_minimal_result_dir(tmp_path / "results"), command="pytest verifier")
    workbook = Path(result["output"])
    manifest_path = Path(result["manifest"])
    original_text = manifest_path.read_text(encoding="utf-8")
    original = json.loads(original_text)

    manifest_path.unlink()
    with pytest.raises(ValueError, match="missing or invalid"):
        verify_run_workbook_export(workbook)

    for invalid_text in ("{not-json", '{"value": Infinity}\n'):
        manifest_path.write_text(invalid_text, encoding="utf-8")
        with pytest.raises(ValueError, match="missing or invalid"):
            verify_run_workbook_export(workbook)

    unsupported_values = (
        [],
        {**original, "artifact_type": "other"},
        {**original, "format_version": 999},
    )
    for value in unsupported_values:
        manifest_path.write_text(json.dumps(value), encoding="utf-8")
        with pytest.raises(ValueError, match="unsupported"):
            verify_run_workbook_export(workbook)

    for integrity in ([], {workbook.name: original["output_integrity"][workbook.name]}):
        manifest_path.write_text(json.dumps({**original, "output_integrity": integrity}), encoding="utf-8")
        with pytest.raises(ValueError, match="complete content set"):
            verify_run_workbook_export(workbook)

    manifest_path.write_text(original_text, encoding="utf-8")
    source_path = provenance.artifact_source_path(workbook)
    source_bytes = source_path.read_bytes()
    source_path.unlink()
    with pytest.raises(ValueError, match="artifact is missing"):
        verify_run_workbook_export(workbook)
    source_path.write_bytes(source_bytes)

    wrong_integrity = json.loads(original_text)
    wrong_integrity["output_integrity"][workbook.name]["size_bytes"] += 1
    manifest_path.write_text(json.dumps(wrong_integrity), encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        verify_run_workbook_export(workbook)

    wrong_generation = json.loads(original_text)
    wrong_generation["generation_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(wrong_generation), encoding="utf-8")
    with pytest.raises(ValueError, match="generation digest"):
        verify_run_workbook_export(workbook)


def test_search_export_row_shapes_cover_full_audit_schema() -> None:
    rich_record = {
        "canonical_id": "PMID:1",
        "title": "=formula-like title",
        "authors": ["A", "B"],
        "year": 2026,
        "pubmed_ids": ["1", "1", ""],
        "doi": "10.1/example",
        "pmcid": None,
        "geo_accessions": ["GSE1"],
        "provider_ids": {"ENA": ["PRJ1", ""]},
        "sources": [
            {"provider": "ENA", "id": "PRJ1", "source_url": "https://example.test/PRJ1"},
            "non-mapping source",
        ],
        "candidates": [
            {"candidate_id": "c1", "provider": "GEO", "verification_state": "verified"},
            "non-mapping candidate",
        ],
        "data_readiness": {"tier": "ready", "priority": 1},
    }
    list_provider_record = {
        "source_unit_id": "U2",
        "accession": "GSE2",
        "provider_ids": ["ArrayExpress:E-MTAB-1", "RAW-ID", ""],
        "deg_input_assessment": "not-a-mapping",
    }
    snapshot = {
        "query": "audit",
        "species": "human",
        "records": [rich_record, list_provider_record],
        "provider_diagnostics": {"pubmed": "ok", "geo": {"status": "partial"}},
    }

    rows = discovery_export.publication_rows(snapshot)
    assert rows[0]["paper_title"] == "=formula-like title"
    assert rows[0]["readiness_tier"] == "ready"
    assert rows[1]["readiness_tier"] == ""
    assert discovery_export.publication_rows({"studies": [rich_record]})[0]["canonical_id"] == "PMID:1"

    identifiers = discovery_export._identifier_rows(snapshot["records"])
    assert {row["identifier_type"] for row in identifiers} >= {"PMID", "DOI", "GEO", "ENA", "ArrayExpress", "provider"}
    assert len([row for row in identifiers if row["identifier_type"] == "PMID"]) == 1
    datasets = discovery_export._dataset_rows(snapshot["records"])
    assert any(row.get("provider") == "ENA" for row in datasets)
    candidates = discovery_export._candidate_rows(snapshot["records"])
    assert candidates == [
        {
            "canonical_id": "PMID:1",
            "candidate_id": "c1",
            "provider": "GEO",
            "name": "",
            "role": "",
            "source_input_type": "",
            "verification_state": "verified",
            "source_url": "",
            "landing_url": "",
            "license": "",
            "reason": "",
        }
    ]
    assert discovery_export._event_rows(snapshot)[0]["provider"] in {"pubmed", "geo"}
    assert discovery_export._event_rows({"provider_events": ["plain event"]}) == [{"detail": "plain event"}]
    assert discovery_export._as_list(None) == []
    assert discovery_export._as_list("one") == ["one"]
    assert discovery_export._safe_cell(None) == ""
    assert discovery_export._safe_cell(True) is True
    assert discovery_export._safe_cell(3.5) == 3.5
    assert discovery_export._safe_cell({"b": 1}) == '{"b": 1}'
    assert discovery_export._safe_cell(["=cmd", "safe"]).startswith("'")
    assert discovery_export._publication_csv([]).startswith("canonical_id,paper_title")
    assert discovery_export.build_publication_search_workbook(snapshot).startswith(b"PK")


def test_set_publication_rejects_ambiguous_targets_and_cleans_preparation_failures(
    tmp_path, monkeypatch
) -> None:
    first = tmp_path / "first.stage"
    second = tmp_path / "second.stage"
    final = tmp_path / "final.txt"
    first.write_text("first", encoding="utf-8")
    second.write_text("second", encoding="utf-8")
    with pytest.raises(ValueError, match="must be unique"):
        provenance.publish_staged_artifacts({first: final, second: final})
    assert first.exists() and second.exists(), "validation must not consume caller staging files"
    provenance.publish_staged_artifacts({})

    final.mkdir()
    with pytest.raises(IsADirectoryError, match="not a file"):
        provenance.publish_staged_artifacts({first: final})
    assert not first.exists()
    assert not list(tmp_path.glob("*.pending"))

    staged = tmp_path / "fallback.stage"
    target = tmp_path / "fallback.txt"
    staged.write_text("new", encoding="utf-8")
    target.write_text("old", encoding="utf-8")
    monkeypatch.setattr(provenance.os, "link", lambda *_args: (_ for _ in ()).throw(OSError("no links")))
    provenance.publish_staged_artifacts({staged: target})
    assert target.read_text(encoding="utf-8") == "new"
    assert not list(tmp_path.glob("*.backup"))


def test_set_publication_rejects_resolved_aliases_and_source_target_overlap(tmp_path) -> None:
    first = tmp_path / "first.stage"
    second = tmp_path / "second.stage"
    final = tmp_path / "final.txt"
    first.write_text("first", encoding="utf-8")
    second.write_text("second", encoding="utf-8")
    relative_final = Path(os.path.relpath(final, Path.cwd()))

    with pytest.raises(ValueError, match="target paths must be unique"):
        provenance.publish_staged_artifacts({first: final, second: relative_final})
    assert first.read_text(encoding="utf-8") == "first"
    assert second.read_text(encoding="utf-8") == "second"
    assert not final.exists()

    same = tmp_path / "same.txt"
    same.write_text("must survive", encoding="utf-8")
    with pytest.raises(ValueError, match="source and target paths must be different"):
        provenance.publish_staged_artifacts({same: relative_final.with_name("same.txt")})
    assert same.read_text(encoding="utf-8") == "must survive"


def test_set_publication_removes_new_targets_and_temporary_files_after_failure(
    tmp_path, monkeypatch
) -> None:
    staged_one = tmp_path / "one.stage"
    staged_two = tmp_path / "two.stage"
    final_one = tmp_path / "one.txt"
    final_two = tmp_path / "two.txt"
    staged_one.write_text("one", encoding="utf-8")
    staged_two.write_text("two", encoding="utf-8")
    real_replace = provenance.os.replace

    def fail_second(source, destination):
        if Path(destination) == final_two and Path(source).suffix == ".pending":
            raise OSError("synthetic second publication failure")
        return real_replace(source, destination)

    monkeypatch.setattr(provenance.os, "replace", fail_second)
    with pytest.raises(OSError, match="second publication"):
        provenance.publish_staged_artifacts({staged_one: final_one, staged_two: final_two})

    assert not final_one.exists() and not final_two.exists()
    assert not staged_one.exists() and not staged_two.exists()
    assert not list(tmp_path.glob("*.pending"))
    assert not list(tmp_path.glob("*.backup"))


def test_serve_authentication_configuration_covers_secure_opt_out_and_network_paths(
    tmp_path, monkeypatch, capsys
) -> None:
    db = tmp_path / "scores.db"
    _minimal_score_db(db)

    with pytest.raises(ValueError, match="must not be empty"):
        api.serve(db, access_token="", quiet=True)
    with pytest.raises(PermissionError, match="without --allow-network"):
        api.serve(db, host="0.0.0.0", quiet=True)

    created: list[dict[str, object]] = []

    class InterruptingServer:
        def __init__(self, address):
            self.server_address = address
            self.db_path = db
            self.closed = False

        def serve_forever(self):
            raise KeyboardInterrupt

        def server_close(self):
            self.closed = True

    def fake_create_server(*_args, **kwargs):
        address = (kwargs["host"], 8765)
        server = InterruptingServer(address)
        created.append({"kwargs": kwargs, "server": server})
        return server

    monkeypatch.setattr(api, "create_server", fake_create_server)
    api.serve(db, host="::1", authenticate_loopback=False, quiet=True)
    assert created[-1]["kwargs"]["access_token"] is None
    assert created[-1]["server"].closed
    assert "http://[::1]:8765" in capsys.readouterr().out

    def browser_failure(_url):
        raise RuntimeError("desktop unavailable")

    monkeypatch.setattr(api.webbrowser, "open", browser_failure)
    api.serve(db, host="0.0.0.0", allow_network=True, open_browser=True, quiet=True)
    generated_token = created[-1]["kwargs"]["access_token"]
    assert isinstance(generated_token, str) and generated_token
    captured = capsys.readouterr()
    assert "#token=" in captured.out
    assert "WARNING: serving on 0.0.0.0" in captured.err
    assert "Could not open a browser automatically" in captured.err
