"""Security regressions for authenticated browser auto-launch."""

from __future__ import annotations

import html
import os
import re
import sqlite3
import threading
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from degora import api


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


def _health(url: str, token: str) -> tuple[int, dict]:
    request = urllib.request.Request(url, headers={"X-DEGORA-Token": token})
    with urllib.request.urlopen(request, timeout=5) as response:
        import json

        return response.status, json.loads(response.read().decode("utf-8"))


@pytest.mark.skipif(os.name == "nt", reason="Linux/WSL owner-only bootstrap contract is POSIX-specific")
def test_private_bootstrap_argv_is_nonsecret_owner_only_and_authenticates(tmp_path, monkeypatch) -> None:
    db = tmp_path / "scores.db"
    _minimal_score_db(db)
    token = "live-capability-for-regression"
    server = api.create_server(db, port=0, quiet=True, access_token=token)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address[:2]
    target = f"http://{host}:{port}#token={urllib.parse.quote(token, safe='')}"
    captured: dict[str, object] = {}

    def fake_open(opened_url: str) -> bool:
        captured["opened_url"] = opened_url
        parsed = urllib.parse.urlsplit(opened_url)
        bootstrap = Path(urllib.request.url2pathname(parsed.path))
        captured["bootstrap"] = bootstrap
        captured["document"] = bootstrap.read_text(encoding="utf-8")
        captured["file_mode"] = bootstrap.stat().st_mode & 0o777
        captured["directory_mode"] = bootstrap.parent.stat().st_mode & 0o777
        return True

    monkeypatch.setattr(api.webbrowser, "open", fake_open)
    cleanup = api._open_browser_with_private_bootstrap(target)
    try:
        opened_url = str(captured["opened_url"])
        parsed_opened = urllib.parse.urlsplit(opened_url)
        assert parsed_opened.scheme == "file"
        assert not parsed_opened.query and not parsed_opened.fragment
        assert token not in opened_url
        assert urllib.parse.quote(token, safe="") not in opened_url
        assert captured["file_mode"] == 0o600
        assert captured["directory_mode"] == 0o700

        document = str(captured["document"])
        match = re.search(r'<meta http-equiv="refresh" content="0;url=([^"]+)">', document)
        assert match is not None
        redirect = html.unescape(match.group(1))
        redirect_parts = urllib.parse.urlsplit(redirect)
        supplied = urllib.parse.parse_qs(redirect_parts.fragment)["token"][0]
        assert supplied == token
        status, health = _health(
            urllib.parse.urlunsplit((redirect_parts.scheme, redirect_parts.netloc, "/api/health", "", "")),
            supplied,
        )
        assert status == 200 and health["status"] == "ok"
    finally:
        bootstrap = Path(captured["bootstrap"])
        cleanup()
        cleanup()
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert not bootstrap.exists()
    assert not bootstrap.parent.exists()


def test_private_bootstrap_has_an_independent_bounded_cleanup_timer(monkeypatch, tmp_path) -> None:
    captured: dict[str, object] = {}

    class FakeTimer:
        daemon = False

        def __init__(self, interval: float, callback) -> None:
            captured["interval"] = interval
            captured["callback"] = callback

        def start(self) -> None:
            captured["started"] = True

        def cancel(self) -> None:
            captured["cancelled"] = True

    def fake_open(opened_url: str) -> bool:
        bootstrap = Path(urllib.request.url2pathname(urllib.parse.urlsplit(opened_url).path))
        captured["bootstrap"] = bootstrap
        return True

    monkeypatch.setattr(api.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(api.threading, "Timer", FakeTimer)
    monkeypatch.setattr(api.webbrowser, "open", fake_open)

    cleanup = api._open_browser_with_private_bootstrap("http://127.0.0.1:8765#token=secret")
    bootstrap = Path(captured["bootstrap"])
    assert captured["interval"] == api.BROWSER_BOOTSTRAP_LIFETIME_SECONDS
    assert captured["started"] is True
    assert bootstrap.exists()

    callback = captured["callback"]
    callback()
    assert not bootstrap.exists()
    assert not bootstrap.parent.exists()
    cleanup()


def test_native_windows_auto_open_fails_closed_without_writing_token_material(
    monkeypatch, tmp_path
) -> None:
    monkeypatch.setattr(api.os, "name", "nt")
    monkeypatch.setattr(api.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(
        api.webbrowser,
        "open",
        lambda _url: pytest.fail("native Windows must not receive an unverified secret bootstrap"),
    )

    with pytest.raises(PermissionError, match="owner-only ACLs cannot be guaranteed"):
        api._open_browser_with_private_bootstrap("http://127.0.0.1:8765#token=secret")

    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    ("has_wslview", "expected_target"),
    [
        (True, "file://wsl.localhost/Ubuntu-24.04/tmp/degora-browser/open-dashboard.html"),
        (False, r"\\wsl.localhost\Ubuntu-24.04\tmp\degora-browser\open-dashboard.html"),
    ],
)
def test_wsl_launcher_receives_a_windows_readable_nonsecret_bootstrap(
    monkeypatch, tmp_path, has_wslview, expected_target
) -> None:
    token = "wsl-live-capability"
    captured: dict[str, object] = {}

    class Conversion:
        stdout = r"\\wsl.localhost\Ubuntu-24.04\tmp\degora-browser\open-dashboard.html" + "\n"

    def fake_which(name: str) -> str | None:
        available = {
            "wslpath": "/usr/bin/wslpath",
            "explorer.exe": "/mnt/c/Windows/explorer.exe",
        }
        if has_wslview:
            available["wslview"] = "/usr/bin/wslview"
        return available.get(name)

    def fake_run(argv, **kwargs):
        captured["conversion_argv"] = argv
        captured["conversion_kwargs"] = kwargs
        return Conversion()

    def fake_popen(argv, **kwargs):
        captured["launch_argv"] = argv
        captured["launch_kwargs"] = kwargs
        return object()

    monkeypatch.setattr(api.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(api, "_running_under_wsl", lambda: True)
    monkeypatch.setattr(api.shutil, "which", fake_which)
    monkeypatch.setattr(api.subprocess, "run", fake_run)
    monkeypatch.setattr(api.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(api.webbrowser, "open", lambda _url: pytest.fail("WSL must use a Windows-readable path"))

    cleanup = api._open_browser_with_private_bootstrap(f"http://127.0.0.1:8765#token={token}")
    try:
        conversion_argv = captured["conversion_argv"]
        launch_argv = captured["launch_argv"]
        assert conversion_argv[:2] == ["/usr/bin/wslpath", "-w"]
        assert launch_argv[1] == expected_target
        assert token not in " ".join(conversion_argv)
        assert token not in " ".join(launch_argv)
        bootstrap = Path(conversion_argv[2])
        assert bootstrap.read_text(encoding="utf-8").count(token) == 2
        assert bootstrap.stat().st_mode & 0o777 == 0o600
        assert bootstrap.parent.stat().st_mode & 0o777 == 0o700
    finally:
        cleanup()


@pytest.mark.parametrize("failure", [False, RuntimeError("desktop unavailable")])
def test_failed_desktop_launch_does_not_leave_token_material(monkeypatch, tmp_path, failure) -> None:
    captured: dict[str, Path] = {}

    def fail_open(opened_url: str) -> bool:
        bootstrap = Path(urllib.request.url2pathname(urllib.parse.urlsplit(opened_url).path))
        captured["bootstrap"] = bootstrap
        if isinstance(failure, Exception):
            raise failure
        return failure

    monkeypatch.setattr(api.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(api, "_running_under_wsl", lambda: False)
    monkeypatch.setattr(api.webbrowser, "open", fail_open)

    with pytest.raises(RuntimeError, match="desktop|declined"):
        api._open_browser_with_private_bootstrap("http://127.0.0.1:8765#token=secret")

    bootstrap = captured["bootstrap"]
    assert not bootstrap.exists()
    assert not bootstrap.parent.exists()
