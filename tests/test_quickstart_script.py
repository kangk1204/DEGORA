from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "degora_quickstart.sh"


def _bash() -> str:
    found = shutil.which("bash")
    if not found:  # pragma: no cover - bash is present on every supported platform.
        pytest.skip("bash is not available")
    return found


def test_quickstart_script_exists_and_is_executable() -> None:
    assert SCRIPT.is_file()
    mode = SCRIPT.stat().st_mode
    assert mode & stat.S_IXUSR, "scripts/degora_quickstart.sh must be executable"
    assert SCRIPT.read_text(encoding="utf-8").startswith("#!/usr/bin/env bash\n")


def test_quickstart_script_parses() -> None:
    result = subprocess.run([_bash(), "-n", str(SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_quickstart_help_lists_every_documented_option() -> None:
    result = subprocess.run([_bash(), str(SCRIPT), "--help"], capture_output=True, text=True, cwd=ROOT)
    assert result.returncode == 0, result.stderr
    for option in ("--port", "--dir", "--config", "--update", "--no-browser", "--no-demo"):
        assert option in result.stdout
    # The help text stops at the end of the comment banner.
    assert "set -euo pipefail" not in result.stdout


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (["--port", "not-a-number"], "whole number"),
        (["--port", "0"], "between 1 and 65535"),
        (["--port"], "needs a value"),
        (["--nonsense"], "unknown option"),
    ],
)
def test_quickstart_rejects_bad_arguments(args: list[str], expected: str) -> None:
    result = subprocess.run([_bash(), str(SCRIPT), *args], capture_output=True, text=True, cwd=ROOT)
    assert result.returncode != 0
    assert expected in result.stderr


def test_quickstart_reports_a_missing_interpreter_instead_of_failing_late(tmp_path: Path) -> None:
    """With no supported interpreter on PATH the script must stop with guidance."""

    empty_bin = tmp_path / "bin"
    empty_bin.mkdir()
    for tool in ("git", "uname", "grep", "awk", "sed", "mktemp", "rm", "mkdir", "cat", "tee"):
        source = shutil.which(tool)
        if source:
            (empty_bin / tool).symlink_to(source)
    env = dict(os.environ, PATH=str(empty_bin))
    result = subprocess.run([_bash(), str(SCRIPT)], capture_output=True, text=True, cwd=tmp_path, env=env)
    assert result.returncode == 1
    assert "no Python 3.10 or newer" in result.stderr


def test_quickstart_avoids_constructs_that_break_on_macos_bash_3_2() -> None:
    """macOS ships bash 3.2, so bash 4+ syntax would fail there but pass in CI."""

    text = SCRIPT.read_text(encoding="utf-8")
    for construct in ("declare -A", "mapfile", "readarray", "${BASH_VERSINFO", "&>>"):
        assert construct not in text, f"{construct} is not available in bash 3.2"
    # `readlink -f` is GNU-only; BSD readlink on macOS does not support it.
    assert "readlink -f" not in text


def test_quickstart_handles_every_supported_platform() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    # Platform detection plus one browser opener per platform.
    assert "PLATFORM=macos" in text and "PLATFORM=wsl" in text and "PLATFORM=linux" in text
    assert "open " in text and "xdg-open" in text and "wslview" in text
    # Debian/Ubuntu splits venv into its own package; the failure must say so.
    assert "-venv" in text


def test_readme_documents_the_quickstart_script() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "scripts/degora_quickstart.sh" in readme
    assert "--no-browser" in readme
