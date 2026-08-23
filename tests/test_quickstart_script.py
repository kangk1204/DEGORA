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
    for option in ("--port", "--dir", "--ref", "--config", "--update", "--no-browser", "--no-demo"):
        assert option in result.stdout
    # The help text stops at the end of the comment banner.
    assert "set -euo pipefail" not in result.stdout


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (["--port", "not-a-number"], "whole number"),
        (["--port", "0"], "between 1 and 65535"),
        (["--port"], "needs a value"),
        (["--ref"], "needs a value"),
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


# --- --ref behaviour -------------------------------------------------------
# `--ref` is the difference between a reviewer running the branch under test and
# a reviewer silently running the default branch, so these exercise real git.


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "git",
            "-c",
            "user.email=test@example.com",
            "-c",
            "user.name=Test",
            "-c",
            "commit.gpgsign=false",
            *args,
        ],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _origin_with_a_branch(tmp_path: Path) -> tuple[Path, str, str]:
    """An origin whose `feature` branch is one commit ahead of `main`."""

    origin = tmp_path / "origin"
    origin.mkdir()
    _git(origin, "init", "--quiet", "--initial-branch=main", ".")
    _write(origin / "pyproject.toml", '[project]\nname = "degora"\nversion = "0"\n')
    _write(origin / "payload.txt", "base\n")
    _git(origin, "add", "-A")
    _git(origin, "commit", "--quiet", "-m", "base")
    base = _git(origin, "rev-parse", "HEAD").stdout.strip()
    _git(origin, "checkout", "--quiet", "-b", "feature")
    _write(origin / "payload.txt", "feature\n")
    _git(origin, "commit", "--quiet", "-am", "feature")
    feature = _git(origin, "rev-parse", "HEAD").stdout.strip()
    _git(origin, "checkout", "--quiet", "main")
    return origin, base, feature


def _checkout_running_the_script(tmp_path: Path, origin: Path) -> Path:
    repo = tmp_path / "checkout"
    subprocess.run(["git", "clone", "--quiet", str(origin), str(repo)], check=True)
    # Untracked, so switching branches never rewrites the script while it runs.
    (repo / "scripts").mkdir(exist_ok=True)
    shutil.copy2(SCRIPT, repo / "scripts" / "degora_quickstart.sh")
    return repo


def _env_that_stops_after_git(tmp_path: Path) -> dict[str, str]:
    """A PATH whose python passes the version gate but cannot build a venv.

    The script then exits right after the git step, which is the part under
    test, without spending a minute on pip.
    """

    bin_dir = tmp_path / "stub-bin"
    bin_dir.mkdir(exist_ok=True)
    for tool in (
        "git",
        "uname",
        "grep",
        "awk",
        "sed",
        "dirname",
        "mktemp",
        "rm",
        "mkdir",
        "cat",
        "tee",
        "bash",
    ):
        source = shutil.which(tool)
        if source and not (bin_dir / tool).exists():
            (bin_dir / tool).symlink_to(source)
    stub = bin_dir / "python3"
    stub.write_text(
        "#!/bin/sh\n"
        'case "$1" in\n'
        '  --version) echo "Python 3.12.0"; exit 0 ;;\n'
        "  -c) exit 0 ;;\n"
        '  -m) [ "$2" = venv ] && { echo "stub: venv disabled" >&2; exit 1; }; exit 0 ;;\n'
        "esac\nexit 0\n",
        encoding="utf-8",
    )
    stub.chmod(0o755)
    return dict(os.environ, PATH=str(bin_dir))


def _run_quickstart(repo: Path, env: dict[str, str], *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [_bash(), str(repo / "scripts" / "degora_quickstart.sh"), *args],
        capture_output=True,
        text=True,
        cwd=repo,
        env=env,
    )


def test_ref_checks_out_the_requested_branch(tmp_path: Path) -> None:
    origin, base, feature = _origin_with_a_branch(tmp_path)
    repo = _checkout_running_the_script(tmp_path, origin)
    assert _git(repo, "rev-parse", "HEAD").stdout.strip() == base

    result = _run_quickstart(repo, _env_that_stops_after_git(tmp_path), "--ref", "feature")

    assert "Checking out feature" in result.stdout, result.stdout + result.stderr
    assert _git(repo, "rev-parse", "HEAD").stdout.strip() == feature
    assert (repo / "payload.txt").read_text(encoding="utf-8") == "feature\n"


def test_ref_fast_forwards_a_stale_local_branch(tmp_path: Path) -> None:
    """A reviewer who ran the branch last week must not silently serve old code."""

    origin, base, feature = _origin_with_a_branch(tmp_path)
    repo = _checkout_running_the_script(tmp_path, origin)
    _git(repo, "checkout", "--quiet", "-b", "feature", base)
    assert _git(repo, "rev-parse", "HEAD").stdout.strip() == base

    result = _run_quickstart(repo, _env_that_stops_after_git(tmp_path), "--ref", "feature")

    assert _git(repo, "rev-parse", "HEAD").stdout.strip() == feature, result.stdout + result.stderr


def test_ref_refuses_to_run_a_diverged_local_branch(tmp_path: Path) -> None:
    origin, base, _feature = _origin_with_a_branch(tmp_path)
    repo = _checkout_running_the_script(tmp_path, origin)
    _git(repo, "checkout", "--quiet", "-b", "feature", base)
    _write(repo / "payload.txt", "local divergence\n")
    _git(repo, "commit", "--quiet", "-am", "local work")
    diverged = _git(repo, "rev-parse", "HEAD").stdout.strip()

    result = _run_quickstart(repo, _env_that_stops_after_git(tmp_path), "--ref", "feature")

    assert result.returncode != 0
    assert "diverged" in result.stderr
    # The local work is still there; nothing was discarded to force the ref.
    assert _git(repo, "rev-parse", "HEAD").stdout.strip() == diverged


def test_ref_reports_an_unknown_name_instead_of_serving_the_wrong_code(tmp_path: Path) -> None:
    origin, base, _feature = _origin_with_a_branch(tmp_path)
    repo = _checkout_running_the_script(tmp_path, origin)

    result = _run_quickstart(repo, _env_that_stops_after_git(tmp_path), "--ref", "no-such-branch")

    assert result.returncode != 0
    assert "could not fetch" in result.stderr
    assert _git(repo, "rev-parse", "HEAD").stdout.strip() == base


def test_quickstart_reuses_an_existing_demo_workspace_instead_of_deleting_it() -> None:
    """"Safe to re-run" has to mean it does not delete the reader's work.

    The script used to `rm -rf degora-demo` before rebuilding, so a config the
    reader had edited there, or results they had kept, went with it - while the
    README called re-running safe.
    """

    script = (Path(__file__).resolve().parents[1] / "scripts" / "degora_quickstart.sh").read_text(encoding="utf-8")

    assert "rm -rf degora-demo" not in script
    assert "rm -rf $DEMO_DIR" not in script
    assert 'rm -rf "$DEMO_DIR"' not in script
    # An existing workspace is reused, and a separate one can be asked for by name.
    assert 'if [ -d "$DEMO_DIR" ]; then' in script
    assert "Reusing the existing demo workspace" in script
    assert "--demo-dir" in script

    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text(encoding="utf-8")
    assert "reused, never\ndeleted" in readme
    assert "--demo-dir NAME" in readme
