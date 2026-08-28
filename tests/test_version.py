from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import degora

ROOT = Path(__file__).resolve().parents[1]


def test_v0436_keeps_the_v13_scoring_and_aggregation_contract() -> None:
    assert degora.__version__ == "0.4.36"
    assert degora.SCORE_VERSION == "degora_score_v1_3_source_unit_mean"


def _place_module_in_checkout(tmp_path: Path, monkeypatch) -> Path:
    checkout = tmp_path / "checkout"
    (checkout / ".git").mkdir(parents=True)
    module_path = checkout / "degora" / "__init__.py"
    module_path.parent.mkdir()
    module_path.touch()
    monkeypatch.setattr(degora, "__file__", str(module_path))
    return checkout


def test_runtime_version_marks_tracked_worktree_changes_dirty(tmp_path, monkeypatch) -> None:
    checkout = _place_module_in_checkout(tmp_path, monkeypatch)
    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append(command)
        if "rev-parse" in command:
            return subprocess.CompletedProcess(command, 0, stdout="abc1234\n", stderr="")
        return subprocess.CompletedProcess(command, 0, stdout=" M degora/api.py\n", stderr="")

    monkeypatch.setattr(degora.subprocess, "run", fake_run)

    info = degora.runtime_version_info()

    assert info["degora_code_revision"] == "abc1234-dirty"
    assert info["degora_code_dirty"] == "true"
    assert degora.format_version_info(info).endswith("(abc1234-dirty)")
    assert commands == [
        [
            "git",
            "-C",
            str(checkout),
            "ls-files",
            "--error-unmatch",
            "--",
            "degora/__init__.py",
        ],
        ["git", "-C", str(checkout), "rev-parse", "--short", "HEAD"],
        ["git", "-C", str(checkout), "status", "--porcelain=v1", "--untracked-files=no"],
    ]


def test_runtime_version_keeps_clean_revision_without_dirty_flag(tmp_path, monkeypatch) -> None:
    _place_module_in_checkout(tmp_path, monkeypatch)

    def fake_run(command, **kwargs):
        stdout = "abc1234\n" if "rev-parse" in command else ""
        return subprocess.CompletedProcess(command, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(degora.subprocess, "run", fake_run)

    info = degora.runtime_version_info()

    assert info["degora_code_revision"] == "abc1234"
    assert "degora_code_dirty" not in info


def test_runtime_version_omits_revision_when_dirty_state_cannot_be_checked(tmp_path, monkeypatch) -> None:
    _place_module_in_checkout(tmp_path, monkeypatch)

    def fake_run(command, **kwargs):
        if "rev-parse" in command:
            return subprocess.CompletedProcess(command, 0, stdout="abc1234\n", stderr="")
        raise subprocess.TimeoutExpired(command, timeout=5)

    monkeypatch.setattr(degora.subprocess, "run", fake_run)

    info = degora.runtime_version_info()

    assert "degora_code_revision" not in info
    assert "degora_code_dirty" not in info


def test_source_checkout_version_is_not_shadowed_by_other_distribution_metadata(tmp_path) -> None:
    fake_site = tmp_path / "fake-site"
    dist_info = fake_site / "degora-9.9.9.dist-info"
    dist_info.mkdir(parents=True)
    (dist_info / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: degora\nVersion: 9.9.9\n",
        encoding="utf-8",
    )
    env = {
        **os.environ,
        "PYTHONPATH": os.pathsep.join([str(fake_site), str(ROOT)]),
    }
    script = """
import json
from importlib import metadata
from pathlib import Path

import degora

print(json.dumps({
    "distribution_version": metadata.version("degora"),
    "module_file": str(Path(degora.__file__).resolve()),
    "runtime": degora.runtime_version_info(),
}))
"""

    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env=env,
        timeout=10,
    )
    payload = json.loads(completed.stdout)

    assert payload["distribution_version"] == "9.9.9"
    assert payload["module_file"] == str((ROOT / "degora" / "__init__.py").resolve())
    assert payload["runtime"]["degora_version"] == degora.__version__


def test_runtime_version_omits_parent_repo_revision_for_untracked_site_package(tmp_path) -> None:
    parent_repo = tmp_path / "parent-repo"
    package_root = parent_repo / "site-packages"
    module_root = package_root / "degora"
    module_root.mkdir(parents=True)
    (module_root / "__init__.py").write_text((ROOT / "degora" / "__init__.py").read_text(encoding="utf-8"))
    subprocess.run(["git", "init"], cwd=parent_repo, check=True, capture_output=True, text=True)
    (parent_repo / "tracked.txt").write_text("tracked\n", encoding="utf-8")
    subprocess.run(["git", "add", "tracked.txt"], cwd=parent_repo, check=True, capture_output=True, text=True)
    subprocess.run(
        ["git", "-c", "user.name=Test", "-c", "user.email=test@example.com", "commit", "-m", "init"],
        cwd=parent_repo,
        check=True,
        capture_output=True,
        text=True,
    )
    script = "import json, degora; print(json.dumps(degora.runtime_version_info(), sort_keys=True))"
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(package_root)},
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    )

    info = json.loads(completed.stdout)
    assert info == {"degora_version": degora.__version__}


def test_runtime_version_omits_revision_when_tracking_check_fails(tmp_path, monkeypatch) -> None:
    _place_module_in_checkout(tmp_path, monkeypatch)

    def fake_run(command, **kwargs):
        if "ls-files" in command:
            raise OSError("git unavailable")
        return subprocess.CompletedProcess(command, 0, stdout="abc1234\n", stderr="")

    monkeypatch.setattr(degora.subprocess, "run", fake_run)

    info = degora.runtime_version_info()

    assert info == {"degora_version": degora.__version__}
