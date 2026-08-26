"""DEGORA thin-slice implementation and local score browser."""

from __future__ import annotations

import subprocess
from pathlib import Path

__all__ = ["SCORE_VERSION", "__version__", "format_version_info", "runtime_version_info"]

__version__ = "0.4.30"
SCORE_VERSION = "degora_score_v1_2_source_unit_mean"


def _is_tracked_by_repo(repo: Path, module_path: Path) -> bool:
    try:
        relative_module = module_path.relative_to(repo)
    except ValueError:
        return False
    try:
        subprocess.run(
            ["git", "-C", str(repo), "ls-files", "--error-unmatch", "--", str(relative_module)],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return True


def _git_revision() -> str:
    module_path = Path(__file__).resolve()
    for parent in module_path.parents:
        if not (parent / ".git").exists():
            continue
        if not _is_tracked_by_repo(parent, module_path):
            continue
        try:
            revision_result = subprocess.run(
                ["git", "-C", str(parent), "rev-parse", "--short", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            )
            status_result = subprocess.run(
                ["git", "-C", str(parent), "status", "--porcelain=v1", "--untracked-files=no"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            )
        except (OSError, subprocess.SubprocessError):
            return ""
        revision = revision_result.stdout.strip()
        if not revision:
            return ""
        if status_result.stdout.strip():
            return f"{revision}-dirty"
        return revision
    return ""


def runtime_version_info() -> dict[str, str]:
    info = {"degora_version": __version__}
    revision = _git_revision()
    if revision:
        info["degora_code_revision"] = revision
        if revision.endswith("-dirty"):
            info["degora_code_dirty"] = "true"
    return info


def format_version_info(info: dict[str, str] | None = None) -> str:
    version_info = info or runtime_version_info()
    text = version_info["degora_version"]
    revision = version_info.get("degora_code_revision", "")
    if revision:
        text = f"{text} ({revision})"
    return text
