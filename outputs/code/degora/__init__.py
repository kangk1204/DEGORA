"""DEGORA thin-slice implementation and local score browser."""

from __future__ import annotations

import subprocess
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

__all__ = ["__version__", "format_version_info", "runtime_version_info"]

__version__ = "0.3.1"


def _installed_version() -> str:
    try:
        return version("degora")
    except PackageNotFoundError:
        return __version__


def _git_revision() -> str:
    for parent in Path(__file__).resolve().parents:
        if not (parent / ".git").exists():
            continue
        try:
            completed = subprocess.run(
                ["git", "-C", str(parent), "rev-parse", "--short", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            )
        except (OSError, subprocess.SubprocessError):
            return ""
        return completed.stdout.strip()
    return ""


def runtime_version_info() -> dict[str, str]:
    info = {"degora_version": _installed_version()}
    revision = _git_revision()
    if revision:
        info["degora_code_revision"] = revision
    return info


def format_version_info(info: dict[str, str] | None = None) -> str:
    version_info = info or runtime_version_info()
    text = version_info["degora_version"]
    revision = version_info.get("degora_code_revision", "")
    if revision:
        text = f"{text} ({revision})"
    return text
