from __future__ import annotations

import re
from pathlib import Path

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised by the Python 3.10 CI job
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]


def _requirement_signature(raw: str) -> tuple[str, tuple[str, ...], tuple[str, ...], str, str]:
    requirement = Requirement(raw)
    return (
        canonicalize_name(requirement.name),
        tuple(sorted(requirement.extras)),
        tuple(sorted(str(specifier) for specifier in requirement.specifier)),
        str(requirement.marker or ""),
        requirement.url or "",
    )


def test_requirements_txt_exactly_matches_project_runtime_dependencies() -> None:
    """The convenience requirements file must not drift from wheel metadata."""

    with (ROOT / "pyproject.toml").open("rb") as handle:
        project_dependencies = tomllib.load(handle)["project"]["dependencies"]
    requirements_dependencies = [
        line.strip()
        for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

    assert [_requirement_signature(item) for item in requirements_dependencies] == [
        _requirement_signature(item) for item in project_dependencies
    ]


def test_manifest_declares_every_release_support_file() -> None:
    manifest_lines = {
        line.strip()
        for line in (ROOT / "MANIFEST.in").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }

    assert {
        "include degora_quickstart.sh",
        "include scripts/degora_quickstart.sh",
        "include requirements.txt",
        "include .gitattributes",
        "include .gitignore",
        "recursive-include tests *.py *.json",
        "include tests/test_live_discovery_regressions.py",
    } <= manifest_lines


def test_every_github_action_is_pinned_to_an_immutable_commit() -> None:
    workflow_path = ROOT / ".github" / "workflows" / "tests.yml"
    if not workflow_path.exists():
        # GitHub metadata is deliberately not part of the installable sdist;
        # the release job also runs this suite from the unpacked archive.
        return
    workflow = workflow_path.read_text(encoding="utf-8")
    uses = re.findall(r"^\s*-\s+uses:\s+([^@\s]+)@([^\s#]+)", workflow, flags=re.MULTILINE)

    assert uses, "the test workflow must use at least one action"
    assert all(re.fullmatch(r"[0-9a-f]{40}", revision) for _action, revision in uses), uses
