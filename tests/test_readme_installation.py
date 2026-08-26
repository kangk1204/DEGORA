import re
import sys
from importlib import metadata
from pathlib import Path

import pytest

import degora

ROOT = Path(__file__).resolve().parents[1]


def test_readme_and_installed_package_metadata_are_consistent() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    workflow_path = ROOT / ".github/workflows/tests.yml"
    workflow = workflow_path.read_text(encoding="utf-8") if workflow_path.exists() else ""

    macos_start = readme.index("macOS (including Apple silicon):")
    linux_start = readme.index("Linux:", macos_start)
    macos_install = readme[macos_start:linux_start]

    assert "python3.12 --version" in macos_install
    assert "python3.12 -m venv .venv" in macos_install
    assert "python3 -m venv .venv" not in macos_install
    assert "requires a different Python" in readme
    assert "Python 3.9.6 not in '>=3.10'" in readme
    assert "degora serve outputs/results/degora-run/degora_scores.db" in readme
    assert "degora_discovery/" in readme
    assert "Run separate Human + Mouse searches" in readme
    assert "not separately exercised on a hosted Windows runner" in readme
    assert '"scipy>=1.11.1,<1.18"' in pyproject
    assert "scipy>=1.11.1,<1.18" in requirements.splitlines()
    if workflow:
        assert "scipy==1.11.1" in workflow
        assert "scipy==1.11.0" not in workflow
        assert "actions/checkout@v" not in workflow
        assert "actions/setup-python@v" not in workflow
        assert re.search(r"actions/checkout@[0-9a-f]{40}\b", workflow)
        assert re.search(r"actions/setup-python@[0-9a-f]{40}\b", workflow)
    assert '"pyarrow>=14.0.1,<25"' in pyproject
    assert "pyarrow>=14.0.1,<25" in requirements.splitlines()

    version_match = re.search(r'^version = "([^"]+)"$', pyproject, flags=re.MULTILINE)
    python_match = re.search(r'^requires-python = ">=([0-9.]+)"$', pyproject, flags=re.MULTILINE)
    assert version_match is not None
    assert python_match is not None
    project_version = version_match.group(1)
    assert degora.__version__ == project_version
    # The history moved out of the README; the guarantee that every released
    # version has an entry moved with it.
    assert f"## {project_version}" in changelog
    assert "CHANGELOG.md" in readme
    try:
        installed_version = metadata.version("degora")
    except metadata.PackageNotFoundError:
        pytest.skip("install DEGORA before checking installed package metadata")
    assert installed_version == project_version
    required_python = tuple(int(part) for part in python_match.group(1).split("."))
    assert sys.version_info[: len(required_python)] >= required_python
    assert f"Python {python_match.group(1)} or newer" in readme
    assert "automated release tests cover Python 3.10-3.13" in readme


def test_xls_support_ships_with_the_package() -> None:
    """read_deg_table accepts .xls, so the reader pandas needs for it must install.

    xlrd sat in the dev extra, so an ordinary `pip install degora` advertised a
    format it could not open: a valid legacy .xls failed with an ImportError
    telling the user to install a package the project already depended on.
    """

    import importlib.util

    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    dependencies = pyproject.split("[project.optional-dependencies]", 1)[0]
    assert '"xlrd>=2.0.1,<3"' in dependencies

    from degora.harmonize import (
        read_deg_table,  # noqa: F401 - the .xls branch lives here
    )

    assert importlib.util.find_spec("xlrd") is not None
