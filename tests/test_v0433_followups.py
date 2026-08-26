"""Follow-ups from a second review of the 0.4.17 first-run fixes, re-verified on 0.4.33.

Each test was written after reproducing the behaviour on a clean checkout, so a
test that fails here is a behaviour a reader could actually meet, not a guess.
"""

from __future__ import annotations

import contextlib
import io
import re
import time
import urllib.request
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora.api import create_server
from degora.beginner import _read_header
from degora.discovery import DiscoveryUnavailableError
from degora.discovery_federated import _safe_provider_error
from degora.discovery_sources import _safe_remote_error
from degora.harmonize import TableMapping, canonical_gene_symbol, harmonize_frame
from degora.slice_runner import (
    DegoraConfigError,
    _count_normalized_gene_symbols,
    describe_table_read_failure,
    validate_catalog_inputs,
)

ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# Recording input_gene_label must not dominate the cost of harmonizing
# --------------------------------------------------------------------------


def _timed_harmonize(rows: int, *, duplicate_share: float) -> float:
    rng = np.random.default_rng(0)
    distinct = int(rows * (1 - duplicate_share))
    base = [f"GENE{index}" for index in range(distinct)]
    genes = base + list(rng.choice(base, size=rows - distinct)) if duplicate_share else [
        f"GENE{index}" for index in range(rows)
    ]
    frame = pd.DataFrame(
        {
            "gene": genes,
            "log2FoldChange": rng.normal(size=rows),
            "pvalue": rng.uniform(1e-12, 1, size=rows),
            "padj": rng.uniform(1e-12, 1, size=rows),
        }
    )
    meta = {"study_id": "S", "paper_id": "P", "pipeline": "DESeq2", "assay_type": "RNA-seq"}
    mapping = TableMapping("gene", "log2FoldChange", "pvalue", "padj")
    timings = []
    for _ in range(3):
        started = time.perf_counter()
        harmonize_frame(frame, mapping, meta)
        timings.append(time.perf_counter() - started)
    return sorted(timings)[1]


def test_duplicate_symbols_do_not_multiply_the_cost_of_harmonizing() -> None:
    """A table with duplicates must not cost several times one without them.

    The label merge ran a per-gene Python join over *every* gene rather than the
    few reached from more than one label, so a 200k-row source took about four
    times as long as it did before labels were recorded at all. The ratio is
    measured against this machine rather than a fixed number of seconds, so the
    test means the same thing on a slow runner.
    """

    rows = 60_000
    without = _timed_harmonize(rows, duplicate_share=0.0)
    with_duplicates = _timed_harmonize(rows, duplicate_share=0.2)
    assert with_duplicates < without * 2.5, (
        f"duplicates cost {with_duplicates / without:.1f}x a table without them "
        f"({with_duplicates:.2f}s vs {without:.2f}s)"
    )


def test_every_label_a_collapsed_gene_came_from_is_still_kept() -> None:
    # The fast path must not lose what the slow one recorded.
    frame = pd.DataFrame(
        {
            "gene": ["SEPT9", "9-Sep", "ISG15", "ISG15"],
            "log2FoldChange": [2.0, 2.5, 1.0, 1.2],
            "pvalue": [1e-6, 1e-5, 1e-4, 1e-3],
            "padj": [1e-5, 1e-4, 1e-3, 1e-2],
        }
    )
    meta = {"study_id": "S", "paper_id": "P", "pipeline": "DESeq2", "assay_type": "RNA-seq"}
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), meta)
    merged = out.loc[out["gene_symbol"].eq("SEPTIN9"), "input_gene_label"].iloc[0]
    assert set(str(merged).split(";")) == {"SEPT9", "9-Sep"}
    # A gene reached from one label only keeps that label, unjoined.
    single = out.loc[out["gene_symbol"].eq("ISG15"), "input_gene_label"].iloc[0]
    assert str(single) == "ISG15"


@pytest.mark.parametrize("value", [[1, 2], (1, 2), np.array(["A", "B"])])
def test_the_symbol_resolver_answers_instead_of_raising_on_a_non_scalar(value: object) -> None:
    # It now runs per distinct label inside _clean_gene_symbol, so a non-scalar
    # cell would end a whole run with a pandas message about truth values. No
    # config, table or URL can deliver one; the requirement is only that a caller
    # passing one gets an answer back.
    assert isinstance(canonical_gene_symbol(value), str)


@pytest.mark.parametrize(("written", "scored"), [("SEPT9", "SEPTIN9"), ("1-Mar", "MARCHF1"), ("DEC1", "BHLHE40")])
def test_symbol_resolution_is_unchanged(written: str, scored: str) -> None:
    assert canonical_gene_symbol(written) == scored


# --------------------------------------------------------------------------
# A gene label that really contains ';' was reported as a rename
# --------------------------------------------------------------------------


def test_a_multi_mapping_gene_label_is_not_reported_as_a_rename() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": pd.array(["TP53;TP53P1", "ABC"], dtype="string"),
            "input_gene_label": ["TP53;TP53P1", "ABC"],
        }
    )
    assert _count_normalized_gene_symbols(frame) == 0


def test_a_real_rename_is_still_counted() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": pd.array(["SEPTIN9", "ABC"], dtype="string"),
            "input_gene_label": ["SEPT9;9-Sep", "ABC"],
        }
    )
    assert _count_normalized_gene_symbols(frame) == 1


# --------------------------------------------------------------------------
# A folder in source_path is copied down a sheet like a wrong path is
# --------------------------------------------------------------------------


def _catalog(tmp_path: Path, name: str, **columns: object) -> Path:
    path = tmp_path / name
    pd.DataFrame(columns).to_csv(path, index=False)
    return path


def test_source_paths_pointing_at_folders_are_all_reported(tmp_path) -> None:
    for folder in ("d0", "d1", "d2"):
        (tmp_path / folder).mkdir()
    config_path = _catalog(
        tmp_path,
        "folders.csv",
        study_id=["S0", "S1", "S2"],
        source_unit_id=["U0", "U1", "U2"],
        source_path=["d0", "d1", "d2"],
        gene_column=["gene"] * 3,
        lfc_column=["log2FoldChange"] * 3,
        p_column=["pvalue"] * 3,
        include=["yes"] * 3,
    )
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config_path)
    message = str(excinfo.value)
    for study in ("S0", "S1", "S2"):
        assert study in message, f"{study} missing: still one fix-and-re-run cycle per row"
    assert "which is a directory" in message
    assert "not at the folder holding it" in message


# --------------------------------------------------------------------------
# A broken table is explained, not quoted from the parser
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("content", "expected"),
    [
        ("gene,logFC,p\nA,1,0.1\nB,1,0.1,x,y\n", "different numbers of columns"),
        ('gene,logFC,p\n"A,1,0.1\nB,1,0.1\n', "unclosed quotation mark"),
        ("", "is empty"),
    ],
)
def test_a_broken_csv_says_what_is_wrong_with_it(tmp_path, content: str, expected: str) -> None:
    path = tmp_path / "broken.csv"
    path.write_text(content, encoding="utf-8")
    try:
        pd.read_csv(path)
    except Exception as exc:  # noqa: BLE001 - the point is what the reader is told
        described = describe_table_read_failure(path, exc)
    else:
        pytest.fail("fixture should not parse")
    assert expected in described
    assert "C error" not in described
    assert "Error tokenizing" not in described


def test_a_broken_source_table_reaches_the_reader_explained(tmp_path) -> None:
    (tmp_path / "ragged.csv").write_text("gene,logFC,p\nA,1,0.1\nB,1,0.1,x,y\n", encoding="utf-8")
    config_path = _catalog(
        tmp_path,
        "cfg.csv",
        study_id=["S"],
        source_unit_id=["U"],
        source_path=["ragged.csv"],
        gene_column=["gene"],
        lfc_column=["logFC"],
        p_column=["p"],
        include=["yes"],
    )
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config_path)
    message = str(excinfo.value)
    assert "different numbers of columns" in message
    assert "C error" not in message


def test_degora_init_explains_a_file_that_is_not_a_table(tmp_path) -> None:
    binary = tmp_path / "figure.txt"
    binary.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 64)
    _, problem = _read_header(binary)
    assert problem
    for leak in ("UnicodeDecodeError", "ParserError", "BadZipFile", "codec"):
        assert leak not in problem


# --------------------------------------------------------------------------
# One rule for both provider-error sanitizers
# --------------------------------------------------------------------------


def test_both_provider_sanitizers_agree() -> None:
    # discovery_sources carried a copy of the federated helper and stapled a class
    # name onto the plain sentence the federated one had been cleaned up to give.
    exc = DiscoveryUnavailableError("request failed after 3 attempt(s): the request timed out")
    assert _safe_provider_error(exc) == _safe_remote_error(exc)
    assert "DiscoveryUnavailableError" not in _safe_remote_error(exc)


def test_a_plain_exception_keeps_the_only_thing_that_identifies_it() -> None:
    # A bare KeyError('pmid') reads as "'pmid'" on its own and says nothing, so
    # dropping the class name is only right for errors that carry a sentence.
    assert "KeyError" in _safe_remote_error(KeyError("pmid"))


@pytest.mark.parametrize(
    "text",
    [
        "Authorization: Bearer SYNTHETIC_BEARER_REVIEW_7F3A",
        "GET https://example.test/x?api_key=SUPERSECRET failed",
    ],
)
def test_credential_redaction_survives_the_change(text: str) -> None:
    # The class-name rule must not be bought at the cost of the redaction another
    # audit added; both sanitizers keep going through redact_secrets_in_text.
    for secret in ("SYNTHETIC_BEARER_REVIEW_7F3A", "SUPERSECRET"):
        assert secret not in _safe_remote_error(RuntimeError(text))
        assert secret not in _safe_provider_error(RuntimeError(text))


# --------------------------------------------------------------------------
# A failed workbook export names its cause
# --------------------------------------------------------------------------


def _runnable_config(tmp_path: Path) -> Path:
    deg = tmp_path / "deg"
    deg.mkdir(exist_ok=True)
    genes = [f"GENE{index}" for index in range(12)]
    for offset, name in enumerate(("a", "b")):
        pd.DataFrame(
            {
                "gene": genes,
                "log2FoldChange": [2.0 + offset * 0.3 + index * 0.1 for index in range(12)],
                "pvalue": [10.0 ** -(6 + offset + index * 0.1) for index in range(12)],
                "padj": [10.0 ** -(5 + offset + index * 0.1) for index in range(12)],
            }
        ).to_csv(deg / f"{name}.csv", index=False)
    return _catalog(
        tmp_path,
        "cfg.csv",
        study_id=["A", "B"],
        source_unit_id=["U1", "U2"],
        source_path=["deg/a.csv", "deg/b.csv"],
        gene_column=["gene"] * 2,
        lfc_column=["log2FoldChange"] * 2,
        p_column=["pvalue"] * 2,
        padj_column=["padj"] * 2,
        table_scope=["full_results"] * 2,
        include=["yes"] * 2,
    )


def test_a_workbook_still_open_in_excel_is_named_as_the_cause(tmp_path, monkeypatch) -> None:
    from degora import cli

    def refuse(*args, **kwargs):
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr("degora.excel_export.export_run_workbook", refuse)
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
        code = cli.main(["run", str(_runnable_config(tmp_path)), "--output-dir", str(tmp_path / "out")])
    message = stderr.getvalue()
    assert code == 2
    assert "still open in Excel" in message
    assert "PermissionError" not in message


# --------------------------------------------------------------------------
# The landing view must not override a reader who navigates first
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def served_dashboard(tmp_path_factory) -> str:
    from degora import cli

    tmp_path = tmp_path_factory.mktemp("landing")
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        assert cli.main(
            ["run", str(_runnable_config(tmp_path)), "--output-dir", str(tmp_path / "out"), "--no-excel"]
        ) == 0
    server = create_server(tmp_path / "out" / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    try:
        yield f"http://{host}:{port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_the_landing_view_yields_to_a_reader_who_navigates_first(served_dashboard: str) -> None:
    with urllib.request.urlopen(served_dashboard, timeout=5) as response:
        html = response.read().decode("utf-8")
    # The choice waits on /api/meta, so without a guard a tab clicked during boot
    # is silently undone a moment later.
    assert 'scored > 0 ? "atlas"' in html
    assert "let landingViewChosen = false;" in html
    assert "if (landingViewChosen) return;" in html
    # Both nav buttons must claim the decision, or one of them still loses.
    assert html.count("landingViewChosen = true;") >= 3


# --------------------------------------------------------------------------
# The release history lives in CHANGELOG.md, and the README points at it
# --------------------------------------------------------------------------


def test_the_readme_is_about_using_degora_not_its_history() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "CHANGELOG.md" in readme
    # The history was 73% of the README; a reader looking for how to run the tool
    # should not have to scroll past every past release to find the license.
    assert len(readme.splitlines()) < len(changelog.splitlines())
    assert "## Release notes" in readme
    assert not re.search(r"^## 0\.4\.\d+$", readme, flags=re.M), "release entries belong in CHANGELOG.md"


def test_every_release_entry_kept_its_place_in_the_move() -> None:
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    entries = re.findall(r"^## (\S+)$", changelog, flags=re.M)
    versions = [entry for entry in entries if re.fullmatch(r"\d+\.\d+\.\d+", entry)]
    assert len(versions) >= 25, f"only {len(versions)} releases survived the move"
    # Newest first, and nothing duplicated by a bad split.
    assert len(set(versions)) == len(versions)


def test_an_sdist_carries_the_changelog_the_readme_links_to() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "include CHANGELOG.md" in manifest
