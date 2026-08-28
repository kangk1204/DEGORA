"""0.4.32: residuals of the 0.4.29 audits, closed.

The validate path was silent about an infinite log2 fold change; the
unusable-row sentence contradicted itself when only infinite rows were set
aside; the browser's analysis job could not be stopped and reported two stages;
a run killed with the server looked finished; the workbook was written in
place; ``--min-studies 1`` said nothing about what it gives up; a table too
small for the scale checks passed in silence; two-replicate groups were not
called out; a curl user learnt the action header from a 400.
"""

from __future__ import annotations

import json
import threading
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd
import pytest

from degora.api import INDEX_HTML
from degora.harmonize import TableMapping, harmonize_frame

META = {"study_id": "TEST001", "paper_id": "PAPER001", "pipeline": "DESeq2", "n_ctrl": 3, "n_treat": 3}


def _write_table(path: Path, lfc: list) -> None:
    pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(len(lfc))],
            "log2FoldChange": lfc,
            "pvalue": [0.001 + 0.01 * i for i in range(len(lfc))],
            "padj": [0.003 + 0.01 * i for i in range(len(lfc))],
        }
    ).to_csv(path, index=False)


def _catalog(path: Path, source: Path, **extra) -> Path:
    row = {
        "study_id": "S1", "paper_id": "P1", "source_path": str(source), "gene_column": "gene",
        "lfc_column": "log2FoldChange", "p_column": "pvalue", "padj_column": "padj", **extra,
    }
    pd.DataFrame([row]).to_csv(path, index=False)
    return path


def test_validate_says_when_an_effect_column_holds_inf(tmp_path: Path) -> None:
    from degora.slice_runner import validate_catalog_inputs

    source = tmp_path / "deg.csv"
    _write_table(source, [2.0, "inf", -1.5, 0.4, 1.1, -0.7, 0.9, 1.8, -2.2, 0.3, 0.6, -0.1])
    report = validate_catalog_inputs(_catalog(tmp_path / "c.csv", source))
    warnings = " ".join(map(str, report.get("warnings", [])))
    assert "holds 1 infinite value(s) (examples: inf)" in warnings
    assert "sets those rows aside as unusable" in warnings


def test_only_infinite_rows_dropped_reads_as_one_sentence() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["TP53", "VEGFA", "HK2", "RPL13A"],
            "log2FoldChange": ["inf", 2.0, -1.5, 0.1],
            "pvalue": [0.001, 1e-8, 1e-4, 0.5],
            "padj": [0.003, 1e-6, 1e-3, 0.8],
        }
    )
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), META)
    warning = str(out["unusable_row_warning"].iloc[0])
    assert warning.startswith("TEST001: 1 row of 4 has an infinite log2 fold change")
    assert "dropped before ranking" not in warning  # no empty generic sentence appended
    assert "The dropped cells were empty" not in warning


def test_the_analysis_reports_its_stages_and_a_stop_rolls_the_folder_back(tmp_path: Path) -> None:
    from degora.discovery_run import run_discovery_analysis

    stages: list[tuple[float, str]] = []

    class Stop(BaseException):
        pass

    def progress(fraction: float, message: str) -> None:
        stages.append((fraction, message))
        if len(stages) == 1:
            raise Stop()

    from test_discovery_run import _prepared_bundle

    (tmp_path / "bundle").mkdir()
    prepared = _prepared_bundle(tmp_path / "bundle")
    output = tmp_path / "run"
    with pytest.raises(Stop):
        run_discovery_analysis(prepared, [{"candidate_id": "x", "mode": "author"}], output, species="human", progress=progress)
    assert not output.exists()  # rolled back like any other failure
    assert stages and stages[0][1].startswith("Deriving contrast 1 of 1")


def test_the_analysis_job_reports_every_stage_through_the_job(tmp_path: Path) -> None:
    from degora.api import DegoraRequestHandler
    from degora.discovery_store import DiscoveryJobManager, DiscoveryStateStore

    store = DiscoveryStateStore(tmp_path / "discovery")
    manager = DiscoveryJobManager(store, max_workers=1)
    handler = object.__new__(DegoraRequestHandler)
    handler.server = type("Server", (), {"discovery_job_manager": manager})()
    seen: list[str] = []

    def fake_analyze(payload, progress=None, before_publish=None):
        progress(0.3, "Deriving contrast 1 of 2.")
        progress(0.62, "Scoring 2 contrast(s) across 2 source unit(s).")
        seen.append("ran")
        assert before_publish is not None
        before_publish()
        return {"run_id": "r1"}

    handler._discovery_analyze = fake_analyze
    started = handler._discovery_analyze_job({"bundle_id": "a" * 16, "species": "human", "selections": []})
    manager.shutdown(wait=True)
    job = store.get_job(started["job_id"])
    assert job["status"] == "completed" and seen == ["ran"]


def test_the_browser_can_stop_an_analysis() -> None:
    assert '<button type="button" class="job-cancel" id="cancelAnalysisJob" hidden>Stop this analysis</button>' in INDEX_HTML
    assert "async function cancelAnalysisJob()" in INDEX_HTML
    # A direct listener: the button sits in the footer, outside the candidate
    # list whose delegated listener handles the other cancels (found live).
    assert '$("cancelAnalysisJob").addEventListener("click", () => { void cancelAnalysisJob(); });' in INDEX_HTML
    assert 'event.target.closest("#cancelAnalysisJob")' not in INDEX_HTML
    assert "state.analysisCancelled = true;" in INDEX_HTML
    assert "Analysis stopped; nothing was written." in INDEX_HTML
    assert "stop.hidden = !state.analyzing;" in INDEX_HTML


def test_a_run_folder_killed_with_the_server_is_labelled(tmp_path: Path) -> None:
    from degora.discovery_run import DISCOVERY_RUN_INTERRUPTED_MARKER, DISCOVERY_RUN_MARKER, mark_unfinished_discovery_runs

    root = tmp_path / "discovery"
    finished = root / "human" / "runs" / "aaaa"
    unfinished = root / "human" / "runs" / "bbbb"
    finished.mkdir(parents=True)
    unfinished.mkdir(parents=True)
    (finished / DISCOVERY_RUN_MARKER).write_text("{}", encoding="utf-8")
    marked = mark_unfinished_discovery_runs(root)
    assert marked == [unfinished]
    note = json.loads((unfinished / DISCOVERY_RUN_INTERRUPTED_MARKER).read_text(encoding="utf-8"))
    assert "run the analysis again" in note["note"]
    assert not (finished / DISCOVERY_RUN_INTERRUPTED_MARKER).exists()
    assert mark_unfinished_discovery_runs(root) == []  # idempotent


def test_the_workbook_set_is_staged_then_published_with_manifest_last() -> None:
    from degora import excel_export

    source = Path(excel_export.__file__).read_text(encoding="utf-8")
    assert "TemporaryDirectory" in source
    assert "publish_staged_artifacts(publication_pairs)" in source
    assert "publication_pairs[staged_manifest] = manifest" in source
    assert source.index("publication_pairs[staged_manifest] = manifest") < source.index(
        "publish_staged_artifacts(publication_pairs)"
    )


def test_min_studies_one_is_named_exploratory(tmp_path: Path) -> None:
    from degora.slice_runner import run_slice

    source = tmp_path / "deg.csv"
    _write_table(source, [2.0, -1.5, 0.4, 1.1, -0.7, 0.9, 1.8, -2.2, 0.3, 0.6, -0.1, 1.4])
    metrics = run_slice(_catalog(tmp_path / "c.csv", source), tmp_path / "out", tmp_path / "h", min_studies=1)
    joined = " ".join(map(str, metrics.get("input_warnings", []) or metrics.get("warnings", [])))
    assert "min_studies=1" in joined and "exploratory prioritisation" in joined


def test_a_table_too_small_for_the_scale_checks_is_said_so(tmp_path: Path) -> None:
    from degora.slice_runner import validate_catalog_inputs

    source = tmp_path / "deg.csv"
    _write_table(source, [2.0, 1.5, 0.4, 1.1, 0.7, 0.9])  # six rows, all positive: unchecked before
    report = validate_catalog_inputs(_catalog(tmp_path / "c.csv", source))
    warnings = " ".join(map(str, report.get("warnings", [])))
    assert "has only 6 numeric value(s), fewer than the 10" in warnings
    report = validate_catalog_inputs(_catalog(tmp_path / "c2.csv", source, lfc_scale="log2"))
    assert "fewer than the 10" not in " ".join(map(str, report.get("warnings", [])))


def test_two_replicate_groups_are_called_exploratory(tmp_path: Path) -> None:
    from test_discovery_run import _fallback_entry, _matrix_candidate

    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, {"c1": [10.0] * 40, "c2": [12.0] * 40, "t1": [30.0] * 40, "t2": [33.0] * 40})
    _row, summary = _fallback_row(
        study=study, candidate=candidate, entry=_fallback_entry(matrix_type="count_matrix", contrast_label="drug vs ctrl"),
        spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
        sequence=1, replay_command="degora",
    )
    assert "compares 2 control against 2 treatment samples" in summary["small_group_note"]


def test_the_search_api_refuses_non_text_over_http(tmp_path: Path) -> None:
    from degora.api import create_server

    db = tmp_path / "degora.db"
    db.touch()
    server = create_server(db, port=0, quiet=True, discovery_root=tmp_path / "discovery")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    try:
        for bad in (["Placenta"], {"q": "x"}, True, 12, None):
            body = json.dumps({"query": bad, "species": "human"}).encode("utf-8")
            request = urllib.request.Request(
                f"http://{host}:{port}/api/discovery/searches", data=body,
                headers={"Content-Type": "application/json", "X-DEGORA-Action": "1"},
            )
            with pytest.raises(urllib.error.HTTPError) as caught:
                urllib.request.urlopen(request, timeout=10)
            assert caught.value.code == 400
            message = json.loads(caught.value.read().decode("utf-8"))["error"]
            assert "must be a string" in message or "at least 2 characters" in message
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_help_and_readme_name_the_boundaries() -> None:
    from degora.cli import build_parser

    parser = build_parser()
    text = parser.format_help()
    assert "discovery-analyze" in text
    readme = Path(__file__).resolve().parents[1].joinpath("README.md").read_text(encoding="utf-8")
    assert "X-DEGORA-Action: 1" in readme
    assert "exploratory prioritisation, not replicated evidence" in readme
    assert "at least 10 numeric values" in readme
    assert "reviewer_attestations" in readme
