"""0.4.31: findings of an extreme multi-study audit of 0.4.29.

F1  an infinite log2 fold change was accepted as the largest effect in a table;
F2  a JSON list in place of the query text was stringified and searched;
O1  the browser's Analyze held one request open for minutes;
O2  a blank search left no lasting message, a file where the output folder
    belongs surfaced as "[Errno 17] File exists", and a finished run left an
    empty lock file that read as active.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from degora.api import INDEX_HTML, DegoraRequestHandler
from degora.harmonize import TableMapping, harmonize_frame

META = {"study_id": "TEST001", "paper_id": "PAPER001", "pipeline": "DESeq2", "n_ctrl": 3, "n_treat": 3}


def test_an_infinite_log2_fold_change_is_set_aside_and_said() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["TP53", "VEGFA", "HK2", "RPL13A"],
            "log2FoldChange": ["inf", 2.0, -1.5, 0.1],
            "pvalue": [0.001, 1e-8, 1e-4, 0.5],
            "padj": [0.003, 1e-6, 1e-3, 0.8],
        }
    )
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), META)
    assert "TP53" not in set(out["gene_symbol"])
    assert int(out["n_rows_dropped_unusable"].iloc[0]) == 1
    warning = str(out["unusable_row_warning"].iloc[0])
    assert "1 row of 4 has an infinite log2 fold change" in warning
    assert "inf" in warning
    # -Inf and a numeric infinity are the same case.
    frame["log2FoldChange"] = [float("-inf"), 2.0, -1.5, 0.1]
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), META)
    assert "TP53" not in set(out["gene_symbol"])
    assert "infinite log2 fold change" in str(out["unusable_row_warning"].iloc[0])


def test_a_table_without_infinities_has_no_such_warning() -> None:
    frame = pd.DataFrame({"gene": ["A", "B", "C"], "log2FoldChange": [2.0, -1.5, 0.1], "pvalue": [1e-8, 1e-4, 0.5], "padj": [1e-6, 1e-3, 0.8]})
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), META)
    assert str(out["unusable_row_warning"].iloc[0]) == ""


def test_the_search_api_refuses_a_query_or_species_that_is_not_text() -> None:
    with pytest.raises(ValueError, match="query must be a string"):
        DegoraRequestHandler._discovery_create_publication_search(None, {"query": ["Placenta"], "species": "human"})
    with pytest.raises(ValueError, match="query must be a string"):
        DegoraRequestHandler._discovery_create_publication_search(None, {"query": 12, "species": "human"})
    with pytest.raises(ValueError, match="species must be a string"):
        DegoraRequestHandler._discovery_create_publication_search(None, {"query": "placenta", "species": ["human"]})


def test_an_analysis_runs_as_a_job_and_returns_the_run_through_it(tmp_path: Path) -> None:
    from degora.discovery_store import DiscoveryJobManager, DiscoveryStateStore

    store = DiscoveryStateStore(tmp_path / "discovery")
    manager = DiscoveryJobManager(store, max_workers=1)
    handler = object.__new__(DegoraRequestHandler)
    handler.server = type("Server", (), {"discovery_job_manager": manager})()
    seen: list[dict] = []
    handler._discovery_analyze = lambda payload, progress=None: (seen.append(payload) or {"run_id": "r1", "n_source_units": 2})

    request = {"bundle_id": "a" * 16, "species": "human", "selections": [], "species_confirmed": True}
    started = handler._discovery_analyze_job(request)
    assert started["status"] == "queued"
    manager.shutdown(wait=True)
    job = store.get_job(started["job_id"])
    assert job is not None and job["status"] == "completed"
    assert job["result"] == {"run_id": "r1", "n_source_units": 2}
    assert seen == [request]
    # The job wrapper keeps the two checks that used to answer at once.
    with pytest.raises(ValueError, match="bundle_id is invalid"):
        handler._discovery_analyze_job({"bundle_id": "nope", "selections": []})
    with pytest.raises(ValueError, match="selections must be a JSON list"):
        handler._discovery_analyze_job({"bundle_id": "a" * 16, "selections": "x"})


def test_the_browser_polls_the_analysis_job_and_explains_a_dropped_connection() -> None:
    assert '"/api/discovery/analyze-jobs"' in INDEX_HTML
    assert "async function pollAnalysisJob(species, requestId, jobId)" in INDEX_HTML
    assert "state.analysisProgress" in INDEX_HTML and "state.analysisMessage" in INDEX_HTML
    assert "function analysisRunningLabel(state)" in INDEX_HTML
    # The stage message is written after the eligibility sentence, or it would be overwritten.
    start = INDEX_HTML.index("function updateAnalysisEligibility() {")
    end = INDEX_HTML.index("\n    }\n", start)
    body = INDEX_HTML[start:end]
    assert body.rindex("review complete.`;") < body.rindex('$("analysisEligibility").textContent = state.analysisMessage;')
    assert "The analysis was interrupted because the local DEGORA server stopped." in INDEX_HTML
    assert "The local DEGORA server did not answer." in INDEX_HTML
    # A blank search leaves a notice that outlives the native validation bubble.
    assert "Nothing was searched: enter a condition, perturbation, disease or pathway in English" in INDEX_HTML


def test_a_file_where_the_output_folder_belongs_is_named_as_such(tmp_path: Path) -> None:
    from degora.provenance import ensure_output_directory, output_directory_lock

    blocker = tmp_path / "results"
    blocker.write_text("not a folder", encoding="utf-8")
    with pytest.raises(NotADirectoryError, match="is a file, not a folder"):
        ensure_output_directory(blocker)
    with pytest.raises(NotADirectoryError, match="is a file, not a folder"):
        with output_directory_lock(blocker):
            pass
    nested = blocker / "deeper"
    with pytest.raises(NotADirectoryError, match="a file stands where a folder is needed"):
        ensure_output_directory(nested)


def test_a_finished_run_leaves_no_lock_file_behind(tmp_path: Path) -> None:
    from degora.provenance import output_directory_lock

    output = tmp_path / "out"
    with output_directory_lock(output):
        assert (output / ".degora-run.lock").exists()
    assert not (output / ".degora-run.lock").exists()
    # ...and the folder is still usable for the next run.
    with output_directory_lock(output):
        pass
    assert not (output / ".degora-run.lock").exists()
