from __future__ import annotations

import csv
import json
from pathlib import Path

from scripts.write_publication_gate_summary import _iter_number_from_name, _latest_existing, build_summary


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_source(path: Path) -> None:
    path.with_suffix(path.suffix + ".source").write_text("test fixture\n")


def test_publication_gate_summary_uses_integrated_iteration11_evidence(tmp_path: Path) -> None:
    project = tmp_path
    results = project / "outputs" / "results" / "iter-11"
    results.mkdir(parents=True)
    (project / "data" / "studies").mkdir(parents=True)
    (project / "data" / "studies" / "curation").mkdir(parents=True)

    _write_csv(
        project / "data" / "studies" / "hypoxia_catalog.csv",
        [
            {
                "study_id": f"HYP{i:03d}",
                "paper_id": f"P{i:03d}",
                "include_in_analysis": "true",
            }
            for i in range(1, 16)
        ],
    )
    _write_csv(
        results / "sensitivity_metrics.csv",
        [
            {
                "variant": "collapse_by_paper_id",
                "recall_at_50": "0.70",
                "recall_at_100": "0.70",
            }
        ],
    )
    _write_source(results / "sensitivity_metrics.csv")
    (results / "slice_metrics.json").write_text(json.dumps({"n_active_catalog_rows": 15}))
    _write_source(results / "slice_metrics.json")

    _write_csv(
        project / "data" / "studies" / "curation" / "iter11_tnfa_rescue_candidate_ledger.csv",
        [{"eligibility_status": "deferred", "paper_id": f"T{i}", "source_family_id": f"T{i}"} for i in range(40)],
    )
    (results / "tnfa_rescue_gate.json").write_text(
        json.dumps(
            {
                "passed": False,
                "scarcity_triggered": True,
                "n_active_rows": 0,
                "n_independent_units": 0,
            }
        )
    )
    _write_source(results / "tnfa_rescue_gate.json")
    (results / "tnfa_rescue_scarcity_report.md").write_text("scarcity\n")
    _write_source(results / "tnfa_rescue_scarcity_report.md")

    _write_csv(
        results / "baselines" / "baseline_parity_matrix.csv",
        [
            {"method_id": "degora_slice", "run_status": "ok"},
            {"method_id": "weighted_stouffer", "run_status": "ok"},
            {"method_id": "robustrankaggreg", "run_status": "blocked"},
        ],
    )
    _write_source(results / "baselines" / "baseline_parity_matrix.csv")

    summary = build_summary(project, 11)

    assert "Claim state: **narrow but defensible paper**" in summary
    assert "Iteration-11 hypoxia slice, diagnostic, and rank-plane artifacts are present" in summary
    assert "scarcity documented; TNF activation failed" in summary
    assert "tier-0 baselines runnable; direct-prior-art blockers remain" in summary
    assert "Lane A/B/C artifacts were not all present" not in summary


def test_latest_existing_uses_iteration_token_not_all_filename_digits(tmp_path: Path) -> None:
    iter9 = tmp_path / "outputs" / "results" / "iter-9" / "Table2_summary.csv"
    iter10 = tmp_path / "outputs" / "results" / "iter-10" / "Table1_summary.csv"
    iter9.parent.mkdir(parents=True)
    iter10.parent.mkdir(parents=True)
    iter9.write_text("old\n")
    iter10.write_text("new\n")

    assert _iter_number_from_name("Table2_summary.csv") == -1
    assert _iter_number_from_name("iter-9") == 9
    assert _latest_existing([iter9, iter10], 10) == iter10
    assert _latest_existing([iter9, iter10], 9) == iter9
