"""Concurrency regressions for public multi-artifact and directory publishers."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import threading
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

import pandas as pd
import pytest

from degora import (
    cli,
    discovery,
    discovery_export,
    discovery_prepare,
    discovery_run,
    excel_export,
    provenance,
    score_db,
    slice_runner,
)


def _thread_call(call: Callable[[], Any]) -> tuple[threading.Thread, list[Any], list[BaseException]]:
    results: list[Any] = []
    failures: list[BaseException] = []

    def invoke() -> None:
        try:
            results.append(call())
        except BaseException as exc:  # noqa: BLE001 - surface worker failures in the test thread.
            failures.append(exc)

    worker = threading.Thread(target=invoke)
    worker.start()
    return worker, results, failures


def _discovery_generation(label: str) -> dict[str, Any]:
    return {
        "query": label,
        "species": {"key": "human", "label": "Human", "scientific_name": "Homo sapiens"},
        "studies": [
            {
                "species": "human",
                "scientific_name": "Homo sapiens",
                "accession": f"GSE_{label}",
                "source_unit_id": f"unit_{label}",
                "title": label,
                "files": [
                    {
                        "source_url": f"https://example.test/{label}.csv",
                        "name": f"{label}.csv",
                        "tier": "strong",
                        "inspection": {
                            "status": "ready_for_review",
                            "mapping": {"gene_column": "gene", "lfc_column": "lfc", "p_column": "p"},
                        },
                    }
                ],
            }
        ],
    }


@pytest.mark.parametrize(
    ("export_name", "committed_name"),
    [
        ("export_search_page", "geo_search_page.json"),
        ("export_discovery_bundle", "discovery_audit.json"),
    ],
)
def test_legacy_exporters_reject_a_racing_generation_before_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    export_name: str,
    committed_name: str,
) -> None:
    output = tmp_path / export_name
    publishing = threading.Event()
    release = threading.Event()
    real_publish = discovery.publish_staged_artifacts

    def paused_publish(artifacts: dict[Path, Path]) -> None:
        publishing.set()
        assert release.wait(timeout=5)
        real_publish(artifacts)

    monkeypatch.setattr(discovery, "publish_staged_artifacts", paused_publish)
    exporter = getattr(discovery, export_name)
    worker, results, failures = _thread_call(
        lambda: exporter(_discovery_generation("generation-A"), output, force=True)
    )
    try:
        assert publishing.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            exporter(_discovery_generation("generation-B"), output, force=True)
    finally:
        release.set()
        worker.join(timeout=5)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(results) == 1
    committed = json.loads((output / committed_name).read_text(encoding="utf-8"))
    assert committed["query"] == "generation-A"
    assert not list(output.glob(".degora-*-export-*"))


def _recognized_bundle(path: Path, label: str) -> None:
    path.mkdir(parents=True)
    (path / discovery.DISCOVERY_BUNDLE_MARKER).write_text(
        json.dumps(
            {
                "artifact_type": discovery.DISCOVERY_BUNDLE_ARTIFACT_TYPE,
                "format_version": discovery.DISCOVERY_BUNDLE_FORMAT_VERSION,
                "species": "human",
            }
        ),
        encoding="utf-8",
    )
    (path / "generation.txt").write_text(label, encoding="utf-8")


def test_preparation_lock_survives_the_old_target_being_moved_to_backup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "prepared"
    _recognized_bundle(target, "OLD")
    old_moved = threading.Event()
    release = threading.Event()
    real_replace = Path.replace

    def replace_then_pause(path: Path, destination: str | Path) -> Path:
        replaced = real_replace(path, destination)
        destination_path = Path(destination)
        if path == target and destination_path.name.startswith(f".{target.name}.backup-"):
            old_moved.set()
            assert release.wait(timeout=5)
        return replaced

    def fake_prepare(*_args: Any, materialize_dir: str | Path | None, query: str, **_kwargs: Any) -> dict[str, Any]:
        assert materialize_dir is not None
        Path(materialize_dir, "generation.txt").write_text(query, encoding="utf-8")
        return _discovery_generation(query)

    monkeypatch.setattr(Path, "replace", replace_then_pause)
    monkeypatch.setattr(discovery, "_prepare_geo_studies_in_place", fake_prepare)
    worker, results, failures = _thread_call(
        lambda: discovery.prepare_geo_studies(
            ["GSE1"], "human", query="generation-A", materialize_dir=target, force=True
        )
    )
    try:
        assert old_moved.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery.prepare_geo_studies(
                ["GSE2"], "human", query="generation-B", materialize_dir=target, force=True
            )
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery.export_discovery_bundle(_discovery_generation("generation-B"), target, force=True)
    finally:
        release.set()
        worker.join(timeout=5)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(results) == 1
    assert (target / "generation.txt").read_text(encoding="utf-8") == "generation-A"
    assert not list(tmp_path.glob(f".{target.name}.backup-*"))


def test_publication_preparation_uses_the_same_stable_target_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "publication-prepared"
    _recognized_bundle(target, "OLD")
    old_moved = threading.Event()
    release = threading.Event()
    real_replace = Path.replace

    def replace_then_pause(path: Path, destination: str | Path) -> Path:
        replaced = real_replace(path, destination)
        destination_path = Path(destination)
        if path == target and destination_path.name.startswith(f".{target.name}.backup-"):
            old_moved.set()
            assert release.wait(timeout=5)
        return replaced

    def fake_prepare(*_args: Any, staging: Path, query: str, **_kwargs: Any) -> dict[str, Any]:
        (staging / "generation.txt").write_text(query, encoding="utf-8")
        return _discovery_generation(query)

    monkeypatch.setattr(Path, "replace", replace_then_pause)
    monkeypatch.setattr(discovery_prepare, "_prepare_into_staging", fake_prepare)
    worker, results, failures = _thread_call(
        lambda: discovery_prepare.prepare_publication_records(
            [], "human", query="generation-A", materialize_dir=target, force=True
        )
    )
    try:
        assert old_moved.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery_prepare.prepare_publication_records(
                [], "human", query="generation-B", materialize_dir=target, force=True
            )
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery_export.export_publication_search({"records": []}, target, force=True)
    finally:
        release.set()
        worker.join(timeout=5)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(results) == 1
    assert (target / "generation.txt").read_text(encoding="utf-8") == "generation-A"
    assert not list(tmp_path.glob(f".{target.name}.backup-*"))


def test_distinct_bundle_and_run_targets_in_one_api_parent_can_overlap(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bundle_entered = threading.Event()
    bundle_release = threading.Event()

    def fake_prepare(
        _accessions: Any,
        _species: str,
        *,
        materialize_dir: str | Path | None,
        query: str,
        **_kwargs: Any,
    ) -> dict[str, Any]:
        if query == "A":
            bundle_entered.set()
            assert bundle_release.wait(timeout=5)
        return {"query": query, "materialize_dir": str(materialize_dir)}

    monkeypatch.setattr(discovery, "_prepare_geo_studies_locked", fake_prepare)
    bundle_parent = tmp_path / "bundles"
    bundle_worker, bundle_results, bundle_failures = _thread_call(
        lambda: discovery.prepare_geo_studies(
            ["GSE1"], "human", query="A", materialize_dir=bundle_parent / "bundle-a"
        )
    )
    try:
        assert bundle_entered.wait(timeout=5)
        # This call must finish while A still owns its sibling target lock.
        bundle_b = discovery.prepare_geo_studies(
            ["GSE2"], "human", query="B", materialize_dir=bundle_parent / "bundle-b"
        )
    finally:
        bundle_release.set()
        bundle_worker.join(timeout=5)
    assert bundle_b["query"] == "B"
    assert not bundle_worker.is_alive()
    assert len(bundle_results) == 1
    assert not bundle_failures, bundle_failures

    run_entered = threading.Event()
    run_release = threading.Event()

    def fake_run(
        _prepared: dict[str, Any],
        _selections: Any,
        output_dir: str | Path,
        **_kwargs: Any,
    ) -> dict[str, Any]:
        if Path(output_dir).name == "run-a":
            run_entered.set()
            assert run_release.wait(timeout=5)
        return {"output_dir": str(output_dir)}

    monkeypatch.setattr(discovery_run, "_run_discovery_analysis_locked", fake_run)
    run_parent = tmp_path / "runs"
    run_worker, run_results, run_failures = _thread_call(
        lambda: discovery_run.run_discovery_analysis(
            {}, [], run_parent / "run-a", species="human"
        )
    )
    try:
        assert run_entered.wait(timeout=5)
        run_b = discovery_run.run_discovery_analysis(
            {}, [], run_parent / "run-b", species="human"
        )
    finally:
        run_release.set()
        run_worker.join(timeout=5)
    assert Path(run_b["output_dir"]).name == "run-b"
    assert not run_worker.is_alive()
    assert len(run_results) == 1
    assert not run_failures, run_failures


@pytest.mark.parametrize(
    ("canonical_name", "alias_name"),
    [("Prepared", "prepared"), ("Caf\u00e9", "Cafe\u0301")],
)
def test_filesystem_equivalent_target_aliases_share_one_lock_identity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    canonical_name: str,
    alias_name: str,
) -> None:
    canonical_target = tmp_path / canonical_name
    canonical_target.mkdir()
    alias_target = tmp_path / alias_name
    try:
        equivalent = alias_target.exists() and alias_target.samefile(canonical_target)
    except OSError:
        equivalent = False
    if not equivalent:
        pytest.skip("filesystem keeps these target spellings distinct")

    assert provenance.publication_target_lock_path(
        canonical_target
    ) == provenance.publication_target_lock_path(alias_target)
    entered = threading.Event()
    release = threading.Event()

    def fake_prepare(
        _accessions: Any,
        _species: str,
        *,
        materialize_dir: str | Path | None,
        query: str,
        **_kwargs: Any,
    ) -> dict[str, Any]:
        if query == "A":
            entered.set()
            assert release.wait(timeout=5)
        return {"query": query, "materialize_dir": str(materialize_dir)}

    monkeypatch.setattr(discovery, "_prepare_geo_studies_locked", fake_prepare)
    worker, results, failures = _thread_call(
        lambda: discovery.prepare_geo_studies(
            ["GSE1"], "human", query="A", materialize_dir=canonical_target, force=True
        )
    )
    try:
        assert entered.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery.prepare_geo_studies(
                ["GSE2"], "human", query="B", materialize_dir=alias_target, force=True
            )
    finally:
        release.set()
        worker.join(timeout=5)

    assert not worker.is_alive()
    assert len(results) == 1
    assert not failures, failures

    # An outer owner using one spelling must be re-entrant when its nested
    # publisher receives another spelling for the same physical target.
    lock_path = provenance.publication_target_lock_path(canonical_target)
    with provenance.output_directory_lock(lock_path):
        nested = discovery.prepare_geo_studies(
            ["GSE3"], "human", query="nested", materialize_dir=alias_target, force=True
        )
    assert nested["query"] == "nested"


def _recognized_run(path: Path, label: str) -> None:
    path.mkdir(parents=True)
    (path / discovery_run.DISCOVERY_RUN_MARKER).write_text(
        json.dumps(
            {
                "artifact_type": discovery_run.DISCOVERY_RUN_ARTIFACT_TYPE,
                "format_version": discovery_run.DISCOVERY_RUN_FORMAT_VERSION,
            }
        ),
        encoding="utf-8",
    )
    (path / "generation.txt").write_text(label, encoding="utf-8")


def test_analysis_lock_prevents_a_successful_generation_from_being_lost(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "analysis"
    _recognized_run(output, "OLD")
    old_moved = threading.Event()
    release = threading.Event()
    real_replace = Path.replace

    def replace_then_pause(path: Path, destination: str | Path) -> Path:
        replaced = real_replace(path, destination)
        destination_path = Path(destination)
        if path == output and destination_path.name.startswith(f".{output.name}.backup-"):
            old_moved.set()
            assert release.wait(timeout=5)
        return replaced

    def fake_execute(
        _prepared: dict[str, Any],
        _selections: Any,
        output_dir: str | Path,
        **_kwargs: Any,
    ) -> dict[str, Any]:
        target = Path(output_dir)
        (target / "generation.txt").write_text("generation-A", encoding="utf-8")
        (target / discovery_run.DISCOVERY_RUN_MARKER).write_text(
            json.dumps(
                {
                    "artifact_type": discovery_run.DISCOVERY_RUN_ARTIFACT_TYPE,
                    "format_version": discovery_run.DISCOVERY_RUN_FORMAT_VERSION,
                }
            ),
            encoding="utf-8",
        )
        return {"status": "complete", "output_dir": str(target)}

    monkeypatch.setattr(Path, "replace", replace_then_pause)
    monkeypatch.setattr(discovery_run, "_execute_discovery_analysis", fake_execute)
    worker, results, failures = _thread_call(
        lambda: discovery_run.run_discovery_analysis({}, [], output, species="human", force=True)
    )
    try:
        assert old_moved.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery_run.run_discovery_analysis({}, [], output, species="human", force=True)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            excel_export.export_run_workbook(output / "results", command="racing direct export")
    finally:
        release.set()
        worker.join(timeout=5)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(results) == 1
    assert (output / "generation.txt").read_text(encoding="utf-8") == "generation-A"
    assert not list(tmp_path.glob(f".{output.name}.backup-*"))


def _minimal_result_dir(path: Path) -> Path:
    path.mkdir()
    with sqlite3.connect(path / "degora_scores.db") as connection:
        pd.DataFrame(
            {
                "gene_symbol": ["ISG15"],
                "quality_weighted_degora_rank": [1],
                "quality_weighted_degora_score": [1.0],
            }
        ).to_sql("genes", connection, index=False)
        pd.DataFrame({"gene_symbol": ["ISG15"], "study_id": ["S1"]}).to_sql(
            "gene_evidence", connection, index=False
        )
        pd.DataFrame({"study_id": ["S1"], "source_unit_id": ["U1"]}).to_sql(
            "studies", connection, index=False
        )
        pd.DataFrame({"key": ["degora_version"], "value": ["test"]}).to_sql(
            "meta", connection, index=False
        )
    (path / "degora_score_metadata.json").write_text("{}\n", encoding="utf-8")
    return path


def _score_generation(label: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": [f"{label}-S1", f"{label}-S2"],
            "paper_id": [f"{label}-P1", f"{label}-P2"],
            "gene_symbol": [f"GENE_{label}", f"GENE_{label}"],
            "lfc": [2.0, 1.8],
            "signed_z": [5.0, 4.5],
            "pvalue": [1e-6, 2e-6],
            "padj": [1e-5, 2e-5],
            "normalized_rank": [0.01, 0.02],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
            "n_genes_in_study": [1000, 1000],
            "pipeline": ["DESeq2", "edgeR"],
            "assay_type": ["RNA-seq", "RNA-seq"],
            "source_input_type": ["author_deg_table", "author_deg_table"],
            "table_scope": ["full_results", "full_results"],
            "species": ["Homo sapiens", "Homo sapiens"],
            "source_path": [f"{label}-source-1.csv", f"{label}-source-2.csv"],
            "source_url": [f"https://example.test/{label}/1", f"https://example.test/{label}/2"],
        }
    )


def _prepared_analysis_bundle(root: Path) -> dict[str, Any]:
    studies: list[dict[str, Any]] = []
    for index, accession in enumerate(("GSE100001", "GSE100002"), start=1):
        source = root / f"{accession}_DESeq2_results.csv"
        pd.DataFrame(
            {
                "gene": ["TP53", "CDKN1A", "VEGFA"],
                "log2FoldChange": [2.0 + index / 10, 1.3, -1.1],
                "pvalue": [0.001, 0.01, 0.03],
                "padj": [0.003, 0.02, 0.04],
            }
        ).to_csv(source, index=False)
        studies.append(
            {
                "species": "human",
                "scientific_name": "Homo sapiens",
                "accession": accession,
                "pubmed_ids": [str(900000 + index)],
                "study_type": "Expression profiling by high throughput sequencing",
                "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "files": [
                    {
                        "candidate_id": f"candidate{index}",
                        "name": source.name,
                        "source_url": f"https://ftp.ncbi.nlm.nih.gov/{source.name}",
                        "role": "deg_table",
                        "inspection": {
                            "status": "ready_for_review",
                            "fetch_scope": "full",
                            "local_path": str(source),
                            "full_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                            "mapping": {
                                "gene_column": "gene",
                                "lfc_column": "log2FoldChange",
                                "p_column": "pvalue",
                                "padj_column": "padj",
                            },
                        },
                    }
                ],
            }
        )
    return {
        "species": {"key": "human", "scientific_name": "Homo sapiens"},
        "materialize_dir": str(root),
        "studies": studies,
    }


def _analysis_selections() -> list[dict[str, Any]]:
    return [
        {
            "candidate_id": f"candidate{index}",
            "mode": "author",
            "contrast_label": "hypoxia versus normoxia",
            "direction_confirmed": True,
            "table_scope": "full_results",
            "n_ctrl": 3,
            "n_treat": 3,
        }
        for index in (1, 2)
    ]


def _slice_catalog(root: Path, label: str) -> Path:
    rows: list[dict[str, Any]] = []
    for index in (1, 2):
        source = root / f"{label}-source-{index}.csv"
        pd.DataFrame(
            {
                "gene_symbol": [f"GENE_{label}", f"ONLY_{label}"],
                "log2FoldChange": [2.0 + index / 10, 1.0],
                "pvalue": [1e-6 * index, 0.01 + index / 1000],
                "padj": [2e-6 * index, 0.02 + index / 1000],
            }
        ).to_csv(source, index=False)
        rows.append(
            {
                "study_id": f"{label}-S{index}",
                "paper_id": f"{label}-P{index}",
                "source_unit_id": f"{label}-U{index}",
                "source_path": str(source),
                "pipeline": "DESeq2",
                "species": "Homo sapiens",
                "gene_column": "gene_symbol",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
                "assay_type": "RNA-seq",
                "source_input_type": "author_deg_table",
                "table_scope": "full_results",
                "include_in_analysis": "yes",
                "n_ctrl": 3,
                "n_treat": 3,
            }
        )
    catalog = root / f"{label}-catalog.csv"
    pd.DataFrame(rows, columns=slice_runner.CATALOG_COLUMNS).to_csv(catalog, index=False)
    return catalog


def test_direct_score_writer_cannot_cross_discovery_commit_barrier(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    output = tmp_path / "analysis"
    harmonized_b = tmp_path / "harmonized-b.csv"
    _score_generation("B").to_csv(harmonized_b, index=False)
    at_commit_barrier = threading.Event()
    release = threading.Event()

    def pause_before_publish() -> None:
        at_commit_barrier.set()
        assert release.wait(timeout=10)

    worker, results, failures = _thread_call(
        lambda: discovery_run.run_discovery_analysis(
            _prepared_analysis_bundle(bundle),
            _analysis_selections(),
            output,
            species="human",
            excel=False,
            before_publish=pause_before_publish,
        )
    )
    try:
        assert at_commit_barrier.wait(timeout=10)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            score_db.write_score_database(
                harmonized_b,
                output / "results",
                min_studies=2,
                extra_metadata={"generation": "B"},
            )
    finally:
        release.set()
        worker.join(timeout=20)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(results) == 1
    assert results[0]["status"] == "complete"
    with sqlite3.connect(results[0]["db_path"]) as connection:
        final_genes = {row[0] for row in connection.execute("SELECT gene_symbol FROM genes")}
    assert "GENE_B" not in final_genes


def test_direct_slice_cannot_cross_discovery_commit_barrier(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    output = tmp_path / "analysis"
    catalog_b = _slice_catalog(tmp_path, "B")
    at_commit_barrier = threading.Event()
    release = threading.Event()

    def pause_before_publish() -> None:
        at_commit_barrier.set()
        assert release.wait(timeout=10)

    worker, results, failures = _thread_call(
        lambda: discovery_run.run_discovery_analysis(
            _prepared_analysis_bundle(bundle),
            _analysis_selections(),
            output,
            species="human",
            excel=False,
            before_publish=pause_before_publish,
        )
    )
    try:
        assert at_commit_barrier.wait(timeout=10)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            slice_runner.run_slice(catalog_b, output / "results", output / "harmonized", 2)
    finally:
        release.set()
        worker.join(timeout=20)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(results) == 1
    final_genes = set(pd.read_csv(output / "results" / "slice_consensus.csv")["gene_symbol"])
    assert "GENE_B" not in final_genes
    standalone_results = tmp_path / "standalone" / "results"
    standalone_metrics = slice_runner.run_slice(
        catalog_b,
        standalone_results,
        tmp_path / "standalone" / "harmonized",
        2,
    )
    assert standalone_metrics["n_consensus_genes"] == 2
    assert set(pd.read_csv(standalone_results / "slice_consensus.csv")["gene_symbol"]) == {
        "GENE_B",
        "ONLY_B",
    }


def test_cli_holds_canonical_results_parent_lock_for_the_whole_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "analysis"
    results_dir = run_dir / "results"
    pipeline_entered = threading.Event()
    release = threading.Event()

    monkeypatch.setattr(cli, "read_excel_settings", lambda _path: {})
    monkeypatch.setattr(cli, "_validate_score_version", lambda _settings: None)
    monkeypatch.setattr(cli, "_validate_for_run", lambda *_args, **_kwargs: {})

    def paused_pipeline(*_args: Any, **_kwargs: Any) -> int:
        pipeline_entered.set()
        assert release.wait(timeout=10)
        return 0

    monkeypatch.setattr(cli, "_run_pipeline", paused_pipeline)
    arguments = SimpleNamespace(
        config=str(tmp_path / "config.xlsx"),
        min_studies=None,
        output_dir=str(results_dir),
        harmonized_dir=str(run_dir / "harmonized"),
        db=None,
        quiet=True,
    )
    worker, results, failures = _thread_call(lambda: cli._run_from_config(arguments))
    try:
        assert pipeline_entered.wait(timeout=10)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            discovery_run.run_discovery_analysis({}, [], run_dir, species="human")
    finally:
        release.set()
        worker.join(timeout=10)

    assert not worker.is_alive()
    assert not failures, failures
    assert results == [0]


def test_distinct_score_outputs_cannot_publish_different_generations_to_one_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    harmonized_a = tmp_path / "harmonized-a.csv"
    harmonized_b = tmp_path / "harmonized-b.csv"
    _score_generation("A").to_csv(harmonized_a, index=False)
    _score_generation("B").to_csv(harmonized_b, index=False)
    output_a = tmp_path / "output-a"
    output_b = tmp_path / "output-b"
    shared_db = tmp_path / "shared" / "scores.db"
    entered = threading.Event()
    release = threading.Event()
    real_write = score_db._write_score_database_locked

    def paused_write(
        harmonized_path: Path,
        output_dir: Path,
        **kwargs: Any,
    ) -> dict[str, Any]:
        if Path(output_dir) == output_a.resolve():
            entered.set()
            assert release.wait(timeout=5)
        return real_write(harmonized_path, output_dir, **kwargs)

    monkeypatch.setattr(score_db, "_write_score_database_locked", paused_write)
    worker, results, failures = _thread_call(
        lambda: score_db.write_score_database(
            harmonized_a,
            output_a,
            db_path=shared_db,
            min_studies=2,
            extra_metadata={"generation": "A"},
        )
    )
    try:
        assert entered.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            score_db.write_score_database(
                harmonized_b,
                output_b,
                db_path=shared_db,
                min_studies=2,
                extra_metadata={"generation": "B"},
            )
    finally:
        release.set()
        worker.join(timeout=10)

    assert not worker.is_alive()
    assert len(results) == 1
    assert not failures, failures
    csv_genes = set(pd.read_csv(output_a / "degora_gene_scores.csv")["gene_symbol"])
    with sqlite3.connect(shared_db) as connection:
        db_genes = {row[0] for row in connection.execute("SELECT gene_symbol FROM genes")}
        db_metadata = dict(connection.execute("SELECT key, value FROM meta"))
    assert csv_genes == db_genes == {"GENE_A"}
    assert db_metadata["generation"] == "A"
    assert not (output_b / "degora_gene_scores.csv").exists()


def test_workbook_export_is_exclusive_and_reentrant_for_the_outer_cli_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    results = _minimal_result_dir(tmp_path / "results")
    output = results / "DEGORA_output.xlsx"
    publishing = threading.Event()
    release = threading.Event()
    real_publish = excel_export.publish_staged_artifacts

    def paused_publish(artifacts: dict[Path, Path]) -> None:
        publishing.set()
        assert release.wait(timeout=5)
        real_publish(artifacts)

    monkeypatch.setattr(excel_export, "publish_staged_artifacts", paused_publish)
    worker, exported, failures = _thread_call(
        lambda: excel_export.export_run_workbook(results, output, command="generation-A")
    )
    try:
        assert publishing.wait(timeout=5)
        with pytest.raises(provenance.OutputDirectoryBusyError):
            excel_export.export_run_workbook(results, output, command="generation-B")
    finally:
        release.set()
        worker.join(timeout=5)

    assert not worker.is_alive()
    assert not failures, failures
    assert len(exported) == 1
    excel_export.verify_run_workbook_export(output)

    # ``degora run`` already owns this exact directory while it calls the
    # workbook exporter.  Same-thread nesting must remain re-entrant.
    monkeypatch.setattr(excel_export, "publish_staged_artifacts", real_publish)
    with provenance.output_directory_lock(results):
        excel_export.export_run_workbook(results, output, command="outer CLI")
    excel_export.verify_run_workbook_export(output)
