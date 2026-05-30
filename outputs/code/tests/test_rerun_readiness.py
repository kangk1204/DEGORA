from __future__ import annotations

from pathlib import Path

from scripts.write_rerun_readiness import build_readiness, render_markdown, write_outputs


def _write_source(path: Path) -> None:
    path.with_suffix(path.suffix + ".source").write_text("make -C outputs/code fixture\n")


def test_rerun_readiness_blocks_until_iteration_harmonized_exists(tmp_path: Path) -> None:
    project = tmp_path
    (project / "data" / "studies").mkdir(parents=True)
    catalog = project / "data" / "studies" / "hypoxia_catalog.csv"
    catalog.write_text("study_id,include_in_analysis\nHYP001,true\n")
    _write_source(catalog)

    readiness = build_readiness(project, 14)

    assert readiness["ready_to_run_full_chain"] is False
    assert "downstream reruns are not executable until slice creates data/deg/harmonized/iter-14_harmonized.csv" in readiness["blockers"]
    markdown = render_markdown(readiness)
    assert "Do not tune DEGORA scoring" in markdown
    assert "make -C outputs/code baseline CORPUS=hypoxia" in markdown


def test_rerun_readiness_writes_source_sidecars(tmp_path: Path) -> None:
    project = tmp_path
    (project / "data" / "studies").mkdir(parents=True)
    (project / "data" / "deg" / "harmonized").mkdir(parents=True)
    catalog = project / "data" / "studies" / "hypoxia_catalog.csv"
    catalog.write_text("study_id,include_in_analysis\nHYP001,true\n")
    harmonized = project / "data" / "deg" / "harmonized" / "iter-14_harmonized.csv"
    harmonized.write_text("study_id,gene_id,symbol,effect,pvalue\n")

    readiness = build_readiness(project, 14)
    json_path, md_path = write_outputs(readiness, project / "outputs" / "results" / "iter-14", "make -C outputs/code rerun-readiness ITER=14")

    assert json_path.exists()
    assert md_path.exists()
    assert json_path.with_suffix(".json.source").read_text().startswith("make -C outputs/code rerun-readiness ITER=14")
    assert md_path.with_suffix(".md.source").read_text().startswith("make -C outputs/code rerun-readiness ITER=14")
