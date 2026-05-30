from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.write_direct_prior_art_no_go import write_no_go_artifacts


def test_direct_prior_art_no_go_marks_subset_candidates_without_superiority_claim(tmp_path: Path) -> None:
    source_ok = tmp_path / "hyp020.csv"
    source_blocked = tmp_path / "other.csv"
    pd.DataFrame(
        columns=["gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    ).to_csv(source_ok, index=False)
    pd.DataFrame(columns=["gene", "logFC", "pvalue"]).to_csv(source_blocked, index=False)
    harmonized = pd.DataFrame(
        {
            "study_id": ["HYP020", "HYP999"],
            "paper_id": ["P1", "P2"],
            "pipeline": ["DESeq2", "limma"],
            "source_path": [str(source_ok), str(source_blocked)],
            "source_url": ["https://example.test/1", "https://example.test/2"],
            "gene_symbol": ["VEGFA", "VEGFA"],
            "lfc": [1.0, 1.1],
            "pvalue": [0.01, 0.02],
        }
    )
    harmonized_path = tmp_path / "harmonized.csv"
    harmonized.to_csv(harmonized_path, index=False)

    summary = write_no_go_artifacts(harmonized_path, tmp_path / "out", corpus="test", command="test command")

    hstouffer = pd.read_csv(tmp_path / "out" / "hstouffer_faithful_input_audit.tsv", sep="\t")
    assert summary["claim_allowed"] is False
    assert "HYP020" in summary["hstouffer"]["compatible_subset_candidates"]
    assert hstouffer.loc[hstouffer["study_id"].eq("HYP020"), "whole_corpus_status"].iloc[0] == "no_go"
    assert (tmp_path / "out" / "direct_prior_art_no_go_report.md.source").exists()
