from __future__ import annotations

import hashlib
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

import degora.discovery_run as discovery_run
from degora.discovery import DiscoveryError
from degora.discovery_run import _author_pipeline, run_discovery_analysis


def test_author_pipeline_recognizes_cuffdiff_gene_exp_diff_filename() -> None:
    assert _author_pipeline("GSE70544_Normoxia_vs_Hypoxia_Gencode_gene_exp.diff.txt.gz") == "Cuffdiff"


def _prepared_bundle(root: Path) -> dict:
    studies = []
    for index, accession in enumerate(("GSE100001", "GSE100002"), start=1):
        path = root / f"{accession}_DESeq2_results.csv"
        pd.DataFrame(
            {
                "gene": ["TP53", "CDKN1A", "VEGFA"],
                "log2FoldChange": [2.0 + index / 10, 1.3, -1.1],
                "pvalue": [0.001, 0.01, 0.03],
                "padj": [0.003, 0.02, 0.04],
            }
        ).to_csv(path, index=False)
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
                        "name": path.name,
                        "source_url": f"https://ftp.ncbi.nlm.nih.gov/{path.name}",
                        "role": "deg_table",
                        "inspection": {
                            "status": "ready_for_review",
                            "fetch_scope": "full",
                            "local_path": str(path),
                            "full_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
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


def _selections() -> list[dict]:
    return [
        {
            "candidate_id": "candidate1",
            "mode": "author",
            "contrast_label": "hypoxia versus normoxia",
            "direction_confirmed": True,
            "table_scope": "full_results",
            "n_ctrl": 3,
            "n_treat": 3,
            "cell_system": "HK-2 renal epithelial cells",
            "duration_h": "24",
            "platform": "NovaSeq 6000",
        },
        {
            "candidate_id": "candidate2",
            "mode": "author",
            "contrast_label": "hypoxia versus normoxia",
            "direction_confirmed": True,
            "table_scope": "full_results",
            "n_ctrl": 4,
            "n_treat": 4,
            "cell_system": "HPTEC",
            "duration_h": "12",
            "platform": "HiSeq 2500",
        },
    ]


def _canonical_no_geo_bundle(root: Path) -> dict:
    prepared = _prepared_bundle(root)
    first, second = prepared["studies"]
    first.pop("accession", None)
    first["pubmed_ids"] = []
    first["doi"] = "https://doi.org/10.1101/2026.01.02.123456"
    first["pmcid"] = "PMC1234567"
    first["source_url"] = ""
    second.pop("accession", None)
    second["pubmed_ids"] = []
    second["canonical_id"] = "AUTHOR-PAPER-002"
    second["provider"] = "author_supplement"
    second["source_url"] = "https://example.org/paper-002"
    return prepared


def test_reviewed_human_candidates_run_end_to_end_without_cross_species_pooling(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    result = run_discovery_analysis(
        _prepared_bundle(bundle),
        _selections(),
        tmp_path / "analysis",
        species="human",
    )
    assert result["status"] == "complete"
    assert result["species"]["key"] == "human"
    assert result["cross_species_pooling"] is False
    assert result["n_source_units"] == 2
    assert Path(result["db_path"]).is_file()
    workbook = result["excel_workbook"]
    assert Path(workbook["output"]) == Path(result["results_dir"]) / "DEGORA_output.xlsx"
    assert Path(workbook["output"]).is_file()
    assert Path(workbook["manifest"]).is_file()
    assert Path(workbook["validation"]).is_file()
    assert Path(str(workbook["output"]) + ".source").is_file()
    assert (Path(result["output_dir"]) / ".degora-discovery-run.json").is_file()
    catalog = pd.read_csv(result["catalog_path"])
    assert set(catalog["species"]) == {"Homo sapiens"}
    assert set(catalog["source_input_type"]) == {"author_deg_table"}
    assert set(catalog["include_in_analysis"]) == {"yes"}
    assert catalog["n_ctrl"].tolist() == [3, 4]
    assert catalog["n_treat"].tolist() == [3, 4]
    assert catalog["cell_system"].tolist() == ["HK-2 renal epithelial cells", "HPTEC"]
    assert catalog["duration_h"].tolist() == [24, 12]
    assert catalog["platform"].tolist() == ["NovaSeq 6000", "HiSeq 2500"]
    with sqlite3.connect(result["db_path"]) as connection:
        assert connection.execute("SELECT COUNT(*) FROM genes").fetchone()[0] == 3
        meta = dict(connection.execute("SELECT key, value FROM meta").fetchall())
        assert meta["discovery_species"] == "human"
        assert meta["discovery_cross_species_pooling"] == "false"
        assert meta["discovery_source_units"] == "PMID:900001,PMID:900002"
        assert connection.execute("SELECT COUNT(DISTINCT source_unit_id) FROM studies").fetchone()[0] == 2


def test_no_geo_author_deg_papers_run_end_to_end_as_two_source_units(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    result = run_discovery_analysis(
        _canonical_no_geo_bundle(bundle),
        _selections(),
        tmp_path / "analysis",
        species="human",
    )

    assert result["status"] == "complete"
    assert result["n_source_units"] == 2
    assert result["source_units"] == ["AUTHOR_SUPPLEMENT:AUTHOR-PAPER-002", "DOI:10.1101/2026.01.02.123456"]
    catalog = pd.read_csv(result["catalog_path"])
    assert catalog["study_id"].str.contains("GSE").sum() == 0
    assert catalog["source_url"].tolist() == [
        "https://ftp.ncbi.nlm.nih.gov/GSE100001_DESeq2_results.csv",
        "https://ftp.ncbi.nlm.nih.gov/GSE100002_DESeq2_results.csv",
    ]
    assert "DOI:10.1101/2026.01.02.123456" in catalog.loc[0, "notes"]
    assert "PMCID:PMC1234567" in catalog.loc[0, "notes"]


def test_shared_no_geo_doi_collapses_to_one_source_unit(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _canonical_no_geo_bundle(bundle)
    for study in prepared["studies"]:
        study.pop("canonical_id", None)
        study.pop("provider", None)
        study["doi"] = "DOI:10.1000/SHARED"

    with pytest.raises(DiscoveryError, match="found 1"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


def test_unstable_no_geo_study_without_identifier_is_rejected(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _canonical_no_geo_bundle(bundle)
    for key in ("doi", "dois", "publication_doi", "canonical_id", "provider_accession", "accession"):
        prepared["studies"][0].pop(key, None)

    with pytest.raises(DiscoveryError, match="stable publication or data identifier"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


def test_mixed_blocked_study_activation_is_rejected(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["mixed_status"] = "mixed_blocked"

    with pytest.raises(DiscoveryError, match="mixed_blocked"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


def test_mixed_rescued_activation_requires_verified_target_species_evidence(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["mixed_status"] = "mixed_rescued"
    prepared["studies"][0]["target_species_verified"] = True

    with pytest.raises(DiscoveryError, match="target_species_verified=true and nonempty evidence"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "missing-evidence", species="human")

    prepared["studies"][0]["target_species_verified"] = False
    prepared["studies"][0]["target_species_evidence"] = "Author table explicitly labels all samples as Homo sapiens."
    with pytest.raises(DiscoveryError, match="target_species_verified=true and nonempty evidence"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "not-verified", species="human")

    prepared["studies"][0]["target_species_verified"] = True
    result = run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")
    assert result["status"] == "complete"
    assert result["n_source_units"] == 2


def test_analysis_requires_explicit_direction_confirmation(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    selections = _selections()
    selections[0]["direction_confirmed"] = False
    with pytest.raises(DiscoveryError, match="direction_confirmed"):
        run_discovery_analysis(prepared, selections, tmp_path / "analysis", species="human")


@pytest.mark.parametrize("field", ["n_ctrl", "n_treat"])
@pytest.mark.parametrize("value", [None, "", 0, -1, True, 3.0, 3.5, "3.0", "three"])
def test_author_candidate_requires_positive_whole_number_group_sizes(
    tmp_path: Path,
    field: str,
    value: object,
) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    selections = _selections()
    selections[0][field] = value
    with pytest.raises(DiscoveryError, match=rf"{field} must be a positive whole number"):
        run_discovery_analysis(_prepared_bundle(bundle), selections, tmp_path / "analysis", species="human")


@pytest.mark.parametrize(
    ("status", "message"),
    [
        ("candidate_header", "pass header review"),
        ("requires_lfc_confirmation", "column_mapping_confirmed"),
        ("requires_pvalue_mapping", "column_mapping_confirmed"),
    ],
)
def test_author_candidate_must_pass_full_value_review(tmp_path: Path, status: str, message: str) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["files"][0]["inspection"]["status"] = status
    with pytest.raises(DiscoveryError, match=message):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


def _multi_cohort_author_bundle(root: Path) -> dict:
    mixed_path = root / "paper_one_supplement.xlsx"
    mixed = pd.DataFrame(
        {
            "dataset": ["COHORT_A", "COHORT_A", "COHORT_A", "COHORT_B", "COHORT_B"],
            "GENE_NAME": ["A", "A", "B", "A", "C"],
            "gene.log2FC": [1.2, 2.5, -0.8, 1.5, -1.1],
            "gene.padj": [0.01, 0.001, 0.04, 0.005, 0.03],
        }
    )
    with pd.ExcelWriter(mixed_path) as writer:
        pd.DataFrame({"note": ["all genes and controls"]}).to_excel(writer, sheet_name="all_data", index=False)
        pd.DataFrame({"note": ["author note row"]}).to_excel(
            writer,
            sheet_name="protein_coding",
            index=False,
            header=False,
        )
        mixed.to_excel(writer, sheet_name="protein_coding", index=False, startrow=1)

    second_path = root / "paper_two_results.csv"
    pd.DataFrame(
        {
            "gene": ["A", "B", "C"],
            "log2FoldChange": [1.0, -0.7, -1.0],
            "pvalue": [0.01, 0.03, 0.02],
            "padj": [0.02, 0.04, 0.03],
        }
    ).to_csv(second_path, index=False)

    return {
        "species": {"key": "human", "scientific_name": "Homo sapiens"},
        "materialize_dir": str(root),
        "studies": [
            {
                "species": "human",
                "scientific_name": "Homo sapiens",
                "canonical_id": "pmid:111",
                "source_unit_id": "PMID:111",
                "pubmed_ids": ["111"],
                "study_type": "publication-linked public table",
                "files": [
                    {
                        "candidate_id": "mixed-author-table",
                        "name": mixed_path.name,
                        "source_url": "https://example.org/paper_one_supplement.xlsx",
                        "role": "unknown_table",
                        "inspection": {
                            "status": "requires_pvalue_mapping",
                            "fetch_scope": "full",
                            "local_path": str(mixed_path),
                            "full_file_sha256": hashlib.sha256(mixed_path.read_bytes()).hexdigest(),
                            "header_row": 2,
                            "sheet_name": "all_data",
                            "mapping": {
                                "gene_column": "GENE_NAME",
                                "lfc_column": "gene.log2FC",
                                "p_column": "",
                                "padj_column": "gene.padj",
                            },
                        },
                    }
                ],
            },
            {
                "species": "human",
                "scientific_name": "Homo sapiens",
                "canonical_id": "pmid:222",
                "source_unit_id": "PMID:222",
                "pubmed_ids": ["222"],
                "study_type": "Expression profiling by high throughput sequencing",
                "files": [
                    {
                        "candidate_id": "second-author-table",
                        "name": second_path.name,
                        "source_url": "https://example.org/paper_two_results.csv",
                        "role": "deg_table",
                        "inspection": {
                            "status": "ready_for_review",
                            "fetch_scope": "full",
                            "local_path": str(second_path),
                            "full_file_sha256": hashlib.sha256(second_path.read_bytes()).hexdigest(),
                            "header_row": 1,
                            "sheet_name": "",
                            "mapping": {
                                "gene_column": "gene",
                                "lfc_column": "log2FoldChange",
                                "p_column": "pvalue",
                                "padj_column": "padj",
                            },
                        },
                    }
                ],
            },
        ],
    }


def _multi_cohort_selections() -> list[dict]:
    common = {
        "candidate_id": "mixed-author-table",
        "mode": "author",
        "direction_confirmed": True,
        "table_scope": "full_results",
        "n_ctrl": 3,
        "n_treat": 3,
        "assay_type": "RNA-seq",
        "pipeline": "DESeq2",
        "sheet_name": "protein_coding",
        "gene_column": "GENE_NAME",
        "lfc_column": "gene.log2FC",
        "p_column": "gene.padj",
        "padj_column": "gene.padj",
        "column_mapping_confirmed": True,
        "adjusted_p_as_pvalue_confirmed": True,
        "row_filter_column": "dataset",
        "row_filter_confirmed": True,
    }
    return [
        {**common, "contrast_label": "case versus control, cohort A", "row_filter_value": "COHORT_A"},
        {**common, "contrast_label": "case versus control, cohort B", "row_filter_value": "COHORT_B"},
        {
            "candidate_id": "second-author-table",
            "mode": "author",
            "contrast_label": "case versus control, second paper",
            "direction_confirmed": True,
            "table_scope": "full_results",
            "n_ctrl": 4,
            "n_treat": 4,
        },
    ]


def test_author_activation_supports_explicit_sheet_mapping_and_repeated_filtered_cohorts(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    result = run_discovery_analysis(
        _multi_cohort_author_bundle(bundle),
        _multi_cohort_selections(),
        tmp_path / "analysis",
        species="human",
    )

    assert result["n_active_contrasts"] == 3
    assert result["n_source_units"] == 2
    catalog = pd.read_csv(result["catalog_path"])
    assert catalog["source_unit_id"].tolist() == ["PMID:111", "PMID:111", "PMID:222"]
    assert catalog["study_id"].is_unique
    assert catalog["assay_type"].tolist()[:2] == ["RNA-seq", "RNA-seq"]
    assert catalog["pipeline"].tolist()[:2] == ["DESeq2", "DESeq2"]
    assert catalog["p_column"].tolist()[:2] == ["pvalue", "pvalue"]
    assert catalog["padj_column"].tolist()[:2] == ["padj", "padj"]
    assert catalog["notes"].str.contains("adjusted-p/FDR column").tolist()[:2] == [True, True]
    assert catalog["notes"].str.contains("Exact source-row subset").tolist()[:2] == [True, True]

    derivations = result["author_derivations"]
    assert [item["row_filter_value"] for item in derivations] == ["COHORT_A", "COHORT_B", ""]
    assert [item["n_rows_after_filter"] for item in derivations[:2]] == [3, 2]
    first = pd.read_csv(catalog.loc[0, "source_path"])
    second = pd.read_csv(catalog.loc[1, "source_path"])
    assert first["gene_symbol"].tolist() == ["A", "A", "B"]
    assert second["gene_symbol"].tolist() == ["A", "C"]
    assert Path(str(catalog.loc[0, "source_path"]) + ".source").is_file()
    assert Path(str(catalog.loc[0, "source_path"]) + ".provenance.json").is_file()


def test_author_activation_can_explicitly_keep_first_duplicate_gene_row(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    selections = _multi_cohort_selections()
    selections[0]["duplicate_gene_policy"] = "keep_first"
    selections[0]["duplicate_gene_policy_confirmed"] = True

    result = run_discovery_analysis(
        _multi_cohort_author_bundle(bundle),
        selections,
        tmp_path / "analysis",
        species="human",
    )

    catalog = pd.read_csv(result["catalog_path"])
    first = pd.read_csv(catalog.loc[0, "source_path"])
    assert first["gene_symbol"].tolist() == ["A", "B"]
    assert first.loc[first["gene_symbol"].eq("A"), "log2FoldChange"].item() == pytest.approx(1.2)
    assert first.loc[first["gene_symbol"].eq("A"), "pvalue"].item() == pytest.approx(0.01)
    derivation = result["author_derivations"][0]
    assert derivation["duplicate_gene_policy"] == "keep_first"
    assert derivation["n_usable_rows_before_duplicate_policy"] == 3
    assert derivation["n_duplicate_gene_rows"] == 2
    assert derivation["n_duplicate_genes"] == 1
    assert derivation["n_usable_output_rows"] == 2
    assert "legacy/manual extraction rule" in catalog.loc[0, "notes"]


def test_author_activation_can_explicitly_clear_detected_optional_padj_column(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    selections = _multi_cohort_selections()
    selections[2]["padj_column"] = ""
    selections[2]["column_mapping_confirmed"] = True

    result = run_discovery_analysis(
        _multi_cohort_author_bundle(bundle),
        selections,
        tmp_path / "analysis",
        species="human",
    )

    catalog = pd.read_csv(result["catalog_path"]).fillna("")
    assert catalog.loc[2, "p_column"] == "pvalue"
    assert catalog.loc[2, "padj_column"] == ""
    third = pd.read_csv(catalog.loc[2, "source_path"])
    assert third.columns.tolist() == ["gene_symbol", "log2FoldChange", "pvalue"]
    assert result["author_derivations"][2]["source_mapping"]["padj_column"] == ""


def test_author_adjusted_p_mapping_and_row_filter_require_explicit_confirmation(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _multi_cohort_author_bundle(bundle)
    selections = _multi_cohort_selections()
    selections[0]["adjusted_p_as_pvalue_confirmed"] = False
    with pytest.raises(DiscoveryError, match="adjusted_p_as_pvalue_confirmed"):
        run_discovery_analysis(prepared, selections, tmp_path / "adjusted-p", species="human")

    selections = _multi_cohort_selections()
    selections[0]["row_filter_confirmed"] = False
    with pytest.raises(DiscoveryError, match="row_filter_confirmed"):
        run_discovery_analysis(prepared, selections, tmp_path / "filter", species="human")

    selections = _multi_cohort_selections()
    selections[0]["duplicate_gene_policy"] = "keep_first"
    with pytest.raises(DiscoveryError, match="duplicate_gene_policy_confirmed"):
        run_discovery_analysis(prepared, selections, tmp_path / "duplicates", species="human")

    selections = _multi_cohort_selections()
    selections[2]["padj_column"] = ""
    with pytest.raises(DiscoveryError, match="column_mapping_confirmed"):
        run_discovery_analysis(prepared, selections, tmp_path / "clear-padj", species="human")


def test_identical_author_candidate_extraction_cannot_be_activated_twice(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    selections = _multi_cohort_selections()
    selections[1]["row_filter_value"] = "COHORT_A"
    selections[1]["duplicate_gene_policy"] = "keep_first"
    selections[1]["duplicate_gene_policy_confirmed"] = True
    with pytest.raises(DiscoveryError, match="same author candidate extraction"):
        run_discovery_analysis(
            _multi_cohort_author_bundle(bundle),
            selections,
            tmp_path / "analysis",
            species="human",
        )


def test_author_candidate_rejects_header_preview_and_lfc_inversion(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["files"][0]["inspection"]["fetch_scope"] = "header_prefix"
    with pytest.raises(DiscoveryError, match="header preview"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "preview", species="human")

    prepared = _prepared_bundle(bundle)
    selections = _selections()
    selections[0]["invert_lfc"] = True
    with pytest.raises(DiscoveryError, match="inversion is not supported"):
        run_discovery_analysis(prepared, selections, tmp_path / "inversion", species="human")


def test_analysis_rejects_species_mismatch(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    with pytest.raises(DiscoveryError, match="species does not match"):
        run_discovery_analysis(_prepared_bundle(bundle), _selections(), tmp_path / "analysis", species="mouse")


def test_same_pubmed_paper_cannot_satisfy_two_source_unit_minimum(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    for study in prepared["studies"]:
        study["pubmed_ids"] = ["123456"]
    with pytest.raises(DiscoveryError, match="found 1"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


def test_analysis_rejects_file_changed_after_preparation(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    Path(prepared["studies"][0]["files"][0]["inspection"]["local_path"]).write_text(
        "gene,log2FoldChange,pvalue\nTP53,-99,0.001\n",
        encoding="utf-8",
    )
    with pytest.raises(DiscoveryError, match="changed after preparation"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


def _normalized_prepared_bundle(root: Path) -> dict:
    studies = []
    for index, accession in enumerate(("GSE200001", "GSE200002"), start=1):
        path = root / f"{accession}_matrix.tsv"
        pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "ctrl_1": [2.0, 5.0, 1.0],
                "ctrl_2": [2.2, 5.2, 1.1],
                "treat_1": [6.0 + index, 5.1, 2.0],
                "treat_2": [6.2 + index, 5.3, 2.1],
            }
        ).to_csv(path, sep="\t", index=False)
        studies.append(
            {
                "species": "human",
                "scientific_name": "Homo sapiens",
                "accession": accession,
                "study_type": "Expression profiling by array",
                "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "files": [
                    {
                        "candidate_id": f"matrix{index}",
                        "name": path.name,
                        "source_url": f"https://ftp.ncbi.nlm.nih.gov/{path.name}",
                        "role": "normalized_expression_matrix",
                        "inspection": {
                            "status": "upstream_matrix_ready_for_contrast",
                            "fetch_scope": "full",
                            "local_path": str(path),
                            "full_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                            "gene_column": "gene",
                            "sample_columns": ["ctrl_1", "ctrl_2", "treat_1", "treat_2"],
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


def _normalized_selections(*, scale: str = "log2") -> list[dict]:
    return [
        {
            "candidate_id": f"matrix{index}",
            "mode": "fallback",
            "contrast_label": "treated versus control",
            "direction_confirmed": True,
            "biological_replicates_confirmed": True,
            "gene_column": "gene",
            "normalized_scale": scale,
            "control_samples": ["ctrl_1", "ctrl_2"],
            "treatment_samples": ["treat_1", "treat_2"],
        }
        for index in (1, 2)
    ]


def test_normalized_fallback_requires_confirmed_scale_and_runs_as_log2_effect(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _normalized_prepared_bundle(bundle)
    missing_scale = _normalized_selections(scale="")
    with pytest.raises(DiscoveryError, match="normalized_scale"):
        run_discovery_analysis(prepared, missing_scale, tmp_path / "missing-scale", species="human")

    result = run_discovery_analysis(
        prepared,
        _normalized_selections(scale="log2"),
        tmp_path / "analysis",
        species="human",
    )
    catalog = pd.read_csv(result["catalog_path"])
    assert set(catalog["normalization"]) == {"public_normalized_matrix_confirmed_log2_scale"}
    assert set(catalog["probe_collapse"]) == {"median_expression"}
    assert all(summary["normalized_scale"] == "log2" for summary in result["fallback_derivations"])
    assert Path(result["analysis_request"]).is_file()
    assert (Path(result["output_dir"]) / "prepared_bundle.json").is_file()


def test_fallback_requires_independent_biological_replicate_attestation(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    selections = _normalized_selections()
    selections[0]["biological_replicates_confirmed"] = False
    output = tmp_path / "analysis"
    with pytest.raises(DiscoveryError, match="biological_replicates_confirmed"):
        run_discovery_analysis(_normalized_prepared_bundle(bundle), selections, output, species="human")
    assert not output.exists()


@pytest.mark.parametrize(
    ("control", "treatment", "message"),
    [
        (["ctrl_1"], ["treat_1", "treat_2"], "at least two"),
        (["ctrl_1", "ctrl_2"], ["ctrl_2", "treat_2"], "disjoint"),
        (["ctrl_1", "missing"], ["treat_1", "treat_2"], "were not found"),
    ],
)
def test_fallback_activation_rejects_invalid_sample_assignments(
    tmp_path: Path,
    control: list[str],
    treatment: list[str],
    message: str,
) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _normalized_prepared_bundle(bundle)
    selections = _normalized_selections()
    selections[0]["control_samples"] = control
    selections[0]["treatment_samples"] = treatment
    with pytest.raises(ValueError, match=message):
        run_discovery_analysis(prepared, selections, tmp_path / "analysis", species="human")


def test_explicit_prepared_source_unit_wins_over_recomputed_pmids(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    for study in prepared["studies"]:
        study["source_unit_id"] = "PMID:SHARED_PREPARED_UNIT"
    output = tmp_path / "analysis"
    with pytest.raises(DiscoveryError, match="found 1"):
        run_discovery_analysis(prepared, _selections(), output, species="human")
    assert not output.exists()


def test_inconsistent_prepared_units_cannot_split_one_shared_pubmed_source(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["pubmed_ids"] = ["111", "999"]
    prepared["studies"][0]["source_unit_id"] = "PMID:111"
    prepared["studies"][1]["pubmed_ids"] = ["222", "999"]
    prepared["studies"][1]["source_unit_id"] = "PMID:222"
    output = tmp_path / "analysis"
    with pytest.raises(DiscoveryError, match="sharing a PubMed ID"):
        run_discovery_analysis(prepared, _selections(), output, species="human")
    assert not output.exists()


def test_failed_late_fallback_activation_leaves_no_partial_run(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _normalized_prepared_bundle(bundle)
    for study in prepared["studies"]:
        study["source_unit_id"] = "PMID:ONE_UNIT"
    output = tmp_path / "analysis"
    with pytest.raises(DiscoveryError, match="found 1"):
        run_discovery_analysis(prepared, _normalized_selections(), output, species="human")
    assert not output.exists()


def test_excel_export_failure_rolls_back_partial_discovery_run(tmp_path: Path, monkeypatch) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    output = tmp_path / "analysis"

    def fail_export(*args, **kwargs):
        raise RuntimeError("synthetic workbook failure")

    monkeypatch.setattr(discovery_run, "export_run_workbook", fail_export)
    with pytest.raises(RuntimeError, match="synthetic workbook failure"):
        run_discovery_analysis(_prepared_bundle(bundle), _selections(), output, species="human")

    assert not output.exists()


def test_failed_forced_replacement_restores_previous_complete_run(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    output = tmp_path / "analysis"
    first = run_discovery_analysis(prepared, _selections(), output, species="human")
    original_request = Path(first["analysis_request"]).read_bytes()
    invalid = _selections()
    invalid[0]["direction_confirmed"] = False

    with pytest.raises(DiscoveryError, match="direction_confirmed"):
        run_discovery_analysis(prepared, invalid, output, species="human", force=True)

    assert Path(first["db_path"]).is_file()
    assert Path(first["analysis_request"]).read_bytes() == original_request


def test_force_refuses_loose_markers_and_preserves_unrelated_directory(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    output = tmp_path / "not-a-degora-run"
    (output / "results").mkdir(parents=True)
    precious = output / "unrelated_keep.txt"
    precious.write_text("do not delete", encoding="utf-8")

    with pytest.raises(DiscoveryError, match="not a recognized DEGORA discovery run"):
        run_discovery_analysis(_prepared_bundle(bundle), _selections(), output, species="human", force=True)

    assert precious.read_text(encoding="utf-8") == "do not delete"


def test_replay_command_records_force_only_when_requested(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    result = run_discovery_analysis(_prepared_bundle(bundle), _selections(), tmp_path / "analysis", species="human")
    source_sidecar = Path(str(result["score_csv"]) + ".source")
    assert source_sidecar.is_file()
    assert "--force" not in source_sidecar.read_text(encoding="utf-8")


def _disjoint_gene_bundle(root: Path) -> dict:
    """Two valid source units whose gene identifier spaces do not overlap."""

    prepared = _prepared_bundle(root)
    genes = (["TP53", "CDKN1A", "VEGFA"], ["ENSG00000141510", "ENSG00000124762", "ENSG00000112715"])
    for study, symbols in zip(prepared["studies"], genes, strict=True):
        candidate = study["files"][0]
        path = Path(candidate["inspection"]["local_path"])
        pd.DataFrame(
            {
                "gene": symbols,
                "log2FoldChange": [2.0, 1.3, -1.1],
                "pvalue": [0.001, 0.01, 0.03],
                "padj": [0.003, 0.02, 0.04],
            }
        ).to_csv(path, index=False)
        candidate["inspection"]["full_file_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return prepared


def test_discovery_analysis_refuses_a_run_that_scored_no_genes(tmp_path: Path) -> None:
    """A run that scores nothing is a failure, not a run with an empty table.

    Two source units that share no gene identifier produce a valid catalog and a
    valid score database with zero rows. `degora run` refuses that; the discovery
    path used to register it as status "complete" with an empty top_genes list
    and a workbook holding nothing, which reads as a finished analysis.
    """

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _disjoint_gene_bundle(bundle)
    output = tmp_path / "run"

    with pytest.raises(DiscoveryError) as excinfo:
        run_discovery_analysis(prepared, _selections(), output, species="human", min_studies=2)

    assert "scored zero genes" in str(excinfo.value)
    # The enclosing transaction must not leave the partial run behind.
    assert not output.exists() or not any(output.iterdir())
