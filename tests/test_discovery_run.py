from __future__ import annotations

import gzip
import hashlib
import json
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

import degora.discovery_run as discovery_run
from degora.discovery import DiscoveryError
from degora.discovery_run import _author_pipeline, run_discovery_analysis
from degora.slice_runner import read_catalog


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


def test_reviewed_author_candidate_reads_gzipped_workbook(tmp_path: Path) -> None:
    import gzip
    import shutil

    from openpyxl import Workbook

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    plain = bundle / "GSE100001_DESeq2_results.xlsx"
    book = Workbook()
    sheet = book.active
    sheet.title = "DEG"
    sheet.append(["Author-supplied differential-expression results"])
    sheet.append(["gene", "log2FoldChange", "pvalue", "padj"])
    sheet.append(["TP53", 2.1, 0.001, 0.003])
    sheet.append(["CDKN1A", 1.3, 0.01, 0.02])
    sheet.append(["VEGFA", -1.1, 0.03, 0.04])
    book.save(plain)
    compressed = bundle / "GSE100001_DESeq2_results.xlsx.gz"
    with plain.open("rb") as source, gzip.open(compressed, "wb") as target:
        shutil.copyfileobj(source, target)

    first_candidate = prepared["studies"][0]["files"][0]
    first_candidate["name"] = compressed.name
    first_candidate["source_url"] = f"https://ftp.ncbi.nlm.nih.gov/{compressed.name}"
    first_candidate["inspection"]["local_path"] = str(compressed)
    first_candidate["inspection"]["full_file_sha256"] = hashlib.sha256(compressed.read_bytes()).hexdigest()
    first_candidate["inspection"]["sheet_name"] = "DEG"
    first_candidate["inspection"]["header_row"] = 2

    selections = _selections()
    selections[0]["sheet_name"] = "DEG"
    result = run_discovery_analysis(prepared, selections, tmp_path / "analysis", species="human")

    assert result["status"] == "complete"
    assert result["author_derivations"][0]["sheet_name"] == "DEG"
    assert result["author_derivations"][0]["header_row"] == 2
    catalog = pd.read_csv(result["catalog_path"])
    materialized = pd.read_csv(catalog.loc[0, "source_path"])
    assert materialized["gene_symbol"].tolist() == ["TP53", "CDKN1A", "VEGFA"]


def test_formula_guarded_author_csv_preserves_raw_gene_in_scoring_and_sqlite(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    for study in prepared["studies"]:
        candidate = study["files"][0]
        path = Path(candidate["inspection"]["local_path"])
        frame = pd.read_csv(path)
        frame.loc[0, "gene"] = "=BAD()"
        frame.to_csv(path, index=False)
        candidate["inspection"]["full_file_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()

    selections = _selections()
    selections[0]["contrast_label"] = "-Dox versus +Dox"
    result = run_discovery_analysis(
        prepared,
        selections,
        tmp_path / "analysis",
        species="human",
    )

    catalog = pd.read_csv(result["catalog_path"])
    assert catalog.loc[0, "hypoxia_modality"] == "'-Dox versus +Dox"
    assert read_catalog(result["catalog_path"]).loc[0, "hypoxia_modality"] == "-Dox versus +Dox"
    catalog_provenance = json.loads(Path(result["catalog_path"] + ".provenance.json").read_text())
    assert catalog_provenance["metadata"]["csv_formula_guard"] == "reversible_apostrophe_prefix_v1"
    published_author_table = pd.read_csv(catalog.loc[0, "source_path"])
    assert published_author_table.loc[0, "gene_symbol"] == "'=BAD()"
    score_csv = pd.read_csv(Path(result["results_dir"]) / "degora_gene_scores.csv")
    assert "'=BAD()" in set(score_csv["gene_symbol"])
    with sqlite3.connect(result["db_path"]) as connection:
        genes = {row[0] for row in connection.execute("SELECT gene_symbol FROM genes")}
    assert "=BAD()" in genes
    assert "'=BAD()" not in genes


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


def test_author_group_sizes_cannot_exceed_known_study_sample_total(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["n_samples"] = 5

    with pytest.raises(DiscoveryError, match=r"n_ctrl\+n_treat=6 exceeds.*n_samples=5.*evidence weight"):
        run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")


@pytest.mark.parametrize("study_sample_total", [6, "unknown"])
def test_author_group_sizes_allow_exact_or_unknown_study_sample_total(
    tmp_path: Path,
    study_sample_total: object,
) -> None:
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    prepared["studies"][0]["n_samples"] = study_sample_total
    prepared["studies"][1]["n_samples"] = 8 if study_sample_total == 6 else "unknown"

    result = run_discovery_analysis(
        prepared,
        _selections(),
        tmp_path / "analysis",
        species="human",
    )

    assert result["status"] == "complete"


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


@pytest.mark.parametrize(
    ("left_identity", "right_identity", "message"),
    [
        ({"accession": "GSE123"}, {"accession": "E-GEOD-123"}, "public study accession"),
        ({"accession": "E-MTAB-456"}, {"accession": "EMTAB456"}, "public study accession"),
        (
            {"pubmed_ids": "PMID:123; PMID:456"},
            {"pubmed_ids": ["PMID-456"]},
            "PubMed ID",
        ),
        ({"pmcid": "PMCID:PMC789"}, {"pmcids": ["PMC-789"]}, "PMC ID"),
        (
            {"doi": "10.1000/Same.DOI"},
            {"publication_doi": "https://doi.org/10.1000/same.doi"},
            "DOI",
        ),
    ],
)
def test_prepared_unit_validation_rejects_public_identifier_aliases_split_across_units(
    left_identity: dict[str, object],
    right_identity: dict[str, object],
    message: str,
) -> None:
    studies = [
        {"source_unit_id": "OPAQUE-UNIT-A", **left_identity},
        {"source_unit_id": "OPAQUE-UNIT-B", **right_identity},
    ]

    with pytest.raises(DiscoveryError, match=message):
        discovery_run._validate_prepared_source_units(studies)


def test_prepared_unit_validation_keeps_truly_distinct_public_units() -> None:
    studies = [
        {"source_unit_id": "OPAQUE-UNIT-A", "accession": "GSE123", "pmid": "PMID:100"},
        {"source_unit_id": "OPAQUE-UNIT-B", "accession": "E-GEOD-124", "pmid": "PMID:101"},
        {"source_unit_id": "OPAQUE-UNIT-C", "accession": "EMTAB456", "pmcid": "PMC200"},
        {"source_unit_id": "OPAQUE-UNIT-D", "accession": "E-MTAB-457", "pmcid": "PMCID:PMC201"},
    ]

    discovery_run._validate_prepared_source_units(studies)


@pytest.mark.parametrize(
    ("study", "expected"),
    [
        ({"source_unit_id": "E-GEOD-123"}, "GSE123"),
        ({"source_unit_id": "EMTAB456"}, "E-MTAB-456"),
        ({"source_unit_id": "PMID-789"}, "PMID:789"),
        ({"source_unit_id": "PMCID:PMC321"}, "PMCID:PMC321"),
        ({"source_unit_id": "PMC-741"}, "PMCID:PMC741"),
        ({"source_unit_id": "12345"}, "PMID:12345"),
        ({"pubmed_ids": "PMID:654; PMID-987"}, "PMID:654"),
    ],
)
def test_paper_source_unit_uses_canonical_public_identifier(study: dict[str, object], expected: str) -> None:
    assert discovery_run._paper_source_unit(study) == expected


def test_bare_and_prefixed_pmid_source_units_collapse_before_counting() -> None:
    studies = [
        {"source_unit_id": "12345"},
        {"source_unit_id": "PMID:12345"},
    ]

    discovery_run._validate_prepared_source_units(studies)

    assert {discovery_run._paper_source_unit(study) for study in studies} == {"PMID:12345"}


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


def test_a_generated_source_unit_id_never_carries_the_list_delimiter() -> None:
    """A DOI may legitimately contain a semicolon; a source_unit_id may not.

    The preflight rejects a hand-written identifier holding the delimiter, but the
    discovery path builds identifiers itself, so it has to be unable to produce one.
    """

    from degora.discovery_run import _paper_source_unit

    assert _paper_source_unit({"doi": "10.1234/abc;def"}) == "DOI:10.1234/abc_def"
    assert _paper_source_unit({"source_unit_id": "U;1"}) == "U_1"
    assert _paper_source_unit({"pubmed_ids": ["12345"]}) == "PMID:12345"


def _matrix_candidate(tmp_path: Path, values: dict[str, list[float]]) -> tuple[dict, dict, Path]:
    """A prepared upstream-matrix candidate whose file sits inside the bundle root."""

    bundle = tmp_path / "bundle"
    bundle.mkdir(exist_ok=True)
    path = bundle / "GSE100002_matrix.csv"
    pd.DataFrame({"gene": [f"G{i}" for i in range(len(next(iter(values.values()))))], **values}).to_csv(path, index=False)
    candidate = {
        "candidate_id": "matrix1",
        "name": path.name,
        "role": "unknown_matrix",
        "inspection": {
            "status": "upstream_matrix_ready_for_contrast",
            "fetch_scope": "full",
            "local_path": str(path),
            "full_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sample_columns": list(values),
            "gene_column": "gene",
            "header_row": 1,
        },
    }
    study = {"accession": "GSE100002", "title": "matrix study", "files": [candidate]}
    return study, candidate, bundle


def _fallback_entry(**overrides) -> dict:
    entry = {
        "candidate_id": "matrix1", "mode": "fallback", "direction_confirmed": True,
        "biological_replicates_confirmed": True, "control_samples": ["c1", "c2"],
        "treatment_samples": ["t1", "t2"], "matrix_type": "normalized_expression_matrix",
        "normalized_scale": "log2", "gene_column": "gene",
    }
    entry.update(overrides)
    return entry


def test_an_unrecognised_matrix_type_is_a_usage_error(tmp_path: Path) -> None:
    """A wrong matrix_type ended in a bare ValueError from the derivation."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, {"c1": [1.0] * 30, "c2": [1.1] * 30, "t1": [2.0] * 30, "t2": [2.2] * 30})

    with pytest.raises(DiscoveryError, match="matrix_type='normalized_expression' is not recognised"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type="normalized_expression"),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )


def test_a_fractional_matrix_is_refused_as_raw_counts(tmp_path: Path) -> None:
    """FPKM selected as count_matrix would be handed to a count model as counts."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(
        tmp_path, {"c1": [1.37] * 30, "c2": [2.51] * 30, "t1": [4.09] * 30, "t2": [3.73] * 30}
    )

    with pytest.raises(DiscoveryError, match="whole numbers"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type="count_matrix"),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )


def test_fallback_rejects_two_columns_from_the_same_geo_sample(tmp_path: Path) -> None:
    """Count and FPKM columns for one GSM are not two biological replicates."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_count": [10.0] * 30,
        "c1_FPKM": [1.0] * 30,
        "t1_count": [20.0] * 30,
        "t1_FPKM": [2.0] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    candidate["inspection"]["sample_labels"] = {
        "c1_count": {"accession": "GSM1001"},
        "c1_FPKM": {"accession": "GSM1001"},
        "t1_count": {"accession": "GSM1002"},
        "t1_FPKM": {"accession": "GSM1002"},
    }
    entry = _fallback_entry(
        control_samples=["c1_count", "c1_FPKM"],
        treatment_samples=["t1_count", "t1_FPKM"],
    )

    with pytest.raises(DiscoveryError, match="same GEO biological sample accession"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_fallback_rejects_paired_measurement_families_without_geo_labels(tmp_path: Path) -> None:
    """Family metadata protects an ambiguous older bundle with no GSM mapping."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_count": [10.0] * 30,
        "c2_count": [12.0] * 30,
        "t1_count": [20.0] * 30,
        "t2_count": [22.0] * 30,
        "c1_FPKM": [1.0] * 30,
        "c2_FPKM": [1.2] * 30,
        "t1_FPKM": [2.0] * 30,
        "t2_FPKM": [2.2] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    candidate["inspection"]["sample_column_families"] = {
        "count": ["c1_count", "c2_count", "t1_count", "t2_count"],
        "normalized": ["c1_FPKM", "c2_FPKM", "t1_FPKM", "t2_FPKM"],
        "paired_base_samples": ["c1", "c2", "t1", "t2"],
        "selected_family": "",
        "selection_basis": "",
    }
    entry = _fallback_entry(
        control_samples=["c1_count", "c1_FPKM"],
        treatment_samples=["t1_count", "t2_count"],
    )

    with pytest.raises(DiscoveryError, match="same biological sample base"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_fallback_rejects_mixed_measurement_families_across_different_samples(tmp_path: Path) -> None:
    """Count controls versus FPKM treatments would compare scales, not biology."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_count": [10.0] * 30,
        "c2_count": [12.0] * 30,
        "t1_FPKM": [2.0] * 30,
        "t2_FPKM": [2.2] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    candidate["inspection"]["sample_column_families"] = {
        "count": ["c1_count", "c2_count"],
        "normalized": ["t1_FPKM", "t2_FPKM"],
        "paired_base_samples": [],
        "selected_family": "",
        "selection_basis": "",
    }
    entry = _fallback_entry(
        control_samples=["c1_count", "c2_count"],
        treatment_samples=["t1_FPKM", "t2_FPKM"],
    )

    with pytest.raises(DiscoveryError, match="mix count and normalized measurement families"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_fallback_rejects_explicit_and_unclassified_measurement_mix(tmp_path: Path) -> None:
    """An unlabeled column cannot silently form a contrast with known FPKM columns."""

    from degora.discovery import inspect_upstream_bytes, normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_FPKM": [1.0] * 30,
        "c2_FPKM": [1.2] * 30,
        "t1_FPKM": [2.0] * 30,
        "t2_FPKM": [2.2] * 30,
        "x1": [100.0] * 30,
        "x2": [120.0] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    source_path = Path(candidate["inspection"]["local_path"])
    inspected = inspect_upstream_bytes(
        source_path.name,
        source_path.read_bytes(),
        declared_role="normalized_expression_matrix",
    )
    candidate["role"] = "normalized_expression_matrix"
    candidate["inspection"] = {**candidate["inspection"], **inspected}
    entry = _fallback_entry(
        control_samples=["c1_FPKM", "c2_FPKM"],
        treatment_samples=["x1", "x2"],
    )

    with pytest.raises(DiscoveryError, match="explicitly labeled measurement columns with unclassified columns"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_fallback_rederives_mixed_families_when_legacy_bundle_has_no_metadata(tmp_path: Path) -> None:
    """Clear suffixes remain fail-closed when an old prepared card lacks family metadata."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_count": [10.0] * 30,
        "c2_count": [12.0] * 30,
        "t1_FPKM": [2.0] * 30,
        "t2_FPKM": [2.2] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    entry = _fallback_entry(
        control_samples=["c1_count", "c2_count"],
        treatment_samples=["t1_FPKM", "t2_FPKM"],
    )

    with pytest.raises(DiscoveryError, match="mix count and normalized measurement families"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


@pytest.mark.parametrize("matrix_type", ["count_matrix", "estimated_count_matrix"])
def test_fallback_rejects_integer_normalized_suffix_as_a_count_role_without_metadata(
    tmp_path: Path, matrix_type: str
) -> None:
    """Integer-looking FPKM values must not pass either count path."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_FPKM": [10.0] * 30,
        "c2_FPKM": [12.0] * 30,
        "t1_FPKM": [20.0] * 30,
        "t2_FPKM": [22.0] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    entry = _fallback_entry(
        control_samples=["c1_FPKM", "c2_FPKM"],
        treatment_samples=["t1_FPKM", "t2_FPKM"],
        matrix_type=matrix_type,
        contrast_label="FPKM treatment vs control",
    )

    with pytest.raises(DiscoveryError, match="normalized-suffix sample columns cannot use a count matrix role"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_fallback_rejects_count_suffix_as_a_normalized_role_without_metadata(tmp_path: Path) -> None:
    """An explicit count suffix overrides an unsafe manual normalized role."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_count": [10.0] * 30,
        "c2_count": [12.0] * 30,
        "t1_count": [20.0] * 30,
        "t2_count": [22.0] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    entry = _fallback_entry(
        control_samples=["c1_count", "c2_count"],
        treatment_samples=["t1_count", "t2_count"],
        matrix_type="normalized_expression_matrix",
        normalized_scale="linear",
        contrast_label="count treatment vs control",
    )

    with pytest.raises(DiscoveryError, match="count-suffix sample columns cannot use a normalized matrix role"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_declared_normalized_counts_filename_allows_count_suffix_with_explicit_scale(tmp_path: Path) -> None:
    """A `_count` suffix does not prove raw counts when the file declares normalization."""

    from degora.discovery import classify_filename, inspect_upstream_bytes, normalize_species
    from degora.discovery_run import _fallback_row

    rows = 40
    values = {
        "c1_count": [1.5 + index * 0.10 for index in range(rows)],
        "c2_count": [1.8 + index * 0.10 for index in range(rows)],
        "t1_count": [3.0 + index * 0.20 for index in range(rows)],
        "t2_count": [3.4 + index * 0.20 for index in range(rows)],
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    old_path = Path(candidate["inspection"]["local_path"])
    source_path = old_path.with_name("normalized_counts.csv")
    old_path.rename(source_path)
    assessment = classify_filename(source_path.name)
    assert assessment["role"] == "normalized_expression_matrix"
    inspected = inspect_upstream_bytes(
        source_path.name,
        source_path.read_bytes(),
        declared_role=assessment["role"],
    )
    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    candidate["name"] = source_path.name
    candidate["role"] = assessment["role"]
    candidate["inspection"] = {
        **candidate["inspection"],
        **inspected,
        "fetch_scope": "full",
        "local_path": str(source_path),
        "full_file_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    }
    row, _summary = _fallback_row(
        study=study,
        candidate=candidate,
        entry=_fallback_entry(
            control_samples=["c1_count", "c2_count"],
            treatment_samples=["t1_count", "t2_count"],
            matrix_type="normalized_expression_matrix",
            normalized_scale="linear",
            contrast_label="normalized-count treatment vs control",
        ),
        spec=normalize_species("human"),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora",
    )

    assert row["source_input_type"] == "normalized_expression_matrix"
    assert Path(row["source_path"]).is_file()


def test_fallback_rejects_inspected_normalized_subtype_mix_across_distinct_bases(tmp_path: Path) -> None:
    """Inspector metadata and the server agree that FPKM-versus-TPM is not biology."""

    from degora.discovery import inspect_upstream_bytes, normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    values = {
        "c1_FPKM": [1.0] * 30,
        "c2_FPKM": [1.2] * 30,
        "x1_FPKM": [1.4] * 30,
        "x2_FPKM": [1.6] * 30,
        "t1_TPM": [2.0] * 30,
        "t2_TPM": [2.2] * 30,
    }
    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    source_path = Path(candidate["inspection"]["local_path"])
    inspected = inspect_upstream_bytes(
        source_path.name,
        source_path.read_bytes(),
        declared_role="normalized_expression_matrix",
    )
    candidate["inspection"] = {**candidate["inspection"], **inspected}
    assert candidate["inspection"]["status"] == "upstream_matrix_ready_for_contrast"
    assert candidate["inspection"]["sample_column_families"]["subtypes_present"]["normalized"] == ["fpkm", "tpm"]
    entry = _fallback_entry(
        control_samples=["c1_FPKM", "c2_FPKM"],
        treatment_samples=["t1_TPM", "t2_TPM"],
    )

    with pytest.raises(DiscoveryError, match="mix normalized measurement subtypes"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


@pytest.mark.parametrize(
    ("family", "values", "control_samples", "treatment_samples", "base_by_column"),
    [
        (
            "normalized",
            {
                "c1_FPKM": [1.0] * 30,
                "c1_TPM": [1.1] * 30,
                "c2_FPKM": [1.2] * 30,
                "t1_FPKM": [2.0] * 30,
                "t2_FPKM": [2.2] * 30,
            },
            ["c1_FPKM", "c1_TPM"],
            ["t1_FPKM", "t2_FPKM"],
            {"c1_FPKM": "c1", "c1_TPM": "c1", "c2_FPKM": "c2", "t1_FPKM": "t1", "t2_FPKM": "t2"},
        ),
        (
            "count",
            {
                "c1_raw_count": [10.0] * 30,
                "c1_count": [11.0] * 30,
                "c2_count": [12.0] * 30,
                "t1_count": [20.0] * 30,
                "t2_count": [22.0] * 30,
            },
            ["c1_raw_count", "c1_count"],
            ["t1_count", "t2_count"],
            {"c1_raw_count": "c1", "c1_count": "c1", "c2_count": "c2", "t1_count": "t1", "t2_count": "t2"},
        ),
    ],
)
def test_fallback_rejects_same_family_subtypes_for_one_base_without_geo_labels(
    tmp_path: Path,
    family: str,
    values: dict[str, list[float]],
    control_samples: list[str],
    treatment_samples: list[str],
    base_by_column: dict[str, str],
) -> None:
    """FPKM+TPM or raw-count+count for one base are not extra replicates."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, values)
    candidate["inspection"]["sample_column_families"] = {
        "count": list(values) if family == "count" else [],
        "normalized": list(values) if family == "normalized" else [],
        "base_by_column": base_by_column,
        "selected_family": family,
        "selection_basis": f"declared_role={family}",
    }
    entry = _fallback_entry(
        control_samples=control_samples,
        treatment_samples=treatment_samples,
        matrix_type="count_matrix" if family == "count" else "normalized_expression_matrix",
    )

    with pytest.raises(DiscoveryError, match="same biological sample base"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_an_r_export_with_row_label_genes_materialises_on_the_author_path(tmp_path: Path) -> None:
    """The inspector said gene_column=row_name; materialisation said row_name was not found."""

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    r_export = bundle / "GSE100001_deseq2_out.csv"
    r_export.write_text(
        '"","baseMean","log2FoldChange","pvalue","padj"\n'
        + "".join(f'"{gene}",{100 + i},{2.0 - i * 0.1},0.001,0.01\n' for i, gene in enumerate(["TP53", "CDKN1A", "VEGFA"])),
        encoding="utf-8",
    )
    first = prepared["studies"][0]["files"][0]
    first["name"] = r_export.name
    first["source_url"] = f"https://ftp.ncbi.nlm.nih.gov/{r_export.name}"
    first["inspection"]["local_path"] = str(r_export)
    first["inspection"]["full_file_sha256"] = hashlib.sha256(r_export.read_bytes()).hexdigest()
    first["inspection"]["mapping"]["gene_column"] = "row_name"
    first["inspection"]["header_row"] = 1

    result = run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")

    assert result["status"] == "complete"
    catalog = pd.read_csv(result["catalog_path"])
    materialized = pd.read_csv(catalog.loc[0, "source_path"])
    assert materialized["gene_symbol"].tolist() == ["TP53", "CDKN1A", "VEGFA"]


def test_a_gzipped_workbook_matrix_reaches_the_welch_derivation(tmp_path: Path) -> None:
    """Preparation opened the .xlsx.gz as a workbook; the fallback read it as CSV and died."""

    import gzip

    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    plain = bundle / "GSE100003_matrix.xlsx"
    pd.DataFrame(
        {"gene": [f"GX{i}" for i in range(40)], "c1": [5.0 + i * 0.4 for i in range(40)], "c2": [5.5 + i * 0.4 for i in range(40)],
         "t1": [9.0 + i * 0.4 for i in range(40)], "t2": [9.5 + i * 0.4 for i in range(40)]}
    ).to_excel(plain, index=False)
    gz = bundle / "GSE100003_matrix.xlsx.gz"
    with plain.open("rb") as src, gzip.open(gz, "wb") as dst:
        dst.write(src.read())
    plain.unlink()
    candidate = {
        "candidate_id": "matrix1", "name": gz.name, "role": "unknown_matrix",
        "inspection": {"status": "upstream_matrix_ready_for_contrast", "fetch_scope": "full", "local_path": str(gz),
                       "full_file_sha256": hashlib.sha256(gz.read_bytes()).hexdigest(),
                       "sample_columns": ["c1", "c2", "t1", "t2"], "gene_column": "gene", "header_row": 1},
    }
    study = {"accession": "GSE100003", "title": "gz workbook", "files": [candidate]}

    row, summary = _fallback_row(
        study=study, candidate=candidate, entry=_fallback_entry(contrast_label="t vs c"), spec=normalize_species("human"),
        bundle_root=bundle, derived_dir=tmp_path / "derived", sequence=1, replay_command="degora",
    )

    derived = pd.read_csv(row["source_path"])
    assert len(derived) == 40, "every gene of the gzipped workbook reached the derived table"
    assert {"gene_symbol", "log2FoldChange", "pvalue", "padj"} <= set(derived.columns)



def test_the_result_warning_list_carries_each_warning_once(tmp_path: Path) -> None:
    """Validation and the run both inspect the tables; each warning arrived from both."""

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    # Make the first author table up-only with small and large values: the
    # effect-scale check reports that shape from validation and from the run,
    # and before deduplication the same sentence appeared twice in the result.
    first = prepared["studies"][0]["files"][0]
    path = Path(first["inspection"]["local_path"])
    genes = [f"GX{i}" for i in range(40)]
    pd.DataFrame({"gene": genes, "log2FoldChange": [0.2 + i * 0.1 for i in range(40)], "pvalue": [0.001] * 40, "padj": [0.01] * 40}).to_csv(path, index=False)
    first["inspection"]["full_file_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    result = run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human", min_studies=1)

    warnings = result["warnings"]
    assert any("no negative values" in w for w in warnings), warnings
    assert len(warnings) == len(set(warnings))


def _fractional_counts(rows: int = 40) -> dict[str, list[float]]:
    """Salmon-style estimated counts: non-negative, fractional, reaching the thousands."""

    return {
        "c1": [round(3.0 + i * 41.37, 2) for i in range(rows)],
        "c2": [round(1.0 + i * 39.11, 2) for i in range(rows)],
        "t1": [round(8.0 + i * 87.53, 2) for i in range(rows)],
        "t2": [round(6.0 + i * 91.27, 2) for i in range(rows)],
    }


def test_estimated_counts_take_the_count_path_once_declared(tmp_path: Path) -> None:
    """Salmon, RSEM and kallisto write fractional counts; the whole-number guard refused them all."""

    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, _fractional_counts())
    row, summary = _fallback_row(
        study=study, candidate=candidate, entry=_fallback_entry(matrix_type="estimated_count_matrix", contrast_label="drug vs ctrl"),
        spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
        sequence=1, replay_command="degora",
    )
    assert row["pipeline"] == "logCPM_Welch_derived_from_public_counts"
    assert "declared estimated counts" in row["notes"]
    assert summary["estimated_counts_note"].startswith("GSE100002: The matrix was declared estimated counts")


def test_the_count_refusal_names_the_estimated_counts_alternative(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, _fractional_counts())
    with pytest.raises(DiscoveryError, match="select matrix_type=estimated_count_matrix"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type="count_matrix"),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )


def test_estimated_counts_still_refuse_what_cannot_be_counts(tmp_path: Path) -> None:
    """The whole-number share says nothing about estimated counts; sign and magnitude still do."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    log_scale = {"c1": [1.37 + i * 0.1 for i in range(40)], "c2": [-2.51 + i * 0.1 for i in range(40)],
                 "t1": [4.09 + i * 0.1 for i in range(40)], "t2": [3.73 + i * 0.1 for i in range(40)]}
    study, candidate, bundle = _matrix_candidate(tmp_path, log_scale)
    with pytest.raises(DiscoveryError, match="never negative"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type="estimated_count_matrix", contrast_label="drug vs ctrl"),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )
    small = {key: [abs(value) for value in values] for key, values in log_scale.items()}
    study, candidate, bundle = _matrix_candidate(tmp_path, small)
    with pytest.raises(DiscoveryError, match="never exceed"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type="estimated_count_matrix", contrast_label="drug vs ctrl"),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )


def test_a_file_labeled_counts_by_the_repository_may_be_declared_estimated(tmp_path: Path) -> None:
    """A declared count_matrix role ignored matrix_type, so a fractional one had no way on."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, _fractional_counts())
    candidate["role"] = "count_matrix"
    with pytest.raises(DiscoveryError, match="whole numbers"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type=""),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )
    row, _summary = _fallback_row(
        study=study, candidate=candidate, entry=_fallback_entry(matrix_type="estimated_count_matrix", contrast_label="drug vs ctrl"),
        spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
        sequence=1, replay_command="degora",
    )
    assert "declared estimated counts" in row["notes"]
    # ...but a declared count file cannot be talked into being a normalized matrix.
    study, candidate, bundle = _matrix_candidate(tmp_path, {"c1": [10.0] * 40, "c2": [12.0] * 40, "t1": [30.0] * 40, "t2": [33.0] * 40})
    candidate["role"] = "count_matrix"
    row, _summary = _fallback_row(
        study=study, candidate=candidate, entry=_fallback_entry(matrix_type="normalized_expression_matrix", normalized_scale="log2", contrast_label="drug vs ctrl"),
        spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
        sequence=1, replay_command="degora",
    )
    assert row["pipeline"] == "logCPM_Welch_derived_from_public_counts"


def _two_arm_matrix(tmp_path: Path):
    """One file, a shared control and two treatment arms (multi-arm design)."""

    values = {"c1": [10.0] * 40, "c2": [12.0] * 40, "t1": [30.0] * 40, "t2": [33.0] * 40, "u1": [50.0] * 40, "u2": [55.0] * 40}
    return _matrix_candidate(tmp_path, values)


def test_one_matrix_file_may_carry_two_contrasts_with_different_groups(tmp_path: Path) -> None:
    """A shared control against two treatment arms was refused as "selected more than once"."""

    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    study, candidate, bundle = _two_arm_matrix(tmp_path)
    spec = normalize_species("human")
    first = _fallback_entry(matrix_type="count_matrix", contrast_label="drug A vs ctrl")
    second = _fallback_entry(matrix_type="count_matrix", contrast_label="drug B vs ctrl", treatment_samples=["u1", "u2"])
    rows = [
        _fallback_row(study=study, candidate=candidate, entry=entry, spec=spec, bundle_root=bundle,
                      derived_dir=tmp_path / "derived", sequence=index, replay_command="degora")[0]
        for index, entry in enumerate((first, second), start=1)
    ]
    assert len(rows) == 2 and rows[0] != rows[1]
    derived = sorted(path.name for path in (tmp_path / "derived").glob("*_welch.csv"))
    assert len(derived) == 2  # the two derivations do not overwrite each other


def test_identical_and_swapped_groups_from_one_file_are_still_refused() -> None:
    from degora.discovery_run import _check_fallback_selection_consistency

    prepared = {"studies": [{"accession": "GSE1", "files": [{"candidate_id": "m1"}]}]}
    same = {"candidate_id": "m1", "mode": "fallback", "control_samples": ["c1", "c2"], "treatment_samples": ["t1", "t2"], "contrast_label": "A"}
    with pytest.raises(Exception, match="selected twice"):
        _check_fallback_selection_consistency([same, {**same, "contrast_label": "B"}], prepared)
    swapped = {**same, "control_samples": ["t1", "t2"], "treatment_samples": ["c1", "c2"], "contrast_label": "B"}
    with pytest.raises(Exception):
        _check_fallback_selection_consistency([same, swapped], prepared)
    # A shared control against another arm is allowed and is not a warning: that
    # is what a multi-arm design looks like. Only a sample that changes role
    # between contrasts is reported.
    other = {**same, "treatment_samples": ["u1", "u2"], "contrast_label": "B"}
    assert _check_fallback_selection_consistency([same, other], prepared) == []
    crossed = {**same, "control_samples": ["t1", "t2"], "treatment_samples": ["u1", "u2"], "contrast_label": "C"}
    warnings = _check_fallback_selection_consistency([same, crossed], prepared)
    assert warnings and "control sample in one selected contrast and a treatment sample in another" in str(warnings[0])


def test_a_declared_group_size_has_a_plausibility_ceiling() -> None:
    from degora.discovery_run import DiscoveryError, _optional_count

    assert _optional_count("3", field="n_ctrl") == 3
    assert _optional_count("10000", field="n_ctrl") == 10000
    with pytest.raises(DiscoveryError, match="not a plausible number of biological replicates"):
        _optional_count("999999", field="n_ctrl")


def test_a_sample_column_is_refused_as_the_gene_column_up_front(tmp_path: Path) -> None:
    """The run used to end in a zero-genes diagnosis long after this could be said."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(tmp_path, {"c1": [10.0] * 40, "c2": [12.0] * 40, "t1": [30.0] * 40, "t2": [33.0] * 40})
    with pytest.raises(DiscoveryError, match="gene_column='c1' is one of the sample columns"):
        _fallback_row(
            study=study, candidate=candidate, entry=_fallback_entry(matrix_type="count_matrix", contrast_label="x", gene_column="c1"),
            spec=normalize_species("human"), bundle_root=bundle, derived_dir=tmp_path / "derived",
            sequence=1, replay_command="degora",
        )


def _matrix_candidate_for_path(
    path: Path,
    *,
    sheet_name: str = "",
    header_row: int = 1,
    gene_column: str = "gene",
) -> tuple[dict, dict, Path]:
    candidate = {
        "candidate_id": "matrix1",
        "name": path.name,
        "role": "normalized_expression_matrix",
        "inspection": {
            "status": "upstream_matrix_ready_for_contrast",
            "fetch_scope": "full",
            "local_path": str(path),
            "full_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sample_columns": ["c1", "c2", "t1", "t2"],
            "gene_column": gene_column,
            "sheet_name": sheet_name,
            "header_row": header_row,
        },
    }
    study = {"accession": "GSE100004", "title": "reader coordinates", "files": [candidate]}
    return study, candidate, path.parent


def test_workbook_preflight_and_derivation_share_the_inspected_sheet_and_header(tmp_path: Path) -> None:
    """A non-first sheet/header used to no-op preflight and then fail derivation."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    path = bundle / "GSE100004_matrix.xlsx"
    real = pd.DataFrame(
        {
            "gene": [f"REAL{i}" for i in range(40)],
            "c1": [100.0 + i for i in range(40)],
            "c2": [110.0 + i for i in range(40)],
            "t1": [500.0 + i for i in range(40)],
            "t2": [520.0 + i for i in range(40)],
        }
    )
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame({"gene": ["DECOY"], "c1": [1.0], "c2": [1.0], "t1": [2.0], "t2": [2.0]}).to_excel(
            writer, sheet_name="decoy", index=False
        )
        real.to_excel(writer, sheet_name="matrix", index=False, startrow=2)
    study, candidate, bundle = _matrix_candidate_for_path(path, sheet_name="matrix", header_row=3)

    with pytest.raises(DiscoveryError, match="normalized_scale=log2"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=_fallback_entry(contrast_label="t vs c", normalized_scale="log2"),
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )

    row, summary = _fallback_row(
        study=study,
        candidate=candidate,
        entry=_fallback_entry(contrast_label="t vs c", normalized_scale="linear"),
        spec=normalize_species("human"),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=2,
        replay_command="degora",
    )
    derived = pd.read_csv(row["source_path"])
    assert set(derived["gene_symbol"]) == set(real["gene"])
    assert summary["sheet_name"] == "matrix"
    assert summary["header_row"] == 3


def test_matrix_preflight_fails_closed_when_inspected_columns_cannot_be_reopened(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(
        tmp_path, {"c1": [1.0] * 30, "c2": [1.1] * 30, "t1": [2.0] * 30, "renamed_t2": [2.1] * 30}
    )
    candidate["inspection"]["sample_columns"] = ["c1", "c2", "t1", "t2"]
    with pytest.raises(DiscoveryError, match=r"preflight could not find.*t2"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=_fallback_entry(contrast_label="t vs c"),
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_matrix_preflight_fails_closed_when_selected_values_are_not_numeric(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    study, candidate, bundle = _matrix_candidate(
        tmp_path, {name: ["not-a-number"] * 30 for name in ("c1", "c2", "t1", "t2")}
    )
    with pytest.raises(DiscoveryError, match="preflight found no numeric values"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=_fallback_entry(contrast_label="t vs c"),
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora",
        )


def test_bom_semicolon_matrix_with_leading_blank_lines_is_read_identically(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    path = bundle / "GSE100004_matrix.csv"
    rows = ["gene;c1;c2;t1;t2"] + [f"GX{i};{i + 1}.1;{i + 1}.2;{i + 5}.1;{i + 5}.2" for i in range(40)]
    path.write_text("\ufeff\n\n" + "\n".join(rows) + "\n", encoding="utf-8")
    study, candidate, bundle = _matrix_candidate_for_path(path)
    row, summary = _fallback_row(
        study=study,
        candidate=candidate,
        entry=_fallback_entry(contrast_label="t vs c", normalized_scale="linear"),
        spec=normalize_species("human"),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora",
    )
    assert summary["n_gene_rows"] == 40
    assert set(pd.read_csv(row["source_path"])["gene_symbol"]) == {f"GX{i}" for i in range(40)}


def test_geo_series_matrix_preamble_is_excluded_by_preflight_and_derivation(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    path = bundle / "GSE100004_series_matrix.txt.gz"
    table = ['"gene"\t"c1"\t"c2"\t"t1"\t"t2"']
    table.extend(f'"GX{i}"\t{i + 1}.1\t{i + 1}.2\t{i + 5}.1\t{i + 5}.2' for i in range(40))
    payload = (
        '!Series_title\t"metadata with tabs"\n\n!series_matrix_table_begin\n'
        + "\n".join(table)
        + "\n!series_matrix_table_end\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(payload)
    study, candidate, bundle = _matrix_candidate_for_path(path)
    row, summary = _fallback_row(
        study=study,
        candidate=candidate,
        entry=_fallback_entry(contrast_label="t vs c", normalized_scale="linear"),
        spec=normalize_species("human"),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora",
    )
    assert summary["n_input_rows"] == 40
    assert set(pd.read_csv(row["source_path"])["gene_symbol"]) == {f"GX{i}" for i in range(40)}


def test_author_table_keeps_header_alignment_with_blank_lines_and_sniffed_separator(tmp_path: Path) -> None:
    """Inspection counted blank lines and sniffed semicolons; activation must do both too."""

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    prepared = _prepared_bundle(bundle)
    first = prepared["studies"][0]["files"][0]
    path = Path(first["inspection"]["local_path"])
    path.write_text(
        "\n\n"
        "gene;log2FoldChange;pvalue;padj\n"
        "TP53;2.1;0.001;0.003\n"
        "CDKN1A;1.3;0.01;0.02\n"
        "VEGFA;-1.1;0.03;0.04\n",
        encoding="utf-8",
    )
    first["inspection"]["header_row"] = 3
    first["inspection"]["full_file_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()

    result = run_discovery_analysis(prepared, _selections(), tmp_path / "analysis", species="human")
    catalog = pd.read_csv(result["catalog_path"])
    materialized = pd.read_csv(catalog.loc[catalog["study_id"].str.contains("GSE100001"), "source_path"].iloc[0])
    assert materialized["gene_symbol"].tolist() == ["TP53", "CDKN1A", "VEGFA"]
    assert "row_name" not in materialized.columns
