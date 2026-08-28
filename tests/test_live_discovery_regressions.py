from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from degora.cli import main
from degora.discovery import classify_filename, inspect_upstream_bytes
from degora.reanalysis import derive_welch_deg, read_matrix_frame, sniff_delimited_separator


@pytest.mark.parametrize(
    ("name", "text", "declared_role", "expected_gene", "expected_samples"),
    [
        (
            "GSE304706_geo_sub_raw_count.csv.gz",
            '\"\",\"A1\",\"A3\",\"A4\",\"A6\"\n'
            '\"Zc3h14\",1028,1236,1106,947\n'
            '\"Troap\",34,37,23,30\n',
            "count_matrix",
            "row_name",
            ["A1", "A3", "A4", "A6"],
        ),
        (
            "GSE317884_all_data_fpkms_human.csv.gz",
            ",NC0_1,NC0_2,NC0_3,NC4_1\n"
            "5S_rRNA,0,0,0.0549,0\n"
            "5_8S_rRNA,0,0,0,0.111\n",
            "normalized_expression_matrix",
            "row_name",
            ["NC0_1", "NC0_2", "NC0_3", "NC4_1"],
        ),
        (
            "GSE337799_rna_counts.tsv.gz",
            '\"HV0021_C_1h_rep1\" \"HV0022_P_1h_rep1\" \"HV0023_HS_1h_rep1\" \"HV0024_HSP_1h_rep1\"\n'
            '\"ENSG00000000003.14\" 1991 2916 2189 2195\n'
            '\"ENSG00000000005.5\" 0 0 0 0\n',
            "count_matrix",
            "row_name",
            ["HV0021_C_1h_rep1", "HV0022_P_1h_rep1", "HV0023_HS_1h_rep1", "HV0024_HSP_1h_rep1"],
        ),
        (
            "GSE318560_raw_counts_matrix.txt.gz",
            '\"HC_148\"\t\"HC_025\"\t\"HC_061\"\t\"HC_128\"\n'
            '\"DDX11L1\"\t6\t0\t0\t2\n'
            '\"WASH7P\"\t51\t44\t18\t72\n',
            "count_matrix",
            "row_name",
            ["HC_148", "HC_025", "HC_061", "HC_128"],
        ),
        (
            "GSE315097_gene_count_mat.txt.gz",
            "EnsemblGene_GeneSymbol\tnormoxia-rep1\tnormoxia-rep2\thypoxia-rep1\thypoxia-rep2\n"
            "ENSG00000000003_TSPAN6\t956\t1192\t1333\t1333\n"
            "ENSG00000000005_TNMD\t0\t0\t0\t0\n",
            "count_matrix",
            "EnsemblGene_GeneSymbol",
            ["normoxia-rep1", "normoxia-rep2", "hypoxia-rep1", "hypoxia-rep2"],
        ),
    ],
)
def test_upstream_inspection_recovers_public_r_row_names(
    name: str,
    text: str,
    declared_role: str,
    expected_gene: str,
    expected_samples: list[str],
) -> None:
    """Common GEO R exports keep gene labels in an unnamed first column."""

    result = inspect_upstream_bytes(name, gzip.compress(text.encode()), declared_role=declared_role)

    assert result["status"] == "upstream_matrix_ready_for_contrast"
    assert result["gene_column"] == expected_gene
    assert result["sample_columns"] == expected_samples


def test_upstream_row_name_recovery_does_not_promote_numeric_row_numbers() -> None:
    """A missing header is not enough: numeric row indices are not gene identifiers."""

    text = (
        '"S1"\t"S2"\t"S3"\t"S4"\n'
        "1\t10\t20\t30\t40\n"
        "2\t11\t21\t31\t41\n"
    )
    result = inspect_upstream_bytes(
        "numeric_row_numbers_raw_counts.tsv.gz",
        gzip.compress(text.encode()),
        declared_role="count_matrix",
    )

    assert result["status"] == "not_upstream_matrix"
    assert result["gene_column"] == ""


@pytest.mark.parametrize(
    "name",
    [
        "GSE315097_gene_count_mat.txt.gz",
        "study_gene_counts_mat.csv.gz",
    ],
)
def test_gene_count_mat_filenames_are_upstream_count_matrices(name: str) -> None:
    classified = classify_filename(name)

    assert classified["role"] == "count_matrix"
    assert classified["tier"] == "upstream"
    assert classified["inspectable"] is True


def test_standalone_count_token_is_classified_without_matching_discount() -> None:
    assert classify_filename("GSE317884_all_data_count_human.csv.gz")["role"] == "count_matrix"
    assert classify_filename("discount_results.csv.gz")["role"] == "unknown_table"
    assert (
        classify_filename("GSE1_Normalized_FPKM_gene_counts_matrix.txt.gz")["role"]
        == "normalized_expression_matrix"
    )


@pytest.mark.parametrize(
    "name",
    [
        "GSE1_counts_per_million.csv.gz",
        "GSE1_CPM_matrix.tsv.gz",
        "GSE1_logCPM.tsv.gz",
        "GSE1_log2CPM.tsv.gz",
        "GSE1_FPKM-UQ.tsv.gz",
        "GSE1_TMM_matrix.tsv.gz",
        "GSE1_voom_expression.tsv.gz",
        "GSE1_normalized_counts.csv.gz",
        "GSE1_normalised_counts.csv.gz",
        "GSE1_VST_counts_matrix.tsv.gz",
        "GSE1_rlog_counts_matrix.tsv.gz",
    ],
)
def test_normalized_count_like_filenames_are_not_classified_as_raw_counts(name: str) -> None:
    assert classify_filename(name)["role"] == "normalized_expression_matrix"


@pytest.mark.parametrize(
    "name",
    [
        "non_normalized_counts.csv",
        "not_normalized_counts.csv",
        "un-normalized_counts.csv",
        "unnormalized_counts.csv",
        "non_normalised_counts.csv",
        "raw_non_normalized_counts.csv",
    ],
)
def test_negated_normalization_tokens_remain_count_matrices(name: str) -> None:
    assert classify_filename(name)["role"] == "count_matrix"


def test_explicit_raw_counts_token_wins_over_conflicting_normalized_hint() -> None:
    assert classify_filename("GSE1_raw_counts_CPM.csv.gz")["role"] == "count_matrix"


def test_integer_like_cpm_stays_on_normalized_fallback_path(tmp_path: Path) -> None:
    name = "GSE1_counts_per_million.csv.gz"
    text = (
        "gene,c1,c2,t1,t2\n"
        "TP53,10,12,20,22\n"
        "VEGFA,30,32,60,64\n"
        "HIF1A,15,16,50,52\n"
        "HMOX1,8,9,40,41\n"
    )
    source = tmp_path / name
    source.write_bytes(gzip.compress(text.encode()))
    role = classify_filename(name)["role"]

    inspected = inspect_upstream_bytes(name, source.read_bytes(), declared_role=role)
    assert role == "normalized_expression_matrix"
    assert inspected["whole_number_share"] == 1.0
    assert inspected["declared_role"] == "normalized_expression_matrix"

    output = tmp_path / "cpm_derived.csv"
    summary = derive_welch_deg(
        source,
        output,
        role=role,
        gene_column="gene",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        normalized_scale="linear",
    )
    assert summary["source_input_type"] == "normalized_expression_matrix"
    assert summary["normalization"] == "public_normalized_matrix_log2x_plus_1"


def test_normalized_deg_results_still_take_the_author_table_path() -> None:
    from degora.discovery import _inspect_preparation_file

    name = "normalized_DEG_results.csv"
    classified = classify_filename(name)
    file_record = {**classified, "source_url": f"https://example.org/{name}"}
    payload = (
        b"gene,log2FoldChange,pvalue,padj\n"
        b"TP53,2.0,0.001,0.01\n"
        b"VEGFA,1.5,0.002,0.02\n"
    )

    inspection = _inspect_preparation_file(
        file_record,
        payload=payload,
        fetch_scope="full",
        accession="GSE1",
        target_dir=None,
    )

    assert classified["role"] == "normalized_expression_matrix"
    assert inspection["status"] == "ready_for_review"
    assert file_record["role"] == "deg_table"
    assert file_record["tier"] == "strong"


def test_deseq2_normalized_counts_stay_on_the_upstream_matrix_path() -> None:
    from degora.discovery import _inspect_preparation_file

    name = "DESeq2_normalized_counts.csv"
    classified = classify_filename(name)
    file_record = {**classified, "source_url": f"https://example.org/{name}"}
    payload = (
        b"gene,ctrl_1,ctrl_2,treat_1,treat_2\n"
        b"TP53,10,12,20,22\n"
        b"VEGFA,30,32,60,64\n"
    )

    inspection = _inspect_preparation_file(
        file_record,
        payload=payload,
        fetch_scope="full",
        accession="GSE1",
        target_dir=None,
    )

    assert classified["role"] == "normalized_expression_matrix"
    assert inspection["status"] == "upstream_matrix_ready_for_contrast"
    assert inspection["declared_role"] == "normalized_expression_matrix"
    assert file_record["role"] == "normalized_expression_matrix"


def test_nonactivatable_author_header_does_not_hide_a_ready_count_matrix() -> None:
    from degora.discovery import _inspect_preparation_file

    name = "raw_counts.csv"
    file_record = {
        **classify_filename(name),
        "source_url": f"https://example.org/{name}",
    }
    payload = (
        b"gene,log2FC,pvalue,c1,c2,t1,t2\n"
        b"TP53,NA,NA,10,12,20,22\n"
        b"VEGFA,NA,NA,30,32,60,64\n"
    )

    inspection = _inspect_preparation_file(
        file_record,
        payload=payload,
        fetch_scope="full",
        accession="GSE1",
        target_dir=None,
    )

    assert file_record["role"] == "count_matrix"
    assert inspection["status"] == "upstream_matrix_ready_for_contrast"
    assert inspection["sample_columns"] == ["c1", "c2", "t1", "t2"]


@pytest.mark.parametrize(
    ("column", "expected"),
    [
        ("sample_counts_per_million", ("normalized", "sample", "cpm")),
        ("sample_CPM", ("normalized", "sample", "cpm")),
        ("sample_logCPM", ("normalized", "sample", "logcpm")),
        ("sample_log2CPM", ("normalized", "sample", "logcpm")),
        ("sample_FPKM-UQ", ("normalized", "sample", "fpkm_uq")),
        ("sample_FPKM_unstranded", ("normalized", "sample", "fpkm")),
        ("sample_TPM_stranded_first", ("normalized", "sample", "tpm")),
        ("sample_TMM", ("normalized", "sample", "tmm")),
        ("sample_voom", ("normalized", "sample", "voom")),
        ("sample_raw_count", ("count", "sample", "count")),
        ("sample_raw_count_stranded_second", ("count", "sample", "count")),
    ],
)
def test_measurement_suffixes_share_one_canonical_vocabulary(
    column: str,
    expected: tuple[str, str, str],
) -> None:
    from degora.discovery import _measurement_column_parts

    assert _measurement_column_parts(column) == expected


def test_qualified_fpkm_and_tpm_columns_do_not_form_one_ready_pool() -> None:
    payload = (
        b"gene,c1_FPKM_unstranded,c2_FPKM_unstranded,t1_TPM_unstranded,t2_TPM_unstranded\n"
        b"TP53,1,1.2,200,220\n"
        b"VEGFA,2,2.2,400,420\n"
    )

    inspected = inspect_upstream_bytes(
        "mixed_normalized_matrix.csv",
        payload,
        declared_role="unknown_matrix",
    )

    assert inspected["status"] == "not_upstream_matrix"
    families = inspected["sample_column_families"]
    assert families["subtypes_present"]["normalized"] == ["fpkm", "tpm"]
    assert families["compatible_pools"] == {"normalized:fpkm": 2, "normalized:tpm": 2}
    assert "multiple measurement subtypes" in inspected["measurement_family_warning"]


def test_runtime_rederives_qualified_measurement_subtypes_from_an_old_bundle(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    source = bundle / "mixed_normalized_matrix.csv"
    sample_columns = [
        "c1_FPKM_unstranded",
        "c2_FPKM_unstranded",
        "t1_TPM_unstranded",
        "t2_TPM_unstranded",
    ]
    pd.DataFrame(
        {
            "gene": [f"GENE{i}" for i in range(30)],
            sample_columns[0]: [1.0] * 30,
            sample_columns[1]: [1.2] * 30,
            sample_columns[2]: [200.0] * 30,
            sample_columns[3]: [220.0] * 30,
        }
    ).to_csv(source, index=False)
    candidate = {
        "candidate_id": "matrix1",
        "name": source.name,
        "role": "normalized_expression_matrix",
        # Simulate a bundle created before qualified scale names were parsed.
        "inspection": {
            "status": "upstream_matrix_ready_for_contrast",
            "fetch_scope": "full",
            "local_path": str(source),
            "full_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "sample_columns": sample_columns,
            "gene_column": "gene",
            "header_row": 1,
        },
    }
    entry = {
        "candidate_id": "matrix1",
        "mode": "fallback",
        "direction_confirmed": True,
        "biological_replicates_confirmed": True,
        "control_samples": sample_columns[:2],
        "treatment_samples": sample_columns[2:],
        "matrix_type": "normalized_expression_matrix",
        "normalized_scale": "linear",
        "gene_column": "gene",
        "contrast_label": "treatment versus control",
    }

    with pytest.raises(DiscoveryError, match="mix normalized measurement subtypes"):
        _fallback_row(
            study={"accession": "GSE1", "title": "qualified scale fixture", "files": [candidate]},
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora discovery-analyze",
        )


def test_logcpm_pool_is_ready_but_warns_about_tpm_columns() -> None:
    payload = (
        b"gene,c1_logCPM,c2_log2CPM,t1_logCPM,t2_log2CPM,x1_TPM,x2_TPM\n"
        b"TP53,1,1.2,2,2.2,3,3.2\n"
        b"VEGFA,2,2.2,4,4.2,6,6.2\n"
    )

    inspected = inspect_upstream_bytes(
        "gene_counts_matrix.csv",
        payload,
        declared_role="normalized_expression_matrix",
    )

    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    families = inspected["sample_column_families"]
    assert families["subtypes_present"]["normalized"] == ["logcpm", "tpm"]
    assert families["compatible_pools"] == {"normalized:logcpm": 4, "normalized:tpm": 2}
    assert "multiple measurement subtypes" in inspected["measurement_family_warning"]


def test_server_rejects_logcpm_and_tpm_mix_even_when_one_pool_was_ready(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    source = bundle / "gene_counts_matrix.csv"
    pd.DataFrame(
        {
            "gene": [f"GENE{i}" for i in range(30)],
            "c1_logCPM": [1.0] * 30,
            "c2_log2CPM": [1.2] * 30,
            "x1_logCPM": [1.4] * 30,
            "x2_log2CPM": [1.6] * 30,
            "t1_TPM": [2.0] * 30,
            "t2_TPM": [2.2] * 30,
        }
    ).to_csv(source, index=False)
    inspected = inspect_upstream_bytes(
        source.name,
        source.read_bytes(),
        declared_role="normalized_expression_matrix",
    )
    candidate = {
        "candidate_id": "matrix1",
        "name": source.name,
        "role": "normalized_expression_matrix",
        "inspection": {
            **inspected,
            "fetch_scope": "full",
            "local_path": str(source),
            "full_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        },
    }
    study = {"accession": "GSE1", "title": "measurement suffix fixture", "files": [candidate]}
    entry = {
        "candidate_id": "matrix1",
        "mode": "fallback",
        "direction_confirmed": True,
        "biological_replicates_confirmed": True,
        "control_samples": ["c1_logCPM", "c2_log2CPM"],
        "treatment_samples": ["t1_TPM", "t2_TPM"],
        "matrix_type": "normalized_expression_matrix",
        "normalized_scale": "log2",
        "gene_column": "gene",
        "contrast_label": "treatment versus control",
    }

    with pytest.raises(DiscoveryError, match="mix normalized measurement subtypes"):
        _fallback_row(
            study=study,
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora discovery-analyze",
        )


def test_server_rejects_explicit_raw_counts_in_a_tampered_normalized_bundle(tmp_path: Path) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    source = bundle / "normalized_counts.csv"
    sample_columns = ["c1_raw_count", "c2_raw_count", "t1_raw_count", "t2_raw_count"]
    pd.DataFrame(
        {
            "gene": [f"GENE{i}" for i in range(30)],
            **{column: list(range(100, 130)) for column in sample_columns},
        }
    ).to_csv(source, index=False)
    inspected = inspect_upstream_bytes(
        source.name,
        source.read_bytes(),
        declared_role="normalized_expression_matrix",
    )
    candidate = {
        "candidate_id": "matrix1",
        "name": source.name,
        "role": "normalized_expression_matrix",
        "inspection": {
            **inspected,
            "fetch_scope": "full",
            "local_path": str(source),
            "full_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        },
    }
    entry = {
        "candidate_id": "matrix1",
        "mode": "fallback",
        "direction_confirmed": True,
        "biological_replicates_confirmed": True,
        "control_samples": sample_columns[:2],
        "treatment_samples": sample_columns[2:],
        "matrix_type": "normalized_expression_matrix",
        "normalized_scale": "linear",
        "gene_column": "gene",
        "contrast_label": "treatment versus control",
    }

    with pytest.raises(DiscoveryError, match="explicit raw-count.*cannot use a normalized matrix role"):
        _fallback_row(
            study={"accession": "GSE1", "title": "tampered role fixture", "files": [candidate]},
            candidate=candidate,
            entry=entry,
            spec=normalize_species("human"),
            bundle_root=bundle,
            derived_dir=tmp_path / "derived",
            sequence=1,
            replay_command="degora discovery-analyze",
        )


def test_explicit_normalized_columns_correct_a_misleading_count_filename_end_to_end(tmp_path: Path) -> None:
    from degora.discovery import _inspect_preparation_file, normalize_species
    from degora.discovery_run import _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    name = "gene_counts_matrix.csv"
    payload = pd.DataFrame(
        {
            "gene": [f"GENE{i}" for i in range(30)],
            "c1_FPKM": [1.0] * 30,
            "c2_FPKM": [1.2] * 30,
            "t1_FPKM": [2.0] * 30,
            "t2_FPKM": [2.2] * 30,
        }
    ).to_csv(index=False).encode()
    file_record = {
        **classify_filename(name),
        "candidate_id": "matrix1",
        "source_url": f"https://example.org/{name}",
    }
    assert file_record["role"] == "count_matrix"
    inspection = _inspect_preparation_file(
        file_record,
        payload=payload,
        fetch_scope="full",
        accession="GSE1",
        target_dir=bundle,
    )
    candidate = {**file_record, "inspection": inspection}
    study = {"accession": "GSE1", "title": "role recovery fixture", "files": [candidate]}
    entry = {
        "candidate_id": "matrix1",
        "mode": "fallback",
        "direction_confirmed": True,
        "biological_replicates_confirmed": True,
        "control_samples": ["c1_FPKM", "c2_FPKM"],
        "treatment_samples": ["t1_FPKM", "t2_FPKM"],
        "matrix_type": "normalized_expression_matrix",
        "normalized_scale": "linear",
        "gene_column": "gene",
        "contrast_label": "treatment versus control",
    }

    row, summary = _fallback_row(
        study=study,
        candidate=candidate,
        entry=entry,
        spec=normalize_species("human"),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora discovery-analyze",
    )

    assert file_record["role"] == "normalized_expression_matrix"
    assert inspection["declared_role"] == "normalized_expression_matrix"
    assert inspection["role_resolution"]["declared_role"] == "count_matrix"
    assert row["source_input_type"] == "normalized_expression_matrix"
    assert summary["normalization"] == "public_normalized_matrix_log2x_plus_1"


@pytest.mark.parametrize(
    ("name", "columns", "expected_role"),
    [
        (
            "normalized_matrix.csv",
            ["c1_raw_count", "c2_raw_count", "t1_raw_count", "t2_raw_count"],
            "count_matrix",
        ),
        (
            "normalized_counts.csv",
            ["c1_count", "c2_count", "t1_count", "t2_count"],
            "normalized_expression_matrix",
        ),
    ],
)
def test_role_recovery_requires_explicit_raw_count_suffixes(
    name: str,
    columns: list[str],
    expected_role: str,
) -> None:
    from degora.discovery import _inspect_preparation_file

    payload = pd.DataFrame(
        {
            "gene": ["TP53", "VEGFA"],
            **{column: [10, 20] for column in columns},
        }
    ).to_csv(index=False).encode()
    file_record = {
        **classify_filename(name),
        "source_url": f"https://example.org/{name}",
    }

    inspection = _inspect_preparation_file(
        file_record,
        payload=payload,
        fetch_scope="full",
        accession="GSE1",
        target_dir=None,
    )

    assert file_record["role"] == expected_role
    assert inspection["declared_role"] == expected_role


@pytest.mark.parametrize(
    ("separator", "description"),
    [
        ("\t", "alpha,beta,gamma,delta,epsilon"),
        (",", "alpha;beta;gamma;delta;epsilon"),
    ],
)
def test_author_delimiter_prefers_distinct_gene_effect_p_columns(
    separator: str,
    description: str,
) -> None:
    from degora.discovery import inspect_candidate_bytes

    text = (
        separator.join(["gene", "log2FoldChange", "pvalue", "description"])
        + "\n"
        + separator.join(["TP53", "2.0", "0.001", description])
        + "\n"
        + separator.join(["VEGFA", "1.5", "0.002", description])
        + "\n"
    )

    inspected = inspect_candidate_bytes("author_results.txt", text.encode())

    assert inspected["status"] == "ready_for_review"
    assert inspected["mapping"] == {
        "gene_column": "gene",
        "lfc_column": "log2FoldChange",
        "p_column": "pvalue",
        "padj_column": "",
    }


def test_tab_author_table_with_comma_rich_annotation_materializes_end_to_end(tmp_path: Path) -> None:
    from degora.discovery import inspect_candidate_bytes
    from degora.discovery_run import _materialize_author_table

    source = tmp_path / "rich_annotation.tsv"
    source.write_text(
        "gene\tlog2FoldChange\tpvalue\tdescription\n"
        "TP53\t2.0\t0.001\talpha,beta,gamma,delta,epsilon\n"
        "VEGFA\t1.5\t0.002\talpha,beta,gamma,delta,epsilon\n",
        encoding="utf-8",
    )
    inspection = inspect_candidate_bytes(source.name, source.read_bytes())

    output, provenance, _ = _materialize_author_table(
        study={"accession": "GSE_RICH"},
        candidate={"candidate_id": "rich1", "name": source.name, "inspection": inspection},
        entry={},
        source_path=source,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora discovery-analyze",
    )

    assert inspection["status"] == "ready_for_review"
    assert pd.read_csv(output)["gene_symbol"].tolist() == ["TP53", "VEGFA"]
    assert provenance["n_usable_output_rows"] == 2


@pytest.mark.parametrize(
    ("rows", "expected_transform", "expected_space", "expected_genes"),
    [
        (
            [
                "ENSG00000000001_TP53,2.0,0.001",
                "ENSG00000000002_EGFR,1.5,0.002",
                "ENSG00000000003_VEGFA,1.2,0.003",
            ],
            "ensembl_gene_symbol_suffix_v1",
            "gene symbol",
            ["TP53", "EGFR", "VEGFA"],
        ),
        (
            [
                "ENSG00000000001_TP53,2.0,0.001",
                "ENSG00000000002_ENSG00000000002,1.5,0.002",
                "ENSG00000000003_NA,1.2,0.003",
                "ENSG00000000004_NULL,1.1,0.004",
            ],
            "ensembl_gene_prefix_v1",
            "Ensembl gene ID",
            [
                "ENSG00000000001",
                "ENSG00000000002",
                "ENSG00000000003",
                "ENSG00000000004",
            ],
        ),
    ],
)
def test_author_composite_gene_ids_materialize_in_one_validated_space(
    tmp_path: Path,
    rows: list[str],
    expected_transform: str,
    expected_space: str,
    expected_genes: list[str],
) -> None:
    from degora.discovery import inspect_candidate_bytes
    from degora.discovery_run import _materialize_author_table

    source = tmp_path / f"{expected_transform}.csv"
    source.write_text(
        "EnsemblGene_GeneSymbol,log2FoldChange,pvalue\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    inspection = inspect_candidate_bytes(source.name, source.read_bytes())

    output, provenance, _ = _materialize_author_table(
        study={"accession": "GSE_COMPOSITE"},
        candidate={"candidate_id": "composite1", "name": source.name, "inspection": inspection},
        entry={},
        source_path=source,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora discovery-analyze",
    )

    assert inspection["status"] == "ready_for_review"
    assert inspection["gene_identifier_space"] == expected_space
    assert inspection["gene_value_transform"]["name"] == expected_transform
    assert pd.read_csv(output)["gene_symbol"].tolist() == expected_genes
    assert provenance["gene_value_transform"] == expected_transform
    assert provenance["gene_value_transform_output_space"] == expected_space


def test_composite_ensembl_symbol_alias_is_narrow() -> None:
    text = (
        "EnsemblGene_GeneScore\tS1\tS2\tS3\tS4\n"
        "ENSG00000000003_TSPAN6\t1\t2\t3\t4\n"
        "ENSG00000000005_TNMD\t2\t3\t4\t5\n"
    )
    result = inspect_upstream_bytes(
        "not_a_gene_symbol_alias_raw_counts.tsv.gz",
        gzip.compress(text.encode()),
        declared_role="count_matrix",
    )

    assert result["status"] == "not_upstream_matrix"
    assert result["gene_column"] == ""


def test_composite_header_without_validated_values_fails_closed() -> None:
    text = (
        "EnsemblGene_GeneSymbol\tS1\tS2\tS3\tS4\n"
        "ENSG00000000003\t1\t2\t3\t4\n"
        "TP53\t2\t3\t4\t5\n"
    )

    result = inspect_upstream_bytes(
        "invalid_composite_raw_counts.tsv.gz",
        gzip.compress(text.encode()),
        declared_role="count_matrix",
    )

    assert result["status"] == "not_upstream_matrix"
    assert "composite Ensembl/gene-symbol identifiers require" in result["reason"]


@pytest.mark.parametrize(
    ("name", "text", "role", "control", "treatment", "scale", "separator"),
    [
        (
            "GSE304706_geo_sub_raw_count.csv.gz",
            ',"A1","A3","A4","A6"\n'
            '"Zc3h14",1028,1236,2106,1947\n'
            '"Troap",340,370,623,730\n'
            '"Trp53",510,530,880,910\n'
            '"Vegfa",220,240,540,560\n',
            "count_matrix",
            ["A1", "A3"],
            ["A4", "A6"],
            None,
            ",",
        ),
        (
            "GSE317884_all_data_fpkms_human.csv.gz",
            ",NC0_1,NC0_2,NC4_1,NC4_2\n"
            "TSPAN6,1.1,1.3,4.0,4.2\n"
            "TNMD,2.0,2.2,5.0,5.1\n"
            "TP53,3.0,3.2,8.0,8.5\n"
            "VEGFA,1.5,1.7,6.0,6.2\n",
            "normalized_expression_matrix",
            ["NC0_1", "NC0_2"],
            ["NC4_1", "NC4_2"],
            "linear",
            ",",
        ),
        (
            "GSE337799_rna_counts.tsv.gz",
            '"HV0021_C_1h_rep1" "HV0037_C_1h_rep2" "HV0023_HS_1h_rep1" "HV0039_HS_1h_rep2"\n'
            '"ENSG00000000003.14" 1991 2916 4189 4195\n'
            '"ENSG00000141510.18" 500 550 1700 1800\n'
            '"ENSG00000146648.19" 320 360 950 1020\n'
            '"ENSG00000171862.12" 710 760 1900 2010\n',
            "count_matrix",
            ["HV0021_C_1h_rep1", "HV0037_C_1h_rep2"],
            ["HV0023_HS_1h_rep1", "HV0039_HS_1h_rep2"],
            None,
            " ",
        ),
        (
            "GSE318560_raw_counts_matrix.txt.gz",
            '"HC_148"\t"HC_025"\t"AD_061"\t"AD_128"\n'
            '"DDX11L1"\t60\t70\t220\t240\n'
            '"WASH7P"\t510\t440\t1810\t1720\n'
            '"TP53"\t330\t370\t980\t1010\n'
            '"VEGFA"\t210\t230\t740\t790\n',
            "count_matrix",
            ["HC_148", "HC_025"],
            ["AD_061", "AD_128"],
            None,
            "\t",
        ),
    ],
)
def test_public_r_row_names_survive_inspection_and_fallback_derivation(
    tmp_path: Path,
    name: str,
    text: str,
    role: str,
    control: list[str],
    treatment: list[str],
    scale: str | None,
    separator: str,
) -> None:
    source = tmp_path / name
    source.write_bytes(gzip.compress(text.encode()))

    inspected = inspect_upstream_bytes(name, source.read_bytes(), declared_role=role)
    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    assert inspected["gene_column"] == "row_name"
    assert sniff_delimited_separator(source, matrix_table=True) == separator
    assert "row_name" in read_matrix_frame(source).columns

    output = tmp_path / f"{name}.derived.csv"
    summary = derive_welch_deg(
        source,
        output,
        role=role,
        gene_column="row_name",
        control_samples=control,
        treatment_samples=treatment,
        normalized_scale=scale,
    )

    derived = pd.read_csv(output)
    assert summary["n_gene_rows"] >= 2
    assert derived["gene_symbol"].notna().all()
    assert derived["log2FoldChange"].notna().all()


@pytest.mark.parametrize("compressed", [False, True])
def test_xlsx_matrix_inspection_matches_fallback_reader(tmp_path: Path, compressed: bool) -> None:
    workbook = tmp_path / "matrix.xlsx"
    pd.DataFrame(
        {
            "gene": ["TP53", "VEGFA", "HIF1A", "HMOX1"],
            "c1": [101, 102, 103, 104],
            "c2": [111, 112, 113, 114],
            "t1": [301, 302, 303, 304],
            "t2": [321, 322, 323, 324],
        }
    ).to_excel(workbook, index=False)
    if compressed:
        source = tmp_path / "matrix.xlsx.gz"
        source.write_bytes(gzip.compress(workbook.read_bytes()))
    else:
        source = workbook

    inspected = inspect_upstream_bytes(source.name, source.read_bytes(), declared_role="count_matrix")
    output = tmp_path / f"{source.name}.derived.csv"
    summary = derive_welch_deg(
        source,
        output,
        role="count_matrix",
        gene_column="gene",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
    )

    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    assert inspected["sample_columns"] == ["c1", "c2", "t1", "t2"]
    assert summary["n_gene_rows"] == 4


@pytest.mark.parametrize("compressed", [False, True])
def test_legacy_xls_matrix_routes_through_bounded_reader(
    monkeypatch,
    compressed: bool,
) -> None:
    frame = pd.DataFrame(
        [
            ["gene", "c1", "c2", "t1", "t2"],
            ["TP53", 101, 111, 301, 321],
            ["VEGFA", 102, 112, 302, 322],
        ]
    )

    def fake_read_excel(*args, **kwargs):
        assert kwargs["sheet_name"] is None
        assert kwargs["header"] is None
        return {"matrix": frame}

    monkeypatch.setattr(pd, "read_excel", fake_read_excel)
    workbook_payload = b"\xd0\xcf\x11\xe0" + b"legacy workbook fixture"
    payload = gzip.compress(workbook_payload) if compressed else workbook_payload
    name = "matrix.xls.gz" if compressed else "matrix.xls"

    inspected = inspect_upstream_bytes(name, payload, declared_role="count_matrix")

    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    assert inspected["gene_column"] == "gene"
    assert inspected["sample_columns"] == ["c1", "c2", "t1", "t2"]


@pytest.mark.parametrize(
    ("prefix", "symbols"),
    [
        ("ENSG", ["TSPAN6", "TNMD", "TP53", "VEGFA"]),
        ("ENSMUSG", ["Trp53", "Vegfa", "Hif1a", "Hmox1"]),
    ],
)
def test_composite_ensembl_symbol_transform_is_applied_and_provenanced(
    tmp_path: Path,
    prefix: str,
    symbols: list[str],
) -> None:
    rows = ["EnsemblGene_GeneSymbol\tc1\tc2\tt1\tt2"]
    for index, symbol in enumerate(symbols, start=1):
        rows.append(f"{prefix}{index:011d}_{symbol}\t{100 + index}\t{110 + index}\t{300 + index}\t{320 + index}")
    text = "\n".join(rows) + "\n"
    source = tmp_path / f"{prefix}_gene_count_mat.txt.gz"
    source.write_bytes(gzip.compress(text.encode()))

    inspected = inspect_upstream_bytes(source.name, source.read_bytes(), declared_role="count_matrix")
    transform = inspected["gene_value_transform"]
    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    assert transform["name"] == "ensembl_gene_symbol_suffix_v1"
    assert transform["match_fraction"] == 1.0

    output = tmp_path / f"{prefix}.derived.csv"
    summary = derive_welch_deg(
        source,
        output,
        role="count_matrix",
        gene_column="EnsemblGene_GeneSymbol",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        gene_value_transform=transform["name"],
    )

    assert set(pd.read_csv(output)["gene_symbol"]) == {symbol.upper() for symbol in symbols}
    assert summary["gene_value_transform"] == "ensembl_gene_symbol_suffix_v1"
    provenance = json.loads(Path(str(output) + ".provenance.json").read_text(encoding="utf-8"))
    assert provenance["metadata"]["gene_value_transform"] == "ensembl_gene_symbol_suffix_v1"


def test_composite_transform_uses_uniform_ensembl_prefix_when_symbols_are_sparse(tmp_path: Path) -> None:
    text = (
        "EnsemblGene_GeneSymbol\tc1\tc2\tt1\tt2\n"
        "ENSG00000000001_TP53\t101\t111\t301\t321\n"
        "ENSG00000000002_ENSG00000000002\t102\t112\t302\t322\n"
        "ENSG00000000003_7SK\t103\t113\t303\t323\n"
        "ENSG00000000004_VEGFA\t104\t114\t304\t324\n"
    )
    source = tmp_path / "composite_counts.tsv.gz"
    source.write_bytes(gzip.compress(text.encode()))
    inspected = inspect_upstream_bytes(source.name, source.read_bytes(), declared_role="count_matrix")

    output = tmp_path / "composite_derived.csv"
    summary = derive_welch_deg(
        source,
        output,
        role="count_matrix",
        gene_column="EnsemblGene_GeneSymbol",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        gene_value_transform=inspected["gene_value_transform"]["name"],
    )

    assert inspected["gene_value_transform"]["name"] == "ensembl_gene_prefix_v1"
    assert set(pd.read_csv(output)["gene_symbol"]) == {
        "ENSG00000000001",
        "ENSG00000000002",
        "ENSG00000000003",
        "ENSG00000000004",
    }
    assert summary["gene_value_transform_output_space"] == "Ensembl gene ID"
    assert summary["gene_value_transform_usable_identifiers"] == 4


def test_composite_missing_symbol_tokens_use_validated_ensembl_prefix(tmp_path: Path) -> None:
    text = (
        "EnsemblGene_GeneSymbol\tc1\tc2\tt1\tt2\n"
        "ENSG00000000001_ENSG00000000001\t101\t111\t301\t321\n"
        "ENSG00000000002_NA\t102\t112\t302\t322\n"
        "ENSG00000000003_NAN\t103\t113\t303\t323\n"
        "ENSG00000000004_NULL\t104\t114\t304\t324\n"
    )
    source = tmp_path / "missing_symbols.tsv.gz"
    source.write_bytes(gzip.compress(text.encode()))
    inspected = inspect_upstream_bytes(source.name, source.read_bytes(), declared_role="count_matrix")

    assert inspected["status"] == "upstream_matrix_ready_for_contrast"
    assert inspected["gene_value_transform"]["name"] == "ensembl_gene_prefix_v1"
    output = tmp_path / "missing_symbols_derived.csv"
    summary = derive_welch_deg(
        source,
        output,
        role="count_matrix",
        gene_column="EnsemblGene_GeneSymbol",
        control_samples=["c1", "c2"],
        treatment_samples=["t1", "t2"],
        gene_value_transform=inspected["gene_value_transform"]["name"],
    )

    assert summary["n_valid_rows"] == 4
    assert set(pd.read_csv(output)["gene_symbol"]) == {
        "ENSG00000000001",
        "ENSG00000000002",
        "ENSG00000000003",
        "ENSG00000000004",
    }


@pytest.mark.parametrize("selected_gene_column", ["EnsemblGene_GeneSymbol", "gene_symbol"])
def test_fallback_applies_composite_transform_only_to_its_source_column(
    tmp_path: Path,
    selected_gene_column: str,
) -> None:
    from degora.discovery import normalize_species
    from degora.discovery_run import _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    source = bundle / "GSE315097_gene_count_mat.csv"
    symbols = [f"GENE{index}" for index in range(1, 31)]
    pd.DataFrame(
        {
            "EnsemblGene_GeneSymbol": [f"ENSG{index:011d}_{symbol}" for index, symbol in enumerate(symbols, 1)],
            "gene_symbol": symbols,
            "c1": range(101, 131),
            "c2": range(111, 141),
            "t1": range(301, 331),
            "t2": range(321, 351),
        }
    ).to_csv(source, index=False)
    candidate = {
        "candidate_id": "matrix1",
        "name": source.name,
        "role": "count_matrix",
        "inspection": {
            "status": "upstream_matrix_ready_for_contrast",
            "fetch_scope": "full",
            "local_path": str(source),
            "full_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "sample_columns": ["c1", "c2", "t1", "t2"],
            "gene_column": "EnsemblGene_GeneSymbol",
            "header_row": 1,
            "gene_value_transform": {
                "name": "ensembl_gene_symbol_suffix_v1",
                "source_header": "EnsemblGene_GeneSymbol",
                "match_fraction": 1.0,
                "sample_rows": 20,
            },
        },
    }
    study = {"accession": "GSE315097", "title": "composite identifier fixture", "files": [candidate]}
    entry = {
        "candidate_id": "matrix1",
        "mode": "fallback",
        "direction_confirmed": True,
        "biological_replicates_confirmed": True,
        "control_samples": ["c1", "c2"],
        "treatment_samples": ["t1", "t2"],
        "matrix_type": "count_matrix",
        "gene_column": selected_gene_column,
        "contrast_label": "treatment versus control",
    }

    row, summary = _fallback_row(
        study=study,
        candidate=candidate,
        entry=entry,
        spec=normalize_species("human"),
        bundle_root=bundle,
        derived_dir=tmp_path / "derived",
        sequence=1,
        replay_command="degora discovery-analyze",
    )

    assert set(pd.read_csv(summary["output_path"])["gene_symbol"]) == set(symbols)
    if selected_gene_column == "EnsemblGene_GeneSymbol":
        assert summary["gene_value_transform"] == "ensembl_gene_symbol_suffix_v1"
        assert "normalized uniformly to gene symbol" in row["notes"]
    else:
        assert "gene_value_transform" not in summary
        assert "normalized uniformly" not in row["notes"]


def test_legacy_geo_prepare_says_when_zero_candidates_are_usable(
    tmp_path, monkeypatch, capsys
) -> None:
    import degora.discovery as discovery

    def fake_prepare(*args, **kwargs):
        output = kwargs["materialize_dir"]
        return {
            "returned_studies": 2,
            "studies": [
                {
                    "accession": "GSE1",
                    "ready_for_review_count": 0,
                    "upstream_matrix_count": 0,
                    "files": [
                        {
                            "inspection": {
                                "status": "not_deg_table",
                                "reason": "missing required columns",
                            }
                        }
                    ],
                },
                {
                    "accession": "GSE2",
                    "ready_for_review_count": 0,
                    "upstream_matrix_count": 0,
                    "files": [],
                },
            ],
            "exports": {"draft_catalog_csv": str(output / "DEGORA_discovery_draft_catalog.csv")},
        }

    monkeypatch.setattr(discovery, "prepare_geo_studies", fake_prepare)

    assert main(
        [
            "discover",
            "hypoxia",
            "--source",
            "geo",
            "--species",
            "human",
            "--select",
            "GSE1",
            "--select",
            "GSE2",
            "--output-dir",
            str(tmp_path / "prepared"),
        ]
    ) == 0

    output = capsys.readouterr().out
    assert "Prepared 2 human GEO studies" in output
    assert "Ready for review: 0 table(s); upstream matrices awaiting contrast: 0" in output
    assert "Usable candidates: 0 - nothing can be activated yet" in output
    assert "GSE1: no usable table" in output
    assert "GSE2: no supplementary files" in output


@pytest.mark.parametrize(
    "status",
    [
        "requires_column_mapping",
        "requires_lfc_confirmation",
        "requires_pvalue_mapping",
    ],
)
def test_legacy_geo_prepare_counts_review_required_author_tables_as_usable(
    tmp_path: Path,
    monkeypatch,
    capsys,
    status: str,
) -> None:
    import degora.discovery as discovery

    def fake_prepare(*args, **kwargs):
        output = kwargs["materialize_dir"]
        return {
            "returned_studies": 1,
            "studies": [
                {
                    "accession": "GSE_REVIEW",
                    "ready_for_review_count": 0,
                    "upstream_matrix_count": 0,
                    "files": [
                        {
                            "inspection": {
                                "status": status,
                                "reason": "an author table needs a reviewer choice",
                            }
                        }
                    ],
                }
            ],
            "exports": {"draft_catalog_csv": str(output / "DEGORA_discovery_draft_catalog.csv")},
        }

    monkeypatch.setattr(discovery, "prepare_geo_studies", fake_prepare)

    assert main(
        [
            "discover",
            "hypoxia",
            "--source",
            "geo",
            "--species",
            "human",
            "--select",
            "GSE_REVIEW",
            "--output-dir",
            str(tmp_path / status),
        ]
    ) == 0

    output = capsys.readouterr().out
    assert "Ready for review: 0 table(s); upstream matrices awaiting contrast: 0" in output
    assert "Author tables awaiting column/scale review: 1." in output
    assert "Usable candidates: 1" in output
    assert "nothing can be activated yet" not in output


def test_legacy_geo_prepare_does_not_call_header_only_candidate_usable(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    from degora import discovery

    def fake_prepare(*args, **kwargs):
        output = kwargs["materialize_dir"]
        return {
            "returned_studies": 1,
            "studies": [
                {
                    "accession": "GSE_HEADER_ONLY",
                    "ready_for_review_count": 0,
                    "upstream_matrix_count": 0,
                    "files": [
                        {
                            "inspection": {
                                "status": "candidate_header",
                                "reason": "header matched but sampled values did not pass",
                            }
                        }
                    ],
                }
            ],
            "exports": {"draft_catalog_csv": str(output / "DEGORA_discovery_draft_catalog.csv")},
        }

    monkeypatch.setattr(discovery, "prepare_geo_studies", fake_prepare)

    assert main(
        [
            "discover",
            "hypoxia",
            "--source",
            "geo",
            "--species",
            "human",
            "--select",
            "GSE_HEADER_ONLY",
            "--output-dir",
            str(tmp_path / "header-only"),
        ]
    ) == 0

    output = capsys.readouterr().out
    assert "Author tables awaiting column/scale review: 0." in output
    assert "Usable candidates: 0 - nothing can be activated yet" in output
    assert "cannot be activated as prepared" in output


def test_prepare_status_prefers_a_usable_matrix_over_a_header_only_candidate() -> None:
    from degora.cli import _prepare_record_lines

    lines = _prepare_record_lines(
        [
            {
                "accession": "GSE_MIXED",
                "ready_for_review_count": 0,
                "upstream_matrix_count": 1,
                "files": [
                    {"inspection": {"status": "candidate_header", "reason": "header only"}},
                    {
                        "inspection": {
                            "status": "upstream_matrix_ready_for_contrast",
                            "reason": "matrix columns were detected",
                        }
                    },
                ],
            }
        ]
    )

    assert len(lines) == 1
    assert "expression matrices were found" in lines[0]
    assert "cannot be activated" not in lines[0]
