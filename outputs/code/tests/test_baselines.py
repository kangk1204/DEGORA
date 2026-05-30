from __future__ import annotations

import json

import degora.baselines as baselines
import pandas as pd
import pytest

from degora.baselines import (
    BASELINE_RESULT_COLUMNS,
    PUBLIC_SUMMARY_TOOL_INPUT_COLUMNS,
    awfisher_adapter,
    awmeta_deg_table_feasibility,
    awmeta_source_input_audit,
    BASELINE_MANIFEST_COLUMNS,
    FAILURE_LEDGER_COLUMNS,
    PARITY_MATRIX_COLUMNS,
    failure_ledger,
    hstouffer_materializer_feasibility,
    metarnaseq_fisher_adapter,
    metarnaseq_invnorm_adapter,
    metavolcanor_adapter,
    public_summary_tool_input_requirements,
    r_preflight_report,
    robustrankaggreg_adapter,
    run_tier0_baselines,
    validate_baseline_result,
    write_baseline_outputs,
    _study_matrix_inputs,
)


def _harmonized() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "gene_symbol": ["VEGFA", "RPL13A", "VEGFA", "RPL13A", "VEGFA", "HK2"],
            "signed_z": [5.0, 0.1, 4.0, -0.2, 3.5, 2.0],
            "lfc": [2.0, 0.1, 1.7, -0.1, 1.5, 1.0],
            "n_ctrl": [3, 3, 4, 4, 5, 5],
            "n_treat": [3, 3, 4, 4, 5, 5],
            "normalized_rank": [0.1, 0.9, 0.1, 0.8, 0.1, 0.2],
            "pvalue": [1e-6, 0.9, 1e-5, 0.8, 1e-4, 0.05],
        }
    )


def test_tier0_baselines_return_uniform_prd_schema() -> None:
    outputs = run_tier0_baselines(_harmonized(), min_studies=2)

    assert set(outputs) == {"degora_slice", "weighted_stouffer", "unweighted_stouffer", "fisher", "rank_product_approx", "sign_vote"}
    for frame in outputs.values():
        assert frame.columns.tolist() == BASELINE_RESULT_COLUMNS
        assert set(frame["status"]) == {"ok"}
        assert frame["rank"].min() == 1
        assert "VEGFA" in set(frame["symbol"])
        validate_baseline_result(frame)


def test_preflight_and_failure_ledger_mark_unavailable_methods_as_open_s1_blockers() -> None:
    preflight = r_preflight_report()
    ledger = failure_ledger("hypoxia", preflight)

    assert {"r_version", "rpy2_available", "interop_status", "packages"}.issubset(preflight)
    assert {"AWFisher", "metaRNASeq", "RankProd", "MetaDE", "DExMA", "MetaIntegrator"}.issubset(preflight["packages"])
    for package_status in preflight["packages"].values():
        assert {"available", "source", "version", "install_path", "failure_log", "message"}.issubset(package_status)
    assert ledger.columns.tolist() == FAILURE_LEDGER_COLUMNS
    assert "awmeta" in set(ledger["method_id"])
    assert {"rankprod_exact", "metade", "dexma", "metaintegrator"}.issubset(set(ledger["method_id"]))
    assert set(ledger["status"]) <= {"open_s1_blocker"}
    assert ledger.loc[ledger["method_id"].eq("awmeta"), "blocker_id"].iloc[0]


def test_awmeta_feasibility_blocks_nonfaithful_variance_imputation() -> None:
    feasibility = awmeta_deg_table_feasibility(_harmonized())

    assert feasibility["faithful"] is False
    assert feasibility["blocker_id"] == "awmeta_variance_inputs_missing"
    assert "deriving SE from p-values" in feasibility["message"]
    assert "derive SE from signed_z" in feasibility["prohibited_approximations"]
    assert "data_acquisition_requirement" in feasibility

    ledger = failure_ledger(
        "hypoxia",
        {"r_version": "R version", "packages": {"metafor": {"available": True}}},
        harmonized=_harmonized(),
    )
    awmeta = ledger.loc[ledger["method_id"].eq("awmeta")].iloc[0]
    assert awmeta["blocker_id"] == "awmeta_variance_inputs_missing"
    assert "non-faithful approximation" in awmeta["message"]


def test_invnorm_matrix_uses_neutral_missing_pvalues() -> None:
    inputs = _study_matrix_inputs(_harmonized(), min_studies=2, missing_pvalue_fill=0.5)

    p_matrix = inputs["p_matrix"]
    assert p_matrix.loc["RPL13A", "S3"] == 0.5


def test_metarnaseq_invnorm_adapter_passes_neutral_missing_pvalues(monkeypatch) -> None:
    captured: dict[str, pd.DataFrame] = {}

    class Completed:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run_rscript(args, *, timeout):
        p_path = args[2]
        out_path = args[3]
        p_matrix = pd.read_csv(p_path)
        captured["p_matrix"] = p_matrix
        pd.DataFrame(
            {
                "gene_symbol": p_matrix["gene_symbol"],
                "pvalue": [0.0] * len(p_matrix),
                "test_statistic": [float("inf")] * len(p_matrix),
            }
        ).to_csv(out_path, index=False)
        return Completed()

    monkeypatch.setattr(baselines, "_run_rscript", fake_run_rscript)
    monkeypatch.setattr(baselines, "_r_package_version", lambda package: "test")

    result = metarnaseq_invnorm_adapter(_harmonized(), min_studies=2)

    p_matrix = captured["p_matrix"].set_index("gene_symbol")
    assert p_matrix.loc["RPL13A", "S3"] == 0.5
    assert result["pvalue"].min() > 0.0
    assert set(result["method_id"]) == {"metarnaseq_invnorm"}


def test_metarnaseq_invnorm_adapter_blocks_all_neutral_sparse_tie(monkeypatch) -> None:
    class Completed:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run_rscript(args, *, timeout):
        p_path = args[2]
        out_path = args[3]
        p_matrix = pd.read_csv(p_path)
        pd.DataFrame(
            {
                "gene_symbol": p_matrix["gene_symbol"],
                "pvalue": [1.0] * len(p_matrix),
                "test_statistic": [0.0] * len(p_matrix),
            }
        ).to_csv(out_path, index=False)
        return Completed()

    monkeypatch.setattr(baselines, "_run_rscript", fake_run_rscript)

    with pytest.raises(RuntimeError, match="metarnaseq_invnorm_uninformative_sparse_public_summary"):
        metarnaseq_invnorm_adapter(_harmonized(), min_studies=2)


def test_awmeta_source_audit_credits_only_explicit_variance_inputs(tmp_path) -> None:
    source = tmp_path / "study.csv"
    pd.DataFrame({"gene": ["VEGFA"], "log2FoldChange": [2.0], "lfcSE": [0.4], "pvalue": [0.01]}).to_csv(source, index=False)
    frame = _harmonized().assign(source_path=str(source), pipeline="DESeq2")

    audit = awmeta_source_input_audit(frame)
    assert {row["study_id"] for row in audit} == {"S1", "S2", "S3"}
    assert all(row["variance_or_weight_columns"] == ["lfcSE"] for row in audit)
    assert all(row["effect_columns"] == ["log2FoldChange"] for row in audit)


def test_awmeta_feasibility_blocks_partial_uniform_variance_inputs() -> None:
    frame = _harmonized().copy()
    frame["lfcSE"] = [0.2, 0.2, None, None, None, None]

    feasibility = awmeta_deg_table_feasibility(frame)

    assert feasibility["faithful"] is False
    assert feasibility["blocker_id"] == "awmeta_incomplete_original_variance_inputs"
    assert feasibility["studies_with_harmonized_variance_candidates"] == ["S1"]
    assert set(feasibility["studies_missing_source_variance_candidates"]) == {"S2", "S3"}


def test_hstouffer_feasibility_blocks_surrogate_deseq2_materialization() -> None:
    report = hstouffer_materializer_feasibility(_harmonized())

    assert report["status"] == "blocked"
    assert report["source_pin"]["commit"] == "306e38c26919f19e7c3dfd6cd646005c502b3310"
    assert report["source_level_method_evidence"]
    assert report["can_materialize_without_imputation_or_filtering"] is False
    assert "lfcSE" in report["missing_original_fields"]
    assert any("signed_z" in blocker for blocker in report["blockers"])
    assert report["source_input_audit"] == []


def test_hstouffer_source_input_audit_records_missing_original_columns(tmp_path) -> None:
    source = tmp_path / "S1.tsv"
    source.write_text("id\tlog2FoldChange\tpvalue\tpadj\nVEGFA\t1.2\t0.01\t0.05\n")
    harmonized = _harmonized().assign(
        paper_id="P1",
        pipeline="DESeq2",
        source_path=str(source),
        source_url="file://S1.tsv",
    )

    report = hstouffer_materializer_feasibility(harmonized)

    assert report["blocked_source_studies"] == ["S1", "S2", "S3"]
    first = report["source_input_audit"][0]
    assert first["local_file_exists"] is True
    assert first["audit_status"] == "blocked_missing_required_columns"
    assert "lfcSE" in first["required_columns_missing"]


def test_public_summary_tool_input_requirements_explains_summary_only_limits() -> None:
    preflight = r_preflight_report()
    hstouffer = hstouffer_materializer_feasibility(_harmonized())
    awmeta = awmeta_deg_table_feasibility(_harmonized())
    ledger = failure_ledger("hypoxia", preflight, hstouffer_feasibility=hstouffer, harmonized=_harmonized())

    table = public_summary_tool_input_requirements(
        corpus="hypoxia",
        harmonized=_harmonized(),
        outputs=run_tier0_baselines(_harmonized(), min_studies=2),
        preflight=preflight,
        ledger=ledger,
        hstouffer_report=hstouffer,
        awmeta_report=awmeta,
    )

    assert table.columns.tolist() == PUBLIC_SUMMARY_TOOL_INPUT_COLUMNS
    assert {
        "awfisher",
        "metarnaseq_fisher",
        "rankprod_exact",
        "dexma",
        "metaintegrator",
        "omicc",
        "imageo",
        "networkanalyst",
        "deet",
        "creeds",
        "crossmeta",
        "generic_pvalue_combiners",
    }.issubset(set(table["method_id"]))
    rankprod = table.loc[table["method_id"].eq("rankprod_exact")].iloc[0]
    assert rankprod["can_run_from_current_public_files"] == "no"
    assert "replicate-level expression" in rankprod["required_inputs"]
    dexma = table.loc[table["method_id"].eq("dexma")].iloc[0]
    assert "phenotype" in dexma["missing_or_nonfaithful_inputs"]
    awmeta_row = table.loc[table["method_id"].eq("awmeta")].iloc[0]
    assert "signed z" in awmeta_row["faithful_adapter_decision"]
    omicc = table.loc[table["method_id"].eq("omicc")].iloc[0]
    assert omicc["current_pipeline_status"] == "workflow_or_resource_prior_art"
    assert "expression matrices" in omicc["missing_or_nonfaithful_inputs"]
    deet = table.loc[table["method_id"].eq("deet")].iloc[0]
    assert deet["manuscript_use"] == "database/resource prior art"
    combiners = table.loc[table["method_id"].eq("generic_pvalue_combiners")].iloc[0]
    assert "Fisher" in combiners["faithful_adapter_decision"]


def test_public_summary_tool_input_requirements_prefers_blocker_over_stale_output() -> None:
    outputs = run_tier0_baselines(_harmonized(), min_studies=2)
    outputs["metarnaseq_invnorm"] = pd.DataFrame(
        {
            "method_id": ["metarnaseq_invnorm"],
            "setting_id": ["default"],
            "gene_id": ["VEGFA"],
            "symbol": ["VEGFA"],
            "rank": [1],
            "score": [1.0],
            "pvalue": [1.0],
            "padj": [1.0],
            "effect": [0.0],
            "direction": ["zero"],
            "n_studies": [2],
            "missingness": [0.0],
            "runtime_s": [0.1],
            "version": ["metaRNASeq 1.0.8"],
            "status": ["ok"],
        }
    )
    ledger = pd.DataFrame(
        {
            "corpus": ["hypoxia"],
            "method_id": ["metarnaseq_invnorm"],
            "setting_id": ["default"],
            "tier": ["tier1"],
            "status": ["open_s1_blocker"],
            "blocker_id": ["metarnaseq_invnorm_uninformative_sparse_public_summary"],
            "message": ["uninformative sparse public-summary result"],
            "resolution": ["report blocked"],
        }
    )

    table = public_summary_tool_input_requirements(
        corpus="hypoxia",
        harmonized=_harmonized(),
        outputs=outputs,
        preflight=r_preflight_report(),
        ledger=ledger,
        hstouffer_report={},
        awmeta_report={},
    )

    invnorm = table.loc[table["method_id"].eq("metarnaseq_invnorm")].iloc[0]
    assert invnorm["current_pipeline_status"] == "blocked:metarnaseq_invnorm_uninformative_sparse_public_summary"
    assert invnorm["can_run_from_current_public_files"] == "blocked_on_current_public_files"
    assert invnorm["manuscript_use"] == "blocked comparator row"


def test_write_baseline_outputs_emits_required_files(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)

    summary = write_baseline_outputs(harmonized_path, tmp_path / "baselines", corpus="hypoxia")

    output_dir = tmp_path / "baselines"
    assert (output_dir / "r_preflight_report.json").exists()
    assert (output_dir / "baseline_failure_ledger.csv").exists()
    assert (output_dir / "baseline_parity_matrix.csv").exists()
    assert (output_dir / "baseline_parity_matrix.csv.source").exists()
    assert (output_dir / "public_summary_tool_input_requirements.csv").exists()
    assert (output_dir / "public_summary_tool_input_requirements.csv.source").exists()
    assert (output_dir / "public_summary_tool_input_requirements.md").exists()
    assert (output_dir / "public_summary_tool_input_requirements.md.source").exists()
    assert (output_dir / "baseline_manifest.csv").exists()
    assert (output_dir / "baseline_manifest.csv.source").exists()
    assert (output_dir / "hstouffer_feasibility_report.json").exists()
    assert (output_dir / "hstouffer_feasibility_report.json.source").exists()
    assert (output_dir / "awmeta_feasibility_report.json").exists()
    assert (output_dir / "awmeta_feasibility_report.json.source").exists()
    assert (output_dir / "hypoxia_weighted_stouffer_default.tsv").exists()
    assert (output_dir / "hypoxia_weighted_stouffer_default.tsv.source").exists()
    assert summary["n_open_s1_blockers"] >= 1
    json.loads((output_dir / "r_preflight_report.json").read_text())
    awmeta_report = json.loads((output_dir / "awmeta_feasibility_report.json").read_text())
    assert awmeta_report["blocker_id"] == "awmeta_variance_inputs_missing"

    parity = pd.read_csv(output_dir / "baseline_parity_matrix.csv")
    manifest = pd.read_csv(output_dir / "baseline_manifest.csv")
    input_requirements = pd.read_csv(output_dir / "public_summary_tool_input_requirements.csv")
    assert parity.columns.tolist() == PARITY_MATRIX_COLUMNS
    assert manifest.columns.tolist() == BASELINE_MANIFEST_COLUMNS
    assert input_requirements.columns.tolist() == PUBLIC_SUMMARY_TOOL_INPUT_COLUMNS
    assert {"robustrankaggreg", "metavolcanor", "awfisher", "metarnaseq_fisher", "metarnaseq_invnorm", "hstouffer", "awmeta"}.issubset(set(parity["method_id"]))
    assert {"rankprod_exact", "metade", "dexma", "metaintegrator", "omicc", "deet", "creeds"}.issubset(set(input_requirements["method_id"]))
    assert {"default", "tuned_min_studies_3"}.issubset(set(parity.loc[parity["method_id"].eq("hstouffer"), "setting_id"]))
    assert parity.loc[parity["method_id"].eq("hstouffer"), "version_or_commit"].str.contains("306e38c26919f19e7c3dfd6cd646005c502b3310").all()
    assert set(parity["claim_allowed"]) <= {"yes", "no"}
    assert {"hstouffer", "awmeta"}.issubset(set(manifest.loc[manifest["artifact_type"].eq("feasibility_report"), "method_id"]))
    assert "public_summary_tool_input_requirements" in set(manifest["artifact_type"])


def test_robustrankaggreg_adapter_emits_direct_prior_art_schema_when_available() -> None:
    preflight = r_preflight_report()
    if not preflight["packages"].get("RobustRankAggreg", {}).get("available"):
        pytest.skip("RobustRankAggreg R package is not available in this environment")

    result = robustrankaggreg_adapter(_harmonized(), min_studies=2)

    assert result.columns.tolist() == BASELINE_RESULT_COLUMNS
    assert set(result["method_id"]) == {"robustrankaggreg"}
    assert set(result["setting_id"]) == {"default"}
    assert set(result["status"]) == {"ok"}
    assert result["rank"].min() == 1
    assert result.iloc[0]["symbol"] == "VEGFA"
    assert result.loc[result["symbol"].eq("VEGFA"), "n_studies"].iloc[0] == 3
    assert result["version"].str.contains("RobustRankAggreg").all()
    validate_baseline_result(result)


def test_metavolcanor_adapter_emits_direct_prior_art_schema_when_available() -> None:
    preflight = r_preflight_report()
    if not preflight["packages"].get("MetaVolcanoR", {}).get("available"):
        pytest.skip("MetaVolcanoR R package is not available in this environment")

    result = metavolcanor_adapter(_harmonized(), min_studies=2)

    assert result.columns.tolist() == BASELINE_RESULT_COLUMNS
    assert set(result["method_id"]) == {"metavolcanor"}
    assert set(result["setting_id"]) == {"default"}
    assert set(result["status"]) == {"ok"}
    assert result["rank"].min() == 1
    assert "VEGFA" in set(result["symbol"])
    assert result.loc[result["symbol"].eq("VEGFA"), "n_studies"].iloc[0] == 3
    assert result.loc[result["symbol"].eq("VEGFA"), "effect"].iloc[0] > 0
    assert result["version"].str.contains("MetaVolcanoR").all()
    validate_baseline_result(result)


def test_awfisher_adapter_emits_direct_prior_art_schema_when_available() -> None:
    preflight = r_preflight_report()
    if not preflight["packages"].get("AWFisher", {}).get("available"):
        pytest.skip("AWFisher R package is not available in this environment")

    result = awfisher_adapter(_harmonized(), min_studies=2)

    assert result.columns.tolist() == BASELINE_RESULT_COLUMNS
    assert set(result["method_id"]) == {"awfisher"}
    assert set(result["setting_id"]) == {"default"}
    assert set(result["status"]) == {"ok"}
    assert result["rank"].min() == 1
    assert "VEGFA" in set(result["symbol"])
    assert result.loc[result["symbol"].eq("VEGFA"), "n_studies"].iloc[0] == 3
    assert result["version"].str.contains("AWFisher").all()
    validate_baseline_result(result)


@pytest.mark.parametrize(
    ("package_method", "adapter"),
    [
        ("metarnaseq_fisher", metarnaseq_fisher_adapter),
        ("metarnaseq_invnorm", metarnaseq_invnorm_adapter),
    ],
)
def test_metarnaseq_adapters_emit_direct_prior_art_schema_when_available(package_method, adapter) -> None:
    preflight = r_preflight_report()
    if not preflight["packages"].get("metaRNASeq", {}).get("available"):
        pytest.skip("metaRNASeq R package is not available in this environment")

    result = adapter(_harmonized(), min_studies=2)

    assert result.columns.tolist() == BASELINE_RESULT_COLUMNS
    assert set(result["method_id"]) == {package_method}
    assert set(result["setting_id"]) == {"default"}
    assert set(result["status"]) == {"ok"}
    assert result["rank"].min() == 1
    assert "VEGFA" in set(result["symbol"])
    assert result.loc[result["symbol"].eq("VEGFA"), "n_studies"].iloc[0] == 3
    assert result["version"].str.contains("metaRNASeq").all()
    validate_baseline_result(result)


def test_write_baseline_outputs_emits_metavolcanor_when_available(tmp_path) -> None:
    preflight = r_preflight_report()
    if not preflight["packages"].get("MetaVolcanoR", {}).get("available"):
        pytest.skip("MetaVolcanoR R package is not available in this environment")
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)

    write_baseline_outputs(harmonized_path, tmp_path / "baselines", corpus="hypoxia")

    output_dir = tmp_path / "baselines"
    result_path = output_dir / "hypoxia_metavolcanor_default.tsv"
    assert result_path.exists()
    assert result_path.with_suffix(result_path.suffix + ".source").exists()
    result = pd.read_csv(result_path, sep="\t")
    assert result.columns.tolist() == BASELINE_RESULT_COLUMNS
    assert set(result["method_id"]) == {"metavolcanor"}
    ledger = pd.read_csv(output_dir / "baseline_failure_ledger.csv")
    assert "metavolcanor_wrapper_missing" not in set(ledger["blocker_id"])
    parity = pd.read_csv(output_dir / "baseline_parity_matrix.csv")
    mv = parity.loc[(parity["method_id"].eq("metavolcanor")) & (parity["setting_id"].eq("default"))]
    assert mv["output_schema_status"].iloc[0] == "uniform_schema_validated"
    assert mv["run_status"].iloc[0] == "ok"
