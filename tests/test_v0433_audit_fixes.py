"""Regression tests for the 0.4.32 audit findings.

Each test names the behaviour that was wrong and asserts the behaviour that
replaced it. Where a finding came from an external audit report the source is
named in the docstring.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from degora.aggregate import time_course_selection_report
from degora.cli import main
from degora.harmonize import normalize_lfc_scale
from degora.score_db import DIRECTION_CONFLICT_RULE, degora_score_table, direction_conflict_warnings
from degora.slice_runner import DegoraConfigError, MAX_DECLARED_GROUP_SIZE, run_slice, validate_catalog_inputs


GENES = ["ISG15", "IFIT1", "MX1", "OAS1", "STAT1", "IRF7", "CXCL10", "DDX58"] + [
    f"G{index:04d}" for index in range(92)
]


def _table(path, seed: int, *, sign: float = 1.0, genes=None) -> None:
    rng = np.random.default_rng(seed)
    labels = list(genes or GENES)
    count = len(labels)
    lfc = rng.normal(0.0, 0.4, count)
    lfc[:8] = sign * np.abs(rng.normal(1.8, 0.3, 8))
    pvalue = np.clip(rng.uniform(0.0, 1.0, count) ** (1 + 2 * np.abs(lfc)), 1e-30, 1.0)
    pd.DataFrame(
        {
            "gene": labels,
            "log2FoldChange": lfc,
            "pvalue": pvalue,
            "padj": np.clip(pvalue * 2, 0, 1),
        }
    ).to_csv(path, index=False)


def _config(tmp_path, rows) -> "object":
    path = tmp_path / "config.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def _two_units(tmp_path, **first_row_overrides) -> "object":
    _table(tmp_path / "a.csv", 1)
    _table(tmp_path / "b.csv", 2)
    rows = [
        {
            "study_id": "A_6h",
            "source_unit_id": "UA",
            "source_path": str(tmp_path / "a.csv"),
            "gene_column": "gene",
            "lfc_column": "log2FoldChange",
            "p_column": "pvalue",
            "padj_column": "padj",
            "time_h": 6,
            "n_ctrl": 3,
            "n_treat": 3,
        },
        {
            "study_id": "B_6h",
            "source_unit_id": "UB",
            "source_path": str(tmp_path / "b.csv"),
            "gene_column": "gene",
            "lfc_column": "log2FoldChange",
            "p_column": "pvalue",
            "padj_column": "padj",
            "time_h": 6,
            "n_ctrl": 3,
            "n_treat": 3,
        },
    ]
    rows[0].update(first_row_overrides)
    return _config(tmp_path, rows)


# --------------------------------------------------------------------------
# F-1: the group-size cap reached only the discovery path
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n_treat", [MAX_DECLARED_GROUP_SIZE + 1, 999_999, 2**31])
def test_catalog_refuses_an_implausible_group_size(tmp_path, n_treat) -> None:
    """A group size above the documented cap was accepted from a workbook.

    `discovery_run.MAX_GROUP_SIZE` refused the same number in the browser's
    review panel, so a typing slip (3 becomes 999999) was refused in one entry
    point and silently accepted in the other, where it saturated the contrast's
    min(sqrt(n_ctrl + n_treat), 4) sample-size weight.
    """

    config = _two_units(tmp_path, n_treat=n_treat)
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config)
    message = str(excinfo.value)
    assert "n_treat" in message
    assert "10,000" in message


def test_catalog_accepts_a_group_size_at_the_cap(tmp_path) -> None:
    config = _two_units(tmp_path, n_treat=MAX_DECLARED_GROUP_SIZE)
    assert validate_catalog_inputs(config)["active_contrasts"] == 2


# --------------------------------------------------------------------------
# F-2: the direction-conflict flag named well-formed source units
# --------------------------------------------------------------------------

def _harmonized_for_direction(n_units: int, reversed_units: set[str]) -> pd.DataFrame:
    """One row per gene per unit, with a shared latent signature.

    Studies of one topic correlate genome-wide, which is what the coherence
    diagnostics measure; independent noise per unit leaves every correlation at
    zero and the flag can never fire either way.
    """

    rng = np.random.default_rng(17)
    genes = [f"G{index:03d}" for index in range(120)]
    latent = rng.normal(0.0, 1.0, len(genes))
    rows = []
    for index in range(n_units):
        unit = f"SU{index + 1}"
        sign = -1.0 if unit in reversed_units else 1.0
        lfc = sign * (latent + rng.normal(0.0, 0.25, len(genes)))
        for gene, value in zip(genes, lfc):
            rows.append(
                {
                    "gene_symbol": gene,
                    "study_id": f"{unit}_1",
                    "source_unit_id": unit,
                    "lfc": float(value),
                    "signed_z": float(value * 2.0),
                    "normalized_rank": 0.5,
                    "n_genes_in_study": len(genes),
                    "pvalue": 0.01,
                    "padj": np.nan,
                    "table_scope": "full_results",
                    "n_ctrl": 3,
                    "n_treat": 3,
                    "source_input_type": "author_deg_table",
                    "assay_type": "RNA-seq",
                    "pipeline": "DESeq2",
                }
            )
    return pd.DataFrame(rows)


def _diagnostics(harmonized: pd.DataFrame) -> pd.DataFrame:
    _scores, _evidence, metadata = degora_score_table(harmonized, min_studies=2)
    return pd.DataFrame(metadata["source_quality_diagnostics"])


def test_only_the_reversed_unit_is_flagged_in_a_three_source_corpus() -> None:
    """Three sources, one reversed: every unit used to be flagged.

    A unit with two pairwise comparisons satisfies "at least half conflict" on a
    single conflicting partner, so the two well-formed units were named as
    disagreeing with the corpus while their own median Spearman sat at 0.00.
    """

    diagnostics = _diagnostics(_harmonized_for_direction(3, {"SU3"}))
    flagged = set(
        diagnostics.loc[diagnostics["source_direction_conflict_flag"], "source_unit_id"]
    )
    assert flagged == {"SU3"}


def test_two_source_corpus_still_flags_both_units() -> None:
    """With two units DEGORA cannot tell which is reversed; both stay flagged."""

    diagnostics = _diagnostics(_harmonized_for_direction(2, {"SU2"}))
    flagged = set(
        diagnostics.loc[diagnostics["source_direction_conflict_flag"], "source_unit_id"]
    )
    assert flagged == {"SU1", "SU2"}


def test_five_source_corpus_with_two_reversed_units_flags_only_those_two() -> None:
    diagnostics = _diagnostics(_harmonized_for_direction(5, {"SU4", "SU5"}))
    flagged = set(
        diagnostics.loc[diagnostics["source_direction_conflict_flag"], "source_unit_id"]
    )
    assert flagged == {"SU4", "SU5"}


def test_direction_conflict_rule_documents_the_median_condition() -> None:
    assert "median pairwise correlation" in DIRECTION_CONFLICT_RULE


def test_a_split_corpus_is_described_as_split(tmp_path) -> None:
    """Four units, two reversed: the corpus is divided, not one source wrong."""

    diagnostics = _diagnostics(_harmonized_for_direction(4, {"SU3", "SU4"}))
    flagged = set(
        diagnostics.loc[diagnostics["source_direction_conflict_flag"], "source_unit_id"]
    )
    # A 2-2 split is genuinely undecidable: DEGORA cannot know which half holds
    # the intended convention, so it flags all four and says the corpus is
    # divided instead of accusing each unit of disagreeing with "the others".
    assert flagged == {"SU1", "SU2", "SU3", "SU4"}
    warnings = direction_conflict_warnings(diagnostics)
    # One warning about the corpus, not four that each contradict themselves.
    assert len(warnings) == 1
    assert "this corpus is split" in warnings[0]
    assert "4 of 4 source units" in warnings[0]
    assert "cannot tell which half" in warnings[0]


def test_a_reversed_majority_flags_the_minority() -> None:
    """Three units, two reversed: the majority convention wins, by definition.

    DEGORA has no way to know which convention the reader intended, so the unit
    that disagrees with the rest is the one named - as it would be if a single
    unit were reversed.
    """

    diagnostics = _diagnostics(_harmonized_for_direction(3, {"SU2", "SU3"}))
    flagged = set(
        diagnostics.loc[diagnostics["source_direction_conflict_flag"], "source_unit_id"]
    )
    assert flagged == {"SU1"}


# --------------------------------------------------------------------------
# F-3: lfc_scale accepted any string
# --------------------------------------------------------------------------

@pytest.mark.parametrize("value", ["log2", "LOG2", "log_2", "log2fc", " log2 "])
def test_lfc_scale_accepts_the_documented_vocabulary(value) -> None:
    assert normalize_lfc_scale(value) == "log2"


@pytest.mark.parametrize("value", ["", None, float("nan"), "nan", "none"])
def test_lfc_scale_treats_an_unset_cell_as_undeclared(value) -> None:
    assert normalize_lfc_scale(value) == ""


@pytest.mark.parametrize("value", ["log10", "ln", "linear", "log", "banana"])
def test_lfc_scale_refuses_anything_else(value) -> None:
    """`log10` and `linear` were silently equivalent to a blank cell.

    `table_scope` and `time_course_mode` both refuse an unknown value; the one
    field whose whole purpose is to record the reader's statement about their
    data did not.
    """

    with pytest.raises(ValueError, match="Unsupported lfc_scale"):
        normalize_lfc_scale(value)


@pytest.mark.parametrize("value", ["log10", "linear"])
def test_validate_refuses_an_unsupported_lfc_scale(tmp_path, value) -> None:
    config = _two_units(tmp_path, lfc_scale=value)
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config)
    assert "lfc_scale" in str(excinfo.value)


# --------------------------------------------------------------------------
# O-1: peak_mean recorded no selection
# --------------------------------------------------------------------------

def test_peak_mean_reports_what_it_selected() -> None:
    """peak_mean drops contrasts and wrote an empty selection report.

    It is the one mode whose selection is per gene, so which contrasts survived
    differs from row to row - the strongest claim on an audit record, not the
    weakest.
    """

    rng = np.random.default_rng(5)
    rows = []
    for unit in ("SU1", "SU2"):
        for index, duration in enumerate((1, 6, 24, 48)):
            for gene_index in range(40):
                rows.append(
                    {
                        "gene_symbol": f"G{gene_index:03d}",
                        "study_id": f"{unit}_{duration}h",
                        "source_unit_id": unit,
                        "duration_h": duration,
                        "time_course_mode": "peak_mean",
                        "lfc": float(rng.normal(0, 1)),
                        "signed_z": float(rng.normal(0, 1) * (index + 1)),
                        "normalized_rank": 0.5,
                    }
                )
    report = time_course_selection_report(pd.DataFrame(rows))
    assert {entry["source_unit_id"] for entry in report} == {"SU1", "SU2"}
    for entry in report:
        assert entry["time_course_mode"] == "peak_mean"
        assert entry["rows_after"] < entry["rows_before"]
        assert 0.0 < entry["row_retention"] < 1.0


def test_mean_selects_nothing_and_is_not_reported() -> None:
    rows = [
        {
            "gene_symbol": "G1",
            "study_id": "SU1_1h",
            "source_unit_id": "SU1",
            "duration_h": 1,
            "time_course_mode": "mean",
            "lfc": 1.0,
            "signed_z": 2.0,
            "normalized_rank": 0.5,
        }
    ]
    assert time_course_selection_report(pd.DataFrame(rows)) == []


# --------------------------------------------------------------------------
# O-2: rank_universe_size_declared lost the catalog's value
# --------------------------------------------------------------------------

def test_declared_rank_universe_survives_the_observed_row_fallback(tmp_path) -> None:
    """The column named `_declared` reported the clamped value.

    Both `_declared` and `_used` read 300 for a catalog that declared 5, so the
    pair that exists to be read against each other carried one number and the
    declaration survived only inside free text.
    """

    config = _two_units(tmp_path, table_scope="deg_only", rank_universe_size=5)
    run_slice(config, tmp_path / "out", tmp_path / "harmonized", min_studies=2)
    harmonized = pd.read_csv(tmp_path / "out" / "slice_harmonized.csv", low_memory=False)
    rows = harmonized[harmonized["study_id"].astype(str).eq("A_6h")]
    declared = pd.to_numeric(rows["rank_universe_size_declared"], errors="coerce").dropna()
    used = pd.to_numeric(rows["rank_universe_size_used"], errors="coerce").dropna()
    assert float(declared.iloc[0]) == 5.0
    assert float(used.iloc[0]) == float(len(GENES))
    assert rows["rank_universe_warning"].astype(str).str.contains("rank_universe_size=5").any()


# --------------------------------------------------------------------------
# O-4: validate said "config OK" for a corpus run then refused
# --------------------------------------------------------------------------

def test_validate_warns_when_source_units_share_no_identifier_space(tmp_path) -> None:
    """A symbols-vs-Ensembl config validated clean and failed at the end of a run.

    The preflight already opens and harmonizes every source table, so it holds
    both identifier sets. It warns rather than refuses: the run-time message is
    computed from the harmonized table and is the better of the two, and a run
    at min_studies=1 is legitimate.
    """

    _table(tmp_path / "sym.csv", 1)
    _table(tmp_path / "ens.csv", 2, genes=[f"ENSG{index:011d}" for index in range(100)])
    config = _config(
        tmp_path,
        [
            {
                "study_id": "A_6h",
                "source_unit_id": "UA",
                "source_path": str(tmp_path / "sym.csv"),
                "gene_column": "gene",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
            },
            {
                "study_id": "B_6h",
                "source_unit_id": "UB",
                "source_path": str(tmp_path / "ens.csv"),
                "gene_column": "gene",
                "lfc_column": "log2FoldChange",
                "p_column": "pvalue",
                "padj_column": "padj",
            },
        ],
    )
    result = validate_catalog_inputs(config)
    warnings = " ".join(result.get("warnings", []))
    assert "share a single gene identifier" in warnings
    assert "UA" in warnings and "UB" in warnings


def test_validate_is_quiet_when_the_units_share_identifiers(tmp_path) -> None:
    config = _two_units(tmp_path)
    result = validate_catalog_inputs(config)
    warnings = " ".join(result.get("warnings", []))
    assert "identifier" not in warnings


# --------------------------------------------------------------------------
# P2 (external audit): --inspection-budget was not a global cap
# --------------------------------------------------------------------------

def test_negative_inspection_budget_is_refused(tmp_path, capsys) -> None:
    """Legacy GEO refused a negative budget; the shared CLI option accepted it.

    Reported by the external v0.4.32 validation as P2, together with the two
    cases below: one root cause in
    `max(1, min(12, inspection_budget // selected_count))`.
    """

    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "discover",
                "hypoxia",
                "--species",
                "human",
                "--output-dir",
                str(tmp_path / "out"),
                "--inspection-budget",
                "-1",
            ]
        )
    assert excinfo.value.code == 2
    assert "non-negative" in capsys.readouterr().err


def test_zero_inspection_budget_keeps_a_valid_per_record_ceiling(tmp_path, monkeypatch) -> None:
    """Global zero skips inspection without violating the backend's 1..12 ceiling."""

    captured = _capture_prepare(monkeypatch, tmp_path)
    assert (
        main(
            [
                "discover",
                "hypoxia",
                "--species",
                "human",
                "--output-dir",
                str(tmp_path / "out"),
                "--select",
                "PMID:900001",
                "--inspection-budget",
                "0",
            ]
        )
        == 0
    )
    assert captured["inspection_budget"] == 0
    assert captured["max_files_per_record"] == 1


def test_a_budget_below_the_selection_count_is_refused(tmp_path, monkeypatch, capsys) -> None:
    """Budget 1 over two selections implied a total of 2, above the cap asked for."""

    _capture_prepare(monkeypatch, tmp_path, records=2)
    exit_code = main(
        [
            "discover",
            "hypoxia",
            "--species",
            "human",
            "--output-dir",
            str(tmp_path / "out"),
            "--select",
            "PMID:900001",
            "--select",
            "PMID:900002",
            "--inspection-budget",
            "1",
        ]
    )
    assert exit_code == 2
    assert "cannot cover" in capsys.readouterr().err


def _capture_prepare(monkeypatch, tmp_path, *, records: int = 1) -> dict:
    """Isolate the federated select path from the network."""

    import sys
    import types

    from degora import discovery_federated

    captured: dict = {}

    def fake_prepare_publication_records(
        selected,
        species,
        *,
        query,
        max_files_per_record,
        materialize_dir,
        force,
        inspection_budget=None,
        before_publish=None,
    ):
        captured.update(
            max_files_per_record=max_files_per_record,
            inspection_budget=inspection_budget,
            selected=len(selected),
        )
        return {
            "returned_records": len(selected),
            "studies": [],
            "exports": {
                "draft_catalog_csv": str(materialize_dir / "DEGORA_discovery_draft_catalog.csv")
            },
        }

    module = types.ModuleType("degora.discovery_prepare")
    module.prepare_publication_records = fake_prepare_publication_records
    monkeypatch.setitem(sys.modules, "degora.discovery_prepare", module)

    snapshot = {
        "records": [
            {
                "canonical_id": f"PMID:90000{index + 1}",
                "pubmed_ids": [f"90000{index + 1}"],
                "species": "human",
                "scientific_name": "Homo sapiens",
                "title": "Synthetic record",
            }
            for index in range(records)
        ]
    }
    monkeypatch.setattr(
        discovery_federated,
        "search_publications",
        lambda query, species, limit: snapshot,
    )
    monkeypatch.setattr(
        discovery_federated,
        "resolve_publication_records",
        lambda selected, species: selected,
    )
    return captured
