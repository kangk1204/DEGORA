"""Regressions for the v0.4.23 fixes: fan-out, identifier spaces, counts, CLI guards."""

from __future__ import annotations

from pathlib import Path


import pandas as pd
import pytest


def test_a_publication_linking_many_series_is_set_aside_with_the_count(tmp_path, monkeypatch) -> None:
    """20 selected publications became 69 series; one consortium paper linked 51."""

    import degora.discovery_prepare as prep
    from degora.discovery import normalize_species

    captured: dict = {}

    def fake_prepare_geo_studies(accessions, species, **kwargs):
        captured["accessions"] = list(accessions)
        return {"studies": [], "excluded_studies": []}

    monkeypatch.setattr(prep, "prepare_geo_studies", fake_prepare_geo_studies)
    mega = {"canonical_id": "pmid:22955616", "title": "consortium", "geo_accessions": [f"GSE{100 + i}" for i in range(8)]}
    normal = {"canonical_id": "accession:GSE343715", "title": "one series", "geo_accessions": ["GSE343715"]}
    excluded: list = []

    prep._prepare_into_staging(
        geo_records=[mega, normal], direct_records=[], spec=normalize_species("human"), query="placenta",
        staging=tmp_path, max_files_per_record=3, transport=None, geo_client=None, excluded=excluded,
    )

    assert captured["accessions"] == ["GSE343715"]
    reasons = [str(item.get("reason") or item.get("exclusion_reason") or item) for item in excluded]
    assert any("links 8 GEO series" in reason and "select them directly" in reason for reason in reasons)


def test_the_inspector_prefers_the_column_whose_values_are_symbols_and_names_the_space() -> None:
    """A column called "Gene Symbol" holding Ensembl IDs joins nothing written in symbols."""

    from degora.discovery import _inspect_rows

    rows = [["GENEID", "GENENAME", "SYMBOL", "logFC", "PValue"]] + [
        [8404 + i, f"protein like {i}", f"SPARCL{i}", 1.0 - i / 100, 0.001] for i in range(30)
    ]
    header = _inspect_rows(rows)
    assert header["mapping"]["gene_column"] == "SYMBOL"
    assert header["gene_identifier_space"] == "gene symbol"

    ensembl_only = [["Gene Symbol", "log2FoldChange", "pvalue"]] + [[f"ENSG{i:011d}", 1.0, 0.001] for i in range(30)]
    header = _inspect_rows(ensembl_only)
    assert header["gene_identifier_space"] == "Ensembl ID"


def test_a_source_unit_in_another_identifier_space_is_flagged_even_with_partial_overlap() -> None:
    """869 genes, all ENSG, reported as a success: overlap alone missed it."""

    from degora.slice_runner import _identifier_space_warnings

    symbols = [f"GENE{i}" for i in range(400)] + [f"ENSG{i:011d}" for i in range(20)]
    ensembl = [f"ENSG{i:011d}" for i in range(450)]
    harmonized = pd.DataFrame(
        {"gene_symbol": symbols + ensembl, "source_unit_id": ["U_symbols"] * len(symbols) + ["U_ensembl"] * len(ensembl)}
    )

    warnings = _identifier_space_warnings(harmonized, min_studies=2)

    assert any("written in Ensembl ID" in w and "gene symbol" in w for w in warnings)


def test_two_symbol_units_with_good_overlap_raise_no_space_warning() -> None:
    from degora.slice_runner import _identifier_space_warnings

    genes = [f"GENE{i}" for i in range(300)]
    harmonized = pd.DataFrame({"gene_symbol": genes + genes[:250], "source_unit_id": ["U1"] * 300 + ["U2"] * 250})

    assert _identifier_space_warnings(harmonized, min_studies=2) == []


def test_the_cli_refuses_a_one_character_query_and_page_zero(tmp_path, capsys) -> None:
    from degora.cli import main

    assert main(["discover", "x", "--species", "human", "--output-dir", str(tmp_path / "a")]) == 2
    assert "at least 2 characters" in capsys.readouterr().err
    with pytest.raises(SystemExit) as excinfo:
        main(["discover", "hypoxia", "--species", "human", "--page", "0", "--output-dir", str(tmp_path / "b")])
    assert excinfo.value.code == 2
    assert "must be 1 or more" in capsys.readouterr().err


def test_serve_without_a_database_says_how_to_get_one(tmp_path) -> None:
    from degora.api import ScoreDatabaseError, _require_degora_score_database

    with pytest.raises(ScoreDatabaseError) as excinfo:
        _require_degora_score_database(tmp_path / "missing.db")

    assert "No run yet?" in str(excinfo.value)
    assert "degora demo" in str(excinfo.value)


def test_a_unit_mixing_identifier_spaces_is_reported_with_the_split() -> None:
    """70% symbols and 30% Ensembl was called "gene symbol" and the 30% joined nothing."""

    from degora.slice_runner import _identifier_space_warnings

    mixed = [f"GENE{i}" for i in range(140)] + [f"ENSG{i:011d}" for i in range(60)]
    other = [f"GENE{i}" for i in range(150)]
    harmonized = pd.DataFrame({"gene_symbol": mixed + other, "source_unit_id": ["U1"] * 200 + ["U2"] * 150})

    warnings = _identifier_space_warnings(harmonized, min_studies=2)

    assert any("mixes identifier spaces" in w and "30%" in w and "Ensembl ID" in w for w in warnings)


def test_a_few_ensembl_fallbacks_in_a_large_symbol_unit_do_not_flip_its_space() -> None:
    """Sampling the head of the sorted list saw only ENSG... and called the unit Ensembl."""

    from degora.slice_runner import _identifier_space_warnings

    big = [f"GENE{i}" for i in range(5000)] + [f"ENSG{i:011d}" for i in range(100)]
    other = [f"GENE{i}" for i in range(4000)]
    harmonized = pd.DataFrame({"gene_symbol": big + other, "source_unit_id": ["U1"] * 5100 + ["U2"] * 4000})

    warnings = _identifier_space_warnings(harmonized, min_studies=2)

    assert not any("written in Ensembl ID" in w for w in warnings)
    assert not any("mixes identifier spaces" in w for w in warnings)  # 2% is below the reporting share


def test_the_zip_member_filter_admits_every_format_the_readers_read() -> None:
    from degora.discovery_prepare import _TABULAR_MEMBER_RE

    for name in ("DEG_results.xlsx.gz", "DEG_results.xls.gz", "DEG_results.xls", "table.csv.gz", "table.tsv"):
        assert _TABULAR_MEMBER_RE.search(name), name
    for name in ("GSE1_RAW.tar", "track.bw", "reads.fastq.gz"):
        assert not _TABULAR_MEMBER_RE.search(name), name


def test_serve_on_a_missing_database_reaches_the_first_run_hint(tmp_path, capsys) -> None:
    from degora.cli import main

    assert main(["serve", str(tmp_path / "missing.db")]) == 2
    err = capsys.readouterr().err
    assert "No run yet?" in err and "degora demo" in err


def test_the_candidate_panels_span_the_grid_and_unusable_studies_are_grouped_last() -> None:
    """Panels fell into the 28px checkbox column; unusable studies were interleaved."""

    from degora.api import INDEX_HTML

    assert ".candidate-row > .candidate-advanced,\n    .candidate-row > .candidate-confirms," in INDEX_HTML
    assert "grid-column: 1 / -1;" in INDEX_HTML
    assert ".candidate-advanced > summary::-webkit-details-marker { display: none; }" in INDEX_HTML
    assert "const analyzable = allStudies.filter((study) => (study.files || []).some(eligibleCandidate));" in INDEX_HTML
    assert 'with no usable table' in INDEX_HTML
    # The analyzable list renders before the grouped remainder.
    assert INDEX_HTML.index("analyzable.map(renderStudy)") < INDEX_HTML.index("unanalyzable.map(renderStudy)")


def _fallback_selection(label, control, treatment, candidate="c1", accession="GSE343715", **extra):
    entry = {
        "accession": accession, "candidate_id": candidate, "mode": "fallback", "contrast_label": label,
        "control_samples": control, "treatment_samples": treatment, "direction_confirmed": True,
        "biological_replicates_confirmed": True, "matrix_type": "normalized_expression_matrix", "normalized_scale": "log2",
    }
    entry.update(extra)
    return entry


def test_the_same_sample_groups_selected_from_two_files_of_one_series_are_refused() -> None:
    """Raw counts, TMM and FPKM of one series are one experiment; all three passed silently."""

    from degora.discovery_run import DiscoveryError, _check_fallback_selection_consistency

    same = [_fallback_selection("FPKM", ["c1", "c2"], ["t1", "t2"], candidate="fpkm"),
            _fallback_selection("TMM", ["c1", "c2"], ["t1", "t2"], candidate="tmm")]
    with pytest.raises(DiscoveryError, match="selected twice") as excinfo:
        _check_fallback_selection_consistency(same)
    assert "raw counts preferred" in str(excinfo.value)


def test_swapped_groups_from_one_series_are_refused_as_contradictory() -> None:
    from degora.discovery_run import DiscoveryError, _check_fallback_selection_consistency

    swapped = [_fallback_selection("A vs C", ["c1", "c2"], ["t1", "t2"], candidate="a"),
               _fallback_selection("C vs A", ["t1", "t2"], ["c1", "c2"], candidate="b")]
    with pytest.raises(DiscoveryError, match="control and treatment swapped"):
        _check_fallback_selection_consistency(swapped)


def test_a_sample_in_both_groups_is_left_to_the_per_selection_check() -> None:
    """The activation step already refuses it in its own words; the guard stays out of the way."""

    from degora.discovery_run import _check_fallback_selection_consistency

    assert _check_fallback_selection_consistency([_fallback_selection("x", ["c1", "s"], ["s", "t1"])]) == []


def test_the_guard_groups_selections_by_the_bundle_series_not_by_sample_names() -> None:
    """ctrl_1/treat_1 from two different studies looked like one experiment selected twice."""

    from degora.discovery_run import _check_fallback_selection_consistency

    prepared = {"studies": [
        {"accession": "GSE1", "files": [{"candidate_id": "a"}]},
        {"accession": "GSE2", "files": [{"candidate_id": "b"}]},
    ]}
    entries = [_fallback_selection("x", ["ctrl_1", "ctrl_2"], ["treat_1", "treat_2"], candidate="a", accession=""),
               _fallback_selection("y", ["ctrl_1", "ctrl_2"], ["treat_1", "treat_2"], candidate="b", accession="")]
    assert _check_fallback_selection_consistency(entries, prepared) == []


def test_option_spellings_labels_and_whitespace_combinations(tmp_path) -> None:
    """Every way a reader can mistype a selection, and what each one gets."""

    from degora.discovery import normalize_species
    from degora.discovery_run import DiscoveryError, _fallback_row

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    path = bundle / "GSE9_matrix.csv"
    pd.DataFrame({"gene": [f"G{i}" for i in range(60)], "c1": [5.0 + i * 0.1 for i in range(60)], "c2": [5.2 + i * 0.1 for i in range(60)],
                  "t1": [7.0 + i * 0.1 for i in range(60)], "t2": [7.3 + i * 0.1 for i in range(60)]}).to_csv(path, index=False)
    import hashlib as _h
    candidate = {"candidate_id": "m", "name": path.name, "role": "unknown_matrix",
                 "inspection": {"status": "upstream_matrix_ready_for_contrast", "fetch_scope": "full", "local_path": str(path),
                                "full_file_sha256": _h.sha256(path.read_bytes()).hexdigest(), "sample_columns": ["c1", "c2", "t1", "t2"],
                                "gene_column": "gene", "header_row": 1}}
    study = {"accession": "GSE9", "title": "t", "files": [candidate]}

    def attempt(**overrides):
        entry = {"candidate_id": "m", "mode": "fallback", "direction_confirmed": True, "biological_replicates_confirmed": True,
                 "control_samples": ["c1", "c2"], "treatment_samples": ["t1", "t2"], "matrix_type": "normalized_expression_matrix",
                 "normalized_scale": "log2", "gene_column": "gene", "contrast_label": "treated vs control"}
        entry.update(overrides)
        try:
            _fallback_row(study=study, candidate=candidate, entry=entry, spec=normalize_species("human"), bundle_root=bundle,
                          derived_dir=tmp_path / f"d{abs(hash(str(overrides)))}", sequence=1, replay_command="degora")
            return "ok"
        except DiscoveryError as exc:
            return str(exc)

    # Tolerated: case and stray whitespace in enumerations, names and labels.
    assert attempt(matrix_type="Normalized-Expression-Matrix") == "ok"
    assert attempt(matrix_type=" NORMALIZED_EXPRESSION_MATRIX ") == "ok"
    assert attempt(normalized_scale="Log2") == "ok"
    assert attempt(normalized_scale=" LOG2 ") == "ok"
    assert attempt(mode="Fallback") in ("ok",)  # mode is normalised by the caller; the row itself does not read it
    assert attempt(control_samples=["c1 ", " c2"], treatment_samples=["t1", "t2 "]) == "ok"
    assert attempt(contrast_label="  Treated VS Control  ") == "ok"
    # Refused, each with the field named.
    assert "matrix_type" in attempt(matrix_type="normalized_expression")
    assert "normalized_scale" in attempt(normalized_scale="log10")
    assert "not found in the inspected matrix" in attempt(control_samples=["C1", "c2"])  # column names are case-sensitive in the file
    assert "contrast_label" in attempt(contrast_label="")
    assert "direction_confirmed" in attempt(direction_confirmed="true")  # a string is not a confirmation
    assert "biological_replicates_confirmed" in attempt(biological_replicates_confirmed=None)
    assert "both" in attempt(control_samples=["c1", "t1"], treatment_samples=["t1", "t2"]).lower() or "disjoint" in attempt(control_samples=["c1", "t1"], treatment_samples=["t1", "t2"])


def test_a_sample_that_changes_role_between_contrasts_is_a_warning_not_a_refusal() -> None:
    """Valid for a time series (T1 vs T0, T2 vs T1); worth a sentence, not a stop."""

    from degora.discovery_run import _check_fallback_selection_consistency

    warnings = _check_fallback_selection_consistency([
        _fallback_selection("T1 vs T0", ["t0a", "t0b"], ["t1a", "t1b"], candidate="a"),
        _fallback_selection("T2 vs T1", ["t1a", "t1b"], ["t2a", "t2b"], candidate="b"),
    ])
    assert len(warnings) == 1 and "t1a, t1b" in warnings[0] and "time series" in warnings[0]


def test_distinct_contrasts_from_different_series_raise_nothing() -> None:
    from degora.discovery_run import _check_fallback_selection_consistency

    assert _check_fallback_selection_consistency([
        _fallback_selection("A", ["c1", "c2"], ["t1", "t2"], candidate="a", accession="GSE1"),
        _fallback_selection("B", ["c1", "c2"], ["t1", "t2"], candidate="b", accession="GSE2"),
    ]) == []


def test_a_linear_matrix_declared_log2_is_refused_before_derivation(tmp_path) -> None:
    from degora.discovery_run import DiscoveryError, _require_plausible_scale

    path = tmp_path / "fpkm.csv"
    pd.DataFrame({"gene": [f"G{i}" for i in range(300)], "c1": [1500.0 + i for i in range(300)], "t1": [2200.0 + i for i in range(300)]}).to_csv(path, index=False)
    with pytest.raises(DiscoveryError, match="looks like a linear matrix"):
        _require_plausible_scale(path, ["c1", "t1"], "log2")
    # The same file declared linear is fine; a log2 file declared linear is not.
    _require_plausible_scale(path, ["c1", "t1"], "linear")
    log2 = tmp_path / "log2.csv"
    pd.DataFrame({"gene": ["A", "B", "C"], "c1": [-1.2, 3.0, 8.5], "t1": [0.5, 2.0, 9.0]}).to_csv(log2, index=False)
    with pytest.raises(DiscoveryError, match="negative values"):
        _require_plausible_scale(log2, ["c1", "t1"], "linear")


def test_the_preferred_file_of_a_series_is_the_least_processed_evidence() -> None:
    """Not the most frequent - a series ships one of each - but the best evidence first."""

    from degora.discovery import annotate_candidate_preference

    def matrix(name):
        return {"name": name, "inspection": {"status": "upstream_matrix_ready_for_contrast", "sample_columns": ["a", "b", "c", "d"]}}

    files = annotate_candidate_preference([
        matrix("GSE343715_Normalized_FPKM_gene_counts_matrix.txt.gz"),
        matrix("GSE343715_Normalized_LOG2_TMM_gene_counts_matrix.txt.gz"),
        matrix("GSE343715_SALMON_tx2gene_counts_matrix.txt.gz"),
    ])
    ranks = {f["name"].split("_", 1)[1][:12]: f["preference_rank"] for f in files}
    assert ranks["Normalized_F"] == 4 and ranks["Normalized_L"] == 3 and ranks["SALMON_tx2ge"] == 2
    assert [f["name"] for f in files if f.get("preferred")] == ["GSE343715_SALMON_tx2gene_counts_matrix.txt.gz"]

    with_author = annotate_candidate_preference([
        matrix("GSE1_counts.txt.gz"),
        {"name": "GSE1_DEG.xlsx", "inspection": {"status": "ready_for_review"}},
    ])
    assert [f["name"] for f in with_author if f.get("preferred")] == ["GSE1_DEG.xlsx"]


def test_the_browser_shows_the_preferred_file_first_and_collapses_the_rest() -> None:
    from degora.api import INDEX_HTML

    assert "function candidatePreferenceRank(candidate)" in INDEX_HTML
    assert ".sort((a, b) => a.rank - b.rank || a.index - b.index)" in INDEX_HTML
    assert '<details class="alternative-candidates">' in INDEX_HTML
    assert "open only if the file above is not the one to use" in INDEX_HTML


def test_search_estimates_the_likely_input_from_file_names_in_the_documented_order() -> None:
    from degora.discovery_federated import likely_input

    assert likely_input({"candidates": [{"name": "GSE2_DESeq2_results_ATRA_vs_ctrl.csv.gz"}]}) == (0, "author DEG table")
    assert likely_input({"candidates": [{"name": "GSE3_raw_counts.txt.gz"}]}) == (2, "raw count matrix")
    assert likely_input({"candidates": [{"name": "GSE4_log2_TMM_matrix.txt.gz"}]}) == (3, "log2 normalised matrix")
    assert likely_input({"candidates": [{"name": "GSE1_RAW.tar"}, {"name": "GSE1_Normalized_FPKM_gene_counts_matrix.txt.gz"}]}) == (4, "normalised matrix")
    assert likely_input({"candidates": [{"name": "GSE5_RAW.tar"}]}) == (9, "no tabular file seen")


def test_within_a_readiness_tier_a_record_naming_a_deg_table_sorts_before_a_count_matrix() -> None:
    from degora.discovery_federated import rank_publication_records

    counts = {"canonical_id": "accession:GSE1", "relevance_rank": 1, "year": 2026,
              "data_readiness": {"priority": 1, "verification_state": "verified_ready", "likely_input_rank": 2, "likely_input": "raw count matrix"}}
    deg = {"canonical_id": "accession:GSE2", "relevance_rank": 2, "year": 2026,
           "data_readiness": {"priority": 1, "verification_state": "verified_ready", "likely_input_rank": 0, "likely_input": "author DEG table"}}
    ranked = rank_publication_records([counts, deg], "data_readiness", "desc")

    assert [row["canonical_id"] for row in ranked] == ["accession:GSE2", "accession:GSE1"]


def test_the_readiness_line_and_readme_state_the_preference_order() -> None:
    from degora.api import INDEX_HTML

    assert "`likely ${detail.likely_input} · `" in INDEX_HTML
    readme = Path("README.md").read_text(encoding="utf-8")
    assert "### What DEGORA prefers, in order" in readme
    for phrase in ("The authors' own results table, covering every gene tested", "A raw count matrix", "A linear normalised matrix"):
        assert phrase in readme
