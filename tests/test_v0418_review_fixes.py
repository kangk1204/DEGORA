"""Regressions for what the v0.4.18 review found."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from degora.formula_safety import restore_formula_text_if_marked
from degora.slice_runner import DegoraConfigError, validate_catalog_inputs


def _write_pair(directory: Path, first: pd.DataFrame, second: pd.DataFrame, **first_row) -> Path:
    first.to_csv(directory / "t1.csv", index=False)
    second.to_csv(directory / "t2.csv", index=False)
    row = {
        "study_id": "S1", "source_unit_id": "U1", "source_path": "t1.csv",
        "gene_column": "gene", "lfc_column": "log2FoldChange", "p_column": "pvalue",
        "include_in_analysis": "yes",
    }
    row.update(first_row)
    other = dict(row, study_id="S2", source_unit_id="U2", source_path="t2.csv", lfc_column="log2FoldChange")
    other.pop("lfc_scale", None)
    config = directory / "cfg.csv"
    pd.DataFrame([row, other]).to_csv(config, index=False)
    return config


def _clean(n: int = 150, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"gene": [f"G{i}" for i in range(n)], "log2FoldChange": rng.normal(0, 1, n), "pvalue": [1e-4] * n})


def test_an_authors_table_with_no_sidecar_is_read_verbatim(tmp_path) -> None:
    """A raw table is not a DEGORA artifact, and its text means what it says.

    An unmapped notes column holding "'=see figure 2" refused the whole table
    and told the reader to restore a .provenance.json sidecar that DEGORA had
    never written. Ambiguity about a guard only exists for a file that carries a
    sidecar which does not vouch for it.
    """

    table = _clean()
    table["notes"] = ["'=see figure 2" if i == 3 else "ok" for i in range(len(table))]
    path = tmp_path / "author.csv"
    table.to_csv(path, index=False)

    frame = pd.read_csv(path)
    # As harmonize calls it: only the mapped columns can make the file ambiguous.
    restored = restore_formula_text_if_marked(frame, path, columns=["gene", "log2FoldChange", "pvalue"])

    assert restored["notes"].tolist() == frame["notes"].tolist()
    # A guard-like value in a column that will be used is still refused.
    with pytest.raises(ValueError, match="no matching DEGORA formula-guard provenance"):
        restore_formula_text_if_marked(frame, path, columns=["notes"])
    # And validate accepts the config built on it.
    config = _write_pair(tmp_path, table, _clean(seed=1))
    assert validate_catalog_inputs(config)["active_contrasts"] == 2


def test_a_sidecar_that_does_not_vouch_still_refuses_guarded_text(tmp_path) -> None:
    table = _clean()
    table["gene"] = ["'=HYPERLINK(\"x\")" if i == 0 else g for i, g in enumerate(table["gene"])]
    path = tmp_path / "out.csv"
    table.to_csv(path, index=False)
    path.with_suffix(".csv.provenance.json").write_text('{"artifact_sha256": "0", "metadata": {}}', encoding="utf-8")

    with pytest.raises(ValueError, match="no matching DEGORA formula-guard provenance"):
        restore_formula_text_if_marked(pd.read_csv(path), path)


def _up_only_log2(n: int = 300) -> pd.DataFrame:
    rng = np.random.default_rng(2)
    return pd.DataFrame({"gene": [f"G{i}" for i in range(n)], "lfc": np.abs(rng.normal(0, 1.5, n)) + 0.05, "pvalue": [1e-4] * n})


def test_an_up_only_log2_table_whose_header_says_nothing_is_refused_with_a_way_out(tmp_path) -> None:
    config = _write_pair(tmp_path, _up_only_log2(), _clean(seed=3), lfc_column="lfc")

    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config)

    message = str(excinfo.value)
    assert "linear fold change" in message
    assert "set lfc_scale to log2" in message


def test_declaring_lfc_scale_log2_turns_the_refusal_into_a_report(tmp_path) -> None:
    """degora init asks whether the values are log2; the answer has to reach validate."""

    config = _write_pair(tmp_path, _up_only_log2(), _clean(seed=3), lfc_column="lfc", lfc_scale="log2")

    result = validate_catalog_inputs(config)

    assert result["active_contrasts"] == 2
    assert any("lfc_scale=log2" in w and "reported, not refused" in w for w in result["warnings"])


def test_lfc_scale_is_case_and_whitespace_tolerant(tmp_path) -> None:
    config = _write_pair(tmp_path, _up_only_log2(), _clean(seed=3), lfc_column="lfc", lfc_scale=" LOG2 ")

    assert validate_catalog_inputs(config)["active_contrasts"] == 2


def test_init_records_the_log2_answer_it_asked_for(tmp_path) -> None:
    from degora.beginner import ContrastAnswers, catalog_row, infer_source_table

    path = tmp_path / "up.csv"
    _up_only_log2().to_csv(path, index=False)
    inference = infer_source_table(path)

    yes = catalog_row(inference, ContrastAnswers(positive_means_up_in_treated=True, effect_is_log2=True, source_unit_id="P1"),
                      study_id="S1", catalog_dir=tmp_path)
    no = catalog_row(inference, ContrastAnswers(positive_means_up_in_treated=True, effect_is_log2=False, source_unit_id="P1"),
                     study_id="S1", catalog_dir=tmp_path)

    assert yes["lfc_scale"] == "log2"
    assert yes["include_in_analysis"] == "yes"
    assert no["lfc_scale"] == ""
    assert no["include_in_analysis"] == "no"


def test_the_template_guide_documents_lfc_scale() -> None:
    from degora.excel_template import _guide_rows

    rows = _guide_rows()
    assert "lfc_scale" in set(rows["column"])


def test_discover_reports_its_stages_while_it_runs(monkeypatch, capsys, tmp_path) -> None:
    """Two live searches were silent for 25 seconds and then timed out.

    The provider layer reports every stage through a callback the browser
    already displays; the CLI simply never passed one.
    """

    import sys
    import types

    from degora.cli import main

    def search_publications(query, species, limit=1000, providers=None, progress=None):
        assert progress is not None, "the CLI must hand the provider a progress callback"
        progress(0.1, "Querying PubMed.")
        progress(0.6, "Resolving linked data.")
        return {"records": [], "total": 0, "provider_status": "complete", "provider_errors": [],
                "species": {"key": "human", "label": "Human"}, "query": query}

    module = types.ModuleType("degora.discovery_federated")
    module.search_publications = search_publications
    module.page_publication_snapshot = lambda *_a, **_k: {"records": [], "total": 0, "page": 1, "page_size": 10}
    module.resolve_publication_records = lambda *_a, **_k: []
    module.filter_publication_records = lambda records, text_filter="": list(records)
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)
    monkeypatch.setitem(sys.modules, "degora.discovery_export", types.SimpleNamespace(
        export_publication_search=lambda *_a, **_k: {"search_csv": str(tmp_path / "x.csv"), "search_json": str(tmp_path / "x.json"), "search_xlsx": str(tmp_path / "x.xlsx")}))

    main(["discover", "anything", "--species", "human", "--limit", "5", "--output-dir", str(tmp_path / "out")])

    err = capsys.readouterr().err
    assert "[ 10%] Querying PubMed." in err
    assert "[ 60%] Resolving linked data." in err


def test_inference_reads_at_most_the_cap_and_defers_the_scope(tmp_path) -> None:
    """A 203 MB interaction file beside the DEG tables stalled init for minutes.

    Only a prefix is read past the cap. A results table is usually sorted by
    p-value, so a prefix would call a full table DEG-only; the scope is left for
    the run, which reads everything.
    """

    from degora.beginner import INFERENCE_MAX_ROWS, describe_inference, infer_source_table

    n = INFERENCE_MAX_ROWS + 1000
    path = tmp_path / "huge.csv"
    pd.DataFrame(
        {"gene": [f"G{i}" for i in range(n)], "log2FoldChange": [1.0] * n, "pvalue": [1e-4] * n}
    ).to_csv(path, index=False)

    inference = infer_source_table(path)

    assert inference.sampled
    assert inference.n_rows == INFERENCE_MAX_ROWS
    assert inference.table_scope == "auto"
    assert "only the first" in inference.table_scope_reason
    assert describe_inference(inference)[0].startswith(f"huge.csv: {INFERENCE_MAX_ROWS:,}+ rows")
    # Column choices are unaffected by the cap.
    assert inference.mapping["gene_column"] == "gene"


def test_a_table_under_the_cap_is_not_marked_sampled(tmp_path) -> None:
    from degora.beginner import infer_source_table

    path = tmp_path / "small.csv"
    pd.DataFrame({"gene": ["A", "B", "C"], "log2FoldChange": [1.0, -1.0, 0.5], "pvalue": [0.01, 0.02, 0.5]}).to_csv(path, index=False)

    inference = infer_source_table(path)

    assert not inference.sampled
    assert inference.n_rows == 3


def test_only_header_like_rows_are_tried_when_locating_a_workbook_table(tmp_path) -> None:
    """Every one of the first ten rows used to cost a full read of the sheet."""

    from degora.beginner import _candidate_header_rows, _read_workbook_raw

    path = tmp_path / "titled.xlsx"
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame([["Supplementary Table 2"], [None], ["gene", "log2FoldChange", "pvalue"], ["TP53", 1.2, 0.01]]).to_excel(
            writer, sheet_name="Table S2", header=False, index=False
        )

    assert _candidate_header_rows(_read_workbook_raw(path)["Table S2"]) == [3]


def test_a_workbook_with_a_readme_sheet_is_still_located(tmp_path) -> None:
    from degora.beginner import infer_source_table

    path = tmp_path / "supp.xlsx"
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame({"README": ["This workbook holds the DEG table on the next sheet."]}).to_excel(writer, sheet_name="README", index=False)
        pd.DataFrame([["Table S1"], ["gene", "log2FoldChange", "pvalue"]] + [[f"G{i}", 1.0, 0.001] for i in range(120)]).to_excel(
            writer, sheet_name="S1", header=False, index=False
        )

    inference = infer_source_table(path)

    assert inference.looks_like_a_deg_table
    assert inference.sheet_name == "S1"
    assert inference.header_row == 2
    assert inference.n_rows == 120


def test_init_names_a_large_file_before_reading_it(tmp_path, monkeypatch) -> None:
    from degora import beginner

    deg = tmp_path / "deg"
    deg.mkdir()
    big = deg / "big.csv"
    pd.DataFrame({"gene": ["A"], "log2FoldChange": [1.0], "pvalue": [0.01]}).to_csv(big, index=False)
    real_stat = Path.stat

    def fake_stat(self, *args, **kwargs):
        result = real_stat(self, *args, **kwargs)
        if self.name == "big.csv":
            return type("S", (), {"st_size": 50_000_000, "st_mode": result.st_mode})()
        return result

    monkeypatch.setattr(Path, "stat", fake_stat)
    lines: list[str] = []
    answers = iter(["human", "unsure", "yes"])
    with pytest.raises(Exception):
        beginner.run_init(tmp_path / "cfg.csv", deg, ask=lambda q, default="": next(answers), echo=lines.append)

    assert any(line.startswith("Reading big.csv (50 MB)") for line in lines)


def test_a_semicolon_delimited_csv_is_read_and_a_comma_in_a_description_does_not_break_it(tmp_path) -> None:
    """R's write.csv2 output is a `.csv` that is semicolon-delimited.

    Read with a comma, a gene description containing one made the row ragged
    and the file failed with "Error tokenizing data" - naming neither the cause
    nor the `sep` column that fixes it.
    """

    from degora.beginner import infer_source_table

    path = tmp_path / "write_csv2.csv"
    path.write_text(
        '"gene";"description";"log2FoldChange";"pvalue"\n'
        + "".join(f'"G{i}";"a protein, part {i}";{1.0 - i * 0.01};0.001\n' for i in range(120)),
        encoding="utf-8",
    )

    inference = infer_source_table(path)

    assert inference.readable
    assert list(inference.columns) == ["gene", "description", "log2FoldChange", "pvalue"]
    assert inference.mapping["gene_column"] == "gene"


def test_a_byte_order_mark_does_not_become_part_of_the_first_header(tmp_path) -> None:
    from degora.beginner import infer_source_table

    path = tmp_path / "bom.csv"
    path.write_bytes(
        b"\xef\xbb\xbfGENEID,SYMBOL,log2FoldChange,pvalue\n"
        + b"".join(f"{7000 + i},GENE{i},{1.0 - i * 0.01},0.001\n".encode() for i in range(120))
    )

    inference = infer_source_table(path)

    assert list(inference.columns)[0] == "GENEID"
    assert inference.mapping["gene_column"] == "SYMBOL"


def test_the_symbol_column_is_preferred_over_a_descriptive_gene_name_column(tmp_path) -> None:
    """GENEID, GENENAME, SYMBOL side by side offered GENENAME: 'SPARC like 1'."""

    from degora.beginner import infer_source_table

    path = tmp_path / "three_names.csv"
    pd.DataFrame(
        {
            "GENEID": [8404 + i for i in range(120)],
            "GENENAME": [f"protein like {i}" for i in range(120)],
            "SYMBOL": [f"SPARCL{i}" for i in range(120)],
            "logFC": [1.0 - i * 0.01 for i in range(120)],
            "PValue": [0.001] * 120,
        }
    ).to_csv(path, index=False)

    inference = infer_source_table(path)

    assert inference.mapping["gene_column"] == "SYMBOL"
    assert inference.identifier_space == "gene symbol"


def test_a_workbook_is_opened_once_during_inference(tmp_path, monkeypatch) -> None:
    """openpyxl parses the shared-string table on every open; the search opened it 24 times."""

    from degora import beginner

    path = tmp_path / "supp.xlsx"
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame({"README": ["see next sheet"]}).to_excel(writer, sheet_name="README", index=False)
        pd.DataFrame([["Table S1"], ["gene", "log2FoldChange", "pvalue"]] + [[f"G{i}", 1.0, 0.001] for i in range(120)]).to_excel(
            writer, sheet_name="S1", header=False, index=False
        )
    opens = []
    real = pd.ExcelFile.__init__

    def counting(self, *args, **kwargs):
        opens.append(1)
        return real(self, *args, **kwargs)

    monkeypatch.setattr(pd.ExcelFile, "__init__", counting)

    inference = beginner.infer_source_table(path)

    assert inference.looks_like_a_deg_table and inference.sheet_name == "S1" and inference.header_row == 2
    assert len(opens) == 1


def _abundance_only_table(tmp_path) -> Path:
    """Counts, FPKM, an FDR and a label - a results file with no effect size."""

    path = tmp_path / "abundance_only.csv"
    pd.DataFrame(
        {
            "#Symbol": [f"G{i}" for i in range(150)],
            "CP1_count": [100 + i for i in range(150)],
            "T1_count": [90 + i for i in range(150)],
            "CP1_FPKM": [10.5 + i for i in range(150)],
            "pvalue": [0.0005 * (i + 1) for i in range(150)],
            "FDR": [0.001 * (i + 1) for i in range(150)],
            "regulated": ["up"] * 150,
        }
    ).to_csv(path, index=False)
    return path


def test_abundance_columns_are_never_offered_as_the_effect_size(tmp_path) -> None:
    """CP1_count was offered as the effect size and taken; a count is not a contrast."""

    from degora.beginner import infer_source_table

    inference = infer_source_table(_abundance_only_table(tmp_path))

    assert inference.plausible.get("lfc_column", ()) == ()


def test_a_table_with_no_effect_size_column_is_skipped_with_that_reason(tmp_path) -> None:
    from degora.beginner import run_init

    deg = tmp_path / "deg"
    deg.mkdir()
    _abundance_only_table(deg)
    lines: list[str] = []
    with pytest.raises(Exception):
        run_init(tmp_path / "cfg.csv", deg, ask=lambda q, default="": "human", echo=lines.append)

    assert any("no effect-size column" in line for line in lines)


def test_init_writes_lfc_scale_into_the_config(tmp_path) -> None:
    from degora.beginner import run_init

    deg = tmp_path / "deg"
    deg.mkdir()
    _up_only_log2().to_csv(deg / "up.csv", index=False)
    answers = iter(["human", "lfc", "yes", "yes", "a vs b", "P1", "3", "3"])
    run_init(tmp_path / "cfg.csv", deg, ask=lambda q, default="": next(answers), echo=lambda _l: None)

    config = pd.read_csv(tmp_path / "cfg.csv")
    assert "lfc_scale" in config.columns
    assert config["lfc_scale"].tolist() == ["log2"]


def test_a_count_column_mapped_as_the_effect_size_is_refused_not_warned_about(tmp_path) -> None:
    """A column of read counts reached scoring with two warnings and exit 0."""

    rng = np.random.default_rng(4)
    counts = pd.DataFrame(
        {"gene": [f"G{i}" for i in range(300)], "CP1_count": rng.integers(0, 2_000_000, 300), "pvalue": [1e-4] * 300}
    )
    config = _write_pair(tmp_path, counts, _clean(seed=5), lfc_column="CP1_count")

    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config)

    assert "not a log2 fold change" in str(excinfo.value)


def test_a_handful_of_absurd_values_is_still_only_a_warning(tmp_path) -> None:
    table = _clean(n=300, seed=6)
    table.loc[0, "log2FoldChange"] = 45.0  # one outlier in three hundred
    config = _write_pair(tmp_path, table, _clean(seed=7))

    result = validate_catalog_inputs(config)

    assert result["active_contrasts"] == 2
    assert any("|log2FC| > 30" in w for w in result["warnings"])


def test_an_entrez_id_column_is_not_offered_as_the_effect_size(tmp_path) -> None:
    """EntrezID was offered as the effect size and taken; validate then refused the run."""

    from degora.beginner import infer_source_table

    path = tmp_path / "entrez.csv"
    pd.DataFrame(
        {
            "Symbol": [f"G{i}" for i in range(200)],
            "EntrezID": [7000 + i * 37 for i in range(200)],
            "FC": [1.5 + i / 100 for i in range(200)],
            "pvalue": [0.001] * 200,
        }
    ).to_csv(path, index=False)

    inference = infer_source_table(path)

    assert "EntrezID" not in inference.plausible["lfc_column"]
    assert "EntrezID" in inference.plausible["gene_column"]
    assert "FC" in inference.plausible["lfc_column"]


@pytest.mark.parametrize("header", ["log2(fc)", "log2 (FC)", "log2.fc", "log2_fc", "Log2(Fold Change)"])
def test_log2_fc_with_any_separator_is_the_effect_column(header: str) -> None:
    """A real table carried `fc` and `log2(fc)` side by side; neither was recognised."""

    from degora.discovery import classify_header

    mapping = classify_header(["gene", "fc", header, "pval"])["mapping"]

    assert mapping["lfc_column"] == header


def test_a_coordinate_column_is_not_offered_as_the_effect_size(tmp_path) -> None:
    """GeneLeft - a base-pair position - was offered on an interaction table."""

    from degora.beginner import infer_source_table

    path = tmp_path / "coords.csv"
    pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(200)],
            "GeneLeft": [1_000_000 + i * 5_000 for i in range(200)],
            "GeneRight": [1_004_000 + i * 5_000 for i in range(200)],
            "log2FoldChange": [1.0 - i / 100 for i in range(200)],
            "pvalue": [0.001] * 200,
        }
    ).to_csv(path, index=False)

    inference = infer_source_table(path)

    assert "GeneLeft" not in inference.plausible["lfc_column"]
    assert "GeneRight" not in inference.plausible["lfc_column"]
    assert "log2FoldChange" in inference.plausible["lfc_column"]


def test_a_query_with_no_english_term_is_refused_before_any_request() -> None:
    """NCBI drops a Korean term silently; the organism filter alone then matches all of GEO."""

    from degora.discovery import DiscoveryError, _query_terms

    with pytest.raises(DiscoveryError, match="no English term"):
        _query_terms("태반")


def test_non_english_terms_are_set_aside_and_named_when_mixed() -> None:
    from degora.discovery import _query_terms, ignored_query_terms

    assert _query_terms("태반 placenta hypoxia") == ["placenta", "hypoxia"]
    assert ignored_query_terms("태반 placenta hypoxia") == ["태반"]


def test_greek_letters_inside_english_terms_are_kept() -> None:
    from degora.discovery import _query_terms, ignored_query_terms

    assert _query_terms("TGF-β EMT") == ["TGF-β", "EMT"]
    assert ignored_query_terms("α-synuclein") == []


def test_discover_names_ignored_terms_and_a_partial_snapshot(monkeypatch, capsys, tmp_path) -> None:
    import sys
    import types

    from degora.cli import main

    def search_publications(query, species, limit=1000, providers=None, progress=None):
        return {
            "records": [], "total": 0, "query": query,
            "species": {"key": "human", "label": "Human"},
            "provider_status": "partial",
            "diagnostics": {"errors": ["ncbi_pubmed: HTTP 429 rate limited"]},
            "ignored_terms": ["태반"],
        }

    module = types.ModuleType("degora.discovery_federated")
    module.search_publications = search_publications
    module.page_publication_snapshot = lambda *_a, **_k: {"records": [], "total": 0, "page": 1, "page_size": 10}
    module.resolve_publication_records = lambda *_a, **_k: []
    module.filter_publication_records = lambda records, text_filter="": list(records)
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)
    monkeypatch.setitem(sys.modules, "degora.discovery_export", types.SimpleNamespace(
        export_publication_search=lambda *_a, **_k: {"search_csv": str(tmp_path / "x.csv"), "search_json": str(tmp_path / "x.json"), "search_xlsx": str(tmp_path / "x.xlsx")}))

    main(["discover", "태반 placenta", "--species", "human", "--limit", "5", "--output-dir", str(tmp_path / "out")])

    err = capsys.readouterr().err
    assert "Ignored (not English): 태반" in err
    assert "WARNING: this snapshot is partial" in err
    assert "HTTP 429" in err


def test_a_single_option_prompt_takes_enter_as_that_option(tmp_path) -> None:
    """The p prompt offered "padj" alone with no default; Enter counted as a failure."""

    from degora.beginner import run_init

    deg = tmp_path / "deg"
    deg.mkdir()
    pd.DataFrame(
        {
            "row_name": [f"G{i}" for i in range(120)],
            "baseMean": [100.0 + i for i in range(120)],
            "log2 fold change": [2.0 - i * 0.01 for i in range(120)],
            "padj": [i / 200.0 for i in range(120)],
        }
    ).to_csv(deg / "padj_only.tsv", sep="\t", index=False)
    other = deg / "other.csv"
    _clean(seed=8).to_csv(other, index=False)
    # Enter (blank) for the p column, then the direction and metadata.
    answers = iter(["human", "yes", "a vs b", "P1", "3", "3", "", "yes", "a vs b", "P2", "3", "3"])
    lines: list[str] = []
    summary = run_init(tmp_path / "cfg.csv", deg, ask=lambda q, default="": next(answers), echo=lines.append)

    assert summary["n_contrasts"] == 2
    assert not any("cannot be used" in line for line in lines)
    config = pd.read_csv(tmp_path / "cfg.csv")
    assert config.set_index("study_id").loc["padj_only", "p_column"] == "padj"


def test_zero_records_come_with_a_next_step(monkeypatch, capsys, tmp_path) -> None:
    import sys
    import types

    from degora.cli import main

    def search_publications(query, species, limit=1000, providers=None, progress=None):
        return {"records": [], "total": 0, "query": query, "species": {"key": "human", "label": "Human"},
                "provider_status": "complete", "diagnostics": {"errors": []}, "ignored_terms": []}

    module = types.ModuleType("degora.discovery_federated")
    module.search_publications = search_publications
    module.page_publication_snapshot = lambda *_a, **_k: {"records": [], "total": 0, "page": 1, "page_size": 10}
    module.resolve_publication_records = lambda *_a, **_k: []
    module.filter_publication_records = lambda records, text_filter="": list(records)
    monkeypatch.setitem(sys.modules, "degora.discovery_federated", module)
    monkeypatch.setitem(sys.modules, "degora.discovery_export", types.SimpleNamespace(
        export_publication_search=lambda *_a, **_k: {"search_csv": str(tmp_path / "x.csv"), "search_json": str(tmp_path / "x.json"), "search_xlsx": str(tmp_path / "x.xlsx")}))

    main(["discover", "qqzzxx nonexistent syndrome", "--species", "human", "--limit", "5", "--output-dir", str(tmp_path / "out")])

    assert "No records matched. Try a broader" in capsys.readouterr().err


def test_an_existing_output_folder_says_what_to_do(tmp_path, capsys) -> None:
    from degora.cli import main

    out = tmp_path / "out"
    out.mkdir()
    (out / "publication_search.csv").write_text("x\n", encoding="utf-8")

    assert main(["discover", "hypoxia", "--species", "human", "--limit", "5", "--output-dir", str(out)]) == 2
    err = capsys.readouterr().err
    assert "pass a new --output-dir" in err
    assert "Discover tab" in err


def test_prepare_summary_names_why_a_record_has_nothing_ready() -> None:
    from degora.cli import _prepare_record_lines

    studies = [
        {"accession": "GSE343715", "ready_for_review_count": 0, "upstream_matrix_count": 3,
         "files": [{"inspection": {"status": "upstream_matrix", "reason": "matrix columns were detected; choose at least two control and two treatment samples"}}]},
        {"accession": "GSE302293", "ready_for_review_count": 0, "upstream_matrix_count": 0,
         "files": [{"inspection": {"status": "requires_pvalue_mapping", "reason": "adjusted significance detected but nominal p-value is missing"}}]},
        {"accession": "GSE1", "ready_for_review_count": 2, "files": []},
        {"accession": "GSE2", "ready_for_review_count": 0, "files": []},
    ]

    lines = _prepare_record_lines(studies)

    assert lines[0].startswith("GSE343715: only expression matrices were found")
    assert lines[1].startswith("GSE302293: a results table with adjusted p-values only")
    assert lines[2].startswith("GSE2: no supplementary files")
    assert len(lines) == 3  # the record with tables ready is not listed


def test_a_guarded_value_in_a_workbook_is_treated_like_one_in_a_csv(tmp_path) -> None:
    """'=ISG15 in an xlsx was re-guarded to ''=ISG15 on every round trip, silently."""

    from degora.harmonize import TableMapping, read_deg_table

    path = tmp_path / "guarded.xlsx"
    pd.DataFrame({"gene": ["'=ISG15", "TP53"], "log2FoldChange": [1.0, 2.0], "pvalue": [0.01, 0.02]}).to_excel(path, index=False)

    # No sidecar and a guard-like value in the gene column: refused, as the CSV path is.
    with pytest.raises(ValueError, match="formula-guard"):
        read_deg_table(path, TableMapping("gene", "log2FoldChange", "pvalue"))
    # The same value in a column that is not used is left alone.
    pd.DataFrame({"gene": ["ISG15"], "notes": ["'=see fig 2"], "log2FoldChange": [1.0], "pvalue": [0.01]}).to_excel(path, index=False)
    frame = read_deg_table(path, TableMapping("gene", "log2FoldChange", "pvalue"))
    assert frame["notes"].tolist() == ["'=see fig 2"]


def test_the_prepare_cap_counts_selected_publications_not_expanded_series(monkeypatch) -> None:
    """20 selected publications linked 22 GEO series; every repository record failed.

    The browser caps the selection at 20 publications. The repository phase
    re-applied the same cap to the series those publications link to, so a
    full selection in which one publication carries three series failed the
    whole repository half with "at most 20 studies can be prepared at once".
    """

    from degora.discovery import DiscoveryError, _validated_accessions

    twenty_two = [f"GSE{100000 + i}" for i in range(22)]
    with pytest.raises(DiscoveryError, match="at most 20"):
        _validated_accessions(twenty_two)
    assert len(_validated_accessions(twenty_two, max_studies=22)) == 22


def test_federated_prepare_hands_the_repository_phase_the_expanded_series_count(monkeypatch) -> None:
    import degora.discovery_prepare as prepare_module

    seen: dict = {}

    def fake_prepare_geo_studies(accessions, species, **kwargs):
        seen["n"] = len(list(accessions))
        seen["max_studies"] = kwargs.get("max_studies")
        return {"studies": [], "excluded_studies": []}

    monkeypatch.setattr(prepare_module, "prepare_geo_studies", fake_prepare_geo_studies)
    source = open(prepare_module.__file__, encoding="utf-8").read()
    assert "max_studies=max(len(accessions), 1)" in source
    # Behaviour is pinned at the call site; the kwarg reaches the cap check.
    from degora.discovery import _validated_accessions

    assert len(_validated_accessions([f"GSE{1 + i}" for i in range(25)], max_studies=25)) == 25


def test_the_readiness_badge_carries_one_phrase_and_the_basis_sits_under_it() -> None:
    from degora.api import INDEX_HTML

    assert "${esc(`${headline} · ${seen}`)}" not in INDEX_HTML
    assert '<span class="dataset-title readiness-basis">${esc(seen)}</span>' in INDEX_HTML
    assert ".readiness-basis { display: block;" in INDEX_HTML


def test_the_species_decision_is_said_once_in_words() -> None:
    from degora.api import INDEX_HTML

    assert 'target_species_verified: "species verified"' in INDEX_HTML
    assert "`species ${String(study.species_decision).replaceAll" not in INDEX_HTML


def test_the_inspect_cell_does_not_clip_its_button() -> None:
    from degora.api import INDEX_HTML

    assert '<td class="inspect-cell"><button class="action-secondary study-inspect"' in INDEX_HTML
    assert "td.inspect-cell { overflow: visible; text-overflow: clip; }" in INDEX_HTML


def test_sort_indicators_are_arrows_not_letters() -> None:
    from degora.api import INDEX_HTML

    assert '"asc" ? "^" : "v"' not in INDEX_HTML
    assert "\u25B4" in INDEX_HTML and "\u25BE" in INDEX_HTML


def test_the_inspector_recognises_an_r_export_whose_gene_names_are_row_labels() -> None:
    """Six DESeq2 result files the guided setup had scored were rejected as not_deg_table."""

    from degora.discovery import _inspect_rows

    rows = [["", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]] + [
        [f"ENSG{i:011d}.1", 100.0 + i, 0.5 - i / 100, 0.1, 1.2, 0.01 + i / 1000, 0.05] for i in range(30)
    ]

    header = _inspect_rows(rows)

    assert header["mapping"]["gene_column"] == "row_name"
    assert header["mapping"]["lfc_column"] == "log2FoldChange"
    assert header["status"] == "ready_for_review"


def test_a_named_first_column_is_left_alone_by_the_row_label_repair() -> None:
    from degora.discovery import _name_row_label_column

    assert _name_row_label_column(["gene", "lfc"], [["A", 1.0], ["B", 2.0]]) == ["gene", "lfc"]
    # An empty header over a numeric column is not a label column.
    assert _name_row_label_column(["", "lfc"], [[1, 1.0], [2, 2.0]]) == ["", "lfc"]
