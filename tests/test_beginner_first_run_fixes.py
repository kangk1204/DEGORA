"""Regressions for the issues a first-time reader hit following the README on v0.4.17.

Each test names the behaviour a beginner actually ran into, not the internal
function that changed, so a future refactor has to keep the user-visible promise
rather than the current implementation.
"""

from __future__ import annotations

import json
import sqlite3
import stat
import subprocess
import threading
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import pytest

from degora.api import ScoreDatabaseError, create_server, serve
from degora.discovery import describe_transport_error, network_failure_message
from degora.excel_export import _curated_lookup, _read_gold_from_config
from degora.harmonize import TableMapping, canonical_gene_symbol, harmonize_frame
from degora.slice_runner import DegoraConfigError, read_catalog, run_slice, validate_catalog_inputs

ROOT = Path(__file__).resolve().parents[1]

# One gene per pair: what a source table might carry, and what DEGORA scores it as.
LEGACY_TO_CURRENT = {
    "SEPT9": "SEPTIN9",
    "9-Sep": "SEPTIN9",
    "MARCH1": "MARCHF1",
    "1-Mar": "MARCHF1",
    "DEC1": "BHLHE40",
    "1-Dec": "BHLHE40",
}


# --------------------------------------------------------------------------
# F-01  Gene symbols resolve the same way everywhere, and leave a trail
# --------------------------------------------------------------------------


@pytest.mark.parametrize(("written", "scored"), sorted(LEGACY_TO_CURRENT.items()))
def test_legacy_and_date_damaged_symbols_resolve_to_the_scored_symbol(written: str, scored: str) -> None:
    assert canonical_gene_symbol(written) == scored


@pytest.mark.parametrize("symbol", ["SEPTIN9", "MARCHF1", "BHLHE40", "ISG15", "TP53"])
def test_current_symbols_are_left_alone(symbol: str) -> None:
    assert canonical_gene_symbol(symbol) == symbol
    assert canonical_gene_symbol(symbol.lower()) == symbol


@pytest.mark.parametrize("value", ["", "   ", None, float("nan"), pd.NA])
def test_unusable_gene_labels_resolve_to_nothing(value: object) -> None:
    assert canonical_gene_symbol(value) == ""


def _deg_table(genes: list[str], *, offset: float = 0.0) -> pd.DataFrame:
    # `offset` keeps two tables with the same gene list from being byte-identical:
    # one result table declared as two independent source units is refused.
    size = len(genes)
    return pd.DataFrame(
        {
            "gene": genes,
            "log2FoldChange": [2.0 + offset + index * 0.1 for index in range(size)],
            "pvalue": [10.0 ** -(6 + index) for index in range(size)],
            "padj": [10.0 ** -(5 + index) for index in range(size)],
        }
    )


def _write_two_source_config(
    tmp_path: Path,
    *,
    first_genes: list[str],
    second_genes: list[str],
    gold: pd.DataFrame | None = None,
) -> Path:
    """Write a runnable two-source-unit Excel config, optionally with a GoldPanel."""

    deg_dir = tmp_path / "deg_tables"
    deg_dir.mkdir(exist_ok=True)
    _deg_table(first_genes).to_csv(deg_dir / "study_a.csv", index=False)
    _deg_table(second_genes, offset=0.05).to_csv(deg_dir / "study_b.csv", index=False)
    contrasts = pd.DataFrame(
        {
            "study_id": ["STUDY_A", "STUDY_B"],
            "source_unit_id": ["PMID_A", "PMID_B"],
            "source_path": ["deg_tables/study_a.csv", "deg_tables/study_b.csv"],
            "gene_column": ["gene", "gene"],
            "lfc_column": ["log2FoldChange"] * 2,
            "p_column": ["pvalue"] * 2,
            "padj_column": ["padj"] * 2,
            "table_scope": ["full_results"] * 2,
            "species": ["Homo sapiens"] * 2,
            "include": ["yes", "yes"],
        }
    )
    config_path = tmp_path / "config.xlsx"
    with pd.ExcelWriter(config_path, engine="openpyxl") as writer:
        contrasts.to_excel(writer, sheet_name="Contrasts", index=False)
        if gold is not None:
            gold.to_excel(writer, sheet_name="GoldPanel", index=False)
    return config_path


def _run(config_path: Path, tmp_path: Path, *, excel: bool = True) -> dict:
    """Run the pipeline the way a reader does, through the CLI, and return its metrics."""

    from degora.cli import main

    argv = [
        "run",
        str(config_path),
        "--output-dir",
        str(tmp_path / "out"),
        "--harmonized-dir",
        str(tmp_path / "harmonized"),
    ]
    if not excel:
        argv.append("--no-excel")
    assert main(argv) == 0
    return json.loads((tmp_path / "out" / "slice_metrics.json").read_text(encoding="utf-8"))


def test_one_gene_written_three_ways_is_scored_once(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "MARCH1", "DEC1", "ISG15"],
        second_genes=["9-Sep", "1-Mar", "1-Dec", "ISG15"],
    )
    metrics = _run(config_path, tmp_path)

    scores = pd.read_csv(tmp_path / "out" / "degora_gene_scores.csv")
    assert set(scores["gene_symbol"]) == {"SEPTIN9", "MARCHF1", "BHLHE40", "ISG15"}
    # Both spellings counted as one gene seen by two independent source units.
    assert set(scores["n_source_units"]) == {2}
    assert metrics["n_gene_symbols_normalized"] == 3


def test_the_label_each_source_carried_survives_into_the_harmonized_table(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "ISG15"],
        second_genes=["9-Sep", "ISG15"],
    )
    _run(config_path, tmp_path)

    harmonized = pd.read_csv(tmp_path / "out" / "slice_harmonized.csv")
    assert "input_gene_label" in harmonized.columns
    repaired = harmonized.loc[harmonized["gene_symbol"].eq("SEPTIN9")]
    labels = dict(zip(repaired["study_id"], repaired["input_gene_label"]))
    assert labels == {"STUDY_A": "SEPT9", "STUDY_B": "9-Sep"}
    # A gene that needed no repair still records the label it arrived with.
    untouched = harmonized.loc[harmonized["gene_symbol"].eq("ISG15"), "input_gene_label"]
    assert set(untouched) == {"ISG15"}


def test_every_label_a_collapsed_gene_came_from_is_kept(tmp_path) -> None:
    # One table spelling the same gene two ways must not lose the other spelling.
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "9-Sep", "ISG15"],
        second_genes=["SEPTIN9", "ISG15"],
    )
    _run(config_path, tmp_path)

    harmonized = pd.read_csv(tmp_path / "out" / "slice_harmonized.csv")
    first = harmonized.loc[
        harmonized["study_id"].eq("STUDY_A") & harmonized["gene_symbol"].eq("SEPTIN9"), "input_gene_label"
    ].iloc[0]
    assert set(str(first).split(";")) == {"SEPT9", "9-Sep"}


def _gold(symbols: list[str], *, locked: str = "yes") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": symbols,
            "expected_direction": ["up"] * len(symbols),
            "role": ["optional_marker"] * len(symbols),
            "evidence_basis": ["regression"] * len(symbols),
            "locked": [locked] * len(symbols),
        }
    )


def test_a_gold_panel_written_with_legacy_symbols_finds_its_genes(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "MARCH1", "DEC1", "ISG15"],
        second_genes=["SEPT9", "MARCH1", "DEC1", "ISG15"],
        gold=_gold(["SEPT9", "MARCH1", "DEC1", "ISG15"]),
    )
    _run(config_path, tmp_path)

    gold, status, reason = _read_gold_from_config(config_path)
    assert status == "locked", reason
    scores = pd.read_csv(tmp_path / "out" / "degora_gene_scores.csv")
    lookup = _curated_lookup(gold, scores)

    assert list(lookup["present_in_degora_output"]) == [True] * 4
    resolved = dict(zip(lookup["gene_symbol"], lookup["resolved_gene_symbol"]))
    assert resolved == {
        "SEPT9": "SEPTIN9",
        "MARCH1": "MARCHF1",
        "DEC1": "BHLHE40",
        "ISG15": "ISG15",
    }


def test_a_gold_panel_gene_absent_from_the_run_is_still_reported_absent(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["ISG15", "IFIT1"],
        second_genes=["ISG15", "IFIT1"],
        gold=_gold(["ISG15", "NOT_A_REAL_GENE"]),
    )
    _run(config_path, tmp_path)

    gold, status, _ = _read_gold_from_config(config_path)
    assert status == "locked"
    scores = pd.read_csv(tmp_path / "out" / "degora_gene_scores.csv")
    lookup = _curated_lookup(gold, scores).set_index("gene_symbol")
    assert bool(lookup.loc["ISG15", "present_in_degora_output"]) is True
    assert bool(lookup.loc["NOT_A_REAL_GENE", "present_in_degora_output"]) is False


def test_run_metrics_gold_panel_resolves_legacy_symbols(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "ISG15"],
        second_genes=["9-Sep", "ISG15"],
        gold=_gold(["SEPT9"]),
    )
    metrics = _run(config_path, tmp_path)
    assert metrics["gold_panel_status"] == "locked"
    assert metrics["gold_panel_gene_count"] == 1
    # The panel said SEPT9; recall is computed against the symbol the run scored.
    assert metrics["recall_at_50"]["recovered"] == ["SEPTIN9"]
    assert metrics["recall_at_50"]["missing"] == []


# --------------------------------------------------------------------------
# F-08  A GoldPanel nobody locked says so
# --------------------------------------------------------------------------


def test_a_gold_panel_with_nothing_locked_reports_why(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["ISG15", "IFIT1"],
        second_genes=["ISG15", "IFIT1"],
        gold=_gold(["ISG15", "IFIT1"], locked="no"),
    )
    metrics = _run(config_path, tmp_path)

    assert metrics["gold_panel_status"] == "not_provided"
    assert "none are locked" in metrics["gold_panel_reason"]
    assert metrics["n_panel_rows"] == 2

    from degora.cli import _run_warning_messages

    assert any("none are locked" in message for message in _run_warning_messages(metrics))

    _, status, reason = _read_gold_from_config(config_path)
    assert status == "not_provided"
    assert "none are locked" in reason


def test_an_empty_gold_panel_is_not_reported_as_a_mistake(tmp_path) -> None:
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["ISG15", "IFIT1"],
        second_genes=["ISG15", "IFIT1"],
    )
    metrics = _run(config_path, tmp_path)
    assert metrics["n_panel_rows"] == 0
    from degora.cli import _run_warning_messages

    assert not any("locked" in message for message in _run_warning_messages(metrics))


def test_the_template_documents_locked_where_a_reader_will_look(tmp_path) -> None:
    from degora.excel_template import write_template

    path = tmp_path / "template.xlsx"
    write_template(path)
    guide = pd.read_excel(path, sheet_name="ColumnGuide", header=1)
    gold_rows = guide.loc[guide["checked_where"].astype(str).str.contains("GoldPanel")]
    assert "locked" in set(gold_rows["column"])
    locked_meaning = gold_rows.loc[gold_rows["column"].eq("locked"), "meaning"].iloc[0]
    assert "no drops it" in str(locked_meaning)

    note = pd.read_excel(path, sheet_name="GoldPanel", header=None).iloc[0, 0]
    assert "locked" in str(note)


# --------------------------------------------------------------------------
# F-06  Duplicate rows that disagree about direction say so
# --------------------------------------------------------------------------


def test_harmonize_warns_when_duplicate_rows_disagree_on_direction() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX", "RPL13A"],
            "log2FoldChange": [2.0, -3.0, 0.1],
            "pvalue": [1e-5, 1e-4, 0.9],
            "padj": [0.001, 0.01, 0.9],
        }
    )
    meta = {"study_id": "CONFLICT", "paper_id": "P1", "pipeline": "DESeq2", "assay_type": "RNA-seq"}
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), meta)

    warning = out.loc[out["gene_symbol"].eq("GENEX"), "gene_symbol_collapse_warning"].iloc[0]
    assert "opposite log2 fold-change signs" in warning
    assert "1 gene(s)" in warning
    # The collapse itself is unchanged: the most significant row still wins.
    assert out.loc[out["gene_symbol"].eq("GENEX"), "lfc"].iloc[0] == pytest.approx(2.0)


def test_duplicate_rows_that_agree_on_direction_stay_quiet() -> None:
    frame = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX", "RPL13A"],
            "log2FoldChange": [2.0, 3.0, 0.1],
            "pvalue": [1e-5, 1e-4, 0.9],
            "padj": [0.001, 0.01, 0.9],
        }
    )
    meta = {"study_id": "AGREE", "paper_id": "P1", "pipeline": "DESeq2", "assay_type": "RNA-seq"}
    out = harmonize_frame(frame, TableMapping("gene", "log2FoldChange", "pvalue", "padj"), meta)
    assert out.loc[out["gene_symbol"].eq("GENEX"), "gene_symbol_collapse_warning"].iloc[0] == ""


# --------------------------------------------------------------------------
# F-02  serve refuses anything that is not a score database
# --------------------------------------------------------------------------


def test_serve_refuses_a_deg_table_instead_of_binding(tmp_path) -> None:
    table = tmp_path / "edger.tsv"
    table.write_text("gene\tlogFC\tPValue\nISG15\t2.0\t0.01\n", encoding="utf-8")
    with pytest.raises(ScoreDatabaseError) as excinfo:
        serve(table, quiet=True)
    message = str(excinfo.value)
    assert "not a DEGORA score database" in message
    assert "degora run" in message


def test_serve_refuses_a_sqlite_file_without_the_score_tables(tmp_path) -> None:
    other = tmp_path / "notes.db"
    with sqlite3.connect(other) as connection:
        connection.execute("CREATE TABLE notes (body TEXT)")
    with pytest.raises(ScoreDatabaseError, match="missing genes, gene_evidence, studies, meta"):
        serve(other, quiet=True)


def test_serve_refuses_a_directory(tmp_path) -> None:
    with pytest.raises(ScoreDatabaseError, match="is a directory"):
        serve(tmp_path, quiet=True)


def test_serve_still_reports_a_missing_file_as_missing(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        serve(tmp_path / "absent.db", quiet=True)


def test_cli_prints_the_serve_error_without_a_traceback(tmp_path, capsys) -> None:
    from degora.cli import main

    table = tmp_path / "edger.tsv"
    table.write_text("gene\tlogFC\tPValue\nISG15\t2.0\t0.01\n", encoding="utf-8")
    assert main(["serve", str(table)]) == 2
    captured = capsys.readouterr()
    assert "not a DEGORA score database" in captured.err
    assert "Traceback" not in captured.err


# --------------------------------------------------------------------------
# F-01  The browser and API answer for the symbol the reader knows
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def served_run(tmp_path_factory) -> dict:
    tmp_path = tmp_path_factory.mktemp("served")
    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "MARCH1", "DEC1", "ISG15"],
        second_genes=["9-Sep", "1-Mar", "1-Dec", "ISG15"],
    )
    _run(config_path, tmp_path, excel=False)
    server = create_server(tmp_path / "out" / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    try:
        yield {"base_url": f"http://{host}:{port}", "tmp_path": tmp_path}
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def _get(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=5) as response:
        return json.loads(response.read().decode("utf-8"))


@pytest.mark.parametrize(("query", "expected"), sorted(LEGACY_TO_CURRENT.items()))
def test_gene_search_finds_the_gene_under_the_name_the_reader_typed(served_run, query, expected) -> None:
    payload = _get(f"{served_run['base_url']}/api/genes?q={urllib.parse.quote(query)}")
    assert [gene["gene_symbol"] for gene in payload["genes"]] == [expected]


def test_partial_gene_search_still_works(served_run) -> None:
    payload = _get(f"{served_run['base_url']}/api/genes?q=SEPT")
    assert [gene["gene_symbol"] for gene in payload["genes"]] == ["SEPTIN9"]


def test_gene_detail_resolves_a_legacy_symbol_and_says_so(served_run) -> None:
    payload = _get(f"{served_run['base_url']}/api/genes/SEPT9")
    assert payload["gene"]["gene_symbol"] == "SEPTIN9"
    assert payload["requested_gene_symbol"] == "SEPT9"
    assert payload["resolved_gene_symbol"] == "SEPTIN9"
    assert payload["evidence"], "evidence must be looked up under the resolved symbol"


def test_gene_detail_for_the_current_symbol_is_unchanged(served_run) -> None:
    payload = _get(f"{served_run['base_url']}/api/genes/SEPTIN9")
    assert payload["gene"]["gene_symbol"] == "SEPTIN9"
    assert "resolved_gene_symbol" not in payload


def test_an_unknown_gene_is_still_not_found(served_run) -> None:
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        _get(f"{served_run['base_url']}/api/genes/NOT_A_REAL_GENE")
    assert excinfo.value.code == 404


# --------------------------------------------------------------------------
# F-03  Every missing source file is reported at once
# --------------------------------------------------------------------------


def _catalog_with_missing_rows(tmp_path: Path) -> Path:
    deg_dir = tmp_path / "deg_tables"
    deg_dir.mkdir()
    _deg_table(["ISG15", "IFIT1"]).to_csv(deg_dir / "real.csv", index=False)
    config_path = tmp_path / "config.csv"
    pd.DataFrame(
        {
            "study_id": ["EXAMPLE_ONE", "EXAMPLE_TWO", "MINE"],
            "source_unit_id": ["EX1", "EX2", "MINE"],
            "source_path": ["missing_one.csv", "missing_two.csv", "deg_tables/real.csv"],
            "gene_column": ["gene"] * 3,
            "lfc_column": ["log2FoldChange"] * 3,
            "p_column": ["pvalue"] * 3,
            "include": ["yes"] * 3,
        }
    ).to_csv(config_path, index=False)
    return config_path


def test_validation_reports_every_missing_source_file_not_just_the_first(tmp_path) -> None:
    config_path = _catalog_with_missing_rows(tmp_path)
    with pytest.raises(DegoraConfigError) as excinfo:
        validate_catalog_inputs(config_path)
    message = str(excinfo.value)
    assert "EXAMPLE_ONE" in message
    assert "EXAMPLE_TWO" in message
    assert "MINE" not in message  # the row that is fine is not blamed


def test_a_run_reports_every_missing_source_file_too(tmp_path) -> None:
    config_path = _catalog_with_missing_rows(tmp_path)
    with pytest.raises(DegoraConfigError) as excinfo:
        run_slice(config_path, tmp_path / "out", tmp_path / "harmonized", 2)
    message = str(excinfo.value)
    assert "EXAMPLE_ONE" in message
    assert "EXAMPLE_TWO" in message


# --------------------------------------------------------------------------
# F-04  The missing-column error does not contradict itself
# --------------------------------------------------------------------------


def test_missing_column_error_lists_only_columns_the_file_had(tmp_path) -> None:
    deg_table = tmp_path / "deseq2.csv"
    _deg_table(["ISG15", "IFIT1"]).to_csv(deg_table, index=False)

    with pytest.raises(DegoraConfigError) as excinfo:
        read_catalog(deg_table)
    message = str(excinfo.value)

    context = message.split("Problems:")[0]
    # The error says source_unit_id is missing; it must not also list it as available.
    assert "source_unit_id" not in context
    assert "hypoxia_modality" not in context
    assert "gene" in context
    assert "Missing required Contrasts column 'source_unit_id (or paper_id)'" in message


def test_passing_a_deg_table_as_the_config_says_what_went_wrong(tmp_path) -> None:
    deg_table = tmp_path / "deseq2.csv"
    _deg_table(["ISG15", "IFIT1"]).to_csv(deg_table, index=False)
    with pytest.raises(DegoraConfigError) as excinfo:
        read_catalog(deg_table)
    message = str(excinfo.value)
    assert "looks like a DEG results table" in message
    assert "degora template" in message


def test_a_real_config_missing_one_column_does_not_get_the_deg_table_hint(tmp_path) -> None:
    config_path = tmp_path / "config.csv"
    pd.DataFrame(
        {
            "study_id": ["A"],
            "source_unit_id": ["U"],
            "source_path": ["a.csv"],
            "lfc_column": ["log2FoldChange"],
            "p_column": ["pvalue"],
        }
    ).to_csv(config_path, index=False)
    with pytest.raises(DegoraConfigError) as excinfo:
        read_catalog(config_path)
    assert "looks like a DEG results table" not in str(excinfo.value)


# --------------------------------------------------------------------------
# F-05  A search that cannot reach the network says why
# --------------------------------------------------------------------------


def test_network_failure_message_names_the_cause_and_the_fix() -> None:
    message = network_failure_message("NCBI could not be reached after 3 attempt(s)", ["the request timed out"])
    assert "Problems:" in message
    assert "- the request timed out" in message
    assert "How to fix:" in message
    assert "proxy" in message
    assert "work offline" in message


def test_transport_errors_are_described_not_named() -> None:
    assert describe_transport_error(TimeoutError()) == "the request timed out"
    described = describe_transport_error(urllib.error.URLError("unreachable"))
    assert "connection could not be established" in described
    assert "URLError" not in described
    assert describe_transport_error(urllib.error.HTTPError("u", 503, "x", {}, None)) == "the server answered HTTP 503"


def test_a_blocked_search_reports_the_provider_reasons() -> None:
    from degora import discovery_federated
    from degora.discovery import DiscoveryUnavailableError

    class DeadProvider:
        name = "ncbi_pubmed"

        def search(self, *args, **kwargs):
            raise DiscoveryUnavailableError("request failed after 3 attempt(s): the request timed out")

    with pytest.raises(DiscoveryUnavailableError) as excinfo:
        discovery_federated.search_publications("hypoxia", "human", providers=[DeadProvider()])
    message = str(excinfo.value)
    assert message.startswith("no publication search provider could be reached")
    assert "the request timed out" in message
    assert "How to fix:" in message


# --------------------------------------------------------------------------
# F-09  The dashboard opens on what the reader just built
# --------------------------------------------------------------------------


def test_dashboard_chooses_its_landing_view_from_the_scored_gene_count(served_run) -> None:
    with urllib.request.urlopen(served_run["base_url"], timeout=5) as response:
        html = response.read().decode("utf-8")
    assert 'showView(Number.isFinite(scored) && scored > 0 ? "atlas" : "discover")' in html
    meta = _get(f"{served_run['base_url']}/api/meta")["meta"]
    assert int(meta["n_gene_scores"]) > 0


# --------------------------------------------------------------------------
# F-11  The run says which rank is the headline one
# --------------------------------------------------------------------------


def test_run_summary_names_the_primary_rank_column(tmp_path, capsys) -> None:
    from degora.cli import main

    config_path = _write_two_source_config(
        tmp_path,
        first_genes=["SEPT9", "ISG15", "IFIT1"],
        second_genes=["9-Sep", "ISG15", "IFIT1"],
    )
    assert main(["run", str(config_path), "--output-dir", str(tmp_path / "out"), "--no-excel"]) == 0
    out = capsys.readouterr().out
    assert "- Primary rank column: quality_weighted_degora_rank" in out
    assert "- Gene symbols normalized: 1" in out
    assert "input_gene_label" in out


# --------------------------------------------------------------------------
# F-07  The quickstart command in the release notes works from a checkout
# --------------------------------------------------------------------------


def test_repository_root_quickstart_forwarder_exists_and_is_executable() -> None:
    forwarder = ROOT / "degora_quickstart.sh"
    assert forwarder.is_file()
    assert forwarder.stat().st_mode & stat.S_IXUSR


def test_repository_root_quickstart_forwarder_forwards() -> None:
    result = subprocess.run(
        ["bash", "degora_quickstart.sh", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert "Inside a checkout" in result.stdout


def test_the_forwarder_says_where_to_get_the_real_script_when_used_alone(tmp_path) -> None:
    solo = tmp_path / "degora_quickstart.sh"
    solo.write_bytes((ROOT / "degora_quickstart.sh").read_bytes())
    result = subprocess.run(
        ["bash", str(solo)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 1
    assert "scripts/degora_quickstart.sh" in result.stderr
    assert "raw.githubusercontent.com" in result.stderr


# --------------------------------------------------------------------------
# Nothing secret is committed
# --------------------------------------------------------------------------


def test_credential_scratch_files_are_ignored() -> None:
    ignored = (ROOT / ".gitignore").read_text(encoding="utf-8")
    for pattern in ("gitignore_git.txt", "*token*.txt", "*.token", ".env"):
        assert pattern in ignored, f"{pattern} must stay out of the repository"
