"""0.4.29: the prepared-evidence card, exercised end to end in node.

The browser's prepare step had ended in "Preparation failed: studies is not
defined" since 0.4.24: a rename inside renderPreparedState left three uses of
the old name, and no test ever executed that function. The card still drew,
so a screenshot looked fine, but the status line, the eligibility check and
the Run button never updated. These tests run the page script under a
permissive fake DOM and call the renderer on a real prepared bundle.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

from degora.api import INDEX_HTML

FIXTURE = Path(__file__).parent / "fixtures" / "prepared_bundle_two_studies.json"
HAS_NODE = subprocess.run(["which", "node"], capture_output=True).returncode == 0


def _script() -> str:
    return re.search(r"<script>(.*?)</script>", INDEX_HTML, re.S).group(1)


def test_page_text_carries_no_control_characters_from_python_escapes() -> None:
    # `content: "\25B8 "` in the page source is an octal escape to Python, not
    # a CSS one: the disclosure triangle rendered as "\x15B8" -> "B8".
    stray = sorted({hex(ord(ch)) for ch in INDEX_HTML if ord(ch) < 32 and ch not in "\n\t"})
    assert stray == [], stray
    assert 'summary::before { content: "\\25B8 "; }' in INDEX_HTML
    assert 'summary::before { content: "\\25BE "; }' in INDEX_HTML


def test_a_conditional_attestation_line_can_actually_hide() -> None:
    # `.confirm-line { display: flex !important }` came after the page-wide
    # `[hidden]` rule at the same specificity, so `hidden` had no effect and
    # every conditional line showed on every author card.
    flex_rule = INDEX_HTML.index(".confirm-line { display: flex !important;")
    hidden_rule = INDEX_HTML.index(".confirm-line[hidden] { display: none !important; }")
    assert hidden_rule > flex_rule


def test_the_renderer_counts_the_study_list_it_built() -> None:
    body = _script()
    start = body.index("function renderPreparedState()")
    end = body.index("\n    function ", start + 10)
    renderer = body[start:end]
    assert "const allStudies = state.prepared.studies || [];" in renderer
    assert "${studies.length" not in renderer
    assert "`${allStudies.length} studies prepared`" in renderer


PRELUDE = r"""
const byId = {}, proxies = {};
function makeEl(id) {
  const store = {
    id, hidden: false, textContent: "", innerHTML: "", value: "", checked: false, disabled: false, open: false,
    dataset: {}, style: {}, children: [], childNodes: [], options: [], files: [],
    parentElement: null, parentNode: null, nextElementSibling: null, previousElementSibling: null,
    firstElementChild: null, lastElementChild: null, selectedIndex: 0, scrollTop: 0, scrollHeight: 0,
    offsetHeight: 0, offsetWidth: 0, clientHeight: 0, clientWidth: 0, tagName: "DIV", type: "", name: "",
    title: "", className: "", classList: { toggle() {}, add() {}, remove() {}, contains() { return false; } },
  };
  if (id) byId[id] = store;
  return new Proxy(function () {}, {
    get(_t, key) {
      if (key === Symbol.toPrimitive) return () => "";
      if (typeof key === "symbol" || key === "then") return undefined;
      if (Object.prototype.hasOwnProperty.call(store, key)) return store[key];
      if (key === "querySelectorAll" || key === "getElementsByClassName" || key === "getElementsByTagName") return () => [];
      if (key === "getAttribute") return () => null;
      if (key === "hasAttribute" || key === "contains" || key === "matches") return () => false;
      if (key === "getBoundingClientRect") return () => ({ top: 0, left: 0, right: 0, bottom: 0, width: 0, height: 0 });
      return makeEl("");
    },
    set(_t, key, value) { store[key] = value; return true; },
    has() { return true; },
    apply() { return makeEl(""); },
  });
}
const elementById = (id) => (proxies[id] = proxies[id] || makeEl(id));
const doc = makeEl("document");
doc.getElementById = elementById;
doc.body = makeEl("body");
doc.documentElement = makeEl("html");
doc.activeElement = makeEl("");
const define = (name, value) => Object.defineProperty(globalThis, name, { value, configurable: true, writable: true });
define("document", doc);
define("window", globalThis);
define("location", { href: "http://127.0.0.1/", origin: "http://127.0.0.1", pathname: "/", search: "", hash: "", protocol: "http:", host: "127.0.0.1" });
define("history", { replaceState() {}, pushState() {}, state: null });
define("navigator", { clipboard: { writeText: async () => {} }, language: "en-US", userAgent: "node-harness" });
const storage = () => ({ getItem: () => null, setItem() {}, removeItem() {}, clear() {} });
define("localStorage", storage());
define("sessionStorage", storage());
define("matchMedia", () => ({ matches: false, media: "", addEventListener() {}, removeEventListener() {}, addListener() {}, removeListener() {} }));
define("getComputedStyle", () => ({ getPropertyValue: () => "" }));
define("innerWidth", 1280);
define("innerHeight", 800);
for (const name of ["scrollTo", "scrollBy", "alert", "addEventListener", "removeEventListener", "cancelAnimationFrame", "clearTimeout", "clearInterval"]) define(name, () => {});
define("confirm", () => false);
define("prompt", () => null);
define("dispatchEvent", () => true);
define("requestAnimationFrame", () => 0);
define("setTimeout", () => 0);
define("setInterval", () => 0);
define("fetch", () => new Promise(() => {}));
"""

PREPARED_EPILOGUE = r"""
const results = {};
function scenario(name, prepared) {
  const state = activeDiscoveryState();
  state.prepared = JSON.parse(JSON.stringify(prepared));
  state.bundleId = prepared.bundle_id || "bundle";
  state.draft = {};
  state.run = null;
  state.preparing = false;
  let error = "";
  try { renderPreparedState(); } catch (e) { error = String((e && e.stack) || e); }
  results[name] = {
    error,
    status: byId.preparedStatus ? byId.preparedStatus.textContent : null,
    cardHidden: byId.preparedCard ? byId.preparedCard.hidden : null,
    html: byId.preparedCandidates ? byId.preparedCandidates.innerHTML : "",
  };
}
scenario("two_studies", FIXTURE);
const oneUsable = JSON.parse(JSON.stringify(FIXTURE));
oneUsable.studies[1].files[0].inspection.status = "not_deg_table";
scenario("one_usable", oneUsable);
// The completion card lives in the same renderer: a finished run must reach it.
{
  const state = activeDiscoveryState();
  state.prepared = JSON.parse(JSON.stringify(FIXTURE));
  state.run = { n_source_units: 2, top_genes: ["ISG15", "IFIT1", "MX1"], warnings: ["unit A: identifiers are Ensembl IDs", "unit A: identifiers are Ensembl IDs", "unit B: DEG-only table"], excel_workbook: "run.xlsx" };
  let error = "";
  try { renderPreparedState(); } catch (e) { error = String((e && e.stack) || e); }
  results.with_run = { error, cardHidden: byId.analysisCompleteCard ? byId.analysisCompleteCard.hidden : null,
    title: byId.analysisCompleteTitle ? byId.analysisCompleteTitle.textContent : null,
    text: byId.analysisCompleteText ? byId.analysisCompleteText.textContent : null,
    warningsHidden: byId.analysisCompleteWarnings ? byId.analysisCompleteWarnings.hidden : null,
    warningsHtml: byId.analysisCompleteWarnings ? byId.analysisCompleteWarnings.innerHTML : null,
    excelDisabled: byId.downloadAnalysisExcel ? byId.downloadAnalysisExcel.disabled : null };
}
console.log(JSON.stringify(results));
"""


FILTER_EPILOGUE = r"""
const state = activeDiscoveryState();
Object.assign(state, { query: "hypoxia", searchId: "abc123", verified: true, loading: false, error: "", studies: [],
  totalHits: 0, totalUnfiltered: 1000, textFilter: "GSE300988", providerStatus: "partial", providerErrors: ["europe_pmc (resolve)"] });
const run = () => { try { renderDiscoveryResults(); return ""; } catch (e) { return String((e && e.stack) || e); } };
const filtered = { error: run(), html: byId.discoveryResults ? byId.discoveryResults.innerHTML : "" };
Object.assign(state, { textFilter: "", totalUnfiltered: 0, totalHits: 0 });
const empty = { error: run(), html: byId.discoveryResults ? byId.discoveryResults.innerHTML : "" };
console.log(JSON.stringify({ filtered, empty }));
"""


COLUMNS_EPILOGUE = r"""
const author = FIXTURE.studies[1].files[0];
const detected = author.inspection.mapping;
const st = activeDiscoveryState();
st.prepared = JSON.parse(JSON.stringify(FIXTURE)); st.bundleId = "bundle"; st.run = null; st.preparing = false;
st.draft = { [author.candidate_id]: { enabled: true, sheetName: author.inspection.sheet_name || "", geneColumn: detected.gene_column,
  lfcColumn: detected.lfc_column, pColumn: detected.p_column, padjColumn: detected.padj_column, tableScope: "auto" } };
const isOpen = () => /class="candidate-advanced" open>/.test(byId.preparedCandidates.innerHTML);
renderPreparedState();
const unedited = isOpen();
st.draft[author.candidate_id].geneColumn = "SYMBOL";
renderPreparedState();
const edited = isOpen();
st.draft[author.candidate_id].geneColumn = detected.gene_column;
st.draft[author.candidate_id].tableScope = "deg_only";
renderPreparedState();
console.log(JSON.stringify({ unedited, edited, scoped: isOpen() }));
"""


def _render_in_node(tmp_path: Path, epilogue: str = PREPARED_EPILOGUE) -> dict:
    harness = PRELUDE + _script() + "\nconst FIXTURE = " + FIXTURE.read_text(encoding="utf-8") + ";\n" + epilogue
    script = tmp_path / "page_harness.js"
    script.write_text(harness, encoding="utf-8")
    result = subprocess.run(["node", str(script)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr[-2000:]
    return json.loads(result.stdout.strip().splitlines()[-1])


@pytest.mark.skipif(not HAS_NODE, reason="node not available")
def test_the_prepared_card_renders_a_real_bundle_without_throwing(tmp_path) -> None:
    out = _render_in_node(tmp_path)["two_studies"]
    assert out["error"] == "", out["error"]
    assert out["status"] == "2 studies prepared"
    assert out["cardHidden"] is False
    html = out["html"]
    # The least processed matrix in front, the other normalizations folded away.
    assert "Shown first: raw counts" in html
    assert html.index("GSE343715_SALMON_tx2gene_counts_matrix.txt.gz") < html.index('<details class="alternative-candidates">')
    assert "<summary>2 other files from this series" in html
    assert "Supplementary Table 7.xlsx" in html
    # ready_for_review with separate p and adjusted-p columns: only the
    # direction attestation is shown; the conditional ones start hidden.
    assert 'data-when="mapping" hidden' in html
    assert 'data-when="padj" hidden' in html
    assert 'data-when="lfc" hidden' in html
    assert "prepared-blocked" not in html
    assert "unanalyzable-group" not in html


@pytest.mark.skipif(not HAS_NODE, reason="node not available")
def test_a_bundle_with_one_usable_study_explains_the_block(tmp_path) -> None:
    out = _render_in_node(tmp_path)["one_usable"]
    assert out["error"] == "", out["error"]
    assert out["status"] == "2 prepared · 1 usable"
    html = out["html"]
    assert "1 of 2 prepared studies produced a usable candidate" in html
    assert "<summary>1 study with no usable table" in html
    # The study with nothing to activate comes after the one that has a candidate.
    assert html.index("GSE343715") < html.index('<details class="unanalyzable-group">')


@pytest.mark.skipif(not HAS_NODE, reason="node not available")
def test_a_filter_that_matches_nothing_keeps_its_box_and_blames_nobody(tmp_path) -> None:
    """Zero matches under a partial snapshot used to drop the filter input and
    show the "data sources did not answer" notice, so the text could not be
    cleared and the reader was told the providers were at fault."""
    out = _render_in_node(tmp_path, FILTER_EPILOGUE)
    filtered, empty = out["filtered"], out["empty"]
    assert filtered["error"] == "", filtered["error"]
    assert 'id="resultFilter"' in filtered["html"]
    assert 'value="GSE300988"' in filtered["html"]
    assert "No record matches \u201cGSE300988\u201d among the 1,000 assessed studies" in filtered["html"]
    assert "did not answer" not in filtered["html"]
    # The count chip beside the box read the honest zero as missing and showed
    # the snapshot total instead ("1,000 of 1,000 match" above "No record matches").
    assert "state.totalHits = Number(data.total ?? data.search?.total ?? 0);" in INDEX_HTML
    # With no filter and no records, the provider notice is still the right one.
    assert empty["error"] == "", empty["error"]
    assert "Some data sources did not answer" in empty["html"]
    assert 'id="resultFilter"' not in empty["html"]


def test_the_species_attestation_counts_records_in_the_plural() -> None:
    # "1 of these record was matched" - the noun after "of these" is plural
    # whatever the count; only the verb changes.
    assert '`${pending.length} of these records ${pending.length === 1 ? "was" : "were"} matched by the `' in INDEX_HTML
    assert "of these record${" not in INDEX_HTML


@pytest.mark.skipif(not HAS_NODE, reason="node not available")
def test_the_columns_panel_opens_only_for_an_edited_mapping(tmp_path) -> None:
    """The draft captures the prefilled detected mapping on every re-render;
    treating any captured value as "set" opened the panel as soon as the
    reader touched anything else on the card."""
    out = _render_in_node(tmp_path, COLUMNS_EPILOGUE)
    assert out == {"unedited": False, "edited": True, "scoped": True}


@pytest.mark.skipif(not HAS_NODE, reason="node not available")
def test_a_finished_run_reaches_the_completion_card(tmp_path) -> None:
    out = _render_in_node(tmp_path)["with_run"]
    assert out["error"] == "", out["error"]
    assert out["cardHidden"] is False
    assert out["title"] == "Human DEGORA analysis complete"
    assert out["text"].startswith("2 independent source units were analyzed separately")
    assert "ISG15, IFIT1, MX1" in out["text"]
    assert out["warningsHidden"] is False
    assert out["warningsHtml"].count("Ensembl IDs") == 1  # deduplicated
    assert "DEG-only table" in out["warningsHtml"]
    assert out["excelDisabled"] is False
