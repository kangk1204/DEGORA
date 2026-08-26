"""The prepared-evidence card, made concise for a first-time reader (v0.4.28)."""

from __future__ import annotations

import re
import subprocess

import pytest

from degora.api import INDEX_HTML


def _script() -> str:
    return re.search(r"<script>(.*?)</script>", INDEX_HTML, re.S).group(1)


def test_author_attestations_beyond_direction_appear_only_when_they_apply() -> None:
    """Six attestation lines stood at equal weight under every table."""

    assert INDEX_HTML.count('class="confirm-line confirm-primary"') == 2  # one per card kind
    for when in ("mapping", "padj", "lfc", "filter", "duplicates"):
        assert f'class="confirm-line confirm-when" data-when="{when}"' in INDEX_HTML
    assert "function refreshConditionalConfirms(row)" in INDEX_HTML
    assert 'const needsMappingConfirm = status !== "ready_for_review";' in INDEX_HTML


def test_a_matrix_has_one_attestation_that_stands_for_both_facts() -> None:
    assert "The Control and Treatment groups assigned in the sample list are right, the comparison is treated minus control, and each column is a separate biological sample, not a repeat measurement of the same one" in INDEX_HTML
    assert 'querySelectorAll(".candidate-row").forEach(refreshConditionalConfirms);' in INDEX_HTML  # re-render keeps the required lines visible
    assert '<input class="biological-replicates-confirmed" type="checkbox" hidden aria-hidden="true"' in INDEX_HTML
    assert 'event.target.classList.contains("direction-confirmed")' in INDEX_HTML


def test_table_scope_lives_in_the_collapsed_columns_panel() -> None:
    """Auto is right for nearly every table; the choice belongs behind the fold."""

    assert "Columns and table scope DEGORA read from the file" in INDEX_HTML
    start = INDEX_HTML.index('<details class="candidate-advanced" ${columnsOpen')
    end = INDEX_HTML.index("</details>", start)
    panel = INDEX_HTML[start:end]
    assert '<label>Table scope<select class="table-scope">' in panel
    assert '<label>Sheet name<input class="sheet-name"' in panel


def test_the_card_offers_a_group_suggestion_from_sample_labels() -> None:
    assert 'class="action-secondary sample-suggest" type="button">Suggest groups from sample labels</button>' in INDEX_HTML
    assert "function sampleSuggestHtml()" in INDEX_HTML
    assert '${sampleSuggestHtml()}${columns.length > 4 ? sampleBulkHtml() : ""}' in INDEX_HTML  # offered on every matrix
    assert "function suggestSampleGroups(row)" in INDEX_HTML
    assert "function suggestedContrastLabel(name)" in INDEX_HTML
    # A Python string carrying a JS regex: the escapes must survive into the page.
    assert "const CONTROL_LABEL_RE = /(?<!anti[ -]?)(?<!non[ -]?)(?:^|[^a-z]|(?<=s[ih]))(control|ctrl|" in INDEX_HTML


@pytest.mark.skipif(subprocess.run(["which", "node"], capture_output=True).returncode != 0, reason="node not available")
def test_group_suggestion_and_label_prefill_behave(tmp_path) -> None:
    """Run the two pure-ish functions in node against a minimal fake row."""

    js = _script()

    def fn(name: str) -> str:
        i = js.index(f"    function {name}(")
        j = js.index("\n    }\n", i) + len("\n    }\n")
        return js[i:j]

    const = js[js.index("    const CONTROL_LABEL_RE") : js.index("\n", js.index("    const CONTROL_LABEL_RE")) + 1]
    harness = const + fn("suggestedContrastLabel") + fn("suggestSampleGroups") + r'''
function updateSampleCounts(){} function capturePreparedDraft(){} function updateAnalysisEligibility(){}
function mkItem(id, title, traits, missing){ const sel={value:""}; return {
  getAttribute:(k)=> k==="title" ? [id, title, traits].filter(Boolean).join(" — ") : null,
  querySelector:(q)=> q===".sample-traits"?(traits?{textContent:traits}:null): q===".sample-label-missing"?(missing?{}:null): q===".sample-label"?{textContent:title}: q===".sample-id"?{textContent:id}: q==="[data-sample]"?sel: null, sel }; }
function mkRow(items){ const note={textContent:""}, label={value:""}; return { querySelectorAll:(q)=> q===".sample-item"?items:[], querySelector:(q)=> q===".sample-suggest-note"?note: q===".contrast-label"?label: null, note, label }; }
const out = [];
let items=[ mkItem("s1","Vehicle_1","treatment: vehicle"), mkItem("s2","Vehicle_2","treatment: vehicle"), mkItem("s3","Drug_1","treatment: rapamycin"), mkItem("s4","Drug_2","treatment: rapamycin") ];
let row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value).join(",") + "|" + row.label.value);
items=[ mkItem("s1","siCtrl_1","transfection: siControl"), mkItem("s2","siCtrl_2","transfection: siControl"), mkItem("s3","siHIF_1","transfection: siHIF1A"), mkItem("s4","siHIF_2","transfection: siHIF1A") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value).join(",") + "|" + row.label.value);
items=[ mkItem("A485_rep1","A485_rep1","treatment: A-485"), mkItem("CPI_rep1","CPI_rep1","treatment: CPI-637"), mkItem("A485_rep2","A485_rep2","treatment: A-485"), mkItem("CPI_rep2","CPI_rep2","treatment: CPI-637") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value).join(",") + "|" + row.note.textContent);
out.push(suggestedContrastLabel("GSE300988_RNAseq_ATRA_vs_cntrl_SKNO1_gene_deseq2_out.txt.gz") + "|" + suggestedContrastLabel("Supplementary Table 7.xlsx"));
// The adversarial cases an independent review raised.
items=[ mkItem("s1","V1","treatment: vehicle"), mkItem("s2","V2","treatment: vehicle"), mkItem("s3","A1","treatment: drugA"), mkItem("s4","A2","treatment: drugA"), mkItem("s5","B1","treatment: drugB"), mkItem("s6","B2","treatment: drugB") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value||"-").join(",") + "|" + row.note.textContent.slice(0, 40));
items=[ mkItem("s1","C1","treatment: control"), mkItem("s2","C2","treatment: control"), mkItem("s3","M1","treatment: mock"), mkItem("s4","M2","treatment: mock") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value||"-").join(",") + "|" + (row.note.textContent.includes("both read as a control") ? "both" : "?"));
items=[ mkItem("s1","I1","antibody: IgG"), mkItem("s2","I2","antibody: IgG"), mkItem("s3","a1","antibody: anti-control antibody"), mkItem("s4","a2","antibody: anti-control antibody") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value||"-").join(",") + "|" + (row.note.textContent.includes("neither reads") ? "neither" : "?"));
items=[ mkItem("s1","Control_1","treatment: vehicle"), mkItem("s2","Control_2","treatment: vehicle"), mkItem("s3","Drug_1","treatment: drug"), mkItem("s4","","",true), mkItem("s5","Drug_2","treatment: drug") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value||"-").join(","));
items=[ mkItem("s1","V1","treatment: vehicle · dose: 0"), mkItem("s2","V2","treatment: vehicle · dose: 0"), mkItem("s3","L1","treatment: drug · dose: 1uM"), mkItem("s4","L2","treatment: drug · dose: 1uM"), mkItem("s5","H1","treatment: drug · dose: 10uM"), mkItem("s6","H2","treatment: drug · dose: 10uM") ];
row=mkRow(items); suggestSampleGroups(row); out.push(row.note.textContent.includes("dose still varies inside a group") ? "dose-noted" : "?");
items=[ mkItem("s1","a","treatment: DMSO"), mkItem("s2","b","treatment: dmso"), mkItem("s3","c","treatment: drug"), mkItem("s4","d","treatment: drug") ];
row=mkRow(items); suggestSampleGroups(row); out.push(items.map(i=>i.sel.value||"-").join(","));
console.log(JSON.stringify(out));
'''
    script = tmp_path / "harness.js"
    script.write_text(harness, encoding="utf-8")
    result = subprocess.run(["node", str(script)], capture_output=True, text=True, check=True)
    lines = __import__("json").loads(result.stdout.strip())
    assert lines[0] == "control,control,treatment,treatment|rapamycin vs vehicle"
    assert lines[1] == "control,control,treatment,treatment|siHIF1A vs siControl"
    assert lines[2].startswith(",,,|Two groups by treatment (A-485 / CPI-637), but neither reads as a control")
    assert lines[3] == "ATRA vs cntrl|"
    assert lines[4].startswith("-,-,-,-,-,-|treatment has 3 values (vehicle / drugA")  # a three-arm design is not split
    assert lines[5] == "-,-,-,-|both"
    assert lines[6] == "-,-,-,-|neither"  # "anti-control" is not a control
    assert lines[7] == "control,control,treatment,-,treatment"  # an unmatched column stays at Ignore
    assert lines[8] == "dose-noted"
    assert lines[9] == "control,control,treatment,treatment"  # DMSO and dmso are one value
