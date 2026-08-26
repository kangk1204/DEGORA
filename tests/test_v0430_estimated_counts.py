"""0.4.30: estimated counts (Salmon, RSEM, kallisto) are named, offered and accepted.

The card put a Salmon `tx2gene_counts` matrix first as "raw counts - the least
processed matrix", and the run then refused it because 32% of its values were
whole numbers. Both were right about the file and wrong about each other.
"""

from __future__ import annotations

from degora.api import INDEX_HTML
from degora.discovery import _inspect_upstream_rows, candidate_preference


def _rows(values: list[list[float]]) -> list[list[str]]:
    header = ["gene", "ctrl_1", "ctrl_2", "drug_1", "drug_2"]
    return [header] + [[f"G{i}", *map(str, row)] for i, row in enumerate(values)]


def test_the_inspector_records_the_whole_number_share_of_sample_values() -> None:
    whole = _inspect_upstream_rows(_rows([[1, 2, 3, 4], [5, 6, 7, 8], [0, 0, 9, 10]]), declared_role="unknown_matrix")
    assert whole["status"] == "upstream_matrix_ready_for_contrast"
    assert whole["whole_number_share"] == 1.0
    fractional = _inspect_upstream_rows(_rows([[1.5, 2, 3.25, 4], [5.5, 6, 7.75, 8], [0.5, 0, 9.25, 10]]), declared_role="unknown_matrix")
    assert fractional["whole_number_share"] == 0.5


def test_a_fractional_count_file_is_called_estimated_counts_not_raw() -> None:
    base = {"name": "GSE1_salmon_gene_counts.txt.gz", "role": "unknown_matrix",
            "inspection": {"status": "upstream_matrix_ready_for_contrast"}}
    rank, reason = candidate_preference({**base, "inspection": {**base["inspection"], "whole_number_share": 0.32}})
    assert rank == 2
    assert reason.startswith("estimated counts (fractional, as Salmon, RSEM and kallisto write them)")
    rank, reason = candidate_preference({**base, "inspection": {**base["inspection"], "whole_number_share": 1.0}})
    assert (rank, reason) == (2, "raw counts - the least processed matrix; DEGORA normalises them itself")
    # A bundle prepared before the share was recorded keeps the old wording.
    assert candidate_preference(base)[1].startswith("raw counts")
    # A repository-declared count file is judged the same way.
    declared = {**base, "role": "count_matrix", "inspection": {**base["inspection"], "whole_number_share": 0.1}}
    assert candidate_preference(declared)[1].startswith("estimated counts")


def test_the_browser_offers_estimated_counts_and_gates_on_the_select_it_shows() -> None:
    assert 'matrixTypeOption(matrixType, "estimated_count_matrix", "Estimated counts (Salmon, RSEM, kallisto; fractional)")' in INDEX_HTML
    # A declared count file with fractional values gets the choice too.
    assert 'const fractionalCounts = role === "count_matrix" && wholeShare !== null && wholeShare < 0.95;' in INDEX_HTML
    assert "const resolvedRole = typeSelect ? typeSelect.value : row.dataset.role;" in INDEX_HTML
    assert "const matrixTypeValid = !matrixType || Boolean(matrixType.value);" in INDEX_HTML
    assert 'if (row.querySelector(".matrix-type")) common.matrix_type = row.querySelector(".matrix-type").value || "";' in INDEX_HTML
    assert 'row.dataset.role === "unknown_matrix" && !row.querySelector(".matrix-type")' not in INDEX_HTML


def test_an_emoji_query_is_refused_for_its_script_not_its_length() -> None:
    import pytest as _pytest

    from degora.api import DegoraRequestHandler

    with _pytest.raises(ValueError, match="no English term"):
        DegoraRequestHandler._discovery_create_publication_search(None, {"query": "\U0001F9EC", "species": "human"})
    with _pytest.raises(ValueError, match="at least 2 characters"):
        DegoraRequestHandler._discovery_create_publication_search(None, {"query": "a", "species": "human"})


def test_a_matrix_card_can_add_another_contrast_from_the_same_file() -> None:
    """A multi-arm design in one file could only activate one arm (D1)."""

    assert "function fallbackCandidateHtml(study, candidate, activationKey, clone)" in INDEX_HTML
    assert 'class="action-secondary clone-candidate" type="button" data-tip="A shared control against another treatment arm' in INDEX_HTML
    assert 'clone-author-candidate clone-candidate' in INDEX_HTML  # one handler for both card kinds
    assert 'event.target.closest(".clone-candidate")' in INDEX_HTML
    assert '+ extra.map((key) => fallbackCandidateHtml(study, candidate, key, true)).join("")' in INDEX_HTML
    assert 'role === "treatment" ? "" : role' in INDEX_HTML  # the clone keeps the control arm only
