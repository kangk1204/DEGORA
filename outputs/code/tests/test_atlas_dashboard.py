from __future__ import annotations

from pathlib import Path

from figures.make_degora_atlas_dashboard import _write_html


def test_atlas_static_html_escapes_script_breakout(tmp_path: Path) -> None:
    output = tmp_path / "atlas.html"

    _write_html(
        output,
        {
            "generated_at": "2026-05-29T00:00:00+00:00",
            "corpora": [{"corpus_id": "x</script><script>alert(1)</script>", "corpus_label": "X"}],
            "methods": [],
            "prior_art": [],
            "gene_columns": ["gene_symbol"],
            "genes": [["CMPK2</script>"]],
        },
    )

    html = output.read_text(encoding="utf-8")

    assert "x</script><script>" not in html
    assert "CMPK2</script>" not in html
    assert "\\u003c/script\\u003e" in html
    assert "const ATLAS =" in html


def test_atlas_static_html_treats_reliability_as_zero_to_100_score(tmp_path: Path) -> None:
    output = tmp_path / "atlas.html"

    _write_html(
        output,
        {
            "generated_at": "2026-05-29T00:00:00+00:00",
            "corpora": [{"corpus_id": "x", "corpus_label": "X"}],
            "methods": [],
            "prior_art": [],
            "gene_columns": [
                "corpus_id",
                "gene_symbol",
                "display_rank",
                "display_direction",
                "display_top_percent",
                "n_source_units",
                "evidence_reliability_score",
            ],
            "genes": [["x", "CMPK2", 1, "up", 0.1, 3, 98.7]],
        },
    )

    html = output.read_text(encoding="utf-8")

    assert 'pct100(geneVal(row, "evidence_reliability_score"), 1)' in html
    assert "fraction100(med(topRows.map(r => geneVal(r, \"evidence_reliability_score\")))" in html
