"""Regression tests locking the behavior-changing fixes from the line-by-line audit.

These cover code paths the pre-existing suite did not exercise (exactly-zero log2FC,
rank-plane ties, missing-column lookup, network path redaction), so the fixes cannot
silently regress.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from degora.aggregate import stouffer_consensus
from degora.api import _looks_like_local_path, _redact_record_paths_for_network
from degora.excel_export import _curated_lookup
from degora.harmonize import TableMapping, harmonize_frame
from degora.rank_plane import rank_plane_gene_summary, rank_plane_points


def test_exactly_zero_lfc_yields_nan_signed_z() -> None:
    # lfc == 0 carries no direction; even with a tiny p-value its signed_z must be
    # NaN (not 0), so it is dropped downstream instead of diluting the Stouffer z.
    frame = pd.DataFrame({"gene": ["A", "B"], "lfc": [2.0, 0.0], "p": [1e-5, 1e-9]})
    out = harmonize_frame(frame, TableMapping("gene", "lfc", "p"), {"study_id": "S1"})
    z_a = out.loc[out["gene_symbol"] == "A", "signed_z"].iloc[0]
    z_b = out.loc[out["gene_symbol"] == "B", "signed_z"].iloc[0]
    assert np.isfinite(z_a)
    assert np.isnan(z_b)


def test_zero_lfc_unit_does_not_dilute_stouffer_or_inflate_n_studies() -> None:
    base = pd.DataFrame(
        {
            "study_id": ["S1", "S2"],
            "paper_id": ["P1", "P2"],
            "source_unit_id": ["P1", "P2"],
            "gene_symbol": ["G", "G"],
            "signed_z": [4.0, 4.0],
            "lfc": [2.0, 2.0],
            "normalized_rank": [0.1, 0.1],
            "n_ctrl": [3, 3],
            "n_treat": [3, 3],
        }
    )
    # harmonize now emits NaN signed_z for a zero-lfc row; such a unit must not change
    # the combined z nor count toward n_studies.
    with_zero = pd.concat(
        [
            base,
            pd.DataFrame(
                {
                    "study_id": ["S3"],
                    "paper_id": ["P3"],
                    "source_unit_id": ["P3"],
                    "gene_symbol": ["G"],
                    "signed_z": [np.nan],
                    "lfc": [0.0],
                    "normalized_rank": [0.9],
                    "n_ctrl": [3],
                    "n_treat": [3],
                }
            ),
        ],
        ignore_index=True,
    )
    z_base = stouffer_consensus(base, min_studies=2).iloc[0]
    z_with = stouffer_consensus(with_zero, min_studies=2).iloc[0]
    assert np.isclose(z_base["stouffer_z"], z_with["stouffer_z"])
    assert int(z_base["n_studies"]) == int(z_with["n_studies"]) == 2


def test_rank_plane_tie_row_is_neutral_not_discordant() -> None:
    # Gene with two up studies and one zero-lfc (tie) study that is the most
    # significant: concordance must be 1.0 (tie excluded), not 0.667 (tie counted
    # against the gene).
    harmonized = pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3"],
            "gene_symbol": ["G", "G", "G"],
            "pvalue": [1e-3, 1e-3, 1e-9],
            "lfc": [2.0, 2.0, 0.0],
        }
    )
    summary = rank_plane_gene_summary(rank_plane_points(harmonized))
    concordance = summary.loc[summary["gene_symbol"] == "G", "effect_sign_concordance"].iloc[0]
    assert np.isclose(concordance, 1.0)


def test_curated_lookup_survives_missing_rank_column() -> None:
    # A legacy/partial score frame without quality_weighted_degora_rank must not crash
    # the export (previously a scalar-vs-Series AttributeError).
    lookup = _curated_lookup(
        pd.DataFrame({"gene_symbol": ["A", "B"]}),
        pd.DataFrame({"gene_symbol": ["A", "C"]}),
    )
    assert "present_in_degora_output" in lookup.columns
    assert (~lookup["present_in_degora_output"]).all()


def test_network_redaction_is_value_based_and_catches_file_urls() -> None:
    record = _redact_record_paths_for_network(
        {
            "source_url": "/Users/x/f.csv",
            "contributing_source_urls": "file:///Users/y",
            "public": "https://example.test",
            "gene_symbol": "VEGFA",
        }
    )
    assert record["source_url"] == "[redacted: local path]"
    assert record["contributing_source_urls"] == "[redacted: local path]"
    assert record["public"] == "https://example.test"
    assert record["gene_symbol"] == "VEGFA"
    assert _looks_like_local_path("file:///Users/x")
    assert not _looks_like_local_path("https://example.test")
