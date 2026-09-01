"""Species-boundary regressions for HGNC retired-symbol resolution."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from degora.ablation import read_gene_list
from degora.harmonize import TableMapping, canonical_gene_symbol, harmonize_frame
from degora.slice_runner import _read_locked_gold_panel, infer_single_species

MOUSE_FALSE_COLLISIONS = [
    ("Arf2", "Arf4"),
    ("Defa2", "Defa1"),
    ("Pabpc2", "Pabpc1"),
    ("Sult1c1", "Sult1c2"),
    ("Tac2", "Tac1"),
]


def _harmonize(genes: list[str], *, species: object = None) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "gene": genes,
            "lfc": [float(index + 1) for index in range(len(genes))],
            "p": [0.01 * (index + 1) for index in range(len(genes))],
        }
    )
    return harmonize_frame(
        frame,
        TableMapping("gene", "lfc", "p"),
        {"study_id": "scope", "species": species},
    )


def _write_gold_panel(path: Path, genes: list[str]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame({"gene_symbol": genes, "locked": ["yes"] * len(genes)}).to_excel(
            writer, sheet_name="GoldPanel", index=False
        )


def test_generic_symbol_resolution_keeps_pre_v0439_behavior() -> None:
    assert canonical_gene_symbol("CTGF") == "CTGF"
    assert canonical_gene_symbol("9-Sep") == "SEPTIN9"
    assert canonical_gene_symbol("ENSG00000141510.16") == "ENSG00000141510"


def test_hgnc_retirements_require_explicit_human_species() -> None:
    assert canonical_gene_symbol("CTGF", species="Homo sapiens") == "CCN2"
    assert canonical_gene_symbol("ctgf", species="human") == "CCN2"
    for species in (None, "", pd.NA, "Mus musculus", "mouse", "unknown"):
        assert canonical_gene_symbol("CTGF", species=species) == "CTGF"


@pytest.mark.parametrize(("mouse_old", "mouse_current"), MOUSE_FALSE_COLLISIONS)
def test_mouse_symbols_that_collide_in_hgnc_remain_distinct(
    mouse_old: str, mouse_current: str
) -> None:
    out = _harmonize([mouse_old, mouse_current], species="Mus musculus")
    assert set(out["gene_symbol"]) == {mouse_old.upper(), mouse_current.upper()}
    assert len(out) == 2


def test_human_harmonization_merges_retired_and_current_symbols() -> None:
    out = _harmonize(["CTGF", "CCN2"], species="Homo sapiens")
    assert out["gene_symbol"].tolist() == ["CCN2"]
    assert set(out["input_gene_label"].iloc[0].split(";")) == {"CTGF", "CCN2"}


def test_active_catalog_species_inference_ignores_excluded_rows() -> None:
    catalog = pd.DataFrame(
        {
            "species": ["human", "Homo sapiens", "Mus musculus"],
            "include_in_analysis": ["yes", "yes", "no"],
        }
    )
    assert infer_single_species(catalog) == "Homo sapiens"


@pytest.mark.parametrize(
    "species",
    [
        ["Homo sapiens", ""],
        ["Homo sapiens", "Mus musculus"],
        ["", ""],
    ],
)
def test_blank_or_mixed_active_catalog_fails_closed(species: list[str]) -> None:
    catalog = pd.DataFrame(
        {"species": species, "include_in_analysis": ["yes"] * len(species)}
    )
    assert infer_single_species(catalog) is None


def test_gold_panel_uses_only_safely_inferred_species(tmp_path: Path) -> None:
    config = tmp_path / "config.xlsx"
    _write_gold_panel(config, ["CTGF"])

    human = pd.DataFrame(
        {"species": ["Homo sapiens"], "include_in_analysis": ["yes"]}
    )
    mixed = pd.DataFrame(
        {
            "species": ["Homo sapiens", "Mus musculus"],
            "include_in_analysis": ["yes", "yes"],
        }
    )
    assert _read_locked_gold_panel(
        config, species=infer_single_species(human)
    )["genes"] == ["CCN2"]
    assert _read_locked_gold_panel(
        config, species=infer_single_species(mixed)
    )["genes"] == ["CTGF"]


@pytest.mark.parametrize(("mouse_old", "mouse_current"), MOUSE_FALSE_COLLISIONS)
def test_mouse_gold_panel_does_not_collapse_distinct_symbols(
    tmp_path: Path, mouse_old: str, mouse_current: str
) -> None:
    config = tmp_path / "config.xlsx"
    _write_gold_panel(config, [mouse_old, mouse_current])
    panel = _read_locked_gold_panel(config, species="Mus musculus")
    assert panel["genes"] == sorted([mouse_old.upper(), mouse_current.upper()])


def test_ablation_gene_list_accepts_explicit_species(tmp_path: Path) -> None:
    genes = tmp_path / "gold.txt"
    genes.write_text("CTGF\nCCN2\n", encoding="utf-8")
    assert read_gene_list(genes) == ["CCN2", "CTGF"]
    assert read_gene_list(genes, species="Homo sapiens") == ["CCN2"]
