"""Component ablation and weight sensitivity for the primary quality-weighted rank.

``score_db.ScoreAblation`` can re-score a run with any subset of the five score
components, with source-quality weighting switched off, or with sample-size
weighting switched off, while holding the gene universe and the source-unit
set fixed. Until now it was reachable only from Python. This module turns one
finished run into the comparison a reader wants to see: how far each variant
moves the primary rank, how much of the top of the list it keeps, and - when a
GoldPanel or gene list is supplied - what happens to recall@k.

Everything here is a re-scoring of the harmonized table a run already wrote
(``slice_harmonized.csv``); no source table is read again and nothing about the
run's own outputs is changed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from .formula_safety import neutralize_formula_text, restore_formula_text_if_marked
from .harmonize import canonical_gene_symbol
from .score_db import PRIMARY_RANK_COLUMN, SCORE_WEIGHTS, ScoreAblation, degora_score_table

DEFAULT_TOP_K = (50, 100)
HARMONIZED_NAME = "slice_harmonized.csv"


def parse_weight_spec(spec: str) -> dict[str, float]:
    """Parse ``support_score=0.4,direction_score=0.3`` into a weight mapping.

    Component names may drop the ``_score``/``_component`` suffix (``support``,
    ``direction``, ``evidence``, ``rank``, ``effect``). Weights are validated by
    ``ScoreAblation``; this only turns text into numbers.
    """

    aliases = {
        "support": "support_score",
        "direction": "direction_score",
        "evidence": "evidence_score",
        "rank": "rank_score_component",
        "effect": "effect_score",
    }
    weights: dict[str, float] = {}
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"weight {part!r} must be written as component=value")
        name, value = (piece.strip() for piece in part.split("=", 1))
        name = aliases.get(name.lower(), name)
        if name not in SCORE_WEIGHTS:
            raise ValueError(f"unknown score component {name!r}; valid components are {sorted(SCORE_WEIGHTS)}")
        try:
            weights[name] = float(value)
        except ValueError as exc:
            raise ValueError(f"weight for {name} must be a number, got {value!r}") from exc
    if not weights:
        raise ValueError("a weight specification names at least one component")
    return weights


def default_ablations() -> list[ScoreAblation]:
    """The shipped ranking, five leave-one-component-out variants, and the two switches."""

    variants = [ScoreAblation(name="full")]
    for component in SCORE_WEIGHTS:
        variants.append(
            ScoreAblation(
                name=f"without_{component}",
                component_weights={name: weight for name, weight in SCORE_WEIGHTS.items() if name != component},
                notes=f"{component} removed; the remaining weights are renormalised",
            )
        )
    variants.append(
        ScoreAblation(name="without_source_quality_weighting", disable_source_quality_weighting=True)
    )
    variants.append(ScoreAblation(name="without_sample_size_weighting", disable_sample_size_weighting=True))
    return variants


def read_gene_list(path: str | Path) -> list[str]:
    """Read gold genes from a text/CSV list or from a DEGORA config's GoldPanel sheet."""

    source = Path(path)
    if source.suffix.lower() in {".xlsx", ".xls"}:
        from .slice_runner import _read_locked_gold_panel

        panel = _read_locked_gold_panel(source)
        if panel["status"] != "locked":
            raise ValueError(f"{source.name}: {panel['reason'] or 'no locked GoldPanel genes'}")
        return list(panel["genes"])
    text = source.read_text(encoding="utf-8-sig")
    genes: list[str] = []
    for line in text.splitlines():
        cell = line.split(",")[0].split("\t")[0].strip()
        if not cell or cell.startswith("#") or cell.lower() in {"gene", "gene_symbol", "symbol"}:
            continue
        symbol = canonical_gene_symbol(cell)
        if symbol:
            genes.append(symbol)
    return sorted(dict.fromkeys(genes))


def _rank_series(scores: pd.DataFrame) -> pd.Series:
    ranked = scores.set_index("gene_symbol")[PRIMARY_RANK_COLUMN]
    return pd.to_numeric(ranked, errors="coerce")


def _top(ranks: pd.Series, k: int) -> set[str]:
    return set(ranks[ranks.le(k)].index.astype(str))


def _recall(ranks: pd.Series, gold: Iterable[str], k: int) -> float | None:
    gold_set = {str(gene).strip().upper() for gene in gold if str(gene).strip()}
    if not gold_set:
        return None
    return len(gold_set & _top(ranks, k)) / len(gold_set)


def run_ablations(
    harmonized: pd.DataFrame,
    *,
    min_studies: int = 2,
    ablations: Iterable[ScoreAblation] | None = None,
    gold_genes: Iterable[str] | None = None,
    top_k: Iterable[int] = DEFAULT_TOP_K,
) -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    """Score every ablation and compare each primary rank with the full one.

    Returns a summary table (one row per variant) and the per-variant rank
    series so a caller can look at individual genes.
    """

    variants = list(ablations) if ablations is not None else default_ablations()
    if not variants or variants[0].name != "full" or not variants[0].is_default:
        variants = [ScoreAblation(name="full"), *[variant for variant in variants if variant.name != "full"]]
    gold = list(gold_genes or [])
    ks = [int(k) for k in top_k]
    ranks: dict[str, pd.Series] = {}
    rows: list[dict[str, Any]] = []
    full_ranks: pd.Series | None = None
    for variant in variants:
        scores, _evidence, _metadata = degora_score_table(
            harmonized, min_studies=min_studies, ablation=variant, include_loo_stability=False
        )
        current = _rank_series(scores)
        ranks[variant.name] = current
        if full_ranks is None:
            full_ranks = current
        aligned = pd.concat([full_ranks.rename("full"), current.rename("variant")], axis=1, join="inner")
        row: dict[str, Any] = {
            "ablation": variant.name,
            "component_weights": json.dumps(variant.weights, sort_keys=True),
            "source_quality_weighting": not variant.disable_source_quality_weighting,
            "sample_size_weighting": not variant.disable_sample_size_weighting,
            "n_genes": int(len(current)),
            "spearman_vs_full": float(aligned["full"].corr(aligned["variant"], method="spearman")) if len(aligned) > 1 else np.nan,
            "median_abs_rank_shift": float((aligned["full"] - aligned["variant"]).abs().median()) if len(aligned) else np.nan,
            "max_abs_rank_shift": float((aligned["full"] - aligned["variant"]).abs().max()) if len(aligned) else np.nan,
        }
        for k in ks:
            full_top = _top(full_ranks, k)
            variant_top = _top(current, k)
            # A corpus smaller than k cannot fill the top k; compare against what
            # exists rather than reporting 12 genes as a 24% overlap of 50.
            denominator = min(k, len(full_ranks))
            row[f"top{k}_overlap_with_full"] = (len(full_top & variant_top) / denominator) if denominator else np.nan
            if gold:
                row[f"recall_at_{k}"] = _recall(current, gold, k)
        row["notes"] = variant.notes
        rows.append(row)
    summary = pd.DataFrame(rows)
    return summary, ranks


def load_harmonized(path: str | Path) -> pd.DataFrame:
    """Read a run's harmonized table from a results folder or a CSV/parquet path."""

    source = Path(path)
    if source.is_dir():
        candidate = source / HARMONIZED_NAME
        if not candidate.exists():
            raise FileNotFoundError(
                f"{source} has no {HARMONIZED_NAME}; pass the results folder written by `degora run`, "
                "or the harmonized CSV itself"
            )
        source = candidate
    if not source.exists():
        raise FileNotFoundError(f"harmonized table does not exist: {source}")
    if source.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(source)
    frame = pd.read_csv(source, low_memory=False)
    return restore_formula_text_if_marked(frame, source)


def write_ablation_report(
    summary: pd.DataFrame,
    ranks: Mapping[str, pd.Series],
    output_dir: str | Path,
) -> dict[str, str]:
    """Write the summary table and the per-gene rank matrix beside each other."""

    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    summary_path = target / "degora_ablation_summary.csv"
    ranks_path = target / "degora_ablation_ranks.csv"
    neutralize_formula_text(summary).to_csv(summary_path, index=False)
    matrix = pd.DataFrame({name: series for name, series in ranks.items()})
    matrix.index.name = "gene_symbol"
    neutralize_formula_text(matrix.reset_index()).to_csv(ranks_path, index=False)
    return {"summary_csv": str(summary_path), "ranks_csv": str(ranks_path)}


def format_summary(summary: pd.DataFrame) -> str:
    """A terminal-friendly view of the summary table."""

    columns = ["ablation", "spearman_vs_full", "median_abs_rank_shift"] + [
        column for column in summary.columns if column.startswith("top") or column.startswith("recall_at_")
    ]
    view = summary[columns].copy()
    for column in columns[1:]:
        view[column] = pd.to_numeric(view[column], errors="coerce").map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    return view.to_string(index=False)
