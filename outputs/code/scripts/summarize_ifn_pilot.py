#!/usr/bin/env python
"""Summarize IFN pilot DEGORA recovery against the locked ISG panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.provenance import write_source_sidecar


def _gold_set(path: Path) -> set[str]:
    frame = pd.read_csv(path)
    return set(frame["gene_symbol"].astype("string").str.upper().dropna())


def _recall_at(score: pd.DataFrame, gold: set[str], k: int) -> dict[str, Any]:
    observed = set(score.head(k)["gene_symbol"].astype(str).str.upper())
    recovered = sorted(observed.intersection(gold))
    return {
        "k": k,
        "n_gold": len(gold),
        "n_recovered": len(recovered),
        "recall": round(len(recovered) / len(gold), 6) if gold else 0.0,
        "recovered": recovered,
        "missing": sorted(gold.difference(observed)),
    }


def summarize(score_csv: Path, gold_csv: Path, output_json: Path, output_tsv: Path, *, command: str) -> dict[str, Any]:
    score = pd.read_csv(score_csv)
    score["gene_symbol"] = score["gene_symbol"].astype("string").str.upper()
    gold = _gold_set(gold_csv)
    gold_ranks = score.loc[score["gene_symbol"].isin(gold)].copy()
    keep = [
        "degora_rank",
        "gene_symbol",
        "evidence_tier",
        "degora_score",
        "top_percent_label",
        "support_label",
        "direction_label",
        "n_source_units",
        "n_contrasts_observed",
        "weighted_lfc",
    ]
    gold_ranks = gold_ranks[[column for column in keep if column in gold_ranks.columns]].sort_values("degora_rank")
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    gold_ranks.to_csv(output_tsv, sep="\t", index=False)

    summary = {
        "corpus": "IFN response derived-count pilot",
        "claim_guardrail": "Derived-count pilot; not an as-published DEG-table validation corpus.",
        "n_scored_genes": int(len(score)),
        "n_gold_genes": len(gold),
        "top_20_genes": score.head(20)["gene_symbol"].tolist(),
        "recall": {str(k): _recall_at(score, gold, k) for k in (10, 20, 50, 100)},
        "gold_rank_table": str(output_tsv),
    }
    output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    for artifact in (output_json, output_tsv):
        write_source_sidecar(
            artifact,
            command,
            inputs=[score_csv, gold_csv],
            metadata={"generator": "ifn-pilot-summary"},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-csv", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/summarize_ifn_pilot.py"
        f" --score-csv {args.score_csv}"
        f" --gold {args.gold}"
        f" --output-json {args.output_json}"
        f" --output-tsv {args.output_tsv}"
    )
    summary = summarize(args.score_csv, args.gold, args.output_json, args.output_tsv, command=command)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
