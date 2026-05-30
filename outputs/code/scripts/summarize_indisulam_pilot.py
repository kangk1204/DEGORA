#!/usr/bin/env python
"""Summarize Indisulam DEGORA anchor recovery with strict claim guardrails."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from degora.provenance import write_source_sidecar


GUARDRAIL_FLAGS = {
    "derived_count_pilot": True,
    "as_published_author_deg_validation": False,
    "primary_min_studies": 2,
    "claim_scope": "exploratory_cross_source_pilot",
}


def _recall_at(score: pd.DataFrame, anchors: set[str], k: int) -> dict[str, Any]:
    observed = set(score.head(k)["gene_symbol"].astype("string").str.upper())
    recovered = sorted(observed.intersection(anchors))
    return {
        "k": k,
        "n_anchor_genes": len(anchors),
        "n_recovered": len(recovered),
        "anchor_recovery": round(len(recovered) / len(anchors), 6) if anchors else 0.0,
        "recovered": recovered,
        "missing_from_top_k": sorted(anchors.difference(observed)),
    }


def summarize(
    score_csv: Path,
    gold_csv: Path,
    output_json: Path,
    output_tsv: Path,
    *,
    command: str,
    corpus: str = "Indisulam derived-count pilot",
    extra_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    score = pd.read_csv(score_csv)
    gold = pd.read_csv(gold_csv)
    score["gene_symbol"] = score["gene_symbol"].astype("string").str.upper()
    gold["gene_symbol"] = gold["gene_symbol"].astype("string").str.upper()
    anchors = set(gold["gene_symbol"].dropna())

    keep = [
        "gene_symbol",
        "degora_rank",
        "rank_label",
        "evidence_tier",
        "degora_score",
        "top_percent_label",
        "support_label",
        "direction_label",
        "n_source_units",
        "n_contrasts_observed",
        "source_units",
        "weighted_lfc",
    ]
    anchor_ranks = gold.merge(score[[column for column in keep if column in score.columns]], on="gene_symbol", how="left")
    anchor_ranks = anchor_ranks.sort_values(["degora_rank", "gene_symbol"], na_position="last")
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    anchor_ranks.to_csv(output_tsv, sep="\t", index=False)

    cross_source = score.loc[pd.to_numeric(score["n_source_units"], errors="coerce").ge(2)].copy()
    guardrails = {**GUARDRAIL_FLAGS, **(extra_metadata or {})}
    summary = {
        "corpus": corpus,
        **guardrails,
        "interpretation_guardrail": (
            "Anchor genes are mechanism/pathway anchors, not guaranteed transcriptional responders, "
            "not validation positives, and not direction claims."
        ),
        "slice_metrics_guardrail": (
            "The generic slice_metrics.json file is produced by the hypoxia-oriented slice runner; "
            "use this Indisulam anchor summary for topic-specific anchor recovery."
        ),
        "n_scored_genes": int(len(score)),
        "n_anchor_genes": len(anchors),
        "n_anchor_genes_scored": int(anchor_ranks["degora_rank"].notna().sum()),
        "n_cross_source_scored_genes": int(len(cross_source)),
        "top_20_genes": score.head(20)["gene_symbol"].tolist(),
        "cross_source_top_20_genes": cross_source.head(20)["gene_symbol"].tolist(),
        "anchor_recovery": {str(k): _recall_at(score, anchors, k) for k in (10, 20, 50, 100)},
        "cross_source_anchor_recovery": {str(k): _recall_at(cross_source, anchors, k) for k in (10, 20, 50, 100)},
        "anchor_rank_table": str(output_tsv),
    }
    output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    for artifact in (output_json, output_tsv):
        write_source_sidecar(
            artifact,
            command,
            inputs=[score_csv, gold_csv],
            metadata={"generator": "indisulam-pilot-summary", **guardrails},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-csv", type=Path, required=True)
    parser.add_argument("--gold", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--corpus", default="Indisulam derived-count pilot")
    parser.add_argument("--extra-metadata", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/summarize_indisulam_pilot.py"
        f" --score-csv {args.score_csv}"
        f" --gold {args.gold}"
        f" --output-json {args.output_json}"
        f" --output-tsv {args.output_tsv}"
        f" --corpus {args.corpus}"
        + (f" --extra-metadata {args.extra_metadata}" if args.extra_metadata else "")
    )
    extra_metadata = json.loads(args.extra_metadata.read_text()) if args.extra_metadata else None
    summary = summarize(
        args.score_csv,
        args.gold,
        args.output_json,
        args.output_tsv,
        command=command,
        corpus=args.corpus,
        extra_metadata=extra_metadata,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
