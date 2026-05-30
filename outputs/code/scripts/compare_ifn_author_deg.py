#!/usr/bin/env python
"""Compare IFN pilot DEGORA results with author-supplied DEG tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from degora.harmonize import TableMapping, assess_table_scope
from degora.provenance import write_source_sidecar


def _safe_neglog10(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").fillna(1.0).clip(lower=0.0, upper=1.0)
    positive = numeric[numeric > 0.0]
    floor = min(float(positive.min()) * 0.1, 1e-300) if not positive.empty else 1e-300
    return -np.log10(numeric.mask(numeric <= 0.0, floor))


def _read_author_ifnb(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name="DESeq2_NHBECells")
    required = {"GeneName", "IFNB_L2FC", "padj_IFNB"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing expected columns: {missing}")

    out = frame.loc[:, ["GeneName", "IFNB_L2FC", "padj_IFNB"]].rename(
        columns={
            "GeneName": "gene_symbol",
            "IFNB_L2FC": "author_ifnb_log2fc",
            "padj_IFNB": "author_ifnb_padj",
        }
    )
    out["gene_symbol"] = out["gene_symbol"].astype("string").str.upper()
    out["author_ifnb_log2fc"] = pd.to_numeric(out["author_ifnb_log2fc"], errors="coerce")
    out["author_ifnb_padj"] = pd.to_numeric(out["author_ifnb_padj"], errors="coerce")
    out = out.dropna(subset=["gene_symbol", "author_ifnb_log2fc", "author_ifnb_padj"])
    out["author_ifnb_significant_padj05"] = out["author_ifnb_padj"] < 0.05
    out["author_ifnb_neglog10_padj"] = _safe_neglog10(out["author_ifnb_padj"])
    out["author_ifnb_abs_log2fc"] = out["author_ifnb_log2fc"].abs()
    out["author_ifnb_direction"] = np.select(
        [out["author_ifnb_log2fc"] > 0, out["author_ifnb_log2fc"] < 0],
        ["up", "down"],
        default="flat",
    )
    out = out.sort_values(["author_ifnb_padj", "author_ifnb_abs_log2fc"], ascending=[True, False])
    out["author_ifnb_padj_rank"] = np.arange(1, len(out) + 1)
    return out


def _author_table_scope(path: Path) -> dict[str, Any]:
    frame = pd.read_excel(path, sheet_name="DESeq2_NHBECells")
    scope = assess_table_scope(
        frame,
        TableMapping(
            gene_column="GeneName",
            lfc_column="IFNB_L2FC",
            p_column="padj_IFNB",
            padj_column="padj_IFNB",
        ),
        declared_scope="auto",
    )
    return {
        "assessment": scope["assessment"],
        "effective_scope": scope["effective_scope"],
        "reason": scope["reason"],
        "value_column": scope["value_column"],
        "n_rows": scope["n_rows"],
        "n_le_0_05": scope["n_le_0_05"],
        "n_gt_0_05": scope["n_gt_0_05"],
        "fraction_le_0_05": scope["fraction_le_0_05"],
        "max_value": scope["max_value"],
    }


def _read_degora(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.upper()
    return frame


def _top_overlap(primary: pd.DataFrame, secondary: pd.DataFrame, k: int) -> dict[str, Any]:
    degora_top = set(primary.nsmallest(k, "degora_rank")["gene_symbol"])
    author_top = set(secondary.nsmallest(k, "author_ifnb_padj_rank")["gene_symbol"])
    overlap = sorted(degora_top.intersection(author_top))
    return {
        "k": k,
        "n_overlap": len(overlap),
        "overlap_fraction_of_k": round(len(overlap) / k, 6),
        "genes": overlap,
    }


def _degora_author_sig(primary: pd.DataFrame, merged: pd.DataFrame, k: int) -> dict[str, Any]:
    top = merged.nsmallest(k, "degora_rank")
    recovered = top.loc[top["author_ifnb_significant_padj05"], "gene_symbol"].tolist()
    total_author_sig = int(primary["author_ifnb_significant_padj05"].sum())
    return {
        "k": k,
        "n_degora_top_with_author_data": int(len(top)),
        "n_author_significant_in_degora_top": int(len(recovered)),
        "fraction_of_degora_top": round(len(recovered) / len(top), 6) if len(top) else 0.0,
        "fraction_of_all_author_significant": round(len(recovered) / total_author_sig, 6) if total_author_sig else 0.0,
        "genes": recovered,
    }


def _derived_comparison(author: pd.DataFrame, path: Path) -> dict[str, Any]:
    frame = pd.read_csv(path)
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.upper()
    frame["log2FoldChange"] = pd.to_numeric(frame["log2FoldChange"], errors="coerce")
    frame["padj"] = pd.to_numeric(frame["padj"], errors="coerce")
    merged = frame.merge(author, on="gene_symbol", how="inner")
    sig = frame.loc[frame["padj"] < 0.05, "gene_symbol"]
    author_sig = set(author.loc[author["author_ifnb_significant_padj05"], "gene_symbol"])
    top20 = set(frame.nsmallest(20, "padj")["gene_symbol"])
    return {
        "derived_deg": str(path),
        "n_derived_genes": int(len(frame)),
        "n_derived_padj05": int((frame["padj"] < 0.05).sum()),
        "n_common_with_author": int(len(merged)),
        "derived_sig_overlap_author_sig": int(sig.isin(author_sig).sum()),
        "derived_top20_overlap_author_sig": sorted(top20.intersection(author_sig)),
        "spearman_log2fc_vs_author_ifnb_log2fc": round(
            float(merged[["log2FoldChange", "author_ifnb_log2fc"]].corr(method="spearman").iloc[0, 1]), 6
        )
        if len(merged) > 1
        else None,
    }


def compare(
    *,
    author_xlsx: Path,
    degora_score_csv: Path,
    derived_degs: list[Path],
    author_output_csv: Path,
    degora_top_output_tsv: Path,
    author_top_output_tsv: Path,
    availability_output_tsv: Path,
    summary_output_json: Path,
    command: str,
) -> dict[str, Any]:
    author = _read_author_ifnb(author_xlsx)
    author_scope = _author_table_scope(author_xlsx)
    degora = _read_degora(degora_score_csv)
    merged = degora.merge(author, on="gene_symbol", how="inner")

    author_output_csv.parent.mkdir(parents=True, exist_ok=True)
    degora_top_output_tsv.parent.mkdir(parents=True, exist_ok=True)
    availability_output_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_output_json.parent.mkdir(parents=True, exist_ok=True)

    author.to_csv(author_output_csv, index=False)

    degora_cols = [
        "degora_rank",
        "rank_label",
        "gene_symbol",
        "degora_score",
        "top_percent_label",
        "support_label",
        "direction_label",
        "weighted_lfc",
        "author_ifnb_log2fc",
        "author_ifnb_padj",
        "author_ifnb_padj_rank",
        "author_ifnb_significant_padj05",
        "author_ifnb_direction",
    ]
    merged.nsmallest(50, "degora_rank").loc[:, degora_cols].to_csv(degora_top_output_tsv, sep="\t", index=False)

    author_top_cols = [
        "author_ifnb_padj_rank",
        "gene_symbol",
        "author_ifnb_log2fc",
        "author_ifnb_padj",
        "author_ifnb_significant_padj05",
        "degora_rank",
        "degora_score",
        "top_percent_label",
        "support_label",
        "direction_label",
    ]
    merged.nsmallest(50, "author_ifnb_padj_rank").loc[:, author_top_cols].to_csv(
        author_top_output_tsv, sep="\t", index=False
    )

    availability = pd.DataFrame(
        [
            {
                "source_unit_id": "GSE147507_NHBE_IFNB",
                "geo_accession": "GSE147507",
                "geo_full_deg_table": "no",
                "paper_full_deg_table": "yes",
                "paper_table_scope": author_scope["effective_scope"],
                "paper_table_scope_assessment": author_scope["assessment"],
                "comparison_status": "compared",
                "notes": (
                    "Cell 2020 supplementary workbook mmc2 contains DESeq2_NHBECells "
                    "IFNB_L2FC/padj_IFNB for hIFNB treatment 4-12hrs; it is not separated "
                    "into the 4h/6h/12h contrasts used in the derived-count pilot."
                ),
            },
            {
                "source_unit_id": "GSE221804_HUHWT_IFNA",
                "geo_accession": "GSE221804",
                "geo_full_deg_table": "no",
                "paper_full_deg_table": "no_full_ranked_table_found",
                "paper_table_scope": "not_available",
                "paper_table_scope_assessment": "not_available",
                "comparison_status": "not_compared",
                "notes": (
                    "GEO provides raw count CSVs for HuhWT IFNa. The linked CMLS paper documents "
                    "DESeq2 LRT and provides integrated gene-list supplements, but no full per-timepoint "
                    "ranked IFNa DEG table matching this pilot."
                ),
            },
        ]
    )
    availability.to_csv(availability_output_tsv, sep="\t", index=False)

    common = merged.copy()
    direction_known = common.loc[common["author_ifnb_direction"].isin(["up", "down"])].copy()
    direction_known["degora_direction"] = np.select(
        [direction_known["weighted_lfc"] > 0, direction_known["weighted_lfc"] < 0],
        ["up", "down"],
        default="flat",
    )
    direction_agree = direction_known["degora_direction"].eq(direction_known["author_ifnb_direction"])

    summary: dict[str, Any] = {
        "claim_guardrail": (
            "GSE147507 comparison uses the authors' aggregate NHBE IFNB 4-12hr DESeq2 table; "
            "it is not a per-timepoint author table."
        ),
        "author_table": str(author_xlsx),
        "author_table_scope": author_scope,
        "normalized_author_ifnb_csv": str(author_output_csv),
        "availability_table": str(availability_output_tsv),
        "n_author_genes": int(len(author)),
        "n_author_ifnb_padj05": int(author["author_ifnb_significant_padj05"].sum()),
        "n_degora_scored_genes": int(len(degora)),
        "n_common_author_degora": int(len(common)),
        "degora_top_author_significant_overlap": {
            str(k): _degora_author_sig(author, common, k) for k in (10, 20, 50, 100)
        },
        "degora_top_vs_author_top_by_padj_overlap": {
            str(k): _top_overlap(common, author, k) for k in (10, 20, 50, 100)
        },
        "spearman_degora_score_vs_author_neglog10_padj": round(
            float(common[["degora_score", "author_ifnb_neglog10_padj"]].corr(method="spearman").iloc[0, 1]), 6
        ),
        "spearman_degora_weighted_lfc_vs_author_log2fc": round(
            float(common[["weighted_lfc", "author_ifnb_log2fc"]].corr(method="spearman").iloc[0, 1]), 6
        ),
        "direction_agreement_common_author_degora": {
            "n_tested": int(len(direction_known)),
            "n_agree": int(direction_agree.sum()),
            "fraction_agree": round(float(direction_agree.mean()), 6) if len(direction_known) else 0.0,
        },
        "derived_timepoint_comparison": [_derived_comparison(author, path) for path in derived_degs],
        "degora_top50_author_annotation": str(degora_top_output_tsv),
        "author_top50_degora_annotation": str(author_top_output_tsv),
    }
    summary_output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    inputs = [author_xlsx, degora_score_csv, *derived_degs]
    metadata = {"generator": "ifn-author-deg-comparison", "author_source": "Blanco-Melo Cell 2020 mmc2"}
    for artifact in (
        author_output_csv,
        degora_top_output_tsv,
        author_top_output_tsv,
        availability_output_tsv,
        summary_output_json,
    ):
        write_source_sidecar(artifact, command, inputs=inputs, metadata=metadata)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--author-xlsx", type=Path, required=True)
    parser.add_argument("--degora-score-csv", type=Path, required=True)
    parser.add_argument("--derived-deg", type=Path, nargs="+", action="append", default=[])
    parser.add_argument("--author-output-csv", type=Path, required=True)
    parser.add_argument("--degora-top-output-tsv", type=Path, required=True)
    parser.add_argument("--author-top-output-tsv", type=Path, required=True)
    parser.add_argument("--availability-output-tsv", type=Path, required=True)
    parser.add_argument("--summary-output-json", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    derived_degs = [path for group in args.derived_deg for path in group]
    command = (
        "PYTHONPATH=outputs/code python outputs/code/scripts/compare_ifn_author_deg.py"
        f" --author-xlsx {args.author_xlsx}"
        f" --degora-score-csv {args.degora_score_csv}"
        + "".join(f" --derived-deg {path}" for path in derived_degs)
        + f" --author-output-csv {args.author_output_csv}"
        + f" --degora-top-output-tsv {args.degora_top_output_tsv}"
        + f" --author-top-output-tsv {args.author_top_output_tsv}"
        + f" --availability-output-tsv {args.availability_output_tsv}"
        + f" --summary-output-json {args.summary_output_json}"
    )
    summary = compare(
        author_xlsx=args.author_xlsx,
        degora_score_csv=args.degora_score_csv,
        derived_degs=derived_degs,
        author_output_csv=args.author_output_csv,
        degora_top_output_tsv=args.degora_top_output_tsv,
        author_top_output_tsv=args.author_top_output_tsv,
        availability_output_tsv=args.availability_output_tsv,
        summary_output_json=args.summary_output_json,
        command=command,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
