#!/usr/bin/env python
"""Derive a gene-level DEG table from a normalized microarray expression matrix.

This is a fallback for already normalized, log-scale microarray matrices. A full
limma/eBayes table remains preferred and can be run on this same normalized
log-scale matrix (lmFit + eBayes need only a design matrix, not raw .CEL files);
the per-gene Welch test here is used for a transparent, dependency-free pipeline,
not because limma is infeasible on these inputs. At small per-group replicate
counts (n=2) the Welch-Satterthwaite p-value has ~1-2 degrees of freedom and is
effectively magnitude-driven -- the regime where limma's moderated variance helps
most -- so such contrasts are best treated as exploratory.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from degora.harmonize import _repair_excel_date_gene_symbol
from degora.provenance import write_source_sidecar


PROBE_COLLAPSE_RULES = ("min_pvalue_max_abs_lfc", "median_expression")


def _read_table(path: Path, *, sep: str | None = None, sheet_name: str | int | None = None) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, sheet_name=0 if sheet_name in (None, "") else sheet_name)
    if sep in (None, ""):
        sep = "\t" if suffixes.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return pd.read_csv(path, sep=sep)


def _parse_samples(value: str) -> list[str]:
    samples = [item.strip() for item in value.split(",") if item.strip()]
    if not samples:
        raise ValueError("sample list is empty")
    return samples


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).clip(0.0, 1.0)
    n = len(values)
    order = np.argsort(values.to_numpy(dtype=float))
    ranked = values.to_numpy(dtype=float)[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(out, index=pvalues.index)


def _clean_symbol(values: pd.Series) -> pd.Series:
    return (
        values.map(_repair_excel_date_gene_symbol)
        .astype("string")
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
        .str.upper()
        .replace({"": pd.NA, "NAN": pd.NA, "NONE": pd.NA})
    )


def _bh_frame(gene_results: pd.DataFrame, collapse_rule: str) -> pd.DataFrame:
    gene_results = gene_results.copy()
    gene_results["padj"] = _bh_adjust(gene_results["pvalue"])
    gene_results["assay_type"] = "microarray"
    gene_results["source_input_type"] = "normalized_expression_matrix"
    gene_results["probe_collapse"] = collapse_rule
    return gene_results


def _probe_level_welch(
    expression: pd.DataFrame,
    genes: pd.Series,
    probes: pd.Series,
    *,
    control_samples: list[str],
    treatment_samples: list[str],
) -> pd.DataFrame:
    control = expression[control_samples]
    treatment = expression[treatment_samples]
    logfc = treatment.mean(axis=1) - control.mean(axis=1)
    test = stats.ttest_ind(
        treatment.to_numpy(dtype=float),
        control.to_numpy(dtype=float),
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )
    pvalue = pd.Series(test.pvalue, index=expression.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    probe_results = pd.DataFrame(
        {
            "gene_symbol": genes.astype(str),
            "probe_id": probes.astype(str),
            "log2FoldChange": logfc,
            "pvalue": pvalue,
        }
    )
    probe_results["abs_lfc"] = probe_results["log2FoldChange"].abs()
    probe_results["n_probes_for_gene"] = probe_results.groupby("gene_symbol")["probe_id"].transform("count")
    return probe_results


def _collapse_best_probe(probe_results: pd.DataFrame) -> pd.DataFrame:
    return (
        probe_results.sort_values(["gene_symbol", "pvalue", "abs_lfc", "probe_id"], ascending=[True, True, False, True])
        .drop_duplicates("gene_symbol", keep="first")
        .drop(columns=["abs_lfc"])
        .sort_values(["pvalue", "gene_symbol"])
        .reset_index(drop=True)
    )


def _collapse_median_expression(
    expression: pd.DataFrame,
    genes: pd.Series,
    probes: pd.Series,
    *,
    control_samples: list[str],
    treatment_samples: list[str],
) -> pd.DataFrame:
    annotated = expression.copy()
    annotated.insert(0, "gene_symbol", genes.astype(str))
    annotated.insert(1, "probe_id", probes.astype(str))
    n_probes = annotated.groupby("gene_symbol")["probe_id"].size().rename("n_probes_for_gene")
    probe_ids = annotated.groupby("gene_symbol")["probe_id"].agg(lambda values: ";".join(sorted({str(value) for value in values if str(value)}))).rename("probe_id")
    gene_expression = annotated.groupby("gene_symbol")[control_samples + treatment_samples].median()
    control = gene_expression[control_samples]
    treatment = gene_expression[treatment_samples]
    logfc = treatment.mean(axis=1) - control.mean(axis=1)
    test = stats.ttest_ind(
        treatment.to_numpy(dtype=float),
        control.to_numpy(dtype=float),
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )
    pvalue = pd.Series(test.pvalue, index=gene_expression.index).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(0.0, 1.0)
    gene_results = pd.DataFrame(
        {
            "gene_symbol": gene_expression.index.astype(str),
            "probe_id": probe_ids.reindex(gene_expression.index).fillna("").to_numpy(),
            "log2FoldChange": logfc.to_numpy(dtype=float),
            "pvalue": pvalue.to_numpy(dtype=float),
            "n_probes_for_gene": n_probes.reindex(gene_expression.index).fillna(0).astype(int).to_numpy(),
        }
    )
    return gene_results.sort_values(["pvalue", "gene_symbol"]).reset_index(drop=True)


def derive_microarray_deg(
    matrix_path: Path,
    output_path: Path,
    *,
    gene_column: str,
    control_samples: list[str],
    treatment_samples: list[str],
    probe_column: str | None = None,
    sep: str | None = None,
    sheet_name: str | int | None = None,
    log2_transform: bool = False,
    collapse_rule: str = "min_pvalue_max_abs_lfc",
    command: str,
    summary_path: Path | None = None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if collapse_rule not in PROBE_COLLAPSE_RULES:
        raise ValueError(f"collapse_rule must be one of: {', '.join(PROBE_COLLAPSE_RULES)}")
    matrix = _read_table(matrix_path, sep=sep, sheet_name=sheet_name)
    if len(control_samples) < 2 or len(treatment_samples) < 2:
        raise ValueError(
            "microarray Welch fallback requires at least two control and two treatment sample columns; "
            "use a validated limma full table or mark the source as exploratory when replicate counts are lower."
        )
    required = [gene_column, *control_samples, *treatment_samples]
    if probe_column:
        required.append(probe_column)
    missing = sorted(set(required).difference(matrix.columns))
    if missing:
        raise ValueError(f"microarray matrix missing columns: {missing}")

    genes = _clean_symbol(matrix[gene_column])
    probes = matrix[probe_column].astype("string").fillna("") if probe_column else pd.Series(matrix.index.astype(str), index=matrix.index)
    expression = matrix[control_samples + treatment_samples].apply(pd.to_numeric, errors="coerce")
    if log2_transform:
        expression = np.log2(expression.clip(lower=0.0) + 1.0)

    # Require finite values: inf/-inf (e.g. from log of a sentinel or a bad cell) is not NaN,
    # so notna() alone would let it through and emit a +/-inf log2 fold-change.
    valid = genes.notna() & np.isfinite(expression).all(axis=1)
    expression = expression.loc[valid].copy()
    genes = genes.loc[valid]
    probes = probes.loc[valid]
    if expression.empty:
        raise ValueError("microarray matrix has no rows with complete gene and expression values after filtering")

    probe_results = _probe_level_welch(
        expression,
        genes,
        probes,
        control_samples=control_samples,
        treatment_samples=treatment_samples,
    )
    if collapse_rule == "min_pvalue_max_abs_lfc":
        gene_results = _collapse_best_probe(probe_results)
    else:
        gene_results = _collapse_median_expression(
            expression,
            genes,
            probes,
            control_samples=control_samples,
            treatment_samples=treatment_samples,
        )
    gene_results = _bh_frame(gene_results, collapse_rule)
    if metadata:
        for key in ["platform", "normalization"]:
            if key in metadata:
                gene_results[key] = metadata[key]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    gene_results.to_csv(output_path, index=False)
    write_source_sidecar(
        output_path,
        command,
        inputs=[matrix_path],
        metadata={"generator": "microarray-normalized-matrix-welch", "collapse_rule": collapse_rule, **(metadata or {})},
    )

    summary = {
        "matrix_path": str(matrix_path),
        "output_path": str(output_path),
        "assay_type": "microarray",
        "source_input_type": "normalized_expression_matrix",
        "pipeline": "welch_microarray_normalized_matrix",
        "preferred_pipeline": "limma/eBayes on this same normalized log matrix, or an author-provided validated full limma table",
        "collapse_rule": collapse_rule,
        "n_input_rows": int(len(matrix)),
        "n_valid_probe_rows": int(len(probe_results)),
        "n_gene_rows": int(len(gene_results)),
        "n_collapsed_probe_rows": int(len(probe_results) - len(gene_results)),
        "top_genes": gene_results.head(20)["gene_symbol"].tolist(),
        **(metadata or {}),
    }
    if summary_path:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        write_source_sidecar(
            summary_path,
            command,
            inputs=[matrix_path, output_path],
            metadata={"generator": "microarray-normalized-matrix-summary", **(metadata or {})},
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--gene-column", required=True)
    parser.add_argument("--probe-column")
    parser.add_argument("--control-samples", required=True, help="Comma-separated control sample columns.")
    parser.add_argument("--treatment-samples", required=True, help="Comma-separated treatment sample columns.")
    parser.add_argument("--sep")
    parser.add_argument("--sheet-name")
    parser.add_argument("--log2-transform", action="store_true")
    parser.add_argument("--collapse-rule", choices=PROBE_COLLAPSE_RULES, default="min_pvalue_max_abs_lfc")
    parser.add_argument("--study-id", default="")
    parser.add_argument("--paper-id", default="")
    parser.add_argument("--source-url", default="")
    parser.add_argument("--platform", default="")
    parser.add_argument("--normalization", default="")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    command_parts = [
        "PYTHONPATH=outputs/code python outputs/code/scripts/derive_microarray_deg.py"
        f" --matrix {args.matrix}",
        f" --output {args.output}",
        f" --gene-column {args.gene_column}",
        f" --control-samples {args.control_samples}",
        f" --treatment-samples {args.treatment_samples}",
        f" --collapse-rule {args.collapse_rule}",
    ]
    if args.summary:
        command_parts.append(f" --summary {args.summary}")
    if args.probe_column:
        command_parts.append(f" --probe-column {args.probe_column}")
    if args.sep:
        command_parts.append(f" --sep {args.sep}")
    if args.sheet_name:
        command_parts.append(f" --sheet-name {args.sheet_name}")
    if args.log2_transform:
        command_parts.append(" --log2-transform")
    command = "".join(command_parts)
    metadata = {
        "study_id": args.study_id,
        "paper_id": args.paper_id,
        "source_url": args.source_url,
        "platform": args.platform,
        "normalization": args.normalization,
    }
    metadata = {key: value for key, value in metadata.items() if value}
    summary = derive_microarray_deg(
        args.matrix,
        args.output,
        gene_column=args.gene_column,
        probe_column=args.probe_column,
        control_samples=_parse_samples(args.control_samples),
        treatment_samples=_parse_samples(args.treatment_samples),
        sep=args.sep,
        sheet_name=args.sheet_name,
        log2_transform=args.log2_transform,
        collapse_rule=args.collapse_rule,
        command=command,
        summary_path=args.summary,
        metadata=metadata,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
