"""Diagnostic tables for post-slice curation and heterogeneity audits."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .aggregate import slice_consensus
from .gold import HIF1A_UP_TARGETS
from .metrics import recall_at_k
from .provenance import write_source_sidecar
from .slice_runner import catalog_include_mask, read_catalog


TARGETS = sorted(HIF1A_UP_TARGETS)


def active_catalog(catalog_path: str | Path) -> pd.DataFrame:
    catalog = read_catalog(catalog_path)
    return catalog.loc[catalog_include_mask(catalog)].copy()


def _uses_padj_surrogate(row: pd.Series) -> bool:
    p_column = str(row.get("p_column", "")).strip().lower()
    padj_column = str(row.get("padj_column", "")).strip().lower()
    notes = str(row.get("notes", "")).lower()
    return (p_column and p_column == padj_column) or "padj" in notes and "surrogate" in notes


def target_support_by_study(harmonized: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    study_targets = pd.MultiIndex.from_product(
        [catalog["study_id"].tolist(), TARGETS],
        names=["study_id", "gene_symbol"],
    ).to_frame(index=False)

    support = (
        harmonized[harmonized["gene_symbol"].isin(TARGETS)]
        .groupby(["study_id", "gene_symbol"], as_index=False)
        .agg(
            lfc=("lfc", "mean"),
            pvalue=("pvalue", "min"),
            padj=("padj", "min"),
            signed_z=("signed_z", "mean"),
            within_study_rank=("within_study_rank", "min"),
            normalized_rank=("normalized_rank", "min"),
            n_genes_in_study=("n_genes_in_study", "max"),
        )
    )
    meta = catalog[
        [
            "study_id",
            "paper_id",
            "pipeline",
            "cell_system",
            "hypoxia_modality",
            "duration_h",
            "p_column",
            "padj_column",
            "notes",
        ]
    ].copy()
    meta["uses_padj_surrogate"] = meta.apply(_uses_padj_surrogate, axis=1)

    out = study_targets.merge(meta, on="study_id", how="left").merge(support, on=["study_id", "gene_symbol"], how="left")
    out["target_present"] = out["lfc"].notna()
    out["observed_direction"] = np.select(
        [out["lfc"].gt(0), out["lfc"].lt(0)],
        ["up_in_hypoxia", "down_in_hypoxia"],
        default="missing_or_zero",
    )
    return out.sort_values(["study_id", "gene_symbol"]).reset_index(drop=True)


def orientation_audit(support: pd.DataFrame) -> pd.DataFrame:
    out = support.copy()
    out["expected_direction"] = "up_in_hypoxia"
    out["orientation_ok"] = np.where(out["target_present"], out["lfc"].gt(0), pd.NA)
    positive_fraction = (
        out.loc[out["target_present"]]
        .groupby("study_id")["orientation_ok"]
        .mean()
        .rename("study_positive_fraction")
        .reset_index()
    )
    out = out.merge(positive_fraction, on="study_id", how="left")
    return out[
        [
            "study_id",
            "paper_id",
            "gene_symbol",
            "expected_direction",
            "observed_direction",
            "lfc",
            "pvalue",
            "padj",
            "orientation_ok",
            "study_positive_fraction",
            "uses_padj_surrogate",
        ]
    ].sort_values(["study_id", "gene_symbol"]).reset_index(drop=True)


def _rerank(consensus: pd.DataFrame) -> pd.DataFrame:
    out = consensus.copy().reset_index(drop=True)
    if not out.empty:
        out["slice_rank"] = np.arange(1, len(out) + 1)
    return out


def _recall(consensus: pd.DataFrame, k: int) -> dict[str, Any]:
    if "gene_symbol" not in consensus.columns:
        return {
            "k": k,
            "n_positives": len(HIF1A_UP_TARGETS),
            "n_recovered": 0,
            "recall": 0.0,
            "recovered": [],
            "missing": sorted(HIF1A_UP_TARGETS),
        }
    return recall_at_k(consensus, HIF1A_UP_TARGETS, k)


def _variant_frame(harmonized: pd.DataFrame, variant: str) -> tuple[pd.DataFrame, str]:
    frame = harmonized.copy()
    notes = ""
    if variant == "active_all":
        notes = "Original iteration-3 active rows."
    elif variant == "exclude_HYP003":
        frame = frame.loc[~frame["study_id"].eq("HYP003")].copy()
        notes = "Drops provisional Bioconductor TFEA.ChIP source."
    elif variant == "exclude_HYP006":
        frame = frame.loc[~frame["study_id"].eq("HYP006")].copy()
        notes = "Drops Bauer acute hypoxia contrast."
    elif variant == "exclude_HYP007":
        frame = frame.loc[~frame["study_id"].eq("HYP007")].copy()
        notes = "Drops Bauer chronic hypoxia contrast."
    elif variant == "exclude_bauer":
        frame = frame.loc[~frame["study_id"].isin(["HYP006", "HYP007"])].copy()
        notes = "Drops both Bauer contrasts."
    elif variant == "exclude_HYP008":
        frame = frame.loc[~frame["study_id"].eq("HYP008")].copy()
        notes = "Drops hDASMC source with sign-flipped LFC."
    elif variant == "exclude_padj_surrogates":
        frame = frame.loc[~frame["study_id"].isin(["HYP001", "HYP006", "HYP007"])].copy()
        notes = "Drops sources where nominal p-values are unavailable and padj is used as p surrogate."
    elif variant == "positive_consensus_only":
        notes = "Post-filters consensus to positive Stouffer z before recall; diagnostic only."
    else:
        raise ValueError(f"Unknown diagnostic variant: {variant}")
    return frame, notes


def sensitivity_metrics(harmonized: pd.DataFrame, min_studies: int = 2) -> pd.DataFrame:
    variants = [
        "active_all",
        "exclude_HYP003",
        "exclude_HYP006",
        "exclude_HYP007",
        "exclude_bauer",
        "exclude_HYP008",
        "exclude_padj_surrogates",
        "positive_consensus_only",
    ]

    records: list[dict[str, Any]] = []
    for variant in variants:
        frame, notes = _variant_frame(harmonized, variant)
        consensus = slice_consensus(frame, min_studies=min_studies)
        if variant == "positive_consensus_only" and not consensus.empty:
            consensus = _rerank(consensus.loc[consensus["stouffer_z"].gt(0)].copy())

        recall50 = _recall(consensus, 50)
        recall100 = _recall(consensus, 100)
        records.append(
            {
                "variant": variant,
                "n_input_rows": int(len(frame)),
                "n_studies": int(frame["study_id"].nunique()) if not frame.empty else 0,
                "n_consensus_genes": int(len(consensus)),
                "recall_at_50": recall50["recall"],
                "n_recovered_at_50": recall50["n_recovered"],
                "recovered_at_50": ";".join(recall50["recovered"]),
                "recall_at_100": recall100["recall"],
                "n_recovered_at_100": recall100["n_recovered"],
                "recovered_at_100": ";".join(recall100["recovered"]),
                "notes": notes,
            }
        )
    return pd.DataFrame.from_records(records)


def write_diagnostics(catalog_path: Path, harmonized_path: Path, output_dir: Path, min_studies: int = 2) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    catalog = active_catalog(catalog_path)
    harmonized = pd.read_csv(harmonized_path, low_memory=False)

    support = target_support_by_study(harmonized, catalog)
    orientation = orientation_audit(support)
    sensitivity = sensitivity_metrics(harmonized, min_studies=min_studies)

    outputs = {
        "target_support_by_study": output_dir / "target_support_by_study.csv",
        "orientation_audit": output_dir / "orientation_audit.csv",
        "sensitivity_metrics": output_dir / "sensitivity_metrics.csv",
    }
    support.to_csv(outputs["target_support_by_study"], index=False)
    orientation.to_csv(outputs["orientation_audit"], index=False)
    sensitivity.to_csv(outputs["sensitivity_metrics"], index=False)

    command = (
        "make -C outputs/code diagnose"
        f" CATALOG={catalog_path.resolve()}"
        f" DIAG_HARMONIZED={harmonized_path.resolve()}"
        f" DIAG_OUTDIR={output_dir.resolve()}"
    )
    for path in outputs.values():
        write_source_sidecar(path, command, inputs=[catalog_path, harmonized_path], metadata={"generator": "diagnostics"})

    summary = {
        "catalog_path": str(catalog_path.resolve()),
        "harmonized_path": str(harmonized_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "n_active_studies": int(len(catalog)),
        "n_target_support_rows": int(len(support)),
        "orientation_positive_fraction_by_study": orientation.groupby("study_id")["orientation_ok"].mean().to_dict(),
        "sensitivity": sensitivity.to_dict(orient="records"),
    }
    summary_path = output_dir / "diagnostic_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(summary_path, command, inputs=[catalog_path, harmonized_path], metadata={"generator": "diagnostics"})
    return summary
