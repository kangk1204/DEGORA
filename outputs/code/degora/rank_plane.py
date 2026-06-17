"""Rank-plane diagnostics over harmonized DEG tables."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .provenance import write_source_sidecar


def _rank_strength(values: pd.Series, *, ascending: bool, universe_size: int | None = None) -> pd.Series:
    """Return percentile-like rank strength where 1.0 is strongest."""

    ranks = values.rank(method="average", ascending=ascending)
    n_observed = int(values.notna().sum())
    n_universe = max(n_observed, int(universe_size or n_observed))
    if n_universe <= 1:
        return pd.Series(np.ones(len(values)), index=values.index, dtype=float)
    return 1.0 - ((ranks - 1.0) / (n_universe - 1.0))


def _collapse_study_gene(harmonized: pd.DataFrame) -> pd.DataFrame:
    required = {"study_id", "gene_symbol", "pvalue", "lfc"}
    missing = sorted(required.difference(harmonized.columns))
    if missing:
        raise ValueError(f"harmonized table is missing required rank-plane columns: {missing}")

    frame = harmonized.dropna(subset=["study_id", "gene_symbol", "pvalue", "lfc"]).copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.upper()
    frame["pvalue"] = pd.to_numeric(frame["pvalue"], errors="coerce")
    frame["lfc"] = pd.to_numeric(frame["lfc"], errors="coerce")
    frame = frame.dropna(subset=["gene_symbol", "pvalue", "lfc"])

    aggregations: dict[str, Any] = {
        "pvalue": ("pvalue", "min"),
        "lfc": ("lfc", "mean"),
    }
    for column in ["paper_id", "pipeline", "cell_system", "source_path"]:
        if column in frame.columns:
            aggregations[column] = (column, "first")
    if "n_genes_in_study" in frame.columns:
        frame["n_genes_in_study"] = pd.to_numeric(frame["n_genes_in_study"], errors="coerce")
        aggregations["n_genes_in_study"] = ("n_genes_in_study", "max")

    return frame.groupby(["study_id", "gene_symbol"], as_index=False).agg(**aggregations)


def rank_plane_points(harmonized: pd.DataFrame) -> pd.DataFrame:
    """Map study-gene rows into p-value rank and signed effect-rank space."""

    collapsed = _collapse_study_gene(harmonized)
    if collapsed.empty:
        return pd.DataFrame(
            columns=[
                "study_id",
                "gene_symbol",
                "pvalue",
                "lfc",
                "p_rank_strength",
                "effect_rank_strength",
                "signed_effect_rank",
                "rank_plane_delta",
                "n_genes_in_rank_universe",
            ]
        )

    pieces = []
    for _, study in collapsed.groupby("study_id", sort=False):
        study = study.copy()
        rank_universe = int(len(study))
        if "n_genes_in_study" in study.columns:
            declared = pd.to_numeric(study["n_genes_in_study"], errors="coerce").dropna()
            if not declared.empty and declared.max() > 0:
                rank_universe = int(max(rank_universe, round(float(declared.max()))))
        study["p_rank_strength"] = _rank_strength(study["pvalue"], ascending=True, universe_size=rank_universe)
        study["effect_rank_strength"] = _rank_strength(study["lfc"].abs(), ascending=False, universe_size=rank_universe)
        study["signed_effect_rank"] = np.sign(study["lfc"]) * study["effect_rank_strength"]
        study["rank_plane_delta"] = study["p_rank_strength"] - study["effect_rank_strength"]
        study["n_genes_in_rank_universe"] = rank_universe
        pieces.append(study)

    out = pd.concat(pieces, ignore_index=True)
    return out.sort_values(
        ["study_id", "p_rank_strength", "effect_rank_strength", "gene_symbol"],
        ascending=[True, False, False, True],
    ).reset_index(drop=True)


def rank_plane_gene_summary(points: pd.DataFrame, *, joint_threshold: float = 0.9) -> pd.DataFrame:
    """Summarize rank-plane support for each gene across studies."""

    if points.empty:
        return pd.DataFrame(
            columns=[
                "gene_symbol",
                "n_studies",
                "median_p_rank_strength",
                "median_effect_rank_strength",
                "median_signed_effect_rank",
                "median_abs_rank_delta",
                "joint_high_fraction",
                "effect_sign_concordance",
                "rank_plane_score",
            ]
        )

    frame = points.copy()
    frame["_joint_high"] = frame["p_rank_strength"].ge(joint_threshold) & frame["effect_rank_strength"].ge(joint_threshold)
    consensus_sign = frame.groupby("gene_symbol")["signed_effect_rank"].median().rename("_median_signed")
    frame = frame.merge(consensus_sign, on="gene_symbol", how="left")
    frame["_consensus_sign"] = np.sign(frame["_median_signed"])
    frame["_sign_concordant"] = np.where(
        frame["_consensus_sign"].ne(0),
        np.sign(frame["signed_effect_rank"]).eq(frame["_consensus_sign"]),
        False,
    )

    out = frame.groupby("gene_symbol", as_index=False).agg(
        n_studies=("study_id", "nunique"),
        median_p_rank_strength=("p_rank_strength", "median"),
        median_effect_rank_strength=("effect_rank_strength", "median"),
        median_signed_effect_rank=("signed_effect_rank", "median"),
        median_abs_rank_delta=("rank_plane_delta", lambda value: float(np.median(np.abs(value)))),
        joint_high_fraction=("_joint_high", "mean"),
        effect_sign_concordance=("_sign_concordant", "mean"),
    )
    out["rank_plane_score"] = (
        out["median_p_rank_strength"]
        * out["median_effect_rank_strength"]
        * out["effect_sign_concordance"]
    )
    return out.sort_values(
        ["rank_plane_score", "joint_high_fraction", "gene_symbol"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def rank_plane_study_summary(points: pd.DataFrame, *, joint_threshold: float = 0.9) -> pd.DataFrame:
    """Summarize within-study p-rank versus effect-rank concordance."""

    if points.empty:
        return pd.DataFrame(columns=["study_id", "n_genes", "spearman_p_vs_effect_rank", "joint_high_count"])

    records = []
    for study_id, study in points.groupby("study_id", sort=False):
        records.append(
            {
                "study_id": study_id,
                "n_genes": int(len(study)),
                "spearman_p_vs_effect_rank": float(
                    study["p_rank_strength"].corr(study["effect_rank_strength"], method="spearman")
                ),
                "joint_high_count": int(
                    (
                        study["p_rank_strength"].ge(joint_threshold)
                        & study["effect_rank_strength"].ge(joint_threshold)
                    ).sum()
                ),
            }
        )
    return pd.DataFrame.from_records(records).sort_values("study_id").reset_index(drop=True)


def write_rank_plane_outputs(
    harmonized_path: Path,
    output_dir: Path,
    *,
    joint_threshold: float = 0.9,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    harmonized = pd.read_csv(harmonized_path, low_memory=False)
    points = rank_plane_points(harmonized)
    gene_summary = rank_plane_gene_summary(points, joint_threshold=joint_threshold)
    study_summary = rank_plane_study_summary(points, joint_threshold=joint_threshold)

    outputs = {
        "rank_plane_points": output_dir / "rank_plane_points.csv",
        "rank_plane_gene_summary": output_dir / "rank_plane_gene_summary.csv",
        "rank_plane_study_summary": output_dir / "rank_plane_study_summary.csv",
    }
    points.to_csv(outputs["rank_plane_points"], index=False)
    gene_summary.to_csv(outputs["rank_plane_gene_summary"], index=False)
    study_summary.to_csv(outputs["rank_plane_study_summary"], index=False)

    command = (
        "make -C outputs/code rank-plane"
        f" HARMONIZED={harmonized_path.resolve()}"
        f" OUTDIR={output_dir.resolve()}"
        f" RANK_PLANE_THRESHOLD={joint_threshold}"
    )
    for path in outputs.values():
        write_source_sidecar(path, command, inputs=[harmonized_path], metadata={"generator": "rank-plane"})

    summary = {
        "harmonized_path": str(harmonized_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "joint_threshold": joint_threshold,
        "n_rank_plane_points": int(len(points)),
        "n_gene_summary_rows": int(len(gene_summary)),
        "n_studies": int(points["study_id"].nunique()) if not points.empty else 0,
        "top_rank_plane_genes": gene_summary.head(20)["gene_symbol"].tolist(),
    }
    summary_path = output_dir / "rank_plane_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_source_sidecar(summary_path, command, inputs=[harmonized_path], metadata={"generator": "rank-plane"})
    return summary
