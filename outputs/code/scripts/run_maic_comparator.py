#!/usr/bin/env python
"""Run MAIC (Meta-Analysis by Information Content) as a same-input comparator.

MAIC (Wang et al., 2022, Bioinformatics; https://github.com/baillielab/maic) aggregates
heterogeneous ranked/unranked gene *lists* by their information content. To compare it
fairly against DEGORA and the other summary-level comparators, this script feeds MAIC the
identical per-study ranked gene lists used by the RobustRankAggreg adapter (study-level,
ordered by within-study normalized rank), parses MAIC's ranked output, and writes a
baseline-format result TSV so MAIC appears alongside the other comparators in the gold
comparator summaries.

MAIC consumes only list membership/rank -- it has no notion of direction or effect size --
so the emitted ``direction`` column is intentionally empty, which makes MAIC's
direction-aware recall zero by construction (a faithful capability difference, not a bug).
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from degora.baselines import BASELINE_RESULT_COLUMNS, _study_level
from degora.provenance import shell_command, write_source_sidecar

MAIC_DIR_DEFAULT = Path("outputs/code/vendor/maic")
MAIC_COMMIT = "fb26061"
MAIC_URL = "https://github.com/baillielab/maic"


def _filter_significant(harmonized: pd.DataFrame, padj_threshold: float = 0.05) -> pd.DataFrame:
    """Keep only per-study significant DEG rows (padj < threshold; fallback pvalue < threshold).

    MAIC is designed for curated/significant gene LISTS; its information-content scoring
    saturates (ties) on full-transcriptome rankings. Restricting each study to its reported
    significant DEGs gives MAIC its native input. See Methods.
    """
    padj = pd.to_numeric(harmonized["padj"], errors="coerce") if "padj" in harmonized.columns else None
    pval = pd.to_numeric(harmonized["pvalue"], errors="coerce") if "pvalue" in harmonized.columns else None
    if padj is None and pval is None:
        raise ValueError("harmonized table lacks both padj and pvalue; cannot define significant lists")
    if padj is None:
        mask = pval < padj_threshold
    elif pval is None:
        mask = padj < padj_threshold
    else:
        mask = (padj < padj_threshold) | (padj.isna() & (pval < padj_threshold))
    return harmonized.loc[mask.fillna(False)].copy()


def _build_maic_input(
    harmonized: pd.DataFrame,
    min_studies: int,
    category: str,
    significant_only: bool = True,
    padj_threshold: float = 0.05,
) -> tuple[str, int, int]:
    """Build a MAIC input file: one RANKED list per study, ordered by normalized rank.

    With ``significant_only`` (default), each study contributes only its significant DEGs
    (padj/pvalue < ``padj_threshold``) -- MAIC's native curated-gene-list regime. The per-study
    rank ordering otherwise mirrors degora.baselines.robustrankaggreg_adapter.
    """
    if significant_only:
        harmonized = _filter_significant(harmonized, padj_threshold)
    by_study = _study_level(harmonized, min_studies)
    required = {"study_id", "gene_symbol", "normalized_rank"}
    missing = required.difference(by_study.columns)
    if missing:
        raise ValueError(f"MAIC input missing required columns: {sorted(missing)}")
    ranked = by_study.loc[:, ["study_id", "gene_symbol", "normalized_rank"]].dropna(
        subset=["study_id", "gene_symbol", "normalized_rank"]
    )
    ranked = ranked.sort_values(["study_id", "normalized_rank", "gene_symbol"]).drop_duplicates(
        ["study_id", "gene_symbol"], keep="first"
    )
    if ranked.empty:
        raise ValueError("MAIC input is empty after study-level filtering")
    lines: list[str] = []
    for study_id, group in ranked.groupby("study_id", sort=True):
        genes = group["gene_symbol"].astype(str).str.strip().str.upper().tolist()
        # MAIC line format: Category <tab> ListLabel <tab> RANKED <tab> NAMED_GENES <tab> gene1 <tab> gene2 ...
        lines.append("\t".join([str(category), str(study_id), "RANKED", "NAMED_GENES", *genes]))
    n_universe = int(ranked["gene_symbol"].nunique())
    n_lists = int(ranked["study_id"].nunique())
    return "\n".join(lines) + "\n", n_universe, n_lists


def _run_maic(maic_dir: Path, input_path: Path, out_dir: Path) -> Path:
    """Run the MAIC CLI and return the path to maic_raw.txt.

    MAIC appends a timestamp to the output folder if it already exists, so the caller must
    pass a non-existent out_dir; we then glob for the produced folder.
    """
    cmd = [sys.executable, "maic.py", "-f", str(input_path.resolve()), "-o", str(out_dir.resolve()), "-q"]
    proc = subprocess.run(cmd, cwd=str(maic_dir), capture_output=True, text=True, timeout=3600)
    candidates = list(out_dir.parent.glob(out_dir.name + "*/maic_raw.txt"))
    if not candidates:
        raise RuntimeError(
            f"MAIC produced no maic_raw.txt (exit {proc.returncode}).\nstderr tail:\n{proc.stderr[-2000:]}"
        )
    return sorted(candidates)[-1]


def _parse_maic(raw_path: Path) -> list[dict[str, object]]:
    """Parse maic_raw.txt; genes are already in descending maic_score order."""
    with open(raw_path) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0]
    score_idx = header.index("maic_score")
    out: list[dict[str, object]] = []
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        out.append({"symbol": row[0].strip().upper(), "score": float(row[score_idx])})
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--min-studies", type=int, default=2)
    parser.add_argument("--maic-dir", type=Path, default=MAIC_DIR_DEFAULT)
    parser.add_argument(
        "--full-lists",
        action="store_true",
        help="Feed MAIC full per-study rankings instead of significant-only DEG lists "
        "(MAIC degenerates on full rankings; default is its native significant-gene-list input).",
    )
    parser.add_argument("--padj-threshold", type=float, default=0.05)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not (args.maic_dir / "maic.py").exists():
        raise FileNotFoundError(
            f"MAIC not found at {args.maic_dir}. Clone it first: "
            f"git clone {MAIC_URL} {args.maic_dir} && (cd {args.maic_dir} && git checkout {MAIC_COMMIT})"
        )
    harmonized = pd.read_csv(args.harmonized, low_memory=False)
    significant_only = not args.full_lists
    text, n_universe, n_lists = _build_maic_input(
        harmonized, args.min_studies, args.corpus, significant_only=significant_only, padj_threshold=args.padj_threshold
    )
    list_mode = "full_rankings" if args.full_lists else f"significant_DEGs_padj<{args.padj_threshold}"

    work = args.baseline_dir / "_maic_work"
    work.mkdir(parents=True, exist_ok=True)
    input_path = work / f"{args.corpus}_maic_input.txt"
    input_path.write_text(text)
    out_dir = work / f"{args.corpus}_maic_out"
    for existing in work.glob(f"{args.corpus}_maic_out*"):
        if existing.is_dir():
            shutil.rmtree(existing)

    start = time.perf_counter()
    raw_path = _run_maic(args.maic_dir, input_path, out_dir)
    runtime = time.perf_counter() - start
    parsed = _parse_maic(raw_path)
    if not parsed:
        raise RuntimeError("MAIC returned no ranked genes")

    frame = pd.DataFrame(
        {
            "method_id": "maic",
            "setting_id": "default",
            "gene_id": [row["symbol"] for row in parsed],
            "symbol": [row["symbol"] for row in parsed],
            "rank": list(range(1, len(parsed) + 1)),
            "score": [row["score"] for row in parsed],
            "pvalue": "",
            "padj": "",
            "effect": "",
            "direction": "",  # MAIC aggregates list membership/rank only; no direction
            "n_studies": "",
            "missingness": 0.0,
            "runtime_s": runtime,
            "version": f"MAIC {MAIC_COMMIT} ({list_mode}; list-aggregation; no direction/effect)",
            "status": "ok",
        }
    )[BASELINE_RESULT_COLUMNS]

    args.baseline_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.baseline_dir / f"{args.corpus}_maic_default.tsv"
    frame.to_csv(out_tsv, sep="\t", index=False)

    # Register the MAIC result in baseline_manifest.csv (the authority used by
    # baseline_result_paths). Upsert so re-runs do not duplicate the row.
    manifest_path = args.baseline_dir / "baseline_manifest.csv"
    manifest_row = {
        "artifact": str(out_tsv.resolve()),
        "method_id": "maic",
        "setting_id": "default",
        "artifact_type": "baseline_result",
        "status": "ok",
        "source_command": "run_maic_comparator.py",
        "input_harmonized": str(args.harmonized.resolve()),
        "rows": len(parsed),
        "notes": "MAIC list-aggregation comparator (no direction); same per-study ranked lists as RobustRankAggreg",
    }
    if manifest_path.exists():
        manifest = pd.read_csv(manifest_path)
        keep = ~(
            manifest["method_id"].astype(str).eq("maic")
            & manifest["artifact_type"].astype(str).eq("baseline_result")
        )
        manifest = manifest.loc[keep]
        manifest = pd.concat([manifest, pd.DataFrame([manifest_row])], ignore_index=True, sort=False)
    else:
        manifest = pd.DataFrame([manifest_row])
    manifest.to_csv(manifest_path, index=False)

    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python",
            "outputs/code/scripts/run_maic_comparator.py",
            "--harmonized",
            str(args.harmonized),
            "--baseline-dir",
            str(args.baseline_dir),
            "--corpus",
            args.corpus,
            "--min-studies",
            str(args.min_studies),
        ]
    )
    write_source_sidecar(
        out_tsv,
        command,
        inputs=[args.harmonized],
        metadata={
            "generator": "maic-comparator",
            "maic_commit": MAIC_COMMIT,
            "maic_source": MAIC_URL,
            "n_lists": n_lists,
            "n_universe": n_universe,
            "list_mode": list_mode,
            "note": "MAIC native input = per-study significant-DEG gene lists; MAIC has no direction so direction-aware recall is 0 by construction.",
        },
    )
    print(
        f"MAIC comparator: ranked {len(parsed)} genes from {n_lists} study lists "
        f"(universe {n_universe}) in {runtime:.1f}s -> {out_tsv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
