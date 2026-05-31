# State

Current iteration: 58
Current phase: statistical reporting enhancement complete
Last updated: 2026-05-29

## Summary of where things stand

The `claude_stat1.txt` statistical enhancement review was accepted where it improved publication rigor without retuning the locked DEGORA score. The main score remains the fixed `degora_quality_weighted_score`; the new statistics are auxiliary reporting, uncertainty, and sensitivity layers.

The repository now has:

- exact binomial recall confidence intervals, top-k precision, hypergeometric enrichment p-values, and BH FDR in comparator summaries and manuscript Table 1;
- background-negative AUROC/AUPRC metrics in deep benchmark summaries;
- source-unit beta-binomial direction posterior fields;
- heterogeneity-aware RE-Stouffer reporting fields;
- beta order-statistic RRA rho/rank fields;
- auxiliary random-effects log2FC reporting fields with CI, tau2, I2, source count, and SE provenance;
- regenerated score DBs and score CSVs for manuscript-facing IFN, ER stress, heat shock, and hypoxia corpora;
- regenerated gold summaries, manuscript tables, figures, Atlas/resource package, and paper-package manifest;
- `environment.yml` and GitHub Actions CI for check plus paper-package regeneration;
- a manuscript package validator that fails on missing required artifacts, reference placeholders, stale hypoxia 0.80-style claims, or missing manuscript headline values.

## Current Benchmark Results

- Table 1 fixed quality-weighted DEGORA recall@10/@50/@100:
  - IFN: 0.35 / 0.85 / 0.90
  - ER stress (RNA-seq + microarray, Table 1 input mode): 0.39 / 0.72 / 0.833; RNA-seq-only primary 0.28 / 0.56 / 0.667
  - Heat shock: 0.38 / 0.69 / 0.75
  - Hypoxia: 0.15 / 0.65 / 0.75
- Table 1 fixed quality-weighted DEGORA recall@100 exact 95% CI / precision@100 / hypergeometric FDR@100:
  - IFN: 0.68-0.99 / 0.18 / 5.21e-35
  - ER stress: 0.59-0.96 / 0.15 / 8.01e-30
  - Heat shock: 0.48-0.93 / 0.12 / 8.09e-23
  - Hypoxia: 0.51-0.91 / 0.15 / 2.35e-34
- Table 2 cross-platform:
  - IFN RNA-only to RNA+microarray: 0.35/0.45/0.85/0.90 to 0.40/0.70/0.80/0.90
  - ER stress RNA-only to RNA+microarray: 0.28/0.39/0.56/0.667 to 0.39/0.50/0.72/0.833 (human-only, three RNA-seq source units)
  - Hypoxia RNA-only to RNA+microarray: 0.15/0.25/0.65/0.75 to 0.15/0.30/0.65/0.75
- Atlas: 6 corpora, 118,633 scored gene rows, 114 method rows, 23 prior-art/resource rows.

## Verification

- Targeted statistical/reporting tests: 10 passed.
- Broader targeted score, comparator, deep-metric, benchmark-stat, and manuscript-package tests: 24 passed.
- `make -C outputs/code paper`: regenerated public-summary/prior-art tables, gold summaries, manuscript tables, Figure 1, Figure 2, cross-platform figure, Atlas HTML/SQLite package, and paper-package manifest.
- Deep point metrics regenerated for IFN, ER stress, heat shock, and hypoxia. IFN/ER/heat-shock bootstrap rows were regenerated with 10 bootstrap repeats; hypoxia point metrics were regenerated with bootstrap disabled because the full corpus remains computationally heavy.
- `make -C outputs/code check`: 128 passed, one known SciPy precision warning in the duplicate-symbol Welch test.
- `git diff --check`: passed.
- `outputs/manuscript/degora_paper_validation.txt`: `failures=0`, `warnings=1`; the warning now fires on the manuscript's own pre-submission packaging language (dataset-level GEO citations + Zenodo DOI + reference-style finalization), not only on a single exact phrase.
- Reproducibility: a deep-metrics/benchmark-recall rename-fallout (`*_piper_advantage_metrics.csv` vs the `*_degora_*` name the figure step expects, and a missing `benchmark_recall_summary.csv`) is fixed, so `make -C outputs/code paper` now runs clean end-to-end from a wiped `outputs/figures/manuscript` in the full working tree. On a bare GitHub clone CI runs package checks only (`make -C outputs/code unit typecheck lint`); the full paper/figure rebuild requires the curated corpus shipped via the Zenodo archive (lean code mirror by design).

## What to do on resume

1. Run the final citation audit, especially dataset-level source-study citations, GEO accession citations, and package citation strings.
2. Prepare public archive/release metadata and DOI.
3. Do final claim-to-artifact and journal-format checks before submission.
4. Add down-regulated gold markers or expanded negative panels if direction-generalization claims become central.
5. Optionally optimize or long-run hypoxia DEGORA bootstrap only if bootstrap intervals become central to the manuscript.

see LAB_NOTEBOOK.md iteration 58 statistical reporting enhancement pass
