# Plan

## Change log

Iteration 58 reviewed `claude_stat1.txt` and implemented the statistically useful recommendations as auxiliary reporting lanes while keeping the locked `degora_quality_weighted_score` unchanged. Comparator summaries, deep benchmark metrics, manuscript Table 1, the supplementary table index, score DBs, manuscript text, figures, and the Atlas package were regenerated.

## Current Conclusion

DEGORA remains manuscript-viable as an auditable public-DEG evidence prioritization and local database tool. The current package is stronger than iteration 57 because benchmark recovery now carries exact uncertainty, null-enrichment/FDR context, background-negative AUROC/AUPRC, and auxiliary per-gene sensitivity lanes without changing the score that generated the manuscript claims.

Fixed claims to preserve:

1. quality-weighted DEGORA is the single manuscript-facing default;
2. primary benchmark outputs are human-only and current with regenerated tables/figures;
3. cross-platform RNA-seq plus human microarray evidence is supported with explicit source-quality and fallback limitations;
4. runnable same-input baselines are compared conservatively, while different-input tools are classified rather than treated as defeated;
5. `metarnaseq_invnorm` is reported as runnable only when the current corpus actually emits informative output; sparse-public-summary blocked cases are blocked rows;
6. high heterogeneity is a review signal, not a disqualifying gate or inferential claim;
7. the local HTML/SQLite Atlas is a concrete resource output and should be emphasized as part of the method's practical value.
8. RE-Stouffer, RRA, direction posterior, and random-effects log2FC outputs are auxiliary reporting lanes, not replacements for the locked score or calibrated meta-analysis.

## Immediate Next Tasks

1. Citation audit:
   - source-study and GEO accession citations for each active benchmark source;
   - package/version citation strings for comparator and plotting packages;
   - final Scientific Reports reference style.
2. Claim-to-artifact consistency audit:
   - Abstract, Results, Table 1, Table 2, Figure 2, Figure 3, Figure 4, source sidecars, and supplementary tables;
   - ensure no stale hypoxia 0.80 claim or stale/deprecated corpus claim remains in manuscript-facing files.
3. Release preparation:
   - public repository cleanup;
   - archive code, catalogs, harmonized tables, generated figures, and source data;
   - DOI and final Data/Code Availability wording.
4. Optional statistical strengthening:
   - add down-regulated gold markers and expanded negative panels if direction-generalization claims are emphasized;
   - optimize or long-run hypoxia DEGORA bootstrap if confidence intervals become part of the central claim.

## Verification Baseline

- `make -C outputs/code paper`: completed in iteration 58; produced gold summaries, manuscript tables, figures, Atlas HTML/SQLite package, and paper-package validation artifacts.
- `make -C outputs/code check`: 128 passed; one known SciPy precision warning. Bare-clone CI runs package checks only (`make unit typecheck lint`); the full `paper`/`figs`/`provenance-check` chain needs the curated corpus + generated intermediates (gitignored here, shipped via Zenodo).
- `git diff --check`: passed.
- Current Table 1 headline values: IFN 0.90, ER stress 0.833 (RNA-seq + microarray input mode; RNA-seq-only primary 0.667), heat shock 0.75, hypoxia 0.75 recall@100 for fixed quality-weighted DEGORA.
- Current Table 1 recall@100 exact 95% CI values: IFN 0.68-0.99, ER stress 0.59-0.96, heat shock 0.48-0.93, hypoxia 0.51-0.91.
- Current Table 1 hypergeometric FDR@100 values: IFN 5.21e-35, ER stress 8.01e-30, heat shock 8.09e-23, hypoxia 2.35e-34.
- Current Table 2 cross-platform values match manuscript text.
- Paper package validation has zero failures and one known pending-citation warning.
