# Biostatinfo Expert Review - 2026-05-28

## Scope

Reviewed the implemented PIPER-DEG scoring, harmonization, derived-count benchmark generation, score DB output, source-quality weighting, and gold-panel comparator/figure outputs.

## Estimand

PIPER-DEG should be framed as a gene-prioritization and evidence-DB method. The target estimand is not a calibrated meta-analysis p-value. The actionable signal is whether a gene shows directionally consistent expression change across independent source units for a biological keyword.

## Findings And Fixes

1. Count-derived RNA-seq inputs lacked an explicit low-count expression filter.
   Fix: added a treatment-label-independent raw-count filter before logCPM/Welch derivation. The default keeps genes with raw count >= 10 in at least 2 selected samples and records filter metadata in each derived table.

2. Derived count benchmarks included non-protein-coding rows that could occupy top ranks despite protein-coding gene-level gold panels.
   Fix: restricted the IFN, ER stress, and heat-shock derived count universes to protein-coding symbols/IDs through HGNC or NCBI annotation before deriving contrasts.

3. Semicolon-joined `table_scope` labels were treated as unknown quality.
   Fix: source-quality parsing now uses the most conservative multiplier among joined scopes. For example, `deg_only;full_results` is scored no better than `deg_only`.

4. Score evidence rows did not expose `table_scope`.
   Fix: `gene_evidence` now includes `table_scope`, making the field used by source-quality weighting auditable in the browser/API database.

5. Gold comparator summaries and recall figures could mix stale baseline TSVs with current outputs.
   Fix: baseline readers now use `baseline_manifest.csv` as the authority when present and raise on malformed/missing manifest artifacts. Comparator summaries and recall figures no longer glob every TSV in the directory when a manifest exists.

## Current Benchmark Readout

| Corpus | PIPER top-10 recall | PIPER top-50 recall | PIPER top-100 recall | Main conclusion |
| --- | ---: | ---: | ---: | --- |
| IFN | 0.35 | 0.85 | 0.90 | Strong top-rank recovery of canonical ISGs; PIPER improves early recall over runnable classical baselines at K=10/50. |
| ER stress primary | 0.39 | 0.72 | 0.72 | Primary-quality ER sources support PIPER as an early-prioritization method. |
| ER stress full | 0.11 | 0.11 | 0.17 | Full noisy sensitivity is unfavorable for the primary score; quality-weighted secondary score recovers to 0.28/0.61/0.61. |
| Heat shock | 0.38 | 0.69 | 0.75 | Heat-shock result is biologically coherent; top genes are dominated by HSPA/HSPH/DNAJ/BAG3 heat-shock markers. |

## Statistical Interpretation

The strongest manuscript claim is not universal superiority over every meta-analysis method. The defensible claim is that PIPER-DEG creates an auditable, source-unit-aware directional evidence database that works well for keyword-level biological signal recovery, especially where users care about consistent up/down direction across heterogeneous DEG tables and count-derived sources.

The ER-full result must be reported as a limitation and as evidence for the quality-weighted secondary ranking. The full set contains a lower-confidence normalized-matrix source that strongly changes the ranking. Hiding that would overstate the method; reporting it strengthens the paper's integrity.

## Remaining Risks

- hStouffer and AWmeta remain blocked because faithful inputs with original variance/SE fields are not available for all datasets. These are blocked comparators, not defeated baselines.
- Derived-count Welch/logCPM tables are pragmatic benchmark inputs, not substitutes for full DESeq2/edgeR re-analysis.
- Gold-panel recall measures positive recovery and direction consistency, not false discovery control.
- ER-full primary performance is weak and should stay in sensitivity/limitation framing.

## Verification

- Targeted tests: `24 passed, 1 warning`.
- Regenerated derived IFN, ER stress, and heat-shock DEG tables with low-count/protein-coding filters.
- Regenerated harmonized slices, score DBs, baseline outputs, gold comparator summaries, source-support reports, ER benchmark report, and heat-shock top-1000 recall figure.
- Comparator/figure stale-output regression tests added.
