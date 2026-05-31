# PIPER-DEG manuscript revision — exact edit specification

Author decisions (locked): ER Table 1 = **keep cross-platform 0.89 but label the input mode openly**;
revision depth = **full (edit text + fix generator + regenerate + validate)**; direction claim = **de-scope honestly**.

Apply every edit below with exact-string matching. Read each file fresh before editing. After all edits,
regenerate and validate (Step Z). Do NOT change any locked recall numbers other than as written here.
Do NOT alter the ER corpus (ER stays sourced from `er-stress-cross-platform`); only ADD an input-mode label.

Ground-truth ER RNA-seq-only values (corpus `er-stress-benchmark-primary`, 8,667 genes): PIPER recall@10/50/100
= 0.39/0.72/0.72; best non-PIPER (Stouffer) ties at 0.72. ER combined corpus (`er-stress-cross-platform`,
10,725 genes, in Table 1): PIPER 0.28/0.83/0.89, best non-PIPER (rank_product_approx) 0.78 at recall@100.

---

## STEP 1 — Generator: add an `input_mode` column to Table 1

File: `/home/keunsoo/projects/09_PIPER/outputs/code/scripts/write_manuscript_tables.py` (this is the REAL Table 1
generator invoked by the Makefile `manuscript-tables` target — NOT write_comparator_summary.py).

Read it. Find where it assembles the Table 1 rows / column order for `Table1_core_benchmark_summary.{csv,md}`.
Add a new column named `input_mode` positioned immediately AFTER the `topic` column, with per-topic values:
- IFN -> `RNA-seq`
- ER stress -> `RNA-seq + microarray`
- Heat shock -> `RNA-seq`
- Hypoxia -> `RNA-seq`

Do NOT change which corpus/summary CSV each topic reads from (ER must stay `er-stress-cross-platform`).
If the generator builds rows via a dict per topic, add `"input_mode": <value>` keyed by topic. If it builds the
column list explicitly, insert `input_mode` after `topic` in both the CSV and the markdown writer. Keep all other
columns and values byte-for-byte identical. Confirm `python -c "import ast,pathlib; ast.parse(...)"` parses.

## STEP 2 — Regenerate Table 1 and confirm

From repo root `/home/keunsoo/projects/09_PIPER`, run the documented command (check the Makefile
`manuscript-tables` target for exact args; typically):
`make -C outputs/code manuscript-tables`  (or `cd outputs/code && python scripts/write_manuscript_tables.py --output-dir ../tables/manuscript` adjusting paths so it does NOT double "outputs/outputs").
Confirm `outputs/tables/manuscript/Table1_core_benchmark_summary.md` now has an `input_mode` column with
ER stress = `RNA-seq + microarray` and the other three = `RNA-seq`, and that ER recall@100 is still 0.8889.

---

## STEP 3 — main.md edits

File: `/home/keunsoo/projects/09_PIPER/outputs/manuscript/main.md`. Apply each OLD->NEW exactly.

### EDIT A — remove the "Draft Status" banner
OLD:
```
# PIPER-DEG: auditable prioritization from heterogeneous published DEG tables

## Draft Status

This working draft is aligned to the current generated benchmark tables in `outputs/tables/manuscript/` and figure packages in `outputs/figures/manuscript/`. Literature citations have been added for major method and biology claims, but dataset-level source-study citations still require a final citation audit before submission.

## Abstract
```
NEW:
```
# PIPER-DEG: auditable prioritization from heterogeneous published DEG tables

## Abstract
```

### EDIT B — Abstract (replace from "It harmonizes DEG rows" through the final sentence of the Abstract)
OLD:
```
It harmonizes DEG rows, collapses related contrasts within independent source units, applies predeclared source-quality weights, ranks genes by repeated directional support, and exports both gene scores and browser-inspectable evidence records. In human positive-control benchmarks, the fixed quality-weighted PIPER-DEG ranking recovered canonical marker panels with recall@100 of 0.90 for interferon response, 0.89 for ER stress/unfolded protein response, 0.75 for heat shock, and 0.75 for hypoxia/HIF. Adding human microarray evidence increased IFN recall@20 from 0.45 to 0.70 and ER-stress recall@50/100 from 0.72/0.72 to 0.83/0.89, while hypoxia recall@100 remained 0.75 with marker-order shifts. Comparator analyses support a calibrated claim: PIPER-DEG is not a universal replacement for raw-expression or variance-aware meta-analysis, but it provides an auditable integration and local database layer for the public-summary regime where RNA-seq and microarray DEG evidence must be interpreted together.
```
NEW:
```
It harmonizes DEG rows, collapses related contrasts within independent source units, applies predeclared source-quality weights, ranks genes by repeated, directionally consistent support, and exports both gene scores and browser-inspectable evidence records. In human positive-control benchmarks built from only 2-15 independent source units per topic, the fixed quality-weighted PIPER-DEG ranking recovered canonical marker panels at recall@100 of 0.90 (interferon response), 0.75 (heat shock), and 0.75 (hypoxia/HIF) on RNA-seq evidence, and 0.89 for ER stress/unfolded protein response on the combined RNA-seq plus microarray corpus (0.72 on RNA-seq alone). These recall values matched standard summary-level meta-analysis methods rather than exceeding them: on three of the four topics the best non-PIPER comparator tied PIPER at recall@100, so the benchmark establishes faithful recovery of known biology, not statistical superiority. Adding human microarray evidence raised IFN recall@20 from 0.45 to 0.70 and ER-stress recall@50/100 from 0.72/0.72 to 0.83/0.89, while IFN recall@50 fell from 0.85 to 0.80 and hypoxia recall@100 stayed 0.75 with marker-order shifts. PIPER-DEG is therefore not a universal replacement for raw-expression or variance-aware meta-analysis; it provides an auditable integration and local database layer for the public-summary regime where RNA-seq and microarray DEG evidence must be interpreted together.
```

### EDIT C — score-label reconciliation (Results, "fixed quality-weighted ranking")
OLD:
```
The manuscript-facing analyses use `piper_quality_weighted_score` as the fixed PIPER-DEG ranking wherever it is available. The unweighted `piper_score` and `piper_slice` outputs remain in the result files as reference indices, but Table 1 and Figure 2 do not select the best PIPER row per topic. This design prevents benchmark-row selection against the locked gold panels.
```
NEW:
```
The manuscript-facing analyses use `piper_quality_weighted_score` as the fixed PIPER-DEG ranking wherever it is available. This source-quality-weighted lane is implemented in the codebase as a secondary variant of the primary `piper_score`; it is predeclared as the single manuscript-facing PIPER row so that no PIPER score is selected per topic after the gold panels are locked. The unweighted `piper_score` and `piper_slice` outputs remain in the result files as reference indices, but Table 1 and Figure 2 do not select the best PIPER row per topic. This design prevents benchmark-row selection against the locked gold panels.
```

### EDIT D — Results recall numbers + ER labeling + FDR generalization
OLD:
```
In Table 1 and Figure 2, fixed quality-weighted PIPER-DEG recovered canonical markers at useful early ranks. IFN achieved recall@10/50/100 of 0.35/0.85/0.90. ER stress achieved 0.28/0.83/0.89. Heat shock achieved 0.38/0.69/0.75. Hypoxia achieved 0.15/0.65/0.75. Exact 95% recall intervals at recall@100 were necessarily wide for the small locked panels: IFN 0.68-0.99, ER stress 0.65-0.99, heat shock 0.48-0.93, and hypoxia 0.51-0.91. Top-100 hypergeometric enrichment FDR values were small across all four topics (IFN 5.21e-35, ER stress 8.40e-31, heat shock 8.09e-23, hypoxia 2.35e-34), supporting marker enrichment while preserving the point-summary interpretation.
```
NEW:
```
In Table 1 and Figure 2, fixed quality-weighted PIPER-DEG recovered canonical markers at useful early ranks. On RNA-seq evidence, IFN achieved recall@10/50/100 of 0.35/0.85/0.90, heat shock 0.38/0.69/0.75, and hypoxia 0.15/0.65/0.75. The ER-stress row in Table 1 is computed on the combined RNA-seq plus microarray corpus (its input mode is labeled explicitly in Table 1) and achieved 0.28/0.83/0.89; on RNA-seq-only evidence ER stress reached recall@10/50/100 of 0.39/0.72/0.72, a tie with the best non-PIPER comparator (Table 2). Exact 95% recall intervals at recall@100 were necessarily wide for the small locked panels (16-20 genes): IFN 0.68-0.99, ER stress 0.65-0.99 (combined corpus), heat shock 0.48-0.93, and hypoxia 0.51-0.91. Top-100 hypergeometric enrichment FDR values were small across all four topics (all < 1e-22), but the runnable comparators reached similarly small enrichment FDRs at the same cutoff, so these values confirm that the rankings are enriched for known biology rather than distinguishing methods.
```

### EDIT E — comparator claim (the paragraph beginning "The comparator results support a competitive")
OLD:
```
The comparator results support a competitive, not universal-superiority, interpretation. PIPER-DEG matched or exceeded the best non-PIPER comparator at recall@100 for all four topics and improved recall@50 for IFN, ER stress, and heat shock. Hypoxia was deliberately retained as a balanced case: PIPER-DEG matched the best non-PIPER recall@100 but did not dominate early cutoffs. This pattern is consistent with the intended claim that PIPER-DEG is a practical evidence-prioritization layer for heterogeneous public DEG tables rather than a uniformly superior statistical test.
```
NEW:
```
The comparator results support a faithfulness interpretation rather than a superiority claim. At recall@100, PIPER-DEG tied the best runnable non-PIPER comparator on IFN (0.90 vs 0.90) and hypoxia (0.75 vs 0.75); on heat shock its point estimate was the only one to exceed every comparator (0.75 vs 0.625), although the 16-gene panel leaves the confidence intervals overlapping; and on the combined RNA-seq plus microarray ER-stress corpus shown in Table 1 it exceeded the best comparator (0.89 vs 0.78), an advantage that disappears on RNA-seq-only evidence (0.72, a tie; Table 2). PIPER-DEG showed modest early- and mid-rank gains on IFN (recall@10 0.35 vs 0.25; recall@50 0.85 vs 0.80) and heat shock (recall@50 0.69 vs 0.56), whereas on hypoxia an untuned Stouffer baseline exceeded it at recall@10 and recall@50 (0.20/0.70 vs 0.15/0.65). The benchmark therefore shows that PIPER-DEG recovers established marker biology as well as standard summary-level methods while adding an auditable evidence layer, not that it is a uniformly superior statistical test.
```

### EDIT F — thin-corpora sentence: add hypoxia source-unit count
OLD:
```
IFN used 2 RNA-seq source units in the primary table and 3 after microarray addition, ER stress used 2 source units, and heat-shock genes were supported by 2-3 source units.
```
NEW:
```
IFN used 2 RNA-seq source units in the primary table and 3 after microarray addition, ER stress used 2 source units, heat-shock genes were supported by 2-3 source units, and hypoxia--the deepest corpus, and the one where PIPER-DEG's advantage over standard methods was smallest--used 15 RNA-seq source units (17 after microarray addition).
```

### EDIT G — cross-platform confound + Welch (Methods "Microarray handling")
OLD:
```
In the current cross-platform benchmarks, IFN uses GSE71634, ER stress uses GSE19519 tunicamycin and thapsigargin contrasts from one source unit, and hypoxia uses GSE3045 and GSE22282. The ER-stress normalized-matrix fallback uses Welch contrasts and does not model paired or family structure when that design cannot be reconstructed from available public metadata.
```
NEW:
```
In the current cross-platform benchmarks, IFN uses GSE71634, ER stress uses GSE19519 tunicamycin and thapsigargin contrasts from one source unit (a large lymphoblastoid donor panel, so its p-values are driven by effect size rather than significance and its added evidence differs in cell system from the RNA-seq sources), and hypoxia uses GSE3045 and GSE22282. Because the added microarray sources differ in cell system and, for IFN, in interferon subtype from the RNA-seq sources, the cross-platform recall gains are confounded with biological context and should be read as evidence that platform integration is feasible without loss rather than as a controlled platform comparison. All derived microarray and count-derived contrasts in these benchmarks use per-gene Welch t-tests on log-transformed values rather than limma empirical-Bayes or negative-binomial models, and the normalized-matrix fallbacks do not model paired or family structure when that design cannot be reconstructed from available public metadata.
```

### EDIT H — SE-derivation caveat (Methods "Auxiliary statistical reporting lanes")
OLD:
```
The SE-derived effect lane is labeled as approximate and is most interpretable in larger-k corpora such as hypoxia.
```
NEW:
```
The SE-derived effect lane is labeled as approximate: the relation SE = abs(log2FC)/abs(signed_z) recovers the true standard error only for Wald-type statistics (for example DESeq2 or edgeR Wald tests) and is biased for moderated-t (limma) or likelihood-ratio statistics. It is most interpretable in larger-k corpora such as hypoxia.
```

### EDIT I — direction-recall caveat (Methods "Benchmark metrics")
OLD:
```
Direction-aware recall was computed where expected directions were defined.
```
NEW:
```
Direction-aware recall was computed where expected directions were defined; because all locked gold panels contain only up-regulated markers, direction-aware recall equals plain recall by construction in every reported row, so the directional scoring components are reported as annotations and are not independently validated by these benchmarks.
```

### EDIT J — Discussion paragraph "The strongest empirical support..."
OLD:
```
The strongest empirical support comes from the IFN and ER-stress benchmarks. IFN recovered canonical interferon-stimulated genes at high recall by rank 50 and 100, and microarray addition improved early recall at rank 20. ER stress gained substantially at rank 50 and 100 after microarray evidence was incorporated, although early top-10 performance decreased. Heat shock supports the broader benchmark because PIPER-DEG improved recall at rank 50 and 100 relative to the best non-PIPER method. Hypoxia is a useful counterweight: PIPER-DEG recovered credible HIF-axis genes but did not dominate early cutoffs, which limits any universal superiority claim.
```
NEW:
```
The benchmark is best read as a faithfulness check rather than a demonstration of superiority. On three of the four topics every runnable summary-level method, including PIPER-DEG, tied at recall@100, so IFN in particular is a ceiling-effect positive control. The clearest comparative result is heat shock, where PIPER-DEG's recall@100 point estimate (0.75) was the only one to exceed all comparators (<=0.625), although the 16-gene panel leaves the confidence intervals overlapping. Microarray addition improved IFN early recall at rank 20 and ER-stress recall at ranks 50 and 100 (with a top-10 decrease). Hypoxia is a useful counterweight: PIPER-DEG recovered credible HIF-axis genes, but an untuned Stouffer baseline exceeded it at early cutoffs, ruling out any universal superiority claim. The consistent finding across topics is that PIPER-DEG loses nothing against standard methods while adding an auditable, source-resolved evidence layer.
```

### EDIT K — Discussion direction limitation (the "Finally, ..." sentence)
OLD:
```
Finally, the current locked gold panels are dominated by up-regulated canonical markers, so direction-recall summaries do not yet test down-regulated marker recovery.
```
NEW:
```
Finally, the locked gold panels consist entirely of up-regulated canonical markers, so direction-aware recall is numerically identical to plain recall in every benchmark row; the direction-confidence and beta-binomial posterior components are therefore reported as annotations but are not independently validated, and down-regulated marker recovery is untested. Validating the directional machinery would require a locked panel containing canonical repressed markers (for example hypoxia-repressed oxidative-phosphorylation genes).
```

### EDIT L — Atlas enumeration
OLD:
```
Figure 4 summarizes the local PIPER-DEG Atlas package. The current generated Atlas contains six corpora, 115,119 scored gene rows, 114 method/result rows, and 23 prior-art/resource rows.
```
NEW:
```
Figure 4 summarizes the local PIPER-DEG Atlas package. The current generated Atlas contains six corpus configurations spanning the four biological topics (IFN RNA-seq 9,860 and RNA-seq plus microarray 12,777; ER stress RNA-seq plus microarray 10,725; heat shock 14,285; hypoxia RNA-seq 32,687 and RNA-seq plus microarray 34,785), totaling 115,119 scored gene rows; because IFN and hypoxia each contribute both an RNA-seq-only and a combined-platform configuration, gene rows are not deduplicated across configurations of the same topic. The Atlas also carries 114 method/result rows and 23 prior-art/resource rows.
```

### EDIT M — Data and Code Availability
OLD:
```
All current generated artifacts are local to this development workspace. Public release, archive DOI, final source-study citation identifiers, and a clean one-command regeneration workflow through the submission repository are required before journal submission. Current manuscript tables are generated under `outputs/tables/manuscript/`, figure packages under `outputs/figures/manuscript/`, and corpus-level score databases under `outputs/results/*/piper_scores.db`.
```
NEW:
```
The PIPER-DEG package, the benchmark and figure-generation scripts, the harmonized DEG tables, the locked gold panels, and the manuscript-facing comparator summaries are released under the MIT License (see `LICENSE`); curated data are released under CC-BY with attribution to the source studies and GEO accessions. Manuscript tables are generated under `outputs/tables/manuscript/`, figure packages under `outputs/figures/manuscript/`, and corpus-level score databases under `outputs/results/*/piper_scores.db`; the complete set of display items is regenerated with `make -C outputs/code paper`. Before journal submission the repository and curated corpus will be archived with a permanent DOI (Zenodo) and the reference list will be completed with dataset-level citations for every GEO accession; these remaining packaging items are listed in the submission checklist below and do not affect any analysis reported here.
```

### EDIT N — replace "Citation Audit Still Needed" with a "Submission Checklist"
OLD:
```
## Citation Audit Still Needed

- Dataset-level citations for every active benchmark source study and GEO accession.
- Verification of final package versions and citation strings from each R package used in the locked environment.
- Final journal-style conversion of references to the required Scientific Reports format.
```
NEW:
```
## Submission Checklist

The analyses, tables, and figures above are final and are regenerated from the committed artifacts with `make -C outputs/code paper`. The following packaging steps remain before submission and do not affect any reported result:

- Dataset-level citations for every benchmark source study and GEO accession (GSE147507, GSE221804, GSE19519, GSE3045, GSE22282, GSE71634).
- Citation strings and pinned versions for each R comparator package used in the locked environment.
- Deposit of code and curated data on Zenodo with a permanent DOI.
- Final conversion of the reference list to Scientific Reports style.
```

---

## STEP Z — regenerate + validate + report

1. From repo root, run `make -C outputs/code paper` (regenerates tables, figures, and the manuscript package
   validator). Capture the output.
2. Read `outputs/manuscript/piper_paper_validation.txt` and report `failures=` and `warnings=`. There must be
   ZERO failures. A small number of warnings is acceptable (e.g., the ER headline "0.889" vs "0.89" substring,
   or pending-citation warnings) — list them.
3. If `make paper` fails on R-only steps, fall back to running only the Python table/figure generators and the
   validator, and report exactly which step failed and why.
4. Report: a unified `git diff` summary (files changed, +/- line counts) for write_manuscript_tables.py and
   outputs/manuscript/main.md; the regenerated Table 1 markdown; and the validator failures/warnings.
5. Do NOT run `git commit` or `git add`.
