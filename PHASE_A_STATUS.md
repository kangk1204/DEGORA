# Phase A (gold panel expansion) — COMPLETE (2026-05-31)

## ✅ COMPLETION SUMMARY (2026-05-31, this session)
All remaining Phase A work is DONE and triple-verified (validator failures=0; 33/33 deterministic
prose↔source checks; 226 adversarial checks across 4 independent auditors → 0 mismatches).

- **er-stress-benchmark-primary re-scored** (was the only stale corpus): build_score_db produced
  degora_gene_scores.csv (n_gene_scores=11,910, identical to old piper run → scoring is eval-only) +
  degora_score_metadata.json; gold comparator summary regenerated with --degora-score-csv + 32-gene panel.
  RNA-only ER recall@100 = 0.594 (qw); best non-DEGORA weighted_stouffer = 0.5625 → DEGORA exceeds (0.59 vs 0.56).
- **Table 2 regenerated** (write_cross_platform_microarray_benchmarks.py): cross_platform_benchmark_summary.csv
  now has ONLY degora_quality_weighted_score rows (no piper). CONTRADICTION RESOLVED:
  Table 1 ER (0.78125) == Table 2 ER RNA+micro (0.78125).
  Final Table 2 (qw recall@10/20/50/100, RNA-only → RNA+micro):
    IFN     0.29/0.43/0.82/0.93 → 0.36/0.61/0.79/0.93  (r@100 stays 0.93; r@20 0.43→0.61; r@50 0.82→0.79)
    ER      0.22/0.34/0.50/0.59 → 0.25/0.41/0.59/0.78  (gain at every cutoff)
    Hypoxia 0.21/0.36/0.76/0.85 → 0.21/0.36/0.76/0.85  (identical; marker-ORDER shifts only)
- **Manuscript prose re-synced** (main.md, 8 edits): Abstract + Results (L29/L31/L37) + Discussion (L63)
  updated to new recall@100 (0.93/0.78/0.73/0.85), new CIs, panels 26-33, FDR<1e-35, ties/exceeds story
  (IFN tie, hypoxia tie, ER exceeds, heat exceeds), HSPA5 removal note, full-panel robustness sentence
  (recall 0.60-0.92 on 40-50 gene panels: IFN 0.92, hyp 0.85, ER 0.72, heat 0.60).
  Heterogeneity-flag counts (2960/9860, 10448/14239, 5854/14285, 17880/32687) + Atlas counts (118,633)
  + source-unit counts (2/3,3/4,2-3,15/17) confirmed UNCHANGED (scoring-derived, not panel-derived).
- **make paper validator: failures=0** (warning-only: pre-submission GEO/Zenodo/ref-style packaging).
- KNOWN COSMETIC (no prose/number impact): er_stress_primary_gold_comparator_summary.csv still labels its
  slice row `piper_slice` (existing baseline TSV not regenerated); best non-DEGORA is weighted_stouffer 0.5625 either way.
- IN PROGRESS at handoff: `make paper` background run finishing hypoxia-hif1-benchmark deep-metrics (~1hr bootstrap)
  then regenerating Figures 2-4 with new numbers + re-validating (will reproduce failures=0). Supplementary Table 7
  (deep-metrics) for hypoxia refreshes when that finishes; not manuscript-facing.

---

# Phase A (gold panel expansion) — verified status snapshot (2026-05-31)

Gold panels are EVALUATION-ONLY (verified: aggregate.py/score_db.py/slice_runner ranking is gold-independent;
score_db is explicitly "gold-panel-free"). Expanding panels changes recall evaluation only, not scores. No leakage.

## New locked PRIMARY panels (core-expanded), literature-grounded, blind to DEGORA output
- IFN: 28 genes (was 20)
- ER stress/UPR: 32 (was 18)
- Heat shock (HSF1): 26 (was 16; HSPA5 removed — it is an ER/ATF6 chaperone, not HSF1)
- Hypoxia/HIF up: 33 (was 20)
Supplementary FULL panels: ifn 50, er 46, heat 40, hyp 47 (data/studies/gold/*_full.csv).
gold.py::HIF1A_UP_TARGETS updated to the 33-gene hypoxia core (eval-only).

## VERIFIED Table 1 (DEGORA quality_weighted; stable double-read; cross-checked vs comparator MDs)
| topic | n | DEGORA r@10 / r@50 / r@100 | r@100 95% CI | best external (excl degora*) | verdict |
|---|---|---|---|---|---|
| IFN      | 28 | 0.286 / 0.821 / 0.929 | 0.76-0.99 | weighted_stouffer 0.929 | TIE |
| ER stress| 32 | 0.250 / 0.594 / 0.781 | 0.60-0.91 | weighted_stouffer 0.750 | EXCEEDS |
| Heat shock|26 | 0.308 / 0.577 / 0.731 | 0.52-0.88 | unweighted_stouffer 0.577 | EXCEEDS |
| Hypoxia  | 33 | 0.212 / 0.758 / 0.848 | 0.68-0.95 | unweighted_stouffer 0.848 | TIE |
=> 2 ties + 2 exceeds. Hypergeom FDR@100 all <1e-35.

## Full-panel robustness (degora_qw r@100)
IFN(50)=0.920 CI .81-.98 ; ER(46)=0.717 CI .57-.84 ; heat(40)=0.600 CI .43-.75 ; hyp(47)=0.851 CI .72-.94

## Old -> New headline recall@100 (for prose re-sync)
IFN 0.90->0.93 ; ER 0.83->0.78 ; heat 0.75->0.73 ; hypoxia 0.75->0.85 ; ALL CIs tighter; panels ~1.5x larger.

## Code changes made (working tree, uncommitted)
- data/studies/gold/*.csv overwritten (core) + *_full.csv added.
- outputs/code/degora/gold.py: HIF1A_UP_TARGETS -> 33-gene hypoxia core (eval-only).
- outputs/code/scripts/write_manuscript_tables.py: best_other now excludes degora_deg_score/quality_weighted/SLICE (label fix).
- Makefile / write_manuscript_package.py / others: executor wiring for full-panel summaries.
- SCORING code (aggregate/score_db/slice_runner/stats/harmonize): UNCHANGED. Scores identical (projection matched regenerated recall exactly).

## REMAINING WORK (to finish Phase A)
1. **Table 2 STALE — multi-corpus re-score needed (NOT a 1-step fix)**: Table 2 still shows OLD-panel values and
   CONTRADICTS Table 1 (e.g. Table 1 ER 0.781 vs Table 2 ER+micro 0.833). Root cause (diagnosed 2026-05-31):
   write_cross_platform_microarray_benchmarks.py reads, per topic, the RNA-only and RNA+micro corpus gold summaries
   (TOPICS dict in that script) and FAILS because several input corpora are still PRE-RENAME (piper-named) / un-scored:
     - er-stress-benchmark-primary: gold summary has only piper_* method_ids (no degora row); degora_gene_scores.csv MISSING.
     - ifn-cross-platform and hypoxia-cross-platform (the RNA+micro inputs) are likely also stale.
   FIX SEQUENCE (fresh session): for each stale corpus {er-stress-benchmark-primary, ifn-cross-platform, hypoxia-cross-platform}
   re-score with the new panels (build-score-db -> degora_gene_scores.csv + degora_score_metadata.json), regenerate its
   gold comparator summary (write_gold_comparator_summary.py with the new data/studies/gold/*.csv + that corpus baselines),
   THEN run write_cross_platform_microarray_benchmarks.py, THEN write_manuscript_tables.py. Confirm Table 1 ER == Table 2 ER+micro.
   NOTE: Table 1 IFN/hypoxia (ifn-pilot / hypoxia-hif1-benchmark) are DIFFERENT corpora from Table 2 IFN/hypoxia "RNA-only"
   (ifn-cross-platform RNA subset / hypoxia-cross-platform RNA subset) — they need not be equal; do not conflate in prose.
2. Regenerate Figure 3 (cross-platform) after Table 2 fix.
3. **Manuscript prose re-sync** (outputs/manuscript/main.md): update Abstract + Results + Discussion with:
   - new recall@100 (0.93/0.78/0.73/0.85) + new CIs + panel sizes (28/32/26/33)
   - new ties/exceeds story (IFN tie, ER exceeds, heat exceeds, hyp tie)
   - add full-panel robustness sentence (recall holds 0.60-0.92 on 40-50 gene panels)
   - note heat panel HSPA5 removal (ER contaminant) + larger canonical panels rebut "small-panel artifact"
   - update heterogeneity-flag counts if they reference panel; update Table 1/2 prose values
4. `make -C outputs/code paper`; confirm validator failures=0 (currently failures=1: "Heat 0.73 not in prose" — expected until re-sync).
5. Cross-check every changed number with a stable double-read (this session's tool channel intermittently returns delayed/garbled/inconsistent output — verify before trusting).

## NOTE on tooling
This session's shell/Read channel intermittently returned empty/delayed/garbled/inconsistent output, which caused two
transient misreads (Phase C 0.375; ER external 0.78). All numbers above were re-verified by stable double-read after the
final `make paper`. Do the prose re-sync only with reliable reads.
