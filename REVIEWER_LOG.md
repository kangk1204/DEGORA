# Adversarial Reviewer Log

(Append-only. One entry per adversarial review pass. The loop populates this in Phase 4 of each iteration.)

## Iteration 1 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on technical soundness and reproducibility.

Major concerns:

1. The vertical slice is useful engineering evidence, but it is not yet a valid S0 study-catalog result. Two of the five entries are public tutorial-derived DESeq2 all-gene outputs for GSE106305, not as-published supplementary DEG tables from the original paper. They are acceptable as a pipeline smoke test only if clearly labeled as such and replaced before full S0.
2. The primary V1 recall target was not met. Recall@50 was 8/20 = 0.40 and recall@100 was 12/20 = 0.60. This is a mixed/weak biological recovery result, not evidence that PIPER works.
3. Ranking noise is likely inflated by non-protein-coding and pseudogene rows in the all-gene tutorial outputs. The top of the consensus contains many noncanonical symbols before several locked HIF targets. A protein-coding-aware filter should be pre-registered and reported as a separate iteration, not silently substituted for iteration 1.
4. GSE76743 mRNA lacks nominal p-values in the downloaded processed table; using adjusted p-values as a conservative p surrogate is defensible for a smoke test but not for final methods claims.
5. The lncRNA-only table from GSE76743 exercises parser heterogeneity but does not satisfy the final hypoxia mRNA catalog intent. Do not count it toward the final S0 pass criterion.

Minor concerns:

- The R environment currently lacks DESeq2 and rpy2, so R-side baselines and DESeqResults loading are not ready. This is not an escape hatch yet because iteration 1 is Python-only, but S1/S4 will require environment work.
- The RRA baseline is explicitly an approximation. This is acceptable for the slice but must be replaced by RobustRankAggreg or a validated implementation before S1 claims.
- The initial slice runtime was longer than expected for 100k rows because aggregation is still Python-loop based. It completed, but vectorization should be considered before scale-up.

Reviewer verdict: continue, but do not scale to full S0/S1 until the slice is debugged. Next iteration should pre-register a minimal ranking-noise fix, retain the iteration-1 metrics as the baseline, and rerun without hiding the failed top-50 result.

## Iteration 2 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on technical soundness and reproducibility.

Major concerns:

1. The pre-registered protein-coding filter did not improve the primary metric. Recall@50 decreased from 0.40 to 0.35, while recall@100 remained 0.60. This falsifies the simple "pseudogene/noncoding noise is the main top-50 problem" hypothesis.
2. The filter removed some useful evidence for VEGFA from source tables where `Gene.type` was available, moving VEGFA from rank 27 to rank 62. This shows that source-specific filtering can change target support and must be handled carefully.
3. The current five-input slice is now clearly an engineering smoke test, not a biological validation set. The most likely limiting factor is input eligibility and heterogeneity: tutorial-derived all-gene outputs, a lncRNA-only parser-test table, and one table using adjusted p as the p-value surrogate.
4. Further tuning on these same five inputs risks overfitting the smoke test to the HIF target list. The next scientifically defensible step is curation replacement, not more scoring tweaks.

Minor concerns:

- The new filter is implemented cleanly and recorded in metrics, so the negative result is reproducible.
- The catalog now has useful optional filter metadata for future studies, but it should remain source-declared rather than target-driven.

Reviewer verdict: do not continue tuning the same slice. Proceed to improved as-published DEG table curation for iteration 3, and keep both iteration-1 and iteration-2 metrics as failed/mixed exploratory evidence.

## Iteration 3 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on technical soundness, data provenance, and overfitting risk.

Major concerns:

1. The curation improved, but the primary V1 metric worsened. Recall@50 is 5/20 = 0.25 and recall@100 is 9/20 = 0.45. This is a null result under the locked pass criterion, even though several top genes are plausible hypoxia markers.
2. Removing the tutorial-derived prostate cancer contrasts made the active catalog more defensible but also removed strong HIF-target support. This shows why the earlier smoke-test recall was not reliable biological validation.
3. Bauer AH and CH are two contrasts from the same paper and cell system. Counting them as independent active rows may overweight one experimental source; iteration 4 needs a source/paper-level sensitivity diagnostic before scale-up.
4. Bauer tables lack nominal p-values, so padj is used as a conservative p surrogate. This is documented, but it changes the effective gene universe and may disadvantage less significant but directionally consistent HIF targets.
5. HYP003 remains only provisionally acceptable because the immediate source is a Bioconductor example object rather than a directly curated article supplement. It should be replaced or provenance-hardened before final S0.
6. The top-ranked genes include SLC2A1, EGLN3, MIR210HG, NDRG1, ADM, P4HA1, SLC2A3, BNIP3L, and VEGFA. The locked target list therefore may be a narrow benchmark rather than a complete biological truth set, but it must not be expanded post hoc.

Minor concerns:

- The new catalog eligibility flag is a clean and auditable way to preserve failed/smoke-test rows without silently dropping them.
- The Mendeley hDASMC sign orientation was checked against canonical HIF genes and sign-flipped in a derived file; this decision must remain traceable.
- The curation ledger is useful, but future iterations should add checksums for downloaded public files before the catalog is scaled.

Reviewer verdict: continue, but only with diagnostics. Do not tune the scoring formula or add the classifier yet. The next iteration should quantify per-study target support, source-level sensitivity, p-value surrogate sensitivity, and heterogeneity so the loop can decide whether to scale S0 or repair a methodological defect.

## Iteration 4 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the null curated slice is a bug, a benchmark limitation, or real cross-source heterogeneity.

Major concerns:

1. The diagnostics do not find a global sign-orientation error. All five active sources have positive-direction fractions above 0.77 for present locked HIF targets, so rerunning the same catalog with sign changes would be unjustified.
2. No sensitivity variant rescues the locked top-50 criterion. The best top-50 diagnostic is excluding HYP003 at 7/20 = 0.35, still below the original smoke-test top-50 and well below any validation threshold.
3. HYP003 remains suspect: dropping it improves top-50 recovery, but the effect is not large enough to justify silent exclusion. It should remain provisional until replaced by better direct supplementary curation.
4. Collapsing Bauer acute/chronic contrasts by paper ID leaves recall@50 unchanged and lowers recall@100 to 0.40. Bauer non-independence is therefore not the main explanation for the null result, but source-level sensitivity should remain part of later validation.
5. The locked target list may be too narrow to reflect all credible hypoxia biology, since top-ranked non-gold genes include strong known hypoxia-response markers. This is not permission to edit the gold list; it argues for additional validation endpoints.
6. The first diagnostic run exposed a scalability problem in the main Python-loop consensus implementation. The vectorized diagnostic path is promising, but it must be equivalence-tested before becoming the main implementation.

Minor concerns:

- Positive-consensus-only filtering improves recall@100 to 0.50 but leaves recall@50 unchanged. This is diagnostic information only, not a pretext for post-hoc score changes.
- Excluding padj-surrogate sources leaves only two active studies and recall@50 0.30, so p-value surrogate use contributes uncertainty but does not explain the full null.
- The new diagnostic artifacts have source files and are reproducible from local iteration-3 outputs.

Reviewer verdict: continue. Before scaling S0 to more studies, make the consensus computation fast enough for repeated catalog runs and prove equivalence to the current method. Then add more eligible studies rather than tuning against the HIF target list.

## Iteration 5 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether a performance refactor changed the method.

Major concerns:

1. The vectorized aggregation is a core method-path change, so equivalence evidence is mandatory. The regression test and iteration-3 versus iteration-5 comparison are adequate for the current thin-slice scope: metrics and top-100 genes match, with only floating-point-scale numeric differences.
2. This iteration must not be interpreted as biological progress. Recall remains 5/20 at top 50 and 9/20 at top 100, exactly matching iteration 3.
3. The implementation now has a stable path for repeated S0 scaling, but the current catalog is still only five active rows. Publication-readiness remains far away from S0's target of at least 20 hypoxia studies.
4. Future curation should add checksums and avoid committing ambiguous candidate data without ledger reasons. Provenance quality now matters more than speed.

Minor concerns:

- The old diagnostic-only vectorized code was removed in favor of the shared aggregate path, reducing divergence risk.
- Returning stable empty-column schemas from aggregation helpers is a slight API improvement, but downstream code should still be tested when empty catalogs or strict `min_studies` filters occur.

Reviewer verdict: continue to S0 scaling. The performance repair is acceptable because equivalence is demonstrated and the null biological result is preserved.

## Iteration 6 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on curation integrity, source independence, and overfitting risk.

Major concerns:

1. The locked V1 metric improves materially after adding GSE132624 melanoma hypoxia contrasts: recall@50 rises to 10/20 = 0.50 and recall@100 to 12/20 = 0.60. This is encouraging but still not a final validation result.
2. The new signal is not independent enough for a strong claim. All three new active rows come from one GEO series and one paper family; the collapse-by-paper diagnostic drops recall@50 to 0.35 and recall@100 to 0.45.
3. Adding the deferred 501mel 12 h table now would amplify the same source-family problem. It should stay deferred until the method has explicit source-family or repeated-timepoint handling.
4. HYP003 remains provisional and now looks increasingly questionable: excluding it improves recall@50 to 0.55 and recall@100 to 0.70. This does not justify silent exclusion, but replacing it with a direct primary table should remain high priority.
5. The padj-surrogate sensitivity also improves recall, suggesting that current conservative p-value substitution may suppress some target recovery. This should be reported as uncertainty, not used as a post-hoc excuse to alter p handling.
6. Several excluded/deferred candidates are not dead ends for the overall project: GSE108676 may become usable after audited transcript-to-gene mapping, and GSE316921 may support V2 raw-count reanalysis. They are correctly excluded from the current S0 as-deposited DEG catalog.

Minor concerns:

- The iteration-6 ledger is stronger than earlier curation because it includes checksums and categorical status before rerun.
- The GSE132624 target orientation fractions are high, reducing concern about sign inversion.
- The top genes include canonical hypoxia-response markers outside the locked HIF list, which supports broader validation endpoints later but does not permit changing the locked V1 list.

Reviewer verdict: continue S0 scaling, but prioritize independent paper/dataset additions over more rows from the same GEO series. Keep the paper-collapse diagnostic as a standing sensitivity check and do not tune scoring while the catalog is still below the target number of independent hypoxia studies.

## Iteration 7 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether independent curation reduces the iteration-6 source-family artifact.

Major concerns:

1. The direction is positive, but the active catalog is still too small: 10 active rows and fewer independent paper families than the S0 target. This remains a development slice, not publication-ready evidence.
2. Paper-collapse recall improved from 0.35/0.45 to 0.40/0.50, but it is still notably below active-all recall of 0.55/0.65. The method is still sensitive to repeated contrasts from the same paper/source family.
3. HYP013 depends on an NCBI ESummary-derived RefSeq-to-symbol mapping. The mapping is auditable and complete for this run, but it is a derived artifact and should remain versioned with its source file; future reruns should not silently remap against a changed NCBI state without noting it.
4. HYP003 remains provisional. Excluding it still improves recall@100 to 0.70, so replacing or hardening that source should remain a priority before final S0.
5. Excluding padj-surrogate sources improves recall@50 to 0.60 and recall@100 to 0.70. This continues to show uncertainty from conservative p-value substitution, not permission to tune p-value handling post hoc.
6. The candidate search found many expression-only and raw-count-only tables. These should not be forced into S0; if the project needs them, they belong in V2 with a separate preregistered raw-count reanalysis protocol.

Minor concerns:

- GSE70544 and GSE108676 both show strong positive orientation for locked HIF targets, reducing concern about contrast inversion.
- The GSE328627 oxygen contrast is tempting but correctly excluded because 1,680 rows is not near-full-universe evidence.
- The iteration-7 ledger is compact and traceable, but future source searches should preserve the broader FTP scan as a reusable artifact if many candidates are screened.

Reviewer verdict: continue S0 scaling with independent primary DEG sources. The improvement is encouraging and less source-family-dependent than iteration 6, but still not robust enough for final validation or method claims. No scoring changes are justified.

## Iteration 8 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether broad curation is improving independent-source robustness.

Major concerns:

1. Active-all recall did not improve after adding HYP014: it remains 0.55 at top 50 and 0.65 at top 100. This is not a new biological performance gain.
2. Paper-collapse recall improved to 0.45/0.55, which is useful, but still materially below active-all. The source-family gap is narrowing, not solved.
3. The broad GEO scan shows many candidate series but few eligible as-deposited full DEG tables. The project should avoid quietly drifting into raw-count reanalysis under the S0 label.
4. HYP003 remains a provenance weakness and has a measurable effect: excluding it raises recall@50 to 0.60 and recall@100 to 0.70. It should be replaced or hardened before final claims.
5. Excluding padj-surrogate sources now raises recall@100 to 0.75. This uncertainty is growing as a reporting issue and should be handled transparently in later methods/results text.
6. GSE95280 is eligible-looking but correctly deferred. Adding more melanoma timecourse rows would improve row count while worsening independence.

Minor concerns:

- HYP014 has excellent locked-target orientation and a full tested sheet, so the inclusion itself is defensible.
- The 167-row GEO scan artifact is valuable and should reduce repeated exploratory searching in later iterations.
- The iAT2 cell system broadens biology beyond vascular/cardiomyocyte/melanoma sources.

Reviewer verdict: continue S0 scaling. The iteration is acceptable because it improves paper-collapsed robustness and preserves negative active-all evidence, but the catalog is still too small and too source-family-sensitive for final validation.

## Iteration 9 — adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the positive iteration-9 jump is technically sound and honestly scoped.

Major concerns:

1. The improvement is substantial: active-all recall rises to 0.70/0.70 and collapse-by-paper recall also rises to 0.70/0.70. This is the first iteration where the paper-family diagnostic no longer lags the active-all diagnostic, but the catalog still has only 13 active rows and 10 paper/source-family units, below the final S0 target.
2. HYP003 is not hardened. The primary GSE89831 audit found raw/sample files but no direct full DEG table, and the TFEA.ChIP example source has pipeline/provenance ambiguity. Excluding HYP003 no longer changes recall, which lowers its immediate influence, but final S0 should still replace it or mark it as provisional in sensitivity tables.
3. HYP016 is useful but not fully method-annotated. The table has LFC and nominal p-value for a direct hypoxia-vs-normoxia contrast, but the pipeline is not resolved and adjusted p-values are absent. This row is acceptable for the current curation slice only if the limitations remain visible.
4. HYP015 has 11 p-value underflow/clipping rows. The clipped rows were reviewed and the registered pipeline still passes, but later methods text must describe how zero p-values from source tables are handled.
5. The positive result should not trigger post-hoc optimism. No scoring changes were made, but the next iteration should test whether the gain persists as independent sources are added, rather than pivoting immediately to manuscript claims.
6. V2 and V7 readiness inventories are good planning artifacts, but they are not validation evidence. No V2 pipeline-perturbation benchmark or V7 classifier accuracy claim exists yet.

Minor concerns:

- HYP015 and HYP016 both show good locked-target orientation, reducing concern about contrast inversion.
- The source-sibling discipline improved because harmonized iteration-9 CSV/parquet outputs now have `.source` files.
- The candidate ledger correctly excludes raw-count/expression-only files from S0 while preserving them for V2 inventory.

Reviewer verdict: continue S0 scaling, with HYP003 replacement/provisional handling and HYP016 pipeline provenance as priority risks. Iteration 9 is a strong curation advance, not publication-ready validation.

## Iteration 10 worker-4 TNF gate audit — 2026-05-25

Reviewer stance: Scientific Reports methods reviewer focused on the cross-domain TNF-alpha lane, source-family independence, and protocol-boundary discipline.

1. The TNF-alpha pass gate is not met in this worker audit: 0 active TNF-alpha rows and 0 independent paper/source-family units were verified. The scarcity branch is therefore the only defensible terminal evidence from this slice, supported by 20 high-signal screened candidates in `data/studies/curation/iter10_tnfa_candidate_ledger.csv`.
2. The hard topic boundary is preserved. IFN/TNF, IL-1β/TNF, TNF/TGF-β, LPS, pathogen, and generic inflammatory-stimulation candidates are not counted as active TNF-only evidence. The project must not auto-switch to IFN or LPS to rescue the TNF gate without an explicit constitution amendment.
3. The new TNF gate checker is appropriately narrow: it counts active rows and independent source-family units from TNF-specific ledger/catalog inputs and accepts the scarcity path only when at least 20 candidates exist with a scarcity report. It does not change PIPER scoring or the locked HIF target list.

Reviewer verdict: the worker-4 evidence is a rigorous negative/blocked TNF-corpus gate result rather than a cross-domain success claim. Continue TNF curation only by resolving full DEG tables for the logged candidates; do not broaden the biology domain silently.

## Iteration 11 final integrated verification - adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the integrated publication upgrade now supports manuscript claims.

Major concerns:

1. The final claim state remains **not ready**. Hypoxia is stable at recall@50/@100 of 0.70/0.70 after paper collapse, but 14 active rows do not meet the PRD's 15-row narrowed-claims fallback or the preferred 20-row S0 target.
2. TNF-alpha is a rigorous scarcity result, not a validation success: 40 candidates were screened, with 0 active rows and 0 independent units. This is acceptable only if reported as a negative limitation.
3. The fallback-domain firewall held. The IFN/LPS/other-domain work is feasibility-only, no amendment exists, no candidate-domain catalog or gold file was activated, and no candidate-domain scoring or baseline output was generated. This reduces domain-shopping risk but does not create validation evidence.
4. Direct-prior-art baseline parity is still blocked. RobustRankAggreg, MetaVolcanoR, hStouffer, and AWmeta have blocked/unresolved rows in the parity matrix, so any claim that PIPER outperforms direct prior art must be removed or deferred.
5. Provenance is materially stronger: iteration-11 result artifacts, including baseline TSV tables, have non-empty `.source` sidecars and `make check` passes. This supports reproducibility but does not overcome the validation and comparator gates.

Minor concerns:

- HYP003/HYP016 sensitivity/provenance concerns remain visible through diagnostics and should stay in limitations until resolved or replaced.
- Classical baseline outputs are useful sanity comparators, but they are not substitutes for direct prior-art tools.
- The alternative-domain feasibility ranking is useful planning material only; any future domain must be locked by amendment before scoring.

Reviewer verdict: do not submit as a strong or narrow methods paper yet. The next defensible milestone is either one more eligible active hypoxia S0 row plus conservative, non-superiority claims, or a fuller S0 target plus resolved direct-prior-art baseline parity. Keep TNF scarcity and baseline blockers explicit in the manuscript trail.

## Iteration 12 Ralph follow-up - adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the new S0 source and baseline parity work genuinely improve manuscript viability.

Major concerns:

1. The claim state improves to **narrow but defensible paper**, but only for conservative non-superiority claims. The active hypoxia count now reaches 15 and paper/source-family units reach 12, with paper-collapse recall@50/@100 preserved at 0.70/0.70.
2. HYP018 is acceptable as an independent near-full GEO processed source, but it is not strong biological support by itself. Its locked HIF-target positive orientation fraction is 0.50, and the source pipeline is recorded as GEO processed contrast unspecified. This must be framed as heterogeneity/robustness stress evidence, not as a clean HIF activation replicate.
3. Direct-prior-art parity is improved but not resolved. RobustRankAggreg now runs and emits a uniform-schema result, and MetaVolcanoR/metafor package availability is verified. However, MetaVolcanoR, hStouffer, and AWmeta/metafor still lack faithful wrappers/adapters. Superiority claims over prior meta-analysis tools remain disallowed.
4. The publication gate appropriately inherits TNF scarcity and alternative-domain firewall evidence, but these are limitations and integrity controls, not positive validation domains.
5. HYP003 remains provisional. The exclude-HYP003 sensitivity is stable at 0.70/0.70, which lowers but does not eliminate the need to mark it carefully in manuscript-facing tables.

Minor concerns:

- The candidate ledger correctly excludes small-row, non-human, replicate-insufficient, no-direct-hypoxia, and counts/expression-only sources rather than inflating S0 with raw reanalysis.
- Provenance is stronger because iter-12 raw, harmonized, baseline, diagnostic, rank-plane, publication-gate, and audit artifacts have `.source` sidecars and `make check` now audits iter-12.
- Rank-plane remains useful as a supplementary heterogeneity diagnostic, not a scoring/tuning mechanism.

Reviewer verdict: proceed toward a narrow, conservative manuscript frame if needed, but do not write a SOTA comparator paper yet. The strongest current contribution is a reproducible, heterogeneity-tolerant DEG evidence-integration workflow with locked recall stability under source-family collapse. The next acceptance-probability gain comes from faithful MetaVolcanoR/hStouffer/AWmeta wrappers or additional independent S0 sources closer to the preferred 20-study target.

## Iteration 13 Team follow-up - adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the added S0 source and comparator parity now support stronger PIPER claims.

Major concerns:

1. The claim state remains **narrow but defensible paper**, not SOTA superiority. PIPER recall@50/@100 is 0.70/0.75, but weighted and unweighted Stouffer match 0.70/0.75, MetaVolcanoR reaches 0.60/0.75, and RobustRankAggreg reaches 0.70/0.70.
2. HYP019 is a defensible S0 addition: it is a full GEO-deposited DESeq2 HK-2 hypoxia-vs-normoxia table, has 3 samples per arm, and has strong locked-target orientation. The decision to defer the paired 786-0 table from the same GSE avoids source-family inflation and should be kept.
3. hStouffer and AWmeta are properly blocked rather than approximated. This is methodologically conservative, but it means direct-prior-art parity is still incomplete under the equal-tuning gate.
4. The apparent PIPER advantage is not higher raw recall versus all baselines. The stronger story is heterogeneity-tolerant, auditable merging of mixed as-deposited DEG tables with stable locked hypoxia recovery under source-family collapse.
5. HYP003/HYP016 caveats remain manuscript limitations. Exclude-HYP003 stability helps, but provenance caveats should remain visible in tables.

Minor concerns:

- The comparator summary is useful and should be cited from the generated CSV/MD rather than recomputed informally.
- Sign vote performs poorly, which is expected for a coarse direction-only method and should not be overused as a straw-man baseline.
- The provenance expansion to iter-13 and 30 passing tests materially improves reproducibility.

Reviewer verdict: the project is now stronger and more honest. It can support a conservative methods/resource manuscript focused on reproducible heterogeneous DEG evidence integration, but not a claim that PIPER outperforms existing meta-analysis tools. The highest-value next work is either more independent S0 sources or a preregistered path for hStouffer/AWmeta using original variance/SE inputs.

## Iteration 14 SOTA-oriented follow-up - adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether failed S0 expansion and faithful-input comparator blockers were handled honestly.

Major concerns:

1. The active S0 evidence base did not grow. Iteration 14 screened 6 additional hypoxia candidates but activated 0 rows, so the project remains at 16 active hypoxia studies and 13 paper/source-family units from iteration 13.
2. This is not a negative mark if framed correctly. The screened candidates were excluded for defensible categorical reasons: gene-list-only, raw-count-only, counts/expression-matrix-only, below full/near-full gene universe, or non-direct drug/HIF-PHI contrasts. The exclusion trail is more valuable than inflating S0 with non-comparable inputs.
3. hStouffer remains blocked for a substantive method reason. Only HYP015 has a compatible original DESeq2-like source header; 15 of 16 active studies lack the required original DESeq2 fields or use non-DESeq2/unspecified pipelines. Running hStouffer by inventing `lfcSE` from p-values or signed-z would be a non-faithful comparator.
4. AWmeta remains blocked for a parallel substantive reason. Only HYP015 has an explicit variance/SE-like source column; 15 of 16 studies lack auditable original variance/SE or documented equivalent weights. Any AWmeta/AW-REM run that silently filters to HYP015 or imputes SE would be methodologically misleading.
5. The iteration strengthens integrity more than performance. It does not justify stronger SOTA claims, but it makes the conservative PIPER claim cleaner: PIPER is useful when heterogeneous DEG tables lack a uniform original variance/SE schema that effect-size meta-analysis tools require.

Minor concerns:

- The iteration-14 rerun-readiness gate correctly blocks downstream reruns while no new `iter-14_harmonized.csv` exists.
- The code review found and removed duplicated AWmeta helper definitions introduced by team merge reconciliation; `make check` now passes with 35 tests.
- GSE167956 raw/count assets should remain clearly labeled as excluded S0 evidence unless a separate preregistered raw-count V2 reanalysis is launched.

Reviewer verdict: keep the manuscript frame narrow and honest. Iteration 14 improves the defense against unfair comparator criticism, but it does not improve empirical S0 coverage or support direct superiority over prior meta-analysis methods.

## Iteration 15 strict S0 expansion - adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the new S0 rows create a real reason to use PIPER instead of existing meta-analysis methods.

Major concerns:

1. The claim state remains **narrow but defensible paper**, not SOTA superiority. PIPER recall@50/@100 is 0.70/0.75, but weighted and unweighted Stouffer match 0.70/0.75, MetaVolcanoR reaches 0.65/0.75, and RobustRankAggreg reaches 0.70/0.70.
2. The new rows materially strengthen the evidence base. HYP020, HYP021, and HYP022 are independent strict S0 additions, raising the active hypoxia set to 19 rows and 16 independent paper/source-family units. Their locked-target orientation fractions are strong enough to reduce concern that the added evidence is biologically inverted.
3. The reason to use PIPER is now clearer but must be phrased precisely. PIPER is defensible as an auditable heterogeneous-DEG integration layer when deposited studies lack uniform original variance/SE/DESeq2-statistic fields. It is not yet defensible as a method that outperforms all existing meta-analysis tools.
4. hStouffer and AWmeta blockers support the use-case rationale only if reported as faithful-input blockers, not as wins. hStouffer is blocked for 18 of 19 active source studies; AWmeta/AW-REM is blocked because uniform original variance/SE or documented equivalent weights are unavailable.
5. The preferred 20 active-row target is still one short. The manuscript can proceed conservatively, but a reviewer may still ask why the search stopped at 19 active rows unless a stopping rule or one more strict source is documented.

Minor concerns:

- The comparator summary title and hStouffer blocker text are now less stale and more manuscript-safe.
- The new PIPER use-case rationale is useful, but it should be treated as a claim-control artifact, not as new validation evidence.
- HYP003 remains provisional and should stay visible in manuscript-facing sensitivity tables.

Reviewer verdict: publication likelihood improved. The defensible paper is a conservative methods/resource manuscript about faithful heterogeneous public DEG integration, with explicit comparator blockers and no superiority claim. The highest-value next step is one more strict independent S0 source or a written stopping rule plus manuscript drafting under the iter-15 rationale.

## Iteration 16 strict S0 target closure - adversarial methods review

Reviewer stance: Scientific Reports methods reviewer focused on whether the new active rows and manuscript skeleton support stronger claims without overreach.

Major concerns:

1. The preferred strict S0 active-row target is now met: 21 active rows and 17 independent paper/source-family units. This removes the iteration-15 reviewer concern that the manuscript stopped one active row short of the target.
2. HYP023 and HYP024 are defensible strict additions because they are direct hypoxia-vs-normoxia RNA-seq source-data tables with adequate replicate counts, large gene universes, and strong locked-target orientation. They must still be treated as one shared paper/source family in collapse claims.
3. The locked recall improvement is real but not uniquely PIPER-specific. PIPER, weighted Stouffer, and unweighted Stouffer all reach 0.70/0.85. This strengthens the biological sanity check but still disallows superiority language.
4. Direct-prior-art blockers remain unresolved. hStouffer is blocked for 20 of 21 active sources, and AWmeta/AW-REM still lacks uniform variance/SE or documented equivalent weights. These are faithful-input blockers, not comparative wins.
5. The manuscript skeleton is appropriately conservative and avoids fabricated references, but it is not submission-ready. Every `[REFERENCE NEEDED]` placeholder must be resolved by citation audit, and every manuscript-facing table needs provenance.

Minor concerns:

- The two Kindrick source-data rows use adjusted p-values as conservative p-value surrogates because nominal p-values are absent. The exclude-padj-surrogate sensitivity is stable, but the limitation should stay visible.
- HYP003 remains provisional and must remain flagged.
- The current manuscript should not imply the full nine-stage methodology is complete; it should distinguish the validated hypoxia vertical slice from future V2/V4/V6/V7 work.

Reviewer verdict: proceed to manuscript hardening. Iteration 16 materially improves publication readiness by closing the S0 active-row target and producing a conservative draft frame, but the paper remains a narrow methods/resource manuscript. The next critical risk is unsupported citation/framing, not additional hypoxia curation.

## Iteration 17 PIPER score DB and local browser - adversarial methods review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new score and web/API layer clarify utility without creating a hidden validation change.

Major concerns:

1. The PIPER score is useful but must be explicitly framed as a prioritization/indexing score, not a calibrated posterior probability, causal score, or new validation endpoint. The metadata warning does this, and the manuscript should preserve that language.
2. The score weights are heuristic. They are defensible for browsing because they expose source support, directional concordance, Stouffer evidence, rank strength, and effect magnitude separately, but they are not tuned or statistically learned. Do not overstate them as an optimized model.
3. The correct unit of repeated evidence is independent source unit, not raw DEG table count. The implementation uses `paper_id`/source-unit collapse for support scoring, which addresses the main inflation risk.
4. The score layer caps non-finite LFC values for finite display/math. This is appropriate for the browser layer and consistent with the methodology's Cuffdiff cap principle, but the log must state that locked validation artifacts were not silently rewritten.
5. A browser/API improves field utility and reproducibility, but it does not strengthen the direct-prior-art superiority claim. Weighted/unweighted Stouffer parity and hStouffer/AWmeta faithful-input blockers remain as before.

Minor concerns:

- The local server is dependency-free and parameterized SQL is used for user-facing queries, reducing deployment and injection risk.
- The SQLite DB is large but reproducible from `make -C outputs/code score-db`; keeping source sidecars for the `.db` artifact is important.
- The UI makes per-gene evidence inspectable, which supports PIPER's provenance claim more directly than a static CSV alone.

Reviewer verdict: accept the score/API/browser iteration as a resource-utility improvement. It supports a stronger "usable evidence atlas" angle for PIPER, but only if manuscript language keeps `piper_score_v1` separate from the locked validation metrics and reports the heuristic nature of the weights.

## Iteration 18 intuitive score presentation - adversarial methods review

Reviewer stance: methods/resource reviewer focused on whether rank, percentile, and tier labels improve interpretability without implying false calibration.

Major concerns:

1. Adding `top_percent`, `percentile`, and tier labels is useful for users, but the labels must not be interpreted as probability of true association. The score metadata and UI warning preserve this distinction.
2. Evidence tiers are heuristic and rank-dependent. They are acceptable as browsing aids because the rules are explicit and stored in metadata, but manuscript text should call them "evidence tiers" or "prioritization tiers," not statistical confidence classes.
3. The labels improve defense against score misuse. `#6 / 34,270`, `top 0.018%`, `16 / 17 source units`, and `100.0% up-concordant` are more interpretable than `94.365929` alone.
4. The tier rules still depend on the current corpus and should not be compared across unrelated corpora without recalculating context-specific total genes and source units.

Minor concerns:

- Keeping iteration-17 artifacts intact and generating iteration-18 artifacts separately preserves provenance.
- API and UI tests now check the new label fields, reducing drift risk.
- The browser still exposes component scores, which is important because aggregate tier alone would hide why a gene ranks highly.

Reviewer verdict: accept. The presentation layer is a meaningful usability improvement and is methodologically safer than showing a bare 0-100 score. Keep the probability caveat visible in any manuscript, README, or UI copy.

## Iteration 32 comparator expansion and benchmark triage - adversarial methods review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new comparison evidence supports stronger claims.

Major concerns:

1. IFN is a favorable benchmark for early prioritization, not a universal win. PIPER-DEG score improves recall@50 over Stouffer/Fisher/MetaVolcanoR/RRA, but Fisher and MetaVolcanoR recover more locked ISGs by top100. The manuscript must report both cutoffs.
2. hStouffer and AWmeta remain blocked for faithful-input reasons. This strengthens the as-published heterogeneous-DEG use-case rationale, but it is not evidence that PIPER outperforms those methods.
3. Indisulam is not a broad benchmark. Recovering `RBM39` and `TYMS` is mechanistically useful, but the broader anchor panel is weakly recovered and should remain a cautionary/mechanistic case.
4. The best-signal source-unit sensitivity is useful for time-course diagnosis, but it must not replace source-unit mean aggregation without permutation/null calibration because it can reintroduce winner-selection bias.
5. ER stress/UPR is a plausible favorable next benchmark, but it must be locked before scoring. Acute/late time windows and stress-overlap exclusions are critical to avoid post-hoc topic shopping.

Minor concerns:

- The generic gold-panel comparator script improves reuse and reduces ad hoc benchmark summaries.
- Deterministic source-unit representative labels improve auditability without changing scores.
- The current no-go dossier is more complete because it audits compatible subsets (`HYP015`, `HYP020`) rather than merely reporting a whole-corpus block.

Reviewer verdict: accept the iteration as a comparative-evidence improvement. The defensible claim is now "PIPER can improve early top-ranked prioritization in favorable heterogeneous DEG settings and remains usable where variance-dependent direct prior-art tools are not faithfully runnable." Do not claim global SOTA superiority.

## Iteration 33 microarray/NutriOmics integration positioning - adversarial methods review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the microarray extension creates a publishable advantage without overstating comparative performance.

Major concerns:

1. The microarray sensitivity case does not support a broad performance-superiority claim. PIPER score recovers only `TYMS` in the top100 anchor panel, while rank-product approximation and RobustRankAggreg recover both `RBM39` and `TYMS`.
2. The stronger claim is a resource/workflow claim: PIPER preserves assay type, platform, normalization, source input type, probe-collapse rule, source-unit support, and direct-vs-mechanism labels in a queryable evidence DB.
3. This is a good fit for a NutriOmics-style database because the user-facing question is usually directional consistency of a nutrient-gene response across independent evidence, not whether a single meta-analysis p-value is formally calibrated.
4. The GSE93829 microarray layer must remain labeled as E7820/RBM39-axis sensitivity evidence, not direct Indisulam treatment evidence. It cannot replace the primary Indisulam RNA-seq result.
5. The current microarray DEG derivation from processed normalized matrices is acceptable as a demonstration path, but limma full tables remain the preferred microarray input when available.

Minor concerns:

- The explicit nutrient-gene query contract is useful and should be copied into README/API documentation later.
- Probe collapse before ranking is an important safeguard; do not allow probe-level duplicate rows into source-unit scoring.
- Rank/percentile and direction-concordance labels are safer user-facing outputs than raw p-values for this DB use case.

Reviewer verdict: accept as a strong positioning improvement. The manuscript should present this as "cross-platform, provenance-preserving nutrient-gene response atlas construction," not as "PIPER beats all meta-analysis methods on microarray data."

## Iteration 34 review-driven pipeline hardening - adversarial methods/code review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the approved implementation fixes remove hidden auditability, provenance, and user-facing interpretation risks.

Major concerns:

1. The evidence DB schema is now safer. Removing ambiguous `pvalue`/`padj` labels from `gene_evidence` and replacing them with `aggregate_pvalue`, `aggregate_padj`, `min_source_pvalue`, and `min_source_padj` prevents readers from confusing a source-unit aggregate evidence value with the best source-row p-value.
2. Source-unit evidence auditability improved materially. Mixed contributing study IDs, pipelines, assay types, platforms, source input types, normalizations, probe-collapse rules, durations, paths, and URLs are preserved instead of being collapsed to the first row.
3. The slice runner no longer contaminates unrelated topics with hypoxia-specific recall. Recall now requires a locked topic-specific Excel `GoldPanel`; otherwise the metric is explicitly marked not applicable.
4. Source provenance is stronger because sidecars now include the actual source DEG input files and the `min_studies` setting, not only the catalog.
5. The microarray fallback guardrail is appropriate. Processed matrices with one sample per group are no longer converted into confirmatory Welch DEG tables; limma full tables or clearly exploratory treatment remain the safer path.
6. Browser/API rendering is safer because DB-derived strings are escaped before insertion into `innerHTML`. Existing SQL parameterization plus output escaping is a stronger local-browser posture.

Minor concerns:

- Existing result DBs generated before iteration 34 are schema-stale; regenerate before producing manuscript tables, screenshots, or public API examples.
- README/API documentation should explicitly explain `GoldPanel`, `SLICE_MIN_STUDIES`, source sidecars, and the distinction between aggregate and minimum source p-values.
- The fixes harden the resource/DB story but do not change the empirical benchmark conclusion: IFN remains favorable, hypoxia remains parity, Indisulam remains mechanistic/cautionary, and microarray remains a cross-platform utility demonstration.

Reviewer verdict: accept the implementation hardening. This materially improves publication readiness for the "provenance-preserving heterogeneous DEG evidence DB" claim, provided all manuscript-facing DB artifacts are regenerated with the new schema before use.

Post-verification addendum: current manuscript-facing DB artifacts for hypoxia (`iter-31`), IFN, Indisulam primary, and Indisulam microarray were regenerated under the hardened schema. A large-corpus performance bottleneck in evidence metadata aggregation was found and fixed with a singleton fast path. Schema smoke checks confirmed the new aggregate/min-source p-value and contributing provenance columns in all four DBs.

## Iteration 35 Tavis-style trust-boundary hardening - adversarial security/code review

Reviewer stance: adversarial local-tool security reviewer focused on reachable parser/API/export surfaces and whether a small testcase proves the behavior.

Major findings and fixes:

1. Local API CORS was too permissive. A running `piperdeg serve` instance returned `Access-Control-Allow-Origin: *`, allowing arbitrary web pages to read local PIPER DB JSON from the user's browser. Removed the wildcard header; the same-origin UI still works.
2. SQLite serving used default read/write connections. The API now opens DBs in read-only URI mode, reducing accidental mutation/journal side effects while browsing.
3. API query text was unbounded. Long `q` strings could force unnecessary scans on large DBs. The API now rejects overlong text parameters and normalizes non-finite numeric params to defaults.
4. Resource-package API capture accepted arbitrary URL schemes/hosts. This could read local `file://` paths or remote attacker-controlled JSON if a bad command was copied into a manuscript export workflow. Capture is now local HTTP(S)-only.
5. Resource `top_n` accepted zero/negative values. This could silently produce empty or surprising top-gene tables. `top_n` now must be at least 1.
6. Provenance `.source` commands could contain newline-delimited injected commands if user-controlled paths were ever placed directly into command strings. Commands are now validated as single-line, and the main user-facing regeneration command builders use shell quoting.

Impact assessment:

- Severity is local-tool moderate, not remote RCE. The highest-impact realistic chain was data exposure from a running local browser/API to a malicious webpage via permissive CORS.
- The resource-package URL issue required the user or automation to pass a malicious `api_url`; it is now blocked by default.
- Provenance command injection required a user to copy/execute a generated `.source` command; single-line validation and quoting reduce that risk in the main workflows.

Reviewer verdict: accept. The fixes remove the main cross-origin local DB exposure and harden adjacent export/provenance surfaces without changing the scientific score algorithm or existing benchmark interpretation.

## Iteration 36 ER stress/UPR benchmark and comparison - adversarial methods review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new ER stress benchmark strengthens the comparison without cherry-picking or overclaiming.

Major concerns:

1. The ER stress gold panel was appropriately locked before scoring. This avoids the most important benchmark-validity failure mode.
2. The primary-quality ER stress tier is favorable to PIPER-DEG, but the claim is narrow. PIPER score improves early recall over Fisher/MetaVolcanoR and weighted Stouffer at top10/top20/top50; at top100 it ties weighted Stouffer and only slightly beats Fisher/MetaVolcanoR.
3. The full sensitivity tier is unfavorable to PIPER-DEG. Adding the lower-confidence normalized-matrix source drops PIPER score sharply while Fisher/MetaVolcanoR remain much higher. This cannot be hidden or relegated to an unreported exploratory check.
4. The primary/full split is defensible only because `GSE103667` was pre-labeled as lower-confidence normalized-matrix evidence with two replicates per arm. The manuscript must describe the tier rule as input-quality based, not performance based.
5. hStouffer and AWmeta remain faithful-input blockers, not wins. The ER stress benchmark does not change that interpretation.
6. The benchmark strengthens a useful PIPER story: early prioritization and provenance-preserving evidence-DB construction under heterogeneous table quality. It does not establish global SOTA superiority.

Minor concerns:

- `GSE84989` was correctly deferred because sample-label metadata was not locked before scoring.
- The current comparator summary is still gene-membership recall, not direction-aware recall. Since the gold panel has expected directions, the next hardening step should count expected-up genes only when recovered in the expected direction.
- A source-support-aware metric would align better with the NutriOmics/evidence-DB use case than p-value-style recall alone.

Reviewer verdict: accept with a required framing caveat. ER stress/UPR is now the best favorable benchmark in hand, but only the primary-quality tier supports that favorable comparison. The full sensitivity result is a genuine limitation and should be used to argue for transparent source-quality tiers rather than suppressed.

## Iteration 37 direction-aware sensitivity reporting - adversarial methods/code review

Reviewer stance: Scientific Reports methods reviewer focused on whether the sensitivity improvement increases interpretability without post-hoc retuning.

Findings:

1. The implementation improves evaluation rather than tuning PIPER. This is methodologically safer: the primary score formula stays `piper_score_v1_2_source_unit_mean`, and the new outputs are comparator/report diagnostics.
2. Direction-aware recall is appropriate because the ER stress gold panel declares expected up-regulation. It prevents simple membership recall from giving full credit to genes recovered in the wrong direction.
3. The new full-sensitivity result is informative: Fisher/MetaVolcanoR still beat PIPER on membership recall, but direction-aware recall exposes a wrong-direction `HSP90B1` recovery. This supports PIPER's direction-consensus value without reversing the full-sensitivity limitation.
4. Source-quality tier effect is now explicit and large for PIPER. Primary-minus-full deltas of +0.28 to +0.56 show that PIPER sensitivity is strongly affected by lower-confidence normalized-matrix evidence. This must be discussed as both a rationale for source-quality tiering and a limitation.
5. The updated tests cover the key failure mode: a method can recover all gold genes by membership while only partially satisfying expected direction.

Reviewer verdict: accept. This strengthens the sensitivity strategy without p-hacking the score. The next defensible metric is source-support-aware recall or a query-focused DB metric, not more score-weight tuning.

## Iteration 38 source-quality weighted sensitivity - adversarial methods/code review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new weighting strategy fixes a real weakness without creating post-hoc overclaiming.

Findings:

1. The implementation correctly preserves the primary score. `piper_score_v1_2_source_unit_mean` remains unchanged and still reports the full-sensitivity limitation. This is important: the unfavorable primary full-set result is not hidden.
2. The source-quality diagnosis is credible because the failing source was already pre-labeled as lower-confidence normalized-matrix evidence. The new diagnostics quantify the issue rather than inventing it after the fact: `GSE103667_THAP` has low static quality, very low median pairwise LFC Spearman, and is marked `recommended_role=sensitivity`.
3. The quality-weighted secondary score is useful but must remain secondary. It recovers many ER stress UPR genes in the noisy full set and improves direction-aware recall over Fisher/MetaVolcanoR at @50/@100, but it is still a heuristic quality-weighted ranking and not a calibrated inferential model.
4. The result strengthens the NutriOmics/evidence-DB positioning. Users often need a robust discovery view when public sources have mixed quality; the new output lets them inspect primary consensus, quality-weighted ranking, and source-level evidence instead of receiving one opaque p-value.
5. The source-coherence guardrail is gold-panel-free, which reduces p-hacking risk. However, it is still data-adaptive and should not be promoted as a universal scoring rule until validated across additional locked topics.
6. The comparator result is now more nuanced:
   - Primary full-set PIPER remains unfavorable: recall@50/@100 = 0.11/0.17.
   - Quality-weighted secondary PIPER improves to 0.61/0.67 and has no wrong-direction ER stress gold hits at @100.
   - Fisher/MetaVolcanoR retain high membership recall but have lower direction-aware recall because of `HSP90B1:down!=up`.

Required manuscript/API wording:

- Call the new score `quality-weighted secondary ranking` or `sensitivity/discovery ranking`.
- Do not call it the main PIPER score unless additional validation promotes it.
- Report source diagnostics and primary/full tier results whenever the secondary score is shown.
- State that normalized-matrix sources can be useful exploratory evidence but should not silently override full-table sources.

Reviewer verdict: accept with labeling constraint. This materially improves the weakness exposed by the ER stress full-sensitivity set and strengthens the evidence-DB story. It does not yet justify a claim of universal superiority or replacement of the primary score.

## Iteration 39 IFN locked-mechanism rerun - adversarial methods review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether replacing Indisulam with IFN reduces benchmark ambiguity without creating a new overclaim.

Findings:

1. The user's objection is valid. Indisulam is a drug perturbation with context-dependent downstream transcriptional effects, so it is weak as a broad method benchmark. It should remain a mechanistic/cautionary case, not a main superiority test.
2. Type-I interferon response is a cleaner benchmark choice. The expected biology is direct and directional: IFN-alpha/beta stimulation should induce canonical ISGs through JAK/STAT signaling. The locked gold panel has explicit `expected_direction=up`, so direction-aware recall is meaningful.
3. The IFN result is favorable but not universal. PIPER primary and quality-weighted secondary recover 0.85 of the locked ISG panel by top50, compared with 0.70 for weighted Stouffer, Fisher, MetaVolcanoR, and the approximate rank-product baseline. This supports early prioritization.
4. The top100 result prevents overclaiming. Fisher and MetaVolcanoR reach 0.95 at top100, while PIPER reaches 0.90. The manuscript must not claim global superiority on IFN; it can claim stronger early ranking and direction-consistent prioritization.
5. The corpus is still a derived-count pilot, not an as-published DEG-table benchmark. The generator and summary correctly label it as such. This is acceptable for a locked-mechanism validation axis, but it should not be used as evidence that PIPER handles heterogeneous published DEG tables better than baselines.
6. The IFN source-quality metadata fix was necessary. Without `source_input_type=derived_count_table`, the new quality diagnostics mislabeled the source units as low quality. After regeneration, both source units are medium-quality primary sources with no outlier flags.
7. hStouffer and AWmeta remain blocked comparators. Their blocked status must remain explicit and cannot be counted as a win.

Reviewer verdict: accept with framing limits. IFN is now the best non-drug locked-mechanism benchmark alongside ER stress/UPR. It strengthens the paper as an evidence-DB prioritization method, but only under the claim that PIPER improves early, direction-aware recovery, not that it dominates every cutoff or every method.

## Iteration 40 IFN top1000 recall curve - adversarial figure/methods review

Reviewer stance: Scientific Reports methods/figure reviewer focused on whether the new curve supports the intended claim without overstating the result.

Findings:

1. The figure is appropriate and more informative than fixed cutoff tables. It shows the whole cumulative recovery path from top1 to top1000, making it harder to cherry-pick a favorable cutoff.
2. The result supports early prioritization. PIPER primary and quality-weighted secondary recover 17/20 locked ISGs by top50, while weighted Stouffer, Fisher, MetaVolcanoR, and rank-product approximation recover 14/20 and RobustRankAggreg recovers 13/20.
3. The same curve limits the claim. By top500/top1000, nearly all plotted methods converge to 19/20. This means the manuscript should not say PIPER has better broad recovery; the claim must be that PIPER ranks canonical IFN genes earlier.
4. Direction-aware recall does not change the IFN result because recovered gold genes are in the expected up direction. This is still worth plotting because it keeps the metric consistent with ER stress/UPR, where wrong-direction recovery mattered.
5. The plotted source data and manifest are present, and the figure package exports PNG/PDF/SVG plus XLSX source data and a DOCX legend. This satisfies the reproducible figure requirement.
6. The limitation remains: IFN is a derived-count pilot, not an as-published DEG-table validation corpus. The figure legend and manuscript text must keep that caveat.

Reviewer verdict: accept. Use this figure for early-prioritization evidence. Do not use it as a global SOTA figure.

## Iteration 41 PIPER score model upgrade - adversarial methods/code review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new score fields improve interpretability without creating a new overclaim or hidden configuration risk.

Findings:

1. The split between `piper_score`, `priority_score`, and `evidence_reliability_score` is conceptually useful. It separates ranking priority from reliability/stability, which is better aligned with a browsing database than a single opaque meta-analysis p-value.
2. The new direction-confidence and leave-one-source-out fields improve auditability. A user can now distinguish a high-ranking gene supported by both source units from a gene whose rank depends heavily on one source.
3. The source reliability shrinkage is acceptable as a secondary discovery weighting, not as a calibrated probability. The metadata and README correctly keep this distinction.
4. Time-course handling is better because it is predeclared by `time_course_mode`. The review found and fixed a silent fallback risk: invalid labels now raise a config error instead of behaving like `mean`.
5. The first IFN DB regeneration revealed a real performance issue in metadata aggregation. The fix is necessary and appropriate: source-unit metadata is reused while gene-source-specific contributor fields remain explicit.
6. The AURC additions make the benchmark harder to cherry-pick. They also constrain the claim: PIPER remains strongest at early IFN prioritization, while Fisher/MetaVolcanoR remain slightly better by broad AURC@1000.
7. Browser/API/resource-package exposure is now consistent with the database fields, so the implementation supports the intended local evidence-DB workflow.

Remaining caveats:

- Leave-one-source-out stability may need another performance pass on much larger corpora with many source units.
- `quality_weighted_piper_score` should remain secondary until validated across additional locked topics.
- None of the new indices should be described as posterior probabilities or inferential p-values.

Reviewer verdict: accept. The upgrade materially improves PIPER-DEG as an evidence prioritization/database tool, provided the manuscript preserves the primary-vs-secondary score distinction and keeps IFN claims focused on early ranking.

## Iteration 42 source-support-aware reporting - adversarial methods/code review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new source-support layer strengthens the evidence-DB claim without post-hoc benchmark inflation.

Findings:

1. The new report addresses the correct estimand for the proposed NutriOmics-style database use case. It asks whether a gene is high-ranked, directionally consistent, and supported by at least two source units, not whether a combined p-value is "significant."
2. The locked validation rows remain clean. `locked_membership_recall`, `locked_direction_recall`, and `locked_source_supported_recall` use only the pre-existing gold panels, so these rows can be used as benchmark diagnostics.
3. The interpretive marker rows are acceptable only because they are explicitly labeled post-output and non-benchmark. This prevents `CMPK2`, `TAP1`, `GBP1`, `CRELD2`, and `PTX3` annotations from quietly expanding the locked panels.
4. IFN is now stronger biologically: 7/10 top genes are locked ISGs and the remaining 3 top10 genes are interpretively supported IFN biology. All top10 annotated genes point up and have full source-unit concordance.
5. ER stress primary remains favorable but not perfect. Seven top10 genes are locked UPR genes and two more are interpretively supported/contextual, but `ADAM19` remains unannotated and should not be described as a textbook UPR marker.
6. ER stress full sensitivity remains negative for the primary score. This is good from an integrity standpoint: the source-support report did not hide the noisy normalized-matrix weakness.
7. Regenerating ER score DBs with the current schema was necessary. The reports now include reliability fields consistently instead of mixing old ER CSVs with new IFN score outputs.
8. The current quality-weighted full sensitivity result is weaker than earlier notes after regeneration: recall@20/@50/@100 is 0.33/0.50/0.61. The manuscript must use the regenerated numbers, not older stronger numbers.

Code/QA review:

- The new script has clear input validation for empty gold panels, missing columns, invalid top N, invalid source-unit thresholds, and invalid concordance bounds.
- Source sidecars and JSON provenance are written for all report artifacts.
- The summary CSV encodes the guardrail in the metric names and criteria field.
- Tests cover the key silent failure mode: membership recovery can pass while direction/source-supported recovery fails.

Remaining caveats:

- The interpretive panels are manually curated and should be citation-backed in the manuscript text if used.
- Source-support recall is PIPER-specific and should not be presented as a symmetric comparator metric against Fisher/Stouffer/RRA, because those methods do not expose per-source directional support in the same output.
- `min_source_units=2` is appropriate for current pilot corpora but should be predeclared per corpus when larger topics are added.

Reviewer verdict: accept. This is a meaningful upgrade for the evidence-DB framing and strengthens publication readiness, provided locked benchmark and interpretive annotation are kept separate in all figures and text.

## Iteration 43 heat-shock/HSF1 benchmark - adversarial methods review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the new benchmark is genuinely IFN-grade and whether the curation choices avoid cherry-picking.

Findings:

1. Heat shock / HSF1 is an appropriate clean-mechanism benchmark. The top recovered genes, `HSPA1A`, `HSPA1B`, `HSPA6`, `DNAJB1`, `BAG3`, and `HSPH1`, are textbook heat-shock response genes and all point in the expected up direction.
2. The locked gold panel is defensible because it was written before PIPER or baseline scoring and uses a compact expected-up chaperone/co-chaperone panel rather than expanding after seeing the output.
3. The active source set is conservative. Excluding `GSE73471` for n=1 control, `GSE130493` for no matched control, and `GSE57397` for no untreated/non-heat arm is methodologically correct.
4. Excluding `GSE123980` from primary use is also reasonable. Its total RNA sample titles indicate heat-shock time points, but total RNA treatment characteristics do not explicitly encode heat exposure, and the clearer TT-seq rows are nascent transcription rather than the main expression-table target.
5. The result supports early prioritization rather than universal superiority. PIPER improves recall@50 and recall@100 over Fisher/MetaVolcanoR/weighted Stouffer, but several methods recover many gold genes by larger K values.
6. PIPER's top10 contains several non-locked but biologically plausible heat response genes (`HSPA7`, `HSPA1L`, `DNAJB4`, `DNAJA4`, `CRYAB`, `ZFAND2A`). These can be discussed as interpretive biology only if separately curated/cited; they must not be counted as locked benchmark hits.
7. The benchmark remains derived-count. It validates the scoring behavior on clean public count matrices, not the complete as-published heterogeneous DEG-table scenario.
8. hStouffer/AWmeta are still blocked faithful-input comparators, so their blocked rows must remain guardrails, not wins.

Reviewer verdict: accept with framing limits. This is a strong second clean-mechanism benchmark alongside IFN and materially improves the paper's validation story. The manuscript should describe it as a locked-mechanism derived-count benchmark supporting early, direction-aware prioritization, not as a general SOTA dominance claim.

## Iteration 44 coding-lord full-code audit - adversarial code/methods review

Reviewer stance: final-integrator code reviewer combining biostatistical unit-of-analysis, debugging, security, data-pipeline, and beginner-UX lenses.

Findings:

1. The accepted critical issue was evidence/score drift under predeclared time-course handling. Before the fix, `collapse_gene_source_units()` could select only the early/late/peak rows for scoring while `study_gene_evidence()` still reported `min_source_pvalue`, contributor IDs, durations, paths, and URLs from all same-source rows. This weakened the SQLite browser's auditability claim.
2. The fix is methodologically appropriate: evidence metadata is now derived from `source_unit_rows_for_aggregation()`, the same selected-row preparation used before source-unit collapse.
3. The replicate-count fix is conservative. When one source unit contains rows with conflicting `n_ctrl`/`n_treat`, quality weighting now uses the minimum observed replicate count instead of whichever row sorted first.
4. The config whitespace fix improves beginner failure modes. A required Excel/CSV cell containing spaces now triggers the direct "empty required values" message rather than a later misleading path or source-column error.
5. The patch avoids broad algorithm changes. It does not alter gold panels, benchmark interpretation, score weights, or manuscript claims.
6. No API SQL-injection or obvious browser XSS issue was found in inspected local SQLite/HTML paths; existing parameterized SQL and escaping remain in place.
7. Citation and novelty expert passes produced no code-actionable finding for this turn because no manuscript or citation claim was changed.

Reviewer verdict: accept. The patch fixes a real auditability bug and a source-quality edge case with narrow blast radius, and the new regressions directly protect both behaviors.

## Iteration 45 biostatinfo expert audit - adversarial methods/results review

Reviewer stance: biostatistics and bioinformatics reviewer focused on unit of analysis, input filtering, result provenance, benchmark contamination, and claim calibration.

Findings:

1. The main estimand remains appropriate only if described as directional evidence prioritization across independent source units, not as calibrated formal meta-analysis significance.
2. Derived-count RNA-seq benchmark inputs needed an explicit low-count expression filter. Adding a label-independent raw-count filter before logCPM/Welch derivation is necessary and statistically preferable to letting low-expression rows enter rank calculations.
3. Restricting count-derived benchmark universes to protein-coding genes is appropriate because the current locked panels and claims are gene-symbol/protein-coding focused. This reduced noncoding artifacts in heat-shock and ER stress derived outputs.
4. Conservative parsing of mixed `table_scope` values is required. A source unit marked `deg_only;full_results` should not receive full-results quality credit.
5. Exposing `table_scope` in score evidence rows improves auditability because users can now inspect the same quality field that affected source weighting.
6. The most important results bug was stale baseline contamination. Comparator and recall-curve scripts previously globbed all TSVs in a baseline directory, so old corpus-name baseline files could silently duplicate or distort benchmark rows. The manifest-authoritative fix is correct.
7. IFN remains strong: PIPER top10/top50/top100 locked recall is 0.35/0.85/0.90 with canonical ISGs in the top ranks.
8. Heat-shock/HSF1 is now cleaner after protein-coding filtering: PIPER top10/top50/top100 recall is 0.38/0.69/0.75, with HSPA/HSPH/DNAJ/BAG3 genes dominating the top ranks.
9. ER primary remains supportive: PIPER top10/top50/top100 recall is 0.39/0.72/0.72 and direction-aware recovery matches membership recovery.
10. ER full remains a limitation. Primary PIPER drops to 0.11/0.11/0.17 at top10/top50/top100, while quality-weighted secondary ranking improves to 0.28/0.61/0.61. This should be reported transparently as evidence for source-quality sensitivity, not hidden.
11. hStouffer and AWmeta remain blocked comparators because faithful variance/SE inputs are missing. They are not defeated baselines.

Reviewer verdict: accept after fixes. The current code/results are much more defensible, but the manuscript claim must stay focused: PIPER-DEG is useful because it builds an auditable directional evidence database across heterogeneous DEG sources. It should not claim universal SOTA statistical superiority.

## Iteration 46 hypoxia/HIF1 benchmark and SOTA comparator gap review

Reviewer stance: Scientific Reports methods/resource reviewer focused on whether the added benchmark is a real gold-standard mechanism and whether the comparator slate anticipates obvious reviewer objections.

Findings:

1. Hypoxia/HIF1 is an appropriate additional benchmark. The biology is externally recognizable, and the locked panel is anchored to expected-up hypoxia/HIF response genes.
2. The benchmark is not a cherry-picked PIPER win. PIPER quality-weighted reaches 0.15/0.35/0.65/0.80 recall at top10/top20/top50/top100, while weighted Stouffer and rank/p-value baselines are competitive or better at early cutoffs.
3. The biological top ranks are still strong: `NDRG1`, `EGLN3`, `AK4`, `ADM`, `HK2`, `VEGFA`, `SLC2A1`, `BNIP3`, and `SLC2A3` are credible hypoxia/HIF-axis signals.
4. This benchmark should be used as balanced validation: PIPER recovers canonical genes and offers auditable source-unit support, but the manuscript must not claim universal superiority over Stouffer/Fisher/rank methods.
5. The current comparator slate is methodologically decent but not package-complete. A reviewer can reasonably ask for AWFisher and exact RankProd because both are named, maintained, and directly relevant to p-value/rank meta-analysis.
6. MetaDE/DExMA should be addressed explicitly even if not fully run, because they are umbrella gene-expression meta-analysis packages that include many method families.
7. MetaIntegrator is a relevant raw-expression/effect-size comparator but is not a faithful baseline for all current heterogeneous DEG-summary inputs unless variance/SE or raw expression matrices are available.
8. hStouffer/AWmeta blocked rows remain valid guardrails. They should stay described as faithful-input blockers, not failed methods or PIPER wins.

Reviewer verdict: accept with comparator follow-up. Hypoxia/HIF1 strengthens benchmark breadth, and the SOTA gap review identifies the next two high-value additions: AWFisher and exact RankProd.

## Iteration 47 comparator coverage and public-summary feasibility table - adversarial methods review

Reviewer stance: skeptical bioinformatics methods reviewer asking whether PIPER-DEG is necessary when existing meta-analysis packages exist, and whether blocked comparator rows are being used as unfair wins.

Findings:

1. Adding AWFisher is a high-value fix. It directly covers an adaptive p-value method a reviewer could reasonably request, and it can be run from public per-study p-value matrices without raw expression.
2. Adding metaRNASeq Fisher is useful but should be framed as p-value-combination parity, not a distinct conceptual win over Fisher/MetaVolcanoR. It mostly confirms the classical p-value-combination lane.
3. metaRNASeq inverse-normal is correctly caveated. Public supplementary DEG tables do not provide uniform one-sided p-values or complete replicate metadata for every source, so the adapter is a documented public-summary adaptation rather than a perfect package-native benchmark.
4. The generated `public_summary_tool_input_requirements.csv` is exactly the right reviewer-facing artifact. It makes the limitation explicit: some packages are not "missing because ignored"; they are blocked because the public DEG-only evidence lacks required raw expression, phenotype, DESeq2-specific, or variance/SE fields.
5. Exact RankProd remains unresolved but honestly separated. Keeping `rank_product_approx` distinct from `rankprod_exact` prevents an overclaim.
6. MetaDE is handled correctly as partial coverage: p-value/rank families are represented by primitive runnable baselines, while effect-size modes require variance or raw-expression-derived inputs.
7. DExMA and MetaIntegrator should not be forced onto the current benchmark. Running them from already-filtered public DEG tables would be methodologically misleading.
8. hStouffer and AWmeta remain faithful-input blockers, not defeated methods. This guardrail must be preserved in the manuscript text.
9. The new results remain balanced. PIPER is strong on IFN and heat-shock top50/top100, but classical methods remain competitive for hypoxia and some early cutoffs. This is a credibility strength if reported honestly.
10. The Indisulam results remain weak across methods and should not be promoted as a positive benchmark. It is better used as a cautionary drug-mechanism/sensitivity example.

Reviewer verdict: accept. The comparator and supplementary-input story is now much more publication-ready. The next writing task is to turn the generated feasibility table into a clean Supplementary Table and ensure the Results section says "runnable where faithful; blocked where required inputs are absent" rather than implying PIPER beats tools that could not be run.

## Iteration 48 deep benchmark metrics and evidence-card package - adversarial methods review

Reviewer stance: skeptical Scientific Reports reviewer asking whether the validation is deep enough, whether PIPER has a visible advantage beyond recall@K, and whether blocked SOTA tools are being handled fairly.

Findings:

1. Adding AURC@1000 improves the benchmark because it uses the full early-prioritization curve rather than a few arbitrary cutoffs.
2. Direction-aware recall is important and should stay in the Results because PIPER's intended use is topic-specific directional evidence, not only membership recovery.
3. Source-unit bootstrap summaries are useful robustness evidence for summary-combination comparators, but the manuscript should not imply these are calibrated PIPER p-value CIs.
4. PIPER's stronger unique metrics are the top100 source-unit count, sign concordance, quality-weighted support, and leave-one-source-unit-out stability. These are the right differentiators against Fisher/Stouffer/rank combiners.
5. IFN remains the cleanest positive benchmark: PIPER quality-weighted top10/top50/top100 is 0.35/0.85/0.90 with canonical ISGs at the top.
6. ER primary is supportive: PIPER quality-weighted top10/top50/top100 is 0.39/0.72/0.72 and AURC@1000 is high.
7. Heat shock is supportive but nuanced: PIPER improves top50/top100 recovery, while Fisher-like methods remain competitive at top10.
8. Hypoxia/HIF1 is balanced, not a PIPER-only win. This is acceptable and useful because it shows the method is not being cherry-picked.
9. ER full noisy data should be presented as sensitivity/limitation evidence. Quality weighting helps, but the noisy corpus still warns that source inclusion matters.
10. The source-level feasibility matrices are a strong supplement. They answer the reviewer question "why not DExMA/MetaIntegrator/RankProd/hStouffer/AWmeta?" with actual missing-input columns rather than hand-waving.
11. Exact RankProd remains unresolved but appropriately blocked. Do not call `rank_product_approx` RankProd in manuscript text.
12. The evidence-card figure is valuable because it shows the database/resource advantage directly: readers can see rank percentile, source count, sign concordance, and per-source log2FC direction for canonical genes.

Reviewer verdict: accept with calibrated claims. The current evidence supports a Scientific Reports-style methods/resource paper if the manuscript says PIPER-DEG is an auditable directional prioritization/database tool for heterogeneous public DEG tables. It should not claim universal SOTA dominance or formal calibrated meta-analysis significance.

## Iteration 49 prior-art comparator/resource coverage - adversarial novelty review

Reviewer stance: skeptical novelty and SOTA-comparator reviewer asking whether PIPER-DEG ignored obvious public expression platforms, raw-expression meta-analysis workflows, signature databases, or generic p-value-combination packages.

Findings:

1. The new comparator/resource coverage closes the main novelty-audit gap. OMiCC, ImaGEO, NetworkAnalyst, crossmeta, DEET, CREEDS, and generic p-value-combiner families are now explicitly represented in generated tables.
2. The classification is fair. OMiCC, ImaGEO, NetworkAnalyst, and crossmeta are relevant prior art, but they require expression matrices, GEO workflow inputs, sample labels, or platform annotation; they are not faithful same-input baselines for heterogeneous supplementary DEG tables alone.
3. DEET and CREEDS are correctly framed as public expression signature resources/databases. They are important context for Related Work, but they do not replace PIPER-DEG's local, topic-specific integration of user-supplied public DEG tables.
4. Generic p-value-combiner packages are not ignored. The method family is covered by Fisher, Stouffer, AWFisher, and metaRNASeq-style rows, with input-faithfulness caveats where one-sided p-values or package-native fields are unavailable.
5. The new `prior_art_coverage_summary.csv` is the right high-level supplementary artifact because it prevents the detailed corpus table from becoming the only place where prior-art coverage is visible.
6. The source-unit feasibility matrices are now stronger reviewer evidence than prose alone: each public source can be inspected for whether a comparator would need missing raw-expression or metadata inputs.
7. The manuscript should still avoid claiming that PIPER-DEG beats methods that could not be run faithfully. The safe claim is public-file compatibility plus auditability, not universal statistical superiority.
8. Exact RankProd, DExMA, MetaIntegrator, hStouffer, AWmeta, OMiCC, ImaGEO, NetworkAnalyst, and crossmeta should remain different-input or blocked rows until raw expression/sample metadata are acquired and a separate raw-expression benchmark is designed.

Reviewer verdict: accept. This resolves the likely "why not existing SOTA tools/resources?" critique at the table and methods-framing level. The remaining risk is rhetorical overclaiming in the manuscript; the generated tables support a strong Scientific Reports-style methods/resource paper if the claims stay calibrated.

## Iteration 50 cross-platform microarray benchmark - adversarial figure/result review

Reviewer stance: skeptical methods/resource reviewer asking whether PIPER-DEG has a visible advantage beyond RNA-seq-only benchmarking and whether microarray integration creates a real, auditable signal rather than source inflation.

Findings:

1. Removing Indisulam from the positive benchmark story is correct. Its anchor recovery remains weak and drug-specific, so it should not carry the manuscript's method validation.
2. The IFN mixed benchmark is supportive. Adding GSE71634 microarray improves top10, top20, and top100 recall, and moves canonical ISGs such as `IFI27`, `IFI44L`, `MX2`, and `IFI6` upward.
3. The ER stress mixed benchmark is supportive at practical broad-prioritization cutoffs. Top50/top100 recall improves from 0.72 to 0.89, and UPR genes such as `EDEM1`, `ATF4`, `ASNS`, `DNAJB9`, and `TRIB3` move upward.
4. ER top10 decreases from 0.39 to 0.22, so the manuscript must not say microarray addition improves every cutoff. The fair claim is sharper recovery by top50/top100 and stronger cross-platform source support.
5. The source-unit treatment is conservative. GSE19519 tunicamycin and thapsigargin rows share one `paper_id`, so the two microarray contrasts do not count as two independent papers.
6. The generated figure package is publication-usable because it includes combined and standalone PNG/PDF/SVG panels, exact source data, a manifest, a DOCX legend, and validation output.
7. The microarray DEG derivation remains a fallback from processed GEO series matrices. The Methods should state that limma full tables are preferred when available, and that these rows are used to demonstrate public processed-matrix integration.

Reviewer verdict: accept with precise wording. The new figure supports PIPER-DEG's cross-platform evidence-database value: public RNA-seq and microarray DEG evidence can be merged without source inflation, and canonical markers become better supported at useful cutoffs. Do not claim universal improvement at top10 or formal calibrated meta-analysis.

## Iteration 51 manuscript figure/table package - adversarial package review

Reviewer stance: Scientific Reports methods/resource reviewer checking whether the manuscript figure, table, and supplementary package is traceable, reproducible, and claim-calibrated.

Findings:

1. The four-figure structure is coherent: Figure 1 explains the method, Figure 2 tests same-input benchmark behavior, Figure 3 demonstrates cross-platform RNA-seq + microarray integration, and Figure 4 shows the database/resource output.
2. The package follows the data-figure requirements. Each main figure has PNG/PDF/SVG outputs, source data, manifest, DOCX legend, validation notes, and source sidecars.
3. Figure 2 is balanced rather than overfit. PIPER-DEG is favorable for IFN, ER primary, and heat shock at useful cutoffs, while hypoxia/HIF1 remains competitive but not a universal win.
4. Figure 3 is manuscript-useful only if the caveat stays explicit: IFN improves at top10/top20/top100, ER stress improves at top50/top100, but ER top10 decreases after adding microarray.
5. Figure 4 strengthens the resource-paper angle because it turns the method into a local searchable atlas and SQLite database, not just a ranked-gene table.
6. The supplementary table index addresses likely reviewer concerns: same-input runnable comparators, different-input raw-expression workflows, resource/database prior art, cross-platform sources, marker rank shifts, and per-source feasibility matrices are all mapped to concrete CSV artifacts.
7. Indisulam remains excluded from the positive figure set, which avoids a weak and drug-specific validation story.
8. The remaining risk is manuscript rhetoric, not the figure package. The paper should avoid saying PIPER-DEG is a formal calibrated meta-analysis significance method or universally superior to every SOTA method.

Reviewer verdict: accept. The figure/table scaffold is ready for manuscript drafting if claims remain centered on public-DEG evidence integration, auditability, directional consistency, cross-platform support, and local DB usability.

## Iteration 52 debugging-expert defensive code review

Reviewer stance: defensive debugging reviewer focused on silent scientific logic errors, HTML/resource boundary failures, malformed SQL edge cases, and invalid benchmark inputs.

Findings:

1. The consensus `weighted_lfc` calculation had a real silent bug for mixed evidence where signed-z/rank were present but LFC was missing. The old denominator used all source-unit weights, so valid LFC evidence could be diluted by LFC-missing sources. The fix uses only valid-LFC weights.
2. The source-unit collapse LFC-valid mask needed to be computed after deterministic sorting. This is a narrow but correct guard against row-order-sensitive missing-LFC handling.
3. The static atlas had a credible script-boundary issue because raw JSON was injected into a `<script>` tag. Public DEG metadata should be treated as untrusted display/input data even in a local HTML file. Escaping JSON control characters is the right minimal fix.
4. Empty top-gene/evidence-card paths should not create SQLite `IN ()` queries. Returning schema-valid empty frames is safer and easier for downstream code to consume.
5. Empty gold panels must fail fast. A benchmark table with no positives is invalid, and silently emitting zero recall would be misleading.
6. Shell-quoted provenance commands are preferable for reproducibility because generated paths/titles may contain spaces.
7. The fixes are narrow and do not alter gold panels, benchmark datasets, method rankings, or the manuscript's calibrated interpretation.

Reviewer verdict: accept. The code is more defensible after this pass. Remaining risk is normal scope risk for a rapidly evolving research codebase, not a known blocking defect: continue to pair manuscript claims with generated source data and validation artifacts.

## Iteration 53 Claude error audit and human-only microarray upgrade

Reviewer stance: adversarial methods/statistics reviewer checking whether the manuscript-facing results are overfit to gold panels, contaminated by species mixing, or dependent on broken comparator implementations.

Findings:

1. The prior benchmark summary had a real optimism risk because it selected between primary and quality-weighted PIPER rows using the locked gold panel. Fixing the manuscript default to quality-weighted PIPER removes this test-set selection problem.
2. The high-N microarray dominance critique was valid. Capping per-contrast source weights at 4 is a pragmatic fix that preserves sample-size awareness while preventing one public matrix from overwhelming independent-source evidence.
3. The `metaRNASeq::invnorm` critique was valid. After fixing neutral missing p-values, p-value underflow still produced p=0 ties, so the final fix needed both a p-value floor and finite tie-breakers. The regenerated IFN/ER/Hypoxia invnorm outputs are no longer p=1 or p=0 alphabetical failures.
4. Human-only primary framing is the right reviewer-facing decision. Mouse ER/hypoxia rows should not be mixed into the main human benchmark by symbol uppercasing. Keeping them as excluded/sensitivity context is acceptable.
5. Replacing stale ER primary manuscript inputs with human ER RNA+microarray cross-platform results is necessary. The old ER primary output still reflected an earlier mixed-species state and should not be used as a primary manuscript benchmark.
6. Adding human hypoxia microarray is scientifically useful even though it does not improve locked-panel recall. It shows the system can integrate human microarray evidence without hiding a neutral/negative result.
7. Figure 3 should now be framed as evidence-support/cross-platform integration, not a universal performance boost. IFN improves early; ER improves at broader cutoffs; hypoxia mixed is slightly worse than RNA-only.
8. The current Table 1 is more honest: PIPER quality-weighted is favorable for IFN, ER stress, and heat shock at useful cutoffs, while hypoxia remains competitive but not a universal win.
9. The Atlas/DB output is stronger after excluding stale ER primary/full corpora from the default dashboard and adding hypoxia cross-platform.

Reviewer verdict: accept with calibrated claims. The corrected pipeline is substantially more defensible for a Scientific Reports-style methods/resource paper. The manuscript should explicitly state that quality-weighted PIPER is the fixed default, human-only active evidence is used for primary figures, processed-matrix microarray rows are fallback public-data integrations, and PIPER scores are prioritization indices rather than calibrated meta-analysis p-values.

## Iteration 54 error1/error2 audit - adversarial methods/code review

Reviewer stance: skeptical reviewer checking whether the response to `claude_error1.txt` and `claude_error2.txt` fixed code and outputs, not only prose.

Findings:

1. The source-unit heterogeneity critique was valid. Reporting only combined Stouffer direction and concordance is insufficient for auditability. Adding `heterogeneity_q`, `heterogeneity_df`, and `heterogeneity_i2` as descriptive fields is an appropriate narrow fix, provided the manuscript does not treat them as calibrated random-effects inference.
2. The Atlas reliability bug was a real generated-output bug. Reliability was already on a 0-100 scale but the HTML/SVG/matplotlib views treated it as 0-1 in places. The regenerated Atlas now displays percent-scale reliability correctly and no longer pushes SVG points off-canvas.
3. The source-unit tool feasibility matrix had a real first-study shortcut. For hStouffer and AWmeta, source-unit compatibility must require all constituent studies to have required original fields. Partial source units are now conservatively blocked as partial.
4. Publication-gate iteration selection was fragile. Explicit `iter-<n>` token parsing is safer than concatenating all digits from filename and parent directory.
5. The ER-stress report's old "use the favorable primary set" language was inappropriate after adopting quality-weighted PIPER as the default. The rewritten report now states the fixed default and preserves primary/full source sets as sensitivity evidence.
6. The README and manuscript had stale aspirational content from the autonomous methodology spec. Rewriting them to separate delivered scope from historical goals reduces reviewer confusion and prevents unimplemented-method claims.
7. Removing current-script and manuscript promotion of the excluded drug benchmark is correct. Drug-specific anchor results should not be part of the positive benchmark story.
8. The regenerated Figure 2 manifest now documents fixed quality-weighted PIPER selection rather than best-row selection, which directly addresses the benchmark cherry-picking concern.
9. Verification is adequate for this patch scope: 120 full tests passed, targeted tests cover the new heterogeneity, feasibility, Atlas scale, and report-parsing behavior, and generated artifact manifests were rechecked.

Reviewer verdict: accept. The remaining paper must still avoid overclaiming formal inferential calibration, but the main code/artifact defects in `error1/error2` have been fixed. The next high-value manuscript work is citation audit and claim-to-figure consistency checking, not another algorithm patch.

## Iteration 55 error3 audit - adversarial methods/code review

Reviewer stance: skeptical reviewer checking whether the response to `claude_error3.txt` fixed remaining comparator, species, heterogeneity, bootstrap, and microarray-design issues in code and generated outputs.

Findings:

1. The stale hypoxia primary output problem was real. The catalog had been made human-only, but generated hypoxia primary files still contained mouse rows. Regenerating the hypoxia slice fixed this: primary and cross-platform hypoxia outputs now contain only `Homo sapiens`.
2. The sparse `metaRNASeq::invnorm` critique was valid. Some public-summary corpora produce all-tie or floor-saturated ranks that should not count as successful same-input comparator runs. The blocker `metarnaseq_invnorm_uninformative_sparse_public_summary` is a conservative fix.
3. Keeping successful `invnorm` runs when they recover markers is also appropriate. IFN and heat shock remain `ok`; ER/hypoxia and IFN mixed cases with uninformative/saturated output are blocked.
4. Adding `heterogeneity_flag` improves usability. The raw I2-style field is useful but hard to scan; a high/moderate/low descriptive label makes the audit field visible without turning it into a score gate.
5. The normalized-matrix microarray limitation is now explicit. GSE19519-derived Welch contrasts are useful fallback public data, but family/pairing structure is not modeled; this must remain in Methods and source notes.
6. Adding PIPER rows to bootstrap output closes the table-structure gap for IFN, ER stress, and heat shock. Hypoxia remains computationally heavy; the deferred hypoxia PIPER bootstrap should be disclosed if bootstrap CIs are discussed.
7. The updated benchmark numbers are more conservative and should replace older text: hypoxia human-only recall@100 is 0.75, not 0.80.
8. Figure and table regeneration was necessary and completed. Figure 2, Figure 3, Figure 4, Table 1, and Table 2 now reflect the corrected hypoxia and comparator statuses.

Reviewer verdict: accept with one caveat. The implementation defects identified in `error3` are fixed or conservatively bounded. The remaining caveat is computational, not conceptual: if source-unit bootstrap CIs become a central manuscript claim, hypoxia needs either an optimized PIPER bootstrap or a clearly reported long-batch result.

## Iteration 56 external review reconsideration

Reviewer stance: publication-readiness reviewer checking whether the current manuscript still hides important caveats after the error3 fixes and paper-writer draft.

Findings:

1. The review's favorable global conclusion is credible: stale benchmark values, quality-weighted cherry-picking risk, invnorm degeneracy, human-only inconsistency, and overbroad methodology claims have been addressed in the current manuscript-facing outputs.
2. The remaining high-value manuscript issue was not another algorithm patch but caveat placement. The paper needed explicit text that high heterogeneity is frequent, especially in ER stress and hypoxia, and that heterogeneity flags are review labels rather than exclusion gates.
3. The positive-control benchmark should not imply population-level performance estimation. IFN, ER stress, and heat shock are deliberately small, biology-anchored corpora, so recall@K should remain a descriptive point summary without confidence intervals unless bootstrap or larger benchmark designs are made central.
4. Baseline comparisons are now honest if described as default/documented same-input adapters. A reviewer could still ask about tuned baselines, so the manuscript should preemptively state that tuned-method dominance is not claimed.
5. Data release, dataset citations, and clean make/CI reproduction are still required before submission. They are pre-submission execution tasks, not current scientific contradictions.

Reviewer verdict: accept with transparency edits. The manuscript is viable as a Scientific Reports-style methods/resource paper if these caveats remain visible and the final package completes citation, release, and reproducibility wiring before submission.

## Iteration 57 error4 follow-up review

Reviewer stance: publication-package reviewer checking whether the favorable `claude_error4.txt` assessment still left concrete pre-submission defects that could undermine reproducibility or supplementary-table honesty.

Findings:

1. The review was correct that publication-critical scientific contradictions were already resolved, but it was also correct that reproduction wiring was still too weak because `make figs` and `make paper` were echo stubs.
2. Replacing those stubs with real table, figure, atlas, and paper-package validation targets materially improves submission readiness.
3. The prior-art/public-summary supplementary tables still had a stale-output risk. A blocked `metarnaseq_invnorm` run could be masked by an older successful TSV unless the failure ledger took precedence.
4. The fix is conservative: ledger blockers now override stale outputs, while truly runnable `metarnaseq_invnorm` corpora remain runnable with the documented weight-proxy caveat.
5. Restricting default public-summary feasibility regeneration to the six current manuscript-facing corpora removes deprecated indisulam/stale ER-primary rows from reviewer-facing summary outputs.
6. Adding `environment.yml` and CI does not solve DOI/release or citation audit, but it closes the practical "how do I rerun the package?" gap.
7. The manuscript package validator is appropriately scoped. It checks required artifacts and known stale manuscript failure modes without claiming to replace final citation or journal-format audits.

Reviewer verdict: accept. The remaining warning is explicit and nonblocking: dataset-level citation audit is still pending before submission. No known code-path or generated-artifact blocker remains from `claude_error4.txt`.

## Iteration 58 statistical enhancement review

Reviewer stance: skeptical biostatistics reviewer checking whether `claude_stat1.txt` should change the method, the benchmark reporting, or the manuscript claims.

Findings:

1. The strongest part of the critique was correct: publication defensibility improves if the paper reports uncertainty and empirical/null enrichment around recall summaries, but the locked PIPER score should not be retuned mid-manuscript.
2. Exact recall intervals and hypergeometric enrichment/FDR are low-risk additions. They expose the small-panel uncertainty while showing that the recovered top-100 marker enrichment is far above random expectation.
3. Background-negative AUROC/AUPRC is a better benchmark supplement than ROC against three housekeeping negatives. It is still panel-dependent, but it avoids the degenerate small-negative setting.
4. The source-unit direction-confidence fix is methodologically useful. A beta-binomial source-unit concordance posterior is more interpretable than treating direction consistency as a vague scalar.
5. RE-Stouffer, RRA rho, and random-effects log2FC belong in auxiliary columns. They should be described as reporting/sensitivity lanes because most current as-published DEG files do not preserve exact standard errors and because source-unit counts are thin for IFN, ER stress, and heat shock.
6. The random-effects log2FC lane is useful but must remain explicitly approximate: current weights are derived from log2FC and two-sided p-values, not from native per-study standard errors in most files.
7. Directional performance remains incompletely tested because the locked gold panels are mostly up-regulated markers. The manuscript should not claim full up/down direction generalization until down-regulated panels are added.
8. The implemented changes strengthen the Scientific Reports-style framing: PIPER-DEG is now an auditable evidence resource with uncertainty/null-enrichment context, not just a ranked-list generator.

Reviewer verdict: accept. The statistical enhancements are correctly implemented as parallel reporting lanes. No new blocker is introduced because the locked score remains unchanged and the manuscript explicitly limits the inferential interpretation of the auxiliary fields.
