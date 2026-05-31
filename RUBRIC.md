# Publication-Readiness Rubric — Scientific Reports tuning

Three gates. All three must pass before the writing phase begins.

**Tuning rationale.** Scientific Reports evaluates on technical soundness and reproducibility rather than subjective importance. Novelty as breakthrough is NOT a gating criterion at this venue; novelty as honest scope of advance IS. The gates below reflect this: Gate 1 (reviewer clearance) and Gate 2 (quantitative checklist) carry most weight; Gate 3 (weighted score) places more weight on reproducibility/clarity/utility than on "novelty score".

---

## Gate 1 — Adversarial reviewer clearance

In Phase 4 of the final iteration, the adversarial reviewer persona produces a critique. Gate 1 passes when:

- No critique is marked "major revision required" with a technical-soundness concern
- No critique flags a core validity issue (label leakage in V7, mismatched filtering in V2, cherry-picked V6 case study, unsupported claim, fabricated reference, broken reproducibility, miscalibrated CI without acknowledgment)
- At most 3 "minor revision" points remain (SR is slightly more tolerant than top-tier journals on minor issues), each addressed in the manuscript's limitations section

The reviewer persona is instantiated fresh for the final gate — no cached judgment from mid-loop. Tell the reviewer explicitly: *"The authors claim this manuscript is ready for submission to Scientific Reports. Evaluate whether the conclusions are supported by the methods and whether the reproducibility expectations of the venue are met."*

Note the framing: not "is this important enough" but "is this supported and reproducible". This is the SR review style.

---

## Gate 2 — Quantitative checklist

Every applicable item must pass or be N/A with reason. Items marked **core** cannot be N/A.

### Core items — methodological soundness

- [ ] **V2 pipeline-perturbation experiment is apples-to-apples.** Matched min-count filtering, matched design formula, matched gene ID system across DESeq2 / edgeR / limma-voom re-analyses. Protocol documented in Methods.
- [ ] **V7 classifier evaluation has no label leakage.** Held-out 30% contains no studies/papers/research groups overlapping the training portion. Documented split.
- [ ] **V4 simulation reports CI coverage probability**, not just point recovery. Empirical coverage of reported 95% CIs is within [0.90, 0.98] for the proposed method (well-calibrated). If outside, this is reported and discussed, not hidden.
- [ ] **Baselines are tuned with the same effort as PIPER.** MetaVolcanoR mode choice, hStouffer cutoff threshold, AWmeta weight initialization documented and run at default + one tuned setting each.
- [ ] **Gold-standard target list (`METHODOLOGY_SPEC.md` §3.2) is locked**, no post-hoc additions or removals.
- [ ] **Per-study reasons for exclusion are categorical and consistent.** Studies excluded after initial catalog entry get a documented exclusion category (curation failure, raw counts unavailable, ambiguous pipeline, mixed perturbation, …). No silent drops.
- [ ] **V6 case study was pre-specified** (or the pre-specified candidate was unworkable and the next on the list was used, with documentation). Not "we tried PIPER on three case studies and report the one with most contamination."

### Core items — reproducibility (heavy weight for SR)

- [ ] **`make slice`, `make benchmark`, `make figs`, `make paper` all run end-to-end from a fresh clone.** CI tests this.
- [ ] **Every figure/table in the manuscript has a `.source` sibling file** with the exact command that produced it.
- [ ] **Raw inputs and harmonized parquets are deposited on Zenodo** with DOI before submission. Catalog CSV under CC-BY.
- [ ] **Code repo has versioned releases** matching the manuscript draft. Tag for submission.
- [ ] **`degpipeline` is PyPI-installable** with versioned releases, CI tests passing, basic README.
- [ ] **R-package versions are pinned** in `environment.yml`. rpy2 bridge tested.

### Core items — honest reporting

- [ ] **Limitations section names actual limitations**, not decorative ones. At minimum: heuristic pipeline detection, SE imputation approximation, two-domain test bed, no full Bayesian hierarchical, GEO publication bias.
- [ ] **No fabricated references** — every citation verified by `ref-finder`.
- [ ] **No fabricated data** — every DEG table is real and traceable to a real paper.
- [ ] **Negative results reported.** If V6 shows little contamination, that is in the paper. If a baseline beats PIPER on a metric, that is in the paper.
- [ ] **Effect sizes alongside p-values** in every consensus signature table. Pooled LFC + CI, not just significance.
- [ ] **Heterogeneity (I²) reported** for every consensus claim. Consensus signatures with I² > 75% flagged as "context-dependent" rather than promoted.

### Bioinformatics-weighted items

- [ ] **Multiple testing correction** documented and applied at each layer (per-study DEG already corrected; consensus signature BH across all tested genes per nutrient).
- [ ] **Pathway redundancy** addressed before reporting top-K (Jaccard ≥ 0.7 clustering).
- [ ] **Gene ID standardization** documented (Ensembl-stable + HGNC-symbol dual labeling).
- [ ] **Cohort composition** reported as a table: per-corpus N studies, species, cell systems, pipelines, time points.

### License / data items

- [ ] **License is MIT or Apache-2.0** for code; CC-BY for data.
- [ ] **No GPL contamination** of the codebase (some R packages are GPL — handled via rpy2 boundary, documented).
- [ ] **All baselines and dependencies properly cited** in Methods.

---

## Gate 3 — Weighted rubric score

Score 0–3 per item (0 = absent, 1 = weak, 2 = solid, 3 = strong). Weighted sum ≥ **0.70** × max → pass.

**SR-specific weighting differs from a Nature Methods-style rubric.** The weights below reflect SR's review priorities:

- Soundness of methodology (weight 3)
- Strength of empirical evidence (weight 3)
- Reproducibility (weight 3)
- Fit between claims and evidence (weight 3)
- Clarity and honesty of writing (weight 2)
- Honest scoping of contribution (weight 2)
- Novelty of the contribution (weight 1) — SR does not heavily weight novelty; PIPER's PARTIAL novelty status as established in `NOVELTY_AUDIT.md` is sufficient at this venue
- Utility for the field (weight 2) — does the method actually enable something a working researcher would use?

Compute after Gates 1 and 2 pass. Pass threshold is 0.70 (vs. 0.75 in the v1 rubric) because SR is more permissive on perceived importance and the audit already softened novelty claims. If Gate 3 fails, iterate on the lowest-scoring axis.

---

# Integrity checklist — runs every iteration

For each item: **pass**, **fail**, or **N/A with reason**. A failed item must be fixed in the same iteration or the iteration is not complete.

- [ ] **No peeking.** I did not look at validation results (V1 recall, V6 contamination, V7 accuracy) before locking the relevant methodology decisions for this iteration.
- [ ] **No silent multiple testing.** If I tried multiple parameter settings or method variants in this iteration, I report all of them, not only the best.
- [ ] **No cherry-picked baselines.** MetaVolcanoR, hStouffer, AWmeta and the classical baselines are what a knowledgeable reader would expect. None were replaced with weaker alternatives.
- [ ] **No unjustified study exclusion.** Any study dropped after initial catalog entry has a documented categorical reason.
- [ ] **No asymmetric tuning.** PIPER hyperparameters were not tuned with more effort than baselines. Tuning logs match across methods.
- [ ] **No leakage in V7.** Held-out evaluation set was never used during prompt iteration or feature engineering.
- [ ] **No leakage in V2.** The pipeline-perturbation re-analyses used matched filtering and design; differences in consensus signatures across pipeline collections reflect pipeline differences, not filtering differences.
- [ ] **No V6 cherry-picking.** Case study was pre-specified; if unworkable, next pre-specified candidate was used with documentation.
- [ ] **No scope creep in claims.** What I am about to log as the result is what was actually tested, not a broader framing.
- [ ] **No fabricated citations or data.** Nothing in notes or writing is invented.
- [ ] **Gold-standard target list unchanged.** No additions or removals to the HIF1α panel after running any analysis.
- [ ] **CI calibration honesty.** If reporting credible intervals for a gene, I either ran V4-style coverage check or explicitly note that calibration was not verified.

When you feel friction against one of these — "but this one time it's fine because..." — that is when the checklist is doing its job. Log the tension in the notebook and choose the cleaner path. Scientific Reports may not gate hard on novelty, but it gates very hard on these.
