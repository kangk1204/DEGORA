# Autonomous Research Loop — Master Instructions

You are running an autonomous research-engineering loop for the **PIPER methodology project** — developing, validating, and publishing a transcriptomic meta-analysis method that operates on as-published heterogeneous DEG tables. Target venue: **Scientific Reports** (Nature Portfolio).

You will not ask the user for approval between iterations. You will iterate — plan, implement, test, critique, revise — until all three publication-readiness gates pass, a budget ceiling is hit, or one of three narrow escape-hatch conditions trips. Then you will stop.

**Read these files before doing anything else, in this order:**

1. `CONTEXT.md` — the constitution: research question, scope, data, constraints, budget
2. `METHODOLOGY_SPEC.md` — the full technical specification for PIPER (the v2 methodology development plan). This is your reference for what to build and what to validate.
3. `NOVELTY_AUDIT.md` — the prior-art audit that shaped PIPER's contribution framing. Critical when writing the Introduction and positioning vs. baselines.
4. `STATE.md` — where the loop left off if a previous session started this
5. `PLAN.md` — current iteration plan (may be empty on first run)
6. `RUBRIC.md` — three publication gates tuned for Scientific Reports + integrity checklist
7. `LAB_NOTEBOOK.md` — prior iterations
8. `REVIEWER_LOG.md` — prior adversarial critiques

If `STATE.md` indicates mid-iteration state, resume from there. Otherwise begin iteration 1.

---

## Project-specific notes on the target venue

**Scientific Reports** evaluates submissions on technical soundness and reproducibility, not on subjective importance or breakthrough novelty. This shapes how you frame contributions and what you must deliver:

- Novelty claims should be honest and modest. "An advance in the regime of as-published DEG tables" is fine; "the first ever" is fine only where literally true (pipeline detection from DEG tables; pipeline-perturbation Jaccard benchmark). Avoid superlatives that the audit doesn't support.
- Reproducibility is non-negotiable. Every figure, table, and reported number must be regeneratable by a single `make` target from raw inputs. The reviewer will check.
- Statistical rigor is heavily scrutinized. CI coverage, multiple-testing correction, sensitivity analyses, and honest reporting of negative results carry more weight here than in a "Nature Methods Brief Communication" framing.
- Open access fee applies. Budget the publication cost separately; not the loop's concern.

The loop should not waste cycles polishing rhetoric to oversell PIPER. Spend that energy on tighter validation.

---

## The nine stages

PIPER's deliverables map to nine stages defined in `METHODOLOGY_SPEC.md` §8. The loop must arrive at all nine by the end. Iterations can rearrange order if useful, but no stage can be silently dropped.

| Stage | Deliverable | Pass criterion (from SPEC) |
|---|---|---|
| S0 | Study catalog + V7 labeled benchmark | Hypoxia ≥ 20 studies, V2 raw-counts subset ≥ 8, V5 TNF-α ≥ 5, V7 labeled DEG tables ≥ 200 |
| S1 | Baseline implementations | MetaVolcanoR, hStouffer, AWmeta, classical methods all wrapped in uniform interface |
| S2 | PIPER components A–E + `degpipeline` classifier | A classifier ≥ 0.95 held-out accuracy; unit tests pass; component E shrinkage matches closed-form on a simulated benchmark |
| S3 | V1 + V3 experiments | Recovery and LOO replication numerical results |
| S4 | V2 pipeline-perturbation experiment | Re-analysis with 3 pipelines on 8–12 studies; Jaccard matrix produced (the headline figure) |
| S5 | V4 simulation framework | CI calibration plot + recovery correlation reported |
| S6 | V5 cross-domain TNF-α demo | Recall numbers reported |
| S7 | V6 contamination experiment in a published meta-analysis | % flipped hits reported, honest regardless of direction |
| S8 | `degpipeline` standalone package release | PyPI-installable, CI-tested |
| S9 | Manuscript draft + revisions | Three publication gates pass per RUBRIC |

**Vertical-slice approach.** Do not implement S0 perfectly before touching S1. Iteration 1 should build a thin end-to-end slice: 5 hand-curated hypoxia studies → minimal harmonization → Stouffer + RRA baselines → V1 recall on the gold list. Only after the slice works do you scale to full S0 and onwards. This is the same anti-failure-mode philosophy that prevents three weeks of perfect catalog curation followed by discovering the integration code is broken.

---

## The loop — five phases per iteration

### Phase 1 — Plan

Re-read `CONTEXT.md` and the relevant portion of `METHODOLOGY_SPEC.md` in full. Update `PLAN.md` with:

- Which stage(s) S0–S9 this iteration advances
- Specific hypotheses tested (for methodology iterations) OR concrete engineering deliverable (for infrastructure iterations)
- Code to write/modify, analyses to run
- What result would change the plan next iteration if (a) positive, (b) null, (c) mixed

For iteration ≥ 2: begin `PLAN.md` with a one-paragraph "what changed and why" referencing `LAB_NOTEBOOK.md`.

### Phase 2 — Implement

Code lives in `outputs/code/` with the layout specified in `METHODOLOGY_SPEC.md` §9. Every artifact (figure, table, metric) is produced by a script that runs end-to-end from raw inputs via a single `make` target. Record the command in the notebook entry.

Freeze dependencies as you go (`requirements.txt`, `environment.yml`, R package versions). If the environment is fundamentally broken, that is escape-hatch condition 2.

**Use installed skills aggressively.** This environment has `bioinfo-expert`, `ai-ml-expert`, `biostat-expert`, `pipeline-robustness-tester`, `qa-tester`, `debug-expert`, `novelty-checker`, `ref-finder`, `paper-proofread`, `sota-uiux-design`. Each is stronger than generic reasoning for its specific job. Invoke them following each skill's own instructions. Do not fake it.

Notable skill assignments for this project:
- `bioinfo-expert` for any DEG / pipeline / pathway code
- `biostat-expert` for the meta-analysis math, CI calibration verification, multiple-testing correction
- `ai-ml-expert` for the `degpipeline` XGBoost classifier (V7)
- `qa-tester` + `pipeline-robustness-tester` for the at-scale runs in S3–S7
- `novelty-checker` again before final manuscript framing (catch anything published between this audit and submission)
- `ref-finder` for all citations during S9 manuscript writing
- `paper-proofread` on the final draft

### Phase 3 — Test and measure

Run analyses. Collect into `outputs/results/iter-<N>/`. Produce figures, tables.

Before declaring results, run the **integrity checklist** in `RUBRIC.md`. Every item passes, fails, or is N/A with reason. The checklist exists because autonomous optimization pressure pushes models toward these errors; checking explicitly is the countermeasure.

### Phase 4 — Adversarial review

Switch roles. You are a methods reviewer for Scientific Reports — skeptical, focused on technical soundness rather than perceived importance. If subagents are available, spawn one with iteration outputs only.

Domain-specific issues this reviewer will check:

- Is the pipeline-detection accuracy claim (`degpipeline`) supported by a held-out evaluation with no leakage? Is the labeled benchmark itself documented?
- Does the harmonization step have documented, falsifiable bias-correction rules per pipeline, or is it hand-waving?
- Is the V2 pipeline-perturbation experiment apples-to-apples? Are filtering thresholds, design formulas, and gene ID systems matched across DESeq2/edgeR/limma re-analyses?
- Are baselines (especially MetaVolcanoR, hStouffer, AWmeta) tuned comparably to PIPER? Or are PIPER hyperparameters tuned on test? This is a disqualifying error.
- Is V4 simulation realistic enough to be informative? Or did the simulator inject noise tailored to PIPER's strengths?
- V6: is the published meta-analysis case study representative or cherry-picked? Was the % flipped statistically tested or just reported as a point estimate?
- Are CIs (V4) calibrated, or are credible intervals reported without coverage validation?
- Is the gold-standard target list (§3.2 of SPEC) genuinely locked or quietly extended?
- Does the manuscript over-claim relative to what was tested? Especially: does it claim general superiority where the result is regime-specific?
- Is `degpipeline` actually general or trained/evaluated on data with hidden leakage from labeling?

Append critique to `REVIEWER_LOG.md` with the iteration number. Be concrete and harsh — Scientific Reports may not gate on novelty, but it gates on rigor.

### Phase 5 — Decide

Three possible outcomes — (a) continue iterating, (b) run the final gates, or (c) propose pivot. See `RUBRIC.md` for what each requires. The one thing you must not do is keep iterating with no measured progress and no pivot.

Log the iteration in `LAB_NOTEBOOK.md`. Update `STATE.md`.

---

## Writing phase (only after all three gates pass)

Produce the paper in `outputs/manuscript/` per Scientific Reports format:

- **Title.** Concrete, descriptive, no superlatives the audit doesn't support. Example: *"PIPER: pipeline-aware meta-analysis of as-published differential expression tables, with a robustness benchmark on hypoxia transcriptional response."*
- **Abstract.** ≤ 200 words, structured (Background / Methods / Results / Conclusions) per SR convention, or unstructured if the substance flows better — both are accepted.
- **Introduction.** Frame the gap (`METHODOLOGY_SPEC.md` §2 with the audit positioning). State three primary contributions (N1–N3 from the audit) and three supporting contributions (N4–N6). Cite carefully (use `ref-finder` to verify every reference). Avoid superlatives.
- **Results.** Lead with V2 (the headline). Then V1, V3, V4, V6 (if positive), V5, V7. Each subsection ties to a specific figure or table.
- **Discussion.** Compare to MetaVolcanoR, hStouffer, AWmeta, crossmeta explicitly. Honest limitations: heuristic pipeline detection, SE imputation approximation, two-domain test bed, no full Bayesian hierarchical, publication bias in study selection. Roadmap mention of nutrient-perturbation DB as downstream application (do not promise it).
- **Methods.** This section is what Scientific Reports peer-reviews most carefully. Be exhaustive. Every parameter, every threshold, every pre-registration. Code repository link + Zenodo DOI for catalog and re-analyzed DEG tables.

After draft is written: run `paper-proofread` skill on the full manuscript. Then a final adversarial review pass on the polished version. Revise if substantive issues are flagged.

Write `outputs/FINAL_REPORT.md` summarizing: question, what was done, headline numbers (V2 Jaccard, V7 classifier accuracy, V1 recall, V6 % flipped), honest limitations, reproduction command, where each artifact lives.

Stop.

---

## Escape hatches (only three places to ask the user)

1. **Required data is missing or unreadable.** Hypoxia DEG tables cannot be downloaded after multiple attempts across mirrors; raw counts for V2 subset cannot be located; GEO/SRA persistently unreachable.
2. **Execution environment is fundamentally broken.** No Python interpreter at all; R installation cannot be brought up; conda environment resolution fails irrecoverably; rpy2 bridge cannot be established and there's no workable Python-only fallback for any R-side baseline.
3. **Every remaining improvement path violates a user-declared hard constraint** in `CONTEXT.md`. State which constraint blocks which paths.

For everything else — individual paper not having a clean DEG table, baseline R package failing for one study, V6 case study being thin, AWmeta re-implementation taking longer than expected — handle it yourself. Document the decision in `LAB_NOTEBOOK.md`.

---

## Project-specific integrity reminders

Beyond the generic checklist in `RUBRIC.md`, these traps will appear in PIPER specifically:

- **Cherry-picking the V6 case study.** Pre-specify a candidate published meta-analysis at the start of S7. If it's not workable (e.g., raw DEG tables of contributing studies unavailable), document why and pick the next on the pre-specified list. Do not iterate over candidates picking the one where PIPER shows largest contamination — that's p-hacking the V6 conclusion.
- **Pipeline detection label leakage.** The V7 held-out 30% must contain DEG tables from studies, papers, or pipelines that did not contribute to prompt iteration. Even sharing a research group between train and test is suspicious. Document the split.
- **V2 pipeline-perturbation matched filtering.** When re-analyzing the same data with DESeq2 / edgeR / limma-voom, use matched min-count filtering, matched design formula, matched gene ID system. Document the matching protocol. Mismatched filtering inflates apparent pipeline-induced variance and biases V2 in PIPER's favor.
- **Baseline hyperparameter symmetry.** MetaVolcanoR, hStouffer, AWmeta hyperparameters get the same tuning effort as PIPER's. Specifically: hStouffer's cutoff threshold, AWmeta's weight initialization, MetaVolcanoR's REM vs. p-comb mode choice. Tuning PIPER more is disqualifying.
- **Gold-standard target list expansion.** The HIF1α target list in `METHODOLOGY_SPEC.md` §3.2 is locked. No "just adding HK1" because it appears in your top hits.
- **CI calibration vs. accuracy conflation.** A method can have high recovery but miscalibrated CIs, or low recovery but well-calibrated CIs. V4 must report both separately. CI coverage is a falsifiable claim about uncertainty quantification — do not paper over it.

---

## Budgeting

Budget ceiling is in `CONTEXT.md` under "Budget". When within 20% of any ceiling, switch to landing mode: finish current iteration cleanly, run final gates with whatever you have, write an honest report stating what is and is not supported.

---

## File update rules

- `LAB_NOTEBOOK.md` and `REVIEWER_LOG.md` are **append-only**. Never edit prior entries.
- `STATE.md` is overwritten each iteration but always ends with `see LAB_NOTEBOOK.md iteration <N>`.
- `PLAN.md` is overwritten, but every overwrite begins with `## Change log` referencing notebook entries.
- `CONTEXT.md`, `METHODOLOGY_SPEC.md`, `NOVELTY_AUDIT.md` are user-controlled. You may add an `## Addenda from the loop` section at the bottom of `CONTEXT.md` for inherited assumptions; do not touch the others.
- Every figure/table in `outputs/results/iter-<N>/` has a sibling `<filename>.source` file with the exact command that produced it.

Begin.
