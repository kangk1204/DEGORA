# Project Context — DEGORA Methodology Development

## Research question

Develop and validate a transcriptomic meta-analysis method that operates directly on as-published heterogeneous DEG tables — supplementary tables from published RNA-seq papers, produced by different analysis pipelines (DESeq2, edgeR, limma-voom, Cuffdiff, NOISeq, sleuth, ballgown, EBSeq, dearseq, …) with different LFC definitions, different shrinkage behavior, and frequently without standard errors. The method (working name: **DEGORA**) auto-detects each table's pipeline, applies bias-aware harmonization, and yields a triplet-concordance consensus signature (magnitude + rank + sign) with per-gene provenance.

The full technical specification is in `METHODOLOGY_SPEC.md`. The prior-art positioning is in `NOVELTY_AUDIT.md`. This file captures only the higher-level project constitution that frames how the loop operates.

## Input mode

Concrete proposal + curated test bed. Methodology paper at *Scientific Reports* tier. Test bed: hypoxia response RNA-seq.

## Domain

Bioinformatics methodology — transcriptomic meta-analysis. Statistical / engineering hybrid contribution. Reference implementation in Python with R interop via rpy2 for established baselines.

## Data

### Test bed corpora

- **Hypoxia core corpus** — 20–30 published hypoxia RNA-seq studies with full DEG tables in supplementary. Hand-curated. Catalog at `data/studies/hypoxia_catalog.csv`. Pipeline diversity required (mixture of DESeq2 / edgeR / limma-voom / older).
- **Hypoxia raw-counts subset for V2** — 8–12 of the above where raw count matrices are publicly accessible on GEO, recount3, or ARCHS4. Used for the pipeline-perturbation experiment (re-analysis with DESeq2 + edgeR + limma-voom under matched protocol).
- **TNF-α cross-domain corpus for V5** — 5–10 published TNF-α stimulation RNA-seq studies, same curation standards. NF-κB target panel as gold standard.
- **Labeled DEG-table benchmark for V7** — ≥ 200 DEG tables with known pipeline labels, used to train + evaluate the `degpipeline` classifier. Collected during S0; 70/30 train/test split with documented anti-leakage protocol.
- **V6 case study** — one previously-published transcriptomic meta-analysis on a clear topic whose contributing DEG tables remain accessible. Pre-specified at the start of S7.

### Supplementary databases

- **MSigDB** (Hallmark, Reactome, KEGG, WikiPathways) — for pathway-level integration in Component D and `HALLMARK_HYPOXIA` validation
- **HIF1α target panel** (locked, `METHODOLOGY_SPEC.md` §3.2) — for V1 positive-control recovery

### Known issues / pitfalls

- Pipeline annotation in some older papers is ambiguous ("DEGs were identified using established Bioconductor tools"). Mark these as `unknown_pipeline` and proceed; do not infer.
- Gene ID systems vary (Ensembl with/without version, Entrez, HGNC symbol, sometimes legacy aliases). All studies must be standardized to Ensembl-stable + HGNC-symbol dual labeling before integration.
- Some DEG tables include only "significant" genes, not the full tested set. These violate the inclusion criterion and are excluded — they are NOT to be imputed.
- A few hypoxia papers report results per cell line or per time point in separate sheets. Treat each contrast as a separate `study_id` with shared `paper_id`.
- Cuffdiff tables can have `±inf` LFC values. Cap at ±10 with `flag_capped=True` recorded in the harmonized table.

## Prior work flagged as relevant

(Full audit in `NOVELTY_AUDIT.md`. Summary here.)

- **MetaVolcanoR** (Prada 2019, Bioconductor) — closest competitor; DEGORA differentiates via pipeline-awareness and triplet concordance
- **hStouffer** (Kim et al. 2026, BMC Bioinformatics; first author D. Kim, senior author J-Y. Lee, KNU) — recent competitor in uniform-pipeline regime; DEGORA addresses rank-integration gap they flagged
- **AWmeta** (Yang et al. 2025, bioRxiv) — adaptive weighting state of the art; DEGORA's weighting is feature-based, not optimization-based
- **crossmeta** (Pickering 2018) — inverse-philosophy comparison (re-process from raw)
- **Robust Rank Aggregation** (Kolde 2012) — component of DEGORA + baseline
- **DerSimonian–Laird REM, REML** — gold-standard baseline
- **metaRNASeq, MetaIntegrator, DExMA, BayesMP, AW-Fisher** — adjacent context

## Hard constraints

- **Test bed remains hypoxia + TNF-α as cross-domain demo.** Do not silently expand to other domains; do not narrow hypoxia to a sub-condition.
- **DEG tables must be as-published.** No silent re-processing from raw counts as a workaround when a paper's table is messy. Either the table is usable or the study is excluded.
- **No fabricated citations.** Every reference in the manuscript verified via `ref-finder`. No exceptions.
- **No fabricated DEG tables.** If a study is excluded, log it with reason in the catalog and the lab notebook. Do not synthesize a "plausible" table.
- **Baseline implementations must use published versions.** MetaVolcanoR from Bioconductor as-published, hStouffer from `CR4CKID/hStouffer` GitHub as-published, AWmeta re-implemented from the bioRxiv preprint with the implementation choices documented. No "improved" baselines that happen to make DEGORA look better.
- **Open-source code.** MIT or Apache-2.0. `degpipeline` separate-releasable.
- **Open data.** Curated catalog + standardized DEG parquets on Zenodo with DOI before submission.
- **Manuscript and code in English.** Internal logs in English for portability.

## Budget ceiling

- **Wall-clock for autonomous agent runtime:** 150 hours (≈ 6 days of agent time). Methodology paper scope is narrower than the DB project's; this should be sufficient.
- **Compute:** CPU-bound work on a Linux workstation. H200 GPU optional and only useful for `degpipeline` XGBoost training (overkill — CPU suffices) and any embedding experiments (none planned in Phase 1).
- **Storage:** ~50 GB free should suffice (no ARCHS4 H5 needed; raw counts for V2 subset are smaller).
- **External API budget:** ≤ 100 USD if any LLM API calls are made during catalog curation. If using local models on H200, GPU time only.
- **Max iterations before pivot trigger:** 5 consecutive iterations without resolving any material reviewer-flagged issue → propose pivot or run final gates.

## Target venue

**Primary: *Scientific Reports* (Nature Portfolio).**

- Open access, broad scope, technical-soundness review (not subjective importance gating)
- Acceptance criteria favor: rigorous methods, reproducibility, honest reporting
- Acceptance criteria do NOT require: breakthrough novelty claims, broad significance pitching
- Word limit: not strict, but aim for 4,500–8,000 words main text
- Display items: up to ~8 main figures/tables; supplementary unlimited
- Manuscript structure: Introduction → Results → Discussion → Methods (Methods last is SR convention, unusual)
- Review focus: validity of conclusions, methodological rigor, code/data availability

**Why this venue fits DEGORA:**

- The audit established DEGORA's novelty is PARTIAL — sharp on pipeline-aware harmonization and the pipeline-perturbation benchmark, incremental elsewhere. Scientific Reports does not gate on novelty intensity, so this is honest framing rather than uphill positioning.
- Methodology + benchmark + reference implementation fits SR's profile well.
- SR's reproducibility expectations align with what DEGORA must deliver anyway.

**Backup tiers if SR rejects on technical grounds rather than scope:** *PLOS ONE*, *F1000Research* (open peer review model), *Briefings in Bioinformatics* if a method-paper angle survives.

**Do NOT target:** BMC Bioinformatics (hStouffer territory overlap, recently published there), or any venue where the user has had recent territory conflicts.

## Domain-specific pitfalls to watch

Distilled from the bioinformatics methodology landscape:

- **Batch effects across studies are the expected heterogeneity here, not a bug.** DEGORA's job is to model them via I², not to "correct" them. Reporting I² openly is required.
- **Pipeline detection is heuristic by definition.** A reviewer will challenge edge cases (DESeq2 with `apeglm` vs. `ashr`; renamed columns; partial outputs). The loop must report classifier accuracy with confidence calibration, not just headline accuracy.
- **Multiple testing.** Each per-study DEG already has padj. DEGORA adds another correction layer; BH-FDR across all tested genes per consensus output. Document.
- **Effect sizes alongside p-values.** Every consensus entry has pooled LFC + 95% CI, not just significance.
- **Cohort composition reported.** Per-corpus table of contributing N studies, species breakdown, cell-system breakdown, pipeline distribution.
- **Pathway redundancy.** Deduplicate top-K Hallmark/Reactome/KEGG by Jaccard ≥ 0.7 before reporting in the manuscript.
- **CI calibration.** Reporting credible intervals (Component E) without verifying coverage probability (V4) is a methods-paper red flag. V4 is where this gets honest.

## Addenda from the loop

(The loop may append assumptions it inherited below this line. Original user content above this line is never edited.)
