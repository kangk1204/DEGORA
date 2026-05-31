# DEGORA Novelty Audit & Contribution Upgrade

> **STATUS (2026-05-31): partly SUPERSEDED — read with the resource framing.**
> The submitted manuscript (`outputs/manuscript/main.md`) was pivoted to a *source-auditable evidence resource/database* contribution. The "genuinely novel" method components claimed below — **N1 (pipeline detection / bias-aware harmonization)** and **N2 (pipeline-perturbation Jaccard benchmark)** — are **NOT implemented** in the codebase and are NOT claimed by the manuscript (pipeline labels come from the user catalog; `harmonize.py` is minimal column-mapping; the only Jaccard is a source-unit ablation). A fresh prior-art search (2026-05-31) reconfirmed PARTIAL novelty and identified **MAIC (Parkinson et al. 2020, Sci Rep, PMID 33339864) as the strongest direct prior art (7/8)** for the quality-weighted repeated-support ranking — it must be cited and differentiated (predeclared vs information-content weights; direction/heterogeneity axis; offline self-contained archive). The defensible contribution is the *integrated as-published-DEG-table resource + offline Evidence Explorer*, not a new method. See `PHASE_A_STATUS.md`.

**Date of search:** 2026-05-25
**Sources covered:** PubMed (E-utilities), Google Scholar (via web_search), bioRxiv, BMC/Springer, journals (BMC Bioinformatics, Bioinformatics, Nature Methods family), GitHub project pages.

---

## Novelty Verdict: **PARTIAL**

DEGORA as originally drafted has five proposed components (A pipeline-aware harmonization, B multi-evidence triplet, C adaptive quality weighting, D hierarchical gene-pathway integration, E EB shrinkage). The audit finds:

- **Components A and the pipeline-perturbation Jaccard benchmark (V2 experiment) are genuinely novel.** No published method auto-detects pipeline from a DEG table or measures meta-analysis robustness to pipeline perturbation.
- **Components B, C, D, E individually are NOT novel enough to load-bear the paper.** Each has a published precedent (MetaVolcanoR, AWmeta, hStouffer, Wei et al., metafor). Their combination is incremental.
- **The positioning "operate on as-published heterogeneous DEG tables, not re-processed counts" is novel and underutilized in the original draft.** This is the actual gap.

The paper survives, but needs a sharper contribution framing. See §3 for the rewrite.

---

## 1. Direct prior art (overlap score ≥ 6 / 8)

Scoring rubric: Method (0–2) + Data (0–2) + Task (0–2) + Claim (0–2). ≥6 = direct prior art. 3–5 = adjacent. ≤2 = background.

### 1.1 MetaVolcanoR (Prada-Medina et al., *Bioconductor* 2019)
- DOI: 10.18129/B9.bioc.MetaVolcanoR
- **What it does:** R package combining DEG results across studies with three strategies: (i) Random Effects Model on LFC+SE, (ii) p-value combining, (iii) vote-counting. Volcano-plot visualization of consensus.
- **Overlap with DEGORA:** Method=1 (same family, different specifics), Data=2 (per-study DEG tables), Task=2 (cross-study consensus DEG), Claim=2 (multi-method consensus). **Score: 7/8**
- **What it doesn't do:** No pipeline detection, no pipeline-aware harmonization, no triplet concordance requirement (it gives three separate strategies as parallel results, not a concordance call), no per-study quality weighting beyond default, no pipeline-perturbation benchmark.
- **DEGORA must:** Cite explicitly, benchmark against, and articulate that DEGORA's contribution is the *pipeline-aware* layer that sits *upstream* of methods like MetaVolcanoR, not a replacement of them.

### 1.2 AWmeta (Yang et al., *bioRxiv* 2025.05.06.650408, currently v5)
- bioRxiv DOI: 10.1101/2025.05.06.650408
- **What it does:** Adaptive-weight transcriptomic meta-analysis. Two modules: AW-Fisher (adaptive p-value aggregation with per-study weights minimizing combined p), AW-REM (weighted REM for effect size).
- **Overlap with DEGORA:** Method=1, Data=2, Task=2, Claim=1 (adaptive weighting overlaps with Component C but mechanism is different). **Score: 6/8**
- **What it doesn't do:** No pipeline awareness, no pipeline-perturbation benchmark, weights are derived from p-value optimization not from study quality features.
- **DEGORA must:** Drop "adaptive weighting" as a standalone novelty claim. Reframe Component C as "auditable per-study quality features driving weights" (different from AWmeta's optimization-based weights) — or fold it into Component A as part of pipeline-aware quality assessment.

### 1.3 hStouffer (Kim et al., *BMC Bioinformatics* 2026, Feb 12; first author D. Kim, senior author J-Y. Lee, KNU)
- DOI: 10.1186/s12859-026-06395-2 — published 3 months before this audit
- **Authors based at:** Kyungpook National University (KNU), Korea — same country as the user; potential collaboration angle rather than purely competitive
- **What it does:** Three-step framework on top of Stouffer's method: (i) p-value capping at 10⁻³ to limit extreme outliers, (ii) dataset-size-dependent dynamic cutoff threshold, (iii) bagging (1000 bootstrap iterations, gene must be DEG in >95%). Plus REM filter for directional consistency.
- **Crucial detail:** They re-analyze 137–500+ ARCHS4 datasets *uniformly with DESeq2 v1.34.0* — so they avoid pipeline heterogeneity by design. They are working in the "uniform pipeline, many studies" regime, not the "heterogeneous published DEG tables" regime.
- **Overlap with DEGORA:** Method=1, Data=2, Task=2, Claim=1. **Score: 6/8**
- **What it doesn't do:** No pipeline-awareness (single-pipeline design), explicitly notes in Limitations: *"the current framework does not incorporate rank-based information, which could provide complementary strengths to *p* value and effect-size–based integration. Incorporating rank-based measures represents a promising direction for future refinement"*. This is a gift — DEGORA's multi-evidence triplet (which includes rank via RRA) directly addresses what hStouffer's authors flag as future work.
- **DEGORA must:** Cite, benchmark against, and explicitly note that DEGORA addresses the rank-information gap hStouffer's authors identified. Position DEGORA as complementary in regime: hStouffer for "uniform re-processed counts, large N", DEGORA for "as-published heterogeneous DEG tables".

---

## 2. Adjacent work (overlap score 3 – 5 / 8)

These are competitive contextually but don't overlap enough to block DEGORA's contribution. Each must be cited and briefly differentiated.

| Tool / Paper | Year | What it does | Overlap | Differentiation in DEGORA's introduction |
|---|---|---|---|---|
| **metaRNASeq** (Rau, Marot, Jaffrézic) | 2014 | Fisher + inverse normal (Stouffer) p-value combination, R package | 5/8 | Foundational; DEGORA addresses pipeline heterogeneity which metaRNASeq does not |
| **MetaIntegrator** (Haynes et al.) | 2017 | DerSimonian-Laird REM per gene + Fisher; requires gene-level effect size + variance | 5/8 | Assumes consistent effect size definition; DEGORA explicitly handles pipeline-specific LFC bias |
| **crossmeta** (Pickering) | 2018 | Re-processes raw data from GEO with consistent pipeline, then meta-analyzes | 4/8 | INVERSE approach (avoids heterogeneity by re-processing); DEGORA accepts published heterogeneity as input |
| **DExMA** (Toro-Domínguez et al.) | 2022 | Handles missing genes across studies via imputation | 3/8 | Different aspect of meta-analysis; complementary |
| **BayesMP / Huo–Song–Tseng** | 2019 | Bayesian latent hierarchical, detects DE in subsets of studies, tight-clustering of meta-patterns | 5/8 | Heterogeneity-aware in a different sense (which studies show signal, not which pipelines); roadmap relation to DEGORA's Phase 2 |
| **AW-Fisher** (Li & Tseng) | 2011 | Foundational adaptive weighting for Fisher's method | 4/8 | Foundational for AWmeta; DEGORA's weighting is feature-based not optimization-based |
| **Wei et al. — Sample-quality weighted RE models** | 2019 | Sample-quality weights in DSL and Bayesian RE models for gene expression | 4/8 | Closest precedent for DEGORA's quality weighting; need to cite and differentiate (per-study vs per-sample) |
| **jGRP** (regulation probability) | 2017 | Transforms expression to "gene regulation profile space" before meta-analysis | 3/8 | Different harmonization mechanism, predates DEGORA's pipeline-detection idea |
| **RankProd / RobustRankAggreg** | 2004, 2012 | Rank-based meta-analysis; foundational for DEGORA's rank axis | 4/8 | Component of DEGORA, not competitor |
| **FedPyDESeq2** | 2024 bioRxiv | Federated DESeq2 across centers; shows federated > meta-analysis in power | 3/8 | Requires raw data access at each site; not applicable to as-published GEO data |
| **arXiv:2601.04122** (Tool Choice Matters) | 2026 May | edgeR vs DESeq2 cross-study generalization benchmark | 3/8 | Motivates DEGORA's problem; not a meta-analysis method itself |
| **"Tool combinations in pipeline versions"** (Merino et al.) *NAR Genom Bioinform* | 2024 | How pipeline tool combinations affect DEG outcomes | 3/8 | Motivates the problem DEGORA solves |

---

## 3. What is genuinely new in DEGORA (revised after audit)

The original draft had five components claimed as novel. After audit:

### 3.1 Genuinely novel (load-bear the paper)

**N1. Pipeline detection and harmonization of as-published DEG tables.**
No published tool auto-detects the analysis pipeline from DEG table structure (column fingerprint + distributional fingerprint) and applies pipeline-specific calibration. crossmeta avoids the problem by re-processing; hStouffer/AWmeta/MetaVolcanoR/MetaIntegrator all assume consistent or convertible effect-size definition across studies. DEGORA explicitly does not assume this.

**N2. Pipeline-perturbation Jaccard benchmark methodology.**
No paper has systematically measured meta-analysis method robustness by applying the method to (a) actual heterogeneous DEG tables and (b) same data re-analyzed with one uniform pipeline, and comparing top-K Jaccard. This is the headline experiment AND a generalizable benchmark methodology that future meta-analysis methods can be evaluated against. It is potentially a separate publishable contribution.

**N3. As-published-DEG-tables operating regime.**
All existing meta-analysis tools require either raw counts (crossmeta, hStouffer, FedPyDESeq2) or DEG tables in a specific harmonized format with consistent effect size definition (MetaVolcanoR, MetaIntegrator, AWmeta). DEGORA operates on whatever a supplementary table provides — DESeq2, edgeR, limma-voom, Cuffdiff, NOISeq, mixed — without re-processing. This is a different problem statement that hasn't been framed as such in the literature.

### 3.2 Defensibly novel as part of the package (but not load-bear alone)

**N4. Triplet concordance call (require all three axes to agree).**
MetaVolcanoR provides three strategies as parallel outputs ("here are three rankings, pick one"). AWmeta provides two modules in sequence (p-value then REM). DEGORA's contribution is the explicit *concordance call*: a gene is high-confidence only if magnitude, rank, and sign axes all agree. This is a subtle but defensible differentiation, and addresses the "future direction" hStouffer's authors explicitly flagged.

**N5. Per-gene provenance graph in the output.**
None of the existing tools output per-gene provenance (which study contributed, which pipeline, which sample size, which lab). This is auditable in a way that other meta-analysis outputs are not, and turns DEGORA's output into something a reviewer can interrogate at the gene level. Easy to implement; high perceived value.

### 3.3 NOT novel — drop or de-emphasize

- ❌ **Adaptive weighting as a standalone contribution.** AWmeta (May 2025) has this. Keep DEGORA's quality-feature-based weighting as a *practical refinement* mentioned in Methods, not a contribution.
- ❌ **Empirical Bayes shrinkage as a contribution.** Well-established; limma's eBayes, DESeq2's apeglm, ashr all use EB. Keep in DEGORA's pipeline as a regularizer; don't claim it as novel.
- ❌ **REM + RRA + Stouffer combination as the headline.** MetaVolcanoR + hStouffer's REM filter cover too much of this. Demote from headline to "components of the triplet axis".
- ⚠️ **Hierarchical gene → pathway integration.** Multiple existing methods do this (MAPE_P, MAPE_G, gene-set meta-analysis literature). Keep as a methodological refinement but not a primary contribution claim.

---

## 4. Five concrete novelty upgrades

These are pivots that, applied to the methodology .md, sharpen the contribution package so it survives a methods reviewer.

### Upgrade 1 — Reposition as the "as-published DEG tables" method

Rewrite the abstract and Introduction to lead with this gap. Concrete language:

> *"Existing transcriptomic meta-analysis methods assume either uniform re-processing from raw data (crossmeta, hStouffer, FedPyDESeq2) or DEG tables with consistent effect-size definitions (MetaVolcanoR, MetaIntegrator, AWmeta). In practice, the supplementary DEG tables a researcher finds in published papers are heterogeneous — produced by DESeq2, edgeR, limma-voom, Cuffdiff, or NOISeq, with different LFC definitions, different shrinkage behavior, and frequently without standard errors. We present DEGORA, the first meta-analysis framework that operates directly on as-published heterogeneous DEG tables, with explicit pipeline detection and bias-aware harmonization."*

This framing makes DEGORA's existence necessary in a way that "yet another adaptive Stouffer variant" does not.

### Upgrade 2 — Promote the pipeline-perturbation benchmark to a standalone methodological contribution

The pipeline-perturbation Jaccard experiment is the most defensible piece of novelty in DEGORA. It deserves more than being one figure in the Results.

- Frame it as a **benchmark methodology** that all future meta-analysis methods should report
- Pre-register the protocol: same studies, three pipelines (DESeq2, edgeR, limma-voom), the four DEG collections setup, Jaccard@50, Jaccard@100, Spearman of signed-z
- Provide a public benchmark script in the DEGORA repo that takes any meta-analysis function and runs it through the pipeline-perturbation evaluation
- This is potentially a separate short paper (*Bioinformatics* Application Note style) if the user wants to split it out

### Upgrade 3 — Add pipeline classifier as a standalone tool

The XGBoost pipeline detector in Component A doesn't need to be buried inside DEGORA. Release as:

- A standalone Python package `degpipeline` with a single function `detect_pipeline(deg_table) → {pipeline: "DESeq2", confidence: 0.97, version_hint: "≥1.20"}`
- Useful to the broader community independent of DEGORA (e.g., GEO curators, paper reviewers, automated meta-analysis pipelines)
- Add this as a fifth concrete contribution claim
- Gives the paper an ML/AI angle that pure statistical meta-analysis papers lack — slight but defensible

### Upgrade 4 — Quantify pipeline-induced contamination in published consensus signatures

This is the experiment that turns DEGORA from "another meta-analysis tool" into "a tool that revealed a real problem in the literature":

- Take a public consensus signature from a high-profile meta-analysis paper (e.g., the MetaVolcanoR psoriasis example, an AWmeta case, or any meta-analysis paper that used heterogeneous published DEG tables)
- Re-run with DEGORA's pipeline-aware harmonization
- Quantify: what fraction of the original consensus hits change when pipeline heterogeneity is corrected for?
- If the answer is "a meaningful fraction" — this is the kind of finding that makes a methods paper feel urgent rather than incremental

Risk: if the answer is "almost no change" then this experiment cuts against DEGORA. Run it early in iteration 2 or 3 and decide whether to include before final framing.

### Upgrade 5 — Add per-gene provenance graph as a deliverable

For every consensus gene, output a structured provenance record:

```json
{
  "gene": "VEGFA",
  "consensus": {"lfc": 2.31, "ci": [1.8, 2.8], "i2": 0.34, "n_studies": 18},
  "axes": {"magnitude_p": 1.2e-9, "rank_rra_p": 3.4e-7, "sign_concordance": 0.94},
  "studies": [
    {"study_id": "HYP001", "pipeline": "DESeq2", "lfc": 2.8, "p": 1e-12, "weight": 0.91, "cell_system": "HeLa", "modality": "1% O2 24h"},
    {"study_id": "HYP002", "pipeline": "edgeR", "logFC": 1.9, "p": 4e-6, "weight": 0.82, "cell_system": "HUVEC", "modality": "CoCl2 100uM 24h"},
    ...
  ]
}
```

This is trivial to implement (it's already in the harmonized parquet), but presenting it as a deliverable changes the perceived value of the output. No competing tool does this.

---

## 5. Revised three-claim contribution sentence

After the upgrades, DEGORA's contribution package becomes:

> *"DEGORA is (1) the first transcriptomic meta-analysis framework designed to operate directly on heterogeneous as-published DEG tables without re-processing, with auto-detection and bias-aware harmonization of the source pipeline; (2) introduces the pipeline-perturbation Jaccard as a generalizable benchmark for evaluating meta-analysis robustness to analytical heterogeneity; and (3) provides a triplet-concordance consensus (magnitude + rank + sign) with auditable per-gene provenance, addressing the rank-information gap flagged as future work by hStouffer (Kim et al., 2026)."*

Auxiliary contribution:
> *"As a byproduct, we release `degpipeline`, a standalone pipeline classifier for GEO supplementary DEG tables, with cross-validated accuracy ≥ 0.95 on a labeled benchmark set."*

This passes a methods reviewer's "what is actually new here?" question without inflating.

---

## 6. Revised V1–V5 validation plan (changes to the methodology .md)

The original validation plan stands, with two additions:

### V6 (NEW) — Pipeline-induced contamination quantification in published meta-analyses (Upgrade 4)

- Reproduce a published meta-analysis from a paper that used heterogeneous DEG inputs (candidates: a psoriasis MetaVolcanoR example, an HCC meta-analysis, any paper with published consensus + raw study list)
- Apply DEGORA and compare consensus hits
- Report: % of original hits that flip status (significant → not significant, or sign reversal), and biological character of the flipped genes

### V7 (NEW) — Pipeline classifier accuracy (Upgrade 3)

- Hand-labeled benchmark: ≥200 DEG tables of known pipeline, with held-out 30% for evaluation
- Report: accuracy, per-pipeline F1, confusion matrix
- Stratify by pipeline version where labels permit

These two additions strengthen the methods-paper claims without changing the test bed (still hypoxia) or the overall structure.

---

## 7. Revised target-venue strategy

After the audit, DEGORA's positioning relative to venues:

| Venue | Fit | Notes |
|---|---|---|
| *Nature Methods* (Brief Comm) | Possible if V2 + V6 are striking | Pipeline-perturbation benchmark + revealed contamination would be the angle |
| *Genome Biology* | Strong fit | Methods + reference implementation + benchmark methodology — their sweet spot |
| *Bioinformatics* | Safest fit | Methodology paper with R/Python implementation |
| *Briefings in Bioinformatics* | Good fit | If framed as benchmark + method |
| *BMC Bioinformatics* | Tertiary | hStouffer just published here three months ago — would need to differentiate hard |

Recommended primary submission: **Genome Biology** (if V6 lands) or **Bioinformatics** (if V6 is mixed). Avoid BMC Bioinformatics as the primary because of hStouffer territory overlap.

---

## 8. Risks introduced by these upgrades

| Risk | Mitigation |
|---|---|
| V6 (pipeline contamination) shows little change in a published consensus | Run early; if negative, demote V6 to a single supplementary figure and reframe contribution toward V2 (pipeline-perturbation benchmark) which doesn't depend on a positive V6 |
| Pipeline classifier accuracy below 0.95 on held-out | Lower the threshold claim to 0.90 or report "high-confidence + low-confidence + unknown" three-class output with separate accuracies; degrade gracefully |
| `degpipeline` standalone is rejected as too small a contribution | Bundle into DEGORA as Component A; do not separate-publish |
| hStouffer authors challenge differentiation in review | Cite favorably, position as complementary regimes (uniform vs heterogeneous input); consider proactive outreach since they are in Korea and may be collaborators not competitors |
| MetaVolcanoR authors note their tool already does multi-method | Pre-empt in Introduction: "MetaVolcanoR provides three strategies as parallel outputs; DEGORA provides a single concordance call requiring all three axes to agree" |
