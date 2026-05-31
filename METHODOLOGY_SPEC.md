# Heterogeneous DEG Integration — Methodology Development Plan v2

> **STATUS (2026-05-31): SUPERSEDED FRAMING — not the submitted contribution.**
> The submitted manuscript (`outputs/manuscript/main.md`) frames DEGORA as a *source-auditable evidence resource/database* built from as-published DEG tables, NOT as a pipeline-aware meta-analysis *method*. The headline method-novelty claims in this document — **N1 (automatic pipeline detection + bias-aware harmonization)** and **N2 (pipeline-perturbation Jaccard benchmark)** — are **NOT implemented** in the current codebase: pipeline labels are read from the user-supplied catalog (`CATALOG_COLUMNS` in `degora/slice_runner.py`; `harmonize.py` does minimal column mapping only, defaulting to `unknown_pipeline`), and the only Jaccard metric is a leave-source-out source-unit ablation (`scripts/write_source_unit_ablation.py`), not a pipeline-perturbation experiment. N1/N2 are retained below only as future-work design notes and are NOT claimed by the manuscript. A prior-art search (2026-05-31) confirmed N1/N2 would be genuinely novel *if built*, but pursuing them requires real new implementation. See `PHASE_A_STATUS.md` and `main.md` for the current, load-bearing contribution.

**Document status:** Revision of v1 (`heterogeneous_deg_meta_methodology_dev.md`) after the novelty audit recorded in `PIPER_NOVELTY_AUDIT.md`. Supersedes v1. The substantive changes from v1 are:

1. **Repositioning** from "robust meta-analysis of heterogeneous studies" to **"the meta-analysis method for as-published heterogeneous DEG tables"** — a narrower, defensible territory.
2. **Component A (pipeline-aware harmonization) promoted to lead** contribution; Components C and E (adaptive weighting, EB shrinkage) demoted to methodological refinements since AWmeta and many tools already cover these.
3. **Pipeline-perturbation Jaccard benchmark elevated** to a standalone methodological contribution.
4. **`degpipeline` standalone classifier** added as a separable deliverable.
5. **Two new validation experiments**: V6 (pipeline-induced contamination quantification in a previously-published consensus) and V7 (classifier accuracy benchmark).
6. **Per-gene provenance graph** added as an output deliverable.
7. **Updated baseline list** to include hStouffer (Kim et al. 2026, BMC Bioinformatics, just 3 months old) and AWmeta (Yang et al. 2025, bioRxiv).
8. **Target venue strategy revised**: Genome Biology / Bioinformatics primary; BMC Bioinformatics demoted due to hStouffer's recent publication in that venue.

Working name remains `DEGORA`. Test bed remains **hypoxia response RNA-seq**.

---

## 0. Three-sentence pitch

DEGORA is the first transcriptomic meta-analysis framework designed to operate directly on the heterogeneous DEG tables published in the supplementary material of RNA-seq papers — produced by different pipelines (DESeq2, edgeR, limma-voom, Cuffdiff, NOISeq, …), with different LFC definitions, different shrinkage behavior, and frequently without standard errors. DEGORA auto-detects the source pipeline of each table, applies bias-aware harmonization to produce a pipeline-invariant representation, and yields a triplet-concordance consensus (magnitude + rank + sign) with auditable per-gene provenance. As a methodological byproduct, DEGORA introduces the *pipeline-perturbation Jaccard* benchmark — a generalizable evaluation protocol that measures any meta-analysis method's robustness to analytical heterogeneity by comparing consensus signatures on the same data analyzed by different pipelines.

---

## 1. Problem statement — precise

### 1.1 The operating regime DEGORA is designed for

A researcher wants to consensus-meta-analyze a biological perturbation (test bed: **hypoxia exposure in mammalian cells**) using studies that have already been published. The realistic data available is:

- For each study $s$: a supplementary DEG table $D_s$ exported by the authors from their analysis software (DESeq2, edgeR, limma-voom, Cuffdiff/Cuffdiff2, NOISeq, sleuth, ballgown, EBSeq, dearseq, or others), in whatever column structure that software produces, with whatever filtering and shrinkage the authors applied.
- Per-paper metadata: sample size, cell system, perturbation specifics — extractable from the paper's methods section.
- Often: no raw count matrix in a usable format. Often: no standard error of LFC. Often: only the "significant" subset (top-K) is given, but a real methods paper requires that we operate only on studies that provide the full DEG table.

This is the actual operating regime of the published literature. DEGORA's premise is that the meta-analysis method should accept this reality rather than require uniform re-processing.

### 1.2 What existing methods require, and why DEGORA's regime is different

| Tool | Input requirement |
|---|---|
| **crossmeta** (Pickering 2018) | Raw data from GEO; performs uniform re-processing across studies |
| **hStouffer** (Kim et al. 2026) | Raw counts from ARCHS4; uniformly DESeq2 re-analyzed |
| **FedPyDESeq2** (2024) | Raw counts at each federated site |
| **MetaIntegrator** (Haynes 2017) | Per-study effect size + variance — assumes consistent definition |
| **MetaVolcanoR** (Prada 2019) | DEG tables with `gene, fold change, p-value` columns — but assumes the FCs are comparable across tools |
| **AWmeta** (Yang 2025) | Per-study p, LFC, within-study variance — assumes consistent meaning |
| **DExMA** (Toro-Domínguez 2022) | DEG tables; handles missing genes |
| **metaRNASeq** (Rau 2014) | DEG p-values; pure p-value combination |
| **DEGORA (this work)** | DEG tables in whatever native column structure their pipeline produced |

This is the gap. Every tool above either re-processes from raw or assumes consistent effect-size definition across studies. None auto-detects the source pipeline or applies pipeline-specific bias correction. DEGORA fills this exact gap.

### 1.3 Output

A consensus signature with four components:

1. **Gene-level table.** For each gene $g$ in $\geq k_{\min}$ studies (default $k_{\min} = 3$): pooled effect $\hat{\theta}_g$ with 95% CI, BH-FDR-adjusted p-value, heterogeneity $I^2_g$, number of contributing studies $N_g$, sign-concordance percentage $\rho_g$, triplet concordance flag $\mathcal{C}_g$.
2. **Pathway-level table.** Hallmark + Reactome + KEGG + WikiPathways with pooled NES, CI, $I^2_P$, redundancy cluster ID.
3. **Per-gene provenance graph (new in v2).** For each consensus gene, a structured record listing every contributing study, its detected pipeline, sample size, cell system, perturbation specifics, contribution weight, and original LFC and p before harmonization. This makes the consensus auditable — any user can trace why a gene is in the top hits.
4. **Diagnostics.** Pipeline detection report (which studies got which pipeline label and at what confidence), pipeline-attributable bias estimates, study-level quality weights, heterogeneity decomposition.

### 1.4 What "robust to pipeline" means operationally

If the same underlying biology is true, the consensus signature should be approximately invariant under:

- The choice of which pipeline was used to produce each contributing study's DEG table
- Replacing a study's DEG table with the output of a different pipeline applied to the same raw data
- Reasonable changes in DEG filtering thresholds within each study

Operational metric: **pipeline-perturbation Jaccard** (the V2 killer experiment) — apply the integration method to (a) the actual heterogeneous DEG tables and (b) re-analyzed homogeneous DEG tables (same studies, one pipeline applied to all), compare top-K Jaccard. Robust methods agree at high Jaccard. Pipeline-naïve methods do not.

---

## 2. Why this is a real gap (positioning vs prior art)

See `PIPER_NOVELTY_AUDIT.md` for the full audit. Summary positioning:

### 2.1 Direct prior art and how DEGORA differs

**MetaVolcanoR (Prada-Medina 2019, Bioconductor).** Provides three meta-analysis *strategies* in parallel: REM, p-value combining, vote-counting. DEGORA differs in (i) pipeline detection and harmonization upstream of these strategies, and (ii) the *concordance call* — a gene is high-confidence only if magnitude, rank, and sign axes all agree, rather than presenting three rankings and letting the user pick.

**AWmeta (Yang 2025 bioRxiv).** Adaptive-weight per-gene per-study weighting on Fisher's and REM. DEGORA's quality weighting is feature-based (sample size, pipeline confidence, DEG table sanity, gene coverage) rather than optimization-based, and is a methodological refinement, not a primary contribution claim.

**hStouffer (Kim et al., *BMC Bioinformatics* 2026; first author D. Kim, senior author J-Y. Lee, KNU).** Three-step stabilization of Stouffer's method (p-value capping + dynamic cutoffs + bagging) plus REM filter. **Key positioning**: hStouffer assumes uniform DESeq2 re-analysis on ARCHS4 counts — it operates on a single pipeline, with the heterogeneity being sample-level not pipeline-level. Kim et al. explicitly state as a limitation: *"the current framework does not incorporate rank-based information, which could provide complementary strengths to p-value and effect-size–based integration."* DEGORA's rank axis directly addresses that limitation, and DEGORA operates in the complementary regime of heterogeneous pipelines.

**crossmeta (Pickering 2018).** Re-processes from raw data with consistent pipeline to avoid heterogeneity. This is the inverse philosophy: crossmeta solves pipeline heterogeneity by eliminating it; DEGORA solves it by modeling it. Both are valid; they cover different use cases (re-processable raw counts vs as-published DEG tables).

### 2.2 The unmatched contributions

After the audit, three contributions are genuinely novel:

- **N1.** Pipeline detection and bias-aware harmonization of as-published DEG tables — no published tool does this.
- **N2.** Pipeline-perturbation Jaccard benchmark methodology — never reported for any meta-analysis method.
- **N3.** Operating regime: working with as-published heterogeneous DEG tables as input, not requiring re-processing or consistent effect-size definitions.

Three supporting contributions strengthen the package without load-bearing alone:

- **N4.** Triplet concordance call — distinct from MetaVolcanoR's parallel-strategies output and AWmeta's two-module sequence.
- **N5.** Per-gene provenance graph as an output deliverable.
- **N6.** `degpipeline` — a standalone pipeline classifier for DEG tables, released as a separable tool.

---

## 3. Test bed — hypoxia response RNA-seq

(Unchanged from v1; reproduced for completeness.)

### 3.1 Why hypoxia

Textbook ground truth (HIF1α target panel is uncontested), multiple pipelines in the literature spanning 2013–present, multiple cell systems (HUVEC, HeLa, MCF7, A549, RCC4, U87, primary fibroblasts), multiple perturbation modalities (physical 1%, 0.5%, 0.1% O₂; chemical CoCl₂, DFO, DMOG; genetic VHL knockout).

### 3.2 Gold-standard targets — locked before any analysis

These define recall at evaluation time. Locked. No post-hoc additions allowed (integrity checklist item).

**Canonical HIF1α direct targets — upregulated:**

- Glycolysis switch: `SLC2A1` (GLUT1), `HK1`, `HK2`, `PFKL`, `ALDOA`, `ENO1`, `PGK1`, `LDHA`, `PDK1`
- Vascular/angiogenic: `VEGFA`, `ANGPT2`, `FLT1`
- Hematopoietic: `EPO` (tissue-dependent)
- pH regulation: `CA9`, `SLC2A3`
- Autophagy/apoptosis: `BNIP3`, `BNIP3L`, `DDIT4` (REDD1)
- HIF feedback: `EGLN1` (PHD2), `EGLN3` (PHD3)

**Canonical downregulated under prolonged hypoxia:** OXPHOS gene set (Hallmark `HALLMARK_OXIDATIVE_PHOSPHORYLATION`), `PPARGC1A` (PGC-1α) in some contexts.

**Negative-control set:** `RPL13A`, `HPRT1`, `TBP`. (`ACTB` and `GAPDH` are borderline — possible context-dependent HIF effects; either include with caveat or exclude.)

**Pathway-level gold standard:** Hallmark `HALLMARK_HYPOXIA` top-1, Reactome `Cellular response to hypoxia` (R-HSA-1234174) and KEGG `HIF-1 signaling pathway` (hsa04066) in top-5.

### 3.3 Study collection protocol

**Target N:** 20–30 hypoxia RNA-seq studies. **Inclusion criteria (locked):**

1. Bulk RNA-seq.
2. Hypoxia vs. normoxia control with ≥ 3 replicates per arm.
3. **Full DEG table in supplementary** — ≥ 10,000 genes reported, not top-K.
4. Pipeline explicitly stated in methods section.
5. Human or mouse.
6. Published 2013 or later.

**Exclusion:** hypoxia + another perturbation without clean hypoxia-only arm; hypoxia as secondary perturbation in a disease cohort; fold change without significance measure; DEG table locked in non-machine-readable format.

**Curation deliverable per study** (CSV `data/studies/catalog.csv`): `study_id, geo_accession, doi, species, cell_system, hypoxia_modality, duration_h, n_ctrl, n_treat, pipeline, pipeline_version, gene_id_system, deg_table_path, n_genes_reported, has_se, has_stat, lfc_column, p_column, padj_column, notes`.

**Time estimate:** 1–2 weeks of manual curation. This is the longest non-automatable bottleneck — but only ~30 studies, fully tractable.

### 3.4 Pipeline re-analysis subset for V2

For the pipeline-perturbation experiment: 8–12 studies with publicly available raw counts (GEO supplementary, recount3, or ARCHS4). Each re-analyzed with DESeq2, edgeR, and limma-voom under matched filtering. ERAPID (the user's tool) is reused as the reproducibility-friendly DESeq2 wrapper here.

---

## 4. Methodology — DEGORA

Five components. Components A and B are the load-bearing contributions; C–E are methodological refinements that complete the framework.

### 4.1 Component A — Pipeline detection and harmonization **(primary contribution)**

#### A.1 Pipeline detection (`piper.detect.pipeline_signature` / `degpipeline.predict`)

Takes a DEG table; returns a pipeline label with confidence. Two signals:

**Column-name fingerprint.** Lookup table over standard pipeline outputs:

- DESeq2 (raw): `{baseMean, log2FoldChange, lfcSE, stat, pvalue, padj}`
- DESeq2 with `lfcShrink` (apeglm): `{baseMean, log2FoldChange, lfcSE, pvalue, padj}` (no `stat`)
- DESeq2 with `lfcShrink` (ashr): same as apeglm but distributional fingerprint differs
- edgeR (glmQLF): `{logFC, logCPM, F, PValue, FDR}` optionally `unshrunk.logFC`
- edgeR (LRT): `{logFC, logCPM, LR, PValue, FDR}`
- limma-voom: `{logFC, AveExpr, t, P.Value, adj.P.Val, B}`
- Cuffdiff/Cuffdiff2: `{test_id, gene_id, gene, locus, sample_1, sample_2, status, value_1, value_2, log2(fold_change), test_stat, p_value, q_value, significant}`
- NOISeq: `{theta, prob, ranking}`
- sleuth: `{target_id, qval, pval, b, se_b, mean_obs, var_obs, tech_var, sigma_sq}`
- ballgown: `{feature, fc, pval, qval}`

**Distributional fingerprint.** When column names are renamed or ambiguous: extract features from the table itself.

- Distribution of LFC at low expression (DESeq2 raw shows inflation; shrunken variants don't)
- Relationship of |LFC| to −log10(p) (different curvature per pipeline)
- Proportion of zero-LFC or exactly-1 fold-change rows
- p-value distribution shape (Cuffdiff has characteristic spikes)
- Presence of extreme values (Cuffdiff can produce ±inf)
- Number of rows reported (some pipelines filter aggressively by default)

An XGBoost classifier trained on a labeled set of DEG tables (≥200 tables, ~20 per pipeline; ≥30% held out for evaluation) gives a probabilistic label `{DESeq2_raw, DESeq2_shrunk_apeglm, DESeq2_shrunk_ashr, edgeR_QLF, edgeR_LRT, limma_voom, limma_trend, cuffdiff, cuffdiff2, NOISeq, sleuth, ballgown, EBSeq, dearseq, unknown}`. Confidence < 0.7 → fall through to `unknown_pipeline` and apply conservative handling.

#### A.2 Harmonization (`piper.harmonize`)

Produces a uniform per-gene per-study representation. For each $(g, s)$:

- $\text{lfc}_s^{(g)}$: raw LFC, preserved for reporting
- $\text{lfc}_s^{*(g)}$: pipeline-bias-corrected LFC. Pipeline-specific calibration rules:
  - DESeq2 raw at low `baseMean`: downweight via `apeglm`-style adaptive shrinkage post-hoc, or fall back to stat column if present
  - DESeq2 shrunk: keep as-is (already corrected)
  - edgeR: keep as-is; mark studies with very low `logCPM` filter as having reduced gene coverage
  - limma-voom: keep as-is (already empirical-Bayes stabilized)
  - Cuffdiff: cap $|\text{LFC}| \leq 10$; drop genes with `status != OK`
  - NOISeq: convert `prob` to pseudo-p
  - sleuth: use `b` (beta) as LFC analog
  - unknown: use as-is with downweight
- $\text{se}_s^{(g)}$: SE if available; otherwise imputed from $|z|$ and LFC: $\text{se} \approx |\text{LFC}| / |z|$ where $z = \Phi^{-1}(1 - p/2) \cdot \text{sign}(\text{LFC})$. Studies with imputed SE are flagged.
- $z_s^{(g)} = \text{sign}(\text{LFC}_s^{(g)}) \cdot \Phi^{-1}(1 - p_s^{(g)}/2)$ — signed-z, the pipeline-invariant magnitude proxy.
- $\text{rank}_s^{(g)}$: within-study rank by $|z_s^{(g)}|$.
- $\text{sign}_s^{(g)} \in \{+1, -1, 0\}$: sign of LFC if padj < 0.1; else 0.

Output: long-format `harmonized.parquet` with per-row provenance.

### 4.2 Component B — Multi-evidence triplet and concordance call **(primary contribution refinement)**

For each gene $g$ in $N_g \geq k_{\min}$ studies, compute three orthogonal scores.

**Magnitude axis.** Random-effects meta-analysis on $(\text{lfc}_s^{*(g)}, \text{se}_s^{(g)})$ via DerSimonian–Laird or REML (`metafor` via rpy2 or `statsmodels.stats.meta_analysis`). Output: pooled $\hat{\theta}_g^{\text{mag}}$, 95% CI, p-value, $I^2_g^{\text{mag}}$.

**Rank axis.** Robust Rank Aggregation (Kolde 2012; `RobustRankAggreg` or Python re-implementation) on within-study rankings. Output: $\hat{p}_g^{\text{rank}}$.

**Sign axis.** Weighted vote: $\text{sign}_g = \sum_s w_s \cdot \text{sign}_s^{(g)} / \sum_s w_s$ with $w_s$ from §4.3. Sign-consistency: fraction of contributing studies agreeing with consensus sign.

**Triplet concordance.** The key differentiation from MetaVolcanoR (which provides three parallel rankings) and AWmeta (which sequences p-value then REM):

$$
\mathcal{C}_g = \mathbb{1}\{\hat{p}_g^{\text{mag, BH}} < \alpha_1\} \cdot \mathbb{1}\{\hat{p}_g^{\text{rank, BH}} < \alpha_2\} \cdot \mathbb{1}\{\rho_g^{\text{sign}} > \alpha_3\}
$$

with defaults $\alpha_1 = \alpha_2 = 0.05$, $\alpha_3 = 0.75$. **High-confidence consensus = all three axes agree.** This is reported as a single boolean in the gene-level output, alongside the three individual axes for auditability.

This implements the rank-information integration that hStouffer's authors flagged as future work.

### 4.3 Component C — Quality-feature-based study weights (refinement, not primary novelty)

Per-study weight:

$$
w_s = \sqrt{n_s} \cdot q_s^{\text{pipeline}} \cdot q_s^{\text{sanity}} \cdot q_s^{\text{coverage}}
$$

- $q_s^{\text{pipeline}} \in [0.5, 1.0]$: known/well-supported (DESeq2, edgeR, limma-voom) → 1.0; legacy (Cuffdiff1, NOISeq) → 0.75; unknown → 0.5.
- $q_s^{\text{sanity}}$: penalty for implausible fraction of significant genes (< 0.5% or > 30%).
- $q_s^{\text{coverage}}$: fraction of a reference gene set present in the table.

Reported per study in the consensus output for auditability. **Not framed as a novel contribution in the paper** — AWmeta has already established adaptive weighting as a domain capability.

### 4.4 Component D — Hierarchical gene → pathway integration (auxiliary)

Two-tier integration. Tier 1 = gene-level (§4.2). Tier 2 = pathway-level: ssGSEA / fgsea NES per study, then random-effects MA across studies, then BH-FDR within pathway library, then Jaccard-based redundancy deduplication.

Surface a per-gene diagnostic: `pathway_consistency ∈ {concordant, gene_only, pathway_only, neither}` indicating whether gene-level and pathway-level evidence agree.

Auxiliary contribution; not primary novelty.

### 4.5 Component E — Empirical-Bayes shrinkage (auxiliary)

$$
\hat{\theta}_g^{\text{shrunk}} = \frac{\tau^2 \hat{\theta}_g^{\text{mag}} + \sigma_g^2 \bar{\theta}_{\mathcal{P}(g)}}{\tau^2 + \sigma_g^2}
$$

Shrinks per-gene effect toward pathway mean weighted by uncertainty. Reported with 95% credible interval. Standard empirical-Bayes; not framed as novel.

### 4.6 Output deliverable — per-gene provenance graph (new in v2)

For every consensus gene, the consensus output includes a structured provenance record (JSON or parquet long format):

```json
{
  "gene": "VEGFA",
  "consensus": {"lfc": 2.31, "ci": [1.80, 2.83], "i2": 0.34, "n_studies": 18},
  "axes": {
    "magnitude_p_bh": 1.2e-9,
    "rank_rra_p_bh": 3.4e-7,
    "sign_concordance": 0.94,
    "triplet_concordant": true
  },
  "studies": [
    {
      "study_id": "HYP001", "geo": "GSEnnnnnn", "doi": "...",
      "pipeline": "DESeq2_raw", "pipeline_confidence": 0.98,
      "lfc_raw": 2.85, "lfc_harmonized": 2.71, "se": 0.21, "p": 1e-12, "padj": 3e-10,
      "weight": 0.91,
      "cell_system": "HeLa", "modality": "1% O2 24h", "n_ctrl": 3, "n_treat": 3
    },
    ...
  ]
}
```

No competing tool outputs anything like this. Cheap to implement; high perceived value for downstream users (a wet-lab biologist can see exactly which papers support a gene's consensus before deciding to follow up).

---

## 5. The `degpipeline` standalone classifier (new deliverable in v2)

The pipeline detector from §4.1 is released as a standalone Python package independent of the rest of DEGORA:

```python
from degpipeline import detect_pipeline

result = detect_pipeline("supplementary_table_3.tsv")
# {pipeline: "DESeq2_raw", confidence: 0.97, version_hint: ">=1.20", flags: []}
```

Use cases beyond DEGORA: GEO curators auditing supplementary data, paper reviewers checking methods consistency, automated meta-analysis pipelines that need to dispatch tool-specific processing, training-data labeling for transcriptomics ML projects.

Released under MIT or Apache-2.0; potentially a separate short publication (*Bioinformatics* Application Note) if labeled benchmark accuracy lands ≥ 0.95.

Implementation: XGBoost trained on a benchmark of ≥ 200 hand-labeled DEG tables (collected during S0 catalog curation), with ≥ 30% strict held-out evaluation. Features = column-name fingerprint (one-hot) + 12 distributional features extracted from the DEG table itself.

---

## 6. Baselines (mandatory comparisons)

All baselines implemented in `piper.baselines` with a uniform interface. The baseline list now reflects the audit:

| Baseline | Implementation | Audit role |
|---|---|---|
| Vote counting (sign concordance threshold) | Custom Python | Floor baseline |
| Fisher's combined p | `scipy.stats.combine_pvalues(method='fisher')` | Classical floor |
| Stouffer's weighted Z (sqrt(N) weights) | `scipy.stats.combine_pvalues(method='stouffer', weights=...)` | Classical floor |
| Random-effects MA on LFC + SE (DerSimonian–Laird) | `metafor::rma` via rpy2 | Statistical gold standard |
| Random-effects MA on LFC + SE (REML) | `metafor::rma` | Statistical gold standard |
| **MetaVolcanoR** (REM, p-comb, vote-counting modes) | `MetaVolcanoR` via rpy2 | Direct prior art — closest in spirit |
| **hStouffer** (p-cap + dynamic cutoff + bagging) | from `CR4CKID/hStouffer` GitHub | Direct prior art — published 3 months ago, must beat or match |
| **AWmeta** (AW-Fisher + AW-REM) | Re-implement from bioRxiv preprint | Direct prior art — adaptive weighting state of the art |
| Robust Rank Aggregation | `RobustRankAggreg` via rpy2 | Component of DEGORA and a baseline |
| MetaIntegrator (DSL REM + Fisher) | `MetaIntegrator` via rpy2 | Established adjacent |
| metaRNASeq (Fisher + inverse normal) | `metaRNASeq` from R-Forge | Foundational |
| RankProd | `RankProd` Bioconductor | Rank-based alternative |
| DExMA | `DExMA` Bioconductor | Handles missing genes — adjacent |
| crossmeta (re-process from raw) | `crossmeta` Bioconductor | Inverse-philosophy comparison; runs only on the V2 subset where raw counts exist |

DEGORA must beat the three direct-prior-art baselines (MetaVolcanoR, hStouffer, AWmeta) on at least one validation criterion. Realistic expectation: DEGORA wins clearly on V2 (pipeline-perturbation Jaccard); competitive on V1; potentially differentiated by V6 (contamination experiment). Honest reporting of where baselines tie or win is required by the integrity checklist.

---

## 7. Validation plan

Seven experiments. All seven committed before any single one is run. Pre-registered hypotheses where applicable.

### 7.1 V1 — Positive-control recovery (unchanged from v1)

**Metric.** Recall@50 and Recall@100 of the gold-standard HIF1α target list (§3.2) in each method's consensus top-K. AUROC over the full gold list (positives) vs. random negatives.

**Expected outcome.** All competent methods recover the gold list at high recall — this is sanity, not differentiation. Methods that fail this are broken.

### 7.2 V2 — Pipeline-perturbation Jaccard (the headline experiment)

**Setup.**
1. 8–12 studies from the catalog with public raw counts.
2. Re-analyze each with DESeq2, edgeR, and limma-voom under matched filtering.
3. Four DEG collections:
   - $D_{\text{DESeq2}}$, $D_{\text{edgeR}}$, $D_{\text{limma}}$ — homogeneous-pipeline collections
   - $D_{\text{mixed}}$ — random pipeline assignment per study, mimicking the realistic heterogeneous input
4. Run each meta-analysis method on each collection.
5. Pairwise Jaccard@50, Jaccard@100, Spearman of signed-z over the four collections.

**Pre-registered hypothesis.** DEGORA achieves mean pairwise Jaccard@50 ≥ 0.85; classical methods (Fisher, Stouffer, REM-without-pipeline-aware-harmonization) achieve ≤ 0.70.

**This is the figure 1 of the paper.** If the hypothesis fails, report honestly and reframe contribution toward V6 if that succeeds.

### 7.3 V3 — Leave-one-study-out replication (unchanged from v1)

Leave each study out; predict its top-100 by signed-z from the consensus of the other $N-1$ studies. Mean Recall@100.

### 7.4 V4 — Simulation with synthetic ground truth (unchanged from v1)

Generate true LFC vector with sparse non-zero entries. Sample 20 "studies" with realistic sample sizes. Inject pipeline-specific noise. Compute recovery correlation, AUROC, **CI coverage probability** (does the reported 95% CI contain $\theta^*_g$ 95% of the time?). The CI coverage check is rarely tested in meta-analysis papers — by itself it's a publishable contribution.

### 7.5 V5 — Cross-domain generalizability (unchanged from v1)

Apply DEGORA as-is to 5–10 **TNF-α stimulation** RNA-seq studies. Canonical NF-κB targets (`NFKBIA, TNFAIP3, CXCL10, CCL2, ICAM1, BIRC3, ...`). Recall@50. Convinces a reviewer that DEGORA is not over-tuned to hypoxia.

### 7.6 V6 — Pipeline-induced contamination in a published meta-analysis (NEW in v2)

**Question.** When a published transcriptomic meta-analysis used heterogeneous DEG inputs, how much of its reported consensus is potentially pipeline-induced artifact?

**Setup.**
1. Identify a previously published meta-analysis that used heterogeneous DEG tables as input (candidates: a MetaVolcanoR-based study, an AWmeta application paper, or any review meta-analysis on a well-defined biological question where the raw DEG tables of contributing studies are still accessible).
2. Reproduce the published consensus using the same DEG tables.
3. Re-run with DEGORA (full pipeline detection + harmonization).
4. Compare top-K hits: how many original hits change significance status? How many change sign? What is the biological character of the flipped genes?

**Metric.** % of original top-100 hits that change status (significant ↔ not significant, or sign reversal). Stratify by pipeline composition of the original study set.

**Hypothesis.** ≥ 10% of original consensus hits change status when pipeline heterogeneity is corrected for. If this holds, it's a publishable empirical finding regardless of DEGORA's methodological standing.

**Risk.** If the answer is "almost no change", V6 cuts against DEGORA's framing. Mitigation: run V6 in iteration 2 or 3, before final paper framing. If negative, demote V6 to a supplementary "no contamination detected in case study X" finding and lean on V2 instead.

### 7.7 V7 — `degpipeline` classifier accuracy (NEW in v2)

**Setup.** Hand-labeled benchmark of ≥ 200 DEG tables across ≥ 10 pipeline classes. 70/30 train/test split, no information leakage (e.g., no two tables from the same study in different splits).

**Metrics.** Overall accuracy, per-pipeline F1, confusion matrix, calibration (does predicted confidence match empirical accuracy?). Stratify by pipeline version where labels permit.

**Pre-registered target.** Overall accuracy ≥ 0.95 on held-out set, ≥ 0.90 per major class (DESeq2, edgeR, limma).

This experiment underwrites the `degpipeline` standalone deliverable.

---

## 8. Implementation stages (revised)

| Stage | Deliverable | Pass criterion | Effort |
|---|---|---|---|
| S0 | Study catalog (20–30 hypoxia + 8–12 raw-count-available + 5–10 TNF-α for V5 + ≥ 200 hand-labeled DEG tables for V7) | All entries downloaded, parsed, labeled | 2–3 weeks (the V7 labeling adds ~1 week to v1's estimate) |
| S1 | Baseline implementations including MetaVolcanoR, hStouffer, AWmeta | Each baseline returns consensus in uniform schema | 1–2 weeks (AWmeta re-implementation is the new effort) |
| S2 | DEGORA components A–E + `degpipeline` classifier | A classifier ≥ 0.95 held-out accuracy (V7 pass); unit tests for each component | 2–3 weeks |
| S3 | V1 + V3 experiments | Numerical results land; integrity checklist passed | 1 week |
| S4 | V2 pipeline-perturbation experiment (the headline) | Re-analysis with 3 pipelines on 8–12 studies; Jaccard matrix; figure 1 produced | 2 weeks |
| S5 | V4 simulation framework | Synthetic data + benchmarking + CI calibration plot | 1 week |
| S6 | V5 cross-domain demo (TNF-α) | Recall numbers in supplementary figure | 3 days |
| S7 | V6 contamination experiment (NEW) | Re-analysis of one published meta-analysis with DEGORA; % flipped reported | 1 week — schedule in iteration 2 so framing can adapt to result |
| S8 | `degpipeline` standalone package release (NEW) | PyPI package, docs, CI, GitHub README; optionally an Application Note draft | 1 week |
| S9 | Manuscript draft + revisions | Three publication gates pass | 3 weeks |

**Total realistic timeline:** 11–14 weeks of focused effort, vs. 8–12 in v1. The increase reflects S0's expanded labeling, S1's added baselines, S7 and S8 new stages.

---

## 9. Code structure (revised)

```
piper/
├── piper/                       # Main importable package
│   ├── catalog/                 # Study catalog + DEG table loading
│   ├── detect/                  # Component A: pipeline detection (uses degpipeline)
│   ├── harmonize/               # Component A: harmonization
│   ├── triplet/                 # Component B: magnitude/rank/sign + concordance call
│   ├── weights/                 # Component C: quality weights
│   ├── pathway/                 # Component D: pathway integration
│   ├── shrink/                  # Component E: EB shrinkage
│   ├── provenance/              # NEW: per-gene provenance graph (§4.6)
│   ├── baselines/               # §6 baseline wrappers (MetaVolcanoR, hStouffer, AWmeta, ...)
│   ├── benchmark/               # V1–V7 experiment runners
│   └── simulate/                # V4 synthetic DEG simulator
├── degpipeline/                 # NEW: standalone pipeline classifier subpackage
│   ├── features/                # Column-name fingerprint + distributional features
│   ├── model/                   # XGBoost classifier and trained weights
│   ├── labeled_benchmark/       # Curated labeled DEG tables for V7
│   └── cli.py                   # Standalone CLI entry
├── data/
│   ├── studies/catalog.csv
│   ├── deg/                     # Standardized per-study DEG parquet
│   ├── reanalysis/              # V2 outputs from DESeq2/edgeR/limma re-analyses
│   ├── pathways/                # MSigDB / Reactome / KEGG / WikiPathways gmt
│   └── v6_case_study/           # Published meta-analysis being audited
├── scripts/                     # One per stage
├── tests/                       # pytest, qa-tester skill
├── results/
├── manuscript/
├── Makefile                     # make benchmark, make figs, make paper, make classifier
├── requirements.txt
├── environment.yml
└── README.md
```

Stack: Python 3.10+, pandas, numpy, scipy, statsmodels, XGBoost, gseapy; R via rpy2 for `metafor`, `RobustRankAggreg`, `MetaVolcanoR`, `MetaIntegrator`, `metaRNASeq`, `RankProd`, `DExMA`, `crossmeta`; hStouffer via its GitHub repo.

---

## 10. Manuscript outline (revised)

### Target venues (priority revised)

1. **Genome Biology** (primary target). Strong fit for "method + reference implementation + benchmark methodology". V2 + V6 together would be their sweet spot.
2. **Bioinformatics** (Oxford). Safest fit; methods paper with reference implementation.
3. **Briefings in Bioinformatics**. If V6 is positive and we frame as benchmark + method.
4. **Nature Methods (Brief Communication)**. Only if V6 is striking (≥ 20% of consensus flips) and V2 hypothesis is met cleanly.
5. **BMC Bioinformatics** (deprioritized). hStouffer just published here three months ago — too much territory overlap; would require sharp differentiation framing that's already harder than choosing a different venue.

Auxiliary: `degpipeline` as a separate *Bioinformatics* Application Note if V7 accuracy is ≥ 0.95.

### Structure

- **Abstract** (~250 words, written last). Lead with the operating-regime gap: "Existing transcriptomic meta-analysis methods assume either re-processing from raw data or DEG tables with consistent effect-size definitions. In practice, the supplementary DEG tables in published RNA-seq papers are heterogeneous in pipeline. We present DEGORA..."
- **Introduction.**
  - Frame the gap (§2) with explicit reference to crossmeta/hStouffer/AWmeta/MetaVolcanoR positioning.
  - State the four contributions (N1–N4 from §2.2) precisely.
  - Confirm via `novelty-checker` skill no closer competitor has emerged since this audit.
- **Results.**
  - 2.1 Motivating analysis — pipeline choice systematically alters DEG outcomes (cite/recapitulate arxiv:2601.04122 evidence + own analysis on a couple of test studies)
  - 2.2 DEGORA architecture overview (Figure 1: pipeline detection → harmonization → triplet → consensus + provenance)
  - 2.3 V2 — Pipeline-perturbation Jaccard (the headline figure; mean Jaccard table across methods × collections)
  - 2.4 V1 + V3 + V4 — Recovery, replication, simulation
  - 2.5 V6 — Pipeline-induced contamination in published meta-analysis (if positive)
  - 2.6 V5 — Cross-domain TNF-α demo (supplementary figure)
  - 2.7 V7 + `degpipeline` — pipeline classifier accuracy + use cases
- **Discussion.** Honest limitations: pipeline detection extensible but heuristic; SE imputation is approximate; test bed is two perturbation domains; full Bayesian hierarchical is future work; nutrient-perturbation database (NutriOmics roadmap) is downstream application.
- **Methods.** Reproducible from manuscript + code repo alone. Component-by-component math.
- **Data and code availability.** GitHub + Zenodo for catalog and re-analyzed DEG tables. `degpipeline` on PyPI.
- **References.** Verify every one with `ref-finder`.

### Anticipated reviewer questions (revised)

| Anticipated question | Pre-empted response |
|---|---|
| "How is this different from MetaVolcanoR?" | MetaVolcanoR provides three strategies as parallel rankings; DEGORA provides a *concordance call* requiring all three axes to agree. DEGORA additionally has pipeline detection and harmonization upstream of the consensus, which MetaVolcanoR does not. |
| "How is this different from AWmeta (Yang 2025)?" | AWmeta's adaptive weighting is optimization-based and operates on consistent effect-size definitions; DEGORA's weighting is feature-based and explicitly handles pipeline-specific bias. AWmeta does not detect or harmonize pipelines. |
| "How is this different from hStouffer (Kim et al. 2026)?" | hStouffer addresses FP inflation at large dataset sizes with uniform DESeq2 re-analysis on ARCHS4; DEGORA addresses pipeline heterogeneity in as-published DEG tables. The regimes are complementary; in fact, hStouffer's authors explicitly flag rank-information integration as future work, which DEGORA's rank axis addresses. |
| "Why not just use crossmeta and re-process everything?" | crossmeta requires raw count availability, which is missing for a substantial fraction of published studies; DEGORA works on the studies for which only the supplementary DEG table is accessible. Both tools are valid for different inputs. |
| "Pipeline detection accuracy?" | V7 reports held-out accuracy with confusion matrix; degrade gracefully via `unknown_pipeline` label when confidence is low. |
| "Test bed is narrow." | V4 simulation + V5 cross-domain TNF-α + V6 case study together address this. |
| "Curated study catalog is small." | Methods paper, not resource paper. 20–30 studies plus 8–12 re-analyzed subset plus 5–10 cross-domain is appropriate scale for benchmarking. The downstream application (nutrient perturbation database) is a separate future deliverable, not part of this paper's claim. |

---

## 11. Risks and mitigations (revised)

| Risk | Likelihood | Mitigation |
|---|---|---|
| Full DEG tables hard to find for 20+ hypoxia studies | Medium | Drop to 15+ with documented exclusions; supplement with recount3 re-analysis if needed |
| V2 pipeline-perturbation Jaccard gap (DEGORA vs baselines) small | Medium-high | Report honestly; if baselines are already at ≥ 0.85, the contribution becomes "first systematic quantification of pipeline-induced variance in meta-analysis methods, and shows existing methods are more robust than expected" — still publishable, different framing |
| V6 contamination experiment shows little flip | Medium | Run early; if negative, demote to supplementary; lean on V2 |
| `degpipeline` classifier accuracy < 0.95 | Low-medium | Lower target to 0.90 with three-class output (high-confidence / low-confidence / unknown); separately publish only if ≥ 0.95 |
| hStouffer authors challenge differentiation in review | Medium | Cite favorably; position as complementary regimes; **consider proactive outreach to the hStouffer group (Kim et al.; senior author J-Y. Lee, KNU) for informal review or co-authorship**, since their stated future-work directly aligns with DEGORA |
| MetaVolcanoR authors note their tool already does multi-method | High | Pre-empt explicitly in Introduction: parallel-strategies vs concordance-call; pipeline detection upstream |
| AWmeta cited as having already solved adaptive weighting | High | Pre-empt: DEGORA's weighting is feature-based (auditable) vs AWmeta's optimization-based; demote weighting from claim to refinement |
| EB shrinkage adds little vs simple weighted mean | Medium | Drop Component E; report negatively as "shrinkage does not materially improve in this regime" — valid scientific finding |
| TNF-α V5 demo fails | Low | Investigate and document; narrow claim to "hypoxia-class transcriptional responses" if needed |
| Reviewer pushback that contribution is "just engineering" | Medium | Frame V2 pipeline-perturbation Jaccard benchmark methodology as the novel scientific contribution; engineering is the implementation |

---

## 12. Outreach consideration — hStouffer authors

Kim et al. (first author Daehee Kim, senior/corresponding author Jun-Yeong Lee) at Kyungpook National University (KNU) published hStouffer in BMC Bioinformatics on 2026-02-12 (DOI 10.1186/s12859-026-06395-2, PMID 41680606). They are in Korea, in the same broad domain, and they explicitly identified rank-based information as a future direction — which is precisely what DEGORA's triplet rank axis addresses.

Three possible postures:

- **Competitive.** Cite, differentiate, submit independently. Lowest risk of friction but loses a potential ally and a domain expert review.
- **Friendly.** Send a preprint draft to Kim/Lee as a "this builds on what you identified as future work; would you be willing to give us informal feedback" message. No commitment, no co-authorship request, but builds rapport.
- **Collaborative.** Approach for a discussion about whether DEGORA + hStouffer could be a joint paper or whether the methodologies are better together. Highest potential payoff (regime complementarity: hStouffer for uniform-pipeline large-n, DEGORA for heterogeneous-pipeline moderate-n) but requires shared work.

This is a strategic call for the user. Recommended default: posture 2 (friendly preprint share) after V2 results are in but before submission. The user can request an email draft when ready.

---

## 13. How to start

S0 (study catalog + DEG table labeling for V7) requires human curation and is the longest non-automatable step. Budget 2–3 weeks before Codex CLI / Claude Code automation takes over.

Then either:

**(a) Manual progression.** Use this v2 document as a reference; build the repo according to §9; advance stage by stage with Codex CLI as a coding assistant on individual scripts.

**(b) Autonomous loop.** Wrap as `INSTRUCTIONS.md` and create the supporting workspace (`CONTEXT.md` adapting this document, plus the empty scaffolds for `STATE.md`, `PLAN.md`, `RUBRIC.md`, `LAB_NOTEBOOK.md`, `REVIEWER_LOG.md`). The user has the `nutriomics_db_phase1` template from a previous session as a reference. Then start Codex CLI with: *"Read INSTRUCTIONS.md and begin the autonomous research loop."*

**Integrity reminders for any execution mode:**

- Gold-standard target list (§3.2) is locked. No additions after running any analysis.
- Inclusion/exclusion criteria (§3.3) are locked.
- Validation experiments V1–V7 are committed. All run regardless of how individual results look.
- No baseline gets dropped to make DEGORA look better. Ties and losses are reported.
- V6 is run in iteration 2 or 3 so the framing can adapt to its outcome before final paper writing — but the experiment itself runs and is reported regardless.

---

## Appendix A — Three-claim contribution sentence (final form)

> *"DEGORA is (1) the first transcriptomic meta-analysis framework designed to operate directly on heterogeneous as-published DEG tables, with auto-detection and bias-aware harmonization of the source pipeline; (2) introduces the pipeline-perturbation Jaccard as a generalizable benchmark for evaluating meta-analysis robustness to analytical heterogeneity; and (3) provides a triplet-concordance consensus (magnitude + rank + sign) with auditable per-gene provenance, addressing the rank-information integration gap flagged as future work by Kim et al. (2026, BMC Bioinformatics)."*

Auxiliary contribution: *"As a byproduct, we release `degpipeline`, a standalone XGBoost-based pipeline classifier for as-published DEG tables, achieving cross-validated accuracy ≥ 0.95 on a labeled benchmark."*

## Appendix B — Quick math summary

Signed-z (pipeline-invariant magnitude proxy):
$$z_s^{(g)} = \text{sign}(\text{LFC}_s^{(g)}) \cdot \Phi^{-1}(1 - p_s^{(g)}/2)$$

DEGORA feature-based study weight (Component C):
$$w_s = \sqrt{n_s} \cdot q_s^{\text{pipeline}} \cdot q_s^{\text{sanity}} \cdot q_s^{\text{coverage}}$$

Random-effects MA pooled effect (Component B magnitude axis):
$$\hat{\theta}_g^{\text{mag}} = \frac{\sum_s w_s^{\text{REM}} \theta_s^{(g)}}{\sum_s w_s^{\text{REM}}}, \quad w_s^{\text{REM}} = \frac{1}{\hat{\tau}_g^2 + \sigma_s^{(g)\,2}}$$

Triplet concordance flag (Component B):
$$\mathcal{C}_g = \mathbb{1}\{\hat{p}_g^{\text{mag, BH}} < \alpha_1\} \cdot \mathbb{1}\{\hat{p}_g^{\text{rank, BH}} < \alpha_2\} \cdot \mathbb{1}\{\rho_g^{\text{sign}} > \alpha_3\}$$

EB shrinkage (Component E):
$$\hat{\theta}_g^{\text{shrunk}} = \frac{\tau^2 \hat{\theta}_g + \sigma_g^2 \bar{\theta}_{\mathcal{P}(g)}}{\tau^2 + \sigma_g^2}$$
