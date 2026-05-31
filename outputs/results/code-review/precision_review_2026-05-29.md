
---

## UPDATE 2026-05-30 — human-only ER-stress corpus rebuilt (NEW BLOCKER found & fixed)

During follow-up verification a NEW integrity issue was found and fixed: the ER-stress
benchmark was **stale and species-contaminated**. The harmonized table had been regenerated
under the human-only species gate, but the score/table artifacts still contained the mouse
source **GSE84450** (Mus musculus, 11,428 rows) plus GSE102505 — so the manuscript's "human
positive-control" claim was inaccurate for ER, and the ER recall numbers (RNA-only 0.72,
cross-platform 0.89) were computed with mouse evidence mixed in.

Fix applied:
- Added two real human ER/UPR bulk RNA-seq sources (verified, tunicamycin vs vehicle):
  **GSE296996** (HepG2, 5v5) and **GSE245918** (OCI-AML3 + HEK293T, 3v3), curated with the
  same logCPM+Welch recipe as GSE102505. Mouse GSE84450/GSE103667 kept excluded.
- Rebuilt ER primary + cross-platform: score-db, R comparators, gold summaries, deep-metrics,
  cross-platform orchestrator, `make paper`. Validator: failures=0, warnings=0.

ER before -> after (mouse removed):
- RNA-seq only recall@10/50/100: 0.39/0.72/0.72 (8,667 genes, 2 units incl. mouse)
  -> **0.28/0.56/0.67** (11,910 genes, 3 human units); best non-PIPER 0.72 (tie) -> 0.61 (PIPER ahead).
- RNA+microarray recall@10/50/100: 0.28/0.83/0.89 (10,725, 2 units)
  -> **0.39/0.72/0.83** (14,239, 4 units); microarray now improves ER at EVERY cutoff
  (incl. recall@10 0.28->0.39; previously @10 had worsened).
- Hallmark broad-panel sensitivity (auxiliary) confirms direction: ER recall@100 0.273->0.330,
  AUROC 0.689->0.739 with microarray.

main.md synchronized (Abstract, Results recall + cross-platform delta, comparator paragraph,
source-unit counts, heterogeneity-flag ER denominator 10,725->14,239, Methods species-gate +
3-human-unit ER baseline, top-gene lists). No stale ER 0.72/0.89/8667/10725/"2 source units"
figures remain. This resolves the "human-only" framing accuracy issue for ER.
