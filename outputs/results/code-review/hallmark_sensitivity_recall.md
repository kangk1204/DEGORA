# Auxiliary Hallmark sensitivity recompute (AUXILIARY — does NOT replace locked compact panels)

**AUXILIARY sensitivity analysis only.** This broad MSigDB Hallmark ground truth is
provided alongside the compact locked panels so the reader can choose which ground
truth to standardize on. It does NOT replace the locked compact gold panels.

- Hallmark source: Enrichr GMT download MSigDB_Hallmark_2020 (https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020)
- sklearn available: True
- Score column: quality_weighted_piper_score DESC, tie-break quality_weighted_piper_rank ASC

## Hallmark set sizes (raw, from source)
- HALLMARK_HYPOXIA: 200 genes
- HALLMARK_INTERFERON_ALPHA_RESPONSE: 97 genes
- HALLMARK_INTERFERON_GAMMA_RESPONSE: 200 genes
- HALLMARK_UNFOLDED_PROTEIN_RESPONSE: 113 genes

## Comparison table

| topic | ground_truth | corpus | panel_total | n_pos_present | universe | recall@100 | recall@200 | recall@500 | AUROC | AUPRC |
|---|---|---|---|---|---|---|---|---|---|---|
| ER | compact | RNA-only | 18 | 18 | 11910 | 0.6667 | 0.8333 | 0.8333 | 0.9396 | 0.3272 |
| ER | compact | +micro | 18 | 18 | 14239 | 0.8333 | 0.8333 | 0.8889 | 0.9848 | 0.5289 |
| ER | compact | DELTA(+micro - RNA-only) | 18 | 18->18 |  | 0.1667 | 0.0000 | 0.0556 | 0.0452 | 0.2018 |
| ER | hallmark | RNA-only | 113 | 110 | 11910 | 0.1091 | 0.2091 | 0.2727 | 0.6891 | 0.0743 |
| ER | hallmark | +micro | 113 | 112 | 14239 | 0.1875 | 0.2321 | 0.3304 | 0.7392 | 0.0951 |
| ER | hallmark | DELTA(+micro - RNA-only) | 113 | 110->112 |  | 0.0784 | 0.0231 | 0.0576 | 0.0501 | 0.0208 |
| Hypoxia | compact | RNA-only | 20 | 20 | 32687 | 0.7500 | 0.8000 | 0.8500 | 0.9463 | 0.2535 |
| Hypoxia | compact | +micro | 20 | 20 | 34785 | 0.7500 | 0.8000 | 0.8500 | 0.9483 | 0.2673 |
| Hypoxia | compact | DELTA(+micro - RNA-only) | 20 | 20->20 |  | 0.0000 | 0.0000 | 0.0000 | 0.0020 | 0.0138 |
| Hypoxia | hallmark | RNA-only | 200 | 200 | 32687 | 0.2350 | 0.3100 | 0.3800 | 0.8257 | 0.2034 |
| Hypoxia | hallmark | +micro | 200 | 200 | 34785 | 0.2400 | 0.3050 | 0.3750 | 0.8281 | 0.2080 |
| Hypoxia | hallmark | DELTA(+micro - RNA-only) | 200 | 200->200 |  | 0.0050 | -0.0050 | -0.0050 | 0.0023 | 0.0046 |
| IFN | compact | RNA-only | 20 | 18 | 9860 | 1.0000 | 1.0000 | 1.0000 | 0.9989 | 0.5585 |
| IFN | compact | +micro | 20 | 19 | 12777 | 0.9474 | 1.0000 | 1.0000 | 0.9991 | 0.6715 |
| IFN | compact | DELTA(+micro - RNA-only) | 20 | 18->19 |  | -0.0526 | 0.0000 | 0.0000 | 0.0002 | 0.1130 |
| IFN | hallmark | RNA-only | 224 | 165 | 9860 | 0.4303 | 0.6000 | 0.6788 | 0.8513 | 0.5244 |
| IFN | hallmark | +micro | 224 | 194 | 12777 | 0.3763 | 0.5464 | 0.6546 | 0.8563 | 0.5031 |
| IFN | hallmark | DELTA(+micro - RNA-only) | 224 | 165->194 |  | -0.0540 | -0.0536 | -0.0241 | 0.0050 | -0.0214 |

## Notes

- recall@k = |panel ∩ top-k| / n_positives_present (n_positives_present = |panel ∩ scored universe|).
- Background-negative AUROC/AUPRC: positives = panel ∩ universe; negatives = all other scored genes; score = quality_weighted_piper_score.
- Hallmark sets are large (~100-200 genes) so they are not fully present in each corpus universe; see n_pos_present vs panel_total.
- Compact panels are the manuscript-locked gold panels and remain the primary ground truth.

