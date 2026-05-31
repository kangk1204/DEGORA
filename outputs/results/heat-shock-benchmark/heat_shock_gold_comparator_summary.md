# Heat shock locked gold-panel comparator summary

Locked positives: `data/studies/gold/heat_shock_hsf1_gold_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/quality_weighted_default | ok | 14285 | 0.31 (0.14-0.52) | 0.38 (0.20-0.59) | 0.58 (0.37-0.77) | 0.73 (0.52-0.88) | 0.19 | 0.00 | 0.73 |  | HSPA6;HSPA1A;HSPA1B;DNAJB1;BAG3;HSPA1L;HSPH1;FOS;DNAJB4;MICB |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 14285 | 0.31 (0.14-0.52) | 0.38 (0.20-0.59) | 0.58 (0.37-0.77) | 0.73 (0.52-0.88) | 0.19 | 0.00 | 0.73 |  | HSPA6;HSPA1A;HSPA1B;DNAJB1;BAG3;HSPA1L;HSPH1;FOS;DNAJB4;MICB |
| degora_slice/locked | ok | 14285 | 0.23 (0.09-0.44) | 0.35 (0.17-0.56) | 0.50 (0.30-0.70) | 0.58 (0.37-0.77) | 0.15 | 0.00 | 0.58 |  | HSPA6;HSPA1A;DNAJB1;BAG3;CCDC121;HSPA1B;MICB;ACTRT3;LCMT2;HSPH1 |
| weighted_stouffer/default | ok | 14285 | 0.23 (0.09-0.44) | 0.35 (0.17-0.56) | 0.50 (0.30-0.70) | 0.58 (0.37-0.77) | 0.15 | 0.00 | 0.58 |  | HSPA6;HSPA1A;DNAJB1;BAG3;CCDC121;HSPA1B;MICB;ACTRT3;LCMT2;HSPH1 |
| unweighted_stouffer/default | ok | 14285 | 0.23 (0.09-0.44) | 0.35 (0.17-0.56) | 0.54 (0.33-0.73) | 0.58 (0.37-0.77) | 0.15 | 0.00 | 0.58 |  | HSPA1A;HSPA6;DNAJB1;HSPA1B;BAG3;CCDC121;LCMT2;MICB;ACTRT3;HSPH1 |
| fisher/default | ok | 14285 | 0.27 (0.12-0.48) | 0.35 (0.17-0.56) | 0.50 (0.30-0.70) | 0.50 (0.30-0.70) | 0.13 | 0.00 | 0.50 |  | HSPA6;CCDC121;HSPA1A;DNAJB1;BAG3;MICB;HSPA1B;ACTRT3;HSPH1;HSP90AB1 |
| rank_product_approx/default | ok | 14285 | 0.19 (0.07-0.39) | 0.31 (0.14-0.52) | 0.42 (0.23-0.63) | 0.42 (0.23-0.63) | 0.11 | 0.00 | 0.42 |  | HSPA6;HSPA1A;PCK1;HSPA1B;NPY4R;CCDC121;BAG3;DNAJB1;RBFOX3;MICB |
| maic/default | ok | 5205 | 0.00 (0.00-0.13) | 0.08 (0.01-0.25) | 0.23 (0.09-0.44) | 0.38 (0.20-0.59) | 0.10 | 0.00 | 0.00 |  | PPME1;YWHAH;TMEM18;GPRC5A;SEH1L;MICB;PWWP3A;LATS2;TLK2;DPYSL3 |
| sign_vote/default | ok | 14285 | 0.00 (0.00-0.13) | 0.00 (0.00-0.13) | 0.00 (0.00-0.13) | 0.04 (0.00-0.20) | 0.01 | 0.17 | 0.04 |  | AAAS;AAGAB;AARS2;AASDHPPT;AATF;ABCA1;ABCA2;ABCA5;ABCA7;ABCB10 |
| metavolcanor/default | ok | 14285 | 0.27 (0.12-0.48) | 0.35 (0.17-0.56) | 0.50 (0.30-0.70) | 0.50 (0.30-0.70) | 0.13 | 0.00 | 0.50 |  | HSPA6;CCDC121;HSPA1A;DNAJB1;BAG3;MICB;HSPA1B;ACTRT3;HSPH1;HSP90AB1 |
| robustrankaggreg/default | ok | 14285 | 0.23 (0.09-0.44) | 0.27 (0.12-0.48) | 0.35 (0.17-0.56) | 0.38 (0.20-0.59) | 0.10 | 0.00 | 0.38 |  | HSPA6;BAG3;HSPA1A;DEDD2;CCDC121;LAG3;PGF;HSPA1L;DNAJB1;HSPA1B |
| awfisher/default | ok | 14285 | 0.27 (0.12-0.48) | 0.31 (0.14-0.52) | 0.46 (0.27-0.67) | 0.50 (0.30-0.70) | 0.13 | 0.00 | 0.50 |  | CCDC121;HSPA6;HSPA1A;DNAJB1;BAG3;ACTRT3;MICB;HSPA1B;HSP90AB1;HSPH1 |
| metarnaseq_fisher/default | ok | 14285 | 0.27 (0.12-0.48) | 0.38 (0.20-0.59) | 0.50 (0.30-0.70) | 0.54 (0.33-0.73) | 0.14 | 0.00 | 0.54 |  | HSPA6;CCDC121;HSPA1A;DNAJB1;BAG3;MICB;HSPA1B;ACTRT3;HSPH1;HSP90AB1 |
| metarnaseq_invnorm/default | ok | 14285 | 0.23 (0.09-0.44) | 0.35 (0.17-0.56) | 0.50 (0.30-0.70) | 0.58 (0.37-0.77) | 0.15 | 0.00 | 0.58 |  | HSPA6;HSPA1A;DNAJB1;BAG3;HSPA1B;CCDC121;MICB;LCMT2;ACTRT3;HSPH1 |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
