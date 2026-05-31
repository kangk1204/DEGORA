# ER stress locked gold-panel comparator summary

Locked positives: `data/studies/gold/er_stress_upr_gold_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/v1_2_source_unit_mean | ok | 14239 | 0.19 (0.07-0.36) | 0.41 (0.24-0.59) | 0.59 (0.41-0.76) | 0.78 (0.60-0.91) | 0.25 | 0.00 | 0.78 |  | DDIT3;HSPA5;TRIB3;HERPUD1;SLC6A9;CCPG1;DNAJB9;SEL1L;CTH;SLC3A2 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 14239 | 0.25 (0.11-0.43) | 0.41 (0.24-0.59) | 0.59 (0.41-0.76) | 0.78 (0.60-0.91) | 0.25 | 0.00 | 0.78 |  | DDIT3;HSPA5;HERPUD1;HYOU1;TRIB3;PDIA4;SEL1L;SLC6A9;MANF;CCPG1 |
| degora_slice/locked | ok | 14239 | 0.16 (0.05-0.33) | 0.25 (0.11-0.43) | 0.56 (0.38-0.74) | 0.75 (0.57-0.89) | 0.24 | 0.00 | 0.75 |  | TRIB3;DDIT3;SEL1L;SLC3A2;CEBPB;NFE2L1;ASNS;PCK2;CTH;HSPA5 |
| weighted_stouffer/default | ok | 14239 | 0.16 (0.05-0.33) | 0.25 (0.11-0.43) | 0.56 (0.38-0.74) | 0.75 (0.57-0.89) | 0.24 | 0.00 | 0.75 |  | TRIB3;DDIT3;SEL1L;SLC3A2;CEBPB;NFE2L1;ASNS;PCK2;CTH;HSPA5 |
| unweighted_stouffer/default | ok | 14239 | 0.16 (0.05-0.33) | 0.31 (0.16-0.50) | 0.56 (0.38-0.74) | 0.75 (0.57-0.89) | 0.24 | 0.00 | 0.75 |  | TRIB3;DDIT3;SEL1L;SLC3A2;HSPA5;PCK2;ASNS;NFE2L1;CEBPB;CTH |
| fisher/default | ok | 14239 | 0.12 (0.04-0.29) | 0.22 (0.09-0.40) | 0.47 (0.29-0.65) | 0.66 (0.47-0.81) | 0.21 | 0.00 | 0.66 |  | TRIB3;CEBPB;NFE2L1;ASNS;CTH;DDIT3;SEL1L;SLC3A2;PCK2;SLC38A2 |
| rank_product_approx/default | ok | 14239 | 0.22 (0.09-0.40) | 0.34 (0.19-0.53) | 0.50 (0.32-0.68) | 0.69 (0.50-0.84) | 0.22 | 0.00 | 0.69 |  | HSPA5;TRIB3;CPB2;SLC3A2;HERPUD1;CALR;HSP90B1;PDIA4;SEL1L;DNAJB11 |
| sign_vote/default | ok | 14239 | 0.00 (0.00-0.11) | 0.00 (0.00-0.11) | 0.00 (0.00-0.11) | 0.00 (0.00-0.11) | 0.00 | 1.00 | 0.00 |  | AACS;AAK1;AAR2;ABCA3;ABCA5;ABCB7;ABCD3;ABHD11;ABHD17B;ABHD2 |
| metavolcanor/default | ok | 14239 | 0.12 (0.04-0.29) | 0.22 (0.09-0.40) | 0.47 (0.29-0.65) | 0.66 (0.47-0.81) | 0.21 | 0.00 | 0.66 |  | TRIB3;CEBPB;NFE2L1;ASNS;CTH;DDIT3;SEL1L;SLC3A2;PCK2;SLC38A2 |
| robustrankaggreg/default | ok | 14239 | 0.19 (0.07-0.36) | 0.31 (0.16-0.50) | 0.50 (0.32-0.68) | 0.69 (0.50-0.84) | 0.22 | 0.00 | 0.69 |  | HSPA5;CALR;HERPUD1;HYOU1;SLC33A1;FDFT1;PDIA4;CRELD2;TTC17;HSP90B1 |
| awfisher/default | ok | 14239 | 0.12 (0.04-0.29) | 0.22 (0.09-0.40) | 0.47 (0.29-0.65) | 0.66 (0.47-0.81) | 0.21 | 0.00 | 0.66 |  | TRIB3;CEBPB;NFE2L1;CTH;ASNS;DDIT3;SEL1L;SLC3A2;SLC38A2;PCK2 |
| metarnaseq_fisher/default | ok | 14239 | 0.12 (0.04-0.29) | 0.19 (0.07-0.36) | 0.47 (0.29-0.65) | 0.66 (0.47-0.81) | 0.21 | 0.00 | 0.66 |  | TRIB3;CEBPB;NFE2L1;ASNS;CTH;DDIT3;SEL1L;SLC3A2;PCK2;SLC38A2 |
| metarnaseq_invnorm/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metarnaseq_invnorm_uninformative_sparse_public_summary |  |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
