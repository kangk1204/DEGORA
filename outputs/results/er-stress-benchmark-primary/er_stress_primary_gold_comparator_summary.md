# ER stress primary locked gold-panel comparator summary

Locked positives: `data/studies/gold/er_stress_upr_gold_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/quality_weighted_default | ok | 11910 | 0.22 (0.09-0.40) | 0.34 (0.19-0.53) | 0.50 (0.32-0.68) | 0.59 (0.41-0.76) | 0.19 | 0.00 | 0.59 |  | HSPA5;KLHDC7B;OSCAR;HSP90B1;PDIA4;HERPUD1;HYOU1;CRELD2;DDIT3;CALR |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 11910 | 0.22 (0.09-0.40) | 0.34 (0.19-0.53) | 0.50 (0.32-0.68) | 0.59 (0.41-0.76) | 0.19 | 0.00 | 0.59 |  | HSPA5;KLHDC7B;OSCAR;HSP90B1;PDIA4;HERPUD1;HYOU1;CRELD2;CALR;DDIT3 |
| weighted_stouffer/default | ok | 11910 | 0.22 (0.09-0.40) | 0.25 (0.11-0.43) | 0.47 (0.29-0.65) | 0.56 (0.38-0.74) | 0.18 | 0.00 | 0.56 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;NUCB2;DNAJB11;CRELD2;HYOU1 |
| unweighted_stouffer/default | ok | 11910 | 0.19 (0.07-0.36) | 0.25 (0.11-0.43) | 0.50 (0.32-0.68) | 0.53 (0.35-0.71) | 0.17 | 0.00 | 0.53 |  | HSPA5;CALR;HSP90B1;PDIA4;HERPUD1;FDFT1;CRELD2;OSCAR;HYOU1;NUCB2 |
| fisher/default | ok | 11910 | 0.19 (0.07-0.36) | 0.28 (0.14-0.47) | 0.47 (0.29-0.65) | 0.56 (0.38-0.74) | 0.18 | 0.00 | 0.56 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| rank_product_approx/default | ok | 11910 | 0.16 (0.05-0.33) | 0.25 (0.11-0.43) | 0.41 (0.24-0.59) | 0.50 (0.32-0.68) | 0.16 | 0.00 | 0.50 |  | HSPA5;CALR;HSP90B1;PDIA4;IGFBP1;DERL3;HHIP;SERPINA5;HERPUD1;OSCAR |
| sign_vote/default | ok | 11910 | 0.00 (0.00-0.11) | 0.00 (0.00-0.11) | 0.00 (0.00-0.11) | 0.00 (0.00-0.11) | 0.00 | 1.00 | 0.00 |  | AACS;AAK1;AAMDC;AAR2;AARS1;AASS;AATK;ABAT;ABCA3;ABCA5 |
| metavolcanor/default | ok | 11910 | 0.19 (0.07-0.36) | 0.28 (0.14-0.47) | 0.47 (0.29-0.65) | 0.56 (0.38-0.74) | 0.18 | 0.00 | 0.56 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| robustrankaggreg/default | ok | 11910 | 0.19 (0.07-0.36) | 0.19 (0.07-0.36) | 0.31 (0.16-0.50) | 0.44 (0.26-0.62) | 0.14 | 0.00 | 0.44 |  | CALR;PDIA4;HSP90B1;HYOU1;HSPA5;HERPUD1;CRELD2;KLHDC7B;SLC33A1;WARS1 |
| awfisher/default | ok | 11910 | 0.19 (0.07-0.36) | 0.28 (0.14-0.47) | 0.47 (0.29-0.65) | 0.56 (0.38-0.74) | 0.18 | 0.00 | 0.56 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| metarnaseq_fisher/default | ok | 11910 | 0.19 (0.07-0.36) | 0.28 (0.14-0.47) | 0.47 (0.29-0.65) | 0.53 (0.35-0.71) | 0.17 | 0.00 | 0.53 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| metarnaseq_invnorm/default | ok | 11910 | 0.22 (0.09-0.40) | 0.25 (0.11-0.43) | 0.47 (0.29-0.65) | 0.56 (0.38-0.74) | 0.18 | 0.00 | 0.56 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;NUCB2;CRELD2;DNAJB11;HYOU1 |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |
| piper_slice/locked | ok | 11910 | 0.22 (0.09-0.40) | 0.25 (0.11-0.43) | 0.47 (0.29-0.65) | 0.56 (0.38-0.74) | 0.18 | 0.00 | 0.56 |  | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;NUCB2;DNAJB11;CRELD2;HYOU1 |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
