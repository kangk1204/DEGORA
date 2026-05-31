# Heat shock deep benchmark metrics

## Point Metrics

| method | status | recall@10 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | AUROC | AUPRC enrichment | AURC | top10 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| awfisher/default | ok | 0.269 (0.116-0.478) | 0.462 (0.266-0.666) | 0.500 (0.299-0.701) | 0.130 | 0.920 | 155.539 | 0.753 | CCDC121;HSPA6;HSPA1A;DNAJB1;BAG3;ACTRT3;MICB;HSPA1B;HSP90AB1;HSPH1 |
| degora_deg_score/v1_2_source_unit_mean | ok | 0.308 (0.143-0.518) | 0.577 (0.369-0.766) | 0.731 (0.522-0.884) | 0.190 | 0.938 | 264.524 | 0.813 | HSPA6;HSPA1A;HSPA1B;DNAJB1;BAG3;HSPA1L;HSPH1;FOS;DNAJB4;MICB |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 0.308 (0.143-0.518) | 0.577 (0.369-0.766) | 0.731 (0.522-0.884) | 0.190 | 0.937 | 263.141 | 0.814 | HSPA6;HSPA1A;HSPA1B;DNAJB1;BAG3;HSPA1L;HSPH1;FOS;DNAJB4;MICB |
| degora_slice/locked | ok | 0.231 (0.090-0.436) | 0.500 (0.299-0.701) | 0.577 (0.369-0.766) | 0.150 | 0.926 | 211.769 | 0.763 | HSPA6;HSPA1A;DNAJB1;BAG3;CCDC121;HSPA1B;MICB;ACTRT3;LCMT2;HSPH1 |
| fisher/default | ok | 0.269 (0.116-0.478) | 0.500 (0.299-0.701) | 0.500 (0.299-0.701) | 0.130 | 0.913 | 182.047 | 0.752 | HSPA6;CCDC121;HSPA1A;DNAJB1;BAG3;MICB;HSPA1B;ACTRT3;HSPH1;HSP90AB1 |
| metarnaseq_fisher/default | ok | 0.269 (0.116-0.478) | 0.500 (0.299-0.701) | 0.538 (0.334-0.734) | 0.140 | 0.917 | 187.194 | 0.762 | HSPA6;CCDC121;HSPA1A;DNAJB1;BAG3;MICB;HSPA1B;ACTRT3;HSPH1;HSP90AB1 |
| metarnaseq_invnorm/default | ok | 0.231 (0.090-0.436) | 0.500 (0.299-0.701) | 0.577 (0.369-0.766) | 0.150 | 0.908 | 205.329 | 0.720 | HSPA6;HSPA1A;DNAJB1;BAG3;HSPA1B;CCDC121;MICB;LCMT2;ACTRT3;HSPH1 |
| metavolcanor/default | ok | 0.269 (0.116-0.478) | 0.500 (0.299-0.701) | 0.500 (0.299-0.701) | 0.130 | 0.913 | 182.047 | 0.752 | HSPA6;CCDC121;HSPA1A;DNAJB1;BAG3;MICB;HSPA1B;ACTRT3;HSPH1;HSP90AB1 |
| rank_product_approx/default | ok | 0.192 (0.066-0.394) | 0.423 (0.234-0.631) | 0.423 (0.234-0.631) | 0.110 | 0.924 | 152.655 | 0.691 | HSPA6;HSPA1A;PCK1;HSPA1B;NPY4R;CCDC121;BAG3;DNAJB1;RBFOX3;MICB |
| robustrankaggreg/default | ok | 0.231 (0.090-0.436) | 0.346 (0.172-0.557) | 0.385 (0.202-0.594) | 0.100 | 0.935 | 133.352 | 0.623 | HSPA6;BAG3;HSPA1A;DEDD2;CCDC121;LAG3;PGF;HSPA1L;DNAJB1;HSPA1B |
| sign_vote/default | ok | 0.000 (0.000-0.132) | 0.000 (0.000-0.132) | 0.038 (0.001-0.196) | 0.010 | 0.834 | 3.718 | 0.069 | AAAS;AAGAB;AARS2;AASDHPPT;AATF;ABCA1;ABCA2;ABCA5;ABCA7;ABCB10 |
| unweighted_stouffer/default | ok | 0.231 (0.090-0.436) | 0.538 (0.334-0.734) | 0.577 (0.369-0.766) | 0.150 | 0.926 | 216.678 | 0.754 | HSPA1A;HSPA6;DNAJB1;HSPA1B;BAG3;CCDC121;LCMT2;MICB;ACTRT3;HSPH1 |
| weighted_stouffer/default | ok | 0.231 (0.090-0.436) | 0.500 (0.299-0.701) | 0.577 (0.369-0.766) | 0.150 | 0.926 | 211.769 | 0.763 | HSPA6;HSPA1A;DNAJB1;BAG3;CCDC121;HSPA1B;MICB;ACTRT3;LCMT2;HSPH1 |
| rankprod_exact/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metade/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| dexma/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metaintegrator/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| hstouffer/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| awmeta/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |

## Source-Unit Bootstrap

| method | metric | mean | 95% CI | n |
| --- | --- | ---: | --- | ---: |
| degora_quality_weighted_score | aurc_at_max_k | 0.737 | 0.445-0.838 | 50 |
| degora_quality_weighted_score | recall_at_100 | 0.589 | 0.269-0.731 | 50 |
| fisher | aurc_at_max_k | 0.612 | 0.263-0.782 | 50 |
| fisher | recall_at_100 | 0.426 | 0.077-0.577 | 50 |
| weighted_stouffer | aurc_at_max_k | 0.631 | 0.263-0.779 | 50 |
| weighted_stouffer | recall_at_100 | 0.462 | 0.077-0.577 | 50 |

## DEGORA Advantage Metrics

| subset | n | median source units | median concordance | median LOO stability |
| --- | ---: | ---: | ---: | ---: |
| all_scored | 14285 | 3.0 | 0.9192926443009034 | 0.940427 |
| top100_primary | 100 | 3.0 | 1.0 | 0.99727 |
| top100_quality_weighted | 100 | 3.0 | 1.0 | 0.997305 |
| locked_gold | 26 | 3.0 | 1.0 | 0.998355 |
