# ER stress cross-platform deep benchmark metrics

## Point Metrics

| method | status | recall@10 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | AUROC | AUPRC enrichment | AURC | top10 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| awfisher/default | ok | 0.125 (0.035-0.290) | 0.469 (0.291-0.653) | 0.656 (0.468-0.814) | 0.210 | 0.976 | 114.672 | 0.818 | TRIB3;CEBPB;NFE2L1;CTH;ASNS;DDIT3;SEL1L;SLC3A2;SLC38A2;PCK2 |
| degora_deg_score/v1_2_source_unit_mean | ok | 0.188 (0.072-0.364) | 0.594 (0.406-0.763) | 0.781 (0.600-0.907) | 0.250 | 0.978 | 218.157 | 0.840 | DDIT3;HSPA5;TRIB3;HERPUD1;SLC6A9;CCPG1;DNAJB9;SEL1L;CTH;SLC3A2 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 0.250 (0.115-0.434) | 0.594 (0.406-0.763) | 0.781 (0.600-0.907) | 0.250 | 0.972 | 237.292 | 0.826 | DDIT3;HSPA5;HERPUD1;HYOU1;TRIB3;PDIA4;SEL1L;SLC6A9;MANF;CCPG1 |
| degora_slice/locked | ok | 0.156 (0.053-0.328) | 0.562 (0.377-0.736) | 0.750 (0.566-0.885) | 0.240 | 0.978 | 163.813 | 0.813 | TRIB3;DDIT3;SEL1L;SLC3A2;CEBPB;NFE2L1;ASNS;PCK2;CTH;HSPA5 |
| fisher/default | ok | 0.125 (0.035-0.290) | 0.469 (0.291-0.653) | 0.656 (0.468-0.814) | 0.210 | 0.977 | 118.746 | 0.823 | TRIB3;CEBPB;NFE2L1;ASNS;CTH;DDIT3;SEL1L;SLC3A2;PCK2;SLC38A2 |
| metarnaseq_fisher/default | ok | 0.125 (0.035-0.290) | 0.469 (0.291-0.653) | 0.656 (0.468-0.814) | 0.210 | 0.978 | 117.859 | 0.822 | TRIB3;CEBPB;NFE2L1;ASNS;CTH;DDIT3;SEL1L;SLC3A2;PCK2;SLC38A2 |
| metavolcanor/default | ok | 0.125 (0.035-0.290) | 0.469 (0.291-0.653) | 0.656 (0.468-0.814) | 0.210 | 0.977 | 118.746 | 0.823 | TRIB3;CEBPB;NFE2L1;ASNS;CTH;DDIT3;SEL1L;SLC3A2;PCK2;SLC38A2 |
| rank_product_approx/default | ok | 0.219 (0.093-0.400) | 0.500 (0.319-0.681) | 0.688 (0.500-0.839) | 0.220 | 0.967 | 174.776 | 0.818 | HSPA5;TRIB3;CPB2;SLC3A2;HERPUD1;CALR;HSP90B1;PDIA4;SEL1L;DNAJB11 |
| robustrankaggreg/default | ok | 0.188 (0.072-0.364) | 0.500 (0.319-0.681) | 0.688 (0.500-0.839) | 0.220 | 0.945 | 151.420 | 0.802 | HSPA5;CALR;HERPUD1;HYOU1;SLC33A1;FDFT1;PDIA4;CRELD2;TTC17;HSP90B1 |
| sign_vote/default | ok | 0.000 (0.000-0.109) | 0.000 (0.000-0.109) | 0.000 (0.000-0.109) | 0.000 | 0.855 | 5.219 | 0.204 | AACS;AAK1;AAR2;ABCA3;ABCA5;ABCB7;ABCD3;ABHD11;ABHD17B;ABHD2 |
| unweighted_stouffer/default | ok | 0.156 (0.053-0.328) | 0.562 (0.377-0.736) | 0.750 (0.566-0.885) | 0.240 | 0.978 | 180.971 | 0.819 | TRIB3;DDIT3;SEL1L;SLC3A2;HSPA5;PCK2;ASNS;NFE2L1;CEBPB;CTH |
| weighted_stouffer/default | ok | 0.156 (0.053-0.328) | 0.562 (0.377-0.736) | 0.750 (0.566-0.885) | 0.240 | 0.978 | 163.813 | 0.813 | TRIB3;DDIT3;SEL1L;SLC3A2;CEBPB;NFE2L1;ASNS;PCK2;CTH;HSPA5 |
| metarnaseq_invnorm/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| rankprod_exact/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metade/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| dexma/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metaintegrator/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| hstouffer/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| awmeta/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |

## Source-Unit Bootstrap

| method | metric | mean | 95% CI | n |
| --- | --- | ---: | --- | ---: |
| degora_quality_weighted_score | aurc_at_max_k | 0.795 | 0.592-0.891 | 50 |
| degora_quality_weighted_score | recall_at_100 | 0.699 | 0.469-0.844 | 50 |
| fisher | aurc_at_max_k | 0.775 | 0.527-0.831 | 50 |
| fisher | recall_at_100 | 0.589 | 0.341-0.688 | 50 |
| weighted_stouffer | aurc_at_max_k | 0.772 | 0.535-0.862 | 50 |
| weighted_stouffer | recall_at_100 | 0.628 | 0.372-0.743 | 50 |

## DEGORA Advantage Metrics

| subset | n | median source units | median concordance | median LOO stability |
| --- | ---: | ---: | ---: | ---: |
| all_scored | 14239 | 4.0 | 0.889477699645202 | 0.964885 |
| top100_primary | 100 | 4.0 | 1.0 | 0.9995615 |
| top100_quality_weighted | 100 | 4.0 | 1.0 | 0.9995259999999999 |
| locked_gold | 32 | 4.0 | 1.0 | 0.999719 |
