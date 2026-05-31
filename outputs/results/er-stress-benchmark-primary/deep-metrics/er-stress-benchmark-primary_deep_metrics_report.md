# ER Stress Primary Deep Benchmark Metrics

## Point Metrics

| method | status | recall@10 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | AUROC | AUPRC enrichment | AURC | top10 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| awfisher/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.611 (0.357-0.827) | 0.110 | 0.979 | 191.808 | 0.844 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| fisher/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.611 (0.357-0.827) | 0.110 | 0.971 | 194.593 | 0.855 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| metarnaseq_fisher/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.556 (0.308-0.785) | 0.100 | 0.972 | 192.089 | 0.847 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| metarnaseq_invnorm/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.611 (0.357-0.827) | 0.110 | 0.948 | 186.818 | 0.787 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;NUCB2;CRELD2;DNAJB11;HYOU1 |
| metavolcanor/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.611 (0.357-0.827) | 0.110 | 0.971 | 194.593 | 0.855 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;DNAJB11;NUCB2;CRELD2;KLHDC7B |
| piper_deg_score/v1_2_source_unit_mean | ok | 0.278 (0.097-0.535) | 0.556 (0.308-0.785) | 0.667 (0.410-0.867) | 0.120 | 0.942 | 219.935 | 0.796 | HSPA5;KLHDC7B;OSCAR;HSP90B1;PDIA4;HERPUD1;HYOU1;CRELD2;DDIT3;CALR |
| piper_quality_weighted_score/quality_weighted_secondary | ok | 0.278 (0.097-0.535) | 0.556 (0.308-0.785) | 0.667 (0.410-0.867) | 0.120 | 0.940 | 216.465 | 0.794 | HSPA5;KLHDC7B;OSCAR;HSP90B1;PDIA4;HERPUD1;HYOU1;CRELD2;CALR;DDIT3 |
| piper_slice/locked | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.611 (0.357-0.827) | 0.110 | 0.921 | 184.609 | 0.769 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;NUCB2;DNAJB11;CRELD2;HYOU1 |
| rank_product_approx/default | ok | 0.278 (0.097-0.535) | 0.500 (0.260-0.740) | 0.556 (0.308-0.785) | 0.100 | 0.977 | 200.415 | 0.819 | HSPA5;CALR;HSP90B1;PDIA4;IGFBP1;DERL3;HHIP;SERPINA5;HERPUD1;OSCAR |
| robustrankaggreg/default | ok | 0.222 (0.064-0.476) | 0.389 (0.173-0.643) | 0.500 (0.260-0.740) | 0.090 | 0.963 | 134.924 | 0.763 | CALR;PDIA4;HSP90B1;HYOU1;HSPA5;HERPUD1;CRELD2;KLHDC7B;SLC33A1;WARS1 |
| sign_vote/default | ok | 0.000 (0.000-0.185) | 0.000 (0.000-0.185) | 0.000 (0.000-0.185) | 0.000 | 0.806 | 3.629 | 0.141 | AACS;AAK1;AAMDC;AAR2;AARS1;AASS;AATK;ABAT;ABCA3;ABCA5 |
| unweighted_stouffer/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.556 (0.308-0.785) | 0.100 | 0.943 | 194.744 | 0.781 | HSPA5;CALR;HSP90B1;PDIA4;HERPUD1;FDFT1;CRELD2;OSCAR;HYOU1;NUCB2 |
| weighted_stouffer/default | ok | 0.222 (0.064-0.476) | 0.500 (0.260-0.740) | 0.611 (0.357-0.827) | 0.110 | 0.921 | 184.609 | 0.769 | HSPA5;CALR;PDIA4;HSP90B1;HERPUD1;FDFT1;NUCB2;DNAJB11;CRELD2;HYOU1 |
| rankprod_exact/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metade/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| dexma/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metaintegrator/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| hstouffer/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| awmeta/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |

## Source-Unit Bootstrap

| method | metric | mean | 95% CI | n |
| --- | --- | ---: | --- | ---: |
| fisher | aurc_at_max_k | 0.763 | 0.499-0.878 | 50 |
| fisher | recall_at_100 | 0.551 | 0.247-0.778 | 50 |
| weighted_stouffer | aurc_at_max_k | 0.731 | 0.497-0.885 | 50 |
| weighted_stouffer | recall_at_100 | 0.543 | 0.247-0.833 | 50 |

## PIPER Advantage Metrics

| subset | n | median source units | median concordance | median LOO stability |
| --- | ---: | ---: | ---: | ---: |
| all_scored | 11910 | 3.0 | 0.8781549444490937 | 0.934005 |
| top100_primary | 100 | 3.0 | 1.0 | 0.998405 |
| top100_quality_weighted | 100 | 3.0 | 1.0 | 0.998405 |
| locked_gold | 18 | 3.0 | 1.0 | 0.999328 |
