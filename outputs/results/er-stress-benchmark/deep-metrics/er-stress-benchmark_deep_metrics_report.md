# ER Stress Full Deep Benchmark Metrics

## Point Metrics

| method | status | recall@10 | recall@50 | recall@100 | AURC | top10 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| awfisher/default | ok | 0.222 | 0.556 | 0.611 | 0.787 | UNC5B;CRELD2;ADAM19;LOXL4;HSP90B1;PPP1R15A;KLHDC8A;SLC1A4;SEL1L;PDIA4 |
| fisher/default | ok | 0.222 | 0.611 | 0.611 | 0.799 | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| metarnaseq_fisher/default | ok | 0.222 | 0.611 | 0.611 | 0.803 | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| metarnaseq_invnorm/default | ok | 0.000 | 0.000 | 0.000 | 0.019 | 0610007P14RIK;0610009B22RIK;0610009O20RIK;0610010F05RIK;0610010K14RIK;0610012G03RIK;0610030E20RIK;0610037L13RIK;1110002L01RIK;1110004E09RIK |
| metavolcanor/default | ok | 0.222 | 0.611 | 0.611 | 0.799 | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| piper_deg_score/v1_2_source_unit_mean | ok | 0.111 | 0.111 | 0.167 | 0.388 | ADAM19;DDIT3;HAS2;PPP1R15A;ACTB;DUSP6;LFNG;LMNB1;MGP;SNAI1 |
| piper_quality_weighted_score/quality_weighted_secondary | ok | 0.278 | 0.611 | 0.611 | 0.811 | DDIT3;PPP1R15A;HAS2;ADAM19;PTX3;LFNG;ATF3;HSPA5;EVA1A;CHAC1 |
| piper_slice/locked | ok | 0.056 | 0.111 | 0.167 | 0.437 | ADAM19;PPP1R15A;ACTB;DUSP6;B4GALT5;FSCN1;POMP;LMNB1;LFNG;SNAI1 |
| rank_product_approx/default | ok | 0.000 | 0.278 | 0.333 | 0.687 | UNC5B;ADAM19;DNAJC15;KLHDC8A;CRELD2;GM6548;PTGS1;ZFP57;MAFB;TXNDC5 |
| robustrankaggreg/default | ok | 0.167 | 0.500 | 0.500 | 0.641 | ADAM19;CRELD2;UNC5B;KLHDC8A;DDIT3;GAS7;DUSP6;HAS2;HSP90B1;HSPA5 |
| sign_vote/default | ok | 0.000 | 0.000 | 0.000 | 0.021 | AAGAB;AAMDC;AAMP;AASDHPPT;ABCC1;ABCD1;ABCD4;ABCE1;ABHD11;ABHD14A |
| unweighted_stouffer/default | ok | 0.056 | 0.111 | 0.111 | 0.379 | ADAM19;ACTB;PPP1R15A;B4GALT5;DUSP6;FSCN1;LFNG;POMP;LMNB1;TUBB6 |
| weighted_stouffer/default | ok | 0.056 | 0.111 | 0.167 | 0.437 | ADAM19;PPP1R15A;ACTB;DUSP6;B4GALT5;FSCN1;POMP;LMNB1;LFNG;SNAI1 |
| rankprod_exact/default | blocked |  |  |  |  |  |
| metade/default | blocked |  |  |  |  |  |
| dexma/default | blocked |  |  |  |  |  |
| metaintegrator/default | blocked |  |  |  |  |  |
| hstouffer/default | blocked |  |  |  |  |  |
| awmeta/default | blocked |  |  |  |  |  |

## Source-Unit Bootstrap

| method | metric | mean | 95% CI | n |
| --- | --- | ---: | --- | ---: |
| fisher | aurc_at_max_k | 0.618 | 0.000-0.813 | 50 |
| fisher | recall_at_100 | 0.372 | 0.000-0.654 | 50 |
| weighted_stouffer | aurc_at_max_k | 0.375 | 0.000-0.840 | 50 |
| weighted_stouffer | recall_at_100 | 0.213 | 0.000-0.654 | 50 |

## PIPER Advantage Metrics

| subset | n | median source units | median concordance | median LOO stability |
| --- | ---: | ---: | ---: | ---: |
| all_scored | 10872 | 3.0 | 0.8907907751451183 | 0.831954 |
| top100_primary | 100 | 3.0 | 1.0 | 0.9844094999999999 |
| top100_quality_weighted | 100 | 3.0 | 1.0 | 0.9862029999999999 |
| locked_gold | 18 | 3.0 | 0.819846437707293 | 0.7413540000000001 |
