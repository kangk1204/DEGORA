# ER-stress-full-gold-comparator-summary

Locked positives: `data/studies/gold/er_stress_upr_gold_panel.csv`.

| method | status | rows | recall@10 | recall@20 | recall@50 | recall@100 | dir@50 | dir@100 | failure | top10 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| piper_deg_score/quality_weighted_default | ok | 10872 | 0.11 | 0.11 | 0.11 | 0.17 | 0.11 | 0.17 |  | ADAM19;DDIT3;HAS2;PPP1R15A;ACTB;DUSP6;LFNG;LMNB1;MGP;SNAI1 |
| piper_quality_weighted_score/quality_weighted_secondary | ok | 10872 | 0.28 | 0.44 | 0.61 | 0.61 | 0.61 | 0.61 |  | DDIT3;PPP1R15A;HAS2;ADAM19;PTX3;LFNG;ATF3;HSPA5;EVA1A;CHAC1 |
| piper_slice/locked | ok | 10872 | 0.06 | 0.11 | 0.11 | 0.17 | 0.11 | 0.17 |  | ADAM19;PPP1R15A;ACTB;DUSP6;B4GALT5;FSCN1;POMP;LMNB1;LFNG;SNAI1 |
| weighted_stouffer/default | ok | 10872 | 0.06 | 0.11 | 0.11 | 0.17 | 0.11 | 0.17 |  | ADAM19;PPP1R15A;ACTB;DUSP6;B4GALT5;FSCN1;POMP;LMNB1;LFNG;SNAI1 |
| unweighted_stouffer/default | ok | 10872 | 0.06 | 0.06 | 0.11 | 0.11 | 0.11 | 0.11 |  | ADAM19;ACTB;PPP1R15A;B4GALT5;DUSP6;FSCN1;LFNG;POMP;LMNB1;TUBB6 |
| fisher/default | ok | 10872 | 0.22 | 0.50 | 0.61 | 0.61 | 0.56 | 0.56 |  | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| rank_product_approx/default | ok | 10872 | 0.00 | 0.06 | 0.28 | 0.33 | 0.28 | 0.33 |  | UNC5B;ADAM19;DNAJC15;KLHDC8A;CRELD2;GM6548;PTGS1;ZFP57;MAFB;TXNDC5 |
| sign_vote/default | ok | 10872 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |  | AAGAB;AAMDC;AAMP;AASDHPPT;ABCC1;ABCD1;ABCD4;ABCE1;ABHD11;ABHD14A |
| metavolcanor/default | ok | 10872 | 0.22 | 0.50 | 0.61 | 0.61 | 0.56 | 0.56 |  | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| robustrankaggreg/default | ok | 10872 | 0.17 | 0.22 | 0.50 | 0.50 | 0.50 | 0.50 |  | ADAM19;CRELD2;UNC5B;KLHDC8A;DDIT3;GAS7;DUSP6;HAS2;HSP90B1;HSPA5 |
| awfisher/default | ok | 10872 | 0.22 | 0.50 | 0.56 | 0.61 | 0.56 | 0.61 |  | UNC5B;CRELD2;ADAM19;LOXL4;HSP90B1;PPP1R15A;KLHDC8A;SLC1A4;SEL1L;PDIA4 |
| metarnaseq_fisher/default | ok | 10872 | 0.22 | 0.50 | 0.61 | 0.61 | 0.56 | 0.56 |  | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| metarnaseq_invnorm/default | ok | 10872 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |  | 0610007P14RIK;0610009B22RIK;0610009O20RIK;0610010F05RIK;0610010K14RIK;0610012G03RIK;0610030E20RIK;0610037L13RIK;1110002L01RIK;1110004E09RIK |
| hstouffer/default | blocked | 0 |  |  |  |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  |  |  |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  |  |  |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  |  |  |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  |  |  |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  |  |  |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
