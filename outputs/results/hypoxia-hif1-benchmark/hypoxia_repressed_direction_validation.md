# Hypoxia repressed-marker direction validation (RNA-seq)

Locked positives: `data/studies/gold/hypoxia_hif1_repressed_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/v1_2_source_unit_mean | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;ADM;ANKRD37;AK4;HK2;PPFIA4 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;AK4;ANKRD37;HK2;ADM;VEGFA |
| degora_slice/locked | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| weighted_stouffer/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| unweighted_stouffer/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | ADM;AK4;ALDOA;BNIP3L;HK2;MXI1;NDRG1;P4HA1;PFKFB4;PGK1 |
| fisher/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| rank_product_approx/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | BNIP3L;TMEM45A;NDRG1;P4HA1;VEGFA;SLC2A1;AK4;MIR210HG;PFKFB4;HIF1A-AS3 |
| sign_vote/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | AEN;AIFM1;AK4;ALDOA;APEH;APOO;ATAD3A;ATAD3B;ATIC;ATP6V0A2 |
| metavolcanor/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| robustrankaggreg/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | AK4;BNIP3L;ALDOA;P4HA1;MXI1;NDRG1;KDM4B;ADM;TMEM45A;KDM3A |
| awfisher/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| metarnaseq_fisher/default | ok | 32687 | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 (0.00-0.21) | 0.00 | 1.00 | 0.00 |  | VEGFA;BNIP3L;PGK1;P4HA1;ADM;HK2;PFKFB4;ALDOA;NDRG1;BHLHE40 |
| metarnaseq_invnorm/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metarnaseq_invnorm_uninformative_sparse_public_summary |  |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
