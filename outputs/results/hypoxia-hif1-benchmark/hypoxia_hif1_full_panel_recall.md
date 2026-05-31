# Hypoxia/HIF FULL-EXPANDED supplementary panel recall

Locked positives: `data/studies/gold/hypoxia_hif1_gold_panel_full.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/v1_2_source_unit_mean | ok | 32687 | 0.19 (0.09-0.33) | 0.36 (0.23-0.51) | 0.64 (0.49-0.77) | 0.85 (0.72-0.94) | 0.40 | 0.00 | 0.85 |  | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;ADM;ANKRD37;AK4;HK2;PPFIA4 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 32687 | 0.21 (0.11-0.36) | 0.32 (0.19-0.47) | 0.64 (0.49-0.77) | 0.85 (0.72-0.94) | 0.40 | 0.00 | 0.85 |  | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;AK4;ANKRD37;HK2;ADM;VEGFA |
| degora_slice/locked | ok | 32687 | 0.21 (0.11-0.36) | 0.38 (0.25-0.54) | 0.74 (0.60-0.86) | 0.83 (0.69-0.92) | 0.39 | 0.00 | 0.83 |  | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| weighted_stouffer/default | ok | 32687 | 0.21 (0.11-0.36) | 0.38 (0.25-0.54) | 0.74 (0.60-0.86) | 0.83 (0.69-0.92) | 0.39 | 0.00 | 0.83 |  | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| unweighted_stouffer/default | ok | 32687 | 0.21 (0.11-0.36) | 0.38 (0.25-0.54) | 0.74 (0.60-0.86) | 0.83 (0.69-0.92) | 0.39 | 0.00 | 0.83 |  | ADM;AK4;ALDOA;BNIP3L;HK2;MXI1;NDRG1;P4HA1;PFKFB4;PGK1 |
| fisher/default | ok | 32687 | 0.19 (0.09-0.33) | 0.32 (0.19-0.47) | 0.64 (0.49-0.77) | 0.81 (0.67-0.91) | 0.38 | 0.00 | 0.81 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| rank_product_approx/default | ok | 32687 | 0.17 (0.08-0.31) | 0.36 (0.23-0.51) | 0.62 (0.46-0.75) | 0.79 (0.64-0.89) | 0.37 | 0.00 | 0.79 |  | BNIP3L;TMEM45A;NDRG1;P4HA1;VEGFA;SLC2A1;AK4;MIR210HG;PFKFB4;HIF1A-AS3 |
| sign_vote/default | ok | 32687 | 0.04 (0.01-0.15) | 0.06 (0.01-0.18) | 0.15 (0.06-0.28) | 0.32 (0.19-0.47) | 0.15 | 0.00 | 0.32 |  | AEN;AIFM1;AK4;ALDOA;APEH;APOO;ATAD3A;ATAD3B;ATIC;ATP6V0A2 |
| metavolcanor/default | ok | 32687 | 0.19 (0.09-0.33) | 0.32 (0.19-0.47) | 0.64 (0.49-0.77) | 0.81 (0.67-0.91) | 0.38 | 0.00 | 0.81 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| robustrankaggreg/default | ok | 32687 | 0.21 (0.11-0.36) | 0.40 (0.26-0.56) | 0.70 (0.55-0.83) | 0.77 (0.62-0.88) | 0.36 | 0.00 | 0.77 |  | AK4;BNIP3L;ALDOA;P4HA1;MXI1;NDRG1;KDM4B;ADM;TMEM45A;KDM3A |
| awfisher/default | ok | 32687 | 0.19 (0.09-0.33) | 0.32 (0.19-0.47) | 0.64 (0.49-0.77) | 0.81 (0.67-0.91) | 0.38 | 0.00 | 0.81 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| metarnaseq_fisher/default | ok | 32687 | 0.21 (0.11-0.36) | 0.38 (0.25-0.54) | 0.64 (0.49-0.77) | 0.81 (0.67-0.91) | 0.38 | 0.00 | 0.81 |  | VEGFA;BNIP3L;PGK1;P4HA1;ADM;HK2;PFKFB4;ALDOA;NDRG1;BHLHE40 |
| metarnaseq_invnorm/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metarnaseq_invnorm_uninformative_sparse_public_summary |  |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
