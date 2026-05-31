# Hypoxia/HIF locked gold-panel comparator summary

Locked positives: `data/studies/gold/hypoxia_hif1_gold_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/v1_2_source_unit_mean | ok | 32687 | 0.18 (0.07-0.35) | 0.39 (0.23-0.58) | 0.76 (0.58-0.89) | 0.85 (0.68-0.95) | 0.28 | 0.00 | 0.85 |  | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;ADM;ANKRD37;AK4;HK2;PPFIA4 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 32687 | 0.21 (0.09-0.39) | 0.36 (0.20-0.55) | 0.76 (0.58-0.89) | 0.85 (0.68-0.95) | 0.28 | 0.00 | 0.85 |  | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;AK4;ANKRD37;HK2;ADM;VEGFA |
| degora_slice/locked | ok | 32687 | 0.24 (0.11-0.42) | 0.48 (0.31-0.66) | 0.82 (0.65-0.93) | 0.85 (0.68-0.95) | 0.28 | 0.00 | 0.85 |  | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| weighted_stouffer/default | ok | 32687 | 0.24 (0.11-0.42) | 0.48 (0.31-0.66) | 0.82 (0.65-0.93) | 0.85 (0.68-0.95) | 0.28 | 0.00 | 0.85 |  | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| unweighted_stouffer/default | ok | 32687 | 0.27 (0.13-0.46) | 0.48 (0.31-0.66) | 0.82 (0.65-0.93) | 0.85 (0.68-0.95) | 0.28 | 0.00 | 0.85 |  | ADM;AK4;ALDOA;BNIP3L;HK2;MXI1;NDRG1;P4HA1;PFKFB4;PGK1 |
| fisher/default | ok | 32687 | 0.21 (0.09-0.39) | 0.33 (0.18-0.52) | 0.73 (0.54-0.87) | 0.82 (0.65-0.93) | 0.27 | 0.00 | 0.82 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| rank_product_approx/default | ok | 32687 | 0.18 (0.07-0.35) | 0.42 (0.25-0.61) | 0.73 (0.54-0.87) | 0.82 (0.65-0.93) | 0.27 | 0.00 | 0.82 |  | BNIP3L;TMEM45A;NDRG1;P4HA1;VEGFA;SLC2A1;AK4;MIR210HG;PFKFB4;HIF1A-AS3 |
| sign_vote/default | ok | 32687 | 0.03 (0.00-0.16) | 0.06 (0.01-0.20) | 0.12 (0.03-0.28) | 0.27 (0.13-0.46) | 0.09 | 0.00 | 0.27 |  | AEN;AIFM1;AK4;ALDOA;APEH;APOO;ATAD3A;ATAD3B;ATIC;ATP6V0A2 |
| metavolcanor/default | ok | 32687 | 0.21 (0.09-0.39) | 0.33 (0.18-0.52) | 0.73 (0.54-0.87) | 0.82 (0.65-0.93) | 0.27 | 0.00 | 0.82 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| robustrankaggreg/default | ok | 32687 | 0.24 (0.11-0.42) | 0.48 (0.31-0.66) | 0.79 (0.61-0.91) | 0.82 (0.65-0.93) | 0.27 | 0.00 | 0.82 |  | AK4;BNIP3L;ALDOA;P4HA1;MXI1;NDRG1;KDM4B;ADM;TMEM45A;KDM3A |
| awfisher/default | ok | 32687 | 0.21 (0.09-0.39) | 0.33 (0.18-0.52) | 0.73 (0.54-0.87) | 0.82 (0.65-0.93) | 0.27 | 0.00 | 0.82 |  | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| metarnaseq_fisher/default | ok | 32687 | 0.30 (0.16-0.49) | 0.52 (0.34-0.69) | 0.73 (0.54-0.87) | 0.82 (0.65-0.93) | 0.27 | 0.00 | 0.82 |  | VEGFA;BNIP3L;PGK1;P4HA1;ADM;HK2;PFKFB4;ALDOA;NDRG1;BHLHE40 |
| metarnaseq_invnorm/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metarnaseq_invnorm_uninformative_sparse_public_summary |  |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
