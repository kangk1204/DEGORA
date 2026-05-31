# IFN FULL-EXPANDED supplementary panel recall

Locked positives: `data/studies/gold/ifn_gold_panel_full.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/v1_2_source_unit_mean | ok | 9860 | 0.18 (0.09-0.31) | 0.30 (0.18-0.45) | 0.70 (0.55-0.82) | 0.92 (0.81-0.98) | 0.46 | 0.00 | 0.92 |  | RSAD2;CMPK2;MX1;TAP1;IFIT1;IFIT3;DDX60;OASL;IFIH1;OAS1 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 9860 | 0.18 (0.09-0.31) | 0.30 (0.18-0.45) | 0.68 (0.53-0.80) | 0.92 (0.81-0.98) | 0.46 | 0.00 | 0.92 |  | RSAD2;CMPK2;MX1;TAP1;IFIT1;IFIT3;DDX60;OASL;IFIH1;OAS1 |
| degora_slice/locked | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.60 (0.45-0.74) | 0.86 (0.73-0.94) | 0.43 | 0.00 | 0.86 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;IFIT3;SAMD9;PLSCR1 |
| weighted_stouffer/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.60 (0.45-0.74) | 0.86 (0.73-0.94) | 0.43 | 0.00 | 0.86 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;IFIT3;SAMD9;PLSCR1 |
| unweighted_stouffer/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.60 (0.45-0.74) | 0.86 (0.73-0.94) | 0.43 | 0.00 | 0.86 |  | TAP1;MX1;RSAD2;CMPK2;TRIM21;DDX60;IFIT1;SAMD9;HERC5;IFIT3 |
| fisher/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.62 (0.47-0.75) | 0.82 (0.69-0.91) | 0.41 | 0.00 | 0.82 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| rank_product_approx/default | ok | 9860 | 0.12 (0.05-0.24) | 0.28 (0.16-0.42) | 0.64 (0.49-0.77) | 0.82 (0.69-0.91) | 0.41 | 0.00 | 0.82 |  | TAP1;MX1;CMPK2;TRIM21;RSAD2;IFIT1;APOL6;HERC5;IFIT3;GBP1 |
| sign_vote/default | ok | 9860 | 0.00 (0.00-0.07) | 0.00 (0.00-0.07) | 0.00 (0.00-0.07) | 0.02 (0.00-0.11) | 0.01 | 0.40 | 0.02 |  | AADAT;AAK1;AARS2;AASDH;AASDHPPT;AATK;ABCB11;ABCB6;ABCB8;ABCC2 |
| metavolcanor/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.62 (0.47-0.75) | 0.82 (0.69-0.91) | 0.41 | 0.00 | 0.82 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| robustrankaggreg/default | ok | 9860 | 0.08 (0.02-0.19) | 0.22 (0.12-0.36) | 0.52 (0.37-0.66) | 0.70 (0.55-0.82) | 0.35 | 0.00 | 0.70 |  | TAP1;TRIM21;MAK;MT1M;TRIM69;HERC5;APOL6;CMPK2;IFIT1;MX1 |
| awfisher/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.58 (0.43-0.72) | 0.80 (0.66-0.90) | 0.40 | 0.00 | 0.80 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;MX2;DDX60;IFIT3;PLSCR1;IFIT1 |
| metarnaseq_fisher/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.62 (0.47-0.75) | 0.82 (0.69-0.91) | 0.41 | 0.00 | 0.82 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| metarnaseq_invnorm/default | ok | 9860 | 0.16 (0.07-0.29) | 0.30 (0.18-0.45) | 0.60 (0.45-0.74) | 0.86 (0.73-0.94) | 0.43 | 0.00 | 0.86 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;SAMD9;IFIT3;PLSCR1 |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
