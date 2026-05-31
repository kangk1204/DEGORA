# IFN locked gold-panel comparator summary

Locked positives: `data/studies/gold/ifn_gold_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/quality_weighted_default | ok | 9860 | 0.29 (0.13-0.49) | 0.43 (0.24-0.63) | 0.82 (0.63-0.94) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | RSAD2;CMPK2;MX1;TAP1;IFIT1;IFIT3;DDX60;OASL;IFIH1;OAS1 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 9860 | 0.29 (0.13-0.49) | 0.43 (0.24-0.63) | 0.82 (0.63-0.94) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | RSAD2;CMPK2;MX1;TAP1;IFIT1;IFIT3;DDX60;OASL;IFIH1;OAS1 |
| degora_slice/locked | ok | 9860 | 0.18 (0.06-0.37) | 0.36 (0.19-0.56) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;IFIT3;SAMD9;PLSCR1 |
| weighted_stouffer/default | ok | 9860 | 0.18 (0.06-0.37) | 0.36 (0.19-0.56) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;IFIT3;SAMD9;PLSCR1 |
| unweighted_stouffer/default | ok | 9860 | 0.21 (0.08-0.41) | 0.36 (0.19-0.56) | 0.75 (0.55-0.89) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;MX1;RSAD2;CMPK2;TRIM21;DDX60;IFIT1;SAMD9;HERC5;IFIT3 |
| fisher/default | ok | 9860 | 0.18 (0.06-0.37) | 0.39 (0.22-0.59) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| rank_product_approx/default | ok | 9860 | 0.21 (0.08-0.41) | 0.43 (0.24-0.63) | 0.79 (0.59-0.92) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | TAP1;MX1;CMPK2;TRIM21;RSAD2;IFIT1;APOL6;HERC5;IFIT3;GBP1 |
| maic/default | ok | 2799 | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 | 1.00 | 0.00 |  | ADAM15;HLA-F;CCDC61;FKBP4;RTP4;ETS1;RPS6KA4;MAN1B1;AK2;INTS3 |
| sign_vote/default | ok | 9860 | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 | 1.00 | 0.00 |  | AADAT;AAK1;AARS2;AASDH;AASDHPPT;AATK;ABCB11;ABCB6;ABCB8;ABCC2 |
| metavolcanor/default | ok | 9860 | 0.18 (0.06-0.37) | 0.39 (0.22-0.59) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| robustrankaggreg/default | ok | 9860 | 0.14 (0.04-0.33) | 0.36 (0.19-0.56) | 0.71 (0.51-0.87) | 0.82 (0.63-0.94) | 0.23 | 0.00 | 0.82 |  | TAP1;TRIM21;MAK;MT1M;TRIM69;HERC5;APOL6;CMPK2;IFIT1;MX1 |
| awfisher/default | ok | 9860 | 0.21 (0.08-0.41) | 0.39 (0.22-0.59) | 0.75 (0.55-0.89) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;MX2;DDX60;IFIT3;PLSCR1;IFIT1 |
| metarnaseq_fisher/default | ok | 9860 | 0.18 (0.06-0.37) | 0.39 (0.22-0.59) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| metarnaseq_invnorm/default | ok | 9860 | 0.18 (0.06-0.37) | 0.36 (0.19-0.56) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;SAMD9;IFIT3;PLSCR1 |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
