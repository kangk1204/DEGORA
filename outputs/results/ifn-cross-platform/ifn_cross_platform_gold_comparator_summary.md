# IFN RNA-seq+microarray comparator summary

Locked positives: `data/studies/gold/ifn_gold_panel.csv`.

| method | status | rows | recall@10 (95% CI) | recall@20 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | FDR@100 | dir@100 | failure | top10 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| degora_deg_score/quality_weighted_default | ok | 12777 | 0.36 (0.19-0.56) | 0.64 (0.44-0.81) | 0.79 (0.59-0.92) | 0.96 (0.82-1.00) | 0.27 | 0.00 | 0.96 |  | CMPK2;IFI27;IFI44L;RSAD2;IFIT1;IFI6;OAS3;OASL;OAS2;USP18 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 12777 | 0.36 (0.19-0.56) | 0.61 (0.41-0.78) | 0.79 (0.59-0.92) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | CMPK2;RSAD2;IFI44L;IFIT1;OASL;IFI27;OAS2;OAS3;IFI6;HERC5 |
| degora_slice/locked | ok | 12777 | 0.32 (0.16-0.52) | 0.50 (0.31-0.69) | 0.82 (0.63-0.94) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | CMPK2;IFI27;IFI44L;USP18;OASL;MX2;OAS2;OAS3;SAMD9;HERC5 |
| weighted_stouffer/default | ok | 12777 | 0.32 (0.16-0.52) | 0.50 (0.31-0.69) | 0.82 (0.63-0.94) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | CMPK2;IFI27;IFI44L;USP18;OASL;MX2;OAS2;OAS3;SAMD9;HERC5 |
| unweighted_stouffer/default | ok | 12777 | 0.32 (0.16-0.52) | 0.54 (0.34-0.72) | 0.82 (0.63-0.94) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | CMPK2;IFI27;IFI44L;USP18;OASL;OAS2;SAMD9;OAS3;HERC5;MX2 |
| fisher/default | ok | 12777 | 0.29 (0.13-0.49) | 0.50 (0.31-0.69) | 0.75 (0.55-0.89) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | CMPK2;IFI27;IFI44L;USP18;MX2;OASL;OAS2;OAS3;SAMD9;DHX58 |
| rank_product_approx/default | ok | 12777 | 0.32 (0.16-0.52) | 0.57 (0.37-0.76) | 0.79 (0.59-0.92) | 0.96 (0.82-1.00) | 0.27 | 0.00 | 0.96 |  | CMPK2;TAP1;MX1;OASL;IFIT1;HERC5;MX2;RSAD2;USP18;OAS3 |
| sign_vote/default | ok | 12777 | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 (0.00-0.12) | 0.00 | 1.00 | 0.00 |  | AADAT;AASDH;AATK;ABCB6;ABCB8;ABCC2;ABCC5;ABCD1;ABCF1;ABHD10 |
| metavolcanor/default | ok | 12777 | 0.29 (0.13-0.49) | 0.50 (0.31-0.69) | 0.75 (0.55-0.89) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | CMPK2;IFI27;IFI44L;USP18;MX2;OASL;OAS2;OAS3;SAMD9;DHX58 |
| robustrankaggreg/default | ok | 12777 | 0.32 (0.16-0.52) | 0.50 (0.31-0.69) | 0.75 (0.55-0.89) | 0.93 (0.76-0.99) | 0.26 | 0.00 | 0.93 |  | CMPK2;HERC5;IFIT1;OAS3;OASL;HELZ2;RSAD2;IFI6;IFIH1;MX2 |
| awfisher/default | ok | 12777 | 0.29 (0.13-0.49) | 0.50 (0.31-0.69) | 0.71 (0.51-0.87) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | CMPK2;IFI27;IFI44L;USP18;MX2;OASL;OAS2;OAS3;SAMD9;DHX58 |
| metarnaseq_fisher/default | ok | 12777 | 0.29 (0.13-0.49) | 0.50 (0.31-0.69) | 0.71 (0.51-0.87) | 0.89 (0.72-0.98) | 0.25 | 0.00 | 0.89 |  | CMPK2;IFI27;IFI44L;USP18;MX2;OASL;OAS2;OAS3;SAMD9;DHX58 |
| metarnaseq_invnorm/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metarnaseq_invnorm_uninformative_sparse_public_summary |  |
| hstouffer/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | hstouffer_deg_table_materializer_blocked |  |
| awmeta/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | awmeta_variance_inputs_missing |  |
| rankprod_exact/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | rankprod_exact_requires_expression_or_origin_labels |  |
| metade/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metade_effect_size_modes_require_variance_or_raw_expression |  |
| dexma/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | dexma_requires_expression_matrices_and_phenotype_metadata |  |
| metaintegrator/default | blocked | 0 |  (-) |  (-) |  (-) |  (-) |  |  |  | metaintegrator_requires_expression_and_phenotype_objects |  |

Interpretation guardrail: hStouffer/AWmeta rows are blocked, not beaten; no superiority claim is allowed from blocked comparators.
