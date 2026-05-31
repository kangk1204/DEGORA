# IFN pilot deep benchmark metrics

## Point Metrics

| method | status | recall@10 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | AUROC | AUPRC enrichment | AURC | top10 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| awfisher/default | ok | 0.214 (0.083-0.410) | 0.750 (0.551-0.893) | 0.893 (0.718-0.977) | 0.250 | 0.998 | 177.921 | 0.900 | TAP1;CMPK2;RSAD2;MX1;TRIM21;MX2;DDX60;IFIT3;PLSCR1;IFIT1 |
| degora_deg_score/v1_2_source_unit_mean | ok | 0.286 (0.132-0.487) | 0.821 (0.631-0.939) | 0.929 (0.765-0.991) | 0.260 | 0.999 | 233.239 | 0.906 | RSAD2;CMPK2;MX1;TAP1;IFIT1;IFIT3;DDX60;OASL;IFIH1;OAS1 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 0.286 (0.132-0.487) | 0.821 (0.631-0.939) | 0.929 (0.765-0.991) | 0.260 | 0.999 | 233.504 | 0.906 | RSAD2;CMPK2;MX1;TAP1;IFIT1;IFIT3;DDX60;OASL;IFIH1;OAS1 |
| degora_slice/locked | ok | 0.179 (0.061-0.369) | 0.786 (0.590-0.917) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 175.698 | 0.901 | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;IFIT3;SAMD9;PLSCR1 |
| fisher/default | ok | 0.179 (0.061-0.369) | 0.786 (0.590-0.917) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 176.535 | 0.901 | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| metarnaseq_fisher/default | ok | 0.179 (0.061-0.369) | 0.786 (0.590-0.917) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 176.535 | 0.901 | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| metarnaseq_invnorm/default | ok | 0.179 (0.061-0.369) | 0.786 (0.590-0.917) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 174.410 | 0.901 | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;SAMD9;IFIT3;PLSCR1 |
| metavolcanor/default | ok | 0.179 (0.061-0.369) | 0.786 (0.590-0.917) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 176.535 | 0.901 | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT3;PLSCR1;IFIT1;SAMD9 |
| rank_product_approx/default | ok | 0.214 (0.083-0.410) | 0.786 (0.590-0.917) | 0.893 (0.718-0.977) | 0.250 | 0.998 | 184.673 | 0.900 | TAP1;MX1;CMPK2;TRIM21;RSAD2;IFIT1;APOL6;HERC5;IFIT3;GBP1 |
| robustrankaggreg/default | ok | 0.143 (0.040-0.327) | 0.714 (0.513-0.868) | 0.821 (0.631-0.939) | 0.230 | 0.997 | 139.982 | 0.886 | TAP1;TRIM21;MAK;MT1M;TRIM69;HERC5;APOL6;CMPK2;IFIT1;MX1 |
| sign_vote/default | ok | 0.000 (0.000-0.123) | 0.000 (0.000-0.123) | 0.000 (0.000-0.123) | 0.000 | 0.760 | 1.869 | 0.025 | AADAT;AAK1;AARS2;AASDH;AASDHPPT;AATK;ABCB11;ABCB6;ABCB8;ABCC2 |
| unweighted_stouffer/default | ok | 0.214 (0.083-0.410) | 0.750 (0.551-0.893) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 174.831 | 0.900 | TAP1;MX1;RSAD2;CMPK2;TRIM21;DDX60;IFIT1;SAMD9;HERC5;IFIT3 |
| weighted_stouffer/default | ok | 0.179 (0.061-0.369) | 0.786 (0.590-0.917) | 0.929 (0.765-0.991) | 0.260 | 0.998 | 175.698 | 0.901 | TAP1;CMPK2;RSAD2;MX1;TRIM21;DDX60;IFIT1;IFIT3;SAMD9;PLSCR1 |
| rankprod_exact/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metade/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| dexma/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metaintegrator/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| hstouffer/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| awmeta/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |

## Source-Unit Bootstrap

| method | metric | mean | 95% CI | n |
| --- | --- | ---: | --- | ---: |
| degora_quality_weighted_score | aurc_at_max_k | 0.908 | 0.864-0.935 | 50 |
| degora_quality_weighted_score | recall_at_100 | 0.917 | 0.857-0.929 | 50 |
| fisher | aurc_at_max_k | 0.895 | 0.849-0.909 | 50 |
| fisher | recall_at_100 | 0.895 | 0.786-0.929 | 50 |
| weighted_stouffer | aurc_at_max_k | 0.895 | 0.849-0.909 | 50 |
| weighted_stouffer | recall_at_100 | 0.895 | 0.786-0.929 | 50 |

## DEGORA Advantage Metrics

| subset | n | median source units | median concordance | median LOO stability |
| --- | ---: | ---: | ---: | ---: |
| all_scored | 9860 | 2.0 | 1.0 | 0.874899 |
| top100_primary | 100 | 2.0 | 1.0 | 0.999239 |
| top100_quality_weighted | 100 | 2.0 | 1.0 | 0.999239 |
| locked_gold | 26 | 2.0 | 1.0 | 0.9994675 |
