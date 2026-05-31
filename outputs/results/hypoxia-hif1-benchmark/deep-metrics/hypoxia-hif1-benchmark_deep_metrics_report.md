# Hypoxia/HIF1 deep benchmark metrics

## Point Metrics

| method | status | recall@10 (95% CI) | recall@50 (95% CI) | recall@100 (95% CI) | precision@100 | AUROC | AUPRC enrichment | AURC | top10 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| awfisher/default | ok | 0.212 (0.090-0.389) | 0.727 (0.545-0.867) | 0.818 (0.645-0.930) | 0.270 | 0.967 | 484.171 | 0.884 | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| degora_deg_score/v1_2_source_unit_mean | ok | 0.182 (0.070-0.355) | 0.758 (0.577-0.889) | 0.848 (0.681-0.949) | 0.280 | 0.968 | 557.262 | 0.871 | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;ADM;ANKRD37;AK4;HK2;PPFIA4 |
| degora_quality_weighted_score/quality_weighted_secondary | ok | 0.212 (0.090-0.389) | 0.758 (0.577-0.889) | 0.848 (0.681-0.949) | 0.280 | 0.967 | 548.340 | 0.871 | EGLN3;NDRG1;ANGPTL4;PFKFB4;TMEM45A;AK4;ANKRD37;HK2;ADM;VEGFA |
| degora_slice/locked | ok | 0.242 (0.111-0.423) | 0.818 (0.645-0.930) | 0.848 (0.681-0.949) | 0.280 | 0.964 | 652.332 | 0.883 | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| fisher/default | ok | 0.212 (0.090-0.389) | 0.727 (0.545-0.867) | 0.818 (0.645-0.930) | 0.270 | 0.965 | 496.096 | 0.882 | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| metarnaseq_fisher/default | ok | 0.303 (0.156-0.487) | 0.727 (0.545-0.867) | 0.818 (0.645-0.930) | 0.270 | 0.968 | 694.436 | 0.886 | VEGFA;BNIP3L;PGK1;P4HA1;ADM;HK2;PFKFB4;ALDOA;NDRG1;BHLHE40 |
| metavolcanor/default | ok | 0.212 (0.090-0.389) | 0.727 (0.545-0.867) | 0.818 (0.645-0.930) | 0.270 | 0.965 | 507.268 | 0.882 | ADM;AK4;ALDOA;ALDOC;ANGPTL4;ANKRD37;ANKZF1;BHLHE40;BNIP3;BNIP3L |
| rank_product_approx/default | ok | 0.182 (0.070-0.355) | 0.727 (0.545-0.867) | 0.818 (0.645-0.930) | 0.270 | 0.949 | 536.589 | 0.859 | BNIP3L;TMEM45A;NDRG1;P4HA1;VEGFA;SLC2A1;AK4;MIR210HG;PFKFB4;HIF1A-AS3 |
| robustrankaggreg/default | ok | 0.242 (0.111-0.423) | 0.788 (0.611-0.910) | 0.818 (0.645-0.930) | 0.270 | 0.960 | 583.859 | 0.846 | AK4;BNIP3L;ALDOA;P4HA1;MXI1;NDRG1;KDM4B;ADM;TMEM45A;KDM3A |
| sign_vote/default | ok | 0.030 (0.001-0.158) | 0.121 (0.034-0.282) | 0.273 (0.133-0.455) | 0.090 | 0.944 | 55.993 | 0.526 | AEN;AIFM1;AK4;ALDOA;APEH;APOO;ATAD3A;ATAD3B;ATIC;ATP6V0A2 |
| unweighted_stouffer/default | ok | 0.273 (0.133-0.455) | 0.818 (0.645-0.930) | 0.848 (0.681-0.949) | 0.280 | 0.966 | 659.834 | 0.883 | ADM;AK4;ALDOA;BNIP3L;HK2;MXI1;NDRG1;P4HA1;PFKFB4;PGK1 |
| weighted_stouffer/default | ok | 0.242 (0.111-0.423) | 0.818 (0.645-0.930) | 0.848 (0.681-0.949) | 0.280 | 0.964 | 652.332 | 0.883 | ADM;AK4;BNIP3L;P4HA1;VEGFA;PGK1;PFKFB4;NDRG1;MXI1;TMEM45A |
| metarnaseq_invnorm/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| rankprod_exact/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metade/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| dexma/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| metaintegrator/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| hstouffer/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |
| awmeta/default | blocked |  (-) |  (-) |  (-) |  |  |  |  |  |

## Source-Unit Bootstrap

| method | metric | mean | 95% CI | n |
| --- | --- | ---: | --- | ---: |
| degora_quality_weighted_score | aurc_at_max_k | 0.859 | 0.827-0.878 | 50 |
| degora_quality_weighted_score | recall_at_100 | 0.823 | 0.788-0.848 | 50 |
| fisher | aurc_at_max_k | 0.849 | 0.783-0.889 | 50 |
| fisher | recall_at_100 | 0.750 | 0.458-0.872 | 50 |
| weighted_stouffer | aurc_at_max_k | 0.860 | 0.813-0.886 | 50 |
| weighted_stouffer | recall_at_100 | 0.825 | 0.727-0.909 | 50 |

## DEGORA Advantage Metrics

| subset | n | median source units | median concordance | median LOO stability |
| --- | ---: | ---: | ---: | ---: |
| all_scored | 32687 | 7.0 | 0.8044363856228389 | 0.995686 |
| top100_primary | 100 | 15.0 | 0.9939505969145692 | 0.999939 |
| top100_quality_weighted | 100 | 15.0 | 0.9939505969145692 | 0.999939 |
| locked_gold | 33 | 15.0 | 0.9939345208113994 | 0.999969 |
