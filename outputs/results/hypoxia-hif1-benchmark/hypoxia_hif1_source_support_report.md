# Hypoxia / HIF1 Source-Support Report

Score table: `outputs/results/hypoxia-hif1-benchmark/piper_gene_scores.csv`.
Locked gold panel: `data/studies/gold/hypoxia_hif1_gold_panel.csv`.

Guardrail: locked recall rows are validation diagnostics. Interpretive marker rows are biology annotation only and must not be reported as pre-locked benchmark performance.

Source-support criterion: `n_source_units >= 2` and `sign_concordance >= 1.0` with expected direction when declared.

## Summary

| metric | K | recovered / denominator | recall | recovered genes |
| --- | ---: | ---: | ---: | --- |
| locked_membership_recall | 10 | 4 / 20 | 0.20 | EGLN3;HK2;SLC2A1;VEGFA |
| locked_direction_recall | 10 | 4 / 20 | 0.20 | EGLN3;HK2;SLC2A1;VEGFA |
| locked_source_supported_recall | 10 | 1 / 20 | 0.05 | VEGFA |
| locked_membership_recall | 20 | 7 / 20 | 0.35 | BNIP3;BNIP3L;EGLN3;HK2;SLC2A1;SLC2A3;VEGFA |
| locked_direction_recall | 20 | 7 / 20 | 0.35 | BNIP3;BNIP3L;EGLN3;HK2;SLC2A1;SLC2A3;VEGFA |
| locked_source_supported_recall | 20 | 2 / 20 | 0.10 | BNIP3L;VEGFA |
| locked_membership_recall | 50 | 13 / 20 | 0.65 | ALDOA;BNIP3;BNIP3L;DDIT4;EGLN1;EGLN3;HK2;LDHA;PDK1;PGK1;SLC2A1;SLC2A3;VEGFA |
| locked_direction_recall | 50 | 13 / 20 | 0.65 | ALDOA;BNIP3;BNIP3L;DDIT4;EGLN1;EGLN3;HK2;LDHA;PDK1;PGK1;SLC2A1;SLC2A3;VEGFA |
| locked_source_supported_recall | 50 | 5 / 20 | 0.25 | ALDOA;BNIP3L;EGLN1;LDHA;VEGFA |
| locked_membership_recall | 100 | 16 / 20 | 0.80 | ALDOA;BNIP3;BNIP3L;CA9;DDIT4;EGLN1;EGLN3;ENO1;HK2;LDHA;PDK1;PFKL;PGK1;SLC2A1;SLC2A3;VEGFA |
| locked_direction_recall | 100 | 16 / 20 | 0.80 | ALDOA;BNIP3;BNIP3L;CA9;DDIT4;EGLN1;EGLN3;ENO1;HK2;LDHA;PDK1;PFKL;PGK1;SLC2A1;SLC2A3;VEGFA |
| locked_source_supported_recall | 100 | 5 / 20 | 0.25 | ALDOA;BNIP3L;EGLN1;LDHA;VEGFA |

## Top Genes

| rank | gene | locked | marker | direction | sources | concordance | reliability | role |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| 1 | NDRG1 | no | no | up | 17.000 | 0.996 | 99.176 |  |
| 2 | EGLN3 | yes | no | up / expected up | 17.000 | 0.996 | 99.180 |  |
| 3 | ANKRD37 | no | no | up | 16.000 | 1.000 | 98.828 |  |
| 4 | AK4 | no | no | up | 17.000 | 1.000 | 99.899 |  |
| 5 | ADM | no | no | up | 16.000 | 0.981 | 98.089 |  |
| 6 | TMEM45A | no | no | up | 16.000 | 1.000 | 98.845 |  |
| 7 | HK2 | yes | no | up / expected up | 17.000 | 0.978 | 99.273 |  |
| 8 | VEGFA | yes | no | up / expected up | 16.000 | 1.000 | 98.833 |  |
| 9 | P4HA1 | no | no | up | 17.000 | 1.000 | 99.325 |  |
| 10 | SLC2A1 | yes | no | up / expected up | 17.000 | 0.978 | 98.495 |  |
| 11 | MXI1 | no | no | up | 17.000 | 0.997 | 99.773 |  |
| 12 | ANGPTL4 | no | no | up | 15.000 | 0.942 | 96.343 |  |
| 13 | BNIP3 | yes | no | up / expected up | 17.000 | 0.985 | 98.124 |  |
| 14 | STC1 | no | no | up | 16.000 | 0.991 | 98.596 |  |
| 15 | BNIP3L | yes | no | up / expected up | 17.000 | 1.000 | 99.899 |  |
| 16 | PPP1R3G | no | no | up | 16.000 | 1.000 | 98.809 |  |
| 17 | SLC2A3 | yes | no | up / expected up | 17.000 | 0.976 | 98.421 |  |
| 18 | PPFIA4 | no | no | up | 16.000 | 0.975 | 97.390 |  |
| 19 | PFKFB4 | no | no | up | 17.000 | 0.935 | 97.738 |  |
| 20 | PFKFB3 | no | no | up | 17.000 | 0.995 | 99.709 |  |
