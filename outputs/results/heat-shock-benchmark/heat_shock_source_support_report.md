# Heat Shock / HSF1 Source-Support Report

Score table: `outputs/results/heat-shock-benchmark/piper_gene_scores.csv`.
Locked gold panel: `data/studies/gold/heat_shock_hsf1_gold_panel.csv`.

Guardrail: locked recall rows are validation diagnostics. Interpretive marker rows are biology annotation only and must not be reported as pre-locked benchmark performance.

Source-support criterion: `n_source_units >= 2` and `sign_concordance >= 1.0` with expected direction when declared.

## Summary

| metric | K | recovered / denominator | recall | recovered genes |
| --- | ---: | ---: | ---: | --- |
| locked_membership_recall | 10 | 6 / 16 | 0.38 | BAG3;DNAJB1;HSPA1A;HSPA1B;HSPA6;HSPH1 |
| locked_direction_recall | 10 | 6 / 16 | 0.38 | BAG3;DNAJB1;HSPA1A;HSPA1B;HSPA6;HSPH1 |
| locked_source_supported_recall | 10 | 6 / 16 | 0.38 | BAG3;DNAJB1;HSPA1A;HSPA1B;HSPA6;HSPH1 |
| locked_membership_recall | 20 | 6 / 16 | 0.38 | BAG3;DNAJB1;HSPA1A;HSPA1B;HSPA6;HSPH1 |
| locked_direction_recall | 20 | 6 / 16 | 0.38 | BAG3;DNAJB1;HSPA1A;HSPA1B;HSPA6;HSPH1 |
| locked_source_supported_recall | 20 | 6 / 16 | 0.38 | BAG3;DNAJB1;HSPA1A;HSPA1B;HSPA6;HSPH1 |
| locked_membership_recall | 50 | 11 / 16 | 0.69 | BAG3;DNAJB1;HSP90AA1;HSP90AB1;HSPA1A;HSPA1B;HSPA6;HSPB1;HSPB8;HSPH1;SERPINH1 |
| locked_direction_recall | 50 | 11 / 16 | 0.69 | BAG3;DNAJB1;HSP90AA1;HSP90AB1;HSPA1A;HSPA1B;HSPA6;HSPB1;HSPB8;HSPH1;SERPINH1 |
| locked_source_supported_recall | 50 | 11 / 16 | 0.69 | BAG3;DNAJB1;HSP90AA1;HSP90AB1;HSPA1A;HSPA1B;HSPA6;HSPB1;HSPB8;HSPH1;SERPINH1 |
| locked_membership_recall | 100 | 12 / 16 | 0.75 | BAG3;DNAJB1;HSP90AA1;HSP90AB1;HSPA1A;HSPA1B;HSPA4L;HSPA6;HSPB1;HSPB8;HSPH1;SERPINH1 |
| locked_direction_recall | 100 | 12 / 16 | 0.75 | BAG3;DNAJB1;HSP90AA1;HSP90AB1;HSPA1A;HSPA1B;HSPA4L;HSPA6;HSPB1;HSPB8;HSPH1;SERPINH1 |
| locked_source_supported_recall | 100 | 12 / 16 | 0.75 | BAG3;DNAJB1;HSP90AA1;HSP90AB1;HSPA1A;HSPA1B;HSPA4L;HSPA6;HSPB1;HSPB8;HSPH1;SERPINH1 |

## Top Genes

| rank | gene | locked | marker | direction | sources | concordance | reliability | role |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| 1 | HSPA6 | yes | no | up / expected up | 3.000 | 1.000 | 99.083 |  |
| 2 | HSPA1A | yes | no | up / expected up | 3.000 | 1.000 | 99.072 |  |
| 3 | HSPA1B | yes | no | up / expected up | 3.000 | 1.000 | 99.039 |  |
| 4 | DNAJB1 | yes | no | up / expected up | 3.000 | 1.000 | 99.062 |  |
| 5 | BAG3 | yes | no | up / expected up | 3.000 | 1.000 | 99.046 |  |
| 6 | HSPA1L | no | no | up | 3.000 | 1.000 | 98.978 |  |
| 7 | HSPH1 | yes | no | up / expected up | 3.000 | 1.000 | 99.000 |  |
| 8 | FOS | no | no | up | 3.000 | 1.000 | 98.976 |  |
| 9 | DNAJB4 | no | no | up | 3.000 | 1.000 | 98.915 |  |
| 10 | MICB | no | no | up | 3.000 | 1.000 | 99.015 |  |
| 11 | ACTRT3 | no | no | up | 3.000 | 1.000 | 98.991 |  |
| 12 | CRYAB | no | no | up | 3.000 | 1.000 | 98.883 |  |
| 13 | DNAJA4 | no | no | up | 3.000 | 1.000 | 98.848 |  |
| 14 | CCDC121 | no | no | up | 3.000 | 1.000 | 99.018 |  |
| 15 | ZFAND2A | no | no | up | 3.000 | 1.000 | 98.884 |  |
| 16 | TENT5A | no | no | up | 3.000 | 1.000 | 98.897 |  |
| 17 | PNLDC1 | no | no | up | 3.000 | 1.000 | 98.802 |  |
| 18 | JMJD6 | no | no | up | 3.000 | 1.000 | 98.917 |  |
| 19 | FKBP4 | no | no | up | 3.000 | 1.000 | 98.905 |  |
| 20 | DEDD2 | no | no | up | 3.000 | 1.000 | 98.760 |  |
