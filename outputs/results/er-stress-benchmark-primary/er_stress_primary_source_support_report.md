# ER Stress/UPR Primary Source-Support Report

Score table: `outputs/results/er-stress-benchmark-primary/piper_gene_scores.csv`.
Locked gold panel: `data/studies/gold/er_stress_upr_gold_panel.csv`.
Interpretive marker panel: `data/studies/interpretive/er_stress_upr_marker_panel.csv`.

Guardrail: locked recall rows are validation diagnostics. Interpretive marker rows are biology annotation only and must not be reported as pre-locked benchmark performance.

Source-support criterion: `n_source_units >= 2` and `sign_concordance >= 1.0` with expected direction when declared.

## Summary

| metric | K | recovered / denominator | recall | recovered genes |
| --- | ---: | ---: | ---: | --- |
| locked_membership_recall | 10 | 7 / 18 | 0.39 | ATF3;CHAC1;DDIT3;HERPUD1;HSP90B1;HSPA5;PPP1R15A |
| locked_direction_recall | 10 | 7 / 18 | 0.39 | ATF3;CHAC1;DDIT3;HERPUD1;HSP90B1;HSPA5;PPP1R15A |
| locked_source_supported_recall | 10 | 7 / 18 | 0.39 | ATF3;CHAC1;DDIT3;HERPUD1;HSP90B1;HSPA5;PPP1R15A |
| known_or_locked_marker_coverage | 10 | 9 / 22 | 0.41 | ATF3;CHAC1;CRELD2;DDIT3;HERPUD1;HSP90B1;HSPA5;PPP1R15A;PTX3 |
| known_or_locked_marker_source_supported | 10 | 9 / 22 | 0.41 | ATF3;CHAC1;CRELD2;DDIT3;HERPUD1;HSP90B1;HSPA5;PPP1R15A;PTX3 |
| locked_membership_recall | 20 | 11 / 18 | 0.61 | ATF3;CHAC1;DDIT3;DERL3;HERPUD1;HSP90B1;HSPA5;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| locked_direction_recall | 20 | 11 / 18 | 0.61 | ATF3;CHAC1;DDIT3;DERL3;HERPUD1;HSP90B1;HSPA5;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| locked_source_supported_recall | 20 | 11 / 18 | 0.61 | ATF3;CHAC1;DDIT3;DERL3;HERPUD1;HSP90B1;HSPA5;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| known_or_locked_marker_coverage | 20 | 14 / 22 | 0.64 | ATF3;CHAC1;CRELD2;DDIT3;DERL3;HERPUD1;HSP90B1;HSPA5;HYOU1;PDIA4;PPP1R15A;PTX3;SEL1L;TRIB3 |
| known_or_locked_marker_source_supported | 20 | 14 / 22 | 0.64 | ATF3;CHAC1;CRELD2;DDIT3;DERL3;HERPUD1;HSP90B1;HSPA5;HYOU1;PDIA4;PPP1R15A;PTX3;SEL1L;TRIB3 |
| locked_membership_recall | 50 | 13 / 18 | 0.72 | ATF3;CHAC1;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;MANF;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| locked_direction_recall | 50 | 13 / 18 | 0.72 | ATF3;CHAC1;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;MANF;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| locked_source_supported_recall | 50 | 13 / 18 | 0.72 | ATF3;CHAC1;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;MANF;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| known_or_locked_marker_coverage | 50 | 16 / 22 | 0.73 | ATF3;CHAC1;CRELD2;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;HYOU1;MANF;PDIA4;PPP1R15A;PTX3;SEL1L;TRIB3 |
| known_or_locked_marker_source_supported | 50 | 16 / 22 | 0.73 | ATF3;CHAC1;CRELD2;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;HYOU1;MANF;PDIA4;PPP1R15A;PTX3;SEL1L;TRIB3 |
| locked_membership_recall | 100 | 13 / 18 | 0.72 | ATF3;CHAC1;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;MANF;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| locked_direction_recall | 100 | 13 / 18 | 0.72 | ATF3;CHAC1;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;MANF;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| locked_source_supported_recall | 100 | 13 / 18 | 0.72 | ATF3;CHAC1;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;MANF;PDIA4;PPP1R15A;SEL1L;TRIB3 |
| known_or_locked_marker_coverage | 100 | 17 / 22 | 0.77 | ATF3;CHAC1;CRELD2;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;HSPB1;HYOU1;MANF;PDIA4;PPP1R15A;PTX3;SEL1L;TRIB3 |
| known_or_locked_marker_source_supported | 100 | 16 / 22 | 0.73 | ATF3;CHAC1;CRELD2;DDIT3;DERL3;DNAJB9;HERPUD1;HSP90B1;HSPA5;HYOU1;MANF;PDIA4;PPP1R15A;PTX3;SEL1L;TRIB3 |

## Top Genes

| rank | gene | locked | marker | direction | sources | concordance | reliability | role |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| 1 | DDIT3 | yes | no | up / expected up | 2.000 | 1.000 | 98.676 |  |
| 2 | PTX3 | no | yes | up / expected up | 2.000 | 1.000 | 98.606 | contextual_er_stress_induced_inflammatory_gene |
| 3 | PPP1R15A | yes | no | up / expected up | 2.000 | 1.000 | 98.584 |  |
| 4 | HSP90B1 | yes | no | up / expected up | 2.000 | 1.000 | 98.585 |  |
| 5 | HSPA5 | yes | no | up / expected up | 2.000 | 1.000 | 98.576 |  |
| 6 | CRELD2 | no | yes | up / expected up | 2.000 | 1.000 | 98.571 | er_stress_inducible_gene |
| 7 | HAS2 | no | no | down | 2.000 | 1.000 | 98.449 |  |
| 8 | HERPUD1 | yes | no | up / expected up | 2.000 | 1.000 | 98.484 |  |
| 9 | ATF3 | yes | no | up / expected up | 2.000 | 1.000 | 98.438 |  |
| 10 | CHAC1 | yes | no | up / expected up | 2.000 | 1.000 | 98.506 |  |
