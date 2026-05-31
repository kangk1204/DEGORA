# ER Stress/UPR Full Source-Support Report

Score table: `outputs/results/er-stress-benchmark/piper_gene_scores.csv`.
Locked gold panel: `data/studies/gold/er_stress_upr_gold_panel.csv`.
Interpretive marker panel: `data/studies/interpretive/er_stress_upr_marker_panel.csv`.

Guardrail: locked recall rows are validation diagnostics. Interpretive marker rows are biology annotation only and must not be reported as pre-locked benchmark performance.

Source-support criterion: `n_source_units >= 2` and `sign_concordance >= 1.0` with expected direction when declared.

## Summary

| metric | K | recovered / denominator | recall | recovered genes |
| --- | ---: | ---: | ---: | --- |
| locked_membership_recall | 10 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| locked_direction_recall | 10 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| locked_source_supported_recall | 10 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| known_or_locked_marker_coverage | 10 | 2 / 22 | 0.09 | DDIT3;PPP1R15A |
| known_or_locked_marker_source_supported | 10 | 2 / 22 | 0.09 | DDIT3;PPP1R15A |
| locked_membership_recall | 20 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| locked_direction_recall | 20 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| locked_source_supported_recall | 20 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| known_or_locked_marker_coverage | 20 | 2 / 22 | 0.09 | DDIT3;PPP1R15A |
| known_or_locked_marker_source_supported | 20 | 2 / 22 | 0.09 | DDIT3;PPP1R15A |
| locked_membership_recall | 50 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| locked_direction_recall | 50 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| locked_source_supported_recall | 50 | 2 / 18 | 0.11 | DDIT3;PPP1R15A |
| known_or_locked_marker_coverage | 50 | 2 / 22 | 0.09 | DDIT3;PPP1R15A |
| known_or_locked_marker_source_supported | 50 | 2 / 22 | 0.09 | DDIT3;PPP1R15A |
| locked_membership_recall | 100 | 3 / 18 | 0.17 | DDIT3;DERL3;PPP1R15A |
| locked_direction_recall | 100 | 3 / 18 | 0.17 | DDIT3;DERL3;PPP1R15A |
| locked_source_supported_recall | 100 | 3 / 18 | 0.17 | DDIT3;DERL3;PPP1R15A |
| known_or_locked_marker_coverage | 100 | 3 / 22 | 0.14 | DDIT3;DERL3;PPP1R15A |
| known_or_locked_marker_source_supported | 100 | 3 / 22 | 0.14 | DDIT3;DERL3;PPP1R15A |

## Top Genes

| rank | gene | locked | marker | direction | sources | concordance | reliability | role |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| 1 | ADAM19 | no | no | down | 3.000 | 1.000 | 98.966 |  |
| 2 | DDIT3 | yes | no | up / expected up | 3.000 | 1.000 | 98.499 |  |
| 3 | HAS2 | no | no | down | 3.000 | 1.000 | 98.668 |  |
| 4 | PPP1R15A | yes | no | up / expected up | 3.000 | 1.000 | 98.764 |  |
| 5 | ACTB | no | no | down | 3.000 | 1.000 | 98.721 |  |
| 6 | DUSP6 | no | no | down | 3.000 | 1.000 | 98.742 |  |
| 7 | LFNG | no | no | down | 3.000 | 1.000 | 98.744 |  |
| 8 | LMNB1 | no | no | down | 3.000 | 1.000 | 98.636 |  |
| 9 | MGP | no | no | down | 3.000 | 1.000 | 98.405 |  |
| 10 | SNAI1 | no | no | down | 3.000 | 1.000 | 98.593 |  |
