# IFN PIPER Source-Support Report

Score table: `outputs/results/ifn-pilot/piper_gene_scores.csv`.
Locked gold panel: `data/studies/gold/ifn_gold_panel.csv`.
Interpretive marker panel: `data/studies/interpretive/ifn_type_i_marker_panel.csv`.

Guardrail: locked recall rows are validation diagnostics. Interpretive marker rows are biology annotation only and must not be reported as pre-locked benchmark performance.

Source-support criterion: `n_source_units >= 2` and `sign_concordance >= 1.0` with expected direction when declared.

## Summary

| metric | K | recovered / denominator | recall | recovered genes |
| --- | ---: | ---: | ---: | --- |
| locked_membership_recall | 10 | 7 / 20 | 0.35 | IFIH1;IFIT1;IFIT3;MX1;OAS1;OASL;RSAD2 |
| locked_direction_recall | 10 | 7 / 20 | 0.35 | IFIH1;IFIT1;IFIT3;MX1;OAS1;OASL;RSAD2 |
| locked_source_supported_recall | 10 | 7 / 20 | 0.35 | IFIH1;IFIT1;IFIT3;MX1;OAS1;OASL;RSAD2 |
| known_or_locked_marker_coverage | 10 | 9 / 23 | 0.39 | CMPK2;IFIH1;IFIT1;IFIT3;MX1;OAS1;OASL;RSAD2;TAP1 |
| known_or_locked_marker_source_supported | 10 | 9 / 23 | 0.39 | CMPK2;IFIH1;IFIT1;IFIT3;MX1;OAS1;OASL;RSAD2;TAP1 |
| locked_membership_recall | 20 | 9 / 20 | 0.45 | IFIH1;IFIT1;IFIT3;MX1;OAS1;OAS2;OAS3;OASL;RSAD2 |
| locked_direction_recall | 20 | 9 / 20 | 0.45 | IFIH1;IFIT1;IFIT3;MX1;OAS1;OAS2;OAS3;OASL;RSAD2 |
| locked_source_supported_recall | 20 | 9 / 20 | 0.45 | IFIH1;IFIT1;IFIT3;MX1;OAS1;OAS2;OAS3;OASL;RSAD2 |
| known_or_locked_marker_coverage | 20 | 12 / 23 | 0.52 | CMPK2;GBP1;IFIH1;IFIT1;IFIT3;MX1;OAS1;OAS2;OAS3;OASL;RSAD2;TAP1 |
| known_or_locked_marker_source_supported | 20 | 12 / 23 | 0.52 | CMPK2;GBP1;IFIH1;IFIT1;IFIT3;MX1;OAS1;OAS2;OAS3;OASL;RSAD2;TAP1 |
| locked_membership_recall | 50 | 17 / 20 | 0.85 | IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1 |
| locked_direction_recall | 50 | 17 / 20 | 0.85 | IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1 |
| locked_source_supported_recall | 50 | 17 / 20 | 0.85 | IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1 |
| known_or_locked_marker_coverage | 50 | 20 / 23 | 0.87 | CMPK2;GBP1;IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;TAP1 |
| known_or_locked_marker_source_supported | 50 | 20 / 23 | 0.87 | CMPK2;GBP1;IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;TAP1 |
| locked_membership_recall | 100 | 18 / 20 | 0.90 | IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;STAT2 |
| locked_direction_recall | 100 | 18 / 20 | 0.90 | IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;STAT2 |
| locked_source_supported_recall | 100 | 18 / 20 | 0.90 | IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;STAT2 |
| known_or_locked_marker_coverage | 100 | 21 / 23 | 0.91 | CMPK2;GBP1;IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;STAT2;TAP1 |
| known_or_locked_marker_source_supported | 100 | 21 / 23 | 0.91 | CMPK2;GBP1;IFI27;IFI44;IFI44L;IFI6;IFIH1;IFIT1;IFIT2;IFIT3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;STAT1;STAT2;TAP1 |

## Top Genes

| rank | gene | locked | marker | direction | sources | concordance | reliability | role |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| 1 | RSAD2 | yes | no | up / expected up | 2.000 | 1.000 | 98.672 |  |
| 2 | CMPK2 | no | yes | up / expected up | 2.000 | 1.000 | 98.658 | antiviral_interferon_stimulated_gene |
| 3 | MX1 | yes | no | up / expected up | 2.000 | 1.000 | 98.652 |  |
| 4 | TAP1 | no | yes | up / expected up | 2.000 | 1.000 | 98.727 | interferon_antigen_processing |
| 5 | IFIT1 | yes | no | up / expected up | 2.000 | 1.000 | 98.590 |  |
| 6 | IFIT3 | yes | no | up / expected up | 2.000 | 1.000 | 98.555 |  |
| 7 | DDX60 | no | no | up | 2.000 | 1.000 | 98.617 |  |
| 8 | OASL | yes | no | up / expected up | 2.000 | 1.000 | 98.539 |  |
| 9 | IFIH1 | yes | no | up / expected up | 2.000 | 1.000 | 98.539 |  |
| 10 | OAS1 | yes | no | up / expected up | 2.000 | 1.000 | 98.547 |  |
