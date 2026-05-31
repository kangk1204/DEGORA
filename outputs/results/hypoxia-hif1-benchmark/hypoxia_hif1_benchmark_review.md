# Hypoxia / HIF1 Benchmark Refresh - 2026-05-28

## Benchmark Definition

Hypoxia / HIF1 was promoted to a refreshed clean-mechanism benchmark because the gold biology is strong and externally recognizable: canonical HIF/hypoxia genes such as `VEGFA`, `CA9`, `EGLN3`, `SLC2A1`, `LDHA`, `PDK1`, and `BNIP3` should move up under low oxygen.

Gold panel:

- `data/studies/gold/hypoxia_hif1_gold_panel.csv`
- 20 expected-up HIF/hypoxia targets.
- Locked before the refreshed `hypoxia-hif1-benchmark` scoring and comparator run.
- External biological anchor: MSigDB `HALLMARK_HYPOXIA` describes genes up-regulated in response to low oxygen levels.

Input corpus:

- `data/studies/hypoxia_catalog.csv`
- 21 active contrasts.
- 17 independent source units.
- 379,898 harmonized rows.
- 33,865 consensus-scored genes.

## Quality-Weighted PIPER Result

| metric | value |
| --- | ---: |
| recall@10 | 0.15 |
| recall@20 | 0.35 |
| recall@50 | 0.65 |
| recall@100 | 0.80 |
| direction recall@100 | 0.80 |

Top 10 quality-weighted genes:

`NDRG1; EGLN3; AK4; ANKRD37; TMEM45A; ADM; HK2; VEGFA; MXI1; ANGPTL4`

Recovered gold genes by top 100:

`ALDOA; BNIP3; BNIP3L; CA9; DDIT4; EGLN1; EGLN3; ENO1; HK2; LDHA; PDK1; PFKL; PGK1; SLC2A1; SLC2A3; VEGFA`

Missing at top 100:

`ANGPT2; EPO; FLT1; HK1`

## Comparator Summary

| method | recall@10 | recall@20 | recall@50 | recall@100 | direction@100 |
| --- | ---: | ---: | ---: | ---: | ---: |
| PIPER quality-weighted | 0.15 | 0.35 | 0.65 | 0.80 | 0.80 |
| Weighted Stouffer | 0.20 | 0.40 | 0.70 | 0.75 | 0.75 |
| Unweighted Stouffer | 0.25 | 0.40 | 0.70 | 0.80 | 0.80 |
| Fisher | 0.15 | 0.35 | 0.70 | 0.75 | 0.75 |
| Rank product approx. | 0.15 | 0.40 | 0.70 | 0.75 | 0.75 |
| MetaVolcanoR | 0.15 | 0.35 | 0.70 | 0.75 | 0.75 |
| RobustRankAggreg | 0.15 | 0.35 | 0.70 | 0.75 | 0.75 |

Interpretation:

- Hypoxia/HIF1 is a strong additional benchmark biologically.
- PIPER quality-weighted is competitive but not dominant in this corpus.
- PIPER recovers more locked genes by top 100 than most non-Stouffer baselines and has highly plausible top genes, but weighted/unweighted Stouffer are at least as strong at early K.
- This benchmark should be used as a balanced validation point, not as a cherry-picked win.
- The strongest manuscript use is as an honest clean-mechanism validation set: canonical HIF genes are recovered, while classical p-value/rank methods remain competitive.

## Generated Artifacts

- `outputs/results/hypoxia-hif1-benchmark/piper_gene_scores.csv`
- `outputs/results/hypoxia-hif1-benchmark/piper_scores.db`
- `outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv`
- `outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_source_support_report.md`
- `outputs/figures/hypoxia-hif1-benchmark/hypoxia_hif1_gold_recall_curve_top1000.png`
- `outputs/results/hypoxia-hif1-benchmark/provenance_audit.json`

## Verification

- Score DB generated successfully.
- Baseline adapters generated manifest-safe outputs.
- Gold comparator summary generated successfully.
- Source-support report generated successfully.
- Top-1000 recall figure package generated successfully.
- Source provenance audit passed.
