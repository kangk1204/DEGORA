# ER Stress / UPR Benchmark Report

ER stress / UPR is reported with the predeclared quality-weighted PIPER-DEG ranking, alongside primary and full sensitivity source sets. The claim is directional evidence integration and auditable cross-platform support, not cherry-picked superiority from a favorable subset.

## Primary Set

| method | status | recall@10 | recall@20 | recall@50 | recall@100 | dir@50 | dir@100 | top10 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| PIPER-DEG score | ok | 0.39 | 0.61 | 0.72 | 0.72 | 0.72 | 0.72 | DDIT3;PTX3;PPP1R15A;HSP90B1;HSPA5;CRELD2;HAS2;HERPUD1;ATF3;CHAC1 |
| PIPER-DEG quality-weighted score | ok | 0.39 | 0.61 | 0.72 | 0.72 | 0.72 | 0.72 | DDIT3;PTX3;PPP1R15A;HSP90B1;HSPA5;CRELD2;HAS2;HERPUD1;ATF3;CHAC1 |
| weighted Stouffer | ok | 0.39 | 0.56 | 0.61 | 0.72 | 0.61 | 0.72 | HSP90B1;PPP1R15A;SEL1L;DDIT3;CHAC1;PTX3;HSPA5;PDIA4;CRELD2;ADAM19 |
| Fisher | ok | 0.39 | 0.56 | 0.61 | 0.67 | 0.61 | 0.67 | HSP90B1;PPP1R15A;SEL1L;LOXL4;DDIT3;UNC5B;CHAC1;PTX3;HSPA5;PDIA4 |
| MetaVolcanoR | ok | 0.39 | 0.56 | 0.61 | 0.67 | 0.61 | 0.67 | HSP90B1;PPP1R15A;SEL1L;LOXL4;DDIT3;UNC5B;CHAC1;PTX3;HSPA5;PDIA4 |
| RobustRankAggreg | ok | 0.22 | 0.44 | 0.50 | 0.50 | 0.50 | 0.50 | DDIT3;ADAM19;GAS7;HAS2;HSP90B1;HYOU1;CRELD2;HSPA5;HERPUD1;ABLIM1 |

## Full Sensitivity

| method | status | recall@10 | recall@20 | recall@50 | recall@100 | dir@50 | dir@100 | top10 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| PIPER-DEG score | ok | 0.11 | 0.11 | 0.11 | 0.17 | 0.11 | 0.17 | ADAM19;DDIT3;HAS2;PPP1R15A;ACTB;DUSP6;LFNG;LMNB1;MGP;SNAI1 |
| PIPER-DEG quality-weighted score | ok | 0.28 | 0.44 | 0.61 | 0.61 | 0.61 | 0.61 | DDIT3;PPP1R15A;HAS2;ADAM19;PTX3;LFNG;ATF3;HSPA5;EVA1A;CHAC1 |
| weighted Stouffer | ok | 0.06 | 0.11 | 0.11 | 0.17 | 0.11 | 0.17 | ADAM19;PPP1R15A;ACTB;DUSP6;B4GALT5;FSCN1;POMP;LMNB1;LFNG;SNAI1 |
| Fisher | ok | 0.22 | 0.50 | 0.61 | 0.61 | 0.56 | 0.56 | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| MetaVolcanoR | ok | 0.22 | 0.50 | 0.61 | 0.61 | 0.56 | 0.56 | UNC5B;CRELD2;ADAM19;HSP90B1;KLHDC8A;SLC1A4;SEL1L;PDIA4;PTX3;PPP1R15A |
| RobustRankAggreg | ok | 0.17 | 0.22 | 0.50 | 0.50 | 0.50 | 0.50 | ADAM19;CRELD2;UNC5B;KLHDC8A;DDIT3;GAS7;DUSP6;HAS2;HSP90B1;HSPA5 |

## Source-Quality Tier Effect

- `piper_recall_primary_minus_full_at_10`: 0.11
- `piper_recall_primary_minus_full_at_20`: 0.17
- `piper_recall_primary_minus_full_at_50`: 0.11
- `piper_recall_primary_minus_full_at_100`: 0.11

## Quality-Weighted Secondary Score Effect

- `quality_weighted_minus_unweighted_score_full_at_10`: 0.17
- `quality_weighted_minus_unweighted_score_full_at_20`: 0.33
- `quality_weighted_minus_unweighted_score_full_at_50`: 0.50
- `quality_weighted_minus_unweighted_score_full_at_100`: 0.44

## Interpretation

- Use the quality-weighted PIPER-DEG row as the fixed manuscript-facing ranking for both primary and sensitivity source sets.
- Report primary and full source sets together; differences are source-quality sensitivity evidence, not a basis for choosing the most favorable subset.
- Keep the unweighted PIPER score as a reference index so the effect of source-quality weighting remains visible.
- The comparison claim is early-rank prioritization and evidence-DB usability, not universal statistical superiority.
- hStouffer and AWmeta remain faithful-input blockers, not defeated baselines.

## Stronger Comparison Strategy

- Pre-lock gold panel and source-quality tiers before scoring.
- Report recall@10, @20, @50, and @100; do not only report the most favorable cutoff.
- Use direction-aware recall when a gold panel declares expected directions; this exposes direction mismatches that simple membership recall hides.
- Use the source-support-aware summaries because PIPER's intended DB use is directional consistency across independent source units.
- Report source-quality diagnostics because normalized-matrix sources can conflict with full-table sources and should not silently flip the primary evidence story.
- Keep primary-quality and full-sensitivity datasets separate instead of hiding lower-confidence sources.
- Preserve per-gene evidence rows so top hits can be audited source by source.
