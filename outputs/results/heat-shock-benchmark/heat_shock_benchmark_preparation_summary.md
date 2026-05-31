# Heat Shock / HSF1 Benchmark Preparation

Gold panel was locked before scoring.

- Active contrasts: 5
- Independent source units: 3
- Gold genes: HSPA1A, HSPA1B, HSPA6, HSPH1, DNAJB1, HSPB1, BAG3, HSPB8, SERPINH1, HSPA4, HSPA4L, HSP90AA1, HSP90AB1, DNAJB6, DNAJC7, HSPA5

## Inclusion Decisions

| source | status | decision | reason |
| --- | --- | --- | --- |
| GSE164834 | active_primary_source | included | SW620 42C heat-shock and recovery arms have 3 control and 3 treatment replicates |
| GSE124609 | active_primary_source | included | WI38 young and senescent RNA-seq heat-shock contrasts have two replicates per arm |
| GSE132447 | active_primary_source | included | primary neuron 3h-recovery RNA-seq arm has complete 3 vs 3 count columns |
| GSE123980 | deferred_candidate | excluded_metadata_ambiguous_for_primary | K562 total RNA-seq sample titles indicate HS15/HS30 but treatment characteristics for total RNA rows do not explicitly encode heat exposure; nascent TT-seq rows are not a primary expression table |
| GSE73471 | deferred_candidate | excluded_low_control_replication | WI38 heat-shock source has only one untreated control in the available matrix |
| GSE130493 | deferred_candidate | excluded_no_matched_control | downloaded RAW archive contains heat-shock samples but no matched control samples |
| GSE57397 | deferred_candidate | excluded_no_untreated_control | HCT116 source contains heat-shock and recovery samples but no matched untreated/non-heat control arm |

## Guardrails

- All active inputs are full NCBI GEO RNA-seq count matrices, not DEG-only lists.
- Same-accession time-point or phenotype rows share a source unit.
- The locked panel asserts expected up-regulation under heat shock.
- The score is a prioritization index, not a calibrated probability.
