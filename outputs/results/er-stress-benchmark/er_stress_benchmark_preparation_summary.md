# ER Stress / UPR Benchmark Preparation

Gold panel was locked before scoring.

- Active contrasts: 8
- Independent source units: 5
- Gold genes: HSPA5, HSP90B1, DNAJB9, HERPUD1, DDIT3, ATF4, ATF3, PPP1R15A, CHAC1, ASNS, PDIA4, EDEM1, SEL1L, DERL3, MANF, TRIB3, DDIT4, XBP1

## Inclusion Decisions

| source | status | decision | reason |
| --- | --- | --- | --- |
| GSE84450 | species_sensitivity_source | excluded_from_manuscript_primary | as-published Cuffdiff full-result DTT ER-stress table; mouse source retained only for non-primary sensitivity after the human-only species gate |
| GSE102505 | active_primary_source | included | public count matrix with clear tunicamycin vs DMSO columns in three cell systems |
| GSE296996 | active_primary_source | included | public HepG2 count matrix with five tunicamycin vs five DMSO replicate columns |
| GSE245918 | active_primary_source | included | public RSEM count matrix with OCI-AML3 and HEK293T tunicamycin vs DMSO 24h arms mapped from the GEO series matrix |
| GSE103667 | species_sensitivity_source | excluded_from_manuscript_primary | normalized thapsigargin vs DMSO matrix with two replicates per arm; mouse source retained only for non-primary sensitivity after the human-only species gate |
| GSE84989 | deferred_candidate | excluded_pending_metadata_lock | sample-label metadata was not locked for this benchmark run |

## Guardrails

- hStouffer/AWmeta failures remain blocker evidence, not wins.
- Same-accession cell-line/time-course rows share a source unit.
- Derived count or normalized-matrix rows are labeled as derived inputs.
- The score is a prioritization index, not a calibrated probability.
