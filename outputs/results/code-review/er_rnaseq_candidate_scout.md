# ER-stress / UPR human bulk RNA-seq candidate scout

Scope: scouting + verification ONLY. No catalog edits, no regeneration, no git.
Goal: add >=1 INDEPENDENT human ER-stress/UPR bulk RNA-seq source unit (paper_id different from
GSE102505, GSE84450, GSE103667, GSE19519) so the human-only RNA-seq ER baseline reaches
min_studies=2 (>=2 independent human source units) when combined with GSE102505.

Date: 2026-05-29. Verification method: NCBI GEO record pages (acc.cgi), NCBI E-utilities
(esearch/esummary on db=gds), and NCBI FTP suppl/ directory listings. All accessions below were
confirmed to exist on GEO with the stated organism/assay.

---

## TASK 1 — In-repo candidate GSE84989: VERDICT = NOT USABLE (wrong species)

Files in repo: `data/deg/raw/er_stress/GSE84989_RAW.tar` + `GSE84989/GSM2255592..GSM2255603_*.cRPKM.txt.gz`
(12 per-sample cRPKM files). Source sidecar:
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84989/suppl/GSE84989_RAW.tar`

Findings:
- ORGANISM: **Mus musculus (mouse)** — NOT human. Confirmed two independent ways:
  1. GEO record (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84989): organism *Mus musculus*,
     murine **intestinal organoids**, Illumina Genome Analyzer RNA-seq.
  2. The per-sample data files themselves use **ENSMUSG** mouse Ensembl gene IDs
     (e.g. `ENSMUSG00000000001  83.37  24089`), which is definitive for mouse.
- DESIGN: 12 samples = CTL#1 x3 (GSM2255592-94), Thapsigargin x3 (GSM2255595-97), CTL#2 x3
  (GSM2255598-600), KRB nutrient-starvation x3 (GSM2255601-603). An ER-stress arm (thapsigargin vs
  control) IS definable, and the "sample labels not locked" concern CAN be resolved from the GEO
  record (the GEO titles map Sample1-12 to the four groups above).
- ASSAY: bulk RNA-seq — OK. But data are cRPKM (normalized), 2 cols/sample (cRPKM, raw read count
  in col 3); raw counts ARE present in col 3, so a count-derived contrast would be technically
  feasible.

WHY IT FAILS: It is mouse. It is excluded by the exact same human-only species gate that already
removed GSE84450 (mouse) and GSE103667 (mouse). It therefore CANNOT contribute to the human-only
RNA-seq ER baseline and does NOT solve the min_studies=2 problem. Do not pursue further for this
purpose. (It could only ever serve a future cross-species analysis, which is out of scope here.)

---

## TASK 2/3 — Ranked human ER-stress/UPR bulk RNA-seq candidates (best first)

All are Homo sapiens, "Expression profiling by high throughput sequencing" (bulk RNA-seq) unless
noted, and have a different paper_id than the four existing ER sources.

### #1 (TOP RECOMMENDATION) — GSE245918
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245918
- Paper: "Unfolded Protein Response governs an Alternative Splicing program conserved from healthy
  to malignant cells." SRA project PRJNA1030384.
- Organism: Homo sapiens. Assay: bulk RNA-seq (n=21 samples).
- Cell systems: OCI-AML3 (AML line) and HEK293T.
- UPR trigger + control: **Tunicamycin vs DMSO** (canonical UPR inducer; clean vehicle control).
  Design per GEO record: OCI-AML3 0h untreated x3; 4h DMSO x3; 4h Tunicamycin x3; 24h DMSO x3;
  24h Tunicamycin x3; HEK293T 24h DMSO x3; HEK293T 24h Tunicamycin x3.
- n_ctrl / n_treat: **3 / 3** per clean contrast (e.g. OCI-AML3 24h DMSO vs 24h TUN, or HEK293T
  24h DMSO vs 24h TUN). Matches the GSE102505 3v3 design exactly.
- Data availability (precise file): `GSE245918_rsem_counts.txt.gz` — a public **RSEM raw count
  matrix**, VERIFIED by direct download + inspection (FTP suppl/, ~998 KB).
  Shape confirmed: **60,671 gene rows x 22 columns** = `Ensembl.id` + 21 sample columns (A1..A21,
  matching the 21 GSM samples). First data row: `ENSG00000000003 0 0 ... 1804 1788 2106 ...`.
  Values are **integer counts** (not floats) — i.e. these are integer expected counts, directly
  usable as the logCPM+Welch count input. Gene IDs are Ensembl gene IDs (ENSG...), so HGNC mapping
  is required (the repo already ships `data/deg/raw/ifn/hgnc_complete_set.txt`, used by the GSE102505
  pipeline). One open item: the A1..A21 column-to-group mapping must be read from the GEO series
  matrix / SOFT to assign which of the 21 columns are the chosen DMSO vs TUN arm.
- Curation path: EXACTLY the GSE102505 recipe — derive logCPM + Welch t-test from the public count
  matrix (`*_derived_logcpm_welch.csv` pattern). RSEM "expected counts" are non-integer but are the
  standard count input and behave identically under logCPM+Welch.
- Inclusion verdict: **INCLUDE (primary-eligible)**. Human, RNA-seq, canonical tunicamycin-vs-DMSO,
  3v3, public count matrix, count-derived (same pipeline as the kept GSE102505 contrasts),
  independent paper_id. Combined with GSE102505 this gives a clean 2 independent human source units
  -> min_studies=2 satisfied.
- Risk: (a) Ensembl gene IDs require ENSG->HGNC mapping (repo already does this for GSE102505 via
  the bundled hgnc set); (b) pick ONE cell-system contrast (e.g. HEK293T 24h DMSO vs 24h TUN, a
  non-cancer line, mirrors NHA) or aggregate within the single paper_id, consistent with the repo's
  time-course/cell-line collapse policy; (c) read the A1..A21 column-to-group mapping from the GEO
  series matrix before scoring.

### #2 (STRONG, dedicated ER-stress study) — GSE296996
- URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296996
- Title: "Gene expression profiling of tunicamycin-induced ER stress in HepG2 hepatoma cells."
- Organism: Homo sapiens. Assay: bulk RNA-seq (n=10).
- Cell system: HepG2 hepatoma. UPR trigger + control: **Tunicamycin vs mock control** (purpose-built
  ER-stress experiment — the cleanest topical match of any candidate).
- n_ctrl / n_treat: **5 / 5**, VERIFIED from the GEO sample table. Control = HepG2 Control rep1-5
  (GSM8981325-GSM8981329); Treatment = HepG2 Tunicamycin rep1-5 (GSM8981330-GSM8981334). Single
  treatment condition, no dose/time confound. This actually EXCEEDS the GSE102505 3v3 depth.
- Data availability (precise file): `GSE296996_count_matrix.txt.gz` (~510 KB) — public **raw count
  matrix**, VERIFIED present on the GEO record (author-supplied; the NCBI-generated rnaseq_counts
  TSV is NOT available for this series — the download endpoint returns HTTP 404). No author DEG
  table (must derive).
- Curation path: derive logCPM + Welch from `GSE296996_count_matrix.txt.gz` (EXACTLY the GSE102505
  recipe).
- Inclusion verdict: **INCLUDE (primary-eligible)** — fully verified: human, bulk RNA-seq, clean
  tunicamycin-vs-mock, 5v5, public author count matrix, count-derivable with the existing pipeline,
  independent paper_id. This is now co-equal with GSE245918 as a top pick and is the most topical.
- Risk: single cell line (HepG2); confirm the count-matrix gene-ID type (symbol vs Ensembl) when
  curating. Otherwise no open blockers.

### #3 (FALLBACK) — GSE281511 (liver) / GSE281194 (breast)
- URLs: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281511 (HCCLM3 liver cancer cells);
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281194 (breast cancer cells).
- Title: "RNA seq of [liver|breast] cancer cells treated with different concentration gradients of
  Thapsigargin (TG), Tunicamycin (TM), and FCCP." Organism: Homo sapiens. Bulk RNA-seq (n=44 each).
- Cell line (GSE281511): HCCLM3 hepatocellular carcinoma (Illumina HiSeq 2000), VERIFIED.
- UPR trigger + control: **Thapsigargin and Tunicamycin vs DMSO** (both canonical UPR inducers),
  plus FCCP (mitochondrial uncoupler — NOT a UPR trigger; exclude FCCP arms) and inhibitor
  combos (melatonin/ISRIB/4u8c/MG132 — exclude).
- n_ctrl / n_treat: **2 / 2 per dose**, VERIFIED. Clean arms exist, e.g. TM_DMSO_1/2 vs TM_2_1/2
  (2 uM TM), and TM 0.5/1 uM; same for TG. Replicate depth is the limiting factor (only 2 per arm,
  vs the >=3 preferred), so this is acceptable-but-flagged per the inclusion criteria.
- Data availability (precise file): `GSE281511_gene_counts.txt.gz` (~1.3 MB) — public **raw gene
  count matrix**, VERIFIED present on the GEO record. No author DEG table. (NCBI-generated counts
  TSV returns HTTP 404 for this series.)
- Curation path: pick ONE standard-dose TM-vs-DMSO (e.g. 2 uM) or TG-vs-DMSO arm; derive
  logCPM + Welch; drop FCCP and all inhibitor-combo arms.
- Inclusion verdict: **FALLBACK / acceptable-with-flag** — verified usable (clean 2v2 TM-vs-DMSO
  and TG-vs-DMSO arms with a public count matrix), but only 2 reps/arm, so flag replicate depth if
  used. Prefer GSE245918 or GSE296996 (both >=3/arm) over this.

### Also-relevant human UPR RNA-seq (lower priority; verify before use)
- GSE323964 / SuperSeries GSE323967 — "IRE1 drives a homeostatic response to reduced protein influx
  into the ER" (Homo sapiens, bulk RNA-seq, n=36/108). UPR-focused but the perturbation is largely
  genetic/IRE1-axis and translation-initiation interference rather than a clean drug-vs-vehicle ER
  stressor; confirm a tunicamycin/thapsigargin/DTT-vs-control arm exists before considering.
- GSE302088 — "Inhibition of Acetyl-CoA-carboxylase induces ER stress ..." (Homo sapiens, n=5).
  ER stress is induced via ACC inhibition (genetic/pharmacologic, not a canonical UPR drug);
  small n. Marginal; not recommended unless a canonical arm is present.
- GSE289189 — TNBC "Bulk RNA-seq" sub-series of the H4K20me3/UPR study (Homo sapiens, n=22).
  UPR-themed but the contrast is epigenetic (KMT5 axis), not a clean ER-stress-drug-vs-vehicle;
  not a direct ER-stress induction contrast.

### Explicitly REJECTED during scouting (do not use)
- GSE84989 — mouse (see Task 1).
- GSE96996 — ChIP-seq (UTX), not RNA-seq; FTP suppl confirms only RAW.tar of BED/WIG.
- GSE112841 — Affymetrix microarray (CEL files), not RNA-seq; also NDRG1/DFO, not ER stress.
- GSE112843 — brain-tissue RPKM, not an ER-stress contrast.
- GSE208030 — ATAC-seq after ER stress, not RNA-seq.
- GSE306046 / GSE305216 / GSE299113 — single-cell ("scRNA-Seq" / toxicology SCTr), excluded (sc).
- GSE275042 — Hi-C / 3D genome (MCF10A DNA damage + ER stress), not expression RNA-seq.
- GSE266904 — CRISPR screen results (TG/TN screen .txt), not an expression count matrix.
- GSE289805 / GSE288477 / GSE270547 — ChIP-seq / ATAC-seq.

---

## RECOMMENDATION

Two candidates are now FULLY verified and co-equal; either alone solves min_studies=2 with GSE102505.

- TOP PICK (most topical, deepest, cleanest): **GSE296996** (paper_id e.g. `GSE296996_HepG2_TUN`).
  Purpose-built human ER-stress study: HepG2, tunicamycin vs mock, **5 vs 5** (exceeds GSE102505's
  3v3), single uncon­founded condition, public author count matrix `GSE296996_count_matrix.txt.gz`
  (verified on the GEO record), count-derivable with the EXACT existing logCPM+Welch recipe.
  Independent paper_id. No open blockers (only confirm gene-ID type at curation time).
- CO-TOP PICK (broadest, two cell systems): **GSE245918** (paper_id e.g. `GSE245918_UPR_TUN`).
  Human bulk RNA-seq, canonical tunicamycin-vs-DMSO, 3v3 per arm, public integer RSEM count matrix
  `GSE245918_rsem_counts.txt.gz` (verified: 60,671 genes x 21 sample cols, ENSG IDs). Pick one
  contrast (e.g. HEK293T 24h DMSO vs TUN) or collapse within paper_id. One curation step: map the
  A1..A21 columns to groups via the series matrix.
- Either GSE296996 OR GSE245918, combined with GSE102505, yields 2 independent human RNA-seq source
  units -> min_studies=2 satisfied. (Adding BOTH would give 3 and is even more robust.)
- FALLBACK: **GSE281511** (HCCLM3 liver, TM/TG-vs-DMSO, verified `GSE281511_gene_counts.txt.gz`) or
  **GSE281194** (breast) — usable but only 2 reps/arm; flag replicate depth, drop FCCP/inhibitor arms.

Anti-fabrication note: every accession above was verified to exist on GEO with the stated
organism/assay via the GEO record and/or E-utilities esummary and/or FTP suppl listing. The two
previously-open confirmations are now CLOSED with direct evidence: (a) `GSE245918_rsem_counts.txt.gz`
was downloaded and inspected (integer counts, 60,671 x 21 samples + Ensembl.id, ENSG gene IDs);
(b) per-arm replicate counts and exact count-file names for GSE296996 (5v5, `GSE296996_count_matrix.txt.gz`)
and GSE281511 (2v2 dose arms, `GSE281511_gene_counts.txt.gz`) were verified from the GEO sample
tables. The NCBI-generated rnaseq_counts TSV endpoint returns HTTP 404 for all of these series, so
the author-supplied count files (named above) are the curation path. No accession or gene data was
invented; remaining items (column-to-group mapping for GSE245918; gene-ID type for GSE296996) are
minor curation-time checks, not existence questions.
