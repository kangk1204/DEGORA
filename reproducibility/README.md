# Reproduce the paper's benchmarks

This folder lets you reproduce the seven benchmark rows from the DEGORA paper. Each
prepared **dataset** is a ready-to-run DEGORA config; the differential-expression (DEG)
tables the configs read are distributed separately as an archived data bundle.

```
reproducibility/
  README.md                       this file
  fetch_reproduction_data.sh      downloads + unpacks the archived DEG tables
  datasets/
    01_ifn_rnaseq/config.xlsx
    02_ifn_rnaseq_array/config.xlsx
    03_er_stress_rnaseq/config.xlsx
    04_er_stress_rnaseq_array/config.xlsx
    05_heat_shock_rnaseq/config.xlsx
    06_hypoxia_hif1_rnaseq/config.xlsx
    07_hypoxia_hif1_rnaseq_array/config.xlsx
```

## The seven datasets

| Dataset | Topic | Arm | Contrasts | Independent source units |
|---|---|---|---|---|
| `01_ifn_rnaseq` | Interferon (IFN) | RNA-seq | 7 | 2 |
| `02_ifn_rnaseq_array` | Interferon (IFN) | RNA-seq + microarray | 12 | 5 |
| `03_er_stress_rnaseq` | ER stress / UPR (tunicamycin) | RNA-seq | 6 | 3 |
| `04_er_stress_rnaseq_array` | ER stress / UPR (tunicamycin) | RNA-seq + microarray | 7 | 4 |
| `05_heat_shock_rnaseq` | Heat shock / HSF1 | RNA-seq | 5 | 3 |
| `06_hypoxia_hif1_rnaseq` | Hypoxia / HIF1 | RNA-seq | 19 | 15 |
| `07_hypoxia_hif1_rnaseq_array` | Hypoxia / HIF1 | RNA-seq + microarray | 26 | 20 |

Each `config.xlsx` is the exact workbook used for the paper. Its `Contrasts` sheet lists
every DEG table and the column mapping; its `GoldPanel` sheet holds the locked canonical
genes used only for the post-run recall lookup (the gold genes are **not** an input to the
score); its `AdvancedSettings` sheet pins the source-unit collapse rule.

## How to run

You only need to fetch the data once.

```bash
# 1. Install DEGORA (see the main repo README).

# 2. From the repo root, download + unpack the archived DEG tables.
#    This writes data/deg/raw/... and data/studies/gold/... at the repo root.
bash reproducibility/fetch_reproduction_data.sh

# 3. Run any dataset. Example: the IFN RNA-seq row.
degora run reproducibility/datasets/01_ifn_rnaseq/config.xlsx

# 4. Open the dashboard for that run (Ctrl+C to stop the server).
degora serve reproducibility/datasets/01_ifn_rnaseq/results/degora_scores.db
```

Run the other datasets the same way - just change the folder name, e.g.
`degora run reproducibility/datasets/05_heat_shock_rnaseq/config.xlsx`.

Each run writes into that dataset's own `results/` subfolder: the ranked gene table
(`degora_gene_scores.csv`), the per-gene/per-source evidence table, the SQLite evidence
database (`degora_scores.db`), harmonized intermediate tables, provenance sidecars, and an
Excel audit workbook. Recall of each topic's locked gold panel reproduces the corresponding
row of Table 1 in the paper.

> If `fetch_reproduction_data.sh` reports that the download URL is still a `PLACEHOLDER`,
> the data archive has not been wired in yet. Download the bundle from the Zenodo record in
> the paper's Data Availability statement and pass it directly:
> `DEGORA_REPRO_DATA_ZIP=/path/to/degora_reproduction_data_v1.zip bash reproducibility/fetch_reproduction_data.sh`

## What's in the data bundle

The archive `degora_reproduction_data_v1.zip` (~95 MB compressed, ~169 MB unpacked)
contains the 50 DEG tables the seven configs reference, plus the five locked gold panels,
laid out under the same `data/deg/raw/...` and `data/studies/gold/...` paths the configs
expect. `MANIFEST.csv` lists every file with its topic, source GEO accession, byte size, and
SHA-256.

## Where the data came from

The Supplementary File **"Benchmark Data Collection and DEGORA Configuration"** (submitted
with the paper) documents, per dataset, every source accession and download URL, how each DEG
table was derived, and how each config was built. If you prefer to regenerate the DEG tables
from the original raw data instead of downloading them, that file lists the exact fetch and
derivation commands for each topic.
