# DEGORA

DEGORA helps you combine many published DEG tables and inspect which genes have repeated, directionally consistent evidence.

You do not need YAML or JSON. The beginner workflow uses one Excel file.

## Very Easy Start

### 1. Install

Open a terminal in this folder and run:

```bash
python -m pip install -e outputs/code
```

Alternative:

```bash
make install
```

### 2. Try the Demo

Before using your own DEG tables, run the tiny demo once:

```bash
degora demo degora-demo
cd degora-demo
degora validate degora_demo_config.xlsx
degora run degora_demo_config.xlsx
degora serve results/degora_scores.db
```

Search for `ISG15`, `IFIT1`, `MX1`, or `OAS1` in the browser.

### 3. Make Your Excel Template

```bash
degora template degora_config.xlsx
```

Alternative:

```bash
make template DEGORA_CONFIG=degora_config.xlsx
```

### 4. Fill the Excel File

Open `degora_config.xlsx` and edit the `Contrasts` sheet.

Fill these columns first:

| Column | What it means |
|---|---|
| `study_id` | A unique name for this DEG contrast, for example `IFN_GSE001_4h` |
| `source_unit_id` | The paper/dataset ID. Use the same value for time points from the same study |
| `source_path` | Path to the DEG table file |
| `gene_column` | Column in the DEG table containing gene names |
| `lfc_column` | Column containing log2 fold change |
| `p_column` | Column containing p-value |
| `padj_column` | Column containing adjusted p-value/FDR, if available |
| `assay_type` | Optional assay label, for example `RNA-seq` or `microarray` |
| `source_input_type` | Optional source type, for example `author_deg_table`, `limma_full_table`, or `normalized_expression_matrix` |
| `platform` | Optional platform label, useful for microarrays, for example `GPL570` |
| `normalization` | Optional source normalization, for example `DESeq2`, `RMA/log2`, or `quantile/log2` |
| `probe_collapse` | For microarrays, how probes were collapsed to gene symbols |
| `table_scope` | Use `full_results` for all tested genes, `deg_only` for significant-gene-only lists, or leave `auto` |
| `rank_universe_size` | For `deg_only` lists, the number of genes originally tested if the paper reports it |
| `condition` | Human-readable contrast label, for example `IFN-beta vs untreated` |
| `time_h` | Time point in hours, if this is time-course data |
| `time_course_mode` | Optional. Use `mean` by default, or `early`, `late`, `peak_mean` for same-source time-course rows |
| `include` | `yes` to use the row, `no` to ignore it |

For time-course data, use one row per time point and keep the same `source_unit_id` for rows from the same paper or dataset. `time_course_mode=mean` is safest for most projects because it summarizes all time points without letting one study count as many independent studies.

### 5. Check the Excel File

```bash
degora validate degora_config.xlsx
```

If something is wrong, DEGORA explains the exact row/column and how to fix it.

### 6. Run and Open Results

```bash
degora run degora_config.xlsx
degora serve outputs/results/degora-run/degora_scores.db
```

Then open the URL printed in the terminal, usually:

```text
http://127.0.0.1:8765
```

## What DEGORA Produces

- `degora_gene_scores.csv`: ranked gene table
- `degora_scores.db`: local SQLite database
- local browser/API: searchable gene evidence
- per-gene evidence rows: which source supported each gene
- source-unit support: avoids overcounting many time points from one paper

The DEGORA score is a transparent prioritization score, not a probability.
Current scoring (`degora_score_v1_2_source_unit_mean`) combines evidence after grouping related contrasts within each independent source unit (`paper_id` when present, otherwise `study_id`). Within a source unit, DEGORA uses mean signed-z/rank/effect summaries rather than choosing the most significant contrast. This prevents multiple time points or repeated probes from one paper from inflating support.

The ranked table also includes interpretation fields for browsing:

- `priority_score`: effect/rank/direction-focused prioritization score.
- `evidence_reliability_score`: support, source quality, direction confidence, and leave-one-source-out rank stability.
- `quality_weighted_degora_rank`: the manuscript-facing default rank, using predeclared source-quality weights.
- `priority_top_percent`: rank as a percent of all scored genes.
- `direction_confidence_index`: direction-consistency index shrunk toward 0.5 when evidence is weak or discordant.
- `loo_rank_stability_score`: whether the rank remains stable after leaving out one source unit at a time.
- `heterogeneity_flag`: descriptive high/moderate/low source-unit heterogeneity review label.

These fields are available in the CSV, SQLite database, JSON API, and local browser.

## Common Fixes

If validation says a DEG column is missing, open the DEG file and copy the column name exactly into `degora_config.xlsx`. Column names are case-sensitive.

If your data are time-course data, do not give every time point a different `source_unit_id` unless they are truly independent studies.

If a paper gives only a significant-DEG list instead of all tested genes, set `table_scope` to `deg_only`. DEGORA will treat genes missing from that table as unreported, not as non-DEG. If the paper reports how many genes were tested, enter that number in `rank_universe_size`; otherwise DEGORA will warn that within-list ranks may be optimistic.

Microarray evidence is supported after assay-aware DEG conversion. Prefer a full limma table when available: set `assay_type=microarray`, `pipeline=limma_microarray`, `source_input_type=limma_full_table`, and record the `platform`, `normalization`, and `probe_collapse` rule. If you only have an already normalized log-scale expression matrix, create a fallback DEG table with:

```bash
PYTHONPATH=outputs/code python outputs/code/scripts/derive_microarray_deg.py \
  --matrix normalized_microarray_matrix.csv \
  --output microarray_deg.csv \
  --summary microarray_deg_summary.json \
  --gene-column gene_symbol \
  --probe-column probe_id \
  --control-samples ctrl1,ctrl2,ctrl3 \
  --treatment-samples drug1,drug2,drug3 \
  --platform GPL570 \
  --normalization RMA_log2
```

Then add `microarray_deg.csv` to the Excel `Contrasts` sheet with `pipeline=welch_microarray_normalized_matrix`. This fallback is exploratory; raw array files or limma full tables are preferred for manuscript-grade microarray analysis. The fallback uses Welch contrasts and does not model paired or family structure unless those designs are handled upstream.

If a row is only a note or placeholder, set `include` to `no`.

## License

Code: MIT (see `LICENSE`; declared in `outputs/code/pyproject.toml`).
