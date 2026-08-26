# DEGORA

DEGORA combines differential-expression result tables into a source-traceable gene ranking and a local browser dashboard. It accepts tabular results from RNA-seq, microarray, and related transcriptomic workflows.

The software can also search public Human or Mouse records, help you inspect available result tables or upstream matrices, and prepare a species-specific analysis. Human and Mouse records are always kept in separate workspaces and separate runs.

## Requirements

- Python 3.10 or newer; automated release tests cover Python 3.10-3.13, and Python 3.12 is recommended
- Ubuntu or macOS; Windows 11 users can run the same Linux workflow in WSL2 Ubuntu
- Internet access only when using public-data search or download features
- Git, unless you use the ZIP download option below

Native Windows PowerShell is not the supported full-workflow environment. On Windows 11, install WSL2 once from an administrator terminal:

```powershell
wsl --install
```

Restart if requested, open Ubuntu, and use the Linux commands below.
The automated release matrix tests native Ubuntu and macOS. WSL2 uses the Ubuntu workflow, but this release is not separately exercised on a hosted Windows runner.

## Install

With Git:

```bash
git clone https://github.com/kangk1204/DEGORA.git
cd DEGORA
```

Without Git, download and unzip either the main-branch ZIP or a tagged release ZIP.
GitHub names those folders `DEGORA-main` and `DEGORA-<version>`, respectively:

```bash
cd DEGORA-main  # main-branch ZIP; a tagged ZIP uses its downloaded folder name
```

Confirm that the interpreter you will use is supported:

```bash
python3 --version
```

The reported version must be 3.10 or newer. Automated release tests cover
Python 3.10-3.13; later versions permitted by the package metadata may not yet
have release-matrix coverage. A virtual environment keeps the
Python version used to create it, so an environment created with Python 3.9
must be removed or renamed and recreated with a supported interpreter.
Do not continue with `python3 -m venv` when the reported version is below 3.10.

macOS (including Apple silicon):

macOS may report Python 3.9 for `/usr/bin/python3` even when a newer Homebrew
Python is installed. Check for the versioned command first:

```bash
python3.12 --version
```

If that command is missing and you use Homebrew, install Python 3.12:

```bash
brew install python@3.12
```

Create the environment with the same versioned command:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python --version
python -m pip install --upgrade pip
python -m pip install -e .
```

Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python --version
python -m pip install --upgrade pip
python -m pip install -e .
```

Windows 11 with WSL Ubuntu:

```bash
python3 -m venv .venv
source .venv/bin/activate
python --version
python -m pip install --upgrade pip
python -m pip install -e .
```

Check the installation:

```bash
degora --version
degora --help
```

### macOS error: `requires a different Python`

If installation reports a message such as `Python 3.9.6 not in '>=3.10'`, the
active virtual environment was created with an unsupported interpreter. From
the repository directory, recreate it with Python 3.12:

```bash
deactivate 2>/dev/null || true
mv .venv .venv-python39-failed
python3.12 -m venv .venv
source .venv/bin/activate
python --version
python -m pip install --upgrade pip
python -m pip install -e .
```

After `degora --version` and the demo below both work, the renamed failed
environment can be deleted.

## One-command quickstart

`scripts/degora_quickstart.sh` performs the whole first run on Ubuntu, macOS
(Intel and Apple silicon), and Windows 11 via WSL2 Ubuntu: it picks a supported
interpreter, creates the virtual environment, installs the package, builds the
demo, chooses a free port, and opens the dashboard.

From a checkout, either path works. `degora_quickstart.sh` in the repository
root forwards to the real script under `scripts/`, so the command is the same
whether or not you remember where the file lives:

```bash
bash scripts/degora_quickstart.sh
bash degora_quickstart.sh          # same thing, from the checkout root
```

Without a checkout, download just that one file and run it; it clones the
repository into `./DEGORA` first:

```bash
curl -fsSLO https://raw.githubusercontent.com/kangk1204/DEGORA/main/scripts/degora_quickstart.sh
bash degora_quickstart.sh
```

Useful options:

```text
--port N        Preferred port (default 8765; the next free port is used when taken)
--dir PATH      Where to place or find the checkout when run standalone
--ref NAME      Branch or tag to run (default: the repository default branch)
--config PATH   Serve your own DEGORA config instead of the bundled demo
--update        git pull an existing checkout before installing
--no-browser    Do not try to open a browser (headless or remote shells)
--no-demo       Skip demo creation and serve an existing database
--demo-dir NAME Demo workspace folder (default: degora-demo)
```

Use `--ref` to review a specific branch or release tag in one command, for
example `bash degora_quickstart.sh --ref main`. It fetches the name from
`origin`, fast-forwards a local copy that is behind, and stops rather than
serving stale code when a local branch of the same name has diverged. `git` is
needed only for cloning, `--ref` and `--update`; an unpacked ZIP folder works
without it. `--config PATH` runs from the config's own folder, so its results
land beside the config rather than inside the checkout, and `--demo-dir` accepts
an absolute path.

The script is safe to re-run: an existing demo workspace is reused, never
deleted, so a config you edited there survives. Use `--demo-dir NAME` if you want
a separate one. It stops with an actionable message when the platform is missing
Python 3.10+, `git`, or the Debian/Ubuntu `python3-venv` package. Press `Ctrl+C` to stop the server. The manual steps below do the same
work one command at a time.

## Run the included demo

```bash
degora demo degora-demo
degora validate degora-demo/degora_demo_config.xlsx
degora run degora-demo/degora_demo_config.xlsx
degora serve degora-demo/results/degora_scores.db
```

The last command starts a local server and prints a browser address, normally `http://127.0.0.1:8765`. Press `Ctrl+C` in the terminal to stop it.

`degora demo degora-demo` writes exactly these starter files:

```text
degora-demo/README.md
degora-demo/degora_demo_config.xlsx
degora-demo/deg_tables/demo_ifn_a_4h.csv
degora-demo/deg_tables/demo_ifn_a_12h.csv
degora-demo/deg_tables/demo_ifn_b_6h.csv
degora-demo/deg_tables/demo_ifn_b_24h.csv
```

`degora run degora-demo/degora_demo_config.xlsx` writes the analysis outputs under `degora-demo/results/`, including:

```text
degora-demo/results/DEGORA_output.xlsx
degora-demo/results/DEGORA_output.manifest.json
degora-demo/results/DEGORA_output.validation.txt
degora-demo/results/degora_gene_scores.csv
degora-demo/results/degora_score_db_summary.json
degora-demo/results/degora_score_metadata.json
degora-demo/results/degora_scores.db
degora-demo/results/degora_source_quality_diagnostics.json
degora-demo/results/degora_source_quality_diagnostics.tsv
degora-demo/results/slice_consensus.csv
degora-demo/results/slice_harmonized.csv
degora-demo/results/slice_metrics.json
```

Each generated result file also has `.source` and `.provenance.json` sidecars that record the command, input checksums, and DEGORA version information. A separate `degora-demo/harmonized/` folder contains the harmonized table copy used by the run.

## Analyze your own tables

If you have a folder of DEG result tables and would rather answer questions than
fill in a spreadsheet:

```bash
degora init my_config.csv --deg-dir path/to/deg_tables
degora validate my_config.csv
degora run my_config.csv
```

`degora init` reads each table and works out what the file can tell it: which
column holds the gene names, the effect size, the p-value and the adjusted
p-value, and whether the table lists every gene tested or only the significant
ones. It shows what it found and what else each column could have been, and asks
only where a file is genuinely ambiguous. Files that are not DEG result tables --
a sample sheet sitting in the same folder -- are recognised and skipped.
Only explicit gene/probe/transcript headers (including recovered R row labels) are
selected automatically. A generic
`ID`, `identifier`, or `feature` column whose values resemble identifiers must be
confirmed, because numeric values alone cannot distinguish Entrez IDs from ranks
or row numbers. Known non-gene labels such as `rank`, `baseMean`, `pathway`, and
`compound` are not offered as gene columns.

Two things are never inferred. For every table it asks whether a positive value
means the gene went **up in the treated samples**. Getting that backwards inverts
every up/down call in the results while leaving them looking entirely reasonable,
so there is nothing later for you to notice. If you answer no, the table is
written into the config **excluded**, with the reason recorded: DEGORA does not
reverse an effect column for you, because that is a correction it cannot verify.
If you are unsure, say so and the table is skipped rather than guessed at. When
the effect column's name does not say it is on a log2 scale (`FoldChange`,
`ratio`), it also asks whether the values **are log2 fold changes**, and shows
their range and how many are negative: a linear fold change (2 = doubling,
0.5 = halving) has no negative values, so DEGORA would read every gene as up.
Answer no and the table is written excluded until it is converted.

Workbooks whose table is not on the first sheet, or whose column names sit
below a title row, are located automatically; the config records `sheet_name`
and `header_row` (the 1-based row holding the column names) so `degora run`
reads the same cells.

Or start from a documented Excel template instead:

```bash
degora template DEGORA_template.xlsx
```

For each contrast, provide at least:

- a stable study and source-unit identifier;
- the result-table path;
- a gene-symbol column;
- a signed effect column such as log2 fold change;
- a nominal p-value column.

Adjusted p-values, sample counts, assay labels, platform details, and source URLs are optional but useful. Validate before running:

```bash
degora validate DEGORA_template.xlsx
degora run DEGORA_template.xlsx
degora serve outputs/results/degora-run/degora_scores.db
```

The output directory contains ranked genes, per-source evidence, provenance sidecars, a SQLite database, and an Excel workbook.

Common exported files are:

```text
DEGORA_output.xlsx
DEGORA_output.manifest.json
DEGORA_output.validation.txt
degora_gene_scores.csv
degora_score_db_summary.json
degora_score_metadata.json
degora_scores.db
degora_source_quality_diagnostics.json
degora_source_quality_diagnostics.tsv
slice_consensus.csv
slice_harmonized.csv
slice_metrics.json
```

### Optional GoldPanel

The Excel template carries a `GoldPanel` sheet for marker genes you already
trust. It changes nothing about the ranking; it adds a recall check to
`DEGORA_output.xlsx` (`Curated_lookup` and the `curated_*` rows of
`Run_summary`) so you can see where your known genes landed.

Two things decide whether the panel is used at all:

- **`locked`** — a row counts when `locked` is `yes` or blank. `locked=no` drops
  that row, and if every row says `no` the whole panel is ignored. The run then
  says so as a non-fatal warning rather than reporting an empty panel silently.
- **`gene_symbol`** — required. Legacy and Excel-damaged symbols are fine:
  `SEPT9`, `MARCH1`, `DEC1`, `9-Sep` and `1-Mar` all resolve to the symbol
  DEGORA scored. `Curated_lookup` shows the symbol you wrote in `gene_symbol`
  and the symbol it resolved to in `resolved_gene_symbol`.

`expected_direction` is recorded beside the observed direction for review and is
not used in scoring. A panel is optional in the ordinary sense: leaving the sheet
empty is a valid, complete run.

## Component ablation and weight sensitivity

`degora ablate` re-scores a finished run with each of the five score components
removed, with source-quality weighting switched off, with sample-size weighting
switched off, and with any weight vectors you name, and reports how far the
primary rank moves under each variant:

```bash
degora ablate degora-demo/results --gold-panel degora-demo/degora_demo_config.xlsx
degora ablate outputs/results/degora-run \
  --weights equal=support=1,direction=1,evidence=1,rank=1,effect=1 \
  --weights effect_heavy=support=0.2,direction=0.2,evidence=0.2,rank=0.1,effect=0.3
```

Every variant scores exactly the same genes over exactly the same source
units: the gene universe and the eligible source-unit set are fixed by the run's
harmonized table, and an ablation only changes how the per-source evidence is
combined. The summary (`ablation/degora_ablation_summary.csv`) reports, per
variant, the Spearman correlation of `quality_weighted_degora_rank` with the
full ranking, the median and maximum rank shift, the top-50 and top-100 overlap
with the full ranking, and recall@k against a GoldPanel or gene list when one is
supplied; `degora_ablation_ranks.csv` holds the rank of every gene under every
variant. Two things to keep in mind when reading it: `evidence_score` and
`rank_score_component` are both derived from the p-value, so removing either
alone moves the ranking less than their joint weight suggests; and
`support_score` is constant in a corpus where every scored gene has the same
number of source units (two source units with `min_studies=2`), where its weight
cannot change the order at all. `without_sample_size_weighting` removes the
per-source-unit sample-size weight; contrasts inside one source unit are still
combined with their sqrt(n) weights, as documented for the collapse rule.

## Search public Human or Mouse records

Search one species at a time:

```bash
degora discover "hypoxia HIF1" \
  --species human \
  --limit 1000 \
  --page 1 \
  --output-dir search-human

degora discover "hypoxia HIF1" \
  --species mouse \
  --limit 1000 \
  --page 1 \
  --output-dir search-mouse
```

The search collects at most 1,000 exact, unique records before sorting and displays 10 rows per page. Detailed file resolution is bounded to the leading 20 records, the first two pages, while later selections are resolved on demand. Rows beyond that report that nothing has been inspected yet - a repository record exists and no file has been opened - until they are prepared. A readiness badge says how many candidate files the estimate rests on, because a record can look ready and hold no usable table: an audit found the top-ranked such candidate had none. `likely_ready` requires a tabular file candidate, and `verified_ready` requires that plus target-species evidence; a bare accession never earns either. Repository records that share a title and have no publication link are marked as a possible single submission, because source units collapse on a shared PubMed ID and an unpublished submission has none. A result is therefore a review queue, not an automatically approved analysis input.

Supplementary tables are inspected as CSV, TSV, TXT and Excel, including gzipped and legacy `.xls` workbooks - the shapes repositories actually serve.

Search exports include JSON, CSV, and Excel snapshots with identifiers, title, authors, journal, year, species evidence, source-unit information, readiness, and provider diagnostics.

Before running an analysis, any selected record that was matched by the species filter rather than a per-record organism check has to be confirmed as that species, and the answer is recorded in the run's metadata rather than assumed.

Species evidence is only as specific as the provider that supplied it. A GEO
record's Series entry lists the organisms of its samples and platforms; DEGORA
reads them, so a record whose samples name exactly the requested species is
`target_species_verified` (and `verified_ready` when it also carries a tabular
file candidate), and a record whose samples span two organisms is quarantined
out of a species-specific preparation at search time. A record found only
through the literature search carries the organism filter that produced the
search and nothing more, so it is labeled `query_constrained` rather than
checked, and mixed-species quarantine cannot apply to it. Confirm the species of
a `query_constrained` record yourself before activating it; the browser asks for
that confirmation only for such records.

You can also use the local browser:

```bash
degora serve path/to/degora_scores.db
```

Open the **Discover** tab, choose **Human** or **Mouse**, search, review the candidates, and select the records to prepare. Author-provided result tables are preferred. Matrix-derived fallback analysis requires explicit group mapping, contrast direction, and biological-replicate confirmation.

### What DEGORA prefers, in order

DEGORA's score is built from published evidence, and some evidence carries more
of it than other evidence. The order below is what the search sorts by within a
readiness tier, what preparation puts first when a series ships several files,
and what a config should aim for. It is an order of evidence quality, not of
frequency: a series usually ships one of each.

1. **The authors' own results table, covering every gene tested** - a
   `full_results` table with a log2 fold change, an unadjusted p-value and, if
   available, an adjusted one, keyed by gene symbol. The authors' statistics
   reflect their design, covariates and model; a full table gives
   `normalized_rank` its true denominator.
2. **The authors' results table for significant genes only** (`deg_only`) -
   the same statistics, but the rank universe is unknown unless
   `rank_universe_size` is given, so its normalised ranks are optimistic and its
   Stouffer p-values are ranking aids rather than genome-wide inference.
3. **A raw count matrix** - the least processed matrix. DEGORA normalises it
   itself and derives a contrast with a Welch test once the reader has assigned
   control and treatment samples and confirmed they are biological replicates.
4. **A log2-normalised matrix** (TMM, VST, rlog) - already normalised; used as
   given after the scale is confirmed.
5. **A linear normalised matrix** (FPKM, TPM, CPM) - transformed with
   log2(x + 1) before testing.

Across all of these: a contrast whose group sizes are known is weighted
`min(sqrt(n_ctrl + n_treat), 4)`, and one whose sizes are unknown is weighted
1.0 - below a two-sample contrast. Gene identifiers join across studies only
when they are written the same way, so symbols join the most, and time points
or normalisations of one series count as one source unit, never as several.

At search time nothing has been opened, so the results table shows this order
as an estimate from file names ("likely author DEG table", "likely raw count
matrix") and sorts by it inside each readiness tier. Preparation opens the
files and ranks what it finds; when one series carries several files of the
same samples, the preferred one is shown first and the rest sit behind one
collapsed line.

Browser discovery analysis uses a fixed `min_studies=2`: a gene must retain
eligible evidence from at least two independent source units. The CLI command
`degora discovery-analyze --min-studies N` exposes the validated threshold when
an intentionally different replication floor is required. `--min-studies 1`
scores genes from a single source unit: the run warns that such a ranking shows
no replication and is exploratory prioritisation, not replicated evidence.

Every `POST` to `/api/discovery/...` must carry the header `X-DEGORA-Action: 1`
(the page adds it; a `curl` call without it is refused with a message naming
the header). It is a cross-site request forgery guard for the local server.

The review panel asks each prepared table only for the confirmations that actually apply to it, and keeps the settings most tables never touch behind a collapsed **Advanced settings** panel that opens by itself if any of them is already set. Where a linked series reports how many samples it holds, that total is shown beside the group-size boxes and the two numbers you enter are checked against it; the split between groups is never guessed, because a results table has one row per gene and the number feeds the source weight directly. A table whose columns were recognised and left alone is asked one thing; a table whose effect column does not say it is log2, or whose adjusted p-value is standing in for a raw one, is asked about that as well. Contrast direction is the exception: it is asked for every table, because reversing it inverts every up/down call while leaving results that look entirely reasonable.

The browser opens on **Evidence atlas** when the database already holds scored genes, which is what `degora serve` after a run is for, and on **Discover** when it does not. Either tab is one click away at any time. Search results are globally ranked before the first page of 10 is shown. **Narrow these results** filters the records already found by title, author, journal or year without running a new search, and reports how many of the total match. The compact table displays publication metadata, linked-data availability, estimated DEG-input readiness, and an **Inspect** action. **Run separate Human + Mouse searches** launches two independent searches; it never pools their records or scores.

While a search runs, the results panel shows the current stage, a completion percentage, and elapsed time, driven by the search job's own progress reporting; `GET /api/discovery/jobs/{job_id}` exposes the same `progress` and `message` fields. At most 20 publications can be selected at a time, and that limit applies across pages: the selection counter reports how many of them are on the page you are looking at, and rows that cannot be selected are shown disabled with the reason.

Starting `degora serve` creates a `degora_discovery/` workspace beside the database. It contains the local search cache (`discovery.sqlite3`), a process lock while the server is running, downloaded candidate files, and any discovery runs. Keep or remove that workspace according to your local data-retention policy after stopping the server.

To reproduce an approved prepared run from saved JSON artifacts:

```bash
degora discovery-analyze BUNDLE_JSON SELECTIONS_JSON \
  --species human \
  --output-dir analysis-human
```

Use a separate command and output directory for Mouse.

## NCBI identification and API key

NCBI identification is optional but recommended for repeated searches:

```bash
export NCBI_EMAIL="you@example.org"
export NCBI_API_KEY="your-key"
```

Create or sign in to an NCBI account and open <https://account.ncbi.nlm.nih.gov/settings/> to generate an API key. Do not place the key in a workbook, command history, source file, or shared archive.

`NCBI_EMAIL` identifies a responsible contact for NCBI requests. `NCBI_API_KEY` raises the allowed request rate for your account. DEGORA reads these values from your shell environment when it contacts NCBI; the synthetic demo does not need them.

## Reproduction boundary

The included synthetic demo is numerically and semantically reproducible from this repository across the supported Python environments. Repeating a run over the same inputs **in the same environment** reproduces `degora_gene_scores.csv`, `degora_scores.db`, and `DEGORA_output.xlsx` byte for byte: the workbook's own timestamps and its archive member timestamps are pinned rather than taken from the clock. Across *different* dependency versions the results are numerically identical but not byte-identical: NumPy and SciPy differ in the last one or two digits of some floating-point fields, which changes the text a float is written as. Ranks, tiers, directions and every displayed value are unaffected. Compare runs from different environments by value at a sane precision, not by checksum; the `.provenance.json` sidecars checksum the *inputs*, which are stable, and the environment that produced each output is recorded beside it. Larger external datasets are not bundled here and must be obtained from their original providers before they can be analyzed.

## Interpretation boundaries

- DEGORA prioritizes genes; its scores are not posterior probabilities.
- The primary output rank is `quality_weighted_degora_rank`. The earlier
  `degora_rank` column is the unweighted audit/reference lane, even when a CSV
  viewer displays it first.
- Gene symbols are resolved to one current symbol before anything is compared.
  Excel date damage is undone (`9-Sep` -> `SEPTIN9`) and legacy symbols are
  updated (`SEPT9` -> `SEPTIN9`, `MARCH1` -> `MARCHF1`, `DEC1` -> `BHLHE40`), so
  two tables that spell the same gene differently still count as one gene. The
  same resolution is applied to the optional GoldPanel and to browser/API
  lookups, so searching `SEPT9` finds the gene DEGORA scored. The label each
  source table actually carried is kept in the `input_gene_label` column of
  `slice_harmonized.csv`, and `degora run` reports how many symbols it changed.
  A version suffix is removed only from accession-shaped identifiers
  (`ENSG00000141510.16`, `NM_000546.5`, an Entrez ID exported as `7157.0`);
  a dotted symbol such as `NKX2.5` is kept as written, and matches a partner
  table only when that table writes it the same way (`NKX2-5` is a different
  label). `DEC1` is also the previous symbol of `DELEC1`; DEGORA maps it to
  `BHLHE40`, and `input_gene_label` keeps the original so the choice is visible.
- `validate` and `run` refuse an effect column whose values have the shape of a
  linear fold change (no negative values, values below 0.5 and above 1), refuse
  a p-value written as a bound (`<1E-16`, `p<0.05`) rather than dropping the
  row, refuse a mapping onto a header the table carries twice, and refuse one
  result table declared under two different source units. Rows dropped for a
  missing gene, effect or p-value are always reported, whatever their share.
- Source units whose log2 fold changes run against the rest of the corpus are
  flagged `source_direction_conflict_flag` in the source-quality diagnostics and
  named in the run warnings as a possibly reversed contrast. A unit qualifies
  when at least half of its pairwise comparisons are below -0.10 and
  significantly negative **and** its own median pairwise correlation is below
  -0.10. The median condition is what keeps a well-formed unit out of it: with
  two comparisons, "at least half" is satisfied by the single comparison against
  a reversed partner, so in a three-source corpus with one reversed contrast
  every unit used to be named. Where the flagged units are not a minority the
  corpus is split rather than one source being wrong, and the warning says so
  instead of accusing each unit in turn. The flag changes no weight and no
  rank.
- Related contrasts are collapsed by source unit before cross-source aggregation.
- `time_course_mode` chooses which contrasts of a source unit are kept before that collapse, and every row sharing a source unit must use the same mode. `mean` keeps all of them; `early` and `late` keep all gene rows at that source unit's globally smallest and largest numeric `duration_h`, and every active row of such a unit must carry `duration_h` as a plain number of hours (`0.5`, `24`; `30min` or a blank cell is refused at validation); `peak_mean` keeps each gene's strongest half by `|signed_z|`, at least two. `peak_mean` selects on statistical strength, not effect size: `signed_z` is derived from the p-value, so a time point with a large fold change but a weak p-value is not the peak.
- A row with `pvalue = 1` or zero effect is neutral evidence and does not contribute directional signed-z support. Values below `1e-300` are floored and reported in the run warnings.
- Stouffer p-values from DEG-only or significance-filtered tables are ranking aids, not calibrated genome-wide inferential p-values or false-discovery rates.
- `heterogeneity_i2` is a sample-size-weighted descriptive dispersion index, not calibrated Higgins I-squared. The heterogeneity-adjusted Stouffer fields are screening aids only.
- Random-effects effect-size intervals are descriptive when few source units contribute; in particular, Hartung-Knapp-Sidik-Jonkman intervals with fewer than three source units are unstable.
- Leave-one-source-out stability is a `priority_rank` diagnostic over global source-unit omission folds. Each fold applies the same `min_studies` eligibility rule and deterministic tie-break as the full priority lane. Median, IQR, and top-N fractions summarize rank-evaluable folds only; the stability score still treats ineligible folds as negative evidence. If no fold keeps a gene eligible, numeric LOO fields are reported as unavailable rather than as zero; the conditional reliability summary then uses the other three mandatory diagnostics with their weights renormalized. These LOO columns are therefore nullable. `evidence_reliability_components_used` records whether three or four diagnostics contributed. Compare reliability values across runs or corpora only when the contributing diagnostic count and LOO eligibility conditions match. The default browser/API ordering remains `quality_weighted_degora_rank`, and reliability does not determine that ordering.
- Missing group sizes receive the documented neutral replicate multiplier of 0.75. `validate` rejects a zero, negative, fractional or non-numeric group size, and one above 10,000 biological replicates for a single contrast — the same bound the browser's review panel applies — so the 0.35 zero-count branch is reachable only through the Python API. The cap matters because the per-contrast weight is `min(sqrt(n_ctrl + n_treat), 4)`: a typing slip that turns 3 into 999999 saturates that weight at its ceiling rather than producing an obviously wrong number.
- `sign_convention` records provenance only. DEGORA does not infer or automatically reverse the input effect direction, so the supplied effect must already represent treatment relative to control.
- `lfc_scale` is optional, and `log2` is the only value it takes. Write it to state that `lfc_column` already holds log2 fold changes when its header does not say so (a column named `lfc`, `beta` or `effect`). Any other value — `log10`, `ln`, `linear` — is refused at validation rather than treated as a blank cell, because DEGORA never converts a scale and a declaration it ignores is worse than none. Without it, a column with no negative values and values on both sides of 1 is refused as a linear fold change; with it, that shape is reported and the run proceeds. `degora init` writes it from the answer to its scale question. DEGORA never converts a scale itself. The shape checks need at least 10 numeric values; a smaller table is taken at its declared or named scale, and the run says so.
- `direction_confirmed`, `biological_replicates_confirmed`, the sample groups of a matrix contrast and `lfc_scale` are the reviewer's statements. DEGORA records them with every run (`reviewer_attestations` in the discovery run summary) and cannot verify them; a wrong attestation produces a ranking that looks ordinary.
- Public-data fallback uses a documented Welch workflow and does not automatically correct study-level batch or condition confounding.
- Result-table semantics, contrast direction, group labels, and species must be reviewed before activation.
- The Search workflow keeps Human and Mouse in separate workspaces and never pools them, and a prepared discovery bundle is refused if its species does not match the run. Scoring itself matches on gene symbol and is **not** species-specific: a hand-written config naming sources from two species produces one pooled ranking, in which those sources can satisfy the `min_studies` replication rule between them. The run warns when it sees more than one `species` value, and every evidence row records the species it came from.

## Development checks

```bash
python -m pip install -e ".[dev]"
make check
make smoke
```

## Release notes

Every release is recorded in [CHANGELOG.md](CHANGELOG.md), newest first, with
what changed and what the previous behaviour was.

## License

MIT
