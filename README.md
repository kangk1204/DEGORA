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
serving stale code when a local branch of the same name has diverged.

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

One thing is never inferred. For every table it asks whether a positive value
means the gene went **up in the treated samples**. Getting that backwards inverts
every up/down call in the results while leaving them looking entirely reasonable,
so there is nothing later for you to notice. If you answer no, the table is
written into the config **excluded**, with the reason recorded: DEGORA does not
reverse an effect column for you, because that is a correction it cannot verify.
If you are unsure, say so and the table is skipped rather than guessed at.

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

Species evidence is only as specific as the provider that supplied it. A record
found through a repository that reports per-sample organisms carries a checked
species label, and a record whose samples span two organisms is quarantined out
of a species-specific preparation. A record found only through the literature
search carries the organism filter that produced the search and nothing more, so
it is labeled `query_constrained` rather than checked, and mixed-species
quarantine cannot apply to it. Confirm the species of a `query_constrained`
record yourself before activating it.

You can also use the local browser:

```bash
degora serve path/to/degora_scores.db
```

Open the **Discover** tab, choose **Human** or **Mouse**, search, review the candidates, and select the records to prepare. Author-provided result tables are preferred. Matrix-derived fallback analysis requires explicit group mapping, contrast direction, and biological-replicate confirmation.

Browser discovery analysis uses a fixed `min_studies=2`: a gene must retain
eligible evidence from at least two independent source units. The CLI command
`degora discovery-analyze --min-studies N` exposes the validated threshold when
an intentionally different replication floor is required.

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
- Related contrasts are collapsed by source unit before cross-source aggregation.
- `time_course_mode` chooses which contrasts of a source unit are kept before that collapse, and every row sharing a source unit must use the same mode. `mean` keeps all of them; `early` and `late` keep all gene rows at that source unit's globally smallest and largest numeric `duration_h`; `peak_mean` keeps each gene's strongest half by `|signed_z|`, at least two. `peak_mean` selects on statistical strength, not effect size: `signed_z` is derived from the p-value, so a time point with a large fold change but a weak p-value is not the peak.
- A row with `pvalue = 1` or zero effect is neutral evidence and does not contribute directional signed-z support. Values below `1e-300` are floored and reported in the run warnings.
- Stouffer p-values from DEG-only or significance-filtered tables are ranking aids, not calibrated genome-wide inferential p-values or false-discovery rates.
- `heterogeneity_i2` is a sample-size-weighted descriptive dispersion index, not calibrated Higgins I-squared. The heterogeneity-adjusted Stouffer fields are screening aids only.
- Random-effects effect-size intervals are descriptive when few source units contribute; in particular, Hartung-Knapp-Sidik-Jonkman intervals with fewer than three source units are unstable.
- Leave-one-source-out stability is a `priority_rank` diagnostic over global source-unit omission folds. Each fold applies the same `min_studies` eligibility rule and deterministic tie-break as the full priority lane. Median, IQR, and top-N fractions summarize rank-evaluable folds only; the stability score still treats ineligible folds as negative evidence. If no fold keeps a gene eligible, numeric LOO fields are reported as unavailable rather than as zero; the conditional reliability summary then uses the other three mandatory diagnostics with their weights renormalized. These LOO columns are therefore nullable. `evidence_reliability_components_used` records whether three or four diagnostics contributed. Compare reliability values across runs or corpora only when the contributing diagnostic count and LOO eligibility conditions match. The default browser/API ordering remains `quality_weighted_degora_rank`, and reliability does not determine that ordering.
- Missing or non-numeric group sizes receive the documented neutral replicate multiplier of 0.75; an explicit zero receives 0.35.
- `sign_convention` records provenance only. DEGORA does not infer or automatically reverse the input effect direction, so the supplied effect must already represent treatment relative to control.
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

### Unreleased

Fixes from a first-run review of v0.4.17 by a reader following this README with
no prior knowledge of the tool.

- Gene-symbol resolution is now applied in every direction, not only on input.
  The optional GoldPanel and browser/API gene lookups resolve legacy and
  Excel-damaged symbols the same way source tables do, so a panel written as
  `SEPT9` no longer reports the gene as absent while the run ranks it as
  `SEPTIN9`, and `/api/genes/SEPT9` answers instead of returning 404. The label
  each source table carried is preserved in the new `input_gene_label` column of
  `slice_harmonized.csv`, and `degora run` reports how many symbols it changed.
- `degora serve` refuses a path that is not a DEGORA score database instead of
  binding and then failing every API call with a 500 behind a dashboard that
  loaded normally.
- Validation reports every active row whose source table is missing, rather than
  stopping at the first, so leaving the template's example rows in place is one
  fix instead of one fix per row.
- The missing-required-column error lists the headers the file actually had.
  Optional catalog defaults are injected before that check, so the message used
  to name `source_unit_id` as available in the same breath as calling it
  missing. A file whose headers look like a DEG results table now says so.
- Search failures name the cause and say what to do about it, in the same
  `Problems:` / `How to fix:` shape as every other DEGORA error. A blocked proxy
  or an offline machine is no longer reported as "retry later", and transport
  failures no longer surface a bare Python exception class name.
- Duplicate rows for one gene whose log2 fold changes disagree in sign now
  populate `gene_symbol_collapse_warning`. Collapsing them is unchanged; leaving
  a directional conflict unreported was the problem.
- A GoldPanel whose rows are all `locked=no` is reported as a non-fatal warning.
  `locked` is now documented in the template's `ColumnGuide` sheet, in the
  `GoldPanel` sheet note, and in this README, instead of only inside the example
  row the template tells you to replace.
- The browser opens on **Evidence atlas** when the database holds scored genes,
  which is what `degora serve` after a run is for. Discover remains the landing
  view for a database with nothing scored yet.
- `degora_quickstart.sh` in the repository root forwards to
  `scripts/degora_quickstart.sh`, so the command works from a fresh checkout's
  root as well.
- `degora run` names the primary rank column in its completion summary, and the
  reproduction boundary now states that results are numerically identical but
  not byte-identical across different NumPy/SciPy versions.
- GitHub Actions now uses the official Node 24 action releases pinned to immutable
  commit SHAs instead of floating Node 20 major tags.
- The browser documentation now states its fixed two-independent-source-unit
  analysis floor. A regression test locks that API policy, and the `run_slice`
  directory-creation comment now matches the actual failure boundary.

### 0.4.17

The primary score contract is unchanged: `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`. Version 0.4.17 supersedes the published
v0.4.16 beginner-input and optional-GoldPanel behavior described below without
changing an unchanged valid config's scores.

Guided initialization no longer promotes value-only `rank`, row-number, count,
statistic, pathway, metabolite, cell-line, or compound columns to gene identifiers.
Generic identifier headers remain recoverable through an explicit column question;
only gene/probe/transcript-specific headers are automatic. Low-level aggregation
also rejects an unknown `time_course_mode` instead of silently treating a typo as
`mean`.

Catalog normalization now fills an empty legacy `paper_id` column from
`source_unit_id` without a pandas 3 dtype error.

CSV/TSV outputs now neutralize spreadsheet formula-like text. Generated files
that are reused as analysis inputs carry a digest-checked provenance marker, so
DEGORA reverses exactly one guard during ingestion and preserves the original
identifier in scoring and SQLite. If that marker is missing or its digest is
invalid, apostrophe-guarded formula-like text is rejected as ambiguous instead of
silently changing a scientific identifier. Early/late time-course runs now record per-source
row/gene retention and warn below the recorded 50% retention threshold; legacy
temporal aliases that fill blank canonical cells are disclosed before validation
and execution.

Custom ablation weights are stored as the normalized floating-point values that
passed validation, so numeric strings cannot survive construction and fail later
inside weighted scoring. The workbook dictionary now labels the legacy
`direction_posterior_mean` field as a shrinkage index rather than a calibrated
posterior probability.

An optional GoldPanel that cannot be read or lacks the required `gene_symbol`
column is no longer treated as though no panel was supplied. Validation and run
commands emit one sanitized, non-fatal warning, while run metrics, manifests,
and workbooks retain the `read_error` or `invalid` status and a bounded reason.
Blank GoldPanel cells do not create a synthetic `NAN` gene.

### 0.4.16

Known issue in the published v0.4.16 tag: guided initialization could accept an
ordinal or other non-gene column such as `rank` as the gene identifier, allowing
two plausible-looking tables to finish successfully with row numbers ranked as
genes. Invalid or unreadable optional GoldPanel content could also be reported as
absent. Version 0.4.17 refuses or explicitly questions the identifier mapping and
surfaces the GoldPanel status; the historical v0.4.16 tag remains unchanged.

The primary default-`mean` score contract is unchanged, so `SCORE_VERSION`
remains `degora_score_v1_2_source_unit_mean`. Configs that explicitly use
`early` or `late` can change: those modes now select the source unit's globally
earliest or latest numeric duration before gene-level collapse, instead of
silently selecting a different earliest/latest row for each gene. One source
unit must declare one normalized temporal mode, and the legacy `temporal_mode`
header is promoted when the canonical column is absent or blank.

Guided config creation now publishes `.csv` and `.xlsx` through an atomic sibling
file, so a forced overwrite that fails cannot corrupt the existing config. Other
output suffixes are refused rather than receiving CSV bytes under a misleading
name. Discovery ordering treats NaN and infinite relevance ranks as missing and
therefore remains independent of provider row order.

Background discovery jobs and their search snapshots now converge to terminal
failure or interruption states even for hard worker failures, server shutdown,
provider failure, and final snapshot-write failure. Invalid custom ablation
weights are rejected before scoring and are snapshotted so caller mutation cannot
change a recorded ablation. The public `discovery-analyze` CLI also validates its
JSON artifacts and `min_studies` without exposing a traceback for reader-correctable
input errors.

### 0.4.15

The score contract is unchanged. This patch release tightens the v0.4.14
beginner, discovery, cancellation, privacy, and provenance paths that reviewers
exercised after the cancellation and guided-init features landed.

`degora init` is stricter about what can be treated as a DEG table. Generic
sample identifiers such as `sample_id` are not accepted as gene columns, binary
significance flags are not accepted as p-values, p-value columns are not offered
as effect-size columns, and common supplementary-table headers such as
`fold-change`, `p-value`, and `q-value` are recognised without role collisions.
If no table is confirmed, the command exits with a clean beginner-facing error
instead of a traceback.

Discovery now separates inspectable result-table candidates from mere linked
data. Search readiness requires a CSV/TSV/Excel/archive-like candidate rather
than only a repository accession, shared-submission warnings are exported with
the publication review tables, and strict assay suffixes such as RNA-seq,
ATAC-seq, CUT&Tag, and ChIP-seq can mark unpublished SuperSeries arms without
collapsing species or time-point suffixes. Persisted search assessments identify
this changed interpretation as assessment version 2, and upstream matrices using
`ensembl_gene_id_version` are recognised.

Cancellation publication and browser status are made consistent: a cancel that
wins prevents a completed search or prepared bundle from being published, while
a cancel that loses to an already completed job adopts the completed result
rather than reporting a false cancellation. Repeated cancel requests now report
the already-cancelled state, and cancelled job bookkeeping is released after
terminal cleanup. If cancellation persistence itself is interrupted, the pending
manager state is cleared and waiting publication work is released rather than
left blocked.

Network-facing API responses and diagnostics redact local paths and
formula-like spreadsheet text more defensively. Source checkouts now report a
dirty code revision, for example `0.4.15 (e2ab3ae-dirty)`, so local patched runs
are not mistaken for the exact public tag. Installed copies cannot inherit the
revision or dirty state of an unrelated parent Git repository.

### 0.4.14

Known issue: the cancellation guarantee described below was not linearizable at
the final search/bundle publication boundary in the published v0.4.14 tag. A
narrow race could expose a completed side effect after the job endpoint reported
`cancelled`. Version 0.4.15 fixes that boundary and adds
deterministic race and failure-path tests.

The score contract is unchanged, and nothing about an existing config or run
behaves differently. This release lets you stop a running job, and answers what
happened when `degora init` was pointed at tables nobody had curated first.

**A search or a preparation can be stopped.** Either issues dozens of paced
requests to public repositories and can take minutes, and there was no way to
end one: closing the tab left the worker downloading, and the only cancellation
that existed stopped every job at once when the server shut down. The Stop
button on the progress card cancels that one job, server-side, at its next
progress report.

What stopping does is stated rather than left to be inferred. Files already
downloaded are kept and a later run reuses them; what is guaranteed is that a
cancelled job never records a result, so partial work cannot be read as a
finished search. A stopped search says it was stopped instead of showing the
empty state, which would have been a claim about the query rather than about
your own action. Pressing Stop a moment too late is answered with "the job had
already finished", and the result stays readable.

**A run that cannot produce anything says so before it is spent finding out.**
Two failures each cost a full run to discover, and both were visible beforehand.

Sources written in different gene identifier spaces share no genes. Two series
downloaded for one topic wrote their gene columns one in Ensembl IDs and one in
symbols; the config validated, the run took its full time and scored zero genes.
`degora init` reads every table anyway, so it now names the identifier space of
each one - gene symbol, Ensembl, RefSeq, Entrez, Affymetrix probe - and warns
before writing the config when the confirmed tables do not agree, naming which
files use which convention. And a config with fewer independent source units
than `min_studies` cannot score a single gene; both `degora init` and
`degora validate` now say that instead of leaving it to the run.

**Questions offer columns that could hold the role.** A table with 43 columns,
32 of them per-sample expression values, asked which one held the gene names and
listed all 43. Nothing in that list was wrong, but a list that long is not a
choice. Candidates are drawn from what the values allow rather than from the
whole header, so that table now offers two. A DESeq2 export carrying only
baseMean, a fold change and padj was offered baseMean as a p-value candidate; a
p-value lies in [0, 1], so only padj is offered now - and the prompt says the
table has no unadjusted p-value, which is worth choosing knowingly.

**R row labels are recognised.** `write.csv` writes gene identifiers as an
unnamed index, which DEGORA recovers under `row_name` - and its own header
classifier had no entry for that name, so the guided setup asked which column
held the gene names for a file whose gene column it had just built itself.

Verified on Ubuntu (Python 3.10-3.13), macOS, a wheel installed without the
development extra, and both the pinned lower-bound and upper-bound dependency
sets.

### 0.4.13

The score contract is unchanged, and nothing about an existing config or run
behaves differently. This release adds one command.

**`degora init` builds a config by asking questions instead of handing you a
spreadsheet.** Point it at a folder of DEG result tables and it reads each one:
which column holds the gene names, the effect size, the p-value and the adjusted
p-value, and whether the table lists every gene tested or only the significant
ones. It shows what it found and what else each column could have been, and asks
only where a file is genuinely ambiguous - two plausible effect columns, or an
effect column whose name does not say it is on a log2 scale. A file that is not a
DEG results table, such as a sample sheet in the same folder, is recognised and
skipped rather than walked through question by question.

The column and scope inference reuses the classifier and the scope assessment the
discovery path already applies to the same kind of file, rather than growing a
second opinion about it.

**Contrast direction is asked for every table and never inferred.** Reversing it
inverts every up/down call in the results while leaving them looking entirely
reasonable, so there is nothing downstream to catch it. Answer that a positive
value means up in the *control* group and the table is written into the config
excluded, with the reason recorded on the row: DEGORA does not reverse an effect
column here any more than it does anywhere else, because that is a correction it
cannot verify. Answer that you are unsure and the table is skipped - pressing
enter never stands in for yes on that question.

### 0.4.12

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces the
same gene scores and evidence as v0.4.11.

**Re-running the quickstart no longer deletes your demo workspace.** The script
removed `degora-demo` before rebuilding it, taking any config you had edited there
and any results you had kept with it, while this page called re-running safe. An
existing workspace is now reused and re-run in place; `--demo-dir NAME` gives you
a separate one.

**`.xls` files open on an ordinary install.** DEG tables ending in `.xls` were
accepted, but the reader pandas needs for them shipped only in the development
extra, so a plain install advertised a format it could not read. `xlrd` is a
runtime dependency now.

**An identifier cannot hold the character that joins identifier lists.**
`contributing_study_ids`, `source_units` and `contributing_source_paths` are
semicolon-joined, so a `study_id` of `A;B` made those lists ambiguous and inflated
`n_contrasts_observed` - three contrasts published as four - while the scores,
which count distinct values rather than splitting text, stayed correct. The
preflight now rejects the delimiter in `study_id`, `paper_id` and
`source_unit_id`, and a source unit derived from a DOI that contains one is
sanitised rather than carried through.

**Opening a different analysis clears the gene filters.** Switching context reset
the rows, the page and the detail pane but not the gene search box, so a filter
typed against one run silently applied to the next: an analysis of 11,886 genes
could open showing the nine that matched, which reads as an analysis that found
nine.

**A ZIP member cannot expand past its cap.** Every archive size limit was
computed from the uncompressed size the archive itself declares, and
`ZipFile.read` decompresses a whole member before it validates the checksum - so a
member declaring 1 KiB and holding a gigabyte of zeros was fully expanded in
memory first. A 1 MB download could force roughly 2 GB of allocation while passing
the member, total, count and depth caps. Members are streamed and stopped at the
cap now.

**A source table that is not a regular file is refused.** `exists()` is true for a
FIFO, and the reader then waited for a writer that never came: no output, no CPU,
no return - indistinguishable from a hang.

**One oversized field no longer ends a preparation with a traceback.** A field past
csv's 128 KiB limit raised an error that none of the preparation handlers caught,
and these files come from public repositories.

**Two runs cannot share one output directory.** The harmonized table is written
seconds into a run and the database tens of seconds later, so two runs pointed at
the same `--output-dir` could interleave and leave one run's contrast table beside
the other's gene scores. Both exited 0, both artifact sets verified against their
own sidecars, and the sidecars were byte-identical because they record only the
command - so nothing downstream could tell the halves came from different runs.
The default output directory is a fixed path, so no flag was needed to hit it. A
run now holds its output directory and a second one is told to wait or use its own.

**A contrast missing a group size is no longer weighted as if it had one.** The
Stouffer weight tested only that the two group sizes summed above zero, so a
contrast with no controls at all drew the weight of a study of its treatment group
- `(0, 5)` earned `sqrt(5)` - and a negative count passed through. Both group sizes
must now be present and positive. The CLI already rejected these values, so this
changes results only for callers using the Python API directly.

**Species are not pooled silently.** The Search workflow keeps Human and Mouse
apart, but scoring matches on gene symbol and is not species-specific: a
hand-written config naming one human and one mouse source produced a pooled
ranking in which those two satisfied the `min_studies` replication rule between
them, with nothing said. Both `validate` and `run` now warn, and the interpretation
boundary above states the real scope.

**`.xls` files open on an ordinary install, and an unreadable one says why.** The
reader pandas needs for legacy Excel shipped only in the development extra, so a
plain install advertised a format it could not read - and reported a perfectly
valid workbook as damaged. A missing engine, a damaged OLE2 file, a renamed text
file and a ZIP that is not a workbook are four different messages now. CI installs
the built wheel with no development extra and runs a real `.xls` through
`validate` and `run`, because the unit tests ran with the extra installed and hid
this for two releases.

**The heterogeneity note no longer claims a bias direction.** `heterogeneity_i2`
was documented as positively biased. Q is not chi-square distributed here and its
scale is arbitrary, so raw `(Q-df)/Q` is frequently negative and clamped to 0,
which makes the reported index effectively bimodal rather than biased in either
direction.

**Three formulas now say what the code does.** The per-contrast Stouffer weight is
published in the score metadata for the first time, including that an unstated
group size scores below a two-sample contrast rather than neutrally. The
direction-confidence tie branch, which credits each source unit one half and
returns 0.5 rather than the 0.25 its formula alone gives, is stated. So is the fact
that the zero-replicate quality branch is unreachable from the CLI, and that an
HKSJ interval at k=2 is built on a t critical value of 12.71 and is uninformative
in practice rather than merely "unstable".

**`peak_mean` says what it is the peak of.** It keeps the strongest half of a
source unit's contrasts by `|signed_z|`, and `signed_z` comes from the p-value, so
a time point with a large fold change but a weak p-value is not the peak. That was
true but undocumented; it is now stated in the score metadata, the catalog help,
the config template, the workbook dictionary and the interpretation boundaries
above, and pinned by a test.

### 0.4.11

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces the
same gene scores and evidence as v0.4.10.

**Species provenance now reaches the prepared bundle and its audit.** v0.4.10 made
the `query_constrained` label correct in the search snapshot, but preparing a
selected publication dropped `species_decision`, `species_evidence` and
`target_species_verified` on the way into the bundle. The species gate reads those
fields off the search record, so preparation kept working while the archived
`discovery_audit.json` -- the document a reviewer opens -- reported none of them,
and the prepare view renders its species provenance line only when the decision
survives. All three are now carried through, and the regression test asserts the
whole chain rather than the helper that copies the fields.

**The unreadable-workbook message matches the extension.** Only `.xlsx` is a ZIP
container; a legacy `.xls` is an OLE2 compound file, so the message named the
wrong container for one of the two extensions it covered.

### 0.4.10

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces the
same gene scores and evidence as v0.4.9. This release makes two things v0.4.9
announced actually reach the data.

**The `query_constrained` species label now survives the pipeline that produces
it.** v0.4.9 added the label and a test that passed, but a live search returned it
for none of 200 records. Normalizing a record's species evidence was not
idempotent: preparation writes the record's display species from the evidence it
has just normalized, and normalizing the same record again folded that copy back
in as a second, independent-looking signal - enough to send every literature-only
record back to `target_species_likely`. Normalization now leaves an
already-normalized record alone, and a record whose only signal is the organism
filter is recognized whichever provider echoed it, not only PubMed. The tests
drive the preparation path rather than a hand-built evidence list, which is why
the gap was invisible.

**Row-loss warnings count rows, not reasons.** The per-reason counts overlap - one
row can be missing its effect and its p-value - and they were summed, so a warning
could report more rows dropped than the table held (one corpus produced "127.4%"),
and a table losing 4% of its rows could be pushed over the 10% threshold and warn
as 12%. The warning now counts the distinct rows the validity mask removed and
keeps the reasons as an explicitly non-exclusive breakdown.

**Duplicate collapse is no longer counted as unusable rows.** The dropped-row count
was taken after duplicate gene symbols were collapsed, so an ordinary probe-level
table appeared to have lost rows its source had in fact supplied. It is settled
from the validity mask before collapse, and rows merged by collapse are reported
separately as `n_rows_merged_by_gene_collapse`.

**An unreadable source table is described, not quoted.** v0.4.9 replaced the
internal pandas message for the config file but not for the DEG tables the config
points at - a renamed CSV or a truncated supplementary download being the commoner
mistake of the two.

### 0.4.9

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces the
same gene scores and evidence as v0.4.8.

**A source that loses rows now says so.** Duplicate gene symbols have always been
reported when they were collapsed. Rows dropped outright -- no gene symbol, no
numeric effect, no numeric p-value -- were not, so a table whose effect column
exported as text lost half its rows between the file and the ranking in silence.
Losses above a tenth of the input now raise a warning that names the responsible
column and shows the values it could not read, and `slice_metrics.json` records
the per-source count.

**`validate` checks that the effect and p-value columns hold numbers.** The
p-value column was checked for being inside [0, 1], which only sees values that
parsed at all. A column of `UP`, `n/a` or a spreadsheet error value passed the
preflight and lost its rows during the run instead.

**Species evidence says what was actually checked.** Mixed-species quarantine
needs two organism labels on one record, and only a provider reporting per-record
organisms can supply them. A record found solely through the literature search
carries the organism filter that produced the search and nothing more, so it is
now labeled `query_constrained` instead of `target_species_likely`. It remains
preparable; only the claim changes. The README says which providers can support a
quarantine decision.

**The config workbooks are reproducible too.** `degora demo` and `degora template`
stamped the clock into the workbooks they write, so two demo runs produced inputs
with different checksums and every provenance sidecar recording an input hash
differed with them. Both now carry the same pinned timestamps the generated
workbook already used.

**`condition` is published alongside `hypoxia_modality`.** The catalog's generic
condition column reaches the SQLite schema, the API and the workbook headers under
a name from one research topic. Both names are now emitted with the same value, so
a reader can move to the neutral one before the old one is ever removed.

Smaller fixes: the preparation summary reports how many tables are ready for
review, not only how many records were prepared; a `.xlsx` that is a valid ZIP but
not a workbook is described as such instead of by an internal pandas option key;
the no-benchmark note in `slice_metrics.json` no longer names one research topic
for every corpus; and CI exercises the top of the dependency ranges as well as the
floor.

### 0.4.8

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces the
same gene scores and evidence as v0.4.7.

**A discovery analysis that scored nothing is now a failure.** `degora run`
already refused a corpus in which no gene had directional evidence from enough
independent source units. The discovery path did not: it registered the run as
complete, with an empty top-gene list and a workbook holding no ranking. It now
raises the same refusal, and the run's own transaction removes the partial
output rather than keeping it as a finished analysis.

**A prepared bundle's audit records the bundle, not the staging directory.** The
audit JSON is written while the bundle is still being staged and is then
published under its final directory. The export paths inside it were captured
before that move, so the archived document pointed at four files that had been
deleted by the time anyone opened it. It now records the published locations,
and those paths are checked for existence by a test.

**Stopping the server reaches a discovery job that is already running.**
Cancelling queued work does not touch a worker that has started, so a
preparation mid-download kept writing into the workspace after the run that
owned it had been recorded as interrupted. Every stage already reports progress,
so that is where a worker now notices the shutdown and unwinds.

**`make check` says which interpreter it is using.** On a machine whose `python3`
is older than 3.10 -- the macOS default -- the target used to select it silently
and fail deep inside the suite on modern syntax. It now names the interpreter and
stops with the command that fixes it.

### 0.4.7

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces a
byte-identical `degora_gene_scores.csv` and `degora_scores.db` to v0.4.6.

**The audit workbook is now reproducible too.** `DEGORA_output.xlsx` was the one
generated artifact whose bytes changed between two identical runs, because the
save time was written into the workbook properties and into every archive member.
Both are pinned now, so the workbook's recorded sha256 verifies the way the CSV
and database ones already did.

**The contrast table no longer loses contrasts.** When a run is scored without a
catalog, the SQLite `studies` table was built from collapsed evidence rows, which
name only the first contributing contrast per gene and source unit. A follow-up
contrast covering fewer genes than its sibling in the same source unit therefore
never appeared, and the reported contrast count disagreed with
`n_contrasts_total` in the same metadata. Both branches now report one row per
contrast with the same columns.

**Count-derived fallbacks use whole-matrix library sizes.** logCPM denominators
were summed after the low-count filter, which made each sample's scaling depend
on how much low-count mass that sample happened to lose. They are taken from the
full count matrix now, as the convention the field expects.

**A malformed workbook explains itself.** A `.xlsx` that is a valid ZIP but not a
workbook raised an engine-selection traceback from the optional settings sheets
before the catalog reader could report it. `degora validate` and `degora run` now
return the same beginner-readable configuration error they return for every other
unreadable config.

Smaller fixes: `degora discover` help and output state the real page size instead
of a stale one; a ZIP member cannot name a Windows drive-relative path; and a
second server that fails to claim a discovery workspace no longer clears the
owning server's in-process record of it.

### 0.4.6

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and the ranking, statistics, and
leave-one-source-out semantics documented above are identical to v0.4.5.
Everything below is discovery, preparation, and the browser.

**A preparation no longer fails whole.** One unreadable supplementary file - a
retired link answering with an error page, a download over the size cap, a URL
shape the checks refuse - used to end the job for every selected publication.
Failures are isolated at the candidate, the record, and the repository step, so
partial results survive. Archive safety violations - path escapes, symlinks,
and the size and count caps - still refuse the whole run, deliberately. A
download that is not the archive it claimed to be is described by what it
actually is, so a moved or login-walled file is distinguishable from a corrupt
one.

**Readiness reports what was seen.** `likely_ready` previously followed from
carrying a repository accession at all, which placed every GEO row in the top
tier of the primary sort key. It now requires a tabular file candidate;
`verified_ready` requires that plus target-species evidence; a bare accession
is `candidate`, and the basis distinguishes "never inspected" from "inspected
and carrying no tabular file". Once files have been opened, each search row
reports the prepared outcome in place of the estimate.

**Records reporting one submission are flagged.** Source units collapse on a
shared PubMed ID, so an unpublished submission deposited as several GEO series
counted as several independent studies - the count the two-source-unit rule
rests on. Records that share a title with no publication link are marked, and
selecting more than one raises a warning before the run. It warns rather than
blocks: only the reader knows whether they are separate experiments.

**Every selection is accounted for.** Publications absorbed into a repository
series another selection already covered appeared in neither the prepared nor
the excluded list, so a selection of twenty reconciled to nothing. They are now
excluded with the series and the covering publication named.

**Group assignment is decidable in place.** The labeled fallback listed bare
GSM accessions. Preparation now resolves matrix columns - including author
matrices headed by the submitter's own column names - to their sample titles
and characteristics, records the mapping in the audit bundle, and offers
filtered bulk assignment that clears the direction attestations it invalidates.
A study that produced nothing lists the files it did find and why each was
unusable.

**Browser.** Ten rows a page, so the page size no longer equals the selection
cap and a selection genuinely spans pages. The cap explains itself before a
click that cannot work. Search and preparation report real progress and lock
the controls they own. The prepared card comes into view when the work starts.
A preparation that cannot reach two source units says so above the review
instead of offering fields that the next preparation would discard.

**Quickstart.** `scripts/degora_quickstart.sh --ref NAME` runs a named branch
or tag, so a reviewer pointed at one no longer runs the default branch instead.

### 0.4.5

Version 0.4.5 preserves the v0.4.4 primary ranking and score contract. It also
ensures that stopping the local server records active fallback discovery work
as interrupted instead of allowing it to appear complete after shutdown. The
LOO nullable-field and rank-evaluable-fold semantics documented above are
unchanged.

## License

MIT
