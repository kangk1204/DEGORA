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
cd DEGORA-main       # main-branch ZIP
# or: cd DEGORA-0.4.6  # v0.4.6 release ZIP
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

From a checkout:

```bash
bash scripts/degora_quickstart.sh
```

Without a checkout, download just that one file and run it; it clones the
repository into `./DEGORA` first:

```bash
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
```

Use `--ref` to review a specific branch or release tag in one command, for
example `bash degora_quickstart.sh --ref v0.4.6`. It fetches the name from
`origin`, fast-forwards a local copy that is behind, and stops rather than
serving stale code when a local branch of the same name has diverged.

The script is safe to re-run and stops with an actionable message when the
platform is missing Python 3.10+, `git`, or the Debian/Ubuntu `python3-venv`
package. Press `Ctrl+C` to stop the server. The manual steps below do the same
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

Create a documented Excel template:

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

The search collects at most 1,000 exact, unique records before sorting and displays 10 rows per page. Detailed file resolution is bounded to the leading 20 records, the first two pages, while later selections are resolved on demand. Rows beyond that report `candidate` readiness - a repository record exists, nothing has been inspected - until they are prepared. `likely_ready` requires a tabular file candidate, and `verified_ready` requires that plus target-species evidence; a bare accession never earns either. Repository records that share a title and have no publication link are marked as a possible single submission, because source units collapse on a shared PubMed ID and an unpublished submission has none. A result is therefore a review queue, not an automatically approved analysis input.

Search exports include JSON, CSV, and Excel snapshots with identifiers, title, authors, journal, year, species evidence, source-unit information, readiness, and provider diagnostics.

You can also use the local browser:

```bash
degora serve path/to/degora_scores.db
```

Open the **Discover** tab, choose **Human** or **Mouse**, search, review the candidates, and select the records to prepare. Author-provided result tables are preferred. Matrix-derived fallback analysis requires explicit group mapping, contrast direction, and biological-replicate confirmation.

The browser opens on **Discover**. Search results are globally ranked before the first 20-row page is shown. The compact table displays publication metadata, linked-data availability, estimated DEG-input readiness, and an **Inspect** action. **Run separate Human + Mouse searches** launches two independent searches; it never pools their records or scores.

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

The included synthetic demo is numerically and semantically reproducible from this repository across the supported Python environments. Repeating a run over the same inputs reproduces `degora_gene_scores.csv`, `degora_scores.db`, and `DEGORA_output.xlsx` byte for byte: the workbook's own timestamps and its archive member timestamps are pinned rather than taken from the clock. Larger external datasets are not bundled here and must be obtained from their original providers before they can be analyzed.

## Interpretation boundaries

- DEGORA prioritizes genes; its scores are not posterior probabilities.
- Related contrasts are collapsed by source unit before cross-source aggregation.
- A row with `pvalue = 1` or zero effect is neutral evidence and does not contribute directional signed-z support. Values below `1e-300` are floored and reported in the run warnings.
- Stouffer p-values from DEG-only or significance-filtered tables are ranking aids, not calibrated genome-wide inferential p-values or false-discovery rates.
- `heterogeneity_i2` is a sample-size-weighted descriptive dispersion index, not calibrated Higgins I-squared. The heterogeneity-adjusted Stouffer fields are screening aids only.
- Random-effects effect-size intervals are descriptive when few source units contribute; in particular, Hartung-Knapp-Sidik-Jonkman intervals with fewer than three source units are unstable.
- Leave-one-source-out stability is a `priority_rank` diagnostic over global source-unit omission folds. Each fold applies the same `min_studies` eligibility rule and deterministic tie-break as the full priority lane. Median, IQR, and top-N fractions summarize rank-evaluable folds only; the stability score still treats ineligible folds as negative evidence. If no fold keeps a gene eligible, numeric LOO fields are reported as unavailable rather than as zero; the conditional reliability summary then uses the other three mandatory diagnostics with their weights renormalized. These LOO columns are therefore nullable. `evidence_reliability_components_used` records whether three or four diagnostics contributed. Compare reliability values across runs or corpora only when the contributing diagnostic count and LOO eligibility conditions match. The default browser/API ordering remains `quality_weighted_degora_rank`, and reliability does not determine that ordering.
- Missing or non-numeric group sizes receive the documented neutral replicate multiplier of 0.75; an explicit zero receives 0.35.
- `sign_convention` records provenance only. DEGORA does not infer or automatically reverse the input effect direction, so the supplied effect must already represent treatment relative to control.
- Public-data fallback uses a documented Welch workflow and does not automatically correct study-level batch or condition confounding.
- Result-table semantics, contrast direction, group labels, and species must be reviewed before activation.
- Human and Mouse evidence is never pooled into one ranking.

## Development checks

```bash
python -m pip install -e ".[dev]"
make check
make smoke
```

## Release notes

### Unreleased

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
