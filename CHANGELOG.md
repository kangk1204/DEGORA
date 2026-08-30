# DEGORA release history

Historical release notes are kept here so the README remains focused on installation and first use.

## 0.4.38

The scoring and aggregation contract remains v1.3. This release changes neither
primary eligibility, score components, source weights, tiers nor primary-rank
semantics.

**Auxiliary RRA ties are deterministic across supported dependency builds.**
`rra_rank` now orders genes by log-rho rounded to 12 decimal places and then by
`gene_symbol`. The reported rho magnitudes still come from the unquantized
log-rho values; primary ranks, score components, eligibility and tiers are
unchanged.

**Manual catalogs canonicalize recognized public source-unit aliases.** The
`validate` and `run` paths now treat GEO/ArrayExpress, PMID/PMCID and DOI spelling
aliases as the same declared independent source unit. Non-numeric opaque user
identifiers and the per-contrast `study_id` fallback keep their original
spelling; bare integral source-unit values are documented as PMID identifiers.
The change closes alias-based `min_studies` inflation without rewriting
prefixed local IDs.

**Source-quality metadata names the actual primary lane.** The persisted rule
now states that the fixed reliability weight belongs to the default primary
quality-weighted ranking. High-quality incoherence and direction-conflict flags
remain advisory because automatically weighting observed agreement would turn a
review diagnostic into outcome-dependent filtering.

## 0.4.37

The scoring and aggregation contract remains v1.3:
`SCORE_VERSION` is still `degora_score_v1_3_source_unit_mean`. This release does
not change eligibility, score components, weights, aggregation, or primary-rank
semantics.

**Source-independence and direction checks are more informative without becoming
outcome filters.** DEGORA now emits bounded, advisory-only warnings when two
declared source units have highly overlapping genes and nearly identical log2 fold
changes and signed z values, or when two contrasts inside one source unit show a
strongly reversed effect pattern. The checks never change weights, scores, or ranks;
when a configured pair bound would be exceeded, DEGORA reports that the advisory was
skipped and directs the user to manual provenance review.

**Small auxiliary random-effects Stouffer probabilities remain nonzero in exported
tables.** The previous 12-decimal-place formatting converted finite tail
probabilities below `5e-13` to exactly zero. DEGORA now preserves the raw finite
floating-point values in CSV, SQLite, and workbook outputs. These screening fields
remain auxiliary and do not enter any score or rank.

**Short auto-classified result tables request explicit scope review.** When fewer
than 3,000 rows are automatically classified as `full_results`, the assessment now
states that a truncated supplementary table cannot be excluded and asks the user to
confirm `table_scope`. The effective scope and all scores remain unchanged.

## 0.4.36

The scoring and aggregation contract is unchanged at v1.3:
`SCORE_VERSION` remains `degora_score_v1_3_source_unit_mean`. This patch is
limited to public-data discovery and input validation; it does not change score
components, weights, aggregation, eligibility, or rank semantics.

**Public source-unit aliases can no longer satisfy independence twice.** GEO
`GSE`/`E-GEOD`, ArrayExpress `E-MTAB`/`EMTAB`, PMID/PMCID prefix variants, and
delimited publication-ID lists are canonicalized before record merging and
prepared-bundle validation. Opaque user-defined source-unit identifiers remain
unchanged, and genuinely distinct public identifiers remain distinct.

**Qualified measurement columns retain their true scale.** Technical suffixes
such as `unstranded`, `stranded_first`, `sense`, and `reverse` no longer hide an
earlier FPKM, TPM, count, or other supported scale label. Inspection and
activation therefore reject mixed normalized subtypes and count/normalized
families even when public exports append strandedness qualifiers.

**Compressed inputs fail closed without discarding valid peers.** Direct gzip
tables and gzip archive members must pass complete-stream and expanded-size
validation before they are marked as fully fetched. Encrypted, unsupported, or
otherwise unreadable ZIP members produce controlled candidate errors, while a
hostile checksum or size violation still refuses the archive.

## 0.4.35

The scoring and aggregation contract is unchanged at v1.3:
`SCORE_VERSION` remains `degora_score_v1_3_source_unit_mean`. This patch changes
only public-data discovery, review, and fallback-input validation; it does not
change score components, weights, aggregation, eligibility, or rank semantics.

**Public matrix activation now keeps measurement scales and biological samples
separate.** Count, FPKM, TPM, CPM, logCPM, FPKM-UQ, TMM, voom, VST, and rlog
suffixes are classified consistently. One contrast cannot mix count and
normalized families, normalized subtypes, explicit and unclassified columns,
or several measurement variants of one base sample or GEO accession. Filename
and prepared-bundle role hints cannot override an explicit incompatible suffix.

**Common public table exports reopen exactly as inspected.** Delimiter choice
prefers a validated gene/effect/p-value or gene/sample header over comma-rich
annotation text. Quoted whitespace matrices, R row names, gzipped XLS/XLSX
workbooks, and validated `EnsemblGene_GeneSymbol` composite identifiers survive
inspection and materialization with bounded readers and transform provenance.

**Independent evidence units fail closed earlier.** `degora init` reuses a
canonical GSE, ArrayExpress, project, PMID, or PMC source-unit default across
related tables, including separator variants, and requires an explicit source
unit when several unlabelled tables are present. Author-table group sizes cannot
exceed a known study sample total.

**Review status is accurate and readable.** The browser applies the server's
group-size ceiling, exposes measurement-family warnings, keeps sample traits and
result actions visible in narrow panels, and presents workflow codes as prose.
The CLI reports exact, review-required, and upstream candidates separately and
does not hide a usable matrix behind a non-activatable header candidate.

## 0.4.34

This release intentionally advances `SCORE_VERSION` to
`degora_score_v1_3_source_unit_mean`. It corrects two boundary cases in the
published formula: a zero component is now truly absorbing in the weighted
geometric mean, and support in a one-source exploratory corpus is normalized to
`1.0` in both scoring lanes. Existing workbooks pinned to the previous score
version fail closed and should be regenerated or explicitly reviewed before
changing their requested version.

**Temporal selection and rank inputs now fail closed.** `early`, `late`, and
`peak_mean` are selected from the raw source-unit rows before score eligibility
is applied, and the same selected frame feeds both audit metadata and scoring.
`peak_mean` is applied exactly once. Every scoring entry point rejects a
non-finite `normalized_rank` or one outside `0 < rank <= 1`; the sole exception
is a retained, verifiably neutral audit row (`lfc = 0` or `pvalue = 1`) whose
`signed_z` and rank are both genuinely missing and which is removed before
scoring. A negative rank can no longer become overwhelming evidence after a
logarithm.

**Score outputs are reproducible and standards-compliant.** Undefined
one-source diagnostics are JSON `null`, and every metadata/provenance writer
uses strict JSON rather than Python's non-standard bare `NaN`/`Infinity`.
Ablation names are non-blank and unique, `full` is reserved for the canonical
configuration, and score components outside `[0, 1]` are refused.

**Concurrent and cancelled work cannot publish a mixed result.** Output locks
are re-entrant only for their owning thread and retain one stable lock inode
across runs. Analysis cancellation is linearized at a pre-publication barrier,
so a cancelled job rolls its complete staged run back. Search and workbook
artifact sets publish with rollback plus a last-member generation manifest
whose SHA-256 entries detect missing, tampered, or mixed files.

**The local service is capability-protected by default.** Loopback serving now
generates a per-run token and places it only in the browser URL fragment. For
`--open-browser`, Python writes an owner-only bootstrap HTML file, passes only
that nonsecret file URL to the desktop launcher, and erases the file within 60
seconds or on shutdown. WSL paths are explicitly converted to Windows-readable
UNC/file locations. DEGORA supplies only the token-free bootstrap location as
its initial process argument; the browser follows the authenticated URL from
inside that page. Native Windows auto-open fails closed because POSIX mode bits
cannot establish an owner-only Windows ACL; WSL2 remains the supported path.
Durable discovery admission is bounded at 64 active jobs by
default and reports HTTP 429 when full. Credential-bearing query, fragment, and
userinfo parts of persisted URLs are recursively redacted while safe accession
and ID parameters remain intact.

**The source release is independently testable.** The sdist now contains both
quickstart entry points, `requirements.txt`, repository metadata, every shipped
test, and its JSON fixtures. CI runs the full suite on Ubuntu Python 3.10-3.14
and macOS, installs Node 24 explicitly, exercises the true dependency floors,
installs and tests the clean sdist, and gates requirements parity, Ruff, strict
incremental mypy, pip-audit, 85% combined branch coverage, check-manifest,
build, and twine. Release tags must match package metadata, be annotated, and
carry a cryptographic signature GitHub verifies.

## 0.4.33

The score and rank calculations are unchanged. For the same inputs, all
previously published gene-score values keep their meanings and column
positions. Audit metadata can differ where 0.4.33 corrects declared
rank-universe reporting, validation or conflict annotations, as described
below. The output schema also expands append-only with three primary-lane audit
statistics: `quality_stouffer_z`, `quality_weighted_lfc` and
`quality_rank_product`. They are added after the existing quality components in
`degora_gene_scores.csv`, the SQLite `genes` table and the workbook, so a reader
can verify the direction, effect and rank evidence used by the primary
quality-weighted lane without changing any score or rank. The release answers
two independent audits of 0.4.32: an external installation and long-run
validation, and a multi-keyword / multi-study / multi-parameter sweep of about
660 `degora run` invocations across 20 topics.

**A group size the browser refuses is refused from a workbook too.** The 10,000
biological-replicate cap lived only in the discovery review panel. From a
catalog, `n_treat=999999` — the exact typing slip 0.4.30 named — passed
`validate` and `run` in silence, and because the per-contrast weight is
`min(sqrt(n_ctrl + n_treat), 4)` it saturated that source at the ceiling instead
of the `sqrt(6) = 2.45` it had earned. On a two-source corpus of 300 genes it
moved 289 ranks, by a median of 15 places and a maximum of 134. `validate` now
applies the same bound it already applied to zero, negative, fractional and
non-numeric sizes.

**The direction-conflict flag stops naming well-formed sources.** The rule was
"at least half of a unit's pairwise comparisons conflict", and half of two is
one — so in a three-source corpus with one reversed contrast, both well-formed
units were flagged too, each with a median Spearman of 0.00 printed inside a
sentence asserting it disagreed with the rest. A unit's own median must now
clear the same -0.10 threshold. Two sources still both flag, because with two
DEGORA cannot tell which is reversed and says so. Where the flagged units are
not a minority the corpus is split rather than one source being wrong, and one
warning now says that instead of four that each contradict their own numbers.

**`lfc_scale` refuses a value it cannot honour.** Only `log2` ever meant
anything; `log10`, `ln`, `linear` and anything else were silently equivalent to
an empty cell. `table_scope` and `time_course_mode` both refuse an unknown
value, and `lfc_scale` is the field whose entire purpose is to record a
statement the reader is making about their data. It is now refused at validation
with the accepted vocabulary named.

**`--inspection-budget` is a global cap again.** Reported by the external
validation as P2: one expression,
`max(1, min(12, inspection_budget // selected_count))`, produced three contract
violations. A negative budget was refused by legacy GEO and accepted by
federated `--select`, which turned it into one file per record; `0` became one
file per record rather than no inspection; and a budget of 1 over two selections
implied a total of 2. The option is validated once at the CLI entry point, `0`
is carried through as zero inspection, and the budget is passed to the
preparation backend as the total it always claimed to be rather than being
re-derived per record. A budget smaller than the selection is refused with the
arithmetic named instead of quietly exceeding itself.

**The search API type-checks `page`.** 0.4.31 gave `query` and `species` a type
contract and 0.4.32 tested it over HTTP; `limit` had one already. `page` did
not, so `0`, `-3`, `1e9`, `1.5` and `"two"` each started a job with 202. The CLI
had always refused them. `POST /api/discovery/searches` now applies the same
bound.

**Smaller.** `peak_mean` records what it selected: `time_course_selection` was
written for `early` and `late` only, and `peak_mean` is the one mode whose
selection is per gene, so which contrasts survived differs from row to row; the
report also carries `row_retention`, because a per-gene subset can keep every
gene while dropping half the evidence. `rank_universe_size_declared` keeps the
catalog's own number when the observed row count overrides it — both columns
used to read the clamped value, so the pair that exists to be read against each
other carried one number. `reviewer_attestations` is written into
`analysis_request.json` and the run marker and printed by
`degora discovery-analyze`; it had existed only as a key of an in-process return
value, which the browser rendered and no artifact recorded. `validate` warns
when source units share no gene identifier space — the commonest way a config
that validates cleanly goes on to score nothing — while `run` keeps its own
later and more precise diagnostic. Federated `publication_search.json` records
`degora_version` and `degora_code_revision`, which the GEO snapshots already
did.

## 0.4.32

The score contract is unchanged. A sweep of every sub-report behind the
0.4.29 audits, checked line by line against 0.4.31, left two defects in the
0.4.31 fixes and a list of smaller residuals; this release closes them.

**Validate says what run says about an infinite effect.** 0.4.31 set aside a
row whose log2 fold change was `Inf`, but only at run time; `degora validate`
still printed "config OK". It now warns with the count and examples. And when
the infinite rows were the only rows set aside, the unusable-row sentence
contradicted itself ("dropped before ranking - (a row can be missing more than
one) ... the dropped cells were empty"); the infinite rows now have their one
sentence and the generic one is reserved for the rest.

**The analysis job reports every stage and can be stopped.** The job reported
two stages; the run now reports each contrast as it is derived, the scoring,
the database and the workbook, and the browser shows them. "Stop this
analysis" ends a running job at the next stage and the output folder is rolled
back, so a stopped run leaves nothing behind; the notice says so. `degora
discovery-analyze` prints the same stages and takes `--no-excel`.

**A run killed with the server is labelled.** A cooperative stop rolls back;
a server killed mid-run left a folder that looked finished. On start, the
server writes `.degora-discovery-run-interrupted.json` into any run folder
that never reached its success marker - nothing is deleted. The Excel workbook
is written beside its target and moved into place, so a kill mid-write no
longer leaves a zero-byte workbook under the final name.

**What a run cannot check is said.** `--min-studies` on `run`, `launch` and
`discovery-analyze` has help text, and a run at `min_studies=1` warns that its
ranking shows no replication. A fold-change column with fewer than 10 numeric
values is taken at its declared or named scale, and the run says so instead of
passing in silence. A matrix contrast with two replicates in a group is called
exploratory. Every discovery run summary carries `reviewer_attestations`,
naming the statements DEGORA records and cannot verify; the README says the
same. The README documents the `X-DEGORA-Action: 1` header a `curl` user
needs, and `degora discover --page N` past the last page says so.

**Smaller.** The search API's type contract is tested over HTTP for a list,
an object, a boolean, a number and null. SciPy's precision-loss warning on a
synthetic fixture is filtered in the test configuration, and a test that left
a file handle open closes it.

## 0.4.31

The score contract is unchanged. This release answers an extreme multi-keyword,
multi-study audit of 0.4.29 (94 adversarial cases, 25 analysis combinations,
a manual browser pass): two defects and three observations, all fixed.

**An infinite log2 fold change is set aside and said.** A row whose effect was
`Inf` or `-Inf` parsed as a number, passed the not-null test and reached the
evidence layer as the largest effect in its table; `validate` and `run` said
nothing. Such rows are now unusable, like a row with no effect at all, and the
run's warnings say how many there were and what they looked like, whatever
their share of the table.

**The browser's analysis is a job.** Running an analysis held one request open
for as long as the derivation and scoring took; on a large pair of matrices a
dropped connection ended in "Failed to fetch" while the server went on and
finished. The Run button now starts a job, polls it as preparation does, and
shows progress and elapsed time; a connection that drops is explained ("reload
this page: the prepared files are kept and a finished run is listed in the
Evidence atlas"), and a run interrupted by a server restart says so.

**Smaller.** The search API refuses a `query` or `species` that is not text
(a list was stringified and searched). A blank or one-character search leaves
a notice that outlives the native validation bubble. A file standing where the
output folder belongs is named as such instead of "[Errno 17] File exists". A
finished run removes its empty `.degora-run.lock`, which read as an active
lock; the lock only means something while a run holds it.

## 0.4.30

The score contract is unchanged. This release closes the last gap a live check
of the prepared-evidence card left open - a matrix the card put first was one
the run then refused - and answers an adversarial review of 0.4.29 (30 cases,
no High finding).

**Estimated counts are named, offered and accepted.** Salmon, RSEM and kallisto
write fractional counts by design. The card put a `tx2gene_counts` matrix first
as "raw counts - the least processed matrix", and the run refused it because
only 32% of its values were whole numbers - both right about the file, each
wrong about the other. The inspector now records the whole-number share of the
sample values it read; a count file whose values are fractional is introduced as
"estimated counts (fractional, as Salmon, RSEM and kallisto write them)"; the
matrix-type choice offers "Estimated counts" beside "Raw counts" and "Normalized
expression", also for a file the repository itself labels as counts; and
`matrix_type=estimated_count_matrix` takes the count path (library-size
normalised, log2, as tximport hands such counts to DESeq2) without the
whole-number test. What cannot be counts is still refused - negative values, or
a matrix that never leaves the log2 range - and the raw-count refusal now names
the estimated-counts alternative. The run's warnings say when a matrix was
taken as estimated counts, with the share of whole numbers it saw.

**One matrix file may carry several contrasts.** A multi-arm design - a shared
control against two or more treatments in one expression matrix - could
activate only one arm: the second activation of the same file was refused as
"selected more than once" whatever its groups. A second contrast with different
groups is now accepted, from the browser ("Add another contrast from this
matrix", which keeps the control arm and starts the treatment arm and the
attestation over) and from a selection file; identical or swapped groups are
still refused, a sample that changes role between contrasts is still reported,
and both contrasts fold into the one source unit the file is, so a series
still casts one vote. The derived tables no longer share a file name.

**Three smaller findings from the same review.** A declared group size is
capped at 10,000 biological replicates, so a typing slip (3 becomes 999,999) is
refused instead of silently accepted; a fallback whose gene column is one of the
sample columns is refused at once ("holds expression values, not gene
identifiers"), not diagnosed as zero genes after the run; and a query with no
Latin letter or digit (an emoji, a word in another script) is refused for its
script, not for its length.

## 0.4.29

The score contract is unchanged. This release repairs the browser's prepare
step, which a live check of the 0.4.28 card showed had been failing since
0.4.24, restores the publication resolvers on NAT64 networks, and fixes five
smaller defects the same check exposed.

**The prepared-evidence card finishes rendering again.** Preparing a selection
in the browser ended with "Preparation failed: studies is not defined". The
bundle had been built and the card drew, so a screenshot looked complete, but
the renderer stopped before the status line, the eligibility check and the Run
button were updated, and the analysis could not be started from the page. A
rename inside the renderer in 0.4.24 left three uses of the old name, and no
test executed that function. A new test now runs the whole page script in
node under a permissive fake DOM and renders a real two-study bundle, and a
second scenario covers the "only one usable study" notice that went through
the same code.

**Conditional attestation lines can hide.** 0.4.28 made five author-table
attestations conditional, but the page's flex rule for those lines outranked
its `[hidden]` rule at equal specificity, so all six lines still showed on
every card. A `.confirm-line[hidden]` rule restores the intended behaviour.

**The disclosure triangle renders.** The collapsed "Columns and table scope"
and "Advanced settings" panels showed "B8" in front of their titles: the CSS
escape `\25B8` was consumed by Python as an octal escape before it reached the
browser. The page text is now asserted to contain no control characters.

**Resolvers work on NAT64 networks.** On an IPv6-only or NAT64 network (many
phone hotspots and campus networks), a DNS64 resolver answers an IPv4-only
host with the well-known prefix `64:ff9b::` plus the IPv4 address (RFC 6052).
Python's `ipaddress` files that prefix under the reserved `::/8` block, so the
public-URL guard refused every Europe PMC, BioStudies, Crossref and DataCite
lookup as "not public" while `curl` on the same machine reached them, and each
search came back as a partial snapshot with no author tables. The guard now
judges such an answer by the IPv4 address it embeds: a public one is allowed,
a loopback or private one is still refused, and the local-use prefix
`64:ff9b:1::/48` stays blocked.

**A filter that matches nothing keeps its box.** Narrowing the results to a
text that no record matched dropped the filter input along with the rows,
under a notice blaming the data sources for an incomplete snapshot. The box
stays, the message says which text matched nothing among how many assessed
studies, the count beside it no longer reads the honest zero as a missing
value and announces "1,000 of 1,000 match", and the provider notice is
reserved for a search that returned no records.

**The columns panel stays closed until a mapping is edited.** The author
card's "Columns and table scope" panel opened by itself as soon as the reader
touched anything else on the card: the draft captures the prefilled detected
mapping on every re-render, and any captured value counted as "set". It now
opens only when a column or sheet differs from what the inspector detected,
or when the table scope is not the default.

**A plural.** The species attestation read "1 of these record was matched";
the noun is plural whatever the count.

## 0.4.28

The score contract is unchanged. This release makes the prepared-evidence card
concise for a first-time reader, without loosening what the analysis requires.

**One attestation that always matters, the rest only when they apply.** Six
attestation lines stood at equal weight under every author table, so the one a
reader must always answer - a positive value means UP in the treated group -
could not be told from the ones that apply to a table with no unadjusted
p-value, a filter, a duplicate-gene policy, an unconfirmed log2 scale, or a
column mapping the inspector was not sure of. Those five appear only when their
condition holds, and follow the reader's edits. The server asks for exactly the
same confirmations it did.

**A matrix asks one question.** "The groups above are right, the comparison is
treated minus control, and each column is a separate biological sample" - one
box for the two facts DEGORA cannot infer. Its hidden twin keeps the contract
(`direction_confirmed` and `biological_replicates_confirmed`) unchanged.

**Groups can be suggested from the sample labels.** GEO's sample
characteristics ("treatment: vehicle" / "treatment: rapamycin", "condition:
normoxia" / "hypoxia", "transfection: siControl" / "siHIF1A") usually already
split the samples in two. "Suggest groups from sample labels" proposes that
split, names the characteristic it came from and which side it took as control
and why, fills the contrast label ("rapamycin vs vehicle"), and leaves every
dropdown editable. Two groups with no recognisable control - two drugs, no
vehicle - are reported as exactly that, and nothing is assigned. A time series
or a three-arm design gets no suggestion.

**What the suggestion refuses to do.** An independent review of the change
tried to make it guess wrong, and each way it found is closed. A design with
three or more arms (vehicle / drugA / drugB, or 0h / 6h / 24h) gets no split -
"a contrast compares two; set the one arm to compare against its control and
leave the others at Ignore" - rather than two arms pooled as "the others". Two
values that both read as a control (control / mock), or neither (IgG /
anti-control antibody), are reported as exactly that; a negating prefix
(anti-, non-) is not a control word. A column GEO could not match stays at
Ignore. `DMSO` and `dmso` are one value. And when a second characteristic - a
dose, a time point, a genotype - still varies inside a proposed group, the note
says which one and that the split pools those levels. The same review found
that re-rendering the card (adding a cohort, switching species) hid a required
mapping confirmation after a column had been edited; the conditional lines
are now recomputed from the restored inputs on every render.

**Table scope is behind the fold.** Auto is right for nearly every table; the
choice sits in the collapsed "Columns DEGORA read" panel with the sheet name.
A DEG-only table without a rank universe still gets the run-time caveat.

**The contrast label is prefilled when the file says the comparison.**
`ATRA_vs_cntrl_SKNO1_gene_deseq2_out.txt.gz` becomes "ATRA vs cntrl" - the
tokens beside the "vs", to edit rather than to type.

## 0.4.27

The score contract is unchanged. This release answers a review of options and
exception handling on v0.4.26 - 32 cases, none wrong - which left one medium
finding and a few small ones, and two independent verification reports.

**The browser's completion card says what the run warned about.** A run that
mixed identifier spaces produced three warnings - a unit written in Ensembl
IDs while the corpus is in symbols, an overlap of 869 of 19,592 identifiers, a
10% minority space - and the CLI printed them, while the browser card printed
the number of source units and eight gene names. The analysis result now
carries one `warnings` list, assembled by the same function the CLI uses, and
the card renders it: identifier-space warnings, rank-universe caveats,
selection warnings and input warnings, in the same words on both surfaces.

**A hand-written catalog whose `p_column` is the adjusted column is told so.**
`degora init` refuses to pick an adjusted column as the p-value by default and
explains why; the browser requires `adjusted_p_as_pvalue_confirmed`; a catalog
with `p_column=padj` passed `validate` and `run` with nothing said. It is a
warning now, in the same words as the other two paths: every p-value read from
that table is already adjusted, which is a usable answer when the table has no
unadjusted one.

**Matrix candidates carry `gene_identifier_space` too.** v0.4.23 said every
candidate carried it; author tables did and matrices did not, so an Ensembl
matrix beside a symbol table was visible only after a run. The inspector
stamps it on matrices from their values.

**Smaller.** The header of the refusal for an invalid `time_course_mode` named
table scope; it names what the check covers. A `sheet_name` on a CSV source
was ignored in silence; it is reported as ignored, with the two things it
usually means, and it no longer lets the same CSV declared twice pass the
duplicate-table check as two files. A readiness tier carried in from upstream
that the record cannot back with a file candidate or an accession is not kept -
the only way a 1959 publication could wear a `likely_ready` badge. The hidden
preparation marker (`.degora-discovery-bundle.json`) handed to
`discovery-analyze` is named for what it is instead of failing on a species
mismatch, and the preparation summary prints the real analysis input path.

**And what an independent review of these changes turned up.** The adjusted-p
warning missed the Seurat and scanpy headers (`p_val_adj`, `pvals_adj`) it was
written for, because a clause required them not to look like a p-value; the
classifier's rule decides alone now. The result's `warnings` list carried each
source-table warning twice, once from validation and once from the run; it is
deduplicated in order. `discovery-analyze` printed only the selection warnings;
it prints the whole list. `gene_identifier_space` was stored on every candidate
and shown on none; the candidate row says it (`identifiers: Ensembl ID`). A
blank `header_row` and `header_row=1` made two identities of one table. A
tautological test assertion was replaced with one that would fail.

A review's proposed year floor for readiness was not adopted: `likely_ready`
already requires a found file candidate, and a cutoff would only hide real
evidence.

## 0.4.26

The score contract is unchanged. This release writes down the order of evidence
DEGORA prefers and makes the search reflect it.

**A "What DEGORA prefers, in order" section in the README.** The authors' own
full-results table first, then a significant-genes-only table, then a raw count
matrix, a log2-normalised matrix, and a linear one - an order of evidence
quality, with the reason for each rank, the group-size weighting, and how
identifiers and source units are counted.

**The search shows and sorts by the likely input.** Nothing is opened at search
time, so each record now carries an estimate from its file names - "likely
author DEG table", "likely raw count matrix" - shown on the readiness line and
used to break ties inside a readiness tier, so that of two records that look
equally ready the one naming a DEG table comes first. This is the same order
preparation applies once it has opened the files, and that the prepared-evidence
card uses when a series ships several.

## 0.4.25

The score contract is unchanged. This release answers one question about the
prepared-evidence card: what happens when a series ships the same samples as
several files, and which of them should a reader use.

**Selections that cannot be independent evidence are refused, with the reason.**
A series often carries raw counts, a log2 TMM matrix and an FPKM matrix of the
same samples. Selecting two or three of them with the same control and
treatment groups - same label or different labels - was accepted without a
word; they collapse into one source unit, so nothing was double-counted, but
nothing said the reader had chosen one experiment three times. Selecting two
with the groups swapped was accepted too, and the two contrasts cancelled
inside the source unit: an ordinary-looking ranking with the wrong genes at the
top. Both are refused now, naming the two selections. A sample that is a
control in one contrast and a treatment in another is reported as a warning,
because a time series or a multi-arm design does exactly that.

**A declared scale is checked against the values.** A linear FPKM matrix
declared `log2` was refused only after derivation, by the harmoniser, with a
message about the derived table. The fallback now checks the selected columns
first: values past 40 are not log2, and negative values are not linear, and
the refusal names the scale to declare instead.

**One file per series in front, the rest behind it.** Each prepared file now
carries `preference_rank` and `preference_reason`, and the file to use is
marked `preferred`. The order is by evidence, not frequency - a series ships
one of each, so there is no majority to take: the authors' own statistics
first, then raw counts (the least processed matrix, which DEGORA normalises
itself), then a log2-normalised matrix, then a linear one that has to be
transformed. The card shows the preferred file and says why; the others sit
under one collapsed line, "usually the same samples in another normalization;
open only if the file above is not the one to use". A file named
`Normalized_FPKM_gene_counts_matrix` is ranked as FPKM, not as counts.

**Selections tolerate the ways they get mistyped.** `Normalized-Expression-Matrix`,
` LOG2 `, `Fallback` and a sample name with a trailing space are read as what
they mean; a wrong value is refused with the field named. Sample column names
stay case-sensitive, because they are the file's own headers. A selection
mistake the derivation used to report as a bare `ValueError` - a sample in both
groups, a non-numeric column - is a usage error with the same words. And the
duplicate-selection guard groups selections by the series they were prepared
from, so `ctrl_1`/`treat_1` in two different studies is not one experiment.

## 0.4.24

The score contract is unchanged. This release closes the remaining findings of
the v0.4.23 audit and two review screens.

**One suffix contract for every reader.** The preparation inspector opened
`.xls` and gzipped workbooks; the Welch fallback read only plain `.xlsx` as a
workbook and the rest as CSV, so a valid `.xlsx.gz` matrix that preparation had
mapped died at derivation on a decode error; and the ZIP member filter admitted
CSV, TSV, TXT and plain `.xlsx` only, so a `DEG_results.xlsx.gz` inside an
archive was silently dropped. `WORKBOOK_SUFFIXES`, `DELIMITED_SUFFIXES` and
`TABULAR_MEMBER_RE` now live in one place and every reader and filter consults
them.

**A source unit that mixes identifier spaces is reported with the split.** The
majority rule called a column that was 70% symbols and 30% Ensembl IDs "gene
symbol", and the 30% then joined nothing in a symbol corpus with nothing said.
A run now reports a minority space at or above 10% of a unit's identifiers -
"about 30% of its gene identifiers are Ensembl ID while the rest are gene
symbol" - without converting or refusing anything. The sample that decides a
unit's space is taken evenly across its sorted identifiers rather than from the
head, where a symbol table's few `ENSG...` fallbacks sort first.

**The prepared-evidence card puts what can be activated first.** Studies with
at least one candidate are listed in their order; studies with no usable table
are grouped under one collapsed heading with their per-file reasons intact,
instead of interleaving three "no usable table" cards between the ones that
matter. The two collapsible panels on a candidate row - the columns DEGORA
read, and the advanced settings - were children of a four-column grid without
a span, so they fell into the 28-pixel checkbox column and rendered one word
per line with their inputs clipped; they span the row now, with one marker
instead of two.

**Smaller.** `degora serve` on a missing database reaches the first-run hint
(`degora demo`, then `degora run`, or `degora launch`) instead of the bare
"does not exist" that answered first. The raw-count preflight says it looked
at the first 2,000 rows and that the derivation checks every row.

## 0.4.23

The score contract is unchanged. This release answers three end-to-end reviews
of the search-to-analysis path on v0.4.21 and v0.4.22, each run on a fresh
clone with live searches.

**A selection cannot fan out into a download of the whole repository.** One
publication can link dozens of GEO series - a consortium publication linked 51 - and a
selection of 20 publications quietly became 69 series to download. A record
linking more than five series is set aside at preparation with the count and a
note to search for the series it wants and select them directly.

**A gene column is judged by its values, everywhere.** A column named "Gene
Symbol" that held Ensembl IDs was accepted on its name, overlapped a symbol
table on the 4% of rows where both sides happened to write Ensembl, and the run
reported success over 869 genes that were all ENSG. One rule, shared by the
guided setup, the preparation inspector and the run, now names the space a
column is written in from its values. The inspector prefers the candidate whose
values join best (symbols, then Ensembl, RefSeq, Entrez) and records
`gene_identifier_space` on every candidate, so the mismatch is visible at
selection. A run warns when a source unit is written in a different space from
the rest of the corpus even when some identifiers overlap; before, that warning
required an overlap of exactly zero.

**A fractional matrix is refused as raw counts.** `matrix_type=count_matrix`
over FPKM, TPM or a log scale would hand a count model something that is not
counts; when fewer than 95% of the selected columns' values are whole numbers,
the selection is refused with the alternative named.

**A gzipped workbook opens on the author path.** Preparation inspected and
mapped an `.xls.gz` author table and said it could be activated; activation
read it as CSV and died on a UTF-8 decode error. The activation reader now
takes the same route as every other reader.

**Smaller.** `degora discover` refuses a one-character query and a `--page`
below 1 instead of searching or silently clamping. `degora serve` on a path
with no database says where one comes from (`degora demo`, then `degora run`,
or `degora launch`) instead of only that the path has none.

## 0.4.22

The score contract is unchanged. Two defects in the path from a search to an
analysis, found by running it end to end on one series that carries author DEG
tables and one that carries expression matrices only.

**An R export's gene column is recognised at preparation.** R writes the gene
names as row labels under an empty header cell. `read_deg_table` and the guided
setup recover that column as `row_name`; the preparation inspector classified
the empty name instead, saw no gene column, and rejected six DESeq2 result
files as `not_deg_table` - files the guided setup had scored the same day. The
inspector now names the column the way the reader does, and the author-table
step restores it too, so the mapping confirmed in the browser is found in the
file it came from. That series went from 0 tables ready to 6.

**A wrong `matrix_type` is a usage error.** An unrecognised value in a fallback
selection ended in a bare traceback from the derivation; the message now names
the two accepted values, `count_matrix` and `normalized_expression_matrix`.

## 0.4.21

The score contract is unchanged. Two screens from a real search, read as a
first-time user would.

**A full selection can be prepared.** The browser caps a selection at 20
publications; the repository phase re-applied the same cap to the GEO series
those publications link to. A selection in which one publication carried three
series expanded to 22, and every repository record failed with "at most 20
studies can be prepared at once" - one usable study and a red banner. The cap
counts what the reader selected; the repository phase takes the expanded list.

**The results table says less, more clearly.** The readiness badge carried a
claim and its caveat in one pill - "data confirmed · nothing inspected yet" -
which read as a contradiction. The badge holds one phrase and the caveat sits
under it, beside the relevance figure. The species note read "species target
species verified"; it says "species verified". The Inspect column showed
"Inspect …", because every cell clips with an ellipsis and this one holds only
a button; it no longer clips. Sort indicators are arrows, not the letters
`^` and `v`.

## 0.4.20

The score contract is unchanged, and a run over unchanged valid inputs produces
the same `degora_gene_scores.csv` as v0.4.19. This release answers a review of
the first-time journey - keyword, selection, preparation, `init`, `validate`,
`run` - on a fresh clone with live searches. What failed loudly was already
clear; these are the places that failed quietly.

**A search with no records says what to try.** "0 globally ranked record(s)"
was a fact with no next step. A query that matched nothing now suggests a
broader English term, a spelling check and the other species; a query where a
source did not answer says so instead, because retrying is the right move
there and rewording is not.

**"0 table(s)" says why, per record.** A record that looked ready in the search
and prepared to nothing had its reason in `discovery_audit.json` and nowhere
else. The preparation summary now prints one line per record with nothing
ready: expression matrices only (choose sample groups in the browser), a
results table with adjusted p-values only (confirm padj in the browser), no
usable table, or no files. The judgement was already right; the reader could
not see it.

**A prompt with one option takes Enter.** The gene prompt offered `[SYMBOL]`
and Enter took it; a padj-only table's p prompt offered `padj` alone with no
default, so Enter counted as a failure and three of them dropped the table -
the most common keystroke on the most common table. One option is now the
default.

**The messages for an existing output folder, and for `--force`, name the
remedy.** Both were true and neither said that each search writes its own
folder, that `--force` replaces a preparation rather than re-running a search,
or that the browser's Discover tab pages through one snapshot without
searching again.

**A guarded value in a workbook is treated like one in a CSV.** The formula
guard restored or refused an apostrophe-guarded value read from a CSV, but the
workbook path did neither, so `'=ISG15` in an `.xlsx` came out as `''=ISG15`
and gained an apostrophe on every round trip with nothing said. Both formats
now follow one rule, scoped to the columns a run uses.

The v0.4.18 note's statement that scores are unchanged is now qualified: an
input carrying a dotted gene symbol changes, in the right direction, because
`CAND1.11` is no longer merged into `CAND1`.

## 0.4.19

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged valid inputs
produces the same `degora_gene_scores.csv` as v0.4.18. What changes is what
`degora init` offers, what `validate` refuses, and how long both take.

Everything here came from pointing the guided setup at three topics' worth of
tables downloaded straight from public repositories, uncurated: 40 files, from
a 200 KB DESeq2 export to a 203 MB interaction table that happened to be in the
same series. That kind of input finds what curated data cannot.

**Columns that can never be an effect size are no longer offered as one.** A
results file of counts, FPKM and an FDR - a real supplementary table with no
fold change in it - was classified as a DEG table and its count column offered
as the effect size, and taken. An Entrez ID column and a base-pair coordinate
column were offered on two other tables, for the same reason: they are numeric.
`degora init` now excludes abundances (count, FPKM, TPM, CPM, expression -
`ReadCount` and `baseMean` as much as `CP1_count`), identifier columns and
columns of large integers from the effect-size candidates, and a table left with no candidate is skipped
with that reason - "no effect-size column" - rather than walked through. The
fallback that kept the p-value columns on offer "so the reader can still pick"
is gone: there was nothing correct to pick.

**`validate` refuses a column that is mostly impossible log2 values.** The
`|log2FC| > 30` check was a warning, so a count column, an Entrez ID column and
a coordinate column each reached scoring with a note and exit 0. When more than
5% of a column's values (and at least twenty of them) are past that line, the
column is refused with the plain statement that it is not a log2 fold change.
One outlier in a small table is still a warning.

**A stated scale is honoured.** `degora init` asks whether an effect column
whose header does not say log2 is on a log2 scale, and the answer went nowhere:
`validate` refused an up-only log2 table named `lfc` as a linear fold change
regardless. The optional `lfc_scale` column carries the answer (`log2`), the
refusal names it as the way out, `init` writes it, and the template's
ColumnGuide documents it.

**`log2(fc)` is recognised.** A table carried `fc` and `log2(fc)` side by side
and the header classifier matched neither, so the reader was offered the linear
one. `log2` followed by `fc` with any separator - `log2(fc)`, `log2 (FC)`,
`log2.fc`, `log2_fc` - is the effect column now.

**The guided setup no longer stalls for minutes on a large file.** Inference
read every table in full; a 203 MB file took minutes with nothing on screen.
Reads are capped at 250,000 rows (any DEG table is far smaller), and a capped
table reports its count as a floor and leaves the scope decision to the run,
which reads everything. A workbook is opened once: openpyxl parses the whole
shared-string table on each open, about three seconds for a 12 MB file, and
locating a titled table on a later sheet opened it twenty-four times. Sheets
are read into memory and candidate header rows are tried by slicing, and only
rows that could head a table are tried at all. The 12 MB workbook went from 23
seconds to 5; the 203 MB file from minutes to under a second. A file over 5 MB
is named before it is read, so the pause has a cause.

**Two common file shapes are read.** A `.csv` written by R's `write.csv2`, or
by Excel where the decimal mark is a comma, is semicolon-delimited; read with a
comma, a gene description containing one made the row ragged and the file
failed with "Error tokenizing data". The header line now decides the delimiter
when the catalog does not. A UTF-8 byte-order mark no longer becomes part of
the first header, which had turned `GENEID` into a column no rule could match.

**A symbol column is preferred over a descriptive name.** `GENEID`, `GENENAME`
and `SYMBOL` side by side offered `GENENAME` - "SPARC like 1" - because its
name matched first. Among several gene-column candidates the default is now the
one whose values join best: symbols, then Ensembl, RefSeq and Entrez IDs. Every
candidate stays on offer.

**An author's table is read as written.** The formula guard added in v0.4.17
refused a raw table because an unmapped notes column held `'=see figure 2`,
and told the reader to restore a provenance sidecar DEGORA had never written.
Only the columns a run will use can make a file ambiguous, and a file with no
sidecar at all is refused only when one of those columns carries guard-like
text - with a message that says what to do about it.

**A search says what it did not search for.** DEGORA searches PubMed and GEO,
which index English text. NCBI drops a term written in another script without
saying so, and a query with nothing left is the organism filter alone - every
Human record in GEO, newest first, presented as a result. A term with no Latin
letter or digit is set aside before any request is sent and named in the output
(`Ignored (not English): ...`); a query with no English term at all is refused
with that explanation. Greek letters inside an English term (`TGF-β`,
`α-synuclein`) are handled by NCBI and are kept. And when a source did not
answer, the CLI now says the snapshot is partial and names the source, instead
of printing a complete-looking list.

**`degora discover` reports its stages.** Two live searches were silent for 98
seconds and then printed everything at once; the provider layer had reported
every stage through a callback the browser displays and the CLI never passed.
It does now: the first line appears in under a second.

**`make smoke` exercises the alias columns.** A blank `paper_id` beside a filled
`source_unit_id`, and a legacy `temporal_mode`, now run through validate and
run in CI under the top of the pandas range - the shapes that had crashed
v0.4.16 without any job noticing.

## 0.4.18

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged valid inputs
produces the same `degora_gene_scores.csv` as v0.4.17 - with one qualification:
an input carrying a dotted gene symbol changes, in the right direction. A
source table that listed `CAND1.11` beside `CAND1` had the two merged into one
gene for seven releases, because `.11` was stripped as a version suffix; kept
apart, `CAND1` moves from up to down and every other rank shifts by at most
one. What changes is which inputs are accepted: a config that passed v0.4.17 can now be refused at
validation when it carries a linear fold-change column, a p-value written as a
bound, a text `duration_h` under `early`/`late`, a mapping onto a repeated
header, or one result table declared as two independent source units. Each of
those was producing a wrong ranking in silence, so the refusal is the fix.

Fixes from an independent code audit of the v0.4.17 branch (three module
audits, edge-case runs, and a live comparison of the GEO organism parser against
the records GEO actually serves), on top of the first-run review below.

Input boundaries that silently changed results:

- A linear fold-change column is refused. Direction is the sign of the effect,
  so a linear ratio (2.5 = up, 0.4 = down) made every gene up with no warning; a
  200-gene table with 97 down-regulated genes called none of them down.
  `validate` and `run` now refuse a column with no negative values and values
  on both sides of 1 unless its name says log2, warn when the name does not say
  log2 and the values could be a signed linear ratio or an up-only list, and
  warn on |log2FC| > 30. `degora init` asks whether the values are log2 when the
  column name does not say so, shows their range, and writes a "no" as an
  excluded row.
- A p-value written as a bound (`<1E-16`, `<0.001`, `p<0.05`) is refused at
  validation instead of being dropped as unreadable. Those rows are the most
  significant genes in a table, and below a 10% share nothing said they were
  gone. Every run now also reports how many rows each source lost to a missing
  or unreadable gene, effect or p-value, whatever the share.
- `early`/`late` time-course units must carry a plain numeric `duration_h` on
  every row. `30min` used to parse as 30 and `4h` as 4, so `early` kept the
  4-hour contrast of a unit whose earliest point was 30 minutes and inverted
  every gene that changed between them; a blank duration was silently dropped.
- Dotted gene symbols are kept. The version-stripping rule that turns
  `ENSG00000141510.16` into `ENSG00000141510` also turned `NKX2.5` and `NKX2.1`
  into one symbol, `NKX2`, that names no gene and matches no partner table. The
  suffix is now removed only from accession-shaped identifiers (Ensembl,
  RefSeq, an Entrez ID exported as `7157.0`). The same rule is one function for
  source tables, the GoldPanel, API lookups and the reanalysis matrices, and it
  treats the literal text `NA`, `<NA>`, `N/A` and `NULL` as missing everywhere.
- A header the source table carries twice cannot be mapped by its bare name.
  pandas renames the second `logFC` to `logFC.1` before anyone sees it, so a
  two-contrast supplementary table bound to the first block in silence; the
  mapping is refused with the pandas name of each copy so a later block can be
  chosen explicitly.
- One result table declared under two different source units is refused when
  the column mappings are identical (the same file, or byte-identical files),
  and warned about when they differ. The same table under U1 and U2 used to
  validate, run, and give every gene "2 / 2 source units" with perfect
  concordance. CI's own `.xls` check had relied on exactly that: its catalog
  declared one workbook under two units, so it now writes two workbooks and
  asserts that every scored gene carries evidence from both.
- A GoldPanel `locked` column typed as 1/0 with a blank cell beside it is read
  correctly. pandas delivers `1.0`, which matched nothing in the flag set, so
  the rows the reader had marked were dropped and the blank rows kept.
  `include_in_analysis` already handled this; the two panel readers now share
  one flag parser.
- A source unit whose log2 fold changes run against the other source units
  (pairwise Spearman below -0.10 and significantly negative for at least half
  of its comparisons) is flagged
  `source_direction_conflict_flag` and named in the run warnings as a possibly
  reversed contrast. The coherence guardrail only ever looked at low-quality
  sources with a near-zero correlation, so a well-documented author table with
  its sign inverted kept full weight and was never mentioned. The flag changes
  no weight and no rank.
- The Welch fallback floors within-group variances at the 1st percentile of
  the matrix's positive within-group variances. Identical replicates gave
  t = infinity and p = 0.0, which the harmonizer floored to p = 1e-300: a
  1.07-fold change with matching duplicates was ranked as the most significant
  gene in a corpus, above a 32-fold change. Genes above the floor get exactly
  scipy's Welch result; the floor and the number of genes it touched are
  recorded in the derived table's provenance.
- A `paper_id`-only catalog (documented as accepted) whose papers use different
  `time_course_mode` values passed `validate` and crashed `run` with a traceback
  naming an empty source unit; with the same mode, the `time_course_selection`
  audit in `slice_metrics.json` reported one unit instead of several. The report
  now resolves source units the way the scorer does.
- The row-loss warning distinguishes empty cells (DESeq2 leaves `pvalue` blank
  for untested genes; no action needed) from text that would not parse.

`degora serve` and the browser:

- The database preflight added for v0.4.17 built its SQLite URI from the raw
  path, so a `#`, `?` or `%XX` in the path made SQLite open a different file,
  read-write, create it when it did not exist, and report the reader's good
  database as not a DEGORA database. The preflight now opens the file exactly
  as every request does, and it checks the `studies` table and the columns the
  dashboard's first requests read, so a database that passes cannot fail
  `/api/health` with a 500 a moment later.
- Ticking "I have confirmed they are Human data" did not enable **Run**, and
  unticking it did not disable it: the listener that recomputes eligibility was
  bound to the candidate list, and the checkbox lives in the card footer. With
  the GEO organism fix below, that box now appears only for records whose
  species could not be checked.
- GEO organism evidence was never read. The parser looked for
  `!Sample_organism_ch1`, which GEO emits only in per-sample records, while
  DEGORA fetches the Series record, which lists `!Series_sample_organism`. Every
  GEO record was therefore `query_constrained`, mixed-species series were not
  quarantined at search time, and `verified_ready` could not be reached by any
  provider. The Series keys are read now, a record whose samples name exactly
  the requested species carries `target_species_verified`, and the test
  fixtures use the keys GEO actually serves.
- Stopping a re-run search left the previous search's rows on screen under the
  new search id, where they could not be prepared. A new search empties the
  previous snapshot's rows when it is submitted.
- `degora demo --species mouse` opens the browser on the Mouse workspace with
  its pre-filled keyword, instead of on Human with an empty box.
- A permissions problem creating the `degora_discovery/` workspace is reported
  as one, not as "port is already in use" after twenty futile bind attempts.
- On Linux and macOS the server can be restarted on the same port straight
  after Ctrl-C; TIME_WAIT connections from the previous run no longer push it
  to the next port with a message about another DEGORA that does not exist.
- `--host ::1` binds as IPv6 and prints a bracketed URL; a loopback alias such
  as `127.0.0.2` is treated as loopback by every check, not only by the
  Host-header check.
- A closed tab or an aborted download no longer prints a nested traceback: the
  handler returns instead of writing a second response to a dead socket, and a
  stalled upload is no longer diagnosed as a filesystem error.
- `Infinity` in a JSON `limit`, a list in `species_confirmation_required_for`,
  and a non-ASCII access token each answer with a 4xx instead of dropping the
  connection.
- Network-mode redaction uses the same path rule as the discovery store, so
  `1%/21% O2` and `ratio (A)/(B)` are no longer replaced by
  `[redacted: local path]`.
- "Narrow these results" filters as you type; the Stop button on the progress
  card keeps its identity across polls so a slow press is not lost; search
  timestamps are ISO-8601 like job timestamps; tuples with non-finite values
  serialise as `null` and a non-finite value can no longer reach the wire as
  the bare token `NaN`.
- A provider error carrying `Authorization: Bearer ...`, `Basic ...`, an
  `X-API-Key` header or a `user:password@host` URL is redacted before it is
  printed or stored.

Preparation, packaging, and the command line:

- One field past csv's 128 KiB limit in a supplementary file no longer ends a
  whole preparation with a traceback; the matrix reader has the guard the
  DEG-table reader already had.
- Every generated file has the permissions a plain write would give it under
  the current umask. Files written through a temporary file (the database, the
  workbook, the sidecars) came out owner-only while the CSV beside them was
  group-readable.
- Two runs with output folders of the same name no longer overwrite each
  other's harmonized copy under `harmonized/`; the second copy gets a suffix
  derived from its output path.
- `degora run` validates before it claims the output directory, so a config
  that fails its preflight leaves no empty results folder or lock file behind;
  a run interrupted by a closed pipe (`| head`) exits 141 instead of 0.
- `degora template` refuses a name that does not end in `.xlsx` before writing
  anything; `degora template <dir>` and `degora demo <file>` say what is wrong
  instead of printing a traceback; a demo keyword starting with `=`, `+`, `-`
  or `@` is stored as text rather than as a live formula.
- `degora init` locates a table on a later sheet or below a title row and
  records `sheet_name` and `header_row`; recognises Seurat `p_val_adj`, scanpy
  `pvals_adj` and R `p.adjust` columns as adjusted p-values and no longer offers
  `pct.1`/`pct.2` as p-values; skips DEGORA's own output tables; and accepts
  only a whole number (or a blank) for a group size and no `;` in a source unit
  id, so what it writes is what `validate` accepts.
- `validate` names a UTF-16 source table and a semicolon-delimited config for
  what they are.
- `scripts/degora_quickstart.sh` needs `git` only to clone, `--ref` or
  `--update`; an unpacked ZIP folder works without it. `--config` runs from the
  config's own folder so its results land beside it, `--demo-dir` accepts an
  absolute path, and the temporary run log is removed on failure. A
  `.gitattributes` file keeps the shell scripts LF-only on Windows checkouts run
  under WSL.
- `--min-studies 0` is reported as a command-line value, not as a spreadsheet
  cell.
- `degora ablate` exposes the component-ablation and weight-sensitivity
  machinery (`score_db.ScoreAblation`) that was reachable only from Python, and
  the RRA and effect-size meta-analysis lanes are computed without a per-gene
  Python loop, which makes a 15,000-gene scoring pass about four times faster
  with results identical to the last digit.

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

## 0.4.17

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

## 0.4.16

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

## 0.4.15

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

## 0.4.14

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

## 0.4.13

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

## 0.4.12

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

## 0.4.11

The score contract is unchanged. `SCORE_VERSION` remains
`degora_score_v1_2_source_unit_mean`, and a run over unchanged inputs produces the
same gene scores and evidence as v0.4.10.

**Species provenance now reaches the prepared bundle and its audit.** v0.4.10 made
the `query_constrained` label correct in the search snapshot, but preparing a
selected publication dropped `species_decision`, `species_evidence` and
`target_species_verified` on the way into the bundle. The species gate reads those
fields off the search record, so preparation kept working while the archived
`discovery_audit.json` -- the document an analyst opens -- reported none of them,
and the prepare view renders its species provenance line only when the decision
survives. All three are now carried through, and the regression test asserts the
whole chain rather than the helper that copies the fields.

**The unreadable-workbook message matches the extension.** Only `.xlsx` is a ZIP
container; a legacy `.xls` is an OLE2 compound file, so the message named the
wrong container for one of the two extensions it covered.

## 0.4.10

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

## 0.4.9

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

## 0.4.8

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

## 0.4.7

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

## 0.4.6

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
or tag, so a user pointed at one no longer runs the default branch instead.

## 0.4.5

Version 0.4.5 preserves the v0.4.4 primary ranking and score contract. It also
ensures that stopping the local server records active fallback discovery work
as interrupted instead of allowing it to appear complete after shutdown. The
LOO nullable-field and rank-evaluable-fold semantics documented above are
unchanged.
