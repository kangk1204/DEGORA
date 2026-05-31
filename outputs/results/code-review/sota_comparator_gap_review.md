# SOTA Comparator Gap Review - 2026-05-28

## Current Runnable Comparator Coverage

Already benchmarked:

- Weighted Stouffer
- Unweighted Stouffer
- Fisher
- Rank product approximation
- MetaVolcanoR-style output
- RobustRankAggreg
- Sign vote
- hStouffer and AWmeta feasibility ledgers as blocked faithful-input comparators

## Reviewer-Relevant Missing Or Partial Comparators

| comparator family | status | why reviewers may ask | recommended action |
| --- | --- | --- | --- |
| Exact RankProd Bioconductor | partially covered | Current implementation is an approximation; RankProd is a long-standing package with meta-analysis support. | Add exact RankProd R adapter if package installation is stable. |
| metaRNASeq / inverse-normal RNA-seq meta-analysis | methodologically covered, package not covered | It is a named RNA-seq meta-analysis package implementing Fisher and inverse-normal methods. | Add package-level no-go/adapter note; exact package results likely redundant with Fisher/Stouffer. |
| MetaDE / MetaOmics | missing umbrella package | MetaDE implements many p-value, effect-size, and rank methods in one microarray meta-analysis suite. | Add a parity row explaining which methods are covered and which require raw expression/effect-size variance. |
| AWFisher | missing | Adaptively weighted Fisher explicitly models study contribution/meta-patterns using p-values. | Add AWFisher as the next high-value p-value-only comparator. |
| AWmeta | blocked but important | Recent adaptively weighted transcriptomic meta-analysis framing overlaps with PIPER's quality-weighted story. | Keep blocked if variance/sample inputs are unavailable; otherwise add on datasets with faithful inputs. |
| MetaIntegrator / random-effects effect-size meta-analysis | partially covered by MetaVolcanoR REM concept, package not covered | Reviewers may ask about effect-size random-effects meta-analysis for multi-cohort signatures. | Add feasibility adapter for cohorts with raw expression or effect-size SE. |
| DExMA / ACAT-based meta-analysis | missing | Provides several meta-analysis methods including ACAT and handles missing genes. | Lower-priority, but useful as supplementary package-level comparator. |
| Gene-set enrichment methods, e.g. GSEA/fgsea/GSVA | orthogonal | Our gold standards are gene sets; reviewers may ask why not pathway enrichment. | Treat as downstream interpretation, not gene-level evidence DB competitor. |

## Current Source Snapshot

- RankProd remains the named package-level rank-product comparator to cover, because Bioconductor describes it as a non-parametric DE method that can combine datasets from different origins for meta-analysis.
- AWFisher is the highest-priority missing p-value-only comparator, because Bioconductor exposes it as an adaptively weighted Fisher implementation with p-value computing, variability index, and meta-pattern output.
- MetaDE is an umbrella microarray meta-analysis package implementing Fisher, Stouffer, adaptively weighted Fisher, min/max/rth ordered p-value, fixed/random effects, rank product/rank sum, and naive rank combinations.
- DExMA is a newer Bioconductor gene-expression meta-analysis package that explicitly handles missing genes and provides QC/GEO/visualization support.
- MetaIntegrator is a cohort-level gene-expression meta-analysis package built around DerSimonian-Laird random-effects plus Fisher sum-of-logs, so it is relevant but requires raw expression or faithful effect-size/variance inputs.
- metaRNASeq is lower priority as a named package because it mainly implements Fisher and inverse-normal p-value combination already represented by current baselines.

## Recommended Next Comparator Additions

1. AWFisher: highest value because it is p-value-only and directly tests an adaptive weighting family.
2. Exact RankProd: replaces the current approximation with a named package result.
3. MetaDE parity/no-go table: prevents reviewers from saying the broad microarray meta-analysis family was ignored.
4. MetaIntegrator feasibility: important only for datasets with raw expression/effect-size variance.

## Positioning

The strongest defense is not that PIPER-DEG beats every statistical meta-analysis method. The defense is:

- PIPER-DEG is built for heterogeneous DEG evidence databases.
- It preserves source-unit support, direction consistency, table scope, source quality, and evidence rows.
- Classical methods usually return only a ranked gene list or combined significance, and they do not provide the same browser/API-ready evidence audit surface.

## External Sources Checked

- RankProd Bioconductor package page.
- metaRNASeq package page.
- MetaDE package documentation.
- AWFisher package documentation.
- MetaVolcanoR package documentation.
- MetaIntegrator package documentation.
- hStouffer GitHub repository.
- AWmeta GitHub repository and preprint surface.
