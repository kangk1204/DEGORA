# Coding Lord Expert Sweep

Task: whole-code ultra review and fixes for PIPER-DEG.

Recorded before implementation. The implementation pass accepts only findings
that are reproducible from local code and testable without changing the
scientific claim scope.

## Expert Roster

1. `env-mixture-expert` - `/home/keunsoo/.codex/skills/env-mixture-expert.backup-20260426112846/SKILL.md`
2. `biostatinfo-expert` - `/tmp/kangs-codex-skills/skills/biostatinfo-expert/SKILL.md`
3. `citing-expert` - `/tmp/kangs-codex-skills/skills/citing-expert/SKILL.md`
4. `coding-expert` - `/tmp/kangs-codex-skills/skills/coding-expert/SKILL.md`
5. `coding-master` - `/tmp/kangs-codex-skills/skills/coding-master/SKILL.md`
6. `debugging-expert` - `/tmp/kangs-codex-skills/skills/debugging-expert/SKILL.md`
7. `demis-expert-agent` - `/tmp/kangs-codex-skills/skills/demis-expert-agent/SKILL.md`
8. `dennis-expert-agent` - `/tmp/kangs-codex-skills/skills/dennis-expert-agent/SKILL.md`
9. `donald-expert-agent` - `/tmp/kangs-codex-skills/skills/donald-expert-agent/SKILL.md`
10. `env-mixture-expert` - `/tmp/kangs-codex-skills/skills/env-mixture-expert/SKILL.md`
11. `gennady-expert-agent` - `/tmp/kangs-codex-skills/skills/gennady-expert-agent/SKILL.md`
12. `geoffrey-expert-agent` - `/tmp/kangs-codex-skills/skills/geoffrey-expert-agent/SKILL.md`
13. `grace-expert-agent` - `/tmp/kangs-codex-skills/skills/grace-expert-agent/SKILL.md`
14. `innovator-expert` - `/tmp/kangs-codex-skills/skills/innovator-expert/SKILL.md`
15. `jeff-expert-agent` - `/tmp/kangs-codex-skills/skills/jeff-expert-agent/SKILL.md`
16. `ken-expert-agent` - `/tmp/kangs-codex-skills/skills/ken-expert-agent/SKILL.md`
17. `linus-expert-agent` - `/tmp/kangs-codex-skills/skills/linus-expert-agent/SKILL.md`
18. `novelty-expert` - `/tmp/kangs-codex-skills/skills/novelty-expert/SKILL.md`
19. `petr-expert-agent` - `/tmp/kangs-codex-skills/skills/petr-expert-agent/SKILL.md`
20. `tavis-expert-agent` - `/tmp/kangs-codex-skills/skills/tavis-expert-agent/SKILL.md`
21. `yann-expert-agent` - `/tmp/kangs-codex-skills/skills/yann-expert-agent/SKILL.md`
22. `yoshua-expert-agent` - `/tmp/kangs-codex-skills/skills/yoshua-expert-agent/SKILL.md`

## Pass Ledger

| Skill | Lens applied | Finding | Recommendation | Verification | Decision |
| --- | --- | --- | --- | --- | --- |
| env-mixture-expert (backup) | Mixture-source bias | No distinct code finding beyond source-quality/evidence consistency. | Use current skill copy for any substantive pass. | No new test needed. | Already covered |
| biostatinfo-expert | Unit of analysis, estimand, multiple testing | `study_gene_evidence()` can report metadata from unselected time-course rows although consensus has already selected early/late/peak rows. | Expose the exact source-unit selected rows and build evidence metadata from those rows. | Add regression where late lower p-values are excluded under `time_course_mode=early`. | Accept |
| citing-expert | Citation/provenance integrity | No citation text was edited; external literature search is not required for this code-only bug fix. | Keep code provenance auditable through sidecars and tests. | Report sidecar and pytest. | No applicable finding |
| coding-expert | Implementation contract | Aggregation selection logic is duplicated implicitly between score and evidence layers. | Factor selected source-unit rows into a reusable helper. | Targeted score DB tests. | Accept |
| coding-master | Cross-skill integration | Accepted findings converge on one invariant: evidence must describe the rows used for scoring. | Implement one small helper and reuse it instead of broader refactor. | Narrow tests then `make check`. | Accept |
| debugging-expert | Failure-path review | Reproducer: early mode with T1/T2/T3 keeps T1 in score but `min_source_pvalue` and contributing IDs include T2/T3. | Build evidence metadata after time-course selection, not before. | Regression asserting selected-only metadata. | Accept |
| demis-expert-agent | Research strategy | Auditability is a paper-critical claim for PIPER-DEG; inconsistent evidence rows weaken the tool's rationale. | Prioritize invariant repair over new features. | Full check suite. | Accept |
| dennis-expert-agent | Interface simplicity | Selection policy is private to `aggregate.py`, forcing other modules to guess. | Add one narrow helper with no new dependency. | Compileall and tests. | Accept |
| donald-expert-agent | Invariant/proof | Invariant should hold: for each evidence row, contributing metadata is a function of selected rows only. | Test early-mode and mixed-replicate edge cases. | Regression tests. | Accept |
| env-mixture-expert (current) | Heterogeneous-source weighting | Source-quality replicate multiplier can use first replicate count from a source-unit, making the weight row-order dependent. | Aggregate replicate counts conservatively with numeric minima. | Regression with mixed replicate counts. | Accept |
| gennady-expert-agent | Edge cases | Whitespace-only config cells can pass required-value validation and fail later as confusing column/file errors. | Strip ordinary whitespace in `_nonempty()` while preserving literal tab separator. | Config validation regression. | Accept |
| geoffrey-expert-agent | Model evaluation | No ML training or representation-learning code path is in scope. | No code change. | Existing benchmark tests remain enough. | No applicable finding |
| grace-expert-agent | Beginner diagnostics | Required config value with spaces should be reported at config-validation time. | Tighten `_nonempty()` and keep actionable `PiperConfigError`. | Config validation regression. | Accept |
| innovator-expert | Transferability/DB use | For NutriOmics-style DB use, browser evidence must reflect selected contrasts exactly. | Repair evidence consistency before adding new benchmark topics. | Score DB regression. | Accept |
| jeff-expert-agent | Data pipeline reliability | Shared selection logic reduces pipeline drift between consensus and SQLite evidence outputs. | Centralize selected-row preparation. | Full `make check`. | Accept |
| ken-expert-agent | Small primitives | Avoid a new scoring abstraction; factor only row preparation. | Add a small helper and reuse local functions. | Compileall. | Accept |
| linus-expert-agent | Regression risk | Large refactor would risk existing benchmark outputs; keep patch narrow. | Edit only aggregate, score DB, validation, and tests. | Full suite. | Accept |
| novelty-expert | Prior-art/claim novelty | No new novelty claim is introduced by this bug fix. | Defer live novelty search to manuscript-framing work. | No new test needed. | No applicable finding |
| petr-expert-agent | Boundary cases | Empty selected frames and missing optional metadata must preserve stable columns. | Helper should return empty selected-row frame with source-unit columns. | Existing empty-frame tests plus targeted tests. | Accept |
| tavis-expert-agent | Trust boundary/security | API uses parameterized SQL and local HTML escaping in inspected paths; config parsing still benefits from stricter blank handling. | Accept whitespace validation; no API security patch found. | Config regression; existing API tests. | Accept |
| yann-expert-agent | Vision-specific evaluation | No image/video model path is in scope. | No code change. | No new test needed. | No applicable finding |
| yoshua-expert-agent | Language-model/NLP lens | No LLM/NLP model path is in scope. | No code change. | No new test needed. | No applicable finding |

## Synthesis

Accepted implementation contract:

- `collapse_gene_source_units()` and `study_gene_evidence()` must use the same
  selected source-unit rows after time-course policy application.
- Evidence metadata, source-quality weights, min source p-values, contributing
  IDs, durations, and URLs must be derived from those selected rows.
- Replicate metadata for a source unit should be conservative when rows disagree.
- Beginner config validation should treat whitespace-only required cells as blank.

## Implemented Fixes

- Added `source_unit_rows_for_aggregation()` in `outputs/code/piper/aggregate.py`
  so consensus and evidence can share the same selected-row invariant.
- Changed `study_gene_evidence()` in `outputs/code/piper/score_db.py` to build
  metadata from the selected source-unit rows after time-course policy is
  applied.
- Changed source-unit replicate metadata to use conservative numeric minima
  instead of row-order-dependent first values.
- Tightened `outputs/code/piper/slice_runner.py` `_nonempty()` so
  whitespace-only config cells are treated as blank while preserving tab
  separator input.
- Added regressions in `outputs/code/tests/test_score_db.py` and
  `outputs/code/tests/test_config_validation.py`.

## Verification

- `PYTHONPATH=outputs/code python3 -m pytest -q outputs/code/tests/test_score_db.py outputs/code/tests/test_config_validation.py outputs/code/tests/test_aggregate_metrics.py`
  -> 26 passed.
- `make -C outputs/code check` -> 96 passed, compileall passed, default
  provenance audit passed; one existing scipy precision warning.
- `git diff --check` -> passed.
