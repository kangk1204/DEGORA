# AGENTS.md

This is an autonomous research-engineering workspace driven by `INSTRUCTIONS.md`. Before doing anything in this directory, read in this order:

1. `INSTRUCTIONS.md` — master prompt defining the loop
2. `CONTEXT.md` — project scope, data, constraints, budget, target venue (Scientific Reports)
3. `METHODOLOGY_SPEC.md` — full technical specification for PIPER (the v2 methodology development plan)
4. `NOVELTY_AUDIT.md` — prior-art audit and contribution positioning (read this before writing any framing or Introduction text)
5. `STATE.md` — current iteration and phase (resume point if mid-loop)
6. `PLAN.md` — current iteration plan
7. `RUBRIC.md` — three publication gates tuned for Scientific Reports + integrity checklist
8. `LAB_NOTEBOOK.md` — prior iteration logs (append-only)
9. `REVIEWER_LOG.md` — prior adversarial critiques (append-only)

After reading, follow the loop in `INSTRUCTIONS.md`. Do not ask the user for approval between iterations. The three escape-hatch conditions are the only places you may pause and ask.

## Conventions for this workspace

- All code in `outputs/code/` per the layout in `METHODOLOGY_SPEC.md` §9
- Use installed skills aggressively (`bioinfo-expert`, `biostat-expert`, `ai-ml-expert`, `pipeline-robustness-tester`, `qa-tester`, `debug-expert`, `novelty-checker`, `ref-finder`, `paper-proofread`)
- Run `make check` (lint + pytest) after every non-trivial code change
- Initialize git in iteration 1 if not already; commit at meaningful checkpoints
- Long-running jobs (V2 pipeline-perturbation re-analyses, V7 classifier training) use `nohup` + log files in `outputs/logs/`
- ARCHS4 / GEO downloads cached under `data/` with existence + checksum checks before re-downloading

## File update rules

- `LAB_NOTEBOOK.md`, `REVIEWER_LOG.md`: append-only, never edit prior entries
- `STATE.md`: overwrite each iteration; last line always references the relevant lab-notebook entry
- `PLAN.md`: overwrite, but begin with `## Change log` referencing the notebook entry that justified the change
- `CONTEXT.md`, `METHODOLOGY_SPEC.md`, `NOVELTY_AUDIT.md`: user-controlled. Loop may add `## Addenda from the loop` at the bottom of `CONTEXT.md` for inherited assumptions; do not touch the other two
- Every figure/table sibling `.source` file with the producing command
