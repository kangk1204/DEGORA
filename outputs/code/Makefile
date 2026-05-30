PYTHON ?= python
PROJECT_ROOT := ../..
CATALOG ?= $(PROJECT_ROOT)/data/studies/hypoxia_catalog.csv
ITER ?= 1
CORPUS ?= hypoxia
OUTDIR ?= $(PROJECT_ROOT)/outputs/results/iter-$(ITER)
HARMONIZED_DIR ?= $(PROJECT_ROOT)/data/deg/harmonized
DIAG_SOURCE_ITER ?= 3
DIAG_HARMONIZED ?= $(PROJECT_ROOT)/data/deg/harmonized/iter-$(DIAG_SOURCE_ITER)_harmonized.csv
DIAG_OUTDIR ?= $(PROJECT_ROOT)/outputs/results/iter-$(ITER)
HARMONIZED ?= $(PROJECT_ROOT)/data/deg/harmonized/iter-$(ITER)_harmonized.csv
SLICE_MIN_STUDIES ?= 2
RANK_PLANE_THRESHOLD ?= 0.9
BASELINE_OUTDIR ?= $(PROJECT_ROOT)/outputs/results/baselines
SCORE_HARMONIZED ?= $(HARMONIZED)
SCORE_DB ?= $(OUTDIR)/degora_scores.db
SCORE_MIN_STUDIES ?= 2
DEGORA_HOST ?= 127.0.0.1
DEGORA_PORT ?= 8765
DEGORA_CONFIG ?= ../../DEGORA_template.xlsx
DEGORA_TEMPLATE ?= ../../DEGORA_template.xlsx
DEGORA_DEMO ?= ../../degora-demo
RESOURCE_SCORE_CSV ?= $(PROJECT_ROOT)/outputs/results/iter-18/degora_gene_scores.csv
RESOURCE_SCORE_DB ?= $(PROJECT_ROOT)/outputs/results/iter-18/degora_scores.db
RESOURCE_OUTDIR ?= $(OUTDIR)
RESOURCE_TOP_N ?= 20
RESOURCE_API_URL ?=
MANUSCRIPT_TABLES_DIR ?= $(PROJECT_ROOT)/outputs/tables/manuscript
MANUSCRIPT_FIGURES_DIR ?= $(PROJECT_ROOT)/outputs/figures/manuscript
MANUSCRIPT_PATH ?= $(PROJECT_ROOT)/outputs/manuscript/main.md
MANUSCRIPT_OUTPUT_DIR ?= $(PROJECT_ROOT)/outputs/manuscript
MICROARRAY_MATRIX ?=
MICROARRAY_DEG_OUTPUT ?= $(PROJECT_ROOT)/data/deg/raw/microarray_derived_deg.csv
MICROARRAY_SUMMARY ?= $(OUTDIR)/microarray_derived_deg_summary.json
MICROARRAY_GENE_COLUMN ?= gene_symbol
MICROARRAY_PROBE_COLUMN ?=
MICROARRAY_CONTROL_SAMPLES ?=
MICROARRAY_TREATMENT_SAMPLES ?=
MICROARRAY_PLATFORM ?=
MICROARRAY_NORMALIZATION ?=
PROVENANCE_AUDIT ?= $(PROJECT_ROOT)/outputs/results/iter-13/provenance_audit.json
PROVENANCE_ITERATIONS ?= --iteration iter-9 --iteration iter-10 --iteration iter-11 --iteration iter-12 --iteration iter-13
PUBLICATION_GATE_SUMMARY ?= $(PROJECT_ROOT)/outputs/results/iter-$(ITER)/publication_gate_summary.md
TNFA_LEDGER ?= $(PROJECT_ROOT)/data/studies/curation/iter10_tnfa_candidate_ledger.csv
TNFA_SCARCITY_REPORT ?= $(PROJECT_ROOT)/outputs/results/iter-10/tnfa_scarcity_report.md
TNFA_GATE_OUTPUT ?=
TNFA_SCARCITY_CANDIDATE_THRESHOLD ?= 20
ROOT_PYTHON := cd $(PROJECT_ROOT) && env PYTHONPATH=outputs/code $(PYTHON)

.PHONY: check slice diagnose rank-plane baseline baseline-preflight score-db resource-package microarray-deg serve degora-template degora-demo degora-validate degora-run tnfa-gate publication-gate rerun-readiness benchmark prior-art-tables gold-summaries manuscript-tables figs paper

check: unit typecheck lint provenance-check

unit:
	PYTHONPATH=. $(PYTHON) -m pytest -q tests

typecheck:
	PYTHONPATH=. $(PYTHON) -m compileall -q degora scripts tests

lint:
	PYTHONPATH=. $(PYTHON) -m compileall -q degora scripts tests

provenance-check:
	PYTHONPATH=. $(PYTHON) scripts/check_source_provenance.py --results-root $(PROJECT_ROOT)/outputs/results $(PROVENANCE_ITERATIONS) --output $(PROVENANCE_AUDIT)

slice:
	PYTHONPATH=. $(PYTHON) scripts/run_slice.py --catalog $(CATALOG) --output-dir $(OUTDIR) --harmonized-dir $(HARMONIZED_DIR) --min-studies $(SLICE_MIN_STUDIES)

diagnose:
	PYTHONPATH=. $(PYTHON) scripts/run_diagnostics.py --catalog $(CATALOG) --harmonized $(DIAG_HARMONIZED) --output-dir $(DIAG_OUTDIR)

rank-plane:
	PYTHONPATH=. $(PYTHON) scripts/run_rank_plane.py --harmonized $(HARMONIZED) --output-dir $(OUTDIR) --joint-threshold $(RANK_PLANE_THRESHOLD)

baseline:
	PYTHONPATH=. $(PYTHON) scripts/run_baselines.py --harmonized $(HARMONIZED) --output-dir $(BASELINE_OUTDIR) --corpus $(CORPUS)

baseline-preflight:
	PYTHONPATH=. $(PYTHON) scripts/run_baseline_preflight.py --output-dir $(BASELINE_OUTDIR) --corpus $(CORPUS)

score-db:
	PYTHONPATH=. $(PYTHON) scripts/build_score_db.py --harmonized $(SCORE_HARMONIZED) --catalog $(CATALOG) --output-dir $(OUTDIR) --db $(SCORE_DB) --min-studies $(SCORE_MIN_STUDIES)

resource-package:
	PYTHONPATH=. $(PYTHON) scripts/export_score_resource_package.py --score-csv $(RESOURCE_SCORE_CSV) --db $(RESOURCE_SCORE_DB) --output-dir $(RESOURCE_OUTDIR) --top-n $(RESOURCE_TOP_N) $(if $(RESOURCE_API_URL),--api-url $(RESOURCE_API_URL),)

microarray-deg:
	PYTHONPATH=. $(PYTHON) scripts/derive_microarray_deg.py --matrix $(MICROARRAY_MATRIX) --output $(MICROARRAY_DEG_OUTPUT) --summary $(MICROARRAY_SUMMARY) --gene-column $(MICROARRAY_GENE_COLUMN) $(if $(MICROARRAY_PROBE_COLUMN),--probe-column $(MICROARRAY_PROBE_COLUMN),) --control-samples $(MICROARRAY_CONTROL_SAMPLES) --treatment-samples $(MICROARRAY_TREATMENT_SAMPLES) $(if $(MICROARRAY_PLATFORM),--platform $(MICROARRAY_PLATFORM),) $(if $(MICROARRAY_NORMALIZATION),--normalization $(MICROARRAY_NORMALIZATION),)

serve:
	PYTHONPATH=. $(PYTHON) scripts/serve_degora_db.py --db $(SCORE_DB) --host $(DEGORA_HOST) --port $(DEGORA_PORT)

degora-template:
	PYTHONPATH=. $(PYTHON) -m degora.cli template $(DEGORA_TEMPLATE)

degora-demo:
	PYTHONPATH=. $(PYTHON) -m degora.cli demo $(DEGORA_DEMO)

degora-validate:
	PYTHONPATH=. $(PYTHON) -m degora.cli validate $(DEGORA_CONFIG)

degora-run:
	PYTHONPATH=. $(PYTHON) -m degora.cli run $(DEGORA_CONFIG)

tnfa-gate:
	PYTHONPATH=. $(PYTHON) scripts/check_tnfa_gate.py --ledger $(TNFA_LEDGER) --scarcity-report $(TNFA_SCARCITY_REPORT) --scarcity-candidate-threshold $(TNFA_SCARCITY_CANDIDATE_THRESHOLD) $(if $(TNFA_GATE_OUTPUT),--output $(TNFA_GATE_OUTPUT),)

publication-gate:
	PYTHONPATH=. $(PYTHON) scripts/write_publication_gate_summary.py --project-root $(PROJECT_ROOT) --iteration $(ITER) --output $(PUBLICATION_GATE_SUMMARY)

rerun-readiness:
	PYTHONPATH=. $(PYTHON) scripts/write_rerun_readiness.py --project-root $(PROJECT_ROOT) --iteration $(ITER)

benchmark: slice

prior-art-tables:
	$(ROOT_PYTHON) outputs/code/scripts/write_public_summary_tool_requirements.py
	$(ROOT_PYTHON) outputs/code/scripts/write_prior_art_coverage_summary.py

gold-summaries:
	$(ROOT_PYTHON) outputs/code/scripts/write_gold_comparator_summary.py --baseline-dir outputs/results/ifn-pilot/baselines --gold data/studies/gold/ifn_gold_panel.csv --degora-score-csv outputs/results/ifn-pilot/degora_gene_scores.csv --title "IFN locked gold-panel comparator summary" --output-csv outputs/results/ifn-pilot/ifn_gold_comparator_summary.csv --output-md outputs/results/ifn-pilot/ifn_gold_comparator_summary.md
	$(ROOT_PYTHON) outputs/code/scripts/write_gold_comparator_summary.py --baseline-dir outputs/results/er-stress-cross-platform/baselines --gold data/studies/gold/er_stress_upr_gold_panel.csv --degora-score-csv outputs/results/er-stress-cross-platform/degora_gene_scores.csv --title "ER stress locked gold-panel comparator summary" --output-csv outputs/results/er-stress-cross-platform/er_stress_cross_platform_gold_comparator_summary.csv --output-md outputs/results/er-stress-cross-platform/er_stress_cross_platform_gold_comparator_summary.md
	$(ROOT_PYTHON) outputs/code/scripts/write_gold_comparator_summary.py --baseline-dir outputs/results/heat-shock-benchmark/baselines --gold data/studies/gold/heat_shock_hsf1_gold_panel.csv --degora-score-csv outputs/results/heat-shock-benchmark/degora_gene_scores.csv --title "Heat shock locked gold-panel comparator summary" --output-csv outputs/results/heat-shock-benchmark/heat_shock_gold_comparator_summary.csv --output-md outputs/results/heat-shock-benchmark/heat_shock_gold_comparator_summary.md
	$(ROOT_PYTHON) outputs/code/scripts/write_gold_comparator_summary.py --baseline-dir outputs/results/hypoxia-hif1-benchmark/baselines --gold data/studies/gold/hypoxia_hif1_gold_panel.csv --degora-score-csv outputs/results/hypoxia-hif1-benchmark/degora_gene_scores.csv --title "Hypoxia/HIF locked gold-panel comparator summary" --output-csv outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.csv --output-md outputs/results/hypoxia-hif1-benchmark/hypoxia_hif1_gold_comparator_summary.md

manuscript-tables: prior-art-tables gold-summaries
	$(ROOT_PYTHON) outputs/code/scripts/write_manuscript_tables.py --output-dir outputs/tables/manuscript

figs: manuscript-tables
	$(ROOT_PYTHON) outputs/code/figures/make_degora_scheme_figure.py --output-dir outputs/figures/manuscript/FIGURE_1_SCHEME
	$(ROOT_PYTHON) outputs/code/figures/make_benchmark_figure.py --output-dir outputs/figures/manuscript/FIGURE_2_BENCHMARK
	$(ROOT_PYTHON) outputs/code/figures/make_cross_platform_microarray_figure.py --results-dir outputs/results/cross-platform-microarray --output-dir outputs/figures/manuscript/DEGORA_CROSS_PLATFORM
	$(ROOT_PYTHON) outputs/code/figures/make_degora_atlas_dashboard.py --output-dir outputs/figures/manuscript/DEGORA_ATLAS --prior-art outputs/results/prior_art_coverage_summary.csv
	$(ROOT_PYTHON) outputs/code/figures/make_degora_evidence_explorer.py --output-dir outputs/figures/manuscript/DEGORA_EVIDENCE_EXPLORER

paper: figs
	$(ROOT_PYTHON) outputs/code/scripts/write_manuscript_package.py --manuscript outputs/manuscript/main.md --tables-dir outputs/tables/manuscript --figures-dir outputs/figures/manuscript --output-dir outputs/manuscript
