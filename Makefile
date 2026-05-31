.PHONY: install template demo validate run serve check unit typecheck lint provenance-check slice diagnose rank-plane baseline baseline-preflight tnfa-gate benchmark figs paper

DEGORA_CONFIG ?= DEGORA_template.xlsx
DEGORA_DB ?= outputs/results/degora-run/degora_scores.db
DEGORA_DEMO ?= degora-demo

install:
	python -m pip install -e outputs/code

template:
	PYTHONPATH=outputs/code python -m degora.cli template $(DEGORA_CONFIG)

demo:
	PYTHONPATH=outputs/code python -m degora.cli demo $(DEGORA_DEMO)

validate:
	PYTHONPATH=outputs/code python -m degora.cli validate $(DEGORA_CONFIG)

run:
	PYTHONPATH=outputs/code python -m degora.cli run $(DEGORA_CONFIG)

serve:
	PYTHONPATH=outputs/code python -m degora.cli serve $(DEGORA_DB)

check:
	$(MAKE) -C outputs/code check

unit:
	$(MAKE) -C outputs/code unit

typecheck:
	$(MAKE) -C outputs/code typecheck

lint:
	$(MAKE) -C outputs/code lint

provenance-check:
	$(MAKE) -C outputs/code provenance-check

slice:
	$(MAKE) -C outputs/code slice

diagnose:
	$(MAKE) -C outputs/code diagnose

rank-plane:
	$(MAKE) -C outputs/code rank-plane

baseline:
	$(MAKE) -C outputs/code baseline

baseline-preflight:
	$(MAKE) -C outputs/code baseline-preflight

tnfa-gate:
	$(MAKE) -C outputs/code tnfa-gate

benchmark:
	$(MAKE) -C outputs/code benchmark

figs:
	$(MAKE) -C outputs/code figs

paper:
	$(MAKE) -C outputs/code paper
