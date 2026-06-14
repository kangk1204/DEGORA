.PHONY: install template demo validate run serve check unit typecheck lint clean

PYTHON ?= python3
DEGORA_CONFIG ?= DEGORA_template.xlsx
DEGORA_DB ?= degora-demo/results/degora_scores.db
DEGORA_DEMO ?= degora-demo

install:
	$(PYTHON) -m pip install -e outputs/code

template:
	PYTHONPATH=outputs/code $(PYTHON) -m degora.cli template $(DEGORA_CONFIG)

demo:
	PYTHONPATH=outputs/code $(PYTHON) -m degora.cli demo $(DEGORA_DEMO)

validate:
	PYTHONPATH=outputs/code $(PYTHON) -m degora.cli validate $(DEGORA_CONFIG)

run:
	PYTHONPATH=outputs/code $(PYTHON) -m degora.cli run $(DEGORA_CONFIG)

serve:
	PYTHONPATH=outputs/code $(PYTHON) -m degora.cli serve $(DEGORA_DB)

check: unit typecheck lint

unit:
	cd outputs/code && PYTHONPATH=. $(PYTHON) -m pytest -q tests

typecheck:
	cd outputs/code && PYTHONPATH=. $(PYTHON) -m compileall -q degora tests

lint:
	cd outputs/code && PYTHONPATH=. $(PYTHON) -m compileall -q degora tests

clean:
	rm -rf .pytest_cache outputs/code/.pytest_cache outputs/code/degora.egg-info build dist degora-demo DEGORA_template.xlsx
