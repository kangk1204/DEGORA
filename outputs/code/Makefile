PYTHON ?= python3

.PHONY: check unit typecheck lint

check: unit typecheck lint

unit:
	PYTHONPATH=. $(PYTHON) -m pytest -q tests

typecheck:
	PYTHONPATH=. $(PYTHON) -m compileall -q degora tests

lint:
	PYTHONPATH=. $(PYTHON) -m compileall -q degora tests
