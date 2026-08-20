PYTHON ?= python3

.PHONY: check unit syntax smoke build audit

check: unit syntax

unit:
	$(PYTHON) -m pytest -q

syntax:
	$(PYTHON) -m compileall -q degora tests

smoke:
	@smoke_root=$$(mktemp -d); \
	trap 'rm -rf -- "$$smoke_root"' EXIT; \
	$(PYTHON) -m degora.cli --version; \
	$(PYTHON) -m degora.cli demo "$$smoke_root/demo"; \
	$(PYTHON) -m degora.cli validate "$$smoke_root/demo/degora_demo_config.xlsx"; \
	$(PYTHON) -m degora.cli run "$$smoke_root/demo/degora_demo_config.xlsx"

build:
	$(PYTHON) -m build

audit:
	$(PYTHON) -m pip_audit
