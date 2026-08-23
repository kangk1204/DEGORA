PYTHON ?= python3

.PHONY: check unit syntax smoke build audit python-version

# macOS ships Python 3.9 as `python3`, so running these targets without
# activating the virtual environment silently selects an unsupported interpreter
# and fails deep inside the test suite on 3.10+ syntax. Name the interpreter in
# use and stop here instead.
python-version:
	@$(PYTHON) -c 'import sys; v = sys.version_info; print("DEGORA: using %s (Python %d.%d.%d)" % (sys.executable, v[0], v[1], v[2])); sys.exit(0 if v[:2] >= (3, 10) else 1)' \
	  || { echo "DEGORA requires Python 3.10 or newer. Try: make PYTHON=.venv/bin/python $(MAKECMDGOALS)"; exit 1; }

check: python-version unit syntax

unit:
	$(PYTHON) -m pytest -q

syntax:
	$(PYTHON) -m compileall -q degora tests

smoke: python-version
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
