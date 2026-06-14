DEGORA
======

DEGORA turns heterogeneous published differential-expression gene (DEG)
tables into a local, inspectable gene-evidence database. It is designed for
users who have DEG result tables and want a reproducible gene ranking plus
source-level evidence that can be reviewed in a browser.

The public repository contains the DEGORA software package, command-line
interface, demo data generator, dashboard server, and package tests.
Development-only logs, draft analyses, and generated publication packages are
not included in this public code snapshot.

.. image:: docs/assets/degora-dashboard-snapshot.png
   :alt: DEGORA dashboard snapshot

Supported Environment
---------------------

Supported paths:

* Linux or WSL/Ubuntu.
* macOS 12 or newer on Apple Silicon, including M1, M2, M3, and M4 Macs.

Native Windows is not the supported path because scientific Python wheels and
shell behavior can differ, especially on Windows-ARM machines.

Recommended Linux/WSL baseline:

* Ubuntu 22.04 or newer, including WSL2 on Windows
* Python 3.10 or 3.11
* ``git``, ``python3-venv``, and ``make``

Recommended macOS baseline:

* macOS 12 or newer on Apple Silicon
* Command Line Tools for Xcode, installed with ``xcode-select --install`` if
  ``git`` or ``make`` is missing
* Python 3.10 or 3.11 from Homebrew, pyenv, or Miniforge
* A native arm64 terminal/Python environment, not a mixed Intel/Rosetta
  environment

Quick Start
-----------

One-time setup:

.. code-block:: bash

   git clone https://github.com/kangk1204/DEGORA.git
   cd DEGORA

   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -e outputs/code

After ``source .venv/bin/activate``, your terminal prompt usually starts with
``(.venv)``. That means the DEGORA environment is active.

Create and run the demo workspace:

.. code-block:: bash

   degora demo degora-demo
   degora validate degora-demo/degora_demo_config.xlsx
   degora run degora-demo/degora_demo_config.xlsx
   degora serve degora-demo/results/degora_scores.db

Open the printed local URL in your browser. The dashboard lets you search
genes, inspect ranked scores, and review the per-source evidence behind each
gene.

Using DEGORA Again Later
------------------------

If you close the terminal after installing DEGORA, you do not need to install
the package again. Open a new terminal, go back to the repository folder, and
reactivate the same environment:

.. code-block:: bash

   cd DEGORA
   source .venv/bin/activate
   degora --help

If ``degora --help`` prints the command help, you are ready to continue. To
reopen the demo dashboard from an existing demo run:

.. code-block:: bash

   cd DEGORA
   source .venv/bin/activate
   degora serve degora-demo/results/degora_scores.db

If the terminal says ``degora: command not found``, the environment is usually
not active. Run ``source .venv/bin/activate`` from inside the DEGORA repository
and try again.

Preparing Your Own Input
------------------------

Create an Excel template:

.. code-block:: bash

   degora template DEGORA_template.xlsx

Each template sheet starts with a short ``#`` note in the first row. Leave
those note rows in place; DEGORA ignores them when reading the workbook. Fill
the ``Contrasts`` sheet with one row per DEG table.

Required fields:

* ``study_id``: contrast identifier, for example
  ``GSE12345_treated_vs_control``.
* ``source_path``: local path to the DEG table file.
* ``gene_column``: column in the source table containing gene symbols.
* ``lfc_column``: column containing signed effect size, usually log2 fold
  change.
* ``p_column``: nominal p-value column. Values must be in ``[0, 1]``.

Optional fields:

* ``padj_column``: adjusted p-value/FDR column. If supplied, values must be in
  ``[0, 1]``.
* ``source_unit_id``: independent source unit, useful when multiple contrasts
  come from one dataset or paper.
* ``paper_id``: publication or dataset group identifier.
* ``n_ctrl`` and ``n_treat``: control and treated/case sample counts.
* ``assay_type``, ``platform``, ``species``, and ``cell_system``: metadata for
  evidence review.
* ``source_url``: source URL for provenance.

Template workbook sheets:

* ``README``: short in-workbook instructions.
* ``Project``: project-level settings such as output folders and the minimum
  number of independent source units required for scoring.
* ``Contrasts``: the main input sheet. This is where you list your DEG tables
  and map each source table's gene, effect-size, and p-value columns.
* ``GoldPanel``: optional known-positive genes for checking recall metrics.
  It is not required for DEGORA scoring, and it does not change the ranked
  gene results. Beginners can leave it empty. If you have a locked marker
  panel for your topic, put gene symbols in ``gene_symbol``; optional columns
  such as ``expected_direction``, ``role``, and ``evidence_basis`` are for your
  own notes. If a ``locked`` column is present, values such as ``yes``,
  ``true``, ``1``, or a blank cell are included in recall calculations.
* ``AdvancedSettings``: optional settings for users who need to change default
  behavior.
* ``ColumnGuide``: examples and explanations for common column names.

Outputs
-------

A DEGORA run writes:

* ``degora_scores.db``: SQLite evidence database used by the dashboard.
* ``degora_gene_scores.csv``: ranked gene table with DEGORA score components.
* ``degora_evidence.csv``: per-source evidence rows behind each gene.
* ``degora_score_metadata.json``: run metadata and input summary.

Development Checks
------------------

Install the test extra before running the package checks:

.. code-block:: bash

   python -m pip install -e 'outputs/code[dev]'
   make check

If you run pytest directly from the repository root instead of using
``make check``, include ``outputs/code`` on ``PYTHONPATH`` so Python can import
the local package:

.. code-block:: bash

   PYTHONPATH=outputs/code python -m pytest -q outputs/code/tests

This runs the public unit tests and byte-compiles the package. The public test
suite is intentionally limited to software/package behavior and does not
require analysis-only data, generated result files, R/Bioconductor packages, or
local artifacts.

License
-------

MIT License. See ``LICENSE``.
