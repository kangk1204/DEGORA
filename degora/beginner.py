"""Guided catalog creation for readers who are not bioinformaticians.

`degora template` hands over a spreadsheet and leaves the reader to fill it in.
That asks them to decide, per DEG table, which column holds the effect, whether it
is on a log2 scale, whether the table lists every gene tested or only the
significant ones, and which direction the comparison runs - judgements that need
the very expertise the tool is meant to stand in for.

This module infers what can be inferred from the file itself, shows the evidence
for each inference, and asks only about what a file cannot tell us. One thing is
never inferred: the contrast direction. Reversing it flips every gene's up/down
call while leaving results that look entirely reasonable, so there is nothing
downstream for a reader to notice. It is asked, and the answer is recorded.
"""

from __future__ import annotations

import os
import re
import tempfile
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .discovery import (
    NON_GENE_IDENTIFIER_RE,
    classify_header,
    is_gene_identifier_header,
)
from .excel_export import _force_formula_like_text
from .harmonize import (
    TableMapping,
    _is_binary_probability_indicator,
    assess_table_scope,
    read_deg_table,
)

# Extensions read_deg_table knows how to open.
SOURCE_TABLE_SUFFIXES = (
    ".csv",
    ".tsv",
    ".txt",
    ".xlsx",
    ".xls",
    ".csv.gz",
    ".tsv.gz",
    ".txt.gz",
    ".xlsx.gz",
    ".xls.gz",
)

# How the reader is asked about contrast direction. The wording avoids "log2FC
# sign" and "reference level": it asks about the experiment, not the file.
DIRECTION_QUESTION = (
    "In this table, does a POSITIVE value mean the gene went UP in the treated "
    "(or disease, or knockout) samples compared with the control samples?"
)
DIRECTION_HELP = (
    "This is the one thing DEGORA cannot work out for you, and getting it backwards "
    "silently inverts every up/down call in the results. If you are unsure, check "
    "the paper or the analysis script for which group was the reference."
)


@dataclass(frozen=True)
class ColumnChoice:
    """One column DEGORA needs, with what it found and what else it could be."""

    role: str
    chosen: str
    alternatives: tuple[str, ...] = ()
    confident: bool = True
    note: str = ""


@dataclass(frozen=True)
class SourceInference:
    """What one DEG table could be read off without asking anybody."""

    path: Path
    n_rows: int
    columns: tuple[str, ...] = ()
    choices: tuple[ColumnChoice, ...] = ()
    table_scope: str = "auto"
    table_scope_reason: str = ""
    plausible: dict[str, tuple[str, ...]] = field(default_factory=dict)
    # Keyed by column name, because the reader may override the gene column and
    # the identifier space has to describe the column actually being used.
    identifier_space_by_column: dict[str, str] = field(default_factory=dict)

    @property
    def identifier_space(self) -> str:
        """The space of the auto-detected gene column, for the walkthrough line."""

        return self.identifier_space_by_column.get(self.mapping.get("gene_column", ""), "")

    def identifier_space_for(self, column: str) -> str:
        return self.identifier_space_by_column.get(column, "")
    readable: bool = True
    problem: str = ""

    @property
    def mapping(self) -> dict[str, str]:
        return {choice.role: choice.chosen for choice in self.choices if choice.chosen}

    @property
    def looks_like_a_deg_table(self) -> bool:
        """Whether this file has anything a DEG table must have.

        A sample sheet or a README sitting in the same folder should be skipped
        rather than walked through question by question.
        """

        if not self.readable:
            return False
        mapping = self.mapping
        classified_roles = sum(
            bool(value)
            for value in (
                mapping.get("gene_column"),
                mapping.get("lfc_column"),
                mapping.get("p_column") or mapping.get("padj_column"),
            )
        )
        if classified_roles < 2:
            return False
        has_gene = bool(mapping.get("gene_column") or self.plausible.get("gene_column"))
        has_effect = bool(mapping.get("lfc_column") or self.plausible.get("lfc_column"))
        has_significance = bool(mapping.get("p_column") or mapping.get("padj_column") or self.plausible.get("p_column"))
        return has_gene and has_effect and has_significance

    @property
    def needs_a_question(self) -> tuple[ColumnChoice, ...]:
        """Column choices a reader has to settle, because the file is ambiguous."""

        return tuple(choice for choice in self.choices if not choice.confident)


@dataclass
class ContrastAnswers:
    """What a reader confirmed about one table. None of this is inferred."""

    positive_means_up_in_treated: bool
    condition: str = ""
    species: str = ""
    source_unit_id: str = ""
    n_ctrl: str = ""
    n_treat: str = ""
    overrides: dict[str, str] = field(default_factory=dict)


def find_source_tables(directory: str | Path) -> list[Path]:
    """Return the DEG-table-shaped files under a directory, deepest name order."""

    root = Path(directory)
    if not root.is_dir():
        return []
    found: list[Path] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name.startswith("."):
            continue
        name = path.name.lower()
        if any(name.endswith(suffix) for suffix in SOURCE_TABLE_SUFFIXES):
            found.append(path)
    return found


def _read_header(path: Path) -> tuple[pd.DataFrame, str]:
    """Read a table for inspection, reporting why if it cannot be read."""

    try:
        frame = read_deg_table(path, TableMapping(gene_column="", lfc_column="", p_column=""))
    except Exception as exc:  # noqa: BLE001 - one unreadable file must not end the walk
        return pd.DataFrame(), f"{type(exc).__name__}: {exc}"
    return frame, ""


# Enough rows to tell a p-value column from a fold change, few enough that a
# 58,000-row table costs nothing to inspect.
PLAUSIBILITY_SAMPLE_ROWS = 2000
# A list longer than this stops being a choice and becomes something to scroll past.
MAX_OPTIONS_SHOWN = 12
# Mirrors the CLI's default Project.min_studies; a config below it scores nothing.
DEFAULT_MIN_SOURCE_UNITS = 2
CATALOG_MARKER_COLUMNS = frozenset({"study_id", "source_path", "gene_column", "lfc_column", "p_column"})


IDENTIFIER_PATTERNS = (
    ("Ensembl ID", re.compile(r"^ENS[A-Z]*[GT]\d{6,}(\.\d+)?$", re.IGNORECASE)),
    ("RefSeq ID", re.compile(r"^[NX][MRP]_\d+(\.\d+)?$", re.IGNORECASE)),
    ("Affymetrix probe ID", re.compile(r"^\d+_[a-z]?_?at$", re.IGNORECASE)),
    ("Entrez ID", re.compile(r"^\d+$")),
    ("gene symbol", re.compile(r"^[A-Za-z][A-Za-z0-9\-.@_]{0,24}$")),
)
IDENTIFIER_SAMPLE_ROWS = 200
UNKNOWN_IDENTIFIER_SPACE = "unrecognised identifiers"
PROBABILITY_COLUMN_NAME_RE = re.compile(
    r"(^|[^A-Za-z0-9])(adj[._ -]?)?p([._ -]?(value|val))($|[^A-Za-z0-9])|"
    r"^p$|^p[._ -](?!vs($|[._ -])).+|"
    r"(^|[^A-Za-z0-9])p[._ -]?adj(ust(ed)?)?($|[^A-Za-z0-9])|"
    r"q[._ -]?value|fdr|false[._ -]?discovery",
    re.IGNORECASE,
)
BARE_IDENTIFIER_RE = re.compile(r"^(?:id|identifier)$", re.IGNORECASE)
SAMPLE_IDENTIFIER_NAME_RE = re.compile(
    r"(?:^|[_ .-])(?:sample|subject|patient|donor|run|file|series)(?:[_ .-]?(?:id|name|accession|identifier))?$|"
    r"^(?:gsm|srr|err|drr)(?:[_ .-]?(?:id|accession|identifier))?$",
    re.IGNORECASE,
)
SAMPLE_ACCESSION_RE = re.compile(r"^(?:GSM|SRR|ERR|DRR)\d+$", re.IGNORECASE)


class BeginnerInitError(ValueError):
    """A reader-correctable reason the guided initializer could not write a config."""


def identifier_space(values: Iterable[Any]) -> str:
    """Name the identifier space a gene column is written in.

    DEGORA matches genes across studies by the identifier itself, so a table of
    Ensembl IDs and a table of symbols have nothing in common even when they
    describe the same genes. Recognising that while the config is being built is
    the difference between a sentence now and a run that scores zero genes later.
    """

    counts: dict[str, int] = {}
    seen = 0
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        seen += 1
        for label, pattern in IDENTIFIER_PATTERNS:
            if pattern.match(text):
                counts[label] = counts.get(label, 0) + 1
                break
        if seen >= IDENTIFIER_SAMPLE_ROWS:
            break
    if not counts or not seen:
        return UNKNOWN_IDENTIFIER_SPACE
    label, hits = max(counts.items(), key=lambda item: item[1])
    # A clear majority, or the column is too mixed to name honestly.
    return label if hits >= seen * 0.7 else UNKNOWN_IDENTIFIER_SPACE


def _can_be_gene_label_column(name: Any, values: Iterable[Any]) -> bool:
    """Conservative value fallback for gene/probe identifier columns."""

    name_text = str(name).strip()
    if SAMPLE_IDENTIFIER_NAME_RE.search(name_text):
        return False
    space = identifier_space(values)
    if space == UNKNOWN_IDENTIFIER_SPACE:
        return is_gene_identifier_header(name, loose=True)
    seen = 0
    sample_accessions = 0
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        seen += 1
        if SAMPLE_ACCESSION_RE.match(text):
            sample_accessions += 1
        if seen >= IDENTIFIER_SAMPLE_ROWS:
            break
    if seen and sample_accessions >= seen * 0.7:
        return False
    if NON_GENE_IDENTIFIER_RE.search(name_text) and not BARE_IDENTIFIER_RE.fullmatch(name_text):
        return False
    return True


def _rescue_gene_mapping_from_values(header: dict[str, Any], plausible: dict[str, tuple[str, ...]]) -> None:
    """Use value-recognised identifiers when a generic header hides the gene column."""

    if header["mapping"].get("gene_column"):
        return
    candidates = tuple(plausible.get("gene_column") or ())
    if len(candidates) != 1:
        return
    candidate = candidates[0]
    header["gene_columns"] = [candidate]
    header["mapping"]["gene_column"] = candidate


def _plausible_columns(frame: pd.DataFrame) -> dict[str, tuple[str, ...]]:
    """Which columns could hold each role, judged by their values, not their names.

    When the header classifier recognises nothing, the reader used to be offered
    every column in the file - 43 of them for one real table, 32 being per-sample
    expression values that could not be a gene name under any reading. A list
    that long is not a choice, it is a wall. Values narrow it honestly: a p-value
    lies in [0, 1], a fold change is numeric, and a gene name is not a number.
    """

    sample = frame.head(PLAUSIBILITY_SAMPLE_ROWS)
    probability: list[str] = []
    numeric: list[str] = []
    labels: list[str] = []
    for name in frame.columns:
        column = pd.to_numeric(sample[name], errors="coerce")
        finite = column[np.isfinite(column)] if len(column) else column
        share_numeric = (len(finite) / len(sample)) if len(sample) else 0.0
        if share_numeric >= 0.5:
            numeric.append(str(name))
            distinct = sample[name].nunique(dropna=True)
            if (
                identifier_space(sample[name]) == "Entrez ID"
                and not NON_GENE_IDENTIFIER_RE.search(str(name).strip())
                and len(finite)
                and distinct / len(finite) >= 0.5
            ):
                labels.append(str(name))
            binary_indicator = _is_binary_probability_indicator(finite, column=str(name))
            if (
                len(finite)
                and not binary_indicator
                and float(finite.min()) >= 0.0
                and float(finite.max()) <= 1.0
            ):
                probability.append(str(name))
        else:
            # A gene column has to distinguish rows; a constant text column cannot.
            distinct = sample[name].astype(str).nunique(dropna=True)
            if distinct > max(1, len(sample) // 100) and _can_be_gene_label_column(name, sample[name]):
                labels.append(str(name))
    lfc_candidates = [name for name in numeric if not PROBABILITY_COLUMN_NAME_RE.search(name)]
    if not lfc_candidates:
        lfc_candidates = numeric
    return {
        "gene_column": tuple(labels),
        "lfc_column": tuple(lfc_candidates),
        "p_column": tuple(probability),
        "padj_column": tuple(probability),
    }


def _mapping_problem(
    mapping: dict[str, str],
    columns: Sequence[str],
    *,
    require_complete: bool = True,
) -> str:
    """Return why a user-provided mapping is not usable, or an empty string."""

    known = set(columns)
    labels = {
        "gene_column": "gene column",
        "lfc_column": "effect-size column",
        "p_column": "p-value column",
        "padj_column": "adjusted p-value column",
    }
    for role in ("gene_column", "lfc_column", "p_column"):
        value = mapping.get(role, "")
        if not value:
            if require_complete:
                return f"the {labels[role]} is required"
            continue
        if value not in known:
            return f"{value!r} is not a column in this table"
    padj = mapping.get("padj_column", "")
    if padj and padj not in known:
        return f"{padj!r} is not a column in this table"

    occupied: dict[str, str] = {}
    for role in ("gene_column", "lfc_column", "p_column", "padj_column"):
        value = mapping.get(role, "")
        if not value:
            continue
        if value in occupied:
            previous = occupied[value]
            if {previous, role} == {"p_column", "padj_column"}:
                continue
            return f"{value!r} cannot be used as both {labels[previous]} and {labels[role]}"
        occupied[value] = role
    return ""


def _ask_column_choice(
    *,
    choice: ColumnChoice,
    inference: SourceInference,
    overrides: dict[str, str],
    ask: Ask,
    echo: Callable[[str], None],
) -> str | None:
    """Ask for one ambiguous mapping, bounded so a bad answer cannot hang init."""

    role_label = choice.role.replace("_", " ")
    options = list(
        dict.fromkeys(
            name
            for name in (choice.chosen, *choice.alternatives, *inference.plausible.get(choice.role, ()))
            if name
        )
    )
    if not options:
        # Nothing matched by name or values. Keeping the bounded full-list fallback
        # gives unusual but real exports one recoverable path.
        options = list(inference.columns)
    echo(f"  DEGORA is not sure which column holds the {role_label}.")
    if choice.note:
        echo(f"  {choice.note}.")
    shown = options[:MAX_OPTIONS_SHOWN]
    hidden = len(options) - len(shown)
    suffix = f", and {hidden} more (type another listed choice)" if hidden > 0 else ""
    echo(f"  Options: {', '.join(shown)}{suffix}")

    last_problem = ""
    for attempt in range(1, 4):
        picked = ask(f"  Which column is the {role_label}?", choice.chosen).strip()
        candidate = dict(inference.mapping)
        candidate.update(overrides)
        if picked:
            candidate[choice.role] = picked
        last_problem = _mapping_problem(candidate, inference.columns, require_complete=False)
        if not last_problem and not picked and not choice.chosen:
            last_problem = f"the {role_label} is required"
        if not last_problem and picked and picked not in options:
            last_problem = f"{picked!r} is not one of the available {role_label} choices"
        if not last_problem:
            return picked if picked else choice.chosen
        remaining = 3 - attempt
        if remaining:
            echo(f"  That answer cannot be used: {last_problem}. Try again ({remaining} left).")
    echo(f"  -> skipped; column mapping was not usable: {last_problem}.")
    return None


def _column_choices(header: dict[str, Any]) -> tuple[ColumnChoice, ...]:
    mapping = header["mapping"]
    choices: list[ColumnChoice] = []
    for role, candidates_key, required in (
        ("gene_column", "gene_columns", True),
        ("lfc_column", "lfc_columns", True),
        ("p_column", "p_columns", True),
        ("padj_column", "padj_columns", False),
    ):
        chosen = mapping.get(role, "")
        candidates = tuple(header.get(candidates_key) or ())
        alternatives = tuple(name for name in candidates if name != chosen)
        note = ""
        confident = True
        if not chosen:
            if required:
                note = "no column in this file looks like it"
                if role == "p_column" and mapping.get("padj_column"):
                    # Worth saying rather than leaving to be discovered: the table
                    # has adjusted p-values only, which is a usable answer and a
                    # consequence the reader should choose knowingly.
                    note = (
                        f"this table has no unadjusted p-value; only "
                        f"{mapping['padj_column']} is available, and using it "
                        f"means every p-value DEGORA reads is already adjusted"
                    )
                confident = False
            # An absent optional column is an answer, not a question: DEGORA runs
            # without an adjusted p-value and asking for one it cannot see wastes
            # the reader's attention on the tables that need none of it.
        elif alternatives:
            note = "more than one column could be this"
            confident = False
        elif role == "lfc_column" and chosen and not header.get("lfc_scale_explicit"):
            note = "the column name does not say the values are on a log2 scale"
            confident = False
        choices.append(
            ColumnChoice(
                role=role,
                chosen=chosen,
                alternatives=alternatives,
                confident=confident,
                note=note,
            )
        )
    return tuple(choices)


def infer_source_table(path: str | Path) -> SourceInference:
    """Read one DEG table and report what it says about itself.

    Pure and side-effect free: everything the guided flow decides is derived from
    this, so the inference can be tested without driving a conversation.
    """

    path = Path(path)
    frame, problem = _read_header(path)
    if problem:
        return SourceInference(path=path, n_rows=0, readable=False, problem=problem)

    literal_columns = {str(name).strip().lower() for name in frame.columns}
    if CATALOG_MARKER_COLUMNS.issubset(literal_columns):
        return SourceInference(
            path=path,
            n_rows=len(frame),
            columns=tuple(str(name) for name in frame.columns),
            problem="this is a DEGORA config catalog, not a source DEG results table",
        )
    if len(frame) == 0:
        return SourceInference(
            path=path,
            n_rows=0,
            columns=tuple(str(name) for name in frame.columns),
            problem="this table has column headers but no data rows",
        )

    plausible = _plausible_columns(frame)
    header = classify_header(frame.columns)
    _rescue_gene_mapping_from_values(header, plausible)
    thresholded_probability_columns: list[str] = []
    for role, candidates_key in (("p_column", "p_columns"), ("padj_column", "padj_columns")):
        candidates = list(header.get(candidates_key) or ())
        thresholded = [
            name
            for name in candidates
            if name in frame.columns and _is_binary_probability_indicator(frame[name], column=name)
        ]
        if not thresholded:
            continue
        thresholded_probability_columns.extend(thresholded)
        header[candidates_key] = [name for name in candidates if name not in thresholded]
        if header["mapping"].get(role) in thresholded:
            header["mapping"][role] = ""
    choices = _column_choices(header)
    mapping = {choice.role: choice.chosen for choice in choices}

    scope_label, scope_reason = "auto", "not enough of a p-value column to tell"
    if mapping.get("p_column"):
        scope = assess_table_scope(
            frame,
            TableMapping(
                gene_column=mapping.get("gene_column", ""),
                lfc_column=mapping.get("lfc_column", ""),
                p_column=mapping["p_column"],
                padj_column=mapping.get("padj_column") or None,
            ),
        )
        scope_label = str(scope["effective_scope"])
        scope_reason = str(scope["reason"])

    probability_problem = ""
    if thresholded_probability_columns and not (
        mapping.get("p_column")
        or mapping.get("padj_column")
        or plausible.get("p_column")
    ):
        probability_problem = (
            f"{', '.join(thresholded_probability_columns)} contains only a binary significance flag; "
            "DEGORA needs unrounded gene-level p-values or adjusted p-values"
        )

    return SourceInference(
        path=path,
        n_rows=len(frame),
        columns=tuple(str(name) for name in frame.columns),
        choices=choices,
        table_scope=scope_label,
        table_scope_reason=scope_reason,
        plausible=plausible,
        identifier_space_by_column={
            str(name): identifier_space(frame[name]) for name in frame.columns
        },
        problem=probability_problem,
    )


def describe_inference(inference: SourceInference) -> list[str]:
    """Render the inference as lines a non-specialist can check."""

    if not inference.readable:
        return [f"{inference.path.name}: could not be read - {inference.problem}"]
    lines = [f"{inference.path.name}: {inference.n_rows:,} rows, {len(inference.columns)} columns"]
    labels = {
        "gene_column": "gene names",
        "lfc_column": "effect size (log2 fold change)",
        "p_column": "p-value",
        "padj_column": "adjusted p-value",
    }
    for choice in inference.choices:
        label = labels[choice.role]
        if not choice.chosen:
            if choice.role == "padj_column":
                lines.append(f"  {label}: not found (optional)")
            else:
                lines.append(f"  {label}: NOT FOUND - {choice.note}")
            continue
        suffix = ""
        if choice.alternatives:
            suffix = f"  (could also be: {', '.join(choice.alternatives)})"
        elif choice.note:
            suffix = f"  ({choice.note})"
        lines.append(f"  {label}: {choice.chosen}{suffix}")
    if inference.table_scope_reason:
        lines.append(f"  table covers: {inference.table_scope} - {inference.table_scope_reason}")
    if inference.identifier_space:
        lines.append(f"  gene names are written as: {inference.identifier_space}")
    if inference.problem:
        lines.append(f"  cannot use this table: {inference.problem}")
    return lines


def catalog_row(
    inference: SourceInference,
    answers: ContrastAnswers,
    *,
    study_id: str,
    catalog_dir: str | Path,
) -> dict[str, Any]:
    """Turn one inference plus one reader's answers into a catalog row.

    The direction answer is not translated into a sign flip: DEGORA never reverses
    an effect column. A table whose positive values mean "up in control" has to be
    corrected at the source, and the row records that it was not usable as-is.
    """

    mapping = dict(inference.mapping)
    mapping.update({role: value for role, value in answers.overrides.items() if value})
    try:
        source_path = inference.path.resolve().relative_to(Path(catalog_dir).resolve()).as_posix()
    except ValueError:
        source_path = str(inference.path.resolve())
    table_scope = inference.table_scope
    if (
        mapping.get("p_column", "") != inference.mapping.get("p_column", "")
        or mapping.get("padj_column", "") != inference.mapping.get("padj_column", "")
    ):
        table_scope = "auto"

    return {
        "study_id": study_id,
        "source_unit_id": answers.source_unit_id or study_id,
        "source_path": source_path,
        "gene_column": mapping.get("gene_column", ""),
        "lfc_column": mapping.get("lfc_column", ""),
        "p_column": mapping.get("p_column", ""),
        "padj_column": mapping.get("padj_column", ""),
        "species": answers.species,
        "condition": answers.condition,
        "n_ctrl": answers.n_ctrl,
        "n_treat": answers.n_treat,
        "table_scope": table_scope,
        "sign_convention": (
            "confirmed_treatment_minus_control"
            if answers.positive_means_up_in_treated
            else "REVERSED - do not analyse until the source table is corrected"
        ),
        "include_in_analysis": "yes" if answers.positive_means_up_in_treated else "no",
        "notes": (
            ""
            if answers.positive_means_up_in_treated
            else "Excluded by degora init: a positive value in this table means up in the control group. "
            "DEGORA never reverses an effect column, so correct the table at source and set "
            "include_in_analysis to yes."
        ),
    }


def build_catalog(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    """Return the rows as a catalog frame in a stable column order."""

    columns = [
        "study_id",
        "source_unit_id",
        "source_path",
        "gene_column",
        "lfc_column",
        "p_column",
        "padj_column",
        "species",
        "condition",
        "n_ctrl",
        "n_treat",
        "table_scope",
        "sign_convention",
        "include_in_analysis",
        "notes",
    ]
    frame = pd.DataFrame(list(rows))
    for column in columns:
        if column not in frame.columns:
            frame[column] = ""
    return frame[columns]


def default_study_id(path: Path, taken: Iterable[str]) -> str:
    """A readable, unique per-contrast id derived from the file name."""

    import re

    stem = re.sub(r"[^A-Za-z0-9]+", "_", path.stem).strip("_") or "contrast"
    stem = stem[:60]
    used = set(taken)
    if stem not in used:
        return stem
    suffix = 2
    while f"{stem}_{suffix}" in used:
        suffix += 1
    return f"{stem}_{suffix}"


Ask = Callable[[str, str], str]


def _prompt(question: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    try:
        answer = input(f"{question}{suffix}: ").strip()
    except EOFError:
        return default
    return answer or default


def _prompt_yes_no(
    question: str,
    help_text: str,
    ask: Ask,
    echo: Callable[[str], None] = print,
) -> bool | None:
    """Return True/False, or None when the reader will not commit.

    There is deliberately no default. The direction question is the one answer a
    reader must actually give, so pressing enter must not stand in for "yes".
    """

    for attempt in range(1, 4):
        answer = ask(f"{question} (yes / no / unsure)", "").strip().lower()
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        if answer in {"u", "unsure", "?", ""}:
            echo(f"  {help_text}")
            if ask("Skip this table for now? (yes/no)", "yes").strip().lower() in {"y", "yes", ""}:
                return None
        remaining = 3 - attempt
        if remaining:
            echo(f"  Please answer yes, no, or unsure ({remaining} tries left).")
    echo("  -> skipped; the contrast direction was not confirmed after three tries.")
    return None


def _write_catalog_atomic(catalog: pd.DataFrame, output: Path) -> None:
    """Publish a complete config without damaging an existing forced target."""

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.",
        suffix=output.suffix,
        dir=output.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        if output.suffix.lower() == ".xlsx":
            with pd.ExcelWriter(temporary, engine="openpyxl") as writer:
                catalog.to_excel(writer, sheet_name="Contrasts", index=False)
                _force_formula_like_text(writer)
        else:
            catalog.to_csv(temporary, index=False)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def run_init(
    output: str | Path,
    deg_dir: str | Path,
    *,
    ask: Ask | None = None,
    echo: Callable[[str], None] = print,
    force: bool = False,
) -> dict[str, Any]:
    """Walk a reader through building a catalog from a folder of DEG tables."""

    ask = ask or _prompt
    output = Path(output)
    output_suffix = output.suffix.lower()
    if output_suffix == ".xls":
        raise BeginnerInitError(
            "degora init cannot write the legacy .xls format; use a .xlsx or .csv output name"
        )
    if output.exists() and output.is_dir():
        raise BeginnerInitError("degora init output must be a file path, not an existing directory")
    if output_suffix not in {".csv", ".xlsx"}:
        raise BeginnerInitError(
            "degora init output must end in .csv or .xlsx so DEGORA and spreadsheet tools "
            "read the generated config in the format it was written"
        )
    if output.exists() and not force:
        raise FileExistsError(f"{output} already exists; pass --force to replace it")

    tables = [path for path in find_source_tables(deg_dir) if path.resolve() != output.resolve()]
    if not tables:
        raise FileNotFoundError(f"no CSV, TSV, TXT or Excel tables were found under {deg_dir}")
    inferences = [(path, infer_source_table(path)) for path in tables]
    empty_tables = [path.name for path, inference in inferences if inference.readable and inference.n_rows == 0]
    if empty_tables and not any(inference.looks_like_a_deg_table for _, inference in inferences):
        detail = ", ".join(empty_tables[:3])
        if len(empty_tables) > 3:
            detail += f", and {len(empty_tables) - 3} more"
        raise BeginnerInitError(f"degora init cannot use table(s) with no data rows: {detail}")

    echo(f"Found {len(tables)} file(s) under {deg_dir}.")
    echo("")
    rows: list[dict[str, Any]] = []
    spaces: dict[str, list[str]] = {}
    active_spaces: dict[str, list[str]] = {}
    skipped: list[dict[str, str]] = []
    taken: list[str] = []

    species = ask("Which species are these tables from? (human / mouse / other)", "human").strip()

    for path, inference in inferences:
        for line in describe_inference(inference):
            echo(line)
        if not inference.looks_like_a_deg_table:
            echo("  -> this does not look like a DEG results table; skipping.")
            echo("")
            skipped.append({"path": path.name, "reason": inference.problem or "not a DEG results table"})
            continue

        overrides: dict[str, str] = {}
        skip_reason = ""
        for choice in inference.needs_a_question:
            picked = _ask_column_choice(
                choice=choice,
                inference=inference,
                overrides=overrides,
                ask=ask,
                echo=echo,
            )
            if picked is None:
                skip_reason = "column mapping was not usable"
                break
            overrides[choice.role] = picked
        if skip_reason:
            echo("")
            skipped.append({"path": path.name, "reason": skip_reason})
            continue
        final_mapping = dict(inference.mapping)
        final_mapping.update({role: value for role, value in overrides.items() if value})
        final_mapping_problem = _mapping_problem(final_mapping, inference.columns)
        if final_mapping_problem:
            echo(f"  -> skipped; column mapping was not usable: {final_mapping_problem}.")
            echo("")
            skipped.append({"path": path.name, "reason": "column mapping was not usable"})
            continue
        if (
            final_mapping.get("p_column", "") != inference.mapping.get("p_column", "")
            or final_mapping.get("padj_column", "") != inference.mapping.get("padj_column", "")
        ):
            echo("  DEGORA will re-check whether this is a full or filtered table using the columns you selected.")

        echo("")
        echo(f"  {DIRECTION_QUESTION}")
        direction = _prompt_yes_no("  Answer", DIRECTION_HELP, ask, echo)
        if direction is None:
            echo("  -> skipped; nothing was guessed about its direction.")
            echo("")
            skipped.append({"path": path.name, "reason": "contrast direction not confirmed"})
            continue

        study_id = default_study_id(path, taken)
        taken.append(study_id)
        answers = ContrastAnswers(
            positive_means_up_in_treated=direction,
            condition=ask("  What was compared? (e.g. hypoxia vs normoxia)", ""),
            species=species,
            source_unit_id=ask("  Which paper or dataset is this from? (same answer groups related tables)", study_id),
            n_ctrl=ask("  How many control samples? (blank if unknown)", ""),
            n_treat=ask("  How many treated samples? (blank if unknown)", ""),
            overrides=overrides,
        )
        row = catalog_row(inference, answers, study_id=study_id, catalog_dir=output.parent)
        rows.append(row)
        chosen_gene_column = overrides.get("gene_column") or inference.mapping.get("gene_column", "")
        chosen_space = inference.identifier_space_for(chosen_gene_column) or UNKNOWN_IDENTIFIER_SPACE
        spaces.setdefault(chosen_space, []).append(path.name)
        if row["include_in_analysis"] == "yes":
            active_spaces.setdefault(chosen_space, []).append(path.name)
        echo("")

    if not rows:
        raise BeginnerInitError(
            "no table was confirmed, so there is nothing to write. Every table was either "
            "unreadable, not a DEG results table, skipped, or had an unusable column mapping."
        )

    catalog = build_catalog(rows)

    # DEGORA matches genes across studies on the identifier itself, so a run mixing
    # Ensembl IDs with symbols scores zero genes. That was only discovered at the end
    # of a full run; every table was read here, so it can be said now, while the
    # reader can still do something about it.
    mixed = [space for space in active_spaces if space != UNKNOWN_IDENTIFIER_SPACE]
    unknown = active_spaces.get(UNKNOWN_IDENTIFIER_SPACE, [])
    identifier_warning = ""
    if len(mixed) > 1:
        detail = "; ".join(
            f"{space}: {', '.join(sorted(active_spaces[space])[:3])}"
            + (f" and {len(active_spaces[space]) - 3} more" if len(active_spaces[space]) > 3 else "")
            for space in sorted(active_spaces)
        )
        identifier_warning = (
            f"These tables do not use one gene identifier space ({detail}). DEGORA matches "
            "genes across studies on the identifier itself, so sources written in different "
            "spaces share no genes and a run scores none. Convert them to one convention - "
            "all symbols, or all Ensembl IDs - before running."
        )
        echo("")
        echo(f"WARNING: {identifier_warning}")
    elif unknown:
        detail = ", ".join(sorted(unknown)[:3])
        if len(unknown) > 3:
            detail += f", and {len(unknown) - 3} more"
        identifier_warning = (
            f"DEGORA classified these columns as unrecognised identifiers and could not verify one "
            f"gene identifier space for: {detail}. They may mix "
            "symbols, Ensembl/Entrez IDs, probes, or non-gene identifiers. Inspect them and convert "
            "all active tables to one convention before relying on cross-study matches."
        )
        echo("")
        echo(f"WARNING: {identifier_warning}")

    # The replication rule needs two independent source units. Five tables from one
    # study is a config that cannot score a gene, and the reader is standing right
    # here - after the run is the expensive place to learn it.
    active_units = {row["source_unit_id"] for row in rows if row["include_in_analysis"] == "yes"}
    replication_warning = ""
    if len(active_units) < DEFAULT_MIN_SOURCE_UNITS:
        replication_warning = (
            f"These tables come from {len(active_units)} independent source unit(s). DEGORA's "
            f"default replication rule needs {DEFAULT_MIN_SOURCE_UNITS}, so a run over this "
            "config scores zero genes. Add a table from another study, or give each independent "
            "study its own answer to the 'which paper or dataset' question."
        )
        echo("")
        echo(f"WARNING: {replication_warning}")

    output.parent.mkdir(parents=True, exist_ok=True)
    _write_catalog_atomic(catalog, output)

    reversed_rows = [row for row in rows if row["include_in_analysis"] == "no"]
    active_rows = [row for row in rows if row["include_in_analysis"] == "yes"]
    return {
        "replication_warning": replication_warning,
        "identifier_spaces": {space: sorted(names) for space, names in spaces.items()},
        "identifier_warning": identifier_warning,
        "config_path": str(output),
        "n_contrasts": len(rows),
        "n_source_units": len({row["source_unit_id"] for row in active_rows}),
        "n_excluded_reversed_direction": len(reversed_rows),
        "skipped": skipped,
    }
