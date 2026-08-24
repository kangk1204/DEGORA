"""Beginner-friendly command line interface for DEGORA."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

from . import format_version_info, runtime_version_info


DEFAULT_OUTPUT_DIR = Path("outputs/results/degora-run")
DEFAULT_HARMONIZED_DIR = Path("data/deg/harmonized")
WARNING_DISPLAY_LIMIT = 8
DISCOVERY_PAGE_SIZE = 10
DISCOVERY_MAX_LIMIT = 1000


class CliUsageError(ValueError):
    """User-facing CLI input error that should print without a traceback."""


def _setting_key(value: Any) -> str:
    return str(value).strip().lower().replace(" ", "_").replace("-", "_")


def _sheet_settings(path: Path, sheet: str, key_column: str) -> dict[str, str]:
    import pandas as pd

    from .excel_io import read_config_sheet

    try:
        frame = read_config_sheet(path, sheet)
    except Exception:  # noqa: BLE001 - a malformed workbook is reported by the catalog reader.
        # A file that is a valid ZIP but not an OOXML workbook raises engine-selection
        # errors from pandas/openpyxl that are neither ValueError nor a DEGORA error.
        # Optional Project/AdvancedSettings sheets must never be the thing that turns a
        # bad config into a traceback; read_catalog raises the beginner-readable error.
        return {}
    if frame.empty:
        return {}
    columns = list(frame.columns)
    if key_column in frame.columns and "value" in frame.columns:
        key_name = key_column
        value_name = "value"
    elif len(columns) >= 2:
        key_name = columns[0]
        value_name = columns[1]
    else:
        return {}
    settings: dict[str, str] = {}
    for _, row in frame.iterrows():
        key = _setting_key(row.get(key_name, ""))
        if not key:
            continue
        value = row.get(value_name)
        if value is None or pd.isna(value):
            continue
        settings[key] = str(value).strip()
    return settings


def read_excel_settings(path: Path) -> dict[str, str]:
    """Read Project and AdvancedSettings sheets when the config is an Excel workbook."""

    if path.suffix.lower() not in {".xlsx", ".xls"}:
        return {}
    settings = _sheet_settings(path, "Project", "field")
    settings.update(_sheet_settings(path, "AdvancedSettings", "setting"))
    return settings


def _validate_score_version(settings: dict[str, str]) -> None:
    """Reject a workbook that explicitly requests an unsupported score contract."""

    configured = settings.get("score_version", "").strip()
    if not configured:
        return
    from . import SCORE_VERSION

    if configured == SCORE_VERSION:
        return
    from .slice_runner import DegoraConfigError

    raise DegoraConfigError(
        "score_version is not supported by this DEGORA installation",
        problems=[f"Workbook requests {configured!r}, but this installation implements {SCORE_VERSION!r}."],
        fixes=[
            f"Set AdvancedSettings.score_version to {SCORE_VERSION!r}, or generate a fresh template with `degora template`."
        ],
    )


def _int_setting(value: Any | None, default: int) -> int:
    if value is None or value == "":
        return default
    try:
        from .aggregate import validate_min_studies

        return validate_min_studies(value)
    except ValueError as exc:
        from .slice_runner import DegoraConfigError

        raise DegoraConfigError(
            "numeric setting is invalid",
            problems=[f"Expected a whole number, but got {value!r}."],
            fixes=["Open the Project or AdvancedSettings sheet and enter a number such as 2."],
        ) from exc


def _bounded_discovery_limit(value: str) -> int:
    try:
        limit = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--limit must be a whole number from 1 to 1000") from exc
    if not 1 <= limit <= DISCOVERY_MAX_LIMIT:
        raise argparse.ArgumentTypeError("--limit must be between 1 and 1000")
    return limit


def _normalize_publication_selection(value: Any) -> str:
    text = str(value or "").strip()
    lowered = text.lower()
    if lowered.startswith("pmid:"):
        return "pmid:" + "".join(ch for ch in text.split(":", 1)[1] if ch.isdigit())
    if lowered.startswith("doi:"):
        return "doi:" + text.split(":", 1)[1].strip().lower().removeprefix("https://doi.org/")
    if lowered.startswith("gse") and text[3:].isdigit():
        return "gse" + text[3:]
    return lowered


def _publication_selection_keys(record: dict[str, Any]) -> set[str]:
    keys: set[str] = set()
    for field in ("canonical_id", "source_unit_id"):
        value = str(record.get(field) or "").strip()
        if value:
            keys.add(_normalize_publication_selection(value))
            keys.add(value.lower())
    for value in _as_list(record.get("pubmed_ids")) + _as_list(record.get("pmid")):
        normalized = _normalize_publication_selection(f"PMID:{value}")
        if normalized != "pmid:":
            keys.add(normalized)
    for value in _as_list(record.get("doi")) + _as_list(record.get("dois")):
        normalized = _normalize_publication_selection(f"DOI:{value}")
        if normalized != "doi:":
            keys.add(normalized)
    for value in _as_list(record.get("geo_accessions")) + _as_list(record.get("accession")):
        text = str(value or "").strip()
        if text:
            keys.add(_normalize_publication_selection(text))
            keys.add(text.lower())
    return {key for key in keys if key}


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if value in (None, ""):
        return []
    return [value]


def _select_publication_records(records: list[dict[str, Any]], selections: list[str]) -> list[dict[str, Any]]:
    normalized_selections = [_normalize_publication_selection(value) for value in selections]
    duplicates = sorted({value for value in normalized_selections if normalized_selections.count(value) > 1})
    if duplicates:
        raise CliUsageError("duplicate --select value(s): " + ", ".join(duplicates))

    matched: list[dict[str, Any]] = []
    for original, normalized in zip(selections, normalized_selections, strict=True):
        candidates = [record for record in records if normalized in _publication_selection_keys(record)]
        if not candidates:
            raise CliUsageError(f"--select did not match any publication/source-unit record: {original}")
        if len(candidates) > 1:
            identifiers = ", ".join(str(record.get("canonical_id") or record.get("source_unit_id") or "?") for record in candidates)
            raise CliUsageError(f"--select matched multiple records for {original}: {identifiers}")
        matched.append(candidates[0])
    return matched


def _print_publication_page(
    records: list[dict[str, Any]],
    *,
    page: int = 1,
    page_size: int = DISCOVERY_PAGE_SIZE,
) -> None:
    """Print one compact, reviewable page of publication rows."""

    start_index = (max(int(page), 1) - 1) * max(int(page_size), 1) + 1
    for index, record in enumerate(records, start=start_index):
        identifier = str(record.get("canonical_id") or record.get("source_unit_id") or "unidentified")
        title = " ".join(str(record.get("paper_title") or record.get("title") or "Untitled publication").split())
        if len(title) > 88:
            title = title[:85].rstrip() + "..."
        readiness = record.get("data_readiness")
        if isinstance(readiness, dict):
            readiness_text = readiness.get("verification_state") or readiness.get("tier") or "unknown"
        else:
            readiness_text = readiness or record.get("readiness") or "unknown"
        year = record.get("year") or "n.d."
        print(f"{index:>2}. {identifier} | {year} | {title} | {str(readiness_text).replace('_', ' ')}")


def _path_setting(value: str | None, default: Path, *, base: Path | None = None) -> Path:
    if value is None or value == "":
        path = default
    else:
        path = Path(value)
    if base is not None and value not in (None, "") and not path.is_absolute():
        return base / path
    return path


def _run_warning_messages(metrics: dict[str, Any]) -> list[str]:
    messages: list[str] = []
    seen: set[str] = set()

    warning_values: list[Any] = []
    for key in ("warnings", "identifier_space_warnings", "rank_universe_warnings"):
        warning_values.extend(metrics.get(key, []) or [])

    try:
        clipped_rows = int(metrics.get("pvalue_clipped_rows", 0) or 0)
    except (TypeError, ValueError):
        clipped_rows = 0
    if clipped_rows:
        text = (
            f"{clipped_rows} row(s) reported pvalue < 1e-300; values were floored "
            "to 1e-300 before signed-z scoring."
        )
        warning_values.append(text)

    for value in warning_values:
        text = str(value).strip()
        if text and text not in seen:
            messages.append(text)
            seen.add(text)

    return messages


def _print_run_warnings(
    metrics: dict[str, Any],
    *,
    metrics_path: Path,
    limit: int = WARNING_DISPLAY_LIMIT,
) -> None:
    messages = _run_warning_messages(metrics)
    if not messages:
        return

    print("", file=sys.stderr)
    print("DEGORA completed with non-fatal input warnings:", file=sys.stderr)
    for message in messages[:limit]:
        print(f"- {message}", file=sys.stderr)
    if len(messages) > limit:
        print(f"- ... {len(messages) - limit} more warning(s); see {metrics_path}", file=sys.stderr)
    else:
        print(f"See details: {metrics_path}", file=sys.stderr)


def _format_validation_items(items: Any) -> str:
    if not items:
        return "(none)"
    return ", ".join(map(str, items))


def _print_validation_summary(validation: dict[str, Any], *, include_excluded: bool = False) -> None:
    print(f"- Active contrasts: {validation['active_contrasts']}")
    if include_excluded:
        print(f"- Excluded contrasts: {validation['excluded_contrasts']}")
    print(f"- Independent source units: {validation['source_units']}")
    print(f"- Required Contrasts columns: {_format_validation_items(validation.get('required_contrasts_columns'))}")
    print(
        "- Required DEG-table mappings: "
        f"{_format_validation_items(validation.get('required_source_table_mappings'))}"
    )
    print(
        "- Optional DEG-table mappings checked when filled: "
        f"{_format_validation_items(validation.get('optional_source_table_mappings'))}"
    )


class _RunProgress:
    """Emit one timestamped line per phase of a run.

    A five-source human corpus takes minutes, and the command printed nothing at
    all until it finished. Lines are flushed so they survive redirection into a
    log, and carry elapsed seconds so a slow phase is identifiable.
    """

    def __init__(self, enabled: bool = True) -> None:
        self.enabled = enabled
        self._started = time.monotonic()
        self._phase_started = self._started
        self._phase = ""
        self._last_update = 0.0

    def start(self, phase: str) -> None:
        if not self.enabled:
            return
        self._phase = phase
        self._phase_started = time.monotonic()
        print(f"[{self._elapsed():>6.1f}s] {phase}...", flush=True)

    def done(self, detail: str = "") -> None:
        if not self.enabled or not self._phase:
            return
        took = time.monotonic() - self._phase_started
        suffix = f" - {detail}" if detail else ""
        print(f"[{self._elapsed():>6.1f}s] {self._phase} done in {took:.1f}s{suffix}", flush=True)
        self._phase = ""

    def update(self, detail: str, *, force: bool = False) -> None:
        if not self.enabled or not self._phase:
            return
        now = time.monotonic()
        if not force and now - self._last_update < 1.0:
            return
        self._last_update = now
        print(f"[{self._elapsed():>6.1f}s] {self._phase}: {detail}", flush=True)

    def _elapsed(self) -> float:
        return time.monotonic() - self._started


def _zero_gene_diagnostic(harmonized_path: Path, min_studies: int) -> str:
    """Explain why nothing scored, using the harmonized table that was written.

    The generic advice listed four possible causes. The commonest one by far -
    two sources that share no gene identifiers because they use different
    identifier spaces - can be stated as a fact instead of a guess.
    """

    try:
        import pandas as pd

        frame = pd.read_csv(harmonized_path, usecols=["gene_symbol", "source_unit_id"])
    except Exception:  # noqa: BLE001 - a diagnostic must never mask the real error.
        return ""
    if frame.empty:
        return " The harmonized table is empty, so no source contributed any usable row."

    by_unit = {
        str(unit): set(group["gene_symbol"].dropna().astype(str))
        for unit, group in frame.groupby("source_unit_id", sort=True)
    }
    if len(by_unit) < min_studies:
        return (
            f" Only {len(by_unit)} source unit(s) produced usable rows, but min_studies is"
            f" {min_studies}. Give each independent study its own source_unit_id, or lower"
            " Project.min_studies."
        )

    counts = ", ".join(f"{unit} ({len(genes):,})" for unit, genes in list(by_unit.items())[:6])
    shared = set.intersection(*by_unit.values()) if by_unit else set()
    if shared:
        return (
            f" {len(shared):,} identifier(s) are shared by all {len(by_unit)} source units, so the"
            " overlap is not the problem; check contrast direction and p-value columns."
            f" Identifiers per source unit: {counts}."
        )

    best_pair = ""
    units = list(by_unit.items())
    for index, (left_name, left) in enumerate(units):
        for right_name, right in units[index + 1 :]:
            if left & right:
                best_pair = (
                    f" The largest overlap between any two units is {len(left & right):,}"
                    f" ({left_name} vs {right_name})."
                )
                break
        if best_pair:
            break
    samples = []
    for unit, genes in list(by_unit.items())[:3]:
        example = sorted(genes)[0] if genes else ""
        samples.append(f"{unit}: {example!r}")
    return (
        f" No identifier is shared by all {len(by_unit)} source units.{best_pair}"
        f" Identifiers per source unit: {counts}. First identifier in each: {'; '.join(samples)}."
        " Sources must use the same identifier space - map symbols, Ensembl IDs and probe IDs"
        " onto one convention before running."
    )


def _run_from_config(args: argparse.Namespace, *, serve_after: bool = False) -> int:
    """Resolve the run's paths and settings, then hold the output directory for it."""

    from .provenance import output_directory_lock

    config = Path(args.config)
    config_base = config.resolve().parent
    settings = read_excel_settings(config)
    _validate_score_version(settings)
    min_studies = _int_setting(
        args.min_studies if args.min_studies is not None else settings.get("min_studies"),
        2,
    )
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else _path_setting(settings.get("output_dir"), DEFAULT_OUTPUT_DIR, base=config_base)
    )
    harmonized_dir = (
        Path(args.harmonized_dir)
        if args.harmonized_dir
        else _path_setting(settings.get("harmonized_dir"), DEFAULT_HARMONIZED_DIR, base=config_base)
    )
    db_path = Path(args.db) if args.db else _path_setting(settings.get("score_db"), output_dir / "degora_scores.db", base=config_base)

    # One claim for the whole pipeline. Harmonization and the database are written
    # tens of seconds apart, so two runs sharing an output directory could each
    # succeed and leave one run's contrast table beside the other's gene scores.
    with output_directory_lock(output_dir):
        return _run_pipeline(args, serve_after, config, settings, min_studies, output_dir, harmonized_dir, db_path)


def _run_pipeline(
    args: argparse.Namespace,
    serve_after: bool,
    config: Path,
    settings: dict[str, str],
    min_studies: int,
    output_dir: Path,
    harmonized_dir: Path,
    db_path: Path,
) -> int:
    from .api import serve as serve_db
    from .excel_export import DEFAULT_WORKBOOK_NAME, export_run_workbook
    from .provenance import shell_command
    from .score_db import write_score_database
    from .slice_runner import run_slice, validate_catalog_inputs

    version_info = runtime_version_info()
    progress = _RunProgress(enabled=not getattr(args, "quiet", False))

    def validation_progress(done: int, total: int, study_id: str) -> None:
        if total <= 0:
            return
        step = max(1, min(1000, total // 100 or 1))
        if done in {0, 1, total} or done % step == 0:
            suffix = f" ({study_id})" if study_id else ""
            progress.update(f"checked {done:,}/{total:,} active row(s){suffix}", force=done in {0, total})

    progress.start("Validating the catalog and source tables")
    validation = validate_catalog_inputs(config, progress=validation_progress)
    progress.done(f"{validation.get('active_contrasts', '?')} contrast(s)")
    print("DEGORA config OK")
    _print_validation_summary(validation)
    source_units = int(validation.get("source_units") or 0)
    if source_units < min_studies:
        raise CliUsageError(
            f"DEGORA cannot run at min_studies={min_studies}: this config has only "
            f"{source_units} independent source unit(s). A run would score zero genes. "
            "Give each independent study its own source_unit_id, add another study, or lower Project.min_studies."
        )

    progress.start("Harmonizing source tables")
    metrics = run_slice(config, output_dir, harmonized_dir, min_studies=min_studies)
    progress.done(f"{int(metrics.get('n_harmonized_rows', 0) or 0):,} harmonized rows")
    _print_run_warnings(metrics, metrics_path=output_dir / "slice_metrics.json")
    harmonized_path = output_dir / "slice_harmonized.csv"
    command = shell_command(
        [
            "degora",
            "run",
            config.resolve(),
            "--output-dir",
            output_dir.resolve(),
            "--harmonized-dir",
            harmonized_dir.resolve(),
            "--db",
            db_path.resolve(),
            "--min-studies",
            min_studies,
        ]
    )
    progress.start("Scoring genes and building the database")
    summary = write_score_database(
        harmonized_path,
        output_dir,
        catalog_path=config,
        db_path=db_path,
        min_studies=min_studies,
        command=command,
        extra_metadata={
            key: settings[key]
            for key in ("demo_search_keyword", "demo_search_species")
            if settings.get(key)
        },
    )
    progress.done(f"{int(summary.get('n_gene_scores', 0) or 0):,} genes scored")
    _print_run_warnings(
        {"warnings": summary.get("significance_warnings", [])},
        metrics_path=output_dir / "degora_score_db_summary.json",
    )
    if int(summary.get("n_gene_scores", 0) or 0) == 0:
        raise CliUsageError(
            f"DEGORA scored zero genes at min_studies={min_studies}. No gene had usable, "
            "directional evidence from enough independent source units."
            + _zero_gene_diagnostic(harmonized_path, min_studies)
        )
    workbook_path = output_dir / DEFAULT_WORKBOOK_NAME
    workbook_summary: dict[str, Any] | None = None
    if not getattr(args, "no_excel", False):
        progress.start(
            f"Writing the audit workbook for {int(summary.get('n_gene_scores', 0) or 0):,} genes"
            " (use --no-excel to skip)"
        )
        try:
            workbook_summary = export_run_workbook(
                output_dir,
                workbook_path,
                config_path=config,
                db_path=db_path,
                command=command,
            )
        except Exception as exc:  # noqa: BLE001 - convert to the stable user-facing CLI error contract
            raise CliUsageError(
                f"Excel workbook export failed ({type(exc).__name__}: {exc}). CSV and SQLite "
                "artifacts were written, but the default run contract is incomplete. Fix the "
                "export error and rerun, or explicitly use --no-excel if a workbook is not needed."
            ) from exc
        workbook_size = workbook_path.stat().st_size if workbook_path.exists() else 0
        progress.done(f"{workbook_size / (1024 * 1024):.1f} MB" if workbook_size else "")

    print("")
    print("DEGORA run complete")
    print(f"- DEGORA version: {format_version_info(version_info)}")
    print(f"- Results folder: {output_dir.resolve()}")
    print(f"- Score table: {Path(summary['score_csv']).resolve()}")
    print(f"- Database: {Path(summary['db_path']).resolve()}")
    if workbook_summary is not None:
        print(f"- Excel workbook: {Path(workbook_summary['output']).resolve()}")
    print(f"- Top genes: {', '.join(summary['top_genes'][:10])}")
    print("")
    print(f"Open browser/API with: degora serve {Path(summary['db_path']).resolve()}")

    if serve_after:
        port = args.port if args.port is not None else _int_setting(settings.get("browser_port"), 8765)
        serve_db(
            Path(summary["db_path"]),
            host=args.host,
            port=port,
            allow_network=args.allow_network,
            access_token=args.token,
        )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="degora",
        description="DEGORA: Excel-first integration of published DEG tables.",
    )
    parser.add_argument("--version", action="version", version=f"DEGORA {format_version_info()}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    template = subparsers.add_parser("template", help="Create an easy Excel config template.")
    template.add_argument("output", nargs="?", default="DEGORA_template.xlsx")
    template.add_argument("--force", action="store_true", help="Overwrite the template if it already exists.")

    init = subparsers.add_parser(
        "init",
        help="Build a config by answering questions about your DEG tables.",
        description=(
            "Read a folder of DEG result tables, work out the column mapping and table scope "
            "from the files themselves, and ask only what a file cannot say. The contrast "
            "direction is always asked and never guessed: reversing it inverts every up/down "
            "call in the results while leaving them looking entirely reasonable."
        ),
    )
    init.add_argument("output", nargs="?", default="degora_config.csv", help="Config file to write.")
    init.add_argument("--deg-dir", default=".", help="Folder holding the DEG tables (default: current folder).")
    init.add_argument("--force", action="store_true", help="Overwrite the config if it already exists.")

    demo = subparsers.add_parser("demo", help="Create a runnable demo workspace with tiny DEG tables.")
    demo.add_argument("output", nargs="?", default="degora-demo")
    demo.add_argument("--force", action="store_true", help="Overwrite demo files if the folder already exists.")
    demo.add_argument(
        "--keyword",
        help='Suggested live GEO search keyword (default: "hypoxia normoxia renal epithelial").',
    )
    demo.add_argument("--species", choices=["human", "mouse"], help="Suggested live GEO search tab (default: human).")

    validate = subparsers.add_parser("validate", help="Check an Excel/CSV config before running analysis.")
    validate.add_argument("config")

    run = subparsers.add_parser("run", help="Run harmonization and build the local score database.")
    run.add_argument("config")
    run.add_argument("--output-dir")
    run.add_argument("--harmonized-dir")
    run.add_argument("--db")
    run.add_argument("--min-studies", type=int)
    run.add_argument("--no-excel", action="store_true", help="Skip the default DEGORA_output.xlsx audit workbook.")

    launch = subparsers.add_parser("launch", help="Run analysis, then optionally start the browser.")
    launch.add_argument("config")
    launch.add_argument("--output-dir")
    launch.add_argument("--harmonized-dir")
    launch.add_argument("--db")
    launch.add_argument("--min-studies", type=int)
    launch.add_argument("--no-excel", action="store_true", help="Skip the default DEGORA_output.xlsx audit workbook.")
    launch.add_argument("--serve", action="store_true", help="Start the browser/API after the run finishes.")
    launch.add_argument("--host", default="127.0.0.1")
    launch.add_argument("--port", type=int)
    launch.add_argument("--allow-network", action="store_true", help="Allow non-loopback browser/API binding with token protection.")
    launch.add_argument("--token", help="Access token for non-loopback browser/API binding; generated when omitted.")

    serve = subparsers.add_parser("serve", help="Open the local browser/API for an existing score DB.")
    serve.add_argument("db")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8765)
    serve.add_argument("--allow-network", action="store_true", help="Allow non-loopback binding with token protection.")
    serve.add_argument("--token", help="Access token for non-loopback binding; generated when omitted.")

    discover = subparsers.add_parser(
        "discover",
        help="Search PubMed plus linked public data, or use explicit legacy GEO mode.",
        description=(
            "Search Human or Mouse publication/source-unit records. Default federated mode searches PubMed "
            "plus linked data, evaluates at most 1000 records, writes a full audit snapshot, and displays "
            f"{DISCOVERY_PAGE_SIZE} rows per globally sorted page. Human and Mouse runs stay separate; searching both species "
            "does not pool evidence. NCBI_EMAIL and NCBI_API_KEY are optional environment variables used "
            "for NCBI requests when configured."
        ),
    )
    discover.add_argument("query", help="Natural-language biological topic or perturbation.")
    discover.add_argument("--species", required=True, choices=["human", "mouse"])
    discover.add_argument(
        "--source",
        choices=["federated", "geo"],
        default="federated",
        help="Search backend: federated PubMed+linked-data workflow (default) or legacy GEO-only mode.",
    )
    discover.add_argument(
        "--limit",
        type=_bounded_discovery_limit,
        default=DISCOVERY_MAX_LIMIT,
        help="Maximum records to evaluate before global sorting; 1..1000 (default: 1000).",
    )
    discover.add_argument(
        "--page",
        type=int,
        default=1,
        help=f"One-based globally sorted display page; {DISCOVERY_PAGE_SIZE} rows per page.",
    )
    discover.add_argument("--output-dir", required=True)
    discover.add_argument(
        "--select",
        action="append",
        default=[],
        metavar="ID",
        help=(
            "Prepare this ID; repeat for multiple selections. Federated accepts PMID:, DOI:, GSE, "
            "canonical_id, or source_unit_id. Legacy GEO mode accepts GSE accessions only."
        ),
    )
    discover.add_argument("--inspection-budget", type=int, default=40)
    discover.add_argument("--force", action="store_true")

    discovery_analyze = subparsers.add_parser(
        "discovery-analyze",
        help="Reproduce a reviewed species-specific discovery run from JSON artifacts.",
    )
    discovery_analyze.add_argument("bundle_json")
    discovery_analyze.add_argument("selections_json")
    discovery_analyze.add_argument("--species", required=True, choices=["human", "mouse"])
    discovery_analyze.add_argument("--output-dir", required=True)
    discovery_analyze.add_argument("--min-studies", type=int, default=2)
    discovery_analyze.add_argument("--force", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "template":
            from .excel_template import write_template

            path = write_template(args.output, force=args.force)
            print(f"Wrote Excel template: {path}")
            print("Next: edit the Contrasts sheet, then run: degora validate <your_config.xlsx>")
            return 0
        if args.command == "init":
            from .beginner import BeginnerInitError, run_init

            try:
                summary = run_init(args.output, args.deg_dir, force=args.force)
            except BeginnerInitError as exc:
                raise CliUsageError(str(exc)) from exc
            print("")
            print(f"Wrote config: {summary['config_path']}")
            print(f"- Contrasts: {summary['n_contrasts']}")
            print(f"- Independent source units: {summary['n_source_units']}")
            if summary["n_excluded_reversed_direction"]:
                print(
                    f"- Excluded (positive means up in the control group): "
                    f"{summary['n_excluded_reversed_direction']} - see the notes column"
                )
            for entry in summary["skipped"]:
                print(f"- Skipped {entry['path']}: {entry['reason']}")
            print("")
            print(f"Next: degora validate {summary['config_path']}")
            return 0
        if args.command == "demo":
            from .demo import write_demo_workspace

            demo = write_demo_workspace(
                args.output,
                force=args.force,
                search_keyword=args.keyword,
                search_species=args.species,
            )
            print(f"Wrote DEGORA demo: {demo['demo_dir']}")
            print("")
            print("Try it:")
            print(f"  cd {demo['demo_dir']}")
            print("  degora validate degora_demo_config.xlsx")
            print("  degora run degora_demo_config.xlsx")
            print("  degora serve results/degora_scores.db")
            print("")
            print(
                f"Live discovery is pre-filled for {demo['search_species']}: "
                f"{demo['search_keyword']} (press Search papers & data yourself)."
            )
            return 0
        if args.command == "validate":
            from .slice_runner import validate_catalog_inputs

            config = Path(args.config)
            if not config.exists():
                validate_catalog_inputs(config)
            settings = read_excel_settings(config)
            _validate_score_version(settings)
            progress = _RunProgress()

            def validation_progress(done: int, total: int, study_id: str) -> None:
                if total <= 0:
                    return
                step = max(1, min(1000, total // 100 or 1))
                if done in {0, 1, total} or done % step == 0:
                    suffix = f" ({study_id})" if study_id else ""
                    progress.update(f"checked {done:,}/{total:,} active row(s){suffix}", force=done in {0, total})

            progress.start("Validating the catalog and source tables")
            validation = validate_catalog_inputs(config, progress=validation_progress)
            progress.done(f"{validation.get('active_contrasts', '?')} contrast(s)")
            print("DEGORA config OK")
            _print_validation_summary(validation, include_excluded=True)
            warnings = [str(message).strip() for message in validation.get("warnings", []) if str(message).strip()]
            # A config with fewer independent source units than the replication rule
            # requires cannot score a single gene. Saying "OK" and letting the run
            # find that out spends the run to deliver a fact already visible here.
            configured_min_studies = _int_setting(settings.get("min_studies"), 2)
            source_units = int(validation.get("source_units") or 0)
            if source_units < configured_min_studies:
                warnings.insert(
                    0,
                    f"This config has {source_units} independent source unit(s) but min_studies is "
                    f"{configured_min_studies}, so a run scores zero genes. Give each independent "
                    f"study its own source_unit_id, add another study, or lower Project.min_studies.",
                )
            if warnings:
                print("", file=sys.stderr)
                print("Non-fatal input warnings:", file=sys.stderr)
                for message in warnings[:WARNING_DISPLAY_LIMIT]:
                    print(f"- {message}", file=sys.stderr)
                if len(warnings) > WARNING_DISPLAY_LIMIT:
                    print(f"- ... {len(warnings) - WARNING_DISPLAY_LIMIT} more warning(s)", file=sys.stderr)
            return 0
        if args.command == "run":
            return _run_from_config(args)
        if args.command == "launch":
            return _run_from_config(args, serve_after=args.serve)
        if args.command == "serve":
            from .api import serve as serve_db

            serve_db(Path(args.db), host=args.host, port=args.port, allow_network=args.allow_network, access_token=args.token)
            return 0
        if args.command == "discover":
            output = Path(args.output_dir).resolve()
            if output.exists() and any(output.iterdir()) and not args.force:
                raise FileExistsError(f"discovery output already exists and is not empty: {output}")
            if args.source == "geo":
                from .discovery import export_search_page, prepare_geo_studies, search_geo

                if args.select and any(not str(value).strip().upper().startswith("GSE") for value in args.select):
                    raise CliUsageError("legacy GEO mode accepts only GSE accessions in --select")
                if args.select:
                    result = prepare_geo_studies(
                        args.select,
                        args.species,
                        query=args.query,
                        inspection_budget=args.inspection_budget,
                        materialize_dir=output,
                        force=args.force,
                    )
                    print(f"Prepared {result['returned_studies']} {args.species} GEO studies: {output}")
                    print(f"Draft catalog (inactive): {result['exports']['draft_catalog_csv']}")
                    print("Review contrast direction and table scope before activation; Human and Mouse remain separate.")
                else:
                    result = search_geo(
                        args.query,
                        args.species,
                        page=args.page,
                        page_size=DISCOVERY_PAGE_SIZE,
                        assess_files=True,
                        global_rank=True,
                        global_limit=args.limit,
                    )
                    exports = export_search_page(result, output, force=args.force)
                    print(
                        f"Found up to {result['total_hits']} pre-filter NCBI hits; "
                        f"wrote legacy GEO globally ranked exact {args.species} page {result['page']} "
                        f"from {result['evaluated_studies']} assessed studies to {exports['search_csv']}"
                    )
                return 0

            from .discovery_export import export_publication_search
            from .discovery_federated import page_publication_snapshot, resolve_publication_records, search_publications

            snapshot = search_publications(args.query, args.species, limit=args.limit)
            records = list(snapshot.get("records", []))
            if args.select:
                try:
                    from .discovery_prepare import prepare_publication_records
                except ModuleNotFoundError as exc:
                    if exc.name and exc.name.endswith("discovery_prepare"):
                        raise CliUsageError(
                            "federated publication preparation backend is unavailable; "
                            "run without --select to export the review snapshot, or use --source geo for legacy GSE preparation"
                        ) from exc
                    raise
                selected = _select_publication_records(records, args.select)
                selected = resolve_publication_records(selected, args.species)
                files_per_record = max(1, min(12, args.inspection_budget // max(1, len(selected))))
                result = prepare_publication_records(
                    selected,
                    args.species,
                    query=args.query,
                    max_files_per_record=files_per_record,
                    materialize_dir=output,
                    force=args.force,
                )
                prepared_count = result.get("returned_records", result.get("returned_studies", len(selected)))
                exports = result.get("exports", {})
                # "Prepared 8 record(s)" alone reads as success even when every one
                # of them resolved to an upstream matrix or no usable table, and the
                # draft catalog is a header with no rows under it.
                ready = sum(
                    int(study.get("ready_for_review_count", 0) or 0)
                    for study in result.get("studies", [])
                    if isinstance(study, dict)
                )
                print(f"Prepared {prepared_count} {args.species} publication/source-unit record(s): {output}")
                print(
                    f"Ready for review: {ready} table(s) across those records"
                    + ("" if ready else " - nothing can be activated yet; see discovery_audit.json for the per-file reason")
                )
                if exports.get("draft_catalog_csv"):
                    print(f"Draft catalog (inactive): {exports['draft_catalog_csv']}")
                print(
                    "Review required before analysis: confirm species, source-unit independence, table scope, "
                    "contrast direction, sample groups, and source provenance. Human and Mouse remain separate."
                )
            else:
                exports = export_publication_search(snapshot, output, force=args.force)
                page = page_publication_snapshot(
                    snapshot,
                    page=args.page,
                    page_size=DISCOVERY_PAGE_SIZE,
                    sort_by="data_readiness",
                    sort_order=None,
                )
                print(
                    f"Wrote federated {args.species} PubMed+linked-data snapshot with "
                    f"{len(records)} globally ranked record(s) evaluated under the {args.limit}-record cap."
                )
                print(
                    f"Display page {page['page']} contains {len(page.get('records', []))} row(s); "
                    f"page size is {DISCOVERY_PAGE_SIZE}."
                )
                _print_publication_page(
                    list(page.get("records", [])),
                    page=int(page["page"]),
                    page_size=DISCOVERY_PAGE_SIZE,
                )
                print(f"Search CSV: {exports['search_csv']}")
                print("Review Human and Mouse snapshots separately; cross-species pooling is not performed.")
            return 0
        if args.command == "discovery-analyze":
            from .aggregate import validate_min_studies
            from .discovery_run import run_discovery_analysis

            bundle_path = Path(args.bundle_json)
            selections_path = Path(args.selections_json)
            try:
                prepared = json.loads(bundle_path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError) as exc:
                raise CliUsageError(f"bundle_json is not readable UTF-8 JSON: {bundle_path}: {exc}") from exc
            try:
                selection_payload = json.loads(selections_path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError) as exc:
                raise CliUsageError(
                    f"selections_json is not readable UTF-8 JSON: {selections_path}: {exc}"
                ) from exc
            selections = selection_payload.get("selections", selection_payload) if isinstance(selection_payload, dict) else selection_payload
            if not isinstance(prepared, dict) or not isinstance(selections, list):
                raise CliUsageError(
                    "bundle_json must contain an object and selections_json must contain a list or {selections: [...]}"
                )
            try:
                min_studies = validate_min_studies(args.min_studies)
            except ValueError as exc:
                raise CliUsageError(str(exc)) from exc
            result = run_discovery_analysis(
                prepared,
                selections,
                args.output_dir,
                species=args.species,
                min_studies=min_studies,
                force=args.force,
            )
            print(f"DEGORA {args.species} discovery run complete: {result['db_path']}")
            print(f"Top genes: {', '.join(result['top_genes'][:10])}")
            return 0
    except (FileExistsError, FileNotFoundError, PermissionError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except Exception as exc:
        if exc.__class__.__name__ == "DegoraConfigError" and exc.__class__.__module__.endswith(".slice_runner"):
            print(str(exc), file=sys.stderr)
            return 2
        if isinstance(exc, CliUsageError) or exc.__class__.__name__ in {
            "DiscoveryUnavailableError",
            "DiscoveryWorkspaceInUseError",
            "OutputDirectoryBusyError",
        } or (
            exc.__class__.__module__.endswith((".discovery", ".discovery_run", ".reanalysis"))
            and isinstance(exc, ValueError)
        ):
            print(str(exc), file=sys.stderr)
            return 2
        raise
    raise AssertionError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
