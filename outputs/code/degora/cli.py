"""Beginner-friendly command line interface for DEGORA."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from . import format_version_info, runtime_version_info


DEFAULT_OUTPUT_DIR = Path("outputs/results/degora-run")
DEFAULT_HARMONIZED_DIR = Path("data/deg/harmonized")
WARNING_DISPLAY_LIMIT = 8


def _setting_key(value: Any) -> str:
    return str(value).strip().lower().replace(" ", "_").replace("-", "_")


def _sheet_settings(path: Path, sheet: str, key_column: str) -> dict[str, str]:
    import pandas as pd

    from .excel_io import read_config_sheet

    try:
        frame = read_config_sheet(path, sheet)
    except ValueError:
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


def _int_setting(value: str | None, default: int) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except ValueError as exc:
        from .slice_runner import DegoraConfigError

        raise DegoraConfigError(
            "numeric setting is invalid",
            problems=[f"Expected a whole number, but got {value!r}."],
            fixes=["Open the Project or AdvancedSettings sheet and enter a number such as 2."],
        ) from exc


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
    for key in ("warnings", "rank_universe_warnings"):
        warning_values.extend(metrics.get(key, []) or [])

    try:
        clipped_rows = int(metrics.get("pvalue_clipped_rows", 0) or 0)
    except (TypeError, ValueError):
        clipped_rows = 0
    if clipped_rows:
        text = (
            f"{clipped_rows} row(s) reported pvalue == 0; values were floored "
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


def _run_from_config(args: argparse.Namespace, *, serve_after: bool = False) -> int:
    from .api import serve as serve_db
    from .excel_export import DEFAULT_WORKBOOK_NAME, export_run_workbook
    from .provenance import shell_command
    from .score_db import write_score_database
    from .slice_runner import run_slice, validate_catalog_inputs

    version_info = runtime_version_info()
    config = Path(args.config)
    config_base = config.resolve().parent
    settings = read_excel_settings(config)
    min_studies = args.min_studies if args.min_studies is not None else _int_setting(settings.get("min_studies"), 2)
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

    validation = validate_catalog_inputs(config)
    print("DEGORA config OK")
    _print_validation_summary(validation)

    metrics = run_slice(config, output_dir, harmonized_dir, min_studies=min_studies)
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
    summary = write_score_database(
        harmonized_path,
        output_dir,
        catalog_path=config,
        db_path=db_path,
        min_studies=min_studies,
        command=command,
    )
    workbook_path = output_dir / DEFAULT_WORKBOOK_NAME
    workbook_summary: dict[str, Any] | None = None
    if not getattr(args, "no_excel", False):
        workbook_summary = export_run_workbook(
            output_dir,
            workbook_path,
            config_path=config,
            db_path=db_path,
            command=command,
        )

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

    demo = subparsers.add_parser("demo", help="Create a runnable demo workspace with tiny DEG tables.")
    demo.add_argument("output", nargs="?", default="degora-demo")
    demo.add_argument("--force", action="store_true", help="Overwrite demo files if the folder already exists.")

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
        if args.command == "demo":
            from .demo import write_demo_workspace

            demo = write_demo_workspace(args.output, force=args.force)
            print(f"Wrote DEGORA demo: {demo['demo_dir']}")
            print("")
            print("Try it:")
            print(f"  cd {demo['demo_dir']}")
            print("  degora validate degora_demo_config.xlsx")
            print("  degora run degora_demo_config.xlsx")
            print("  degora serve results/degora_scores.db")
            return 0
        if args.command == "validate":
            from .slice_runner import validate_catalog_inputs

            validation = validate_catalog_inputs(Path(args.config))
            print("DEGORA config OK")
            _print_validation_summary(validation, include_excluded=True)
            warnings = [str(message).strip() for message in validation.get("warnings", []) if str(message).strip()]
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
    except (FileExistsError, FileNotFoundError, PermissionError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except Exception as exc:
        if exc.__class__.__name__ == "DegoraConfigError" and exc.__class__.__module__.endswith(".slice_runner"):
            print(str(exc), file=sys.stderr)
            return 2
        raise
    raise AssertionError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
