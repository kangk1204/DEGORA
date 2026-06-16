#!/usr/bin/env python
"""Export one DEGORA run folder to an Excel audit workbook."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from degora.excel_export import DEFAULT_WORKBOOK_NAME, export_run_workbook
from degora.provenance import shell_command


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--db", type=Path)
    args = parser.parse_args(argv)
    output = args.output or (args.result_dir / DEFAULT_WORKBOOK_NAME)
    command = shell_command(
        [
            "env",
            "PYTHONPATH=outputs/code",
            "python3",
            "outputs/code/scripts/export_degora_run_workbook.py",
            "--result-dir",
            args.result_dir.as_posix(),
            "--output",
            output.as_posix(),
            *(["--config", args.config.as_posix()] if args.config else []),
            *(["--db", args.db.as_posix()] if args.db else []),
        ]
    )
    print(
        json.dumps(
            export_run_workbook(
                args.result_dir,
                output,
                config_path=args.config,
                db_path=args.db,
                command=command,
            ),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
