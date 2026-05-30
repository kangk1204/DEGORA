#!/usr/bin/env python
"""Serve the local DEGORA score database through JSON API and browser UI."""

from __future__ import annotations

import argparse
from pathlib import Path

from degora.api import serve


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    serve(args.db, host=args.host, port=args.port, quiet=args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
