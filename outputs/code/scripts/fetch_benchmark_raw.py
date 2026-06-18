#!/usr/bin/env python
"""Download the raw count matrices and annotations that the RNA-seq derivation
scripts consume, so the derivation is reproducible from scratch without any
manual download.

write_ifn_derived_deg.py, write_er_stress_benchmark.py, and
write_heat_shock_benchmark.py each hold a deterministic logCPM-Welch derivation
plus a RAW_INPUTS table of {url, path}, but they expect the raw inputs to already
be present and raise FileNotFoundError otherwise. This fetcher reads those same
RAW_INPUTS tables and downloads each missing file to its declared path (atomic,
validated, never overwriting an existing file), closing the one gap between the
recorded source URLs and a fully agent-free reproduction.

    PYTHONPATH=outputs/code python outputs/code/scripts/fetch_benchmark_raw.py --topic heat
    make -C outputs/code fetch-rnaseq            # fetch raw for all three RNA-seq topics
"""

from __future__ import annotations

import argparse
import importlib
import urllib.request
from pathlib import Path
from typing import Any

# topic -> derivation module that declares RAW_INPUTS = {key: {"url", "path", ...}}
SCRIPTS = {
    "ifn": "scripts.write_ifn_derived_deg",
    "er": "scripts.write_er_stress_benchmark",
    "heat": "scripts.write_heat_shock_benchmark",
}
DOWNLOAD_TIMEOUT_SECONDS = 300
MAX_BYTES = 3 * 1024 * 1024 * 1024  # 3 GiB cap (NCBI annotation files are large)


def _resolve_url(url: Any) -> str:
    """RAW_INPUTS urls are sometimes written as a multi-line string tuple."""

    if isinstance(url, (tuple, list)):
        return "".join(url)
    return str(url or "")


def _download(url: str, path: Path) -> tuple[bool, str]:
    tmp = path.with_name(path.name + ".part")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(url, headers={"User-Agent": "DEGORA-fetch/1.0"})
        total = 0
        with urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response, open(tmp, "wb") as handle:
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_BYTES:
                    tmp.unlink(missing_ok=True)
                    return False, f"exceeds {MAX_BYTES // (1024 ** 3)} GiB cap"
                handle.write(chunk)
        if total == 0:
            tmp.unlink(missing_ok=True)
            return False, "empty response"
        tmp.replace(path)
        return True, f"{total:,} bytes"
    except Exception as exc:  # noqa: BLE001 - report any fetch failure, never crash the run
        tmp.unlink(missing_ok=True)
        return False, f"{type(exc).__name__}: {exc}"


def fetch_topic(topic: str) -> dict[str, int]:
    module = importlib.import_module(SCRIPTS[topic])
    raw_inputs: dict[str, dict[str, Any]] = getattr(module, "RAW_INPUTS", {})
    counts = {"present": 0, "downloaded": 0, "failed": 0}
    for entry in raw_inputs.values():
        path = Path(str(entry["path"]))
        url = _resolve_url(entry.get("url", ""))
        if path.exists():
            counts["present"] += 1
            print(f"  present   {path.name}")
            continue
        if not url.lower().startswith(("http://", "https://")):
            counts["failed"] += 1
            print(f"  NO-URL    {path.name}")
            continue
        print(f"  fetching  {path.name} ...", flush=True)
        ok, note = _download(url, path)
        counts["downloaded" if ok else "failed"] += 1
        print(f"    {'OK' if ok else 'FAIL'}: {note}")
    return counts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--topic", choices=[*SCRIPTS, "all"], default="all")
    args = parser.parse_args(argv)
    topics = list(SCRIPTS) if args.topic == "all" else [args.topic]
    totals = {"present": 0, "downloaded": 0, "failed": 0}
    for topic in topics:
        print(f"[{topic}]")
        for key, value in fetch_topic(topic).items():
            totals[key] += value
    print(f"\nraw inputs: present={totals['present']} downloaded={totals['downloaded']} failed={totals['failed']}")
    return 0 if totals["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
