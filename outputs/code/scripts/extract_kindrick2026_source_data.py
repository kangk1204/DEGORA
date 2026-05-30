#!/usr/bin/env python
"""Extract Kindrick 2026 hypoxia RNA-seq source-data tables."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


TABLES = {
    "pc3": {
        "sheet": "Figure 4A&C source data ",
        "output": "GSE296192_PC3_hypoxia_vs_normoxia.csv",
        "expected_rows_min": 10000,
    },
    "hct116": {
        "sheet": "Figure 4B&D source data",
        "output": "GSE296192_HCT116_hypoxia_vs_normoxia.csv",
        "expected_rows_min": 10000,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_table(workbook: Path, sheet: str) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name=sheet, header=None)
    header = raw.iloc[1].tolist()
    first_columns = [idx for idx, value in enumerate(header) if pd.notna(value)]
    if not first_columns:
        raise ValueError(f"No header row detected in sheet {sheet!r}")
    # The source-data sheets place the gene identifier in the unlabeled first column,
    # followed by replicate values, foldchange, and Adjusted_P_Value.
    last_required = max(i for i, value in enumerate(header) if value in {"foldchange", "Adjusted_P_Value"})
    frame = raw.iloc[2:, : last_required + 1].copy()
    columns = header[: last_required + 1]
    columns[0] = "gene_symbol"
    frame.columns = columns
    frame = frame.loc[frame["gene_symbol"].notna()].copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype("string").str.strip()
    frame["hypoxia_log2FoldChange"] = pd.to_numeric(frame["foldchange"], errors="coerce")
    frame["pvalue"] = pd.to_numeric(frame["Adjusted_P_Value"], errors="coerce")
    frame["padj"] = pd.to_numeric(frame["Adjusted_P_Value"], errors="coerce")
    frame = frame.loc[
        frame["gene_symbol"].notna()
        & frame["hypoxia_log2FoldChange"].notna()
        & frame["pvalue"].notna()
    ].copy()
    return frame


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    workbook = args.workbook.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_hash = sha256(workbook)

    for key, spec in TABLES.items():
        frame = extract_table(workbook, str(spec["sheet"]))
        if len(frame) < int(spec["expected_rows_min"]):
            raise ValueError(f"{key} extracted only {len(frame)} rows")
        output = output_dir / str(spec["output"])
        frame.to_csv(output, index=False)
        command = (
            "PYTHONPATH=outputs/code python outputs/code/scripts/extract_kindrick2026_source_data.py "
            f"--workbook {workbook} --output-dir {output_dir}"
        )
        output.with_suffix(output.suffix + ".source").write_text(
            command
            + "\n"
            + f"raw_workbook_sha256={raw_hash}\n"
            + f"sheet={spec['sheet']}\n"
        )
        print(f"{output}: rows={len(frame)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
