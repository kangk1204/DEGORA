#!/usr/bin/env python
"""Generate CORE-EXPANDED (primary) and FULL-EXPANDED (supplementary) gold panels.

Each topic's primary locked CSV is OVERWRITTEN with the core-expanded canonical
gene list, and a sibling ``*_full.csv`` supplementary CSV is created with the
full-expanded list (core + extended-only). The existing column schema of every
primary CSV is preserved exactly; the supplementary CSV reuses the matching
primary schema. All rows have ``expected_direction = up``; the per-topic tier
column (``panel_role`` or ``evidence_class``) is set to ``core_canonical`` for
core genes and ``extended_canonical`` for extended-only genes.

These panels are EVALUATION-ONLY (recall / prioritization diagnostics). They are
locked blind to DEGORA output and never feed the scoring pipeline.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
GOLD_DIR = PROJECT_ROOT / "data/studies/gold"

EVIDENCE_BASIS = (
    "Expanded canonical panel (MSigDB Hallmark + topic-specific target/ChIP "
    "literature); locked blind to DEGORA output."
)

CORE_TIER = "core_canonical"
EXTENDED_TIER = "extended_canonical"

# --- Gene lists (order preserved as provided) -------------------------------

IFN_CORE = [
    "ISG15", "RSAD2", "CMPK2", "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3", "OAS1",
    "OAS2", "OAS3", "OASL", "IFI44", "IFI44L", "IFI6", "IFITM1", "USP18",
    "DDX58", "IFIH1", "HERC5", "STAT1", "IRF7", "IFI27", "BST2", "ISG20",
    "IRF9", "STAT2", "EIF2AK2",
]
IFN_EXTENDED = [
    "DDX60", "HERC6", "IFI35", "IFITM3", "LY6E", "RTP4", "PARP9", "PARP12",
    "DTX3L", "LGALS3BP", "EPSTI1", "SAMD9", "SAMD9L", "PLSCR1", "XAF1",
    "TRIM22", "SP100", "ADAR", "SAMHD1", "UBE2L6", "UBA7", "TRIM5",
]

ER_CORE = [
    "HSPA5", "DDIT3", "ATF4", "ATF3", "HERPUD1", "TRIB3", "DDIT4", "ASNS",
    "CHAC1", "PPP1R15A", "ATF6", "XBP1", "DNAJB9", "DNAJB11", "DNAJC3",
    "PDIA4", "PDIA6", "HSP90B1", "HYOU1", "MANF", "CRELD2", "SDF2L1", "EDEM1",
    "SEL1L", "SERP1", "WFS1", "ERO1A", "STC2", "CEBPG", "PSAT1", "MTHFD2",
    "DERL3",
]
ER_EXTENDED = [
    "ERN1", "EIF2AK3", "PDIA5", "ERO1B", "VEGFA", "CEBPB", "PHGDH", "SLC7A11",
    "SLC7A5", "SLC1A4", "SLC3A2", "CALR", "CANX", "WIPI1",
]

HEAT_CORE = [
    "HSPA1A", "HSPA1B", "HSPA1L", "HSPA6", "HSPA4L", "HSPH1", "HSPD1", "HSPE1",
    "HSP90AA1", "HSP90AB1", "HSPB1", "HSPB8", "SERPINH1", "DNAJB1", "DNAJA1",
    "DNAJB4", "DNAJB6", "DNAJA4", "BAG3", "BAG2", "STIP1", "FKBP4", "CHORDC1",
    "AHSA1", "HSPA4", "DNAJC7",
]
HEAT_EXTENDED = [
    "HSPA8", "PTGES3", "CRYAB", "BAG1", "UBB", "UBC", "ZFAND2A", "CDC37",
    "HSPA2", "HSPB2", "DNAJA2", "HDAC6", "DEDD2", "HSF1",
]

HYPOXIA_CORE = [
    "VEGFA", "SLC2A1", "SLC2A3", "HK2", "PGK1", "LDHA", "ALDOA", "ALDOC",
    "ENO1", "PFKFB3", "PFKFB4", "PDK1", "CA9", "ADM", "EPO", "ANGPTL4",
    "EGLN3", "EGLN1", "NDRG1", "BNIP3", "BNIP3L", "DDIT4", "P4HA1", "P4HA2",
    "PLOD2", "KDM3A", "KDM4B", "MXI1", "BHLHE40", "FLT1", "HK1", "PFKL",
    "ANGPT2",
]
HYPOXIA_EXTENDED = [
    "ENO2", "LOX", "VLDLR", "AK4", "ANKRD37", "TMEM45A", "PGM1", "GAPDH",
    "TPI1", "STC2", "GBE1", "ADORA2B", "SERPINE1", "INSIG2",
]


def _note(topic: str, tier: str) -> str:
    role = "core canonical" if tier == CORE_TIER else "extended canonical"
    return (
        f"{topic} {role} target; expanded gold panel locked blind to DEGORA "
        "output; used only for recall/prioritization diagnostics."
    )


# topic name, primary filename, full filename, column schema, tier column,
# lock column + value, core list, extended list.
PANELS = [
    {
        "topic": "Type I/II interferon (ISG)",
        "primary": "ifn_gold_panel.csv",
        "full": "ifn_gold_panel_full.csv",
        "columns": ["gene_symbol", "expected_direction", "panel_role", "locked_at", "evidence_basis", "notes"],
        "tier_column": "panel_role",
        "lock_column": "locked_at",
        "lock_value": "2026-05-30",
        "core": IFN_CORE,
        "extended": IFN_EXTENDED,
    },
    {
        "topic": "ER stress / UPR",
        "primary": "er_stress_upr_gold_panel.csv",
        "full": "er_stress_upr_gold_panel_full.csv",
        "columns": ["gene_symbol", "expected_direction", "evidence_class", "locked", "evidence_basis", "notes"],
        "tier_column": "evidence_class",
        "lock_column": "locked",
        "lock_value": "yes",
        "core": ER_CORE,
        "extended": ER_EXTENDED,
    },
    {
        "topic": "Heat shock / HSF1",
        "primary": "heat_shock_hsf1_gold_panel.csv",
        "full": "heat_shock_hsf1_gold_panel_full.csv",
        "columns": ["gene_symbol", "expected_direction", "evidence_class", "locked", "evidence_basis", "notes"],
        "tier_column": "evidence_class",
        "lock_column": "locked",
        "lock_value": "yes",
        "core": HEAT_CORE,
        "extended": HEAT_EXTENDED,
    },
    {
        "topic": "Hypoxia / HIF",
        "primary": "hypoxia_hif1_gold_panel.csv",
        "full": "hypoxia_hif1_gold_panel_full.csv",
        "columns": ["gene_symbol", "panel_role", "expected_direction", "locked_at_utc", "evidence_basis", "notes"],
        "tier_column": "panel_role",
        "lock_column": "locked_at_utc",
        "lock_value": "2026-05-30T00:00:00Z",
        "core": HYPOXIA_CORE,
        "extended": HYPOXIA_EXTENDED,
    },
]


def _build_frame(spec: dict, genes_with_tier: list[tuple[str, str]]) -> pd.DataFrame:
    records = []
    for gene, tier in genes_with_tier:
        record = {
            "gene_symbol": gene,
            "expected_direction": "up",
            spec["tier_column"]: tier,
            spec["lock_column"]: spec["lock_value"],
            "evidence_basis": EVIDENCE_BASIS,
            "notes": _note(spec["topic"], tier),
        }
        records.append(record)
    frame = pd.DataFrame.from_records(records)
    return frame[spec["columns"]]


def _write_sidecars(csv_path: Path, command: str, generator: str) -> None:
    source_path = Path(str(csv_path) + ".source")
    source_path.write_text(command + "\n")
    data = csv_path.read_bytes()
    rel = csv_path.relative_to(PROJECT_ROOT).as_posix()
    prov = {
        "artifact_path": rel,
        "artifact_sha256": hashlib.sha256(data).hexdigest(),
        "artifact_size_bytes": len(data),
        "command": command,
        "inputs": [],
        "metadata": {"generator": generator},
    }
    prov_path = Path(str(csv_path) + ".provenance.json")
    prov_path.write_text(json.dumps(prov, indent=2) + "\n")


def main() -> None:
    command = "make -C outputs/code expanded-gold-panels"
    for spec in PANELS:
        core_set = list(dict.fromkeys(spec["core"]))
        extended_set = list(dict.fromkeys(spec["extended"]))

        # Primary: CORE-EXPANDED (all core_canonical).
        primary_frame = _build_frame(spec, [(g, CORE_TIER) for g in core_set])
        primary_path = GOLD_DIR / spec["primary"]
        primary_frame.to_csv(primary_path, index=False)
        _write_sidecars(primary_path, command, f"expanded-core-gold-panel:{spec['primary']}")

        # Supplementary: FULL-EXPANDED (core_canonical + extended_canonical).
        full_rows = [(g, CORE_TIER) for g in core_set] + [(g, EXTENDED_TIER) for g in extended_set]
        full_frame = _build_frame(spec, full_rows)
        full_path = GOLD_DIR / spec["full"]
        full_frame.to_csv(full_path, index=False)
        _write_sidecars(full_path, command, f"expanded-full-gold-panel:{spec['full']}")

        print(
            f"{spec['topic']}: core={len(core_set)} "
            f"full={len(core_set) + len(extended_set)} "
            f"-> {spec['primary']} / {spec['full']}"
        )


if __name__ == "__main__":
    main()
