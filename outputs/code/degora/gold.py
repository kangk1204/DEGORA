"""Locked validation gene panels from the DEGORA methodology specification.

HIF1A_UP_TARGETS holds the CORE-EXPANDED hypoxia/HIF canonical panel that
matches data/studies/gold/hypoxia_hif1_gold_panel.csv. It is evaluation-only
(recall / prioritization diagnostics) and is locked blind to DEGORA output.
"""

HIF1A_UP_TARGETS = frozenset(
    {
        "VEGFA",
        "SLC2A1",
        "SLC2A3",
        "HK2",
        "PGK1",
        "LDHA",
        "ALDOA",
        "ALDOC",
        "ENO1",
        "PFKFB3",
        "PFKFB4",
        "PDK1",
        "CA9",
        "ADM",
        "EPO",
        "ANGPTL4",
        "EGLN3",
        "EGLN1",
        "NDRG1",
        "BNIP3",
        "BNIP3L",
        "DDIT4",
        "P4HA1",
        "P4HA2",
        "PLOD2",
        "KDM3A",
        "KDM4B",
        "MXI1",
        "BHLHE40",
        "FLT1",
        "HK1",
        "PFKL",
        "ANGPT2",
    }
)

NEGATIVE_CONTROL_TARGETS = frozenset({"RPL13A", "HPRT1", "TBP"})
