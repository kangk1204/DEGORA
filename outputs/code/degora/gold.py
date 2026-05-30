"""Locked validation gene panels from METHODOLOGY_SPEC.md section 3.2."""

HIF1A_UP_TARGETS = frozenset(
    {
        "SLC2A1",
        "HK1",
        "HK2",
        "PFKL",
        "ALDOA",
        "ENO1",
        "PGK1",
        "LDHA",
        "PDK1",
        "VEGFA",
        "ANGPT2",
        "FLT1",
        "EPO",
        "CA9",
        "SLC2A3",
        "BNIP3",
        "BNIP3L",
        "DDIT4",
        "EGLN1",
        "EGLN3",
    }
)

NEGATIVE_CONTROL_TARGETS = frozenset({"RPL13A", "HPRT1", "TBP"})
