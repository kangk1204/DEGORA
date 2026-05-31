#!/usr/bin/env python3
"""Auxiliary broad ground-truth (MSigDB Hallmark) sensitivity recompute.

Compares marker recovery for RNA-seq-only vs RNA-seq+microarray PIPER corpora
using (G1) the existing compact locked gold panels and (G2) large MSigDB Hallmark
gene sets, against the quality-weighted PIPER ranking.

Outputs:
  outputs/results/code-review/hallmark_sensitivity_recall.csv
  outputs/results/code-review/hallmark_sensitivity_recall.md
  outputs/results/code-review/hallmark_sensitivity_recall.csv.source
"""
import json
import os
import sys
import subprocess
import urllib.request

import pandas as pd

ROOT = "/home/keunsoo/projects/09_PIPER"
OUT = os.path.join(ROOT, "outputs/results/code-review")
os.makedirs(OUT, exist_ok=True)

try:
    from sklearn.metrics import roc_auc_score, average_precision_score
    HAVE_SK = True
except Exception:
    HAVE_SK = False


def manual_auroc(scores, labels):
    # Mann-Whitney U based AUROC
    import numpy as np
    s = np.asarray(scores, float)
    y = np.asarray(labels, int)
    order = s.argsort()
    ranks = np.empty(len(s), float)
    ranks[order] = np.arange(1, len(s) + 1)
    # average ties
    # simple tie handling
    df = pd.DataFrame({"s": s, "r": ranks})
    df["r"] = df.groupby("s")["r"].transform("mean")
    ranks = df["r"].values
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    sum_pos = ranks[y == 1].sum()
    return (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def manual_auprc(scores, labels):
    import numpy as np
    s = np.asarray(scores, float)
    y = np.asarray(labels, int)
    order = (-s).argsort()
    y = y[order]
    tp = 0
    fp = 0
    n_pos = y.sum()
    if n_pos == 0:
        return float("nan")
    prev_recall = 0.0
    ap = 0.0
    for i in range(len(y)):
        if y[i] == 1:
            tp += 1
        else:
            fp += 1
        precision = tp / (tp + fp)
        recall = tp / n_pos
        if y[i] == 1:
            ap += precision * (recall - prev_recall)
            prev_recall = recall
    return ap


def auroc(scores, labels):
    if HAVE_SK:
        return float(roc_auc_score(labels, scores))
    return manual_auroc(scores, labels)


def auprc(scores, labels):
    if HAVE_SK:
        return float(average_precision_score(labels, scores))
    return manual_auprc(scores, labels)


# ---------------------------------------------------------------------------
# 1) Obtain MSigDB Hallmark gene sets
# ---------------------------------------------------------------------------
def load_hallmark():
    source = None
    lib = None
    # (i) gseapy
    try:
        import gseapy
        lib = gseapy.get_library("MSigDB_Hallmark_2020")
        source = "gseapy.get_library('MSigDB_Hallmark_2020')"
        return lib, source
    except Exception as e:
        sys.stderr.write(f"[gseapy failed] {e}\n")
    # (ii) Enrichr GMT download
    gmt_path = "/tmp/hallmark.gmt"
    try:
        url = ("https://maayanlab.cloud/Enrichr/geneSetLibrary?"
               "mode=text&libraryName=MSigDB_Hallmark_2020")
        req = urllib.request.Request(url, headers={"User-Agent": "curl/7.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = r.read().decode("utf-8", "replace")
        if not data.strip():
            raise ValueError("Enrichr returned an empty body (no network / blocked)")
        with open(gmt_path, "w") as f:
            f.write(data)
        lib = parse_gmt(gmt_path)
        # Enrichr MSigDB_Hallmark_2020 uses friendly set names (e.g. "Hypoxia",
        # "Interferon Alpha Response") rather than the HALLMARK_* prefix; accept
        # the library as long as the expected friendly sets are present.
        if lib and ("Hypoxia" in lib or "Unfolded Protein Response" in lib):
            source = ("Enrichr GMT download MSigDB_Hallmark_2020 "
                      "(https://maayanlab.cloud/Enrichr/geneSetLibrary"
                      "?mode=text&libraryName=MSigDB_Hallmark_2020)")
            return lib, source
    except Exception as e:
        sys.stderr.write(f"[Enrichr download failed] {e}\n")
    # (iii) local gmt search
    try:
        res = subprocess.run(
            ["bash", "-lc",
             "find / -iname '*hallmark*.gmt' 2>/dev/null; "
             "find /opt /usr /home -iname '*.gmt' 2>/dev/null | head -50"],
            capture_output=True, text=True, timeout=120)
        candidates = [l for l in res.stdout.splitlines() if l.strip()]
        for c in candidates:
            lib = parse_gmt(c)
            if lib and any("HALLMARK" in k.upper() for k in lib):
                source = f"local gmt file {c}"
                return lib, source
    except Exception as e:
        sys.stderr.write(f"[local gmt search failed] {e}\n")
    return None, None


def parse_gmt(path):
    out = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                name = parts[0]
                genes = [g.strip() for g in parts[2:] if g.strip()]
                # strip trailing weights like GENE,1.0
                genes = [g.split(",")[0] for g in genes]
                out[name] = genes
    except Exception:
        return {}
    return out


HALLMARK, HALLMARK_SOURCE = load_hallmark()
if not HALLMARK:
    sys.stderr.write("FATAL: could not obtain real Hallmark sets. Aborting Task B.\n")
    print("HALLMARK_UNAVAILABLE")
    sys.exit(2)

# normalize keys
HK = {k.upper(): v for k, v in HALLMARK.items()}

# Map prescribed HALLMARK_* names to Enrichr MSigDB_Hallmark_2020 friendly names.
HALLMARK_ALIASES = {
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": [
        "INTERFERON ALPHA RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE"],
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": [
        "INTERFERON GAMMA RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"],
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": [
        "UNFOLDED PROTEIN RESPONSE", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
    "HALLMARK_HYPOXIA": ["HYPOXIA", "HALLMARK_HYPOXIA"],
}


def hk_get(name):
    name = name.upper()
    for alias in HALLMARK_ALIASES.get(name, [name]):
        if alias.upper() in HK:
            return HK[alias.upper()]
    if name in HK:
        return HK[name]
    return None


# ---------------------------------------------------------------------------
# 2) Config: corpora and panels
# ---------------------------------------------------------------------------
CORPORA = {
    "IFN": {
        "rna": "outputs/results/ifn-pilot/piper_gene_scores.csv",
        "micro": "outputs/results/ifn-cross-platform/piper_gene_scores.csv",
        "compact": "data/studies/gold/ifn_gold_panel.csv",
        "hallmark": ["HALLMARK_INTERFERON_ALPHA_RESPONSE",
                     "HALLMARK_INTERFERON_GAMMA_RESPONSE"],
    },
    "ER": {
        "rna": "outputs/results/er-stress-benchmark-primary/piper_gene_scores.csv",
        "micro": "outputs/results/er-stress-cross-platform/piper_gene_scores.csv",
        "compact": "data/studies/gold/er_stress_upr_gold_panel.csv",
        "hallmark": ["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
    },
    "Hypoxia": {
        "rna": "outputs/results/hypoxia-hif1-benchmark/piper_gene_scores.csv",
        "micro": "outputs/results/hypoxia-cross-platform/piper_gene_scores.csv",
        "compact": "data/studies/gold/hypoxia_hif1_gold_panel.csv",
        "hallmark": ["HALLMARK_HYPOXIA"],
    },
}


def load_panel_csv(path):
    df = pd.read_csv(os.path.join(ROOT, path))
    col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
    return set(str(g).strip().upper() for g in df[col].dropna())


def load_ranking(path):
    df = pd.read_csv(os.path.join(ROOT, path))
    # gene symbol col
    gcol = "gene_symbol" if "gene_symbol" in df.columns else None
    if gcol is None:
        for c in df.columns:
            if "gene_symbol" in c:
                gcol = c
                break
    # score col
    scol = None
    rcol = None
    if "quality_weighted_piper_score" in df.columns:
        scol = "quality_weighted_piper_score"
    if "quality_weighted_piper_rank" in df.columns:
        rcol = "quality_weighted_piper_rank"
    used = scol
    if scol is None:
        # fall back: any quality_weighted score-like column
        cands = [c for c in df.columns if "quality_weighted" in c
                 and "rank" not in c]
        if cands:
            scol = cands[0]
            used = scol
        else:
            scol = "piper_score"
            used = "piper_score (no quality_weighted column present)"
    df["_gene"] = df[gcol].astype(str).str.strip().str.upper()
    # rank ordering: score DESC, tie-break rank ASC
    if rcol is not None:
        df = df.sort_values([scol, rcol], ascending=[False, True])
    else:
        df = df.sort_values([scol], ascending=False)
    df = df.reset_index(drop=True)
    return df, "_gene", scol, used


def compute_metrics(df, gcol, scol, panel):
    universe = set(df[gcol])
    present = panel & universe
    n_present = len(present)
    res = {"panel_total": len(panel), "n_positives_present": n_present}
    if n_present == 0:
        for k in ("recall@100", "recall@200", "recall@500", "auroc", "auprc"):
            res[k] = float("nan")
        return res
    genes_in_order = df[gcol].tolist()
    for k in (100, 200, 500):
        topk = set(genes_in_order[:k])
        res[f"recall@{k}"] = len(present & topk) / n_present
    labels = df[gcol].isin(present).astype(int).values
    scores = df[scol].values
    res["auroc"] = auroc(scores, labels)
    res["auprc"] = auprc(scores, labels)
    return res


# ---------------------------------------------------------------------------
# 3-4) Compute and assemble comparison
# ---------------------------------------------------------------------------
rows = []
panel_sizes = {}
hallmark_set_sizes = {}

for topic, cfg in CORPORA.items():
    compact_panel = load_panel_csv(cfg["compact"])
    panel_sizes[(topic, "compact")] = len(compact_panel)
    # union of hallmark sets present
    hall_union = set()
    for hname in cfg["hallmark"]:
        g = hk_get(hname)
        if g:
            gs = set(x.strip().upper() for x in g)
            hallmark_set_sizes[hname] = len(gs)
            hall_union |= gs
    panel_sizes[(topic, "hallmark")] = len(hall_union)

    panels = {"compact": compact_panel, "hallmark": hall_union}

    metrics_store = {}
    for corpus in ("rna", "micro"):
        df, gcol, scol, used = load_ranking(cfg[corpus])
        for gtname, panel in panels.items():
            m = compute_metrics(df, gcol, scol, panel)
            m.update({
                "topic": topic,
                "ground_truth": gtname,
                "corpus": "RNA-only" if corpus == "rna" else "+micro",
                "score_column_used": used,
                "universe_size": len(set(df[gcol])),
            })
            metrics_store[(gtname, corpus)] = m
            rows.append(m)

    # deltas RNA-only -> +micro per ground truth
    for gtname in ("compact", "hallmark"):
        r = metrics_store[(gtname, "rna")]
        mi = metrics_store[(gtname, "micro")]
        delta = {
            "topic": topic,
            "ground_truth": gtname,
            "corpus": "DELTA(+micro - RNA-only)",
            "score_column_used": "",
            "universe_size": "",
            "panel_total": r["panel_total"],
            "n_positives_present": f"{r['n_positives_present']}->{mi['n_positives_present']}",
        }
        for k in ("recall@100", "recall@200", "recall@500", "auroc", "auprc"):
            try:
                delta[k] = mi[k] - r[k]
            except Exception:
                delta[k] = float("nan")
        rows.append(delta)

# build dataframe with stable column order
cols = ["topic", "ground_truth", "corpus", "panel_total", "n_positives_present",
        "universe_size", "recall@100", "recall@200", "recall@500",
        "auroc", "auprc", "score_column_used"]
out_df = pd.DataFrame(rows)
for c in cols:
    if c not in out_df.columns:
        out_df[c] = ""
out_df = out_df[cols]

# sort: topic, ground_truth, then RNA-only/+micro/DELTA
corder = {"RNA-only": 0, "+micro": 1, "DELTA(+micro - RNA-only)": 2}
out_df["_co"] = out_df["corpus"].map(corder).fillna(9)
out_df = out_df.sort_values(["topic", "ground_truth", "_co"]).drop(columns="_co")

csv_path = os.path.join(OUT, "hallmark_sensitivity_recall.csv")
out_df.to_csv(csv_path, index=False, float_format="%.4f")

# markdown
def fmt(v):
    if isinstance(v, float):
        if v != v:
            return "NA"
        return f"{v:.4f}"
    return str(v)

md = []
md.append("# Auxiliary Hallmark sensitivity recompute (AUXILIARY — does NOT replace locked compact panels)")
md.append("")
md.append("**AUXILIARY sensitivity analysis only.** This broad MSigDB Hallmark ground truth is")
md.append("provided alongside the compact locked panels so the reader can choose which ground")
md.append("truth to standardize on. It does NOT replace the locked compact gold panels.")
md.append("")
md.append(f"- Hallmark source: {HALLMARK_SOURCE}")
md.append(f"- sklearn available: {HAVE_SK}")
md.append("- Score column: quality_weighted_piper_score DESC, tie-break quality_weighted_piper_rank ASC")
md.append("")
md.append("## Hallmark set sizes (raw, from source)")
for k, v in sorted(hallmark_set_sizes.items()):
    md.append(f"- {k}: {v} genes")
md.append("")
md.append("## Comparison table")
md.append("")
header = "| topic | ground_truth | corpus | panel_total | n_pos_present | universe | recall@100 | recall@200 | recall@500 | AUROC | AUPRC |"
md.append(header)
md.append("|" + "---|" * 11)
for _, r in out_df.iterrows():
    md.append("| " + " | ".join([
        fmt(r["topic"]), fmt(r["ground_truth"]), fmt(r["corpus"]),
        fmt(r["panel_total"]), fmt(r["n_positives_present"]),
        fmt(r["universe_size"]),
        fmt(r["recall@100"]), fmt(r["recall@200"]), fmt(r["recall@500"]),
        fmt(r["auroc"]), fmt(r["auprc"]),
    ]) + " |")
md.append("")
md.append("## Notes")
md.append("")
md.append("- recall@k = |panel ∩ top-k| / n_positives_present (n_positives_present = |panel ∩ scored universe|).")
md.append("- Background-negative AUROC/AUPRC: positives = panel ∩ universe; negatives = all other scored genes; score = quality_weighted_piper_score.")
md.append("- Hallmark sets are large (~100-200 genes) so they are not fully present in each corpus universe; see n_pos_present vs panel_total.")
md.append("- Compact panels are the manuscript-locked gold panels and remain the primary ground truth.")
md.append("")
md_path = os.path.join(OUT, "hallmark_sensitivity_recall.md")
with open(md_path, "w") as f:
    f.write("\n".join(md) + "\n")

# source sidecar
src = []
src.append("AUXILIARY MSigDB Hallmark sensitivity recompute (does NOT replace locked compact panels).")
src.append("Command: PYTHONPATH=outputs/code python3 outputs/results/code-review/_hallmark_sensitivity.py")
src.append(f"Hallmark source/version: {HALLMARK_SOURCE}")
src.append(f"Hallmark sets used: {json.dumps(hallmark_set_sizes)}")
src.append(f"sklearn_available: {HAVE_SK}")
src.append("Compact gold panels: data/studies/gold/{ifn_gold_panel,er_stress_upr_gold_panel,hypoxia_hif1_gold_panel}.csv")
src.append("Ranking: quality_weighted_piper_score DESC, tie-break quality_weighted_piper_rank ASC.")
with open(csv_path + ".source", "w") as f:
    f.write("\n".join(src) + "\n")
with open(md_path + ".source", "w") as f:
    f.write("\n".join(src) + "\n")

# print machine-readable summary for the agent
print("HALLMARK_SOURCE=" + str(HALLMARK_SOURCE))
print("HALLMARK_SET_SIZES=" + json.dumps(hallmark_set_sizes))
print("PANEL_SIZES=" + json.dumps({f"{t}/{g}": s for (t, g), s in panel_sizes.items()}))
print("SKLEARN=" + str(HAVE_SK))
print("=== TABLE ===")
print(out_df.to_string(index=False))
print("=== END ===")
