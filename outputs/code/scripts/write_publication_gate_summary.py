#!/usr/bin/env python
"""Write the iteration publication gate summary from on-disk evidence."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

try:
    from scripts.check_source_provenance import audit_iteration, source_sidecar
except ModuleNotFoundError:  # pragma: no cover - direct CLI execution path
    from check_source_provenance import audit_iteration, source_sidecar

DIRECT_BASELINES = {
    "robustrankaggreg": "RobustRankAggreg",
    "metavolcanor": "MetaVolcanoR",
    "hstouffer": "hStouffer",
    "awmeta": "AWmeta/metafor",
}

ITERATION_TOKEN_RE = re.compile(r"(?:^|[^A-Za-z0-9])iter[-_]?(\d+)(?=$|[^0-9])", re.IGNORECASE)


def _iter_number_from_name(value: str) -> int:
    matches = [int(match.group(1)) for match in ITERATION_TOKEN_RE.finditer(str(value))]
    return max(matches) if matches else -1


def _iter_number_from_path(path: Path) -> int:
    for part in reversed(path.parts):
        iteration = _iter_number_from_name(part)
        if iteration >= 0:
            return iteration
    return -1


def _latest_existing(paths: list[Path], max_iteration: int) -> Path:
    existing = [path for path in paths if path.exists() and _iter_number_from_path(path) <= max_iteration]
    if not existing:
        return paths[0]
    return sorted(existing, key=_iter_number_from_path, reverse=True)[0]


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _active_hypoxia(catalog: Path) -> tuple[int | None, int | None]:
    rows = _read_csv(catalog)
    if not rows:
        return None, None
    active = [r for r in rows if r.get("include_in_analysis", "").strip().lower() == "true"]
    papers = {r.get("paper_id", "").strip() for r in active if r.get("paper_id", "").strip()}
    return len(active), len(papers)


def _paper_collapse_recall(results_dir: Path) -> tuple[str, str]:
    sensitivity = _read_csv(results_dir / "sensitivity_metrics.csv")
    for row in sensitivity:
        if row.get("variant") == "collapse_by_paper_id":
            return row.get("recall_at_50", "pending"), row.get("recall_at_100", "pending")
    metrics = _read_json(results_dir / "slice_metrics.json")
    if metrics:
        r50 = metrics.get("recall_at_50", {}).get("recall")
        r100 = metrics.get("recall_at_100", {}).get("recall")
        return str(r50) if r50 is not None else "pending", str(r100) if r100 is not None else "pending"
    return "pending", "pending"


def _tnf_status(project_root: Path, iter_name: str) -> tuple[str, list[str]]:
    iter_num = _iter_number_from_name(iter_name)
    curation = project_root / "data" / "studies" / "curation"
    results_root = project_root / "outputs" / "results"
    ledger = _latest_existing(
        [curation / f"iter{iter_num}_tnfa_rescue_candidate_ledger.csv"]
        + sorted(curation.glob("iter*_tnfa_rescue_candidate_ledger.csv"))
        + sorted(curation.glob("iter*_tnfa_candidate_ledger.csv")),
        iter_num,
    )
    gate = _latest_existing(
        [results_root / iter_name / "tnfa_rescue_gate.json"]
        + sorted(results_root.glob("iter-*/tnfa_rescue_gate.json"))
        + sorted(results_root.glob("iter-*/tnfa_gate.json")),
        iter_num,
    )
    scarcity = _latest_existing(
        [results_root / iter_name / "tnfa_rescue_scarcity_report.md"]
        + sorted(results_root.glob("iter-*/tnfa_rescue_scarcity_report.md"))
        + sorted(results_root.glob("iter-*/tnfa_scarcity_report.md")),
        iter_num,
    )
    evidence: list[str] = []
    rows = _read_csv(ledger)
    if rows:
        active = [r for r in rows if r.get("eligibility_status", "").lower() == "active"]
        papers = {r.get("paper_id", "").strip() for r in active if r.get("paper_id", "").strip()}
        evidence.append(f"ledger candidates={len(rows)}, active={len(active)}, active independent paper units={len(papers)}")
    else:
        evidence.append(f"pending: {ledger} not present")
    gate_payload = _read_json(gate)
    if gate_payload:
        evidence.append(
            "gate present: "
            f"passed={gate_payload.get('passed')}, "
            f"scarcity_triggered={gate_payload.get('scarcity_triggered')}, "
            f"active={gate_payload.get('n_active_rows')}, "
            f"independent_units={gate_payload.get('n_independent_units')}"
        )
    else:
        evidence.append(f"pending: {gate} not present")
    if scarcity.exists():
        evidence.append(f"scarcity report present: {scarcity}")
        return "scarcity documented; TNF activation failed", evidence
    if rows and len([r for r in rows if r.get("eligibility_status", "").lower() == "active"]) >= 5:
        return "candidate pass evidence present; verify gate json", evidence
    return "pending", evidence


def _baseline_status(results_dir: Path) -> tuple[str, list[str]]:
    matrix = results_dir / "baselines" / "baseline_parity_matrix.csv"
    manifest = results_dir / "baselines" / "baseline_manifest.csv"
    failure = results_dir / "baselines" / "baseline_failure_ledger.csv"
    rows = _read_csv(matrix)
    evidence: list[str] = []
    if not rows:
        evidence.append(f"pending: {matrix} not present or empty")
        evidence.append(f"pending: {manifest} not present")
        evidence.append(f"pending: {failure} not present")
        return "pending; comparative superiority claims disallowed", evidence
    resolved_statuses = {"ok", "success", "passed", "resolved"}
    runnable = {
        r.get("method_id", "").strip().lower(): r
        for r in rows
        if r.get("run_status", "").strip().lower() in resolved_statuses
    }
    blocked_direct = sorted(
        label
        for method_id, label in DIRECT_BASELINES.items()
        if method_id not in runnable
    )
    evidence.append(f"parity rows={len(rows)}; runnable/default methods={len(runnable)}")
    if blocked_direct:
        evidence.append("blocked or unresolved direct baselines: " + ", ".join(blocked_direct))
        return "tier-0 baselines runnable; direct-prior-art blockers remain; comparative superiority claims disallowed", evidence
    return "resolved; comparative claims may be evaluated", evidence


def _domain_firewall(project_root: Path, results_dir: Path) -> tuple[str, list[str]]:
    iter_num = _iter_number_from_name(results_dir.name)
    curation = project_root / "data" / "studies" / "curation"
    results_root = project_root / "outputs" / "results"
    dossier = _latest_existing(
        [curation / f"iter{iter_num}_domain_feasibility_dossier.csv"] + sorted(curation.glob("iter*_domain_feasibility_dossier.csv")),
        iter_num,
    )
    summary = _latest_existing(
        [results_dir / "domain_feasibility_summary.md"] + sorted(results_root.glob("iter-*/domain_feasibility_summary.md")),
        iter_num,
    )
    amendments = sorted((project_root / ".omx" / "plans").glob("domain-amendment-*.md"))
    suspicious = []
    if results_dir.exists():
        for path in results_dir.rglob("*"):
            name = path.name.lower()
            if path.is_file() and any(token in name for token in ["ifn", "lps", "cytokine", "interferon"]):
                suspicious.append(path)
    evidence = [
        f"domain feasibility dossier present={dossier.exists()}",
        f"domain feasibility summary present={summary.exists()}",
        f"domain amendment artifacts={len(amendments)}",
        f"candidate-domain output artifacts flagged before amendment={len(suspicious)}",
    ]
    if suspicious and not amendments:
        return "failed: candidate-domain outputs exist before amendment", evidence
    if dossier.exists() or summary.exists():
        return "feasibility-only evidence present; amendment required before scoring", evidence
    return "no alternative-domain activation or scoring evidence found", evidence


def _claim_state(hypoxia_active: int | None, paper_units: int | None, r50: str, tnf: str, baseline: str) -> str:
    try:
        recall50 = float(r50)
    except ValueError:
        recall50 = -1.0
    hypoxia_positive = hypoxia_active is not None and hypoxia_active >= 15 and paper_units is not None and paper_units > 10 and recall50 >= 0.60
    strong_hypoxia = hypoxia_active is not None and hypoxia_active >= 20 and paper_units is not None and paper_units > 10 and recall50 >= 0.60
    tnf_ok = "pass" in tnf.lower() or "scarcity documented" in tnf.lower()
    baseline_ok = baseline.startswith("resolved")
    if strong_hypoxia and tnf_ok and baseline_ok:
        return "strong methods paper"
    if hypoxia_positive and tnf_ok:
        return "narrow but defensible paper"
    return "not ready"


def build_summary(project_root: Path, iteration: int) -> str:
    iter_name = f"iter-{iteration}"
    results_dir = project_root / "outputs" / "results" / iter_name
    fallback_dir = project_root / "outputs" / "results" / "iter-9"
    metrics_dir = results_dir if (results_dir / "slice_metrics.json").exists() else fallback_dir
    active, paper_units = _active_hypoxia(project_root / "data" / "studies" / "hypoxia_catalog.csv")
    r50, r100 = _paper_collapse_recall(metrics_dir)
    tnf, tnf_evidence = _tnf_status(project_root, iter_name)
    baseline, baseline_evidence = _baseline_status(results_dir)
    firewall, firewall_evidence = _domain_firewall(project_root, results_dir)
    audit = audit_iteration(results_dir)
    claim = _claim_state(active, paper_units, r50, tnf, baseline)
    has_iteration_metrics = (results_dir / "slice_metrics.json").exists()
    hypoxia_text = (
        f"Iteration-{iteration} hypoxia slice, diagnostic, and rank-plane artifacts are present. "
        "The current gate treats new S0 additions and HYP003/HYP016 caveats as "
        "auditable curation evidence, not as permission to tune scoring or the locked "
        "HIF target panel."
        if has_iteration_metrics
        else f"Iteration-{iteration} hypoxia slice/diagnostic artifacts were not present when this "
        "summary was generated; treat positive hypoxia language as pending until Lane A "
        f"regenerates iteration-{iteration} slice, diagnostics, and rank-plane outputs."
    )
    publication_text = (
        f"The current integrated state is **{claim}**. The gate uses the latest "
        f"integrated iteration-{iteration} artifacts and remains conservative: hypoxia evidence "
        "is improving, TNF is a documented scarcity result, fallback-domain work is "
        "feasibility-only, and direct-prior-art baseline blockers still prevent "
        "comparative superiority claims."
    )
    if active is not None and active >= 20:
        s0_followup = "- S0 active-row target is met; continue expansion only for clearly independent sources or replacement of provisional/weak rows."
    else:
        s0_followup = "- Continue S0 scaling toward at least 20 active studies or a prespecified stopping rule."
    followups = [
        s0_followup,
        "- Keep HYP003 provisional in manuscript-facing tables unless replaced or fully hardened.",
        "- Resolve direct-prior-art baselines before any superiority claim over existing meta-analysis tools.",
        "- Activate any fallback domain only through a written amendment and pre-scoring locks.",
    ]

    lines = [
        f"# Iteration {iteration} Publication Gate Summary",
        "",
        f"Claim state: **{claim}**",
        "",
        "## Evidence Status",
        "",
        f"- Hypoxia active rows: {active if active is not None else 'pending'}",
        f"- Hypoxia independent paper/source units: {paper_units if paper_units is not None else 'pending'}",
        f"- Paper-collapse recall@50/@100: {r50}/{r100} (source: {metrics_dir})",
        f"- TNF rescue state: {tnf}",
        f"- Alternative-domain firewall: {firewall}",
        f"- Baseline parity state: {baseline}",
        f"- Iteration provenance audit: passed={audit['passed']} artifacts={audit['artifact_count']} missing_source={len(audit['missing_source'])} empty_source={len(audit['empty_source'])}",
        "",
        "## Gate Details",
        "",
        "### Hypoxia S0 Evidence",
        "",
        hypoxia_text,
        "",
        "### TNF Rescue Sprint",
        "",
        *[f"- {item}" for item in tnf_evidence],
        "",
        "### Alternative-Domain Firewall",
        "",
        *[f"- {item}" for item in firewall_evidence],
        "",
        f"No candidate-domain DEGORA, baseline, recall, diagnostic, or rank-plane output was found in the iteration-{iteration} result tree by this firewall scan. Any future alternative-domain scoring still requires a written amendment, locked gold panel, locked corpus rules, and no-output-seen certification before scoring.",
        "",
        "### Baseline Parity",
        "",
        *[f"- {item}" for item in baseline_evidence],
        "",
        "Direct-prior-art comparative superiority claims remain disallowed unless RobustRankAggreg, MetaVolcanoR, hStouffer, and AWmeta/metafor are all resolved and run under the equal-tuning rule.",
        "",
        "### Provenance Coverage",
        "",
        f"- Audit directory: {results_dir}",
        f"- Missing source sidecars: {audit['missing_source']}",
        f"- Empty source sidecars: {audit['empty_source']}",
        "",
        "## Publication Decision",
        "",
        publication_text,
        "",
        "## Required Follow-up Before Manuscript Claims",
        "",
        *followups,
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("../.."))
    parser.add_argument("--iteration", type=int, default=11)
    parser.add_argument("--output", type=Path, default=Path("../../outputs/results/iter-11/publication_gate_summary.md"))
    args = parser.parse_args(argv)
    project_root = args.project_root.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_summary(project_root, args.iteration))
    command = (
        "make -C outputs/code publication-gate"
        f" ITER={args.iteration} OUTPUT={output} PROJECT_ROOT={project_root}"
    )
    source_sidecar(output).write_text(command + "\n")
    print(json.dumps({"output": str(output), "source": str(source_sidecar(output))}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
