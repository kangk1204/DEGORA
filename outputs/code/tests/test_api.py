from __future__ import annotations

import json
import threading
import urllib.request
from urllib.error import HTTPError

import pandas as pd

from degora.api import create_server
from degora.score_db import write_score_database


def _harmonized() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S1", "S2"],
            "paper_id": ["P1", "P2", "P1", "P2"],
            "gene_symbol": ["VEGFA", "VEGFA", "RPL13A", "RPL13A"],
            "lfc": [2.0, 1.8, 0.1, -0.1],
            "signed_z": [5.0, 4.5, 0.1, -0.1],
            "pvalue": [1e-7, 1e-6, 0.9, 0.8],
            "padj": [1e-5, 1e-4, 0.9, 0.9],
            "normalized_rank": [0.02, 0.03, 0.9, 0.8],
            "n_ctrl": [3, 4, 3, 4],
            "n_treat": [3, 4, 3, 4],
            "n_genes_in_study": [1000] * 4,
            "pipeline": ["DESeq2", "edgeR", "DESeq2", "edgeR"],
            "assay_type": ["RNA-seq", "microarray", "RNA-seq", "microarray"],
            "source_input_type": ["author_deg_table", "limma_full_table", "author_deg_table", "limma_full_table"],
            "platform": ["", "GPL570", "", "GPL570"],
            "normalization": ["", "RMA/log2", "", "RMA/log2"],
            "probe_collapse": ["", "author_gene_level", "", "author_gene_level"],
            "species": ["Homo sapiens"] * 4,
            "cell_system": ["A", "B", "A", "B"],
            "hypoxia_modality": ["1% O2"] * 4,
            "duration_h": ["24"] * 4,
            "source_path": ["source.csv"] * 4,
            "source_url": ["https://example.test"] * 4,
        }
    )


def _get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=5) as response:
        return json.loads(response.read().decode("utf-8"))


def test_local_api_serves_health_gene_list_and_detail(tmp_path) -> None:
    harmonized_path = tmp_path / "harmonized.csv"
    _harmonized().to_csv(harmonized_path, index=False)
    write_score_database(harmonized_path, tmp_path, db_path=tmp_path / "degora_scores.db")

    server = create_server(tmp_path / "degora_scores.db", port=0, quiet=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"

    try:
        with urllib.request.urlopen(base_url, timeout=5) as response:
            html = response.read().decode("utf-8")
        with urllib.request.urlopen(f"{base_url}/api/health", timeout=5) as response:
            assert response.headers.get("Access-Control-Allow-Origin") is None
        health = _get_json(f"{base_url}/api/health")
        genes = _get_json(f"{base_url}/api/genes?q=VEGF&limit=5")
        wildcard = _get_json(f"{base_url}/api/genes?q=_&limit=5")
        detail = _get_json(f"{base_url}/api/genes/VEGFA")
        try:
            _get_json(f"{base_url}/api/genes?q={'A' * 129}")
        except HTTPError as exc:
            assert exc.code == 400
        else:  # pragma: no cover - assertion branch
            raise AssertionError("long gene query should be rejected")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert "<title>DEGORA</title>" in html
    assert "const esc =" in html
    assert health["status"] == "ok"
    assert health["gene_count"] == 2
    assert genes["genes"][0]["gene_symbol"] == "VEGFA"
    assert wildcard["count"] == 0
    assert genes["genes"][0]["rank_label"] == "#1 / 2"
    assert "top_percent_label" in genes["genes"][0]
    assert "evidence_tier" in genes["genes"][0]
    assert detail["gene"]["gene_symbol"] == "VEGFA"
    assert detail["gene"]["support_label"] == "2 / 2 source units"
    assert detail["gene"]["direction_label"] == "100.0% up-concordant"
    assert len(detail["evidence"]) == 2
    assert {row["assay_type"] for row in detail["evidence"]} == {"RNA-seq", "microarray"}
    assert all("contributing_study_ids" in row for row in detail["evidence"])
