from __future__ import annotations

import hashlib
import io
import json
import urllib.error
import zipfile
from pathlib import Path

import pytest

from degora.discovery import DiscoveryError, DiscoveryUnavailableError
from degora.discovery_prepare import prepare_publication_records
from degora.discovery_run import run_discovery_analysis


class FakeTransport:
    def __init__(self, payloads: dict[str, bytes]) -> None:
        self.payloads = payloads
        self.urls: list[str] = []

    def validate_url(self, url: str) -> str:
        if not url.startswith("https://"):
            raise DiscoveryError("public-source URLs must use HTTPS")
        return url

    def get_bytes(self, url: str, *, max_bytes: int) -> bytes:
        self.urls.append(url)
        payload = self.payloads[url]
        if len(payload) > max_bytes:
            raise DiscoveryError("remote response exceeds the safety cap")
        return payload


def _deg_table() -> bytes:
    return b"gene,log2FoldChange,pvalue,padj\nTP53,2.1,0.001,0.003\nVEGFA,-1.2,0.02,0.04\n"


def _matrix_table() -> bytes:
    return b"gene,ctrl1,ctrl2,treat1,treat2\nTP53,1,1.1,3,3.2\nVEGFA,6,6.2,2,2.3\n"


def _record(**overrides) -> dict:
    record = {
        "canonical_id": "pmid:1",
        "source_unit_id": "PMID:1",
        "pmid": "1",
        "title": "Publication",
        "species_decision": "target_species_verified",
        "target_species_verified": True,
        "target_species_evidence": "The downloadable table is explicitly Homo sapiens.",
        "direct_file_candidates": [
            {
                "url": "https://zenodo.org/files/deg.csv",
                "name": "author_DESeq2_results.csv",
                "role": "deg_table",
            }
        ],
    }
    record.update(overrides)
    return record


def _zip(entries: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, payload in entries.items():
            archive.writestr(name, payload)
    return buffer.getvalue()


def test_no_geo_author_deg_table_materializes_full_bundle(tmp_path: Path) -> None:
    result = prepare_publication_records(
        [_record()],
        "human",
        query="hypoxia",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
    )

    assert result["species"]["key"] == "human"
    assert result["query"] == "hypoxia"
    assert result["returned_studies"] == 1
    candidate = result["studies"][0]["files"][0]
    assert candidate["candidate_id"]
    assert candidate["role"] == "deg_table"
    assert candidate["inspection"]["status"] == "ready_for_review"
    assert candidate["inspection"]["fetch_scope"] == "full"
    assert candidate["inspection"]["full_file_sha256"] == hashlib.sha256(_deg_table()).hexdigest()
    assert Path(candidate["inspection"]["local_path"]).is_file()
    assert (tmp_path / "prepared" / ".degora-discovery-bundle.json").is_file()
    audit = json.loads(Path(result["exports"]["audit_json"]).read_text(encoding="utf-8"))
    assert audit["studies"][0]["files"][0]["inspection"]["local_path"] == candidate["inspection"]["local_path"]


def test_upstream_fallback_matrix_is_inspected_as_full_candidate(tmp_path: Path) -> None:
    record = _record(
        canonical_id="doi:10.1/matrix",
        source_unit_id="DOI:10.1/matrix",
        doi="10.1/matrix",
        direct_file_candidates=[
            {
                "url": "https://zenodo.org/files/expression_matrix.csv",
                "name": "normalized_expression_matrix.csv",
                "role": "normalized_expression_matrix",
            }
        ],
    )
    result = prepare_publication_records(
        [record],
        "Homo sapiens",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/expression_matrix.csv": _matrix_table()}),
    )

    candidate = result["studies"][0]["files"][0]
    assert candidate["role"] == "normalized_expression_matrix"
    assert candidate["inspection"]["status"] == "upstream_matrix_ready_for_contrast"
    assert candidate["inspection"]["sample_columns"] == ["ctrl1", "ctrl2", "treat1", "treat2"]


def test_direct_preparation_prioritizes_author_deg_before_bounded_matrix_fallback(tmp_path: Path) -> None:
    matrix_candidates = [
        {
            "url": f"https://zenodo.org/files/matrix_{index}.csv",
            "name": f"normalized_expression_matrix_{index}.csv",
            "role": "normalized_expression_matrix",
        }
        for index in range(4)
    ]
    author_url = "https://zenodo.org/files/z_author_DESeq2_results.csv"
    record = _record(
        direct_file_candidates=[
            *matrix_candidates,
            {"url": author_url, "name": "z_author_DESeq2_results.csv", "role": "deg_table"},
        ]
    )
    transport = FakeTransport(
        {
            **{candidate["url"]: _matrix_table() for candidate in matrix_candidates},
            author_url: _deg_table(),
        }
    )

    result = prepare_publication_records(
        [record],
        "human",
        materialize_dir=tmp_path / "prepared",
        max_files_per_record=1,
        transport=transport,
    )

    assert transport.urls == [author_url]
    assert result["studies"][0]["files"][0]["role"] == "deg_table"


def test_transient_candidate_failure_does_not_cancel_other_candidates_or_records(tmp_path: Path) -> None:
    failed_url = "https://zenodo.org/files/a_author_DEG.csv"
    recovered_url = "https://zenodo.org/files/b_author_DEG.csv"
    other_url = "https://zenodo.org/files/c_author_DEG.csv"

    class TransientTransport(FakeTransport):
        def get_bytes(self, url: str, *, max_bytes: int) -> bytes:
            self.urls.append(url)
            if url == failed_url:
                raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
            payload = self.payloads[url]
            assert len(payload) <= max_bytes
            return payload

    first = _record(
        direct_file_candidates=[
            {"url": failed_url, "name": "a_author_DEG.csv", "role": "deg_table"},
            {"url": recovered_url, "name": "b_author_DEG.csv", "role": "deg_table"},
        ]
    )
    second = _record(
        canonical_id="pmid:2",
        source_unit_id="PMID:2",
        pmid="2",
        direct_file_candidates=[{"url": other_url, "name": "c_author_DEG.csv", "role": "deg_table"}],
    )
    result = prepare_publication_records(
        [first, second],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=TransientTransport({recovered_url: _deg_table(), other_url: _deg_table()}),
    )

    assert result["returned_studies"] == 2
    first_study = next(study for study in result["studies"] if study["canonical_id"] == "pmid:1")
    assert first_study["files"][0]["name"].endswith("b_author_DEG.csv")
    assert first_study["candidate_errors"][0]["status"] == "unavailable"


def test_safe_nested_zip_extracts_only_tabular_members(tmp_path: Path) -> None:
    inner = _zip({"nested/author_DEG.tsv": b"gene\tlog2FC\tpvalue\nTP53\t2\t0.01\nVEGFA\t-1\t0.02\n"})
    outer = _zip({"notes/readme.md": b"ignore", "inner.zip": inner})
    result = prepare_publication_records(
        [_record(direct_file_candidates=[{"url": "https://zenodo.org/files/bundle.zip", "name": "bundle.zip"}])],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/bundle.zip": outer}),
    )

    files = result["studies"][0]["files"]
    assert len(files) == 1
    assert files[0]["name"] == "author_DEG.tsv"
    assert files[0]["inspection"]["status"] == "ready_for_review"
    assert Path(files[0]["inspection"]["local_path"]).is_file()
    assert not any(path.name == "readme.md" for path in (tmp_path / "prepared").rglob("*"))


def test_archive_member_name_sanitization_never_overwrites_a_peer(tmp_path: Path) -> None:
    first = b"gene,log2FoldChange,pvalue,padj\nTP53,2.1,0.001,0.003\n"
    second = b"gene,log2FoldChange,pvalue,padj\nVEGFA,-1.2,0.02,0.04\n"
    bundle = _zip({"a/b.csv": first, "a_b.csv": second})
    result = prepare_publication_records(
        [_record(direct_file_candidates=[{"url": "https://zenodo.org/files/collision.zip", "name": "collision.zip"}])],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/collision.zip": bundle}),
    )

    files = result["studies"][0]["files"]
    assert len(files) == 2
    local_paths = [Path(item["inspection"]["local_path"]) for item in files]
    assert len(set(local_paths)) == 2
    expected = {"b.csv": first, "a_b.csv": second}
    for item, local_path in zip(files, local_paths, strict=True):
        assert local_path.read_bytes() == expected[item["name"]]
        assert item["inspection"]["full_file_sha256"] == hashlib.sha256(
            expected[item["name"]]
        ).hexdigest()


def test_zip_slip_and_bomb_rejection_rolls_back(tmp_path: Path) -> None:
    target = tmp_path / "prepared"
    slip = _zip({"../escape.csv": _deg_table()})
    with pytest.raises(DiscoveryError, match="unsafe member path"):
        prepare_publication_records(
            [_record(direct_file_candidates=[{"url": "https://zenodo.org/files/slip.zip", "name": "slip.zip"}])],
            "human",
            materialize_dir=target,
            transport=FakeTransport({"https://zenodo.org/files/slip.zip": slip}),
        )
    assert not target.exists()

    huge = _zip({"huge.csv": b"x" * (51 * 1024 * 1024)})
    with pytest.raises(DiscoveryError, match="oversized member"):
        prepare_publication_records(
            [_record(direct_file_candidates=[{"url": "https://zenodo.org/files/huge.zip", "name": "huge.zip"}])],
            "human",
            materialize_dir=target,
            transport=FakeTransport({"https://zenodo.org/files/huge.zip": huge}),
        )
    assert not target.exists()


def test_mixed_quarantine_and_rescue_policy(tmp_path: Path) -> None:
    blocked = _record(canonical_id="pmid:blocked", mixed_quarantined=True)
    unverified = _record(canonical_id="pmid:unverified", mixed_rescued=True, target_species_verified=True, target_species_evidence="")
    rescued = _record(
        canonical_id="pmid:rescued",
        source_unit_id="PMID:rescued",
        mixed_rescued=True,
        target_species_verified=True,
        target_species_evidence="Supplement labels this sheet as Homo sapiens only.",
    )
    result = prepare_publication_records(
        [blocked, unverified, rescued],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
    )

    assert result["returned_studies"] == 1
    assert result["studies"][0]["canonical_id"] == "pmid:rescued"
    reasons = {item["canonical_id"]: item["reason"] for item in result["excluded_studies"]}
    assert "mixed_quarantined" in reasons["pmid:blocked"]
    assert "mixed_rescued" in reasons["pmid:unverified"]


def test_publication_only_record_is_excluded_not_crashing(tmp_path: Path) -> None:
    result = prepare_publication_records(
        [_record(direct_file_candidates=[])],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({}),
    )

    assert result["returned_studies"] == 0
    assert result["excluded_studies"][0]["reason"] == "publication record has no GEO accession or direct file candidate"
    assert (tmp_path / "prepared" / "discovery_audit.json").is_file()


def test_transaction_rollback_preserves_existing_bundle_on_failed_force(tmp_path: Path) -> None:
    target = tmp_path / "prepared"
    prepare_publication_records(
        [_record()],
        "human",
        materialize_dir=target,
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
    )
    before = {path.relative_to(target): path.read_bytes() for path in target.rglob("*") if path.is_file()}

    with pytest.raises(KeyError):
        prepare_publication_records(
            [_record(direct_file_candidates=[{"url": "https://zenodo.org/files/missing.csv", "name": "missing.csv"}])],
            "human",
            materialize_dir=target,
            transport=FakeTransport({}),
            force=True,
        )

    after = {path.relative_to(target): path.read_bytes() for path in target.rglob("*") if path.is_file()}
    assert after == before


def test_prepared_author_candidate_is_compatible_with_run_discovery_analysis(tmp_path: Path) -> None:
    prepared = prepare_publication_records(
        [_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
    )
    candidate_id = prepared["studies"][0]["files"][0]["candidate_id"]
    result = run_discovery_analysis(
        prepared,
        [
            {
                "candidate_id": candidate_id,
                "mode": "author",
                "contrast_label": "hypoxia versus normoxia",
                "direction_confirmed": True,
                "table_scope": "full_results",
                "n_ctrl": 3,
                "n_treat": 3,
                "cell_system": "HK-2",
                "duration_h": "24",
                "platform": "public table",
            }
        ],
        tmp_path / "analysis",
        species="human",
        min_studies=1,
    )

    assert result["status"] == "complete"
    assert Path(result["db_path"]).is_file()
