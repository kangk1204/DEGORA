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


def test_prepare_publication_before_publish_runs_once_immediately_before_publish(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    events: list[str] = []
    import degora.discovery_prepare as prepare_module

    real_publish = prepare_module._publish_prepared_bundle

    def wrapped_publish(staging: Path, target: Path, *, force: bool) -> None:
        events.append("publish")
        real_publish(staging, target, force=force)

    monkeypatch.setattr(prepare_module, "_publish_prepared_bundle", wrapped_publish)

    prepare_publication_records(
        [_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
        before_publish=lambda: events.append("before"),
    )

    assert events == ["before", "publish"]


def test_prepare_publication_before_publish_exception_aborts_publish(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    events: list[str] = []
    import degora.discovery_prepare as prepare_module

    def wrapped_publish(staging: Path, target: Path, *, force: bool) -> None:
        events.append("publish")

    def fail_before_publish() -> None:
        events.append("before")
        raise RuntimeError("commit barrier failed")

    monkeypatch.setattr(prepare_module, "_publish_prepared_bundle", wrapped_publish)

    with pytest.raises(RuntimeError, match="commit barrier failed"):
        prepare_publication_records(
            [_record()],
            "human",
            materialize_dir=tmp_path / "prepared",
            transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
            before_publish=fail_before_publish,
        )

    assert events == ["before"]
    assert not (tmp_path / "prepared").exists()


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


# --- one bad candidate must not cost the whole selection --------------------
# A researcher who selects 20 studies and loses all of them to a single
# retired supplementary file has to start the download over from nothing.


class FlakyTransport(FakeTransport):
    """Serves payloads, and raises whatever a real public source might raise."""

    def __init__(self, payloads: dict[str, bytes], failures: dict[str, Exception] | None = None) -> None:
        super().__init__(payloads)
        self.failures = failures or {}

    def get_bytes(self, url: str, *, max_bytes: int) -> bytes:
        self.urls.append(url)
        failure = self.failures.get(url)
        if failure is not None:
            raise failure
        return super().get_bytes(url, max_bytes=max_bytes)


def _html_error_page() -> bytes:
    return b"<!DOCTYPE html>\n<html><head><title>404 Not Found</title></head><body>Gone</body></html>\n"


def _zip_candidate(url: str, name: str = "supplement.zip") -> dict:
    return {"url": url, "name": name, "role": "deg_table"}


def _healthy_record() -> dict:
    return _record(canonical_id="pmid:2", source_unit_id="PMID:2", pmid="2", title="Healthy study")


def test_an_unreadable_archive_excludes_its_own_study_and_spares_the_rest(tmp_path: Path) -> None:
    broken = _record(
        canonical_id="pmid:broken",
        source_unit_id="PMID:broken",
        pmid="broken",
        title="Study whose supplement was retired",
        direct_file_candidates=[_zip_candidate("https://zenodo.org/files/supplement.zip")],
    )
    result = prepare_publication_records(
        [broken, _healthy_record()],
        "human",
        query="hypoxia renal epithelial",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport(
            {
                "https://zenodo.org/files/supplement.zip": _html_error_page(),
                "https://zenodo.org/files/deg.csv": _deg_table(),
            }
        ),
    )

    assert result["returned_studies"] == 1
    assert result["studies"][0]["canonical_id"] == "pmid:2"
    excluded = {item["canonical_id"]: item for item in result["excluded_studies"]}
    reason = excluded["pmid:broken"]["reason"]
    assert "not a valid ZIP" in reason
    # The reader is told what actually arrived, not just that a ZIP was expected.
    assert "web page" in reason
    assert excluded["pmid:broken"]["candidate_errors"][0]["status"] == "rejected"


@pytest.mark.parametrize(
    ("failure", "expected_status", "expected_fragment"),
    [
        (DiscoveryError("remote response exceeds the 26214400-byte safety cap"), "rejected", "safety cap"),
        (DiscoveryUnavailableError("public source did not respond"), "unavailable", "temporarily unavailable"),
    ],
)
def test_every_public_source_failure_class_is_survivable(
    tmp_path: Path, failure: Exception, expected_status: str, expected_fragment: str
) -> None:
    """DiscoveryError and DiscoveryUnavailableError are siblings, so both must be caught."""

    failing = _record(
        canonical_id="pmid:huge",
        source_unit_id="PMID:huge",
        pmid="huge",
        title="Study with an unusable candidate",
        direct_file_candidates=[
            {"url": "https://zenodo.org/files/huge.csv", "name": "huge_DESeq2.csv", "role": "deg_table"}
        ],
    )
    result = prepare_publication_records(
        [failing, _healthy_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FlakyTransport(
            {"https://zenodo.org/files/deg.csv": _deg_table()},
            {"https://zenodo.org/files/huge.csv": failure},
        ),
    )

    assert result["returned_studies"] == 1
    excluded = {item["canonical_id"]: item for item in result["excluded_studies"]}
    assert excluded["pmid:huge"]["candidate_errors"][0]["status"] == expected_status
    assert expected_fragment in excluded["pmid:huge"]["reason"]


def test_a_study_keeps_the_candidates_that_did_work(tmp_path: Path) -> None:
    mixed = _record(
        canonical_id="pmid:mixed",
        source_unit_id="PMID:mixed",
        pmid="mixed",
        title="Study with one good and one broken candidate",
        direct_file_candidates=[
            _zip_candidate("https://zenodo.org/files/broken.zip", name="broken.zip"),
            {"url": "https://zenodo.org/files/good.csv", "name": "author_DESeq2_results.csv", "role": "deg_table"},
        ],
    )
    result = prepare_publication_records(
        [mixed],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport(
            {
                "https://zenodo.org/files/broken.zip": b"\x1f\x8bnot really a zip",
                "https://zenodo.org/files/good.csv": _deg_table(),
            }
        ),
    )

    study = result["studies"][0]
    names = [item["name"] for item in study["files"]]
    assert len(names) == 1 and names[0].endswith("author_DESeq2_results.csv"), names
    assert study["candidate_errors"][0]["status"] == "rejected"
    assert "gzip" in study["candidate_errors"][0]["error"]


def test_a_corrupt_member_does_not_discard_the_tables_beside_it(tmp_path: Path) -> None:
    outer = _zip({"author_DESeq2_results.csv": _deg_table(), "extra.zip": b"not a zip at all"})
    record = _record(
        canonical_id="pmid:nested",
        source_unit_id="PMID:nested",
        pmid="nested",
        title="Study with a corrupt nested archive",
        direct_file_candidates=[_zip_candidate("https://zenodo.org/files/bundle.zip", name="bundle.zip")],
    )
    result = prepare_publication_records(
        [record],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/bundle.zip": outer}),
    )

    study = result["studies"][0]
    names = [item["name"] for item in study["files"]]
    assert len(names) == 1 and names[0].endswith("author_DESeq2_results.csv"), names
    assert any("extra.zip" in item["error"] for item in study["candidate_errors"])


def test_an_archive_safety_violation_still_refuses_the_whole_run(tmp_path: Path) -> None:
    """A malformed download is one study's problem; a hostile archive is not.

    Making unreadable files survivable must not quietly downgrade zip-slip and
    zip-bomb rejection into a skipped row nobody reads.
    """

    target = tmp_path / "prepared"
    slip = _zip({"../escape.csv": _deg_table()})
    with pytest.raises(DiscoveryError, match="unsafe member path"):
        prepare_publication_records(
            [
                _record(direct_file_candidates=[_zip_candidate("https://zenodo.org/files/slip.zip", name="slip.zip")]),
                _healthy_record(),
            ],
            "human",
            materialize_dir=target,
            transport=FakeTransport(
                {
                    "https://zenodo.org/files/slip.zip": slip,
                    "https://zenodo.org/files/deg.csv": _deg_table(),
                }
            ),
        )
    assert not target.exists()
    assert not (tmp_path / "escape.csv").exists()


def test_a_rejected_candidate_leaves_no_downloaded_file_behind(tmp_path: Path) -> None:
    prepared = tmp_path / "prepared"
    record = _record(
        canonical_id="pmid:rejected",
        source_unit_id="PMID:rejected",
        pmid="rejected",
        title="Study whose only candidate is unreadable",
        direct_file_candidates=[_zip_candidate("https://zenodo.org/files/broken.zip", name="broken.zip")],
    )
    result = prepare_publication_records(
        [record, _healthy_record()],
        "human",
        materialize_dir=prepared,
        transport=FakeTransport(
            {
                "https://zenodo.org/files/broken.zip": _html_error_page(),
                "https://zenodo.org/files/deg.csv": _deg_table(),
            }
        ),
    )

    assert result["returned_studies"] == 1
    assert not list(prepared.rglob("*broken.zip")), "the unusable download was shipped in the bundle"


# --- the reader has to be able to tell the arms apart -----------------------


class LabelledGeoClient:
    """A GEO series whose matrix columns are bare GSM accessions."""

    SAMPLE_SOFT = "\n".join(
        [
            "^SAMPLE = GSM320836",
            "!Sample_title = Normoxia replicate 1",
            "!Sample_source_name_ch1 = renal proximal tubule epithelial cells",
            "!Sample_characteristics_ch1 = treatment: normoxia",
            "!Sample_characteristics_ch1 = time: 24 h",
            "^SAMPLE = GSM320853",
            "!Sample_title = Hypoxia 1% O2 replicate 1",
            "!Sample_characteristics_ch1 = treatment: 1 percent oxygen",
            "^SAMPLE = GSM320854",
            "!Sample_title = Normoxia replicate 2",
            "^SAMPLE = GSM320855",
            "!Sample_title = Hypoxia 1% O2 replicate 2",
        ]
    )

    def __init__(self, *, with_labels: bool = True) -> None:
        self.with_labels = with_labels
        self.sample_calls: list[str] = []

    def accession_summaries(self, accessions, species):
        return [
            {
                "accession": accession,
                "taxon": "Homo sapiens",
                "title": "Renal epithelial hypoxia series",
                "gdstype": "Expression profiling by high throughput sequencing",
                "pubmedids": [],
                "pdat": "2015/01/01",
            }
            for accession in accessions
        ]

    def publication_summaries(self, pmids):
        return {}

    def fetch_geo_soft(self, accession):
        return "\n".join(
            [
                f"^SERIES = {accession}",
                "!Series_title = Renal epithelial hypoxia series",
                "!Series_overall_design = hypoxia versus normoxia",
                "!Series_sample_organism = Homo sapiens",
                f"!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/{accession}/suppl/{accession}_series_matrix.txt.gz",
            ]
        )

    def fetch_geo_sample_soft(self, accession):
        self.sample_calls.append(accession)
        if not self.with_labels:
            raise DiscoveryUnavailableError("GEO did not answer the sample request")
        return self.SAMPLE_SOFT

    def fetch_candidate(self, url, *, full):
        import gzip

        payload = gzip.compress(
            b"ID_REF\tGSM320836\tGSM320853\tGSM320854\tGSM320855\n"
            b"TP53\t10\t40\t11\t44\nCDKN1A\t4\t30\t5\t28\n"
        )
        return payload, "full" if full else "header_prefix"


def _geo_record() -> dict:
    return _record(
        canonical_id="pmid:geo",
        source_unit_id="PMID:geo",
        pmid="geo",
        title="Renal epithelial hypoxia",
        geo_accessions=["GSE100001"],
    )


def test_group_assignment_shows_the_submitter_labels_not_only_accessions(tmp_path: Path) -> None:
    """Choosing control and treatment from bare GSM ids means guessing."""

    client = LabelledGeoClient()
    result = prepare_publication_records(
        [_geo_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({}),
        geo_client=client,
    )

    study = result["studies"][0]
    labels = study["sample_labels"]
    assert labels["GSM320836"]["title"] == "Normoxia replicate 1"
    assert labels["GSM320853"]["title"] == "Hypoxia 1% O2 replicate 1"
    assert "treatment: normoxia" in labels["GSM320836"]["characteristics"]
    # Every column the reader is asked to assign has something to read.
    matrix = next(item for item in study["files"] if item.get("inspection", {}).get("sample_columns"))
    for column in matrix["inspection"]["sample_columns"]:
        assert labels.get(column.upper(), {}).get("title"), column


def test_a_label_lookup_failure_never_fails_the_preparation(tmp_path: Path) -> None:
    result = prepare_publication_records(
        [_geo_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({}),
        geo_client=LabelledGeoClient(with_labels=False),
    )

    assert result["returned_studies"] == 1
    assert result["studies"][0]["sample_labels"] == {}


def test_labels_are_not_fetched_when_nobody_will_assign_groups(tmp_path: Path) -> None:
    """The author-DEG path needs no group assignment, so it pays no extra request."""

    client = LabelledGeoClient()
    prepare_publication_records(
        [_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
        geo_client=client,
    )

    assert client.sample_calls == []


def test_publications_sharing_a_series_are_explained_not_dropped(tmp_path: Path) -> None:
    """A selection of N must reconcile to prepared + excluded, with nothing missing.

    Two publications reporting one GEO series are one dataset, so they stay one
    study - but the absorbed publication used to appear in neither list, which
    is how a selection of 20 quietly became 11 with no account of the rest.
    """

    first = _geo_record()
    second = _record(
        canonical_id="pmid:geo-2",
        source_unit_id="PMID:geo-2",
        pmid="geo-2",
        title="A second paper on the same series",
        geo_accessions=["GSE100001"],
    )
    result = prepare_publication_records(
        [first, second, _healthy_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
        geo_client=LabelledGeoClient(),
    )

    prepared = {study.get("canonical_id") for study in result["studies"]}
    excluded = {item.get("canonical_id"): item for item in result["excluded_studies"]}
    assert {"pmid:geo", "pmid:geo-2", "pmid:2"} <= prepared | set(excluded)
    assert "pmid:geo-2" in excluded
    assert "GSE100001" in excluded["pmid:geo-2"]["reason"]
    assert "already prepared" in excluded["pmid:geo-2"]["reason"]


def test_the_repository_phase_says_how_many_selections_it_covers(tmp_path: Path) -> None:
    """"12 repository record(s)" against a selection of 20 reads like a miscount."""

    messages: list[str] = []
    prepare_publication_records(
        [_geo_record(), _healthy_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
        geo_client=LabelledGeoClient(),
        progress=lambda fraction, message: messages.append(message),
    )

    repository = next(text for text in messages if "repository series" in text)
    assert "1 of 2 selected publications" in repository, repository


class AuthorMatrixGeoClient(LabelledGeoClient):
    """A series whose matrix is headed by the submitter's own column names."""

    SAMPLE_SOFT = "\n".join(
        [
            "^SAMPLE = GSM6072341",
            "!Sample_title = 4641CERM6M24M",
            "!Sample_characteristics_ch1 = transgene: induced",
            "^SAMPLE = GSM6072342",
            "!Sample_title = 4709CERM6m24M",
            "!Sample_characteristics_ch1 = transgene: induced",
            "^SAMPLE = GSM6072343",
            "!Sample_title = 4728CERM12M24M",
            "!Sample_characteristics_ch1 = transgene: uninduced",
            "^SAMPLE = GSM6072344",
            "!Sample_title = 4754CERM12M24M",
            "!Sample_characteristics_ch1 = transgene: uninduced",
        ]
    )

    def fetch_geo_soft(self, accession):
        return "\n".join(
            [
                f"^SERIES = {accession}",
                "!Series_title = Mouse mammary gland series",
                "!Series_overall_design = induced versus uninduced",
                "!Series_sample_organism = Homo sapiens",
                f"!Series_supplementary_file = https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/{accession}/suppl/{accession}_TPM_matrix.txt.gz",
            ]
        )

    def fetch_candidate(self, url, *, full):
        import gzip

        payload = gzip.compress(
            b"Gene\t4641CERM6M24M_S2\t4709CERM6m24M_S2\t4728CERM12M24M_S3\t4754CERM12M24M_S3\n"
            b"Stat5a\t4\t5\t6\t7\nEsr1\t9\t8\t7\t6\n"
        )
        return payload, "full" if full else "header_prefix"


def test_an_author_matrix_carries_a_label_for_every_column(tmp_path: Path) -> None:
    result = prepare_publication_records(
        [_geo_record()],
        "human",
        materialize_dir=tmp_path / "prepared",
        transport=FakeTransport({}),
        geo_client=AuthorMatrixGeoClient(),
    )

    matrix = next(
        item
        for item in result["studies"][0]["files"]
        if item.get("inspection", {}).get("sample_columns")
    )
    resolved = matrix["inspection"]["sample_labels"]
    columns = matrix["inspection"]["sample_columns"]
    assert set(resolved) == set(columns), (sorted(resolved), sorted(columns))
    assert resolved["4728CERM12M24M_S3"]["characteristics"] == ["transgene: uninduced"]
    # The mapping is recorded in the bundle, not only drawn on screen.
    audit = json.loads(Path(result["exports"]["audit_json"]).read_text(encoding="utf-8"))
    audit_matrix = next(
        item
        for item in audit["studies"][0]["files"]
        if item.get("inspection", {}).get("sample_labels")
    )
    assert audit_matrix["inspection"]["sample_labels"]["4641CERM6M24M_S2"]["accession"] == "GSM6072341"


def test_the_persisted_audit_points_at_the_published_bundle(tmp_path: Path) -> None:
    """Every path the archived audit records has to resolve for the reader.

    The audit JSON is written into a staging directory and then published under
    the target, and the staging directory is removed immediately afterwards. The
    export paths were captured before that move, so the persisted document -- the
    one kept for provenance -- described four files that no longer existed, while
    only the in-memory return value carried the real ones.
    """

    target = tmp_path / "prepared"
    result = prepare_publication_records(
        [_geo_record(), _healthy_record()],
        "human",
        materialize_dir=target,
        transport=FakeTransport({"https://zenodo.org/files/deg.csv": _deg_table()}),
        geo_client=LabelledGeoClient(),
    )

    persisted = json.loads((target / "discovery_audit.json").read_text(encoding="utf-8"))
    exports = persisted["exports"]
    assert set(exports) == {"output_dir", "audit_json", "candidates_csv", "draft_catalog_csv"}
    for key, value in exports.items():
        assert Path(value).exists(), f"persisted audit records a missing {key}: {value}"
        assert ".prepare-" not in value, f"persisted audit records a staging path for {key}: {value}"

    # The archived document and the returned object must agree about the bundle.
    assert exports == result["exports"]
    assert Path(exports["output_dir"]).resolve() == target.resolve()


def test_species_provenance_reaches_the_prepared_bundle_and_its_audit(tmp_path: Path) -> None:
    """The audit a reviewer opens has to say which species claim a study carries.

    The species gate reads these fields off the search record, so preparation kept
    working while dropping them: the returned study and `discovery_audit.json`
    both reported null, and the prepare UI renders its species provenance line
    only when `species_decision` survives. Assert the whole chain, not the helper
    -- a helper-level check is what missed this.
    """

    from degora.discovery_federated import _prepare_record, normalize_species

    url = "https://zenodo.org/files/deg.csv"
    record = _prepare_record(
        {
            "provider": "pubmed",
            "pmid": "41932308",
            "paper_title": "A literature-only publication",
            "species_evidence": [{"species": "Homo sapiens", "basis": "PubMed organism-constrained query"}],
            "direct_file_candidates": [{"name": "deg.csv", "source_url": url, "role": "deg_table"}],
        },
        normalize_species("human"),
    )
    assert record["species_decision"] == "query_constrained"

    target = tmp_path / "prepared"
    result = prepare_publication_records(
        [record],
        "human",
        materialize_dir=target,
        transport=FakeTransport({url: _deg_table()}),
    )

    study = result["studies"][0]
    assert study["species_decision"] == "query_constrained"
    assert study["species_evidence"] == [
        {"species": "Human", "basis": "PubMed organism-constrained query"}
    ]
    assert study["target_species_verified"] is False

    persisted = json.loads((target / "discovery_audit.json").read_text(encoding="utf-8"))["studies"][0]
    assert persisted["species_decision"] == study["species_decision"]
    assert persisted["target_species_verified"] == study["target_species_verified"]
    assert [item["basis"] for item in persisted["species_evidence"]] == [
        "PubMed organism-constrained query"
    ]
