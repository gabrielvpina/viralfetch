"""Tests for the NCBI E-utilities client — parsing and transport behaviour.

No real network: a fake session serves responses from an in-memory registry.
"""

import json
from pathlib import Path

import pytest

from viralfetch.cache import Cache
from viralfetch.config import Config
from viralfetch.ncbi import (
    MAX_BATCH,
    NCBIClient,
    NCBIError,
    RateLimiter,
    _meta_from_json,
    _parse_elink,
    _parse_esummary,
    _parse_taxonomy_xml,
    _split_fasta,
    _split_genbank,
)

FIXTURES = Path(__file__).parent / "fixtures"


def noop_sleep(_seconds):
    pass


# -- parsing (frozen fixtures) --------------------------------------------

def test_parse_esummary_fixture():
    body = (FIXTURES / "esummary_nuccore.json").read_text()
    found = _parse_esummary(body)
    assert set(found) == {"MN908947", "AY274119"}
    meta = _meta_from_json(found["MN908947"])
    assert meta.accession == "MN908947.3"
    assert meta.length == 29903
    assert meta.sourcedb == "insd"


def test_split_fasta_fixture():
    body = (FIXTURES / "efetch_nuccore.fasta").read_text()
    records = _split_fasta(body)
    assert set(records) == {"MN908947", "AY274119"}
    assert records["MN908947"].startswith(">MN908947.3")
    assert records["AY274119"].rstrip().endswith("CTCAC")


def test_split_genbank_fixture():
    body = (FIXTURES / "efetch_nuccore.gb").read_text()
    records = _split_genbank(body)
    assert set(records) == {"MN908947", "AY274119"}
    assert records["MN908947"].startswith("LOCUS")
    assert records["MN908947"].rstrip().endswith("//")


# -- rate limiter ---------------------------------------------------------

def test_rate_limiter_spaces_requests():
    slept = []
    clock = iter([0.0, 0.0, 0.0, 0.4, 0.4])  # monotonic returns

    def fake_clock():
        return next(clock)

    rl = RateLimiter(per_second=3, sleep=slept.append, clock=fake_clock)
    rl.wait()  # first call, last=0 -> gap = 1/3, but now=0 so sleeps ~0.333
    rl.wait()
    assert slept  # it slept at least once to honour the interval


# -- fake session ---------------------------------------------------------

class FakeResponse:
    def __init__(self, text, status_code=200):
        self.text = text
        self.status_code = status_code


class FakeSession:
    """Serves esummary/efetch/elink from an in-memory registry.

    ``registry`` maps base-accession -> record dict (with a ``uid``).
    ``protein_links`` maps a source UID -> list of protein UIDs (for elink).
    ``taxonomy`` maps a source UID -> taxid; ``tax_xml`` maps taxid -> xml.
    """

    def __init__(self, registry: dict, fail_times: int = 0, fail_status: int = 503,
                 protein_links: dict | None = None, taxonomy: dict | None = None,
                 tax_xml: dict | None = None):
        self.registry = registry
        self.calls = []
        self.fail_times = fail_times
        self.fail_status = fail_status
        self.protein_links = protein_links or {}
        self.taxonomy = taxonomy or {}
        self.tax_xml = tax_xml or {}

    def post(self, url, data=None, timeout=None):
        self.calls.append(data)
        if self.fail_times > 0:
            self.fail_times -= 1
            return FakeResponse("", self.fail_status)
        ids = data["id"].split(",")
        if url.endswith("esummary.fcgi"):
            result = {"uids": []}
            for acc in ids:
                rec = self.registry.get(acc.split(".")[0])
                if rec:
                    result["uids"].append(rec["uid"])
                    result[rec["uid"]] = rec
            return FakeResponse(json.dumps({"result": result}))
        if url.endswith("elink.fcgi"):
            links = []
            for uid in ids:
                if data["db"] == "protein":
                    links += self.protein_links.get(uid, [])
                elif data["db"] == "taxonomy" and uid in self.taxonomy:
                    links.append(self.taxonomy[uid])
            return FakeResponse(json.dumps({"linksets": [
                {"linksetdbs": [{"links": links}]}
            ]}))
        if url.endswith("efetch.fcgi"):
            if data.get("db") == "taxonomy":
                return FakeResponse(self.tax_xml.get(ids[0], "<TaxaSet></TaxaSet>"))
            by_uid = {r["uid"]: r for r in self.registry.values()}
            chunks = []
            for ident in ids:
                rec = self.registry.get(ident.split(".")[0]) or by_uid.get(ident)
                if rec:
                    chunks.append(f">{rec['accessionversion']} {rec['organism']}\nACGT\n")
            return FakeResponse("".join(chunks))
        return FakeResponse("", 404)


REGISTRY = {
    "MN908947": {
        "uid": "1", "caption": "MN908947", "accessionversion": "MN908947.3",
        "organism": "SARS-CoV-2", "slen": "29903", "moltype": "rna",
        "biomol": "genomic", "topology": "linear", "completeness": "complete",
        "sourcedb": "insd", "updatedate": "2020/03/30",
    },
    "AY274119": {
        "uid": "2", "caption": "AY274119", "accessionversion": "AY274119.3",
        "organism": "SARS-CoV Tor2", "slen": "29751", "moltype": "rna",
        "biomol": "genomic", "topology": "linear", "completeness": "complete",
        "sourcedb": "refseq", "updatedate": "2018/08/13",
    },
}


def make_client(session, tmp_path, cache_enabled=True):
    cfg = Config(email="tester@example.com")
    cache = Cache(tmp_path, enabled=cache_enabled)
    return NCBIClient(cfg, cache=cache, session=session, sleep=noop_sleep)


def test_missing_email_raises():
    with pytest.raises(Exception):
        NCBIClient(Config(email=None), session=FakeSession(REGISTRY))


def test_esummary_partial_failure_reported(tmp_path):
    session = FakeSession(REGISTRY)
    client = make_client(session, tmp_path)
    result = client.esummary_nuccore(["MN908947", "AY274119", "ZZ000000"])
    got = {r.accession for r in result.records}
    assert got == {"MN908947.3", "AY274119.3"}
    assert result.missing == ["ZZ000000"]  # reported, not raised


def test_esummary_uses_cache_second_call(tmp_path):
    session = FakeSession(REGISTRY)
    client = make_client(session, tmp_path)
    client.esummary_nuccore(["MN908947"])
    n_calls = len(session.calls)
    client.esummary_nuccore(["MN908947"])  # should be fully cached
    assert len(session.calls) == n_calls  # no new request


def test_batching_splits_at_200(tmp_path):
    # 201 fake accessions -> two batches.
    reg = {f"ACC{i:05d}": {
        "uid": str(i), "caption": f"ACC{i:05d}",
        "accessionversion": f"ACC{i:05d}.1", "organism": "x", "slen": "1",
        "moltype": "rna", "biomol": "genomic", "topology": "linear",
        "completeness": "complete", "sourcedb": "insd", "updatedate": "2020",
    } for i in range(201)}
    session = FakeSession(reg)
    client = make_client(session, tmp_path)
    result = client.esummary_nuccore(list(reg))
    assert len(result.records) == 201
    assert len(session.calls) == 2
    assert MAX_BATCH == 200


def test_retry_on_5xx_then_success(tmp_path):
    session = FakeSession(REGISTRY, fail_times=2, fail_status=503)
    client = make_client(session, tmp_path)
    result = client.esummary_nuccore(["MN908947"])
    assert {r.accession for r in result.records} == {"MN908947.3"}
    assert len(session.calls) == 3  # 2 failures + 1 success


def test_retry_exhausted_raises(tmp_path):
    session = FakeSession(REGISTRY, fail_times=99, fail_status=503)
    client = make_client(session, tmp_path)
    with pytest.raises(NCBIError):
        client.esummary_nuccore(["MN908947"])


def test_efetch_fasta_missing_and_cache(tmp_path):
    session = FakeSession(REGISTRY)
    client = make_client(session, tmp_path)
    result = client.efetch_nuccore(["MN908947", "ZZ000000"], "fasta")
    assert result.returned == ["MN908947"]
    assert result.missing == ["ZZ000000"]
    assert ">MN908947.3" in result.text


# -- elink / taxonomy parsing (frozen fixtures) ---------------------------

def test_parse_elink_fixture():
    body = (FIXTURES / "elink_protein.json").read_text()
    assert _parse_elink(body) == ["1798172433", "1798172434", "1798172435"]


def test_parse_taxonomy_xml_fixture():
    lineage = _parse_taxonomy_xml((FIXTURES / "taxonomy.xml").read_text())
    assert lineage.taxid == "2697049"
    assert lineage.rank == "no rank"
    ranks = dict(lineage.lineage)
    assert ranks["realm"] == "Riboviria"
    assert ranks["genus"] == "Betacoronavirus"
    # The taxon itself is appended last.
    assert lineage.lineage[-1][1] == lineage.name


# -- elink-backed client paths (nuccore UID resolution first) -------------

def test_nuccore_uids_resolves_via_esummary(tmp_path):
    session = FakeSession(REGISTRY)
    client = make_client(session, tmp_path)
    assert client.nuccore_uids(["MN908947", "AY274119"]) == ["1", "2"]


def test_protein_uids_for_links(tmp_path):
    session = FakeSession(REGISTRY, protein_links={"1": ["101", "102"]})
    client = make_client(session, tmp_path)
    assert client.protein_uids_for(["MN908947"]) == ["101", "102"]


def test_protein_uids_empty_not_cached(tmp_path):
    session = FakeSession(REGISTRY, protein_links={})  # no links
    client = make_client(session, tmp_path)
    assert client.protein_uids_for(["MN908947"]) == []
    # An empty result must not be cached (it would poison later runs).
    key_calls = len(session.calls)
    client.protein_uids_for(["MN908947"])
    assert len(session.calls) > key_calls  # refetched, not served from cache


def test_taxid_for_accession(tmp_path):
    session = FakeSession(REGISTRY, taxonomy={"1": "2697049"})
    client = make_client(session, tmp_path)
    assert client.taxid_for_accession("MN908947") == "2697049"


def test_efetch_taxonomy_parses(tmp_path):
    xml = (FIXTURES / "taxonomy.xml").read_text()
    session = FakeSession(REGISTRY, tax_xml={"2697049": xml})
    client = make_client(session, tmp_path)
    lineage = client.efetch_taxonomy("2697049")
    assert lineage.name.startswith("Severe acute")


def test_efetch_all_returns_every_record(tmp_path):
    # Protein-style: ids are UIDs; efetch_all takes all split records.
    reg = {"YP_009": {
        "uid": "9", "caption": "YP_009", "accessionversion": "YP_009.1",
        "organism": "spike", "slen": "1", "moltype": "aa", "biomol": "peptide",
        "topology": "linear", "completeness": "complete", "sourcedb": "refseq",
        "updatedate": "2020",
    }}
    session = FakeSession(reg)
    client = make_client(session, tmp_path)
    result = client.efetch_all("protein", ["9"], "fasta")
    assert result.returned == ["YP_009"]
    assert ">YP_009.1" in result.text
