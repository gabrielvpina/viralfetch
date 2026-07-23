"""Tests for the ICTV-vs-NCBI lineage comparison service."""

import pytest

from viralfetch import compare
from viralfetch.ncbi import NcbiLineage
from viralfetch.queries import TaxonNotFound


class StubNcbi:
    """Minimal stand-in exposing only what compare_ncbi calls."""

    def __init__(self, taxid, lineage):
        self._taxid = taxid
        self._lineage = lineage
        self.seen_accession = None

    def taxid_for_accession(self, accession):
        self.seen_accession = accession
        return self._taxid

    def efetch_taxonomy(self, taxid):
        return self._lineage


class StubSearchNcbi:
    """Stand-in for the name->taxid->lineage path (tax_ncbi/family_via_ncbi)."""

    def __init__(self, taxid, lineage):
        self._taxid = taxid
        self._lineage = lineage
        self.seen_name = None

    def esearch_taxid(self, name):
        self.seen_name = name
        return self._taxid

    def efetch_taxonomy(self, taxid):
        return self._lineage


def test_compare_builds_both_lineages(vmr):
    ncbi_lineage = NcbiLineage(
        taxid="999", name="Alphavirus one", rank="species",
        lineage=[("realm", "Testviria"), ("genus", "Alphavirus"), ("species", "Alphavirus one")],
    )
    stub = StubNcbi("999", ncbi_lineage)
    result = compare.compare_ncbi(vmr, stub, "Alphavirus one")

    # Representative accession is the exemplar's first accession.
    assert result.representative_accession == "NC_000001"
    assert stub.seen_accession == "NC_000001"
    # ICTV lineage comes from the VMR taxon.
    ictv = dict(result.ictv)
    assert ictv["family"] == "Alphaviridae"
    assert result.ncbi.taxid == "999"


def test_compare_not_found(vmr):
    stub = StubNcbi("1", None)
    with pytest.raises(TaxonNotFound):
        compare.compare_ncbi(vmr, stub, "Nope")


def test_compare_no_taxid(vmr):
    stub = StubNcbi(None, None)  # NCBI returns no taxid
    with pytest.raises(compare.NoNcbiTaxid):
        compare.compare_ncbi(vmr, stub, "Alphavirus one")


# -- tax --ncbi (direct NCBI lineage, no VMR) -----------------------------

def _ncbi_lineage():
    return NcbiLineage(
        taxid="2697049", name="SARS-CoV-2", rank="species",
        lineage=[("realm", "Riboviria"), ("family", "Coronaviridae"),
                 ("species", "SARS-CoV-2")],
    )


def test_tax_ncbi_returns_lineage():
    stub = StubSearchNcbi("2697049", _ncbi_lineage())
    lineage = compare.tax_ncbi(stub, "SARS-CoV-2")
    assert stub.seen_name == "SARS-CoV-2"
    assert lineage.taxid == "2697049"
    assert lineage.name == "SARS-CoV-2"


def test_tax_ncbi_not_found():
    stub = StubSearchNcbi(None, None)  # esearch finds nothing
    with pytest.raises(compare.NcbiTaxonNotFound):
        compare.tax_ncbi(stub, "Nope")


# -- family_via_ncbi (fallback for `text`) --------------------------------

def test_family_via_ncbi_extracts_family():
    stub = StubSearchNcbi("2697049", _ncbi_lineage())
    assert compare.family_via_ncbi(stub, "SARS-CoV-2") == "Coronaviridae"


def test_family_via_ncbi_no_taxon_returns_none():
    stub = StubSearchNcbi(None, None)
    assert compare.family_via_ncbi(stub, "Nope") is None


def test_family_via_ncbi_no_family_rank_returns_none():
    lineage = NcbiLineage(
        taxid="10", name="Riboviria", rank="realm",
        lineage=[("realm", "Riboviria")],  # no family rank
    )
    stub = StubSearchNcbi("10", lineage)
    assert compare.family_via_ncbi(stub, "Riboviria") is None
