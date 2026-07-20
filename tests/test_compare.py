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
