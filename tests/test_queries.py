"""Tests for the taxonomy query service layer."""

import pytest

from viralfetch import queries
from viralfetch.queries import InvalidRank, TaxonNotFound


def test_vmr_indices(vmr):
    assert len(vmr.species) == 4
    assert set(vmr.species) == {
        "Alphavirus one", "Alphavirus two", "Betavirus one", "Gammavirus one"
    }
    assert vmr.empty_accession_rows == 1
    assert vmr.unparsed_rows == []


def test_tax_family_lineage(vmr):
    view = queries.tax(vmr, "Alphaviridae")
    assert view.taxon.rank == "family"
    assert view.taxon.lineage["realm"] == "Testviria"
    assert view.taxon.lineage["family"] == "Alphaviridae"
    assert "genus" not in view.taxon.lineage  # family node stops at family
    assert view.isolate_summary is None


def test_tax_case_insensitive(vmr):
    assert queries.tax(vmr, "alphaviridae").taxon.name == "Alphaviridae"


def test_tax_species_isolate_summary(vmr):
    view = queries.tax(vmr, "Alphavirus one")
    s = view.isolate_summary
    assert s is not None
    assert s.total == 2
    assert s.exemplars == 1
    assert s.additional == 1
    assert s.accessions == 2
    assert s.genome_compositions == {"ssRNA(+)": 2}


def test_tax_species_segmented_accession_count(vmr):
    # "Alphavirus two" exemplar has 2 segment accessions -> counted separately.
    view = queries.tax(vmr, "Alphavirus two")
    assert view.isolate_summary.accessions == 2


def test_tax_not_found_suggests(vmr):
    with pytest.raises(TaxonNotFound) as exc:
        queries.tax(vmr, "AlphaviridaX")
    assert "Alphaviridae" in exc.value.suggestions


def test_members_at_rank_with_species_counts(vmr):
    view = queries.members(vmr, "Alphaviridae", rank="genus")
    names = {m.name: m.species_count for m in view.members}
    assert names == {"Alphavirus": 2, "Betavirus": 1}


def test_members_family_under_order(vmr):
    view = queries.members(vmr, "Testord", rank="family")
    assert {m.name for m in view.members} == {"Alphaviridae", "Gammaviridae"}


def test_members_count_only(vmr):
    view = queries.members(vmr, "Alphaviridae", rank="genus", count=True)
    assert view.count_only is True
    assert view.breakdown == {"genus": 2}


def test_members_breakdown_no_rank(vmr):
    view = queries.members(vmr, "Testord")
    assert view.rank is None
    assert view.breakdown["family"] == 2
    assert view.breakdown["genus"] == 3
    assert view.breakdown["species"] == 4


def test_members_invalid_rank(vmr):
    # "realm" is above family -> invalid as a member rank of a family.
    with pytest.raises(InvalidRank):
        queries.members(vmr, "Alphaviridae", rank="realm")


def test_members_not_found(vmr):
    with pytest.raises(TaxonNotFound):
        queries.members(vmr, "Nope")
