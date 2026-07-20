"""Tests for the sequence service (local accession resolution)."""

import pytest

from viralfetch import sequences
from viralfetch.queries import TaxonNotFound


def test_accessions_for_species(vmr):
    name, accs = sequences.accessions_for_species(vmr, "Alphavirus two")
    assert name == "Alphavirus two"
    # Segmented exemplar contributes both segment accessions.
    assert accs == ["NC_000003", "NC_000004"]


def test_accessions_dedupe_and_order(vmr):
    name, accs = sequences.accessions_for_species(vmr, "Alphavirus one")
    assert accs == ["NC_000001", "NC_000002"]


def test_accessions_case_insensitive(vmr):
    name, _ = sequences.accessions_for_species(vmr, "alphavirus one")
    assert name == "Alphavirus one"


def test_species_not_found_suggests(vmr):
    with pytest.raises(TaxonNotFound):
        sequences.accessions_for_species(vmr, "Nope virus")


def test_higher_taxon_rejected(vmr):
    with pytest.raises(sequences.NotASpecies) as exc:
        sequences.accessions_for_species(vmr, "Alphaviridae")
    assert exc.value.rank == "family"
