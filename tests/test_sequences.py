"""Tests for the sequence service (local accession resolution)."""

import pytest

from viralfetch import sequences
from viralfetch.ncbi import SeqMeta
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


# -- taxon-level (Phase 4) -------------------------------------------------

def test_accessions_for_taxon_pools_family(vmr):
    name, rank, accs = sequences.accessions_for_taxon(vmr, "Alphaviridae")
    assert name == "Alphaviridae"
    assert rank == "family"
    # All accessions across Alphavirus (4) + Betavirus (1).
    assert set(accs) == {"NC_000001", "NC_000002", "NC_000003", "NC_000004", "NC_000005"}


def test_taxon_aggregate_counts(vmr):
    agg = sequences.taxon_aggregate(vmr, "Alphaviridae")
    assert agg.species == 3  # Alphavirus one/two + Betavirus one
    assert agg.isolates == 4
    assert agg.accessions == 5
    # NC_ accessions are RefSeq (underscore prefix).
    assert agg.refseq == 5
    assert agg.moltype_breakdown["ssRNA(+)"] == 4
    assert agg.moltype_breakdown["dsDNA"] == 1


def _meta(acc, moltype, biomol):
    return SeqMeta(acc, "org", 1, moltype, biomol, "linear", "complete", "insd", "2020")


def test_filter_records_moltype_normalised():
    recs = [_meta("A", "ss-RNA", "genomic"), _meta("B", "dna", "genomic")]
    kept = sequences.filter_records(recs, moltype="ssRNA", biomol=None)
    assert [r.accession for r in kept] == ["A"]  # ss-RNA matches ssRNA


def test_filter_records_biomol():
    recs = [_meta("A", "rna", "genomic"), _meta("B", "rna", "mRNA")]
    kept = sequences.filter_records(recs, moltype=None, biomol="mRNA")
    assert [r.accession for r in kept] == ["B"]
