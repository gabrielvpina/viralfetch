"""Parser tests against the hard cases from SPEC section 6 plus real VMR shapes."""

from viralfetch.accessions import is_parseable, parse_accessions
from viralfetch.models import Accession


def test_simple():
    assert parse_accessions("NC_045512") == [Accession("NC_045512", None)]


def test_simple_wgs_long_prefix():
    # WGS accessions have a 4-6 letter prefix; must not be rejected.
    assert parse_accessions("CAJDJZ010000002") == [
        Accession("CAJDJZ010000002", None)
    ]


def test_segmented_labelled_rna():
    assert parse_accessions("RNA1: NC_003615; RNA2: NC_003616") == [
        Accession("NC_003615", "RNA1"),
        Accession("NC_003616", "RNA2"),
    ]


def test_geminivirus_dna_segments():
    assert parse_accessions("DNA-A: X15656; DNA-B: X15657") == [
        Accession("X15656", "DNA-A"),
        Accession("X15657", "DNA-B"),
    ]


def test_multiple_unlabelled():
    assert parse_accessions("KT732816; KT732815; KT732817") == [
        Accession("KT732816", None),
        Accession("KT732815", None),
        Accession("KT732817", None),
    ]


def test_coordinate_annotation_stripped():
    assert parse_accessions("AE006468 (2844298.2877981)") == [
        Accession("AE006468", None)
    ]


def test_textual_annotation_stripped():
    assert parse_accessions("AB012345 (partial)") == [Accession("AB012345", None)]


def test_partial_as_label():
    # Real VMR shape: "partial:" used as a segment label.
    assert parse_accessions("partial: KF360970; partial: KF360971") == [
        Accession("KF360970", "partial"),
        Accession("KF360971", "partial"),
    ]


def test_versioned_accession_kept():
    assert parse_accessions("MN908947.3") == [Accession("MN908947.3", None)]


def test_empty_and_none():
    assert parse_accessions("") == []
    assert parse_accessions(None) == []
    assert parse_accessions("   ") == []


def test_whitespace_between_segments():
    assert parse_accessions("A: EU623082;   B: EU623083") == [
        Accession("EU623082", "A"),
        Accession("EU623083", "B"),
    ]


def test_is_parseable():
    assert is_parseable("NC_045512") is True
    assert is_parseable("") is False
    assert is_parseable(None) is False
    assert is_parseable("not-an-accession") is False
