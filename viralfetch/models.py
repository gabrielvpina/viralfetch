"""Core data structures shared across viralfetch.

These are plain, presentation-free dataclasses. Business logic returns these
(or dicts of them); only the ``render`` layer turns them into output.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# The ICTV taxonomic ranks, ordered from most inclusive to least. These map
# directly onto columns of the VMR TSV.
RANKS: tuple[str, ...] = (
    "realm",
    "subrealm",
    "kingdom",
    "subkingdom",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus",
    "subgenus",
    "species",
)

# English plurals for rank labels (used only for display text).
RANK_PLURAL: dict[str, str] = {
    "realm": "realms",
    "subrealm": "subrealms",
    "kingdom": "kingdoms",
    "subkingdom": "subkingdoms",
    "phylum": "phyla",
    "subphylum": "subphyla",
    "class": "classes",
    "subclass": "subclasses",
    "order": "orders",
    "suborder": "suborders",
    "family": "families",
    "subfamily": "subfamilies",
    "genus": "genera",
    "subgenus": "subgenera",
    "species": "species",
}


def plural(rank: str) -> str:
    """Human-readable plural of a rank label."""
    return RANK_PLURAL.get(rank, rank + "s")


@dataclass(frozen=True)
class Accession:
    """A single GenBank/RefSeq accession parsed from the VMR free-text field.

    ``segment`` is the segment/molecule label when the source virus is
    segmented (e.g. ``"RNA1"``, ``"DNA-A"``), otherwise ``None``.
    """

    accession: str
    segment: str | None = None


@dataclass
class Isolate:
    """An exemplar or additional isolate of a species, as listed in the VMR."""

    isolate_id: str
    species: str
    exemplar: bool  # True for the exemplar ("E"), False for additional ("A")
    virus_names: str
    abbreviations: str
    designation: str
    genome_composition: str
    host_source: str
    raw_accession: str
    accessions: list[Accession] = field(default_factory=list)


@dataclass
class Taxon:
    """A node in the ICTV taxonomy, with its full lineage.

    ``lineage`` maps rank name -> taxon name for every populated rank at or
    above this node (missing ranks are simply absent from the dict).
    """

    name: str
    rank: str
    lineage: dict[str, str] = field(default_factory=dict)

    def parent_rank_value(self, rank: str) -> str | None:
        return self.lineage.get(rank)


@dataclass
class Sequence:
    """A biological sequence fetched from NCBI (nt or aa)."""

    accession: str
    organism: str
    length: int
    moltype: str
    definition: str
    data: str  # FASTA/GenBank text, or empty when only metadata was fetched


@dataclass
class Chapter:
    """A parsed ICTV Report chapter rendered to Markdown."""

    slug: str
    title: str
    markdown: str
    authors: str | None = None
    citation: str | None = None
    doi: str | None = None
    url: str | None = None
