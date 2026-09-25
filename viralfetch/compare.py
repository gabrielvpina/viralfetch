"""ICTV-vs-NCBI lineage comparison (SPEC section 5.2).

Flow: take the local ICTV lineage from the VMR, find a representative accession
for the taxon, elink nuccore -> taxonomy for the NCBI taxid, efetch the NCBI
lineage, and return both side by side. Divergences are the *product* of this
command, not an error — NCBI commonly lags ICTV.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .models import RANKS, Taxon
from .ncbi import NcbiLineage, NCBIClient
from .queries import TaxonNotFound, descendant_isolates
from .vmr import VMR


class NoRepresentativeAccession(Exception):
    """Raised when a taxon has no accession to anchor the NCBI lookup."""

    def __init__(self, name: str):
        self.name = name
        super().__init__(f"no accession available for {name!r}")


class NoNcbiTaxid(Exception):
    """Raised when NCBI has no taxonomy link for the representative accession."""

    def __init__(self, accession: str):
        self.accession = accession
        super().__init__(f"NCBI returned no taxid for {accession!r}")


class NcbiTaxonNotFound(Exception):
    """Raised when NCBI's taxonomy database has no match for a queried name."""

    def __init__(self, name: str):
        self.name = name
        super().__init__(f"NCBI taxonomy has no match for {name!r}")


@dataclass
class LineageComparison:
    taxon: str
    representative_accession: str
    ictv: list[tuple[str, str]] = field(default_factory=list)  # (rank, name)
    ncbi: NcbiLineage | None = None


def _representative_accession(vmr: VMR, name: str) -> str:
    _taxon, isolates = descendant_isolates(vmr, name)
    # Prefer an exemplar isolate's first accession; fall back to any.
    for iso in sorted(isolates, key=lambda i: not i.exemplar):
        if iso.accessions:
            return iso.accessions[0].accession
    raise NoRepresentativeAccession(name)


def compare_ncbi(vmr: VMR, ncbi: NCBIClient, name: str) -> LineageComparison:
    taxon = vmr.find(name)
    if taxon is None:
        raise TaxonNotFound(name, vmr.suggest(name))

    accession = _representative_accession(vmr, name)
    taxid = ncbi.taxid_for_accession(accession)
    if not taxid:
        raise NoNcbiTaxid(accession)

    ncbi_lineage = ncbi.efetch_taxonomy(taxid)
    ictv = [(rank, taxon.lineage[rank]) for rank in RANKS if rank in taxon.lineage]
    return LineageComparison(
        taxon=taxon.name,
        representative_accession=accession,
        ictv=ictv,
        ncbi=ncbi_lineage,
    )


def tax_ncbi(ncbi: NCBIClient, name: str) -> NcbiLineage:
    """Look a taxon's lineage up directly in NCBI's taxonomy (no VMR).

    Resolves ``name`` to a taxid with esearch, then efetches its lineage.
    Raises :class:`NcbiTaxonNotFound` when NCBI knows no such taxon.
    """
    taxid = ncbi.esearch_taxid(name)
    if not taxid:
        raise NcbiTaxonNotFound(name)
    return ncbi.efetch_taxonomy(taxid)


def lineage_via_ncbi(ncbi: NCBIClient, name: str) -> NcbiLineage | None:
    """The NCBI lineage of ``name``, or ``None`` if NCBI knows no such taxon.

    The soft-failing sibling of :func:`tax_ncbi`, used as a fallback for names
    the local VMR does not know.
    """
    taxid = ncbi.esearch_taxid(name)
    if not taxid:
        return None
    return ncbi.efetch_taxonomy(taxid)


def family_via_ncbi(ncbi: NCBIClient, name: str) -> str | None:
    """Resolve the family of ``name`` through NCBI taxonomy, or ``None``.

    A fallback for names the local VMR does not know: search NCBI for the name
    and return the ``family`` rank of its lineage (virus family names are shared
    between NCBI and ICTV). Returns ``None`` if NCBI knows no such taxon or its
    lineage has no family rank.
    """
    lineage = lineage_via_ncbi(ncbi, name)
    if lineage is None:
        return None
    for rank, taxon_name in lineage.lineage:
        if rank == "family":
            return taxon_name
    return None


def nearest_vmr_taxon(
    vmr: VMR, lineage: list[tuple[str, str]], broadest: str = "family"
) -> Taxon | None:
    """The deepest taxon of an NCBI lineage that the VMR knows, or ``None``.

    Used to redirect a name the VMR doesn't know (an NCBI-only name, an old
    synonym, a strain) to the closest VMR taxon. Matches broader than
    ``broadest`` are ignored — "members of Riboviria" is rarely what was meant.
    """
    limit = RANKS.index(broadest)
    for _rank, name in reversed(lineage):
        taxon = vmr.find(name)
        if taxon is not None:
            return taxon if RANKS.index(taxon.rank) >= limit else None
    return None
