"""Sequence service: resolve a species to accessions (local VMR), then fetch
metadata or records from NCBI.

Like the taxonomy services, these functions return plain result objects and
never print. Accession resolution is purely local; only the fetch touches the
network.
"""

from __future__ import annotations

from .ncbi import MetaResult, NCBIClient, RecordsResult
from .queries import TaxonNotFound
from .vmr import VMR


def accessions_for_species(vmr: VMR, name: str) -> tuple[str, list[str]]:
    """Return ``(canonical_species_name, accessions)`` for a species.

    Raises :class:`TaxonNotFound` (with suggestions) if the name is unknown, or
    :class:`NotASpecies` if it resolves to a higher-rank taxon.
    """
    taxon = vmr.find(name)
    if taxon is None:
        raise TaxonNotFound(name, vmr.suggest(name))
    if taxon.rank != "species":
        raise NotASpecies(taxon.name, taxon.rank)

    isolates = vmr.isolates_by_species.get(taxon.name.casefold(), [])
    seen: set[str] = set()
    accessions: list[str] = []
    for iso in isolates:
        for acc in iso.accessions:
            if acc.accession not in seen:
                seen.add(acc.accession)
                accessions.append(acc.accession)
    return taxon.name, accessions


class NotASpecies(Exception):
    """Raised when ``seq <name>`` gets a non-species taxon (use --taxon)."""

    def __init__(self, name: str, rank: str):
        self.name = name
        self.rank = rank
        super().__init__(f"{name!r} is a {rank}, not a species")


def seq_meta(vmr: VMR, ncbi: NCBIClient, name: str) -> tuple[str, MetaResult]:
    species, accessions = accessions_for_species(vmr, name)
    return species, ncbi.esummary_nuccore(accessions)


def seq_records(
    vmr: VMR, ncbi: NCBIClient, name: str, rettype: str
) -> tuple[str, RecordsResult]:
    species, accessions = accessions_for_species(vmr, name)
    return species, ncbi.efetch_nuccore(accessions, rettype)
