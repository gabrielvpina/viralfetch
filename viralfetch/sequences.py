"""Sequence service: resolve a species to accessions (local VMR), then fetch
metadata or records from NCBI.

Like the taxonomy services, these functions return plain result objects and
never print. Accession resolution is purely local; only the fetch touches the
network.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field

from .ncbi import MetaResult, NCBIClient, RecordsResult, SeqMeta
from .queries import TaxonNotFound, descendant_isolates
from .vmr import VMR


def _is_refseq(accession: str) -> bool:
    """RefSeq accessions carry a ``XX_`` prefix (e.g. NC_, NP_, YP_)."""
    return "_" in accession


def _dedupe(accessions):
    seen: set[str] = set()
    out: list[str] = []
    for acc in accessions:
        if acc not in seen:
            seen.add(acc)
            out.append(acc)
    return out


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


def accessions_for_taxon(vmr: VMR, name: str) -> tuple[str, str, list[str]]:
    """Return ``(canonical_name, rank, accessions)`` for any taxon (inclusive).

    Unlike :func:`accessions_for_species`, this accepts any rank and pools the
    accessions of every isolate beneath it.
    """
    taxon, isolates = descendant_isolates(vmr, name)
    accs = [a.accession for iso in isolates for a in iso.accessions]
    return taxon.name, taxon.rank, _dedupe(accs)


@dataclass
class TaxonAggregate:
    name: str
    rank: str
    species: int
    isolates: int
    accessions: int
    moltype_breakdown: dict[str, int] = field(default_factory=dict)
    refseq: int = 0


def taxon_aggregate(vmr: VMR, name: str) -> TaxonAggregate:
    """A purely-local summary of what a taxon contains, to gauge a download.

    ``moltype_breakdown`` uses the VMR genome-composition field (local); a true
    per-record moltype needs esummary and is not worth fetching just to decide.
    """
    taxon, isolates = descendant_isolates(vmr, name)
    species = {iso.species for iso in isolates if iso.species}
    compositions: Counter[str] = Counter()
    refseq = 0
    total_acc = 0
    for iso in isolates:
        for acc in iso.accessions:
            total_acc += 1
            compositions[iso.genome_composition or "(unknown)"] += 1
            if _is_refseq(acc.accession):
                refseq += 1
    return TaxonAggregate(
        name=taxon.name,
        rank=taxon.rank,
        species=len(species),
        isolates=len(isolates),
        accessions=total_acc,
        moltype_breakdown=dict(compositions),
        refseq=refseq,
    )


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def filter_records(
    records: list[SeqMeta], moltype: str | None, biomol: str | None
) -> list[SeqMeta]:
    """Filter esummary records by moltype/biomol (SPEC section 5.6).

    Matching is normalised and substring-based, so ``--moltype ssRNA`` matches
    ``ss-RNA``. These are nuccore fields; the filter runs locally over the
    esummary result — it is never a query-side parameter.
    """
    mt = _norm(moltype) if moltype else None
    bm = _norm(biomol) if biomol else None
    out = []
    for r in records:
        if mt and mt not in _norm(r.moltype):
            continue
        if bm and bm not in _norm(r.biomol):
            continue
        out.append(r)
    return out


# -- unified fetch entry points -------------------------------------------

def resolve_accessions(vmr: VMR, *, species: str | None, taxon: str | None) -> tuple[str, list[str]]:
    """Resolve the target (species or taxon) to a label and accession list."""
    if taxon is not None:
        name, _rank, accs = accessions_for_taxon(vmr, taxon)
        return name, accs
    name, accs = accessions_for_species(vmr, species)
    return name, accs


def meta(
    ncbi: NCBIClient,
    accessions: list[str],
    *,
    moltype: str | None = None,
    biomol: str | None = None,
) -> MetaResult:
    """nuccore metadata, optionally filtered by moltype/biomol."""
    result = ncbi.esummary_nuccore(accessions)
    if moltype or biomol:
        result.records = filter_records(result.records, moltype, biomol)
    return result


def records(
    ncbi: NCBIClient,
    accessions: list[str],
    rettype: str,
    *,
    moltype: str | None = None,
    biomol: str | None = None,
) -> RecordsResult:
    """nuccore fasta/gb records. Filters require an esummary pass first."""
    if moltype or biomol:
        kept = filter_records(ncbi.esummary_nuccore(accessions).records, moltype, biomol)
        accessions = [r.accession for r in kept]
    return ncbi.efetch_nuccore(accessions, rettype)


def protein_records(ncbi: NCBIClient, accessions: list[str], rettype: str) -> RecordsResult:
    """Protein path: elink nuccore -> protein, then efetch db=protein.

    Protein UIDs come from elink; efetch returns records keyed by protein
    accession, so we fetch them all rather than matching per requested UID.
    """
    uids = ncbi.protein_uids_for(accessions)
    return ncbi.efetch_all("protein", uids, rettype)


def protein_meta(ncbi: NCBIClient, accessions: list[str]) -> MetaResult:
    uids = ncbi.protein_uids_for(accessions)
    return ncbi.esummary_all("protein", uids)


# Backwards-compatible species-only helpers (Phase 3).

def seq_meta(vmr: VMR, ncbi: NCBIClient, name: str) -> tuple[str, MetaResult]:
    species, accessions = accessions_for_species(vmr, name)
    return species, ncbi.esummary_nuccore(accessions)


def seq_records(
    vmr: VMR, ncbi: NCBIClient, name: str, rettype: str
) -> tuple[str, RecordsResult]:
    species, accessions = accessions_for_species(vmr, name)
    return species, ncbi.efetch_nuccore(accessions, rettype)
