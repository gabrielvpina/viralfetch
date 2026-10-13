"""Annotate a list of names with their full taxonomy (``viralfetch annotate``).

One input line = one name (a species or any higher taxon; blank lines and
``#`` comments are skipped). Each becomes one table row with the ICTV lineage
spread over one column per rank, plus what links it to NCBI: the taxid, the
exemplar GenBank/RefSeq accession(s) of a species, and links to both.

Resolution, per name:

1. the local VMR (``status = "ictv"``);
2. else NCBI taxonomy (needs an email; skipped with ``ncbi=None``) — the
   nearest VMR taxon in NCBI's lineage gives the ICTV columns
   (``status = "ncbi"``), or NCBI's own ranks when no VMR taxon is close
   enough (``status = "ncbi_only"``);
3. else ``status = "not_found"``, with "did you mean" suggestions in ``note``.

Repeated names are resolved once. NCBI failures never stop the run: the row
keeps its local data and the error goes to ``note``.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, field

from . import compare
from .models import RANKS, Taxon
from .ncbi import NCBIClient, NCBIError
from .queries import descendant_isolates
from .vmr import VMR

NUCCORE_URL = "https://www.ncbi.nlm.nih.gov/nuccore/{ids}"
TAXONOMY_URL = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}"

COLUMNS = (
    ["query", "status", "name", "rank"]
    + list(RANKS)
    + ["genome", "host_source", "species_count", "ncbi_taxid", "ncbi_taxonomy_url",
       "accession", "accession_url", "note"]
)

Progress = Callable[[int, int, "Annotation"], None]


@dataclass
class Annotation:
    query: str
    status: str                          # ictv | ncbi | ncbi_only | not_found
    name: str | None = None
    rank: str | None = None
    lineage: dict[str, str] = field(default_factory=dict)
    genome: str | None = None
    host_source: str | None = None
    species_count: int | None = None
    ncbi_taxid: str | None = None
    accessions: list[str] = field(default_factory=list)
    note: str | None = None

    def row(self) -> dict[str, str | int | None]:
        """The table row, keyed by :data:`COLUMNS`."""
        row: dict[str, str | int | None] = {
            "query": self.query, "status": self.status,
            "name": self.name, "rank": self.rank,
        }
        row.update({rank: self.lineage.get(rank) for rank in RANKS})
        row.update({
            "genome": self.genome,
            "host_source": self.host_source,
            "species_count": self.species_count,
            "ncbi_taxid": self.ncbi_taxid,
            "ncbi_taxonomy_url": TAXONOMY_URL.format(taxid=self.ncbi_taxid) if self.ncbi_taxid else None,
            "accession": ";".join(self.accessions) or None,
            "accession_url": NUCCORE_URL.format(ids=",".join(self.accessions)) if self.accessions else None,
            "note": self.note,
        })
        return row


def read_names(lines: Iterable[str]) -> list[str]:
    """Names from an input file: one per line, blanks and ``#`` comments skipped."""
    names = []
    for line in lines:
        name = line.split("#", 1)[0].strip().strip('"').strip()
        if name:
            names.append(name)
    return names


def _vmr_fields(vmr: VMR, taxon: Taxon) -> dict:
    _t, isolates = descendant_isolates(vmr, taxon.name)
    genomes = sorted({i.genome_composition for i in isolates if i.genome_composition})
    hosts: dict[str, None] = {}
    for iso in isolates:
        for h in iso.host_source.split(","):
            if h.strip():
                hosts.setdefault(h.strip(), None)
    accessions: list[str] = []
    if taxon.rank == "species":
        for iso in isolates:
            if iso.exemplar:
                accessions.extend(a.accession for a in iso.accessions)
    return {
        "lineage": dict(taxon.lineage),
        "genome": "; ".join(genomes) or None,
        "host_source": "; ".join(hosts) or None,
        "species_count": len({i.species for i in isolates}),
        "accessions": list(dict.fromkeys(accessions)),
    }


def annotate_one(vmr: VMR, name: str, ncbi: NCBIClient | None) -> Annotation:
    taxon = vmr.find(name)
    if taxon is not None:
        ann = Annotation(query=name, status="ictv", name=taxon.name, rank=taxon.rank,
                         **_vmr_fields(vmr, taxon))
        if ncbi is not None:
            try:
                ann.ncbi_taxid = ncbi.esearch_taxid(taxon.name)
            except NCBIError as exc:
                ann.note = f"NCBI taxid lookup failed: {exc}"
        return ann

    lineage = None
    if ncbi is not None:
        try:
            lineage = compare.lineage_via_ncbi(ncbi, name)
        except NCBIError as exc:
            return Annotation(query=name, status="not_found",
                              note=f"not in the VMR; NCBI lookup failed: {exc}")
    if lineage is None:
        suggestions = vmr.suggest(name)
        note = "not in the VMR" + ("" if ncbi is not None else " (NCBI not consulted)")
        if suggestions:
            note += "; did you mean: " + ", ".join(suggestions)
        return Annotation(query=name, status="not_found", note=note)

    nearest = compare.nearest_vmr_taxon(vmr, lineage.lineage)
    if nearest is not None:
        ann = Annotation(query=name, status="ncbi", name=nearest.name, rank=nearest.rank,
                         **_vmr_fields(vmr, nearest))
        ann.ncbi_taxid = lineage.taxid
        ann.note = (f"NCBI: {lineage.name} ({lineage.rank or 'no rank'}, taxid {lineage.taxid}); "
                    f"ICTV columns from its {nearest.rank} {nearest.name}")
        return ann
    return Annotation(
        query=name, status="ncbi_only", name=lineage.name, rank=lineage.rank or "no rank",
        lineage={r: n for r, n in lineage.lineage if r in RANKS},
        ncbi_taxid=lineage.taxid,
        note="NCBI only: no ICTV taxon down to family; lineage columns are NCBI's",
    )


def annotate(
    vmr: VMR,
    names: list[str],
    *,
    ncbi: NCBIClient | None = None,
    progress: Progress | None = None,
) -> list[Annotation]:
    """One :class:`Annotation` per name, in input order."""
    cache: dict[str, Annotation] = {}
    out: list[Annotation] = []
    for i, name in enumerate(names, 1):
        key = name.casefold()
        if key not in cache:
            cache[key] = annotate_one(vmr, name, ncbi)
        ann = cache[key]
        if ann.query != name:  # same name, other spelling: keep the input as typed
            ann = Annotation(**{**ann.__dict__, "query": name})
        out.append(ann)
        if progress is not None:
            progress(i, len(names), ann)
    return out
