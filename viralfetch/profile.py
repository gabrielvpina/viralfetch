"""The default command: a one-stop profile of a taxon (``viralfetch <name>``).

It gathers what each source is best at, side by side:

- **ICTV (VMR)** — the official lineage, how many species / isolates sit under
  the taxon, their genome compositions and host sources;
- **NCBI** — the taxid, NCBI's lineage compared rank by rank with ICTV's, the
  ICTV status NCBI records for the name, and how many nucleotide records,
  complete genomes, RefSeq genomes and assemblies it holds (NCBI Datasets);
  above species, the genome sizes and proteins of its RefSeq genomes
  (:mod:`viralfetch.genomics`);
- **ViralZone** — the curated factsheet closest to the taxon (host, geography,
  diseases, transmission, vaccines, antivirals, ...), credited to ViralZone.

The VMR and ViralZone are local; NCBI is optional. Without an email (or with
``local=True``) the profile is built from the local sources alone and says so.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

from . import compare, genomics
from . import trees as trees_mod
from . import viralzone as vz_mod
from .datasets import VIRUS_PORTAL, DatasetsClient, SequenceStats, TaxonReport, sorted_genomes
from .genomics import GenomeSummary
from .models import RANKS, Taxon
from .ncbi import NcbiLineage, NCBIClient, NCBIError
from .queries import TaxonNotFound, descendant_isolates
from .vmr import VMR

TOP_HOSTS = 5
REFSEQ_LISTED = 25


@dataclass
class RankRow:
    """One rank of the ICTV-vs-NCBI comparison. ``agree`` is ``None`` when only
    one side has the rank."""

    rank: str
    ictv: str | None
    ncbi: str | None
    agree: bool | None


@dataclass
class VmrSummary:
    species: int
    isolates: int
    exemplars: int
    genome_compositions: dict[str, int] = field(default_factory=dict)
    host_sources: dict[str, int] = field(default_factory=dict)  # top TOP_HOSTS
    exemplar_accessions: list[str] = field(default_factory=list)  # species only


@dataclass
class ViralZoneSection:
    factsheet: dict
    related: list[dict] = field(default_factory=list)  # other matches, broader
    host_groups: list[str] = field(default_factory=list)
    note: str | None = None  # e.g. "no Orthoflavivirus denguei factsheet; ..."
    credit: dict = field(default_factory=lambda: dict(vz_mod.CREDIT))


@dataclass
class Profile:
    query: str
    name: str                      # the taxon profiled (VMR name when known)
    rank: str
    ictv: Taxon | None = None
    ncbi: NcbiLineage | None = None
    ncbi_report: TaxonReport | None = None
    comparison: list[RankRow] = field(default_factory=list)
    vmr: VmrSummary | None = None
    sequences: SequenceStats | None = None
    genomes: GenomeSummary | None = None   # above species: RefSeq genome features
    viralzone: ViralZoneSection | None = None
    family: str | None = None
    has_trees: bool = False
    notes: list[str] = field(default_factory=list)     # how the name was resolved
    warnings: list[str] = field(default_factory=list)  # sources that were skipped / failed


# -- pieces ---------------------------------------------------------------------

def compare_lineages(ictv: Taxon | None, ncbi: NcbiLineage | None) -> list[RankRow]:
    """ICTV and NCBI lineages aligned on the ICTV ranks, root first."""
    ictv_ranks = ictv.lineage if ictv else {}
    ncbi_ranks: dict[str, str] = {}
    for rank, name in (ncbi.lineage if ncbi else []):
        if rank in RANKS:
            ncbi_ranks.setdefault(rank, name)
    rows: list[RankRow] = []
    for rank in RANKS:
        a, b = ictv_ranks.get(rank), ncbi_ranks.get(rank)
        if a is None and b is None:
            continue
        agree = None if a is None or b is None else a.casefold() == b.casefold()
        rows.append(RankRow(rank=rank, ictv=a, ncbi=b, agree=agree))
    return rows


def vmr_summary(vmr: VMR, taxon: Taxon) -> VmrSummary:
    _taxon, isolates = descendant_isolates(vmr, taxon.name)
    compositions = Counter(i.genome_composition for i in isolates if i.genome_composition)
    hosts: Counter[str] = Counter()
    for iso in isolates:
        for host in (h.strip() for h in iso.host_source.split(";")):
            if host:
                hosts[host] += 1
    exemplar_accessions: list[str] = []
    if taxon.rank == "species":
        for iso in isolates:
            if iso.exemplar:
                exemplar_accessions.extend(
                    f"{a.segment}: {a.accession}" if a.segment else a.accession
                    for a in iso.accessions
                )
    return VmrSummary(
        species=len({i.species for i in isolates}),
        isolates=len(isolates),
        exemplars=sum(1 for i in isolates if i.exemplar),
        genome_compositions=dict(compositions.most_common()),
        host_sources=dict(hosts.most_common(TOP_HOSTS)),
        exemplar_accessions=exemplar_accessions,
    )


def viralzone_section(
    vz: vz_mod.ViralZone,
    *,
    query: str,
    taxon: Taxon | None,
    ncbi: NcbiLineage | None,
) -> ViralZoneSection | None:
    """The most specific ViralZone factsheet for the taxon, plus broader ones."""
    names = [query]
    if taxon is not None:
        names.append(taxon.name)
    if ncbi is not None:
        names.append(ncbi.name)
        names.extend(name for _rank, name in reversed(ncbi.lineage))
    matches = vz.match(
        taxids=[ncbi.taxid] if ncbi else [],
        names=names,
        lineage=taxon.lineage if taxon else None,
    )
    # Keep factsheets on the taxon's own lineage (or listing it as a species):
    # a name shared with an unrelated clade must not pull in its factsheet.
    # Their own lineage is not used for this — it may lag or carry typos.
    allowed = {n.casefold() for n in names}
    if taxon is not None:
        allowed.update(v.casefold() for v in taxon.lineage.values())
    matches = [
        m for m in matches
        if m["name"].casefold() in allowed
        or (ncbi is not None and str(m.get("taxid")) == str(ncbi.taxid))
        or any(sp["name"].casefold() in allowed for sp in m.get("species", []))
    ]
    if not matches:
        return None
    primary, related = matches[0], matches[1:]

    subject = taxon.name if taxon else query
    note = None
    if primary["name"].casefold() not in {subject.casefold(), query.casefold()} and not (
        ncbi and str(primary.get("taxid")) == str(ncbi.taxid)
    ):
        note = (f"ViralZone has no {subject} factsheet; showing the closest one, "
                f"{primary['name']} ({primary['rank']}).")

    host_groups: set[str] = set()
    for name in [subject, primary["name"]]:
        t = vz.taxon(name)
        if t is not None:
            host_groups.update(t["hosts"])
    if not host_groups and primary["rank"] in {"family", "subfamily"}:
        key = "family" if primary["rank"] == "family" else "subfamily"
        for t in vz.taxa:
            if t["rank"] == "genus" and (t.get(key) or "").casefold() == primary["name"].casefold():
                host_groups.update(t["hosts"])
    return ViralZoneSection(
        factsheet=primary, related=related, host_groups=sorted(host_groups), note=note,
    )


def ncbi_stats(
    ncbi: NCBIClient,
    datasets: DatasetsClient | None,
    taxid: str,
    warnings: list[str],
    *,
    genomes: bool = False,
) -> tuple[TaxonReport | None, SequenceStats | None, GenomeSummary | None]:
    """Record counts (E-utilities), the RefSeq listing / taxonomy report
    (Datasets) and, with ``genomes``, the genome summary of the RefSeq
    sequences for ``taxid`` — the queries run concurrently.

    Each API keeps its own queue, so its rate limit holds. A failed query
    leaves its field empty and adds a warning; nothing here raises.
    """
    def counts() -> tuple[int | None, int | None]:
        total = complete = None
        try:
            total = ncbi.nuccore_count(f"txid{taxid}[Organism:exp]")
            complete = ncbi.nuccore_count(f'txid{taxid}[Organism:exp] AND "complete genome"[Title]')
        except NCBIError as exc:
            warnings.append(f"NCBI record count failed: {exc}")
        return total, complete

    # With the genome summary, the full RefSeq dataset report replaces the
    # short listing (the listing is its first records).
    def datasets_part():
        report = refseq = records = None
        if datasets is None:
            return report, refseq, records
        try:
            report = datasets.taxonomy_report(taxid)
            if genomes:
                records = datasets.refseq_reports(taxid)
                refseq = records[0], sorted_genomes(records[1])[:REFSEQ_LISTED]
            else:
                refseq = datasets.refseq_genomes(taxid)
        except NCBIError as exc:
            warnings.append(f"NCBI Datasets request failed: {exc}")
        return report, refseq, records

    def annotations():
        try:
            return datasets.refseq_reports(taxid, "annotation_report")[1]
        except NCBIError as exc:
            warnings.append(f"NCBI Datasets annotation request failed: {exc}")
            return None

    with ThreadPoolExecutor(max_workers=3) as pool:
        counts_future = pool.submit(counts)
        datasets_future = pool.submit(datasets_part)
        annotation_future = pool.submit(annotations) if genomes and datasets else None
        total, complete = counts_future.result()
        report, refseq, records = datasets_future.result()
        annotated = annotation_future.result() if annotation_future else None

    summary = None
    if records is not None:
        summary = genomics.summarize(taxid, records[1], annotated or [], total=records[0])
    if summary is not None:  # topology and genetic code: the records' nucleotide summaries
        try:
            meta = ncbi.esummary_nuccore(genomics.metadata_accessions(summary))
            genomics.apply_metadata(summary, meta.records)
        except NCBIError as exc:
            warnings.append(f"NCBI topology / genetic code lookup failed: {exc}")

    if total is None and complete is None and refseq is None:
        return report, None, summary
    refseq_total, listed = refseq if refseq is not None else (None, [])
    stats = SequenceStats(
        taxid=int(taxid), total=total, complete=complete, refseq=refseq_total,
        refseq_genomes=listed, portal_url=VIRUS_PORTAL.format(taxid=taxid),
    )
    return report, stats, summary


def _ncbi_rank(ncbi: NcbiLineage | None, rank: str) -> str | None:
    for r, name in (ncbi.lineage if ncbi else []):
        if r == rank:
            return name
    return None


def _has_trees(family: str | None) -> bool:
    if not family:
        return False
    for entry in trees_mod.index().get("families", []):
        if entry.get("family", "").casefold() == family.casefold():
            return entry.get("status") == "included"
    return False


# -- the profile --------------------------------------------------------------------

def build(
    vmr: VMR,
    query: str,
    *,
    ncbi: NCBIClient | None = None,
    datasets: DatasetsClient | None = None,
    viralzone: vz_mod.ViralZone | None = None,
) -> Profile:
    """Profile ``query`` from every available source.

    ``ncbi`` / ``datasets`` are optional: pass ``None`` to stay local. NCBI
    failures are recorded in ``warnings`` rather than raised, so a flaky network
    still yields the local half. Raises :class:`TaxonNotFound` when neither the
    VMR nor NCBI knows the name.
    """
    taxon = vmr.find(query)
    notes: list[str] = []
    warnings: list[str] = []
    ncbi_lineage: NcbiLineage | None = None

    if ncbi is not None:
        try:
            ncbi_lineage, note = _ncbi_lineage(ncbi, vmr, query, taxon)
            if note:
                notes.append(note)
        except NCBIError as exc:
            warnings.append(f"NCBI taxonomy lookup failed: {exc}")

    in_vmr = taxon is not None
    if taxon is None:
        if ncbi_lineage is None:
            raise TaxonNotFound(query, vmr.suggest(query))
        nearest = compare.nearest_vmr_taxon(vmr, ncbi_lineage.lineage)
        if nearest is not None and nearest.name.casefold() != ncbi_lineage.name.casefold():
            notes.append(
                f"{query!r} is not in the local VMR; NCBI places it in {nearest.rank} "
                f"{nearest.name}, whose ICTV lineage is shown."
            )
        elif nearest is None:
            notes.append(f"{query!r} is not in the local VMR; only NCBI knows it.")
        taxon = nearest

    # The subject is what was asked for: the VMR taxon, or NCBI's record of a
    # name the VMR doesn't know (whose nearest VMR taxon still gives the lineage).
    subject_is_vmr = in_vmr or ncbi_lineage is None
    profile = Profile(
        query=query,
        name=taxon.name if subject_is_vmr else ncbi_lineage.name,
        rank=taxon.rank if subject_is_vmr else (ncbi_lineage.rank or "no rank"),
        ictv=taxon,
        ncbi=ncbi_lineage,
        notes=notes,
        warnings=warnings,
    )
    profile.comparison = compare_lineages(taxon, ncbi_lineage)
    if taxon is not None:
        profile.vmr = vmr_summary(vmr, taxon)
    profile.family = (taxon.lineage.get("family") if taxon else None) or _ncbi_rank(ncbi_lineage, "family")
    profile.has_trees = _has_trees(profile.family)

    if ncbi is not None and ncbi_lineage is not None:
        profile.ncbi_report, profile.sequences, profile.genomes = ncbi_stats(
            ncbi, datasets, ncbi_lineage.taxid, warnings,
            genomes=profile.rank in RANKS and profile.rank != "species",
        )

    if viralzone is not None:
        profile.viralzone = viralzone_section(
            viralzone, query=query, taxon=taxon, ncbi=ncbi_lineage,
        )
    return profile


def _ncbi_lineage(
    ncbi: NCBIClient, vmr: VMR, query: str, taxon: Taxon | None,
) -> tuple[NcbiLineage | None, str | None]:
    """``(NCBI record, note)`` for the taxon.

    By name first (ICTV and NCBI share most names). A VMR species NCBI files
    under another name is reached through its exemplar accession; that gives
    the isolate's own NCBI record, which is reported as such. Higher ranks get
    no such fallback — an isolate's counts would misrepresent a genus.
    """
    lineage = compare.lineage_via_ncbi(ncbi, query)
    if lineage is not None or taxon is None:
        return lineage, None
    if taxon.name.casefold() != query.casefold():
        lineage = compare.lineage_via_ncbi(ncbi, taxon.name)
        if lineage is not None:
            return lineage, None
    if taxon.rank != "species":
        return None, f"NCBI taxonomy has no {taxon.rank} named {taxon.name}."
    try:
        accession = compare.representative_accession(vmr, taxon.name)
    except compare.NoRepresentativeAccession:
        return None, f"NCBI taxonomy has no species named {taxon.name}."
    taxid = ncbi.taxid_for_accession(accession)
    if not taxid:
        return None, f"NCBI taxonomy has no species named {taxon.name}."
    leaf = ncbi.efetch_taxonomy(taxid)
    return leaf, (f"NCBI taxonomy has no species named {taxon.name}; NCBI data are for "
                  f"{leaf.name} (taxid {leaf.taxid}), the record of its exemplar {accession}.")
