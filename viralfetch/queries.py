"""Taxonomy query services over an indexed :class:`~viralfetch.vmr.VMR`.

This is the service layer the CLI calls. Every function returns plain view
dataclasses (or raises :class:`TaxonNotFound` / :class:`InvalidRank`);
nothing here prints.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .models import RANKS, Isolate, Taxon
from .vmr import VMR


class TaxonNotFound(Exception):
    """Raised when a taxon name has no exact match. Carries suggestions."""

    def __init__(self, name: str, suggestions: list[str]):
        self.name = name
        self.suggestions = suggestions
        super().__init__(f"taxon not found: {name!r}")


class InvalidRank(Exception):
    """Raised when a requested member rank is not below the target taxon."""

    def __init__(self, rank: str, taxon: Taxon, valid: list[str]):
        self.rank = rank
        self.taxon = taxon
        self.valid = valid
        super().__init__(f"invalid rank {rank!r} for {taxon.rank} {taxon.name!r}")


@dataclass
class IsolateSummary:
    total: int
    exemplars: int
    additional: int
    accessions: int
    genome_compositions: dict[str, int]
    example_accessions: list[str]


@dataclass
class TaxonView:
    taxon: Taxon
    isolate_summary: IsolateSummary | None = None


@dataclass
class MemberEntry:
    name: str
    rank: str
    species_count: int


@dataclass
class MembersView:
    parent: Taxon
    rank: str | None  # requested rank filter (None => breakdown overview)
    count_only: bool
    members: list[MemberEntry] = field(default_factory=list)
    breakdown: dict[str, int] = field(default_factory=dict)  # rank -> distinct count


@dataclass
class TaxonTreeNode:
    name: str
    rank: str
    children: list["TaxonTreeNode"] = field(default_factory=list)


@dataclass
class TreeView:
    root: TaxonTreeNode
    total: int  # number of descendant taxa (excludes the root)


def _resolve(vmr: VMR, name: str) -> Taxon:
    taxon = vmr.find(name)
    if taxon is None:
        raise TaxonNotFound(name, vmr.suggest(name))
    return taxon


def _lineage_of(vmr: VMR, iso: Isolate) -> dict[str, str]:
    """The full lineage dict for an isolate, via its species taxon."""
    sp = vmr.taxa.get(iso.species.casefold())
    return sp.lineage if sp else {}


def _rows_under(vmr: VMR, target: Taxon) -> list[Isolate]:
    """Isolate rows whose lineage carries ``target`` at ``target.rank``."""
    return [
        iso
        for iso in vmr.isolates
        if _lineage_of(vmr, iso).get(target.rank) == target.name
    ]


def descendant_isolates(vmr: VMR, name: str) -> tuple[Taxon, list[Isolate]]:
    """Resolve a taxon and return every isolate row beneath it (inclusive)."""
    taxon = _resolve(vmr, name)
    return taxon, _rows_under(vmr, taxon)


def report_target(vmr: VMR, name: str) -> tuple[str, str | None]:
    """Resolve which ICTV Report chapter (a family) describes ``name``.

    The ICTV Report is organised by family: a subfamily/genus/subgenus/species
    maps to its family's chapter, and a family maps to itself. Returns
    ``(family_name, note)`` where ``note`` explains a redirect (``None`` when
    the name was already a family, or a rank at/above family with no family in
    its lineage — in which case the name is returned unchanged to try as-is).

    Raises :class:`TaxonNotFound` (with suggestions) if the name is unknown.
    """
    taxon = _resolve(vmr, name)
    if taxon.rank == "family":
        return taxon.name, None
    family = taxon.lineage.get("family")
    if family:
        note = (
            f"{taxon.name!r} is a {taxon.rank}; the ICTV Report has no "
            f"{taxon.rank} chapter — showing its family, {family}."
        )
        return family, note
    return taxon.name, None


def tax(vmr: VMR, name: str) -> TaxonView:
    """Full lineage of a taxon, plus an isolate summary when it is a species."""
    taxon = _resolve(vmr, name)
    summary = None
    if taxon.rank == "species":
        summary = _isolate_summary(vmr, taxon.name)
    return TaxonView(taxon=taxon, isolate_summary=summary)


def _isolate_summary(vmr: VMR, species: str) -> IsolateSummary:
    isolates = vmr.isolates_by_species.get(species.casefold(), [])
    exemplars = sum(1 for i in isolates if i.exemplar)
    compositions: Counter[str] = Counter()
    n_accessions = 0
    examples: list[str] = []
    for iso in isolates:
        if iso.genome_composition:
            compositions[iso.genome_composition] += 1
        n_accessions += len(iso.accessions)
        for acc in iso.accessions:
            if len(examples) < 5:
                examples.append(acc.accession)
    return IsolateSummary(
        total=len(isolates),
        exemplars=exemplars,
        additional=len(isolates) - exemplars,
        accessions=n_accessions,
        genome_compositions=dict(compositions),
        example_accessions=examples,
    )


def members(
    vmr: VMR,
    name: str,
    rank: str | None = None,
    count: bool = False,
) -> MembersView:
    """Child taxa of ``name`` at (or aggregated below) a given rank.

    - ``rank`` set: distinct taxa at that rank descending from the target.
    - ``rank`` unset: a per-rank breakdown of distinct descendant counts.
    - ``count`` with a rank: return only the tally for that rank.
    """
    target = _resolve(vmr, name)
    target_idx = RANKS.index(target.rank)
    ranks_below = list(RANKS[target_idx + 1:])

    if rank is not None:
        rank = rank.lower()
        if rank not in ranks_below:
            raise InvalidRank(rank, target, ranks_below)

    rows = _rows_under(vmr, target)

    if rank is not None:
        entries = _members_at_rank(vmr, rows, rank)
        if count:
            return MembersView(
                parent=target,
                rank=rank,
                count_only=True,
                breakdown={rank: len(entries)},
            )
        return MembersView(parent=target, rank=rank, count_only=False, members=entries)

    # No specific rank: per-rank breakdown of distinct descendant counts.
    breakdown: dict[str, int] = {}
    for r in ranks_below:
        distinct = {
            v for iso in rows if (v := _lineage_of(vmr, iso).get(r))
        }
        if distinct:
            breakdown[r] = len(distinct)
    return MembersView(parent=target, rank=None, count_only=count, breakdown=breakdown)


def members_tree(vmr: VMR, name: str) -> TreeView:
    """Full descendant subtree of ``name``, nested by populated ranks below it.

    Optional ranks that are empty for a given lineage are skipped, so a genus
    can sit directly under a family when no subfamily is assigned — matching
    how the VMR actually records taxonomy.
    """
    target = _resolve(vmr, name)
    target_idx = RANKS.index(target.rank)
    ranks_below = RANKS[target_idx + 1:]
    rows = _rows_under(vmr, target)

    # Nested mutable tree: {(rank, name): {"rank", "name", "children": {...}}}
    root_children: dict[tuple[str, str], dict] = {}
    for iso in rows:
        lin = _lineage_of(vmr, iso)
        cursor = root_children
        for r in ranks_below:
            value = lin.get(r)
            if not value:
                continue
            key = (r, value)
            node = cursor.get(key)
            if node is None:
                node = cursor[key] = {"rank": r, "name": value, "children": {}}
            cursor = node["children"]

    count = 0

    def build(children: dict[tuple[str, str], dict]) -> list[TaxonTreeNode]:
        nonlocal count
        nodes: list[TaxonTreeNode] = []
        for data in children.values():
            count += 1
            nodes.append(
                TaxonTreeNode(
                    name=data["name"],
                    rank=data["rank"],
                    children=build(data["children"]),
                )
            )
        # Group by rank order, then alphabetically within a rank.
        nodes.sort(key=lambda n: (RANKS.index(n.rank), n.name.casefold()))
        return nodes

    root = TaxonTreeNode(
        name=target.name, rank=target.rank, children=build(root_children)
    )
    return TreeView(root=root, total=count)


def _members_at_rank(vmr: VMR, rows: list[Isolate], rank: str) -> list[MemberEntry]:
    """Distinct taxa at ``rank`` among rows, each with its species count."""
    species_by_taxon: dict[str, set[str]] = defaultdict(set)
    for iso in rows:
        lin = _lineage_of(vmr, iso)
        value = lin.get(rank)
        if not value:
            continue
        species_by_taxon[value].add(lin.get("species", ""))
    entries = [
        MemberEntry(name=name, rank=rank, species_count=len(sp - {""}))
        for name, sp in species_by_taxon.items()
    ]
    entries.sort(key=lambda e: e.name.casefold())
    return entries
