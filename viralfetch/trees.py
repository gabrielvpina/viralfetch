"""ICTV phylogenetic trees: locate, load, and resolve a taxon to its family tree.

The ``ictv-trees`` data set ships one directory per family (slug = family name
lower-cased), each with a ``chapter.md`` and one or more ``trees/tree<N>/``
folders holding a Newick tree (``tree.nwk``), its metadata (``tree.json``) and a
``members.tsv`` mapping each tip to a virus and its full ICTV lineage.

Design note (SPEC section 4): nothing here prints. Functions return view
dataclasses (or raise :class:`TreesNotFound`); only the render layer emits text.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

from . import queries
from .vmr import VMR


# -- data location --------------------------------------------------------

def _root() -> Path:
    """Filesystem path to the bundled ``ictv-trees`` data set."""
    return Path(resources.files("viralfetch").joinpath("ictv-trees"))


def _slug(family: str) -> str:
    return family.strip().lower()


# -- errors ---------------------------------------------------------------

class TreesNotFound(Exception):
    """The name is unknown to the VMR and to every tree's members."""

    def __init__(self, name: str, suggestions: list[str]):
        self.name = name
        self.suggestions = suggestions
        super().__init__(f"no tree found for {name!r}")


# -- Newick ---------------------------------------------------------------

@dataclass
class NewickNode:
    """A node of a parsed Newick tree."""

    name: str | None = None          # tip label (leaves only)
    support: str | None = None       # internal-node support value, if any
    length: float | None = None      # branch length to the parent
    children: list["NewickNode"] = field(default_factory=list)

    @property
    def is_tip(self) -> bool:
        return not self.children

    def tip_labels(self) -> list[str]:
        if self.is_tip:
            return [self.name] if self.name else []
        out: list[str] = []
        for c in self.children:
            out.extend(c.tip_labels())
        return out


def parse_newick(text: str) -> NewickNode:
    """Parse a Newick string into a :class:`NewickNode` tree.

    Handles the common flavour used by these files: unquoted tip labels,
    ``:branch_length`` suffixes, and internal-node support values after ``)``.
    """
    s = text.strip()

    def read_token(i: int) -> tuple[str, int]:
        j = i
        while j < len(s) and s[j] not in ",():;":
            j += 1
        return s[i:j].strip(), j

    def parse_node(i: int) -> tuple[NewickNode, int]:
        node = NewickNode()
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                child, i = parse_node(i)
                node.children.append(child)
                if i < len(s) and s[i] == ",":
                    i += 1
                    continue
                if i < len(s) and s[i] == ")":
                    i += 1
                    break
                break  # malformed, but stop rather than loop forever
        label, i = read_token(i)
        length = None
        if i < len(s) and s[i] == ":":
            length, i = read_token(i + 1)
        if node.children:
            node.support = label or None
        else:
            node.name = label or None
        try:
            node.length = float(length) if length else None
        except ValueError:
            node.length = None
        return node, i

    root, _ = parse_node(0)
    return root


# -- view dataclasses -----------------------------------------------------

@dataclass
class TreeDoc:
    """One family tree: its Newick, metadata, and tip → virus mapping."""

    tree_id: str                         # e.g. "tree1"
    figure_label: str | None
    caption: str | None
    molecule: str | None
    method: str | None
    region: str | None
    n_tips: int
    newick: str
    root: NewickNode
    tip_rows: dict[str, dict[str, str]]  # tip_label -> members.tsv row
    align_path: Path | None = None       # alignment.fasta beside the tree, if any
    matched: set[str] = field(default_factory=set)  # tip_labels hit by the query

    def display_name(self, tip_label: str) -> str:
        row = self.tip_rows.get(tip_label)
        if row and row.get("name"):
            return row["name"]
        return _prettify_tip(tip_label)


@dataclass
class TreesResult:
    """Resolution outcome for ``viralfetch tree <name>``."""

    query: str
    family: str
    slug: str
    trees: list[TreeDoc]
    note: str | None = None       # redirect note (genus/species -> family)
    source: str = "vmr"           # "vmr" | "member" | "ncbi"
    chapter_path: Path | None = None
    matched_rank: str | None = None   # rank the highlight fell back to (NCBI path)
    matched_value: str | None = None  # its value, e.g. the genus that was hit

    @property
    def has_trees(self) -> bool:
        return bool(self.trees)


def _prettify_tip(tip_label: str) -> str:
    label = tip_label.replace("_ACCESSION_NOT_ON_SPREADSHEET", "")
    return label.replace("_", " ").strip()


# -- loading --------------------------------------------------------------

def _load_members(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tip = row.get("tip_label")
            if tip:
                rows[tip] = row
    return rows


def _load_tree(tree_dir: Path) -> TreeDoc | None:
    nwk = tree_dir / "tree.nwk"
    if not nwk.is_file():
        return None
    newick = nwk.read_text(encoding="utf-8").strip()
    meta: dict = {}
    meta_path = tree_dir / "tree.json"
    if meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}
    members_path = tree_dir / "members.tsv"
    tip_rows = _load_members(members_path) if members_path.is_file() else {}
    root = parse_newick(newick)
    return TreeDoc(
        tree_id=tree_dir.name,
        figure_label=meta.get("figure_label"),
        caption=meta.get("caption"),
        molecule=meta.get("molecule"),
        method=meta.get("method"),
        region=meta.get("region"),
        n_tips=meta.get("n_tips") or len(root.tip_labels()),
        newick=newick,
        root=root,
        tip_rows=tip_rows,
        align_path=(tree_dir / "alignment.fasta") if (tree_dir / "alignment.fasta").is_file() else None,
    )


def _load_family_trees(slug: str) -> list[TreeDoc]:
    trees_dir = _root() / slug / "trees"
    if not trees_dir.is_dir():
        return []
    docs = []
    for tree_dir in sorted(trees_dir.iterdir()):
        if tree_dir.is_dir():
            doc = _load_tree(tree_dir)
            if doc is not None:
                docs.append(doc)
    return docs


def _chapter_path(slug: str) -> Path | None:
    path = _root() / slug / "chapter.md"
    return path if path.is_file() else None


# -- matching -------------------------------------------------------------

def _mark_matches(docs: list[TreeDoc], *, rank: str | None, value: str, name: str) -> int:
    """Flag the tips each tree that the query points at; return the hit count.

    A resolved taxon matches tips whose ``rank`` column equals its ``value``
    (so a species highlights its own tip, a genus highlights its whole clade).
    A free-text member name matches tips by their ``name`` or ``species``.
    """
    target = value.casefold()
    raw = name.casefold()
    hits = 0
    for doc in docs:
        for tip, row in doc.tip_rows.items():
            if rank and row.get(rank, "").casefold() == target:
                doc.matched.add(tip)
                hits += 1
            elif not rank and raw in (row.get("name", "").casefold(),
                                      row.get("species", "").casefold()):
                doc.matched.add(tip)
                hits += 1
    return hits


# -- resolution -----------------------------------------------------------

def resolve(
    vmr: VMR,
    name: str,
    *,
    ncbi_lineage: Callable[[str], Sequence[tuple[str, str]] | None] | None = None,
) -> TreesResult:
    """Resolve ``name`` to its family tree(s).

    First tries the VMR: any taxon (species/genus/…/family) maps to its family,
    whose tree(s) are then loaded and the query's tips highlighted. If the name
    is unknown to the VMR, falls back to searching every tree's members for a
    matching virus, and then — when ``ncbi_lineage`` is given — to NCBI taxonomy:
    the name's NCBI lineage points at a family, and its nearest rank that the
    trees do carry (species → subgenus → genus → subfamily) highlights the
    related tips. Raises :class:`TreesNotFound` if nothing matches.

    ``ncbi_lineage`` is injected (a callable returning ``[(rank, name)]`` or
    ``None``) so this module stays offline and testable; the CLI passes a
    best-effort NCBI lookup that swallows config/network errors.
    """
    taxon = None
    try:
        taxon = queries.tax(vmr, name).taxon
    except queries.TaxonNotFound:
        taxon = None

    if taxon is not None:
        family = taxon.name if taxon.rank == "family" else taxon.lineage.get("family")
        if family:
            slug = _slug(family)
            docs = _load_family_trees(slug)
            # A family query shows the whole tree — no single tip to highlight.
            if taxon.rank != "family":
                _mark_matches(docs, rank=taxon.rank, value=taxon.name, name=name)
            note = None
            if taxon.rank != "family":
                note = (f"{taxon.name!r} is a {taxon.rank}; showing the tree of "
                        f"its family, {family}.")
            return TreesResult(query=name, family=family, slug=slug, trees=docs,
                               note=note, source="vmr", chapter_path=_chapter_path(slug))

    hit = _search_members(name)
    if hit is not None:
        slug, family = hit
        docs = _load_family_trees(slug)
        _mark_matches(docs, rank=None, value=name, name=name)
        note = (f"{name!r} was found as a member of the {family} tree "
                f"(not an ICTV taxon name).")
        return TreesResult(query=name, family=family, slug=slug, trees=docs,
                           note=note, source="member", chapter_path=_chapter_path(slug))

    if ncbi_lineage is not None:
        result = _resolve_via_ncbi(name, ncbi_lineage)
        if result is not None:
            return result

    raise TreesNotFound(name, vmr.suggest(name))


# -- NCBI fallback --------------------------------------------------------

# members.tsv ranks worth highlighting, deepest first. ``family`` is left out:
# the whole tree already is the family, so marking it would highlight every tip.
_RELATED_RANKS = ("species", "subgenus", "genus", "subfamily")


def _resolve_via_ncbi(
    name: str,
    ncbi_lineage: Callable[[str], Sequence[tuple[str, str]] | None],
) -> TreesResult | None:
    """Resolve ``name`` through its NCBI lineage; ``None`` if that leads nowhere.

    NCBI knows strains, synonyms and taxa the VMR release predates. Its lineage
    shares the ICTV family names, so the family points at a bundled tree, and
    the deepest shared rank highlights the query's closest relatives on it.
    """
    lineage = ncbi_lineage(name)
    if not lineage:
        return None

    ranks = {rank: value for rank, value in lineage if rank and rank != "no rank"}
    ncbi_name = lineage[-1][1]

    family = ranks.get("family")
    via_family = family is not None
    if family:
        slug = _slug(family)
    else:
        # No family rank (an unclassified taxon): find one via a lower rank that
        # some tree's members do record.
        hit = _search_ranks(ranks)
        if hit is None:
            return None
        slug, family = hit

    docs = _load_family_trees(slug)
    rank, value = _mark_related(docs, ranks, ncbi_name)
    return TreesResult(
        query=name,
        family=family,
        slug=slug,
        trees=docs,
        note=_ncbi_note(name, ncbi_name, family, rank, value, docs, via_family),
        source="ncbi",
        chapter_path=_chapter_path(slug),
        matched_rank=rank,
        matched_value=value,
    )


def _mark_related(
    docs: list[TreeDoc], ranks: dict[str, str], ncbi_name: str
) -> tuple[str | None, str | None]:
    """Highlight the tips closest to the query; return the ``(rank, value)`` used.

    Tries the taxon's own name first (tips carry virus names as well as species),
    then walks up the lineage until a rank the trees record produces a hit —
    so a strain NCBI knows but ICTV does not still lands on its genus' clade.
    """
    if _mark_matches(docs, rank=None, value=ncbi_name, name=ncbi_name):
        return None, ncbi_name
    for rank in _RELATED_RANKS:
        value = ranks.get(rank)
        if value and _mark_matches(docs, rank=rank, value=value, name=ncbi_name):
            return rank, value
    return None, None


def _ncbi_note(
    name: str,
    ncbi_name: str,
    family: str,
    rank: str | None,
    value: str | None,
    docs: list[TreeDoc],
    via_family: bool,
) -> str:
    """Explain what NCBI resolved and how faithful the highlight is."""
    if via_family:
        head = f"{name!r} is not in the local VMR; NCBI places it in family {family}"
    else:
        head = (f"{name!r} is not in the local VMR and NCBI's lineage for it names "
                f"no family; its ranks lead to the {family} tree")
    if not docs:
        return f"{head}."
    if value is None:
        return (f"{head} — showing that tree, but no tip of it matches "
                f"{ncbi_name!r} or any rank of its lineage.")
    if rank is None:
        return f"{head} — showing that tree, highlighting {value!r}."
    return (f"{head} — no tip matches {ncbi_name!r} itself, so its {rank} "
            f"{value!r} is highlighted instead (closest relatives on the tree).")


def _search_ranks(ranks: dict[str, str]) -> tuple[str, str] | None:
    """Find a family whose members carry one of ``ranks``; deepest rank wins.

    The fallback-of-the-fallback, for NCBI lineages with no family rank.
    """
    root = _root()
    if not root.is_dir():
        return None
    for rank in _RELATED_RANKS:
        value = ranks.get(rank)
        if not value:
            continue
        needle = value.casefold()
        for family_dir in sorted(root.iterdir()):
            for members in (family_dir / "trees").glob("*/members.tsv"):
                if needle not in members.read_text(encoding="utf-8").casefold():
                    continue
                for row in _load_members(members).values():
                    if row.get(rank, "").casefold() == needle:
                        family = row.get("family") or family_dir.name.capitalize()
                        return family_dir.name, family
    return None


def _search_members(name: str) -> tuple[str, str] | None:
    """Scan every tree's members for ``name``; return ``(slug, family)`` or None.

    Only used as a fallback for names the VMR does not know, so the linear scan
    over the (small) TSVs is acceptable.
    """
    needle = name.strip().casefold()
    root = _root()
    if not root.is_dir():
        return None
    for family_dir in sorted(root.iterdir()):
        trees_dir = family_dir / "trees"
        if not trees_dir.is_dir():
            continue
        for members in trees_dir.glob("*/members.tsv"):
            text = members.read_text(encoding="utf-8")
            if needle not in text.casefold():
                continue
            rows = _load_members(members)
            for row in rows.values():
                if needle in (row.get("name", "").casefold(),
                              row.get("species", "").casefold()):
                    family = row.get("family") or family_dir.name.capitalize()
                    return family_dir.name, family
    return None


def index() -> dict:
    """The bundled ``_index.json`` (families, counts), or ``{}`` if absent."""
    path = _root() / "_index.json"
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
