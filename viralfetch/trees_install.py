"""Rebuild the ICTV trees/alignments data set from the ICTV Report.

Port of the ``viraltree`` pipeline (github.com/gabrielvpina/viraltree), which
produced the bundled ``ictv-trees``. For every VMR family with a published
Report chapter it reads the chapter's *Resources* page, downloads the official
alignment and tree of each figure, parses their members (accession + name),
molecule, method and region, joins the members to the active VMR for their
lineage, and writes the same layout the bundled data set uses::

    ictv-trees/
      _index.json
      <family>/chapter.md  resources.json  provenance.json
      <family>/trees/tree<N>/alignment.fasta  tree.nwk  tree.json  members.tsv

The result lands in the user data dir, never in the installed package; it is
built in a staging directory and only swapped in once it is complete, so a
failed or interrupted run never leaves the tool without trees.

Parsing is deterministic and dependency-free (no Biopython): Newick and NEXUS
trees, FASTA / Pearson FASTA and MEGA alignments.

Design note (SPEC section 4): nothing in this module prints. It returns data;
progress is reported through an optional callback.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from . import trees as trees_mod
from . import vmr as vmr_mod
from .ictv import BASE_URL, CHAPTER_PATH, ChapterNotFound, ICTVClient, _slug

TOOL_VERSION = f"viralfetch/{__version__}"

# Acceptance gate for the per-figure ``quality`` flag in resources.json. A
# figure failing it is still kept (its files are the reference data); the flag
# only tells consumers how far to trust the parsed members.
MIN_TIPS = 3
MIN_ALIGNMENT_COVERAGE = 0.80


class TreesInstallError(Exception):
    """The rebuilt data set is unusable; the current one was left in place."""


# -- member / figure records ------------------------------------------------

# ICTV lineage ranks, as members.tsv columns and tree.json lineage keys.
_RANKS = (
    "realm", "subrealm", "kingdom", "subkingdom", "phylum", "subphylum",
    "class", "subclass", "order", "suborder", "family", "subfamily",
    "genus", "subgenus", "species",
)
_MEMBERS_COLUMNS = (
    ["tip_label", "accession", "name", "molecule", "resolved_via",
     "ncbi_taxid", "vmr_match"] + list(_RANKS) + ["ictv_id"]
)


@dataclass
class Member:
    """A tree tip / alignment sequence and what was resolved for it.

    The NCBI fields are kept (always unset) so tree.json keeps the schema of
    the bundled data set, which viraltree could enrich from NCBI.
    """

    tip_label: str
    accession: str | None = None
    name: str | None = None
    abbreviation: str | None = None
    resolved_via: str = "none"
    ncbi_db: str | None = None
    ncbi_taxid: str | None = None
    raw_seq_fetched: bool = False
    vmr_match: bool = False
    lineage: dict[str, str] = field(default_factory=dict)
    ictv_id: str | None = None


@dataclass
class TreeFigure:
    """A parsed reference phylogeny, serialised as ``tree.json``."""

    family: str
    family_slug: str
    figure_label: str
    source: dict = field(default_factory=dict)
    caption: str | None = None
    molecule: str = "unknown"      # AA | nt | unknown
    method: str | None = None
    region: str | None = None
    region_kind: str = "unknown"   # gene | protein | whole_genome | region | unknown
    region_source: str = "none"    # filename | caption | headers | none
    n_members: int = 0
    n_tips: int = 0
    alignment_coverage: float = 0.0
    tip_alignment_consistent: bool = True
    members: list[Member] = field(default_factory=list)
    parsed_at: str | None = None
    tool_version: str | None = None


# -- Newick / NEXUS -----------------------------------------------------------

# Same token grammar as Biopython's Newick reader, so labels come out the same
# as in the bundled data set (which viraltree parsed with Bio.Phylo).
_NEWICK_TOKEN_RE = re.compile(
    r"(\(|\)|[^\s\(\)\[\]\'\:\;\,]+|\:\ ?[+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?"
    r"|\,|\[(?:\\.|[^\]])*\]|\'(?:\\.|[^\'])*\'|\;|\n)"
)


class _NewickError(ValueError):
    pass


class _Clade:
    __slots__ = ("name", "children", "parent")

    def __init__(self, parent: "_Clade | None" = None):
        self.name: str | None = None
        self.children: list[_Clade] = []
        self.parent = parent


def _newick_tips(text: str) -> list[str]:
    """Tip labels of one Newick tree statement, in tree order."""
    root = current = _Clade()
    opened = closed = 0

    def close(clade: _Clade) -> _Clade | None:
        parent = clade.parent
        if parent is not None:
            parent.children.append(clade)
            clade.parent = None
        return parent

    tokens = _NEWICK_TOKEN_RE.finditer(text.strip())
    for match in tokens:
        token = match.group()
        if token.startswith("'"):
            # Two adjacent quoted chunks are an escaped quote inside one label.
            current.name = token[1:-1] if not current.name else current.name + token[:-1]
        elif token.startswith("[") or token.startswith(":") or token == "\n":
            continue  # comment / branch length
        elif token == "(":
            current = _Clade(current)
            opened += 1
        elif token == ",":
            if current is root:  # no outer parentheses: grow a new root
                root = _Clade()
                current.parent = root
            current = _Clade(close(current))
        elif token == ")":
            parent = close(current)
            if parent is None:
                raise _NewickError("parenthesis mismatch")
            current = parent
            closed += 1
        elif token == ";":
            break
        else:
            current.name = token
    if opened != closed:
        raise _NewickError(f"{opened} open vs {closed} close parentheses")
    if next(tokens, None) is not None:
        raise _NewickError("text after the closing semicolon")
    close(current)
    close(root)

    tips: list[str] = []
    stack = [root]
    while stack:
        clade = stack.pop()
        if clade.children:
            stack.extend(reversed(clade.children))
        elif clade.name:
            tips.append(clade.name)
    return tips


def _newick_statements(text: str) -> list[str]:
    """Split a Newick file into tree statements (one per ``;``-ended line run)."""
    out: list[str] = []
    buf = ""
    for line in text.splitlines():
        buf += line.rstrip()
        if buf.endswith(";"):
            out.append(buf)
            buf = ""
    if buf:
        out.append(buf)
    return out


_NEXUS_TREES_RE = re.compile(r"begin\s+trees\s*;(.*?)\bend(?:block)?\s*;", re.I | re.S)
_NEXUS_TRANSLATE_RE = re.compile(r"\btranslate\b(.*?);", re.I | re.S)
_NEXUS_TREE_RE = re.compile(r"\btree\s+[^=]*=\s*(?:\[[^\]]*\]\s*)*(.*?;)", re.I | re.S)


def _nexus_trees(text: str) -> list[list[str]]:
    """Tip lists of every tree in a NEXUS file's TREES block(s)."""
    trees: list[list[str]] = []
    for block in _NEXUS_TREES_RE.findall(text):
        translate: dict[str, str] = {}
        m = _NEXUS_TRANSLATE_RE.search(block)
        if m:
            for pair in m.group(1).split(","):
                parts = pair.strip().split(None, 1)
                if len(parts) == 2:
                    translate[parts[0]] = parts[1].strip().strip("'\"")
            block = block[:m.start()] + block[m.end():]
        for statement in _NEXUS_TREE_RE.findall(block):
            tips = _newick_tips(" ".join(statement.split()))
            trees.append([translate.get(t, t) for t in tips])
    return trees


def parse_tree_tips(text: str) -> list[str]:
    """Tip labels of an ICTV tree file, in tree order (``[]`` if unparseable).

    Handles the two irregularities seen across ICTV tree files: some are NEXUS
    rather than raw Newick, and some hold several trees (bootstrap replicates,
    or a stray 1-tip artifact after the real tree) — the tree with the most
    tips wins, since replicates share a tip set.
    """
    trees: list[list[str]] = []
    if text.lstrip()[:6].upper() == "#NEXUS":
        try:
            trees = _nexus_trees(text)
        except _NewickError:
            trees = []  # a malformed NEXUS may still hold a usable Newick line
    if not trees:
        try:
            trees = [_newick_tips(s) for s in _newick_statements(text)]
        except _NewickError:
            return []
    if not trees:
        return []
    # Labels with spaces (quoted, or NEXUS) are normalised to underscores so
    # they line up with the alignment headers.
    return [t.replace(" ", "_") for t in max(trees, key=len)]


# -- alignments -------------------------------------------------------------

def parse_alignment(text: str) -> list[tuple[str, str]]:
    """``(label, sequence)`` pairs of an ICTV alignment file.

    FASTA (comment lines before the first record are skipped) or MEGA
    (``#mega``). Returns ``[]`` for an empty, binary or unrecognised file.
    """
    stripped = text.lstrip()
    if stripped[:5].lower() == "#mega":
        return _parse_mega(text)
    if not stripped or "\x00" in stripped[:200]:
        return []
    out: list[tuple[str, str]] = []
    label: str | None = None
    seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if label is not None:
                out.append((label, "".join(seq)))
            label, seq = line[1:].strip(), []
        elif label is not None:
            seq.append("".join(line.split()))
    if label is not None:
        out.append((label, "".join(seq)))
    return out


def _parse_mega(text: str) -> list[tuple[str, str]]:
    """Parse a MEGA alignment (``#Name`` records, ``!command;`` lines)."""
    out: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("!") or s.lower().startswith("#mega"):
            continue
        if s.startswith("#"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = s[1:].strip(), []
        elif name is not None:
            seq.append(s)
    if name is not None:
        out.append((name, "".join(seq)))
    return out


_AA_ONLY = set("EFILPQZXBJOU*")  # residues that never occur in DNA/RNA
_NT_LETTERS = set("ACGTUN")


def detect_molecule(sequences: Iterable[str]) -> str:
    """``"nt"`` or ``"AA"`` from the alignment alphabet (``"unknown"`` if empty)."""
    letters = aa_hits = nt_hits = 0
    for seq in sequences:
        for ch in seq.upper():
            if ch in "-.?*XN ":  # gaps / ambiguity — uninformative
                continue
            letters += 1
            if ch in _AA_ONLY:
                aa_hits += 1
            elif ch in _NT_LETTERS:
                nt_hits += 1
    if letters == 0:
        return "unknown"
    if aa_hits / letters > 0.05:
        return "AA"
    if nt_hits / letters > 0.9:
        return "nt"
    return "AA"


# -- member labels ----------------------------------------------------------

# Accession token: 1–6 letters, optional RefSeq underscore, ≥4 digits, optional
# version. Labels put it first (``ACC_Name``), last (``Abbrev_ACC``) or in the
# middle (``ABBREV_ACC_688_bp``).
_LEADING_ACC_RE = re.compile(r"^([A-Z]{1,6}_?\d{4,}(?:\.\d+)?)(?:[_ |]|$)")
_ANY_ACC_RE = re.compile(r"(?:^|[_ |])([A-Z]{1,6}_?\d{4,}(?:\.\d+)?)(?=[_ |]|$)")
_CONTROL_RE = re.compile(r"[\x00-\x1f\x7f]")
# A curation note ICTV appends when the accession isn't in the VMR sheet.
_NOT_ON_SHEET_RE = re.compile(r"[_ ]?ACCESSION[_ ]NOT[_ ]ON[_ ]SPREADSHEET", re.I)


def _canonical_label(label: str) -> str:
    """Separator/quote/case-insensitive key matching a tip to its header."""
    text = _CONTROL_RE.sub("", label).strip().strip("'\"")
    text = _NOT_ON_SHEET_RE.sub("", text)
    return re.sub(r"[\s_]+", " ", text).strip().lower()


def parse_member_label(label: str) -> Member:
    """Split a tip/header label into accession + name."""
    cleaned = _NOT_ON_SHEET_RE.sub("", _CONTROL_RE.sub("", label).strip().strip("'\"").strip())
    cleaned = cleaned.strip("_ ")
    accession = None
    remainder = cleaned
    m = _LEADING_ACC_RE.match(cleaned)
    if m:
        accession = m.group(1)
        remainder = cleaned[m.end():]
    else:
        m = _ANY_ACC_RE.search(cleaned)
        if m:
            accession = m.group(1)
            remainder = cleaned[:m.start()] + " " + cleaned[m.end():]
    name = re.sub(r"[_\s]+", " ", remainder).strip() or None
    return Member(tip_label=label, accession=accession, name=name)


def _member_key(member: Member) -> str:
    """Identity key: base accession if any (absorbs name typos), else the label."""
    if member.accession:
        return "acc:" + member.accession.split(".", 1)[0].upper()
    return "name:" + _canonical_label(member.tip_label)


# -- method / region --------------------------------------------------------

_METHODS = [
    (re.compile(r"maximum[\s-]*likelihood|\bML\b|RAxML|IQ-?TREE|PhyML", re.I), "maximum likelihood"),
    (re.compile(r"neigh-?bou?r[\s-]*joining|\bNJ\b", re.I), "neighbor-joining"),
    (re.compile(r"bayesian|MrBayes|BEAST|posterior", re.I), "bayesian"),
    (re.compile(r"maximum[\s-]*parsimony|parsimony", re.I), "maximum parsimony"),
    (re.compile(r"minimum[\s-]*evolution", re.I), "minimum evolution"),
]

# Abbreviations use letter-lookarounds, not \b: file names delimit tokens with
# underscores, a word character. First match wins, so specific terms go first.
_REGIONS = [
    (re.compile(r"RdRp|RNA[\s-]*dependent[\s-]*RNA[\s-]*polymerase", re.I), "RdRp", "protein"),
    (re.compile(r"helicase|(?<![A-Za-z])Hel(?![A-Za-z])", re.I), "helicase", "protein"),
    (re.compile(r"replicase|(?<![A-Za-z])Rep(?![A-Za-z])", re.I), "replicase", "protein"),
    (re.compile(r"major capsid protein|(?<![A-Za-z])MCP(?![A-Za-z])", re.I),
     "major capsid protein", "protein"),
    (re.compile(r"nucleo[\s-]*capsid|nucleoprotein|(?<![A-Za-z])NP(?![A-Za-z])", re.I),
     "nucleocapsid", "protein"),
    (re.compile(r"glyco[\s-]*protein|(?<![A-Za-z])GP(?![A-Za-z])", re.I), "glycoprotein", "protein"),
    (re.compile(r"capsid|coat protein|(?<![A-Za-z])CP(?![A-Za-z])", re.I), "capsid", "protein"),
    (re.compile(r"poly[\s-]*protein", re.I), "polyprotein", "protein"),
    (re.compile(r"terminase|(?<![A-Za-z])TerL(?![A-Za-z])", re.I), "terminase", "protein"),
    (re.compile(r"movement protein|(?<![A-Za-z])MP(?![A-Za-z])", re.I), "movement protein", "protein"),
    # "L" (large) protein = the polymerase of mononegaviruses/bunyaviruses; a
    # protein/gene qualifier is required so a bare "L" cannot match.
    (re.compile(r"large protein|(?<![A-Za-z])L[\s-]*(?:protein|gene|sequences?)", re.I),
     "L (polymerase)", "protein"),
    (re.compile(r"polymerase", re.I), "polymerase", "protein"),
    # Whole-genome intent outranks a "DNA-A" sub-clause in the same caption.
    (re.compile(r"whole[\s-]*genome|full[\s-]*genome|complete genome|genomes?\b", re.I),
     "whole genome", "whole_genome"),
    (re.compile(r"DNA-?A", re.I), "DNA-A", "region"),
]


def parse_method(text: str) -> str | None:
    return next((label for pattern, label in _METHODS if pattern.search(text)), None)


def parse_region(text: str) -> tuple[str | None, str]:
    for pattern, label, kind in _REGIONS:
        if pattern.search(text):
            return label, kind
    return None, "unknown"


# -- VMR join ---------------------------------------------------------------

_ACCESSION_RE = re.compile(r"\b[A-Z]{1,6}_?\d{4,}(?:\.\d+)?\b")
_RANK_COLUMNS = {rank: rank.capitalize() for rank in _RANKS}


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip().lower())


def _base_accession(acc: str) -> str:
    return acc.split(".", 1)[0].upper()


class VMRIndex:
    """VMR rows indexed for member -> lineage resolution.

    Accession cells are often multi-segment (``"A: EU623082; B: EU623083"``),
    so every accession token in a cell points at its row. Names are matched
    against species, virus names and abbreviations.
    """

    def __init__(self, rows: Iterable[dict[str, str]]):
        self.by_accession: dict[str, dict[str, str]] = {}
        self.by_name: list[dict[str, dict[str, str]]] = [{}, {}, {}]  # species, name, abbrev
        for row in rows:
            for acc in _ACCESSION_RE.findall((row.get("Virus GENBANK accession") or "").upper()):
                self.by_accession.setdefault(_base_accession(acc), row)
            species = _norm(row.get("Species") or "")
            if species:
                self.by_name[0].setdefault(species, row)
            for i, col in ((1, "Virus name(s)"), (2, "Virus name abbreviation(s)")):
                for name in re.split(r"[;,]", row.get(col) or ""):
                    if name.strip():
                        self.by_name[i].setdefault(_norm(name), row)

    @classmethod
    def load(cls, path: Path) -> "VMRIndex":
        with path.open(newline="", encoding="utf-8") as fh:
            return cls(csv.DictReader(fh, delimiter="\t"))

    def resolve(self, member: Member) -> dict[str, str] | None:
        """The member's VMR row: by accession first, then by name."""
        if member.accession:
            row = self.by_accession.get(_base_accession(member.accession))
            if row is not None:
                return row
        if member.name:
            key = _norm(member.name)
            for index in self.by_name:
                if key in index:
                    return index[key]
        return None

    def annotate(self, members: Iterable[Member]) -> None:
        """Fill each member's lineage and ICTV id, in place."""
        for member in members:
            row = self.resolve(member)
            member.vmr_match = row is not None
            if row is not None:
                member.lineage = {r: (row.get(c) or "").strip() for r, c in _RANK_COLUMNS.items()}
                member.ictv_id = (row.get("ICTV_ID") or "").strip() or None


# -- figure assembly --------------------------------------------------------

def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _filename(url: str | None) -> str:
    return url.rsplit("/", 1)[-1] if url else ""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def build_figure(
    *,
    family: str,
    slug: str,
    figure_label: str,
    alignment_url: str | None,
    tree_url: str | None,
    tree_text: str,
    alignment_text: str,
    caption: str | None,
    chapter_url: str,
) -> TreeFigure:
    """Parse one figure's tree + alignment into a :class:`TreeFigure`."""
    tips = parse_tree_tips(tree_text)
    records = parse_alignment(alignment_text)
    header_labels = [label for label, _ in records]

    filenames = f"{_filename(alignment_url)} {_filename(tree_url)}"
    # The region (which protein/gene the tree is built on) is read from the
    # most specific place that states it: the per-figure file names, then the
    # caption (often shared by panels A/B), then the alignment headers.
    region, region_kind, region_source = None, "unknown", "none"
    for source, blob in (("filename", filenames), ("caption", caption or ""),
                         ("headers", " ".join(header_labels))):
        region, region_kind = parse_region(blob)
        if region is not None:
            region_source = source
            break

    # Members are the tree's tips; the alignment headers stand in only when
    # there is no parseable tree.
    tip_members = [parse_member_label(t) for t in tips]
    header_members = [parse_member_label(h) for h in header_labels]
    tip_keys = {_member_key(m) for m in tip_members}
    header_keys = {_member_key(m) for m in header_members}
    members: list[Member] = []
    seen: set[str] = set()
    for member in tip_members or header_members:
        key = _member_key(member)
        if key not in seen:
            seen.add(key)
            members.append(member)

    return TreeFigure(
        family=family,
        family_slug=slug,
        figure_label=figure_label,
        source={
            "chapter_url": chapter_url,
            "resources_url": chapter_url + "/resources",
            "alignment_url": alignment_url,
            "tree_url": tree_url,
            "alignment_sha256": _sha256(alignment_text),
            "tree_sha256": _sha256(tree_text),
        },
        caption=caption,
        molecule=detect_molecule(seq for _, seq in records),
        method=parse_method(f"{caption or ''} {filenames}"),
        region=region,
        region_kind=region_kind,
        region_source=region_source,
        n_members=len(members),
        n_tips=len(tips),
        alignment_coverage=len(tip_keys & header_keys) / len(tip_keys) if tip_keys else 0.0,
        tip_alignment_consistent=bool(tip_keys) and tip_keys == header_keys,
        members=members,
        parsed_at=_now(),
        tool_version=TOOL_VERSION,
    )


def figure_quality(figure: TreeFigure) -> str:
    """``"ok"`` or why the parsed members of a figure are not trustworthy."""
    if figure.molecule == "unknown":
        return "molecule_indeterminate"
    if figure.n_tips < MIN_TIPS:
        return "too_few_members"
    if figure.alignment_coverage < MIN_ALIGNMENT_COVERAGE:
        return "low_alignment_coverage"
    return "ok"


# -- writing ----------------------------------------------------------------

def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_figure(tree_dir: Path, figure: TreeFigure, tree_text: str, alignment_text: str) -> None:
    """Write a figure's files; the ICTV originals are kept verbatim."""
    tree_dir.mkdir(parents=True, exist_ok=True)
    if tree_text.strip():
        (tree_dir / "tree.nwk").write_text(tree_text, encoding="utf-8")
    if alignment_text.strip():
        (tree_dir / "alignment.fasta").write_text(alignment_text, encoding="utf-8")
    _write_json(tree_dir / "tree.json", asdict(figure))
    lines = ["\t".join(_MEMBERS_COLUMNS)]
    for m in figure.members:
        row = (
            [m.tip_label, m.accession or "", m.name or "", figure.molecule,
             m.resolved_via, m.ncbi_taxid or "", str(m.vmr_match)]
            + [m.lineage.get(rank, "") for rank in _RANKS]
            + [m.ictv_id or ""]
        )
        lines.append("\t".join(c.replace("\t", " ") for c in row))
    (tree_dir / "members.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _figure_number(label: str) -> str | None:
    m = re.search(r"fig(?:ure)?[_\s]*([0-9]+[A-Za-z]*)", label, re.I)
    return m.group(1).upper() if m else None


def _download(client: ICTVClient, url: str | None) -> str:
    """A resource file's text, or ``""`` if missing or failed to download.

    A single bad file must not sink its figure (the other file may be fine).
    """
    if not url:
        return ""
    try:
        return client.download_file(url)
    except Exception:
        return ""


def build_family(
    client: ICTVClient,
    vmr: VMRIndex,
    family: str,
    out_dir: Path,
    *,
    vmr_file: str,
    fresh: bool = True,
) -> dict:
    """Fetch, parse and write one family under ``out_dir``; return its index entry.

    A figure is kept whenever its alignment or its tree came through, even if
    the members could not be parsed — the ICTV files are the reference data.
    """
    slug = _slug(family)
    entry: dict = {"family": family, "slug": slug}
    try:
        resources = client.fetch_resources(family, fresh=fresh)
    except ChapterNotFound:
        entry.update(status="omitted", reason="no_resources", n_trees=0)
        return entry

    chapter_url = BASE_URL + CHAPTER_PATH.format(slug=slug)
    try:
        markdown = client.fetch_chapter(family, fresh=fresh).markdown
        captions = client.fetch_figure_captions(family)  # the page just fetched
    except ChapterNotFound:
        markdown, captions = "", {}

    kept: list[tuple[TreeFigure, str, str]] = []
    manifest: list[dict] = []
    for res in resources:
        alignment_text = _download(client, res.alignment_url)
        tree_text = _download(client, res.tree_url)
        caption = None
        num = _figure_number(res.figure_label)
        if num:
            # Panels ("5A") often share their figure's caption ("5").
            caption = captions.get(num) or captions.get(re.sub(r"[A-Za-z]+$", "", num))
        figure = build_figure(
            family=family, slug=slug, figure_label=res.figure_label,
            alignment_url=res.alignment_url, tree_url=res.tree_url,
            tree_text=tree_text, alignment_text=alignment_text,
            caption=caption, chapter_url=chapter_url,
        )
        vmr.annotate(figure.members)
        keep = bool(alignment_text.strip()) or figure.n_tips > 0
        manifest.append({
            "figure_label": res.figure_label, "caption": None,
            "alignment_url": res.alignment_url, "tree_url": res.tree_url,
            "molecule": figure.molecule, "n_members": figure.n_members,
            "n_tips": figure.n_tips,
            "alignment_coverage": round(figure.alignment_coverage, 3),
            "tip_alignment_consistent": figure.tip_alignment_consistent,
            "quality": figure_quality(figure), "kept": keep,
        })
        if keep:
            kept.append((figure, tree_text, alignment_text))

    if not kept:
        entry.update(status="omitted", reason="no_resources", n_trees=0)
        return entry

    family_dir = out_dir / slug
    # tree1, tree2, … in page order; the ICTV label stays in tree.json.
    for i, (figure, tree_text, alignment_text) in enumerate(kept, 1):
        _write_figure(family_dir / "trees" / f"tree{i}", figure, tree_text, alignment_text)
    if markdown:
        (family_dir / "chapter.md").write_text(markdown, encoding="utf-8")
    _write_json(family_dir / "resources.json", manifest)
    _write_json(family_dir / "provenance.json", {
        "family": family, "slug": slug, "chapter_url": chapter_url,
        "vmr_file": vmr_file, "tool_version": TOOL_VERSION, "generated_at": _now(),
    })
    entry.update(
        status="included",
        n_trees=len(kept),
        n_members=sum(f.n_members for f, _, _ in kept),
        molecules=sorted({f.molecule for f, _, _ in kept}),
    )
    return entry


# -- build + install --------------------------------------------------------

@dataclass
class InstalledTrees:
    path: Path
    included: int
    omitted: int
    errors: list[dict]       # families that failed (their previous trees are kept)
    families: list[str]      # families rebuilt in this run
    partial: bool            # only some families were rebuilt


Progress = Callable[[int, int, dict], None]


def _families_with_chapters(client: ICTVClient, fresh: bool) -> list[str]:
    """VMR families that have a published Report chapter."""
    chapters = set(client.list_chapters(fresh=fresh))
    return [f for f in vmr_mod.load().by_rank.get("family", []) if _slug(f) in chapters]


def _match_families(requested: Iterable[str]) -> list[str]:
    """Map user-typed family names to VMR family names (case-insensitive)."""
    known = {f.casefold(): f for f in vmr_mod.load().by_rank.get("family", [])}
    out, unknown = [], []
    for name in requested:
        family = known.get(name.strip().casefold())
        (out if family else unknown).append(family or name)
    if unknown:
        raise TreesInstallError(f"not a family in the active VMR: {', '.join(unknown)}")
    return out


def install(
    client: ICTVClient,
    *,
    families: Iterable[str] | None = None,
    fresh: bool = True,
    progress: Progress | None = None,
) -> InstalledTrees:
    """Rebuild the trees data set from ICTV and install it in the user data dir.

    With ``families``, only those are rebuilt and merged into the data set in
    use (bundled or installed); otherwise every family with a published chapter
    is. The new data set replaces the installed one only once it is complete
    and at least one family made it in.
    """
    requested = list(families or [])
    partial = bool(requested)
    names = _match_families(requested) if partial else _families_with_chapters(client, fresh)

    target = installed_root()
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = target.with_name(target.name + ".partial")
    shutil.rmtree(staging, ignore_errors=True)

    vmr_path = vmr_mod.active_path()
    vmr = VMRIndex.load(vmr_path)
    current = trees_mod.root()
    previous = {e.get("slug"): e for e in trees_mod.index().get("families", [])}
    _ignore = shutil.ignore_patterns(".DS_Store", "__pycache__")
    try:
        if partial:
            # Start from the data set in use, then replace the requested families.
            shutil.copytree(current, staging, ignore=_ignore)
            for name in names:
                shutil.rmtree(staging / _slug(name), ignore_errors=True)
        else:
            staging.mkdir(parents=True)
            previous = {s: e for s, e in previous.items() if s in {_slug(n) for n in names}}

        entries: list[dict] = []   # this run's outcome per family (reported)
        indexed: list[dict] = []   # what _index.json records
        for i, name in enumerate(names, 1):
            slug = _slug(name)
            try:
                entry = build_family(client, vmr, name, staging,
                                     vmr_file=vmr_path.name, fresh=fresh)
                indexed.append(entry)
            except Exception as exc:  # one bad family must not abort the run
                shutil.rmtree(staging / slug, ignore_errors=True)
                entry = {"family": name, "slug": slug, "status": "error",
                         "error": f"{type(exc).__name__}: {exc}"}
                # Keep the family's current trees rather than losing them.
                if (current / slug).is_dir() and slug in previous:
                    shutil.copytree(current / slug, staging / slug, ignore=_ignore)
                    indexed.append(previous[slug])
                else:
                    indexed.append(entry)
            entries.append(entry)
            if progress is not None:
                progress(i, len(names), entry)

        rebuilt = {e["slug"] for e in indexed}
        merged = indexed + (
            [e for s, e in previous.items() if s not in rebuilt] if partial else []
        )
        merged.sort(key=lambda e: e.get("family", "").casefold())
        counts = {
            status: sum(1 for e in merged if e.get("status") == status)
            for status in ("included", "omitted", "error")
        }
        if not counts["included"]:
            raise TreesInstallError(
                "no family could be built — the ICTV site layout may have changed"
            )
        _write_json(staging / "_index.json", {
            "tool_version": TOOL_VERSION,
            "vmr_file": vmr_path.name,
            "generated_at": _now(),
            "counts": {"total": len(merged), **counts},
            "families": merged,
        })
        _swap_in(staging, target)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    return InstalledTrees(
        path=target,
        included=counts["included"],
        omitted=counts["omitted"],
        errors=[e for e in entries if e.get("status") == "error"],
        families=names,
        partial=partial,
    )


def _swap_in(staging: Path, target: Path) -> None:
    """Replace ``target`` with ``staging`` with no window where neither exists."""
    old = target.with_name(target.name + ".old")
    shutil.rmtree(old, ignore_errors=True)
    if target.exists():
        os.replace(target, old)
    os.replace(staging, target)
    shutil.rmtree(old, ignore_errors=True)


def installed_root() -> Path:
    """Where ``viralfetch update --trees`` installs the rebuilt data set."""
    return trees_mod.installed_root()


def reset() -> bool:
    """Remove the installed data set, reverting to the bundled one.

    Returns whether anything was removed.
    """
    target = installed_root()
    if not target.exists():
        return False
    shutil.rmtree(target)
    return True
