"""Rich renderer — the default human-facing output.

Rich auto-detects TTYs and disables colour when stdout is piped, satisfying the
"no decoration in pipes" requirement (SPEC section 7).
"""

from __future__ import annotations

from rich.console import Console, Group
from rich.markdown import Markdown
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.tree import Tree

import io
import os
import re
import sys

from ..cache import NamespaceInfo
from ..compare import LineageComparison
from ..ictv import VMRUpdate
from ..models import RANKS, Chapter, plural
from ..ncbi import MetaResult, RecordsResult
from ..queries import Diagnostics, MembersView, TaxonTreeNode, TaxonView, TreeView
from ..sequences import TaxonAggregate

_out = Console()
_err = Console(stderr=True)


def tax(view: TaxonView) -> None:
    taxon = view.taxon
    tree = Tree(Text("lineage", style="bold"))
    node = tree
    for rank in RANKS:
        value = taxon.lineage.get(rank)
        if not value:
            continue
        is_self = value == taxon.name and rank == taxon.rank
        label = Text()
        label.append(f"{rank}: ", style="dim")
        label.append(value, style="bold cyan" if is_self else "white")
        node = node.add(label)

    header = Text()
    header.append(taxon.name, style="bold")
    header.append(f"  ({taxon.rank})", style="dim")
    _out.print(Panel(tree, title=header, title_align="left", expand=False))

    summary = view.isolate_summary
    if summary is not None:
        table = Table(show_header=False, box=None, pad_edge=False)
        table.add_column(style="dim")
        table.add_column()
        table.add_row("isolates", f"{summary.total} ({summary.exemplars} exemplar, {summary.additional} additional)")
        table.add_row("accessions", str(summary.accessions))
        if summary.genome_compositions:
            comp = ", ".join(f"{k} ×{v}" for k, v in summary.genome_compositions.items())
            table.add_row("genome", comp)
        if summary.example_accessions:
            table.add_row("examples", ", ".join(summary.example_accessions))
        _out.print(Panel(table, title="isolates", title_align="left", expand=False))


def members(view: MembersView) -> None:
    parent = view.parent
    if view.rank is not None and not view.count_only:
        table = Table(title=f"{plural(view.rank)} in {parent.rank} {parent.name}")
        table.add_column(view.rank, style="cyan")
        table.add_column("species", justify="right", style="green")
        for m in view.members:
            table.add_row(m.name, str(m.species_count))
        table.caption = f"{len(view.members)} {plural(view.rank)}"
        _out.print(table)
    elif view.rank is not None and view.count_only:
        n = view.breakdown.get(view.rank, 0)
        _out.print(f"[green]{n}[/] {plural(view.rank)} in {parent.rank} [bold]{parent.name}[/]")
    else:
        table = Table(title=f"members of {parent.rank} {parent.name} by rank")
        table.add_column("rank", style="cyan")
        table.add_column("count", justify="right", style="green")
        for rank, n in view.breakdown.items():
            table.add_row(rank, str(n))
        _out.print(table)
        _out.print(
            f"[dim]Tip: add [/][bold]--tree[/][dim] to list every member of "
            f"{parent.name} as a hierarchy.[/]"
        )


def _add_tree_node(parent: Tree, node: TaxonTreeNode) -> None:
    label = Text()
    label.append(node.name, style="cyan")
    label.append(f"  ({node.rank})", style="dim")
    branch = parent.add(label)
    for child in node.children:
        _add_tree_node(branch, child)


def members_tree(view: TreeView) -> None:
    root_label = Text()
    root_label.append(view.root.name, style="bold")
    root_label.append(f"  ({view.root.rank})", style="dim")
    tree = Tree(root_label)
    for child in view.root.children:
        _add_tree_node(tree, child)
    _out.print(tree)
    _out.print(f"[dim]{view.total} descendant taxa[/]")


def seq_meta(species: str, result: MetaResult) -> None:
    table = Table(title=f"nuccore metadata — {species}")
    for col in ("accession", "organism", "len", "moltype", "biomol", "topology", "completeness", "source", "updated"):
        justify = "right" if col == "len" else "left"
        table.add_column(col, justify=justify, style="cyan" if col == "accession" else None)
    for r in result.records:
        table.add_row(
            r.accession, r.organism, str(r.length or "-"), r.moltype, r.biomol,
            r.topology, r.completeness, r.sourcedb, r.updatedate,
        )
    _out.print(table)
    _out.print(f"[dim]{len(result.records)} record(s)[/]")
    if result.missing:
        warn(f"{len(result.missing)} accession(s) not returned by NCBI: {', '.join(result.missing)}")


def seq_records(species: str, result: RecordsResult, output: str | None) -> None:
    if output:
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(result.text)
        _err.print(
            f"[green]Wrote[/] {len(result.returned)} {result.rettype} record(s) "
            f"for [bold]{species}[/] to [bold]{output}[/]"
        )
    else:
        # Raw records to stdout (pipeable); no Rich decoration on sequence data.
        sys.stdout.write(result.text)
        _err.print(
            f"[dim]{len(result.returned)} {result.rettype} record(s) for {species}[/]"
        )
    if result.missing:
        warn(f"{len(result.missing)} accession(s) not returned by NCBI: {', '.join(result.missing)}")


def seq_aggregate(agg: TaxonAggregate) -> None:
    table = Table(show_header=False, box=None, pad_edge=False)
    table.add_column(style="dim")
    table.add_column(justify="right")
    table.add_row("species", str(agg.species))
    table.add_row("isolates", str(agg.isolates))
    table.add_row("accessions", str(agg.accessions))
    table.add_row("RefSeq", str(agg.refseq))
    _out.print(Panel(table, title=f"{agg.name} ({agg.rank}) — download estimate", title_align="left", expand=False))

    if agg.moltype_breakdown:
        comp = Table(title="by genome composition")
        comp.add_column("composition", style="cyan")
        comp.add_column("accessions", justify="right", style="green")
        for k, v in sorted(agg.moltype_breakdown.items(), key=lambda kv: -kv[1]):
            comp.add_row(k, str(v))
        _out.print(comp)
    _out.print(
        f"[dim]Fetch with [/][bold]--fasta[/][dim] or [/][bold]--gb[/][dim] "
        f"(use --moltype/--biomol to narrow, --yes to skip confirmation).[/]"
    )


def compare(cmp: LineageComparison) -> None:
    ncbi = cmp.ncbi
    ncbi_lineage = ncbi.lineage if ncbi else []
    ictv_names = {n.casefold() for _, n in cmp.ictv}
    ncbi_names = {n.casefold() for _, n in ncbi_lineage}

    def fmt(rank: str, name: str, other: set[str]) -> Text:
        t = Text()
        t.append(f"{rank}: ", style="dim")
        t.append(name, style="white" if name.casefold() in other else "bold yellow")
        return t

    table = Table(title=f"ICTV vs NCBI lineage — {cmp.taxon}", caption=f"rep. accession {cmp.representative_accession}")
    table.add_column("ICTV (VMR) 2026", style="cyan")
    table.add_column("NCBI (taxid " + (ncbi.taxid if ncbi else "?") + ")", style="magenta")

    ictv_lines = [fmt(r, n, ncbi_names) for r, n in cmp.ictv]
    ncbi_lines = [fmt(r, n, ictv_names) for r, n in ncbi_lineage]
    for i in range(max(len(ictv_lines), len(ncbi_lines))):
        left = ictv_lines[i] if i < len(ictv_lines) else Text("")
        right = ncbi_lines[i] if i < len(ncbi_lines) else Text("")
        table.add_row(left, right)
    _out.print(table)
    _out.print("[yellow]Highlighted[/] = present in one lineage but not the other (NCBI often lags ICTV).")


def text(
    chapter: Chapter,
    markdown: str,
    figures: dict[str, bytes] | None = None,
    fig_width: int | None = None,
) -> None:
    """Render a chapter's Markdown with headings, tables and italics.

    In an interactive terminal the output is paged (like ``man``/``git``):
    ``less -F`` prints short chapters inline and opens a scrollable view for
    long ones, starting at the top. In a pipe or redirect it prints directly.

    ``figures`` maps image URL to PNG bytes; when given, each figure is drawn as
    terminal graphics **in its original place** in the text (matched by URL to
    its ``![alt](url)`` reference). ``fig_width`` caps figure width in character
    cells (default: the full terminal width).
    """
    body = _with_figures(chapter, markdown, figures, fig_width)
    if _out.is_terminal:
        # -F: quit if it fits one screen; -R: keep colours/italics; -X: don't
        # wipe the screen on exit. Respect a user's own $LESS if they set one.
        os.environ.setdefault("LESS", "FRX")
        with _out.pager(styles=True):
            _out.print(body)
    else:
        _out.print(body)


def figures_supported() -> bool:
    """Whether terminal graphics can be drawn (needs Pillow)."""
    try:
        import PIL  # noqa: F401
        return True
    except ImportError:
        return False


# A standalone image block: a whole line that is just ``![alt](url)``. The
# parser emits chapter figures this way, so each can be swapped for graphics.
_IMG_BLOCK = re.compile(r"(?m)^[ \t]*!\[(?P<alt>[^\]]*)\]\((?P<url>[^)\s]+)\)[ \t]*$")


def _with_figures(chapter: Chapter, markdown: str, figures, fig_width):
    """Return a renderable for the chapter, drawing figures inline when asked.

    ``figures is None`` means the feature is off — render plain Markdown (image
    references stay as text). Otherwise split the Markdown at each standalone
    image block and swap in half-block graphics for the URLs we have bytes for;
    unknown or undecodable images fall back to their Markdown reference.
    """
    if figures is None:
        return Markdown(markdown)

    if not figures_supported():
        note = Text("Figures need Pillow, which appears to be missing — "
                    "reinstall viralfetch (or run: pip install pillow).", style="yellow")
        return Group(Markdown(markdown), note) if chapter.images else Markdown(markdown)

    width = min(fig_width or _out.width, _out.width)
    parts: list = []
    last = 0
    for m in _IMG_BLOCK.finditer(markdown):
        art = _blockart(figures[m["url"]], width) if m["url"] in figures else None
        if art is None:
            continue  # leave this reference in the surrounding Markdown text
        pre = markdown[last:m.start()]
        if pre.strip():
            parts.append(Markdown(pre))
        if m["alt"].strip():
            parts.append(Text(m["alt"].strip(), style="italic dim"))
        parts.append(art)
        last = m.end()
    tail = markdown[last:]
    if tail.strip() or not parts:
        parts.append(Markdown(tail))
    return Group(*parts)


# Block-Elements glyphs (U+2580–U+259F) indexed by a 4-bit coverage mask over
# the cell's 2x2 sub-pixel grid — bit 3=top-left, 2=top-right, 1=bottom-left,
# 0=bottom-right. A "1" bit is painted in the foreground colour, a "0" bit in
# the background. This whole block enjoys near-universal font support, so the
# picture renders the same on any truecolor terminal — no graphics protocol.
_GLYPHS = " ▗▖▄▝▐▞▟▘▚▌▙▀▜▛█"

# For each mask, the sub-pixel indices (0=TL,1=TR,2=BL,3=BR) that are "on"/"off".
_MASKS = [([i for i in range(4) if p >> (3 - i) & 1],
          [i for i in range(4) if not (p >> (3 - i) & 1)]) for p in range(16)]

# Non-space glyphs, for tests/callers that need to spot an image row.
_BLOCK_GLYPHS = _GLYPHS[1:]

# sRGB <-> linear-light lookup tables. Downscaling averages colours; doing that
# in linear light (not gamma-encoded sRGB) keeps brightness and hue faithful
# instead of muddy/darkened.
def _srgb_to_linear(c):
    c /= 255
    return (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4) * 255


def _linear_to_srgb(c):
    c /= 255
    return (c * 12.92 if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055) * 255


_SRGB2LIN = [round(_srgb_to_linear(i)) for i in range(256)]
_LIN2SRGB = [round(_linear_to_srgb(i)) for i in range(256)]

# Floyd–Steinberg error-diffusion weights (dx, dy, weight/16) — spreads each
# sub-pixel's quantisation error to its not-yet-drawn neighbours, so smooth
# gradients dither instead of banding under the two-colours-per-cell limit.
_FS = ((1, 0, 7 / 16), (-1, 1, 3 / 16), (0, 1, 5 / 16), (1, 1, 1 / 16))


def _blockart(data: bytes, cols: int):
    """Convert PNG/JPEG bytes to a Rich ``Text`` of truecolor block-element art.

    Each character cell packs a 2x2 grid of sub-pixels — double the horizontal
    density of a half-block — by choosing, per cell, the block glyph and the two
    colours (foreground/background) that best fit its four sub-pixels. The image
    is downscaled in linear light (LANCZOS) and error-diffusion dithered, so
    gradients stay smooth. Returns ``None`` if the bytes cannot be decoded.
    Width is capped at ``cols`` cells; height follows the aspect ratio.
    """
    from PIL import Image

    try:
        img = Image.open(io.BytesIO(data)).convert("RGB")
    except Exception:
        return None

    src_w, src_h = img.size
    cols = max(1, min(cols, src_w // 2 or 1))  # 2 sub-pixels per cell, no upscaling
    sub_w = 2 * cols
    # Cells are ~twice as tall as wide; sampling 2 sub-px per cell in each axis,
    # this height keeps the picture's aspect ratio on screen.
    sub_h = max(2, round(cols * src_h / src_w))
    if sub_h % 2:  # need an even count to pair rows into cells
        sub_h += 1
    # Gamma-correct downscale: sRGB -> linear, resize (average in linear), -> sRGB.
    # RGB point() wants one 256-entry table per band, hence the ``* 3``.
    img = img.point(_SRGB2LIN * 3).resize((sub_w, sub_h), Image.LANCZOS).point(_LIN2SRGB * 3)

    # Mutable float buffer so error diffusion can bleed into later sub-pixels.
    px = img.load()
    buf = [[list(px[x, y]) for x in range(sub_w)] for y in range(sub_h)]

    text = Text(no_wrap=True, overflow="crop")
    for cy in range(0, sub_h, 2):
        for cx in range(0, sub_w, 2):
            cell = (buf[cy][cx], buf[cy][cx + 1], buf[cy + 1][cx], buf[cy + 1][cx + 1])
            mask, fg, bg = _best_cell(cell)
            _diffuse(buf, cx, cy, cell, mask, fg, bg, sub_w, sub_h)
            text.append(_GLYPHS[mask], style=f"#{_hex(fg)} on #{_hex(bg)}")
        if cy + 2 < sub_h:
            text.append("\n")
    return text


def _best_cell(cell):
    """Pick the block mask + (fg, bg) colours that best fit four sub-pixels.

    ``cell`` is ``(TL, TR, BL, BR)`` RGB tuples. Every 2-colour split of the
    cell corresponds to exactly one block glyph; the split with the lowest
    squared-error reconstruction wins (this includes the plain half-block).
    """
    best_err = None
    best = None
    for mask, (on, off) in enumerate(_MASKS):
        fg = _mean(cell, on) if on else _mean(cell, off)
        bg = _mean(cell, off) if off else fg
        err = sum(_dist2(cell[i], fg) for i in on) + sum(_dist2(cell[i], bg) for i in off)
        if best_err is None or err < best_err:
            best_err = err
            best = (mask, fg, bg)
    return best


def _diffuse(buf, cx, cy, cell, mask, fg, bg, sub_w, sub_h):
    """Push each sub-pixel's (colour − assigned) error onto later neighbours."""
    positions = ((cx, cy), (cx + 1, cy), (cx, cy + 1), (cx + 1, cy + 1))
    for k, (x, y) in enumerate(positions):
        assigned = fg if (mask >> (3 - k)) & 1 else bg
        err = (cell[k][0] - assigned[0], cell[k][1] - assigned[1], cell[k][2] - assigned[2])
        for dx, dy, w in _FS:
            nx, ny = x + dx, y + dy
            if not (0 <= nx < sub_w and 0 <= ny < sub_h):
                continue
            if cx <= nx <= cx + 1 and cy <= ny <= cy + 1:
                continue  # inside the current cell — already decided
            target = buf[ny][nx]
            target[0] += err[0] * w
            target[1] += err[1] * w
            target[2] += err[2] * w


def _mean(cell, idxs):
    n = len(idxs)
    return (sum(cell[i][0] for i in idxs) / n,
            sum(cell[i][1] for i in idxs) / n,
            sum(cell[i][2] for i in idxs) / n)


def _dist2(c, m):
    return (c[0] - m[0]) ** 2 + (c[1] - m[1]) ** 2 + (c[2] - m[2]) ** 2


def _hex(c):
    return f"{_clamp(c[0]):02x}{_clamp(c[1]):02x}{_clamp(c[2]):02x}"


def _clamp(v):
    return 0 if v < 0 else 255 if v > 255 else round(v)


def text_raw(markdown: str) -> None:
    """Emit the raw Markdown to stdout (for redirection to a file)."""
    sys.stdout.write(markdown if markdown.endswith("\n") else markdown + "\n")


# -- phylogenetic trees ---------------------------------------------------

def tree(result, doc, others) -> None:
    """Render a family tree as an indented cladogram, highlighting the query.

    ``doc`` is the selected tree; ``others`` is ``[(n, TreeDoc)]`` for the
    family's remaining trees (shown as a footer hint). Long trees are paged.
    """
    header = Text(" · ".join(_tree_header_bits(result, doc)), style="bold")
    renderables = [header, _ascii_tree(doc, _out.width)]
    if doc.caption:
        renderables.append(Text(f"\nCaption: {doc.caption}", style="dim"))
    if others:
        labels = ", ".join(f"--tree {n} ({d.region or d.tree_id})" for n, d in others)
        renderables.append(Text(f"\nOther trees for {result.family}: {labels}", style="dim"))
    body = Group(*renderables)
    if _out.is_terminal:
        os.environ.setdefault("LESS", "FRX")
        with _out.pager(styles=True):
            _out.print(body)
    else:
        _out.print(body)


def _tree_header_bits(result, doc) -> list[str]:
    bits = [result.family]
    if doc.figure_label:
        bits.append(doc.figure_label.replace("_", " "))
    if doc.region:
        bits.append(f"{doc.region} ({doc.molecule})" if doc.molecule else doc.region)
    elif doc.molecule:
        bits.append(doc.molecule)
    if doc.method:
        bits.append(doc.method)
    bits.append(f"{doc.n_tips} tips")
    if len(doc.matched) == 1:
        bits.append("1 match")
    elif doc.matched:
        bits.append(f"{len(doc.matched)} matches")
    return bits


# Combine an incoming branch from the left with whatever box glyph is already
# at a junction cell, so branches meet cleanly.
_LEFT_ARM = {" ": "─", "│": "┤", "┌": "┬", "└": "┴", "├": "┼"}


def _terminals(node) -> list:
    if node.is_tip:
        return [node]
    out = []
    for child in node.children:
        out.extend(_terminals(child))
    return out


def _ascii_tree(doc, width: int) -> Text:
    """Draw the tree left-to-right: tips one per line, branches in box-drawing.

    This handles deep, ladder-like phylogenies (common here) that an indented
    outline would push far off-screen. Column x is proportional to a node's
    distance from the root (branch lengths, or levels if none), scaled to fit;
    the queried tips are highlighted.
    """
    root = doc.root
    tips = _terminals(root)
    if not tips:
        return Text("(empty tree)")

    labels = {id(t): (doc.display_name(t.name) if t.name else "?") for t in tips}
    matched = {id(t) for t in tips if t.name in doc.matched}
    max_label = max(len(v) for v in labels.values())
    draw_w = max(12, min(width, 100) - max_label - 3)
    height = 2 * len(tips) - 1

    row: dict[int, int] = {id(t): 2 * i for i, t in enumerate(tips)}

    def calc_row(node) -> int:
        if node.is_tip:
            return row[id(node)]
        child_rows = [calc_row(c) for c in node.children]
        r = (child_rows[0] + child_rows[-1]) // 2
        row[id(node)] = r
        return r

    calc_row(root)

    depth: dict[int, float] = {}

    def calc_depth(node, acc: float) -> None:
        depth[id(node)] = acc
        for c in node.children:
            calc_depth(c, acc + (c.length or 0.0))

    calc_depth(root, 0.0)
    if max(depth.values()) <= 0:  # no usable branch lengths — fall back to levels
        depth.clear()

        def calc_level(node, lvl: int) -> None:
            depth[id(node)] = lvl
            for c in node.children:
                calc_level(c, lvl + 1)

        calc_level(root, 0)
    maxd = max(depth.values()) or 1
    col = {k: int(round(v / maxd * (draw_w - 1))) for k, v in depth.items()}

    grid = [[" "] * draw_w for _ in range(height)]

    def draw(node, startcol: int) -> None:
        r, c = row[id(node)], col[id(node)]
        for x in range(startcol, c):
            if grid[r][x] == " ":
                grid[r][x] = "─"
        if node.children:
            top, bot = row[id(node.children[0])], row[id(node.children[-1])]
            for y in range(top, bot + 1):
                if grid[y][c] == " ":
                    grid[y][c] = "│"
            for child in node.children:
                grid[row[id(child)]][c] = "├"
            grid[top][c] = "┌"
            grid[bot][c] = "└"
            grid[r][c] = _LEFT_ARM.get(grid[r][c], grid[r][c])
            for child in node.children:
                draw(child, c + 1)

    draw(root, 0)

    tip_by_row = {row[id(t)]: t for t in tips}
    out = Text()
    for y in range(height):
        out.append("".join(grid[y]).rstrip(), style="dim")
        tip = tip_by_row.get(y)
        if tip is not None:
            if id(tip) in matched:
                out.append(f" {labels[id(tip)]}  ← match", style="bold yellow")
            else:
                out.append(f" {labels[id(tip)]}")
        if y != height - 1:
            out.append("\n")
    return out


def tree_newick(newick: str) -> None:
    """Emit the raw Newick string to stdout (for piping to other tools)."""
    sys.stdout.write(newick if newick.endswith("\n") else newick + "\n")


def tree_chapter(markdown: str) -> None:
    """Render a family's bundled chapter text (Markdown), paged like ``text``."""
    md = Markdown(markdown)
    if _out.is_terminal:
        os.environ.setdefault("LESS", "FRX")
        with _out.pager(styles=True):
            _out.print(md)
    else:
        _out.print(md)


# -- multiple sequence alignments -----------------------------------------

def msa(alignment, start: int, end: int, show_consensus: bool) -> None:
    """Render columns ``start..end`` (1-based) of an alignment, coloured by
    residue and wrapped into blocks with a column ruler (via ``alv``).

    Matched records are marked ``▶``; an optional consensus row leads. ``alv``
    windows and block-wraps itself, so no horizontal scrolling is needed.
    """
    import contextlib
    import io as _io

    from alv import get_alv_objects
    from alv.alignmentterminal import AlignmentTerminal
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from ..msa import consensus as _consensus

    # Cap names so alv's left margin stays narrow enough to keep a name and its
    # sequence on one line (alv sizes the margin to the longest id).
    def _rid(name: str, matched: bool) -> str:
        mark = "▶ " if matched else ""
        return (mark + name)[:22]

    records = []
    if show_consensus:
        records.append(SeqRecord(Seq(_consensus(alignment.rows)), id="consensus", description=""))
    for row in alignment.rows:
        records.append(SeqRecord(Seq(row.seq), id=_rid(row.name, row.matched), description=""))

    header = Text(" · ".join(_msa_header_bits(alignment, start, end)), style="bold")
    if not records:
        _out.print(Group(header, Text("(empty alignment)", style="dim")))
        return

    width = _out.width if _out.is_terminal else 100
    # al_end is exclusive, so `end` (1-based inclusive) maps straight across.
    alignment_obj, painter = get_alv_objects(
        MultipleSeqAlignment(records), "guess", "guess", 1, start - 1, end
    )
    terminal = AlignmentTerminal(width)
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        # Pass chosen_width=0 so alv sizes each block to the terminal *after*
        # the name margin — otherwise it assumes a tiny margin and the lines
        # overflow and wrap, splitting names from their sequences.
        terminal.output_alignment(alignment_obj, painter, 0)
    body = Group(header, Text.from_ansi(buf.getvalue()))
    if _out.is_terminal:
        os.environ.setdefault("LESS", "RX")
        with _out.pager(styles=True):
            _out.print(body)
    else:
        _out.print(body)


def _msa_header_bits(a, start: int, end: int) -> list[str]:
    bits = [a.family, a.tree_id]
    if a.molecule:
        bits.append(a.molecule)
    bits.append(f"cols {start}–{end} of {a.total_cols}")
    bits.append(f"{a.n_rows} seqs")
    if len(a.matched_names) == 1:
        bits.append("1 match")
    elif a.matched_names:
        bits.append(f"{len(a.matched_names)} matches")
    return bits


def msa_fasta(alignment, start: int, end: int) -> None:
    """Emit the windowed alignment as FASTA to stdout, for piping onward."""
    for row in alignment.rows:
        sys.stdout.write(f">{row.name}\n{row.seq[start - 1:end]}\n")


def not_found(name: str, suggestions: list[str]) -> None:
    msg = Text()
    msg.append("No taxon named ", style="red")
    msg.append(repr(name), style="bold red")
    msg.append(".")
    if suggestions:
        msg.append("\nDid you mean: ", style="yellow")
        msg.append(", ".join(suggestions), style="bold")
    _err.print(msg)


def error(message: str) -> None:
    _err.print(Text(message, style="bold red"))


def warn(message: str) -> None:
    _err.print(Text(message, style="yellow"))


# -- utilities (cache / config / update / diagnose) -----------------------

def _human_bytes(n: int) -> str:
    size = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"


def cache_info(infos: list[NamespaceInfo]) -> None:
    table = Table(title="cache", title_style="bold")
    table.add_column("namespace")
    table.add_column("entries", justify="right")
    table.add_column("size", justify="right")
    for info in infos:
        table.add_row(info.namespace, str(info.entries), _human_bytes(info.bytes))
    _out.print(table)


def cache_cleared(removed: int, *, texts: bool, seqs: bool, images: bool) -> None:
    picked = [n for n, f in (("texts", texts), ("seqs", seqs), ("images", images)) if f]
    scope = picked[0] if len(picked) == 1 else "all"
    _out.print(Text(f"Cleared {removed} cached entr{'y' if removed == 1 else 'ies'} ({scope}).",
                    style="green"))


def config_view(view: dict) -> None:
    lines = Text()
    lines.append("NCBI email:   ", style="dim")
    lines.append(f"{view['email'] or '(not set)'}\n")
    lines.append("NCBI api key: ", style="dim")
    lines.append(f"{view['api_key']}\n")
    lines.append("cache dir:    ", style="dim")
    lines.append(f"{view['cache_dir']}\n")
    lines.append("config file:  ", style="dim")
    exists = "" if view["config_file_exists"] else "  (not created yet)"
    lines.append(f"{view['config_file']}{exists}")
    _out.print(Panel(lines, title="viralfetch config", expand=False))


def update_status(u: VMRUpdate) -> None:
    if u.up_to_date:
        _out.print(Text(f"VMR is up to date ({u.current}).", style="green"))
        return
    msg = Text()
    msg.append("A newer VMR is available.\n", style="bold yellow")
    msg.append("  current: ", style="dim")
    msg.append(f"{u.current}\n")
    msg.append("  latest:  ", style="dim")
    msg.append(f"{u.latest}\n")
    msg.append("  download: ", style="dim")
    msg.append(f"{u.latest_url}")
    _out.print(msg)


def diagnose(d: Diagnostics) -> None:
    summary = Table.grid(padding=(0, 2))
    summary.add_row("isolates", str(d.isolates))
    summary.add_row("accessions", str(d.accessions))
    summary.add_row("empty-accession rows", str(d.empty_accession_rows))
    summary.add_row("unparsed rows", str(len(d.unparsed)))
    _out.print(Panel(summary, title="VMR accession-parser diagnostics", expand=False))
    if d.unparsed:
        table = Table(title=f"rows that yielded zero accessions ({len(d.unparsed)})")
        table.add_column("species")
        table.add_column("raw accession field")
        for species, raw in d.unparsed:
            table.add_row(species, raw)
        _out.print(table)
