"""Rich renderer — the default human-facing output.

Rich auto-detects TTYs and disables colour when stdout is piped, satisfying the
"no decoration in pipes" requirement (SPEC section 7).
"""

from __future__ import annotations

from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.tree import Tree

import os
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


def text(chapter: Chapter, markdown: str) -> None:
    """Render a chapter's Markdown with headings, tables and italics.

    In an interactive terminal the output is paged (like ``man``/``git``):
    ``less -F`` prints short chapters inline and opens a scrollable view for
    long ones, starting at the top. In a pipe or redirect it prints directly.
    """
    md = Markdown(markdown)
    if _out.is_terminal:
        # -F: quit if it fits one screen; -R: keep colours/italics; -X: don't
        # wipe the screen on exit. Respect a user's own $LESS if they set one.
        os.environ.setdefault("LESS", "FRX")
        with _out.pager(styles=True):
            _out.print(md)
    else:
        _out.print(md)


def text_raw(markdown: str) -> None:
    """Emit the raw Markdown to stdout (for redirection to a file)."""
    sys.stdout.write(markdown if markdown.endswith("\n") else markdown + "\n")


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


def cache_cleared(removed: int, *, texts: bool, seqs: bool) -> None:
    scope = "texts" if texts and not seqs else "seqs" if seqs and not texts else "all"
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
