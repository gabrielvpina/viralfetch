"""Rich renderer — the default human-facing output.

Rich auto-detects TTYs and disables colour when stdout is piped, satisfying the
"no decoration in pipes" requirement (SPEC section 7).
"""

from __future__ import annotations

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.tree import Tree

import sys

from ..models import RANKS, plural
from ..ncbi import MetaResult, RecordsResult
from ..queries import MembersView, TaxonTreeNode, TaxonView, TreeView

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
