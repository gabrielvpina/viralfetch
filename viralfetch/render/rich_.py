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

from ..models import RANKS, plural
from ..queries import MembersView, TaxonView

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
