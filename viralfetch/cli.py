"""Command-line interface.

Thin by design (SPEC section 4): commands resolve config, call a service in
:mod:`viralfetch.queries`, and hand the resulting view dataclass to the
:mod:`viralfetch.render` layer. No business logic or printing lives here.
"""

from __future__ import annotations

import typer

from . import config as config_mod
from . import queries
from . import render
from .vmr import load

app = typer.Typer(
    name="viralfetch",
    help="Query and download viral taxonomy, metadata and sequences.",
    no_args_is_help=True,
    add_completion=True,
)


@app.callback()
def main(
    ctx: typer.Context,
    json_out: bool = typer.Option(False, "--json", help="Emit pure JSON on stdout (for jq)."),
    no_cache: bool = typer.Option(False, "--no-cache", help="Force refetch, ignore cache."),
    verbose: bool = typer.Option(False, "--verbose", help="Verbose diagnostics on stderr."),
    email: str = typer.Option(None, "--email", help="NCBI email (overrides $NCBI_EMAIL)."),
    api_key: str = typer.Option(None, "--api-key", help="NCBI API key (overrides $NCBI_API_KEY)."),
) -> None:
    """Resolve global configuration and stash it on the context."""
    ctx.obj = config_mod.resolve(
        email=email,
        api_key=api_key,
        fmt="json" if json_out else "rich",
        verbose=verbose,
        no_cache=no_cache,
    )


@app.command()
def tax(ctx: typer.Context, name: str = typer.Argument(..., help="Taxon name (any rank).")) -> None:
    """Show the full ICTV lineage of a taxon (realm -> species)."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    vmr = load()
    try:
        view = queries.tax(vmr, name)
    except queries.TaxonNotFound as exc:
        out.not_found(exc.name, exc.suggestions)
        raise typer.Exit(1)
    out.tax(view)


@app.command()
def members(
    ctx: typer.Context,
    taxon: str = typer.Argument(..., help="Parent taxon name."),
    rank: str = typer.Option(None, "--rank", help="Restrict to a rank below the parent (e.g. genus)."),
    count: bool = typer.Option(False, "--count", help="Show aggregated counts only."),
    tree: bool = typer.Option(False, "--tree", help="List the full descendant subtree as a hierarchy."),
) -> None:
    """List child taxa of a taxon at any rank below it (local, no network).

    With --tree, render the entire descendant hierarchy (subfamily -> genus ->
    species) rooted at the taxon. Without flags, show a per-rank breakdown.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    vmr = load()
    try:
        if tree:
            out.members_tree(queries.members_tree(vmr, taxon))
            return
        view = queries.members(vmr, taxon, rank=rank, count=count)
    except queries.TaxonNotFound as exc:
        out.not_found(exc.name, exc.suggestions)
        raise typer.Exit(1)
    except queries.InvalidRank as exc:
        out.error(
            f"Rank {exc.rank!r} is not below {exc.taxon.rank} {exc.taxon.name!r}. "
            f"Valid ranks: {', '.join(exc.valid)}."
        )
        raise typer.Exit(2)
    out.members(view)


if __name__ == "__main__":  # pragma: no cover
    app()
