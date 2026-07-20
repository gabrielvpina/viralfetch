"""Command-line interface.

Thin by design (SPEC section 4): commands resolve config, call a service in
:mod:`viralfetch.queries`, and hand the resulting view dataclass to the
:mod:`viralfetch.render` layer. No business logic or printing lives here.
"""

from __future__ import annotations

import sys

import typer

from . import compare
from . import config as config_mod
from . import queries
from . import render
from . import sequences
from .cache import Cache
from .ncbi import NCBIClient, NCBIError
from .vmr import load

LARGE_DOWNLOAD = 500  # accessions above which a fetch asks for confirmation


def _make_client(cfg: config_mod.Config, out) -> NCBIClient:
    """Build an NCBI client or exit(3) with a helpful message if no email."""
    try:
        return NCBIClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
    except config_mod.ConfigError as exc:
        out.error(str(exc))
        raise typer.Exit(3)

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
def tax(
    ctx: typer.Context,
    name: str = typer.Argument(..., help="Taxon name (any rank)."),
    compare_ncbi: bool = typer.Option(
        False, "--compare-ncbi", help="Show the ICTV lineage beside NCBI's, highlighting divergences."
    ),
) -> None:
    """Show the full ICTV lineage of a taxon (realm -> species).

    With --compare-ncbi, fetch the NCBI taxonomy lineage for a representative
    accession and render both side by side. Divergences are expected and are
    the product of the command — NCBI commonly lags ICTV.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    vmr = load()

    if compare_ncbi:
        try:
            client = _make_client(cfg, out)
            result = compare.compare_ncbi(vmr, client, name)
        except queries.TaxonNotFound as exc:
            out.not_found(exc.name, exc.suggestions)
            raise typer.Exit(1)
        except (compare.NoRepresentativeAccession, compare.NoNcbiTaxid) as exc:
            out.error(str(exc))
            raise typer.Exit(4)
        except NCBIError as exc:
            out.error(f"NCBI request failed: {exc}")
            raise typer.Exit(4)
        out.compare(result)
        return

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


@app.command()
def seq(
    ctx: typer.Context,
    species: str = typer.Argument(None, help="Species name (VMR). Omit when using --taxon."),
    taxon: str = typer.Option(None, "--taxon", help="Operate on a whole taxon (any rank) instead of one species."),
    meta: bool = typer.Option(False, "--meta", help="Sequence metadata via esummary (default)."),
    fasta: bool = typer.Option(False, "--fasta", help="FASTA sequences via efetch."),
    gb: bool = typer.Option(False, "--gb", help="Full GenBank records via efetch."),
    moltype: str = typer.Option(None, "--moltype", help="Filter nuccore results by moltype (e.g. ssRNA); matches ss-RNA etc."),
    biomol: str = typer.Option(None, "--biomol", help="Filter nuccore results by biomol (e.g. genomic, mRNA, cRNA)."),
    protein: bool = typer.Option(False, "--protein", help="Fetch PROTEINS via elink (nuccore->protein). Not a nuccore filter."),
    output: str = typer.Option(None, "-o", "--output", help="Write records to a file instead of stdout."),
    yes: bool = typer.Option(False, "--yes", help="Skip the confirmation prompt for large downloads."),
) -> None:
    """Fetch NCBI sequence data for a species or a whole taxon.

    Accessions are resolved locally from the VMR. Formats are mutually
    exclusive; --meta is the default.

    Molecule filters: --moltype and --biomol are db=nuccore fields, filtered
    locally over the esummary result. --protein is NOT a nuccore filter:
    proteins live in db=protein, reached via elink (nuccore->protein).

    On a whole taxon, --meta shows a local aggregate (species/isolates/
    accessions, RefSeq count, composition breakdown) so you can decide whether
    to download. Fetches of more than %d accessions ask for confirmation
    unless --yes is given.
    """ % LARGE_DOWNLOAD
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)

    if (species is None) == (taxon is None):
        out.error("Provide exactly one of a species argument or --taxon.")
        raise typer.Exit(2)

    chosen = [name for name, on in (("meta", meta), ("fasta", fasta), ("gb", gb)) if on]
    if len(chosen) > 1:
        out.error("Choose only one of --meta, --fasta, --gb.")
        raise typer.Exit(2)
    mode = chosen[0] if chosen else "meta"

    if protein and (moltype or biomol):
        out.warn("--moltype/--biomol are nuccore fields and are ignored with --protein.")

    vmr = load()
    try:
        # Taxon --meta with no record-level need => cheap local aggregate.
        if taxon is not None and mode == "meta" and not protein:
            out.seq_aggregate(sequences.taxon_aggregate(vmr, taxon))
            return

        name, accessions = sequences.resolve_accessions(vmr, species=species, taxon=taxon)
        client = _make_client(cfg, out)

        if protein:
            if mode == "meta":
                out.seq_meta(name, sequences.protein_meta(client, accessions))
            else:
                rettype = "fasta" if mode == "fasta" else "gb"
                _confirm_or_exit(client.protein_uids_for(accessions), yes, cfg, out)
                out.seq_records(name, sequences.protein_records(client, accessions, rettype), output)
            return

        if mode == "meta":
            result = sequences.meta(client, accessions, moltype=moltype, biomol=biomol)
            out.seq_meta(name, result)
        else:
            rettype = "fasta" if mode == "fasta" else "gb"
            _confirm_or_exit(accessions, yes, cfg, out)
            result = sequences.records(client, accessions, rettype, moltype=moltype, biomol=biomol)
            out.seq_records(name, result, output)
    except queries.TaxonNotFound as exc:
        out.not_found(exc.name, exc.suggestions)
        raise typer.Exit(1)
    except sequences.NotASpecies as exc:
        out.error(f"{exc.name!r} is a {exc.rank}, not a species. Use --taxon to fetch a whole taxon.")
        raise typer.Exit(2)
    except NCBIError as exc:
        out.error(f"NCBI request failed: {exc}")
        raise typer.Exit(4)


def _confirm_or_exit(items: list[str], yes: bool, cfg: config_mod.Config, out) -> None:
    """Guard large downloads. Above the threshold, require an explicit yes."""
    n = len(items)
    if n <= LARGE_DOWNLOAD or yes:
        return
    # Never block on a prompt in JSON mode or a non-interactive stdin.
    if cfg.format == "json" or not sys.stdin.isatty():
        out.error(f"{n} records to download (> {LARGE_DOWNLOAD}). Re-run with --yes to proceed.")
        raise typer.Exit(2)
    if not typer.confirm(f"About to download {n} records. Continue?"):
        raise typer.Exit(0)


if __name__ == "__main__":  # pragma: no cover
    app()
