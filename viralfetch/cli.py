"""Command-line interface.

Thin by design (SPEC section 4): commands resolve config, call a service in
:mod:`viralfetch.queries`, and hand the resulting view dataclass to the
:mod:`viralfetch.render` layer. No business logic or printing lives here.
"""

from __future__ import annotations

import sys

import typer

from . import __version__
from . import compare
from . import config as config_mod
from . import ictv
from . import queries
from . import render
from . import sequences
from .cache import Cache
from .ictv import ICTVClient
from .ncbi import NCBIClient, NCBIError
from .vmr import VMR_FILENAME, load

LARGE_DOWNLOAD = 500  # accessions above which a fetch asks for confirmation


def _make_client(cfg: config_mod.Config, out) -> NCBIClient:
    """Build an NCBI client or exit(3) with a helpful message if no email."""
    try:
        return NCBIClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
    except config_mod.ConfigError as exc:
        out.error(str(exc))
        raise typer.Exit(3)


def _make_ictv_client(cfg: config_mod.Config, out) -> ICTVClient:
    """Build an ICTV client or exit(3) if no email (needed for the User-Agent)."""
    try:
        return ICTVClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
    except config_mod.ConfigError as exc:
        out.error(str(exc))
        raise typer.Exit(3)


def _complete_taxon(incomplete: str) -> list[str]:
    """Shell-completion source: taxon names starting with the typed prefix."""
    try:
        vmr = load()
    except Exception:
        return []
    needle = incomplete.casefold()
    names = sorted({t.name for t in vmr.taxa.values()}, key=str.casefold)
    return [n for n in names if n.casefold().startswith(needle)][:40]


def _mask(key: str | None) -> str:
    """Mask an API key for display, revealing only the last four characters."""
    if not key:
        return "(not set)"
    return ("…" + key[-4:]) if len(key) > 4 else "****"

HELP = """Query and download viral taxonomy, metadata and sequences.

Combines the ICTV VMR (local, embedded), NCBI E-utilities (remote) and the
ICTV Report (remote). The VMR is the local index; everything else is fetched
on demand and cached.

Global options below apply to every command and go [bold]before[/] it, e.g.
[cyan]viralfetch --json tax Coronaviridae[/]. Each command has its own arguments
and options — run [bold]viralfetch COMMAND --help[/] to see them in their own boxes.
"""

# Help-panel titles group each command's options into their own boxes.
_FORMAT = "Output format"
_SELECT = "Molecule selection"
_TARGET = "Target & output"

# Command-group panels split the main --help listing into two sections.
_PANEL_QUERY = "Query & retrieval"
_PANEL_CONFIG = "Configuration & maintenance"

app = typer.Typer(
    name="viralfetch",
    help=HELP,
    no_args_is_help=True,
    add_completion=True,
    rich_markup_mode="rich",
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"viralfetch {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    ctx: typer.Context,
    version: bool = typer.Option(
        None, "--version", callback=_version_callback, is_eager=True,
        help="Show the version and exit.",
    ),
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


@app.command(rich_help_panel=_PANEL_QUERY)
def tax(
    ctx: typer.Context,
    name: str = typer.Argument(..., help="Taxon name (any rank).", autocompletion=_complete_taxon),
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


@app.command(rich_help_panel=_PANEL_QUERY)
def members(
    ctx: typer.Context,
    taxon: str = typer.Argument(..., help="Parent taxon name.", autocompletion=_complete_taxon),
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


@app.command(rich_help_panel=_PANEL_QUERY)
def seq(
    ctx: typer.Context,
    species: str = typer.Argument(None, help="Species name (VMR). Omit when using --taxon.", autocompletion=_complete_taxon),
    taxon: str = typer.Option(None, "--taxon", help="Operate on a whole taxon (any rank) instead of one species.", rich_help_panel=_TARGET, autocompletion=_complete_taxon),
    meta: bool = typer.Option(False, "--meta", help="Metadata via esummary (default).", rich_help_panel=_FORMAT),
    fasta: bool = typer.Option(False, "--fasta", help="FASTA sequences via efetch.", rich_help_panel=_FORMAT),
    gb: bool = typer.Option(False, "--gb", help="Full GenBank records via efetch.", rich_help_panel=_FORMAT),
    moltype: str = typer.Option(None, "--moltype", help="Filter nuccore results by moltype (e.g. ssRNA); matches ss-RNA etc.", rich_help_panel=_SELECT),
    biomol: str = typer.Option(None, "--biomol", help="Filter nuccore results by biomol (e.g. genomic, mRNA, cRNA).", rich_help_panel=_SELECT),
    protein: bool = typer.Option(False, "--protein", help="Fetch PROTEINS via elink (nuccore->protein). Not a nuccore filter.", rich_help_panel=_SELECT),
    output: str = typer.Option(None, "-o", "--output", help="Write records to a file instead of stdout.", rich_help_panel=_TARGET),
    yes: bool = typer.Option(False, "--yes", help="Skip the confirmation prompt for large downloads.", rich_help_panel=_TARGET),
) -> None:
    """Fetch NCBI sequence data for a species or a whole taxon (accessions come
    from the VMR). Output formats are mutually exclusive; --meta is the default.

    Note: --moltype/--biomol filter nuccore records locally, while --protein is
    a separate path (elink nuccore->protein), not a nuccore filter.
    """
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


@app.command(rich_help_panel=_PANEL_QUERY)
def text(
    ctx: typer.Context,
    name: str = typer.Argument(..., help="Family name (ICTV Report chapter).", autocompletion=_complete_taxon),
    section: str = typer.Option(None, "--section", help="Show only a section by heading (e.g. summary)."),
    raw: bool = typer.Option(False, "--raw", help="Emit raw Markdown to stdout (for redirecting to a file)."),
) -> None:
    """Fetch an ICTV Report chapter and render it (headings, tables, italics).

    The original page URL and the chapter's references/attribution are shown at
    the top, and the content is CC BY 4.0. Images are omitted.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    client = _make_ictv_client(cfg, out)

    # The ICTV Report is organised by family: map a genus/species/etc. to its
    # family chapter. A name unknown to the VMR is tried verbatim, and its
    # suggestions are kept in case that fetch 404s.
    vmr = load()
    suggestions: list[str] = []
    try:
        target, note = queries.report_target(vmr, name)
    except queries.TaxonNotFound as exc:
        target, note, suggestions = name, None, exc.suggestions
    if note:
        out.warn(note)  # to stderr, so --raw/--json stdout stays clean

    try:
        chapter = client.fetch_chapter(target)
        markdown = ictv.section_markdown(chapter, section) if section else chapter.markdown
    except ictv.ChapterNotFound as exc:
        if suggestions:
            out.not_found(name, suggestions)
        else:
            out.error(f"No ICTV Report chapter found for {name!r} (tried {exc.url}).")
        raise typer.Exit(1)
    except ictv.SectionNotFound as exc:
        out.error(f"Section {exc.section!r} not found. Available: {', '.join(exc.available)}.")
        raise typer.Exit(2)
    except ictv.ChapterParseError as exc:
        out.error(f"Could not parse the ICTV chapter: {exc}")
        raise typer.Exit(4)
    except ictv.ICTVError as exc:
        out.error(f"ICTV request failed: {exc}")
        raise typer.Exit(4)

    if raw and cfg.format != "json":
        out.text_raw(markdown)
    else:
        out.text(chapter, markdown)


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


@app.command(rich_help_panel=_PANEL_CONFIG)
def diagnose(ctx: typer.Context) -> None:
    """Report VMR accession-parser quality (rows that yielded zero accessions).

    The empty/unparsed counts are the parser's quality indicator (SPEC section
    6); a spike means the free-text accession field grew an unhandled shape.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    out.diagnose(queries.diagnostics(load()))


@app.command(rich_help_panel=_PANEL_CONFIG)
def update(ctx: typer.Context) -> None:
    """Check whether a newer VMR is published on ictv.global/vmr.

    The VMR ships embedded in the package; this only reports whether a newer
    release exists (and where to download it), it does not replace the file.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    client = _make_ictv_client(cfg, out)
    try:
        status = client.check_vmr_update(VMR_FILENAME)
    except ictv.ICTVError as exc:
        out.error(f"ICTV request failed: {exc}")
        raise typer.Exit(4)
    out.update_status(status)


@app.command(rich_help_panel=_PANEL_CONFIG)
def config(
    ctx: typer.Context,
    store_ncbi_email: str = typer.Option(None, "--store-ncbi-email", help="Persist an NCBI email to the config file."),
    store_ncbi_apikey: str = typer.Option(None, "--store-ncbi-apikey", help="Persist an NCBI API key to the config file."),
) -> None:
    """Show the effective NCBI email/API key (masked) and cache/config paths.

    With --store-ncbi-email / --store-ncbi-apikey, persist those values first.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    if store_ncbi_email is not None or store_ncbi_apikey is not None:
        config_mod.store(email=store_ncbi_email, api_key=store_ncbi_apikey)
        cfg = config_mod.resolve(fmt=cfg.format)  # refresh from env + file
    out.config_view({
        "email": cfg.email,
        "api_key": _mask(cfg.api_key),
        "cache_dir": str(config_mod.CACHE_DIR),
        "config_file": str(config_mod.CONFIG_FILE),
        "config_file_exists": config_mod.CONFIG_FILE.is_file(),
    })


cache_app = typer.Typer(help="Inspect or clear the on-disk cache.")
app.add_typer(cache_app, name="cache", rich_help_panel=_PANEL_CONFIG)


@cache_app.command("info")
def cache_info_cmd(ctx: typer.Context) -> None:
    """Show per-namespace entry counts and total size."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    out.cache_info(Cache(config_mod.CACHE_DIR).info())


@cache_app.command("clear")
def cache_clear_cmd(
    ctx: typer.Context,
    texts: bool = typer.Option(False, "--texts", help="Clear only ICTV chapter text (30-day TTL namespace)."),
    seqs: bool = typer.Option(False, "--seqs", help="Clear only sequences/metadata (permanent namespace)."),
) -> None:
    """Remove cached entries. With neither flag, clear everything."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    removed = Cache(config_mod.CACHE_DIR).clear(texts=texts, seqs=seqs)
    out.cache_cleared(removed, texts=texts, seqs=seqs)


if __name__ == "__main__":  # pragma: no cover
    app()
