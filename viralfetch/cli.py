"""Command-line interface.

Thin by design (SPEC section 4): commands resolve config, call a service in
:mod:`viralfetch.queries`, and hand the resulting view dataclass to the
:mod:`viralfetch.render` layer. No business logic or printing lives here.
"""

from __future__ import annotations

import difflib
import shutil
import sys

import typer
from typer.core import TyperGroup

from . import __version__
from . import annotate as annotate_mod
from . import compare
from . import config as config_mod
from . import profile as profile_mod
from . import ictv
from . import queries
from . import render
from . import msa as msa_mod
from . import sequences
from . import trees as trees_mod
from . import trees_install
from . import viralzone as viralzone_mod
from . import viralzone_install
from . import vmr as vmr_mod
from . import vmr_install
from .cache import Cache
from .datasets import DatasetsClient
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


def _mask(key: str | None) -> str:
    """Mask an API key for display, revealing only the last four characters."""
    if not key:
        return "(not set)"
    return ("…" + key[-4:]) if len(key) > 4 else "****"

HELP = """Query viral taxonomy (ICTV VMR), sequences (NCBI), the ICTV Report and ViralZone.

[bold]viralfetch NAME[/] profiles any taxon or virus name: ICTV vs NCBI taxonomy,
NCBI genome counts and the ViralZone factsheet ([cyan]viralfetch Flaviviridae[/]).
Global options go [bold]before[/] the name or command: [cyan]viralfetch --csv tax Coronaviridae[/].
Run [bold]viralfetch COMMAND --help[/] for command options.
"""

# Help-panel titles group each command's options into their own boxes.
_FORMAT = "Output format"
_SELECT = "Molecule selection"
_TARGET = "Target & output"

# Command-group panels split the main --help listing into two sections.
_PANEL_QUERY = "Query & retrieval"
_PANEL_CONFIG = "Configuration & maintenance"

class _DefaultToInfo(TyperGroup):
    """Route ``viralfetch NAME ...`` to the ``info`` command when NAME is not a
    command, so the profile is the default mode."""

    def resolve_command(self, ctx, args):
        if args and args[0] not in self.commands and not args[0].startswith("-"):
            return "info", self.commands["info"], args
        return super().resolve_command(ctx, args)


app = typer.Typer(
    name="viralfetch",
    cls=_DefaultToInfo,
    help=HELP,
    no_args_is_help=True,
    add_completion=False,
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
    json_out: bool = typer.Option(False, "--json", help="JSON output."),
    csv_out: bool = typer.Option(False, "--csv", help="CSV table output."),
    tsv_out: bool = typer.Option(False, "--tsv", help="TSV table output."),
    no_cache: bool = typer.Option(False, "--no-cache", help="Ignore the cache."),
    verbose: bool = typer.Option(False, "--verbose", help="Verbose output."),
    email: str = typer.Option(None, "--email", help="NCBI email (overrides $NCBI_EMAIL)."),
    api_key: str = typer.Option(None, "--api-key", help="NCBI API key (overrides $NCBI_API_KEY)."),
) -> None:
    """Resolve global configuration and stash it on the context."""
    chosen = [f for f, on in (("json", json_out), ("csv", csv_out), ("tsv", tsv_out)) if on]
    if len(chosen) > 1:
        render.get("rich").error("Choose only one of --json, --csv and --tsv.")
        raise typer.Exit(2)
    ctx.obj = config_mod.resolve(
        email=email,
        api_key=api_key,
        fmt=chosen[0] if chosen else "rich",
        verbose=verbose,
        no_cache=no_cache,
    )
    # Remind on every call until an email is configured. `config` has its own,
    # more detailed warning (and may be storing the email right now).
    if not ctx.obj.email and ctx.invoked_subcommand != "config":
        render.get(ctx.obj.format).warn(
            "NCBI email is not configured yet — commands that reach NCBI (the profile, seq, text, "
            f"tax fallback/--ncbi/--compare-ncbi, tree/msa fallback, update) won't work. "
            f"{config_mod.EMAIL_HOWTO}"
        )


@app.command(rich_help_panel=_PANEL_QUERY)
def info(
    ctx: typer.Context,
    name: list[str] = typer.Argument(..., help="Taxon or virus name (quotes optional)."),
    local: bool = typer.Option(False, "--local", help="Skip NCBI: VMR and ViralZone only."),
    full: bool = typer.Option(
        False, "--full",
        help="Also the rank-by-rank lineage table, RefSeq list, every protein and ViralZone biology.",
    ),
) -> None:
    """Profile a taxon: ICTV vs NCBI taxonomy, NCBI genomes, ViralZone factsheet.

    Above species, also the genome sizes and proteins of its RefSeq genomes.

    The default command: [cyan]viralfetch Flaviviridae[/] is [cyan]viralfetch info Flaviviridae[/].
    ViralZone content is credited to the SIB Swiss Institute of Bioinformatics.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    query = " ".join(name).strip()

    ncbi_client = datasets_client = None
    skipped = None
    if local:
        skipped = "NCBI was skipped (--local)."
    else:
        cache = Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache)
        try:
            ncbi_client = NCBIClient(cfg, cache=cache)
            datasets_client = DatasetsClient(cfg, cache=cache)
        except config_mod.ConfigError:
            skipped = "NCBI was skipped: no NCBI email is configured."

    try:
        result = profile_mod.build(
            load(), query, ncbi=ncbi_client, datasets=datasets_client,
            viralzone=viralzone_mod.load(),
        )
    except queries.TaxonNotFound as exc:
        out.not_found(exc.name, exc.suggestions)
        command = difflib.get_close_matches(query, list(ctx.parent.command.commands), n=1)
        if command and " " not in query:
            out.warn(f"Did you mean the command `viralfetch {command[0]}`?")
        raise typer.Exit(1)
    if skipped:
        result.warnings.insert(0, skipped)
    out.profile(result, full=full)


@app.command(rich_help_panel=_PANEL_QUERY)
def annotate(
    ctx: typer.Context,
    file: str = typer.Argument(..., help="Text file with one taxon name per line ('-' reads stdin)."),
    output: str = typer.Option(
        None, "--output", "-o",
        help="Write the table here; the format follows --csv/--tsv/--json, else the extension "
             "(.csv, .tsv, .json; CSV by default).",
    ),
    local: bool = typer.Option(False, "--local", help="Skip NCBI: no taxids, no NCBI fallback."),
) -> None:
    """Annotate a list of names (species or higher taxa) with their full taxonomy.

    One row per input line: the ICTV lineage (one column per rank), genome and host,
    the NCBI taxid and the exemplar accession(s) with links. Names the VMR doesn't
    know are placed through NCBI taxonomy.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    try:
        if file == "-":
            names = annotate_mod.read_names(sys.stdin)
        else:
            with open(file, encoding="utf-8") as fh:
                names = annotate_mod.read_names(fh)
    except OSError as exc:
        out.error(f"Cannot read {file}: {exc.strerror or exc}")
        raise typer.Exit(2)
    if not names:
        out.error(f"No names in {file} (one per line; blank lines and # comments are skipped).")
        raise typer.Exit(2)

    client = None
    if not local:
        try:
            client = NCBIClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
        except config_mod.ConfigError:
            out.warn("No NCBI email configured: annotating from the local VMR only "
                     "(no taxids, no NCBI fallback).")

    annotations = annotate_mod.annotate(load(), names, ncbi=client, progress=out.annotate_progress)

    if cfg.machine_readable:
        writer = out
    else:  # the table is the product: CSV, or the format the -o extension names
        ext = (output or "").lower().rsplit(".", 1)[-1] if output and "." in output else ""
        writer = render.get({"tsv": "tsv", "tab": "tsv", "json": "json"}.get(ext, "csv"))
    writer.annotate(annotations, output)
    out.annotate_summary(annotations, output)
    if all(a.status == "not_found" for a in annotations):
        raise typer.Exit(1)


@app.command(rich_help_panel=_PANEL_QUERY)
def tax(
    ctx: typer.Context,
    name: str = typer.Argument(..., help="Taxon name (any rank)."),
    compare_ncbi: bool = typer.Option(
        False, "--compare-ncbi", help="Compare with the NCBI lineage."
    ),
    ncbi: bool = typer.Option(
        False, "--ncbi", help="Use NCBI taxonomy instead of the VMR."
    ),
) -> None:
    """Show the ICTV lineage of a taxon.

    Names missing from the VMR fall back to NCBI taxonomy.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)

    if compare_ncbi and ncbi:
        out.error("Choose only one of --ncbi and --compare-ncbi.")
        raise typer.Exit(2)

    if ncbi:
        client = _make_client(cfg, out)
        try:
            lineage = compare.tax_ncbi(client, name)
        except compare.NcbiTaxonNotFound as exc:
            out.not_found(exc.name, [])
            raise typer.Exit(1)
        except NCBIError as exc:
            out.error(f"NCBI request failed: {exc}")
            raise typer.Exit(4)
        out.tax_ncbi(lineage)
        return

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
        # Not in the local VMR: fall back to NCBI taxonomy before giving up.
        lineage = _tax_via_ncbi(cfg, out, name)
        if lineage is None:
            out.not_found(exc.name, exc.suggestions)
            raise typer.Exit(1)
        out.warn(f"{name!r} is not in the local VMR; showing its NCBI taxonomy lineage instead.")
        out.tax_ncbi(lineage)
        return
    out.tax(view)


def _tax_via_ncbi(cfg: config_mod.Config, out, name: str):
    """Best-effort NCBI lineage for a name the VMR doesn't know; never fatal.

    Returns ``None`` when NCBI knows no such taxon, when no email is configured
    (with a hint on how to set one), or when the request fails.
    """
    try:
        client = NCBIClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
    except config_mod.ConfigError:
        out.warn(
            f"{name!r} is not in the local VMR and no NCBI email is configured, so "
            "NCBI taxonomy was not searched."
        )
        return None
    try:
        return compare.lineage_via_ncbi(client, name)
    except NCBIError as exc:
        out.warn(f"NCBI taxonomy fallback failed: {exc}")
        return None


@app.command(rich_help_panel=_PANEL_QUERY)
def members(
    ctx: typer.Context,
    taxon: str = typer.Argument(..., help="Parent taxon name."),
    rank: str = typer.Option(None, "--rank", help="Only this rank (e.g. genus)."),
    count: bool = typer.Option(False, "--count", help="Counts only."),
    tree: bool = typer.Option(False, "--tree", help="Show the full descendant hierarchy."),
    host: str = typer.Option(
        None, "--host",
        help="Only species with this host: bacteria, archaea, vertebrates, invertebrates, plants, "
             "fungi, protists, algae, soil, marine, ... (VMR); human (ViralZone, per genus).",
    ),
    genome: str = typer.Option(
        None, "--genome",
        help="Only species with this genome: dsDNA, ssDNA, dsRNA, ssRNA (any polarity), "
             "ssRNA(+), ssRNA(-), ssRNA-RT, dsDNA-RT, or RNA / DNA / RT.",
    ),
) -> None:
    """List the taxa below a taxon, optionally filtered by host and genome.

    Names missing from the VMR are mapped via NCBI taxonomy.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    vmr = load()
    taxon = _redirect_via_ncbi(cfg, out, vmr, taxon)
    filters = {"host": host, "genome": genome,
               "viralzone": viralzone_mod.load() if host else None}
    try:
        if tree:
            out.members_tree(queries.members_tree(vmr, taxon, **filters))
            return
        view = queries.members(vmr, taxon, rank=rank, count=count, **filters)
    except queries.TaxonNotFound as exc:
        out.not_found(exc.name, exc.suggestions)
        raise typer.Exit(1)
    except queries.InvalidFilter as exc:
        out.error(f"Unknown --{exc.option} {exc.value!r}. Valid: {', '.join(exc.valid)}.")
        raise typer.Exit(2)
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
    species: str = typer.Argument(None, help="Species name (VMR). Omit when using --taxon."),
    taxon: str = typer.Option(None, "--taxon", help="Use a whole taxon (any rank).", rich_help_panel=_TARGET),
    meta: bool = typer.Option(False, "--meta", help="Metadata via esummary (default).", rich_help_panel=_FORMAT),
    fasta: bool = typer.Option(False, "--fasta", help="FASTA sequences via efetch.", rich_help_panel=_FORMAT),
    gb: bool = typer.Option(False, "--gb", help="Full GenBank records via efetch.", rich_help_panel=_FORMAT),
    moltype: str = typer.Option(None, "--moltype", help="Filter by moltype (e.g. ssRNA).", rich_help_panel=_SELECT),
    biomol: str = typer.Option(None, "--biomol", help="Filter by biomol (e.g. genomic, mRNA).", rich_help_panel=_SELECT),
    protein: bool = typer.Option(False, "--protein", help="Fetch the linked proteins instead.", rich_help_panel=_SELECT),
    output: str = typer.Option(None, "-o", "--output", help="Write to a file.", rich_help_panel=_TARGET),
    yes: bool = typer.Option(False, "--yes", help="Don't ask before large downloads.", rich_help_panel=_TARGET),
) -> None:
    """Fetch NCBI sequences for a species or taxon.

    Accessions come from the VMR. Choose one of --meta (default), --fasta, --gb.
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
    if taxon is not None:
        taxon = _redirect_via_ncbi(cfg, out, vmr, taxon)
    else:
        species = _redirect_via_ncbi(cfg, out, vmr, species)
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
    name: str = typer.Argument(..., help="Family name (ICTV Report chapter)."),
    section: str = typer.Option(None, "--section", help="Only this section (e.g. summary)."),
    raw: bool = typer.Option(False, "--raw", help="Print raw Markdown."),
    images: bool = typer.Option(False, "--images", help="Draw the figures in the terminal."),
    fig_width: int = typer.Option(None, "--fig-width", help="Figure width in columns.", min=8),
) -> None:
    """Show a family's ICTV Report chapter (CC BY 4.0)."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    client = _make_ictv_client(cfg, out)

    # The ICTV Report is organised by family: map a genus/species/etc. to its
    # family chapter. A name unknown to the VMR is looked up in NCBI taxonomy
    # (its family, if any), then tried verbatim; suggestions are kept for a 404.
    vmr = load()
    suggestions: list[str] = []
    try:
        target, note = queries.report_target(vmr, name)
    except queries.TaxonNotFound as exc:
        target, note, suggestions = name, None, exc.suggestions
        family = _family_via_ncbi(cfg, out, name)
        if family:
            target, note, suggestions = family, (
                f"{name!r} is not in the local VMR; NCBI places it in family "
                f"{family} — showing that chapter."
            ), []
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
        figures = None
        if _want_figures(images, raw, cfg):
            # Pillow is a core dependency, so this normally fetches; the guard
            # only skips the (large) downloads in a broken env missing Pillow,
            # where the empty list lets the renderer say so in-band.
            figures = _fetch_figures(client, chapter) if out.figures_supported() else {}
        out.text(chapter, markdown, figures=figures, fig_width=fig_width)


def _family_via_ncbi(cfg: config_mod.Config, out, name: str) -> str | None:
    """Best-effort NCBI lookup of a name's family; never fatal for `text`.

    Used as a fallback when the VMR doesn't know the name. A missing email or an
    NCBI hiccup just yields ``None`` so the caller tries the name verbatim.
    """
    try:
        client = NCBIClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
        return compare.family_via_ncbi(client, name)
    except (config_mod.ConfigError, NCBIError):
        return None


def _redirect_via_ncbi(cfg: config_mod.Config, out, vmr, name: str) -> str:
    """Map a name the VMR doesn't know to the closest VMR taxon via NCBI.

    Returns ``name`` unchanged when the VMR knows it, or when NCBI can't place
    it (no email, request failure, or no VMR taxon down to family) — the caller
    then reports "not found" with the usual suggestions.
    """
    if vmr.find(name) is not None:
        return name
    lineage = _ncbi_lineage_lookup(cfg)(name)
    taxon = compare.nearest_vmr_taxon(vmr, lineage) if lineage else None
    if taxon is None:
        return name
    out.warn(
        f"{name!r} is not in the local VMR; NCBI places it in {taxon.rank} "
        f"{taxon.name} — showing that {taxon.rank} instead."
    )
    return taxon.name


def _ncbi_lineage_lookup(cfg: config_mod.Config):
    """A best-effort NCBI lineage lookup for `tree`/`msa`; never fatal.

    Returns a callable yielding ``[(rank, name)]`` for a taxon, or ``None`` when
    NCBI knows no such name — or when there is no email configured / the request
    fails, so the local-only behaviour is unchanged offline.
    """
    def lookup(name: str):
        try:
            client = NCBIClient(cfg, cache=Cache(config_mod.CACHE_DIR, enabled=not cfg.no_cache))
            lineage = compare.lineage_via_ncbi(client, name)
        except (config_mod.ConfigError, NCBIError):
            return None
        return lineage.lineage if lineage else None

    return lookup


def _want_figures(images: bool, raw: bool, cfg: config_mod.Config) -> bool:
    """Draw figures only for the interactive Rich view — not JSON, raw, or pipes."""
    return images and not raw and not cfg.machine_readable and sys.stdout.isatty()


def _fetch_figures(client: ICTVClient, chapter) -> dict[str, bytes]:
    """Fetch each figure's bytes (cached, polite), keyed by URL for inline
    placement; unavailable ones are skipped."""
    figures: dict[str, bytes] = {}
    for img in chapter.images:
        data = client.fetch_image(img.url)
        if data is not None:
            figures[img.url] = data
    return figures


@app.command(rich_help_panel=_PANEL_QUERY)
def tree(
    ctx: typer.Context,
    name: str = typer.Argument(..., help="Virus/taxon name, or a member of a tree."),
    tree_n: int = typer.Option(None, "--tree", help="Pick a tree when a family has several (1-based).", min=1),
    newick: bool = typer.Option(False, "--newick", help="Print the Newick string."),
    chapter: bool = typer.Option(False, "--chapter", help="Show the family's chapter instead."),
) -> None:
    """Show the ICTV tree of a taxon's family, highlighting the taxon.

    Works offline; names unknown locally fall back to NCBI taxonomy.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    vmr = load()

    try:
        result = trees_mod.resolve(vmr, name, ncbi_lineage=_ncbi_lineage_lookup(cfg))
    except trees_mod.TreesNotFound as exc:
        out.not_found(name, exc.suggestions)
        raise typer.Exit(1)

    if result.note:
        out.warn(result.note)  # to stderr, so --newick/--json stdout stays clean

    if chapter:
        if result.chapter_path is None:
            out.error(f"No chapter text is bundled for {result.family}.")
            raise typer.Exit(1)
        out.tree_chapter(result.chapter_path.read_text(encoding="utf-8"))
        return

    if not result.has_trees:
        out.error(f"No published tree is bundled for family {result.family}.")
        raise typer.Exit(1)

    doc = _select_tree(result, tree_n, out)
    others = [(i + 1, d) for i, d in enumerate(result.trees) if d is not doc]

    if newick and cfg.format != "json":
        out.tree_newick(doc.newick)
    else:
        out.tree(result, doc, others)


def _select_tree(result, tree_n: int | None, out):
    """Pick a tree from a resolution result (shared by `tree` and `msa`)."""
    docs = result.trees
    default_idx = next((i for i, d in enumerate(docs) if d.matched), 0)
    idx = (tree_n - 1) if tree_n is not None else default_idx
    if not 0 <= idx < len(docs):
        out.error(f"--tree {tree_n} is out of range (family {result.family} has {len(docs)}).")
        raise typer.Exit(2)
    return docs[idx]


@app.command(rich_help_panel=_PANEL_QUERY)
def msa(
    ctx: typer.Context,
    name: str = typer.Argument(..., help="Virus/taxon name, or a member of a tree."),
    tree_n: int = typer.Option(None, "--tree", help="Pick a tree when a family has several (1-based).", min=1),
    col_range: str = typer.Option(None, "--range", help="Column window, 1-based inclusive (e.g. 100:180)."),
    consensus: bool = typer.Option(False, "--consensus", help="Prepend a per-column consensus row."),
    fasta: bool = typer.Option(False, "--fasta", help="Print as FASTA."),
) -> None:
    """Show the ICTV alignment of a taxon's family.

    Shows the columns that fit the terminal; use --range for others.
    The query's sequences are marked ▶.
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    vmr = load()

    try:
        result = trees_mod.resolve(vmr, name, ncbi_lineage=_ncbi_lineage_lookup(cfg))
    except trees_mod.TreesNotFound as exc:
        out.not_found(name, exc.suggestions)
        raise typer.Exit(1)

    if result.note:
        out.warn(result.note)
    if not result.has_trees:
        out.error(f"No published tree is bundled for family {result.family}.")
        raise typer.Exit(1)

    doc = _select_tree(result, tree_n, out)
    if doc.align_path is None:
        out.error(f"No alignment is bundled for {result.family} {doc.tree_id}.")
        raise typer.Exit(1)

    alignment = msa_mod.load_alignment(doc, result.family)
    if col_range:
        try:
            start, end = msa_mod.parse_range(col_range, alignment.total_cols)
        except ValueError:
            out.error(f"Invalid --range {col_range!r} (use e.g. 100:180).")
            raise typer.Exit(2)
    elif fasta:
        # FASTA is for files and pipelines: export every column by default.
        start, end = 1, alignment.total_cols
    else:
        # Default viewport: a leading window that fits the terminal (alv wraps
        # it into blocks); widen or move it with --range.
        term_cols = shutil.get_terminal_size().columns
        start, end = 1, min(alignment.total_cols, max(40, (term_cols - 22) * 3))

    if fasta and cfg.format != "json":
        out.msa_fasta(alignment, start, end)
    else:
        out.msa(alignment, start, end, consensus)


def _confirm_or_exit(items: list[str], yes: bool, cfg: config_mod.Config, out) -> None:
    """Guard large downloads. Above the threshold, require an explicit yes."""
    n = len(items)
    if n <= LARGE_DOWNLOAD or yes:
        return
    # Never block on a prompt in JSON / CSV / TSV mode or a non-interactive stdin.
    if cfg.machine_readable or not sys.stdin.isatty():
        out.error(f"{n} records to download (> {LARGE_DOWNLOAD}). Re-run with --yes to proceed.")
        raise typer.Exit(2)
    if not typer.confirm(f"About to download {n} records. Continue?"):
        raise typer.Exit(0)


@app.command(rich_help_panel=_PANEL_CONFIG)
def diagnose(ctx: typer.Context) -> None:
    """Report VMR rows with no parseable accession."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    out.diagnose(queries.diagnostics(load()))


@app.command(rich_help_panel=_PANEL_CONFIG)
def update(
    ctx: typer.Context,
    yes: bool = typer.Option(False, "--yes", help="Install a newer VMR without asking."),
    reset: bool = typer.Option(False, "--reset", help="Go back to the bundled VMR (or trees / ViralZone data, with --trees / --viralzone)."),
    trees: bool = typer.Option(False, "--trees", help="Rebuild the trees and alignments from the ICTV Report."),
    family: list[str] = typer.Option(None, "--family", "-f", help="With --trees: only this family (repeatable)."),
    viralzone: bool = typer.Option(False, "--viralzone", help="Rebuild the ViralZone factsheets data set."),
) -> None:
    """Update the VMR, the trees and alignments (--trees), or the ViralZone data (--viralzone).

    --trees takes 10-15 minutes for all families; update the VMR first.
    --viralzone takes about 20 minutes (one request per second).
    """
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)

    if family and not trees:
        out.error("--family only applies to --trees.")
        raise typer.Exit(2)
    if trees and viralzone:
        out.error("--trees and --viralzone cannot be combined; run them one at a time.")
        raise typer.Exit(2)
    if viralzone and yes:
        out.error("--yes does not apply to --viralzone.")
        raise typer.Exit(2)
    if reset and (yes or family):
        out.error("--reset cannot be combined with --yes or --family.")
        raise typer.Exit(2)

    if viralzone:
        if reset:
            out.viralzone_reset(viralzone_install.reset())
            return
        _update_viralzone(cfg, out)
        return

    if trees:
        if reset:
            out.trees_reset(trees_install.reset())
            return
        _update_trees(_make_ictv_client(cfg, out), family, out)
        return

    if reset:
        out.vmr_reset(vmr_install.reset(), VMR_FILENAME)
        return

    client = _make_ictv_client(cfg, out)
    try:
        status = client.check_vmr_update(vmr_mod.active_filename())
    except ictv.ICTVError as exc:
        out.error(f"ICTV request failed: {exc}")
        raise typer.Exit(4)

    # JSON / CSV / TSV: data on stdout, so never prompt; install only on --yes.
    if cfg.machine_readable or status.up_to_date:
        installed = _install_vmr(client, status, out) if yes and not status.up_to_date else None
        out.update_status(status, installed)
        return

    out.update_status(status)
    if not yes:
        if not _interactive():
            out.warn("Re-run with --yes to download and install it.")
            return
        if not typer.confirm("Download and install it?", default=False):
            return
    out.vmr_installed(_install_vmr(client, status, out))


def _update_trees(client: ICTVClient, families: list[str] | None, out) -> None:
    """Rebuild and install the trees data set, or exit(4) on failure."""
    try:
        result = trees_install.install(client, families=families, progress=out.trees_progress)
    except (ictv.ICTVError, trees_install.TreesInstallError) as exc:
        out.error(f"Could not rebuild the trees: {exc} (the current trees are unchanged).")
        raise typer.Exit(4)
    out.trees_installed(result)


def _update_viralzone(cfg: config_mod.Config, out) -> None:
    """Rebuild and install the ViralZone data set, or exit(4) on failure."""
    try:
        client = viralzone_install.ViralZoneClient(cfg)
    except config_mod.ConfigError as exc:
        out.error(str(exc))
        raise typer.Exit(3)
    try:
        result = viralzone_install.install(client, progress=out.viralzone_progress)
    except (viralzone_install.ViralZoneError, viralzone_install.ViralZoneInstallError) as exc:
        out.error(f"Could not rebuild the ViralZone data: {exc} (the current data are unchanged).")
        raise typer.Exit(4)
    out.viralzone_installed(result)


def _interactive() -> bool:
    """Whether stdin is a terminal we can prompt on."""
    return sys.stdin.isatty()


def _install_vmr(client: ICTVClient, status, out) -> vmr_install.InstalledVMR:
    """Download, convert and install ``status.latest``, or exit(4) on failure."""
    try:
        data = client.download_vmr(status.latest_url)
        return vmr_install.install(data, status.latest)
    except (ictv.ICTVError, vmr_install.VMRInstallError) as exc:
        out.error(f"Could not install {status.latest}: {exc} (the current VMR is unchanged).")
        raise typer.Exit(4)


@app.command(rich_help_panel=_PANEL_CONFIG)
def config(
    ctx: typer.Context,
    store_ncbi_email: str = typer.Option(None, "--store-ncbi-email", help="Save an NCBI email."),
    store_ncbi_apikey: str = typer.Option(None, "--store-ncbi-apikey", help="Save an NCBI API key."),
) -> None:
    """Show or store the NCBI email and API key."""
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

    # Nudge the user to persist an NCBI email: a session-scoped env var or flag
    # won't survive a new terminal, and NCBI-bound commands need one.
    if not config_mod.stored().get("email"):
        if cfg.email:
            out.warn(
                "NCBI email is set for this session only (via $NCBI_EMAIL/--email) and is "
                "not persisted. Store it with "
                "`viralfetch config --store-ncbi-email you@example.com`."
            )
        else:
            out.warn(
                "No NCBI email is stored. Commands that reach NCBI (seq, text, "
                "tax fallback/--ncbi/--compare-ncbi, tree/msa fallback, update) "
                f"need one. {config_mod.EMAIL_HOWTO}"
            )


cache_app = typer.Typer(help="Show or clear the cache.")
app.add_typer(cache_app, name="cache", rich_help_panel=_PANEL_CONFIG)


@cache_app.command("info")
def cache_info_cmd(ctx: typer.Context) -> None:
    """Show cache size."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    out.cache_info(Cache(config_mod.CACHE_DIR).info())


@cache_app.command("clear")
def cache_clear_cmd(
    ctx: typer.Context,
    texts: bool = typer.Option(False, "--texts", help="Only ICTV chapters and NCBI counts."),
    seqs: bool = typer.Option(False, "--seqs", help="Only sequences, metadata and ICTV files."),
    images: bool = typer.Option(False, "--images", help="Only chapter figures."),
) -> None:
    """Clear the cache (all of it with no flag)."""
    cfg: config_mod.Config = ctx.obj
    out = render.get(cfg.format)
    removed = Cache(config_mod.CACHE_DIR).clear(texts=texts, seqs=seqs, images=images)
    out.cache_cleared(removed, texts=texts, seqs=seqs, images=images)


if __name__ == "__main__":  # pragma: no cover
    app()
