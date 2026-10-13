"""CSV / TSV renderer — one delimited table on stdout, ready for R, pandas or awk.

Every command maps onto a tidy table: one row per observation, a header row,
no decoration. Multi-line values are joined with " | " so each record stays on
one line; booleans are ``true``/``false``; missing values are empty. Raw
formats the user asked for explicitly (``--newick``, ``--fasta``, ``--raw``,
FASTA/GenBank records) are written as they are. Warnings and errors go to
stderr as plain text.
"""

from __future__ import annotations

import csv
import re
import sys
from collections.abc import Iterable
from dataclasses import asdict

from ..models import RANKS
from ..viralzone import CREDIT


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple, set)):
        return ";".join(_cell(v) for v in value)
    text = str(value)
    return re.sub(r"\s*[\r\n]+\s*", " | ", text).replace("\t", " ").strip()


def _flatten(d: dict, prefix: str = "") -> dict:
    out: dict = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, key + "."))
        else:
            out[key] = v
    return out


def _sections(markdown: str) -> list[tuple[str, int, str]]:
    """``(heading, level, text)`` per Markdown section (text before the first
    heading is section "", level 0)."""
    out: list[tuple[str, int, str]] = []
    heading, level, buf = "", 0, []
    for line in markdown.splitlines():
        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            if heading or "".join(buf).strip():
                out.append((heading, level, "\n".join(buf).strip()))
            heading, level, buf = m.group(2).strip(), len(m.group(1)), []
        else:
            buf.append(line)
    if heading or "".join(buf).strip():
        out.append((heading, level, "\n".join(buf).strip()))
    return out


class TableRenderer:
    def __init__(self, delimiter: str):
        self.delimiter = delimiter

    # -- plumbing -------------------------------------------------------------

    def _write(self, header: list[str], rows: Iterable, output: str | None = None) -> None:
        """Write a table; rows are sequences or dicts keyed by ``header``."""
        fh = open(output, "w", newline="", encoding="utf-8") if output else sys.stdout
        try:
            writer = csv.writer(fh, delimiter=self.delimiter, lineterminator="\n")
            writer.writerow(header)
            for row in rows:
                if isinstance(row, dict):
                    row = [row.get(col) for col in header]
                writer.writerow([_cell(v) for v in row])
        finally:
            if output:
                fh.close()

    def _record(self, d: dict) -> None:
        """A single flat record as a one-row table."""
        flat = _flatten(d)
        self._write(list(flat), [flat])

    @staticmethod
    def _stderr(prefix: str, message: str) -> None:
        print(f"{prefix}{message}", file=sys.stderr)

    # -- messages -------------------------------------------------------------

    def not_found(self, name: str, suggestions: list[str]) -> None:
        msg = f"No taxon named {name!r}."
        if suggestions:
            msg += " Did you mean: " + ", ".join(suggestions)
        self._stderr("error: ", msg)

    def error(self, message: str) -> None:
        self._stderr("error: ", message)

    def warn(self, message: str) -> None:
        self._stderr("warning: ", message)

    def figures_supported(self) -> bool:
        return False

    # -- taxonomy -------------------------------------------------------------

    def tax(self, view) -> None:
        taxon = view.taxon
        rows = [(r, taxon.lineage[r]) for r in RANKS if taxon.lineage.get(r)]
        self._write(["rank", "name"], rows)

    def tax_ncbi(self, lineage) -> None:
        rows = [(r or "no rank", n, lineage.taxid if n == lineage.name else None)
                for r, n in lineage.lineage]
        self._write(["rank", "name", "taxid"], rows)

    def compare(self, cmp) -> None:
        rows = [("ictv", r, n) for r, n in cmp.ictv]
        if cmp.ncbi is not None:
            rows += [("ncbi", r or "no rank", n) for r, n in cmp.ncbi.lineage]
        self._write(["source", "rank", "name"], rows)

    def members(self, view) -> None:
        filters = [view.filters.get("host"), view.filters.get("genome")]
        if view.rank is not None and not view.count_only:
            self._write(["name", "rank", "species_count", "parent", "host_filter", "genome_filter"],
                        [(m.name, m.rank, m.species_count, view.parent.name, *filters)
                         for m in view.members])
        else:
            self._write(["rank", "count", "parent", "host_filter", "genome_filter"],
                        [(r, n, view.parent.name, *filters) for r, n in view.breakdown.items()])

    def members_tree(self, view) -> None:
        rows: list = []

        def walk(node, parent: str | None, path: list[str], depth: int) -> None:
            path = path + [node.name]
            if depth:
                rows.append((depth, node.rank, node.name, parent, " > ".join(path)))
            for child in node.children:
                walk(child, node.name, path, depth + 1)

        walk(view.root, None, [], 0)
        self._write(["depth", "rank", "name", "parent", "path"], rows)

    # -- sequences ------------------------------------------------------------

    def seq_meta(self, species: str, result) -> None:
        header = ["accession", "organism", "length", "moltype", "biomol", "topology",
                  "completeness", "sourcedb", "updatedate"]
        self._write(header, [asdict(r) for r in result.records])
        if result.missing:
            self.warn(f"{len(result.missing)} accession(s) not returned by NCBI: "
                      + ", ".join(result.missing))

    def seq_records(self, species: str, result, output: str | None) -> None:
        # FASTA / GenBank are formats of their own: written as they are.
        if output:
            with open(output, "w", encoding="utf-8") as fh:
                fh.write(result.text)
        else:
            sys.stdout.write(result.text)
        if result.missing:
            self.warn(f"{len(result.missing)} accession(s) not returned by NCBI: "
                      + ", ".join(result.missing))

    def seq_aggregate(self, agg) -> None:
        row = {"name": agg.name, "rank": agg.rank, "species": agg.species,
               "isolates": agg.isolates, "accessions": agg.accessions, "refseq": agg.refseq}
        row.update({f"accessions_{k}": v for k, v in agg.moltype_breakdown.items()})
        self._write(list(row), [row])

    # -- ICTV Report ----------------------------------------------------------

    def text(self, chapter, markdown: str, figures=None, fig_width=None) -> None:
        self._write(["chapter", "url", "section", "level", "text"],
                    [(chapter.title, chapter.url, h, lvl, body)
                     for h, lvl, body in _sections(markdown)])

    def text_raw(self, markdown: str) -> None:
        sys.stdout.write(markdown if markdown.endswith("\n") else markdown + "\n")

    def tree(self, result, doc, others) -> None:
        rows = []
        for tip in doc.root.tip_labels():
            info = doc.tip_rows.get(tip, {})
            rows.append((result.family, doc.tree_id, tip, doc.display_name(tip),
                         info.get("accession"), info.get("species"), info.get("genus"),
                         tip in doc.matched))
        self._write(["family", "tree", "tip_label", "name", "accession", "species", "genus",
                     "matched"], rows)

    def tree_newick(self, newick: str) -> None:
        sys.stdout.write(newick.strip() + "\n")

    def tree_chapter(self, markdown: str) -> None:
        self._write(["section", "level", "text"], _sections(markdown))

    def msa(self, alignment, start: int, end: int, show_consensus: bool) -> None:
        from ..msa import window
        view = window(alignment, start, end, add_consensus=show_consensus)
        rows = []
        if view.consensus:
            rows.append(("consensus", view.consensus, False, view.start, view.n_cols))
        rows += [(r.name, r.seq, r.matched, view.start, view.n_cols) for r in view.rows]
        self._write(["name", "sequence", "matched", "start", "n_cols"], rows)

    def msa_fasta(self, alignment, start: int, end: int) -> None:
        sys.stdout.write("".join(f">{r.name}\n{r.seq[start - 1:end]}\n" for r in alignment.rows))

    # -- profile & annotate ---------------------------------------------------

    def profile(self, p, full: bool = False) -> None:
        self._write(list(profile_row(p)), [profile_row(p)])

    def annotate(self, annotations, output: str | None = None) -> None:
        from ..annotate import COLUMNS
        self._write(list(COLUMNS), [a.row() for a in annotations], output)

    def annotate_progress(self, i: int, n: int, ann) -> None:
        pass  # stdout carries only the table

    def annotate_summary(self, annotations, output: str | None) -> None:
        pass

    # -- utilities ------------------------------------------------------------

    def cache_info(self, infos) -> None:
        self._write(["namespace", "entries", "bytes"],
                    [(i.namespace, i.entries, i.bytes) for i in infos])

    def cache_cleared(self, removed: int, *, texts: bool, seqs: bool, images: bool) -> None:
        picked = [n for n, f in (("texts", texts), ("seqs", seqs), ("images", images)) if f]
        self._record({"cleared": removed, "scope": picked[0] if len(picked) == 1 else "all"})

    def config_view(self, view: dict) -> None:
        self._write(["key", "value"], list(_flatten(view).items()))

    def update_status(self, u, installed=None) -> None:
        row = {"current": u.current, "latest": u.latest, "latest_url": u.latest_url,
               "up_to_date": u.up_to_date}
        if installed is not None:
            row.update({"installed": installed.filename, "isolates": installed.isolates,
                        "species": installed.species})
        self._record(row)

    def vmr_installed(self, v) -> None:
        self._record({"installed": v.filename, "path": str(v.path), "isolates": v.isolates,
                      "species": v.species})

    def vmr_reset(self, removed: list[str], bundled: str) -> None:
        self._record({"removed": removed, "active": bundled})

    def trees_progress(self, i: int, n: int, entry: dict) -> None:
        pass

    def trees_installed(self, t) -> None:
        self._record({"path": str(t.path), "included": t.included, "omitted": t.omitted,
                      "rebuilt": t.families, "errors": len(t.errors)})

    def trees_reset(self, removed: bool) -> None:
        self._record({"removed": removed, "active": "bundled"})

    def viralzone_progress(self, i: int, n: int, entry: dict) -> None:
        pass

    def viralzone_installed(self, v) -> None:
        self._record({"path": str(v.path), "factsheets": v.factsheets, "taxa": v.taxa,
                      "failed": len(v.failed)})

    def viralzone_reset(self, removed: bool) -> None:
        self._record({"removed": removed, "active": "bundled"})

    def diagnose(self, d) -> None:
        self._record({"isolates": d.isolates, "accessions": d.accessions,
                      "empty_accession_rows": d.empty_accession_rows,
                      "unparsed_rows": len(d.unparsed)})


def profile_row(p) -> dict:
    """The taxon profile flattened to one row (shared by the CSV/TSV output)."""
    row: dict = {"query": p.query, "name": p.name, "rank": p.rank}
    lineage = p.ictv.lineage if p.ictv else {}
    row.update({rank: lineage.get(rank) for rank in RANKS})
    rep, seqs, vmr = p.ncbi_report, p.sequences, p.vmr
    row.update({
        "ncbi_taxid": p.ncbi.taxid if p.ncbi else None,
        "ncbi_name": p.ncbi.name if p.ncbi else None,
        "ictv_ncbi_agree": (all(r.agree is not False for r in p.comparison)
                            if p.ncbi else None),
        "ictv_status": rep.ictv_status if rep else None,
        "molecule": rep.genomic_moltype if rep else None,
        "species_count": vmr.species if vmr else None,
        "isolates": vmr.isolates if vmr else None,
        "genome": "; ".join(vmr.genome_compositions) if vmr else None,
        "host_source": "; ".join(vmr.host_sources) if vmr else None,
        "nucleotide_records": seqs.total if seqs else None,
        "complete_genomes": seqs.complete if seqs else None,
        "refseq_sequences": seqs.refseq if seqs else None,
        "assemblies": rep.assemblies if rep else None,
    })
    z = p.viralzone
    fs = z.factsheet if z else {}
    row.update({
        "viralzone_factsheet": fs.get("name"),
        "viralzone_url": fs.get("url"),
        "viralzone_host_groups": z.host_groups if z else None,
        "viralzone_host": fs.get("host"),
        "viralzone_geography": fs.get("geography"),
        "viralzone_diseases": fs.get("diseases"),
        "viralzone_transmission": fs.get("transmission"),
        "viralzone_vaccine": fs.get("vaccine"),
        "viralzone_antivirals": fs.get("antiviral_drugs"),
        "viralzone_replication": fs.get("replication_site"),
        "viralzone_credit": CREDIT["source"] if z else None,
    })
    return row
