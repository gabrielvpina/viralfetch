"""Multiple sequence alignments: load, window, and summarise a family's MSA.

Each family tree ships an aligned FASTA (``alignment.fasta``) beside its Newick.
This module loads one, maps its records to readable virus names, marks the
records the query points at, and windows the (often huge) column range down to a
viewport. Rendering — the colouring — is the render layer's job.

Design note (SPEC section 4): nothing here prints.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

GAP = "-"


@dataclass
class Row:
    """One aligned record: its display name, the gapped sequence, matched flag."""

    name: str
    seq: str
    matched: bool = False


@dataclass
class Alignment:
    """A family's multiple sequence alignment (one tree's ``alignment.fasta``)."""

    family: str
    tree_id: str
    molecule: str | None
    rows: list[Row]
    n_cols: int
    start: int = 1                       # 1-based first column of this window
    consensus: str | None = None
    total_cols: int = 0                  # columns before windowing
    matched_names: list[str] = field(default_factory=list)

    @property
    def n_rows(self) -> int:
        return len(self.rows)


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into ``[(header, sequence)]``, preserving order."""
    records: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            records.append((line[1:].strip(), []))
        elif records and line.strip():
            records[-1][1].append(line.strip())
    return [(h, "".join(parts)) for h, parts in records]


def _header_key(header: str) -> str:
    """Normalise a FASTA header to a tip_label (some are single-quoted)."""
    return header.strip().strip("'\"")


def load_alignment(doc, family: str) -> Alignment:
    """Build an :class:`Alignment` from a tree's ``alignment.fasta``.

    Records are labelled with their virus name (from the tree's ``members.tsv``)
    when known, and flagged when they are among the query's matched tips.
    """
    text = doc.align_path.read_text(encoding="utf-8")
    rows: list[Row] = []
    matched_names: list[str] = []
    for header, seq in parse_fasta(text):
        key = _header_key(header)
        row = doc.tip_rows.get(key)
        name = (row.get("name") if row else None) or _prettify(key)
        is_match = key in doc.matched
        rows.append(Row(name=name, seq=seq, matched=is_match))
        if is_match:
            matched_names.append(name)
    n_cols = max((len(r.seq) for r in rows), default=0)
    return Alignment(
        family=family,
        tree_id=doc.tree_id,
        molecule=doc.molecule,
        rows=rows,
        n_cols=n_cols,
        total_cols=n_cols,
        matched_names=matched_names,
    )


def _prettify(label: str) -> str:
    return label.replace("_ACCESSION_NOT_ON_SPREADSHEET", "").replace("_", " ").strip()


def consensus(rows: list[Row]) -> str:
    """Per-column majority residue (ties broken by first seen); gap if all gaps."""
    if not rows:
        return ""
    width = max(len(r.seq) for r in rows)
    out: list[str] = []
    for col in range(width):
        counts: Counter[str] = Counter()
        for r in rows:
            ch = r.seq[col] if col < len(r.seq) else GAP
            if ch != GAP:
                counts[ch] += 1
        out.append(counts.most_common(1)[0][0] if counts else GAP)
    return "".join(out)


def parse_range(spec: str, n_cols: int) -> tuple[int, int]:
    """Parse a 1-based inclusive column range like ``"100:180"`` or ``"100-180"``.

    Either bound may be omitted (``":180"``, ``"100:"``). Returns ``(start, end)``
    clamped to ``1..n_cols``. Raises ``ValueError`` on a malformed spec.
    """
    sep = ":" if ":" in spec else "-"
    lo, _, hi = spec.partition(sep)
    start = int(lo) if lo.strip() else 1
    end = int(hi) if hi.strip() else n_cols
    if start < 1 or end < start:
        raise ValueError(f"bad range {spec!r}")
    return max(1, start), min(n_cols, end)


def window(alignment: Alignment, start: int, end: int, *, add_consensus: bool) -> Alignment:
    """Return a copy of ``alignment`` sliced to columns ``start..end`` (1-based).

    ``add_consensus`` prepends a consensus row computed over the *window*.
    """
    lo, hi = start - 1, end  # to 0-based half-open
    rows = [Row(name=r.name, seq=r.seq[lo:hi], matched=r.matched) for r in alignment.rows]
    cons = consensus(rows) if add_consensus else None
    return Alignment(
        family=alignment.family,
        tree_id=alignment.tree_id,
        molecule=alignment.molecule,
        rows=rows,
        n_cols=end - start + 1,
        start=start,
        consensus=cons,
        total_cols=alignment.total_cols,
        matched_names=alignment.matched_names,
    )
