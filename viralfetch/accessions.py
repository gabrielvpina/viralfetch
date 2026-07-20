"""Parser for the VMR "Virus GENBANK accession" free-text field.

The column is not a clean identifier. Observed shapes (VMR MSL41):

    NC_045512                              -> simple
    RNA1: NC_003615; RNA2: NC_003616       -> segmented, labelled
    DNA-A: X15656; DNA-B: X15657           -> geminivirus
    KT732816; KT732815; KT732817           -> multiple, unlabelled
    AE006468 (2844298.2877981)             -> with coordinate annotation
    AB012345 (partial)                     -> with textual annotation
    CAJDJZ010000002                        -> WGS accession (long prefix)

The parser normalises every field into ``list[Accession]`` — one entry per
accession, which is what lets the rest of the tool work at accession
granularity (essential for segmented viruses, where each segment is a
separate GenBank record).
"""

from __future__ import annotations

import re

from .models import Accession

# A GenBank / RefSeq / WGS accession token. Covers:
#   - classic:  X15656, AB012345, EU623082
#   - RefSeq:   NC_045512, AC_123456
#   - WGS:      CAJDJZ010000002 (up to 6 letters + long digit run)
# Optionally versioned with a trailing ".N".
_ACCESSION_RE = re.compile(r"^[A-Za-z]{1,6}_?\d{5,}(?:\.\d+)?$")

# Parenthetical annotations to strip before extracting the token, e.g.
# "(2844298.2877981)" or "(partial)".
_PAREN_RE = re.compile(r"\([^)]*\)")


def parse_accessions(raw: str | None) -> list[Accession]:
    """Parse a raw VMR accession field into a list of :class:`Accession`.

    Unparseable chunks are dropped silently here; callers that need a
    quality signal should compare ``len(parse_accessions(raw))`` against the
    presence of a non-empty ``raw`` (see :func:`count_unparsed`).
    """
    if not raw:
        return []

    out: list[Accession] = []
    for chunk in raw.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue

        segment: str | None = None
        # A leading "LABEL:" marks a segment. Labels are arbitrary free text
        # ("RNA1", "DNA-A", "S", "partial", ...); split on the first colon.
        if ":" in chunk:
            label, _, rest = chunk.partition(":")
            segment = label.strip() or None
            chunk = rest.strip()

        # Drop coordinate / textual annotations before tokenising, then trim
        # trailing punctuation (VMR occasionally has a dangling "." with no
        # version, e.g. "PP173676.").
        token = _PAREN_RE.sub("", chunk).strip().rstrip(".")
        if not token:
            continue

        if _ACCESSION_RE.match(token):
            out.append(Accession(accession=token, segment=segment))

    return out


def is_parseable(raw: str | None) -> bool:
    """True if a non-empty field yielded at least one accession."""
    if not (raw and raw.strip()):
        return False
    return len(parse_accessions(raw)) > 0
