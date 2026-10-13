"""Host and genome filters for taxon listings (``members --host / --genome``).

Both work on VMR isolate rows, so every count downstream (species per genus,
per-rank breakdown, the descendant tree) reflects the filter.

- **Genome** matches the VMR "Genome" column. A bare molecule (``ssRNA``,
  ``dsDNA``) covers its polarities (``ssRNA(+)``, ``ssRNA(-)``, ``ssRNA(+/-)``)
  but not the reverse-transcribing groups (``ssRNA-RT``), which are named on
  their own; ``RNA`` / ``DNA`` / ``RT`` are broad classes.
- **Host** matches the VMR "Host source" column, per species (bacteria,
  vertebrates, plants, soil, ...). The VMR has no "human" (or "eukaryotic
  microorganism") category; those come from ViralZone's host groups, which are
  recorded per genus — ``--host human`` keeps the genera with human-infecting
  members.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable

from .models import Isolate


class InvalidFilter(Exception):
    """Raised for a --host / --genome value no data set knows."""

    def __init__(self, option: str, value: str, valid: list[str]):
        self.option = option
        self.value = value
        self.valid = valid
        super().__init__(f"unknown {option} {value!r}")

# Host categories of the VMR "Host source" column ("(S)" = sampled source).
VMR_HOSTS = (
    "bacteria", "archaea", "vertebrates", "invertebrates", "plants", "fungi",
    "protists", "algae", "soil", "freshwater", "marine", "sewage", "air",
    "phytobiome", "unknown",
)
# Host groups only ViralZone records (per genus).
VIRALZONE_HOSTS = ("human", "micro")
_HOST_ALIASES = {
    "humans": "human", "homo sapiens": "human",
    "bacterium": "bacteria", "bacterial": "bacteria", "phage": "bacteria", "phages": "bacteria",
    "archaeon": "archaea", "archaeal": "archaea",
    "fungus": "fungi", "fungal": "fungi",
    "protist": "protists", "protozoa": "protists", "protozoan": "protists",
    "alga": "algae",
    "eukaryotic microorganisms": "micro", "microorganisms": "micro",
}
_BROAD_GENOMES = {"rna", "dna", "rt"}


def normalize_host(value: str) -> str:
    """The canonical host term for ``value``, or :class:`InvalidFilter`."""
    term = " ".join(value.casefold().split())
    term = _HOST_ALIASES.get(term, term)
    if term not in VMR_HOSTS and term + "s" in VMR_HOSTS:  # "plant" -> "plants"
        term += "s"
    if term not in VMR_HOSTS and term not in VIRALZONE_HOSTS:
        raise InvalidFilter("host", value, list(VMR_HOSTS + VIRALZONE_HOSTS))
    return term


def isolate_hosts(iso: Isolate) -> set[str]:
    """The VMR host terms of an isolate ("invertebrates, vertebrates (S)")."""
    return {
        re.sub(r"\s*\(s\)$", "", t.strip().casefold())
        for t in iso.host_source.split(",") if t.strip()
    }


def _genome_key(value: str) -> str:
    key = re.sub(r"\s+", "", value).casefold()
    for sign in ("+/-", "+", "-"):  # "ssRNA+" / "+ssRNA" -> "ssrna(+)"
        if key.endswith(sign) and not key.endswith("-rt") and "(" not in key:
            return f"{key[:-len(sign)]}({sign})"
        if key.startswith(sign) and "(" not in key:
            return f"{key[len(sign):]}({sign})"
    return key


def genome_matcher(value: str, known: Iterable[str]) -> Callable[[str], bool]:
    """A predicate over VMR genome compositions (which may list several parts,
    "ssRNA(-); ssRNA(+/-)"); :class:`InvalidFilter` if it matches none."""
    query = _genome_key(value)

    def part_matches(part: str) -> bool:
        part = _genome_key(part)
        if query in _BROAD_GENOMES:
            return part.endswith("-rt") if query == "rt" else query in part
        if part == query:
            return True
        return "(" not in query and part.split("(")[0] == query

    def matches(composition: str) -> bool:
        return any(part_matches(p) for p in composition.split(";") if p.strip())

    known = sorted({k for k in known if k})
    if not any(matches(k) for k in known):
        parts = sorted({p.strip() for k in known for p in k.split(";")})
        raise InvalidFilter("genome", value, parts + ["RNA", "DNA", "RT"])
    return matches


def filter_rows(
    rows: list[Isolate],
    lineage_of: Callable[[Isolate], dict[str, str]],
    *,
    host: str | None = None,
    genome: str | None = None,
    known_genomes: Iterable[str] = (),
    viralzone=None,
) -> list[Isolate]:
    """``rows`` restricted to a host and/or genome (see the module docstring).

    ``viralzone`` (a :class:`viralfetch.viralzone.ViralZone`) is needed only
    for the ViralZone-only host groups.
    """
    if genome:
        matches = genome_matcher(genome, known_genomes)
        rows = [r for r in rows if matches(r.genome_composition)]
    if host:
        term = normalize_host(host)
        if term in VIRALZONE_HOSTS:
            if viralzone is None:
                return []
            genera: dict[str, bool] = {}

            def keep(iso: Isolate) -> bool:
                genus = lineage_of(iso).get("genus", "")
                if genus not in genera:
                    taxon = viralzone.taxon(genus) if genus else None
                    genera[genus] = bool(taxon and term in taxon["hosts"])
                return genera[genus]

            rows = [r for r in rows if keep(r)]
        else:
            rows = [r for r in rows if term in isolate_hosts(r)]
    return rows
