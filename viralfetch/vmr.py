"""Load and index the VMR TSV.

The VMR is the local source of truth for taxonomy. One ships embedded in the
package; `viralfetch update` can install a newer release into the user data
dir, which is then preferred while it is newer than the embedded one. The TSV
is read-only and never modified; all normalisation happens here into in-memory
structures.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from pathlib import Path

from . import config as config_mod
from .accessions import is_parseable, parse_accessions
from .models import RANKS, Isolate, Taxon

# Basename of the embedded VMR release. Bump when a newer VMR is vendored.
VMR_FILENAME = "VMR_MSL41.v1.20260320.tsv"

# Mapping from our lowercase rank keys to the VMR TSV column headers.
_RANK_COLUMN = {
    "realm": "Realm",
    "subrealm": "Subrealm",
    "kingdom": "Kingdom",
    "subkingdom": "Subkingdom",
    "phylum": "Phylum",
    "subphylum": "Subphylum",
    "class": "Class",
    "subclass": "Subclass",
    "order": "Order",
    "suborder": "Suborder",
    "family": "Family",
    "subfamily": "Subfamily",
    "genus": "Genus",
    "subgenus": "Subgenus",
    "species": "Species",
}


@dataclass
class VMR:
    """An indexed, in-memory view of the VMR.

    Indices are keyed by ``name.casefold()`` for case-insensitive lookup.
    """

    isolates: list[Isolate]
    # name (casefolded) -> Taxon
    taxa: dict[str, Taxon]
    # rank -> list of taxon names (original casing), sorted, unique
    by_rank: dict[str, list[str]]
    # species name (casefolded) -> its isolates
    isolates_by_species: dict[str, list[Isolate]]
    # rows whose accession field produced zero accessions (quality signal)
    unparsed_rows: list[Isolate] = field(default_factory=list)
    empty_accession_rows: int = 0

    @property
    def species(self) -> list[str]:
        """Distinct species names present in the VMR."""
        return self.by_rank.get("species", [])

    def find(self, name: str) -> Taxon | None:
        return self.taxa.get(name.casefold())

    def suggest(self, name: str, limit: int = 5) -> list[str]:
        """Suggestions for a near-miss taxon name.

        Substring matches first (best for partial input), then fuzzy matches
        via ``difflib`` so trailing typos like "CoronaviridaX" still resolve.
        """
        import difflib

        needle = name.casefold()
        substring = [t.name for key, t in self.taxa.items() if needle in key]
        substring.sort(key=lambda n: (not n.casefold().startswith(needle), len(n)))

        out: list[str] = []
        seen: set[str] = set()
        for candidate in substring:
            key = candidate.casefold()
            if key not in seen:
                seen.add(key)
                out.append(candidate)
            if len(out) >= limit:
                return out

        for key in difflib.get_close_matches(needle, self.taxa.keys(), n=limit, cutoff=0.7):
            name_ = self.taxa[key].name
            if key not in seen:
                seen.add(key)
                out.append(name_)
        return out[:limit]


def release_key(name: str) -> tuple[int, int, int]:
    """Sortable (MSL, version, date) key parsed from a VMR filename.

    Filenames vary over releases — ``VMR_MSL39.v1_20240912``,
    ``VMR_MSL41.v1.20260320``, older ``VMR_MSL38_v1`` with no date — so each
    component is extracted independently and defaults to 0 when absent.
    """
    msl = re.search(r"MSL(\d+)", name)
    ver = re.search(r"[._]v(\d+)", name, re.I)
    date = re.search(r"(\d{8})", name)
    return (
        int(msl.group(1)) if msl else 0,
        int(ver.group(1)) if ver else 0,
        int(date.group(1)) if date else 0,
    )


def installed_dir() -> Path:
    """Where `viralfetch update` installs a downloaded VMR."""
    return config_mod.DATA_DIR / "vmr"


def installed_path() -> Path | None:
    """The VMR installed by `update`, if any (the newest, should several exist)."""
    candidates = list(installed_dir().glob("VMR_*.tsv"))
    return max(candidates, key=lambda p: release_key(p.name)) if candidates else None


def _data_path() -> Path:
    """Absolute path to the embedded VMR TSV."""
    return Path(resources.files("viralfetch").joinpath("data", VMR_FILENAME))


def active_path() -> Path:
    """The VMR in use: an installed one while it is newer than the embedded one.

    A package upgrade that vendors a newer VMR thus wins over a stale install.
    """
    installed = installed_path()
    if installed is not None and release_key(installed.name) > release_key(VMR_FILENAME):
        return installed
    return _data_path()


def active_filename() -> str:
    """Basename of the VMR in use (see :func:`active_path`)."""
    return active_path().name


def _lineage_from_row(row: dict[str, str]) -> dict[str, str]:
    lineage: dict[str, str] = {}
    for rank, col in _RANK_COLUMN.items():
        val = (row.get(col) or "").strip()
        if val:
            lineage[rank] = val
    return lineage


def _detect_rank(name: str, lineage: dict[str, str]) -> str:
    """The rank at which ``name`` sits in a lineage (defaults to lowest)."""
    for rank in reversed(RANKS):
        if lineage.get(rank) == name:
            return rank
    return "species"


def _load(path: Path) -> VMR:
    isolates: list[Isolate] = []
    taxa: dict[str, Taxon] = {}
    by_rank: dict[str, set[str]] = defaultdict(set)
    isolates_by_species: dict[str, list[Isolate]] = defaultdict(list)
    unparsed_rows: list[Isolate] = []
    empty_accession_rows = 0

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            lineage = _lineage_from_row(row)
            species = lineage.get("species", "")

            # Register every populated rank in this row as a Taxon.
            for rank, value in lineage.items():
                key = value.casefold()
                if key not in taxa:
                    # A taxon's lineage is everything at or above its own rank.
                    idx = RANKS.index(rank)
                    ancestors = {
                        r: lineage[r]
                        for r in RANKS[: idx + 1]
                        if r in lineage
                    }
                    taxa[key] = Taxon(name=value, rank=rank, lineage=ancestors)
                by_rank[rank].add(value)

            raw_acc = (row.get("Virus GENBANK accession") or "").strip()
            isolate = Isolate(
                isolate_id=(row.get("Isolate ID") or "").strip(),
                species=species,
                exemplar=(row.get("Exemplar or additional isolate") or "").strip() == "E",
                virus_names=(row.get("Virus name(s)") or "").strip(),
                abbreviations=(row.get("Virus name abbreviation(s)") or "").strip(),
                designation=(row.get("Virus isolate designation") or "").strip(),
                genome_composition=(row.get("Genome") or "").strip(),
                host_source=(row.get("Host source") or "").strip(),
                raw_accession=raw_acc,
                accessions=parse_accessions(raw_acc),
            )
            isolates.append(isolate)
            if species:
                isolates_by_species[species.casefold()].append(isolate)

            if not raw_acc:
                empty_accession_rows += 1
            elif not is_parseable(raw_acc):
                unparsed_rows.append(isolate)

    by_rank_sorted = {
        rank: sorted(names, key=str.casefold) for rank, names in by_rank.items()
    }

    return VMR(
        isolates=isolates,
        taxa=taxa,
        by_rank=by_rank_sorted,
        isolates_by_species=dict(isolates_by_species),
        unparsed_rows=unparsed_rows,
        empty_accession_rows=empty_accession_rows,
    )


@lru_cache(maxsize=1)
def load() -> VMR:
    """Load and index the active VMR (cached for the process lifetime)."""
    return _load(active_path())


def load_from(path: str | Path) -> VMR:
    """Load an arbitrary VMR TSV (used by tests with fixtures)."""
    return _load(Path(path))


if __name__ == "__main__":  # pragma: no cover - manual smoke check
    vmr = load()
    print(f"species: {len(vmr.species)}", file=sys.stderr)
    print(f"isolates: {len(vmr.isolates)}", file=sys.stderr)
    print(f"empty-accession rows: {vmr.empty_accession_rows}", file=sys.stderr)
    print(f"unparsed rows: {len(vmr.unparsed_rows)}", file=sys.stderr)
