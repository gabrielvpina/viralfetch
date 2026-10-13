"""Read the ViralZone data set.

ViralZone (https://viralzone.expasy.org) is the SIB Swiss Institute of
Bioinformatics' knowledge base of curated virus factsheets. A snapshot of its
factsheets ships with the package (``data/viralzone.json``); ``viralfetch
update --viralzone`` can install a fresher one into the user data dir, which is
then preferred while it is at least as new as the bundled one (the same rule as
the VMR and the trees).

Factsheets are keyed by ViralZone page id and carry the NCBI taxid of their
taxon, which is the most reliable join key: ViralZone may follow an older ICTV
release than the VMR in use, so names can differ.

Anything shown from this data set must credit ViralZone (:data:`CREDIT`).

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from pathlib import Path

from . import config as config_mod

BASE_URL = "https://viralzone.expasy.org"
DATA_FILENAME = "viralzone.json"

CREDIT = {
    "source": "ViralZone, SIB Swiss Institute of Bioinformatics",
    "source_url": BASE_URL,
    "citation": (
        "De Castro E, Hulo C, Masson P, Auchincloss A, Bridge A, Le Mercier P. "
        "ViralZone 2024 provides higher-resolution images and advanced "
        "virus-specific resources. Nucleic Acids Res. 2024;52(D1):D817-D821. "
        "doi:10.1093/nar/gkad946"
    ),
    "license": (
        "ViralZone images are licensed CC BY 4.0 "
        "(https://creativecommons.org/licenses/by/4.0/); data are provided by "
        "ViralZone for research, educational and informational purposes, not "
        "for medical use."
    ),
}

# One-line credit for space-constrained displays (the full citation is in CREDIT).
CREDIT_SHORT = "ViralZone (SIB) · CC BY 4.0 · doi:10.1093/nar/gkad946"

# Most specific first: the factsheet closest to a query wins.
_RANK_ORDER = ("species", "virus", "genus", "subfamily", "family")


def bundled_path() -> Path:
    return Path(resources.files("viralfetch").joinpath("data", DATA_FILENAME))


def installed_path() -> Path:
    """Where ``viralfetch update --viralzone`` installs a rebuilt data set."""
    return config_mod.DATA_DIR / DATA_FILENAME


def _generated_at(path: Path) -> str:
    try:
        with path.open(encoding="utf-8") as fh:
            return json.load(fh).get("generated_at", "")
    except (OSError, json.JSONDecodeError):
        return ""


def active_path() -> Path:
    installed = installed_path()
    stamp = _generated_at(installed)
    if stamp and stamp >= _generated_at(bundled_path()):
        return installed
    return bundled_path()


@dataclass
class ViralZone:
    """An indexed, in-memory view of the ViralZone data set."""

    generated_at: str
    factsheets: dict[int, dict]
    taxa: list[dict]
    by_taxid: dict[int, dict] = field(default_factory=dict)
    by_name: dict[str, dict] = field(default_factory=dict)
    # species / virus names listed on a factsheet (casefolded) -> factsheets
    by_species: dict[str, list[dict]] = field(default_factory=dict)
    # genus/subfamily/family names from the realm pages (casefolded) -> taxon
    taxa_by_name: dict[str, dict] = field(default_factory=dict)

    def get(self, page_id: int) -> dict | None:
        return self.factsheets.get(page_id)

    def find(self, name: str) -> dict | None:
        """The factsheet titled ``name`` (case-insensitive), else ``None``."""
        return self.by_name.get(name.casefold())

    def for_taxid(self, taxid: int | str) -> dict | None:
        try:
            return self.by_taxid.get(int(taxid))
        except (TypeError, ValueError):
            return None

    def taxon(self, name: str) -> dict | None:
        """A taxon as laid out on ViralZone's realm pages (with host groups),
        whether or not it has a factsheet."""
        return self.taxa_by_name.get(name.casefold())

    def match(
        self,
        *,
        taxids: list[int | str] | None = None,
        names: list[str] | None = None,
        lineage: dict[str, str] | None = None,
    ) -> list[dict]:
        """Factsheets relevant to a taxon, most specific first, without repeats.

        ``taxids`` (the taxon and its NCBI ancestors) are tried first, then
        ``names`` (factsheet titles or species they list), then the ICTV
        ``lineage`` (``{"species": ..., "genus": ..., "family": ...}``).
        """
        found: dict[int, dict] = {}

        def add(fs: dict | None) -> None:
            if fs is not None:
                found.setdefault(fs["id"], fs)

        for t in taxids or []:
            add(self.for_taxid(t))
        for n in names or []:
            add(self.find(n))
            for fs in self.by_species.get(n.casefold(), []):
                add(fs)
        for rank in ("species", "genus", "subfamily", "family"):
            value = (lineage or {}).get(rank)
            if not value:
                continue
            add(self.find(value))
            if rank == "species":
                for fs in self.by_species.get(value.casefold(), []):
                    add(fs)
        order = {r: i for i, r in enumerate(_RANK_ORDER)}
        return sorted(found.values(), key=lambda fs: order.get(fs.get("rank"), len(order)))


def load_from(path: Path) -> ViralZone:
    with Path(path).open(encoding="utf-8") as fh:
        data = json.load(fh)
    factsheets = {int(k): v for k, v in data.get("factsheets", {}).items()}
    vz = ViralZone(
        generated_at=data.get("generated_at", ""),
        factsheets=factsheets,
        taxa=data.get("taxa", []),
    )
    for fs in factsheets.values():
        if fs.get("taxid") is not None:
            vz.by_taxid.setdefault(int(fs["taxid"]), fs)
        vz.by_name.setdefault(fs["name"].casefold(), fs)
        for sp in fs.get("species", []):
            vz.by_species.setdefault(sp["name"].casefold(), []).append(fs)
    for taxon in vz.taxa:
        vz.taxa_by_name.setdefault(taxon["name"].casefold(), taxon)
    return vz


@lru_cache(maxsize=1)
def load() -> ViralZone | None:
    """The active data set, or ``None`` when none is bundled or installed."""
    path = active_path()
    if not path.is_file():
        return None
    return load_from(path)
