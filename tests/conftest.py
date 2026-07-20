"""Shared test fixtures: a tiny in-memory VMR built from concise rows."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from viralfetch.vmr import VMR, load_from

# Full VMR column header (order matters for DictReader round-trips).
HEADER = [
    "Isolate ID", "Species Sort", "Isolate Sort", "Realm", "Subrealm",
    "Kingdom", "Subkingdom", "Phylum", "Subphylum", "Class", "Subclass",
    "Order", "Suborder", "Family", "Subfamily", "Genus", "Subgenus",
    "Species", "ICTV_ID", "Exemplar or additional isolate", "Virus name(s)",
    "Virus name abbreviation(s)", "Virus isolate designation",
    "Virus GENBANK accession", "Genome coverage", "Genome", "Host source",
    "Accessions Link",
]


def _row(**cells: str) -> dict[str, str]:
    row = {col: "" for col in HEADER}
    row.update(cells)
    return row


# A minimal but structurally realistic taxonomy.
_ROWS = [
    _row(**{
        "Isolate ID": "I1", "Realm": "Testviria", "Kingdom": "Testvirae",
        "Phylum": "Testphy", "Class": "Testcla", "Order": "Testord",
        "Family": "Alphaviridae", "Genus": "Alphavirus",
        "Species": "Alphavirus one",
        "Exemplar or additional isolate": "E",
        "Virus GENBANK accession": "NC_000001", "Genome": "ssRNA(+)",
    }),
    _row(**{
        "Isolate ID": "I2", "Realm": "Testviria", "Kingdom": "Testvirae",
        "Phylum": "Testphy", "Class": "Testcla", "Order": "Testord",
        "Family": "Alphaviridae", "Genus": "Alphavirus",
        "Species": "Alphavirus one",
        "Exemplar or additional isolate": "A",
        "Virus GENBANK accession": "NC_000002", "Genome": "ssRNA(+)",
    }),
    _row(**{
        "Isolate ID": "I3", "Realm": "Testviria", "Kingdom": "Testvirae",
        "Phylum": "Testphy", "Class": "Testcla", "Order": "Testord",
        "Family": "Alphaviridae", "Genus": "Alphavirus",
        "Species": "Alphavirus two",
        "Exemplar or additional isolate": "E",
        "Virus GENBANK accession": "RNA1: NC_000003; RNA2: NC_000004",
        "Genome": "ssRNA(+)",
    }),
    _row(**{
        "Isolate ID": "I4", "Realm": "Testviria", "Kingdom": "Testvirae",
        "Phylum": "Testphy", "Class": "Testcla", "Order": "Testord",
        "Family": "Alphaviridae", "Genus": "Betavirus",
        "Species": "Betavirus one",
        "Exemplar or additional isolate": "E",
        "Virus GENBANK accession": "NC_000005", "Genome": "dsDNA",
    }),
    _row(**{
        "Isolate ID": "I5", "Realm": "Testviria", "Kingdom": "Testvirae",
        "Phylum": "Testphy", "Class": "Testcla", "Order": "Testord",
        "Family": "Gammaviridae", "Genus": "Gammavirus",
        "Species": "Gammavirus one",
        "Exemplar or additional isolate": "E",
        "Virus GENBANK accession": "",  # empty accession row
        "Genome": "ssRNA(-)",
    }),
]


@pytest.fixture
def mini_vmr_path(tmp_path: Path) -> Path:
    path = tmp_path / "mini_vmr.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(_ROWS)
    return path


@pytest.fixture
def vmr(mini_vmr_path: Path) -> VMR:
    return load_from(mini_vmr_path)
