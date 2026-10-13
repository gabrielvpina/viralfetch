"""Output layer — the *only* place that writes to stdout.

Business logic returns view dataclasses; the CLI hands them here. Selecting a
renderer by format is what makes ``--json`` work without duplicating logic:

    renderer = render.get(config.format)
    renderer.tax(view)

All backends (:mod:`rich_`, :mod:`json_` and the CSV / TSV
:class:`table_.TableRenderer`) expose the same function names.
"""

from __future__ import annotations

from . import json_, rich_
from .table_ import TableRenderer

_TABLES = {"csv": TableRenderer(","), "tsv": TableRenderer("\t")}


def get(fmt: str):
    """The renderer for ``fmt``: "rich", "json", "csv" or "tsv"."""
    if fmt == "json":
        return json_
    return _TABLES.get(fmt, rich_)
