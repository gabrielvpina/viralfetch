"""Output layer — the *only* place that writes to stdout.

Business logic returns view dataclasses; the CLI hands them here. Selecting a
renderer by format is what makes ``--json`` work without duplicating logic:

    renderer = render.get(config.format)
    renderer.tax(view)

Both backends (:mod:`rich_` and :mod:`json_`) expose the same function names.
"""

from __future__ import annotations

from . import json_, rich_


def get(fmt: str):
    """Return the renderer module for ``fmt`` ("json" or "rich")."""
    return json_ if fmt == "json" else rich_
