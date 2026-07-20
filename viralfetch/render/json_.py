"""JSON renderer — pure JSON on stdout, no decoration, ready for ``jq``.

Warnings and progress must go to stderr in this mode; only the payload prints
to stdout. Errors also go to stderr (see :func:`error`).
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict

from ..queries import MembersView, TaxonTreeNode, TaxonView, TreeView


def _emit(payload) -> None:
    json.dump(payload, sys.stdout, ensure_ascii=False, indent=2)
    sys.stdout.write("\n")


def tax(view: TaxonView) -> None:
    payload = {
        "name": view.taxon.name,
        "rank": view.taxon.rank,
        "lineage": view.taxon.lineage,
    }
    if view.isolate_summary is not None:
        payload["isolates"] = asdict(view.isolate_summary)
    _emit(payload)


def members(view: MembersView) -> None:
    payload: dict = {
        "parent": {"name": view.parent.name, "rank": view.parent.rank},
    }
    if view.rank is not None and not view.count_only:
        payload["rank"] = view.rank
        payload["members"] = [
            {"name": m.name, "rank": m.rank, "species_count": m.species_count}
            for m in view.members
        ]
    else:
        payload["breakdown"] = view.breakdown
    _emit(payload)


def _node_to_dict(node: TaxonTreeNode) -> dict:
    return {
        "name": node.name,
        "rank": node.rank,
        "children": [_node_to_dict(c) for c in node.children],
    }


def members_tree(view: TreeView) -> None:
    payload = {"total": view.total, "tree": _node_to_dict(view.root)}
    _emit(payload)


def not_found(name: str, suggestions: list[str]) -> None:
    print(
        json.dumps(
            {"error": "taxon_not_found", "name": name, "suggestions": suggestions},
            ensure_ascii=False,
        ),
        file=sys.stderr,
    )


def error(message: str) -> None:
    print(json.dumps({"error": message}, ensure_ascii=False), file=sys.stderr)


def warn(message: str) -> None:
    print(json.dumps({"warning": message}, ensure_ascii=False), file=sys.stderr)
