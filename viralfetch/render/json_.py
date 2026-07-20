"""JSON renderer — pure JSON on stdout, no decoration, ready for ``jq``.

Warnings and progress must go to stderr in this mode; only the payload prints
to stdout. Errors also go to stderr (see :func:`error`).
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict

from ..compare import LineageComparison
from ..models import Chapter
from ..ncbi import MetaResult, RecordsResult
from ..queries import MembersView, TaxonTreeNode, TaxonView, TreeView
from ..sequences import TaxonAggregate


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


def seq_meta(species: str, result: MetaResult) -> None:
    _emit({
        "species": species,
        "records": [asdict(r) for r in result.records],
        "missing": result.missing,
    })


def seq_records(species: str, result: RecordsResult, output: str | None) -> None:
    if output:
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(result.text)
    else:
        sys.stdout.write(result.text)
    # Summary always to stderr so stdout stays pure data (file or records).
    summary = {
        "species": species,
        "rettype": result.rettype,
        "returned": len(result.returned),
        "missing": result.missing,
    }
    if output:
        summary["written"] = output
    print(json.dumps(summary, ensure_ascii=False), file=sys.stderr)


def seq_aggregate(agg: TaxonAggregate) -> None:
    _emit({
        "name": agg.name,
        "rank": agg.rank,
        "species": agg.species,
        "isolates": agg.isolates,
        "accessions": agg.accessions,
        "refseq": agg.refseq,
        "moltype_breakdown": agg.moltype_breakdown,
    })


def compare(cmp: LineageComparison) -> None:
    ncbi = cmp.ncbi
    _emit({
        "taxon": cmp.taxon,
        "representative_accession": cmp.representative_accession,
        "ictv": [{"rank": r, "name": n} for r, n in cmp.ictv],
        "ncbi": {
            "taxid": ncbi.taxid,
            "name": ncbi.name,
            "rank": ncbi.rank,
            "lineage": [{"rank": r, "name": n} for r, n in ncbi.lineage],
        } if ncbi else None,
    })


def text(chapter: Chapter, markdown: str) -> None:
    _emit({
        "slug": chapter.slug,
        "title": chapter.title,
        "url": chapter.url,
        "doi": chapter.doi,
        "markdown": markdown,
    })


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
