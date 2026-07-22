"""Tests for the `tree` feature: Newick parsing, resolution, and rendering."""

import json

import pytest
from typer.testing import CliRunner

from viralfetch import trees
from viralfetch.cli import app
from viralfetch.render import rich_
from viralfetch.vmr import load

runner = CliRunner()


@pytest.fixture(scope="module")
def vmr():
    return load()


# -- Newick parser --------------------------------------------------------

def test_parse_newick_structure_and_values():
    root = trees.parse_newick("((A:0.1,B:0.2)0.9:0.3,C:0.4);")
    assert not root.is_tip
    assert root.tip_labels() == ["A", "B", "C"]
    internal = root.children[0]
    assert internal.support == "0.9"          # support kept on internal nodes
    assert internal.length == 0.3
    assert internal.children[0].name == "A"
    assert internal.children[0].length == 0.1


def test_parse_newick_tolerates_labels_with_dots_and_underscores():
    root = trees.parse_newick("(AF086833_EBOV:0.30,FJ217161_BDBV:0.25);")
    assert root.tip_labels() == ["AF086833_EBOV", "FJ217161_BDBV"]


# -- resolution -----------------------------------------------------------

def test_resolve_species_highlights_its_tips(vmr):
    r = trees.resolve(vmr, "Betacoronavirus pandemicum")
    assert r.family == "Coronaviridae"
    assert r.source == "vmr"
    assert r.has_trees
    doc = r.trees[0]
    assert doc.matched
    for tip in doc.matched:
        assert doc.tip_rows[tip]["species"] == "Betacoronavirus pandemicum"


def test_resolve_genus_highlights_whole_clade(vmr):
    doc = trees.resolve(vmr, "Betacoronavirus").trees[0]
    assert doc.matched
    assert all(doc.tip_rows[t]["genus"] == "Betacoronavirus" for t in doc.matched)


def test_resolve_family_shows_tree_without_highlighting(vmr):
    r = trees.resolve(vmr, "Coronaviridae")
    assert r.note is None
    assert all(not doc.matched for doc in r.trees)


def test_member_search_fallback_for_non_taxon_name(vmr):
    # A tip's virus *name* that is not an ICTV taxon: the VMR misses it, so the
    # member scan finds it.
    r = trees.resolve(vmr, "porcine epidemic diarrhea virus")
    assert r.source == "member"
    assert r.family == "Coronaviridae"
    assert any(doc.matched for doc in r.trees)


def test_unknown_name_raises_with_suggestions(vmr):
    with pytest.raises(trees.TreesNotFound) as exc:
        trees.resolve(vmr, "Notarealvirusxyz")
    assert isinstance(exc.value.suggestions, list)


def test_family_without_bundled_tree_has_no_trees(vmr):
    r = trees.resolve(vmr, "Ahmunviridae")  # in VMR/index but omitted (no resources)
    assert r.family == "Ahmunviridae"
    assert not r.has_trees


def test_family_may_carry_several_trees(vmr):
    r = trees.resolve(vmr, "Coronaviridae")
    assert len(r.trees) >= 2  # RdRp + helicase


# -- rendering ------------------------------------------------------------

def _render(renderable, width=100):
    import io
    from rich.console import Console
    buf = io.StringIO()
    Console(file=buf, width=width, force_terminal=False).print(renderable)
    return buf.getvalue()


def test_ascii_tree_draws_branches_and_all_tips(vmr):
    doc = trees.resolve(vmr, "Coronaviridae").trees[0]
    out = _render(rich_._ascii_tree(doc, 100))
    assert any(ch in out for ch in "─│├┌└")          # box-drawing branches
    # every tip's display name is present, one per line
    for tip in doc.tip_rows.values():
        assert tip["name"] in out


def test_ascii_tree_marks_the_query(vmr):
    doc = trees.resolve(vmr, "Betacoronavirus pandemicum").trees[0]
    out = _render(rich_._ascii_tree(doc, 100))
    assert "← match" in out


# -- CLI ------------------------------------------------------------------

def test_tree_cli_json_payload():
    result = runner.invoke(app, ["--json", "tree", "Betacoronavirus pandemicum"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["family"] == "Coronaviridae"
    assert payload["tree"]["newick"].startswith("(")
    assert payload["tree"]["matched"]
    assert payload["other_trees"]  # a second tree is listed


def test_tree_cli_newick_is_raw_stdout():
    result = runner.invoke(app, ["tree", "Filoviridae", "--newick"])
    assert result.exit_code == 0
    assert result.stdout.strip().startswith("(")
    assert result.stdout.strip().endswith(";")


def test_tree_cli_unknown_exits_1():
    result = runner.invoke(app, ["tree", "Notarealvirusxyz"])
    assert result.exit_code == 1


def test_tree_cli_out_of_range_exits_2():
    result = runner.invoke(app, ["tree", "Filoviridae", "--tree", "9"])
    assert result.exit_code == 2
