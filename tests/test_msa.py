"""Tests for the `msa` feature: alignment loading, windowing, and the command."""

import json

import pytest
from typer.testing import CliRunner

from viralfetch import msa, trees
from viralfetch.cli import app
from viralfetch.vmr import load

runner = CliRunner()


@pytest.fixture(scope="module")
def vmr():
    return load()


# -- parsing / helpers ----------------------------------------------------

def test_parse_fasta_multiline_and_order():
    text = ">a\nAC-\nGT\n>b\nACGGT\n"
    assert msa.parse_fasta(text) == [("a", "AC-GT"), ("b", "ACGGT")]


def test_consensus_majority_and_gaps():
    rows = [msa.Row("x", "AAC"), msa.Row("y", "A-C"), msa.Row("z", "TTC")]
    # col0: A,A,T -> A ; col1: A,-,T -> A ; col2: C,C,C -> C
    assert msa.consensus(rows) == "AAC"
    assert msa.consensus([msa.Row("x", "-"), msa.Row("y", "-")]) == "-"


def test_parse_range_variants():
    assert msa.parse_range("100:180", 316) == (100, 180)
    assert msa.parse_range("100-180", 316) == (100, 180)
    assert msa.parse_range(":50", 316) == (1, 50)
    assert msa.parse_range("300:", 316) == (300, 316)
    assert msa.parse_range("100:999", 316) == (100, 316)  # clamped
    for bad in ("5:2", "0:10", "abc"):
        with pytest.raises(ValueError):
            msa.parse_range(bad, 316)


# -- loading / windowing on real data -------------------------------------

def _corona_doc(vmr):
    return trees.resolve(vmr, "Betacoronavirus pandemicum").trees[0]


def test_load_alignment_maps_names_and_marks_matches(vmr):
    doc = _corona_doc(vmr)
    aln = msa.load_alignment(doc, "Coronaviridae")
    assert aln.n_rows == 55
    assert aln.total_cols == 316
    assert aln.matched_names  # the SARS coronaviruses
    assert any(r.matched for r in aln.rows)
    # names come from members.tsv, not the raw tip labels
    assert all("_ACCESSION_NOT_ON_SPREADSHEET" not in r.name for r in aln.rows)


def test_window_slices_columns_and_adds_consensus(vmr):
    aln = msa.load_alignment(_corona_doc(vmr), "Coronaviridae")
    view = msa.window(aln, 10, 25, add_consensus=True)
    assert view.start == 10
    assert view.n_cols == 16
    assert all(len(r.seq) == 16 for r in view.rows)
    assert view.consensus is not None and len(view.consensus) == 16
    assert view.total_cols == 316  # remembers the full width


# -- CLI ------------------------------------------------------------------

def test_msa_cli_json_windowed():
    result = runner.invoke(app, ["--json", "msa", "Betacoronavirus pandemicum", "--range", "1:12"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["family"] == "Coronaviridae"
    assert payload["n_cols"] == 12
    assert payload["total_cols"] == 316
    assert payload["matched"]
    assert all(len(row["seq"]) == 12 for row in payload["rows"])


def test_msa_cli_rich_renders_via_alv():
    result = runner.invoke(app, ["msa", "Coronaviridae", "--range", "1:40", "--consensus"])
    assert result.exit_code == 0
    assert "cols 1–40 of 316" in result.stdout   # header
    assert "consensus" in result.stdout          # consensus row drawn by alv


def test_msa_cli_fasta_is_raw_stdout():
    result = runner.invoke(app, ["msa", "Coronaviridae", "--range", "1:8", "--fasta"])
    assert result.exit_code == 0
    assert result.stdout.startswith(">")
    # every record is truncated to the window
    seqs = [ln for ln in result.stdout.splitlines() if not ln.startswith(">")]
    assert seqs and all(len(s) == 8 for s in seqs)


def test_msa_cli_bad_range_exits_2():
    result = runner.invoke(app, ["msa", "Coronaviridae", "--range", "9:2"])
    assert result.exit_code == 2


def test_msa_cli_family_without_alignment_exits_1():
    # Pleolipoviridae has a tree but no bundled alignment.fasta.
    result = runner.invoke(app, ["msa", "Pleolipoviridae"])
    assert result.exit_code == 1


def test_msa_cli_unknown_exits_1():
    result = runner.invoke(app, ["msa", "Notarealvirusxyz"])
    assert result.exit_code == 1
