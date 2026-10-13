"""members --host/--genome, `annotate`, and the CSV/TSV output of every command."""

import csv
import io
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from viralfetch import annotate as annotate_mod
from viralfetch import compare, filters, queries
from viralfetch import viralzone as vz_mod
from viralfetch.cli import app
from viralfetch.ncbi import NcbiLineage, NCBIError
from viralfetch.render.table_ import TableRenderer, _cell, _sections
from viralfetch.vmr import load_from

from .conftest import HEADER, _row

runner = CliRunner()

_BASE = {"Realm": "Testviria", "Kingdom": "Testvirae", "Phylum": "Testphy",
         "Class": "Testcla", "Order": "Testord", "Exemplar or additional isolate": "E"}
_ROWS = [
    _row(**_BASE, **{"Isolate ID": "1", "Family": "Alphaviridae", "Genus": "Alphavirus",
                     "Species": "Alphavirus one", "Genome": "ssRNA(+)",
                     "Host source": "invertebrates, vertebrates",
                     "Virus GENBANK accession": "NC_000001"}),
    _row(**_BASE, **{"Isolate ID": "2", "Family": "Alphaviridae", "Genus": "Alphavirus",
                     "Species": "Alphavirus two", "Genome": "ssRNA(-); ssRNA(+/-)",
                     "Host source": "plants (S)",
                     "Virus GENBANK accession": "RNA1: NC_000003; RNA2: NC_000004"}),
    _row(**_BASE, **{"Isolate ID": "3", "Family": "Alphaviridae", "Genus": "Betavirus",
                     "Species": "Betavirus one", "Genome": "ssRNA-RT",
                     "Host source": "vertebrates", "Virus GENBANK accession": "NC_000005"}),
    _row(**_BASE, **{"Isolate ID": "4", "Family": "Gammaviridae", "Genus": "Gammavirus",
                     "Species": "Gammavirus one", "Genome": "dsDNA",
                     "Host source": "bacteria", "Virus GENBANK accession": "NC_000006"}),
]


@pytest.fixture
def hvmr(tmp_path):
    path = tmp_path / "vmr.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(_ROWS)
    return load_from(path)


def _vz(**hosts_by_genus) -> vz_mod.ViralZone:
    vz = vz_mod.ViralZone(generated_at="x", factsheets={}, taxa=[
        {"name": g, "rank": "genus", "id": None, "family": None, "subfamily": None,
         "hosts": h, "page": 1} for g, h in hosts_by_genus.items()])
    vz.taxa_by_name = {t["name"].casefold(): t for t in vz.taxa}
    return vz


def _genera(view) -> dict[str, int]:
    return {m.name: m.species_count for m in view.members}


# -- filters ----------------------------------------------------------------------

@pytest.mark.parametrize("value, expected", [
    ("vertebrates", {"Alphavirus": 1, "Betavirus": 1}),
    ("Vertebrate", {"Alphavirus": 1, "Betavirus": 1}),   # singular, any case
    ("plant", {"Alphavirus": 1}),                          # "(S)" source still counts
    ("invertebrates", {"Alphavirus": 1}),
    ("bacteria", {}),
])
def test_host_filter_vmr_terms(hvmr, value, expected):
    assert _genera(queries.members(hvmr, "Alphaviridae", "genus", host=value)) == expected


def test_host_filter_human_uses_viralzone_genera(hvmr):
    vz = _vz(Alphavirus=["human", "vertebrate"], Betavirus=["vertebrate"])
    view = queries.members(hvmr, "Alphaviridae", "genus", host="humans", viralzone=vz)
    assert _genera(view) == {"Alphavirus": 2}
    assert view.filters == {"host": "humans"}
    assert queries.members(hvmr, "Alphaviridae", "genus", host="human").members == []  # no data


@pytest.mark.parametrize("value, expected", [
    ("ssRNA", {"Alphavirus": 2}),             # every polarity, not ssRNA-RT
    ("ssRNA(+)", {"Alphavirus": 1}),
    ("+ssRNA", {"Alphavirus": 1}),
    ("ssRNA(+/-)", {"Alphavirus": 1}),        # one part of a multi-part composition
    ("ssRNA-RT", {"Betavirus": 1}),
    ("RT", {"Betavirus": 1}),
    ("RNA", {"Alphavirus": 2, "Betavirus": 1}),
])
def test_genome_filter(hvmr, value, expected):
    assert _genera(queries.members(hvmr, "Alphaviridae", "genus", genome=value)) == expected


def test_combined_filters_and_breakdown(hvmr):
    view = queries.members(hvmr, "Testviria", host="vertebrates", genome="ssRNA")
    assert view.breakdown == {"kingdom": 1, "phylum": 1, "class": 1, "order": 1,
                              "family": 1, "genus": 1, "species": 1}
    tree = queries.members_tree(hvmr, "Alphaviridae", genome="dsDNA")
    assert tree.total == 0 and tree.filters == {"genome": "dsDNA"}


@pytest.mark.parametrize("kwargs, option", [({"host": "martians"}, "host"),
                                            ({"genome": "xRNA"}, "genome")])
def test_invalid_filters(hvmr, kwargs, option):
    with pytest.raises(filters.InvalidFilter) as exc:
        queries.members(hvmr, "Alphaviridae", "genus", **kwargs)
    assert exc.value.option == option and exc.value.valid


def test_members_filter_cli():
    result = runner.invoke(app, ["--json", "members", "Flaviviridae", "--rank", "genus",
                                 "--host", "human"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.stdout)
    assert payload["filters"] == {"host": "human"}
    assert "Orthoflavivirus" in {m["name"] for m in payload["members"]}
    bad = runner.invoke(app, ["members", "Flaviviridae", "--genome", "xRNA"])
    assert bad.exit_code == 2 and "Unknown --genome" in bad.stderr


# -- annotate -------------------------------------------------------------------

class FakeNCBI:
    def __init__(self, fail: bool = False):
        self.fail = fail

    def esearch_taxid(self, name):
        if self.fail:
            raise NCBIError("down")
        return {"alphavirus one": "101", "alphaviridae": "100"}.get(name.casefold())


NCBI = {
    "old alpha name": NcbiLineage("555", "Old alpha name", "no rank", [
        ("family", "Alphaviridae"), ("genus", "Alphavirus"), ("species", "Alphavirus two"),
        ("no rank", "Old alpha name")]),
    "far away virus": NcbiLineage("777", "Far away virus", "species", [
        ("realm", "Otherviria"), ("family", "Farviridae"), ("species", "Far away virus")]),
}


@pytest.fixture
def ncbi_lookup(monkeypatch):
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda c, n: NCBI.get(n.casefold()))


def test_read_names_skips_blanks_and_comments():
    lines = ["Alphavirus one\n", "  \n", "# header\n", '"Betavirus" # quoted, commented\n']
    assert annotate_mod.read_names(lines) == ["Alphavirus one", "Betavirus"]


def test_annotate_rows(hvmr, ncbi_lookup):
    names = ["alphavirus one", "Alphavirus two", "Alphaviridae", "Old alpha name",
             "Far away virus", "Alphavirus oen", "ALPHAVIRUS ONE"]
    anns = annotate_mod.annotate(hvmr, names, ncbi=FakeNCBI())
    rows = [a.row() for a in anns]
    assert [r["status"] for r in rows] == ["ictv", "ictv", "ictv", "ncbi", "ncbi_only",
                                           "not_found", "ictv"]
    sp = rows[0]
    assert (sp["name"], sp["rank"], sp["family"], sp["genus"], sp["species"]) == (
        "Alphavirus one", "species", "Alphaviridae", "Alphavirus", "Alphavirus one")
    assert sp["ncbi_taxid"] == "101" and sp["ncbi_taxonomy_url"].endswith("id=101")
    assert sp["accession"] == "NC_000001"
    assert sp["accession_url"] == "https://www.ncbi.nlm.nih.gov/nuccore/NC_000001"
    seg = rows[1]
    assert seg["accession"] == "NC_000003;NC_000004"
    assert seg["accession_url"].endswith("/nuccore/NC_000003,NC_000004")
    fam = rows[2]
    assert fam["species"] is None and fam["species_count"] == 3 and fam["accession"] is None
    assert fam["genome"] == "ssRNA(+); ssRNA(-); ssRNA(+/-); ssRNA-RT"
    old = rows[3]  # old name placed through NCBI onto the VMR species
    assert (old["name"], old["ncbi_taxid"], old["accession"]) == ("Alphavirus two", "555", "NC_000003;NC_000004")
    assert "Old alpha name" in old["note"]
    far = rows[4]
    assert (far["family"], far["species"], far["ncbi_taxid"]) == ("Farviridae", "Far away virus", "777")
    assert "did you mean" in rows[5]["note"]
    assert rows[6]["query"] == "ALPHAVIRUS ONE"  # input kept as typed
    assert list(rows[0]) == annotate_mod.COLUMNS


def test_annotate_local_and_ncbi_failure(hvmr, ncbi_lookup):
    local = annotate_mod.annotate(hvmr, ["Alphavirus one", "Old alpha name"])
    assert local[0].ncbi_taxid is None
    assert local[1].status == "not_found" and "NCBI not consulted" in local[1].note
    failing = annotate_mod.annotate(hvmr, ["Alphavirus one"], ncbi=FakeNCBI(fail=True))
    assert failing[0].status == "ictv" and "lookup failed" in failing[0].note


def test_annotate_cli_writes_csv(tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("Coronaviridae\nBetacoronavirus pandemicum\nNotavirusviridae\n")
    out = tmp_path / "out.csv"
    result = runner.invoke(app, ["annotate", str(names), "-o", str(out), "--local"])
    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(out.open()))
    assert [r["status"] for r in rows] == ["ictv", "ictv", "not_found"]
    assert rows[1]["genus"] == "Betacoronavirus" and rows[1]["accession"]
    assert "3 names annotated" in result.stderr


def test_annotate_cli_formats(tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("Coronaviridae\n")
    tsv = runner.invoke(app, ["annotate", str(names), "-o", str(tmp_path / "o.tsv"), "--local"])
    assert tsv.exit_code == 0
    assert (tmp_path / "o.tsv").read_text().split("\n")[0].startswith("query\tstatus\t")
    js = runner.invoke(app, ["--json", "annotate", str(names), "--local"])
    assert json.loads(js.stdout)[0]["family"] == "Coronaviridae"
    stdout_csv = runner.invoke(app, ["annotate", str(names), "--local"])
    assert stdout_csv.stdout.startswith("query,status,name,rank,realm")


def test_annotate_cli_bad_input(tmp_path):
    assert runner.invoke(app, ["annotate", str(tmp_path / "missing.txt")]).exit_code == 2
    empty = tmp_path / "empty.txt"
    empty.write_text("# nothing\n\n")
    assert runner.invoke(app, ["annotate", str(empty), "--local"]).exit_code == 2
    nothing = tmp_path / "nothing.txt"
    nothing.write_text("Nopeviridae\n")
    assert runner.invoke(app, ["annotate", str(nothing), "--local"]).exit_code == 1


# -- CSV / TSV output ---------------------------------------------------------------

def test_cell_normalisation():
    assert _cell(None) == "" and _cell(True) == "true" and _cell(["a", "b"]) == "a;b"
    assert _cell("line one\nline two\tend") == "line one | line two end"


def test_markdown_sections():
    md = "intro\n# Title\ntext\n## Sub\nmore\n"
    assert _sections(md) == [("", 0, "intro"), ("Title", 1, "text"), ("Sub", 2, "more")]


def test_global_format_flags_are_exclusive():
    result = runner.invoke(app, ["--csv", "--tsv", "tax", "Coronaviridae"])
    assert result.exit_code == 2


@pytest.mark.parametrize("args, header", [
    (["tax", "Coronaviridae"], ["rank", "name"]),
    (["members", "Coronaviridae", "--rank", "genus"],
     ["name", "rank", "species_count", "parent", "host_filter", "genome_filter"]),
    (["members", "Coronaviridae"], ["rank", "count", "parent", "host_filter", "genome_filter"]),
    (["members", "Coronaviridae", "--tree"], ["depth", "rank", "name", "parent", "path"]),
    (["seq", "--taxon", "Coronaviridae", "--meta"],
     ["name", "rank", "species", "isolates", "accessions", "refseq", "accessions_ssRNA(+)"]),
    (["tree", "Coronaviridae"],
     ["family", "tree", "tip_label", "name", "accession", "species", "genus", "matched"]),
    (["msa", "Coronaviridae", "--range", "1:20"], ["name", "sequence", "matched", "start", "n_cols"]),
    (["diagnose"], ["isolates", "accessions", "empty_accession_rows", "unparsed_rows"]),
    (["cache", "info"], ["namespace", "entries", "bytes"]),
])
@pytest.mark.parametrize("flag, delimiter", [("--csv", ","), ("--tsv", "\t")])
def test_commands_emit_tables(args, header, flag, delimiter):
    result = runner.invoke(app, [flag, *args])
    assert result.exit_code == 0, result.output
    rows = list(csv.reader(io.StringIO(result.stdout), delimiter=delimiter))
    assert rows[0] == header and len(rows) > 1
    assert all(len(r) == len(header) for r in rows)


def test_profile_csv_is_one_wide_row():
    result = runner.invoke(app, ["--csv", "Coronaviridae", "--local"])
    assert result.exit_code == 0, result.output
    header, row = list(csv.reader(io.StringIO(result.stdout)))
    assert len(header) == len(set(header))
    rec = dict(zip(header, row))
    assert rec["family"] == "Coronaviridae" and rec["species"] == ""
    assert int(rec["species_count"]) > 0
    assert rec["viralzone_credit"].startswith("ViralZone")


def test_explicit_raw_format_wins_over_csv():
    result = runner.invoke(app, ["--csv", "tree", "Coronaviridae", "--newick"])
    assert result.exit_code == 0
    assert result.stdout.strip().endswith(";") and "," in result.stdout


def test_table_renderer_writes_file(tmp_path):
    out = tmp_path / "t.tsv"
    TableRenderer("\t")._write(["a", "b"], [{"a": 1, "b": "x\ny"}], str(out))
    assert out.read_text() == "a\tb\n1\tx | y\n"
