"""Rebuilding the ICTV trees data set (`viralfetch update --trees`)."""

import json

import pytest
from typer.testing import CliRunner

from viralfetch import config as config_mod
from viralfetch import ictv
from viralfetch import trees as trees_mod
from viralfetch import trees_install as ti
from viralfetch import vmr as vmr_mod
from viralfetch.cli import app
from viralfetch.ictv import ChapterNotFound, ICTVClient, ICTVError
from viralfetch.models import Chapter, Resource

runner = CliRunner()

TREE = "((X15656_Tomato_yellow_leaf_curl_virus:0.1,FM877473_African_cassava_mosaic_virus:0.2)90:0.1,Unknown_isolate_1:0.3);\n"
ALIGNMENT = (
    ">X15656 Tomato yellow leaf curl virus\nACGT-ACGT\n"
    ">FM877473 African cassava mosaic virus\nACGTTACGT\n"
    ">Unknown isolate 1\nACGA-ACGT\n"
)


# -- parsing ----------------------------------------------------------------

def test_newick_tips_in_tree_order():
    assert ti.parse_tree_tips(TREE) == [
        "X15656_Tomato_yellow_leaf_curl_virus",
        "FM877473_African_cassava_mosaic_virus",
        "Unknown_isolate_1",
    ]


def test_newick_quoted_labels_and_comments():
    text = "(('A virus':1[&c=1],B_virus:2)[&x]0.9:1,C);"
    assert ti.parse_tree_tips(text) == ["A_virus", "B_virus", "C"]


def test_newick_picks_the_tree_with_most_tips():
    text = "(A,B);\n(A,B,C,D);\n(Z);\n"
    assert ti.parse_tree_tips(text) == ["A", "B", "C", "D"]


def test_newick_malformed_yields_nothing():
    assert ti.parse_tree_tips("((A,B);") == []
    assert ti.parse_tree_tips("") == []


def test_nexus_applies_translate_table():
    text = """#NEXUS
begin trees;
    translate
        1 MN908947,
        2 'Bat coronavirus',
        3 NC_004718
    ;
    tree con_50 = [&U] (1:0.1[&prob=1],(2:0.2,3:0.3)1.00:0.1);
end;
"""
    assert ti.parse_tree_tips(text) == ["MN908947", "Bat_coronavirus", "NC_004718"]


def test_alignment_fasta_skips_leading_comments():
    text = "; exported by some tool\n>a x\nAC GT\nAA\n>b\nTTTT\n"
    assert ti.parse_alignment(text) == [("a x", "ACGTAA"), ("b", "TTTT")]


def test_alignment_mega():
    text = "#mega\n!Title demo;\n#seq1\nACGT\nAC\n#seq2\nTTTT\n"
    assert ti.parse_alignment(text) == [("seq1", "ACGTAC"), ("seq2", "TTTT")]


def test_alignment_binary_or_empty():
    assert ti.parse_alignment("") == []
    assert ti.parse_alignment("\x00\x01binary") == []


def test_molecule_detection():
    assert ti.detect_molecule(["ACGT-ACGTN"]) == "nt"
    assert ti.detect_molecule(["MKLVEFPQ"]) == "AA"
    assert ti.detect_molecule(["----"]) == "unknown"


@pytest.mark.parametrize("label, accession, name", [
    ("X15656_Tomato_yellow_leaf_curl_virus", "X15656", "Tomato yellow leaf curl virus"),
    ("TYLCV_NC_004005", "NC_004005", "TYLCV"),
    ("ABC_KX123456_688_bp", "KX123456", "ABC 688 bp"),
    ("Some_virus_ACCESSION_NOT_ON_SPREADSHEET", None, "Some virus"),
])
def test_member_labels(label, accession, name):
    member = ti.parse_member_label(label)
    assert (member.accession, member.name) == (accession, name)


def test_method_and_region():
    assert ti.parse_method("Maximum likelihood tree (IQ-TREE)") == "maximum likelihood"
    assert ti.parse_region("Fam_RdRp_aa.fasta") == ("RdRp", "protein")
    assert ti.parse_region("complete genome sequences") == ("whole genome", "whole_genome")
    assert ti.parse_region("nothing here") == (None, "unknown")


# -- ICTV pages -------------------------------------------------------------

def test_parse_chapter_index():
    html = ('<a href="/report/chapter/geminiviridae">x</a>'
            '<a href="/report/chapter/Coronaviridae/">y</a><a href="/other">z</a>')
    assert ictv.parse_chapter_index(html) == ["coronaviridae", "geminiviridae"]


def test_parse_resources_groups_links_by_figure():
    html = """<div class="field--name-field-mt-srv-body">
      <h2>Sequence alignments and tree files</h2>
      <h3>Figure 3. Geminiviridae:</h3>
      <a href=" /sites/fig3.fasta ">Alignment</a> <a href="/sites/fig3.nwk">Tree</a>
      <p><strong>Figure 4A</strong></p>
      <a href="https://ictv.global/sites/fig4a_RdRp.fas">download</a>
      <h2>Other section</h2><a href="/sites/ignored.nwk">Tree</a>
    </div>"""
    got = ictv.parse_resources(html, url="https://ictv.global/x/resources")
    assert got == [
        Resource("Figure_3", "https://ictv.global/sites/fig3.fasta", "https://ictv.global/sites/fig3.nwk"),
        Resource("Figure_4A", "https://ictv.global/sites/fig4a_RdRp.fas", None),
    ]


def test_parse_resources_without_section():
    html = '<div class="field--name-field-mt-srv-body"><h2>Nothing</h2></div>'
    assert ictv.parse_resources(html, url="u") == []


# -- install ----------------------------------------------------------------

BASE = "https://ictv.global/sites/"


class FakeClient(ICTVClient):
    """Serves canned ICTV pages; ``fail`` names families whose pages error."""

    def __init__(self, resources: dict[str, list[Resource]], fail=()):
        self.cache = None
        self.resources = resources
        self.fail = set(fail)
        self.files = {BASE + "gem.nwk": TREE, BASE + "gem.fasta": ALIGNMENT}

    def list_chapters(self, *, fresh=False):
        return ["geminiviridae", "coronaviridae", "riboviria"]

    def fetch_resources(self, name, *, fresh=False):
        if name in self.fail:
            raise ICTVError("HTTP 500")
        if name not in self.resources:
            raise ChapterNotFound(name, "url")
        return self.resources[name]

    def fetch_chapter(self, name, *, fresh=False):
        return Chapter(slug=name.lower(), title=name, markdown=f"# {name}\n")

    def fetch_figure_captions(self, name, *, fresh=False):
        return {"1": "Figure 1. Maximum likelihood phylogeny of complete genome sequences."}

    def download_file(self, url):
        return self.files[url]


GEM = {"Geminiviridae": [Resource("Figure_1", BASE + "gem.fasta", BASE + "gem.nwk")]}


@pytest.fixture
def data_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(config_mod, "DATA_DIR", tmp_path / "data")
    return tmp_path / "data"


def test_install_full_rebuild(data_dir):
    seen = []
    result = ti.install(FakeClient(GEM), progress=lambda i, n, e: seen.append((i, n, e["status"])))

    assert result.path == data_dir / "ictv-trees"
    assert (result.included, result.omitted, result.errors) == (1, 1, [])
    assert seen == [(1, 2, "omitted"), (2, 2, "included")]  # VMR order; realms skipped
    assert trees_mod.root() == result.path

    tree_dir = result.path / "geminiviridae" / "trees" / "tree1"
    assert (tree_dir / "tree.nwk").read_text() == TREE          # verbatim
    assert (tree_dir / "alignment.fasta").read_text() == ALIGNMENT
    meta = json.loads((tree_dir / "tree.json").read_text())
    assert meta["molecule"] == "nt"
    assert meta["method"] == "maximum likelihood"
    assert (meta["region"], meta["region_source"]) == ("whole genome", "caption")
    assert meta["n_tips"] == 3 and meta["alignment_coverage"] == 1.0

    rows = {r["accession"]: r for r in trees_mod._load_members(tree_dir / "members.tsv").values()}
    assert rows["X15656"]["genus"] == "Begomovirus"             # joined to the VMR
    assert rows["X15656"]["vmr_match"] == "True"
    assert rows[""]["vmr_match"] == "False"                     # unknown isolate kept

    index = json.loads((result.path / "_index.json").read_text())
    assert index["counts"] == {"total": 2, "included": 1, "omitted": 1, "error": 0}
    assert not list(data_dir.glob("*.partial"))

    # The rebuilt data set drives `tree` resolution.
    resolved = trees_mod.resolve(vmr_mod.load(), "Begomovirus")
    assert resolved.slug == "geminiviridae" and resolved.trees[0].matched


def test_install_partial_merges_into_current(data_dir):
    result = ti.install(FakeClient(GEM), families=["geminiviridae"])
    assert result.partial and result.families == ["Geminiviridae"]
    # Every other bundled family is carried over untouched.
    assert (result.path / "coronaviridae" / "trees").is_dir()
    index = json.loads((result.path / "_index.json").read_text())
    assert index["counts"] == trees_mod.index()["counts"]
    assert len(index["families"]) == len(json.loads(
        (trees_mod.bundled_root() / "_index.json").read_text())["families"])


def test_failed_family_keeps_its_previous_trees(data_dir):
    before = sorted(p.name for p in (trees_mod.bundled_root() / "coronaviridae" / "trees").iterdir())
    result = ti.install(FakeClient(GEM, fail={"Coronaviridae"}), families=["Coronaviridae"])
    assert [e["family"] for e in result.errors] == ["Coronaviridae"]
    assert sorted(p.name for p in (result.path / "coronaviridae" / "trees").iterdir()) == before
    entry = next(e for e in trees_mod.index()["families"] if e["slug"] == "coronaviridae")
    assert entry["status"] == "included"


def test_nothing_built_leaves_current_data_set(data_dir):
    with pytest.raises(ti.TreesInstallError):
        ti.install(FakeClient({}))
    assert trees_mod.root() == trees_mod.bundled_root()
    assert not (data_dir / "ictv-trees").exists()


def test_unknown_family_is_rejected(data_dir):
    with pytest.raises(ti.TreesInstallError, match="Notaviridae"):
        ti.install(FakeClient(GEM), families=["Notaviridae"])


def test_stale_install_loses_to_newer_bundle(data_dir):
    result = ti.install(FakeClient(GEM))
    index_path = result.path / "_index.json"
    index = json.loads(index_path.read_text())
    index["generated_at"] = "2000-01-01T00:00:00+00:00"
    index_path.write_text(json.dumps(index))
    assert trees_mod.root() == trees_mod.bundled_root()


def test_reset(data_dir):
    assert ti.reset() is False
    ti.install(FakeClient(GEM))
    assert ti.reset() is True
    assert trees_mod.root() == trees_mod.bundled_root()


# -- CLI --------------------------------------------------------------------

def test_cli_update_trees(monkeypatch, data_dir):
    monkeypatch.setattr("viralfetch.cli._make_ictv_client", lambda cfg, out: FakeClient(GEM))
    result = runner.invoke(app, ["--json", "update", "--trees"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.stdout)["installed"]
    assert (payload["included"], payload["omitted"]) == (1, 1)

    result = runner.invoke(app, ["--json", "update", "--trees", "--reset"])
    assert json.loads(result.stdout) == {"removed": True, "active": "bundled"}


def test_cli_update_trees_failure_exits_4(monkeypatch, data_dir):
    monkeypatch.setattr("viralfetch.cli._make_ictv_client", lambda cfg, out: FakeClient({}))
    result = runner.invoke(app, ["update", "--trees"])
    assert result.exit_code == 4


@pytest.mark.parametrize("args", [
    ["update", "--family", "Geminiviridae"],
    ["update", "--trees", "--reset", "--family", "Geminiviridae"],
])
def test_cli_update_trees_bad_flags(args):
    assert runner.invoke(app, args).exit_code == 2
