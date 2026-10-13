"""The ViralZone data set: parsing factsheets / realm pages, building, reading."""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from viralfetch import config as config_mod
from viralfetch import viralzone as vz_mod
from viralfetch import viralzone_install as vi
from viralfetch.cli import app

FIXTURES = Path(__file__).parent / "fixtures" / "viralzone"
runner = CliRunner()


def _fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


# -- links -------------------------------------------------------------------

@pytest.mark.parametrize("href, expected", [
    ("/23", 23),
    ("/preview%5Fby%5Fspecies/207", 207),
    ("/preview_by_protein/10557", 10557),
    ("/by%5Fspecies/664", 664),
    ("https://viralzone.expasy.org/176", 176),
    ("https://www.uniprot.org/taxonomy/11266", None),
    ("/resources/Filovirus_virion.svg", None),
    ("#tab2", None),
    (None, None),
])
def test_page_id(href, expected):
    assert vi.page_id(href) == expected


# -- realm pages ---------------------------------------------------------------

def test_taxonomy_index_lists_realm_pages():
    pages = vi.parse_taxonomy_index(_fixture("home.html"))
    ids = [p["id"] for p in pages]
    assert 237 in ids and 11365 in ids
    assert {"id": 237, "label": "Negarnaviricota"} in pages
    assert len(ids) == len(set(ids))


def test_realm_page_hierarchy_and_hosts():
    taxa = vi.parse_realm_page(_fixture("realm_237.html"), 237)
    by_name = {t["name"]: t for t in taxa}
    assert by_name["Filoviridae"] == {
        "name": "Filoviridae", "rank": "family", "id": 23, "family": None,
        "subfamily": None, "hosts": [], "page": 237,
    }
    ebola = by_name["Orthoebolavirus"]
    assert (ebola["rank"], ebola["id"], ebola["family"], ebola["hosts"]) == (
        "genus", 207, "Filoviridae", ["human"])
    # Folded genus without a factsheet, multi-host tag, "!" marker stripped.
    assert by_name["Acridvirus"]["id"] is None
    assert by_name["Acridvirus"]["family"] == "Lispiviridae"
    assert by_name["Dianlovirus"]["hosts"] == ["unknown", "vertebrate"]
    assert by_name["Thamnovirus"]["family"] == "Filoviridae"
    assert len({(t["rank"], t["name"], t["family"]) for t in taxa}) == len(taxa)


def test_factsheet_candidates_skip_buttons_and_navigation():
    ids = vi.factsheet_candidates(_fixture("realm_237.html"))
    assert 23 in ids and 207 in ids
    assert 10179 not in ids  # "Mononega-type capping" button
    assert 678 not in ids    # "Human virus" in the site menu


# -- factsheets ----------------------------------------------------------------

def test_family_factsheet():
    fs = vi.parse_factsheet(_fixture("23.html"), 23)
    assert (fs["name"], fs["taxid"], fs["rank"]) == ("Filoviridae", 11266, "family")
    assert fs["url"] == "https://viralzone.expasy.org/23"
    assert fs["baltimore"] == "Group V: Negative sense ssRNA viruses"
    assert fs["lineage"]["order"] == "Mononegavirales"
    assert fs["lineage"]["family"] == "Filoviridae"
    assert "genus" not in fs["lineage"]
    assert "Orthoebolavirus" in fs["children"]["genus"]
    assert {"name": "Zaire ebolavirus", "highlighted": True} in fs["species"]
    assert fs["summary"]["genome"].startswith("Negative-stranded RNA linear genome")
    assert fs["summary"]["virion"].startswith("Filamentous")
    assert fs["replication_site"] == "CYTOPLASMIC"
    assert not fs["summary"]["replication"].startswith("CYTOPLASMIC")
    assert fs["host"] == "Reservoir: Bats\nOccasional: Human and primates"
    assert fs["diseases"].startswith("Hemorragic fever")
    assert fs["transmission"] == "Zoonosis, Contact with body fluids"
    assert fs["vaccine"] == "None"
    assert fs["antiviral_drugs"] is None
    assert [i["title"] for i in fs["host_virus_interactions"]][0] == (
        "Inhibition of host Interferon induction")
    assert fs["reference_strains"][0]["name"] == "Zaire ebolavirus"
    assert fs["reference_strains"][0]["sequences"] == [{"accession": "NC_002549", "segment": None}]
    assert fs["reference_strains"][0]["assemblies"] == ["GCF_000848505"]
    assert fs["links"]["ictv_taxnode_id"] == "202201560"
    assert fs["links"]["ictv_report"].endswith("/filoviridae/filoviridae")
    assert fs["images"]["genome"].endswith("Filovirus%5Fgenome.svg")


def test_genus_factsheet_with_variant_box_titles():
    fs = vi.parse_factsheet(_fixture("178.html"), 178)
    assert (fs["name"], fs["rank"]) == ("Simplexvirus", "genus")
    assert fs["lineage"]["subfamily"] == "Alphaherpesvirinae"
    assert fs["host"] == "Human, mammals"        # NATURAL HOSTS
    assert fs["tropism"].startswith("Primary site of infection")  # TROPISM
    assert fs["etymology"].startswith("Simplex virus")
    assert "Acyclovir" in fs["antiviral_drugs"]
    assert len(fs["species"]) > 1
    assert [s["name"] for s in fs["reference_strains"]] == [
        "Human herpesvirus 1", "Human herpesvirus 2"]


def test_segmented_reference_strain():
    fs = vi.parse_factsheet(_fixture("130.html"), 130)
    strain = fs["reference_strains"][0]
    assert strain["name"].startswith("Influenza A virus A/Goose/Guangdong/1/96")
    assert [s["segment"] for s in strain["sequences"]][:2] == ["RNA1", "RNA2"]
    assert len(strain["sequences"]) == 8
    assert fs["rank"] == "virus"


def test_lineage_row_listing_siblings_keeps_own_name():
    fs = vi.parse_factsheet(_fixture("7777.html"), 7777)
    assert fs["lineage"]["genus"] == "Betasatellite"
    assert fs["reference_strains"][0]["sequences"][0]["accession"] == "NC_003403"


def _record(**overrides) -> dict:
    record = {
        "name": "Testviridae", "lineage": {}, "children": {}, "summary": {},
        "replication_site": None, "other_fields": [], "etymology": None,
    }
    record.update(overrides)
    return record


@pytest.mark.parametrize("first_line, site", [
    ("CYTOPLASMIC", "CYTOPLASMIC"),
    ("NUCLEUS", "NUCLEAR"),
    ("CYTOPLASMIC/NUCLEAR", "CYTOPLASMIC AND NUCLEAR"),
    ("NUCLEO-CYTOPLASMIC", "CYTOPLASMIC AND NUCLEAR"),
])
def test_normalize_replication_site(first_line, site):
    rec = vi.normalize(_record(summary={"replication": f"{first_line}\n- Attachment."}))
    assert rec["replication_site"] == site
    assert rec["summary"]["replication"] == "- Attachment."


def test_normalize_leaves_unstated_site_alone():
    rec = vi.normalize(_record(summary={"replication": "- Virus penetrates into the cell."}))
    assert rec["replication_site"] is None
    assert rec["summary"]["replication"] == "- Virus penetrates into the cell."


def test_normalize_rank_labels_and_boxes():
    rec = vi.normalize(_record(
        name="Phasivirus",
        lineage={"family": "Phenuiviridae", "genera": "Phasivirus", "sub-family": "X"},
        other_fields=[{"title": "ORIGIN OF NAME", "text": "From ..."},
                      {"title": "GENERAL", "text": "Back to adenoviridae page"}],
    ))
    assert rec["lineage"] == {"family": "Phenuiviridae", "genus": "Phasivirus", "subfamily": "X"}
    assert rec["rank"] == "genus"
    assert rec["etymology"] == "From ..."
    assert rec["other_fields"] == []


def test_normalize_strips_textile_links():
    rec = vi.normalize(_record(
        antiviral_drugs='Ribavirin\n[YFV]"Favipiravir":http://en.wikipedia.org/wiki/Favipiravir (Experimental)',
        summary={"genome": '"Myndus":https://x.org/a vector'},
    ))
    assert rec["antiviral_drugs"] == "Ribavirin\n[YFV]Favipiravir (Experimental)"
    assert rec["summary"]["genome"] == "Myndus vector"


def test_page_without_taxid_is_not_a_factsheet():
    assert vi.parse_factsheet(_fixture("realm_237.html"), 237) is None


# -- build / install -----------------------------------------------------------

class FakeClient:
    def __init__(self, pages: dict[str, str]):
        self.pages = pages
        self.requested: list[str] = []

    def get(self, path: str):
        self.requested.append(path)
        return self.pages.get(path)


def _site() -> dict[str, str]:
    pages = {"/": _fixture("home.html"), "/237": _fixture("realm_237.html")}
    for pid in (23, 178, 130, 7777):
        pages[f"/{pid}"] = _fixture(f"{pid}.html")
    return pages


def test_build_refuses_a_broken_site(monkeypatch):
    monkeypatch.setattr(vi, "MIN_FACTSHEETS", 1000)
    with pytest.raises(vi.ViralZoneInstallError):
        vi.build(FakeClient(_site()))


def test_install_and_load(tmp_path, monkeypatch):
    monkeypatch.setattr(vi, "MIN_FACTSHEETS", 1)
    monkeypatch.setattr(config_mod, "DATA_DIR", tmp_path)
    progress = []
    result = vi.install(FakeClient(_site()), progress=lambda i, n, e: progress.append(e))
    assert result.path == tmp_path / vz_mod.DATA_FILENAME
    # Realm 237 links Filoviridae and H5N1; Simplexvirus is not reachable from it.
    assert result.factsheets == 2
    assert any(e["status"] == "factsheet" for e in progress)

    data = json.loads(result.path.read_text(encoding="utf-8"))
    assert data["source"] == vz_mod.CREDIT["source"]
    assert "gkad946" in data["citation"]
    assert data["factsheets"]["23"]["name"] == "Filoviridae"

    vz_mod.load.cache_clear()
    loaded = vz_mod.load()
    assert loaded is not None and loaded.find("filoviridae")["taxid"] == 11266
    assert vi.reset() is True
    assert vi.reset() is False


# -- reader --------------------------------------------------------------------

@pytest.fixture
def dataset(tmp_path) -> vz_mod.ViralZone:
    factsheets = {
        str(pid): vi.parse_factsheet(_fixture(f"{pid}.html"), pid)
        for pid in (23, 178, 130)
    }
    path = tmp_path / "vz.json"
    path.write_text(json.dumps({
        "generated_at": "2026-10-12T00:00:00+00:00",
        "taxa": vi.parse_realm_page(_fixture("realm_237.html"), 237),
        "factsheets": factsheets,
    }), encoding="utf-8")
    return vz_mod.load_from(path)


def test_reader_lookups(dataset):
    assert dataset.find("FILOVIRIDAE")["id"] == 23
    assert dataset.for_taxid("10294")["name"] == "Simplexvirus"
    assert dataset.for_taxid("nope") is None
    assert dataset.taxon("Acridvirus")["family"] == "Lispiviridae"


def test_match_orders_most_specific_first(dataset):
    hits = dataset.match(
        taxids=[11266],
        lineage={"family": "Filoviridae", "genus": "Simplexvirus"},
    )
    assert [h["name"] for h in hits] == ["Simplexvirus", "Filoviridae"]


def test_match_by_listed_species(dataset):
    hits = dataset.match(names=["Zaire ebolavirus"])
    assert [h["name"] for h in hits] == ["Filoviridae"]


def test_active_path_prefers_newer_install(tmp_path, monkeypatch):
    monkeypatch.setattr(config_mod, "DATA_DIR", tmp_path)
    bundled = tmp_path / "bundled.json"
    bundled.write_text(json.dumps({"generated_at": "2026-01-01"}), encoding="utf-8")
    monkeypatch.setattr(vz_mod, "bundled_path", lambda: bundled)
    installed = tmp_path / vz_mod.DATA_FILENAME
    installed.write_text(json.dumps({"generated_at": "9999-01-01"}), encoding="utf-8")
    assert vz_mod.active_path() == installed
    installed.write_text(json.dumps({"generated_at": "1999-01-01"}), encoding="utf-8")
    assert vz_mod.active_path() == bundled


# -- CLI -----------------------------------------------------------------------

def test_update_viralzone_flag_conflicts():
    assert runner.invoke(app, ["update", "--viralzone", "--trees"]).exit_code == 2
    assert runner.invoke(app, ["update", "--viralzone", "--yes"]).exit_code == 2


def test_update_viralzone_reset(tmp_path, monkeypatch):
    monkeypatch.setattr(config_mod, "DATA_DIR", tmp_path)
    result = runner.invoke(app, ["--json", "update", "--viralzone", "--reset"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.output) == {"removed": False, "active": "bundled"}
