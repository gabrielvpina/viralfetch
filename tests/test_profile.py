"""The default command: taxon profile (`viralfetch NAME` / `viralfetch info NAME`)."""

import json

import pytest
from typer.testing import CliRunner

from viralfetch import compare
from viralfetch import config as config_mod
from viralfetch import profile
from viralfetch import viralzone as vz_mod
from viralfetch.cache import Cache
from viralfetch.cli import app
from viralfetch.datasets import DatasetsClient, GenomeRecord, TaxonReport
from viralfetch.ncbi import MetaResult, NcbiLineage, NCBIError, SeqMeta
from viralfetch.queries import TaxonNotFound

runner = CliRunner()

NCBI_ALPHAVIRUS = NcbiLineage(
    taxid="111", name="Alphavirus", rank="genus",
    lineage=[("no rank", "Viruses"), ("realm", "Testviria"), ("kingdom", "Testvirae"),
             ("phylum", "Testphy"), ("class", "Testcla"), ("order", "Otherord"),
             ("family", "Alphaviridae"), ("genus", "Alphavirus")],
)
NCBI_STRAIN = NcbiLineage(
    taxid="999", name="Alpha virus strain X", rank="no rank",
    lineage=NCBI_ALPHAVIRUS.lineage + [("species", "Alphavirus one"), ("no rank", "Alpha virus strain X")],
)


class FakeNCBI:
    """Stands in for NCBIClient; lineage lookups are stubbed via compare."""

    def __init__(self, fail_counts: bool = False):
        self.fail_counts = fail_counts
        self.terms: list[str] = []

    def nuccore_count(self, term):
        self.terms.append(term)
        if self.fail_counts:
            raise NCBIError("count boom")
        return 30 if "complete genome" in term else 120

    def taxid_for_accession(self, accession):
        return None

    def esummary_nuccore(self, accessions):
        self.summarised = list(accessions)
        return MetaResult(records=[
            SeqMeta(a, "org", 1, "rna", "genomic", "circular" if a == "NC_000003.1" else "linear",
                    "complete", "refseq", "2020", geneticcode="1")
            for a in accessions
        ])

    def efetch_taxonomy(self, taxid):
        raise AssertionError("not expected")


class FakeDatasets:
    def __init__(self, fail: bool = False):
        self.fail = fail
        self.asked: list = []

    def taxonomy_report(self, taxid):
        self.asked.append(taxid)
        if self.fail:
            raise NCBIError("boom")
        return TaxonReport(taxid=int(taxid), name="Alphavirus", rank="genus",
                           ictv_status="Name is currently accepted by the ICTV.",
                           genomic_moltype="ssRNA(+)", assemblies=7)

    def refseq_genomes(self, taxid):
        return 1, [GenomeRecord("NC_000001.1", "alphavirus one", 111, 9000, "complete", None,
                                "RefSeq", "2001-01-01", "2020-01-01")]

    def refseq_reports(self, taxid, kind="dataset_report"):
        if kind == "annotation_report":
            return len(GENOME_ANNOTATIONS), GENOME_ANNOTATIONS
        return len(GENOME_RECORDS), GENOME_RECORDS


def _genome(acc, taxid, length, *proteins, completeness="COMPLETE"):
    record = {"accession": acc, "virus": {"tax_id": taxid, "organism_name": f"virus {taxid}"},
              "length": length, "completeness": completeness, "source_database": "RefSeq",
              "host": {"organism_name": "Homo sapiens"}}
    cds, start = [], 1
    for name, aa in proteins:
        end = start + aa * 3 + 2
        cds.append({"name": name, "nucleotide": {"range": [{"begin": str(start), "end": str(end)}]},
                    "protein": {"range": [{"begin": "1", "end": str(aa)}]}})
        start = end + 1
    return record, {"accession": acc, "genes": [{"cds": [c]} for c in cds]}


_GENOMES = [
    _genome("NC_000001.1", 1, 9000, ("nonstructural polyprotein", 2500), ("capsid protein", 260)),
    _genome("NC_000002.1", 2, 9500, ("polyprotein", 2500), ("coat protein", 250)),
    _genome("NC_000003.1", 3, 10000, ("putative capsid protein", 270), ("hypothetical protein", 90)),
    _genome("NC_000004.1", 4, 1200, ("NS5 protein", 345)),                  # a single-gene record
    _genome("NC_000005.1", 5, 9800, ("capsid protein", 250), completeness="PARTIAL"),
]
GENOME_RECORDS = [g[0] for g in _GENOMES]
GENOME_ANNOTATIONS = [g[1] for g in _GENOMES]


def _viralzone() -> vz_mod.ViralZone:
    fam = {"id": 1, "name": "Alphaviridae", "rank": "family", "taxid": 1000, "url": "u1",
           "lineage": {"family": "Alphaviridae"}, "species": [], "host": "Vertebrates"}
    gen = {"id": 2, "name": "Alphavirus", "rank": "genus", "taxid": 111, "url": "u2",
           "lineage": {"family": "Alphavridae", "genus": "Alphavirus"},  # typo, as on the site
           "species": [{"name": "Alpha virus strain X", "highlighted": True}],
           "host": "Human", "diseases": "Fever"}
    other = {"id": 3, "name": "Gammavirus", "rank": "genus", "taxid": 333, "url": "u3",
             "lineage": {"family": "Gammaviridae"}, "species": [], "host": "Plants"}
    vz = vz_mod.ViralZone(generated_at="x", factsheets={1: fam, 2: gen, 3: other}, taxa=[
        {"name": "Alphavirus", "rank": "genus", "id": 2, "family": "Alphaviridae",
         "subfamily": None, "hosts": ["human", "vertebrate"], "page": 9},
        {"name": "Betavirus", "rank": "genus", "id": None, "family": "Alphaviridae",
         "subfamily": None, "hosts": ["invertebrate"], "page": 9},
    ])
    for fs in vz.factsheets.values():
        vz.by_taxid[fs["taxid"]] = fs
        vz.by_name[fs["name"].casefold()] = fs
        for sp in fs["species"]:
            vz.by_species.setdefault(sp["name"].casefold(), []).append(fs)
    for t in vz.taxa:
        vz.taxa_by_name[t["name"].casefold()] = t
    return vz


@pytest.fixture
def ncbi_knows(monkeypatch):
    table = {"alphavirus": NCBI_ALPHAVIRUS, "alpha virus strain x": NCBI_STRAIN}
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: table.get(name.casefold()))


# -- service ----------------------------------------------------------------------

def test_profile_of_a_vmr_genus(vmr, ncbi_knows):
    datasets, ncbi = FakeDatasets(), FakeNCBI()
    p = profile.build(vmr, "alphavirus", ncbi=ncbi, datasets=datasets, viralzone=_viralzone())
    assert (p.name, p.rank, p.family) == ("Alphavirus", "genus", "Alphaviridae")
    assert p.notes == [] and p.warnings == []
    rows = {r.rank: r for r in p.comparison}
    assert rows["family"].agree is True
    assert (rows["order"].ictv, rows["order"].ncbi, rows["order"].agree) == ("Testord", "Otherord", False)
    assert "no rank" not in rows
    assert p.vmr.species == 2 and p.vmr.isolates == 3 and p.vmr.exemplars == 2
    assert p.vmr.genome_compositions == {"ssRNA(+)": 3}
    assert datasets.asked == ["111"]
    assert ncbi.terms == ["txid111[Organism:exp]",
                          'txid111[Organism:exp] AND "complete genome"[Title]']
    assert (p.sequences.total, p.sequences.complete, p.sequences.refseq) == (120, 30, 5)
    assert p.sequences.refseq_genomes[0].accession == "NC_000001.1"
    assert p.sequences.portal_url.endswith("taxid:111")
    assert p.ncbi_report.assemblies == 7
    # Genome summary: one partial and one single-gene record set aside.
    g = p.genomes
    assert (g.genomes, g.partial, g.short) == (3, 1, 1)
    assert (g.size_min, g.size_median, g.size_max) == (9000, 9500, 10000)
    assert g.proteins[0].name == "capsid protein" and g.proteins[0].share == 1.0
    assert g.proteins[0].median_aa == 260
    assert g.polyprotein_genomes == 2 and g.hosts == {"Homo sapiens": 3}
    # The genus factsheet wins (despite its own lineage's typo); the family is related.
    assert p.viralzone.factsheet["name"] == "Alphavirus"
    assert [r["name"] for r in p.viralzone.related] == ["Alphaviridae"]
    assert p.viralzone.note is None
    assert p.viralzone.host_groups == ["human", "vertebrate"]
    assert p.viralzone.credit["source"].startswith("ViralZone")


def test_profile_of_a_name_only_ncbi_knows(vmr, ncbi_knows):
    p = profile.build(vmr, "Alpha virus strain X", ncbi=FakeNCBI(), datasets=FakeDatasets(),
                      viralzone=_viralzone())
    assert (p.name, p.rank) == ("Alpha virus strain X", "no rank")
    assert p.ictv.name == "Alphavirus one"  # nearest VMR taxon gives the ICTV side
    assert "NCBI places it in species Alphavirus one" in p.notes[0]
    # Listed as a species on the Alphavirus factsheet.
    assert p.viralzone.factsheet["name"] == "Alphavirus"


def test_profile_local_only(vmr):
    p = profile.build(vmr, "Alphavirus one", viralzone=_viralzone())
    assert p.ncbi is None and p.sequences is None
    assert all(r.ncbi is None for r in p.comparison)
    assert p.vmr.exemplar_accessions == ["NC_000001"]
    assert p.viralzone.factsheet["name"] == "Alphavirus"
    assert "no Alphavirus one factsheet" in p.viralzone.note


def test_profile_ignores_unrelated_factsheets(vmr):
    p = profile.build(vmr, "Gammavirus", viralzone=_viralzone())
    assert p.viralzone.factsheet["name"] == "Gammavirus"
    p = profile.build(vmr, "Betavirus", viralzone=_viralzone())
    assert p.viralzone.factsheet["name"] == "Alphaviridae"
    assert p.viralzone.related == []


def test_profile_without_viralzone_match(vmr):
    vz = vz_mod.ViralZone(generated_at="x", factsheets={}, taxa=[])
    assert profile.build(vmr, "Gammaviridae", viralzone=vz).viralzone is None


def test_profile_datasets_failure_keeps_counts(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(fail=True),
                      viralzone=_viralzone())
    assert (p.sequences.total, p.sequences.refseq, p.sequences.refseq_genomes) == (120, None, [])
    assert p.ncbi_report is None
    assert p.warnings == ["NCBI Datasets request failed: boom"]


def test_profile_all_ncbi_stats_failing(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(fail_counts=True),
                      datasets=FakeDatasets(fail=True), viralzone=_viralzone())
    assert p.sequences is None
    assert sorted(p.warnings) == ["NCBI Datasets request failed: boom",
                                  "NCBI record count failed: count boom"]


def test_profile_unknown_name(vmr, monkeypatch):
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: None)
    with pytest.raises(TaxonNotFound):
        profile.build(vmr, "Nopeviridae", ncbi=FakeNCBI(), viralzone=_viralzone())


def test_higher_rank_gets_no_accession_fallback(vmr, monkeypatch):
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: None)
    p = profile.build(vmr, "Alphaviridae", ncbi=FakeNCBI(), datasets=FakeDatasets())
    assert p.ncbi is None and p.sequences is None
    assert p.notes == ["NCBI taxonomy has no family named Alphaviridae."]


# -- Datasets client ------------------------------------------------------------------

class FakeResponse:
    def __init__(self, status: int, payload: dict):
        self.status_code = status
        self._payload = payload
        self.text = json.dumps(payload)

    def json(self):
        return self._payload


class FakeSession:
    def __init__(self, routes):
        self.routes = routes  # (substring, response) pairs, first match wins
        self.urls: list[str] = []

    def get(self, url, headers=None, timeout=None):
        self.urls.append(url)
        for needle, response in self.routes:
            if needle in url:
                return response
        return FakeResponse(404, {})


def _client(session, tmp_path=None) -> DatasetsClient:
    cfg = config_mod.Config(email="a@b.c", api_key=None)
    cache = Cache(tmp_path) if tmp_path else None
    return DatasetsClient(cfg, cache=cache, session=session, sleep=lambda s: None)


REFSEQ = {"total_count": 2, "reports": [
    {"accession": "NC_2.1", "virus": {"tax_id": 5, "organism_name": "b virus"}, "length": 100,
     "completeness": "COMPLETE", "segment": "S", "source_database": "RefSeq",
     "release_date": "2001-01-01T00:00:00Z", "update_date": "2020-02-02T00:00:00Z"},
    {"accession": "NC_1.1", "virus": {"tax_id": 4, "organism_name": "A virus"}, "length": 200,
     "completeness": "COMPLETE"},
]}


def test_refseq_genomes_sorted_by_organism():
    session = FakeSession([("refseq_only", FakeResponse(200, REFSEQ))])
    total, genomes = _client(session).refseq_genomes(11050)
    assert total == 2
    assert [g.accession for g in genomes] == ["NC_1.1", "NC_2.1"]
    assert genomes[1].segment == "S"
    assert genomes[1].update_date == "2020-02-02"


def test_datasets_retries_then_raises():
    session = FakeSession([("dataset_report", FakeResponse(500, {}))])
    with pytest.raises(NCBIError):
        _client(session).refseq_genomes(1)
    assert len(session.urls) == 2  # one retry only on virus reports


def test_taxonomy_report_parsing_and_cache(tmp_path):
    payload = {"reports": [{"taxonomy": {
        "tax_id": 11050, "rank": "FAMILY",
        "current_scientific_name": {"name": "Flaviviridae", "notes": [
            {"name": "ICTV Status", "note": "Name is currently accepted by the ICTV."}]},
        "counts": [{"type": "COUNT_TYPE_GENE", "count": 3}, {"type": "COUNT_TYPE_ASSEMBLY", "count": 994}],
        "genomic_moltype": "ssRNA(+)", "children": [1, 2],
    }}]}
    session = FakeSession([("taxonomy/taxon/11050", FakeResponse(200, payload))])
    client = _client(session, tmp_path)
    report = client.taxonomy_report(11050)
    assert report == TaxonReport(taxid=11050, name="Flaviviridae", rank="family",
                                 ictv_status="Name is currently accepted by the ICTV.",
                                 genomic_moltype="ssRNA(+)", assemblies=994, children=2)
    client.taxonomy_report(11050)
    assert len(session.urls) == 1  # second call served from the cache


# -- CLI -------------------------------------------------------------------------------

def test_name_alone_runs_the_profile():
    result = runner.invoke(app, ["--json", "Coronaviridae", "--local"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.stdout)
    assert (payload["name"], payload["rank"]) == ("Coronaviridae", "family")
    assert payload["ncbi"] is None
    assert payload["warnings"] == ["NCBI was skipped (--local)."]
    assert payload["viralzone"]["factsheet"]["name"] == "Coronaviridae"
    assert payload["viralzone"]["credit"]["citation"].endswith("gkad946")


def test_unquoted_multiword_name():
    result = runner.invoke(app, ["--json", "Betacoronavirus", "pandemicum", "--local"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["name"] == "Betacoronavirus pandemicum"


def test_info_rich_output_credits_viralzone():
    result = runner.invoke(app, ["info", "Coronaviridae", "--local"])
    assert result.exit_code == 0, result.output
    assert "ViralZone · Coronaviridae" in result.stdout
    assert "ViralZone (SIB) · CC BY 4.0 · doi:10.1093/nar/gkad946" in result.stdout
    assert "https://viralzone.expasy.org/30" in result.stdout  # piped: the factsheet URL is spelt out


def test_existing_commands_still_route():
    result = runner.invoke(app, ["--json", "tax", "Coronaviridae"])
    assert result.exit_code == 0
    assert json.loads(result.stdout)["name"] == "Coronaviridae"


def test_unknown_name_suggests_command(monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "test@example.com")
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: None)
    result = runner.invoke(app, ["membrs"])
    assert result.exit_code == 1
    assert "viralfetch members" in result.stderr


def test_no_email_profile_stays_local(tmp_path, monkeypatch):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    result = runner.invoke(app, ["--json", "Coronaviridae"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.stdout)
    assert payload["ncbi"] is None
    assert "no NCBI email" in payload["warnings"][0]


# -- rich rendering ------------------------------------------------------------------

def _render(p, full=False, width=100) -> str:
    import io
    from rich.console import Console
    from viralfetch.render import rich_
    buf = io.StringIO()
    saved = rich_._out
    rich_._out = Console(file=buf, width=width, color_system=None)
    try:
        rich_.profile(p, full=full)
    finally:
        rich_._out = saved
    return buf.getvalue()


def test_render_compact_profile(vmr, ncbi_knows, monkeypatch):
    monkeypatch.setattr(compare, "lineage_via_ncbi",
                        lambda client, name: NCBI_ALPHAVIRUS if name.casefold() == "alphavirus" else None)
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    p.comparison = [r for r in p.comparison if r.agree is not False]  # make the lineages agree
    out = _render(p)
    assert "Testviria › Testvirae" in out and "Alphavirus" in out
    assert "✓ ICTV and NCBI agree" in out
    assert "lineage — ICTV vs NCBI" not in out          # table only on disagreement / --full
    assert "RefSeq sequence" in out and "NC_000001.1" not in out  # RefSeq list only with --full
    assert "● human" in out and "ViralZone · Alphavirus" in out
    assert "--full:" in out


def test_render_disagreement_shows_the_table(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    out = _render(p)
    assert "≠ ICTV and NCBI differ at order" in out
    assert "lineage — ICTV vs NCBI" in out and "Otherord" in out


def test_render_full_adds_details(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    out = _render(p, full=True)
    assert "NC_000001.1" in out
    assert "also: Alphaviridae (family)" in out
    assert "--full:" not in out


def test_render_narrow_terminal_stacks_cards(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    out = _render(p, width=70)
    assert "╭─ ICTV (VMR)" in out and "╭─ NCBI" in out
    assert max(len(line) for line in out.splitlines()) <= 70


def test_field_text_truncates_and_bolds_qualifiers():
    from viralfetch.render import rich_
    long = "DENV: " + " ".join(["word"] * 80)
    compact = rich_._field_text(long, full=False)
    assert compact.plain.endswith(" …") and len(compact.plain) < 260
    assert compact.spans and compact.plain.startswith("DENV:")
    assert rich_._field_text(long, full=True).plain == long
