"""Genome summary of a taxon's RefSeq genomes (the profile's "Genomes" section)."""

import pytest

from viralfetch import config as config_mod
from viralfetch import genomics
from viralfetch.datasets import DatasetsClient
from viralfetch.render.table_ import profile_row

from .test_profile import (
    FakeDatasets, FakeNCBI, FakeResponse, _render, _viralzone, ncbi_knows, profile,  # noqa: F401
)


def _record(acc, taxid, length, segment=None, completeness="COMPLETE", host=None):
    r = {"accession": acc, "virus": {"tax_id": taxid, "organism_name": f"virus {taxid}"},
         "length": length, "completeness": completeness}
    if segment:
        r["segment"] = segment
    if host:
        r["host"] = {"organism_name": host}
    return r


def _cds(name, begin, end, peptides=()):
    cds = {"name": name,
           "nucleotide": {"range": [{"begin": str(begin), "end": str(end)}]},
           "protein": {"range": [{"begin": "1", "end": str((end - begin + 1) // 3 - 1)}]}}
    if peptides:
        cds["mature_peptide"] = [
            {"name": n, "nucleotide": {"range": [{"begin": str(b), "end": str(e)}]}}
            for n, b, e in peptides
        ]
    return cds


def _annotation(acc, *cds):
    return {"accession": acc, "genes": [{"cds": [c]} for c in cds]}


# -- protein names ---------------------------------------------------------------------

@pytest.mark.parametrize("name, key", [
    ("spike protein", "spike glycoprotein"),
    ("Spike glycoprotein", "spike glycoprotein"),
    ("putative coat protein", "capsid protein"),
    ("replication-associated protein", "replication-associated protein (Rep)"),
    ("Rep", "replication-associated protein (Rep)"),
    ("nucleocapsid phosphoprotein", "nucleocapsid protein"),
    ("major capsid protein; VP5", "major capsid protein"),
    ("CDS: dUTPase; the Herpesviridae are ...", "dutpase"),
    ("NS5-like protein", "ns5"),
    ("ORF1ab polyprotein", "ORF1ab polyprotein"),
    ("hypothetical protein", None),
    ("CDS: ", None),
    ("ORF", None),
])
def test_protein_key(name, key):
    assert genomics.protein_key(name) == key


# -- genomes ---------------------------------------------------------------------------

def test_segments_form_one_genome_unlabelled_records_stay_apart():
    records = [
        _record("NC_1", 1, 2700, "DNA-A"), _record("NC_2", 1, 2600, "DNA-B"),
        _record("NC_3", 1, 2710, "DNA-A"),                  # a second RefSeq of the virus
        _record("NC_4", 2, 28000), _record("NC_5", 2, 27900),  # two complete RefSeqs
    ]
    groups = genomics.group_genomes(records)
    sizes = sorted([sorted(r["accession"] for r in g) for g in groups.values()])
    assert sizes == [["NC_1", "NC_2"], ["NC_3"], ["NC_4"], ["NC_5"]]


def test_summary_statistics_and_set_asides():
    records = [_record(f"NC_{i}", i, n, host="Homo sapiens" if i < 3 else None)
               for i, n in enumerate([9000, 9500, 10000, 10500, 11000], 1)]
    records += [_record("NC_short", 6, 1000),                       # single-gene record
                _record("NC_part", 7, 9800, completeness="PARTIAL")]
    s = genomics.summarize(42, records, [], total=10)
    assert (s.genomes, s.short, s.partial) == (5, 1, 1)
    assert (s.size_min, s.size_q1, s.size_median, s.size_q3, s.size_max) == (9000, 9500, 10000, 10500, 11000)
    assert s.segments == {1: 5}
    assert sum(b.count for b in s.histogram) == 5 and not s.log_scale
    assert s.smallest.accessions == ["NC_1"] and s.largest.size == 11000
    assert s.hosts == {"Homo sapiens": 2} and s.hosts_recorded == 2
    assert s.sampled and (s.records_used, s.refseq_records) == (7, 10)
    assert s.proteins == [] and s.proteins_per_genome is None


def test_too_few_genomes_is_no_summary():
    assert genomics.summarize(1, [_record("NC_1", 1, 9000), _record("NC_2", 2, 9100)], []) is None
    assert genomics.summarize(1, [], []) is None


def test_wide_size_range_uses_a_log_histogram():
    records = [_record(f"NC_{i}", i, n) for i, n in enumerate([20_000, 40_000, 60_000, 300_000], 1)]
    s = genomics.summarize(1, records, [])
    assert s.log_scale and len(s.histogram) == genomics.HISTOGRAM_BINS
    assert s.histogram[0].low == 20_000 and s.histogram[-1].high == 300_000


def test_proteins_counted_once_per_genome_with_synonyms():
    records = [_record(f"NC_{i}", i, 10_000) for i in (1, 2, 3, 4)]
    annotations = [
        _annotation("NC_1", _cds("spike protein", 1, 3000), _cds("membrane protein", 3001, 3600)),
        _annotation("NC_2", _cds("spike glycoprotein", 1, 3300), _cds("hypothetical protein", 4000, 4300)),
        _annotation("NC_3", _cds("S protein", 1, 3000), _cds("membrane glycoprotein", 3001, 3300)),
        _annotation("NC_4", _cds("polyprotein", 1, 9000)),
    ]
    s = genomics.summarize(1, records, annotations)
    by_name = {p.name: p for p in s.proteins}
    assert by_name["spike glycoprotein"].genomes == 2          # NC_3's "S protein" is not merged
    assert by_name["spike glycoprotein"].share == 0.5
    assert by_name["spike glycoprotein"].median_aa == 1049      # longest copy per genome
    assert "hypothetical protein" not in by_name
    assert by_name["membrane protein"].genomes == 2
    assert [p.name for p in s.proteins[:2]] == ["membrane protein", "spike glycoprotein"]
    assert s.proteins_per_genome == (1, 2, 2)
    assert s.polyprotein_genomes == 1
    assert s.coding_density == pytest.approx(0.36, abs=0.01)


def test_mature_peptides_in_genome_order():
    records = [_record(f"NC_{i}", i, 10_800) for i in (1, 2, 3)]
    peptides = [("nonstructural protein NS1", 2400, 3450), ("capsid protein C", 100, 450),
                ("envelope protein E", 900, 2400)]
    annotations = [_annotation(f"NC_{i}", _cds("polyprotein", 100, 10_300, peptides)) for i in (1, 2, 3)]
    s = genomics.summarize(1, records, annotations)
    assert [p.name for p in s.mature_peptides] == [
        "capsid protein C", "envelope protein E", "nonstructural protein NS1"]
    assert s.mature_peptide_genomes == 3
    assert s.mature_peptides[0].median_aa == 117               # from the coding span
    assert s.polyprotein_genomes == 3


# -- topology and genetic code ----------------------------------------------------------

def _meta(acc, topology, code):
    from viralfetch.ncbi import SeqMeta
    return SeqMeta(acc, "org", 1, "dna", "genomic", topology, "complete", "refseq", "2020",
                   geneticcode=code)


def test_apply_metadata_counts_genomes():
    records = [_record("NC_1", 1, 2700, "DNA-A"), _record("NC_2", 1, 2600, "DNA-B"),
               _record("NC_3", 2, 2800), _record("NC_4", 3, 2750), _record("NC_5", 4, 2900)]
    s = genomics.summarize(1, records, [])
    assert sorted(genomics.metadata_accessions(s)) == ["NC_1", "NC_2", "NC_3", "NC_4", "NC_5"]
    genomics.apply_metadata(s, [
        _meta("NC_1.1", "circular", "11"), _meta("NC_2.1", "linear", "11"),   # one genome, mixed
        _meta("NC_3.1", "circular", "11"), _meta("NC_4.2", "not-set", "1"),
        # NC_5: no summary -> not counted
    ])
    assert s.metadata_genomes == 3
    assert s.topology == {"mixed": 1, "circular": 1, "not recorded": 1}
    assert s.genetic_codes == {"11": 2, "1": 1}
    assert genomics.genetic_code_label("11") == "11 (bacterial, archaeal, plastid)"
    assert genomics.genetic_code_label("99") == "99"


def test_profile_genomes_topology_and_code(vmr, ncbi_knows):
    ncbi = FakeNCBI()
    p = profile.build(vmr, "Alphavirus", ncbi=ncbi, datasets=FakeDatasets(), viralzone=_viralzone())
    assert ncbi.summarised == ["NC_000001.1", "NC_000002.1", "NC_000003.1"]  # full-length only
    assert p.genomes.topology == {"linear": 2, "circular": 1}
    assert p.genomes.genetic_codes == {"1": 3}
    out = _render(p)
    assert "linear 67% · circular 33%" in out
    assert "1 (standard) in every genome" in out
    row = profile_row(p)
    assert (row["topology"], row["genetic_codes"]) == ("linear:2; circular:1", "1:3")


def test_metadata_failure_is_a_warning(vmr, ncbi_knows):
    from viralfetch.ncbi import NCBIError

    class NCBI(FakeNCBI):
        def esummary_nuccore(self, accessions):
            raise NCBIError("boom")

    p = profile.build(vmr, "Alphavirus", ncbi=NCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    assert p.genomes.size_median == 9500 and p.genomes.topology == {}
    assert any("genetic code" in w for w in p.warnings)


def test_json_profile_omits_genome_accessions(vmr, ncbi_knows, capsys):
    from viralfetch.render import json_
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    json_.profile(p)
    out = capsys.readouterr().out
    assert '"genetic_codes"' in out and "genome_accessions" not in out


# -- Datasets paging -------------------------------------------------------------------

class PagedSession:
    def __init__(self, pages):
        self.pages = pages
        self.urls: list[str] = []

    def get(self, url, headers=None, timeout=None):
        self.urls.append(url)
        page = 1 if "page_token=t1" in url else 2 if "page_token=t2" in url else 0
        return FakeResponse(200, self.pages[page])


def _client(session) -> DatasetsClient:
    return DatasetsClient(config_mod.Config(email="a@b.c", api_key=None), session=session,
                          sleep=lambda s: None)


def test_refseq_reports_follows_pages_up_to_the_cap():
    pages = [
        {"total_count": 5, "reports": [{"accession": "a"}, {"accession": "b"}], "next_page_token": "t1"},
        {"total_count": 5, "reports": [{"accession": "c"}, {"accession": "d"}], "next_page_token": "t2"},
        {"total_count": 5, "reports": [{"accession": "e"}]},
    ]
    session = PagedSession(pages)
    total, reports = _client(session).refseq_reports(11050, "annotation_report")
    assert total == 5 and [r["accession"] for r in reports] == list("abcde")
    assert all("/annotation_report?" in u and "filter.refseq_only=true" in u for u in session.urls)

    session = PagedSession(pages)
    total, reports = _client(session).refseq_reports(11050, max_records=3)
    assert total == 5 and len(reports) == 4 and len(session.urls) == 2
    assert "page_size=1" in session.urls[1]                      # only what the cap leaves


# -- profile ---------------------------------------------------------------------------

def test_species_profile_has_no_genome_summary(vmr, ncbi_knows):
    class Datasets(FakeDatasets):
        def refseq_reports(self, taxid, kind="dataset_report"):
            raise AssertionError("not for a species")

    p = profile.build(vmr, "Alphavirus one", ncbi=FakeNCBI(), datasets=Datasets(),
                      viralzone=_viralzone())
    assert p.genomes is None


def test_annotation_failure_keeps_the_sizes(vmr, ncbi_knows):
    from viralfetch.ncbi import NCBIError

    class Datasets(FakeDatasets):
        def refseq_reports(self, taxid, kind="dataset_report"):
            if kind == "annotation_report":
                raise NCBIError("boom")
            return super().refseq_reports(taxid, kind)

    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=Datasets(), viralzone=_viralzone())
    assert p.genomes.size_median == 9500 and p.genomes.proteins == []
    assert any("annotation" in w for w in p.warnings)


def test_render_genomes_panel(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    out = _render(p)
    assert "Genomes · 3 RefSeq genomes" in out
    assert "9.50 kb median" in out
    assert "capsid protein" in out and "100%" in out
    assert "set aside: 1 partial, 1 shorter than ⅓ of the median" in out
    full = _render(p, full=True)
    assert "smallest" in full and "NC_000001.1" in full


def test_profile_row_has_genome_columns(vmr, ncbi_knows):
    p = profile.build(vmr, "Alphavirus", ncbi=FakeNCBI(), datasets=FakeDatasets(), viralzone=_viralzone())
    row = profile_row(p)
    assert (row["refseq_genomes"], row["genome_size_median"]) == (3, 9500)
    assert row["common_proteins"].startswith("capsid protein (100%)")
