"""Genome features of a taxon, summarised from its RefSeq genomes.

The profile of a genus, family or any higher taxon (``viralfetch Flaviviridae``)
gets a "Genomes" section built from two NCBI Datasets reports over the RefSeq
sequences under the taxon: the dataset report (length, completeness, segment,
host) and the annotation report (CDS and mature peptides).

The unit is a *genome*: one RefSeq record, or the segments of a segmented
virus summed (records of one organism labelled with distinct segments; NCBI
also keeps several complete RefSeqs of one organism, which stay separate
genomes). A genome with a partial segment is set aside, and so is one shorter than a third of the
median — RefSeq also holds single-gene records flagged "complete" (an NS5
fragment of a flavivirus, say) that would otherwise skew every number.

From the full-length genomes:

- size distribution (quartiles, range, a histogram) and segment counts;
- proteins per genome and coding density (share of the genome inside a CDS);
- how many genomes express a polyprotein;
- the proteins most genomes encode, their names normalised so that
  "spike protein" and "spike glycoprotein" count together, with median length;
- the mature peptides (cleavage products of polyproteins), when annotated;
- the hosts recorded on the RefSeq entries;
- topology (circular / linear) and genetic code (NCBI translation table),
  added by :func:`apply_metadata` from the records' nucleotide summaries.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

import math
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field

MIN_GENOMES = 3            # fewer full-length genomes: no distribution to speak of
SHORT_FRACTION = 1 / 3     # below this share of the median size: set aside as partial
PROTEIN_MIN_SHARE = 0.05   # proteins kept in the summary
PROTEIN_MAX = 30
MATURE_PEPTIDE_MIN_SHARE = 0.1  # of the genomes, to list mature peptides at all
HISTOGRAM_BINS = 40
POLYPROTEIN_SPAN = 0.6     # a CDS covering this share of its segment...
POLYPROTEIN_AA = 1000      # ...and this long, counts as a polyprotein
TOP_HOSTS = 5
METADATA_GENOMES = 1000    # topology / genetic code: from at most this many genomes

# NCBI translation tables (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
GENETIC_CODES = {
    "1": "standard", "2": "vertebrate mitochondrial", "3": "yeast mitochondrial",
    "4": "mold/protozoan mitochondrial, Mycoplasma", "5": "invertebrate mitochondrial",
    "6": "ciliate nuclear", "9": "echinoderm mitochondrial", "10": "euplotid nuclear",
    "11": "bacterial, archaeal, plastid", "12": "alternative yeast nuclear",
    "13": "ascidian mitochondrial", "14": "alternative flatworm mitochondrial",
    "15": "Blepharisma nuclear", "16": "chlorophycean mitochondrial",
    "25": "SR1 and Gracilibacteria", "26": "Pachysolen nuclear", "27": "karyorelict nuclear",
    "28": "Condylostoma nuclear", "29": "Mesodinium nuclear", "30": "peritrich nuclear",
    "31": "Blastocrithidia nuclear", "33": "Cephalodiscidae mitochondrial",
}


@dataclass
class ProteinFreq:
    name: str
    genomes: int          # full-length genomes encoding it
    share: float          # of the full-length genomes
    median_aa: int | None
    position: float | None = None   # median start, as a share of the genome length


@dataclass
class SizeBin:
    low: int
    high: int
    count: int


@dataclass
class GenomeExtreme:
    organism: str
    size: int
    accessions: list[str]


@dataclass
class GenomeSummary:
    """Genome features over the full-length RefSeq genomes under a taxon."""

    taxid: int
    refseq_records: int            # RefSeq sequences NCBI holds under the taxon
    records_used: int              # of those, the ones fetched (capped on big taxa)
    genomes: int                   # full-length genomes the numbers below describe
    partial: int = 0               # genomes set aside: a segment flagged partial
    short: int = 0                 # genomes set aside: < SHORT_FRACTION of the median
    size_min: int | None = None
    size_q1: int | None = None
    size_median: int | None = None
    size_q3: int | None = None
    size_max: int | None = None
    log_scale: bool = False        # histogram bins on a log scale
    histogram: list[SizeBin] = field(default_factory=list)
    smallest: GenomeExtreme | None = None
    largest: GenomeExtreme | None = None
    segments: dict[int, int] = field(default_factory=dict)          # segments -> genomes
    proteins_per_genome: tuple[int, int, int] | None = None         # min, median, max
    coding_density: float | None = None                              # median share
    polyprotein_genomes: int = 0
    proteins: list[ProteinFreq] = field(default_factory=list)
    mature_peptides: list[ProteinFreq] = field(default_factory=list)
    mature_peptide_genomes: int = 0
    hosts: dict[str, int] = field(default_factory=dict)              # top TOP_HOSTS
    hosts_recorded: int = 0                                          # genomes with a host
    topology: dict[str, int] = field(default_factory=dict)           # circular/linear/mixed/not recorded
    genetic_codes: dict[str, int] = field(default_factory=dict)      # table number -> genomes
    metadata_genomes: int = 0                                         # genomes the two above cover
    # Accessions of each full-length genome (input to apply_metadata; not output).
    genome_accessions: list[list[str]] = field(default_factory=list, repr=False)

    @property
    def sampled(self) -> bool:
        return self.records_used < self.refseq_records


# -- protein names --------------------------------------------------------------------

_PREFIX = re.compile(r"^(putative|probable|predicted|possible|truncated|similar to|"
                     r"conserved|uncharacterised|uncharacterized)\s+")
_SUFFIX = re.compile(r"\s+(protein|gene product|polypeptide|precursor)$")
_JUNK = re.compile(r"^(hypothetical|uncharacteri[sz]ed|unknown|unnamed|orf|cds|"
                   r"novel|predicted)\b|^\W*$")

# Normalised name -> the label it is counted and shown under. Only names that
# mean one thing across families: no single letters ("M" is the matrix of a
# paramyxovirus but the membrane protein of a coronavirus).
_SYNONYMS: dict[str, str] = {}
for label, names in {
    "capsid protein": ["capsid", "coat", "cp", "coat protein", "capsid protein vp1"],
    "major capsid protein": ["major capsid", "mcp"],
    "spike glycoprotein": ["spike", "spike glycoprotein", "s glycoprotein", "surface glycoprotein s"],
    "membrane protein": ["membrane", "membrane glycoprotein"],
    "envelope protein": ["envelope", "small envelope", "env", "envelope glycoprotein",
                         "envelope polyprotein", "env polyprotein", "envelope glycoprotein e"],
    "nucleocapsid protein": ["nucleocapsid", "nucleoprotein", "nucleocapsid phosphoprotein"],
    "matrix protein": ["matrix"],
    "phosphoprotein": ["phosphoprotein"],
    "fusion protein": ["fusion", "fusion glycoprotein"],
    "attachment glycoprotein": ["attachment glycoprotein", "attachment", "cell attachment",
                                "glycoprotein g"],
    "hemagglutinin-neuraminidase": ["hemagglutinin-neuraminidase", "haemagglutinin-neuraminidase",
                                    "hn", "hn glycoprotein"],
    "RNA-dependent RNA polymerase": ["rna-dependent rna polymerase", "rna dependent rna polymerase",
                                     "rna-directed rna polymerase", "rdrp", "rna polymerase",
                                     "large polymerase", "large structural", "polymerase",
                                     "l protein", "rna-dependent rna polymerase l"],
    "replication-associated protein (Rep)": ["replication-associated", "replication associated",
                                             "rep", "replication initiator",
                                             "replication initiation", "replication-associated rep"],
    "replication enhancer (REn)": ["replication enhancer", "replication enhancement", "ren"],
    "transcriptional activator (TrAP)": ["transcriptional activator", "transcription activator",
                                         "transactivator", "trap"],
    "pre-coat protein": ["pre-coat", "precoat"],
    "movement protein": ["movement", "mp", "cell-to-cell movement"],
    "nuclear shuttle protein": ["nuclear shuttle"],
    "Gag": ["gag", "gag polyprotein"],
    "Pol": ["pol", "pol polyprotein"],
    "Gag-Pol": ["gag-pol", "gag-pol polyprotein", "gag/pol", "gag-pol fusion"],
    "ORF1ab polyprotein": ["orf1ab polyprotein", "orf1ab", "polyprotein 1ab", "pp1ab",
                           "replicase polyprotein 1ab", "1ab polyprotein"],
    "ORF1a polyprotein": ["orf1a polyprotein", "orf1a", "polyprotein 1a", "pp1a",
                          "replicase polyprotein 1a", "1a polyprotein"],
    "polyprotein": ["polyprotein", "genome polyprotein", "flavivirus polyprotein",
                    "polyprotein precursor"],
    "terminase large subunit": ["terminase large subunit", "large terminase subunit",
                                "terminase large", "large terminase", "terl"],
    "terminase small subunit": ["terminase small subunit", "small terminase subunit",
                                "terminase small", "small terminase", "ters"],
    "portal protein": ["portal", "portal protein", "capsid portal"],
    "endolysin": ["endolysin", "lysin"],
}.items():
    for n in names:
        _SYNONYMS[n] = label


def protein_key(name: str) -> str | None:
    """The label a protein name is counted under; ``None`` for names that say
    nothing ("hypothetical protein", "ORF", an empty "CDS:")."""
    # "major capsid protein; VP5", "CDS: dUTPase; the Herpesviridae are ...":
    # the name is what comes before the first ";".
    n = re.sub(r"^cds:\s*", "", name.strip(), flags=re.IGNORECASE).split(";")[0]
    n = " ".join(n.casefold().replace("_", " ").split()).strip(" :;,.")
    previous = None
    while n != previous:
        previous = n
        n = _PREFIX.sub("", n)
        n = _SUFFIX.sub("", n)
        n = re.sub(r"-like$", "", n).strip()
    if _JUNK.search(n):
        return None
    return _SYNONYMS.get(n, n)


def _display_names(spellings: dict[str, Counter]) -> dict[str, str]:
    """Each key's label: the synonym label, else its most frequent spelling."""
    labels = set(_SYNONYMS.values())
    out = {}
    for key, counter in spellings.items():
        if key in labels:
            out[key] = key
        else:
            name = counter.most_common(1)[0][0]
            name = re.sub(r"^cds:\s*", "", name.strip(), flags=re.IGNORECASE).split(";")[0].strip()
            out[key] = _PREFIX.sub("", name) or name
    return out


# -- the summary ----------------------------------------------------------------------

def _int(value) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _ranges(location: dict | None) -> list[tuple[int, int]]:
    out = []
    for r in (location or {}).get("range") or []:
        a, b = _int(r.get("begin")), _int(r.get("end"))
        if a is not None and b is not None:
            out.append((min(a, b), max(a, b)))
    return out


def _covered(intervals: list[tuple[int, int]]) -> int:
    total, end = 0, 0
    for a, b in sorted(intervals):
        if b <= end:
            continue
        total += b - max(a, end + 1) + 1
        end = b
    return total


def _aa_length(cds: dict) -> int | None:
    """Protein length; from the coding span when the protein has no range."""
    spans = _ranges(cds.get("protein"))
    if spans:
        return sum(b - a + 1 for a, b in spans)
    return sum(b - a + 1 for a, b in _ranges(cds.get("nucleotide"))) // 3 or None


def _start(location: dict | None) -> int | None:
    return min((a for a, _b in _ranges(location)), default=None)


Seen = dict[str, tuple[int | None, float | None]]   # key -> (aa length, relative position)


def _note(seen: Seen, key: str, aa: int | None, start: int | None, whole: int | None) -> None:
    """Record one protein of a genome, keeping its longest copy."""
    position = round((start - 1) / whole, 4) if start and whole else None
    if key not in seen or (aa or 0) > (seen[key][0] or 0):
        seen[key] = (aa, position)


def _histogram(sizes: list[int]) -> tuple[bool, list[SizeBin]]:
    lo, hi = min(sizes), max(sizes)
    if lo == hi:
        return False, [SizeBin(lo, hi, len(sizes))]
    log = hi / lo > 4
    f = math.log10 if log else float
    inv = (lambda x: 10 ** x) if log else (lambda x: x)
    a, b = f(lo), f(hi)
    width = (b - a) / HISTOGRAM_BINS
    counts = [0] * HISTOGRAM_BINS
    for s in sizes:
        counts[min(int((f(s) - a) / width), HISTOGRAM_BINS - 1)] += 1
    bins = [SizeBin(round(inv(a + i * width)), round(inv(a + (i + 1) * width)), c)
            for i, c in enumerate(counts)]
    return log, bins


def _freqs(found: dict[str, list[tuple]], names: dict[str, str], n: int) -> list[ProteinFreq]:
    """``found``: key -> one ``(aa, position)`` per genome encoding it."""
    out = []
    for key, entries in found.items():
        share = len(entries) / n
        if share < PROTEIN_MIN_SHARE:
            continue
        aa = [x for x, _p in entries if x]
        pos = [p for _x, p in entries if p is not None]
        out.append(ProteinFreq(names[key], len(entries), round(share, 4),
                               round(statistics.median(aa)) if aa else None,
                               round(statistics.median(pos), 4) if pos else None))
    out.sort(key=lambda p: (-p.genomes, p.name.casefold()))
    return out[:PROTEIN_MAX]


def group_genomes(records: list[dict]) -> dict[str, list[dict]]:
    """Records grouped into genomes: a record without a segment label is a
    genome; labelled records of one organism form one, a repeated label
    starting another (two RefSeqs of the same segmented virus)."""
    genomes: dict[str, list[dict]] = {}
    by_organism: dict[str, list[list[dict]]] = defaultdict(list)
    for r in records:
        virus = r.get("virus") or {}
        organism = str(virus.get("tax_id") or virus.get("organism_name") or "")
        segment = (r.get("segment") or "").strip().casefold()
        if not segment or not organism:
            genomes[r.get("accession") or str(len(genomes))] = [r]
            continue
        sets = by_organism[organism]
        target = next((g for g in sets
                       if segment not in {(x.get("segment") or "").strip().casefold() for x in g}), None)
        if target is None:
            target = []
            sets.append(target)
        target.append(r)
    for organism, sets in by_organism.items():
        for i, g in enumerate(sets):
            genomes[f"{organism}#{i}"] = g
    return genomes


def summarize(
    taxid: int | str,
    records: list[dict],
    annotations: list[dict],
    *,
    total: int | None = None,
) -> GenomeSummary | None:
    """The genome summary of a taxon from its Datasets RefSeq reports.

    ``records`` are dataset-report entries, ``annotations`` annotation-report
    entries (matched by accession), ``total`` how many RefSeq records NCBI
    holds when only part was fetched. ``None`` when fewer than
    :data:`MIN_GENOMES` full-length genomes remain.
    """
    by_accession = {a.get("accession"): a for a in annotations}
    organisms = group_genomes(records)

    complete = {k: rs for k, rs in organisms.items()
                if all((r.get("completeness") or "").upper() == "COMPLETE" for r in rs)}
    sizes = {k: sum(_int(r.get("length")) or 0 for r in rs) for k, rs in complete.items()}
    sizes = {k: s for k, s in sizes.items() if s > 0}
    summary = GenomeSummary(taxid=int(taxid), refseq_records=max(total or 0, len(records)),
                            records_used=len(records), genomes=0,
                            partial=len(organisms) - len(complete))
    if not sizes:
        return None
    floor = statistics.median(sizes.values()) * SHORT_FRACTION
    full = {k: s for k, s in sizes.items() if s >= floor}
    summary.short = len(sizes) - len(full)
    summary.genomes = n = len(full)
    if n < MIN_GENOMES:
        return None

    values = sorted(full.values())
    q1, median, q3 = statistics.quantiles(values, n=4, method="inclusive")
    summary.size_min, summary.size_max = values[0], values[-1]
    summary.size_q1, summary.size_median, summary.size_q3 = round(q1), round(median), round(q3)
    summary.log_scale, summary.histogram = _histogram(values)
    for attr, pick in (("smallest", min), ("largest", max)):
        key = pick(full, key=lambda k: full[k])
        rs = organisms[key]
        setattr(summary, attr, GenomeExtreme(
            organism=(rs[0].get("virus") or {}).get("organism_name", ""),
            size=full[key], accessions=[r.get("accession", "") for r in rs]))
    summary.segments = dict(sorted(Counter(len(organisms[k]) for k in full).items()))
    summary.genome_accessions = [[r.get("accession", "") for r in organisms[k]] for k in full]

    proteins: dict[str, list[tuple]] = defaultdict(list)
    peptides: dict[str, list[tuple]] = defaultdict(list)
    spellings: dict[str, Counter] = defaultdict(Counter)
    per_genome, densities = [], []
    hosts: Counter[str] = Counter()
    for key in full:
        seen: Seen = {}
        seen_mp: Seen = {}
        cds_count = coded = 0
        polyprotein = False
        host = None
        for r in organisms[key]:
            host = host or ((r.get("host") or {}).get("organism_name"))
            length = _int(r.get("length")) or 0
            spans: list[tuple[int, int]] = []
            for gene in (by_accession.get(r.get("accession")) or {}).get("genes") or []:
                for cds in gene.get("cds") or []:
                    cds_count += 1
                    nt = _ranges(cds.get("nucleotide"))
                    spans += nt
                    aa = _aa_length(cds)
                    name = cds.get("name") or ""
                    span = sum(b - a + 1 for a, b in nt)
                    if "polyprotein" in name.casefold() or (
                        length and span >= POLYPROTEIN_SPAN * length and (aa or 0) >= POLYPROTEIN_AA
                    ):
                        polyprotein = True
                    k = protein_key(name)
                    if k:
                        spellings[k][name] += 1
                        _note(seen, k, aa, _start(cds.get("nucleotide")), length)
                    for mp in cds.get("mature_peptide") or []:
                        mk = protein_key(mp.get("name") or "")
                        if mk:
                            spellings[mk][mp.get("name")] += 1
                            _note(seen_mp, mk, _aa_length(mp), _start(mp.get("nucleotide")), length)
            coded += _covered(spans)
        for k, entry in seen.items():
            proteins[k].append(entry)
        for k, entry in seen_mp.items():
            peptides[k].append(entry)
        if seen_mp:
            summary.mature_peptide_genomes += 1
        if cds_count:
            per_genome.append(cds_count)
            densities.append(min(coded / full[key], 1.0))
        summary.polyprotein_genomes += polyprotein
        if host:
            hosts[host] += 1

    names = _display_names(spellings)
    summary.proteins = _freqs(proteins, names, n)
    if summary.mature_peptide_genomes >= max(MIN_GENOMES, MATURE_PEPTIDE_MIN_SHARE * n):
        summary.mature_peptides = sorted(  # in polyprotein order: C, prM, E, NS1, ...
            _freqs(peptides, names, summary.mature_peptide_genomes),
            key=lambda p: (p.position is None, p.position or 0))
    if per_genome:
        summary.proteins_per_genome = (min(per_genome), round(statistics.median(per_genome)),
                                       max(per_genome))
        summary.coding_density = round(statistics.median(densities), 4)
    summary.hosts = dict(hosts.most_common(TOP_HOSTS))
    summary.hosts_recorded = sum(hosts.values())
    return summary


def _base(accession: str) -> str:
    return accession.split(".")[0]


def apply_metadata(summary: GenomeSummary, metadata) -> None:
    """Fill topology and genetic code from nucleotide summaries.

    ``metadata`` is an iterable of records with ``accession``, ``topology``
    and ``geneticcode`` (:class:`viralfetch.ncbi.SeqMeta`). A genome is
    counted once its records are all there; a segmented genome mixing
    circular and linear segments is "mixed", and takes its segments' most
    frequent genetic code.
    """
    by_accession = {_base(m.accession): m for m in metadata}
    topology: Counter[str] = Counter()
    codes: Counter[str] = Counter()
    covered = 0
    for accessions in summary.genome_accessions[:METADATA_GENOMES]:
        metas = [by_accession.get(_base(a)) for a in accessions]
        if not metas or any(m is None for m in metas):
            continue
        covered += 1
        shapes = {(m.topology or "").lower() for m in metas} - {""}
        shapes = {"not recorded" if t == "not-set" else t for t in shapes}
        if shapes:
            topology[shapes.pop() if len(shapes) == 1 else "mixed"] += 1
        genome_codes = Counter(m.geneticcode for m in metas if m.geneticcode)
        if genome_codes:
            codes[genome_codes.most_common(1)[0][0]] += 1
    summary.topology = dict(topology.most_common())
    summary.genetic_codes = dict(codes.most_common())
    summary.metadata_genomes = covered


def metadata_accessions(summary: GenomeSummary) -> list[str]:
    """The accessions :func:`apply_metadata` needs, in genome order."""
    return [a for g in summary.genome_accessions[:METADATA_GENOMES] for a in g]


def genetic_code_label(code: str) -> str:
    name = GENETIC_CODES.get(code)
    return f"{code} ({name})" if name else code
