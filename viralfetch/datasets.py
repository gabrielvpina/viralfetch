"""Minimal NCBI Datasets v2 client: RefSeq listings and taxonomy reports.

E-utilities (:mod:`viralfetch.ncbi`) resolve names and lineages and count
records (``esearch`` stays fast even for millions of them). Datasets adds the
RefSeq sequences of a taxon with their metadata and annotation (the genome
summary of :mod:`viralfetch.genomics`), and the taxonomy report's ICTV
name status, genomic molecule type and assembly count. Its own record counts
are not used: on large taxa they take tens of seconds or fail outright.

Responses are cached in the TEXTS namespace for a week: these are live
summaries that grow with every GenBank release, unlike accession records.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from urllib.parse import urlencode

import requests

from . import __version__
from .cache import COUNTS_TTL, TEXTS, Cache
from .config import Config
from .ictv import REPO_URL
from .ncbi import NCBIError, RateLimiter

BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
VIRUS_PORTAL = "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:{taxid}"
_RETRY_STATUS = {429, 500, 502, 503, 504}
REPORT_PAGE = 1000
MAX_REPORT_RECORDS = 3000  # big taxa (an order, a class) are summarised from a sample


@dataclass
class GenomeRecord:
    """One NCBI Virus nucleotide record (a RefSeq genome, typically)."""

    accession: str
    organism: str
    taxid: int | None
    length: int | None
    completeness: str
    segment: str | None
    source_database: str
    release_date: str
    update_date: str


@dataclass
class TaxonReport:
    """The parts of a Datasets taxonomy report the profile shows."""

    taxid: int
    name: str
    rank: str
    ictv_status: str | None = None     # e.g. "Name is currently accepted by the ICTV."
    genomic_moltype: str | None = None
    assemblies: int | None = None
    children: int = 0


@dataclass
class SequenceStats:
    """NCBI nucleotide records under a taxid (descendants included).

    ``complete`` counts records titled "complete genome" (Entrez has no
    completeness field; NCBI Virus' own flag agrees within a few percent but
    cannot be counted quickly on large taxa). ``None`` = that query failed.
    """

    taxid: int
    total: int | None
    complete: int | None
    refseq: int | None
    refseq_genomes: list[GenomeRecord] = field(default_factory=list)
    portal_url: str = ""


class DatasetsClient:
    def __init__(
        self,
        config: Config,
        cache: Cache | None = None,
        session: requests.Session | None = None,
        sleep=time.sleep,
        max_retries: int = 4,
    ):
        self.api_key = config.api_key
        self.cache = cache
        self.session = session or requests.Session()
        self.max_retries = max_retries
        self._limiter = RateLimiter(config.rate_limit, sleep=sleep)
        self._sleep = sleep
        self.user_agent = f"viralfetch/{__version__} (+{REPO_URL})"

    def _get(self, path: str, params: dict | None = None, *, retries: int | None = None) -> dict:
        query = urlencode(sorted((params or {}).items()))
        url = f"{BASE_URL}/{path}" + (f"?{query}" if query else "")
        if self.cache:
            cached = self.cache.get(TEXTS, f"datasets:{url}", ttl=COUNTS_TTL)
            if cached is not None:
                return json.loads(cached)
        headers = {"User-Agent": self.user_agent, "Accept": "application/json"}
        if self.api_key:
            headers["api-key"] = self.api_key
        last_exc: Exception | None = None
        for attempt in range(retries or self.max_retries):
            self._limiter.wait()
            try:
                resp = self.session.get(url, headers=headers, timeout=60)
            except requests.RequestException as exc:
                last_exc = exc
                self._sleep(0.5 * (2 ** attempt))
                continue
            if resp.status_code in _RETRY_STATUS:
                last_exc = NCBIError(f"HTTP {resp.status_code} from Datasets {path}")
                self._sleep(0.5 * (2 ** attempt))
                continue
            if resp.status_code != 200:
                raise NCBIError(f"HTTP {resp.status_code} from Datasets {path}: {resp.text[:200]}")
            data = resp.json()
            if self.cache:
                self.cache.set(TEXTS, f"datasets:{url}", json.dumps(data))
            return data
        raise NCBIError(f"Datasets {path} failed after {retries or self.max_retries} attempts: {last_exc}")

    # -- virus ----------------------------------------------------------------

    def _virus_report(self, taxid: int | str, page_size: int, **filters: bool) -> dict:
        params: dict = {"page_size": page_size}
        for name, on in filters.items():
            if on:
                params[f"filter.{name}"] = "true"
        return self._get(f"virus/taxon/{taxid}/dataset_report", params, retries=2)

    def refseq_genomes(self, taxid: int | str, limit: int = 25) -> tuple[int, list[GenomeRecord]]:
        """``(total, first records)`` of the RefSeq sequences under ``taxid``,
        sorted by organism (Datasets' own order is arbitrary)."""
        data = self._virus_report(taxid, limit, refseq_only=True)
        return int(data.get("total_count") or 0), sorted_genomes(data.get("reports", []))

    def refseq_reports(
        self, taxid: int | str, kind: str = "dataset_report", max_records: int = MAX_REPORT_RECORDS,
    ) -> tuple[int, list[dict]]:
        """``(total, raw reports)`` of every RefSeq sequence under ``taxid``, up
        to ``max_records``; ``kind`` is ``dataset_report`` (length, completeness,
        host, ...) or ``annotation_report`` (genes, CDS, mature peptides)."""
        reports: list[dict] = []
        total, token = 0, None
        while len(reports) < max_records:
            params: dict = {"page_size": min(REPORT_PAGE, max_records - len(reports)),
                            "filter.refseq_only": "true"}
            if token:
                params["page_token"] = token
            data = self._get(f"virus/taxon/{taxid}/{kind}", params, retries=2)
            reports += data.get("reports") or []
            total = max(total, int(data.get("total_count") or 0))
            token = data.get("next_page_token")
            if not token or not data.get("reports"):
                break
        return max(total, len(reports)), reports

    # -- taxonomy -------------------------------------------------------------

    def taxonomy_report(self, taxid: int | str) -> TaxonReport | None:
        data = self._get(f"taxonomy/taxon/{taxid}/dataset_report")
        reports = data.get("reports") or []
        if not reports or "taxonomy" not in reports[0]:
            return None
        return _taxon_report(reports[0]["taxonomy"])


def _int(value) -> int | None:
    try:
        return int(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def sorted_genomes(reports: list[dict]) -> list[GenomeRecord]:
    """Dataset-report entries as :class:`GenomeRecord`, sorted by organism."""
    records = [_genome_record(r) for r in reports]
    records.sort(key=lambda g: (g.organism.casefold(), g.segment or "", g.accession))
    return records


def _genome_record(r: dict) -> GenomeRecord:
    virus = r.get("virus") or {}
    return GenomeRecord(
        accession=r.get("accession", ""),
        organism=virus.get("organism_name", ""),
        taxid=_int(virus.get("tax_id")),
        length=_int(r.get("length")),
        completeness=(r.get("completeness") or "").lower(),
        segment=r.get("segment") or None,
        source_database=r.get("source_database", ""),
        release_date=(r.get("release_date") or "")[:10],
        update_date=(r.get("update_date") or "")[:10],
    )


def _taxon_report(tax: dict) -> TaxonReport:
    name_info = tax.get("current_scientific_name") or {}
    status = None
    for note in name_info.get("notes") or []:
        if note.get("name") == "ICTV Status":
            status = note.get("note")
    assemblies = None
    for count in tax.get("counts") or []:
        if count.get("type") == "COUNT_TYPE_ASSEMBLY":
            assemblies = _int(count.get("count"))
    return TaxonReport(
        taxid=_int(tax.get("tax_id")) or 0,
        name=name_info.get("name", ""),
        rank=(tax.get("rank") or "").lower(),
        ictv_status=status,
        genomic_moltype=tax.get("genomic_moltype"),
        assemblies=assemblies,
        children=len(tax.get("children") or []),
    )
