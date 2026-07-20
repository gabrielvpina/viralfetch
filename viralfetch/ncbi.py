"""NCBI E-utilities client.

Hard rules enforced here (SPEC section 3):

- Only the E-utilities endpoints are used — never NCBI HTML scraping.
- Every request carries ``tool=viralfetch`` and a real ``email`` (no default;
  a missing email is a hard error, raised earlier by :meth:`Config.require_email`).
- One central rate limiter (3 req/s without an API key, 10 with one) — no
  ``sleep()`` scattered around.
- Requests are POSTed (accession lists blow past URL length limits).
- Accessions are sent in batches of at most 200.
- Partial failure is normal: missing accessions are reported, never swallowed
  and never raised.

Immutable data (accession metadata and sequences) is cached per-accession, so
overlapping requests reuse prior fetches.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field

import requests

from .cache import SEQS, Cache
from .config import Config

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_BATCH = 200
_RETRY_STATUS = {429, 500, 502, 503, 504}


class NCBIError(Exception):
    """A non-recoverable error talking to E-utilities."""


@dataclass
class SeqMeta:
    accession: str
    organism: str
    length: int | None
    moltype: str
    biomol: str
    topology: str
    completeness: str
    sourcedb: str
    updatedate: str


@dataclass
class MetaResult:
    records: list[SeqMeta] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)


@dataclass
class RecordsResult:
    rettype: str
    text: str = ""
    returned: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)


def _base_accession(acc: str) -> str:
    """Accession without a version suffix (``MN908947.3`` -> ``MN908947``)."""
    return acc.split(".", 1)[0]


class RateLimiter:
    """Enforces a minimum interval between requests, process-wide.

    Single central instance per client; there are no ad-hoc sleeps elsewhere.
    """

    def __init__(self, per_second: int, sleep=time.sleep, clock=time.monotonic):
        self._min_interval = 1.0 / per_second
        self._sleep = sleep
        self._clock = clock
        self._last = 0.0

    def wait(self) -> None:
        now = self._clock()
        gap = self._min_interval - (now - self._last)
        if gap > 0:
            self._sleep(gap)
        self._last = self._clock()


class NCBIClient:
    def __init__(
        self,
        config: Config,
        cache: Cache | None = None,
        session: requests.Session | None = None,
        sleep=time.sleep,
        max_retries: int = 4,
    ):
        self.email = config.require_email()  # fails loudly if unset
        self.api_key = config.api_key
        self.tool = "viralfetch"
        self.cache = cache
        self.session = session or requests.Session()
        self.max_retries = max_retries
        self._limiter = RateLimiter(config.rate_limit, sleep=sleep)
        self._sleep = sleep

    # -- transport ---------------------------------------------------------

    def _params(self, extra: dict) -> dict:
        params = {"tool": self.tool, "email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(extra)
        return params

    def _post(self, endpoint: str, params: dict) -> str:
        """Rate-limited, retried POST returning the response body as text."""
        url = f"{BASE_URL}/{endpoint}"
        data = self._params(params)
        last_exc: Exception | None = None
        for attempt in range(self.max_retries):
            self._limiter.wait()
            try:
                resp = self.session.post(url, data=data, timeout=60)
            except requests.RequestException as exc:
                last_exc = exc
                self._backoff(attempt)
                continue
            if resp.status_code in _RETRY_STATUS:
                last_exc = NCBIError(f"HTTP {resp.status_code} from {endpoint}")
                self._backoff(attempt)
                continue
            if resp.status_code != 200:
                raise NCBIError(f"HTTP {resp.status_code} from {endpoint}: {resp.text[:200]}")
            return resp.text
        raise NCBIError(f"{endpoint} failed after {self.max_retries} attempts: {last_exc}")

    def _backoff(self, attempt: int) -> None:
        self._sleep(0.5 * (2 ** attempt))

    @staticmethod
    def _batches(items: list[str]):
        for i in range(0, len(items), MAX_BATCH):
            yield items[i:i + MAX_BATCH]

    # -- esummary (metadata) ----------------------------------------------

    def esummary_nuccore(self, accessions: list[str]) -> MetaResult:
        """Fetch nuccore metadata for accessions, caching each permanently."""
        result = MetaResult()
        to_fetch: list[str] = []
        for acc in accessions:
            cached = self.cache.get(SEQS, f"summary:{acc}") if self.cache else None
            if cached is not None:
                result.records.append(_meta_from_json(json.loads(cached)))
            else:
                to_fetch.append(acc)

        for batch in self._batches(to_fetch):
            body = self._post(
                "esummary.fcgi",
                {"db": "nuccore", "id": ",".join(batch), "retmode": "json"},
            )
            found = _parse_esummary(body)  # base_accession -> record dict
            for acc in batch:
                record = found.get(_base_accession(acc))
                if record is None:
                    result.missing.append(acc)
                    continue
                if self.cache:
                    self.cache.set(SEQS, f"summary:{acc}", json.dumps(record))
                result.records.append(_meta_from_json(record))
        return result

    # -- efetch (fasta / gb) ----------------------------------------------

    def efetch_nuccore(self, accessions: list[str], rettype: str) -> RecordsResult:
        """Fetch fasta/gb records for accessions, caching each permanently."""
        result = RecordsResult(rettype=rettype)
        parts: list[str] = []
        to_fetch: list[str] = []
        for acc in accessions:
            cached = self.cache.get(SEQS, f"{rettype}:{acc}") if self.cache else None
            if cached is not None:
                parts.append(cached)
                result.returned.append(acc)
            else:
                to_fetch.append(acc)

        for batch in self._batches(to_fetch):
            body = self._post(
                "efetch.fcgi",
                {"db": "nuccore", "id": ",".join(batch), "rettype": rettype, "retmode": "text"},
            )
            records = _split_records(body, rettype)  # base_accession -> record text
            for acc in batch:
                record = records.get(_base_accession(acc))
                if record is None:
                    result.missing.append(acc)
                    continue
                if self.cache:
                    self.cache.set(SEQS, f"{rettype}:{acc}", record)
                parts.append(record)
                result.returned.append(acc)

        result.text = "".join(parts)
        return result


# -- parsing helpers -------------------------------------------------------

def _parse_esummary(body: str) -> dict[str, dict]:
    """Map base-accession -> esummary record from a JSON esummary body."""
    data = json.loads(body)
    result = data.get("result", {})
    out: dict[str, dict] = {}
    for uid in result.get("uids", []):
        record = result.get(uid)
        if not record:
            continue
        caption = record.get("caption") or _base_accession(record.get("accessionversion", ""))
        if caption:
            out[caption] = record
    return out


def _meta_from_json(record: dict) -> SeqMeta:
    slen = record.get("slen")
    try:
        length = int(slen) if slen not in (None, "") else None
    except (TypeError, ValueError):
        length = None
    return SeqMeta(
        accession=record.get("accessionversion") or record.get("caption", ""),
        organism=record.get("organism", ""),
        length=length,
        moltype=record.get("moltype", ""),
        biomol=record.get("biomol", ""),
        topology=record.get("topology", ""),
        completeness=record.get("completeness", ""),
        sourcedb=record.get("sourcedb", ""),
        updatedate=record.get("updatedate", ""),
    )


_FASTA_HEADER = re.compile(r"^>(\S+)")
_GB_VERSION = re.compile(r"^VERSION\s+(\S+)", re.MULTILINE)
_GB_ACCESSION = re.compile(r"^ACCESSION\s+(\S+)", re.MULTILINE)


def _split_records(body: str, rettype: str) -> dict[str, dict]:
    """Split a concatenated efetch body into per-accession record text."""
    if rettype == "fasta":
        return _split_fasta(body)
    return _split_genbank(body)


def _split_fasta(body: str) -> dict[str, str]:
    out: dict[str, str] = {}
    current_key: str | None = None
    lines: list[str] = []
    for line in body.splitlines(keepends=True):
        if line.startswith(">"):
            if current_key is not None:
                out[current_key] = "".join(lines)
            m = _FASTA_HEADER.match(line)
            current_key = _base_accession(m.group(1)) if m else None
            lines = [line]
        elif current_key is not None:
            lines.append(line)
    if current_key is not None:
        out[current_key] = "".join(lines)
    return out


def _split_genbank(body: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for chunk in body.split("//\n"):
        if not chunk.strip():
            continue
        record = chunk + "//\n"
        m = _GB_VERSION.search(record) or _GB_ACCESSION.search(record)
        if m:
            out[_base_accession(m.group(1))] = record
    return out
