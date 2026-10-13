"""Build the ViralZone data set (``viralfetch/data/viralzone.json``).

ViralZone (https://viralzone.expasy.org, SIB Swiss Institute of Bioinformatics)
publishes one hand-curated *factsheet* per family, subfamily, genus and for a
few notable viruses. This module crawls them and keeps the descriptive,
complementary content — virion, genome, replication, host, geography, diseases,
transmission, vaccines, antivirals, host-virus interactions, reference strains
and cross-references — in a single JSON file. The per-protein UniProt listings
on each factsheet are left out: UniProt is their primary source.

The crawl starts at the taxonomy diagram on the home page, which links one page
per realm/phylum. Those pages lay out the whole tree (subphylum, class, order,
family, subfamily, genus, each genus tagged with its host groups) and link the
factsheets. A linked page is a factsheet when its ``<h1>`` carries an NCBI
taxid (``data-about="taxid:N"``), which is also the join key with NCBI and the
VMR, since ViralZone may lag one ICTV release behind.

The result lands in the user data dir, never in the installed package (the
bundled copy is refreshed by the maintainer with ``--out``); it is written to a
temporary file and only moved into place once complete and loadable.

Design note (SPEC section 4): nothing in this module prints. It returns data;
progress is reported through an optional callback.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.robotparser
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qs, unquote, urljoin, urlparse

import requests
from bs4 import BeautifulSoup, NavigableString, Tag

from . import __version__
from . import viralzone as vz_mod
from .config import Config
from .ictv import REPO_URL
from .ncbi import RateLimiter

TOOL_VERSION = f"viralfetch/{__version__}"
BASE_URL = vz_mod.BASE_URL

_RETRY_STATUS = {429, 500, 502, 503, 504}
# Below this many factsheets the site layout has most likely changed.
MIN_FACTSHEETS = 100

Progress = Callable[[int, int, dict], None]


class ViralZoneError(Exception):
    """A ViralZone page could not be fetched."""


class ViralZoneInstallError(Exception):
    """The rebuilt data set is unusable; the current one was left in place."""


# -- transport --------------------------------------------------------------

class ViralZoneClient:
    """Polite HTTP client for ViralZone: >= 1 s between requests, robots.txt
    honoured, identified User-Agent (the same policy as the ICTV client)."""

    def __init__(
        self,
        config: Config,
        session: requests.Session | None = None,
        sleep=time.sleep,
        clock=time.monotonic,
        max_retries: int = 3,
    ):
        email = config.require_email()  # fails loudly if unset
        self.session = session or requests.Session()
        self.max_retries = max_retries
        self._limiter = RateLimiter(1, sleep=sleep, clock=clock)
        self._sleep = sleep
        self._robots: urllib.robotparser.RobotFileParser | None = None
        self.user_agent = f"viralfetch/{__version__} (+{REPO_URL}; {email})"

    def _check_robots(self, url: str) -> None:
        if self._robots is None:
            rp = urllib.robotparser.RobotFileParser()
            try:
                body = self.session.get(
                    f"{BASE_URL}/robots.txt",
                    headers={"User-Agent": self.user_agent}, timeout=30,
                ).text
                rp.parse(body.splitlines())
            except requests.RequestException:
                rp.allow_all = True  # unreachable: do not invent rules
            self._robots = rp
        if not self._robots.can_fetch(self.user_agent, url):
            raise ViralZoneError(f"robots.txt disallows {url}")

    def get(self, path: str) -> str | None:
        """A page's HTML, or ``None`` on 404. Retries on 5xx / 429."""
        url = urljoin(BASE_URL + "/", path.lstrip("/"))
        self._check_robots(url)
        last_exc: Exception | None = None
        for attempt in range(self.max_retries):
            self._limiter.wait()
            try:
                resp = self.session.get(url, headers={"User-Agent": self.user_agent}, timeout=60)
            except requests.RequestException as exc:
                last_exc = exc
                self._sleep(0.5 * (2 ** attempt))
                continue
            if resp.status_code in _RETRY_STATUS:
                last_exc = ViralZoneError(f"HTTP {resp.status_code} from {url}")
                self._sleep(0.5 * (2 ** attempt))
                continue
            if resp.status_code == 404:
                return None
            if resp.status_code != 200:
                raise ViralZoneError(f"HTTP {resp.status_code} from {url}")
            return resp.content.decode("utf-8", errors="replace")
        raise ViralZoneError(f"{url} failed after {self.max_retries} attempts: {last_exc}")


# -- helpers ----------------------------------------------------------------

# Factsheet links come in several spellings that all redirect to ``/<id>``.
_PAGE_ID_RE = re.compile(r"^/(?:(?:preview_)?by_(?:species|protein)/)?(\d+)/?$")
_TAXID_RE = re.compile(r"taxid:(\d+)")
_WS_RE = re.compile(r"[ \t\r\f\v\xa0]+")
_BLOCKS = {"p", "div", "li", "ul", "ol", "tr", "table", "h1", "h2", "h3", "h4", "fieldset"}


def page_id(href: str | None) -> int | None:
    """The ViralZone page id an internal link points to, if any."""
    if not href:
        return None
    parsed = urlparse(href)
    if parsed.netloc and parsed.netloc != urlparse(BASE_URL).netloc:
        return None
    m = _PAGE_ID_RE.match(unquote(parsed.path))
    return int(m.group(1)) if m else None


def _soup(html: str) -> BeautifulSoup:
    return BeautifulSoup(html, "lxml")


def _clean_lines(text: str) -> str:
    lines = (_WS_RE.sub(" ", line).strip() for line in text.split("\n"))
    return "\n".join(line for line in lines if line)


def _text(node, skip: Callable[[Tag], bool] | None = None) -> str:
    """Readable text of a node: ``<br>`` and block elements become newlines,
    list items get a bullet, images and scripts are dropped."""
    parts: list[str] = []

    def walk(n) -> None:
        if isinstance(n, NavigableString):
            if type(n) is NavigableString:  # not a Comment / CData / ...
                parts.append(str(n).replace("\n", " "))
            return
        if not isinstance(n, Tag) or n.name in {"script", "style", "img", "svg", "button", "legend"}:
            return
        if skip is not None and skip(n):
            return
        if n.name == "br":
            parts.append("\n")
            return
        block = n.name in _BLOCKS
        if block:
            parts.append("\n")
        if n.name == "li":
            parts.append("- ")
        for child in n.children:
            walk(child)
        if block:
            parts.append("\n")

    walk(node)
    return _clean_lines("".join(parts))


def _key(title: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", title.casefold()).strip("_")


# -- realm / phylum pages ---------------------------------------------------

def parse_taxonomy_index(html: str) -> list[dict]:
    """``[{"id", "label"}]`` of the realm/phylum pages in the home-page diagram."""
    soup = _soup(html)
    diagram = soup.find(id="diagram")
    if diagram is None:
        return []
    out: list[dict] = []
    seen: set[int] = set()
    for a in diagram.find_all("a", href=True):
        pid = page_id(a["href"])
        if pid is None or pid in seen:
            continue
        seen.add(pid)
        out.append({"id": pid, "label": " ".join(a.get_text(" ").split())})
    return out


_TAXON_NAME_RE = re.compile(r"^[A-Z][a-z]+(?:-[a-z]+)?$")
_RANK_SUFFIXES = (
    ("viridae", "family"), ("satellitidae", "family"), ("viriformidae", "family"),
    ("viroidae", "family"), ("virinae", "subfamily"), ("satellitinae", "subfamily"),
    ("virales", "order"), ("viricetes", "class"), ("viricotina", "subphylum"),
    ("viricota", "phylum"),
)


def _rank_from_name(name: str) -> str | None:
    lowered = name.casefold()
    for suffix, rank in _RANK_SUFFIXES:
        if lowered.endswith(suffix):
            return rank
    if lowered.endswith(("virus", "viroid", "satellite", "viriform")):
        return "genus"
    return None


def _own_link(li: Tag) -> Tag | None:
    """The item's own link (before any nested list), if it has one."""
    for child in li.children:
        if isinstance(child, Tag) and child.name in {"ul", "ol", "li"}:
            return None
        if isinstance(child, Tag):
            a = child if child.name == "a" else child.find("a")
            if a is not None and a.get("href"):
                return a
    return None


def _own_label(li: Tag) -> str:
    label = ""
    for child in li.children:
        if isinstance(child, Tag) and child.name in {"ul", "ol", "li", "button"}:
            break
        label += child.get_text(" ") if isinstance(child, Tag) else str(child)
    return " ".join(label.replace("!", " ").split())


def parse_realm_page(html: str, page: int) -> list[dict]:
    """Taxa laid out on a realm/phylum page, in page order, without repeats.

    Families and subfamilies are the bold / underlined list items, genera the
    ``<li data-host="...">`` items (whose host groups are recorded); folded
    lists of families without a factsheet are tagged like genera, so the rank
    is read from the name's suffix. The ``family`` of a genus is the last family
    seen above it. Entries without a factsheet keep ``id = None``.
    """
    soup = _soup(html)
    out: list[dict] = []
    seen: set[tuple] = set()
    family: str | None = None
    subfamily: str | None = None
    for li in soup.find_all("li"):
        hosts = li.get("data-host")
        if hosts is None:
            head = li.find(["strong", "ins"])
            if head is None or head.find_parent("li") is not li:
                continue
            name = " ".join(head.get_text(" ").split())
        else:
            name = (_own_label(li).split(" ") or [""])[0]
        if not _TAXON_NAME_RE.match(name):
            continue
        rank = _rank_from_name(name)
        if rank not in {"family", "subfamily", "genus"}:
            continue
        if rank == "family":
            family, subfamily = name, None
        elif rank == "subfamily":
            subfamily = name
        key = (rank, name, family if rank != "family" else None)
        if key in seen:
            continue
        seen.add(key)
        link = _own_link(li)
        out.append({
            "name": name, "rank": rank,
            "id": page_id(link["href"]) if link is not None else None,
            "family": family if rank != "family" else None,
            "subfamily": subfamily if rank == "genus" else None,
            "hosts": sorted({h.strip() for h in (hosts or "").split(",") if h.strip()}),
            "page": page,
        })
    return out


def factsheet_candidates(html: str) -> list[int]:
    """Ids of every internal page linked from a realm page's taxonomy layout
    (description paragraphs and the "PFAM ..." buttons are not taxa)."""
    soup = _soup(html)
    ids: list[int] = []
    seen: set[int] = set()
    for a in soup.find_all("a", href=True):
        if "btn" in (a.get("class") or []) or a.find_parent(id="diagram") is not None:
            continue
        if a.find_parent(["header", "nav", "footer"]) is not None:
            continue
        pid = page_id(a["href"])
        if pid is not None and pid not in seen:
            seen.add(pid)
            ids.append(pid)
    return ids


# -- factsheets -------------------------------------------------------------

# Sidebar box titles vary between factsheets; map them onto stable keys.
_FIELD_ALIASES = {
    "host": "host", "natural_host": "host", "natural_hosts": "host", "hosts": "host",
    "tropism": "tropism", "cell_tropism": "tropism",
    "cell_receptor_s": "receptors", "cell_receptors": "receptors", "receptor": "receptors",
    "receptors": "receptors",
    "geography": "geography",
    "associated_diseases": "diseases", "associated_disease": "diseases", "diseases": "diseases",
    "transmission": "transmission",
    "vaccine": "vaccine", "vaccines": "vaccine",
    "antiviral_drugs": "antiviral_drugs", "antiviral_drug": "antiviral_drugs",
    "etymology": "etymology",
}
_SKIP_FIELDS = {"db_links", "taxonomy", "virus", "species", "interactions"}
_TEXT_FIELDS = ("host", "tropism", "receptors", "geography", "diseases",
                "transmission", "vaccine", "antiviral_drugs", "etymology")
_SECTION_ALIASES = {"gene_expression": "gene_expression", "virion": "virion",
                    "genome": "genome", "enzymes": "enzymes", "replication": "replication"}
# The replication box opens with the cellular site, spelt in several ways
# ("CYTOPLASMIC", "NUCLEUS", "CYTOPLASMIC/NUCLEAR", ...).
_SITE_WORD = r"(?:CYTOPLASM(?:IC)?|NUCLE(?:AR|US)|PERIPLASM(?:IC)?)"
_SITE_LINE_RE = re.compile(rf"^{_SITE_WORD}(?:\s*(?:/|,|AND|&)\s*{_SITE_WORD})*$")
_SITE_NAMES = {"CYTOPLASM": "CYTOPLASMIC", "NUCLEUS": "NUCLEAR", "PERIPLASM": "PERIPLASMIC"}
# Rank labels on the TAXONOMY box are not always singular / unhyphenated.
_RANK_LABELS = {"genera": "genus", "families": "family", "subfamilies": "subfamily",
                "sub_family": "subfamily", "sub_genus": "subgenus", "sub_order": "suborder"}
_FIELD_ALIASES_EXTRA = {"origin_of_name": "etymology"}
# A few factsheets carry raw Textile links ("label":http://...) in their text.
_TEXTILE_LINK_RE = re.compile(r'"([^"\n]+)":(?:https?|ftp)://[^\s\]\)]+')
_NONE_VALUES = {"", "none", "n/a", "unknown", "-"}


@dataclass
class _Box:
    title: str
    node: Tag


def _boxes(scope: Tag) -> list[_Box]:
    """Every ``<fieldset><legend>`` box under ``scope``, innermost legend wins."""
    out: list[_Box] = []
    for legend in scope.find_all("legend"):
        fs = legend.find_parent("fieldset")
        if fs is None:
            continue
        title = " ".join(legend.get_text(" ").split())
        if title:
            out.append(_Box(title, fs))
    return out


def _box_text(box: _Box) -> str:
    # A nested box (the t2 sidebar wraps several) must not leak into its parent.
    return _text(box.node, skip=lambda t: t.name == "fieldset" and t is not box.node)


def _sections(main: Tag) -> list[tuple[str, str, list[Tag]]]:
    """``(id, title, nodes)`` for each ``<h4>`` section of a tab's main column."""
    out: list[tuple[str, str, list[Tag]]] = []
    for h4 in main.find_all("h4"):
        nodes: list[Tag] = []
        for sib in h4.next_siblings:
            if isinstance(sib, Tag) and (sib.name in {"h4", "h2"} or sib.find("h4")):
                break
            nodes.append(sib)
        title = " ".join(h4.get_text(" ").split())
        out.append((h4.get("id") or _key(title), title, nodes))
    return out


def _nodes_text(nodes: list) -> str:
    wrapper = BeautifulSoup("<div></div>", "lxml").div
    for n in nodes:
        wrapper.append(_copy(n))
    return _text(wrapper)


def _copy(node):
    if isinstance(node, NavigableString):
        return NavigableString(str(node))
    return BeautifulSoup(str(node), "lxml").find(node.name) or NavigableString("")


def _images(nodes: list) -> list[str]:
    urls: list[str] = []
    for n in nodes:
        if not isinstance(n, Tag):
            continue
        imgs = [n] if n.name == "img" else n.find_all("img")
        for img in imgs:
            if "preview_image" in (img.get("class") or []):
                continue  # the PNG twin of an inline SVG
            src = img.get("src")
            if src:
                urls.append(urljoin(BASE_URL + "/", src))
    return urls


def _lineage(box: _Box | None, name: str) -> tuple[dict, dict, str | None]:
    """``(lineage, children, baltimore)`` from the TAXONOMY box.

    A rank row with several names (e.g. the genera of a family page) lists
    children, unless the page's own taxon is among them.
    """
    lineage: dict[str, str] = {}
    children: dict[str, list[str]] = {}
    baltimore = None
    if box is None:
        return lineage, children, baltimore
    first = box.node.find("a")
    group = box.node.find(string=re.compile(r"Group\s+[IVX]+"))
    if first is not None and re.match(r"Group\s+[IVX]+", first.get_text(" ").strip()):
        baltimore = " ".join(first.get_text(" ").split())
    elif group is not None:
        baltimore = " ".join(str(group).split())
    for tr in box.node.find_all("tr"):
        cells = tr.find_all("td")
        if len(cells) < 2:
            continue
        rank = cells[0].get_text(" ").strip().rstrip(":").strip().casefold()
        if not rank:
            continue
        # Names are links / bold, separated by (sometimes malformed) <br>s.
        names = [" ".join(t.get_text(" ").split()) for t in cells[1].find_all(["a", "strong", "b"])]
        names = [n for n in names if n] or [
            n for n in (x.lstrip("- ").strip() for x in _text(cells[1]).split("\n")) if n
        ]
        if len(names) == 1:
            lineage[rank] = names[0]
        elif name in names:  # the page's own rank, listed with its siblings
            lineage[rank] = name
        elif names:
            children[rank] = names
    return lineage, children, baltimore


def _species(box: _Box | None) -> list[dict]:
    if box is None:
        return []
    out: list[dict] = []
    seen: set[str] = set()
    for em in box.node.find_all(["em", "i"]):
        name = " ".join(em.get_text(" ").split())
        if not name or name in seen:
            continue
        seen.add(name)
        out.append({"name": name, "highlighted": em.find_parent("strong") is not None})
    if not out:
        for line in _box_text(box).split("\n"):
            name = line.lstrip("- ").strip()
            if name and name not in seen:
                seen.add(name)
                out.append({"name": name, "highlighted": False})
    return out


_NUCCORE_RE = re.compile(r"\b([A-Z]{1,6}_?\d{5,}(?:\.\d+)?)\b")


def _is_strain_box(key: str) -> bool:
    """REFERENCE STRAIN(S) boxes, including prefixed ones ("HIV-1 REFERENCE
    STRAIN", "GROUP M REFERENCE STRAINS")."""
    return key.endswith(("reference_strain", "reference_strains"))


def _segments(node: Tag) -> list[tuple[str, list[Tag]]]:
    """A box's content as ``(text, links)`` lines, split at ``<br>`` and block
    elements (strain names and their links sit on separate lines)."""
    out: list[tuple[str, list[Tag]]] = []
    text: list[str] = []
    links: list[Tag] = []

    def flush() -> None:
        line = " ".join("".join(text).split())
        if line or links:
            out.append((line, list(links)))
        text.clear()
        links.clear()

    def walk(n) -> None:
        if isinstance(n, NavigableString):
            if type(n) is NavigableString:
                text.append(str(n))
            return
        if not isinstance(n, Tag) or n.name in {"legend", "script", "style", "img"}:
            return
        if n.name == "br":
            flush()
            return
        if n.name == "a" and n.get("href"):
            links.append(n)
            text.append(n.get_text(" "))
            return
        block = n.name in _BLOCKS
        if block:
            flush()
        for child in n.children:
            walk(child)
        if block:
            flush()

    for child in node.children:
        walk(child)
    flush()
    return out


def _sequence_accession(href: str) -> str | None:
    """The nucleotide accession a "Sequence" link points to, if any."""
    parsed = urlparse(href)
    q = parse_qs(parsed.query)
    if "nucleotide" in href or "/nuccore/" in href:
        acc = (q.get("val") or q.get("term") or [""])[0] or unquote(parsed.path)
    elif "expasyfetch" in href:
        acc = unquote(parsed.query)
    else:
        return None
    m = _NUCCORE_RE.search(acc)
    return m.group(1) if m else None


def _reference_strains(boxes: list[_Box]) -> list[dict]:
    """Strains listed in the REFERENCE STRAIN(S) boxes: a name line followed by
    sequence / assembly / proteome links. Lines with links but no name before
    them (box-level notes) are dropped."""
    strains: list[dict] = []
    seen: set[tuple] = set()
    for box in boxes:
        current: dict | None = None
        for line, links in _segments(box.node):
            if not links:
                current = {"name": line, "sequences": [], "assemblies": [], "proteome_url": None}
                strains.append(current)
                continue
            if current is None:
                continue
            for a in links:
                href = a["href"]
                label = " ".join(a.get_text(" ").split())
                acc = _sequence_accession(href)
                if acc:
                    segment = label if label.casefold() != "sequence" else None
                    current["sequences"].append({"accession": acc, "segment": segment})
                elif "/datasets/genome/" in href:
                    acc = unquote(href.rstrip("/").rsplit("/", 1)[-1])
                    if acc not in current["assemblies"]:
                        current["assemblies"].append(acc)
                elif "uniprot.org" in href and current["proteome_url"] is None:
                    current["proteome_url"] = href
    out: list[dict] = []
    for strain in strains:
        if not (strain["sequences"] or strain["assemblies"]):
            continue
        key = (strain["name"], tuple(x["accession"] for x in strain["sequences"]))
        if key not in seen:
            seen.add(key)
            out.append(strain)
    return out


def _db_links(box: _Box | None) -> dict:
    links: dict = {}
    other: list[dict] = []
    if box is None:
        return links
    for a in box.node.find_all("a", href=True):
        href = a["href"]
        parsed = urlparse(href)
        q = parse_qs(parsed.query)
        host = parsed.netloc
        if "Taxonomy/Browser" in href:
            links["ncbi_taxonomy"] = href
        elif "ncbi.nlm.nih.gov/labs/virus" in href:
            links["ncbi_virus"] = href
        elif "uniprot.org" in host:
            links["uniprot"] = href
        elif "ictv.global" in host and "taxondetails" in href:
            links["ictv_taxonomy"] = href
            if q.get("taxnode_id"):
                links["ictv_taxnode_id"] = q["taxnode_id"][0]
        elif "ictv.global" in host and "/report/" in href:
            links["ictv_report"] = href
        elif "bv-brc.org" in host:
            links["bv_brc"] = href
        elif "bfvd" in host:
            links["bfvd"] = href
        else:
            label = " ".join(a.get_text(" ").split()) or (a.get("title") or host)
            other.append({"label": label, "url": href})
    if other:
        links["other"] = other
    return links


def _name_and_taxid(soup: BeautifulSoup) -> tuple[str | None, int | None]:
    h1 = soup.find("h1", attrs={"data-about": True})
    if h1 is None:
        return None, None
    m = _TAXID_RE.search(h1.get("data-about", ""))
    if not m:
        return None, None
    name = None
    for child in h1.children:
        if isinstance(child, NavigableString) and str(child).strip():
            name = " ".join(str(child).split())
            break
    if not name:
        content = h1.get("content", "")
        name = content.removesuffix(" ViralZone page").strip() or None
    return name, int(m.group(1))


def _rank(name: str, lineage: dict) -> str:
    for rank, value in lineage.items():
        if value == name:
            return rank
    return _rank_from_name(name) or "virus"


def parse_factsheet(html: str, page: int) -> dict | None:
    """A factsheet's structured content, or ``None`` if the page is not one
    (no taxid on its title)."""
    soup = _soup(html)
    name, taxid = _name_and_taxid(soup)
    if name is None or taxid is None:
        return None
    tabs = {div.get("id"): div for div in soup.select("#wiki div[id^=t]")}
    main_tab = tabs.get("t0") or soup.find(id="wiki") or soup
    boxes = _boxes(main_tab)
    by_key: dict[str, list[_Box]] = {}
    for box in boxes:
        by_key.setdefault(_key(box.title), []).append(box)

    def first(*keys: str) -> _Box | None:
        for k in keys:
            if by_key.get(k):
                return by_key[k][0]
        return None

    lineage, children, baltimore = _lineage(first("taxonomy"), name)
    record: dict = {
        "id": page,
        "url": f"{BASE_URL}/{page}",
        "name": name,
        "taxid": taxid,
        "rank": _rank(name, lineage),
        "baltimore": baltimore,
        "lineage": lineage,
        "children": children,
        "species": _species(first("virus", "species")),
    }

    # Main-column sections (VIRION, GENOME, ...), main tab only.
    main_col = main_tab.find("div", class_="main") or main_tab
    summary: dict[str, str] = {}
    extra_sections: list[dict] = []
    images: dict[str, str] = {}
    for sid, title, nodes in _sections(main_col):
        key = _SECTION_ALIASES.get(_key(sid)) or _SECTION_ALIASES.get(_key(title))
        text = _nodes_text(nodes)
        if key in {"virion", "genome"}:
            imgs = _images(nodes)
            if imgs:
                images[key] = imgs[0]
        if not text:
            continue
        if key:
            summary[key] = text
        else:
            extra_sections.append({"title": title, "text": text})
    record["summary"] = summary
    record["replication_site"] = None
    record["other_sections"] = extra_sections

    # Sidebar boxes: host, ecology & disease, therapy.
    fields: dict[str, str] = {}
    other_fields: list[dict] = []
    all_boxes = boxes + [b for tid, tab in tabs.items() if tid != "t0" for b in _boxes(tab)]
    for box in all_boxes:
        k = _key(box.title)
        if k in _SKIP_FIELDS or _is_strain_box(k):
            continue
        text = _box_text(box)
        if not text or (text.casefold() in _NONE_VALUES
                        and k not in {"vaccine", "vaccines", "antiviral_drugs"}):
            continue
        alias = _FIELD_ALIASES.get(k)
        if alias:
            fields.setdefault(alias, text)
        elif box.node.find_parent(id="t0") is not None or tabs.get("t0") is None:
            if text and not any(o["title"] == box.title for o in other_fields):
                other_fields.append({"title": box.title, "text": text})
    for f in _TEXT_FIELDS:
        record[f] = fields.get(f)
    record["other_fields"] = other_fields

    # Host-virus interaction tab: one entry per <h4> mechanism.
    interactions: list[dict] = []
    for tid, tab in tabs.items():
        if (tab.get("title") or "").casefold() != "host-virus interaction":
            continue
        col = tab.find("div", class_="main") or tab
        for _sid, title, nodes in _sections(col):
            text = _nodes_text(nodes)
            if title or text:
                interactions.append({"title": title, "text": text})
    record["host_virus_interactions"] = interactions

    record["reference_strains"] = _reference_strains(
        [b for k, bs in by_key.items() if _is_strain_box(k) for b in bs]
    )
    record["links"] = _db_links(first("db_links"))
    record["images"] = images
    return normalize(record)


def normalize(record: dict) -> dict:
    """Smooth over spelling variants between factsheets, in place.

    Works on an already-parsed record (so it can be re-applied to a stored
    data set): rank labels, the replication site, alternative box titles and
    navigation-only boxes, raw Textile links.
    """
    for key in _TEXT_FIELDS:
        if record.get(key):
            record[key] = _TEXTILE_LINK_RE.sub(r"\1", record[key])
    record["summary"] = {k: _TEXTILE_LINK_RE.sub(r"\1", v) for k, v in record["summary"].items()}
    lineage = {}
    for rank, value in record["lineage"].items():
        k = _key(rank)
        lineage[_RANK_LABELS.get(k, k)] = value
    record["lineage"] = lineage
    record["children"] = {_RANK_LABELS.get(_key(r), _key(r)): v for r, v in record["children"].items()}
    record["rank"] = _rank(record["name"], lineage)

    replication = record["summary"].get("replication")
    if replication and not record.get("replication_site"):
        lines = replication.split("\n")
        for i, line in enumerate(lines[:3]):
            upper = line.strip().upper().replace("NUCLEO-CYTOPLASMIC", "NUCLEAR/CYTOPLASMIC")
            if _SITE_LINE_RE.match(upper):
                words = re.findall(_SITE_WORD, upper)
                sites = sorted({_SITE_NAMES.get(w, w) for w in words})
                record["replication_site"] = " AND ".join(sites)
                rest = "\n".join(lines[:i] + lines[i + 1:]).strip()
                if rest:
                    record["summary"]["replication"] = rest
                else:
                    del record["summary"]["replication"]
                break

    kept = []
    for box in record["other_fields"]:
        alias = _FIELD_ALIASES_EXTRA.get(_key(box["title"]))
        if alias:
            if not record.get(alias):
                record[alias] = box["text"]
        elif not box["text"].startswith("Back to "):
            kept.append(box)
    record["other_fields"] = kept
    return record


# -- build / install --------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass
class InstalledViralZone:
    path: Path
    factsheets: int
    taxa: int
    failed: list[dict] = field(default_factory=list)


def build(client: ViralZoneClient, *, progress: Progress | None = None) -> dict:
    """Crawl ViralZone and return the data set (see :mod:`viralfetch.viralzone`)."""
    home = client.get("/")
    if home is None:
        raise ViralZoneInstallError("the ViralZone home page is unavailable")
    realm_pages = parse_taxonomy_index(home)
    if not realm_pages:
        raise ViralZoneInstallError(
            "no taxonomy diagram on the ViralZone home page — the site layout may have changed"
        )

    taxa: list[dict] = []
    candidates: list[int] = []
    seen_ids: set[int] = set()
    realm_ids = {p["id"] for p in realm_pages}
    for rp in realm_pages:
        html = client.get(f"/{rp['id']}")
        if html is None:
            continue
        taxa.extend(parse_realm_page(html, rp["id"]))
        for pid in factsheet_candidates(html):
            if pid not in seen_ids and pid not in realm_ids:
                seen_ids.add(pid)
                candidates.append(pid)

    factsheets: dict[str, dict] = {}
    failed: list[dict] = []
    for i, pid in enumerate(candidates, 1):
        entry: dict = {"id": pid}
        try:
            html = client.get(f"/{pid}")
            record = parse_factsheet(html, pid) if html else None
            if record is not None:
                factsheets[str(pid)] = record
                entry.update(status="factsheet", name=record["name"], rank=record["rank"])
            else:
                entry.update(status="skipped")
        except Exception as exc:  # one bad page must not abort the crawl
            entry.update(status="error", error=f"{type(exc).__name__}: {exc}")
            failed.append(entry)
        if progress is not None:
            progress(i, len(candidates), entry)

    if len(factsheets) < MIN_FACTSHEETS:
        raise ViralZoneInstallError(
            f"only {len(factsheets)} factsheets parsed — the ViralZone layout may have changed"
        )
    return {
        **vz_mod.CREDIT,
        "tool_version": TOOL_VERSION,
        "generated_at": _now(),
        "realm_pages": realm_pages,
        "counts": {"factsheets": len(factsheets), "taxa": len(taxa), "failed": len(failed)},
        "failed": failed,
        "taxa": taxa,
        "factsheets": dict(sorted(factsheets.items(), key=lambda kv: int(kv[0]))),
    }


def write(data: dict, dest: Path) -> None:
    """Write the data set atomically, after checking it loads."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".{dest.name}.partial")
    try:
        tmp.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        loaded = vz_mod.load_from(tmp)
        if not loaded.factsheets:
            raise ViralZoneInstallError("the rebuilt data set has no factsheets")
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)
    vz_mod.load.cache_clear()


def install(
    client: ViralZoneClient,
    *,
    dest: Path | None = None,
    progress: Progress | None = None,
) -> InstalledViralZone:
    """Rebuild the data set and install it (user data dir unless ``dest``)."""
    data = build(client, progress=progress)
    target = dest or vz_mod.installed_path()
    write(data, target)
    return InstalledViralZone(
        path=target, factsheets=data["counts"]["factsheets"],
        taxa=data["counts"]["taxa"], failed=data["failed"],
    )


def reset() -> bool:
    """Remove the installed data set, reverting to the bundled one."""
    path = vz_mod.installed_path()
    if not path.exists():
        return False
    path.unlink()
    vz_mod.load.cache_clear()
    return True



if __name__ == "__main__":  # maintainer entry point: refresh the bundled copy
    import sys

    from . import config as config_mod

    _dest = Path(sys.argv[1]) if len(sys.argv) > 1 else vz_mod.bundled_path()
    _result = install(
        ViralZoneClient(config_mod.resolve()),
        dest=_dest,
        progress=lambda i, n, e: print(f"[{i}/{n}] {e}", file=sys.stderr),
    )
    print(f"wrote {_result.path}: {_result.factsheets} factsheets, "
          f"{_result.taxa} taxa, {len(_result.failed)} failed", file=sys.stderr)
