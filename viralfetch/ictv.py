"""ICTV Report chapter client: fetch, parse, and convert to Markdown.

The ICTV Report (``ictv.global/report``) publishes a descriptive chapter per
family. This module fetches a chapter page, scopes to its main content
container, and turns that HTML into Markdown — preserving section headings,
the italics of scientific names (mandatory in viral nomenclature),
characteristic tables, and figure references (as Markdown ``![alt](url)`` with
absolute URLs; the images themselves are not downloaded here).

Politeness (SPEC section 3): a descriptive ``User-Agent`` carrying the contact
email, a minimum 1-second gap between requests, ``robots.txt`` honoured, and
chapter HTML cached with a 30-day TTL. Content is CC BY 4.0, so the original
page URL and the chapter's references/attribution are kept at the top of the
output.

Parsing fails **loudly** when an expected structural landmark is missing (the
content container, the title heading, or the ``Summary`` section). A silent
``None`` would mask a breaking change in the ICTV site theme (SPEC section 5.7).
"""

from __future__ import annotations

import re
import time
import urllib.robotparser
from dataclasses import dataclass
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup

from . import __version__
from .cache import IMAGES, TEXT_TTL, TEXTS, Cache
from .config import Config
from .models import Chapter, ChapterImage
from .ncbi import RateLimiter

BASE_URL = "https://ictv.global"
REPO_URL = "https://github.com/gabrielvpina/viralfetch"
CHAPTER_PATH = "/report/chapter/{slug}/{slug}"
CONTENT_SELECTOR = ".field--name-field-mt-srv-body"
LICENSE_NOTE = "Content is licensed CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)."

_RETRY_STATUS = {429, 500, 502, 503, 504}
_DOI_RE = re.compile(r"10\.\d{4,}/[^\s\"'<)\]]+")

# HTML tags we treat as block-level, as transparent containers to descend into,
# and as content to drop entirely (SPEC 5.7: never scrape unscoped). Images are
# captured as Markdown ``![alt](url)`` with an absolute URL; only genuinely
# non-content chrome (scripts, menus, forms) is dropped.
_BLOCK_TAGS = {"h1", "h2", "h3", "h4", "h5", "h6", "p", "ul", "ol", "table", "blockquote", "figure"}
_TRANSPARENT = {"div", "section", "article", "picture"}
_SKIP = {"script", "style", "nav", "form", "button", "iframe"}


# -- errors ---------------------------------------------------------------

class ICTVError(Exception):
    """Base class for ICTV fetch/parse failures."""


class ChapterNotFound(ICTVError):
    def __init__(self, name: str, url: str):
        self.name = name
        self.url = url
        super().__init__(f"no ICTV Report chapter for {name!r} at {url}")


class ChapterParseError(ICTVError):
    """An expected structural landmark was missing — the site layout changed."""


class SectionNotFound(ICTVError):
    def __init__(self, section: str, available: list[str]):
        self.section = section
        self.available = available
        super().__init__(f"section {section!r} not found")


class RobotsDisallowed(ICTVError):
    """robots.txt forbids fetching the requested path."""


# -- client ---------------------------------------------------------------

def _slug(name: str) -> str:
    """Derive the URL slug from a taxon name (families are single words)."""
    return re.sub(r"\s+", "-", name.strip().lower())


class ICTVClient:
    """Polite HTTP client for ICTV Report chapters."""

    def __init__(
        self,
        config: Config,
        cache: Cache | None = None,
        session: requests.Session | None = None,
        sleep=time.sleep,
        clock=time.monotonic,
        max_retries: int = 3,
    ):
        self.email = config.require_email()  # fails loudly if unset
        self.cache = cache
        self.session = session or requests.Session()
        self.max_retries = max_retries
        # >= 1 second between requests (SPEC section 3).
        self._limiter = RateLimiter(1, sleep=sleep, clock=clock)
        self._sleep = sleep
        self._robots: urllib.robotparser.RobotFileParser | None = None
        self.user_agent = f"viralfetch/{__version__} (+{REPO_URL}; {self.email})"

    # -- transport --------------------------------------------------------

    def _robots_parser(self) -> urllib.robotparser.RobotFileParser | None:
        if self._robots is None:
            rp = urllib.robotparser.RobotFileParser()
            try:
                body = self.session.get(
                    f"{BASE_URL}/robots.txt",
                    headers={"User-Agent": self.user_agent},
                    timeout=30,
                ).text
                rp.parse(body.splitlines())
                self._robots = rp
            except requests.RequestException:
                # If robots.txt is unreachable, do not invent rules; proceed.
                self._robots = urllib.robotparser.RobotFileParser()
                self._robots.allow_all = True
        return self._robots

    def _check_robots(self, path: str) -> None:
        rp = self._robots_parser()
        if rp is not None and not rp.can_fetch(self.user_agent, BASE_URL + path):
            raise RobotsDisallowed(f"robots.txt disallows {path}")

    def _send(self, url: str, path: str) -> requests.Response:
        """Rate-limited, robots-checked GET returning the final response.

        Retries on 5xx; any other status (including 404) is returned to the
        caller, which decides what a missing resource means.
        """
        self._check_robots(path)
        headers = {"User-Agent": self.user_agent}
        last_exc: Exception | None = None
        for attempt in range(self.max_retries):
            self._limiter.wait()
            try:
                resp = self.session.get(url, headers=headers, timeout=60)
            except requests.RequestException as exc:
                last_exc = exc
                self._backoff(attempt)
                continue
            if resp.status_code in _RETRY_STATUS:
                last_exc = ICTVError(f"HTTP {resp.status_code} from {url}")
                self._backoff(attempt)
                continue
            return resp
        raise ICTVError(f"{url} failed after {self.max_retries} attempts: {last_exc}")

    def _fetch(self, path: str) -> tuple[int, str]:
        """GET a path under the ICTV domain, returning ``(status_code, text)``."""
        resp = self._send(BASE_URL + path, path)
        return resp.status_code, resp.text

    def _get(self, path: str, name: str) -> str:
        status, text = self._fetch(path)
        if status == 404:
            raise ChapterNotFound(name, BASE_URL + path)
        if status != 200:
            raise ICTVError(f"HTTP {status} from {BASE_URL + path}: {text[:200]}")
        return text

    def _backoff(self, attempt: int) -> None:
        self._sleep(0.5 * (2 ** attempt))

    # -- public API -------------------------------------------------------

    def fetch_chapter(self, name: str) -> Chapter:
        """Fetch and parse the ICTV Report chapter for ``name``."""
        slug = _slug(name)
        path = CHAPTER_PATH.format(slug=slug)
        url = BASE_URL + path
        cache_key = f"chapter:{slug}"
        html = self.cache.get(TEXTS, cache_key, ttl=TEXT_TTL) if self.cache else None
        if html is None:
            html = self._get(path, name)
            if self.cache:
                self.cache.set(TEXTS, cache_key, html)
        return parse_chapter(html, slug=slug, url=url)

    def fetch_image(self, url: str) -> bytes | None:
        """Fetch a chapter figure's bytes, cached permanently.

        Returns ``None`` when the figure is unavailable — a broken or missing
        image must never sink a chapter, whose text is the primary content.
        Only figures hosted on the ICTV domain are fetched; an off-site ``src``
        is refused rather than followed (SPEC section 3).
        """
        cache_key = f"image:{url}"
        if self.cache:
            cached = self.cache.get_bytes(IMAGES, cache_key)  # permanent, no TTL
            if cached is not None:
                return cached
        parsed = urlparse(url)
        if parsed.netloc != urlparse(BASE_URL).netloc:
            return None
        try:
            resp = self._send(url, parsed.path)
        except ICTVError:
            return None
        if resp.status_code != 200:
            return None
        data = resp.content
        if self.cache:
            self.cache.set_bytes(IMAGES, cache_key, data)
        return data

    def check_vmr_update(self, current_filename: str) -> "VMRUpdate":
        """Compare the embedded VMR against the newest one on ictv.global/vmr."""
        status, html = self._fetch(VMR_PATH)
        if status != 200:
            raise ICTVError(f"HTTP {status} from {BASE_URL + VMR_PATH}")
        return _compare_vmr(html, current_filename)


# -- VMR version check ----------------------------------------------------

VMR_PATH = "/vmr"


@dataclass
class VMRUpdate:
    current: str
    latest: str | None
    latest_url: str | None
    up_to_date: bool


def _vmr_key(name: str) -> tuple[int, int, int]:
    """Sortable (MSL, version, date) key parsed from a VMR filename.

    Filenames vary over releases — ``VMR_MSL39.v1_20240912``,
    ``VMR_MSL41.v1.20260320``, older ``VMR_MSL38_v1`` with no date — so each
    component is extracted independently and defaults to 0 when absent.
    """
    msl = re.search(r"MSL(\d+)", name)
    ver = re.search(r"[._]v(\d+)", name, re.I)
    date = re.search(r"(\d{8})", name)
    return (
        int(msl.group(1)) if msl else 0,
        int(ver.group(1)) if ver else 0,
        int(date.group(1)) if date else 0,
    )


def _compare_vmr(html: str, current_filename: str) -> VMRUpdate:
    soup = BeautifulSoup(html, "lxml")
    candidates: dict[str, str] = {}
    for a in soup.find_all("a", href=True):
        m = re.search(r"VMR_MSL\d[^\"'/]*\.xlsx", a["href"], re.I)
        if m:
            fname = m.group(0)
            href = a["href"]
            candidates[fname] = href if href.startswith("http") else BASE_URL + href
    if not candidates:
        raise ICTVError("no VMR downloads found on the ICTV page — the layout may have changed")
    latest = max(candidates, key=_vmr_key)
    up_to_date = _vmr_key(current_filename) >= _vmr_key(latest)
    return VMRUpdate(
        current=current_filename,
        latest=latest,
        latest_url=candidates[latest],
        up_to_date=up_to_date,
    )


# -- parsing --------------------------------------------------------------

def parse_chapter(html: str, *, slug: str, url: str) -> Chapter:
    """Parse chapter HTML into a :class:`Chapter` with Markdown content.

    Raises :class:`ChapterParseError` if the content container, the title
    heading, or the ``Summary`` section is missing.
    """
    soup = BeautifulSoup(html, "lxml")
    body = soup.select_one(CONTENT_SELECTOR)
    if body is None:
        raise ChapterParseError(
            f"content container {CONTENT_SELECTOR!r} not found — the ICTV page layout may have changed"
        )

    blocks = list(_iter_blocks(body))
    title_el = next((b for b in blocks if b.name in ("h1", "h2")), None)
    if title_el is None:
        raise ChapterParseError("no chapter title heading found")
    title = _collapse_ws(title_el.get_text(" ", strip=True)).strip()

    summary_idx = next(
        (i for i, b in enumerate(blocks)
         if b.name == "h2" and b.get_text(strip=True).lower().startswith("summary")),
        None,
    )
    if summary_idx is None:
        raise ChapterParseError(
            "no 'Summary' section found — the ICTV page layout may have changed"
        )

    title_idx = blocks.index(title_el)
    preamble = blocks[title_idx + 1:summary_idx]
    content = blocks[summary_idx:]

    references_md = "\n\n".join(_block_md(b) for b in preamble if _block_md(b)).strip()
    content_md = "\n\n".join(_block_md(b) for b in content if _block_md(b)).strip()
    doi_match = _DOI_RE.search(references_md)
    doi = doi_match.group(0).rstrip(".,;") if doi_match else None

    markdown = _assemble(title, url, references_md, content_md)
    return Chapter(
        slug=slug,
        title=title,
        markdown=markdown,
        authors=None,
        citation=references_md or None,
        doi=doi,
        url=url,
        images=_collect_images(content),
    )


def _collect_images(blocks) -> list[ChapterImage]:
    """Gather chapter figures (in document order) from the content blocks.

    Duplicate URLs are dropped so a figure that a page repeats is listed once.
    """
    seen: set[str] = set()
    out: list[ChapterImage] = []
    for block in blocks:
        imgs = [block] if block.name == "img" else block.find_all("img")
        for img in imgs:
            src = img.get("src")
            if not src:
                continue
            url = _resolve_url(src)
            if url in seen:
                continue
            seen.add(url)
            out.append(ChapterImage(url=url, alt=_collapse_ws(img.get("alt") or "").strip()))
    return out


def _resolve_url(src: str) -> str:
    """Resolve an image ``src`` to an absolute URL against the ICTV domain."""
    return urljoin(BASE_URL + "/", src.strip())


def _iter_blocks(node):
    """Yield block-level elements in document order, descending transparent
    wrappers and dropping menus/scripts. Standalone images are yielded as
    blocks so a figure outside any paragraph is not lost."""
    for child in node.children:
        name = getattr(child, "name", None)
        if name is None or name in _SKIP:
            continue
        if name in _BLOCK_TAGS or name == "img":
            yield child
        elif name in _TRANSPARENT:
            yield from _iter_blocks(child)
        else:
            yield from _iter_blocks(child)


def _block_md(el) -> str:
    name = el.name
    if name in ("h1", "h2", "h3", "h4", "h5", "h6"):
        return f"{'#' * int(name[1])} {_inline(el).strip()}"
    if name in ("ul", "ol"):
        return _list_md(el, ordered=(name == "ol"))
    if name == "table":
        # ICTV wraps each figure in a one-column table (image row + caption
        # row). Unwrap those into standalone blocks so the figure sits inline
        # in document order; only genuine data tables become Markdown tables.
        if el.find("img") is not None:
            return _figure_table_md(el)
        return _table_md(el)
    if name == "blockquote":
        return "> " + _inline(el).strip()
    if name == "img":
        return _img_md(el)
    if name == "figure":
        return _figure_md(el)
    return _inline(el).strip()  # p and anything else


def _inline(node) -> str:
    """Convert inline content to Markdown, preserving italics and bold."""
    parts: list[str] = []
    for child in node.children:
        name = getattr(child, "name", None)
        if name is None:
            parts.append(str(child))
        elif name in _SKIP:
            continue
        elif name in ("em", "i"):
            parts.append(_emph(_inline(child), "*"))
        elif name in ("strong", "b"):
            parts.append(_emph(_inline(child), "**"))
        elif name == "img":
            parts.append(_img_md(child))
        elif name == "br":
            parts.append(" ")
        else:  # a, span, sup, sub, and other inline wrappers: keep text only
            parts.append(_inline(child))
    return _collapse_ws("".join(parts))


def _img_md(el) -> str:
    """Render an ``<img>`` as Markdown ``![alt](url)`` with an absolute URL."""
    src = el.get("src")
    if not src:
        return ""
    alt = _collapse_ws(el.get("alt") or "").strip()
    return f"![{alt}]({_resolve_url(src)})"


def _figure_table_md(table) -> str:
    """Unwrap a figure-wrapper table into standalone blocks (one per cell).

    Each cell becomes its own block: an image cell yields ``![alt](url)`` on its
    own line (so the renderer can draw it in place), a caption cell yields a
    paragraph. Empty cells are dropped.
    """
    blocks: list[str] = []
    for tr in table.find_all("tr"):
        for cell in tr.find_all(["th", "td"], recursive=False):
            text = _inline(cell).strip()
            if text:
                blocks.append(text)
    return "\n\n".join(blocks)


def _figure_md(el) -> str:
    """Render a ``<figure>`` as its image followed by its caption in italics."""
    parts: list[str] = []
    img = el.find("img")
    if img is not None:
        md = _img_md(img)
        if md:
            parts.append(md)
    caption = el.find("figcaption")
    if caption is not None:
        text = _inline(caption).strip()
        if text:
            parts.append(f"*{text}*")
    return "\n\n".join(parts)


def _emph(text: str, marker: str) -> str:
    """Wrap ``text`` in an emphasis ``marker`` without swallowing edge spaces
    (Markdown emphasis breaks if the marker hugs whitespace)."""
    stripped = text.strip()
    if not stripped:
        return text
    lead = text[:len(text) - len(text.lstrip())]
    trail = text[len(text.rstrip()):]
    return f"{lead}{marker}{stripped}{marker}{trail}"


def _list_md(el, ordered: bool) -> str:
    lines = []
    for i, li in enumerate(el.find_all("li", recursive=False), start=1):
        marker = f"{i}." if ordered else "-"
        text = _inline(li).strip()
        if text:
            lines.append(f"{marker} {text}")
    return "\n".join(lines)


def _table_md(table) -> str:
    grid: list[list[str]] = []
    for tr in table.find_all("tr"):
        cells = tr.find_all(["th", "td"], recursive=False)
        row = [(_inline(c).strip() or " ").replace("|", "\\|") for c in cells]
        if row:
            grid.append(row)
    if not grid:
        return ""
    ncol = max(len(r) for r in grid)
    for r in grid:
        r.extend([" "] * (ncol - len(r)))
    out = ["| " + " | ".join(grid[0]) + " |",
           "| " + " | ".join(["---"] * ncol) + " |"]
    out.extend("| " + " | ".join(r) + " |" for r in grid[1:])
    return "\n".join(out)


def _collapse_ws(text: str) -> str:
    text = text.replace("​", "").replace("\xa0", " ")
    return re.sub(r"[ \t\r\n\f\v]+", " ", text)


def _assemble(title: str, url: str, references_md: str, content_md: str) -> str:
    header = [f"# {title}", "", f"*Source: {url}*"]
    if references_md:
        header += ["", references_md]
    header += ["", f"*{LICENSE_NOTE}*", "", "---"]
    return "\n".join(header) + "\n\n" + content_md + "\n"


# -- section selection ----------------------------------------------------

_SPLIT_RE = re.compile(r"(?m)^(?=##\s)")


def section_markdown(chapter: Chapter, section: str) -> str:
    """Reduce a chapter's Markdown to the section(s) whose heading matches.

    The top header (source URL, references, licence) is always kept. Matching
    is a case-insensitive substring test on ``##``-level headings.
    """
    marker = "\n---\n"
    md = chapter.markdown
    if marker not in md:
        return md
    head, body = md.split(marker, 1)
    chunks = [c for c in _SPLIT_RE.split(body) if c.strip()]
    wanted = section.strip().lower()
    picked = [c for c in chunks if wanted in _heading_text(c).lower()]
    if not picked:
        raise SectionNotFound(section, [_heading_text(c) for c in chunks if _heading_text(c)])
    return head + marker + "\n" + "\n\n".join(c.strip() for c in picked) + "\n"


def _heading_text(chunk: str) -> str:
    first = chunk.lstrip().splitlines()[0] if chunk.strip() else ""
    return first.lstrip("#").strip() if first.startswith("##") else ""
