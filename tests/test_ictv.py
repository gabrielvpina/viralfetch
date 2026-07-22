"""Tests for the ICTV Report chapter client and HTML->Markdown parser.

No real network: five frozen chapter HTMLs under ``fixtures/chapters/`` drive
the parser (a regression guard — if the ICTV theme changes, these break
loudly), and a fake session serves them for the client transport tests.
"""

from pathlib import Path

import pytest

from viralfetch.config import Config
from viralfetch.cache import Cache
from viralfetch import ictv
from viralfetch.ictv import (
    ChapterNotFound,
    ChapterParseError,
    ICTVClient,
    RobotsDisallowed,
    SectionNotFound,
    _compare_vmr,
    _slug,
    _vmr_key,
    parse_chapter,
    section_markdown,
)

CHAPTERS = Path(__file__).parent / "fixtures" / "chapters"
SLUGS = ["coronaviridae", "geminiviridae", "rhabdoviridae", "filoviridae", "poxviridae"]


def _load(slug):
    html = (CHAPTERS / f"{slug}.html").read_text(encoding="utf-8")
    return parse_chapter(html, slug=slug, url=f"https://ictv.global/report/chapter/{slug}/{slug}")


def noop_sleep(_seconds):
    pass


# -- parser over frozen fixtures ------------------------------------------

@pytest.mark.parametrize("slug", SLUGS)
def test_parses_landmarks(slug):
    ch = _load(slug)
    assert ch.title.lower().startswith("family:")
    assert slug in ch.title.lower()
    md = ch.markdown
    # Original URL + attribution are kept at the very top (before the rule).
    head = md.split("\n---\n", 1)[0]
    assert f"*Source: https://ictv.global/report/chapter/{slug}/{slug}*" in head
    assert "CC BY 4.0" in head
    # Section headings survive.
    assert "\n## Summary" in md


@pytest.mark.parametrize("slug", SLUGS)
def test_preserves_italic_scientific_names(slug):
    # The family name itself appears italicised somewhere in the body.
    md = _load(slug).markdown
    assert f"*{slug.capitalize()}*" in md


def test_characteristics_table_becomes_markdown():
    md = _load("coronaviridae").markdown
    assert "| --- | --- |" in md  # a Markdown table separator row
    assert "Characteristic" in md and "Description" in md
    # Italics are preserved inside table cells.
    assert "*Betacoronavirus muris*" in md


@pytest.mark.parametrize("slug", SLUGS)
def test_figures_captured_as_images(slug):
    ch = _load(slug)
    # Every chapter carries figures; they are collected structurally...
    assert ch.images, f"no figures captured for {slug}"
    for img in ch.images:
        assert img.url.startswith("https://ictv.global/")  # absolute, resolved
    # ...and rendered inline as Markdown pointing at the same absolute URLs.
    for img in ch.images:
        assert f"]({img.url})" in ch.markdown


def test_figure_urls_are_unique_and_exclude_boilerplate_logos():
    ch = _load("coronaviridae")
    urls = [img.url for img in ch.images]
    assert len(urls) == len(set(urls))  # de-duplicated
    # The ICTV header/logo images live outside the content container.
    assert not any("ictvLogo" in u or "ICTV%20Report%20Header" in u for u in urls)


@pytest.mark.parametrize("slug", SLUGS)
def test_figures_are_standalone_blocks_not_table_rows(slug):
    # ICTV wraps figures in layout tables; the parser unwraps them so each
    # image sits on its own line (as a block) and can be drawn in place.
    import re
    md = _load(slug).markdown
    assert "| ![" not in md  # no image left inside a Markdown table cell
    for img in _load(slug).images:
        # The reference appears as its own block line, surrounded by blank lines.
        assert re.search(rf"(?m)^!\[[^\]]*\]\({re.escape(img.url)}\)\s*$", md)


def test_doi_extracted_when_present_in_preamble():
    assert _load("poxviridae").doi == "10.1099/jgv.0.001849"


def test_doi_absent_is_none_not_fabricated():
    # Corona's DOI lives in the reference list, not the top block: stays None.
    assert _load("coronaviridae").doi is None


# -- loud failure on layout change ----------------------------------------

def test_missing_container_raises():
    with pytest.raises(ChapterParseError):
        parse_chapter("<html><body><p>hi</p></body></html>", slug="x", url="u")


def test_missing_summary_raises():
    html = (
        '<div class="field--name-field-mt-srv-body">'
        "<h2>Family: Nowhere</h2><p>authors</p></div>"
    )
    with pytest.raises(ChapterParseError):
        parse_chapter(html, slug="nowhere", url="u")


# -- section selection ----------------------------------------------------

def test_section_markdown_keeps_only_matching_section():
    ch = _load("coronaviridae")
    md = section_markdown(ch, "summary")
    assert "## Summary" in md
    assert "## Virion" not in md
    # Header (source + references) is still present.
    assert "*Source:" in md


def test_section_markdown_unknown_raises():
    ch = _load("coronaviridae")
    with pytest.raises(SectionNotFound) as exc:
        section_markdown(ch, "nonexistent-section")
    assert exc.value.available  # lists the real sections


# -- slug -----------------------------------------------------------------

def test_slug():
    assert _slug("Coronaviridae") == "coronaviridae"
    assert _slug("  Geminiviridae ") == "geminiviridae"


# -- client transport (fake session) --------------------------------------

class FakeResponse:
    def __init__(self, text="", status_code=200, content=b""):
        self.text = text
        self.status_code = status_code
        self.content = content


class FakeSession:
    def __init__(self, pages, robots="User-agent: *\nDisallow: /admin/\n", fail_times=0, images=None):
        self.pages = pages          # url -> html
        self.images = images or {}  # url -> bytes
        self.robots = robots
        self.fail_times = fail_times
        self.calls = []

    def get(self, url, headers=None, timeout=None):
        self.calls.append(url)
        if url.endswith("/robots.txt"):
            return FakeResponse(self.robots)
        if self.fail_times > 0:
            self.fail_times -= 1
            return FakeResponse("", 503)
        if url in self.pages:
            return FakeResponse(self.pages[url])
        if url in self.images:
            return FakeResponse(content=self.images[url])
        return FakeResponse("", 404)


def _client(session, tmp_path, cache_enabled=True):
    cfg = Config(email="tester@example.com")
    cache = Cache(tmp_path, enabled=cache_enabled)
    return ICTVClient(cfg, cache=cache, session=session, sleep=noop_sleep)


def _url(slug):
    return f"https://ictv.global/report/chapter/{slug}/{slug}"


def test_missing_email_raises(tmp_path):
    with pytest.raises(Exception):
        ICTVClient(Config(email=None), session=FakeSession({}))


def test_user_agent_carries_email_and_repo(tmp_path):
    client = _client(FakeSession({}), tmp_path)
    assert client.user_agent.startswith("viralfetch/")
    assert "tester@example.com" in client.user_agent
    assert "github.com/gabrielvpina/viralfetch" in client.user_agent


def test_fetch_chapter_parses(tmp_path):
    html = (CHAPTERS / "coronaviridae.html").read_text(encoding="utf-8")
    session = FakeSession({_url("coronaviridae"): html})
    client = _client(session, tmp_path)
    ch = client.fetch_chapter("Coronaviridae")
    assert ch.title == "Family: Coronaviridae"
    assert "## Summary" in ch.markdown


def test_fetch_chapter_uses_cache(tmp_path):
    html = (CHAPTERS / "coronaviridae.html").read_text(encoding="utf-8")
    session = FakeSession({_url("coronaviridae"): html})
    client = _client(session, tmp_path)
    client.fetch_chapter("Coronaviridae")
    n = len(session.calls)
    client.fetch_chapter("Coronaviridae")  # served from the 30-day TTL cache
    assert len(session.calls) == n  # no new request


def test_fetch_chapter_404_raises(tmp_path):
    session = FakeSession({})  # nothing registered -> 404
    client = _client(session, tmp_path)
    with pytest.raises(ChapterNotFound):
        client.fetch_chapter("Nosuchviridae")


def test_robots_disallow_blocks_fetch(tmp_path):
    html = (CHAPTERS / "coronaviridae.html").read_text(encoding="utf-8")
    session = FakeSession({_url("coronaviridae"): html},
                          robots="User-agent: *\nDisallow: /report/\n")
    client = _client(session, tmp_path)
    with pytest.raises(RobotsDisallowed):
        client.fetch_chapter("Coronaviridae")


def test_retry_on_5xx_then_success(tmp_path):
    html = (CHAPTERS / "coronaviridae.html").read_text(encoding="utf-8")
    session = FakeSession({_url("coronaviridae"): html}, fail_times=2)
    client = _client(session, tmp_path)
    ch = client.fetch_chapter("Coronaviridae")
    assert ch.title == "Family: Coronaviridae"


# -- image fetch ----------------------------------------------------------

IMG_URL = "https://ictv.global/system/files/inline-images/OPSR.Corona.Fig1_.v1.png"


def test_fetch_image_returns_bytes_and_caches(tmp_path):
    session = FakeSession({}, images={IMG_URL: b"\x89PNG\r\n\x1a\n-figure-bytes"})
    client = _client(session, tmp_path)
    data = client.fetch_image(IMG_URL)
    assert data == b"\x89PNG\r\n\x1a\n-figure-bytes"
    n = len(session.calls)
    again = client.fetch_image(IMG_URL)  # served from the permanent cache
    assert again == data
    assert len(session.calls) == n  # no new request


def test_fetch_image_offsite_src_refused(tmp_path):
    off = "https://evil.example.com/tracker.png"
    session = FakeSession({}, images={off: b"nope"})
    client = _client(session, tmp_path)
    assert client.fetch_image(off) is None
    assert off not in session.calls  # never requested


def test_fetch_image_missing_returns_none_not_raise(tmp_path):
    session = FakeSession({})  # unknown url -> 404
    client = _client(session, tmp_path)
    assert client.fetch_image(IMG_URL) is None


def test_fetch_image_robots_disallow_skips_image(tmp_path):
    # robots disallowing the image path is respected but non-fatal: the figure
    # is skipped (None), never fetched, and the chapter text is unaffected.
    session = FakeSession({}, images={IMG_URL: b"data"},
                          robots="User-agent: *\nDisallow: /system/\n")
    client = _client(session, tmp_path)
    assert client.fetch_image(IMG_URL) is None
    assert IMG_URL not in session.calls  # never requested


# -- VMR update check -----------------------------------------------------

def test_vmr_key_orders_releases():
    assert _vmr_key("VMR_MSL41.v1.20260320.tsv") > _vmr_key("VMR_MSL40.v2.20260223.xlsx")
    assert _vmr_key("VMR_MSL40.v2.20260223") > _vmr_key("VMR_MSL40.v1.20250307")
    assert _vmr_key("VMR_MSL38_v1") > _vmr_key("VMR_MSL37")  # missing date -> 0


_VMR_PAGE = """
<html><body>
  <a href="/sites/default/files/VMR/VMR_MSL40.v2.20260223.xlsx">old</a>
  <a href="/sites/default/files/VMR/VMR_MSL41.v1.20260320.xlsx">latest</a>
</body></html>
"""


def test_compare_vmr_up_to_date():
    result = _compare_vmr(_VMR_PAGE, "VMR_MSL41.v1.20260320.tsv")
    assert result.up_to_date is True
    assert result.latest == "VMR_MSL41.v1.20260320.xlsx"
    assert result.latest_url.startswith("https://ictv.global/")


def test_compare_vmr_newer_available():
    result = _compare_vmr(_VMR_PAGE, "VMR_MSL40.v2.20260223.tsv")
    assert result.up_to_date is False
    assert result.latest == "VMR_MSL41.v1.20260320.xlsx"


def test_compare_vmr_no_downloads_raises():
    with pytest.raises(ictv.ICTVError):
        _compare_vmr("<html><body>nothing</body></html>", "VMR_MSL41.v1.20260320.tsv")


def test_check_vmr_update_via_client(tmp_path):
    session = FakeSession({"https://ictv.global/vmr": _VMR_PAGE})
    client = _client(session, tmp_path)
    result = client.check_vmr_update("VMR_MSL41.v1.20260320.tsv")
    assert result.up_to_date is True
