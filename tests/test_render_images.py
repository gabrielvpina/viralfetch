"""Terminal-graphics rendering of ICTV chapter figures."""

import io

import pytest

from viralfetch import config as config_mod
from viralfetch.cli import _fetch_figures, _want_figures
from viralfetch.models import Chapter, ChapterImage
from viralfetch.render import rich_


# -- CLI gating: when should figures be drawn? ----------------------------

def _cfg(fmt):
    return config_mod.Config(email="t@example.com", format=fmt)


@pytest.mark.parametrize(
    "images,raw,fmt,tty,expected",
    [
        (True, False, "rich", True, True),    # the one case that draws
        (False, False, "rich", True, False),  # flag off
        (True, True, "rich", True, False),    # --raw stays plain
        (True, False, "json", True, False),   # JSON never draws
        (True, False, "rich", False, False),  # piped: no graphics
    ],
)
def test_want_figures_matrix(monkeypatch, images, raw, fmt, tty, expected):
    monkeypatch.setattr("sys.stdout.isatty", lambda: tty)
    assert _want_figures(images, raw, _cfg(fmt)) is expected


# -- CLI fetch: unavailable figures are skipped, not fatal ----------------

class _FakeClient:
    def __init__(self, table):
        self.table = table  # url -> bytes | None

    def fetch_image(self, url):
        return self.table.get(url)


def test_fetch_figures_skips_unavailable():
    chapter = Chapter(
        slug="x", title="X", markdown="",
        images=[
            ChapterImage(url="https://ictv.global/a.png", alt="A"),
            ChapterImage(url="https://ictv.global/b.png", alt="B"),  # returns None
            ChapterImage(url="https://ictv.global/c.png", alt="C"),
        ],
    )
    client = _FakeClient({
        "https://ictv.global/a.png": b"aa",
        "https://ictv.global/b.png": None,
        "https://ictv.global/c.png": b"cc",
    })
    figures = _fetch_figures(client, chapter)
    # A dict keyed by URL (for inline placement); the unavailable B is dropped.
    assert figures == {
        "https://ictv.global/a.png": b"aa",
        "https://ictv.global/c.png": b"cc",
    }


# -- block-element conversion (Pillow, a core dependency) -----------------

def _png(w, h):
    """A gradient PNG: varied enough that cells need real block glyphs."""
    Image = pytest.importorskip("PIL.Image")
    img = Image.new("RGB", (w, h))
    px = img.load()
    for y in range(h):
        for x in range(w):
            px[x, y] = ((x * 255) // max(1, w - 1), (y * 255) // max(1, h - 1), 128)
    buf = io.BytesIO()
    img.save(buf, "PNG")
    return buf.getvalue()


def _has_block(line):
    return any(ch in rich_._BLOCK_GLYPHS for ch in line)


def test_blockart_preserves_aspect_and_width():
    data = _png(40, 20)  # 2:1 image
    art = rich_._blockart(data, cols=20)
    assert art is not None
    lines = art.plain.split("\n")
    assert all(len(line) == 20 for line in lines)   # 20 cells wide (40 sub-px)
    # 20 cells over a 2:1 image => 10 sub-rows => 5 cell rows.
    assert len(lines) == 5
    assert set(art.plain) <= set(rich_._GLYPHS + "\n")  # only block-element glyphs


def test_blockart_doubles_horizontal_density():
    # Two sub-pixels per cell: cols is capped at half the source width.
    data = _png(8, 8)
    art = rich_._blockart(data, cols=200)  # asked wider than the image
    assert max(len(l) for l in art.plain.split("\n")) == 4  # src_w // 2


def test_blockart_bad_bytes_returns_none():
    pytest.importorskip("PIL.Image")
    assert rich_._blockart(b"not-an-image", cols=20) is None


def test_clamp_keeps_dithered_colours_in_gamut():
    # Error diffusion can push a channel past 0..255; _clamp must reel it in.
    assert rich_._clamp(-40) == 0
    assert rich_._clamp(321) == 255
    assert rich_._clamp(127.6) == 128
    # _hex always yields a valid 6-digit triple, even for out-of-range input.
    assert rich_._hex((-5, 300, 127.4)) == "00ff7f"


def _chapter_with_images(n=1):
    imgs = [ChapterImage(url=f"https://ictv.global/{i}.png", alt=f"fig{i}") for i in range(n)]
    return Chapter(slug="x", title="X", markdown="", images=imgs)


def _render(renderable, width=100):
    from rich.console import Console
    buf = io.StringIO()
    Console(file=buf, width=width, force_terminal=False).print(renderable)
    return buf.getvalue()


# -- inline figure placement ----------------------------------------------

_MD = (
    "intro paragraph\n\n"
    "![virion morphology](https://ictv.global/a.png)\n\n"
    "Figure 1. the caption\n\n"
    "closing paragraph\n"
)


def test_with_figures_none_renders_plain_markdown():
    # Feature off: the image reference stays as Markdown text, no graphics.
    out = _render(rich_._with_figures(_chapter_with_images(), _MD, None, None))
    assert not _has_block(out)
    assert "intro paragraph" in out and "closing paragraph" in out


def test_with_figures_draws_image_in_place():
    pytest.importorskip("PIL.Image")
    figs = {"https://ictv.global/a.png": _png(20, 10)}
    out = _render(rich_._with_figures(_chapter_with_images(), _MD, figs, None))
    lines = out.splitlines()
    img_rows = [i for i, ln in enumerate(lines) if _has_block(ln)]
    intro = next(i for i, ln in enumerate(lines) if "intro paragraph" in ln)
    closing = next(i for i, ln in enumerate(lines) if "closing paragraph" in ln)
    assert img_rows, "figure was not drawn"
    # The graphic sits between the intro and the closing text — in its place.
    assert intro < min(img_rows) <= max(img_rows) < closing
    assert "virion morphology" in out  # alt shown as a caption


def test_with_figures_unknown_url_falls_back_to_reference():
    pytest.importorskip("PIL.Image")
    out = _render(rich_._with_figures(_chapter_with_images(), _MD, {}, None))
    assert not _has_block(out)             # no bytes for the URL -> no graphic
    assert "virion morphology" in out      # the Markdown reference survives


def test_fig_width_caps_figure_columns():
    pytest.importorskip("PIL.Image")
    md = "![x](https://ictv.global/a.png)\n"
    figs = {"https://ictv.global/a.png": _png(200, 100)}
    out = _render(rich_._with_figures(_chapter_with_images(), md, figs, 30), width=120)
    block_lines = [ln for ln in out.splitlines() if _has_block(ln)]
    assert block_lines and max(len(ln) for ln in block_lines) == 30


def _hide_pillow(monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "PIL" or name.startswith("PIL."):
            raise ImportError("no PIL")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)


def test_without_pillow_shows_in_band_hint_and_text(monkeypatch):
    _hide_pillow(monkeypatch)
    # Text still renders; a visible in-band note explains the missing Pillow.
    out = _render(rich_._with_figures(_chapter_with_images(1), _MD, {}, None))
    assert "intro paragraph" in out
    assert "Pillow" in out
    assert not _has_block(out)


def test_without_pillow_no_hint_when_chapter_has_no_figures(monkeypatch):
    _hide_pillow(monkeypatch)
    md = "just text\n"
    out = _render(rich_._with_figures(Chapter(slug="x", title="X", markdown=""), md, {}, None))
    assert "Pillow" not in out
    assert "just text" in out


def test_figures_supported_reflects_pillow(monkeypatch):
    _hide_pillow(monkeypatch)
    assert rich_.figures_supported() is False
