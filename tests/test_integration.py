"""Live NCBI integration tests. Excluded from the default run.

Run explicitly with a configured email:
    NCBI_EMAIL=you@example.com pytest -m network
"""

import os

import pytest

from viralfetch.cache import Cache
from viralfetch.config import CACHE_DIR, Config
from viralfetch.ictv import ICTVClient
from viralfetch.ncbi import NCBIClient
from viralfetch.sequences import seq_meta
from viralfetch.vmr import load

pytestmark = pytest.mark.network


@pytest.mark.skipif(not os.environ.get("NCBI_EMAIL"), reason="NCBI_EMAIL not set")
def test_seq_meta_returns_real_data():
    cfg = Config(email=os.environ["NCBI_EMAIL"])
    client = NCBIClient(cfg, cache=Cache(CACHE_DIR))
    species, result = seq_meta(load(), client, "Betacoronavirus pandemicum")
    assert species == "Betacoronavirus pandemicum"
    assert result.records  # at least one real record came back
    assert all(r.moltype for r in result.records)


@pytest.mark.skipif(not os.environ.get("NCBI_EMAIL"), reason="NCBI_EMAIL not set")
def test_ictv_chapter_returns_real_text():
    cfg = Config(email=os.environ["NCBI_EMAIL"])
    client = ICTVClient(cfg, cache=Cache(CACHE_DIR))
    chapter = client.fetch_chapter("Coronaviridae")
    assert chapter.title == "Family: Coronaviridae"
    assert "## Summary" in chapter.markdown
