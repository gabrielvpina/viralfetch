"""Tests for the on-disk TTL cache."""

import os
import time

from viralfetch.cache import IMAGES, SEQS, TEXTS, Cache


def test_set_get_roundtrip(tmp_path):
    c = Cache(tmp_path)
    assert c.get(SEQS, "k") is None
    c.set(SEQS, "k", "value")
    assert c.get(SEQS, "k") == "value"


def test_bytes_roundtrip(tmp_path):
    c = Cache(tmp_path)
    assert c.get_bytes(IMAGES, "fig") is None
    c.set_bytes(IMAGES, "fig", b"\x89PNG\x00\xff")
    assert c.get_bytes(IMAGES, "fig") == b"\x89PNG\x00\xff"


def test_bytes_disabled_cache_is_noop(tmp_path):
    c = Cache(tmp_path, enabled=False)
    c.set_bytes(IMAGES, "fig", b"data")
    assert c.get_bytes(IMAGES, "fig") is None


def test_disabled_cache_is_noop(tmp_path):
    c = Cache(tmp_path, enabled=False)
    c.set(SEQS, "k", "value")
    assert c.get(SEQS, "k") is None


def test_ttl_expiry(tmp_path):
    c = Cache(tmp_path)
    c.set(TEXTS, "k", "value")
    path = c._path(TEXTS, "k")
    old = time.time() - 100
    os.utime(path, (old, old))
    assert c.get(TEXTS, "k", ttl=50) is None  # expired
    assert c.get(TEXTS, "k", ttl=200) == "value"  # still fresh
    assert c.get(TEXTS, "k") == "value"  # no ttl => permanent


def test_clear_selective(tmp_path):
    c = Cache(tmp_path)
    c.set(SEQS, "a", "1")
    c.set(TEXTS, "b", "2")
    removed = c.clear(texts=True)
    assert removed == 1
    assert c.get(SEQS, "a") == "1"
    assert c.get(TEXTS, "b") is None


def test_clear_selective_images(tmp_path):
    c = Cache(tmp_path)
    c.set(TEXTS, "b", "2")
    c.set_bytes(IMAGES, "fig", b"data")
    removed = c.clear(images=True)
    assert removed == 1
    assert c.get(TEXTS, "b") == "2"
    assert c.get_bytes(IMAGES, "fig") is None


def test_clear_all(tmp_path):
    c = Cache(tmp_path)
    c.set(SEQS, "a", "1")
    c.set(TEXTS, "b", "2")
    c.set_bytes(IMAGES, "fig", b"3")
    assert c.clear() == 3


def test_info(tmp_path):
    c = Cache(tmp_path)
    c.set(SEQS, "a", "hello")
    c.set_bytes(IMAGES, "fig", b"1234")
    info = {i.namespace: i for i in c.info()}
    assert info[SEQS].entries == 1
    assert info[SEQS].bytes == len("hello")
    assert info[TEXTS].entries == 0
    assert info[IMAGES].entries == 1
    assert info[IMAGES].bytes == 4
