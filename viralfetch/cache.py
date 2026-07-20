"""On-disk cache with a simple TTL.

Two namespaces (SPEC section 3):

- ``SEQS``  — sequences and accession metadata. Accessions are immutable, so
  these entries are cached **permanently** (no TTL).
- ``TEXTS`` — ICTV chapter text. Cached with a **30-day TTL**.

The store is intentionally dumb: one file per entry, keyed by a hash of the
logical key; expiry is decided from the file mtime. No "smart" invalidation
(SPEC section 10).
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass
from pathlib import Path

SEQS = "seqs"
TEXTS = "texts"

TEXT_TTL = 30 * 24 * 60 * 60  # 30 days, in seconds


def _hash(key: str) -> str:
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


@dataclass
class NamespaceInfo:
    namespace: str
    entries: int
    bytes: int


class Cache:
    """A file-backed key/value cache. Disable with ``enabled=False``."""

    def __init__(self, base_dir: Path, enabled: bool = True):
        self.base_dir = Path(base_dir)
        self.enabled = enabled

    def _path(self, namespace: str, key: str) -> Path:
        return self.base_dir / namespace / _hash(key)

    def get(self, namespace: str, key: str, ttl: int | None = None) -> str | None:
        """Return cached text, or ``None`` if absent, expired, or disabled."""
        if not self.enabled:
            return None
        path = self._path(namespace, key)
        if not path.is_file():
            return None
        if ttl is not None and (time.time() - path.stat().st_mtime) > ttl:
            return None
        try:
            return path.read_text(encoding="utf-8")
        except OSError:
            return None

    def set(self, namespace: str, key: str, value: str) -> None:
        """Store ``value`` under ``key``. No-op when disabled."""
        if not self.enabled:
            return
        path = self._path(namespace, key)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value, encoding="utf-8")

    def clear(self, *, texts: bool = False, seqs: bool = False) -> int:
        """Remove cached entries. With neither flag, clear everything.

        Returns the number of entries removed.
        """
        both = not texts and not seqs
        targets = []
        if texts or both:
            targets.append(TEXTS)
        if seqs or both:
            targets.append(SEQS)

        removed = 0
        for ns in targets:
            ns_dir = self.base_dir / ns
            if not ns_dir.is_dir():
                continue
            for entry in ns_dir.iterdir():
                if entry.is_file():
                    entry.unlink()
                    removed += 1
        return removed

    def info(self) -> list[NamespaceInfo]:
        """Per-namespace entry counts and total bytes."""
        out: list[NamespaceInfo] = []
        for ns in (SEQS, TEXTS):
            ns_dir = self.base_dir / ns
            entries = 0
            total = 0
            if ns_dir.is_dir():
                for entry in ns_dir.iterdir():
                    if entry.is_file():
                        entries += 1
                        total += entry.stat().st_size
            out.append(NamespaceInfo(namespace=ns, entries=entries, bytes=total))
        return out
