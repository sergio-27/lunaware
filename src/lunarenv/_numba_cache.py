"""Discard numba's on-disk cache when any compiled source in the package changes.

numba's ``cache=True`` keys a cached kernel on the source file that defines it, not on
the files of the functions it calls. A kernel in ``thermal/solver.py`` that inlines
``heat_capacity`` from ``thermal/properties.py`` therefore keeps the old compiled code
after ``properties.py`` changes, silently. To prevent that, a hash of every module that
uses numba is kept next to the cache; when it differs, the package's cache files are
deleted and numba recompiles on first use.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

_PACKAGE = Path(__file__).resolve().parent
_STAMP = _PACKAGE / "__pycache__" / "numba_sources.sha256"


def _sources_digest() -> str:
    digest = hashlib.sha256()
    for path in sorted(_PACKAGE.rglob("*.py")):
        data = path.read_bytes()
        if b"numba" in data:
            digest.update(str(path.relative_to(_PACKAGE)).encode())
            digest.update(data)
    return digest.hexdigest()


def _cache_dirs() -> list[Path]:
    dirs = [_PACKAGE]
    custom = os.environ.get("NUMBA_CACHE_DIR")
    if custom:
        dirs.append(Path(custom))
    return dirs


def invalidate_if_stale() -> bool:
    """Delete stale numba cache files of this package; True if anything was deleted."""
    try:
        current = _sources_digest()
        if _STAMP.is_file() and _STAMP.read_text().strip() == current:
            return False
        removed = False
        for root in _cache_dirs():
            for pattern in ("*.nbi", "*.nbc"):
                for path in root.rglob(pattern):
                    if root == _PACKAGE or "lunarenv" in str(path):
                        path.unlink(missing_ok=True)
                        removed = True
        _STAMP.parent.mkdir(exist_ok=True)
        _STAMP.write_text(current + "\n")
        return removed
    except OSError:  # read-only installation: numba then caches elsewhere or not at all
        return False
