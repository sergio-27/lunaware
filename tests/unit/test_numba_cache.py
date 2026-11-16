from __future__ import annotations

from lunarenv import _numba_cache


def test_stale_cache_files_are_removed(tmp_path, monkeypatch):
    package = tmp_path / "pkg"
    (package / "__pycache__").mkdir(parents=True)
    source = package / "kernels.py"
    source.write_text("import numba\n")
    stale = package / "__pycache__" / "kernels.f-1.py312.nbi"
    stale.write_text("x")
    monkeypatch.setattr(_numba_cache, "_PACKAGE", package)
    monkeypatch.setattr(_numba_cache, "_STAMP", package / "__pycache__" / "stamp")
    monkeypatch.delenv("NUMBA_CACHE_DIR", raising=False)
    assert _numba_cache.invalidate_if_stale() is True and not stale.exists()
    stale.write_text("x")
    assert _numba_cache.invalidate_if_stale() is False and stale.exists()  # sources unchanged
    source.write_text("import numba\nX = 1\n")
    assert _numba_cache.invalidate_if_stale() is True and not stale.exists()
