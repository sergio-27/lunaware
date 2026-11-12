"""Shared fixtures.

Tests marked ``data`` read real files from the local data cache. They are skipped, not
failed, when a file is absent, so the unit suite runs on a fresh checkout without
network access. Set ``LUNARENV_DATA_DIR`` to point at a populated cache; a ``.data``
directory at the repository root is used when the variable is unset.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from hypothesis import settings

from lunarenv.data.store import ENV_DATA_DIR, DataStore

# Compiled kernels are built on first call; a per-example time limit would turn that
# one-off cost into spurious failures.
settings.register_profile("lunarenv", deadline=None)
settings.load_profile("lunarenv")

_REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def store() -> DataStore:
    root = os.environ.get(ENV_DATA_DIR)
    return DataStore(Path(root) if root else _REPO_ROOT / ".data")


@pytest.fixture(scope="session")
def require(store: DataStore):
    """Return a function that skips the test unless the given catalog ids are cached."""

    def _require(*entry_ids: str) -> None:
        missing = [e for e in entry_ids if not store.is_available(e)]
        if missing:
            pytest.skip(f"not in the data cache: {', '.join(missing)}")

    return _require
