"""Identity of a run and of its results.

``run_id`` answers "what was asked, of which code, with which data": it is a hash of the
resolved scenario, the checksums of every input file, and the versions of the package
and of each model. Two runs with the same ``run_id`` were given identical inputs.

``result_checksum`` answers "what came out": a hash of the numerical arrays of the
result (names, shapes, types and bytes), independent of file format and timestamps.

Guarantee: on the same platform with the locked dependency set (``uv.lock``), equal
``run_id`` implies equal ``result_checksum``; the test-suite checks this. Across
platforms or library versions the arrays may differ in the last bits, which the
checksum will show; the recorded environment tells such cases apart.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np


def canonical_json(value: Any) -> str:
    """Deterministic JSON: sorted keys, no insignificant whitespace, no NaN."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_hex(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def compute_run_id(
    scenario: Mapping[str, Any],
    data_checksums: Mapping[str, str],
    model_versions: Mapping[str, str],
    package_version: str,
) -> str:
    """Identifier of a run: first 16 hex digits of the hash of everything it depends on."""
    payload = {
        "scenario": scenario,
        "data": dict(sorted(data_checksums.items())),
        "models": dict(sorted(model_versions.items())),
        "lunarenv": package_version,
    }
    return sha256_hex(canonical_json(payload))[:16]


def arrays_checksum(arrays: Mapping[str, np.ndarray]) -> str:
    """Hash of named arrays: order-independent, sensitive to dtype, shape and every bit."""
    digest = hashlib.sha256()
    for name in sorted(arrays):
        array = np.ascontiguousarray(arrays[name])
        if array.dtype.kind in "OU":
            array = np.array([str(v) for v in array.ravel()], dtype="S").reshape(array.shape)
        digest.update(name.encode("utf-8"))
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.tobytes())
    return digest.hexdigest()


def source_revision() -> str | None:
    """Git commit of the source tree the package runs from, if it is a checkout.

    ``None`` for an installed distribution (the package version then identifies the
    code) or when git is unavailable. A trailing ``+dirty`` marks uncommitted changes.
    """
    root = Path(__file__).resolve().parents[3]
    if not (root / ".git").exists():
        return None
    try:
        head = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5
        )
        if head.returncode != 0:
            return None
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        return head.stdout.strip() + ("+dirty" if status.stdout.strip() else "")
    except (OSError, subprocess.SubprocessError):
        return None


def environment() -> dict[str, str | None]:
    """Software environment recorded with every run."""
    import numba
    import spiceypy
    import xarray

    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": np.__version__,
        "numba": numba.__version__,
        "xarray": xarray.__version__,
        "spiceypy": spiceypy.__version__,
        "spice_toolkit": str(spiceypy.tkvrsn("TOOLKIT")),
        "source_revision": source_revision(),
    }
