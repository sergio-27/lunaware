"""Extract Diviner Polar Cumulative Product (PCP) bins around the Shackleton rim test point.

Reads the 96 summer diurnal PCP tables (``pcp_avg_tbol_pols_sum_slonNN_240.tab``, PDS
bundle ``urn-nasa-pds-lro_diviner_derived1``, collection ``data_derived_pcp``) with HTTP
range requests: the tables are fixed-length ASCII records sorted by x then y, so a
bisection on x finds the three grid columns of interest without downloading the
~190 MB files.

Output: ``src/lunarenv/validation/reference/diviner_pcp_shackleton.json`` with, for each
of the 3 x 3 grid cells centred on the cell that contains the test point and each
subsolar-longitude bin, the bolometric brightness temperature, plus provenance.

Run:  uv run python benchmarks/extract_diviner_pcp.py
      uv run python benchmarks/extract_diviner_pcp.py --sites   (TH-V3e blocks)
Prints only counts and coordinates, not temperatures.

With ``--sites`` the blocks are read from the selection already recorded in
``diviner_pcp_sites.json`` (``benchmarks/select_th_v3e_blocks.py``) and the values are
added to that file.
"""

from __future__ import annotations

import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src/lunarenv/validation/reference/diviner_pcp_shackleton.json"
SITES = ROOT / "src/lunarenv/validation/reference/diviner_pcp_sites.json"
BASE = (
    "https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/"
    "data_derived_pcp/diurnal/slon/pols/"
)
NAME = "pcp_avg_tbol_pols_sum_slon{:02d}_240.tab"
RECORD = 59  # bytes, from the PDS3 label (RECORD_BYTES); the first record is the header
RADIUS_M = 1737400.0
PIXEL_M = 240.0
STEP = PIXEL_M / RADIUS_M  # grid spacing in normalized units
TEST_POINT_M = (-2822.5, -9062.5)
CENTRE = (round(TEST_POINT_M[0] / PIXEL_M), round(TEST_POINT_M[1] / PIXEL_M))  # (-12, -38)
CELLS = [(CENTRE[0] + i, CENTRE[1] + j) for i in (-1, 0, 1) for j in (-1, 0, 1)]


def _curl(*args: str) -> bytes:
    """curl with TLS verification against the system store (the bundled Python may lack it)."""
    for attempt in range(5):
        done = subprocess.run(
            ["curl", "-sfL", "--max-time", "120", *args], capture_output=True, check=False
        )
        if done.returncode == 0:
            return done.stdout
        time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"curl failed: {' '.join(args)}")


def _get(url: str, start: int, stop: int) -> bytes:
    data = _curl("-r", f"{start}-{stop - 1}", url)
    if len(data) != stop - start:
        raise RuntimeError(f"short read of {url} [{start}, {stop}): {len(data)} bytes")
    return data


def _size(url: str) -> int:
    for line in _curl("-I", url).decode("ascii", "replace").splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":", 1)[1])
    raise RuntimeError(f"no Content-Length for {url}")


def _x_at(url: str, row: int) -> float:
    return float(_get(url, RECORD * (row + 1), RECORD * (row + 1) + 9))


def _first_row_with_x_at_least(url: str, rows: int, x: float) -> int:
    lo, hi = 0, rows
    while lo < hi:
        mid = (lo + hi) // 2
        if _x_at(url, mid) < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def extract(slon: int, cells: list[tuple[int, int]] = CELLS) -> dict[str, object]:
    url = BASE + NAME.format(slon)
    size = _size(url)
    rows = size // RECORD - 1
    kx = sorted({c[0] for c in cells})
    ky = sorted({c[1] for c in cells})
    x_lo = (kx[0] - 0.5) * STEP
    x_hi = (kx[-1] + 0.5) * STEP
    start = _first_row_with_x_at_least(url, rows, x_lo)
    stop = _first_row_with_x_at_least(url, rows, x_hi)
    chunk = _get(url, RECORD * (start + 1), RECORD * (stop + 1)).decode("ascii")
    values: dict[str, float] = {}
    for line in chunk.splitlines():
        x, y, _, _, tbol = (float(v) for v in line.split(","))
        cell = (round(x / STEP), round(y / STEP))
        if cell in cells and ky[0] <= cell[1] <= ky[-1]:
            values[f"{cell[0]},{cell[1]}"] = tbol
    return {
        "file": NAME.format(slon),
        "bytes": size,
        "rows_read": [start, stop],
        "sunlon_deg": [3.75 * (slon - 1), 3.75 * slon],
        "tbol_k": values,
    }


def block_cells(centre: tuple[int, int]) -> list[tuple[int, int]]:
    return [(centre[0] + i, centre[1] + j) for i in (-1, 0, 1) for j in (-1, 0, 1)]


def main_sites() -> int:
    data = json.loads(SITES.read_text(encoding="utf-8"))
    data["source"] = (
        "Diviner PCP (diurnal, summer), LRO-L-DLRE-5-GCP-V1.0 / "
        "urn:nasa:pds:lro_diviner_derived1:data_derived_pcp"
    )
    data["url"] = BASE
    data["sites"] = {}
    for label, block in data["selection"].items():
        cells = block_cells(tuple(block["centre_cell"]))
        with ThreadPoolExecutor(6) as pool:
            bins = list(pool.map(lambda n, c=cells: extract(n, c), range(1, 97)))
        counts = [len(entry["tbol_k"]) for entry in bins]  # type: ignore[arg-type]
        print(
            f"{label}: bins with all 9 cells {sum(c == 9 for c in counts)}, any {sum(c > 0 for c in counts)}"
        )
        data["sites"][label] = {
            "cells": [
                {"key": f"{i},{j}", "x_m": i * PIXEL_M, "y_m": j * PIXEL_M} for i, j in cells
            ],
            "bins": bins,
        }
    data["retrieved_utc"] = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    SITES.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {SITES.relative_to(ROOT)}")
    return 0


def main() -> int:
    with ThreadPoolExecutor(6) as pool:  # six files at a time; each file is read serially
        bins = list(pool.map(extract, range(1, 97)))
    for slon, entry in enumerate(bins, start=1):
        print(f"slon{slon:02d}: {len(entry['tbol_k'])} of {len(CELLS)} cells")  # type: ignore[arg-type]
    out = {
        "source": "Diviner PCP (diurnal, summer), LRO-L-DLRE-5-GCP-V1.0 / "
        "urn:nasa:pds:lro_diviner_derived1:data_derived_pcp",
        "url": BASE,
        "label": "Nadir observations 2009-07-05 to 2019-02-17 with subsolar latitude below "
        "the equator (summer); 240 m polar stereographic bins (radius 1737.4 km); "
        "3.75 deg of subsolar longitude per file; empty bins omitted. Radiances are "
        "binned per channel and the bolometric brightness temperature is derived from "
        "the binned radiances (Diviner RDR SIS, section 9.1).",
        "retrieved_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pixel_m": PIXEL_M,
        "radius_m": RADIUS_M,
        "test_point_m": list(TEST_POINT_M),
        "cells": [{"key": f"{i},{j}", "x_m": i * PIXEL_M, "y_m": j * PIXEL_M} for i, j in CELLS],
        "bins": bins,
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    import sys

    raise SystemExit(main_sites() if "--sites" in sys.argv[1:] else main())
