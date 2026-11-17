"""Extract the times and positions of Diviner RDR observations over the five TH-V3 blocks.

For TH-V3h (observation timing). Queries the PDS Geosciences Node Orbital Data Explorer
(ODE) Diviner RDR point service (``oderest.rsl.wustl.edu/livegds``, the service behind
the "DIVINER RDR Query V2" page) for channel 7 records in a latitude/longitude box around
each block, 2009-07-05 to 2019-02-17, and keeps the records the summer diurnal Polar
Cumulative Products are built from (Diviner RDR SIS sections 2.2.2 and 9.1): activity
flag 110 (on Moon, standard nadir), subsolar latitude below the equator, brightness
temperature 10-450 K.

Each kept record is one observation (one detector, one 0.128 s integration). Stored per
observation: Julian date (UTC), polar stereographic x and y of the footprint centre (m),
subsolar longitude (deg), orbit. Radiances and temperatures are not stored and not printed: the
timing is extracted blind to the measured values.

Output: ``src/lunarenv/validation/reference/diviner_rdr_times.json.gz``.

Run:  uv run python benchmarks/extract_diviner_rdr_times.py
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from lunarenv.geometry.projection import SouthPolarStereographic

ROOT = Path(__file__).resolve().parents[1]
SITES = ROOT / "src/lunarenv/validation/reference/diviner_pcp_sites.json"
SHACKLETON = ROOT / "src/lunarenv/validation/reference/diviner_pcp_shackleton.json"
OUT = ROOT / "src/lunarenv/validation/reference/diviner_rdr_times.json.gz"
SERVICE = "https://oderest.rsl.wustl.edu/livegds/"
START, STOP = "2009-07-05", "2019-02-17"
PIXEL_M = 240.0
#: kept: centres within the 3 x 3 block grown by this margin (m)
MARGIN_M = 240.0
FIELDS = (
    "date,utc,jdate,orbit,sundist,sunlat,sunlon,sclk,sclat,sclon,scrad,scalt,el_cmd,az_cmd,"
    "af,vert_lat,vert_lon,c,det,vlookx,vlooky,vlookz,radiance,tb,clat,clon,cemis,csunzen,"
    "csunazi,cloctime,qca,qge,qmi"
).split(",")
PROJ = SouthPolarStereographic()


def _curl(url: str) -> bytes:
    """curl with TLS verification against the system store."""
    for attempt in range(5):
        done = subprocess.run(
            ["curl", "-sfL", "-A", "lunarenv-validation", "--max-time", "300", url],
            capture_output=True,
            check=False,
        )
        if done.returncode == 0:
            return done.stdout
        time.sleep(5.0 * (attempt + 1))
    raise RuntimeError(f"curl failed: {url}")


def _query(box: tuple[float, float, float, float]) -> tuple[str, bytes]:
    west, east, south, north = box
    url = (
        f"{SERVICE}?query=divinerrdr&source=GDSWebDLRERDR&odemetadb=moon&o=j&results=cv"
        f"&minobtime={START}&maxobtime={STOP}&channel=fffffftff&async=t"
        f"&westlon={west:.4f}&eastlon={east:.4f}&minlat={south:.5f}&maxlat={north:.5f}"
    )
    job = json.loads(_curl(url))["GDSResults"]
    if job.get("Status") != "Success" and "JobId" not in job:
        raise RuntimeError(job.get("Error", "query refused"))
    status_url = f"{SERVICE}?query=status&source=GDSWebDLRERDR&o=j&jobid={job['JobId']}"
    for _ in range(90):
        state = json.loads(_curl(status_url))["GDSResults"]
        if state.get("Status", "").lower() == "error":
            raise RuntimeError(str(state))
        if state.get("StateSummary", {}).get("State") == "Finished":
            files = state["ResultFiles"]["ResultFile"]
            files = files if isinstance(files, list) else [files]
            tab = next(f["URL"] for f in files if f["URL"].endswith(".tab"))
            return url, _curl(tab)
        time.sleep(20.0)
    raise RuntimeError("ODE job did not finish")


def _box(xc: float, yc: float, half: float) -> tuple[float, float, float, float]:
    t = np.linspace(-half, half, 21)
    xx, yy = np.meshgrid(xc + t, yc + t)
    lon, lat = PROJ.inverse(xx.ravel(), yy.ravel())
    lon = np.degrees(lon) % 360.0
    lat = np.degrees(lat)
    return float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())


def _parse(table: bytes) -> dict[str, np.ndarray]:
    lines = (line for line in io.StringIO(table.decode("ascii")) if not line.startswith("#"))
    rows = [r for r in csv.reader(lines) if r]
    cols = {k: [r[i].strip().strip('"') for r in rows] for i, k in enumerate(FIELDS)}
    out = {k: np.asarray(cols[k]) for k in ("utc", "date")}
    for k in (
        "jdate",
        "orbit",
        "sunlat",
        "sunlon",
        "scalt",
        "af",
        "c",
        "det",
        "tb",
        "clat",
        "clon",
        "cemis",
    ):
        out[k] = np.asarray(cols[k], dtype=np.float64)
    return out


def _blocks() -> dict[str, tuple[float, float]]:
    """Block centres (m): the middle cell of each site's 3 x 3 Diviner cells."""
    out = {}
    shack = json.loads(SHACKLETON.read_text(encoding="utf-8"))
    xs = [c["x_m"] for c in shack["cells"]]
    ys = [c["y_m"] for c in shack["cells"]]
    out["Shackleton"] = (float(np.mean(xs)), float(np.mean(ys)))
    sites = json.loads(SITES.read_text(encoding="utf-8"))["sites"]
    for name, site in sites.items():
        xs = [c["x_m"] for c in site["cells"]]
        ys = [c["y_m"] for c in site["cells"]]
        out[name] = (float(np.mean(xs)), float(np.mean(ys)))
    return out


def main() -> int:
    result: dict[str, object] = {
        "source": "LRO Diviner RDR (LRO-L-DLRE-4-RDR-V1.0), channel 7, via the PDS Geosciences "
        "Node ODE Diviner RDR point service",
        "service": SERVICE,
        "selection": "activity flag 110 (on Moon, standard nadir), subsolar latitude < 0, "
        f"10 K <= tb <= 450 K, {START} to {STOP}; centres within the 3 x 3 block grown by "
        f"{MARGIN_M:.0f} m",
        "fields": "jd (UTC), x_m, y_m (polar stereographic, R = 1737.4 km), sunlon_deg, orbit",
        "retrieved_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sites": {},
    }
    sites: dict[str, object] = {}
    half = 1.5 * PIXEL_M + MARGIN_M
    for name, (xc, yc) in _blocks().items():
        url, table = _query(_box(xc, yc, half))
        d = _parse(table)
        lon = np.radians(d["clon"])
        lat = np.radians(d["clat"])
        x, y = PROJ.forward(lon, lat)
        keep = (
            (d["af"] == 110)
            & (d["c"] == 7)
            & (d["sunlat"] < 0.0)
            & (d["tb"] >= 10.0)
            & (d["tb"] <= 450.0)
            & (np.abs(x - xc) <= half)
            & (np.abs(y - yc) <= half)
        )
        order = np.argsort(d["jdate"][keep], kind="stable")
        sites[name] = {
            "query": url,
            "table_sha256": hashlib.sha256(table).hexdigest(),
            "records_returned": int(d["jdate"].size),
            "jd": np.round(d["jdate"][keep][order], 7).tolist(),
            "x_m": np.round(np.asarray(x)[keep][order], 1).tolist(),
            "y_m": np.round(np.asarray(y)[keep][order], 1).tolist(),
            "sunlon_deg": np.round(d["sunlon"][keep][order], 4).tolist(),
            "orbit": d["orbit"][keep][order].astype(int).tolist(),
        }
        years = 2000 + (np.asarray(sites[name]["jd"]) - 2451544.5) / 365.25  # type: ignore[arg-type]
        print(
            f"{name}: {d['jdate'].size} returned, {int(keep.sum())} kept; per year 2009-2018:",
            np.histogram(years, bins=np.arange(2009, 2020))[0].tolist(),
        )
    result["sites"] = sites
    payload = (json.dumps(result, separators=(",", ":")) + "\n").encode("utf-8")
    OUT.write_bytes(gzip.compress(payload, mtime=0))
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
