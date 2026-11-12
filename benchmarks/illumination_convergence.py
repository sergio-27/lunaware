"""Convergence and sensitivity studies of the illumination model.

Each study varies one numerical parameter around a baseline and compares against the
finest setting of that study, at three observers of different character in the Site 4
tile. Metrics: horizon profile difference, mean solar illumination over 2026, longest
continuous shadow and illumination periods.

Run:  LUNARENV_DATA_DIR=.data uv run python benchmarks/illumination_convergence.py
Output: docs/convergence/illumination.json (raw) and illumination_tables.md.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from lunarenv.data.dem import Bounds, load_dem
from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.illumination.disc import CHORD, PROFILE
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_series
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack, coarsen

OUT = Path(__file__).resolve().parents[1] / "docs" / "convergence"
SITES = {
    "A (rim crest, best lit)": (-2847.5, -9052.5),
    "B (rim crest, north)": (-2647.5, -5952.5),
    "C (inner flank, poorly lit)": (-1000.0, -7000.0),
}
HEIGHT_ABOVE_SURFACE_M = 2.0
ARCSEC = np.degrees(1.0) * 3600.0


@dataclass(frozen=True)
class Config:
    n_azimuth: int = 1440
    step_fraction: float = 0.5
    near_factor: int = 1  # block-averaging of the 5 m near DEM
    near_radius_m: float = 6000.0
    mid: str = "80"  # "80" or "40" (m/pix)
    mid_factor: int = 1
    mid_radius_m: float = 100_000.0
    far_factor: int = 3  # block-averaging of the 80 m DEM for the far level
    n_strips: int = 32
    step_minutes: int = 60
    mode: int = PROFILE
    limb_u: float = 0.6


class Study:
    def __init__(self) -> None:
        self.store = DataStore()
        self.near = load_dem(self.store, "pgda78.site04.surf")
        self.regional = load_dem(self.store, "pgda90.ldem_80s_80mpp")
        self.fine: object | None = None
        self.xy = np.array(list(SITES.values()))
        self.surface = self.near.sample(self.xy[:, 0], self.xy[:, 1])
        self._ephemeris: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        with spice_session(self.store) as ephemeris:
            self.sun_radius = ephemeris.body_radius_m("SUN")
            for minutes in (240, 120, 60, 30, 15, 5, 1):
                utc = np.arange(
                    np.datetime64("2026-01-01T00:00"),
                    np.datetime64("2027-01-01T00:00"),
                    np.timedelta64(minutes, "m"),
                )
                et = ephemeris.utc_to_et(utc)
                self._ephemeris[minutes] = (et, ephemeris.position_from_center("SUN", et))
        self._coarse: dict[tuple[str, int], object] = {}

    def level_dem(self, which: str, factor: int):
        key = (which, factor)
        if key not in self._coarse:
            if which == "near":
                base = self.near
            elif which == "80":
                base = self.regional
            else:
                if self.fine is None:
                    window = Bounds.around(-3000.0, -7500.0, 130_000.0)
                    self.fine = load_dem(self.store, "pgda90.ldem_80s_40mpp", window)
                base = self.fine
            self._coarse[key] = coarsen(base, factor)
        return self._coarse[key]

    def stack(self, c: Config) -> TerrainStack:
        near = self.level_dem("near", c.near_factor)
        mid = self.level_dem(c.mid, c.mid_factor)
        far = self.level_dem("80", c.far_factor)
        reach = 860_000.0
        return TerrainStack(
            (
                TerrainLevel(near, 0.0, c.near_radius_m, c.step_fraction * near.grid.pixel_size_m),
                TerrainLevel(
                    mid, c.near_radius_m, c.mid_radius_m, c.step_fraction * mid.grid.pixel_size_m
                ),
                TerrainLevel(far, c.mid_radius_m, reach, c.step_fraction * far.grid.pixel_size_m),
            )
        )

    def run(self, c: Config) -> dict:
        start = time.perf_counter()
        horizons = compute_horizons(
            self.stack(c),
            self.xy[:, 0],
            self.xy[:, 1],
            height_above_surface_m=HEIGHT_ABOVE_SURFACE_M,
            n_azimuth=c.n_azimuth,
            surface_height_m=self.surface,
        )
        seconds = time.perf_counter() - start
        et, sun = self._ephemeris[c.step_minutes]
        sites = []
        for i in range(len(SITES)):
            series = illumination_series(
                horizons,
                et,
                sun,
                self.sun_radius,
                observer=i,
                n_strips=c.n_strips,
                limb_u=c.limb_u,
                mode=c.mode,
            )
            sites.append(
                {
                    "mean_pct": 100.0 * series.mean_visible_fraction,
                    "mean_flux_pct": 100.0 * float(series.flux_fraction.mean()),
                    "lcsp_d": series.longest_shadow().duration_s / 86400.0,
                    "lcip_d": series.longest_illumination().duration_s / 86400.0,
                }
            )
        return {
            "horizon": horizons.elevation_rad.astype(np.float64),
            "complete": float(horizons.complete.mean()),
            "seconds_per_observer": seconds / len(SITES),
            "sites": sites,
        }


def resample(profile: np.ndarray, n: int) -> np.ndarray:
    """Linear periodic interpolation of a horizon profile onto ``n`` azimuths."""
    m = profile.shape[-1]
    position = np.arange(n) * m / n
    i0 = np.floor(position).astype(int) % m
    fraction = position - np.floor(position)
    return profile[..., i0] * (1 - fraction) + profile[..., (i0 + 1) % m] * fraction


def compare(result: dict, reference: dict) -> dict:
    n = max(result["horizon"].shape[1], reference["horizon"].shape[1])
    difference = resample(result["horizon"], n) - resample(reference["horizon"], n)
    return {
        "horizon_rms_arcsec": [float(np.sqrt(np.mean(d**2)) * ARCSEC) for d in difference],
        "horizon_max_arcsec": [float(np.abs(d).max() * ARCSEC) for d in difference],
        "horizon_bias_arcsec": [float(d.mean() * ARCSEC) for d in difference],
        "d_mean_pp": [
            a["mean_pct"] - b["mean_pct"] for a, b in zip(result["sites"], reference["sites"])
        ],
        "d_mean_flux_pp": [
            a["mean_flux_pct"] - b["mean_flux_pct"]
            for a, b in zip(result["sites"], reference["sites"])
        ],
        "d_lcsp_h": [
            24 * (a["lcsp_d"] - b["lcsp_d"]) for a, b in zip(result["sites"], reference["sites"])
        ],
        "d_lcip_h": [
            24 * (a["lcip_d"] - b["lcip_d"]) for a, b in zip(result["sites"], reference["sites"])
        ],
    }


STUDIES = {
    "azimuths": ("n_azimuth", [180, 360, 720, 1440, 2880, 5760, 11520]),
    "ray step (pixels)": ("step_fraction", [2.0, 1.0, 0.5, 0.25, 0.125]),
    "near DEM block-averaging (x 5 m)": ("near_factor", [8, 4, 2, 1]),
    "near radius (m)": ("near_radius_m", [500.0, 1000.0, 2000.0, 4000.0, 6000.0, 8000.0]),
    "mid DEM block-averaging (x 80 m)": ("mid_factor", [3, 2, 1]),
    "mid DEM product (m/pix)": ("mid", ["80", "40"]),
    "mid radius (m)": ("mid_radius_m", [20_000.0, 50_000.0, 100_000.0, 200_000.0, 400_000.0]),
    "far DEM block-averaging (x 80 m)": ("far_factor", [6, 3, 2, 1]),
    "disc strips": ("n_strips", [4, 8, 16, 32, 64, 128, 512]),
    "time step (min)": ("step_minutes", [240, 120, 60, 30, 15, 5, 1]),
}
MODEL_FORM = {
    "chord disc (Mazarico2011)": {"mode": CHORD},
    "uniform disc (no limb darkening)": {"limb_u": 0.0},
}


def main() -> int:
    only = set(sys.argv[1:])
    study = Study()
    baseline = Config()
    OUT.mkdir(parents=True, exist_ok=True)
    out_json = OUT / "illumination.json"
    results = json.loads(out_json.read_text()) if out_json.exists() else {}
    base = study.run(baseline)
    results["baseline"] = {
        "config": baseline.__dict__,
        "sites": base["sites"],
        "complete": base["complete"],
    }
    for name, (field, values) in STUDIES.items():
        if only and name not in only:
            continue
        if field == "mid" and not study.store.is_available("pgda90.ldem_80s_40mpp"):
            print(f"skip {name}: 40 m product not cached")
            continue
        runs = [study.run(replace(baseline, **{field: v})) for v in values]
        reference = runs[-1]
        rows = []
        for value, run in zip(values, runs):
            rows.append(
                {
                    "value": value,
                    "seconds_per_observer": run["seconds_per_observer"],
                    "complete": run["complete"],
                    "sites": run["sites"],
                    **compare(run, reference),
                }
            )
        results[name] = {"field": field, "rows": rows}
        print(f"done: {name}", flush=True)
    if not only:
        rows = []
        for name, change in MODEL_FORM.items():
            run = study.run(replace(baseline, **change))
            rows.append({"value": name, "sites": run["sites"], **compare(run, base)})
        results["model form (difference from baseline)"] = {"field": "model", "rows": rows}
    out_json.write_text(json.dumps(results, indent=1))
    write_tables(results)
    return 0


def write_tables(results: dict) -> None:
    names = list(SITES)
    lines = [
        "<!-- generated by benchmarks/illumination_convergence.py; do not edit -->",
        "",
        "Sites: "
        + "; ".join(
            f"**{n[0]}** = {n[3:-1]} at ({x:.1f}, {y:.1f}) m" for n, (x, y) in SITES.items()
        )
        + f". Observer {HEIGHT_ABOVE_SURFACE_M:g} m above the surface. Period: 2026.",
        "",
        "Baseline illumination: "
        + "; ".join(
            f"{n[0]}: mean {s['mean_pct']:.2f} %, longest shadow {s['lcsp_d']:.2f} d, "
            f"longest illumination {s['lcip_d']:.1f} d"
            for n, s in zip(names, results["baseline"]["sites"])
        ),
        "",
    ]
    for name, block in results.items():
        if name == "baseline":
            continue
        lines += [f"### {name}", ""]
        lines.append(
            "| value | horizon RMS (arcsec) A / B / C | horizon max (arcsec) A / B / C "
            "| d mean illumination (pp) A / B / C | d mean flux-weighted (pp) A / B / C "
            "| d longest shadow (h) A / B / C | s per observer |"
        )
        lines.append("|---|---|---|---|---|---|---|")
        for row in block["rows"]:
            fmt = lambda key, spec: " / ".join(format(v, spec) for v in row[key])  # noqa: E731
            lines.append(
                f"| {row['value']} | {fmt('horizon_rms_arcsec', '.2f')} | {fmt('horizon_max_arcsec', '.1f')} "
                f"| {fmt('d_mean_pp', '+.3f')} | {fmt('d_mean_flux_pp', '+.3f')} | {fmt('d_lcsp_h', '+.1f')} "
                f"| {row.get('seconds_per_observer', float('nan')):.3f} |"
            )
        lines.append("")
    (OUT / "illumination_tables.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
