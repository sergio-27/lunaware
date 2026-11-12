"""Mean-illumination map of a window on the Shackleton rim, with timings.

Computes horizons for every 5 m pixel of a square window and the 2026 statistics
(hourly). Writes docs/benchmarks/site_map.json (timings) and .data/site_map_2026.npz.

Run:  LUNARENV_DATA_DIR=.data uv run python benchmarks/site_map.py [half_width_m] [n_azimuth]
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numba
import numpy as np

from lunarenv.data.dem import Bounds, load_dem
from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_statistics
from lunarenv.illumination.terrain import build_stack

CENTER = (-2847.5, -9052.5)
ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    half_width = float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0
    n_azimuth = int(sys.argv[2]) if len(sys.argv) > 2 else 1440
    store = DataStore()
    near = load_dem(store, "pgda78.site04.surf")
    regional = load_dem(store, "pgda90.ldem_80s_80mpp")
    stack = build_stack(near, regional, near_radius_m=8000.0, mid_radius_m=100_000.0)
    window = near.crop(Bounds.around(CENTER[0], CENTER[1], half_width))
    xx, yy = np.meshgrid(window.grid.x_centers(), window.grid.y_centers())

    compute_horizons(stack, xx.ravel()[:2], yy.ravel()[:2], n_azimuth=16)  # compile
    start = time.perf_counter()
    horizons = compute_horizons(
        stack,
        xx.ravel(),
        yy.ravel(),
        height_above_surface_m=2.0,
        n_azimuth=n_azimuth,
        diagnostics=False,
    )
    t_horizon = time.perf_counter() - start

    with spice_session(store) as ephemeris:
        utc = np.arange(
            np.datetime64("2026-01-01T00:00"),
            np.datetime64("2027-01-01T00:00"),
            np.timedelta64(60, "m"),
        )
        et = ephemeris.utc_to_et(utc)
        sun = ephemeris.position_from_center("SUN", et)
        radius = ephemeris.body_radius_m("SUN")
    illumination_statistics(
        type(horizons)(
            horizons.x_m[:2],
            horizons.y_m[:2],
            horizons.observer_height_m[:2],
            horizons.elevation_rad[:2],
            horizons.complete[:2],
            None,
            None,
        ),
        et[:4],
        sun[:4],
        radius,
    )  # compile
    start = time.perf_counter()
    stats = illumination_statistics(horizons, et, sun, radius)
    t_stats = time.perf_counter() - start

    shape = xx.shape
    np.savez_compressed(
        store.root / "site_map_2026.npz",
        x=window.grid.x_centers(),
        y=window.grid.y_centers(),
        height=window.heights_m,
        mean=stats.mean_visible_fraction.reshape(shape).astype(np.float32),
        mean_flux=stats.mean_flux_fraction.reshape(shape).astype(np.float32),
        lcsp_d=(stats.longest_shadow_s / 86400).reshape(shape).astype(np.float32),
        lcip_d=(stats.longest_illumination_s / 86400).reshape(shape).astype(np.float32),
        complete=horizons.complete.mean(axis=1).reshape(shape).astype(np.float32),
    )
    timing = {
        "observers": int(xx.size),
        "n_azimuth": n_azimuth,
        "epochs": int(et.size),
        "threads": numba.get_num_threads(),
        "cpu_count": os.cpu_count(),
        "horizon_s": t_horizon,
        "statistics_s": t_stats,
        "horizon_ms_per_observer": 1e3 * t_horizon / xx.size,
        "rays_per_s": xx.size * n_azimuth / t_horizon,
        "observer_epochs_per_s": xx.size * et.size / t_stats,
    }
    out = ROOT / "docs" / "benchmarks"
    out.mkdir(parents=True, exist_ok=True)
    name = f"site_map_{int(2 * half_width)}m_{n_azimuth}az.json"
    (out / name).write_text(json.dumps(timing, indent=1) + "\n")
    print(json.dumps(timing, indent=1))
    best = int(np.argmax(stats.mean_visible_fraction))
    print(
        f"best pixel x={xx.ravel()[best]:.1f} y={yy.ravel()[best]:.1f} "
        f"mean={100 * stats.mean_visible_fraction[best]:.2f}% "
        f"lcsp={stats.longest_shadow_s[best] / 86400:.2f} d"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
