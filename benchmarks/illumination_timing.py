"""Timing and scaling of the illumination kernels.

Measures (after compilation):
  * horizon tracing versus number of observers, azimuths and ray step;
  * time-series statistics versus number of epochs (1, 5 and 18.6 years, hourly);
  * the cost of one point with its 100-member terrain ensemble (extrapolated from one member).

Run:  LUNARENV_DATA_DIR=.data uv run python benchmarks/illumination_timing.py
Output: docs/benchmarks/illumination_timing.json
"""

from __future__ import annotations

import json
import os
import platform
import time
from pathlib import Path

import numba
import numpy as np

from lunarenv.data.dem import load_dem
from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_statistics
from lunarenv.illumination.terrain import build_stack

ROOT = Path(__file__).resolve().parents[1]
CENTER = (-2847.5, -9052.5)


def observers(n: int) -> tuple[np.ndarray, np.ndarray]:
    """n observers on a square lattice of 5 m pitch around the rim site."""
    side = int(np.ceil(np.sqrt(n)))
    offsets = (np.arange(side) - side / 2) * 5.0
    xx, yy = np.meshgrid(CENTER[0] + offsets, CENTER[1] + offsets)
    return xx.ravel()[:n], yy.ravel()[:n]


def timed(function, *args, **kwargs):
    """Result, wall seconds and CPU seconds (all threads) of a call.

    CPU time is the robust figure: wall time depends on what else the machine is doing.
    """
    wall, cpu = time.perf_counter(), time.process_time()
    result = function(*args, **kwargs)
    return result, time.perf_counter() - wall, time.process_time() - cpu


def main() -> int:
    store = DataStore()
    near = load_dem(store, "pgda78.site04.surf")
    regional = load_dem(store, "pgda90.ldem_80s_80mpp")
    out: dict = {
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "numba_threads": numba.get_num_threads(),
        "horizon": [],
        "statistics": [],
    }

    def stack(step_fraction: float):
        return build_stack(
            near,
            regional,
            near_radius_m=8000.0,
            mid_radius_m=100_000.0,
            step_fraction=step_fraction,
        )

    base = stack(0.5)
    samples = base.sample_plan()[0].size
    compute_horizons(base, *observers(2), n_azimuth=16)  # compile
    for n_obs, n_az, step in [
        (100, 1440, 0.5),
        (1000, 1440, 0.5),
        (4000, 1440, 0.5),
        (1000, 360, 0.5),
        (1000, 720, 0.5),
        (1000, 2880, 0.5),
        (1000, 1440, 1.0),
        (1000, 1440, 0.25),
    ]:
        x, y = observers(n_obs)
        _, seconds, cpu = timed(
            compute_horizons,
            stack(step),
            x,
            y,
            height_above_surface_m=2.0,
            n_azimuth=n_az,
            diagnostics=False,
        )
        out["horizon"].append(
            {
                "observers": n_obs,
                "n_azimuth": n_az,
                "step_fraction": step,
                "wall_s": seconds,
                "cpu_s": cpu,
                "cpu_ms_per_observer": 1e3 * cpu / n_obs,
                "cpu_us_per_ray": 1e6 * cpu / (n_obs * n_az),
                "effective_cores": cpu / seconds,
            }
        )
        print(out["horizon"][-1], flush=True)
    out["plan_samples_per_ray_max"] = int(samples)

    x, y = observers(1000)
    horizons = compute_horizons(
        base, x, y, height_above_surface_m=2.0, n_azimuth=1440, diagnostics=False
    )
    with spice_session(store) as ephemeris:
        radius = ephemeris.body_radius_m("SUN")
        for years, stop in [(1.0, "2027-01-01"), (5.0, "2031-01-01"), (18.6, "2044-08-08")]:
            utc = np.arange(
                np.datetime64("2026-01-01T00:00"),
                np.datetime64(f"{stop}T00:00"),
                np.timedelta64(60, "m"),
            )
            et, t_time, _ = timed(ephemeris.utc_to_et, utc)
            sun, t_ephemeris, _ = timed(ephemeris.position_from_center, "SUN", et)
            illumination_statistics(horizons, et[:8], sun[:8], radius)  # compile
            _, seconds, cpu = timed(illumination_statistics, horizons, et, sun, radius)
            out["statistics"].append(
                {
                    "years": years,
                    "epochs": int(et.size),
                    "observers": 1000,
                    "wall_s": seconds,
                    "cpu_s": cpu,
                    "cpu_ns_per_observer_epoch": 1e9 * cpu / (1000 * et.size),
                    "ephemeris_wall_s": t_time + t_ephemeris,
                }
            )
            print(out["statistics"][-1], flush=True)

    # one point, one terrain member: load a 41 MB DEM, trace, one year of statistics
    x1, y1 = observers(1)
    wall, cpu = time.perf_counter(), time.process_time()
    member = load_dem(store, "pgda78.site04.clone.0001")
    h1 = compute_horizons(
        base.with_near(member),
        x1,
        y1,
        height_above_surface_m=2.0,
        n_azimuth=1440,
        diagnostics=False,
    )
    illumination_statistics(h1, et[:8760], sun[:8760], radius)
    out["one_point_one_member"] = {
        "wall_s": time.perf_counter() - wall,
        "cpu_s": time.process_time() - cpu,
    }
    print("one member:", out["one_point_one_member"])
    out["load_average_at_end"] = os.getloadavg()

    target = ROOT / "docs" / "benchmarks"
    target.mkdir(parents=True, exist_ok=True)
    (target / "illumination_timing.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
