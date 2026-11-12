"""Terrain-uncertainty ensemble for the three study sites, and its convergence with size.

Uses every Site 4 error clone present in the data cache. Writes
docs/convergence/illumination_ensemble.json.

Run:  LUNARENV_DATA_DIR=.data uv run python benchmarks/illumination_ensemble.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.illumination.ensemble import QUANTILES, cached_clones, terrain_ensemble
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_statistics
from lunarenv.illumination.settings import IlluminationSettings, load_site_stack

ROOT = Path(__file__).resolve().parents[1]
SITES = {"A": (-2847.5, -9052.5), "B": (-2647.5, -5952.5), "C": (-1000.0, -7000.0)}
HEIGHTS_M = (2.0, 10.0)


def main() -> int:
    store = DataStore()
    settings = IlluminationSettings()
    stack = load_site_stack(store, settings)
    clone_ids = [e.id for e in store.catalog.group("site04-clones") if store.is_available(e.id)]
    xy = np.array(list(SITES.values()))
    with spice_session(store) as ephemeris:
        utc = np.arange(
            np.datetime64("2026-01-01T00:00"),
            np.datetime64("2027-01-01T00:00"),
            np.timedelta64(60, "m"),
        )
        et = ephemeris.utc_to_et(utc)
        sun = ephemeris.position_from_center("SUN", et)
        radius = ephemeris.body_radius_m("SUN")

    out: dict = {"members": len(clone_ids), "quantiles": QUANTILES, "sites": SITES, "heights": {}}
    for height in HEIGHTS_M:
        nominal = illumination_statistics(
            compute_horizons(
                stack,
                xy[:, 0],
                xy[:, 1],
                height_above_surface_m=height,
                n_azimuth=settings.n_azimuth,
                diagnostics=False,
            ),
            et,
            sun,
            radius,
        )
        ensemble = terrain_ensemble(
            stack,
            cached_clones(store, clone_ids),
            xy[:, 0],
            xy[:, 1],
            et,
            sun,
            radius,
            height_above_surface_m=height,
            settings=settings,
        )
        mean = 100.0 * ensemble.mean_visible_fraction
        lcsp = ensemble.longest_shadow_s / 86400.0
        block = {
            "nominal_mean_pct": (100 * nominal.mean_visible_fraction).tolist(),
            "nominal_lcsp_d": (nominal.longest_shadow_s / 86400).tolist(),
            "ensemble_mean_pct": {
                "mean": mean.mean(0).tolist(),
                "std": mean.std(0, ddof=1).tolist(),
                "quantiles": np.quantile(mean, QUANTILES, axis=0).tolist(),
            },
            "ensemble_lcsp_d": {
                "mean": lcsp.mean(0).tolist(),
                "std": lcsp.std(0, ddof=1).tolist(),
                "quantiles": np.quantile(lcsp, QUANTILES, axis=0).tolist(),
            },
            # running estimates: how the ensemble mean and spread settle as members are added
            "running": [
                {
                    "n": n,
                    "mean_pct": mean[:n].mean(0).tolist(),
                    "std_pct": mean[:n].std(0, ddof=1).tolist(),
                }
                for n in range(4, len(clone_ids) + 1)
                if n in (4, 8, 13, 16, 25, 32, 50, 64, 75, 100) or n == len(clone_ids)
            ],
        }
        out["heights"][f"{height:g}"] = block
        for i, name in enumerate(SITES):
            print(
                f"dz={height:>4g} m site {name}: nominal {block['nominal_mean_pct'][i]:6.2f} % | "
                f"ensemble {mean[:, i].mean():6.2f} +/- {mean[:, i].std(ddof=1):.2f} % "
                f"[{mean[:, i].min():.2f}, {mean[:, i].max():.2f}] | LCSP nominal "
                f"{block['nominal_lcsp_d'][i]:.2f} d, ensemble {lcsp[:, i].mean():.2f} +/- {lcsp[:, i].std(ddof=1):.2f} d"
            )
    target = ROOT / "docs" / "convergence"
    (target / "illumination_ensemble.json").write_text(json.dumps(out, indent=1) + "\n")
    print("members:", len(clone_ids))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
