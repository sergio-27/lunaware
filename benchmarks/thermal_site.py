"""Thermal model at the Shackleton rim test point: results, convergence, uncertainty.

Period 2026 (spin-up cycle: 2025). The surface element is the DEM surface at the best-lit
pixel of the 1 x 1 km window, with the DEM's own slope.

Run:  LUNARENV_DATA_DIR=.data uv run python benchmarks/thermal_site.py
Output: docs/convergence/thermal.json, docs/figures/thermal_site_2026.png
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lunarenv.data.dem import load_dem
from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.illumination.eclipse import earth_occultation_factor
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_series
from lunarenv.illumination.settings import IlluminationSettings, load_site_stack
from lunarenv.thermal.model import ThermalSettings, run_thermal, spin_up, summarize
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.solver import solve_conduction
from lunarenv.thermal.surface import (
    absorbed_solar_flux,
    neglected_terrain_flux_bound,
    slope_and_facing,
    surface_normal,
)
from lunarenv.thermal.uq import DEFAULT_RANGES, rank_correlation, sample_properties

ROOT = Path(__file__).resolve().parents[1]
SITE = (-2822.5, -9062.5)
PROJ = SouthPolarStereographic()
YEAR_MIN = 365 * 24 * 60


class Site:
    def __init__(self) -> None:
        self.store = DataStore()
        self.settings = IlluminationSettings()
        self.stack = load_site_stack(self.store, self.settings)
        self.ephemeris: dict[int, dict] = {}
        with spice_session(self.store) as eph:
            self.r_sun, self.r_earth = eph.body_radius_m("SUN"), eph.body_radius_m("EARTH")
            for minutes in (60, 30, 15, 5):
                utc = np.datetime64("2025-01-01T00:00") + np.arange(
                    2 * YEAR_MIN // minutes + 1
                ) * np.timedelta64(minutes, "m")
                et = eph.utc_to_et(utc)
                self.ephemeris[minutes] = {
                    "et": et,
                    "sun": eph.position_from_center("SUN", et),
                    "earth": eph.position_from_center("EARTH", et),
                }

    def forcing(
        self, minutes=60, dem=None, properties=None, *, height=0.0, baseline=None, eclipse=True
    ):
        """Absorbed flux series for a terrain member; also shadow mask, normal and horizon."""
        properties = properties or RegolithProperties()
        stack = self.stack if dem is None else self.stack.with_near(dem)
        near = stack.levels[0].dem
        horizons = compute_horizons(
            stack, *SITE, height_above_surface_m=height, n_azimuth=self.settings.n_azimuth
        )
        eph = self.ephemeris[minutes]
        series = illumination_series(horizons, eph["et"], eph["sun"], self.r_sun)
        normal = surface_normal(near, *SITE, baseline_m=baseline)
        factor = None
        if eclipse:
            position = PROJ.to_cartesian(*SITE, float(horizons.observer_height_m[0]))
            factor = earth_occultation_factor(
                eph["sun"], eph["earth"], position, self.r_sun, self.r_earth
            )
        flux = absorbed_solar_flux(series, normal, properties, eclipse_factor=factor)
        return {
            "flux": flux,
            "series": series,
            "normal": normal,
            "horizons": horizons,
            "factor": factor,
        }


def metrics(result, reference=None, stride=1):
    surface = result.surface_temperature_k[::stride]
    out = {
        "surface_min": float(surface.min()),
        "surface_max": float(surface.max()),
        "surface_mean": float(surface.mean()),
    }
    for key, depth in (("mean_0.3m", 0.3), ("mean_1m", 1.0)):
        if depth <= result.depth_m[-1]:
            out[key] = float(result.at_depth(depth)[::stride].mean())
    if reference is not None:
        out["max_abs_surface_diff"] = float(np.abs(surface - reference["surface"]).max())
        out["rms_surface_diff"] = float(np.sqrt(np.mean((surface - reference["surface"]) ** 2)))
        for key in ("surface_min", "surface_max", "surface_mean", "mean_0.3m", "mean_1m"):
            if key in out:
                out["d_" + key] = out[key] - reference[key]
    return out


def main() -> int:
    site = Site()
    props = RegolithProperties()
    base = ThermalSettings()
    out: dict = {"site_xy_m": SITE, "period": "2026 (spin-up cycle 2025)"}
    lon, lat = PROJ.inverse(*SITE)

    # ------------------------------------------------------------------ nominal
    nominal = site.forcing(60)
    spin = 365 * 24
    start = time.perf_counter()
    cpu = time.process_time()
    result = run_thermal(
        nominal["flux"],
        3600.0,
        spin,
        props,
        base,
        latitude_rad=float(lat),
        shadow=nominal["series"].visible_fraction[spin:] <= 0.0,
    )
    out["nominal_wall_s"] = time.perf_counter() - start
    out["nominal_cpu_s"] = time.process_time() - cpu
    summary = summarize(result)
    slope, facing = slope_and_facing(nominal["normal"])
    bound = neglected_terrain_flux_bound(nominal["horizons"].elevation_rad[0], props)
    out["nominal"] = {
        "latitude_deg": float(np.degrees(lat)),
        "longitude_deg": float(np.degrees(lon) % 360),
        "slope_deg": float(np.degrees(slope)),
        "facing_grid_azimuth_deg": float(np.degrees(facing)),
        "mean_illumination_pct_2026": 100 * float(nominal["series"].visible_fraction[spin:].mean()),
        "mean_absorbed_flux_w_m2": float(result.absorbed_flux_w_m2.mean()),
        "max_absorbed_flux_w_m2": float(result.absorbed_flux_w_m2.max()),
        "summary": {
            k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in asdict(summary).items()
        },
        "spinup_cycles": result.spinup.cycles,
        "spinup_residual_k": result.spinup.residual_k,
        "energy_residual_relative": abs(result.energy_residual_j_m2)
        / float(result.absorbed_flux_w_m2.sum() * 3600.0),
        "warnings": list(result.warnings),
        "cells": base.grid().n_cells,
        "neglected_flux": asdict(bound),
        "eclipse_hours_2026": float(np.sum(nominal["factor"][spin:] < 1.0)),
    }
    print(json.dumps(out["nominal"], indent=1))

    # ------------------------------------------------- time integration (forcing 60 min)
    def run(flux, step, settings, spinup_samples):
        t0 = time.process_time()
        r = run_thermal(flux, step, spinup_samples, props, settings)
        return r, time.process_time() - t0

    reference_run, _ = run(
        nominal["flux"], 3600.0, replace_settings(base, scheme="bdf2", substeps=64), spin
    )
    reference = {"surface": reference_run.surface_temperature_k, **metrics(reference_run)}
    rows = []
    for scheme in ("explicit", "backward_euler", "crank_nicolson", "bdf2"):
        for substeps in (None,) if scheme == "explicit" else (1, 2, 4, 8, 16):
            settings = replace_settings(base, scheme=scheme, substeps=substeps or 1)
            if scheme == "explicit":
                grid = settings.grid()
                from lunarenv.thermal.solver import explicit_stable_step_s

                substeps_used = int(
                    np.ceil(3600.0 / (0.8 * explicit_stable_step_s(grid, props, 400.0)))
                )
                settings = replace_settings(base, scheme=scheme, substeps=substeps_used)
            r, seconds = run(nominal["flux"], 3600.0, settings, spin)
            rows.append(
                {
                    "scheme": scheme,
                    "substeps": settings.substeps,
                    "step_s": 3600.0 / settings.substeps,
                    "cpu_s": seconds,
                    **metrics(r, reference),
                }
            )
            print(
                rows[-1]["scheme"],
                rows[-1]["substeps"],
                f"{rows[-1]['max_abs_surface_diff']:.4f} K",
                f"{seconds:.2f} s",
                flush=True,
            )
    out["time_integration"] = rows

    # --------------------------------------------------- forcing sampling (solver step 5 min)
    rows = []
    reference5 = None
    for minutes in (5, 15, 30, 60):
        f = site.forcing(minutes)
        r, seconds = run(
            f["flux"],
            60.0 * minutes,
            replace_settings(base, substeps=minutes // 5),
            YEAR_MIN // minutes,
        )
        stride = 60 // minutes
        if reference5 is None:
            reference5 = {"surface": r.surface_temperature_k[::stride], **metrics(r, stride=stride)}
        rows.append(
            {
                "forcing_minutes": minutes,
                "cpu_s": seconds,
                **metrics(r, reference5, stride=stride),
                "mean_flux": float(f["flux"][YEAR_MIN // minutes :].mean()),
            }
        )
        print("forcing", minutes, f"{rows[-1]['max_abs_surface_diff']:.3f} K", flush=True)
    out["forcing_sampling"] = rows

    # -------------------------------------------------------------------------- grid
    fine_settings = replace_settings(
        base, first_cell_m=0.00025, growth=1.05, depth_m=4.0, substeps=4
    )
    fine_run, _ = run(nominal["flux"], 3600.0, fine_settings, spin)
    fine = {"surface": fine_run.surface_temperature_k, **metrics(fine_run)}
    rows = []
    variations = (
        [("first_cell_m", v) for v in (0.008, 0.004, 0.002, 0.001, 0.0005)]
        + [("growth", v) for v in (1.3, 1.2, 1.12, 1.06)]
        + [("depth_m", v) for v in (0.3, 0.5, 1.0, 2.0, 4.0)]
    )
    for name, value in variations:
        settings = replace_settings(base, substeps=4, **{name: value})
        r, seconds = run(nominal["flux"], 3600.0, settings, spin)
        row = {
            "parameter": name,
            "value": value,
            "cells": settings.grid().n_cells,
            "cpu_s": seconds,
        }
        m = metrics(r, fine)
        rows.append({**row, **m})
        print(name, value, f"{m['max_abs_surface_diff']:.4f} K", flush=True)
    out["grid"] = rows

    # ----------------------------------------------------------------------- spin-up
    cold = spin_up(nominal["flux"][: spin + 1], 3600.0, props, base, initial_temperature_k=60.0)
    hot = spin_up(nominal["flux"][: spin + 1], 3600.0, props, base, initial_temperature_k=320.0)
    grid = base.grid()
    after = [
        solve_conduction(
            grid, props, s.temperature_k, nominal["flux"][spin:], 3600.0, scheme="bdf2"
        )
        for s in (cold, hot)
    ]
    out["spinup"] = {
        "cycles": [cold.cycles, hot.cycles],
        "converged": [cold.converged, hot.converged],
        "max_state_difference_k": float(np.abs(cold.temperature_k - hot.temperature_k).max()),
        "max_surface_difference_in_period_k": float(
            np.abs(after[0].surface_temperature_k - after[1].surface_temperature_k).max()
        ),
        "max_any_depth_difference_in_period_k": float(
            np.abs(after[0].temperature_k - after[1].temperature_k).max()
        ),
    }
    print("spin-up", out["spinup"])

    # -------------------------------------------------------------------- model form
    def variant(label, **kwargs):
        p = kwargs.pop("properties", props)
        f = site.forcing(60, properties=p, **kwargs)
        r = run_thermal(f["flux"], 3600.0, spin, p, base)
        return {
            "variant": label,
            **metrics(r, {"surface": result.surface_temperature_k, **metrics(result)}),
            "mean_flux": float(f["flux"][spin:].mean()),
        }

    out["model_form"] = [
        variant(
            "constant albedo (a = b = 0)", properties=RegolithProperties(albedo_a=0.0, albedo_b=0.0)
        ),
        variant("no lunar eclipses", eclipse=False),
        variant("slope over 20 m baseline", baseline=20.0),
        variant("illumination seen from 2 m above the surface", height=2.0),
    ]
    for row in out["model_form"]:
        print(row["variant"], {k: round(v, 2) for k, v in row.items() if k.startswith("d_")})

    # ------------------------------------------------------------------- uncertainty
    clone_ids = [
        e.id for e in site.store.catalog.group("site04-clones") if site.store.is_available(e.id)
    ]
    members = [nominal] + [site.forcing(60, dem=load_dem(site.store, c)) for c in clone_ids]
    slopes = [np.degrees(slope_and_facing(m["normal"])[0]) for m in members]
    n = 128
    properties, values = sample_properties(n, seed=20261112)
    rng = np.random.default_rng(20261112)
    terrain_index = rng.integers(1, len(members), size=n) if len(members) > 1 else np.zeros(n, int)
    start = time.process_time()
    outputs = []
    surfaces = np.empty((n, YEAR_MIN // 60 + 1), dtype=np.float32)
    for i, p in enumerate(properties):
        member = members[int(terrain_index[i])]
        # the albedo of this member changes the absorbed fraction of the same illumination
        flux = absorbed_solar_flux(
            member["series"], member["normal"], p, eclipse_factor=member["factor"]
        )
        r = run_thermal(flux, 3600.0, spin, p, base)
        s = summarize(r)
        surfaces[i] = r.surface_temperature_k
        outputs.append(
            [
                s.surface_min_k,
                s.surface_max_k,
                s.surface_mean_k,
                s.lunation_amplitude_mean_k,
                float(r.at_depth(0.3).mean()),
                float(r.at_depth(1.0).mean()),
                s.max_heating_rate_k_per_h,
            ]
        )
    uq_cpu = time.process_time() - start
    outputs = np.array(outputs)
    names = [
        "surface_min",
        "surface_max",
        "surface_mean",
        "lunation_amplitude_mean",
        "mean_0.3m",
        "mean_1m",
        "max_heating_rate_k_per_h",
    ]
    quantiles = (0.05, 0.25, 0.5, 0.75, 0.95)
    inputs = np.column_stack(
        [
            values,
            [slopes[int(k)] for k in terrain_index],
            [
                100 * float(members[int(k)]["series"].visible_fraction[spin:].mean())
                for k in terrain_index
            ],
        ]
    )
    input_names = [r.field for r in DEFAULT_RANGES] + [
        "terrain_slope_deg",
        "terrain_mean_illumination_pct",
    ]
    out["uncertainty"] = {
        "members": n,
        "terrain_members": len(members) - 1,
        "cpu_s": uq_cpu,
        "quantile_levels": quantiles,
        "terrain_slope_deg_range": [float(min(slopes)), float(max(slopes))],
        "outputs": {
            name: {
                "nominal": float(nom),
                "quantiles": np.quantile(outputs[:, j], quantiles).tolist(),
                "quantiles_first_64": np.quantile(outputs[:64, j], quantiles).tolist(),
                "quantiles_first_32": np.quantile(outputs[:32, j], quantiles).tolist(),
                "rank_correlation": dict(
                    zip(input_names, np.round(rank_correlation(inputs, outputs[:, j]), 2).tolist())
                ),
            }
            for j, (name, nom) in enumerate(
                zip(
                    names,
                    [
                        summary.surface_min_k,
                        summary.surface_max_k,
                        summary.surface_mean_k,
                        summary.lunation_amplitude_mean_k,
                        float(result.at_depth(0.3).mean()),
                        float(result.at_depth(1.0).mean()),
                        summary.max_heating_rate_k_per_h,
                    ],
                )
            )
        },
    }
    for name, block in out["uncertainty"]["outputs"].items():
        print(
            f"{name:>26s}: nominal {block['nominal']:7.2f} | 5-50-95 % {block['quantiles'][0]:7.2f} {block['quantiles'][2]:7.2f} {block['quantiles'][4]:7.2f}"
        )

    (ROOT / "docs" / "convergence" / "thermal.json").write_text(json.dumps(out, indent=1) + "\n")

    # ------------------------------------------------------------------------ figure
    days = result.time_s / 86400.0
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True, constrained_layout=True)
    axes[0].plot(days, result.absorbed_flux_w_m2, lw=0.6, color="tab:orange")
    axes[0].set_ylabel("absorbed solar flux (W m$^{-2}$)")
    axes[0].set_title(
        f"Shackleton rim test point ({SITE[0]:.1f}, {SITE[1]:.1f}) m, slope {np.degrees(slope):.1f} deg: thermal history, 2026"
    )
    low, high = np.quantile(surfaces, [0.05, 0.95], axis=0)
    axes[1].fill_between(
        days, low, high, color="tab:red", alpha=0.25, lw=0, label="surface, 5-95 % of the ensemble"
    )
    axes[1].plot(
        days,
        result.surface_temperature_k,
        lw=0.7,
        color="tab:red",
        label="surface, nominal parameters",
    )
    axes[1].set_ylabel("temperature (K)")
    axes[1].legend(loc="lower right", fontsize=8)
    for depth, color in (
        (0.05, "tab:purple"),
        (0.1, "tab:blue"),
        (0.3, "tab:green"),
        (1.0, "black"),
    ):
        axes[2].plot(days, result.at_depth(depth), lw=0.9, color=color, label=f"{depth:g} m")
    axes[2].set_ylabel("temperature at depth (K)")
    axes[2].set_xlabel("days since 2026-01-01 00:00 UTC")
    axes[2].legend(loc="lower right", fontsize=8, ncol=4)
    (ROOT / "docs" / "figures").mkdir(exist_ok=True)
    fig.savefig(ROOT / "docs" / "figures" / "thermal_site_2026.png", dpi=130)
    return 0


def replace_settings(settings: ThermalSettings, **changes) -> ThermalSettings:
    return settings.model_copy(update=changes)


if __name__ == "__main__":
    raise SystemExit(main())
