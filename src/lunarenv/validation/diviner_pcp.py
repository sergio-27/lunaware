"""Forward model of the Diviner summer Polar Cumulative Product at the Shackleton rim.

Pre-registered in ``docs/review/preregistration-th-v3-2026-11-15.md``. For every 5 m
DEM pixel sampled inside the 3 x 3 Diviner cells of 240 m around the test point, the
engine computes the horizon at ground level, the direct solar flux on the DEM facet
(with Earth occultation) and the surface temperature, hourly from 2008-07-05 (one-year
spin-up) to 2019-02-18. Samples in the southern summer (subsolar latitude < 0) are
grouped by subsolar longitude in the 3.75 deg bins of the product, and the emitted
radiance ``eps * T^4`` is averaged over pixels and samples, as Diviner averages
radiances before converting them to a bolometric brightness temperature.

The pixels are independent, so they are spread over worker processes.
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from importlib import resources
from multiprocessing import get_context

import numba
import numpy as np

from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import spice_session
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.illumination.eclipse import earth_occultation_factor
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_series
from lunarenv.illumination.settings import IlluminationSettings, load_site_stack
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack, coarsen
from lunarenv.thermal.albedo_laws import BURGER2024, FOOTE2020
from lunarenv.thermal.model import ThermalSettings, run_thermal
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.surface import absorbed_solar_flux, surface_normal

PROJ = SouthPolarStereographic()
REFERENCE = "diviner_pcp_shackleton.json"
N_BINS = 96
BIN_DEG = 3.75
SPINUP_START = np.datetime64("2008-07-05T00:00")
PERIOD_START = np.datetime64("2009-07-05T00:00")
PERIOD_STOP = np.datetime64("2019-02-18T00:00")
#: Albedo laws compared: the engine default (Hayne2017 Table A1), a constant 0.12 and the
#: two bounds of the albedo envelope (lunarenv.thermal.albedo_laws).
LAWS: dict[str, RegolithProperties] = {
    "hayne": RegolithProperties(),
    "constant": RegolithProperties(albedo_a=0.0, albedo_b=0.0),
    # envelope bounds, added for TH-V3c (pre-registered 2026-11-16)
    "burger2024": BURGER2024.apply(RegolithProperties()),
    "foote2020": FOOTE2020.apply(RegolithProperties()),
}


def load_reference() -> dict:  # type: ignore[type-arg]
    path = resources.files("lunarenv.validation") / "reference" / REFERENCE
    return json.loads(path.read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def cell_pixels(
    centre_x: float, centre_y: float, spacing_m: float
) -> tuple[np.ndarray, np.ndarray]:
    """Pixel centres every ``spacing_m`` in a 240 m cell, half a spacing in from its edges."""
    offsets = np.arange(-120.0 + spacing_m / 2, 120.0, spacing_m)
    xx, yy = np.meshgrid(centre_x + offsets, centre_y + offsets)
    return xx.ravel(), yy.ravel()


# ------------------------------------------------------------------- worker processes
_W: dict = {}  # type: ignore[type-arg]


def _horizon_stack(
    stack: TerrainStack, factor: int, settings: IlluminationSettings
) -> TerrainStack:
    """The stack with its near (5 m) level block-averaged by ``factor``, for horizons only."""
    if factor == 1:
        return stack
    near = stack.levels[0]
    dem = coarsen(near.dem, factor)
    level = TerrainLevel(
        dem, near.start_m, near.end_m, settings.step_fraction * dem.grid.pixel_size_m
    )
    return TerrainStack((level, *stack.levels[1:]))


def _init_worker(data_root: str, horizon_coarsening: int = 1) -> None:
    # one numba thread per worker; set at run time (changing NUMBA_NUM_THREADS after numba
    # is imported fails as soon as a kernel is compiled in the worker)
    numba.set_num_threads(1)
    store = DataStore(data_root)
    settings = IlluminationSettings()
    utc = np.arange(SPINUP_START, PERIOD_STOP, np.timedelta64(1, "h"))
    with spice_session(store) as eph:
        et = eph.utc_to_et(utc)
        sun = eph.position_from_center("SUN", et)
        earth = eph.position_from_center("EARTH", et)
        r_sun, r_earth = eph.body_radius_m("SUN"), eph.body_radius_m("EARTH")
    spin = int(np.searchsorted(utc, PERIOD_START))
    period_sun = sun[spin:]
    latitude = np.arcsin(period_sun[:, 2] / np.linalg.norm(period_sun, axis=1))
    longitude = np.degrees(np.arctan2(period_sun[:, 1], period_sun[:, 0])) % 360.0
    summer = latitude < 0.0
    bins = np.minimum((longitude // BIN_DEG).astype(np.int64), N_BINS - 1)
    # label each southern summer (consecutive run of samples with subsolar latitude < 0)
    starts = np.flatnonzero(summer & ~np.concatenate([[False], summer[:-1]]))
    season = np.cumsum(np.isin(np.arange(summer.size), starts)) - 1
    stack = load_site_stack(store, settings)
    _W.update(
        store=store,
        settings=settings,
        stack=stack,
        horizon_stack=_horizon_stack(stack, horizon_coarsening, settings),
        et=et,
        sun=sun,
        earth=earth,
        r_sun=r_sun,
        r_earth=r_earth,
        spin=spin,
        summer=summer,
        bins=bins,
        season=season,
        n_seasons=int(season[summer].max()) + 1,
    )


def _pixel(xy: tuple[float, float]) -> dict[str, np.ndarray]:
    """Per-bin sums for one pixel: radiance per law (all and per summer), visible fraction."""
    w = _W
    x, y = xy
    s: IlluminationSettings = w["settings"]
    horizons = compute_horizons(
        w["horizon_stack"], np.array([x]), np.array([y]), n_azimuth=s.n_azimuth, diagnostics=False
    )
    series = illumination_series(
        horizons,
        w["et"],
        w["sun"],
        w["r_sun"],
        n_strips=s.n_strips,
        limb_u=s.limb_darkening,
        mode=s.disc_mode_code,
    )
    position = PROJ.to_cartesian(x, y, float(horizons.observer_height_m[0]))
    factor = earth_occultation_factor(w["sun"], w["earth"], position, w["r_sun"], w["r_earth"])
    normal = surface_normal(w["stack"].levels[0].dem, x, y)
    spin, summer, bins, season = w["spin"], w["summer"], w["bins"], w["season"]
    index = bins[summer]
    out: dict[str, np.ndarray] = {
        "visible": np.bincount(
            index, weights=series.visible_fraction[spin:][summer], minlength=N_BINS
        ),
        "count": np.bincount(index, minlength=N_BINS).astype(np.float64),
    }
    seasonal_index = season[summer] * N_BINS + index
    n_seasonal = w["n_seasons"] * N_BINS
    out["season_count"] = np.bincount(seasonal_index, minlength=n_seasonal).astype(np.float64)
    for name, props in LAWS.items():
        flux = absorbed_solar_flux(series, normal, props, eclipse_factor=factor)
        result = run_thermal(flux, 3600.0, spin, props, ThermalSettings())
        radiance = props.emissivity * result.surface_temperature_k[summer] ** 4
        out[name] = np.bincount(index, weights=radiance, minlength=N_BINS)
        out[name + "_season"] = np.bincount(seasonal_index, weights=radiance, minlength=n_seasonal)
    return out


@dataclass(frozen=True, slots=True)
class CellSums:
    """Sums over the pixels of one cell (arrays per subsolar-longitude bin)."""

    n_pixels: int
    sums: dict[str, np.ndarray]


def forward(
    store: DataStore,
    cells: list[tuple[float, float]],
    spacing_m: float,
    horizon_coarsening: int = 1,
) -> list[CellSums]:
    """Sum the per-pixel bin statistics over the pixels of each cell.

    ``horizon_coarsening`` > 1 computes the horizons on the site DEM block-averaged by that
    factor (the surface normal stays on the 5 m DEM): a sensitivity test, not TH-V3.
    """
    tasks: list[tuple[int, tuple[float, float]]] = []
    for c, (cx, cy) in enumerate(cells):
        xs, ys = cell_pixels(cx, cy, spacing_m)
        tasks += [(c, (float(x), float(y))) for x, y in zip(xs, ys, strict=True)]
    workers = max(1, min(8, (os.cpu_count() or 2) - 1))
    totals: list[dict[str, np.ndarray]] = [{} for _ in cells]
    counts = [0] * len(cells)
    with ProcessPoolExecutor(
        workers,
        mp_context=get_context("spawn"),
        initializer=_init_worker,
        initargs=(str(store.root), horizon_coarsening),
    ) as pool:
        for (c, _), result in zip(
            tasks, pool.map(_pixel, [t[1] for t in tasks], chunksize=8), strict=True
        ):
            counts[c] += 1
            for key, value in result.items():
                totals[c][key] = totals[c].get(key, 0.0) + value
    return [CellSums(n, t) for n, t in zip(counts, totals, strict=True)]
