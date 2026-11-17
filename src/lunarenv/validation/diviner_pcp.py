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
from lunarenv.illumination.products import illumination_series, illumination_statistics
from lunarenv.illumination.settings import SITE04_NEAR, IlluminationSettings, load_site_stack
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack, coarsen
from lunarenv.thermal import terrain_radiation as tr
from lunarenv.thermal.albedo_laws import BURGER2024, FOOTE2020
from lunarenv.thermal.model import ThermalSettings, run_thermal
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.roughness import rough_radiance, sun_vectors
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


def _init_worker(
    data_root: str,
    horizon_coarsening: int = 1,
    near_id: str = SITE04_NEAR,
    laws: tuple[str, ...] | None = None,
) -> None:
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
    stack = load_site_stack(store, settings, near_id=near_id)
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
        laws=tuple(LAWS) if laws is None else laws,
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
    for name in w["laws"]:
        props = LAWS[name]
        flux = absorbed_solar_flux(series, normal, props, eclipse_factor=factor)
        result = run_thermal(flux, 3600.0, spin, props, ThermalSettings())
        radiance = props.emissivity * result.surface_temperature_k[summer] ** 4
        out[name] = np.bincount(index, weights=radiance, minlength=N_BINS)
        out[name + "_season"] = np.bincount(seasonal_index, weights=radiance, minlength=n_seasonal)
    return out


#: Laws for which TH-V3d adds the diffuse terms (pre-registered 2026-11-16).
DIFFUSE_LAWS = ("hayne", "foote2020")


def _receiver(xy: tuple[float, float]) -> dict[str, np.ndarray]:
    """TH-V3d: per-bin change of the summer radiance from terrain radiation and Earthshine
    at one receiver, per law, and from the low-temperature conductivity (Hayne2017 law)."""
    from lunarenv.runner import _diffuse_absorbed

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
    site = tr.site_radiation(
        w["stack"],
        x,
        y,
        {"ground": 0.0},
        w["et"],
        w["sun"],
        w["r_sun"],
        eclipse_factor=factor,
        n_strips=s.n_strips,
        limb_u=s.limb_darkening,
        mode=s.disc_mode_code,
    )
    relevant = site.relevant([("ground", normal)])
    earth = tr.earthshine(horizons, w["et"], w["earth"], w["sun"], w["r_earth"], projection=PROJ)
    e_vis, e_ir, e_inc = earth.on(normal)
    zeros = np.zeros_like(e_vis)
    spin, summer, bins = w["spin"], w["summer"], w["bins"]
    index = bins[summer]

    def radiance(flux: np.ndarray, props: RegolithProperties) -> np.ndarray:
        result = run_thermal(flux, 3600.0, spin, props, ThermalSettings())
        rad = props.emissivity * result.surface_temperature_k[summer] ** 4
        return np.bincount(index, weights=rad, minlength=N_BINS)

    out: dict[str, np.ndarray] = {"count": np.bincount(index, minlength=N_BINS).astype(float)}
    for name in DIFFUSE_LAWS:
        props = LAWS[name]
        site.compute(
            name,
            props,
            relevant,
            step_s=3600.0,
            spinup_samples=spin,
            settings=ThermalSettings(),
            workers=1,
        )
        ts, ti = site.irradiance("ground", name, normal)
        terrain = _diffuse_absorbed(
            {
                "terrain_scattered": ts,
                "terrain_infrared": ti,
                "earthshine_reflected": zeros,
                "earthshine_infrared": zeros,
            },
            None,
            props,
        )
        earthshine = _diffuse_absorbed(
            {
                "terrain_scattered": zeros,
                "terrain_infrared": zeros,
                "earthshine_reflected": e_vis,
                "earthshine_infrared": e_ir,
            },
            e_inc,
            props,
        )
        direct = absorbed_solar_flux(series, normal, props, eclipse_factor=factor)
        base = radiance(direct, props)
        with_diffuse = radiance(direct + terrain + earthshine, props)
        out[name + "_change"] = with_diffuse - base
        out[name + "_terrain_abs"] = np.array([terrain[spin:][summer].mean()])
        out[name + "_earth_abs"] = np.array([earthshine[spin:][summer].mean()])
        if name == "hayne":
            low_t = props.model_copy(update={"conductivity_model": "woods_robinson2019"})
            out["hayne_lowt_change"] = radiance(direct + terrain + earthshine, low_t) - with_diffuse
    out["left_out_bound"] = np.array([site.left_out_bound_w_m2("ground", "hayne", normal)])
    return out


def forward_diffuse(
    store: DataStore, cells: list[tuple[float, float]], spacing_m: float = 120.0
) -> list[list[dict[str, np.ndarray]]]:
    """TH-V3d receivers (every ``spacing_m`` in each cell): per-receiver results, per cell."""
    tasks: list[tuple[int, tuple[float, float]]] = []
    for c, (cx, cy) in enumerate(cells):
        xs, ys = cell_pixels(cx, cy, spacing_m)
        tasks += [(c, (float(x), float(y))) for x, y in zip(xs, ys, strict=True)]
    workers = max(1, min(8, (os.cpu_count() or 2) - 1))
    out: list[list[dict[str, np.ndarray]]] = [[] for _ in cells]
    with ProcessPoolExecutor(
        workers,
        mp_context=get_context("spawn"),
        initializer=_init_worker,
        initargs=(str(store.root), 1),
    ) as pool:
        for (c, _), result in zip(
            tasks, pool.map(_receiver, [t[1] for t in tasks], chunksize=1), strict=True
        ):
            out[c].append(result)
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
    near_id: str = SITE04_NEAR,
    laws: tuple[str, ...] | None = None,
) -> list[CellSums]:
    """Sum the per-pixel bin statistics over the pixels of each cell.

    ``horizon_coarsening`` > 1 computes the horizons on the site DEM block-averaged by that
    factor (the surface normal stays on the 5 m DEM): a sensitivity test, not TH-V3.
    ``laws`` limits the albedo laws computed (default: all of ``LAWS``).
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
        initargs=(str(store.root), horizon_coarsening, near_id, laws),
    ) as pool:
        for (c, _), result in zip(
            tasks, pool.map(_pixel, [t[1] for t in tasks], chunksize=8), strict=True
        ):
            counts[c] += 1
            for key, value in result.items():
                totals[c][key] = totals[c].get(key, 0.0) + value
    return [CellSums(n, t) for n, t in zip(counts, totals, strict=True)]


# ------------------------------------------------------- TH-V3e: block selection
#: Sites of TH-V3e (pre-registered 2026-11-17): label, site DEM.
TH_V3E_SITES = (
    ("Site01", "pgda78.site01.surf"),
    ("Site07", "pgda78.site07.surf"),
    ("Site11", "pgda78.site11.surf"),
    ("Site23", "pgda78.site23.surf"),
)
#: Cell indices (x / 240 m, y / 240 m) of the TH-V3 block, excluded from the selection.
TH_V3_CELLS = frozenset((i, j) for i in (-13, -12, -11) for j in (-39, -38, -37))
CELL_M = 240.0
SELECTION_MARGIN_M = 1000.0


def summer_axis(store: DataStore) -> tuple[list[np.ndarray], list[np.ndarray], float]:
    """Hourly epochs and Sun positions of each southern summer of the TH-V3 period.

    Each summer is a contiguous run of samples with subsolar latitude < 0, so it has a
    uniform step. Returns (epochs per summer, Sun positions per summer, Sun radius).
    """
    utc = np.arange(PERIOD_START, PERIOD_STOP, np.timedelta64(1, "h"))
    with spice_session(store) as eph:
        et = eph.utc_to_et(utc)
        sun = eph.position_from_center("SUN", et)
        r_sun = eph.body_radius_m("SUN")
    summer = sun[:, 2] < 0.0
    edges = np.flatnonzero(np.diff(np.concatenate([[0], summer.astype(np.int8), [0]])))
    runs = list(zip(edges[::2], edges[1::2], strict=True))
    return [et[a:b] for a, b in runs], [sun[a:b] for a, b in runs], r_sun


def select_block(
    store: DataStore, near_id: str, exclude: frozenset[tuple[int, int]] = TH_V3_CELLS
) -> dict[str, object]:
    """The 3 x 3 Diviner block of a site tile with the highest summer visible fraction.

    Pre-registered rule (docs/review/preregistration-th-v3e-2026-11-17.md): candidate
    cells are the 240 m cells whose centres lie at least 1 km inside the tile; each cell
    centre gets its ground-level horizon and mean visible solar disc over the summer
    samples; the block with the highest mean of its nine values wins, blocks touching
    ``exclude`` are skipped, ties go to the first block in row-major order (largest y
    first, then smallest x). Uses model geometry only.
    """
    settings = IlluminationSettings()
    stack = load_site_stack(store, settings, near_id=near_id)
    bounds = stack.levels[0].dem.grid.bounds
    i_range = np.arange(
        int(np.ceil((bounds.x_min + SELECTION_MARGIN_M) / CELL_M)),
        int(np.floor((bounds.x_max - SELECTION_MARGIN_M) / CELL_M)) + 1,
    )
    j_range = np.arange(
        int(np.floor((bounds.y_max - SELECTION_MARGIN_M) / CELL_M)),
        int(np.ceil((bounds.y_min + SELECTION_MARGIN_M) / CELL_M)) - 1,
        -1,
    )  # rows from largest y down
    jj, ii = np.meshgrid(j_range, i_range, indexing="ij")
    horizons = compute_horizons(
        stack,
        (ii * CELL_M).ravel().astype(np.float64),
        (jj * CELL_M).ravel().astype(np.float64),
        n_azimuth=settings.n_azimuth,
        diagnostics=False,
    )
    ets, suns, r_sun = summer_axis(store)
    total = np.zeros(ii.size)
    count = 0
    for et, sun in zip(ets, suns, strict=True):
        stats = illumination_statistics(
            horizons,
            et,
            sun,
            r_sun,
            n_strips=settings.n_strips,
            limb_u=settings.limb_darkening,
            mode=settings.disc_mode_code,
        )
        total += stats.mean_visible_fraction * et.size
        count += et.size
    visible = (total / count).reshape(ii.shape)
    best: tuple[float, int, int] | None = None
    for r in range(1, visible.shape[0] - 1):
        for c in range(1, visible.shape[1] - 1):
            cells = {
                (int(ii[r + a, c + b]), int(jj[r + a, c + b]))
                for a in (-1, 0, 1)
                for b in (-1, 0, 1)
            }
            if cells & exclude:
                continue
            value = float(visible[r - 1 : r + 2, c - 1 : c + 2].mean())
            if best is None or value > best[0]:
                best = (value, r, c)
    if best is None:
        raise ValueError(f"no admissible block in {near_id}")
    value, r, c = best
    centre = (int(ii[r, c]), int(jj[r, c]))
    lon, lat = PROJ.inverse(centre[0] * CELL_M, centre[1] * CELL_M)
    return {
        "near_id": near_id,
        "centre_cell": list(centre),
        "centre_m": [centre[0] * CELL_M, centre[1] * CELL_M],
        "latitude_deg": float(np.degrees(lat)),
        "east_longitude_deg": float(np.degrees(lon) % 360.0),
        "block_mean_visible_fraction": value,
        "candidate_cells": int(ii.size),
        "summer_samples": count,
    }


# ------------------------------------------------------- TH-V3f: sub-DEM roughness
#: Laws and RMS slopes of TH-V3f (pre-registered 2026-11-17).
ROUGH_LAWS = ("hayne", "foote2020")
ROUGH_SLOPES_DEG = (20.0, 35.0)


def rough_key(law: str, rms_deg: float) -> str:
    return f"{law}_rough{rms_deg:g}"


def _rough_receiver(task: tuple[float, float, int, int, bool]) -> dict[str, np.ndarray]:
    """TH-V3f: per-bin summer radiance sums of the smooth DEM facet and of the rough
    surface at one receiver, per law and RMS slope. ``task`` = (x, y, slope nodes,
    azimuths, Hayne2017 at the first slope only)."""
    w = _W
    x, y, n_nodes, n_azimuth, single = task
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
    irradiance = series.normal_irradiance_w_m2 * factor
    sun = sun_vectors(series.elevation_rad, series.azimuth_rad)
    spin, summer, bins = w["spin"], w["summer"], w["bins"]
    index = bins[summer]
    out: dict[str, np.ndarray] = {
        "count": np.bincount(index, minlength=N_BINS).astype(float),
        "xy": np.array([x, y]),
    }
    laws = ("hayne",) if single else ROUGH_LAWS
    slopes = ROUGH_SLOPES_DEG[:1] if single else ROUGH_SLOPES_DEG
    for name in laws:
        result = rough_radiance(
            irradiance,
            sun,
            normal,
            LAWS[name],
            slopes,
            3600.0,
            spin,
            ThermalSettings(),
            n_nodes=n_nodes,
            n_azimuth=n_azimuth,
        )
        out[name] = np.bincount(index, weights=result.smooth[summer], minlength=N_BINS)
        for rms, rad in result.rough.items():
            out[rough_key(name, rms)] = np.bincount(index, weights=rad[summer], minlength=N_BINS)
    return out


def forward_rough(
    store: DataStore,
    cells: list[tuple[float, float]],
    near_id: str = SITE04_NEAR,
    spacing_m: float = 120.0,
    convergence: tuple[float, float] | None = None,
) -> tuple[list[list[dict[str, np.ndarray]]], dict[str, np.ndarray] | None]:
    """TH-V3f receivers (every ``spacing_m`` in each cell, 3 slope nodes x 12 azimuths):
    per-receiver results per cell, and, if ``convergence`` gives a point, its result with
    5 nodes x 24 azimuths (Hayne2017, first RMS slope)."""
    tasks: list[tuple[int, tuple[float, float, int, int, bool]]] = []
    for c, (cx, cy) in enumerate(cells):
        xs, ys = cell_pixels(cx, cy, spacing_m)
        tasks += [(c, (float(x), float(y), 3, 12, False)) for x, y in zip(xs, ys, strict=True)]
    if convergence is not None:
        tasks.append((-1, (float(convergence[0]), float(convergence[1]), 5, 24, True)))
    workers = max(1, min(8, (os.cpu_count() or 2) - 1))
    out: list[list[dict[str, np.ndarray]]] = [[] for _ in cells]
    fine: dict[str, np.ndarray] | None = None
    with ProcessPoolExecutor(
        workers,
        mp_context=get_context("spawn"),
        initializer=_init_worker,
        initargs=(str(store.root), 1, near_id),
    ) as pool:
        # the long convergence task first, so it does not finish last on its own
        order = sorted(range(len(tasks)), key=lambda i: tasks[i][0] >= 0)
        results = pool.map(_rough_receiver, [tasks[i][1] for i in order], chunksize=1)
        for i, result in zip(order, results, strict=True):
            c = tasks[i][0]
            if c < 0:
                fine = result
            else:
                out[c].append(result)
    return out, fine
