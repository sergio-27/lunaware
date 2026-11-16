"""Single-bounce terrain radiation and Earthshine on a receiver (review change 3).

Energy reaching a surface element besides direct sunlight:

* **Terrain-scattered sunlight** and **terrain thermal emission**, single bounce: each
  visible terrain element is lit by the Sun only (no light it receives from other
  terrain), scatters ``A(theta) E cos(theta)`` and emits ``epsilon sigma T^4``, both
  Lambertian.
* **Earthshine**: sunlight reflected by the Earth and the Earth's own infrared, with the
  peak values of NASA HLS-UG-001 section 11.1.

Method (pre-registered in docs/review/preregistration-th-a3-2026-11-16.md):

1. *Visible terrain.* Rays leave the receiver over the terrain stack with the geometry
   and sample plan of the horizon tracer (:mod:`lunarenv.illumination.horizon`). A sample
   is visible when its elevation exceeds that of every closer sample on the ray; it then
   fills the elevation interval from that previous maximum up to its own elevation (the
   first sample fills from -90 degrees).
2. *Patches.* Visible samples are grouped by azimuth sector and distance band around the
   site. Each patch is represented by its visible sample of largest solid angle; that
   point gets its own horizon, DEM normal, direct sunlight and 1-D temperature, with the
   same models as the receiver.
3. *Receiver irradiance.* Each sample contributes its patch's exitance M times its
   cosine-weighted solid angle over pi, integrated exactly in elevation (clipped at the
   receiver plane) and at the ray azimuth over the azimuth step. Sky and terrain together
   give a weight of 1 for a full hemisphere.

What is left out: multiple scattering between terrain elements, roughness below the DEM
resolution, non-Lambertian scattering and emission, shadowing of the ground by the asset
itself, and the variation of sunlight within a patch (one representative point).
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from multiprocessing import get_context

import numpy as np
from numpy.typing import ArrayLike, NDArray

from lunarenv.core.constants import SOLAR_CONSTANT_W_M2, STEFAN_BOLTZMANN_W_M2_K4
from lunarenv.geometry.projection import FloatArray, SouthPolarStereographic
from lunarenv.illumination.disc import EDDINGTON_LIMB_DARKENING, PROFILE
from lunarenv.illumination.horizon import HorizonSet, compute_horizons
from lunarenv.illumination.products import IlluminationSeries, illumination_series
from lunarenv.illumination.terrain import TerrainStack
from lunarenv.thermal.model import ThermalSettings, spin_up
from lunarenv.thermal.properties import RegolithProperties, albedo
from lunarenv.thermal.solver import solve_conduction
from lunarenv.thermal.surface import incidence_cosine, surface_normal

#: Default patch layout for scenarios: 36 sectors of 10 degrees and distance bands
#: growing by a factor sqrt(2) from 25 m (the last band runs to the end of the stack).
#: Chosen on the analytic bowl (TH-A3, informational): radial resolution matters most.
DEFAULT_SECTORS = 36
DEFAULT_BAND_EDGES_M: tuple[float, ...] = (
    0.0,
    *(25.0 * 2.0 ** (k / 2.0) for k in range(31)),
    float("inf"),
)
#: Rays from the receiver and azimuths of each patch's horizon.
DEFAULT_RECEIVER_RAYS = 360
DEFAULT_PATCH_AZIMUTHS = 360
#: Representative points per patch (stratified by solid angle along distance).
DEFAULT_POINTS_PER_PATCH = 3
_RAY_CHUNK = 24

#: Earthshine at full Earth (illuminated fraction 1), normal incidence: reflected sunlight
#: and Earth infrared. HLS-UG-001 (NASA, 4 Jan 2021) section 11.1, citing Glenar et al.
#: (2019): "the infrared contribution is estimated to be 56 mW/m2 while the reflected
#: solar is 100 mW/m2" at peak illumination.
EARTHSHINE_REFLECTED_W_M2 = 0.100
EARTHSHINE_INFRARED_W_M2 = 0.056
#: Earth-Moon distance the guidebook values are taken to refer to (an assumption: the
#: guidebook does not state it). Earthshine scales as the inverse square of the distance.
EARTH_MOON_MEAN_DISTANCE_M = 3.844e8


# ------------------------------------------------------------------- visible terrain
@dataclass(frozen=True, slots=True)
class VisibleTerrain:
    """Terrain samples visible from one receiver position.

    Sample ``i`` lies on ray ``ray[i]`` (grid azimuth ``2 pi ray / n_rays``) at ground
    distance ``distance_m[i]`` and fills the elevation interval
    ``[elevation_low_rad[i], elevation_high_rad[i]]`` of the receiver's view.
    """

    x_m: float
    y_m: float
    observer_height_m: float
    n_rays: int
    ray: NDArray[np.int64]
    distance_m: FloatArray
    sample_x_m: FloatArray
    sample_y_m: FloatArray
    sample_height_m: FloatArray
    level: NDArray[np.int64]
    elevation_low_rad: FloatArray
    elevation_high_rad: FloatArray

    @property
    def azimuth_rad(self) -> FloatArray:
        return np.asarray(2.0 * np.pi * self.ray / self.n_rays)

    @property
    def solid_angle_sr(self) -> FloatArray:
        step = 2.0 * np.pi / self.n_rays
        return np.asarray(step * (np.sin(self.elevation_high_rad) - np.sin(self.elevation_low_rad)))


def _sample_levels(
    stack: TerrainStack, x: FloatArray, y: FloatArray, start_level: NDArray[np.int64] | int
) -> tuple[FloatArray, NDArray[np.int64]]:
    """Bilinear height from the first level at or after ``start_level`` covering each
    point (as the horizon tracer does); NaN and level -1 where no level covers it."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    start = np.broadcast_to(np.asarray(start_level), x.shape)
    height = np.full(x.shape, np.nan)
    level_of = np.full(x.shape, -1, dtype=np.int64)
    for index, level in enumerate(stack.levels):
        todo = np.isnan(height) & (start <= index)
        if not np.any(todo):
            continue
        grid = level.dem.grid
        col = (x[todo] - grid.x_min) / grid.pixel_size_m - 0.5
        row = (grid.y_max - y[todo]) / grid.pixel_size_m - 0.5
        n_rows, n_cols = grid.shape
        inside = (col >= 0.0) & (row >= 0.0) & (col <= n_cols - 1) & (row <= n_rows - 1)
        if not np.any(inside):
            continue
        col, row = col[inside], row[inside]
        c0 = np.minimum(col.astype(np.int64), n_cols - 2)
        r0 = np.minimum(row.astype(np.int64), n_rows - 2)
        fc, fr = col - c0, row - r0
        z = level.dem.heights_m
        top = z[r0, c0] * (1.0 - fc) + z[r0, c0 + 1] * fc
        bottom = z[r0 + 1, c0] * (1.0 - fc) + z[r0 + 1, c0 + 1] * fc
        target = np.flatnonzero(todo)[inside]
        height.flat[target] = top * (1.0 - fr) + bottom * fr
        level_of.flat[target] = index
    return height, level_of


def visible_terrain(
    stack: TerrainStack,
    x_m: float,
    y_m: float,
    *,
    height_above_surface_m: float = 0.0,
    n_rays: int = DEFAULT_RECEIVER_RAYS,
) -> VisibleTerrain:
    """Terrain samples visible from a receiver ``height_above_surface_m`` above the
    surface at map point ``(x_m, y_m)``."""
    if n_rays < 8:
        raise ValueError("need at least 8 rays")
    projection = stack.levels[0].dem.grid.projection
    radius = stack.radius_m
    surface, _ = _sample_levels(stack, np.array([x_m]), np.array([y_m]), 0)
    if not np.isfinite(surface[0]):
        raise ValueError("receiver lies outside the terrain stack")
    r_obs = radius + float(surface[0]) + height_above_surface_m
    distance, plan_level = stack.sample_plan()
    c, s = np.cos(distance / radius), np.sin(distance / radius)
    e_x, e_y, up = (np.asarray(v) for v in projection.local_basis(x_m, y_m))
    parts: dict[str, list[NDArray[np.generic]]] = {
        k: [] for k in ("ray", "j", "x", "y", "h", "lv", "lo", "hi")
    }
    for first in range(0, n_rays, _RAY_CHUNK):
        rays = np.arange(first, min(first + _RAY_CHUNK, n_rays))
        azimuth = 2.0 * np.pi * rays / n_rays
        t = np.sin(azimuth)[:, None] * e_x + np.cos(azimuth)[:, None] * e_y  # (k, 3)
        u = up[None, None, :] * c[None, :, None] + t[:, None, :] * s[None, :, None]
        scale = 2.0 * radius / (1.0 - u[..., 2])
        x = u[..., 1] * scale
        y = u[..., 0] * scale
        height, level = _sample_levels(stack, x, y, plan_level[None, :])
        # a ray ends where it first runs out of terrain data
        missing = np.isnan(height)
        ends = np.where(missing.any(axis=1), missing.argmax(axis=1), distance.size)
        valid = np.arange(distance.size)[None, :] < ends[:, None]
        r_terrain = radius + np.where(valid, height, 0.0)
        tangent = np.where(valid, (r_terrain * c - r_obs) / (r_terrain * s), -np.inf)
        previous = np.maximum.accumulate(
            np.concatenate([np.full((rays.size, 1), -np.inf), tangent[:, :-1]], axis=1), axis=1
        )
        seen = valid & (tangent > previous)
        k, j = np.nonzero(seen)
        parts["ray"].append(rays[k])
        parts["j"].append(j)
        parts["x"].append(x[k, j])
        parts["y"].append(y[k, j])
        parts["h"].append(height[k, j])
        parts["lv"].append(level[k, j])
        parts["lo"].append(np.arctan(previous[k, j]))  # arctan(-inf) = -pi/2
        parts["hi"].append(np.arctan(tangent[k, j]))
    cat = {k: np.concatenate(v) for k, v in parts.items()}
    return VisibleTerrain(
        x_m=float(x_m),
        y_m=float(y_m),
        observer_height_m=r_obs - radius,
        n_rays=n_rays,
        ray=cat["ray"].astype(np.int64),
        distance_m=distance[cat["j"]],
        sample_x_m=cat["x"].astype(np.float64),
        sample_y_m=cat["y"].astype(np.float64),
        sample_height_m=cat["h"].astype(np.float64),
        level=cat["lv"].astype(np.int64),
        elevation_low_rad=cat["lo"].astype(np.float64),
        elevation_high_rad=cat["hi"].astype(np.float64),
    )


def cosine_weights(view: VisibleTerrain, normal: ArrayLike) -> FloatArray:
    """Weight of each sample for a receiver with unit ``normal`` (receiver's local grid
    frame): (1/pi) * integral of max(0, n.d) cos(e) de dphi over the sample's interval.

    Irradiance from a Lambertian element of exitance M is M times this weight; a full
    hemisphere above the receiver plane has weight 1.
    """
    n = np.asarray(normal, dtype=np.float64)
    phi = view.azimuth_rad
    a = n[0] * np.sin(phi) + n[1] * np.cos(phi)  # horizontal part of n along the ray
    b = n[2]
    low, high = view.elevation_low_rad.copy(), view.elevation_high_rad.copy()
    # n.d = a cos(e) + b sin(e) changes sign at tan(e) = -a / b
    if b > 0.0:
        low = np.maximum(low, np.arctan(-a / b))
    elif b < 0.0:
        high = np.minimum(high, np.arctan(-a / b))
    else:
        high = np.where(a > 0.0, high, low)

    def primitive(e: FloatArray) -> FloatArray:
        return np.asarray(a * (e / 2.0 + np.sin(2.0 * e) / 4.0) + b * np.sin(e) ** 2 / 2.0)

    integral = np.where(high > low, primitive(high) - primitive(low), 0.0)
    return np.asarray(np.maximum(integral, 0.0) * (2.0 * np.pi / view.n_rays) / np.pi)


# --------------------------------------------------------------------------- patches
@dataclass(frozen=True, slots=True)
class Patches:
    """Terrain patches seen by receivers at one map point.

    Each patch is represented by up to ``points_per_patch`` points (arrays indexed by
    point, sorted by ``owner``, the patch of each point); the patch exitance is the mean
    over its points. ``solid_angle_sr`` is per patch: the largest over the receivers.
    """

    x_m: FloatArray
    y_m: FloatArray
    height_m: FloatArray
    level: NDArray[np.int64]
    owner: NDArray[np.int64]
    solid_angle_sr: FloatArray

    @property
    def size(self) -> int:
        """Number of patches."""
        return int(self.solid_angle_sr.size)

    @property
    def n_points(self) -> int:
        return int(self.x_m.size)

    def mean_over_points(self, values: FloatArray) -> FloatArray:
        """Per-patch mean of a per-point array (first axis)."""
        starts = np.flatnonzero(np.r_[True, self.owner[1:] != self.owner[:-1]])
        counts = np.diff(np.r_[starts, self.owner.size])
        sums = np.add.reduceat(np.asarray(values), starts, axis=0)
        shape = (-1,) + (1,) * (sums.ndim - 1)
        return np.asarray(sums / counts.reshape(shape))


def build_patches(
    views: Sequence[VisibleTerrain],
    *,
    n_sectors: int = DEFAULT_SECTORS,
    band_edges_m: Sequence[float] = DEFAULT_BAND_EDGES_M,
    points_per_patch: int = DEFAULT_POINTS_PER_PATCH,
) -> tuple[Patches, list[NDArray[np.int64]]]:
    """Group the visible samples of receivers at one map point into patches.

    The representative points of a patch are stratified: its samples (of all views,
    weighted by solid angle) are ordered by distance and azimuth, and the sample at the
    middle of each of ``points_per_patch`` equal shares of the weight is taken.
    Returns the patches and, for each view, the patch index of every sample.
    """
    if not views:
        raise ValueError("need at least one view")
    if any((v.x_m, v.y_m) != (views[0].x_m, views[0].y_m) for v in views):
        raise ValueError("sectors and bands are relative to one map point: views must share it")
    if points_per_patch < 1:
        raise ValueError("need at least one point per patch")
    edges = np.asarray(band_edges_m, dtype=np.float64)
    if edges[0] != 0.0 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("band edges must start at 0 and increase")
    n_bands = edges.size - 1
    keys = []
    for view in views:
        sector = np.minimum(
            (view.azimuth_rad / (2.0 * np.pi) * n_sectors).astype(np.int64), n_sectors - 1
        )
        band = np.clip(np.searchsorted(edges, view.distance_m, side="right") - 1, 0, n_bands - 1)
        keys.append(sector * n_bands + band)
    unique = np.unique(np.concatenate(keys))
    n = unique.size
    indices = [np.searchsorted(unique, key) for key in keys]
    omega = np.zeros(n)
    for view, index in zip(views, indices, strict=True):
        omega = np.maximum(omega, np.bincount(index, view.solid_angle_sr, minlength=n))

    def cat(name: str) -> NDArray[np.generic]:
        return np.concatenate([getattr(v, name) for v in views])

    owner = np.concatenate(indices)
    weight = np.concatenate([v.solid_angle_sr for v in views])
    distance, azimuth = cat("distance_m"), cat("azimuth_rad")
    order = np.lexsort((azimuth, distance, owner))
    owner, weight = owner[order], weight[order]
    cumulative = np.cumsum(weight)
    starts = np.flatnonzero(np.r_[True, owner[1:] != owner[:-1]])
    before = np.r_[0.0, cumulative][starts]
    total = np.add.reduceat(weight, starts)
    group = np.repeat(np.arange(starts.size), np.diff(np.r_[starts, owner.size]))
    # position of each sample's centre within its patch's weight, in [0, 1], plus the
    # patch number so that the sequence increases over all patches
    position = (cumulative - 0.5 * weight - before[group]) / np.where(total > 0, total, 1.0)[group]
    position = group + np.clip(position, 0.0, 1.0 - 1e-12)
    targets = (
        np.arange(starts.size)[:, None]
        + (np.arange(points_per_patch)[None, :] + 0.5) / points_per_patch
    ).ravel()
    pick = np.clip(np.searchsorted(position, targets), 0, owner.size - 1)
    # stay inside the patch (searchsorted may step into the next one)
    pick = np.clip(
        pick,
        starts[targets.astype(np.int64)],
        np.r_[starts[1:], owner.size][targets.astype(np.int64)] - 1,
    )
    pick = np.unique(pick)
    chosen = order[pick]
    return (
        Patches(
            cat("sample_x_m")[chosen].astype(np.float64),
            cat("sample_y_m")[chosen].astype(np.float64),
            cat("sample_height_m")[chosen].astype(np.float64),
            cat("level")[chosen].astype(np.int64),
            owner[pick],
            omega,
        ),
        indices,
    )


# -------------------------------------------------------------------- patch sunlight
@dataclass(frozen=True, slots=True)
class PatchSunlight:
    """Direct sunlight on every representative point over the whole time axis (spin-up
    included), shape ``(n_points, n_times)``."""

    patches: Patches
    cosine: FloatArray  # clipped at 0
    irradiance_w_m2: FloatArray  # normal irradiance times eclipse factor
    horizons: HorizonSet
    normals: FloatArray  # (n_points, 3), each in its point's local grid frame


def _patch_normal(stack: TerrainStack, x: float, y: float, level: int) -> FloatArray:
    for lv in stack.levels[max(level, 0) :]:
        try:
            return surface_normal(lv.dem, x, y)
        except ValueError:
            continue
    return np.array([0.0, 0.0, 1.0])


def patch_sunlight(
    stack: TerrainStack,
    patches: Patches,
    et_s: FloatArray,
    sun_position_m: FloatArray,
    sun_radius_m: float,
    *,
    eclipse_factor: ArrayLike | None = None,
    n_azimuth: int = DEFAULT_PATCH_AZIMUTHS,
    n_strips: int = 32,
    limb_u: float = EDDINGTON_LIMB_DARKENING,
    mode: int = PROFILE,
) -> PatchSunlight:
    """Horizon, normal and direct sunlight of every patch's representative point.

    The Earth-eclipse factor of the site is applied to every patch (the patches lie
    within a few hundred kilometres of the site).
    """
    # height from the finest level covering the point, as the tracer samples near terrain
    height, _ = _sample_levels(stack, patches.x_m, patches.y_m, 0)
    height = np.where(np.isfinite(height), height, patches.height_m)
    horizons = compute_horizons(
        stack,
        patches.x_m,
        patches.y_m,
        n_azimuth=n_azimuth,
        diagnostics=False,
        surface_height_m=height,
    )
    projection: SouthPolarStereographic = stack.levels[0].dem.grid.projection
    n_times = np.asarray(et_s).size
    cosine = np.empty((patches.n_points, n_times))
    irradiance = np.empty((patches.n_points, n_times))
    normals = np.empty((patches.n_points, 3))
    factor = None if eclipse_factor is None else np.asarray(eclipse_factor, dtype=np.float64)
    for i in range(patches.n_points):
        series: IlluminationSeries = illumination_series(
            horizons,
            et_s,
            sun_position_m,
            sun_radius_m,
            observer=i,
            projection=projection,
            n_strips=n_strips,
            limb_u=limb_u,
            mode=mode,
        )
        normal = _patch_normal(
            stack, float(patches.x_m[i]), float(patches.y_m[i]), int(patches.level[i])
        )
        normals[i] = normal
        cosine[i] = np.clip(
            incidence_cosine(normal, series.elevation_rad, series.azimuth_rad), 0.0, 1.0
        )
        e = series.normal_irradiance_w_m2
        irradiance[i] = e if factor is None else e * factor
    return PatchSunlight(patches, cosine, irradiance, horizons, normals)


# -------------------------------------------------------------------- patch exitance
@dataclass(frozen=True, slots=True)
class PatchExitance:
    """Exitance (W m-2) of every patch over the whole time axis (mean over its points):
    scattered sunlight and thermal emission, shape ``(n_patches, n_times)``."""

    reflected_w_m2: FloatArray
    emitted_w_m2: FloatArray


def _temperatures(
    fluxes: FloatArray,
    step_s: float,
    spinup_samples: int,
    properties: RegolithProperties,
    settings: ThermalSettings,
) -> FloatArray:
    """Surface temperature over the whole axis: periodic state of the spin-up cycle, then
    the whole series from it (the period part equals :func:`run_thermal`'s)."""
    grid = settings.grid()
    out = np.empty_like(fluxes)
    for i, flux in enumerate(fluxes):
        start = spin_up(flux[: spinup_samples + 1], step_s, properties, settings)
        result = solve_conduction(
            grid,
            properties,
            start.temperature_k,
            flux,
            step_s,
            scheme=settings.scheme,
            substeps=settings.substeps,
        )
        out[i] = result.surface_temperature_k
    return out


def _init_worker() -> None:
    import numba

    numba.set_num_threads(1)


def patch_exitance(
    sunlight: PatchSunlight,
    properties: RegolithProperties,
    *,
    step_s: float,
    spinup_samples: int,
    settings: ThermalSettings | None = None,
    equilibrium: bool = False,
    only: NDArray[np.bool_] | None = None,
    workers: int | None = None,
) -> PatchExitance:
    """Scattered and emitted exitance of each patch under one set of properties.

    ``equilibrium`` replaces the 1-D thermal model by instantaneous radiative equilibrium
    (emitted = absorbed), as in the analytic bowl case. ``only`` (per patch) restricts
    the computation to those patches; the others get zero exitance.
    """
    incident = sunlight.irradiance_w_m2 * sunlight.cosine
    if only is not None:
        incident = np.where(np.asarray(only)[sunlight.patches.owner][:, None], incident, 0.0)
    reflectance = albedo(np.arccos(sunlight.cosine), properties)
    reflected = reflectance * incident
    absorbed = incident - reflected
    patches = sunlight.patches
    if equilibrium:
        return PatchExitance(
            patches.mean_over_points(reflected), patches.mean_over_points(absorbed)
        )
    settings = settings or ThermalSettings()
    todo = np.arange(absorbed.shape[0])
    if only is not None:
        todo = todo[np.asarray(only)[sunlight.patches.owner]]
    n = todo.size
    temperature = np.zeros_like(absorbed)
    if workers is None:
        workers = max(1, min(8, (os.cpu_count() or 2) - 1))
    if n == 0:
        pass
    elif workers <= 1 or n < 2 * workers:
        temperature[todo] = _temperatures(
            absorbed[todo], step_s, spinup_samples, properties, settings
        )
    else:
        chunks = np.array_split(todo, workers * 3)
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=get_context("spawn"), initializer=_init_worker
        ) as pool:
            futures = [
                (
                    chunk,
                    pool.submit(
                        _temperatures, absorbed[chunk], step_s, spinup_samples, properties, settings
                    ),
                )
                for chunk in chunks
                if chunk.size
            ]
            for chunk, future in futures:
                temperature[chunk] = future.result()
    emitted = properties.emissivity * STEFAN_BOLTZMANN_W_M2_K4 * temperature**4
    if only is not None:
        emitted = np.where(np.asarray(only)[sunlight.patches.owner][:, None], emitted, 0.0)
    return PatchExitance(patches.mean_over_points(reflected), patches.mean_over_points(emitted))


def receiver_irradiance(
    view: VisibleTerrain,
    patch_index: NDArray[np.int64],
    exitance: PatchExitance,
    normal: ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    """Terrain irradiance (W m-2) on a receiver: scattered sunlight and infrared."""
    weights = np.bincount(
        patch_index, cosine_weights(view, normal), minlength=exitance.reflected_w_m2.shape[0]
    )
    return weights @ exitance.reflected_w_m2, weights @ exitance.emitted_w_m2


def diffuse_albedo(properties: RegolithProperties, n: int = 2000) -> float:
    """Albedo of the law under isotropic illumination: integral of A(theta) sin(2 theta)
    over [0, pi/2] (midpoint rule)."""
    theta = (np.arange(n) + 0.5) * (np.pi / 2.0) / n
    return float(np.sum(albedo(theta, properties) * np.sin(2.0 * theta)) * (np.pi / 2.0) / n)


# ----------------------------------------------------------------- site bookkeeping
#: Patches that cannot send more than this to any receiver of interest are not computed.
RELEVANCE_W_M2 = 1e-4
#: Largest exitance a patch can have: the solar irradiance at the smallest Sun distance
#: (DE440, 2024-2044: 0.9806 au).
MAX_EXITANCE_W_M2 = SOLAR_CONSTANT_W_M2 / 0.98**2


def patch_layout(refinement: int = 1) -> tuple[int, tuple[float, ...]]:
    """Sectors and band edges: the defaults, or ``refinement`` times as many of each."""
    if refinement < 1:
        raise ValueError("refinement must be at least 1")
    step = 0.5 / refinement
    edges = (0.0, *(25.0 * 2.0 ** (k * step) for k in range(30 * refinement + 1)), float("inf"))
    return DEFAULT_SECTORS * refinement, edges


@dataclass
class SiteRadiation:
    """Visible terrain, patches and patch exitance for receivers at one map point.

    ``views`` maps a receiver name (e.g. ground, asset) to what it sees. Exitance is
    computed per albedo case with :meth:`compute`, for the patches that matter.
    """

    views: dict[str, VisibleTerrain]
    indices: dict[str, NDArray[np.int64]]
    sunlight: PatchSunlight
    exitance: dict[str, PatchExitance] = field(default_factory=dict)
    computed: dict[str, NDArray[np.bool_]] = field(default_factory=dict)

    @property
    def patches(self) -> Patches:
        return self.sunlight.patches

    def weights(self, view: str, normal: ArrayLike) -> FloatArray:
        """Per-patch cosine weight for a receiver of the named view."""
        return np.bincount(
            self.indices[view],
            cosine_weights(self.views[view], normal),
            minlength=self.patches.size,
        )

    def relevant(self, receivers: Sequence[tuple[str, ArrayLike]]) -> NDArray[np.bool_]:
        """Patches that can send more than :data:`RELEVANCE_W_M2` to any receiver."""
        mask = np.zeros(self.patches.size, dtype=bool)
        for view, normal in receivers:
            mask |= self.weights(view, normal) * MAX_EXITANCE_W_M2 > RELEVANCE_W_M2
        return mask

    def compute(
        self,
        case: str,
        properties: RegolithProperties,
        only: NDArray[np.bool_],
        **kwargs: object,
    ) -> None:
        self.exitance[case] = patch_exitance(self.sunlight, properties, only=only, **kwargs)  # type: ignore[arg-type]
        self.computed[case] = np.asarray(only, dtype=bool)

    def irradiance(self, view: str, case: str, normal: ArrayLike) -> tuple[FloatArray, FloatArray]:
        """Scattered and infrared terrain irradiance over the whole time axis."""
        return receiver_irradiance(
            self.views[view], self.indices[view], self.exitance[case], normal
        )

    def left_out_bound_w_m2(self, view: str, case: str, normal: ArrayLike) -> float:
        """Upper bound of what the patches not computed could send to the receiver."""
        weights = self.weights(view, normal)
        return float(weights[~self.computed[case]].sum() * MAX_EXITANCE_W_M2)


def site_radiation(
    stack: TerrainStack,
    x_m: float,
    y_m: float,
    heights_above_surface_m: dict[str, float],
    et_s: FloatArray,
    sun_position_m: FloatArray,
    sun_radius_m: float,
    *,
    eclipse_factor: ArrayLike | None = None,
    refinement: int = 1,
    n_rays: int = DEFAULT_RECEIVER_RAYS,
    points_per_patch: int = DEFAULT_POINTS_PER_PATCH,
    n_strips: int = 32,
    limb_u: float = EDDINGTON_LIMB_DARKENING,
    mode: int = PROFILE,
) -> SiteRadiation:
    """Views of the named receivers, their common patches and the patches' sunlight."""
    views = {
        name: visible_terrain(stack, x_m, y_m, height_above_surface_m=h, n_rays=n_rays)
        for name, h in heights_above_surface_m.items()
    }
    sectors, edges = patch_layout(refinement)
    patches, indices = build_patches(
        list(views.values()),
        n_sectors=sectors,
        band_edges_m=edges,
        points_per_patch=points_per_patch,
    )
    light = patch_sunlight(
        stack,
        patches,
        et_s,
        sun_position_m,
        sun_radius_m,
        eclipse_factor=eclipse_factor,
        n_strips=n_strips,
        limb_u=limb_u,
        mode=mode,
    )
    return SiteRadiation(views, dict(zip(views, indices, strict=True)), light)


# ------------------------------------------------------------------------- Earthshine
@dataclass(frozen=True, slots=True)
class Earthshine:
    """Earthshine at a receiver position: normal irradiance of each component (zero when
    the Earth is below the terrain horizon) and the Earth's apparent direction."""

    elevation_rad: FloatArray
    azimuth_rad: FloatArray
    visible_fraction: FloatArray
    illuminated_fraction: FloatArray
    reflected_normal_w_m2: FloatArray
    infrared_normal_w_m2: FloatArray

    def on(self, normal: ArrayLike) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Reflected and infrared irradiance on a surface, and the incidence angle."""
        cosine = np.clip(incidence_cosine(normal, self.elevation_rad, self.azimuth_rad), 0.0, 1.0)
        return (
            np.asarray(self.reflected_normal_w_m2 * cosine),
            np.asarray(self.infrared_normal_w_m2 * cosine),
            np.asarray(np.arccos(cosine)),
        )


def earthshine(
    horizons: HorizonSet,
    et_s: FloatArray,
    earth_position_m: FloatArray,
    sun_position_m: FloatArray,
    earth_radius_m: float,
    *,
    observer: int = 0,
    projection: SouthPolarStereographic | None = None,
) -> Earthshine:
    """Earthshine from the HLS-UG-001 section 11.1 approximation.

    Normal irradiance ``(0.100 gamma + 0.056) (d_mean / d)^2`` W m-2 times the visible
    fraction of the Earth's disc, with ``gamma = (1 + cos alpha) / 2`` the illuminated
    fraction of the disc at Sun-Earth-observer phase angle ``alpha``.
    """
    projection = projection or SouthPolarStereographic()
    series = illumination_series(
        horizons,
        et_s,
        earth_position_m,
        earth_radius_m,
        observer=observer,
        projection=projection,
        limb_u=0.0,
    )
    site = projection.to_cartesian(
        horizons.x_m[observer], horizons.y_m[observer], horizons.observer_height_m[observer]
    )
    earth = np.asarray(earth_position_m, dtype=np.float64)
    to_sun = np.asarray(sun_position_m, dtype=np.float64) - earth
    to_site = np.asarray(site)[None, :] - earth
    distance = np.linalg.norm(to_site, axis=1)
    cos_alpha = np.sum(to_sun * to_site, axis=1) / (np.linalg.norm(to_sun, axis=1) * distance)
    gamma = 0.5 * (1.0 + cos_alpha)
    scale = (EARTH_MOON_MEAN_DISTANCE_M / distance) ** 2 * series.visible_fraction
    return Earthshine(
        series.elevation_rad,
        series.azimuth_rad,
        series.visible_fraction,
        gamma,
        np.asarray(EARTHSHINE_REFLECTED_W_M2 * gamma * scale),
        np.asarray(EARTHSHINE_INFRARED_W_M2 * scale),
    )
