"""Horizon elevation seen from points on (or above) the terrain.

Method (PLAN.md decisions D1-D4). For each observer and each of ``n_azimuth`` equally
spaced directions, a ray is marched along the great circle leaving the observer in that
direction. At each sample the terrain height is interpolated bilinearly and the
elevation of that terrain point above the observer's horizontal plane is computed from
exact spherical geometry:

    tan(e) = ((R + h) cos(psi) - (R + h_obs)) / ((R + h) sin(psi)),   psi = d / R.

The horizon is the maximum over the ray. Lunar curvature is therefore built in rather
than added as a correction.

Stopping rule. With ``h_cap`` the highest terrain in the stack, no terrain farther than
the current sample can appear above the elevation that ``h_cap`` would have at the
current sample (that elevation decreases with distance once past the tangent point).
When this bound drops below the horizon already found, the ray is finished and marked
*complete*. A ray that runs out of terrain data before that is marked *incomplete*: the
data cannot prove that nothing farther away rises higher. Whether an incomplete ray
matters for the Sun is decided by :func:`sun_relevant_incomplete`.

Azimuth convention: *grid azimuth*, measured from map +y towards map +x (see
``lunarenv.geometry.projection``). Direction ``k`` has azimuth ``2 pi k / n_azimuth``.

Assumptions and limits:

* the terrain is the DEM surface, bilinear between pixel centres; relief below the DEM
  resolution (boulders, small craters) does not exist for the tracer;
* the observer sits ``height_above_surface_m`` above the interpolated surface, along
  the local vertical; at zero height the horizon includes the local slope;
* azimuth sampling is discrete; features narrower than the azimuth step can be missed
  between rays (quantified in the convergence study).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import get_parallel_chunksize, njit, prange, set_parallel_chunksize
from numpy.typing import ArrayLike, NDArray

from lunarenv.geometry.projection import FloatArray
from lunarenv.geometry.visibility import (
    TERRAIN_HEIGHT_CAP_M,
    lowest_solar_limb_elevation,
    occluding_height,
)
from lunarenv.illumination.terrain import TerrainStack

_STATUS_BOUND = 2  # stopped by the bound: complete
_CHECK_EVERY = 8  # samples between evaluations of the stopping rule
#: Observers handed to a thread at a time. The cost per observer varies by an order of
#: magnitude (rays end early in enclosed terrain), so the default split into one block
#: per thread leaves most threads idle; small blocks keep all of them busy.
PARALLEL_CHUNK = 16


@dataclass(frozen=True, slots=True)
class HorizonSet:
    """Horizons of a set of observers.

    ``elevation_rad[i, k]`` is the horizon elevation of observer ``i`` at grid azimuth
    ``2 pi k / n_azimuth``. ``distance_m`` is the ground distance of the terrain point
    that defines it; ``stop_distance_m`` is where the ray ended; ``complete`` says
    whether the stopping rule proved the horizon final. The two distance arrays are
    ``None`` when diagnostics were not requested (large maps).
    """

    x_m: FloatArray
    y_m: FloatArray
    observer_height_m: FloatArray
    elevation_rad: NDArray[np.float32]
    complete: NDArray[np.bool_]
    distance_m: NDArray[np.float32] | None
    stop_distance_m: NDArray[np.float32] | None

    @property
    def n_azimuth(self) -> int:
        return int(self.elevation_rad.shape[1])

    @property
    def azimuth_rad(self) -> FloatArray:
        return 2.0 * np.pi * np.arange(self.n_azimuth) / self.n_azimuth


@njit(cache=True, inline="always")
def _bilinear(dem, x_min, y_max, inverse_pixel, x, y):  # type: ignore[no-untyped-def]
    """Interpolated height, or (False, 0) outside the rectangle of pixel centres."""
    col = (x - x_min) * inverse_pixel - 0.5
    row = (y_max - y) * inverse_pixel - 0.5
    n_rows, n_cols = dem.shape
    if col < 0.0 or row < 0.0 or col > n_cols - 1 or row > n_rows - 1:
        return False, 0.0
    c0 = min(int(col), n_cols - 2)
    r0 = min(int(row), n_rows - 2)
    fc = col - c0
    fr = row - r0
    top = dem[r0, c0] * (1.0 - fc) + dem[r0, c0 + 1] * fc
    bottom = dem[r0 + 1, c0] * (1.0 - fc) + dem[r0 + 1, c0 + 1] * fc
    return True, top * (1.0 - fr) + bottom * fr


@njit(cache=True, parallel=True)
def _trace(  # type: ignore[no-untyped-def]
    dems,
    grid_meta,
    plan_cos,
    plan_sin,
    plan_level,
    obs_x,
    obs_y,
    obs_h,
    n_azimuth,
    radius,
    h_cap,
    out_tan,
    out_best,
    out_stop,
    out_complete,
    diagnostics,
):
    n_levels = len(dems)
    n_samples = plan_cos.size
    r_cap = radius + h_cap
    for i in prange(obs_x.size):
        x0 = obs_x[i]
        y0 = obs_y[i]
        r_obs = radius + obs_h[i]
        # grid-aligned local basis (same algebra as SouthPolarStereographic.local_basis)
        t2 = (x0 * x0 + y0 * y0) / (4.0 * radius * radius)
        s = 1.0 / (1.0 + t2)
        q = s / (2.0 * radius * radius)
        up0, up1, up2 = y0 * s / radius, x0 * s / radius, -(1.0 - t2) * s
        ex0, ex1, ex2 = -x0 * y0 * q, 1.0 - x0 * x0 * q, x0 * s / radius
        ey0, ey1, ey2 = 1.0 - y0 * y0 * q, -x0 * y0 * q, y0 * s / radius
        # past the tangent point of the cap sphere the bound decreases with distance
        cos_tangent = min(1.0, r_cap / r_obs)
        for k in range(n_azimuth):
            azimuth = 2.0 * math.pi * k / n_azimuth
            sin_a = math.sin(azimuth)
            cos_a = math.cos(azimuth)
            t0 = sin_a * ex0 + cos_a * ey0
            t1 = sin_a * ex1 + cos_a * ey1
            t2_ = sin_a * ex2 + cos_a * ey2
            best = -math.inf
            best_index = -1
            stop_index = n_samples - 1
            complete = False
            for j in range(n_samples):
                c = plan_cos[j]
                sn = plan_sin[j]
                # unit vector of the ground point, then its map coordinates
                u0 = up0 * c + t0 * sn
                u1 = up1 * c + t1 * sn
                u2 = up2 * c + t2_ * sn
                scale = 2.0 * radius / (1.0 - u2)
                x = u1 * scale
                y = u0 * scale
                found = False
                height = 0.0
                for level in range(plan_level[j], n_levels):
                    found, height = _bilinear(
                        dems[level],
                        grid_meta[level, 0],
                        grid_meta[level, 1],
                        grid_meta[level, 2],
                        x,
                        y,
                    )
                    if found:
                        break
                if not found:  # ran out of terrain data
                    stop_index = j
                    bound = (r_cap * c - r_obs) / (r_cap * sn)
                    complete = bound < best and c <= cos_tangent
                    break
                r_terrain = radius + height
                tangent = (r_terrain * c - r_obs) / (r_terrain * sn)
                if tangent > best:
                    best = tangent
                    best_index = j
                if j % _CHECK_EVERY == 0 or j == n_samples - 1:
                    bound = (r_cap * c - r_obs) / (r_cap * sn)
                    if bound < best and c <= cos_tangent:
                        stop_index = j
                        complete = True
                        break
            out_tan[i, k] = best
            out_complete[i, k] = complete
            if diagnostics:
                out_best[i, k] = best_index
                out_stop[i, k] = stop_index


def _run_trace(  # type: ignore[no-untyped-def]
    dems,
    grid_meta,
    psi,
    level,
    x,
    y,
    observer_height,
    n_azimuth,
    radius,
    h_cap,
    tangent,
    best_index,
    stop_index,
    complete,
    diagnostics,
) -> None:
    _trace(
        dems,
        grid_meta,
        np.cos(psi),
        np.sin(psi),
        level,
        x,
        y,
        observer_height,
        n_azimuth,
        radius,
        h_cap,
        tangent,
        best_index,
        stop_index,
        complete,
        diagnostics,
    )


def compute_horizons(
    stack: TerrainStack,
    x_m: ArrayLike,
    y_m: ArrayLike,
    *,
    height_above_surface_m: float = 0.0,
    n_azimuth: int = 720,
    diagnostics: bool = True,
    surface_height_m: ArrayLike | None = None,
) -> HorizonSet:
    """Trace the horizon of each observer over the terrain stack.

    Observers are given by map coordinates and must lie inside the interpolation domain
    of the nearest level, which supplies their surface height unless
    ``surface_height_m`` is given (used to hold the observer fixed while the terrain
    model is varied in resolution studies).
    """
    if n_azimuth < 8:
        raise ValueError("n_azimuth must be at least 8")
    if height_above_surface_m < 0.0:
        raise ValueError("the observer cannot be below the surface")
    x = np.atleast_1d(np.asarray(x_m, dtype=np.float64)).ravel()
    y = np.atleast_1d(np.asarray(y_m, dtype=np.float64)).ravel()
    if x.shape != y.shape:
        raise ValueError("x and y must have the same number of elements")
    if surface_height_m is None:
        surface = stack.levels[0].dem.sample(x, y)
    else:
        surface = np.broadcast_to(np.asarray(surface_height_m, dtype=np.float64), x.shape)
    observer_height = np.ascontiguousarray(surface + height_above_surface_m)

    radius = stack.radius_m
    distance, level = stack.sample_plan()
    psi = distance / radius
    dems = tuple(np.ascontiguousarray(lv.dem.heights_m, dtype=np.float32) for lv in stack.levels)
    grid_meta = np.array(
        [
            [lv.dem.grid.x_min, lv.dem.grid.y_max, 1.0 / lv.dem.grid.pixel_size_m]
            for lv in stack.levels
        ],
        dtype=np.float64,
    )
    shape = (x.size, n_azimuth)
    tangent = np.empty(shape, dtype=np.float64)
    complete = np.empty(shape, dtype=np.bool_)
    diag_shape = shape if diagnostics else (1, 1)
    best_index = np.empty(diag_shape, dtype=np.int32)
    stop_index = np.empty(diag_shape, dtype=np.int32)
    previous_chunk = get_parallel_chunksize()
    set_parallel_chunksize(PARALLEL_CHUNK)
    try:
        _run_trace(
            dems,
            grid_meta,
            psi,
            level,
            x,
            y,
            observer_height,
            n_azimuth,
            radius,
            stack.max_height_m,
            tangent,
            best_index,
            stop_index,
            complete,
            diagnostics,
        )
    finally:
        set_parallel_chunksize(previous_chunk)
    if not np.all(np.isfinite(tangent)):
        raise RuntimeError("a ray found no terrain at all; the stack does not cover the observers")
    return HorizonSet(
        x_m=x,
        y_m=y,
        observer_height_m=observer_height,
        elevation_rad=np.arctan(tangent).astype(np.float32),
        complete=complete,
        distance_m=distance[best_index].astype(np.float32) if diagnostics else None,
        stop_distance_m=distance[stop_index].astype(np.float32) if diagnostics else None,
    )


def sun_relevant_incomplete(horizons: HorizonSet, stack: TerrainStack) -> NDArray[np.bool_]:
    """Rays whose incompleteness could affect solar illumination.

    An incomplete ray matters only if terrain beyond its end could rise above both the
    horizon already found and the lowest elevation the solar disc ever reaches at the
    site. The test assumes terrain up to :data:`TERRAIN_HEIGHT_CAP_M` may exist beyond
    the data. Needs a horizon set computed with diagnostics.
    """
    if horizons.stop_distance_m is None:
        raise ValueError("this horizon set was computed without diagnostics")
    projection = stack.levels[0].dem.grid.projection
    _, latitude = projection.inverse(horizons.x_m, horizons.y_m)
    threshold = np.maximum(
        horizons.elevation_rad.astype(np.float64),
        lowest_solar_limb_elevation(latitude)[:, None],
    )
    stop = horizons.stop_distance_m.astype(np.float64)
    needed = occluding_height(
        stop, horizons.observer_height_m[:, None], threshold, projection.radius_m
    )
    past_grazing = stop / projection.radius_m + threshold >= 0.0
    could_matter = ~((needed > TERRAIN_HEIGHT_CAP_M) & past_grazing)
    return np.asarray(~horizons.complete & could_matter)
