"""Illumination products: time series at a point and statistics over a map.

Definitions follow Barker2021 so that results are comparable:

* *illumination* at an instant is the visible fraction of the solar disc;
* *average illumination* is its mean over the time samples of a period;
* the *longest continuous shadow period* is the longest run of samples in which no part
  of the disc is visible; the *longest continuous illumination period* is the longest
  run in which any part is visible.

Run lengths are counted in whole time steps (``n * dt``), so they carry a resolution
of one step at each end. The time-step convergence study quantifies the effect.

The apparent direction of the body is the Moon-centred apparent position minus the
observer's body-fixed position (accuracy: see ``lunarenv.geometry.ephemeris``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import get_parallel_chunksize, njit, prange, set_parallel_chunksize
from numpy.typing import NDArray

from lunarenv.core.constants import ASTRONOMICAL_UNIT_M, SOLAR_CONSTANT_W_M2
from lunarenv.geometry.projection import FloatArray, SouthPolarStereographic
from lunarenv.illumination.disc import EDDINGTON_LIMB_DARKENING, PROFILE, disc_fractions
from lunarenv.illumination.horizon import PARALLEL_CHUNK, HorizonSet


@dataclass(frozen=True, slots=True)
class Run:
    """A run of consecutive time samples: ``count`` samples starting at ``start_index``."""

    count: int
    start_index: int
    duration_s: float


@dataclass(frozen=True, slots=True)
class IlluminationSeries:
    """Apparent geometry and visibility of a body from one observer over time.

    Angles in radians (grid azimuth); ``visible_fraction`` is geometric,
    ``flux_fraction`` is weighted by limb darkening. ``normal_irradiance_w_m2`` is the
    direct solar irradiance on a surface facing the Sun (meaningful for the Sun only).
    """

    et_s: FloatArray
    step_s: float
    elevation_rad: FloatArray
    azimuth_rad: FloatArray
    angular_radius_rad: FloatArray
    distance_m: FloatArray
    visible_fraction: FloatArray
    flux_fraction: FloatArray

    @property
    def normal_irradiance_w_m2(self) -> FloatArray:
        scale = (ASTRONOMICAL_UNIT_M / self.distance_m) ** 2
        return np.asarray(SOLAR_CONSTANT_W_M2 * scale * self.flux_fraction)

    @property
    def mean_visible_fraction(self) -> float:
        return float(self.visible_fraction.mean())

    def longest_shadow(self) -> Run:
        return longest_run(self.visible_fraction <= 0.0, self.step_s)

    def longest_illumination(self) -> Run:
        return longest_run(self.visible_fraction > 0.0, self.step_s)


def longest_run(mask: NDArray[np.bool_], step_s: float) -> Run:
    """Longest run of ``True`` in a boolean series (the first one if several tie)."""
    flags = np.asarray(mask, dtype=np.bool_)
    if flags.ndim != 1:
        raise ValueError("mask must be one-dimensional")
    if not flags.any():
        return Run(0, -1, 0.0)
    padded = np.concatenate(([False], flags, [False])).astype(np.int8)
    edges = np.diff(padded)
    starts = np.flatnonzero(edges == 1)
    lengths = np.flatnonzero(edges == -1) - starts
    best = int(np.argmax(lengths))
    return Run(int(lengths[best]), int(starts[best]), float(lengths[best] * step_s))


@njit(cache=True)
def _apparent(bx, by, bz, x0, y0, h0, radius, body_radius):  # type: ignore[no-untyped-def]
    """Elevation, grid azimuth, angular radius and distance of a body from an observer.

    The local basis is the grid-aligned one of ``SouthPolarStereographic.local_basis``,
    written out in scalars for the compiler.
    """
    t2 = (x0 * x0 + y0 * y0) / (4.0 * radius * radius)
    s = 1.0 / (1.0 + t2)
    q = s / (2.0 * radius * radius)
    up0, up1, up2 = y0 * s / radius, x0 * s / radius, -(1.0 - t2) * s
    r = radius + h0
    dx = bx - up0 * r
    dy = by - up1 * r
    dz = bz - up2 * r
    distance = math.sqrt(dx * dx + dy * dy + dz * dz)
    vertical = dx * up0 + dy * up1 + dz * up2
    along_x = dx * (-x0 * y0 * q) + dy * (1.0 - x0 * x0 * q) + dz * (x0 * s / radius)
    along_y = dx * (1.0 - y0 * y0 * q) + dy * (-x0 * y0 * q) + dz * (y0 * s / radius)
    elevation = math.atan2(vertical, math.sqrt(along_x * along_x + along_y * along_y))
    azimuth = math.atan2(along_x, along_y)
    if azimuth < 0.0:
        azimuth += 2.0 * math.pi
    return elevation, azimuth, math.asin(body_radius / distance), distance


@njit(cache=True, parallel=True)
def _series_kernel(  # type: ignore[no-untyped-def]
    body, x0, y0, h0, radius, body_radius, horizon, n_strips, limb_u, mode, out
):
    for i in prange(body.shape[0]):
        elevation, azimuth, rho, distance = _apparent(
            body[i, 0], body[i, 1], body[i, 2], x0, y0, h0, radius, body_radius
        )
        geometric, flux = disc_fractions(elevation, azimuth, rho, horizon, n_strips, limb_u, mode)
        out[i, 0] = elevation
        out[i, 1] = azimuth
        out[i, 2] = rho
        out[i, 3] = distance
        out[i, 4] = geometric
        out[i, 5] = flux


@njit(cache=True, parallel=True)
def _map_kernel(  # type: ignore[no-untyped-def]
    body, obs_x, obs_y, obs_h, radius, body_radius, horizons, n_strips, limb_u, mode, out
):
    n_times = body.shape[0]
    for i in prange(obs_x.size):
        horizon = horizons[i].astype(np.float64)
        sum_geometric = 0.0
        sum_flux = 0.0
        dark = 0
        lit = 0
        longest_dark = 0
        longest_lit = 0
        for t in range(n_times):
            elevation, azimuth, rho, _ = _apparent(
                body[t, 0],
                body[t, 1],
                body[t, 2],
                obs_x[i],
                obs_y[i],
                obs_h[i],
                radius,
                body_radius,
            )
            geometric, flux = disc_fractions(
                elevation, azimuth, rho, horizon, n_strips, limb_u, mode
            )
            sum_geometric += geometric
            sum_flux += flux
            if geometric > 0.0:
                lit += 1
                dark = 0
                if lit > longest_lit:
                    longest_lit = lit
            else:
                dark += 1
                lit = 0
                if dark > longest_dark:
                    longest_dark = dark
        out[i, 0] = sum_geometric / n_times
        out[i, 1] = sum_flux / n_times
        out[i, 2] = longest_dark
        out[i, 3] = longest_lit


#: Largest departure from a uniform step accepted on a time axis, seconds. A grid that
#: is uniform in UTC is not exactly uniform in TDB: periodic relativistic terms move it
#: by under 2 ms, and a leap second inside the period shifts it by 1 s.
_STEP_TOLERANCE_S = 1.5


def _uniform_step(et_s: FloatArray) -> float:
    if et_s.ndim != 1 or et_s.size < 2:
        raise ValueError("need a one-dimensional time axis with at least two samples")
    steps = np.diff(et_s)
    step = float(np.median(steps))
    if step <= 0.0 or float(np.abs(steps - step).max()) > _STEP_TOLERANCE_S:
        raise ValueError("the time axis must be uniformly spaced and increasing")
    return step


def illumination_series(
    horizons: HorizonSet,
    et_s: FloatArray,
    body_position_m: FloatArray,
    body_radius_m: float,
    *,
    observer: int = 0,
    projection: SouthPolarStereographic | None = None,
    n_strips: int = 32,
    limb_u: float = EDDINGTON_LIMB_DARKENING,
    mode: int = PROFILE,
) -> IlluminationSeries:
    """Visibility of a body over time from one observer of a horizon set.

    ``body_position_m`` is the Moon-centred apparent position of the body in the
    body-fixed frame at each epoch, shape ``(n, 3)``
    (:meth:`lunarenv.geometry.ephemeris.Ephemeris.position_from_center`).
    """
    projection = projection or SouthPolarStereographic()
    et = np.asarray(et_s, dtype=np.float64)
    step = _uniform_step(et)
    body = np.ascontiguousarray(body_position_m, dtype=np.float64)
    if body.shape != (et.size, 3):
        raise ValueError("body positions must have shape (n_times, 3)")
    out = np.empty((et.size, 6), dtype=np.float64)
    _series_kernel(
        body,
        float(horizons.x_m[observer]),
        float(horizons.y_m[observer]),
        float(horizons.observer_height_m[observer]),
        projection.radius_m,
        body_radius_m,
        np.ascontiguousarray(horizons.elevation_rad[observer], dtype=np.float64),
        n_strips,
        limb_u,
        mode,
        out,
    )
    return IlluminationSeries(et, step, *(np.ascontiguousarray(out[:, c]) for c in range(6)))


@dataclass(frozen=True, slots=True)
class IlluminationStatistics:
    """Per-observer statistics over a period (same order as the horizon set)."""

    mean_visible_fraction: FloatArray
    mean_flux_fraction: FloatArray
    longest_shadow_s: FloatArray
    longest_illumination_s: FloatArray
    step_s: float
    n_times: int


def illumination_statistics(
    horizons: HorizonSet,
    et_s: FloatArray,
    body_position_m: FloatArray,
    body_radius_m: float,
    *,
    projection: SouthPolarStereographic | None = None,
    n_strips: int = 32,
    limb_u: float = EDDINGTON_LIMB_DARKENING,
    mode: int = PROFILE,
) -> IlluminationStatistics:
    """Mean illumination and longest shadow/illumination runs for every observer."""
    projection = projection or SouthPolarStereographic()
    et = np.asarray(et_s, dtype=np.float64)
    step = _uniform_step(et)
    body = np.ascontiguousarray(body_position_m, dtype=np.float64)
    if body.shape != (et.size, 3):
        raise ValueError("body positions must have shape (n_times, 3)")
    out = np.empty((horizons.x_m.size, 4), dtype=np.float64)
    previous_chunk = get_parallel_chunksize()
    set_parallel_chunksize(PARALLEL_CHUNK)  # observers differ in cost; keep threads busy
    try:
        _map_kernel(
            body,
            horizons.x_m,
            horizons.y_m,
            horizons.observer_height_m,
            projection.radius_m,
            body_radius_m,
            np.ascontiguousarray(horizons.elevation_rad),
            n_strips,
            limb_u,
            mode,
            out,
        )
    finally:
        set_parallel_chunksize(previous_chunk)
    return IlluminationStatistics(
        mean_visible_fraction=np.ascontiguousarray(out[:, 0]),
        mean_flux_fraction=np.ascontiguousarray(out[:, 1]),
        longest_shadow_s=out[:, 2] * step,
        longest_illumination_s=out[:, 3] * step,
        step_s=step,
        n_times=int(et.size),
    )
