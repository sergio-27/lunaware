"""Energy input at the surface: absorbed sunlight on a tilted surface element.

    Q_abs = (1 - A(theta)) * E_sun * cos(theta) * f_eclipse

with ``E_sun`` the direct solar irradiance on a plane facing the Sun (solar constant
scaled by distance and by the limb-darkened visible fraction of the disc), ``theta`` the
angle between the Sun's centre and the surface normal, ``A`` the incidence-dependent
albedo of Hayne2017 and ``f_eclipse`` the Earth-occultation factor.

What is left out, and what it does (PLAN.md, decision 1 of section 11):

* **Sunlight scattered by surrounding terrain and thermal emission from it.** Zero
  effect on a surface in full sunlight compared with the direct beam; the only energy
  source of a surface in shadow, where the model is therefore too cold.
  :func:`neglected_terrain_flux_bound` gives an upper bound from the horizon profile:
  terrain cannot send more than it receives, so the flux from terrain is at most the
  solar irradiance times the fraction of the receiver's view occupied by terrain.
* **Partial discs**: the incidence angle is that of the disc centre.
* **Sub-pixel roughness**: the normal is that of the DEM at its sampling baseline. The
  albedo law is an empirical fit that already absorbs roughness effects on level
  equatorial ground; applying it to a tilted polar facet is an extrapolation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from lunarenv.core.constants import SOLAR_CONSTANT_W_M2, STEFAN_BOLTZMANN_W_M2_K4
from lunarenv.data.dem import Dem
from lunarenv.geometry.projection import FloatArray
from lunarenv.illumination.products import IlluminationSeries
from lunarenv.thermal.properties import RegolithProperties, albedo

#: Smallest Sun-Moon distance in au (DE440, 2024-2044: 0.9806), for flux upper bounds.
_MIN_SUN_DISTANCE_AU = 0.98


def normal_from_slope(slope_rad: float, facing_azimuth_rad: float) -> FloatArray:
    """Unit normal, in the local grid frame (e_x, e_y, up), of a plane tilted by
    ``slope_rad`` and facing (downhill towards) grid azimuth ``facing_azimuth_rad``."""
    if not 0.0 <= slope_rad < np.pi / 2:
        raise ValueError("slope must be in [0, 90) degrees")
    s = np.sin(slope_rad)
    return np.array(
        [s * np.sin(facing_azimuth_rad), s * np.cos(facing_azimuth_rad), np.cos(slope_rad)]
    )


def surface_normal(
    dem: Dem, x_m: float, y_m: float, *, baseline_m: float | None = None
) -> FloatArray:
    """Unit surface normal of the DEM at a point, in the local grid frame.

    Central differences of the interpolated surface over ``baseline_m`` on either side
    (default: one pixel). Map distances are converted to ground distances at the height
    of the surface with the projection scale factor.
    """
    grid = dem.grid
    step = baseline_m if baseline_m is not None else grid.pixel_size_m
    if step <= 0.0:
        raise ValueError("baseline must be positive")
    height = float(dem.sample(x_m, y_m))
    radius = grid.projection.radius_m
    # ground metres per map metre at the surface height
    ground = (radius + height) / radius / float(grid.projection.scale_factor(x_m, y_m))
    dh_dx = float(dem.sample(x_m + step, y_m) - dem.sample(x_m - step, y_m)) / (2 * step * ground)
    dh_dy = float(dem.sample(x_m, y_m + step) - dem.sample(x_m, y_m - step)) / (2 * step * ground)
    normal = np.array([-dh_dx, -dh_dy, 1.0])
    return np.asarray(normal / np.linalg.norm(normal))


def slope_and_facing(normal: ArrayLike) -> tuple[float, float]:
    """Slope (radians from horizontal) and facing grid azimuth of a unit normal."""
    n = np.asarray(normal, dtype=np.float64)
    slope = float(np.arctan2(np.hypot(n[0], n[1]), n[2]))
    return slope, float(np.mod(np.arctan2(n[0], n[1]), 2 * np.pi))


def incidence_cosine(
    normal: ArrayLike, elevation_rad: ArrayLike, azimuth_rad: ArrayLike
) -> FloatArray:
    """Cosine of the angle between the surface normal and a direction in the sky."""
    n = np.asarray(normal, dtype=np.float64)
    elevation = np.asarray(elevation_rad, dtype=np.float64)
    azimuth = np.asarray(azimuth_rad, dtype=np.float64)
    horizontal = np.cos(elevation)
    return np.asarray(
        n[0] * horizontal * np.sin(azimuth)
        + n[1] * horizontal * np.cos(azimuth)
        + n[2] * np.sin(elevation)
    )


def absorbed_solar_flux(
    series: IlluminationSeries,
    normal: ArrayLike,
    properties: RegolithProperties,
    *,
    eclipse_factor: ArrayLike | None = None,
) -> FloatArray:
    """Absorbed direct solar flux (W m-2) on a surface element at every epoch."""
    cosine = np.clip(incidence_cosine(normal, series.elevation_rad, series.azimuth_rad), 0.0, 1.0)
    irradiance = series.normal_irradiance_w_m2
    if eclipse_factor is not None:
        irradiance = irradiance * np.asarray(eclipse_factor, dtype=np.float64)
    absorptance = 1.0 - albedo(np.arccos(cosine), properties)
    return np.asarray(absorptance * irradiance * cosine)


def analytic_absorbed_flux(
    hour_angle_rad: ArrayLike,
    latitude_rad: float,
    properties: RegolithProperties,
    *,
    declination_rad: float = 0.0,
    sun_distance_au: float = 1.0,
) -> FloatArray:
    """Absorbed flux on level, unobstructed ground (Hayne2017, A9-A11).

    Used for the equatorial and Apollo-site verification cases, where no terrain model
    is involved. Hour angle 0 is local noon.
    """
    hour = np.asarray(hour_angle_rad, dtype=np.float64)
    cosine = np.sin(latitude_rad) * np.sin(declination_rad) + np.cos(latitude_rad) * np.cos(
        declination_rad
    ) * np.cos(hour)
    cosine = np.clip(cosine, 0.0, 1.0)
    irradiance = SOLAR_CONSTANT_W_M2 / sun_distance_au**2
    return np.asarray((1.0 - albedo(np.arccos(cosine), properties)) * irradiance * cosine)


@dataclass(frozen=True, slots=True)
class NeglectedFluxBound:
    """Upper bound of the energy the model does not include (terrain to surface).

    ``terrain_view_factor`` is the fraction of a level receiver's hemisphere-weighted
    view filled by terrain; ``flux_w_m2`` bounds scattered sunlight plus terrain thermal
    emission reaching it; ``shadow_temperature_floor_k`` is the temperature that flux
    alone would sustain.
    """

    terrain_view_factor: float
    flux_w_m2: float
    shadow_temperature_floor_k: float


def neglected_terrain_flux_bound(
    horizon_elevation_rad: ArrayLike, properties: RegolithProperties
) -> NeglectedFluxBound:
    """Bound on the terrain-to-surface flux from an observer's horizon profile.

    For a level Lambertian receiver, terrain below elevation ``e_h(azimuth)`` occupies
    the view factor ``mean(sin^2(max(e_h, 0)))``. A terrain element cannot radiate
    (reflected plus thermal) more than the solar irradiance it could receive, so the
    flux is at most that irradiance, at the smallest Sun distance, times the view
    factor. The bound is loose by design: most visible terrain is obliquely lit.
    """
    elevation = np.clip(np.asarray(horizon_elevation_rad, dtype=np.float64), 0.0, None)
    view = float(np.mean(np.sin(elevation) ** 2))
    flux = view * SOLAR_CONSTANT_W_M2 / _MIN_SUN_DISTANCE_AU**2
    floor = (flux / (properties.emissivity * STEFAN_BOLTZMANN_W_M2_K4)) ** 0.25
    return NeglectedFluxBound(view, float(flux), float(floor))


def plane_normal(tilt_rad: float, facing_azimuth_rad: float) -> FloatArray:
    """Unit normal of a plane tilted by ``tilt_rad`` from horizontal (0 to 90 degrees
    inclusive: a vertical wall is allowed), facing grid azimuth ``facing_azimuth_rad``."""
    if not 0.0 <= tilt_rad <= np.pi / 2 + 1e-12:
        raise ValueError("tilt must be in [0, 90] degrees")
    s = np.sin(tilt_rad)
    return np.array(
        [s * np.sin(facing_azimuth_rad), s * np.cos(facing_azimuth_rad), np.cos(tilt_rad)]
    )


def direct_irradiance(
    series: IlluminationSeries, normal: ArrayLike, *, eclipse_factor: ArrayLike | None = None
) -> FloatArray:
    """Direct solar irradiance (W m-2) on a plane, before any absorption.

    Limb-darkened visible fraction of the disc above the terrain horizon, Earth
    occultation, and the cosine of the incidence on the plane (zero from behind).
    """
    cosine = np.clip(incidence_cosine(normal, series.elevation_rad, series.azimuth_rad), 0.0, 1.0)
    irradiance = series.normal_irradiance_w_m2
    if eclipse_factor is not None:
        irradiance = irradiance * np.asarray(eclipse_factor, dtype=np.float64)
    return np.asarray(irradiance * cosine)


@dataclass(frozen=True, slots=True)
class ViewFactors:
    """Fractions of a plane's diffuse (cosine-weighted) view filled by sky and by ground.

    ``sky + ground = 1``. ``ground`` is the terrain above the local horizontal plus the
    ground below it that faces the plane; it is the factor that scales terrain-scattered
    sunlight and terrain infrared reaching the plane.
    """

    sky: float
    ground: float


def view_factors(
    horizon_elevation_rad: ArrayLike,
    azimuth_rad: ArrayLike,
    normal: ArrayLike,
    *,
    n_elevation: int = 360,
) -> ViewFactors:
    """Sky and ground view factors of a plane from the horizon profile at its location.

    The sky is every direction above the terrain horizon ``e_h(azimuth)`` (which may be
    negative for an elevated observer). The view factor integrates ``max(cos, 0) / pi``
    of the angle to the normal over that solid angle, by the midpoint rule in elevation
    and the horizon's own azimuth grid. For a level plane and non-negative horizons it
    equals ``1 - mean(sin^2 e_h)``.
    """
    horizon = np.asarray(horizon_elevation_rad, dtype=np.float64)
    azimuth = np.asarray(azimuth_rad, dtype=np.float64)
    if horizon.shape != azimuth.shape or horizon.ndim != 1:
        raise ValueError("horizon and azimuth must be 1-D arrays of the same length")
    d_azimuth = 2.0 * np.pi / azimuth.size
    # midpoints between the horizon and the zenith, per azimuth
    u = (np.arange(n_elevation) + 0.5) / n_elevation
    elevation = horizon[:, None] + (np.pi / 2 - horizon[:, None]) * u[None, :]
    d_elevation = (np.pi / 2 - horizon) / n_elevation
    cosine = np.clip(incidence_cosine(normal, elevation, azimuth[:, None]), 0.0, None)
    sky = float(np.sum(cosine * np.cos(elevation) * d_elevation[:, None]) * d_azimuth / np.pi)
    sky = min(max(sky, 0.0), 1.0)
    return ViewFactors(sky, 1.0 - sky)


__all__ = [
    "NeglectedFluxBound",
    "ViewFactors",
    "absorbed_solar_flux",
    "analytic_absorbed_flux",
    "direct_irradiance",
    "incidence_cosine",
    "neglected_terrain_flux_bound",
    "normal_from_slope",
    "plane_normal",
    "slope_and_facing",
    "surface_normal",
    "view_factors",
]
