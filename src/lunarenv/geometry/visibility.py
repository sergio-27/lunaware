"""Line-of-sight geometry on a sphere: how far away terrain can matter.

These closed-form relations answer the data question "how large must the terrain model
be?" (PLAN.md decision D4) and later give the horizon tracer its exact stopping rule.

Setting: reference sphere of radius ``R``; observer at height ``h_o`` and a target at
height ``h_t`` above it, separated by a ground distance ``d`` measured along the
reference sphere (central angle ``psi = d / R``). Elevation is measured from the
observer's local horizontal plane. Heights may be negative.

Criterion for the extent of the far-field DEM (Sun). Terrain can hide part of the solar
disc only if it rises above the lowest elevation the disc ever reaches at the site,

    e_low = -(beta_max + colatitude + rho_max + parallax),

where ``beta_max`` bounds the sub-solar latitude (the lunar obliquity to the ecliptic),
``rho_max`` bounds the Sun's angular radius, and ``parallax`` bounds the difference
between the Sun's direction from the site and from the Moon's centre. A DEM reaching a
ground distance ``d`` from the observer is sufficient when no terrain beyond ``d`` can
appear above ``e_low``, i.e. when :func:`occluding_height` at ``d`` exceeds the highest
terrain that exists (it increases with distance from there on).

The two astronomical bounds below are verified against DE440 over 2024-2044, a full
18.6-year nodal cycle, in ``tests/data/test_visibility_bounds.py``.
"""

from __future__ import annotations

from typing import Final

import numpy as np
from numpy.typing import ArrayLike, NDArray

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M
from lunarenv.geometry.projection import FloatArray

#: Upper bound of the absolute sub-solar latitude on the Moon, radians (1.60 deg).
MAX_SUBSOLAR_LATITUDE_RAD: Final = float(np.radians(1.60))

#: Upper bound of the Sun's angular radius seen from the Moon, radians (0.2720 deg).
MAX_SOLAR_ANGULAR_RADIUS_RAD: Final = float(np.radians(0.2720))

#: Upper bound of the angle between the Sun seen from the surface and from the centre.
MAX_SOLAR_PARALLAX_RAD: Final = 1.2e-5

#: Height that exceeds all lunar terrain relative to the 1737.4 km sphere. The highest
#: point in the LOLA global model is about 10.8 km (Smith2010, to be verified); the cap
#: only needs to be an upper bound.
TERRAIN_HEIGHT_CAP_M: Final = 11_000.0


def apparent_elevation(
    ground_distance_m: ArrayLike,
    target_height_m: ArrayLike,
    observer_height_m: ArrayLike,
    radius_m: float = MOON_REFERENCE_RADIUS_M,
) -> FloatArray:
    """Elevation angle (radians) of a target seen by the observer, on a sphere."""
    psi = np.asarray(ground_distance_m, dtype=np.float64) / radius_m
    r_target = radius_m + np.asarray(target_height_m, dtype=np.float64)
    r_observer = radius_m + np.asarray(observer_height_m, dtype=np.float64)
    return np.asarray(
        np.arctan2(r_target * np.cos(psi) - r_observer, r_target * np.sin(psi)), dtype=np.float64
    )


def occluding_height(
    ground_distance_m: ArrayLike,
    observer_height_m: ArrayLike,
    elevation_rad: ArrayLike,
    radius_m: float = MOON_REFERENCE_RADIUS_M,
) -> FloatArray:
    """Height a target at a given distance needs in order to appear at ``elevation_rad``.

    Inverse of :func:`apparent_elevation` in the target height:
    ``R + h_t = (R + h_o) cos(e) / cos(psi + e)``. Returns ``inf`` where
    ``psi + e >= pi/2``: no finite height reaches that elevation from there.
    """
    psi = np.asarray(ground_distance_m, dtype=np.float64) / radius_m
    elevation = np.asarray(elevation_rad, dtype=np.float64)
    r_observer = radius_m + np.asarray(observer_height_m, dtype=np.float64)
    cosine = np.cos(psi + elevation)
    with np.errstate(divide="ignore", invalid="ignore"):
        height = r_observer * np.cos(elevation) / cosine - radius_m
    return np.asarray(np.where(cosine > 0.0, height, np.inf), dtype=np.float64)


def horizon_dip(
    observer_height_m: ArrayLike, radius_m: float = MOON_REFERENCE_RADIUS_M
) -> FloatArray:
    """Elevation (negative, radians) of the horizon of a smooth sphere of radius ``radius_m``.

    ``observer_height_m`` is the height above that sphere and must not be negative.
    """
    height = np.asarray(observer_height_m, dtype=np.float64)
    if np.any(height < 0.0):
        raise ValueError("observer must be on or above the sphere")
    return np.asarray(-np.arccos(radius_m / (radius_m + height)), dtype=np.float64)


def lowest_solar_limb_elevation(latitude_rad: ArrayLike) -> FloatArray:
    """Lower bound (radians) of the elevation of any point of the solar disc at a site.

    Valid for any epoch. Flat-horizon geometry; the site's own terrain does not enter.
    """
    colatitude = np.pi / 2.0 - np.abs(np.asarray(latitude_rad, dtype=np.float64))
    return np.asarray(
        -(
            MAX_SUBSOLAR_LATITUDE_RAD
            + colatitude
            + MAX_SOLAR_ANGULAR_RADIUS_RAD
            + MAX_SOLAR_PARALLAX_RAD
        ),
        dtype=np.float64,
    )


def dem_extent_is_sufficient_for_sun(
    distance_to_dem_edge_m: ArrayLike,
    observer_height_m: ArrayLike,
    latitude_rad: ArrayLike,
    *,
    terrain_height_cap_m: float = TERRAIN_HEIGHT_CAP_M,
    radius_m: float = MOON_REFERENCE_RADIUS_M,
) -> NDArray[np.bool_]:
    """True where terrain beyond the DEM edge provably cannot hide any part of the Sun.

    ``distance_to_dem_edge_m`` is the shortest ground distance from the observer to the
    boundary of the terrain model. The test is conservative: it assumes terrain as high
    as ``terrain_height_cap_m`` could sit right at the edge.

    The height needed to reach a negative elevation ``e`` first decreases with distance
    and then increases beyond the grazing point ``psi = -e`` (about 64 km for
    ``e = -2.1 deg``). The edge value bounds everything farther out only past that
    point, so observers whose DEM edge is closer are reported as not covered.
    """
    limb = lowest_solar_limb_elevation(latitude_rad)
    edge = np.asarray(distance_to_dem_edge_m, dtype=np.float64)
    needed = occluding_height(edge, observer_height_m, limb, radius_m)
    beyond_grazing_point = edge / radius_m + limb >= 0.0
    return np.asarray((needed > terrain_height_cap_m) & beyond_grazing_point)
