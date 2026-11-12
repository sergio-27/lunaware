"""Synthetic terrain and independent reference computations for verification cases."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.data.dem import Dem, GridGeometry
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.geometry.visibility import apparent_elevation

PROJ = SouthPolarStereographic()


def analytic_dem(
    height: Callable[[np.ndarray, np.ndarray], np.ndarray],
    *,
    center: tuple[float, float],
    half_width_m: float,
    pixel_m: float,
) -> Dem:
    """DEM sampled from ``height(x, y)`` at pixel centres."""
    n = round(2 * half_width_m / pixel_m)
    grid = GridGeometry(
        x_min=center[0] - half_width_m,
        y_max=center[1] + half_width_m,
        pixel_size_m=pixel_m,
        width=n,
        height=n,
        projection=PROJ,
    )
    xx, yy = np.meshgrid(grid.x_centers(), grid.y_centers())
    return Dem(np.asarray(height(xx, yy), dtype=np.float32), grid)


def reference_horizon(
    height: Callable[[np.ndarray, np.ndarray], np.ndarray],
    observer_xy: tuple[float, float],
    observer_height_m: float,
    grid_azimuth_rad: np.ndarray,
    *,
    max_distance_m: float,
    step_m: float,
) -> np.ndarray:
    """Horizon elevation computed without the tracer's code path.

    Ground points are generated with the spherical direct geodesic formula in
    longitude/latitude, from the *true* azimuth (grid azimuth minus the longitude), and
    only then projected. The terrain height is the analytic function, not a raster.
    """
    lon0, lat0 = PROJ.inverse(*observer_xy)
    psi = np.arange(step_m, max_distance_m, step_m) / R
    result = np.empty(grid_azimuth_rad.size)
    for i, grid_azimuth in enumerate(grid_azimuth_rad):
        true_azimuth = grid_azimuth - lon0
        sin_lat = np.sin(lat0) * np.cos(psi) + np.cos(lat0) * np.sin(psi) * np.cos(true_azimuth)
        lat = np.arcsin(np.clip(sin_lat, -1.0, 1.0))
        lon = lon0 + np.arctan2(
            np.sin(true_azimuth) * np.sin(psi) * np.cos(lat0),
            np.cos(psi) - np.sin(lat0) * sin_lat,
        )
        x, y = PROJ.forward(lon, lat)
        result[i] = apparent_elevation(psi * R, height(x, y), observer_height_m).max()
    return result
