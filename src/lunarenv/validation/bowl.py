"""Spherical bowl crater with a single-bounce analytic solution (TH-A3).

A bowl of rim radius ``r`` and depth ``H`` is a cap of a sphere of radius
``R = (r^2 + H^2) / (2 H)``. Inside a sphere the view factor from any interior element
to another is ``dA / (4 pi R^2)``, so a shadowed interior point receives the power that
enters the opening, ``F sin(e) pi r^2``, scattered or re-emitted once, divided by
``4 pi R^2`` (Ingersoll et al. 1992, single-scattering term):

    E_vis = A F sin(e) r^2 / (4 R^2),   E_IR = (1 - A) F sin(e) r^2 / (4 R^2)

the second with the walls in instantaneous radiative equilibrium and a constant albedo.
Pre-registered in docs/review/preregistration-th-a3-2026-11-16.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lunarenv.core.constants import ASTRONOMICAL_UNIT_M, SOLAR_CONSTANT_W_M2
from lunarenv.data.dem import Dem, GridGeometry
from lunarenv.geometry.projection import FloatArray, SouthPolarStereographic
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_series
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack
from lunarenv.thermal import terrain_radiation as tr
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.surface import surface_normal

RIM_RADIUS_M = 500.0
DEPTH_M = 200.0
PIXEL_M = 2.5
HALF_WIDTH_M = 800.0
CENTRE_XY_M = (0.0, -60000.0)
SUN_ELEVATION_DEG = 10.0
ALBEDO = 0.12
RECEIVER_OFFSETS_M = (0.0, -150.0, 150.0)  # along map x, across the Sun direction
SUN_RADIUS_M = 6.957e8


def sphere_radius_m(rim_radius_m: float = RIM_RADIUS_M, depth_m: float = DEPTH_M) -> float:
    return (rim_radius_m**2 + depth_m**2) / (2.0 * depth_m)


def bowl_stack(
    rim_radius_m: float = RIM_RADIUS_M, depth_m: float = DEPTH_M, pixel_m: float = PIXEL_M
) -> TerrainStack:
    """One-level terrain stack: level ground with the bowl at :data:`CENTRE_XY_M`."""
    projection = SouthPolarStereographic()
    n = round(2 * HALF_WIDTH_M / pixel_m)
    x0, y0 = CENTRE_XY_M
    grid = GridGeometry(x0 - HALF_WIDTH_M, y0 + HALF_WIDTH_M, pixel_m, n, n, projection)
    xc = grid.x_centers() - x0
    yc = grid.y_centers() - y0
    rho2 = xc[None, :] ** 2 + yc[:, None] ** 2
    radius = sphere_radius_m(rim_radius_m, depth_m)
    centre_height = radius - depth_m  # sphere centre above the floor, relative to the rim
    inside = rho2 < rim_radius_m**2
    height = np.where(inside, centre_height - np.sqrt(np.maximum(radius**2 - rho2, 0.0)), 0.0)
    dem = Dem(height.astype(np.float32), grid, "synthetic-bowl")
    return TerrainStack((TerrainLevel(dem, 0.0, HALF_WIDTH_M * 1.5, pixel_m / 2.0),))


def sun_positions(n_times: int = 3, elevation_deg: float = SUN_ELEVATION_DEG) -> FloatArray:
    """Fixed Sun at 1 au, ``elevation_deg`` above the horizontal of the bowl centre,
    grid azimuth 0, repeated ``n_times``."""
    projection = SouthPolarStereographic()
    _, e_y, up = projection.local_basis(*CENTRE_XY_M)
    e = np.radians(elevation_deg)
    centre = projection.to_cartesian(*CENTRE_XY_M, 0.0)
    position = centre + ASTRONOMICAL_UNIT_M * (np.cos(e) * e_y + np.sin(e) * up)
    return np.tile(position, (n_times, 1))


def analytic_irradiance(
    elevation_deg: float = SUN_ELEVATION_DEG, albedo_value: float = ALBEDO
) -> tuple[float, float]:
    """Single-bounce (E_vis, E_IR) at any shadowed interior point."""
    radius = sphere_radius_m()
    base = (
        SOLAR_CONSTANT_W_M2 * np.sin(np.radians(elevation_deg)) * RIM_RADIUS_M**2 / (4 * radius**2)
    )
    return albedo_value * base, (1.0 - albedo_value) * base


@dataclass(frozen=True, slots=True)
class BowlResult:
    offsets_m: tuple[float, ...]
    in_shadow: tuple[bool, ...]
    visible_w_m2: tuple[float, ...]
    infrared_w_m2: tuple[float, ...]
    reference_visible_w_m2: float
    reference_infrared_w_m2: float
    n_patches: int


def run_bowl(
    *,
    n_sectors: int,
    band_edges_m: tuple[float, ...],
    n_rays: int,
    patch_azimuths: int,
    points_per_patch: int = tr.DEFAULT_POINTS_PER_PATCH,
) -> BowlResult:
    """Single-bounce irradiance at the receivers with the given patch settings."""
    stack = bowl_stack()
    projection = stack.levels[0].dem.grid.projection
    sun = sun_positions()
    et = np.arange(3, dtype=np.float64) * 3600.0
    properties = RegolithProperties(albedo_normal=ALBEDO, albedo_a=0.0, albedo_b=0.0)
    x0, y0 = CENTRE_XY_M
    points = [(x0 + dx, y0) for dx in RECEIVER_OFFSETS_M]
    receivers = compute_horizons(
        stack,
        [p[0] for p in points],
        [p[1] for p in points],
        n_azimuth=patch_azimuths,
        diagnostics=False,
    )
    shadow, vis, ir, n_patches = [], [], [], 0
    for i, (x, y) in enumerate(points):
        # sectors and bands are relative to the receiver: one patch set per receiver
        view = tr.visible_terrain(stack, x, y, n_rays=n_rays)
        patches, (index,) = tr.build_patches(
            [view],
            n_sectors=n_sectors,
            band_edges_m=band_edges_m,
            points_per_patch=points_per_patch,
        )
        light = tr.patch_sunlight(stack, patches, et, sun, SUN_RADIUS_M, n_azimuth=patch_azimuths)
        exitance = tr.patch_exitance(
            light, properties, step_s=3600.0, spinup_samples=1, equilibrium=True
        )
        series = illumination_series(
            receivers, et, sun, SUN_RADIUS_M, observer=i, projection=projection
        )
        shadow.append(bool(series.visible_fraction[0] <= 0.0))
        normal = surface_normal(stack.levels[0].dem, x, y)
        e_vis, e_ir = tr.receiver_irradiance(view, index, exitance, normal)
        vis.append(float(e_vis[0]))
        ir.append(float(e_ir[0]))
        n_patches = max(n_patches, patches.size)
    ref_vis, ref_ir = analytic_irradiance()
    return BowlResult(
        RECEIVER_OFFSETS_M, tuple(shadow), tuple(vis), tuple(ir), ref_vis, ref_ir, n_patches
    )
