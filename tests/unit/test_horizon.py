"""Horizon tracer against closed-form and independently computed horizons.

ILL-A1: smooth sphere (horizon dip). ILL-A2: analytic hill, off-pole observer.
"""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.geometry.visibility import horizon_dip
from lunarenv.illumination.horizon import compute_horizons, sun_relevant_incomplete
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack, build_stack, coarsen
from tests.helpers import analytic_dem, reference_horizon


def flat(level_m: float):
    return lambda x, y: np.full(np.shape(x), level_m)


def single_level(dem, reach_m, step_m):
    return TerrainStack((TerrainLevel(dem, 0.0, reach_m, step_m),))


@pytest.mark.parametrize(("terrain_m", "above_m"), [(0.0, 2.0), (1500.0, 10.0), (-2500.0, 50.0)])
def test_ill_a1_smooth_sphere_gives_the_horizon_dip(terrain_m, above_m):
    """Flat terrain at constant height: the horizon is the dip of a sphere, every azimuth."""
    center = (4000.0, -6000.0)
    dem = analytic_dem(flat(terrain_m), center=center, half_width_m=20_000.0, pixel_m=20.0)
    stack = single_level(dem, 19_000.0, 10.0)
    horizons = compute_horizons(
        stack, center[0], center[1], height_above_surface_m=above_m, n_azimuth=36
    )
    expected = float(horizon_dip(above_m, R + terrain_m))
    assert np.max(np.abs(horizons.elevation_rad[0] - expected)) < 1e-7
    assert horizons.complete.all()
    # ground distances are measured on the reference sphere of radius R
    tangent_distance = R * np.arccos((R + terrain_m) / (R + terrain_m + above_m))
    assert np.max(np.abs(horizons.distance_m[0] - tangent_distance)) < 10.0  # one step


def hill(x, y):
    return 400.0 * np.exp(-((x - 3800.0) ** 2 + (y + 4500.0) ** 2) / (2.0 * 600.0**2))


def test_ill_a2_analytic_hill_matches_independent_reference():
    """Off-pole observer next to a Gaussian hill; reference via true-azimuth geodesics."""
    observer = (3000.0, -5000.0)
    dem = analytic_dem(hill, center=(3500.0, -4500.0), half_width_m=6_000.0, pixel_m=5.0)
    stack = single_level(dem, 5_000.0, 2.5)
    n_azimuth = 72
    horizons = compute_horizons(stack, *observer, height_above_surface_m=2.0, n_azimuth=n_azimuth)
    reference = reference_horizon(
        hill,
        observer,
        float(horizons.observer_height_m[0]),
        horizons.azimuth_rad,
        max_distance_m=5_000.0,
        step_m=0.5,
    )
    error = np.abs(horizons.elevation_rad[0] - reference)
    # bilinear interpolation of a 5 m raster: height error <= px^2 |z''| / 8 = 3.5 mm
    assert error.max() < 2e-5
    # the hill is in the direction of +x, +y: that is where the horizon is highest
    peak_azimuth = horizons.azimuth_rad[np.argmax(horizons.elevation_rad[0])]
    assert peak_azimuth == pytest.approx(np.arctan2(800.0, 500.0), abs=np.radians(5.0))


def test_local_slope_defines_the_horizon_of_a_surface_observer():
    slope = 0.1
    dem = analytic_dem(lambda x, y: slope * x, center=(0.0, 0.0), half_width_m=2_000.0, pixel_m=5.0)
    stack = single_level(dem, 1_500.0, 2.5)
    horizons = compute_horizons(stack, 2.5, 2.5, height_above_surface_m=0.0, n_azimuth=8)
    elevation = horizons.elevation_rad[0]
    assert elevation[2] == pytest.approx(np.arctan(slope), abs=1e-5)  # grid azimuth 90 deg = +x
    assert elevation[6] == pytest.approx(-np.arctan(slope), abs=1e-5)  # downslope
    assert elevation[0] == pytest.approx(0.0, abs=1e-5)  # along the contour


def test_nested_levels_agree_with_a_single_fine_level():
    observer = (3000.0, -5000.0)
    fine = analytic_dem(hill, center=(3500.0, -4500.0), half_width_m=6_000.0, pixel_m=5.0)
    reference = compute_horizons(
        single_level(fine, 5_000.0, 2.5), *observer, n_azimuth=72, height_above_surface_m=2.0
    )
    near = fine.crop(type(fine.grid.bounds)(2500.0, -5500.0, 3500.0, -4500.0))
    stack = TerrainStack(
        (
            TerrainLevel(near, 0.0, 400.0, 2.5),
            TerrainLevel(coarsen(fine, 4), 400.0, 5_000.0, 10.0),
        )
    )
    nested = compute_horizons(stack, *observer, n_azimuth=72, height_above_surface_m=2.0)
    assert np.max(np.abs(nested.elevation_rad - reference.elevation_rad)) < 2e-4


def test_rays_that_run_out_of_data_are_reported_incomplete():
    """A 1 km patch cannot contain the horizon of an observer 50 m up (tangent at 13 km)."""
    dem = analytic_dem(flat(0.0), center=(0.0, 0.0), half_width_m=500.0, pixel_m=5.0)
    stack = single_level(dem, 5_000.0, 2.5)
    horizons = compute_horizons(stack, 0.0, 0.0, height_above_surface_m=50.0, n_azimuth=16)
    assert not horizons.complete.any()
    assert np.all(horizons.stop_distance_m < 800.0)
    assert sun_relevant_incomplete(horizons, stack).all()


def test_build_stack_layout_and_sample_plan():
    near = analytic_dem(flat(0.0), center=(0.0, 0.0), half_width_m=1_000.0, pixel_m=5.0)
    regional = analytic_dem(flat(0.0), center=(0.0, 0.0), half_width_m=48_000.0, pixel_m=80.0)
    stack = build_stack(
        near, regional, near_radius_m=800.0, mid_radius_m=20_000.0, far_coarsening=3
    )
    assert [level.dem.grid.pixel_size_m for level in stack.levels] == [5.0, 80.0, 240.0]
    distance, level = stack.sample_plan()
    assert np.all(np.diff(distance) > 0) and distance[0] == pytest.approx(1.25)
    assert list(np.unique(level)) == [0, 1, 2]
    assert distance[level == 1].min() > 800.0 and distance[level == 0].max() < 800.0
    with pytest.raises(ValueError):
        build_stack(near, regional, near_radius_m=800.0, mid_radius_m=500.0)


def test_coarsen_is_a_block_mean():
    dem = analytic_dem(lambda x, y: x + 2 * y, center=(0.0, 0.0), half_width_m=40.0, pixel_m=5.0)
    coarse = coarsen(dem, 4)
    assert coarse.grid.pixel_size_m == 20.0 and coarse.grid.shape == (4, 4)
    xx, yy = np.meshgrid(coarse.grid.x_centers(), coarse.grid.y_centers())
    assert np.allclose(coarse.heights_m, xx + 2 * yy)  # the mean of a plane is its centre value


def test_invalid_requests_are_rejected():
    dem = analytic_dem(flat(0.0), center=(0.0, 0.0), half_width_m=500.0, pixel_m=5.0)
    stack = single_level(dem, 400.0, 2.5)
    with pytest.raises(ValueError):
        compute_horizons(stack, 0.0, 0.0, n_azimuth=4)
    with pytest.raises(ValueError):
        compute_horizons(stack, 0.0, 0.0, height_above_surface_m=-1.0)
    with pytest.raises(ValueError, match="outside"):
        compute_horizons(stack, 5_000.0, 0.0)


def test_surface_height_override_holds_the_observer_fixed():
    """With the same absolute observer height, raising the whole terrain lowers nothing nearby."""
    dem = analytic_dem(flat(100.0), center=(0.0, 0.0), half_width_m=20_000.0, pixel_m=20.0)
    stack = single_level(dem, 19_000.0, 10.0)
    default = compute_horizons(stack, 0.0, 0.0, height_above_surface_m=5.0, n_azimuth=8)
    forced = compute_horizons(
        stack, 0.0, 0.0, height_above_surface_m=5.0, n_azimuth=8, surface_height_m=100.0
    )
    assert np.array_equal(default.elevation_rad, forced.elevation_rad)
    raised = compute_horizons(
        stack, 0.0, 0.0, height_above_surface_m=5.0, n_azimuth=8, surface_height_m=150.0
    )
    assert raised.observer_height_m[0] == 155.0
    assert np.all(raised.elevation_rad < default.elevation_rad)  # deeper dip from higher up
