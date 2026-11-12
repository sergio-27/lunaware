"""Illumination products on synthetic geometry."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.core.constants import ASTRONOMICAL_UNIT_M, SOLAR_CONSTANT_W_M2
from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.geometry.visibility import horizon_dip
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import (
    illumination_series,
    illumination_statistics,
    longest_run,
)
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack
from tests.helpers import PROJ, analytic_dem

SUN_RADIUS_M = 695_700_000.0


def test_longest_run():
    step = 3600.0
    mask = np.array([1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1], dtype=bool)
    run = longest_run(mask, step)
    assert (run.count, run.start_index, run.duration_s) == (3, 3, 3 * step)  # first of the ties
    assert longest_run(~mask, step).count == 2
    assert longest_run(np.zeros(5, bool), step) == longest_run(np.zeros(1, bool), step)
    assert longest_run(np.zeros(5, bool), step).count == 0
    assert longest_run(np.ones(5, bool), step).count == 5


def body_at(observer_xyz, basis, elevation, azimuth, distance):
    e_x, e_y, up = basis
    horizontal = np.sin(azimuth)[:, None] * e_x + np.cos(azimuth)[:, None] * e_y
    direction = np.cos(elevation)[:, None] * horizontal + np.sin(elevation)[:, None] * up
    return observer_xyz + distance[:, None] * direction


@pytest.fixture(scope="module")
def flat_world():
    """Two observers 2 m and 100 m above a smooth sphere."""
    dem = analytic_dem(
        lambda x, y: np.zeros(np.shape(x)), center=(0.0, 0.0), half_width_m=30_000.0, pixel_m=20.0
    )
    stack = TerrainStack((TerrainLevel(dem, 0.0, 25_000.0, 10.0),))
    low = compute_horizons(stack, 1000.0, -2000.0, height_above_surface_m=2.0, n_azimuth=90)
    high = compute_horizons(stack, 1000.0, -2000.0, height_above_surface_m=100.0, n_azimuth=90)
    return low, high


def test_series_recovers_the_geometry_it_was_given(flat_world):
    low, _ = flat_world
    xy = (float(low.x_m[0]), float(low.y_m[0]))
    observer = PROJ.to_cartesian(*xy, low.observer_height_m[0])
    n = 400
    elevation = np.linspace(-0.02, 0.02, n)
    azimuth = np.linspace(0.1, 6.1, n)
    distance = np.full(n, 1.01 * ASTRONOMICAL_UNIT_M)
    sun = body_at(observer, PROJ.local_basis(*xy), elevation, azimuth, distance)
    et = 1000.0 + 3600.0 * np.arange(n)
    series = illumination_series(low, et, sun, SUN_RADIUS_M, n_strips=512)
    assert np.max(np.abs(series.elevation_rad - elevation)) < 1e-12
    assert np.max(np.abs(series.azimuth_rad - azimuth)) < 1e-10
    assert np.allclose(series.distance_m, distance, rtol=1e-12)
    rho = np.arcsin(SUN_RADIUS_M / distance)
    assert np.allclose(series.angular_radius_rad, rho)

    dip = float(horizon_dip(2.0, R))
    clear = elevation - rho > dip + 1e-6
    hidden = elevation + rho < dip - 1e-6
    assert np.allclose(series.visible_fraction[clear], 1.0) and np.all(
        series.visible_fraction[hidden] == 0.0
    )
    partial = ~clear & ~hidden
    y0 = (dip - elevation[partial]) / rho[partial]
    expected = (np.arccos(y0) - y0 * np.sqrt(1 - y0**2)) / np.pi
    assert np.max(np.abs(series.visible_fraction[partial] - expected)) < 5e-6
    # irradiance: inverse-square law times the flux fraction
    assert np.allclose(
        series.normal_irradiance_w_m2, SOLAR_CONSTANT_W_M2 / 1.01**2 * series.flux_fraction
    )
    assert series.step_s == 3600.0


def test_statistics_agree_with_the_series(flat_world):
    low, high = flat_world
    xy = (float(low.x_m[0]), float(low.y_m[0]))
    observer = PROJ.to_cartesian(*xy, low.observer_height_m[0])
    rng = np.random.default_rng(3)
    n = 500
    # a Sun wandering around the horizon so that lit, partial and dark samples all occur
    elevation = 0.012 * np.sin(np.linspace(0, 9 * np.pi, n)) + rng.normal(0, 5e-4, n) - 0.004
    sun = body_at(
        observer,
        PROJ.local_basis(*xy),
        elevation,
        np.linspace(0, 12, n) % (2 * np.pi),
        np.full(n, ASTRONOMICAL_UNIT_M),
    )
    et = 7200.0 * np.arange(n)
    for horizons in (low, high):
        series = illumination_series(horizons, et, sun, SUN_RADIUS_M)
        stats = illumination_statistics(horizons, et, sun, SUN_RADIUS_M)
        assert stats.mean_visible_fraction[0] == pytest.approx(
            series.mean_visible_fraction, abs=1e-14
        )
        assert stats.mean_flux_fraction[0] == pytest.approx(series.flux_fraction.mean(), abs=1e-14)
        assert stats.longest_shadow_s[0] == series.longest_shadow().duration_s
        assert stats.longest_illumination_s[0] == series.longest_illumination().duration_s
        assert 0.0 < series.mean_visible_fraction < 1.0 and series.longest_shadow().count > 1
    # a higher observer sees a lower horizon: never less illumination
    assert (
        illumination_statistics(high, et, sun, SUN_RADIUS_M).mean_visible_fraction[0]
        > illumination_statistics(low, et, sun, SUN_RADIUS_M).mean_visible_fraction[0]
    )


def test_time_axis_must_be_uniform(flat_world):
    low, _ = flat_world
    body = np.tile([0.0, 0.0, -ASTRONOMICAL_UNIT_M], (3, 1))
    with pytest.raises(ValueError, match="uniform"):
        illumination_series(low, np.array([0.0, 3600.0, 9000.0]), body, SUN_RADIUS_M)
    with pytest.raises(ValueError, match="shape"):
        illumination_series(low, np.array([0.0, 3600.0, 7200.0]), body[:2], SUN_RADIUS_M)
    # millisecond non-uniformity of a UTC grid expressed in TDB is accepted
    illumination_series(low, np.array([0.0, 3600.001, 7200.0]), body, SUN_RADIUS_M)
