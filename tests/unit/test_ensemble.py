"""Terrain-ensemble propagation and settings, on synthetic terrain."""

from __future__ import annotations

import numpy as np
import pytest
from pydantic import ValidationError

from lunarenv.core.constants import ASTRONOMICAL_UNIT_M
from lunarenv.data.dem import Dem
from lunarenv.illumination.disc import CHORD, PROFILE
from lunarenv.illumination.ensemble import QUANTILES, terrain_ensemble
from lunarenv.illumination.settings import IlluminationSettings
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack
from tests.helpers import PROJ, analytic_dem

SUN_RADIUS_M = 695_700_000.0


def ridge(height):
    """A wall of the given height 300 m north (map +y) of the observer at the origin."""
    return lambda x, y: height * np.exp(-((y - 300.0) ** 2) / (2.0 * 40.0**2))


def test_settings_defaults_validation_and_modes():
    settings = IlluminationSettings()
    assert settings.disc_mode_code == PROFILE
    assert IlluminationSettings(disc_mode="chord").disc_mode_code == CHORD
    with pytest.raises(ValidationError):
        IlluminationSettings(n_azimuth=2)
    with pytest.raises(ValidationError):
        IlluminationSettings(unknown=1)
    with pytest.raises(ValidationError):
        settings.n_azimuth = 10  # frozen


def test_higher_terrain_gives_less_illumination_and_quantiles_are_ordered():
    nominal = analytic_dem(ridge(20.0), center=(0.0, 0.0), half_width_m=1_000.0, pixel_m=5.0)
    stack = TerrainStack((TerrainLevel(nominal, 0.0, 900.0, 2.5),))
    heights = [10.0, 20.0, 30.0, 40.0]
    members = [
        (
            f"h{h:.0f}",
            Dem(
                analytic_dem(
                    ridge(h), center=(0.0, 0.0), half_width_m=1_000.0, pixel_m=5.0
                ).heights_m,
                nominal.grid,
            ),
        )
        for h in heights
    ]
    # the Sun sweeps up and down behind the ridge (grid azimuth 0 = map +y)
    n = 240
    elevation = np.linspace(0.0, 0.15, n)
    observer = PROJ.to_cartesian(0.0, 0.0, 2.0)
    _, e_y, up = PROJ.local_basis(0.0, 0.0)
    direction = np.cos(elevation)[:, None] * e_y + np.sin(elevation)[:, None] * up
    sun = observer + ASTRONOMICAL_UNIT_M * direction
    et = 3600.0 * np.arange(n)
    settings = IlluminationSettings(n_azimuth=360)
    ensemble = terrain_ensemble(
        stack,
        members,
        0.0,
        0.0,
        et,
        sun,
        SUN_RADIUS_M,
        height_above_surface_m=2.0,
        settings=settings,
    )
    assert ensemble.member_ids == ("h10", "h20", "h30", "h40")
    mean = ensemble.mean_visible_fraction[:, 0]
    assert np.all(np.diff(mean) < 0.0)  # taller ridge, less Sun
    assert np.all(np.diff(ensemble.longest_shadow_s[:, 0]) > 0.0)
    # ridge top seen from 300 m: elevation ~ (h - 2) / 300; the Sun clears it later
    expected_dark_fraction = np.array([(h - 2.0) / 300.0 for h in heights]) / 0.15
    assert np.allclose(1.0 - mean, expected_dark_fraction, atol=0.03)
    quantiles = ensemble.quantiles("mean_visible_fraction")
    assert quantiles.shape == (len(QUANTILES), 1) and np.all(np.diff(quantiles[:, 0]) >= 0.0)


def test_empty_ensemble_and_wrong_grid_are_rejected():
    nominal = analytic_dem(ridge(20.0), center=(0.0, 0.0), half_width_m=500.0, pixel_m=5.0)
    stack = TerrainStack((TerrainLevel(nominal, 0.0, 400.0, 2.5),))
    body = np.tile([0.0, 0.0, -ASTRONOMICAL_UNIT_M], (2, 1))
    et = np.array([0.0, 3600.0])
    with pytest.raises(ValueError, match="at least one"):
        terrain_ensemble(stack, [], 0.0, 0.0, et, body, SUN_RADIUS_M)
    other = analytic_dem(ridge(20.0), center=(0.0, 0.0), half_width_m=400.0, pixel_m=5.0)
    with pytest.raises(ValueError, match="same grid"):
        terrain_ensemble(stack, [("x", other)], 0.0, 0.0, et, body, SUN_RADIUS_M)


def test_regional_error_stacks_bound_the_horizon():
    from lunarenv.illumination.ensemble import regional_error_stacks, shifted
    from lunarenv.illumination.horizon import compute_horizons
    from lunarenv.illumination.terrain import build_stack

    def hills(x, y):
        return 300.0 * np.exp(-((x - 6_000.0) ** 2 + (y - 2_000.0) ** 2) / (2.0 * 1_500.0**2))

    def zero(x, y):
        return np.zeros(np.shape(x))

    near = analytic_dem(zero, center=(0.0, 0.0), half_width_m=600.0, pixel_m=5.0)
    regional = analytic_dem(hills, center=(0.0, 0.0), half_width_m=24_000.0, pixel_m=80.0)
    rng = np.random.default_rng(1)
    error = Dem(rng.uniform(0.0, 20.0, regional.heights_m.shape).astype(np.float32), regional.grid)
    radii = {"near_radius_m": 500.0, "mid_radius_m": 10_000.0}
    settings = IlluminationSettings(n_azimuth=64, **radii)
    nominal = build_stack(near, regional, **radii)
    lowered, raised = regional_error_stacks(near, regional, error, 2.0, settings)
    # any realization within +/- 2 errors, here a random one, lies between the bounds
    noise = rng.uniform(-2.0, 2.0, regional.heights_m.shape) * error.heights_m
    random = build_stack(near, shifted(regional, noise, 1.0), **radii)

    def horizon(stack):
        return compute_horizons(stack, 0.0, 0.0, n_azimuth=64, diagnostics=False).elevation_rad[0]

    low, mid, high, rnd = (horizon(s) for s in (lowered, nominal, raised, random))
    assert np.all(low <= mid) and np.all(mid <= high)
    assert np.all(low <= rnd) and np.all(rnd <= high)
    assert np.any(high > low)
    with pytest.raises(ValueError, match="grid"):
        regional_error_stacks(near, regional, near, 2.0, settings)
    with pytest.raises(ValueError, match="positive"):
        regional_error_stacks(near, regional, error, 0.0, settings)
