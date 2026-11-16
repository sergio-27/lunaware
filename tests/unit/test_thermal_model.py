"""Surface energy input, eclipse factor, spin-up and summaries."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.core.constants import ASTRONOMICAL_UNIT_M, SOLAR_CONSTANT_W_M2
from lunarenv.illumination.eclipse import disc_overlap_area, earth_occultation_factor
from lunarenv.illumination.products import IlluminationSeries
from lunarenv.thermal.model import (
    SYNODIC_MONTH_S,
    ThermalSettings,
    resistance_to_cells,
    run_thermal,
    spin_up,
    summarize,
)
from lunarenv.thermal.properties import RegolithProperties, kirchhoff
from lunarenv.thermal.surface import (
    absorbed_solar_flux,
    analytic_absorbed_flux,
    incidence_cosine,
    neglected_terrain_flux_bound,
    normal_from_slope,
    slope_and_facing,
    surface_normal,
)
from tests.helpers import analytic_dem

PROPS = RegolithProperties()
PER = 360
STEP = SYNODIC_MONTH_S / PER
HOUR = 2 * np.pi * np.arange(2 * PER + 1) / PER
EQUATOR = analytic_absorbed_flux(HOUR, 0.0, PROPS)


# --------------------------------------------------------------------------- eclipse
def test_disc_overlap_limits_and_closed_form():
    assert disc_overlap_area(3.0, 1.0, 1.5) == 0.0  # apart
    assert disc_overlap_area(0.2, 1.0, 3.0) == pytest.approx(np.pi)  # small disc inside
    # two unit discs one radius apart: lens area 2 pi / 3 - sqrt(3) / 2
    assert disc_overlap_area(1.0, 1.0, 1.0) == pytest.approx(2 * np.pi / 3 - np.sqrt(3) / 2)
    separation = np.linspace(0.0, 5.0, 200)
    area = disc_overlap_area(separation, 1.0, 3.5)
    assert np.all(np.diff(area) <= 1e-12) and area[0] == pytest.approx(np.pi) and area[-1] == 0.0


def test_earth_occultation_factor_geometry():
    observer = np.zeros(3)
    sun = np.array([[ASTRONOMICAL_UNIT_M, 0.0, 0.0]] * 3)
    r_sun, r_earth, d_earth = 6.957e8, 6.371e6, 3.844e8
    offsets = np.array([0.0, np.arcsin(r_earth / d_earth), 0.05])  # centred, limb on centre, clear
    earth = d_earth * np.column_stack([np.cos(offsets), np.sin(offsets), np.zeros(3)])
    factor = earth_occultation_factor(sun, earth, observer, r_sun, r_earth)
    assert factor[0] == 0.0  # Earth's disc is 3.6 times wider than the Sun's
    assert 0.4 < factor[1] < 0.6  # limb through the Sun's centre: about half
    assert factor[2] == 1.0


# --------------------------------------------------------------------------- surface
def test_normals_slopes_and_incidence():
    normal = normal_from_slope(np.radians(20.0), np.radians(90.0))  # tilted towards map +x
    assert np.allclose(normal, [np.sin(np.radians(20.0)), 0.0, np.cos(np.radians(20.0))])
    slope, facing = slope_and_facing(normal)
    assert slope == pytest.approx(np.radians(20.0)) and facing == pytest.approx(np.radians(90.0))
    # Sun at the zenith: incidence equals the slope
    assert incidence_cosine(normal, np.pi / 2, 0.0) == pytest.approx(np.cos(np.radians(20.0)))
    # Sun on the horizon in the direction the surface faces
    assert incidence_cosine(normal, 0.0, np.radians(90.0)) == pytest.approx(
        np.sin(np.radians(20.0))
    )
    assert incidence_cosine(normal, 0.0, np.radians(270.0)) < 0.0  # behind the slope
    with pytest.raises(ValueError):
        normal_from_slope(np.radians(95.0), 0.0)


def test_surface_normal_of_a_tilted_plane():
    dem = analytic_dem(
        lambda x, y: 0.1 * x - 0.05 * y, center=(0.0, 0.0), half_width_m=200.0, pixel_m=5.0
    )
    expected = np.array([-0.1, 0.05, 1.0]) / np.linalg.norm([-0.1, 0.05, 1.0])
    assert np.allclose(surface_normal(dem, 12.5, -7.5), expected, atol=1e-6)
    assert np.allclose(surface_normal(dem, 12.5, -7.5, baseline_m=20.0), expected, atol=1e-6)


def series_for(elevation, azimuth, flux_fraction=1.0, distance_au=1.0):
    n = np.size(elevation)
    full = np.full(n, 1.0)
    return IlluminationSeries(
        et_s=np.arange(n) * 3600.0,
        step_s=3600.0,
        elevation_rad=np.asarray(elevation, float),
        azimuth_rad=np.asarray(azimuth, float),
        angular_radius_rad=0.00465 * full,
        distance_m=distance_au * ASTRONOMICAL_UNIT_M * full,
        visible_fraction=flux_fraction * full,
        flux_fraction=flux_fraction * full,
    )


def test_absorbed_flux_on_level_and_tilted_ground():
    level = np.array([0.0, 0.0, 1.0])
    overhead = series_for([np.pi / 2], [0.0])
    assert absorbed_solar_flux(overhead, level, PROPS)[0] == pytest.approx(
        0.88 * SOLAR_CONSTANT_W_M2
    )
    # grazing Sun: level ground absorbs almost nothing, a slope facing the Sun a lot
    low = series_for([np.radians(1.0)], [np.radians(90.0)])
    facing = normal_from_slope(np.radians(30.0), np.radians(90.0))
    away = normal_from_slope(np.radians(30.0), np.radians(270.0))
    q_level = absorbed_solar_flux(low, level, PROPS)[0]
    q_facing = absorbed_solar_flux(low, facing, PROPS)[0]
    assert 0.0 < q_level < 5.0 and q_facing > 300.0
    assert absorbed_solar_flux(low, away, PROPS)[0] == 0.0
    # distance, partial disc and eclipse scale the flux
    half = series_for([np.radians(1.0)], [np.radians(90.0)], flux_fraction=0.5, distance_au=2.0)
    assert absorbed_solar_flux(half, facing, PROPS)[0] == pytest.approx(q_facing / 8.0)
    assert absorbed_solar_flux(low, facing, PROPS, eclipse_factor=[0.25])[0] == pytest.approx(
        q_facing / 4
    )


def test_analytic_flux_and_neglected_flux_bound():
    assert EQUATOR[0] == pytest.approx(0.88 * SOLAR_CONSTANT_W_M2)
    assert np.all(EQUATOR[PER // 4 : 3 * PER // 4] < 1e-9)  # night (round-off at the terminator)
    assert neglected_terrain_flux_bound(np.full(360, -0.01), PROPS).flux_w_m2 == 0.0
    walls = neglected_terrain_flux_bound(np.full(360, np.radians(30.0)), PROPS)
    assert walls.terrain_view_factor == pytest.approx(0.25)
    assert walls.flux_w_m2 == pytest.approx(0.25 * SOLAR_CONSTANT_W_M2 / 0.98**2)
    assert 250.0 < walls.shadow_temperature_floor_k < 300.0


# ------------------------------------------------------------------------------ model
def test_spin_up_does_not_depend_on_the_initial_state():
    settings = ThermalSettings()
    cold = spin_up(EQUATOR[: PER + 1], STEP, PROPS, settings, initial_temperature_k=120.0)
    hot = spin_up(EQUATOR[: PER + 1], STEP, PROPS, settings, initial_temperature_k=380.0)
    assert cold.converged and hot.converged and max(cold.cycles, hot.cycles) < 30
    assert np.max(np.abs(cold.temperature_k - hot.temperature_k)) < 0.01


def test_periodic_state_closes_and_satisfies_the_mean_identity():
    """End of a cycle equals its start; mean psi at depth = mean psi at the surface + Q R."""
    settings = ThermalSettings()
    result = run_thermal(EQUATOR, STEP, PER, PROPS, settings)
    assert result.warnings == ()
    assert np.max(np.abs(result.temperature_k[-1] - result.temperature_k[0])) < 5e-3
    chi = PROPS.radiative_parameter
    surface = kirchhoff(result.surface_temperature_k[1:], 1.0, chi).mean()
    cells = kirchhoff(result.temperature_k[1:], 1.0, chi).mean(axis=0)
    expected = surface + PROPS.heat_flow_w_m2 * resistance_to_cells(settings.grid(), PROPS)
    assert np.max(np.abs(cells - expected) / expected) < 2e-5
    # the heat pumped down by the T^3 conductivity makes the interior warmer than the mean surface
    assert result.temperature_k[:, -1].mean() - result.surface_temperature_k[:-1].mean() > 35.0


def test_equatorial_cycle_reproduces_published_temperatures():
    """Pre-check of validation case TH-V1 (Hayne2017 Table A2, after Vasavada2012):
    about 385 K at noon, 101 K at midnight, 95 K before sunrise, for A0 = 0.12."""
    fine = 2 * np.pi * np.arange(2 * 720 + 1) / 720
    result = run_thermal(
        analytic_absorbed_flux(fine, 0.0, PROPS), SYNODIC_MONTH_S / 720, 720, PROPS
    )
    surface = result.surface_temperature_k
    assert surface.max() == pytest.approx(385.0, abs=1.0)
    assert surface[360] == pytest.approx(101.0, abs=1.5)
    assert surface.min() == pytest.approx(95.0, abs=2.0)


def test_summary_and_depth_interpolation():
    result = run_thermal(EQUATOR, STEP, PER, PROPS)
    summary = summarize(result)
    surface = result.surface_temperature_k
    assert summary.surface_max_k == surface.max() and summary.surface_min_k == surface.min()
    assert summary.lunation_amplitude_max_k == pytest.approx(surface.max() - surface.min(), abs=1.0)
    assert summary.max_heating_rate_k_per_h > summary.max_cooling_rate_k_per_h > 0.0
    assert np.all(np.diff(summary.depth_max_k - summary.depth_min_k) < 0.0)  # damped with depth
    assert np.array_equal(result.at_depth(0.0), surface)
    centre = result.depth_m[5]
    assert np.allclose(result.at_depth(centre), result.temperature_k[:, 5])
    with pytest.raises(ValueError):
        result.at_depth(99.0)


def test_validity_warnings():
    dim = 0.02 * EQUATOR  # a surface that never gets warm
    shadow = dim[PER:] == 0.0
    result = run_thermal(dim, STEP, PER, PROPS, latitude_rad=np.radians(-89.7), shadow=shadow)
    text = " | ".join(result.warnings)
    assert "heat-capacity fit" in text and "extrapolation" in text and "no direct sunlight" in text
    quick = ThermalSettings(spinup_max_cycles=1)
    assert any("spin-up" in w for w in run_thermal(EQUATOR, STEP, PER, PROPS, quick).warnings)
    with pytest.raises(ValueError):
        run_thermal(EQUATOR, STEP, 1, PROPS)


def test_woods_robinson_heat_capacity_reaches_the_solver():
    # cold ground, where the two c_p laws differ; also guards against a stale numba cache
    # (compiled solver kernels do not recompile when properties.py changes)
    dim = 0.02 * EQUATOR
    hayne = run_thermal(dim, STEP, PER, PROPS)
    wr_props = PROPS.model_copy(update={"heat_capacity_model": "woods_robinson2019"})
    wr = run_thermal(dim, STEP, PER, wr_props)
    assert hayne.surface_temperature_k.min() < 90.0
    assert abs(wr.surface_temperature_k.min() - hayne.surface_temperature_k.min()) > 0.01
