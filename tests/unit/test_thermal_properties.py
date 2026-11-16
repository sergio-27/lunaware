"""Regolith property laws: published values, internal consistency, grid."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from pydantic import ValidationError

from lunarenv.thermal.grid import geometric_grid, uniform_grid
from lunarenv.thermal.properties import (
    RegolithProperties,
    albedo,
    conductivity,
    conductivity_dt,
    contact_conductivity,
    density,
    enthalpy,
    heat_capacity,
    kirchhoff,
    radiative_equilibrium_temperature,
    temperature_from_enthalpy,
)

PROPS = RegolithProperties()
P = PROPS.as_array()
temperature = st.floats(20.0, 450.0)


def test_defaults_are_hayne2017_table_a1():
    assert (PROPS.surface_conductivity_w_m_k, PROPS.deep_conductivity_w_m_k) == (7.4e-4, 3.4e-3)
    assert (PROPS.surface_density_kg_m3, PROPS.deep_density_kg_m3) == (1100.0, 1800.0)
    assert (PROPS.scale_height_m, PROPS.radiative_parameter, PROPS.emissivity) == (0.06, 2.7, 0.95)
    assert PROPS.heat_capacity_coefficients == (-3.6125, 2.7431, 2.3616e-3, -1.2340e-5, 8.9093e-9)
    assert (PROPS.albedo_normal, PROPS.albedo_a, PROPS.albedo_b) == (0.12, 0.06, 0.25)
    assert PROPS.heat_flow_w_m2 == 0.018


def test_density_and_contact_conductivity_limits():
    assert density(0.0, P) == pytest.approx(1100.0)
    assert density(10.0, P) == pytest.approx(1800.0)
    assert density(0.06, P) == pytest.approx(1800.0 - 700.0 / np.e)
    assert contact_conductivity(1100.0, P) == pytest.approx(7.4e-4)
    assert contact_conductivity(1800.0, P) == pytest.approx(3.4e-3)
    assert contact_conductivity(1450.0, P) == pytest.approx(0.5 * (7.4e-4 + 3.4e-3))


def test_conductivity_at_reference_temperature():
    # at 350 K the radiative term equals chi times the contact term (Hayne2017, A4)
    assert conductivity(350.0, 7.4e-4, 2.7) == pytest.approx(7.4e-4 * 3.7)
    assert conductivity(0.0, 7.4e-4, 2.7) == 7.4e-4


def test_heat_capacity_polynomial_by_hand():
    t = 300.0
    by_hand = -3.6125 + 2.7431 * t + 2.3616e-3 * t**2 - 1.2340e-5 * t**3 + 8.9093e-9 * t**4
    assert heat_capacity(t, P) == pytest.approx(by_hand, rel=1e-13)
    assert 740.0 < by_hand < 780.0  # J kg-1 K-1, the familiar lunar-soil value near 300 K
    assert heat_capacity(1.0, P) < 0.0 < heat_capacity(2.0, P)  # negative below 1.3 K, as stated


@given(temperature)
def test_enthalpy_and_kirchhoff_are_the_integrals_they_claim(t):
    h = 1e-3
    d_enthalpy = (enthalpy(t + h, P) - enthalpy(t - h, P)) / (2 * h)
    assert d_enthalpy == pytest.approx(heat_capacity(t, P), rel=1e-7)
    d_phi = (kirchhoff(t + h, 1e-3, 2.7) - kirchhoff(t - h, 1e-3, 2.7)) / (2 * h)
    assert d_phi == pytest.approx(conductivity(t, 1e-3, 2.7), rel=1e-7)
    d_k = (conductivity(t + h, 1e-3, 2.7) - conductivity(t - h, 1e-3, 2.7)) / (2 * h)
    assert d_k == pytest.approx(conductivity_dt(t, 1e-3, 2.7), rel=1e-6)


@given(temperature, st.floats(20.0, 450.0))
def test_temperature_from_enthalpy_round_trip(t, guess):
    assert temperature_from_enthalpy(enthalpy(t, P), guess, P) == pytest.approx(t, abs=1e-8)


def test_albedo_law_values():
    assert albedo(0.0, PROPS) == pytest.approx(0.12)
    assert albedo(np.pi / 4, PROPS) == pytest.approx(0.12 + 0.06 + 0.25 / 256)
    assert albedo(np.pi / 2, PROPS) == pytest.approx(0.12 + 0.06 * 8 + 0.25)  # 0.85
    assert albedo(2.0, PROPS) == albedo(np.pi / 2, PROPS)  # clipped beyond grazing
    bright = RegolithProperties(albedo_normal=0.5)
    assert albedo(np.pi / 2, bright) == 1.0  # capped


def test_radiative_equilibrium_of_the_subsolar_point():
    # (1 - 0.12) * 1361 W m-2 at emissivity 0.95: the well-known ~387 K
    assert radiative_equilibrium_temperature(0.88 * 1361.0, 0.95) == pytest.approx(386.5, abs=0.5)
    assert radiative_equilibrium_temperature(-5.0, 0.95) == 0.0


def test_invalid_properties_are_rejected():
    with pytest.raises(ValidationError):
        RegolithProperties(surface_density_kg_m3=2000.0)  # denser at the top than at depth
    with pytest.raises(ValidationError):
        RegolithProperties(emissivity=1.2)
    with pytest.raises(ValidationError):
        RegolithProperties(bogus=1.0)


def test_grids():
    grid = geometric_grid(0.002, 1.2, 1.0)
    assert grid.edges_m[0] == 0.0 and grid.depth_m >= 1.0
    assert np.allclose(grid.thickness_m[1:] / grid.thickness_m[:-1], 1.2)
    assert np.allclose(grid.centers_m, grid.edges_m[:-1] + grid.thickness_m / 2)
    uniform = uniform_grid(0.01, 0.5)
    assert uniform.n_cells == 50 and np.allclose(uniform.thickness_m, 0.01)
    for bad in ((0.0, 1.2, 1.0), (0.01, 0.9, 1.0), (2.0, 1.2, 1.0)):
        with pytest.raises(ValueError):
            geometric_grid(*bad)


def test_woods_robinson_heat_capacity_table_and_enthalpy():
    wr = RegolithProperties(heat_capacity_model="woods_robinson2019").as_array()
    # DSNE (SLS-SPEC-159 Rev I) Table 3.4.6.2.2-1 rows, and the midpoint of 40-50 K
    assert heat_capacity(40.0, wr) == pytest.approx(89.0)
    assert heat_capacity(91.0, wr) == pytest.approx(258.0)
    assert heat_capacity(45.0, wr) == pytest.approx(0.5 * (89.0 + 122.0))
    # polynomial above 300 K, and the default law is unchanged
    assert heat_capacity(350.0, wr) == pytest.approx(heat_capacity(350.0, P), rel=1e-13)
    assert heat_capacity(43.0, P) == pytest.approx(117.8, abs=0.1)
    # enthalpy is the integral of c_p (continuous across 300 K) and inverts exactly
    for t in (15.0, 43.0, 95.0, 299.5, 300.5, 380.0):
        h = 1e-3
        derivative = (enthalpy(t + h, wr) - enthalpy(t - h, wr)) / (2 * h)
        assert derivative == pytest.approx(heat_capacity(t, wr), rel=1e-6)
        assert temperature_from_enthalpy(enthalpy(t, wr), 0.7 * t, wr) == pytest.approx(t)


def test_low_temperature_conductivity_branches() -> None:
    from lunarenv.thermal.properties import WR_K_LOW_K, WR_K_TOP_K, low_temperature_factor

    assert low_temperature_factor(300.0) == 1.0
    assert low_temperature_factor(WR_K_TOP_K - 1e-9) == pytest.approx(1.0, abs=1e-8)
    below = low_temperature_factor(WR_K_LOW_K - 1e-9)
    assert below == pytest.approx(low_temperature_factor(WR_K_LOW_K), abs=1e-8)
    # values fixed by the pre-registration (docs/review/preregistration-f09-2026-11-16.md)
    assert low_temperature_factor(20.0) == pytest.approx(0.152, abs=5e-4)
    assert low_temperature_factor(95.0) == pytest.approx(0.689, abs=5e-4)
    assert low_temperature_factor(0.0) == 0.0


@pytest.mark.parametrize("t", [5.0, 19.0, 21.0, 43.0, 100.0, 149.0, 151.0, 300.0])
def test_low_temperature_kirchhoff_is_the_integral(t: float) -> None:
    h = 1e-4
    d_phi = (kirchhoff(t + h, 1e-3, 2.7, 1.0) - kirchhoff(t - h, 1e-3, 2.7, 1.0)) / (2 * h)
    assert d_phi == pytest.approx(conductivity(t, 1e-3, 2.7, 1.0), rel=1e-7)
    d_k = (conductivity(t + h, 1e-3, 2.7, 1.0) - conductivity(t - h, 1e-3, 2.7, 1.0)) / (2 * h)
    assert d_k == pytest.approx(conductivity_dt(t, 1e-3, 2.7, 1.0), rel=1e-5)


def test_low_temperature_option_leaves_warm_conductivity_and_default_unchanged() -> None:
    for t in (150.0, 250.0, 400.0):
        assert conductivity(t, 1e-3, 2.7, 1.0) == conductivity(t, 1e-3, 2.7)
    assert conductivity(43.0, 1e-3, 2.7, 1.0) < 0.35 * conductivity(43.0, 1e-3, 2.7)
    p = RegolithProperties(conductivity_model="woods_robinson2019").as_array()
    assert p[-1] == 1.0 and RegolithProperties().as_array()[-1] == 0.0
