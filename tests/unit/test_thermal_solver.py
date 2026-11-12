"""Conduction solver: analytic solutions (TH-A1, TH-A2) and physical properties (TH-P1)."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4 as SIGMA
from lunarenv.thermal.grid import geometric_grid, uniform_grid
from lunarenv.thermal.properties import RegolithProperties, kirchhoff
from lunarenv.thermal.solver import SolverError, explicit_stable_step_s, solve_conduction

SCHEMES = ("explicit", "backward_euler", "crank_nicolson", "bdf2")

# constant-property medium for the analytic cases
K, RHO, CP = 1.0e-3, 1500.0, 600.0
KAPPA = K / (RHO * CP)
CONSTANT = RegolithProperties(
    surface_conductivity_w_m_k=K,
    deep_conductivity_w_m_k=K,
    surface_density_kg_m3=RHO,
    deep_density_kg_m3=RHO,
    radiative_parameter=0.0,
    heat_capacity_coefficients=(CP, 0.0, 0.0, 0.0, 0.0),
    heat_flow_w_m2=0.0,
)
PERIOD = 2.55024e6  # synodic month, s
OMEGA = 2 * np.pi / PERIOD
SKIN = np.sqrt(2 * KAPPA / OMEGA)  # 3.0 cm
MEAN, AMPLITUDE = 250.0, 100.0


def wave(z, t):
    """Carslaw & Jaeger: damped, phase-lagged wave under a periodic surface temperature."""
    return MEAN + AMPLITUDE * np.exp(-z / SKIN) * np.cos(OMEGA * t - z / SKIN)


def wave_error(scheme, cells_per_skin, steps_per_period, substeps=None):
    grid = uniform_grid(SKIN / cells_per_skin, 12 * SKIN)
    times = np.arange(steps_per_period + 1) * PERIOD / steps_per_period
    result = solve_conduction(
        grid,
        CONSTANT,
        wave(grid.centers_m, 0.0),
        wave(0.0, times),
        PERIOD / steps_per_period,
        scheme=scheme,
        boundary="temperature",
        substeps=substeps,
    )
    exact = wave(grid.centers_m[None, :], times[:, None])
    return float(np.abs(result.temperature_k - exact).max())


@pytest.mark.parametrize("scheme", SCHEMES)
def test_th_a1_periodic_wave_is_reproduced(scheme):
    assert wave_error(scheme, 20, 2000) < 0.25  # K, of a 100 K amplitude


def test_th_a1_spatial_convergence_is_second_order():
    errors = [wave_error("bdf2", cells, 20000) for cells in (5, 10, 20)]
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(orders > 1.8), (errors, orders)


@pytest.mark.parametrize(
    ("scheme", "expected"), [("backward_euler", 1.0), ("crank_nicolson", 2.0), ("bdf2", 2.0)]
)
def test_th_a1_temporal_convergence_order(scheme, expected):
    """Error against a spatially identical reference with a much smaller time step."""
    grid = uniform_grid(SKIN / 10, 12 * SKIN)

    def run(steps, name=scheme):
        times = np.arange(steps + 1) * PERIOD / steps
        result = solve_conduction(
            grid,
            CONSTANT,
            wave(grid.centers_m, 0.0),
            wave(0.0, times),
            PERIOD / steps,
            scheme=name,
            boundary="temperature",
        )
        return result.temperature_k[-1]

    reference = run(12800, "bdf2")
    errors = [float(np.abs(run(steps) - reference).max()) for steps in (100, 200, 400)]
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all(np.abs(orders - expected) < 0.25), (errors, orders)


def test_th_a2_steady_state_with_constant_properties_is_exact():
    """eps sigma Ts^4 = Q_abs + Q and a linear profile of gradient Q / K."""
    q_abs, heat_flow = 60.0, 0.5
    props = CONSTANT.model_copy(update={"heat_flow_w_m2": heat_flow})
    grid = geometric_grid(0.004, 1.15, 0.6)
    result = solve_conduction(
        grid,
        props,
        np.full(grid.n_cells, 150.0),
        np.full(400, q_abs),
        1.0e9,
        scheme="backward_euler",
    )
    surface = ((q_abs + heat_flow) / (props.emissivity * SIGMA)) ** 0.25
    assert result.surface_temperature_k[-1] == pytest.approx(surface, abs=1e-7)
    expected = surface + heat_flow * grid.centers_m / K
    assert np.max(np.abs(result.temperature_k[-1] - expected)) < 1e-6


def test_th_a2_steady_state_with_temperature_dependent_conductivity_is_exact():
    """With K(T) at uniform density the Kirchhoff potential is linear in depth.

    The fluxes are written in the Kirchhoff variable, so the discrete steady state
    reproduces this on any grid, not just in the limit of fine cells.
    """
    q_abs, heat_flow = 60.0, 2.0
    props = RegolithProperties(
        surface_conductivity_w_m_k=K,
        deep_conductivity_w_m_k=K,
        surface_density_kg_m3=RHO,
        deep_density_kg_m3=RHO,
        radiative_parameter=2.7,
        heat_flow_w_m2=heat_flow,
    )
    surface = ((q_abs + heat_flow) / (props.emissivity * SIGMA)) ** 0.25
    for grid in (uniform_grid(0.05, 0.4), geometric_grid(0.004, 1.3, 0.4)):
        result = solve_conduction(
            grid,
            props,
            np.full(grid.n_cells, 200.0),
            np.full(400, q_abs),
            1.0e9,
            scheme="backward_euler",
        )
        assert result.surface_temperature_k[-1] == pytest.approx(surface, abs=1e-6)
        potential = kirchhoff(result.temperature_k[-1], K, 2.7)
        expected = kirchhoff(surface, K, 2.7) + heat_flow * grid.centers_m
        assert np.max(np.abs(potential - expected)) < 1e-9
        assert result.temperature_k[-1, -1] - surface > 100.0  # a strongly nonlinear case


def test_th_a2_steady_state_with_the_full_property_model_converges_at_second_order():
    """Depth-dependent density: psi(z) = psi_s + Q * integral of dz / K_c(z).

    The only discretization error left is the sampling of K_c at cell centres.
    """
    heat_flow, q_abs = 0.5, 30.0
    props = RegolithProperties(heat_flow_w_m2=heat_flow)
    p = props.as_array()
    surface = ((q_abs + heat_flow) / (props.emissivity * SIGMA)) ** 0.25

    def kc(z):
        rho = 1800.0 - 700.0 * np.exp(-z / 0.06)
        return 3.4e-3 - (3.4e-3 - 7.4e-4) * (1800.0 - rho) / 700.0

    errors = []
    for cell in (0.01, 0.005, 0.0025):
        grid = uniform_grid(cell, 0.4)
        result = solve_conduction(
            grid,
            props,
            np.full(grid.n_cells, 150.0),
            np.full(600, q_abs),
            1.0e9,
            scheme="backward_euler",
        )
        fine = np.linspace(0.0, 0.4, 400001)
        resistance = np.concatenate(
            ([0.0], np.cumsum(0.5 * (1 / kc(fine[1:]) + 1 / kc(fine[:-1])) * np.diff(fine)))
        )
        exact_psi = kirchhoff(surface, 1.0, p[5]) + heat_flow * np.interp(
            grid.centers_m, fine, resistance
        )
        errors.append(
            np.max(np.abs(kirchhoff(result.temperature_k[-1], 1.0, p[5]) - exact_psi) / exact_psi)
        )
    # the density scale height is 6 cm, so cells of 1 cm are barely in the asymptotic
    # range: the observed order rises towards 2 as the grid is refined
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert errors[-1] < 3e-4
    assert orders[0] > 1.6 and orders[1] > 1.8 and orders[1] > orders[0]


def test_th_a2_vanishing_conductivity_gives_radiative_equilibrium():
    props = CONSTANT.model_copy(
        update={"surface_conductivity_w_m_k": 1e-12, "deep_conductivity_w_m_k": 1e-12}
    )
    grid = uniform_grid(0.01, 0.2)
    flux = 50.0 + 400.0 * np.sin(np.linspace(0, np.pi, 50)) ** 2
    result = solve_conduction(
        grid, props, np.full(grid.n_cells, 200.0), flux, 3600.0, scheme="bdf2"
    )
    assert np.allclose(result.surface_temperature_k, (flux / (0.95 * SIGMA)) ** 0.25, atol=1e-4)
    assert np.allclose(result.temperature_k, 200.0, atol=1e-4)  # nothing reaches the interior


# ------------------------------------------------------------------ realistic regolith
HAYNE = RegolithProperties()
GRID = geometric_grid(0.003, 1.15, 1.2)


def day_night_flux(n_cycles=3, per_cycle=240, peak=1100.0):
    phase = 2 * np.pi * np.arange(n_cycles * per_cycle + 1) / per_cycle
    return peak * np.clip(np.cos(phase), 0.0, None), PERIOD / per_cycle


@pytest.mark.parametrize("scheme", SCHEMES)
def test_th_p1_energy_is_conserved(scheme):
    """Stored enthalpy changes by the time-integrated net boundary flux."""
    flux, step = day_night_flux()
    result = solve_conduction(GRID, HAYNE, np.full(GRID.n_cells, 240.0), flux, step, scheme=scheme)
    throughput = float(np.sum(flux) * step)  # J m-2 absorbed over the run
    assert abs(result.energy_residual_j_m2) < 1e-9 * throughput
    assert abs(result.stored_energy_change_j_m2) > 1e-4 * throughput  # the test is not vacuous


@pytest.mark.parametrize("scheme", ["explicit", "backward_euler", "bdf2"])
def test_th_p1_temperatures_stay_within_physical_bounds(scheme):
    flux, step = day_night_flux()
    initial = 240.0
    result = solve_conduction(
        GRID, HAYNE, np.full(GRID.n_cells, initial), flux, step, scheme=scheme
    )
    hottest = (flux.max() / (HAYNE.emissivity * SIGMA)) ** 0.25
    # with the 18 mW m-2 of interior heat the bound moves by a negligible amount
    assert result.surface_temperature_k.max() <= hottest + 0.01
    assert result.temperature_k.max() <= max(hottest, initial) + 0.01
    assert result.temperature_k.min() > 40.0 and result.surface_temperature_k.min() > 40.0
    # the surface cools monotonically through the night
    night = result.surface_temperature_k[61:180]
    assert np.all(np.diff(night) < 0.0)


def test_crank_nicolson_overshoots_after_an_abrupt_change():
    """Why Crank-Nicolson is not the default (PLAN.md decision D9).

    The run starts with the Sun suddenly overhead on 240 K ground. Crank-Nicolson is
    not L-stable: it does not damp the stiff modes that this excites, and the first
    cells overshoot the radiative-equilibrium temperature, which no physical solution
    can exceed. Abrupt forcing is the norm at the poles (sunrise over a ridge), so
    this is disqualifying there. BDF2, also second order, stays within the bound.
    """
    flux, step = day_night_flux()
    initial = np.full(GRID.n_cells, 240.0)
    hottest = (flux.max() / (HAYNE.emissivity * SIGMA)) ** 0.25
    cn = solve_conduction(GRID, HAYNE, initial, flux, step, scheme="crank_nicolson")
    bdf2 = solve_conduction(GRID, HAYNE, initial, flux, step, scheme="bdf2")
    assert cn.temperature_k.max() > hottest + 5.0
    assert bdf2.temperature_k.max() < hottest


@pytest.mark.parametrize("scheme", ["explicit", "backward_euler"])
def test_th_p1_more_absorbed_flux_never_cools_anything(scheme):
    """Comparison principle, which the monotone schemes inherit from the equation."""
    flux, step = day_night_flux(n_cycles=2)
    initial = np.full(GRID.n_cells, 240.0)
    substeps = 40 if scheme == "explicit" else 1  # identical discretization in both runs
    base = solve_conduction(GRID, HAYNE, initial, flux, step, scheme=scheme, substeps=substeps)
    extra = flux + 20.0 * (np.arange(flux.size) % 50 < 10)
    more = solve_conduction(GRID, HAYNE, initial, extra, step, scheme=scheme, substeps=substeps)
    assert np.all(more.temperature_k >= base.temperature_k - 1e-9)
    assert np.all(more.surface_temperature_k >= base.surface_temperature_k - 1e-9)
    assert more.temperature_k.max() > base.temperature_k.max()


def test_all_schemes_converge_to_the_same_solution():
    flux, step = day_night_flux(n_cycles=2, per_cycle=480)
    initial = np.full(GRID.n_cells, 240.0)
    results = {
        "explicit": solve_conduction(GRID, HAYNE, initial, flux, step, scheme="explicit"),
        "backward_euler": solve_conduction(
            GRID, HAYNE, initial, flux, step, scheme="backward_euler", substeps=64
        ),
        "crank_nicolson": solve_conduction(
            GRID, HAYNE, initial, flux, step, scheme="crank_nicolson", substeps=8
        ),
        "bdf2": solve_conduction(GRID, HAYNE, initial, flux, step, scheme="bdf2", substeps=8),
    }
    reference = results["bdf2"].surface_temperature_k
    for name, result in results.items():
        assert np.max(np.abs(result.surface_temperature_k - reference)) < 0.3, name


def test_explicit_substeps_follow_the_stability_limit():
    flux, step = day_night_flux(n_cycles=1)
    result = solve_conduction(
        GRID, HAYNE, np.full(GRID.n_cells, 240.0), flux, step, scheme="explicit"
    )
    limit = explicit_stable_step_s(GRID, HAYNE, 400.0)
    assert step / result.substeps <= limit
    with pytest.raises(SolverError, match="stability limit"):
        solve_conduction(
            GRID, HAYNE, np.full(GRID.n_cells, 240.0), flux, step, scheme="explicit", substeps=1
        )


def test_invalid_inputs_are_rejected():
    flux, step = day_night_flux(n_cycles=1)
    with pytest.raises(ValueError):
        solve_conduction(GRID, HAYNE, np.full(3, 240.0), flux, step)
    with pytest.raises(ValueError):
        solve_conduction(GRID, HAYNE, np.full(GRID.n_cells, 240.0), flux[:1], step)
    with pytest.raises(ValueError):
        solve_conduction(GRID, HAYNE, np.full(GRID.n_cells, 1.0), flux, step)
    with pytest.raises(ValueError):
        solve_conduction(GRID, HAYNE, np.full(GRID.n_cells, 240.0), flux, -1.0)
