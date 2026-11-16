"""Thermal history of the regolith at a point: spin-up, integration and summary.

Initial condition. Temperatures at depth remember decades of surface history (the
conduction time across 2 m of regolith is of order 25 years), so the state at the start
of the period cannot be guessed. It is obtained by spin-up:

1. the forcing of the year that precedes the period is applied cyclically;
2. after each cycle the mean temperature profile is set to the one the periodic state
   must have. Because the conductivity is separable (see ``solver``), the time-mean of
   the Kirchhoff variable ``psi`` at any depth equals its time-mean at the surface plus
   the interior heat flow times the thermal resistance in between. The surface
   temperature depends only weakly on what lies below, so every cell is shifted, in
   ``psi``, by the difference between that target and its own mean over the cycle;
3. cycles repeat until the end-of-cycle state stops changing.

Step 2 replaces decades of simulated relaxation of the deep layers by an exact
property of the solution; it does not alter the converged state, which satisfies the
identity by itself.

The year before the period is not exactly periodic (a year is 12.37 lunations and the
illumination geometry drifts over 18.6 years). The upper layers forget this within a
few lunations, and the deep profile depends only on yearly means; the residual
sensitivity to the initial state is measured in ``docs/convergence/thermal.md``.

Validity flags raised with the result:

* temperatures below 90 K, where the heat-capacity fit is at or beyond its stated range;
* sites poleward of 70 degrees, where the Diviner-derived properties are extrapolated;
* epochs in shadow, where neglected terrain radiation is the only energy source.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from lunarenv.geometry.projection import FloatArray
from lunarenv.thermal.grid import DepthGrid, geometric_grid
from lunarenv.thermal.properties import (
    CP_POSITIVE_MIN_K,
    CP_VALID_MIN_K,
    WOODS_ROBINSON_MIN_K,
    RegolithProperties,
    conductivity,
    contact_conductivity,
    density,
    kirchhoff,
    radiative_equilibrium_temperature,
)
from lunarenv.thermal.solver import SchemeName, solve_conduction

SYNODIC_MONTH_S = 2.55024e6  # Hayne2017 Table A1 (29.5306 days)

#: Numerical error of the surface temperature at the default settings with hourly
#: forcing, from docs/convergence/thermal.md (Shackleton rim test point, 2026).
#: Largest instantaneous differences against the finest setting tried; errors of
#: extremes and means are 10 to 100 times smaller.
NUMERICAL_ERROR_ESTIMATE = {
    "time_integration_k": 0.02,
    "forcing_sampling_60min_k": 0.12,
    "grid_k": 0.05,
    "spinup_k": 0.0002,
    "surface_temperature_total_k": 0.2,
}
JULIAN_YEAR_S = 365.25 * 86400.0


class ThermalSettings(BaseModel):
    """Numerical parameters of the thermal model (defaults: docs/convergence/thermal.md)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    first_cell_m: float = Field(0.002, gt=0.0, description="Thickness of the top cell.")
    growth: float = Field(1.12, ge=1.0, description="Thickness ratio of consecutive cells.")
    depth_m: float = Field(2.0, gt=0.0, description="Depth of the model domain.")
    scheme: SchemeName = Field("bdf2", description="Time integrator.")
    substeps: int = Field(2, ge=1, description="Solver steps per forcing sample.")
    spinup_tolerance_k: float = Field(1e-3, gt=0.0)
    spinup_max_cycles: int = Field(60, ge=1)

    def grid(self) -> DepthGrid:
        return geometric_grid(self.first_cell_m, self.growth, self.depth_m)


@dataclass(frozen=True, slots=True)
class SpinUp:
    """State reached by the spin-up and how well it converged."""

    temperature_k: FloatArray
    cycles: int
    residual_k: float
    converged: bool


@dataclass(frozen=True, slots=True)
class ThermalResult:
    """Temperatures over the period at every forcing time and cell."""

    time_s: FloatArray
    depth_m: FloatArray
    surface_temperature_k: FloatArray
    temperature_k: FloatArray
    absorbed_flux_w_m2: FloatArray
    spinup: SpinUp
    energy_residual_j_m2: float
    warnings: tuple[str, ...] = field(default_factory=tuple)

    def at_depth(self, depth_m: float) -> FloatArray:
        """Temperature series at a depth, linear between the surface and cell centres."""
        if not 0.0 <= depth_m <= self.depth_m[-1]:
            raise ValueError("depth outside the model domain")
        nodes = np.concatenate(([0.0], self.depth_m))
        values = np.column_stack([self.surface_temperature_k, self.temperature_k])
        index = int(np.clip(np.searchsorted(nodes, depth_m) - 1, 0, nodes.size - 2))
        weight = (depth_m - nodes[index]) / (nodes[index + 1] - nodes[index])
        return np.asarray(values[:, index] * (1.0 - weight) + values[:, index + 1] * weight)


def _psi(temperature: FloatArray | float, chi: float) -> FloatArray:
    return np.asarray(kirchhoff(np.asarray(temperature, dtype=np.float64), 1.0, chi))


def _psi_inverse(psi: FloatArray, guess: FloatArray, chi: float) -> FloatArray:
    """Temperature whose Kirchhoff variable is ``psi`` (Newton; psi is monotonic)."""
    t = np.array(guess, dtype=np.float64)
    for _ in range(60):
        step = (_psi(t, chi) - psi) / conductivity(t, 1.0, chi)
        t = np.maximum(t - step, CP_POSITIVE_MIN_K)
        if float(np.abs(step).max()) < 1e-11:
            break
    return t


def resistance_to_cells(grid: DepthGrid, properties: RegolithProperties) -> FloatArray:
    """Thermal resistance in the Kirchhoff variable from the surface to each cell centre."""
    p = properties.as_array()
    dz = grid.thickness_m
    kc = np.array([contact_conductivity(density(z, p), p) for z in grid.centers_m])
    half = 0.5 * dz / kc
    return np.asarray(np.cumsum(2.0 * half) - half)


def spin_up(
    absorbed_flux_w_m2: FloatArray,
    step_s: float,
    properties: RegolithProperties,
    settings: ThermalSettings,
    *,
    initial_temperature_k: float | None = None,
) -> SpinUp:
    """Periodic steady state under one cycle of forcing (normally the preceding year)."""
    grid = settings.grid()
    chi = properties.radiative_parameter
    flux = np.asarray(absorbed_flux_w_m2, dtype=np.float64)
    if initial_temperature_k is None:
        equilibrium = radiative_equilibrium_temperature(flux, properties.emissivity)
        initial_temperature_k = float(max(equilibrium.mean(), 40.0))
    state = np.full(grid.n_cells, float(initial_temperature_k))
    offset = properties.heat_flow_w_m2 * resistance_to_cells(grid, properties)
    residual = np.inf
    for cycle in range(1, settings.spinup_max_cycles + 1):
        result = solve_conduction(
            grid,
            properties,
            state,
            flux,
            step_s,
            scheme=settings.scheme,
            substeps=settings.substeps,
        )
        # one full cycle: the last sample repeats the phase of the first, so it is left out
        surface_mean = float(_psi(result.surface_temperature_k[1:], chi).mean())
        cell_mean = _psi(result.temperature_k[1:], chi).mean(axis=0)
        end = result.temperature_k[-1]
        target = _psi(end, chi) + (surface_mean + offset - cell_mean)
        new_state = _psi_inverse(target, end, chi)
        residual = float(np.abs(new_state - state).max())
        state = new_state
        if residual < settings.spinup_tolerance_k:
            return SpinUp(state, cycle, residual, True)
    return SpinUp(state, settings.spinup_max_cycles, residual, False)


def run_thermal(
    absorbed_flux_w_m2: FloatArray,
    step_s: float,
    spinup_samples: int,
    properties: RegolithProperties | None = None,
    settings: ThermalSettings | None = None,
    *,
    latitude_rad: float | None = None,
    shadow: np.ndarray | None = None,
) -> ThermalResult:
    """Thermal history over a period preceded by its spin-up cycle.

    ``absorbed_flux_w_m2`` is sampled every ``step_s``; its first ``spinup_samples``
    intervals are the spin-up cycle and the rest is the period, so the period starts at
    index ``spinup_samples``. ``shadow`` optionally marks period samples with no direct
    sunlight (for the validity flag).
    """
    properties = properties or RegolithProperties()
    settings = settings or ThermalSettings()
    flux = np.asarray(absorbed_flux_w_m2, dtype=np.float64)
    if not 2 <= spinup_samples < flux.size - 1:
        raise ValueError("need a spin-up cycle followed by at least two period samples")
    grid = settings.grid()
    start = spin_up(flux[: spinup_samples + 1], step_s, properties, settings)
    period_flux = flux[spinup_samples:]
    result = solve_conduction(
        grid,
        properties,
        start.temperature_k,
        period_flux,
        step_s,
        scheme=settings.scheme,
        substeps=settings.substeps,
    )
    warnings: list[str] = []
    if not start.converged:
        warnings.append(
            f"spin-up did not reach {settings.spinup_tolerance_k} K in {start.cycles} cycles "
            f"(residual {start.residual_k:.3g} K)"
        )
    coldest = float(min(result.surface_temperature_k.min(), result.temperature_k.min()))
    cp_valid_min = (
        CP_VALID_MIN_K if properties.heat_capacity_model == "hayne2017" else WOODS_ROBINSON_MIN_K
    )
    if coldest < cp_valid_min:
        warnings.append(
            f"temperatures reach {coldest:.1f} K, below the {cp_valid_min:.0f} K down to which "
            "the heat-capacity fit is stated to be valid"
        )
    if latitude_rad is not None and abs(latitude_rad) > np.radians(70.0):
        warnings.append(
            "regolith properties were derived from data within 70 degrees of the equator; "
            "their use at this latitude is an extrapolation"
        )
    if shadow is not None and bool(np.any(shadow)):
        warnings.append(
            f"{100.0 * float(np.mean(shadow)):.1f} % of the period has no direct sunlight; "
            "terrain-scattered and terrain-emitted radiation are neglected, so temperatures "
            "in shadow are lower bounds"
        )
    return ThermalResult(
        time_s=np.arange(period_flux.size, dtype=np.float64) * step_s,
        depth_m=grid.centers_m,
        surface_temperature_k=result.surface_temperature_k,
        temperature_k=result.temperature_k,
        absorbed_flux_w_m2=period_flux,
        spinup=start,
        energy_residual_j_m2=result.energy_residual_j_m2,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class ThermalSummary:
    """Extremes and cycle amplitudes over the period (temperatures in K)."""

    surface_min_k: float
    surface_max_k: float
    surface_mean_k: float
    lunation_amplitude_mean_k: float
    lunation_amplitude_max_k: float
    max_heating_rate_k_per_h: float
    max_cooling_rate_k_per_h: float
    depth_m: FloatArray
    depth_min_k: FloatArray
    depth_max_k: FloatArray
    depth_mean_k: FloatArray


def summarize(
    result: ThermalResult, depths_m: tuple[float, ...] = (0.05, 0.1, 0.3, 1.0)
) -> ThermalSummary:
    """Extremes, mean, cycle amplitudes and fastest temperature changes.

    The amplitude of a thermal cycle is the surface maximum minus minimum within each
    synodic month of the period (whole months only).
    """
    surface = result.surface_temperature_k
    step = float(result.time_s[1] - result.time_s[0])
    per_lunation = max(2, round(SYNODIC_MONTH_S / step))
    count = surface.size // per_lunation
    if count:
        windows = surface[: count * per_lunation].reshape(count, per_lunation)
        amplitude = windows.max(axis=1) - windows.min(axis=1)
    else:
        amplitude = np.array([surface.max() - surface.min()])
    rate = np.diff(surface) / step * 3600.0
    depths = np.array([d for d in depths_m if d <= result.depth_m[-1]])
    series = [result.at_depth(float(d)) for d in depths]
    return ThermalSummary(
        surface_min_k=float(surface.min()),
        surface_max_k=float(surface.max()),
        surface_mean_k=float(surface.mean()),
        lunation_amplitude_mean_k=float(amplitude.mean()),
        lunation_amplitude_max_k=float(amplitude.max()),
        max_heating_rate_k_per_h=float(rate.max()),
        max_cooling_rate_k_per_h=float(-rate.min()),
        depth_m=depths,
        depth_min_k=np.array([s.min() for s in series]),
        depth_max_k=np.array([s.max() for s in series]),
        depth_mean_k=np.array([s.mean() for s in series]),
    )
