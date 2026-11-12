"""One-dimensional heat conduction with a radiating surface: discretization and solvers.

Equation:  rho(z) c_p(T) dT/dt = d/dz ( K(T, rho) dT/dz ),  z positive downwards.

Boundary conditions
    surface (z = 0):  Q_abs + K dT/dz = eps sigma T_s^4      (Hayne2017, A7)
    base (z = z*):    K dT/dz = Q                             (A12)

Spatial discretization (PLAN.md decision D8): conservative finite volumes.

* Unknowns are cell temperatures plus a massless surface temperature ``T_s``.
* The state variable advanced in time is the specific enthalpy ``e(T)``, the integral of
  ``c_p``; with a temperature-dependent heat capacity this is what makes the stored
  energy change by exactly the net boundary flux.
* The conductivity is separable, ``K = K_c(z) k(T)`` with ``k(T) = 1 + chi (T/350)^3``,
  so the heat flux is ``-K_c(z) d(psi)/dz`` with the Kirchhoff variable
  ``psi(T) = T + chi T^4 / (4 * 350^3)``, which is continuous wherever temperature is.
  Fluxes are therefore written as differences of ``psi`` across series resistances
  ``dz / (2 K_c)``. This is exact for steady conduction whatever the temperature
  contrast, including the factor-of-several change of conductivity across the top
  millimetres at sunrise, and no rule for averaging conductivities is needed.
* A consequence used by the spin-up: in a periodic state the time-mean of ``psi``
  at depth equals its time-mean at the surface plus the interior heat flow times the
  resistance in between.

Time integration (decision D9); the default is chosen by measurement in
``docs/convergence/thermal.md``:

* ``EXPLICIT``: forward Euler. Conditionally stable; sub-steps are set from the
  diffusion limit. This is the class of scheme used by Hayne2017.
* ``BACKWARD_EULER``: first order, L-stable.
* ``CRANK_NICOLSON``: second order, A-stable but not L-stable.
* ``BDF2``: second order, L-stable; first step by backward Euler.

Implicit steps solve the full nonlinear system (conductivity, heat capacity and the
T^4 surface term) by Newton iteration with the exact tridiagonal Jacobian.

Energy accounting. For the one-step schemes the change of stored enthalpy equals the
time-integrated net boundary flux to solver tolerance; the solver returns both sides.
For BDF2 the corresponding two-step identity is accumulated instead.

Assumptions: one-dimensional (no lateral conduction), properties depend on depth and
temperature only, no rocks, no phase change or volatile transport.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numba import njit

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4
from lunarenv.geometry.projection import FloatArray
from lunarenv.thermal.grid import DepthGrid
from lunarenv.thermal.properties import (
    CP_POSITIVE_MIN_K,
    P_CHI,
    P_EMIS,
    P_Q,
    RegolithProperties,
    conductivity,
    contact_conductivity,
    density,
    enthalpy,
    heat_capacity,
    kirchhoff,
    temperature_from_enthalpy,
)

EXPLICIT = 0
BACKWARD_EULER = 1
CRANK_NICOLSON = 2
BDF2 = 3
SCHEMES = {
    "explicit": EXPLICIT,
    "backward_euler": BACKWARD_EULER,
    "crank_nicolson": CRANK_NICOLSON,
    "bdf2": BDF2,
}
SchemeName = Literal["explicit", "backward_euler", "crank_nicolson", "bdf2"]

BC_FLUX = 0  # radiating surface driven by an absorbed flux
BC_TEMPERATURE = 1  # prescribed surface temperature
BoundaryName = Literal["flux", "temperature"]

_SIGMA = STEFAN_BOLTZMANN_W_M2_K4
_NEWTON_TOLERANCE_K = 1e-9
_MAX_NEWTON = 40

STATUS_OK = 0
STATUS_NOT_CONVERGED = 1
STATUS_NOT_FINITE = 2


@njit(cache=True)
def _surface_balance(q_abs, t0, kc0, chi, half_dz, emissivity, guess):  # type: ignore[no-untyped-def]
    """Surface temperature from Q_abs = eps sigma T^4 + G(T, T_0).

    The residual decreases monotonically in T, so Newton is safeguarded by bisection on
    a bracket that always contains the root.
    """
    phi0 = kirchhoff(t0, kc0, chi)
    low = CP_POSITIVE_MIN_K
    high = max(2.0 * guess, 2.0 * t0, 100.0)
    while q_abs - emissivity * _SIGMA * high**4 - (kirchhoff(high, kc0, chi) - phi0) / half_dz > 0:
        high *= 2.0
    t = min(max(guess, low), high)
    for _ in range(100):
        f = q_abs - emissivity * _SIGMA * t**4 - (kirchhoff(t, kc0, chi) - phi0) / half_dz
        if f > 0.0:
            low = t
        else:
            high = t
        slope = -4.0 * emissivity * _SIGMA * t**3 - conductivity(t, kc0, chi) / half_dz
        candidate = t - f / slope
        if not (low < candidate < high):
            candidate = 0.5 * (low + high)
        if abs(candidate - t) < _NEWTON_TOLERANCE_K:
            return candidate
        t = candidate
    return t


@njit(cache=True)
def _psi(temperature, chi):  # type: ignore[no-untyped-def]
    """Kirchhoff variable: integral of k(T) = 1 + chi (T/350)^3."""
    return kirchhoff(temperature, 1.0, chi)


@njit(cache=True)
def _conductances(dz, kc):  # type: ignore[no-untyped-def]
    """Face conductances in the Kirchhoff variable (W m-2 per unit psi).

    ``g[0]`` links the surface to the first cell centre; ``g[i + 1]`` links cells ``i``
    and ``i + 1``; the base face carries the prescribed heat flow and has no entry.
    """
    n = dz.size
    g = np.empty(n)
    g[0] = kc[0] / (0.5 * dz[0])
    for i in range(n - 1):
        g[i + 1] = 1.0 / (0.5 * dz[i] / kc[i] + 0.5 * dz[i + 1] / kc[i + 1])
    return g


@njit(cache=True)
def _fluxes(t_surface, temp, g, chi, q_base, flux):  # type: ignore[no-untyped-def]
    """Downward heat flux at every cell face; ``flux[i]`` is the top face of cell ``i``."""
    n = temp.size
    flux[0] = g[0] * (_psi(t_surface, chi) - _psi(temp[0], chi))
    for i in range(n - 1):
        flux[i + 1] = g[i + 1] * (_psi(temp[i], chi) - _psi(temp[i + 1], chi))
    flux[n] = -q_base


@njit(cache=True)
def _implicit_step(  # type: ignore[no-untyped-def]
    t_surface, temp, mass, g, p, forcing, boundary, dt, theta, a0, history, old_divergence
):
    """Advance one implicit step in place. Returns (surface temperature, status).

    Solves, for every cell,
        mass (a0 e(T) + history) / dt + theta (F_bottom - F_top) + (1 - theta) old_div = 0
    together with the surface equation, by Newton iteration.
    """
    n = temp.size
    chi = p[P_CHI]
    emissivity = p[P_EMIS]
    q_base = p[P_Q]
    lower = np.empty(n + 1)
    diagonal = np.empty(n + 1)
    upper = np.empty(n + 1)
    residual = np.empty(n + 1)
    delta = np.empty(n + 1)
    ts = t_surface
    for _ in range(_MAX_NEWTON):
        # ---- surface row (index 0)
        k_above = conductivity(ts, 1.0, chi)
        psi_above = _psi(ts, chi)
        if boundary == BC_FLUX:
            flux_top = g[0] * (psi_above - _psi(temp[0], chi))
            residual[0] = forcing - emissivity * _SIGMA * ts**4 - flux_top
            diagonal[0] = -4.0 * emissivity * _SIGMA * ts**3 - g[0] * k_above
            upper[0] = g[0] * conductivity(temp[0], 1.0, chi)
        else:
            residual[0] = ts - forcing
            diagonal[0] = 1.0
            upper[0] = 0.0
        lower[0] = 0.0
        # ---- cells (index i + 1); the face above cell i has conductance g[i]
        for i in range(n):
            k_self = conductivity(temp[i], 1.0, chi)
            psi_self = _psi(temp[i], chi)
            flux_top = g[i] * (psi_above - psi_self)
            if i < n - 1:
                flux_bottom = g[i + 1] * (psi_self - _psi(temp[i + 1], chi))
                d_bottom_self = g[i + 1] * k_self
                d_bottom_below = -g[i + 1] * conductivity(temp[i + 1], 1.0, chi)
            else:
                flux_bottom = -q_base
                d_bottom_self = 0.0
                d_bottom_below = 0.0
            row = i + 1
            residual[row] = (
                mass[i] * (a0 * enthalpy(temp[i], p) + history[i]) / dt
                + theta * (flux_bottom - flux_top)
                + (1.0 - theta) * old_divergence[i]
            )
            lower[row] = -theta * g[i] * k_above
            diagonal[row] = mass[i] * a0 * heat_capacity(temp[i], p) / dt + theta * (
                d_bottom_self + g[i] * k_self
            )
            upper[row] = theta * d_bottom_below
            k_above = k_self
            psi_above = psi_self
        # ---- tridiagonal solve of J delta = -residual (Thomas algorithm)
        for row in range(1, n + 1):
            factor = lower[row] / diagonal[row - 1]
            diagonal[row] -= factor * upper[row - 1]
            residual[row] -= factor * residual[row - 1]
        delta[n] = -residual[n] / diagonal[n]
        for row in range(n - 1, -1, -1):
            delta[row] = (-residual[row] - upper[row] * delta[row + 1]) / diagonal[row]
        ts = max(ts + delta[0], CP_POSITIVE_MIN_K)
        largest = abs(delta[0])
        for i in range(n):
            temp[i] = max(temp[i] + delta[i + 1], CP_POSITIVE_MIN_K)
            largest = max(largest, abs(delta[i + 1]))
        if not math.isfinite(largest):
            return ts, STATUS_NOT_FINITE
        if largest < _NEWTON_TOLERANCE_K:
            return ts, STATUS_OK
    return ts, STATUS_NOT_CONVERGED


@njit(cache=True)
def _cell_arrays(edges, p):  # type: ignore[no-untyped-def]
    """Cell thickness, areal mass and contact conductivity."""
    n = edges.size - 1
    dz = edges[1:] - edges[:-1]
    mass = np.empty(n)
    kc = np.empty(n)
    for i in range(n):
        rho = density(0.5 * (edges[i] + edges[i + 1]), p)
        mass[i] = rho * dz[i]
        kc[i] = contact_conductivity(rho, p)
    return dz, mass, kc


@njit(cache=True)
def _integrate(  # type: ignore[no-untyped-def]
    edges,
    p,
    initial,
    forcing,
    dt_forcing,
    substeps,
    scheme,
    boundary,
    out_surface,
    out_temp,
    balance,
):
    """Integrate over the forcing series; outputs are stored at every forcing time."""
    n = initial.size
    chi = p[P_CHI]
    emissivity = p[P_EMIS]
    q_base = p[P_Q]
    dz, mass, kc = _cell_arrays(edges, p)
    g = _conductances(dz, kc)
    half = 0.5 * dz[0]

    temp = initial.copy()
    e_now = np.empty(n)
    for i in range(n):
        e_now[i] = enthalpy(temp[i], p)
    e_previous = e_now.copy()
    if boundary == BC_FLUX:
        ts = _surface_balance(forcing[0], temp[0], kc[0], chi, half, emissivity, temp[0])
    else:
        ts = forcing[0]
    flux = np.empty(n + 1)
    history = np.empty(n)
    old_divergence = np.zeros(n)
    out_surface[0] = ts
    out_temp[0, :] = temp
    dt = dt_forcing / substeps
    stored_change = 0.0  # left-hand side of the scheme's energy identity, J m-2
    boundary_input = 0.0  # right-hand side: time-integrated net boundary flux, J m-2
    step_count = 0
    status = STATUS_OK

    for k in range(forcing.size - 1):
        for s in range(substeps):
            weight = (s + 1.0) / substeps
            forcing_new = forcing[k] * (1.0 - weight) + forcing[k + 1] * weight
            if scheme == EXPLICIT:
                _fluxes(ts, temp, g, chi, q_base, flux)
                for i in range(n):
                    change = -dt * (flux[i + 1] - flux[i]) / mass[i]
                    e_now[i] += change
                    stored_change += mass[i] * change
                    temp[i] = temperature_from_enthalpy(e_now[i], temp[i], p)
                boundary_input += dt * (flux[0] - flux[n])
                if boundary == BC_FLUX:
                    ts = _surface_balance(forcing_new, temp[0], kc[0], chi, half, emissivity, ts)
                else:
                    ts = forcing_new
            else:
                use_bdf2 = scheme == BDF2 and step_count > 0
                if use_bdf2:
                    theta, a0 = 1.0, 1.5
                    for i in range(n):
                        history[i] = -2.0 * e_now[i] + 0.5 * e_previous[i]
                else:
                    theta = 0.5 if scheme == CRANK_NICOLSON else 1.0
                    a0 = 1.0
                    for i in range(n):
                        history[i] = -e_now[i]
                if theta < 1.0:
                    _fluxes(ts, temp, g, chi, q_base, flux)
                    for i in range(n):
                        old_divergence[i] = flux[i + 1] - flux[i]
                    old_net = flux[0] - flux[n]
                else:
                    old_net = 0.0
                ts, status = _implicit_step(
                    ts,
                    temp,
                    mass,
                    g,
                    p,
                    forcing_new,
                    boundary,
                    dt,
                    theta,
                    a0,
                    history,
                    old_divergence,
                )
                if status != STATUS_OK:
                    return status
                _fluxes(ts, temp, g, chi, q_base, flux)
                boundary_input += dt * (theta * (flux[0] - flux[n]) + (1.0 - theta) * old_net)
                for i in range(n):
                    e_new = enthalpy(temp[i], p)
                    if use_bdf2:
                        stored_change += mass[i] * (
                            1.5 * e_new - 2.0 * e_now[i] + 0.5 * e_previous[i]
                        )
                    else:
                        stored_change += mass[i] * (e_new - e_now[i])
                    e_previous[i] = e_now[i]
                    e_now[i] = e_new
            step_count += 1
        if not math.isfinite(ts):
            return STATUS_NOT_FINITE
        out_surface[k + 1] = ts
        out_temp[k + 1, :] = temp
    balance[0] = stored_change
    balance[1] = boundary_input
    return status


@njit(cache=True)
def _explicit_limit(edges, p, t_hot):  # type: ignore[no-untyped-def]
    """Largest stable forward-Euler step (s) if every cell were at ``t_hot``.

    Diffusion number 1/2: dt <= mass c_p / (2 k(T) (g_above + g_below)). Hot is the
    restrictive case because the conductivity grows faster than the heat capacity.
    """
    dz, mass, kc = _cell_arrays(edges, p)
    g = _conductances(dz, kc)
    n = dz.size
    k_hot = conductivity(t_hot, 1.0, p[P_CHI])
    limit = math.inf
    for i in range(n):
        below = g[i + 1] if i < n - 1 else 0.0
        limit = min(limit, mass[i] * heat_capacity(t_hot, p) / (2.0 * k_hot * (g[i] + below)))
    return limit


@dataclass(frozen=True, slots=True)
class ConductionResult:
    """Temperatures at every forcing time.

    ``temperature_k[k, i]`` is cell ``i`` at forcing time ``k``; ``surface_temperature_k``
    is the massless surface node. ``stored_energy_change_j_m2`` and
    ``boundary_energy_input_j_m2`` are the two sides of the scheme's energy identity;
    their difference measures how well the solver conserved energy.
    """

    surface_temperature_k: FloatArray
    temperature_k: FloatArray
    stored_energy_change_j_m2: float
    boundary_energy_input_j_m2: float
    substeps: int
    scheme: str

    @property
    def energy_residual_j_m2(self) -> float:
        return self.stored_energy_change_j_m2 - self.boundary_energy_input_j_m2


class SolverError(RuntimeError):
    """The conduction solver failed (non-convergence or non-finite temperatures)."""


def explicit_stable_step_s(
    grid: DepthGrid, properties: RegolithProperties, hottest_k: float
) -> float:
    """Largest stable forward-Euler step for temperatures up to ``hottest_k``."""
    return float(_explicit_limit(grid.edges_m, properties.as_array(), hottest_k))


def solve_conduction(
    grid: DepthGrid,
    properties: RegolithProperties,
    initial_temperature_k: FloatArray,
    forcing: FloatArray,
    forcing_step_s: float,
    *,
    scheme: SchemeName = "bdf2",
    substeps: int | None = None,
    boundary: BoundaryName = "flux",
) -> ConductionResult:
    """Integrate the conduction equation over a uniformly sampled forcing series.

    ``forcing`` is the absorbed surface flux in W m-2 (``boundary="flux"``) or the
    surface temperature in K (``boundary="temperature"``), sampled every
    ``forcing_step_s`` and interpolated linearly in between. The solver takes
    ``substeps`` steps per forcing interval; for the explicit scheme the default is the
    smallest number that satisfies the stability limit with a safety factor of 0.8.
    """
    forcing = np.ascontiguousarray(forcing, dtype=np.float64)
    initial = np.ascontiguousarray(initial_temperature_k, dtype=np.float64)
    if forcing.ndim != 1 or forcing.size < 2:
        raise ValueError("forcing must be a series of at least two samples")
    if initial.shape != (grid.n_cells,):
        raise ValueError("initial temperature must have one value per cell")
    if forcing_step_s <= 0.0:
        raise ValueError("forcing step must be positive")
    if np.any(initial < CP_POSITIVE_MIN_K) or not np.all(np.isfinite(forcing)):
        raise ValueError("initial temperatures or forcing are out of range")
    code = SCHEMES[scheme]
    bc = BC_FLUX if boundary == "flux" else BC_TEMPERATURE
    p = properties.as_array()
    if substeps is None:
        if code == EXPLICIT:
            if bc == BC_FLUX:
                hottest = (max(float(forcing.max()), 0.0) / (p[P_EMIS] * _SIGMA)) ** 0.25
            else:
                hottest = float(forcing.max())
            hottest = max(hottest, float(initial.max()))
            limit = _explicit_limit(grid.edges_m, p, hottest)
            substeps = max(1, math.ceil(forcing_step_s / (0.8 * limit)))
        else:
            substeps = 1
    if substeps < 1:
        raise ValueError("substeps must be at least 1")
    if code == EXPLICIT:
        # An unstable explicit run does not necessarily produce non-finite numbers, so
        # the stability condition is checked up front instead of trusting the output.
        if bc == BC_FLUX:
            hottest = (max(float(forcing.max()), 0.0) / (p[P_EMIS] * _SIGMA)) ** 0.25
        else:
            hottest = float(forcing.max())
        limit = _explicit_limit(grid.edges_m, p, max(hottest, float(initial.max())))
        if forcing_step_s / substeps > limit:
            raise SolverError(
                f"explicit step {forcing_step_s / substeps:.3g} s exceeds the stability limit "
                f"{limit:.3g} s; use at least {math.ceil(forcing_step_s / limit)} substeps"
            )
    out_surface = np.empty(forcing.size)
    out_temp = np.empty((forcing.size, grid.n_cells))
    balance = np.zeros(2)
    status = _integrate(
        grid.edges_m,
        p,
        initial,
        forcing,
        float(forcing_step_s),
        int(substeps),
        code,
        bc,
        out_surface,
        out_temp,
        balance,
    )
    if status == STATUS_NOT_CONVERGED:
        raise SolverError("Newton iteration did not converge; reduce the time step")
    if status == STATUS_NOT_FINITE or not np.all(np.isfinite(out_temp)):
        raise SolverError("non-finite temperatures; the time step is unstable")
    return ConductionResult(
        out_surface, out_temp, float(balance[0]), float(balance[1]), int(substeps), scheme
    )
