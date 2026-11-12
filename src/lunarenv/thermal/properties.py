"""Thermophysical properties of lunar regolith after Hayne et al. (2017), Appendix A.

Equations (numbers refer to Hayne2017):

* density (A2):        rho(z) = rho_d - (rho_d - rho_s) exp(-z / H)
* conductivity (A4-5): K(T, rho) = K_c(rho) [1 + chi (T / 350)^3],
                       K_c = K_d - (K_d - K_s) (rho_d - rho) / (rho_d - rho_s)
* heat capacity (A6):  c_p(T) = c0 + c1 T + c2 T^2 + c3 T^3 + c4 T^4
* albedo (A8):         A(theta) = A0 + a (theta / (pi/4))^3 + b (theta / (pi/2))^8

Default values are Table A1 of the paper. Two derived quantities are added because the
solver needs them:

* specific enthalpy ``e(T)``, the integral of ``c_p`` from 0 to T, so that energy is
  conserved exactly while the heat capacity varies with temperature;
* the Kirchhoff potential ``phi(T)``, the integral of ``K`` over temperature, which
  gives the exact steady conductive flux across a layer of uniform density.

Range of validity (stated by the source):

* the property fit comes from Diviner data at latitudes within 70 degrees of the
  equator; use poleward of that is an extrapolation;
* the heat-capacity polynomial is described as valid from below 90 K to above 400 K
  and becomes negative below 1.3 K. Temperatures under :data:`CP_VALID_MIN_K` are
  flagged by the model; the solver itself only requires ``T > CP_POSITIVE_MIN_K``.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit
from pydantic import BaseModel, ConfigDict, Field, model_validator

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4

#: Lowest temperature at which the source describes the c_p fit as valid.
CP_VALID_MIN_K = 90.0
#: The solver refuses temperatures below this (c_p is positive above 1.3 K).
CP_POSITIVE_MIN_K = 5.0
_REFERENCE_T_K = 350.0


class RegolithProperties(BaseModel):
    """Parameters of the regolith model (defaults: Hayne2017, Table A1)."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    surface_conductivity_w_m_k: float = Field(7.4e-4, gt=0.0, description="K_s")
    deep_conductivity_w_m_k: float = Field(3.4e-3, gt=0.0, description="K_d")
    surface_density_kg_m3: float = Field(1100.0, gt=0.0, description="rho_s")
    deep_density_kg_m3: float = Field(1800.0, gt=0.0, description="rho_d")
    scale_height_m: float = Field(0.06, gt=0.0, description="H (0.068 m global mode)")
    radiative_parameter: float = Field(2.7, ge=0.0, description="chi")
    heat_capacity_coefficients: tuple[float, float, float, float, float] = Field(
        (-3.6125, 2.7431, 2.3616e-3, -1.2340e-5, 8.9093e-9),
        description="c0..c4 of c_p(T) in J kg-1 K-1",
    )
    emissivity: float = Field(0.95, gt=0.0, le=1.0)
    albedo_normal: float = Field(0.12, ge=0.0, lt=1.0, description="A0")
    albedo_a: float = Field(0.06, ge=0.0)
    albedo_b: float = Field(0.25, ge=0.0)
    heat_flow_w_m2: float = Field(0.018, ge=0.0, description="Q, interior heat flow")

    @model_validator(mode="after")
    def _check(self) -> RegolithProperties:
        if self.deep_density_kg_m3 < self.surface_density_kg_m3:
            raise ValueError("deep density must not be below surface density")
        return self

    def as_array(self) -> np.ndarray:
        """Parameters packed for the compiled kernels (see ``P_*`` indices)."""
        return np.array(
            [
                self.surface_conductivity_w_m_k,
                self.deep_conductivity_w_m_k,
                self.surface_density_kg_m3,
                self.deep_density_kg_m3,
                self.scale_height_m,
                self.radiative_parameter,
                *self.heat_capacity_coefficients,
                self.emissivity,
                self.heat_flow_w_m2,
            ],
            dtype=np.float64,
        )


# indices into RegolithProperties.as_array()
P_KS, P_KD, P_RHOS, P_RHOD, P_H, P_CHI, P_C0, P_C1, P_C2, P_C3, P_C4, P_EMIS, P_Q = range(13)


@njit(cache=True)
def density(z, p):  # type: ignore[no-untyped-def]
    """Bulk density (kg m-3) at depth z (m)."""
    return p[P_RHOD] - (p[P_RHOD] - p[P_RHOS]) * math.exp(-z / p[P_H])


@njit(cache=True)
def contact_conductivity(rho, p):  # type: ignore[no-untyped-def]
    """Phonon (contact) conductivity K_c (W m-1 K-1), linear in density."""
    span = p[P_RHOD] - p[P_RHOS]
    if span <= 0.0:
        return p[P_KD]
    return p[P_KD] - (p[P_KD] - p[P_KS]) * (p[P_RHOD] - rho) / span


@njit(cache=True)
def conductivity(temperature, k_contact, chi):  # type: ignore[no-untyped-def]
    """Thermal conductivity (W m-1 K-1) including the radiative T^3 term."""
    ratio = temperature / _REFERENCE_T_K
    return k_contact * (1.0 + chi * ratio * ratio * ratio)


@njit(cache=True)
def conductivity_dt(temperature, k_contact, chi):  # type: ignore[no-untyped-def]
    """Derivative of the conductivity with respect to temperature."""
    return k_contact * chi * 3.0 * temperature * temperature / _REFERENCE_T_K**3


@njit(cache=True)
def kirchhoff(temperature, k_contact, chi):  # type: ignore[no-untyped-def]
    """Kirchhoff potential: integral of K from 0 to T (W m-1)."""
    return k_contact * (temperature + chi * temperature**4 / (4.0 * _REFERENCE_T_K**3))


@njit(cache=True)
def heat_capacity(temperature, p):  # type: ignore[no-untyped-def]
    """Specific heat capacity c_p (J kg-1 K-1)."""
    t = temperature
    return p[P_C0] + t * (p[P_C1] + t * (p[P_C2] + t * (p[P_C3] + t * p[P_C4])))


@njit(cache=True)
def enthalpy(temperature, p):  # type: ignore[no-untyped-def]
    """Specific enthalpy (J kg-1): integral of c_p from 0 to T."""
    t = temperature
    return t * (
        p[P_C0]
        + t * (p[P_C1] / 2.0 + t * (p[P_C2] / 3.0 + t * (p[P_C3] / 4.0 + t * p[P_C4] / 5.0)))
    )


@njit(cache=True)
def temperature_from_enthalpy(e, guess, p):  # type: ignore[no-untyped-def]
    """Invert :func:`enthalpy` by Newton iteration (monotonic above 1.3 K)."""
    t = guess
    for _ in range(50):
        step = (enthalpy(t, p) - e) / heat_capacity(t, p)
        t -= step
        if t < CP_POSITIVE_MIN_K:
            t = CP_POSITIVE_MIN_K
        if abs(step) < 1e-10 * max(t, 1.0):
            break
    return t


def albedo(incidence_rad: np.ndarray | float, properties: RegolithProperties) -> np.ndarray:
    """Bolometric Bond albedo at solar incidence angle theta (clipped to [0, pi/2]).

    The law is an empirical fit to equatorial brightness temperatures (Keihm1984,
    Vasavada2012, constants from Hayne2017). At grazing incidence it rises steeply
    (A0 + 8a + b at 90 degrees, 0.85 with the defaults), which matters at the poles,
    where sunlight is always grazing on level ground. The result is capped at 1.
    """
    theta = np.clip(np.asarray(incidence_rad, dtype=np.float64), 0.0, np.pi / 2.0)
    value = (
        properties.albedo_normal
        + properties.albedo_a * (theta / (np.pi / 4.0)) ** 3
        + properties.albedo_b * (theta / (np.pi / 2.0)) ** 8
    )
    return np.asarray(np.minimum(value, 1.0))


def radiative_equilibrium_temperature(
    absorbed_flux_w_m2: np.ndarray | float, emissivity: float
) -> np.ndarray:
    """Temperature (K) at which a surface emits what it absorbs."""
    flux = np.maximum(np.asarray(absorbed_flux_w_m2, dtype=np.float64), 0.0)
    return np.asarray((flux / (emissivity * STEFAN_BOLTZMANN_W_M2_K4)) ** 0.25)
