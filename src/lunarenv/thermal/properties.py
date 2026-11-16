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

Alternative heat capacity (``heat_capacity_model="woods_robinson2019"``): the
Woods-Robinson et al. (2019) estimate as tabulated by DSNE (SLS-SPEC-159 Rev I, Table
3.4.6.2.2-1, 10-300 K, stated valid 10-400 K), interpolated linearly in temperature,
linear to zero below 10 K, and the Hayne polynomial above 300 K (0.3 % below the table
there). It differs from the polynomial below about 100 K, where the polynomial falls
towards zero; above 100 K the two agree within 2 %.

Low-temperature conductivity (``conductivity_model="woods_robinson2019"``, review finding
F09, pre-registered in docs/review/preregistration-f09-2026-11-16.md): the contact term
gets the temperature shape of Woods-Robinson et al. (2019, section 4.3) below 150 K,

    K = K_c(rho) [g(T) + chi (T / 350)^3],

with g = 1 above 150 K, g = k_am(T) / k_am(150 K) between 20 and 150 K (k_am: their fit to
amorphous SiO2) and linear to zero below 20 K, where the source calls its model
nonphysical. The sign of G_am is taken as negative: as printed (positive) the fit gives
6.4 W m-1 K-1 at 300 K, against 1-2 in the source's own Figure 6c; with the negative sign
it gives 1.38. This is a sensitivity option: Hayne2017's K_s, K_d and chi, fitted with a
temperature-independent contact term, are not refitted.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numba import njit, vectorize
from pydantic import BaseModel, ConfigDict, Field, model_validator

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4

#: Lowest temperature at which the source describes the c_p fit as valid.
CP_VALID_MIN_K = 90.0
#: The solver refuses temperatures below this (c_p is positive above 1.3 K).
CP_POSITIVE_MIN_K = 5.0
_REFERENCE_T_K = 350.0
#: Lowest temperature of the Woods-Robinson 2019 estimate as stated by DSNE.
WOODS_ROBINSON_MIN_K = 10.0

#: DSNE (SLS-SPEC-159 Rev I) Table 3.4.6.2.2-1: Woods-Robinson et al. (2019) c_p.
WOODS_ROBINSON_T_K = np.array(
    [
        0.0,
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        80,
        91,
        100,
        120,
        140,
        160,
        181,
        200,
        220,
        240,
        260,
        280,
        300,
    ],
    dtype=np.float64,
)
WOODS_ROBINSON_CP = 1000.0 * np.array(
    [
        0.0,
        0.009,
        0.027,
        0.053,
        0.089,
        0.122,
        0.158,
        0.191,
        0.227,
        0.258,
        0.287,
        0.347,
        0.404,
        0.458,
        0.513,
        0.564,
        0.609,
        0.653,
        0.698,
        0.733,
        0.773,
    ],
    dtype=np.float64,
)
# enthalpy at each table temperature (trapezoids: exact for the linear interpolant)
_WR_E = np.concatenate(
    (
        [0.0],
        np.cumsum(
            np.diff(WOODS_ROBINSON_T_K) * 0.5 * (WOODS_ROBINSON_CP[1:] + WOODS_ROBINSON_CP[:-1])
        ),
    )
)
_WR_TOP = 300.0


#: Woods-Robinson et al. (2019) k_am(T) coefficients A..I of
#: A + B T^-4 + C T^-3 + D T^-2 + E T^-1 + F T + G T^2 + H T^3 + I T^4 (W m-1 K-1);
#: G with the sign inferred from the source's Figure 6c (see module docstring).
WR_KAM = np.array(
    [
        -2.03297e-1,
        -11.472,
        22.5793,
        -14.3084,
        3.41742,
        0.01101,
        -2.80491e-5,
        3.35837e-8,
        -1.40021e-11,
    ],
    dtype=np.float64,
)
#: Above this temperature the low-temperature option equals the default.
WR_K_TOP_K = 150.0
#: Below this temperature g(T) is taken linear to zero.
WR_K_LOW_K = 20.0


@njit(cache=True)
def _kam(t):  # type: ignore[no-untyped-def]
    c = WR_KAM
    return (
        c[0]
        + c[1] / t**4
        + c[2] / t**3
        + c[3] / t**2
        + c[4] / t
        + t * (c[5] + t * (c[6] + t * (c[7] + t * c[8])))
    )


@njit(cache=True)
def _kam_dt(t):  # type: ignore[no-untyped-def]
    c = WR_KAM
    return (
        -4.0 * c[1] / t**5
        - 3.0 * c[2] / t**4
        - 2.0 * c[3] / t**3
        - c[4] / t**2
        + c[5]
        + t * (2.0 * c[6] + t * (3.0 * c[7] + t * 4.0 * c[8]))
    )


@njit(cache=True)
def _kam_integral(t):  # type: ignore[no-untyped-def]
    c = WR_KAM
    return (
        c[0] * t
        - c[1] / (3.0 * t**3)
        - c[2] / (2.0 * t**2)
        - c[3] / t
        + c[4] * math.log(t)
        + t * t * (c[5] / 2.0 + t * (c[6] / 3.0 + t * (c[7] / 4.0 + t * c[8] / 5.0)))
    )


@njit(cache=True)
def low_temperature_factor(temperature):  # type: ignore[no-untyped-def]
    """g(T): Woods-Robinson 2019 shape of the contact conductivity, 1 above 150 K."""
    t = temperature
    if t >= WR_K_TOP_K:
        return 1.0
    top = _kam(WR_K_TOP_K)
    if t >= WR_K_LOW_K:
        return _kam(t) / top
    return _kam(WR_K_LOW_K) / top * t / WR_K_LOW_K


@njit(cache=True)
def _low_temperature_factor_dt(t):  # type: ignore[no-untyped-def]
    if t >= WR_K_TOP_K:
        return 0.0
    top = _kam(WR_K_TOP_K)
    if t >= WR_K_LOW_K:
        return _kam_dt(t) / top
    return _kam(WR_K_LOW_K) / top / WR_K_LOW_K


@njit(cache=True)
def _low_temperature_factor_integral(t):  # type: ignore[no-untyped-def]
    """Integral of g from 0 to T."""
    top = _kam(WR_K_TOP_K)
    g_low = _kam(WR_K_LOW_K) / top
    if t < WR_K_LOW_K:
        return 0.5 * g_low * t * t / WR_K_LOW_K
    below = 0.5 * g_low * WR_K_LOW_K
    if t < WR_K_TOP_K:
        return below + (_kam_integral(t) - _kam_integral(WR_K_LOW_K)) / top
    return below + (_kam_integral(WR_K_TOP_K) - _kam_integral(WR_K_LOW_K)) / top + t - WR_K_TOP_K


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
    albedo_quadratic_deg: tuple[float, float, float] | None = Field(
        None,
        description="c0, c1, c2 of an albedo quadratic in incidence angle in degrees, "
        "A = c0 + c1 i + c2 i^2 (Foote2020 form); replaces A0, a and b when set",
    )
    heat_flow_w_m2: float = Field(0.018, ge=0.0, description="Q, interior heat flow")
    heat_capacity_model: Literal["hayne2017", "woods_robinson2019"] = Field(
        "hayne2017",
        description="c_p law: Hayne2017 polynomial (heat_capacity_coefficients) or the "
        "Woods-Robinson 2019 table of DSNE, polynomial above 300 K",
    )
    conductivity_model: Literal["hayne2017", "woods_robinson2019"] = Field(
        "hayne2017",
        description="contact-conductivity temperature law: Hayne2017 (constant) or the "
        "Woods-Robinson 2019 shape below 150 K (sensitivity option, review finding F09)",
    )

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
                1.0 if self.heat_capacity_model == "woods_robinson2019" else 0.0,
                1.0 if self.conductivity_model == "woods_robinson2019" else 0.0,
            ],
            dtype=np.float64,
        )


# indices into RegolithProperties.as_array()
(P_KS, P_KD, P_RHOS, P_RHOD, P_H, P_CHI, P_C0, P_C1, P_C2, P_C3, P_C4, P_EMIS, P_Q, P_CPM, P_KM) = (
    range(15)
)


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


@vectorize(["float64(float64, float64)"], cache=True)  # type: ignore[no-untyped-call]
def _contact(temperature, low_t):  # type: ignore[no-untyped-def]
    return low_temperature_factor(temperature) if low_t > 0.5 else 1.0


@vectorize(["float64(float64, float64)"], cache=True)  # type: ignore[no-untyped-call]
def _contact_dt(temperature, low_t):  # type: ignore[no-untyped-def]
    return _low_temperature_factor_dt(temperature) if low_t > 0.5 else 0.0


@vectorize(["float64(float64, float64)"], cache=True)  # type: ignore[no-untyped-call]
def _contact_integral(temperature, low_t):  # type: ignore[no-untyped-def]
    return _low_temperature_factor_integral(temperature) if low_t > 0.5 else temperature


@njit(cache=True)
def conductivity(temperature, k_contact, chi, low_t=0.0):  # type: ignore[no-untyped-def]
    """Thermal conductivity (W m-1 K-1) including the radiative T^3 term.

    ``low_t`` > 0.5 selects the Woods-Robinson 2019 low-temperature contact term.
    """
    ratio = temperature / _REFERENCE_T_K
    contact = _contact(temperature, low_t)
    return k_contact * (contact + chi * ratio * ratio * ratio)


@njit(cache=True)
def conductivity_dt(temperature, k_contact, chi, low_t=0.0):  # type: ignore[no-untyped-def]
    """Derivative of the conductivity with respect to temperature."""
    contact = _contact_dt(temperature, low_t)
    return k_contact * (contact + chi * 3.0 * temperature * temperature / _REFERENCE_T_K**3)


@njit(cache=True)
def kirchhoff(temperature, k_contact, chi, low_t=0.0):  # type: ignore[no-untyped-def]
    """Kirchhoff potential: integral of K from 0 to T (W m-1)."""
    contact = _contact_integral(temperature, low_t)
    return k_contact * (contact + chi * temperature**4 / (4.0 * _REFERENCE_T_K**3))


@njit(cache=True)
def _polynomial_cp(t, p):  # type: ignore[no-untyped-def]
    return p[P_C0] + t * (p[P_C1] + t * (p[P_C2] + t * (p[P_C3] + t * p[P_C4])))


@njit(cache=True)
def _polynomial_e(t, p):  # type: ignore[no-untyped-def]
    return t * (
        p[P_C0]
        + t * (p[P_C1] / 2.0 + t * (p[P_C2] / 3.0 + t * (p[P_C3] / 4.0 + t * p[P_C4] / 5.0)))
    )


@njit(cache=True)
def _table_index(t):  # type: ignore[no-untyped-def]
    i = np.searchsorted(WOODS_ROBINSON_T_K, t) - 1
    return min(max(i, 0), WOODS_ROBINSON_T_K.size - 2)


@njit(cache=True)
def heat_capacity(temperature, p):  # type: ignore[no-untyped-def]
    """Specific heat capacity c_p (J kg-1 K-1)."""
    t = temperature
    if p[P_CPM] < 0.5 or t >= _WR_TOP:
        return _polynomial_cp(t, p)
    i = _table_index(t)
    t0 = WOODS_ROBINSON_T_K[i]
    slope = (WOODS_ROBINSON_CP[i + 1] - WOODS_ROBINSON_CP[i]) / (WOODS_ROBINSON_T_K[i + 1] - t0)
    return WOODS_ROBINSON_CP[i] + slope * (t - t0)


@njit(cache=True)
def enthalpy(temperature, p):  # type: ignore[no-untyped-def]
    """Specific enthalpy (J kg-1): integral of c_p from 0 to T."""
    t = temperature
    if p[P_CPM] < 0.5:
        return _polynomial_e(t, p)
    if t >= _WR_TOP:
        return _WR_E[-1] + _polynomial_e(t, p) - _polynomial_e(_WR_TOP, p)
    i = _table_index(t)
    t0 = WOODS_ROBINSON_T_K[i]
    return _WR_E[i] + 0.5 * (t - t0) * (WOODS_ROBINSON_CP[i] + heat_capacity(t, p))


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

    With ``albedo_quadratic_deg`` set, the law is instead a quadratic in the incidence
    angle in degrees (the form of Foote2020's laboratory fits), clipped to [0, 1].
    """
    theta = np.clip(np.asarray(incidence_rad, dtype=np.float64), 0.0, np.pi / 2.0)
    if properties.albedo_quadratic_deg is not None:
        c0, c1, c2 = properties.albedo_quadratic_deg
        degrees = np.degrees(theta)
        return np.asarray(np.clip(c0 + degrees * (c1 + degrees * c2), 0.0, 1.0))
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
