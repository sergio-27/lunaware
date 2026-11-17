"""Sub-DEM roughness: Gaussian slope facets with Smith (1967) cast shadows.

Validation only (TH-V3f, pre-registered in
``docs/review/preregistration-th-v3f-2026-11-17.md``); the engine does not use it.

The model follows Bandfield et al. (2015, Icarus 248, 357, section 2): below the DEM
scale the surface is a Gaussian random surface with RMS slope ``theta_rms``
(``s = tan(theta_rms)``, ``s^2 = <|grad h|^2>``, random azimuths). Each facet runs the
1-D thermal model with its own direct flux. A front-facing facet is lit with the Smith
(1967) probability ``G1``, the same for every facet orientation, so while it faces a Sun
above the mean plane each facet is split into a lit population (its own thermal history)
and a cast-shadow population (a run with no direct sunlight); otherwise it follows its
own history. Sloped facets also receive sunlight and infrared
from adjacent surfaces in proportion to the fraction ``(1 - cos beta) / 2`` of their
hemisphere filled, taken as a surface parallel to the DEM plane with the smooth facet's
incident sunlight and temperature. The view is nadir along the DEM normal, so facets
are weighted by their probability per unit area of the mean plane.

Not included: lateral conduction between facets, correlation between shadowing and
orientation, multiple scattering beyond the adjacent-surface term, oblique views.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4
from lunarenv.geometry.projection import FloatArray
from lunarenv.thermal.model import ThermalSettings, run_thermal
from lunarenv.thermal.properties import RegolithProperties, albedo
from lunarenv.thermal.terrain_radiation import diffuse_albedo


@dataclass(frozen=True, slots=True)
class FacetSet:
    """Quadrature of the slope distribution: gradients in the mean-plane frame.

    ``gradient`` (k, 2), ``weight`` (k,) summing to 1, ``node`` (k,) index of the slope
    node of each facet (facets of one node share a slope and differ in azimuth).
    """

    gradient: FloatArray
    weight: FloatArray
    node: np.ndarray
    n_nodes: int


def facet_set(rms_slope_deg: float, n_nodes: int = 3, n_azimuth: int = 12) -> FacetSet:
    """Gauss-Laguerre nodes in ``|grad h|^2 / s^2`` (exponentially distributed for an
    isotropic Gaussian surface) times ``n_azimuth`` equally spaced azimuths."""
    if rms_slope_deg <= 0.0:
        raise ValueError("RMS slope must be positive")
    s = float(np.tan(np.radians(rms_slope_deg)))
    x, w = np.polynomial.laguerre.laggauss(n_nodes)
    azimuth = (np.arange(n_azimuth) + 0.5) * 2.0 * np.pi / n_azimuth
    magnitude = s * np.sqrt(x)
    gradient = np.stack(
        [
            np.outer(magnitude, np.cos(azimuth)).ravel(),
            np.outer(magnitude, np.sin(azimuth)).ravel(),
        ],
        axis=1,
    )
    weight = np.repeat(w / n_azimuth, n_azimuth)
    node = np.repeat(np.arange(n_nodes), n_azimuth)
    return FacetSet(gradient, weight / weight.sum(), node, n_nodes)


def smith_lambda(nu: ArrayLike) -> FloatArray:
    """Smith (1967) shadowing term for a Gaussian surface, ``nu = cot(theta) / s``
    with ``s`` the RMS of the total gradient (``sqrt(2)`` times one component's)."""
    v = np.asarray(nu, dtype=np.float64)
    erfc = np.vectorize(math.erfc, otypes=[np.float64])(v)
    return np.asarray(0.5 * (np.exp(-(v**2)) / (v * np.sqrt(np.pi)) - erfc))


def lit_probability(cos_zenith: ArrayLike, rms_slope_deg: float) -> FloatArray:
    """Probability that a front-facing facet is outside cast shadows (Smith 1967),
    the Sun at zenith angle ``arccos(cos_zenith)`` from the mean-plane normal."""
    c = np.asarray(cos_zenith, dtype=np.float64)
    s = float(np.tan(np.radians(rms_slope_deg)))
    out = np.zeros_like(c)
    up = c > 0.0
    cz = np.minimum(c[up], 1.0)
    with np.errstate(divide="ignore"):
        nu = cz / np.sqrt(np.maximum(1.0 - cz**2, 0.0)) / s
    lam = np.where(np.isfinite(nu), smith_lambda(np.where(np.isfinite(nu), nu, 1.0)), 0.0)
    out[up] = 1.0 / (1.0 + lam)
    return out


def facet_normals(normal0: ArrayLike, gradient: FloatArray) -> FloatArray:
    """Unit normals (k, 3) of facets with ``gradient`` relative to the plane normal to
    ``normal0`` (grid frame: x, y, up)."""
    n0 = np.asarray(normal0, dtype=np.float64)
    n0 = n0 / np.linalg.norm(n0)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, n0)
    u /= np.linalg.norm(u)
    v = np.cross(n0, u)
    n = n0[None, :] - gradient[:, :1] * u[None, :] - gradient[:, 1:2] * v[None, :]
    return np.asarray(n / np.linalg.norm(n, axis=1, keepdims=True))


def sun_vectors(elevation_rad: ArrayLike, azimuth_rad: ArrayLike) -> FloatArray:
    """Unit Sun directions (t, 3) in the grid frame (as ``surface.incidence_cosine``)."""
    el = np.asarray(elevation_rad, dtype=np.float64)
    az = np.asarray(azimuth_rad, dtype=np.float64)
    return np.stack([np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)], axis=1)


@dataclass(frozen=True, slots=True)
class RoughResult:
    """Period radiance ``eps * T^4`` (K^4) of the smooth facet and of each rough surface."""

    smooth: FloatArray
    rough: dict[float, FloatArray]


def rough_radiance(
    irradiance_w_m2: ArrayLike,
    sun: FloatArray,
    normal0: ArrayLike,
    properties: RegolithProperties,
    rms_slopes_deg: Sequence[float],
    step_s: float,
    spinup_samples: int,
    settings: ThermalSettings | None = None,
    *,
    n_nodes: int = 3,
    n_azimuth: int = 12,
) -> RoughResult:
    """Period radiance of the smooth DEM facet and of the rough surface at each RMS slope.

    ``irradiance_w_m2`` is the normal irradiance of the visible solar disc (including any
    eclipse factor) and ``sun`` the unit Sun directions, both over spin-up and period.
    The adjacent-surface infrared uses the smooth facet's temperature over the period and,
    over the spin-up (where the solver does not return temperatures), its radiative
    equilibrium with the absorbed flux.
    """
    settings = settings or ThermalSettings()
    eps = properties.emissivity
    e = np.asarray(irradiance_w_m2, dtype=np.float64)
    n0 = np.asarray(normal0, dtype=np.float64)
    n0 = n0 / np.linalg.norm(n0)
    cos0 = np.clip(sun @ n0, 0.0, 1.0)
    incident0 = e * cos0
    flux0 = (1.0 - albedo(np.arccos(cos0), properties)) * incident0
    smooth = run_thermal(flux0, step_s, spinup_samples, properties, settings)
    t0 = smooth.surface_temperature_k
    smooth_rad = eps * t0**4
    emitted0 = np.concatenate([flux0[:spinup_samples], eps * STEFAN_BOLTZMANN_W_M2_K4 * t0**4])
    reflected0 = albedo(np.arccos(cos0), properties) * incident0
    neighbour = (1.0 - diffuse_albedo(properties)) * reflected0 + eps * emitted0
    period = slice(spinup_samples, None)
    rough: dict[float, FloatArray] = {}
    for rms in rms_slopes_deg:
        facets = facet_set(rms, n_nodes, n_azimuth)
        g1 = lit_probability(sun[period] @ n0, rms)
        normals = facet_normals(n0, facets.gradient)
        slope = np.arctan(np.hypot(facets.gradient[:, 0], facets.gradient[:, 1]))
        shadow_rad = []
        diffuse = []
        for j in range(facets.n_nodes):
            beta = slope[facets.node == j][0]
            d = 0.5 * (1.0 - np.cos(beta)) * neighbour
            diffuse.append(d)
            t = run_thermal(d, step_s, spinup_samples, properties, settings).surface_temperature_k
            shadow_rad.append(eps * t**4)
        total = np.zeros_like(smooth_rad)
        for k in range(facets.weight.size):
            cos_k = np.clip(sun @ normals[k], 0.0, 1.0)
            flux = (1.0 - albedo(np.arccos(cos_k), properties)) * e * cos_k
            flux = flux + diffuse[facets.node[k]]
            facet = run_thermal(flux, step_s, spinup_samples, properties, settings)
            lit = eps * facet.surface_temperature_k**4
            # the lit / cast-shadow split applies while the facet faces a Sun above the mean
            # plane; otherwise no part of it is directly lit and it follows its own history
            front = (cos_k[period] > 0.0) & (g1 > 0.0)
            mixed = g1 * lit + (1.0 - g1) * shadow_rad[facets.node[k]]
            total += facets.weight[k] * np.where(front, mixed, lit)
        rough[float(rms)] = total
    return RoughResult(smooth_rad, rough)
