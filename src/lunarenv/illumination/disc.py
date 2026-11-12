"""Visible fraction of a luminous disc (Sun or Earth) above a terrain horizon.

The disc is cut into vertical strips. Strip ``j`` sits at horizontal angular offset
``xi_j = rho sin(theta_j)`` from the disc centre, with ``theta_j`` equally spaced in
(-pi/2, pi/2); its half-height is ``c_j = rho cos(theta_j)``. In each strip the terrain
hides everything below the horizon elevation at that strip's azimuth. Integrating in
``theta`` removes the square-root end-point behaviour of the chord length, so the
midpoint rule is exact for an unobstructed disc and converges quickly otherwise.

Two horizon modes (PLAN.md decision D5):

* ``PROFILE``: each strip uses the horizon interpolated at its own azimuth. This
  resolves terrain structure across the 0.53 deg width of the Sun.
* ``CHORD``: every strip uses the horizon at the azimuth of the disc centre, i.e. the
  disc is cut by one horizontal line. This is the approximation of Mazarico2011 and
  Barker2021 and is offered to reproduce their numbers.

Limb darkening. Intensity follows the linear law ``I(mu)/I(1) = 1 - u (1 - mu)`` with
``mu`` the cosine of the emergent angle. ``u = 0.6`` is the Eddington approximation of
a grey atmosphere, ``I(mu) = (2 + 3 mu) / 5``, appropriate for the wavelength-integrated
Sun; ``u = 0`` gives a uniform disc. With this law every strip integral has a closed
form, so the flux fraction costs the same as the geometric one.

Approximations: the sky around the disc is treated as flat (errors of order rho^2,
about 2e-5), and the horizon is linear between azimuth samples.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit
from numpy.typing import NDArray

PROFILE = 0
CHORD = 1

#: Linear limb-darkening coefficient of the Eddington grey atmosphere.
EDDINGTON_LIMB_DARKENING = 0.6


@njit(cache=True, inline="always")
def interpolate_horizon(horizon, azimuth):  # type: ignore[no-untyped-def]
    """Horizon elevation at a grid azimuth, linear and periodic between samples."""
    n = horizon.size
    position = azimuth * n / (2.0 * math.pi)
    position -= math.floor(position / n) * n
    i0 = int(position)
    if i0 >= n:  # guards position == n after round-off
        i0 = 0
        position = 0.0
    fraction = position - i0
    i1 = i0 + 1
    if i1 == n:
        i1 = 0
    return horizon[i0] * (1.0 - fraction) + horizon[i1] * fraction


@njit(cache=True, inline="always")
def _chord_primitive(y, c):  # type: ignore[no-untyped-def]
    """Primitive of sqrt(c^2 - y^2) for |y| <= c."""
    ratio = y / c
    if ratio > 1.0:
        ratio = 1.0
    elif ratio < -1.0:
        ratio = -1.0
    return 0.5 * (y * math.sqrt(max(c * c - y * y, 0.0)) + c * c * math.asin(ratio))


@njit(cache=True)
def disc_fractions(elevation, azimuth, radius, horizon, n_strips, limb_u, mode):  # type: ignore[no-untyped-def]
    """Geometric and flux-weighted visible fractions of a disc.

    ``elevation`` and ``azimuth`` locate the disc centre (radians, grid azimuth),
    ``radius`` is its angular radius, ``horizon`` the horizon profile of the observer.
    Returns ``(geometric, flux)``, both in [0, 1].
    """
    cos_elevation = max(math.cos(elevation), 1e-9)
    area = 0.0
    flux = 0.0
    full_area = 0.0
    full_flux = 0.0
    d_theta = math.pi / n_strips
    for j in range(n_strips):
        theta = -0.5 * math.pi + (j + 0.5) * d_theta
        c = radius * math.cos(theta)
        weight = c * d_theta
        # the same quadrature applied to the unobstructed strip gives the normalization,
        # so a disc clear of the horizon yields exactly 1 and quadrature errors cancel
        full_area += weight * 2.0 * c
        full_flux += weight * ((1.0 - limb_u) * 2.0 * c + limb_u * 0.5 * math.pi * c * c / radius)
        if mode == PROFILE:
            strip_azimuth = azimuth + radius * math.sin(theta) / cos_elevation
        else:
            strip_azimuth = azimuth
        y_low = interpolate_horizon(horizon, strip_azimuth) - elevation
        if y_low >= c:
            continue  # strip fully hidden
        if y_low < -c:
            y_low = -c
        visible = c - y_low
        area += weight * visible
        if limb_u != 0.0:
            mu_integral = (0.25 * math.pi * c * c - _chord_primitive(y_low, c)) / radius
            flux += weight * ((1.0 - limb_u) * visible + limb_u * mu_integral)
    geometric = min(max(area / full_area, 0.0), 1.0)
    if limb_u == 0.0:
        return geometric, geometric
    return geometric, min(max(flux / full_flux, 0.0), 1.0)


def disc_fractions_array(
    elevation: NDArray[np.float64],
    azimuth: NDArray[np.float64],
    radius: NDArray[np.float64],
    horizon: NDArray[np.float32] | NDArray[np.float64],
    *,
    n_strips: int = 32,
    limb_u: float = EDDINGTON_LIMB_DARKENING,
    mode: int = PROFILE,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Vectorized convenience wrapper around :func:`disc_fractions` (not for hot loops)."""
    elevation, azimuth, radius = np.broadcast_arrays(elevation, azimuth, radius)
    geometric = np.empty(elevation.shape)
    flux = np.empty(elevation.shape)
    profile = np.ascontiguousarray(horizon, dtype=np.float64)
    for index in np.ndindex(elevation.shape):
        geometric[index], flux[index] = disc_fractions(
            float(elevation[index]),
            float(azimuth[index]),
            float(radius[index]),
            profile,
            n_strips,
            limb_u,
            mode,
        )
    return geometric, flux
