"""Occultation of the Sun by the Earth as seen from the Moon (lunar eclipses).

During a lunar eclipse the Earth's disc covers part or all of the Sun for an observer on
the Moon. The terrain horizon knows nothing about this, so the factor computed here
multiplies the solar flux obtained from the horizon model.

Model: two uniform circular discs. The visible fraction of the Sun is one minus the
area of overlap divided by the area of the solar disc.

Simplifications and their effect:

* no limb darkening of the Sun in the overlap (a few percent of the flux during
  partial phases);
* the Earth is a sphere of its mean radius with a sharp limb: no atmosphere. Sunlight
  refracted through the atmosphere, which gives the eclipsed Moon its red glow, is
  ignored (about 1e-4 of full sunlight, negligible for the energy budget), and the
  atmosphere's enlargement of the shadow (about 1-2 %, Danjon's rule in eclipse
  almanacs) is not applied, so contact times differ from almanac values by a minute
  or two.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from lunarenv.geometry.projection import FloatArray


def disc_overlap_area(
    separation: ArrayLike, radius_a: ArrayLike, radius_b: ArrayLike
) -> FloatArray:
    """Area common to two discs of angular radii ``a`` and ``b`` whose centres are
    ``separation`` apart (flat-sky geometry; all angles in the same unit)."""
    d, a, b = np.broadcast_arrays(
        np.asarray(separation, dtype=np.float64),
        np.asarray(radius_a, dtype=np.float64),
        np.asarray(radius_b, dtype=np.float64),
    )
    small = np.minimum(a, b)
    area = np.zeros(d.shape, dtype=np.float64)
    contained = d <= np.abs(a - b)
    area[contained] = np.pi * small[contained] ** 2
    partial = (~contained) & (d < a + b)
    dp, ap, bp = d[partial], a[partial], b[partial]
    alpha = np.arccos(np.clip((dp**2 + ap**2 - bp**2) / (2 * dp * ap), -1.0, 1.0))
    beta = np.arccos(np.clip((dp**2 + bp**2 - ap**2) / (2 * dp * bp), -1.0, 1.0))
    area[partial] = ap**2 * (alpha - np.sin(2 * alpha) / 2) + bp**2 * (beta - np.sin(2 * beta) / 2)
    return area


def earth_occultation_factor(
    sun_position_m: ArrayLike,
    earth_position_m: ArrayLike,
    observer_position_m: ArrayLike,
    sun_radius_m: float,
    earth_radius_m: float,
) -> FloatArray:
    """Fraction of the solar disc not covered by the Earth, per epoch (1 = no eclipse).

    Positions are Moon-centred, body-fixed, shape ``(n, 3)`` for the bodies and ``(3,)``
    for the observer.
    """
    observer = np.asarray(observer_position_m, dtype=np.float64)
    to_sun = np.asarray(sun_position_m, dtype=np.float64) - observer
    to_earth = np.asarray(earth_position_m, dtype=np.float64) - observer
    d_sun = np.linalg.norm(to_sun, axis=-1)
    d_earth = np.linalg.norm(to_earth, axis=-1)
    cross = np.linalg.norm(np.cross(to_sun, to_earth), axis=-1)
    separation = np.arctan2(cross, np.sum(to_sun * to_earth, axis=-1))
    rho_sun = np.arcsin(sun_radius_m / d_sun)
    rho_earth = np.arcsin(earth_radius_m / d_earth)
    covered = disc_overlap_area(separation, rho_sun, rho_earth) / (np.pi * rho_sun**2)
    return np.asarray(np.clip(1.0 - covered, 0.0, 1.0))
