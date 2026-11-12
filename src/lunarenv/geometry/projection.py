"""South polar stereographic projection on a sphere, and its link to 3-D coordinates.

Conventions (those of the LOLA polar products, verified against PROJ in the tests):

* sphere of radius ``R``; projection plane tangent at the south pole (scale 1 there);
* longitude 0 points along map ``+y``, longitude 90 E along map ``+x``;
* east longitude, planetocentric latitude, radians;
* height ``h`` is radial distance above the sphere.

Body-fixed Cartesian axes follow the mean-Earth/polar-axis frame: ``+Z`` north pole,
``+X`` at longitude 0, ``+Y`` at longitude 90 E.

With ``rho = hypot(x, y)`` and ``t = rho / (2R)`` (the tangent of half the angular
distance from the pole), every relation below is algebraic in ``x, y``. Nothing is
singular at the pole, which lies inside the MVP test tile.

Simplification: the Moon is treated as a sphere, as in the LOLA products themselves.
Planetocentric and planetodetic latitude therefore coincide and the local vertical is
radial.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class SouthPolarStereographic:
    """South polar stereographic projection on a sphere, true scale at the pole."""

    radius_m: float = MOON_REFERENCE_RADIUS_M

    def __post_init__(self) -> None:
        if not self.radius_m > 0.0:
            raise ValueError("radius must be positive")

    # ------------------------------------------------------------ lon/lat <-> map
    def forward(self, lon_rad: ArrayLike, lat_rad: ArrayLike) -> tuple[FloatArray, FloatArray]:
        """Longitude/latitude to map coordinates in metres.

        The north pole is the projection's point at infinity and must not be passed.
        """
        lon = np.asarray(lon_rad, dtype=np.float64)
        lat = np.asarray(lat_rad, dtype=np.float64)
        rho = 2.0 * self.radius_m * np.tan(np.pi / 4.0 + lat / 2.0)
        return rho * np.sin(lon), rho * np.cos(lon)

    def inverse(self, x_m: ArrayLike, y_m: ArrayLike) -> tuple[FloatArray, FloatArray]:
        """Map coordinates to longitude in (-pi, pi] and latitude.

        At the pole the longitude is undefined; 0 is returned.
        """
        x = np.asarray(x_m, dtype=np.float64)
        y = np.asarray(y_m, dtype=np.float64)
        rho = np.hypot(x, y)
        lat = 2.0 * np.arctan(rho / (2.0 * self.radius_m)) - np.pi / 2.0
        return np.arctan2(x, y), lat

    def scale_factor(self, x_m: ArrayLike, y_m: ArrayLike) -> FloatArray:
        """Map distance per unit distance on the sphere, ``k = 1 + t**2``.

        The projection is conformal, so ``k`` is the same in every direction.
        """
        x = np.asarray(x_m, dtype=np.float64)
        y = np.asarray(y_m, dtype=np.float64)
        return np.asarray(1.0 + (x * x + y * y) / (4.0 * self.radius_m**2))

    # ---------------------------------------------------------- map <-> Cartesian
    def to_cartesian(self, x_m: ArrayLike, y_m: ArrayLike, h_m: ArrayLike = 0.0) -> FloatArray:
        """Map coordinates and height to body-fixed Cartesian metres, shape ``(..., 3)``."""
        x, y, h = np.broadcast_arrays(
            np.asarray(x_m, dtype=np.float64),
            np.asarray(y_m, dtype=np.float64),
            np.asarray(h_m, dtype=np.float64),
        )
        radius = self.radius_m
        t2 = (x * x + y * y) / (4.0 * radius * radius)
        s = 1.0 / (1.0 + t2)
        r_over_radius = (radius + h) / radius
        out = np.empty((*x.shape, 3), dtype=np.float64)
        out[..., 0] = r_over_radius * y * s
        out[..., 1] = r_over_radius * x * s
        out[..., 2] = -(radius + h) * (1.0 - t2) * s
        return out

    def from_cartesian(self, position_m: ArrayLike) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Body-fixed Cartesian metres to map coordinates and height.

        Valid everywhere except at the north pole, the projection's point at infinity.
        """
        p = np.asarray(position_m, dtype=np.float64)
        r = np.linalg.norm(p, axis=-1)
        # t = tan(c/2) with c the angular distance from the south pole;
        # tan(c/2) = sin(c) / (1 + cos(c)), and -Z/r = cos(c).
        scale = 2.0 * self.radius_m / (r - p[..., 2])
        return p[..., 1] * scale, p[..., 0] * scale, r - self.radius_m

    # --------------------------------------------------------------- local frames
    def local_basis(
        self, x_m: ArrayLike, y_m: ArrayLike
    ) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Orthonormal basis ``(e_x, e_y, up)`` of the *grid-aligned* local frame.

        ``up`` is the outward radial direction; ``e_x`` and ``e_y`` are the horizontal
        unit vectors pointing along map ``+x`` and ``+y``. The triad is right-handed
        (``e_x x e_y = up``) and, unlike east/north, continuous through the pole.
        Each array has shape ``(..., 3)`` in body-fixed axes.
        """
        x, y = np.broadcast_arrays(
            np.asarray(x_m, dtype=np.float64), np.asarray(y_m, dtype=np.float64)
        )
        radius = self.radius_m
        t2 = (x * x + y * y) / (4.0 * radius * radius)
        s = 1.0 / (1.0 + t2)
        q = s / (2.0 * radius * radius)

        up = np.empty((*x.shape, 3), dtype=np.float64)
        up[..., 0] = y * s / radius
        up[..., 1] = x * s / radius
        up[..., 2] = -(1.0 - t2) * s

        e_x = np.empty_like(up)
        e_x[..., 0] = -x * y * q
        e_x[..., 1] = 1.0 - x * x * q
        e_x[..., 2] = x * s / radius

        e_y = np.empty_like(up)
        e_y[..., 0] = 1.0 - y * y * q
        e_y[..., 1] = -x * y * q
        e_y[..., 2] = y * s / radius
        return e_x, e_y, up

    def grid_convergence(self, x_m: ArrayLike, y_m: ArrayLike) -> FloatArray:
        """Angle to add to a true azimuth to obtain the grid azimuth.

        True azimuth is measured from local north towards east; grid azimuth from map
        ``+y`` towards map ``+x``. Their difference equals the longitude. Returns 0 at
        the pole, where true azimuth is undefined.
        """
        return np.arctan2(np.asarray(x_m, dtype=np.float64), np.asarray(y_m, dtype=np.float64))


def elevation_azimuth(
    direction: ArrayLike, e_x: ArrayLike, e_y: ArrayLike, up: ArrayLike
) -> tuple[FloatArray, FloatArray]:
    """Elevation and grid azimuth of a direction seen from a local frame.

    ``direction`` need not be normalized. Elevation is the angle above the local
    horizontal plane in [-pi/2, pi/2]; grid azimuth is measured from ``e_y`` towards
    ``e_x`` in [0, 2 pi).
    """
    d = np.asarray(direction, dtype=np.float64)
    d = d / np.linalg.norm(d, axis=-1, keepdims=True)
    vertical = np.sum(d * np.asarray(up), axis=-1)
    along_x = np.sum(d * np.asarray(e_x), axis=-1)
    along_y = np.sum(d * np.asarray(e_y), axis=-1)
    # atan2 keeps full precision near the horizon and near the zenith alike
    elevation = np.arctan2(vertical, np.hypot(along_x, along_y))
    azimuth = np.mod(np.arctan2(along_x, along_y), 2.0 * np.pi)
    return elevation, azimuth
