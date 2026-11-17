"""Diviner footprint and binning applied to a model radiance field (TH-V3g).

Pre-registered in ``docs/review/preregistration-th-v3g-2026-11-17.md``. Validation only.

A Diviner observation measures the mean radiance over its footprint: the geometric IFOV
rectangle (6.7 x 3.4 mrad, Diviner RDR SIS Table 1) smeared in track by the motion during
the 0.128 s integration. The gridded products assign each observation to every bin that
one of nine points of its IFOV rectangle falls in (centre, four edge midpoints, four
corners; SIS section 5.4), each point counting once. With observation centres uniform
over the area, a bin's value is the mean over the nine points ``p`` of the mean
observation radiance over centres in the bin shifted by ``-d_p``.

The field is a regular grid in polar stereographic metres (``x`` east, ``y`` north of
the projection), one value per subsolar-longitude bin.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M

#: Standard lunar GM (m^3 s^-2), for the ground speed only.
GM_MOON = 4.9028e12
INTEGRATION_S = 0.128
IFOV_IN_TRACK_RAD = 6.7e-3
IFOV_CROSS_TRACK_RAD = 3.4e-3


def ground_speed(altitude_m: float) -> float:
    """Speed of the sub-spacecraft point on a circular orbit at ``altitude_m``."""
    r = MOON_REFERENCE_RADIUS_M + altitude_m
    return float(np.sqrt(GM_MOON / r) * MOON_REFERENCE_RADIUS_M / r)


@dataclass(frozen=True, slots=True)
class Footprint:
    """Geometric IFOV rectangle (m) and in-track smear during one integration (m)."""

    in_track_m: float
    cross_track_m: float
    smear_m: float

    @classmethod
    def at_altitude(cls, altitude_m: float = 50e3) -> Footprint:
        return cls(
            IFOV_IN_TRACK_RAD * altitude_m,
            IFOV_CROSS_TRACK_RAD * altitude_m,
            INTEGRATION_S * ground_speed(altitude_m),
        )

    def scaled(self, factor: float) -> Footprint:
        return Footprint(
            self.in_track_m * factor, self.cross_track_m * factor, self.smear_m * factor
        )


def _axes(theta_rad: float) -> tuple[np.ndarray, np.ndarray]:
    u = np.array([np.cos(theta_rad), np.sin(theta_rad)])  # in track
    return u, np.array([-u[1], u[0]])  # cross track


def _box_samples(length_m: float, step_m: float) -> np.ndarray:
    if length_m <= 0.0:
        return np.zeros(1)
    n = max(1, int(np.ceil(length_m / step_m)))
    return (np.arange(n) + 0.5) * (length_m / n) - length_m / 2


def response(
    footprint: Footprint, theta_rad: float, step_m: float = 10.0
) -> tuple[np.ndarray, np.ndarray]:
    """Offsets (k, 2) and normalised weights (k,) of the footprint response.

    Uniform cross track; in track the IFOV box convolved with the smear box (a
    trapezoid), sampled at the midpoints of about ``step_m`` intervals.
    """
    length, smear = footprint.in_track_m, footprint.smear_m
    s = _box_samples(length + smear, step_m)
    if smear > 0.0 and length > 0.0:
        w_in = np.clip(
            np.minimum(length / 2, s + smear / 2) - np.maximum(-length / 2, s - smear / 2),
            0.0,
            None,
        )
    else:
        w_in = np.ones_like(s)
    c = _box_samples(footprint.cross_track_m, step_m)
    u, v = _axes(theta_rad)
    offsets = (s[:, None, None] * u + c[None, :, None] * v).reshape(-1, 2)
    weights = np.repeat(w_in, c.size).astype(np.float64)
    return offsets, weights / weights.sum()


def bin_points(footprint: Footprint, theta_rad: float) -> np.ndarray:
    """The nine binning points (9, 2) of the geometric IFOV rectangle (no smear)."""
    u, v = _axes(theta_rad)
    a, b = footprint.in_track_m / 2, footprint.cross_track_m / 2
    return np.array([i * a * u + j * b * v for i in (0, -1, 1) for j in (0, -1, 1)])


@dataclass(frozen=True, slots=True)
class Field:
    """Values (ny, nx, nb) at nodes ``x0 + step * i``, ``y0 + step * j``."""

    x0: float
    y0: float
    step: float
    values: np.ndarray

    def sample(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Bilinear interpolation (n, nb); points must lie inside the node grid."""
        ny, nx = self.values.shape[:2]
        fx = (np.asarray(x) - self.x0) / self.step
        fy = (np.asarray(y) - self.y0) / self.step
        tol = 1e-9
        if fx.min() < -tol or fy.min() < -tol or fx.max() > nx - 1 + tol or fy.max() > ny - 1 + tol:
            raise ValueError("footprint reaches outside the model field")
        fx = np.clip(fx, 0.0, nx - 1)
        fy = np.clip(fy, 0.0, ny - 1)
        i0 = np.minimum(fx.astype(np.int64), nx - 2)
        j0 = np.minimum(fy.astype(np.int64), ny - 2)
        tx = (fx - i0)[:, None]
        ty = (fy - j0)[:, None]
        v = self.values
        return np.asarray(
            (1 - ty) * ((1 - tx) * v[j0, i0] + tx * v[j0, i0 + 1])
            + ty * ((1 - tx) * v[j0 + 1, i0] + tx * v[j0 + 1, i0 + 1])
        )


def cell_radiance(
    field: Field,
    centres_x: np.ndarray,
    centres_y: np.ndarray,
    footprint: Footprint,
    theta_rad: float,
    *,
    binning: Footprint | None = None,
    nine_points: bool = True,
    step_m: float = 10.0,
) -> np.ndarray:
    """Model bin value (nb,) of one cell whose observation centres are sampled at
    ``centres_x``, ``centres_y`` (uniform over the cell).

    ``footprint`` sets the response; ``binning`` (default: the same) sets the nine
    binning points; ``nine_points=False`` bins by the centre only.
    """
    offsets, weights = response(footprint, theta_rad, step_m)
    points = bin_points(binning or footprint, theta_rad) if nine_points else np.zeros((1, 2))
    total = np.zeros(field.values.shape[2])
    for d in points:
        x = (centres_x - d[0])[:, None] + offsets[None, :, 0]
        y = (centres_y - d[1])[:, None] + offsets[None, :, 1]
        values = field.sample(x.ravel(), y.ravel()).reshape(x.shape[0], x.shape[1], -1)
        total += np.einsum("ikb,k->b", values, weights) / x.shape[0]
    return total / len(points)


@dataclass(frozen=True, slots=True)
class Layout:
    """Model sampling around a block: 20 m nodes inside, coarser pixels in a ring."""

    x0: float  # first node
    y0: float
    n: int  # nodes per side
    inner_m: float
    ring_m: float
    block_first: int  # index of the first node inside the block
    block_nodes: int

    @property
    def factor(self) -> int:
        return round(self.ring_m / self.inner_m)

    def inner_points(self) -> tuple[np.ndarray, np.ndarray]:
        idx = np.arange(self.block_first, self.block_first + self.block_nodes)
        xx, yy = np.meshgrid(self.x0 + self.inner_m * idx, self.y0 + self.inner_m * idx)
        return xx.ravel(), yy.ravel()

    def ring_points(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Centres and ring-pixel indices (jx, jy) of the ring pixels outside the block."""
        f = self.factor
        m = self.n // f
        j = np.arange(m)
        jx, jy = np.meshgrid(j, j)
        lo, hi = self.block_first // f, (self.block_first + self.block_nodes) // f
        outside = ~((jx >= lo) & (jx < hi) & (jy >= lo) & (jy < hi))
        jx, jy = jx[outside], jy[outside]
        origin_x = self.x0 - self.inner_m / 2
        origin_y = self.y0 - self.inner_m / 2
        return (
            origin_x + self.ring_m * (jx + 0.5),
            origin_y + self.ring_m * (jy + 0.5),
            jx,
            jy,
        )

    def assemble(self, inner: np.ndarray, ring: np.ndarray) -> Field:
        """Field from inner values (block_nodes^2, nb) and ring values (n_ring, nb)."""
        nb = inner.shape[1]
        values = np.zeros((self.n, self.n, nb))
        covered = np.zeros((self.n, self.n), dtype=bool)
        _, _, jx, jy = self.ring_points()
        f = self.factor
        for a in range(f):
            for b in range(f):
                values[jy * f + b, jx * f + a] = ring
                covered[jy * f + b, jx * f + a] = True
        s = slice(self.block_first, self.block_first + self.block_nodes)
        values[s, s] = inner.reshape(self.block_nodes, self.block_nodes, nb)
        covered[s, s] = True
        if not covered.all():
            raise ValueError("field not covered")
        return Field(self.x0, self.y0, self.inner_m, values)


def block_layout(
    xmin_m: float,
    ymin_m: float,
    block_m: float = 720.0,
    inner_m: float = 20.0,
    ring_m: float = 60.0,
    margin_m: float = 660.0,
) -> Layout:
    """Layout for a square block with lower-left corner (``xmin_m``, ``ymin_m``)."""
    f = ring_m / inner_m
    for value in (f, block_m / ring_m, margin_m / ring_m):
        if abs(value - round(value)) > 1e-9:
            raise ValueError("block, margin and ring spacing must be multiples")
    n = round((block_m + 2 * margin_m) / inner_m)
    return Layout(
        xmin_m - margin_m + inner_m / 2,
        ymin_m - margin_m + inner_m / 2,
        n,
        inner_m,
        ring_m,
        round(margin_m / inner_m),
        round(block_m / inner_m),
    )
