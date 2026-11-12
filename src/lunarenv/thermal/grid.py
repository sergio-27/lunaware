"""Depth grid for the one-dimensional conduction solver.

Finite-volume layout: ``n`` cells stacked from the surface downwards. Cell ``i`` spans
``edges[i]`` to ``edges[i + 1]`` and carries one temperature at its centre. Cell
thickness grows geometrically with depth, because the temperature wave is sharp near
the surface and smooth below. The surface itself (z = 0) is a massless node handled by
the solver.

The three parameters (first-cell thickness, growth ratio, total depth) are chosen by
the convergence study in ``docs/convergence/thermal.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lunarenv.geometry.projection import FloatArray


@dataclass(frozen=True, slots=True)
class DepthGrid:
    """Cell edges (m, increasing, starting at 0)."""

    edges_m: FloatArray

    def __post_init__(self) -> None:
        edges = self.edges_m
        if edges.ndim != 1 or edges.size < 3:
            raise ValueError("need at least two cells")
        if edges[0] != 0.0 or np.any(np.diff(edges) <= 0.0):
            raise ValueError("edges must start at 0 and increase")

    @property
    def n_cells(self) -> int:
        return int(self.edges_m.size - 1)

    @property
    def thickness_m(self) -> FloatArray:
        return np.diff(self.edges_m)

    @property
    def centers_m(self) -> FloatArray:
        return 0.5 * (self.edges_m[:-1] + self.edges_m[1:])

    @property
    def depth_m(self) -> float:
        return float(self.edges_m[-1])


def geometric_grid(first_cell_m: float, growth: float, depth_m: float) -> DepthGrid:
    """Cells of thickness ``first_cell_m * growth**i`` down to at least ``depth_m``."""
    if first_cell_m <= 0.0 or depth_m <= first_cell_m:
        raise ValueError("need 0 < first_cell_m < depth_m")
    if growth < 1.0:
        raise ValueError("growth must be at least 1")
    edges = [0.0]
    thickness = first_cell_m
    while edges[-1] < depth_m:
        edges.append(edges[-1] + thickness)
        thickness *= growth
    return DepthGrid(np.asarray(edges, dtype=np.float64))


def uniform_grid(cell_m: float, depth_m: float) -> DepthGrid:
    """Cells of equal thickness (used by the analytic verification cases)."""
    count = round(depth_m / cell_m)
    if count < 2:
        raise ValueError("need at least two cells")
    return DepthGrid(np.linspace(0.0, count * cell_m, count + 1))
