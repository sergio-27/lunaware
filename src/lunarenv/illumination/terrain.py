"""Nested terrain models for horizon tracing.

A horizon ray needs fine terrain close to the observer and can use coarser terrain far
away, where a pixel subtends a small angle. A :class:`TerrainStack` is an ordered list
of DEM levels, each assigned a range of ground distances and a sampling step. The
ranges and steps are parameters whose values are chosen by the convergence study in
``docs/convergence/illumination.md``, not constants of the method.

Coarser levels are produced by block averaging, which is what a coarser LOLA product
is (mean height per bin). Averaging lowers crests, so a coarse level slightly
under-estimates the horizon; the study quantifies this.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from lunarenv.data.dem import Dem, GridGeometry


@dataclass(frozen=True, slots=True)
class TerrainLevel:
    """A DEM used for ground distances in ``[start_m, end_m)`` with a fixed step."""

    dem: Dem
    start_m: float
    end_m: float
    step_m: float

    def __post_init__(self) -> None:
        if not (0.0 <= self.start_m < self.end_m and self.step_m > 0.0):
            raise ValueError("level needs 0 <= start < end and a positive step")
        if min(self.dem.grid.shape) < 2:
            raise ValueError("a terrain level needs at least 2 x 2 pixels")


@dataclass(frozen=True, slots=True)
class TerrainStack:
    """Terrain levels ordered from near to far, covering contiguous distance ranges."""

    levels: tuple[TerrainLevel, ...]

    def __post_init__(self) -> None:
        if not self.levels:
            raise ValueError("a terrain stack needs at least one level")
        if self.levels[0].start_m != 0.0:
            raise ValueError("the first level must start at distance 0")
        for inner, outer in zip(self.levels, self.levels[1:], strict=False):
            if inner.end_m != outer.start_m:
                raise ValueError("levels must cover contiguous distance ranges")
            if inner.dem.grid.projection != outer.dem.grid.projection:
                raise ValueError("levels must share one projection")

    @property
    def radius_m(self) -> float:
        return self.levels[0].dem.grid.projection.radius_m

    @property
    def max_height_m(self) -> float:
        """Highest terrain in any level: the cap used by the ray stopping rule."""
        return max(float(level.dem.heights_m.max()) for level in self.levels)

    def sample_plan(self) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
        """Ground distances of the samples along every ray, and the level of each.

        Within a level, samples sit at the centres of consecutive steps, so the first
        sample is half a step from the observer and no sample falls on a level boundary.
        """
        distances: list[NDArray[np.float64]] = []
        levels: list[NDArray[np.int64]] = []
        for index, level in enumerate(self.levels):
            count = max(1, round((level.end_m - level.start_m) / level.step_m))
            step = (level.end_m - level.start_m) / count
            offsets = (np.arange(count, dtype=np.float64) + 0.5) * step
            distances.append(np.asarray(level.start_m + offsets, dtype=np.float64))
            levels.append(np.full(count, index, dtype=np.int64))
        return np.concatenate(distances), np.concatenate(levels)

    def with_near(self, dem: Dem) -> TerrainStack:
        """Same stack with the nearest level's terrain replaced (error clones)."""
        first = self.levels[0]
        if dem.grid != first.dem.grid:
            raise ValueError("replacement DEM must be on the same grid as the near level")
        return TerrainStack(
            (TerrainLevel(dem, first.start_m, first.end_m, first.step_m), *self.levels[1:])
        )


def coarsen(dem: Dem, factor: int) -> Dem:
    """Block-average a DEM by an integer factor.

    Rows and columns that do not fill a whole block at the south and east edges are
    dropped, so the result stays aligned with the top-left corner.
    """
    if factor < 1:
        raise ValueError("factor must be a positive integer")
    if factor == 1:
        return dem
    grid = dem.grid
    rows, cols = grid.height // factor, grid.width // factor
    if rows < 2 or cols < 2:
        raise ValueError("DEM too small for this coarsening factor")
    block = dem.heights_m[: rows * factor, : cols * factor].astype(np.float64)
    mean = block.reshape(rows, factor, cols, factor).mean(axis=(1, 3))
    coarse_grid = GridGeometry(
        x_min=grid.x_min,
        y_max=grid.y_max,
        pixel_size_m=grid.pixel_size_m * factor,
        width=cols,
        height=rows,
        projection=grid.projection,
    )
    return Dem(mean.astype(np.float32), coarse_grid, dem.source_id, dem.source_sha256)


def build_stack(
    near: Dem,
    regional: Dem,
    *,
    near_radius_m: float,
    mid_radius_m: float,
    step_fraction: float = 0.5,
    far_coarsening: int = 3,
    wide: Dem | None = None,
) -> TerrainStack:
    """Three-level stack: site DEM, regional DEM, coarsened regional DEM.

    * near: ``near`` at its own resolution out to ``near_radius_m``;
    * mid: ``regional`` at its own resolution out to ``mid_radius_m``;
    * far: ``regional`` block-averaged by ``far_coarsening`` to the edge of the data.

    Every level is sampled every ``step_fraction`` pixels. Where a ray leaves the near
    DEM before ``near_radius_m`` the tracer falls through to the next level, so the
    near DEM does not have to be centred on the observers.

    With ``wide``, a coarser DEM that reaches beyond ``regional``, a fourth level is
    added. In the same way a ray that leaves ``regional`` falls through to it, so the
    horizon no longer ends at the edge of the regional product. Inside ``regional`` the
    terrain sampled is the same with and without it.
    """
    if not 0.0 < near_radius_m < mid_radius_m:
        raise ValueError("need 0 < near_radius_m < mid_radius_m")
    if not 0.0 < step_fraction <= 2.0:
        raise ValueError("step_fraction must be in (0, 2]")
    far = coarsen(regional, far_coarsening)
    bounds = regional.grid.bounds
    # no observer-to-data distance can exceed the diagonal of the regional product
    reach = float(np.hypot(bounds.x_max - bounds.x_min, bounds.y_max - bounds.y_min))
    if mid_radius_m >= reach:
        raise ValueError("mid_radius_m lies beyond the regional DEM")
    levels = [
        TerrainLevel(near, 0.0, near_radius_m, step_fraction * near.grid.pixel_size_m),
        TerrainLevel(
            regional, near_radius_m, mid_radius_m, step_fraction * regional.grid.pixel_size_m
        ),
        TerrainLevel(far, mid_radius_m, reach, step_fraction * far.grid.pixel_size_m),
    ]
    if wide is not None:
        if not wide.grid.bounds.contains(bounds):
            raise ValueError("the wide DEM must cover the regional DEM")
        wide_bounds = wide.grid.bounds
        wide_reach = float(
            np.hypot(wide_bounds.x_max - wide_bounds.x_min, wide_bounds.y_max - wide_bounds.y_min)
        )
        if wide_reach <= reach:
            raise ValueError("the wide DEM must reach beyond the regional DEM")
        levels.append(TerrainLevel(wide, reach, wide_reach, step_fraction * wide.grid.pixel_size_m))
    return TerrainStack(tuple(levels))


def regional_edge_distance_m(regional: GridGeometry, x_m: float, y_m: float) -> float:
    """Lower bound of the ground distance from a map point to the edge of a DEM.

    The shortest map distance to the edge, divided by the largest scale factor of the
    projection inside the DEM (map distances are longer than ground distances away from
    the pole). Negative if the point is outside.
    """
    bounds = regional.bounds
    map_distance = min(
        x_m - bounds.x_min, bounds.x_max - x_m, y_m - bounds.y_min, bounds.y_max - y_m
    )
    corner_x = max(abs(bounds.x_min), abs(bounds.x_max))
    corner_y = max(abs(bounds.y_min), abs(bounds.y_max))
    return float(map_distance / regional.projection.scale_factor(corner_x, corner_y))
