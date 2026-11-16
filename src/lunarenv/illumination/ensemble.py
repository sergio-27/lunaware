"""Propagation of terrain uncertainty to illumination with DEM error clones.

The PGDA site products come with 100 "clones": DEM realizations whose differences from
the nominal DEM have the statistics of its height error (orbit, pointing and
interpolation errors; Barker2021). Running the model on each clone turns terrain
uncertainty into a distribution of illumination results without any assumption of
linearity. This matters: raising the horizon takes one high pixel along a ray, lowering
it takes many, so height errors bias illumination downwards (Barker2021, section 5).

Scope: the clones cover the site tile only. Terrain beyond it is the same in every
member (Barker2021 argue the far field matters less: it is much farther away and
averaged over many more altimeter shots). The regional DEM's error is not available as
clones, only as a per-pixel rms height error map (Barker2023, section 3), so
:func:`regional_error_stacks` turns it into a bound instead: the regional terrain raised
and lowered everywhere by ``k`` times its error. The horizon along a ray is the maximum
of an increasing function of each sample's height, so any regional error with
``|dh| <= k * error`` pixel by pixel, whatever its correlation, gives a horizon between
the two, and every illumination statistic that is monotone in the horizon (visible
fraction, longest shadow, longest illumination) is bounded with it.

The observer follows the terrain: in each member it stands on that member's surface.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from lunarenv.data.dem import Dem, load_dem
from lunarenv.data.store import DataStore
from lunarenv.geometry.projection import FloatArray
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_statistics
from lunarenv.illumination.settings import IlluminationSettings
from lunarenv.illumination.terrain import TerrainStack, build_stack

QUANTILES = (0.01, 0.05, 0.50, 0.95, 0.99)


@dataclass(frozen=True, slots=True)
class IlluminationEnsemble:
    """Per-member statistics; arrays have shape ``(n_members, n_observers)``."""

    member_ids: tuple[str, ...]
    mean_visible_fraction: FloatArray
    mean_flux_fraction: FloatArray
    longest_shadow_s: FloatArray
    longest_illumination_s: FloatArray

    def quantiles(self, field: str) -> FloatArray:
        """Quantiles :data:`QUANTILES` of a field over the members, per observer."""
        values = getattr(self, field)
        return np.asarray(np.quantile(values, QUANTILES, axis=0))


def cached_clones(store: DataStore, clone_ids: Sequence[str]) -> Iterator[tuple[str, Dem]]:
    """Yield ``(id, DEM)`` for clones in the data cache, one at a time (41 MB each)."""
    for clone_id in clone_ids:
        yield clone_id, load_dem(store, clone_id)


def terrain_ensemble(
    stack: TerrainStack,
    members: Iterable[tuple[str, Dem]],
    x_m: ArrayLike,
    y_m: ArrayLike,
    et_s: FloatArray,
    sun_position_m: FloatArray,
    sun_radius_m: float,
    *,
    height_above_surface_m: float = 0.0,
    settings: IlluminationSettings | None = None,
) -> IlluminationEnsemble:
    """Illumination statistics of each observer on each near-field terrain realization.

    ``members`` supplies alternative DEMs for the nearest level of ``stack``; they must
    be on the same grid as the nominal one.
    """
    settings = settings or IlluminationSettings()
    ids: list[str] = []
    columns: list[list[FloatArray]] = [[], [], [], []]
    for member_id, dem in members:
        horizons = compute_horizons(
            stack.with_near(dem),
            x_m,
            y_m,
            height_above_surface_m=height_above_surface_m,
            n_azimuth=settings.n_azimuth,
            diagnostics=False,
        )
        stats = illumination_statistics(
            horizons,
            et_s,
            sun_position_m,
            sun_radius_m,
            n_strips=settings.n_strips,
            limb_u=settings.limb_darkening,
            mode=settings.disc_mode_code,
        )
        ids.append(member_id)
        for column, values in zip(
            columns,
            (
                stats.mean_visible_fraction,
                stats.mean_flux_fraction,
                stats.longest_shadow_s,
                stats.longest_illumination_s,
            ),
            strict=True,
        ):
            column.append(values)
    if not ids:
        raise ValueError("need at least one ensemble member")
    stacked = [np.stack(column) for column in columns]
    return IlluminationEnsemble(tuple(ids), *stacked)


def shifted(dem: Dem, offset_m: ArrayLike, factor: float) -> Dem:
    """``dem`` with ``factor * offset_m`` added to its heights (same grid and provenance)."""
    offset = np.asarray(offset_m, dtype=np.float64)
    if offset.shape != dem.heights_m.shape:
        raise ValueError("offset must have the DEM's shape")
    heights = dem.heights_m.astype(np.float64) + factor * offset
    return Dem(heights.astype(np.float32), dem.grid, dem.source_id, dem.source_sha256)


def regional_error_stacks(
    near: Dem,
    regional: Dem,
    error: Dem,
    k: float,
    settings: IlluminationSettings | None = None,
) -> tuple[TerrainStack, TerrainStack]:
    """Terrain stacks with the regional DEM lowered and raised by ``k`` times its error.

    ``error`` is the rms height error map of ``regional`` (same grid). The near level and
    the observer are unchanged; the far level is block-averaged from the shifted DEM.
    Returns ``(lowered, raised)``: horizons on them bound the horizon of any regional
    terrain within ``k`` errors of the nominal one (see the module docstring).
    """
    if error.grid != regional.grid:
        raise ValueError("the error map must be on the regional DEM's grid")
    if not k > 0.0:
        raise ValueError("k must be positive")
    if float(error.heights_m.min()) < 0.0:
        raise ValueError("height errors must be non-negative")
    settings = settings or IlluminationSettings()

    def stack(sign: float) -> TerrainStack:
        return build_stack(
            near,
            shifted(regional, error.heights_m, sign * k),
            near_radius_m=settings.near_radius_m,
            mid_radius_m=settings.mid_radius_m,
            step_fraction=settings.step_fraction,
            far_coarsening=settings.far_coarsening,
        )

    return stack(-1.0), stack(1.0)
