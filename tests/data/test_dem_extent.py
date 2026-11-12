"""PLAN.md decision D4 applied to the test tile: is the 80S regional DEM large enough?"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio

from lunarenv.data.dem import load_dem
from lunarenv.geometry.visibility import dem_extent_is_sufficient_for_sun

pytestmark = pytest.mark.data


def test_products_are_pixel_registered(store, require):
    """GMT ``node_offset = 1`` means pixel registration, as the PGDA README states."""
    require("pgda78.site04.surf")
    with rasterio.open(store.path("pgda78.site04.surf")) as dataset:
        assert dataset.tags()["NC_GLOBAL#node_offset"] == "1"
        assert dataset.tags()["AREA_OR_POINT"] == "Area"


def test_regional_dem_covers_every_sun_occluder_for_the_site04_tile(store, require):
    """For which observers in the tile is terrain beyond the 80S product provably irrelevant?

    Conservative on every count: shortest distance to the product edge, terrain at the
    height cap assumed right at the edge, lowest possible solar limb.
    """
    require("pgda78.site04.surf")
    site = load_dem(store, "pgda78.site04.surf")
    half_width = 304_000.0  # LDEM_80S_80MPP_ADJ extends to +/-304 km (asserted elsewhere)
    xx, yy = np.meshgrid(site.grid.x_centers()[::8], site.grid.y_centers()[::8])
    heights = site.heights_m[::8, ::8].astype(np.float64)
    # shortest map distance to the product edge, converted to ground distance with the
    # largest scale factor on the way (k <= 1.0077 at the edge): conservative
    edge_distance = (half_width - np.maximum(np.abs(xx), np.abs(yy))) / 1.0077
    _, latitude = site.grid.projection.inverse(xx, yy)
    sufficient = dem_extent_is_sufficient_for_sun(edge_distance, heights, latitude)
    fraction = float(sufficient.mean())
    highest_uncovered = float(heights[~sufficient].max()) if fraction < 1.0 else float("-inf")
    print(
        f"80S product provably sufficient for {100 * fraction:.2f}% of the tile; "
        f"guaranteed for all observers above {highest_uncovered:.0f} m"
    )
    # every observer on or above the datum (rim, flanks, plains) is covered
    assert np.all(sufficient[heights >= 0.0])
    assert fraction > 0.7
