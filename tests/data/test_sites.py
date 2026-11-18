"""A scenario runs at any terrain set whose site DEM is in the cache, by name alone."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.api import run_scenario
from lunarenv.data.dem import load_dem
from lunarenv.geometry.ephemeris import DE440
from lunarenv.sites import TERRAIN_SETS
from tests.helpers import PROJ

pytestmark = pytest.mark.data


@pytest.mark.parametrize("name", sorted(TERRAIN_SETS))
def test_run_at_the_centre_of_the_site_tile(name, store):
    terrain = TERRAIN_SETS[name]
    required = (*terrain.required_ids, *DE440.kernel_ids)
    missing = [i for i in required if not store.is_available(i)]
    if missing:
        pytest.skip(f"not in the data cache: {', '.join(missing)}")
    bounds = load_dem(store, terrain.near_id).grid.bounds
    x, y = 0.5 * (bounds.x_min + bounds.x_max), 0.5 * (bounds.y_min + bounds.y_max)
    longitude, latitude = (float(np.degrees(v)) for v in PROJ.inverse(x, y))

    result = run_scenario(
        {
            "name": f"centre-{name}",
            "site": {"terrain": name, "longitude_deg": longitude % 360.0, "latitude_deg": latitude},
            "period": {
                "start_utc": "2026-03-01T00:00:00",
                "end_utc": "2026-03-03T00:00:00",
                "step_minutes": 120,
            },
            "illumination": {"settings": {"n_azimuth": 360}},
            "thermal": {"enabled": False},
            "uncertainty": {"members": 0, "terrain_clones": 0, "regional_error_k": 0.0},
        },
        store=store,
    )
    site = result.metadata["site"]
    assert site["terrain_set"] == name
    np.testing.assert_allclose([site["x_m"], site["y_m"]], [x, y], atol=1e-3)
    dataset = result.dataset
    assert dataset.sizes["time"] == 25
    assert np.isfinite(dataset["horizon_elevation"].values).all()
    # the Sun stays within the lunar obliquity plus the colatitude of the site
    assert np.abs(dataset["sun_elevation"].values).max() < 90.0 + latitude + 1.6 + 1.0
    assert 0.0 <= float(dataset["sun_visible_fraction"].mean()) <= 1.0
    # with the wide DEM in the cache no ray is left without terrain that matters
    wide = site["wide_regional_dem"]
    if store.is_available(terrain.wide_id):
        assert wide["used"] == wide["needed"]
        assert not [w for w in result.metadata["validity"]["warnings"] if "too small" in w]
        assert (terrain.wide_id in {d["id"] for d in result.metadata["data"]}) == wide["used"]
    else:
        assert not wide["used"]


def test_the_reference_point_does_not_need_the_wide_dem(store):
    terrain = TERRAIN_SETS["site04"]
    missing = [i for i in (*terrain.required_ids, *DE440.kernel_ids) if not store.is_available(i)]
    if missing:
        pytest.skip(f"not in the data cache: {', '.join(missing)}")
    result = run_scenario(
        {
            "name": "rim",
            "site": {"x_m": -2822.5, "y_m": -9062.5},
            "period": {
                "start_utc": "2026-03-01T00:00:00",
                "end_utc": "2026-03-03T00:00:00",
                "step_minutes": 120,
            },
            "illumination": {"settings": {"n_azimuth": 360}},
            "thermal": {"enabled": False},
            "uncertainty": {"members": 0, "terrain_clones": 0, "regional_error_k": 0.0},
        },
        store=store,
    )
    wide = result.metadata["site"]["wide_regional_dem"]
    assert not wide["needed"] and not wide["used"]
    assert terrain.wide_id not in {d["id"] for d in result.metadata["data"]}
