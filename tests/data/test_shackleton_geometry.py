"""Independent check of georeferencing against published Shackleton crater geometry.

Zuber2012 (Nature 486, 378), Table 1, derived from LOLA altimetry:
rim centre at latitude -89.655 deg, longitude 129.174 deg E; lunar radius at the floor
centre 1734.63 km; mean rim diameter 21 km; mean rim height above datum 1.3 km; mean
depth rim to floor 4.1 +/- 0.05 km; floor topography range 0.210 km.

These numbers were not used to build the reader. An axis swap, a sign error, a wrong
central meridian or a wrong height unit would each move the crater by kilometres or
change heights by orders of magnitude. Tolerances are set by the precision of the
published values and by the natural irregularity of the rim, not tuned.
"""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M
from lunarenv.data.dem import Bounds, load_dem
from lunarenv.geometry.projection import SouthPolarStereographic

pytestmark = pytest.mark.data

PROJ = SouthPolarStereographic()
CENTER_LAT_DEG, CENTER_LON_DEG = -89.655, 129.174
FLOOR_RADIUS_M = 1_734_630.0
RIM_DIAMETER_M = 21_000.0
RIM_HEIGHT_M = 1_300.0
DEPTH_M = 4_100.0


@pytest.fixture(scope="module")
def crater(store):
    if not store.is_available("pgda90.ldem_80s_80mpp"):
        pytest.skip("not in the data cache: pgda90.ldem_80s_80mpp")
    cx, cy = (
        float(v) for v in PROJ.forward(np.radians(CENTER_LON_DEG), np.radians(CENTER_LAT_DEG))
    )
    dem = load_dem(store, "pgda90.ldem_80s_80mpp", Bounds.around(cx, cy, 20_000.0))
    return cx, cy, dem


def radial_profile(cx, cy, dem, radii):
    azimuth = np.radians(np.arange(0.0, 360.0, 1.0))
    x = cx + radii[:, None] * np.sin(azimuth)
    y = cy + radii[:, None] * np.cos(azimuth)
    return dem.sample(x, y)  # (radius, azimuth)


def test_floor_height_at_published_centre(crater):
    cx, cy, dem = crater
    floor = float(np.median(radial_profile(cx, cy, dem, np.array([0.0, 500.0, 1000.0]))))
    expected = FLOOR_RADIUS_M - MOON_REFERENCE_RADIUS_M  # -2770 m
    assert floor == pytest.approx(expected, abs=210.0)  # published floor relief


def test_rim_diameter_height_and_depth(crater):
    cx, cy, dem = crater
    radii = np.arange(8_000.0, 13_000.0, 80.0)
    heights = radial_profile(cx, cy, dem, radii)
    rim_radius = float(radii[np.argmax(heights.mean(axis=1))])
    rim_height = float(heights.max(axis=0).mean())  # crest height per azimuth, averaged
    floor = float(np.median(radial_profile(cx, cy, dem, np.array([0.0, 500.0, 1000.0]))))
    print(
        f"rim radius {rim_radius:.0f} m, mean rim-crest height {rim_height:.0f} m, "
        f"depth {rim_height - floor:.0f} m"
    )
    assert 2.0 * rim_radius == pytest.approx(RIM_DIAMETER_M, abs=1_000.0)  # "21 km"
    assert rim_height == pytest.approx(RIM_HEIGHT_M, abs=150.0)  # "1.3 km"
    assert rim_height - floor == pytest.approx(DEPTH_M, abs=250.0)


def test_crater_is_centred_where_published(crater):
    """The rim-crest circle fitted to the DEM is centred on the published coordinates."""
    cx, cy, dem = crater
    radii = np.arange(8_000.0, 13_000.0, 80.0)
    heights = radial_profile(cx, cy, dem, radii)
    crest_radius = radii[np.argmax(heights, axis=0)]  # per azimuth
    azimuth = np.radians(np.arange(0.0, 360.0, 1.0))
    px, py = crest_radius * np.sin(azimuth), crest_radius * np.cos(azimuth)
    # algebraic least-squares circle: x^2 + y^2 = 2 a x + 2 b y + c
    design = np.column_stack([2 * px, 2 * py, np.ones_like(px)])
    (a, b, _), *_ = np.linalg.lstsq(design, px**2 + py**2, rcond=None)
    offset = float(np.hypot(a, b))
    print(f"fitted rim-circle centre offset from published centre: {offset:.0f} m")
    # 0.001 deg of latitude is 30 m; the rim is irregular at the few-hundred-metre level
    assert offset < 500.0
