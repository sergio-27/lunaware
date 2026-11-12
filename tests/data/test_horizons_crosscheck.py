"""ILL-X1: apparent Sun and Earth directions at the site against JPL Horizons.

The reference tables are verbatim responses of the JPL Horizons API (retrieved
2026-11-12), stored in ``lunarenv/validation/reference``; the comparison itself lives in
``lunarenv.validation.cases_illumination`` and is exercised here.

Tolerance, fixed before the first comparison: 2e-6 rad (0.4 arcsec). Expected
differences are the DE440/DE441 and frame-realization differences (D6 measured 3e-7 rad
between DE440 and DE421) and the 1e-6 deg rounding of the site coordinates sent to
Horizons (3 cm). The Sun's angular radius is 4.65e-3 rad.
"""

from __future__ import annotations

import pytest

from lunarenv.geometry.ephemeris import DE440
from lunarenv.validation.cases_illumination import horizons_separation, read_horizons_table

pytestmark = pytest.mark.data
TOLERANCE_RAD = 2e-6


@pytest.mark.parametrize(
    ("body", "fixture"),
    [("SUN", "horizons_10_site04_rim.txt"), ("EARTH", "horizons_399_site04_rim.txt")],
)
def test_ill_x1_direction_matches_jpl_horizons(store, require, body, fixture):
    require(*DE440.kernel_ids)
    assert read_horizons_table(fixture)[0].size > 250
    exact, approximate = horizons_separation(store, body, fixture)
    print(
        f"{body} vs Horizons: surface observer max {exact.max():.2e} rad; "
        f"centre approximation max {approximate.max():.2e} rad"
    )
    assert exact.max() < TOLERANCE_RAD
    # the fast path used for maps adds at most its documented aberration term
    assert approximate.max() < TOLERANCE_RAD + 5.2e-7
