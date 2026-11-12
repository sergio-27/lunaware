"""Earth-occultation factor against a published lunar eclipse.

Total lunar eclipse of 2026 March 3 (F. Espenak, NASA GSFC, eclipse.gsfc.nasa.gov,
figure LE2026Mar03T): penumbral contacts P1 08:44:22 and P4 14:22:59 UT, umbral
contacts U1 09:50:00 and U4 13:17:10, totality U2 11:04:26 to U3 12:02:45, greatest
eclipse 11:33:37 UT. Those contacts are for the Moon's limb entering Earth's shadow as
enlarged by Danjon's rule; a single site sees its own contacts inside that window.
"""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.geometry.ephemeris import DE440, spice_session
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.illumination.eclipse import earth_occultation_factor

pytestmark = pytest.mark.data


def test_total_lunar_eclipse_of_2026_march_3(store, require):
    require(*DE440.kernel_ids)
    site = SouthPolarStereographic().to_cartesian(-2822.5, -9062.5, 1733.0)
    utc = np.arange(
        np.datetime64("2026-03-03T06:00"), np.datetime64("2026-03-03T17:00"), np.timedelta64(1, "m")
    )
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(utc)
        factor = earth_occultation_factor(
            ephemeris.position_from_center("SUN", et),
            ephemeris.position_from_center("EARTH", et),
            site,
            ephemeris.body_radius_m("SUN"),
            ephemeris.body_radius_m("EARTH"),
        )
    minute = lambda text: int(np.flatnonzero(utc == np.datetime64(f"2026-03-03T{text}"))[0])  # noqa: E731
    eclipsed = np.flatnonzero(factor < 1.0)
    dark = np.flatnonzero(factor == 0.0)
    print(
        f"site sees the Sun dimmed {utc[eclipsed[0]]} to {utc[eclipsed[-1]]}, "
        f"fully hidden {utc[dark[0]]} to {utc[dark[-1]]}"
    )
    assert factor[minute("11:34")] == 0.0  # greatest eclipse
    # nothing outside the published penumbral window (3 min margin: no Danjon enlargement here)
    assert eclipsed[0] >= minute("08:44") - 3 and eclipsed[-1] <= minute("14:23") + 3
    # the site's total phase lies within the published umbral phase and lasts a similar time
    assert dark[0] >= minute("09:50") and dark[-1] <= minute("13:17")
    assert 45 <= dark.size <= 215
    assert np.all(np.diff(factor[: dark[0]]) <= 1e-12) and np.all(
        np.diff(factor[dark[-1] :]) >= -1e-12
    )
