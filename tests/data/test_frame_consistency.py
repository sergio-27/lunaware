"""Frame-consistency check between the DE440 chain used by lunarenv and DE421."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.data.dem import load_dem
from lunarenv.geometry.ephemeris import DE421, DE440, spice_session
from lunarenv.geometry.projection import SouthPolarStereographic

pytestmark = pytest.mark.data

PROJ = SouthPolarStereographic()
# Test point on the Shackleton rim inside the Site 4 tile (map metres); the height is
# read from the DEM where available and is irrelevant at the tolerances used otherwise.
SITE_XY = (-1000.0, -7000.0)


def angle_between(a, b):
    cross = np.linalg.norm(np.cross(a, b), axis=-1)
    return np.arctan2(cross, np.sum(a * b, axis=-1))


@pytest.fixture(scope="module")
def epochs_utc():
    # 2024-2034, every 73 h: samples all lunar phases, seasons and half a nodal cycle
    return np.arange(
        np.datetime64("2024-01-01T00:00"),
        np.datetime64("2034-01-01T00:00"),
        np.timedelta64(73, "h"),
    )


@pytest.fixture(scope="module")
def site_position(store):
    height = 0.0
    if store.is_available("pgda78.site04.surf"):
        height = float(load_dem(store, "pgda78.site04.surf").sample(*SITE_XY))
    return PROJ.to_cartesian(*SITE_XY, height)


def test_d6_de440_frame_agrees_with_de421(
    store, require, epochs_utc, site_position, record_property
):
    """PLAN.md decision D6.

    The LOLA products live in the DE421 mean-Earth frame. lunarenv uses DE440 with the
    frame MOON_ME_DE440_ME421. The direction of the Sun and the Earth seen from the site
    must agree with a pure DE421 computation far below the scales that matter:
    the Sun's angular radius is 4.65e-3 rad and one 5 m pixel at 5 km subtends 1e-3 rad.
    Tolerance fixed before the first run: 1e-5 rad (2 arcsec).
    """
    require(*DE421.kernel_ids, *DE440.kernel_ids)
    in_range = epochs_utc[epochs_utc < np.datetime64("2034-01-01")]
    directions = {}
    for kernel_set in (DE440, DE421):
        with spice_session(store, kernel_set) as ephemeris:
            et = ephemeris.utc_to_et(in_range)
            directions[kernel_set.name] = {
                body: ephemeris.position_from_center(body, et) - site_position
                for body in ("SUN", "EARTH")
            }
    for body in ("SUN", "EARTH"):
        difference = angle_between(directions[DE440.name][body], directions[DE421.name][body])
        record_property(f"d6_max_rad_{body}", float(difference.max()))
        print(
            f"D6 {body}: max {difference.max():.3e} rad, "
            f"rms {np.sqrt(np.mean(difference**2)):.3e} rad"
        )
        assert difference.max() < 1e-5
