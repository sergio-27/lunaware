"""Astronomical bounds used by the DEM-extent criterion, verified against DE440."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.geometry.ephemeris import DE440, spice_session
from lunarenv.geometry.visibility import (
    MAX_SOLAR_ANGULAR_RADIUS_RAD,
    MAX_SOLAR_PARALLAX_RAD,
    MAX_SUBSOLAR_LATITUDE_RAD,
)

pytestmark = pytest.mark.data


def test_solar_bounds_hold_over_a_nodal_cycle(store, require):
    """2024-2044 every 6 h: more than one 18.6-year cycle, all seasons and perihelia."""
    require(*DE440.kernel_ids)
    utc = np.arange(
        np.datetime64("2024-01-01T00:00"), np.datetime64("2044-01-01T00:00"), np.timedelta64(6, "h")
    )
    with spice_session(store, DE440) as ephemeris:
        sun = ephemeris.position_from_center("SUN", ephemeris.utc_to_et(utc))
        solar_radius = ephemeris.body_radius_m("SUN")
        moon_radius = ephemeris.body_radius_m("MOON")
    distance = np.linalg.norm(sun, axis=-1)
    subsolar_latitude = np.abs(np.arcsin(sun[:, 2] / distance)).max()
    angular_radius = np.arcsin(solar_radius / distance).max()
    parallax = np.arcsin(moon_radius / distance).max()
    print(
        f"max |subsolar latitude| {np.degrees(subsolar_latitude):.4f} deg, "
        f"max solar radius {np.degrees(angular_radius):.4f} deg, max parallax {parallax:.3e} rad"
    )
    # bounds must hold, and must not be so loose that they hide a mistake
    assert 0.95 * MAX_SUBSOLAR_LATITUDE_RAD < subsolar_latitude < MAX_SUBSOLAR_LATITUDE_RAD
    assert 0.99 * MAX_SOLAR_ANGULAR_RADIUS_RAD < angular_radius < MAX_SOLAR_ANGULAR_RADIUS_RAD
    assert 0.9 * MAX_SOLAR_PARALLAX_RAD < parallax < MAX_SOLAR_PARALLAX_RAD
