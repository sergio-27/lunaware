"""SPICE chain: time conversion, frames and apparent positions."""

from __future__ import annotations

import numpy as np
import pytest
import spiceypy

from lunarenv.core.constants import (
    ASTRONOMICAL_UNIT_M,
    MOON_REFERENCE_RADIUS_M,
    SPEED_OF_LIGHT_M_S,
)
from lunarenv.data.dem import load_dem
from lunarenv.geometry.ephemeris import DE421, DE440, SpiceSession, spice_session
from lunarenv.geometry.projection import SouthPolarStereographic, elevation_azimuth

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


@pytest.fixture(scope="module")
def de440(store, require_module):
    require_module(*DE440.kernel_ids)
    with spice_session(store, DE440) as ephemeris:
        yield ephemeris


@pytest.fixture(scope="module")
def require_module(store):
    def _require(*entry_ids):
        missing = [e for e in entry_ids if not store.is_available(e)]
        if missing:
            pytest.skip(f"not in the data cache: {', '.join(missing)}")

    return _require


def test_utc_to_et_reference_value_and_round_trip(de440):
    # TDB - UTC at J2000 = 32 leap seconds + 32.184 s, minus the periodic term
    assert de440.utc_to_et("2000-01-01T12:00:00") == pytest.approx(64.183927285, abs=1e-6)
    stamps = np.array(["2025-06-30T23:59:59.500", "2026-11-12T16:00:00.000"])
    assert list(de440.et_to_utc(de440.utc_to_et(stamps))) == list(stamps)
    as_datetime = stamps.astype("datetime64[ms]")
    assert np.array_equal(de440.utc_to_et(as_datetime), de440.utc_to_et(stamps))


def test_body_radii_are_consistent_with_the_dem_datum(de440):
    assert de440.body_radius_m("SUN") == 695_700_000.0  # IAU 2015 nominal
    _, moon_km = spiceypy.bodvrd("MOON", "RADII", 3)
    assert np.all(np.asarray(moon_km) * 1000.0 == MOON_REFERENCE_RADIUS_M)


def test_sun_distance_and_subsolar_latitude_are_physical(de440, epochs_utc):
    sun = de440.position_from_center("SUN", de440.utc_to_et(epochs_utc))
    distance_au = np.linalg.norm(sun, axis=-1) / ASTRONOMICAL_UNIT_M
    assert 0.980 < distance_au.min() < 0.986 and 1.014 < distance_au.max() < 1.020
    subsolar_latitude = np.degrees(np.arcsin(sun[:, 2] / np.linalg.norm(sun, axis=-1)))
    # the lunar equator is inclined 1.54 deg to the ecliptic
    assert 1.50 < np.abs(subsolar_latitude).max() < 1.60


def test_sub_earth_point_stays_within_libration_limits(de440, epochs_utc):
    earth = de440.position_from_center("EARTH", de440.utc_to_et(epochs_utc))
    unit = earth / np.linalg.norm(earth, axis=-1, keepdims=True)
    longitude = np.degrees(np.arctan2(unit[:, 1], unit[:, 0]))
    latitude = np.degrees(np.arcsin(unit[:, 2]))
    assert np.abs(longitude).max() < 8.3 and np.abs(latitude).max() < 7.0
    # mean-Earth frame: the sub-Earth point averages to the origin
    assert abs(longitude.mean()) < 0.3 and abs(latitude.mean()) < 0.3
    distance_km = np.linalg.norm(earth, axis=-1) / 1e3
    assert 356_000 < distance_km.min() and distance_km.max() < 407_000


def test_sun_elevation_at_the_polar_site_is_grazing(de440, epochs_utc, site_position):
    """At 89.8 S the Sun never rises more than obliquity + colatitude above the horizon."""
    et = de440.utc_to_et(epochs_utc)
    direction = de440.position_from_center("SUN", et) - site_position
    elevation, azimuth = elevation_azimuth(direction, *PROJ.local_basis(*SITE_XY))
    _, latitude = PROJ.inverse(*SITE_XY)
    limit = np.radians(1.60) + (np.pi / 2 + latitude)
    assert np.abs(elevation).max() < limit
    assert np.abs(elevation).max() > limit - np.radians(0.3)
    # the Sun circles the horizon: every azimuth octant is visited
    assert len(np.unique((azimuth // (np.pi / 4)).astype(int))) == 8


def test_center_observer_approximation_error(de440, epochs_utc, site_position, record_property):
    """Subtracting the site vector from the Moon-centred apparent position is accurate.

    Reference: SPICE with the observer placed at the surface point. The dominant error
    is stellar aberration applied to the wrong vector, bounded by
    (lunar radius / minimum distance) * (v / c), with v the observer's barycentric speed
    (Earth orbital speed at perihelion plus lunar orbital speed, 31.4 km/s).
    """
    et = de440.utc_to_et(epochs_utc[::10])
    v_over_c = 31.4e3 / SPEED_OF_LIGHT_M_S
    for body, minimum_distance_m in (("SUN", 0.983 * ASTRONOMICAL_UNIT_M), ("EARTH", 3.56e8)):
        approximate = de440.position_from_center(body, et) - site_position
        exact = de440.position_from_site(body, et, site_position)
        error = float(angle_between(approximate, exact).max())
        bound = MOON_REFERENCE_RADIUS_M / minimum_distance_m * v_over_c
        record_property(f"center_approximation_max_rad_{body}", error)
        print(f"centre-observer approximation, {body}: max {error:.3e} rad (bound {bound:.3e})")
        assert error < 1.05 * bound


def test_sessions_are_exclusive(store, de440):
    with pytest.raises(RuntimeError, match="already active"):
        SpiceSession(store, DE421).__enter__()
