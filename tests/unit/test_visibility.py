"""Closed-form line-of-sight relations on a sphere."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.geometry.projection import SouthPolarStereographic, elevation_azimuth
from lunarenv.geometry.visibility import (
    apparent_elevation,
    dem_extent_is_sufficient_for_sun,
    horizon_dip,
    lowest_solar_limb_elevation,
    occluding_height,
)

PROJ = SouthPolarStereographic()
distance = st.floats(10.0, 6.0e5)
height = st.floats(-9.0e3, 1.1e4)


@given(distance, height, height)
def test_apparent_elevation_matches_3d_vectors(d, h_target, h_observer):
    """Same quantity computed from Cartesian positions through the projection code."""
    observer_xy = (2000.0, -7000.0)
    lon, lat = PROJ.inverse(*observer_xy)
    # target at ground distance d along the meridian towards the equator
    target_xy = PROJ.forward(lon, lat + d / R)
    line = PROJ.to_cartesian(*target_xy, h_target) - PROJ.to_cartesian(*observer_xy, h_observer)
    elevation, _ = elevation_azimuth(line, *PROJ.local_basis(*observer_xy))
    assert abs(elevation - apparent_elevation(d, h_target, h_observer)) < 1e-9


@given(distance, height, st.floats(-0.1, 0.5))
def test_occluding_height_inverts_apparent_elevation(d, h_observer, elevation):
    needed = occluding_height(d, h_observer, elevation)
    assert np.isfinite(needed)
    assert abs(apparent_elevation(d, needed, h_observer) - elevation) < 1e-9


def test_occluding_height_is_infinite_beyond_reach():
    assert np.isinf(occluding_height(0.5 * np.pi * R, 0.0, 0.01))


@given(st.floats(0.0, 2.0e4))
def test_horizon_dip_is_the_grazing_ray(h):
    """At the dip angle, the required target height is minimal and equals zero."""
    dip = float(horizon_dip(h))
    tangent_distance = R * np.arccos(R / (R + h))
    assert abs(occluding_height(tangent_distance, h, dip)) < 1e-6
    for d in (0.5 * tangent_distance + 1.0, 1.5 * tangent_distance + 1.0):
        assert occluding_height(d, h, dip) >= -1e-6


def test_horizon_dip_rejects_observer_below_sphere():
    with pytest.raises(ValueError):
        horizon_dip(-1.0)


def test_lowest_solar_limb_at_the_pole_and_off_pole():
    at_pole = np.degrees(lowest_solar_limb_elevation(-np.pi / 2))
    assert at_pole == pytest.approx(-(1.60 + 0.2720) - np.degrees(1.2e-5), abs=1e-9)
    off_pole = np.degrees(lowest_solar_limb_elevation(np.radians(-89.0)))
    assert off_pole == pytest.approx(at_pole - 1.0, abs=1e-9)


def test_dem_extent_criterion_is_monotonic_in_distance_and_height():
    latitude = np.radians(-89.8)
    distances = np.array([50e3, 150e3, 289e3, 600e3])
    sufficient = dem_extent_is_sufficient_for_sun(distances, 1300.0, latitude)
    assert list(sufficient) == sorted(sufficient)  # once sufficient, stays sufficient
    assert sufficient[-1] and not sufficient[0]
    # a lower observer looks up at distant terrain: the requirement gets stricter
    low = dem_extent_is_sufficient_for_sun(distances, -2800.0, latitude)
    assert np.all(low <= sufficient)
    # closer than the grazing point the edge value is not a bound: never reported sufficient
    assert not dem_extent_is_sufficient_for_sun(30e3, 1300.0, latitude, terrain_height_cap_m=-1e9)
