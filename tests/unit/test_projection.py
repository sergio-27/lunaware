"""Projection and local-frame geometry: identities, PROJ agreement, differential checks."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from rasterio.crs import CRS
from rasterio.warp import transform as proj_transform

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.geometry.projection import SouthPolarStereographic, elevation_azimuth

PROJ = SouthPolarStereographic()

# Region of use: within 30 degrees of the south pole (the 60S products are the widest).
MAP_LIMIT_M = 9.0e5
coordinate = st.floats(-MAP_LIMIT_M, MAP_LIMIT_M, allow_nan=False)
height = st.floats(-9.0e3, 1.1e4, allow_nan=False)


def spherical_to_cartesian(lon, lat, r):
    return np.stack(
        [r * np.cos(lat) * np.cos(lon), r * np.cos(lat) * np.sin(lon), r * np.sin(lat)], axis=-1
    )


def test_pole_maps_to_origin_and_axes_follow_convention():
    x, y = PROJ.forward(0.0, -np.pi / 2)
    assert abs(x) < 1e-6 and abs(y) < 1e-6
    # longitude 0 is +y, longitude 90E is +x
    x0, y0 = PROJ.forward(0.0, np.radians(-89.0))
    x90, y90 = PROJ.forward(np.pi / 2, np.radians(-89.0))
    assert y0 > 0 and abs(x0) < 1e-6
    assert x90 > 0 and abs(y90) < 1e-6


def test_forward_matches_proj():
    rng = np.random.default_rng(1)
    lon = rng.uniform(-180.0, 180.0, 500)
    lat = rng.uniform(-90.0, -60.0, 500)
    geographic = CRS.from_proj4(f"+proj=longlat +R={R} +no_defs")
    stereo = CRS.from_proj4(f"+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +R={R} +units=m +no_defs")
    px, py = proj_transform(geographic, stereo, lon.tolist(), lat.tolist())
    x, y = PROJ.forward(np.radians(lon), np.radians(lat))
    assert np.max(np.abs(x - np.array(px))) < 1e-6
    assert np.max(np.abs(y - np.array(py))) < 1e-6


@given(coordinate, coordinate)
def test_inverse_then_forward_is_identity(x, y):
    lon, lat = PROJ.inverse(x, y)
    x2, y2 = PROJ.forward(lon, lat)
    assert abs(x2 - x) < 1e-6 and abs(y2 - y) < 1e-6


@given(coordinate, coordinate, height)
def test_cartesian_matches_spherical_formula(x, y, h):
    lon, lat = PROJ.inverse(x, y)
    expected = spherical_to_cartesian(lon, lat, R + h)
    assert np.max(np.abs(PROJ.to_cartesian(x, y, h) - expected)) < 1e-6


@given(coordinate, coordinate, height)
def test_cartesian_round_trip(x, y, h):
    x2, y2, h2 = PROJ.from_cartesian(PROJ.to_cartesian(x, y, h))
    assert abs(x2 - x) < 1e-6 and abs(y2 - y) < 1e-6 and abs(h2 - h) < 1e-6


@given(coordinate, coordinate)
def test_local_basis_is_orthonormal_and_right_handed(x, y):
    e_x, e_y, up = PROJ.local_basis(x, y)
    basis = np.stack([e_x, e_y, up])
    assert np.max(np.abs(basis @ basis.T - np.eye(3))) < 1e-12
    assert np.max(np.abs(np.cross(e_x, e_y) - up)) < 1e-12


def test_local_basis_at_pole():
    e_x, e_y, up = PROJ.local_basis(0.0, 0.0)
    assert (
        np.allclose(e_x, [0, 1, 0]) and np.allclose(e_y, [1, 0, 0]) and np.allclose(up, [0, 0, -1])
    )


@given(coordinate, coordinate)
@settings(max_examples=200)
def test_basis_vectors_are_tangent_to_map_axes(x, y):
    """e_x, e_y equal the normalized derivative of the surface point along map x, y."""
    step = 10.0
    e_x, e_y, up = PROJ.local_basis(x, y)
    dx = PROJ.to_cartesian(x + step, y) - PROJ.to_cartesian(x - step, y)
    dy = PROJ.to_cartesian(x, y + step) - PROJ.to_cartesian(x, y - step)
    assert np.max(np.abs(dx / np.linalg.norm(dx) - e_x)) < 1e-9
    assert np.max(np.abs(dy / np.linalg.norm(dy) - e_y)) < 1e-9
    assert np.max(np.abs(PROJ.to_cartesian(x, y) / R - up)) < 1e-12


@given(coordinate, coordinate)
@settings(max_examples=200)
def test_scale_factor_is_map_over_ground_distance(x, y):
    step = 10.0
    ground = np.linalg.norm(PROJ.to_cartesian(x + step, y) - PROJ.to_cartesian(x - step, y))
    assert abs(2 * step / ground - PROJ.scale_factor(x, y)) < 1e-9


@pytest.mark.parametrize("lon_deg", [-170.0, -90.0, -30.0, 0.0, 45.0, 90.0, 179.0])
def test_grid_convergence_links_true_and_grid_azimuth(lon_deg):
    """Local north seen in the grid frame has grid azimuth equal to the longitude."""
    lon, lat = np.radians(lon_deg), np.radians(-88.0)
    x, y = PROJ.forward(lon, lat)
    north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    basis = PROJ.local_basis(x, y)
    for true_azimuth, vector in ((0.0, north), (np.pi / 2, east)):
        elevation, grid_azimuth = elevation_azimuth(vector, *basis)
        expected = np.mod(true_azimuth + PROJ.grid_convergence(x, y), 2 * np.pi)
        assert abs(elevation) < 1e-12
        assert abs(np.angle(np.exp(1j * (grid_azimuth - expected)))) < 1e-12


def test_elevation_azimuth_of_known_directions():
    e_x, e_y, up = np.eye(3)
    elevation, azimuth = elevation_azimuth(
        [[0, 1, 0], [1, 0, 0], [0, 0, 2], [0, 1, 1]], e_x, e_y, up
    )
    assert np.allclose(elevation, [0, 0, np.pi / 2, np.pi / 4])
    assert np.allclose(azimuth[[0, 1, 3]], [0, np.pi / 2, 0])


def test_radius_must_be_positive():
    with pytest.raises(ValueError):
        SouthPolarStereographic(0.0)
