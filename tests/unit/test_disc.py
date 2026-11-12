"""Visible fraction of a disc: closed forms (ILL-A3), brute-force reference, properties."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from lunarenv.illumination.disc import (
    CHORD,
    PROFILE,
    disc_fractions,
    interpolate_horizon,
)

RHO = 0.00465  # solar angular radius, rad
N_AZ = 1440


def flat_horizon(elevation: float) -> np.ndarray:
    return np.full(N_AZ, elevation)


def segment_fraction(y0: float) -> float:
    """Area fraction of a unit disc above the horizontal line at height y0 in [-1, 1]."""
    return float((np.arccos(y0) - y0 * np.sqrt(1.0 - y0 * y0)) / np.pi)


def row_integrated_flux_fraction(y0: float, u: float) -> float:
    """Limb-darkened flux fraction of a unit disc above the line y = y0, in closed form.

    Integrates the linear limb-darkening law along horizontal rows (the tracer's kernel
    integrates along vertical strips): the row at height y contributes
    2 (1 - u) sqrt(1 - y^2) + u (pi / 2) (1 - y^2).
    """
    root = np.sqrt(1.0 - y0 * y0)
    uniform = 2.0 * (np.pi / 4.0 - 0.5 * (y0 * root + np.arcsin(y0)))
    darkened = 0.5 * np.pi * (2.0 / 3.0 - (y0 - y0**3 / 3.0))
    return float(((1.0 - u) * uniform + u * darkened) / (np.pi * (1.0 - u / 3.0)))


@pytest.mark.parametrize("y0", [-0.95, -0.5, -0.1, 0.0, 0.3, 0.8, 0.99])
@pytest.mark.parametrize("mode", [PROFILE, CHORD])
def test_ill_a3_straight_horizon_gives_the_circular_segment(y0, mode):
    elevation = 0.02
    horizon = flat_horizon(elevation + y0 * RHO)
    geometric, _ = disc_fractions(elevation, 1.0, RHO, horizon, 512, 0.6, mode)
    assert geometric == pytest.approx(segment_fraction(y0), abs=2e-6)


def test_strip_quadrature_converges_at_second_order_or_better():
    y0 = 0.37
    horizon = flat_horizon(0.02 + y0 * RHO)
    errors = [
        abs(disc_fractions(0.02, 1.0, RHO, horizon, n, 0.0, PROFILE)[0] - segment_fraction(y0))
        for n in (8, 16, 32, 64, 128)
    ]
    ratios = np.array(errors[:-1]) / np.array(errors[1:])
    assert errors[-1] < 1e-5
    assert np.exp(np.mean(np.log(ratios))) > 2.8  # at least ~h^1.5 on average; kink-limited


@pytest.mark.parametrize("y0", [-0.9, -0.7, 0.0, 0.45, 0.97])
@pytest.mark.parametrize("u", [0.3, 0.6, 1.0])
def test_ill_a3_limb_darkened_flux_matches_closed_form(y0, u):
    horizon = flat_horizon(0.02 + y0 * RHO)
    _, flux = disc_fractions(0.02, 1.0, RHO, horizon, 512, u, PROFILE)
    assert flux == pytest.approx(row_integrated_flux_fraction(y0, u), abs=2e-6)
    assert row_integrated_flux_fraction(-1.0, u) == pytest.approx(1.0, abs=1e-12)


def test_uniform_disc_has_equal_geometric_and_flux_fractions():
    horizon = flat_horizon(0.02 + 0.2 * RHO)
    geometric, flux = disc_fractions(0.02, 1.0, RHO, horizon, 64, 0.0, PROFILE)
    assert geometric == flux


def test_limb_darkening_weights_the_centre():
    """With the limb hidden, the visible share of the flux exceeds the share of the area."""
    horizon = flat_horizon(0.02 - 0.6 * RHO)  # lower 20 % of the disc hidden
    geometric, flux = disc_fractions(0.02, 1.0, RHO, horizon, 256, 0.6, PROFILE)
    assert flux > geometric


def test_profile_mode_resolves_structure_across_the_disc():
    """A cliff edge at the Sun's azimuth hides half of the disc; a chord cannot see that."""
    azimuth = np.pi
    n = 20 * N_AZ  # sampled finely enough that the edge is sharp on the scale of the disc
    fine = np.where(np.arange(n) * 2 * np.pi / n < azimuth, 0.5, -0.5)
    geometric, _ = disc_fractions(0.0, azimuth, RHO, fine, 256, 0.0, PROFILE)
    assert geometric == pytest.approx(0.5, abs=0.02)
    chord, _ = disc_fractions(0.0, azimuth, RHO, fine, 256, 0.0, CHORD)
    assert chord in (0.0, 1.0)


def test_horizon_interpolation_is_periodic_and_linear():
    horizon = np.array([0.0, 1.0, 2.0, 3.0])
    step = 2 * np.pi / 4
    assert interpolate_horizon(horizon, 0.5 * step) == pytest.approx(0.5)
    assert interpolate_horizon(horizon, 3.5 * step) == pytest.approx(1.5)  # wraps 3 -> 0
    assert interpolate_horizon(horizon, -0.5 * step) == pytest.approx(1.5)
    assert interpolate_horizon(horizon, 2 * np.pi + step) == pytest.approx(1.0)
    assert interpolate_horizon(horizon, 2 * np.pi) == pytest.approx(0.0)


profiles = st.lists(st.floats(-0.02, 0.06), min_size=32, max_size=32).map(np.array)


@given(
    profiles,
    st.floats(0.0, 0.04),
    st.floats(0.0, 2 * np.pi),
    st.floats(0.0, 0.02),
    st.floats(0.0, 1.0),
)
@settings(max_examples=300, deadline=None)
def test_fraction_bounds_and_monotonicity(profile, elevation, azimuth, lift, u):
    geometric, flux = disc_fractions(elevation, azimuth, RHO, profile, 32, u, PROFILE)
    assert 0.0 <= geometric <= 1.0 and 0.0 <= flux <= 1.0
    # raising the terrain everywhere can never reveal more of the disc
    lifted_geometric, lifted_flux = disc_fractions(
        elevation, azimuth, RHO, profile + lift, 32, u, PROFILE
    )
    assert lifted_geometric <= geometric + 1e-12 and lifted_flux <= flux + 1e-12
    # visible area and visible flux vanish together
    assert (geometric == 0.0) == (flux == 0.0)


@given(st.floats(-0.05, 0.05), st.floats(0.0, 2 * np.pi), st.floats(0.0, 1.0))
def test_disc_clear_of_the_horizon_is_fully_visible_or_fully_hidden(elevation, azimuth, u):
    above = disc_fractions(
        elevation, azimuth, RHO, flat_horizon(elevation - 1.001 * RHO), 32, u, PROFILE
    )
    below = disc_fractions(
        elevation, azimuth, RHO, flat_horizon(elevation + 1.001 * RHO), 32, u, PROFILE
    )
    assert above == pytest.approx((1.0, 1.0), abs=1e-12)
    assert below == (0.0, 0.0)


@given(profiles, st.floats(0.0, 0.04), st.floats(0.0, 2 * np.pi))
@settings(max_examples=100, deadline=None)
def test_mirror_symmetry(profile, elevation, azimuth):
    """Reflecting the sky about azimuth 0 leaves the visible fraction unchanged."""
    mirrored = np.roll(profile[::-1], 1)  # sample k -> sample -k (mod n)
    direct = disc_fractions(elevation, azimuth, RHO, profile, 32, 0.6, PROFILE)
    reflected = disc_fractions(elevation, 2 * np.pi - azimuth, RHO, mirrored, 32, 0.6, PROFILE)
    assert direct == pytest.approx(reflected, abs=1e-9)
