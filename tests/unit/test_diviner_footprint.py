import numpy as np
import pytest

from lunarenv.validation.diviner_footprint import (
    Field,
    Footprint,
    bin_points,
    block_layout,
    cell_radiance,
    ground_speed,
    response,
)
from lunarenv.validation.diviner_pcp import cell_pixels


def test_nominal_footprint_at_50_km() -> None:
    fp = Footprint.at_altitude(50e3)
    assert fp.in_track_m == pytest.approx(335.0)
    assert fp.cross_track_m == pytest.approx(170.0)
    assert ground_speed(50e3) == pytest.approx(1610.0, abs=2.0)
    assert fp.smear_m == pytest.approx(206.0, abs=0.5)


@pytest.mark.parametrize("theta", [0.0, 0.4, 2.0])
def test_response_is_a_normalised_symmetric_trapezoid(theta: float) -> None:
    fp = Footprint(335.0, 170.0, 206.0)
    off, w = response(fp, theta, 5.0)
    assert w.sum() == pytest.approx(1.0)
    assert np.allclose(np.sum(w[:, None] * off, axis=0), 0.0, atol=1e-9)
    u = np.array([np.cos(theta), np.sin(theta)])
    s = off @ u
    assert np.abs(s).max() <= (335.0 + 206.0) / 2
    # in-track variance of a box (L) convolved with a box (S): (L^2 + S^2) / 12
    assert np.sum(w * s**2) == pytest.approx((335.0**2 + 206.0**2) / 12, rel=2e-3)


def test_bin_points_are_the_ifov_corners_and_edges() -> None:
    fp = Footprint(300.0, 100.0, 50.0)
    p = bin_points(fp, 0.0)
    assert p.shape == (9, 2)
    assert np.allclose(p[0], 0.0)
    assert np.abs(p[:, 0]).max() == pytest.approx(150.0)
    assert np.abs(p[:, 1]).max() == pytest.approx(50.0)


def _linear_field(a: float, b: float, c: float) -> Field:
    xs = np.arange(0.0, 2000.0 + 1e-9, 20.0)
    xx, yy = np.meshgrid(xs, xs)
    values = (a + b * xx + c * yy)[:, :, None] * np.array([1.0, 2.0])
    return Field(0.0, 0.0, 20.0, values)


@pytest.mark.parametrize("theta", [0.0, 0.7])
def test_linear_field_is_unchanged_by_a_symmetric_footprint(theta: float) -> None:
    field = _linear_field(5.0, 0.01, -0.02)
    cx, cy = cell_pixels(1000.0, 1000.0, 20.0)
    point = cell_radiance(field, cx, cy, Footprint(0.0, 0.0, 0.0), theta, nine_points=False)
    full = cell_radiance(field, cx, cy, Footprint.at_altitude(), theta)
    assert np.allclose(point, [5.0 + 10.0 - 20.0, 2 * (5.0 + 10.0 - 20.0)])
    assert np.allclose(full, point)


def test_point_footprint_reproduces_the_cell_mean() -> None:
    rng = np.random.default_rng(1)
    values = rng.random((60, 60, 3))
    field = Field(0.0, 0.0, 20.0, values)
    cx, cy = cell_pixels(590.0, 590.0, 20.0)  # nodes 24..35 in both directions
    got = cell_radiance(field, cx, cy, Footprint(0.0, 0.0, 0.0), 0.3, nine_points=False)
    assert np.allclose(got, values[24:36, 24:36].reshape(-1, 3).mean(axis=0))


def test_bright_cell_in_dark_surroundings_is_diluted() -> None:
    values = np.zeros((60, 60, 1))
    values[24:36, 24:36] = 1.0
    field = Field(0.0, 0.0, 20.0, values)
    cx, cy = cell_pixels(590.0, 590.0, 20.0)
    got = cell_radiance(field, cx, cy, Footprint.at_altitude(), 0.0)[0]
    assert 0.0 < got < 0.8


def test_footprint_outside_the_field_is_refused() -> None:
    field = _linear_field(0.0, 0.0, 0.0)
    cx, cy = cell_pixels(130.0, 130.0, 20.0)
    with pytest.raises(ValueError):
        cell_radiance(field, cx, cy, Footprint.at_altitude(), 0.0)


def test_layout_covers_the_ring_and_reassembles() -> None:
    lay = block_layout(-1000.0, 2000.0)
    ix, _ = lay.inner_points()
    rx, ry, _, _ = lay.ring_points()
    assert ix.size == 36 * 36
    assert lay.n == 102
    # ring pixels tile the 2040 m square outside the 720 m block
    assert rx.size * 60.0**2 + ix.size * 20.0**2 == pytest.approx(2040.0**2)
    assert ix.min() == pytest.approx(-990.0) and ix.max() == pytest.approx(-290.0)
    inside = (rx > -1000) & (rx < -280) & (ry > 2000) & (ry < 2720)
    assert not inside.any()
    field = lay.assemble(np.c_[ix + 0.0], np.c_[rx + 0.0])
    # node values equal x within half a ring pixel
    xs = lay.x0 + 20.0 * np.arange(lay.n)
    assert np.abs(field.values[:, :, 0] - xs[None, :]).max() <= 20.0 + 1e-9
    cx, cy = cell_pixels(-880.0, 2120.0, 20.0)
    nominal = Footprint.at_altitude()
    got = cell_radiance(field, cx, cy, nominal.scaled(1.5), 1.0, binning=nominal)
    assert np.isfinite(got).all()
