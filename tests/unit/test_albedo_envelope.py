from __future__ import annotations

import numpy as np
import pytest

from lunarenv.thermal.albedo_laws import BURGER2024, CONSTANT, ENVELOPE, FOOTE2020, HAYNE2017
from lunarenv.thermal.properties import RegolithProperties, albedo
from lunarenv.thermal.surface import plane_normal, view_factors

BASE = RegolithProperties()


@pytest.mark.parametrize(
    ("law", "expected"),
    [(BURGER2024, 0.860), (HAYNE2017, 0.784), (FOOTE2020, 0.238), (CONSTANT, 0.12)],
)
def test_albedo_at_the_test_point_incidence(law, expected):
    # values at 88.2 deg quoted by the 2026-11-13 review for the DSNE Rev J envelope
    assert float(albedo(np.radians(88.2), law.apply(BASE))) == pytest.approx(expected, abs=5e-4)


def test_foote_quadratic_by_hand_and_envelope_order():
    i = 60.0
    by_hand = 1.3299e-5 * i**2 - 1.0248e-4 * i + 0.14382
    assert float(albedo(np.radians(i), FOOTE2020.apply(BASE))) == pytest.approx(by_hand)
    theta = np.radians(np.linspace(80.0, 90.0, 11))
    cold, hot = (albedo(theta, ENVELOPE[c].apply(BASE)) for c in ("cold", "hot"))
    assert np.all(cold > albedo(theta, BASE)) and np.all(albedo(theta, BASE) > hot)


def test_law_replaces_every_albedo_field():
    sampled = BASE.model_copy(update={"albedo_normal": 0.2, "emissivity": 0.97})
    hot = FOOTE2020.apply(sampled)
    assert hot.emissivity == 0.97 and hot.albedo_quadratic_deg is not None
    back = HAYNE2017.apply(hot)
    assert back.albedo_quadratic_deg is None and back.albedo_normal == 0.12


def test_view_factors_of_level_and_vertical_planes():
    azimuth = np.linspace(0.0, 2 * np.pi, 720, endpoint=False)
    horizon = np.radians(4.0) * (1 + np.cos(azimuth))
    level = view_factors(horizon, azimuth, plane_normal(0.0, 0.0))
    assert level.sky == pytest.approx(1 - np.mean(np.sin(horizon) ** 2), abs=1e-5)
    assert level.sky + level.ground == 1.0
    wall = view_factors(np.zeros(720), azimuth, plane_normal(np.pi / 2, 1.0))
    assert wall.sky == pytest.approx(0.5, abs=1e-4)
    # a horizon 2 degrees below the horizontal adds a band in front of the wall
    lower = view_factors(np.full(720, np.radians(-2.0)), azimuth, plane_normal(np.pi / 2, 0.0))
    band = (np.radians(2.0) / 2 + np.sin(2 * np.radians(2.0)) / 4) * 2 / np.pi
    assert lower.sky == pytest.approx(0.5 + band, abs=1e-4)
    with pytest.raises(ValueError):
        plane_normal(np.radians(91.0), 0.0)
