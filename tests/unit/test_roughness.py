import numpy as np
import pytest

from lunarenv.thermal.model import ThermalSettings
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.roughness import (
    facet_normals,
    facet_set,
    lit_probability,
    rough_radiance,
    smith_lambda,
    sun_vectors,
)


@pytest.mark.parametrize("rms", [10.0, 20.0, 35.0])
def test_facet_set_reproduces_the_slope_moments(rms: float) -> None:
    f = facet_set(rms)
    s2 = np.tan(np.radians(rms)) ** 2
    assert f.weight.sum() == pytest.approx(1.0)
    assert np.sum(f.weight * np.sum(f.gradient**2, axis=1)) == pytest.approx(s2, rel=1e-12)
    assert np.allclose(np.sum(f.weight[:, None] * f.gradient, axis=0), 0.0, atol=1e-12)
    # each component carries half the variance
    assert np.sum(f.weight * f.gradient[:, 0] ** 2) == pytest.approx(s2 / 2, rel=1e-12)


@pytest.mark.parametrize("nu", [0.2, 0.7, 1.5, 3.0])
def test_smith_lambda_matches_its_definition(nu: float) -> None:
    # Smith 1967: Lambda(mu) = (1/mu) * int_mu^inf (q - mu) p(q) dq, p the pdf of the
    # gradient component along the Sun azimuth (std s / sqrt(2)); nu = mu / s.
    s = 0.4
    sigma_c = s / np.sqrt(2.0)
    mu = nu * s
    q = np.linspace(mu, mu + 12 * sigma_c, 200_001)
    pdf = np.exp(-0.5 * (q / sigma_c) ** 2) / (sigma_c * np.sqrt(2 * np.pi))
    integral = np.trapezoid((q - mu) * pdf, q)
    assert smith_lambda(nu) == pytest.approx(integral / mu, rel=1e-6)


def test_lit_probability_limits() -> None:
    cz = np.array([1.0, 0.5, 0.05, 0.0, -0.2])
    g = lit_probability(cz, 20.0)
    assert g[0] == 1.0
    assert g[3] == 0.0 and g[4] == 0.0
    assert np.all(np.diff(g[:4]) < 0)
    # rougher surfaces shadow more
    assert lit_probability(np.array([0.1]), 35.0)[0] < lit_probability(np.array([0.1]), 10.0)[0]


def test_facet_normals_are_unit_and_tilted_about_the_mean_normal() -> None:
    n0 = np.array([0.1, -0.05, 1.0])
    n0 /= np.linalg.norm(n0)
    f = facet_set(20.0)
    n = facet_normals(n0, f.gradient)
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)
    slope = np.arctan(np.hypot(f.gradient[:, 0], f.gradient[:, 1]))
    assert np.allclose(np.arccos(n @ n0), slope)
    assert np.allclose(facet_normals(n0, np.zeros((1, 2)))[0], n0)


def test_rough_radiance_tends_to_smooth_for_tiny_roughness() -> None:
    step = 3600.0
    hours = np.arange(int(29.53 * 24 * 3))
    hour_angle = 2 * np.pi * hours / (29.53 * 24)
    lat = np.radians(45.0)
    elevation = np.arcsin(np.cos(lat) * np.cos(hour_angle))
    azimuth = np.where(np.sin(hour_angle) > 0, np.pi / 2, -np.pi / 2)  # east-west only
    sun = sun_vectors(elevation, azimuth)
    irradiance = np.where(elevation > 0, 1361.0, 0.0)
    spin = int(29.53 * 24 * 2)
    props = RegolithProperties()
    out = rough_radiance(
        irradiance, sun, [0.0, 0.0, 1.0], props, [0.01, 20.0], step, spin, ThermalSettings()
    )
    lit = elevation[spin:] > 0.2
    t_smooth = out.smooth**0.25
    assert np.allclose(out.rough[0.01][lit] ** 0.25, t_smooth[lit], atol=0.05)
    assert not np.allclose(out.rough[20.0][lit] ** 0.25, t_smooth[lit], atol=0.5)
    # at night no facet is lit: each follows its own (daytime-heated) history, so the
    # rough surface stays close to the smooth one rather than to a never-lit population
    night = elevation[spin:] < -0.3
    assert np.allclose(out.rough[20.0][night] ** 0.25, t_smooth[night], atol=3.0)
