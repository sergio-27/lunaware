"""Sampling of regolith parameters."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.uq import (
    DEFAULT_RANGES,
    ParameterRange,
    latin_hypercube,
    rank_correlation,
    sample_properties,
)


def test_latin_hypercube_stratifies_every_dimension_and_is_reproducible():
    sample = latin_hypercube(50, 4, seed=7)
    assert sample.shape == (50, 4) and np.all((sample >= 0) & (sample < 1))
    for d in range(4):  # exactly one point in each of the 50 strata
        assert sorted(np.floor(sample[:, d] * 50).astype(int)) == list(range(50))
    assert np.array_equal(sample, latin_hypercube(50, 4, seed=7))
    assert not np.array_equal(sample, latin_hypercube(50, 4, seed=8))


def test_default_ranges_contain_the_nominal_values_and_name_real_fields():
    nominal = RegolithProperties()
    for parameter in DEFAULT_RANGES:
        value = getattr(nominal, parameter.field)
        assert parameter.low <= value <= parameter.high, parameter.field
        assert parameter.basis
    assert len({p.field for p in DEFAULT_RANGES}) == len(DEFAULT_RANGES)


def test_sampled_properties_lie_in_their_intervals():
    members, values = sample_properties(40, seed=1)
    assert len(members) == 40 and values.shape == (40, len(DEFAULT_RANGES))
    for d, parameter in enumerate(DEFAULT_RANGES):
        column = np.array([getattr(m, parameter.field) for m in members])
        assert np.array_equal(column, values[:, d])
        assert column.min() >= parameter.low and column.max() <= parameter.high
        assert column.max() - column.min() > 0.9 * (parameter.high - parameter.low)
    assert all(m.radiative_parameter == 2.7 for m in members)  # not varied
    with pytest.raises(ValueError):
        sample_properties(4, 0, ranges=(ParameterRange("emissivity", 0.9, 0.9, "x"),))


def test_rank_correlation_recovers_monotonic_dependence():
    rng = np.random.default_rng(0)
    x = rng.random((200, 3))
    y = np.exp(3 * x[:, 0]) - 0.2 * x[:, 1]
    rho = rank_correlation(x, y)
    assert rho[0] > 0.95 and rho[1] < 0.0 and abs(rho[2]) < 0.2
