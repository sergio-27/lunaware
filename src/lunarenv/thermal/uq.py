"""Uncertainty of the regolith parameters and its propagation by sampling.

Each parameter is given an interval and sampled uniformly within it by Latin hypercube
sampling (McKay1979): with ``n`` samples every parameter visits each of its ``n``
equal-probability strata exactly once, which estimates output quantiles with fewer
runs than plain random sampling.

The intervals are *ranges of plausible values*, not measured distributions. Where the
literature gives a spread it is used; where it does not, the interval is an assumption
and is marked as such. They are inputs the user can and should override for a
specific site.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from lunarenv.geometry.projection import FloatArray
from lunarenv.thermal.properties import RegolithProperties


@dataclass(frozen=True, slots=True)
class ParameterRange:
    """Interval sampled for one field of :class:`RegolithProperties`."""

    field: str
    low: float
    high: float
    basis: str


#: Default intervals. "Hayne2017" refers to the text and tables of that paper as read;
#: values it quotes from other authors are attributed in the basis string.
DEFAULT_RANGES: tuple[ParameterRange, ...] = (
    ParameterRange(
        "scale_height_m",
        0.047,
        0.089,
        "Hayne2017: global mode 0.068 m, standard deviation 0.007 m; interval is +/-3 sigma",
    ),
    ParameterRange(
        "surface_conductivity_w_m_k",
        6.0e-4,
        9.0e-4,
        "Hayne2017 Table A2: 0.6e-3 (Vasavada2012) to 0.9e-3 (lower end of Keihm1973 in situ)",
    ),
    ParameterRange(
        "deep_conductivity_w_m_k",
        3.4e-3,
        7.0e-3,
        "Hayne2017: 3.4e-3 adopted; Table A2 quotes 7.0e-3 at ~1 m from Vasavada2012",
    ),
    ParameterRange(
        "surface_density_kg_m3",
        1000.0,
        1300.0,
        "ASSUMPTION: +/-10-20 % around the 1100 kg m-3 of Hayne2017 Table A1",
    ),
    ParameterRange(
        "deep_density_kg_m3", 1800.0, 1900.0, "Hayne2017 Table A2: 1800-1900 kg m-3 at 30-60 cm"
    ),
    ParameterRange("emissivity", 0.95, 0.98, "Hayne2017: 0.95 adopted, 0.98 was the earlier value"),
    ParameterRange(
        "albedo_normal",
        0.12,
        0.25,
        "Hayne2017: albedo ~5 % (maria) to ~25 % (highlands); 0.12 is the lunar average; "
        "the site is highland terrain",
    ),
    ParameterRange(
        "heat_flow_w_m2",
        0.012,
        0.024,
        "ASSUMPTION: +/-33 % around the 0.018 W m-2 of Hayne2017 Table A1; polar value unknown",
    ),
)
#: Parameters held fixed: the radiative parameter chi (no published spread was accessed)
#: and the heat-capacity polynomial.


def latin_hypercube(n_samples: int, n_dimensions: int, seed: int) -> FloatArray:
    """Latin hypercube sample on the unit cube, shape ``(n_samples, n_dimensions)``.

    Deterministic for a given seed (``numpy.random.Generator`` with PCG64).
    """
    if n_samples < 1 or n_dimensions < 1:
        raise ValueError("need at least one sample and one dimension")
    rng = np.random.default_rng(seed)
    strata = (np.arange(n_samples)[:, None] + rng.random((n_samples, n_dimensions))) / n_samples
    for d in range(n_dimensions):
        strata[:, d] = strata[rng.permutation(n_samples), d]
    return np.asarray(strata)


def sample_properties(
    n_samples: int,
    seed: int,
    *,
    base: RegolithProperties | None = None,
    ranges: tuple[ParameterRange, ...] = DEFAULT_RANGES,
) -> tuple[list[RegolithProperties], FloatArray]:
    """Property sets drawn from the parameter intervals, and the sampled values.

    Returns the list of property objects and an array ``(n_samples, len(ranges))`` of
    the sampled values in the order of ``ranges``.
    """
    base = base or RegolithProperties()
    unit = latin_hypercube(n_samples, len(ranges), seed)
    values = np.empty_like(unit)
    for d, parameter in enumerate(ranges):
        if not parameter.high > parameter.low:
            raise ValueError(f"empty interval for {parameter.field}")
        values[:, d] = parameter.low + unit[:, d] * (parameter.high - parameter.low)
    members = [
        base.model_copy(update={p.field: float(v) for p, v in zip(ranges, row, strict=True)})
        for row in values
    ]
    return [RegolithProperties.model_validate(m.model_dump()) for m in members], values


def rank_correlation(inputs: FloatArray, output: FloatArray) -> FloatArray:
    """Spearman rank correlation of each input column with an output (sensitivity ranking)."""

    def ranks(a: FloatArray) -> FloatArray:
        return np.asarray(np.argsort(np.argsort(a, axis=0), axis=0), dtype=np.float64)

    x = ranks(np.asarray(inputs, dtype=np.float64))
    y = ranks(np.asarray(output, dtype=np.float64).reshape(-1, 1))[:, 0]
    x = x - x.mean(axis=0)
    y = y - y.mean()
    denominator = np.sqrt((x**2).sum(axis=0) * (y**2).sum())
    return np.asarray((x * y[:, None]).sum(axis=0) / np.where(denominator > 0, denominator, 1.0))
