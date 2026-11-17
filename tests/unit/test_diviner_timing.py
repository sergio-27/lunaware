import numpy as np
import pytest

from lunarenv.validation.cases_thermal import _timed_cells
from lunarenv.validation.diviner_pcp import N_BINS, PERIOD_START, at_hours, hour_index


def test_hour_index_of_the_period_start_and_one_day_later() -> None:
    # 2009-07-05T00:00 UTC is JD 2455017.5
    assert PERIOD_START == np.datetime64("2009-07-05T00:00")
    assert hour_index(np.array([2455017.5, 2455018.5, 2455017.5 + 1.5 / 24])) == pytest.approx(
        [0.0, 24.0, 1.5]
    )


def test_at_hours_interpolates_linearly_between_samples() -> None:
    series = np.array([0.0, 10.0, 30.0, 60.0])
    assert at_hours(series, np.array([0.0, 0.5, 1.25, 3.0])) == pytest.approx(
        [0.0, 5.0, 15.0, 60.0]
    )


def test_at_hours_refuses_times_outside_the_series() -> None:
    with pytest.raises(ValueError):
        at_hours(np.arange(4.0), np.array([3.5]))
    with pytest.raises(ValueError):
        at_hours(np.arange(4.0), np.array([-0.5]))


def test_timed_cells_weights_observations_by_membership() -> None:
    values = np.array([[1.0, 2.0, 4.0], [10.0, 20.0, 40.0]])
    member = np.array([[1.0, 1.0, 0.0], [0.0, 2.0, 1.0]])
    bins = np.array([5, 5, 5])
    out = _timed_cells(values, member, bins, 2)
    assert out.shape == (2, N_BINS)
    assert out[0, 5] == pytest.approx(1.5)
    assert out[1, 5] == pytest.approx((2 * 20.0 + 40.0) / 3)
    assert np.isnan(out[:, 4]).all()
