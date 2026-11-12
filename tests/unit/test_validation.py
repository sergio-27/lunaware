"""Validation framework: verdict logic, runner behaviour, stored results, report."""

from __future__ import annotations

import numpy as np
import pytest
import rasterio

from lunarenv.data.dem import read_raster
from lunarenv.data.store import DataStore
from lunarenv.validation import core
from lunarenv.validation.cases_illumination import BARKER_TABLE2, _region
from lunarenv.validation.core import (
    Case,
    Check,
    Outcome,
    Status,
    get_validation_status,
    load_report,
    registry,
    render_markdown,
    run_cases,
    save_report,
)


def test_check_verdicts():
    assert Check(quantity="q", reference=10.0, measured=10.4, tolerance=0.5).passed is True
    assert Check(quantity="q", reference=10.0, measured=10.6, tolerance=0.5).passed is False
    assert (
        Check(quantity="q", reference=10.0, measured=99.0).passed is None
    )  # no tolerance: information
    assert (
        Check(quantity="q", reference=10.0, measured=99.0, tolerance=1.0, informational=True).passed
        is None
    )
    assert Check(quantity="q", measured=1.0, tolerance=1.0).passed is None


def test_registry_is_consistent():
    cases = registry()
    ids = [c.id for c in cases]
    assert len(ids) == len(set(ids))
    for case in cases:
        assert case.source and case.metric and case.tolerance
        assert (case.run is None) == bool(case.pending_reason), case.id  # pending cases explain why
    # at least one published comparison per model, as the project definition requires
    published = {c.id[:2] for c in cases if c.kind == "published" and c.run is not None}
    assert published == {"IL", "TH"}


def _fake_registry(calls):
    def good(_store):
        calls.append("good")
        return Outcome([Check(quantity="x", reference=1.0, measured=1.0, tolerance=0.1)], "fine")

    def bad(_store):
        calls.append("bad")
        return Outcome(
            [Check(quantity="x", reference=1.0, measured=2.0, tolerance=0.1)], "off", ("why",)
        )

    def slow(_store):
        calls.append("slow")
        return Outcome(
            [Check(quantity="x", reference=1.0, measured=1.0, informational=True)], "context"
        )

    return (
        Case("T-1", "good", "analytic", "s", "m", "t", run=good),
        Case("T-2", "bad", "published", "s", "m", "t", run=bad),
        Case("T-3", "slow", "published", "s", "m", "t", slow=True, run=slow),
        Case("T-4", "needs data", "published", "s", "m", "t", requires=("naif.lsk",), run=good),
        Case("T-5", "pending", "published", "s", "m", "t", pending_reason="paper not accessed"),
    )


def test_runner_statuses_selection_and_retention(tmp_path, monkeypatch):
    calls: list[str] = []
    monkeypatch.setattr(core, "registry", lambda: _fake_registry(calls))
    store = DataStore(tmp_path / "empty-cache")
    quick = run_cases(store, include_slow=False)
    status = {r.id: r.status for r in quick.results}
    assert status == {
        "T-1": Status.PASS,
        "T-2": Status.FAIL,
        "T-3": Status.NOT_RUN,
        "T-4": Status.NOT_RUN,
        "T-5": Status.PENDING,
    }
    assert "naif.lsk" in quick.by_id()["T-4"].summary and quick.by_id()["T-2"].notes == ["why"]
    assert calls == ["good", "bad"]

    full = run_cases(store, ["T-3"], previous=quick)  # only the slow one is re-run
    assert calls[-1] == "slow" and len(calls) == 3
    assert full.by_id()["T-3"].status is Status.INFO
    assert full.by_id()["T-2"] == quick.by_id()["T-2"]  # kept from the previous report

    path = tmp_path / "results.json"
    save_report(full, path)
    assert load_report(path) == full
    summary = get_validation_status(path)
    assert summary["counts"] == {"PASS": 1, "FAIL": 1, "INFO": 1, "NOT RUN": 1, "PENDING": 1}
    assert [c["id"] for c in summary["cases"]] == ["T-1", "T-2", "T-3", "T-4", "T-5"]

    text = render_markdown(full)
    assert "| T-2 | bad | s | m | t | off | **FAIL** |" in text
    assert "| x | 1 | 2 | 1 | 0.1 | FAIL |" in text and "- why" in text
    assert "paper not accessed" in text


def test_table_cells_escape_pipes():
    """A metric such as max |a - b| must not break the Markdown table."""
    report = core.ValidationReport(
        results=[
            core.CaseResult(
                id="T-9",
                title="t",
                kind="analytic",
                source="s",
                metric="max |a - b|",
                tolerance="1",
                status=Status.PASS,
                summary="ok",
                checks=[Check(quantity="|x|", reference=0.0, measured=0.0, tolerance=1.0)],
            )
        ]
    )
    text = render_markdown(report)
    row = next(line for line in text.splitlines() if line.startswith("| T-9"))
    assert "max \\|a - b\\|" in row and row.replace("\\|", "").count("|") == 8
    assert "| \\|x\\| |" in text


def test_status_without_stored_results(tmp_path):
    summary = get_validation_status(tmp_path / "missing.json")
    assert summary["generated_utc"] is None
    assert all(c["status"] == "NOT RUN" for c in summary["cases"])


def test_barker_table_is_transcribed_consistently():
    """Internal consistency of the transcribed table: 7 regions, centroids on pixel centres."""
    assert [row[0] for row in BARKER_TABLE2] == list(range(1, 8))
    for _, x_km, y_km, area, one_metre, five_metres in BARKER_TABLE2:
        for coordinate in (x_km, y_km):
            assert round(coordinate * 1000 % 5, 6) == 2.5  # centres of the 5 m lattice
        assert area % 25 == 0  # whole pixels
        assert (
            five_metres[0] > one_metre[0] and five_metres[2] < one_metre[2]
        )  # higher is better lit


def test_region_growing():
    mean = np.zeros((6, 6))
    mean[1:3, 1:3] = 0.9
    mean[3, 3] = 0.8  # touches the block diagonally
    mean[5, 5] = 0.95  # separate
    region = _region(mean, (1, 1), 0.7)
    assert region.sum() == 5 and region[3, 3] and not region[5, 5]
    assert not _region(mean, (0, 0), 0.7).any()


def test_read_raster_applies_scale_and_nodata(tmp_path):
    """Packed integer rasters (the PGDA illumination maps) are returned in physical units."""
    path = tmp_path / "packed.tif"
    crs = rasterio.crs.CRS.from_proj4(
        "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +x_0=0 +y_0=0 +R=1737400 +units=m +no_defs"
    )
    data = np.array([[0, 25000], [-32768, 12500]], dtype=np.int16)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=2,
        height=2,
        count=1,
        dtype="int16",
        crs=crs,
        transform=rasterio.transform.from_origin(0.0, 120.0, 60.0, 60.0),
        nodata=-32768,
    ) as dataset:
        dataset.write(data, 1)
        dataset.scales = (4e-5,)
    values, grid = read_raster(path, allow_nodata=True)
    assert grid.pixel_size_m == 60.0
    assert values[0, 1] == pytest.approx(1.0) and values[1, 1] == pytest.approx(0.5)
    assert values[0, 0] == 0.0 and np.isnan(values[1, 0])
