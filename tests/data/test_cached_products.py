"""Regression and consistency checks on the real cached products."""

from __future__ import annotations

import numpy as np
import pytest

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M
from lunarenv.data.dem import Bounds, ingest_crop, load_dem, read_dem, read_raster

pytestmark = pytest.mark.data

SITE04 = "pgda78.site04.surf"
REGIONAL = "pgda90.ldem_80s_80mpp"


def test_cached_files_match_manifest_and_pins(store):
    status = store.verify()
    if not status:
        pytest.skip("the data cache is empty")
    assert {k: v for k, v in status.items() if v != "ok"} == {}


def test_site04_grid_and_height_statistics(store, require):
    """Pins the unit (metres) and extent of the Site 4 product (PGDA-78)."""
    require(SITE04)
    dem = load_dem(store, SITE04)
    assert dem.grid.shape == (3200, 3200) and dem.grid.pixel_size_m == 5.0
    assert dem.grid.bounds == Bounds(-9000.0, -15000.0, 7000.0, 1000.0)
    assert dem.grid.projection.radius_m == MOON_REFERENCE_RADIUS_M
    assert dem.source_sha256 == store.record(SITE04).sha256
    assert float(dem.heights_m.min()) == pytest.approx(-2847.877, abs=1e-3)
    assert float(dem.heights_m.max()) == pytest.approx(1805.073, abs=1e-3)


def test_regional_grid_and_height_statistics(store, require):
    require(REGIONAL)
    dem = load_dem(store, REGIONAL)
    assert dem.grid.shape == (7600, 7600) and dem.grid.pixel_size_m == 80.0
    assert dem.grid.bounds == Bounds(-304000.0, -304000.0, 304000.0, 304000.0)
    assert float(dem.heights_m.min()) == pytest.approx(-7295.517, abs=1e-2)
    assert float(dem.heights_m.max()) == pytest.approx(7026.331, abs=1e-2)


def test_site_and_regional_products_are_co_registered(store, require):
    """Independent products must describe the same terrain in the same place.

    The 5 m tile is block-averaged to 80 m and compared with the 80 m product. If either
    reader mis-handled pixel registration or axis order, the misfit would be minimized
    at a non-zero shift (or be enormous: the tile has 4.6 km of relief).
    """
    require(SITE04, REGIONAL)
    # The 80 m lattice is offset by 40 m from the tile edges, so compare on the largest
    # area made of whole 80 m cells inside the tile (199 x 199 cells of 16 x 16 pixels).
    common = Bounds(-8960.0, -14960.0, 6960.0, 960.0)
    site = load_dem(store, SITE04, common)
    regional = load_dem(store, REGIONAL, common)
    assert regional.grid.bounds == site.grid.bounds == common
    assert regional.grid.shape == (199, 199) and site.grid.shape == (3184, 3184)
    coarse = site.heights_m.astype(np.float64).reshape(199, 16, 199, 16).mean(axis=(1, 3))
    reference = regional.heights_m.astype(np.float64)

    def rms(shift_rows: int, shift_cols: int) -> float:
        a = coarse[5:-5, 5:-5]
        b = reference[5 + shift_rows : 194 + shift_rows, 5 + shift_cols : 194 + shift_cols]
        return float(np.sqrt(np.mean((a - b) ** 2)))

    shifts = {(dr, dc): rms(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)}
    assert min(shifts, key=shifts.get) == (0, 0)
    print(
        f"site vs regional at 80 m: rms {shifts[(0, 0)]:.2f} m; shifted by one cell: "
        f"{min(v for k, v in shifts.items() if k != (0, 0)):.1f} m"
    )
    assert shifts[(0, 0)] < 2.0  # metres; a one-pixel shift gives tens of metres
    assert min(v for k, v in shifts.items() if k != (0, 0)) > 5.0 * shifts[(0, 0)]


def test_error_clone_is_a_full_dem_consistent_with_the_error_map(store, require):
    """Clones are heights (not perturbations) differing from the nominal DEM by ~toterr."""
    require(SITE04, "pgda78.site04.clone.0001", "pgda78.site04.toterr")
    nominal = load_dem(store, SITE04)
    clone = load_dem(store, "pgda78.site04.clone.0001")
    error_map, grid = read_raster(store.path("pgda78.site04.toterr"))
    assert clone.grid == nominal.grid == grid
    difference = clone.heights_m.astype(np.float64) - nominal.heights_m
    assert abs(difference.mean()) < 0.01
    rms = float(np.sqrt(np.mean(difference**2)))
    formal = float(np.sqrt(np.mean(error_map.astype(np.float64) ** 2)))
    assert 0.3 < rms < 0.6  # Barker2021: median RMS height error 0.30-0.50 m
    assert 0.5 < rms / formal < 2.0


def test_crop_ingest_registers_provenance(store, require, tmp_path):
    require(SITE04)
    from lunarenv.data.store import DataStore

    # a scratch store sharing the real raw file, so the real manifest is not modified
    scratch = DataStore(tmp_path, catalog=store.catalog)
    target = scratch._raw_path(store.catalog.get(SITE04))
    target.parent.mkdir(parents=True)
    target.symlink_to(store.path(SITE04))
    scratch.fetch(SITE04)  # adopted: content equals the pin, no download

    window = Bounds.around(-1000.0, -7000.0, 500.0)
    record = ingest_crop(scratch, SITE04, window)
    assert record.parent_sha256 == store.record(SITE04).sha256
    crop = read_dem(scratch.derived_path(record.id))
    expected = load_dem(store, SITE04, window)
    assert crop.grid == expected.grid and np.array_equal(crop.heights_m, expected.heights_m)
