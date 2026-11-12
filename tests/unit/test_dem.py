"""DEM grid conventions, strict georeferencing checks, interpolation and cropping."""

from __future__ import annotations

import numpy as np
import pytest
import rasterio
from hypothesis import given
from hypothesis import strategies as st

from lunarenv.data.dem import (
    Bounds,
    Dem,
    GeoreferencingError,
    GridGeometry,
    read_dem,
    read_raster,
    write_dem,
)
from lunarenv.geometry.projection import SouthPolarStereographic

PROJ = SouthPolarStereographic()
GRID = GridGeometry(
    x_min=-100.0, y_max=250.0, pixel_size_m=5.0, width=40, height=30, projection=PROJ
)


def plane_dem(a=0.3, b=-0.7, c=12.0, grid=GRID) -> Dem:
    xx, yy = np.meshgrid(grid.x_centers(), grid.y_centers())
    return Dem((a * xx + b * yy + c).astype(np.float32), grid)


def test_pixel_centres_and_row_order():
    assert GRID.x_centers()[0] == -97.5 and GRID.x_centers()[-1] == 97.5
    assert GRID.y_centers()[0] == 247.5 and GRID.y_centers()[-1] == 102.5
    assert GRID.bounds == Bounds(-100.0, 100.0, 100.0, 250.0)
    row, col = GRID.fractional_index(-97.5, 247.5)
    assert row == 0.0 and col == 0.0


def test_window_snaps_outwards_to_whole_pixels():
    assert GRID.window(Bounds(-91.0, 201.0, -79.0, 239.0)) == (2, 10, 1, 5)
    # bounds already on pixel edges are not enlarged
    assert GRID.window(Bounds(-90.0, 200.0, -80.0, 240.0)) == (2, 10, 2, 4)
    assert GRID.window(GRID.bounds) == (0, 30, 0, 40)


def test_window_outside_grid_raises():
    with pytest.raises(ValueError, match="beyond"):
        GRID.window(Bounds(-200.0, 100.0, 0.0, 250.0))


def test_crop_preserves_georeferencing():
    dem = plane_dem()
    crop = dem.crop(Bounds(-90.0, 200.0, -80.0, 240.0))
    assert crop.grid.shape == (8, 2)
    assert crop.grid.x_centers()[0] == -87.5 and crop.grid.y_centers()[0] == 237.5
    assert crop.sample(-87.5, 237.5) == dem.sample(-87.5, 237.5)


@given(st.floats(-97.5, 97.5), st.floats(102.5, 247.5), st.floats(-2, 2), st.floats(-2, 2))
def test_bilinear_sampling_is_exact_for_planes(x, y, a, b):
    dem = plane_dem(a, b, 5.0)
    assert abs(dem.sample(x, y) - (a * x + b * y + 5.0)) < 2e-3  # float32 storage


def test_sampling_reproduces_pixel_values_and_rejects_outside():
    rng = np.random.default_rng(0)
    dem = Dem(rng.normal(size=GRID.shape).astype(np.float32), GRID)
    xx, yy = np.meshgrid(GRID.x_centers(), GRID.y_centers())
    assert np.array_equal(dem.sample(xx, yy).astype(np.float32), dem.heights_m)
    with pytest.raises(ValueError, match="outside"):
        dem.sample(-99.0, 200.0)  # between the grid edge and the first pixel centre


def test_write_then_read_round_trip(tmp_path):
    dem = plane_dem()
    path = tmp_path / "dem.tif"
    write_dem(dem, path)
    back = read_dem(path, source_id="x", source_sha256="y")
    assert back.grid == dem.grid and np.array_equal(back.heights_m, dem.heights_m)
    assert (back.source_id, back.source_sha256) == ("x", "y")
    part = read_dem(path, Bounds(-91.0, 201.0, -79.0, 239.0))
    assert np.array_equal(part.heights_m, dem.heights_m[2:10, 1:5])


def _write(path, *, proj4, transform, data=None, nodata=None):
    data = np.zeros((4, 4), np.float32) if data is None else data
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype="float32",
        crs=rasterio.crs.CRS.from_proj4(proj4),
        transform=transform,
        nodata=nodata,
    ) as dataset:
        dataset.write(data, 1)


GOOD = "+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +x_0=0 +y_0=0 +R=1737400 +units=m +no_defs"
NORTH_UP = rasterio.transform.from_origin(0.0, 20.0, 5.0, 5.0)


@pytest.mark.parametrize(
    ("proj4", "transform", "message"),
    [
        (GOOD.replace("R=1737400", "R=1737000"), NORTH_UP, "radius"),
        (GOOD.replace("lat_0=-90 +lat_ts=-90", "lat_0=90 +lat_ts=90"), NORTH_UP, "lat_0"),
        (GOOD.replace("lon_0=0", "lon_0=30"), NORTH_UP, "lon_0"),
        (GOOD.replace("lat_ts=-90", "lat_ts=-71"), NORTH_UP, "lat_ts"),
        (GOOD.replace("x_0=0", "x_0=100"), NORTH_UP, "x_0"),
        ("+proj=eqc +R=1737400 +units=m +no_defs", NORTH_UP, "projection"),
        (GOOD, rasterio.transform.from_origin(0.0, 20.0, 5.0, 4.0), "square"),
        (GOOD, rasterio.Affine(5.0, 0.5, 0.0, 0.0, -5.0, 20.0), "rotated"),
    ],
)
def test_unexpected_georeferencing_is_refused(tmp_path, proj4, transform, message):
    path = tmp_path / "bad.tif"
    _write(path, proj4=proj4, transform=transform)
    with pytest.raises(GeoreferencingError, match=message):
        read_dem(path)


def test_nodata_and_implausible_heights_are_refused(tmp_path):
    holes = np.zeros((4, 4), np.float32)
    holes[1, 1] = -9999.0
    path = tmp_path / "holes.tif"
    _write(path, proj4=GOOD, transform=NORTH_UP, data=holes, nodata=-9999.0)
    with pytest.raises(ValueError, match="no-data"):
        read_dem(path)
    values, _ = read_raster(path, allow_nodata=True)
    assert np.isnan(values[1, 1]) and np.isnan(values).sum() == 1

    path = tmp_path / "huge.tif"
    _write(path, proj4=GOOD, transform=NORTH_UP, data=np.full((4, 4), 5.0e4, np.float32))
    with pytest.raises(ValueError, match="plausible"):
        read_dem(path)
