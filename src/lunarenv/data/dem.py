"""Digital elevation models: strict reading, cropping and sampling.

A DEM is only useful to the physics if its georeferencing is exactly what the code
assumes. :func:`read_dem` therefore refuses any raster that is not a north-up, square
pixel, south polar stereographic grid on the expected sphere, instead of reprojecting or
guessing.

Grid convention: GDAL "pixel is area". ``x_min``/``y_max`` are the outer edges of the
top-left pixel; pixel centres sit half a pixel inside. Row 0 is the northernmost
(largest ``y``) row. This matches the "pixel-registered" statement in the PGDA README.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from numpy.typing import ArrayLike, NDArray
from rasterio.windows import Window

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M
from lunarenv.data.store import DataStore, DerivedRecord
from lunarenv.geometry.projection import FloatArray, SouthPolarStereographic

# Heights outside this range cannot be lunar topography relative to the 1737.4 km sphere
# (the global extremes are about -9.1 km and +10.8 km); they indicate a unit mix-up.
_PLAUSIBLE_HEIGHT_RANGE_M = (-10_000.0, 11_000.0)
_GEOREF_TOLERANCE = 1e-9


class GeoreferencingError(ValueError):
    """A raster's georeferencing is not what lunarenv requires."""


@dataclass(frozen=True, slots=True)
class Bounds:
    """Axis-aligned rectangle in map coordinates (metres)."""

    x_min: float
    y_min: float
    x_max: float
    y_max: float

    def __post_init__(self) -> None:
        if not (self.x_max > self.x_min and self.y_max > self.y_min):
            raise ValueError("bounds must have positive width and height")

    @classmethod
    def around(cls, x_m: float, y_m: float, half_width_m: float) -> Bounds:
        """Square of side ``2 * half_width_m`` centred on a point."""
        return cls(x_m - half_width_m, y_m - half_width_m, x_m + half_width_m, y_m + half_width_m)

    def contains(self, other: Bounds) -> bool:
        return (
            self.x_min <= other.x_min
            and self.y_min <= other.y_min
            and self.x_max >= other.x_max
            and self.y_max >= other.y_max
        )


@dataclass(frozen=True, slots=True)
class GridGeometry:
    """Regular north-up grid in a south polar stereographic projection."""

    x_min: float
    y_max: float
    pixel_size_m: float
    width: int
    height: int
    projection: SouthPolarStereographic

    def __post_init__(self) -> None:
        if self.pixel_size_m <= 0.0 or self.width <= 0 or self.height <= 0:
            raise ValueError("grid must have positive pixel size and shape")

    @property
    def shape(self) -> tuple[int, int]:
        return (self.height, self.width)

    @property
    def bounds(self) -> Bounds:
        return Bounds(
            self.x_min,
            self.y_max - self.height * self.pixel_size_m,
            self.x_min + self.width * self.pixel_size_m,
            self.y_max,
        )

    def x_centers(self) -> FloatArray:
        return self.x_min + (np.arange(self.width, dtype=np.float64) + 0.5) * self.pixel_size_m

    def y_centers(self) -> FloatArray:
        """Pixel-centre ``y`` per row; decreasing, because row 0 is the northernmost."""
        return self.y_max - (np.arange(self.height, dtype=np.float64) + 0.5) * self.pixel_size_m

    def fractional_index(self, x_m: ArrayLike, y_m: ArrayLike) -> tuple[FloatArray, FloatArray]:
        """Continuous ``(row, col)`` such that integer values fall on pixel centres."""
        col = (np.asarray(x_m, dtype=np.float64) - self.x_min) / self.pixel_size_m - 0.5
        row = (self.y_max - np.asarray(y_m, dtype=np.float64)) / self.pixel_size_m - 0.5
        return row, col

    def window(self, bounds: Bounds) -> tuple[int, int, int, int]:
        """Smallest block of whole pixels covering ``bounds``: ``(row0, row1, col0, col1)``.

        Half-open ranges. Raises if ``bounds`` is not fully inside the grid; silently
        clipping would hide a terrain model that is too small for the request.
        """
        if not self.bounds.contains(bounds):
            raise ValueError(f"requested {bounds} extends beyond the grid {self.bounds}")
        size = self.pixel_size_m
        col0 = math.floor((bounds.x_min - self.x_min) / size + _GEOREF_TOLERANCE)
        col1 = math.ceil((bounds.x_max - self.x_min) / size - _GEOREF_TOLERANCE)
        row0 = math.floor((self.y_max - bounds.y_max) / size + _GEOREF_TOLERANCE)
        row1 = math.ceil((self.y_max - bounds.y_min) / size - _GEOREF_TOLERANCE)
        return row0, row1, col0, col1

    def subgrid(self, row0: int, row1: int, col0: int, col1: int) -> GridGeometry:
        if not (0 <= row0 < row1 <= self.height and 0 <= col0 < col1 <= self.width):
            raise ValueError("window is empty or outside the grid")
        return GridGeometry(
            x_min=self.x_min + col0 * self.pixel_size_m,
            y_max=self.y_max - row0 * self.pixel_size_m,
            pixel_size_m=self.pixel_size_m,
            width=col1 - col0,
            height=row1 - row0,
            projection=self.projection,
        )


@dataclass(frozen=True, slots=True)
class Dem:
    """Heights above the reference sphere on a regular grid.

    ``heights_m`` has shape ``(height, width)``; it is stored as read (float32 in the
    LOLA products, i.e. about 0.1 mm resolution at kilometre heights).
    """

    heights_m: NDArray[np.float32]
    grid: GridGeometry
    source_id: str = ""
    source_sha256: str = ""

    def __post_init__(self) -> None:
        if self.heights_m.shape != self.grid.shape:
            raise ValueError("height array shape does not match the grid")

    def crop(self, bounds: Bounds) -> Dem:
        """Sub-DEM of whole pixels covering ``bounds`` (a view, not a copy)."""
        row0, row1, col0, col1 = self.grid.window(bounds)
        return Dem(
            self.heights_m[row0:row1, col0:col1],
            self.grid.subgrid(row0, row1, col0, col1),
            self.source_id,
            self.source_sha256,
        )

    def sample(self, x_m: ArrayLike, y_m: ArrayLike) -> FloatArray:
        """Bilinearly interpolated height at map coordinates.

        The interpolant is defined on the rectangle spanned by pixel centres. Points
        outside it raise: extrapolating terrain is never the right answer here.
        """
        row, col = self.grid.fractional_index(x_m, y_m)
        row_max, col_max = self.grid.height - 1, self.grid.width - 1
        slack = 1e-9  # absorbs round-off for points exactly on the outer centres
        if (
            np.any(row < -slack)
            or np.any(col < -slack)
            or np.any(row > row_max + slack)
            or np.any(col > col_max + slack)
        ):
            raise ValueError("sample point outside the DEM interpolation domain")
        row = np.clip(row, 0.0, row_max)
        col = np.clip(col, 0.0, col_max)
        r0 = np.minimum(np.floor(row).astype(np.intp), max(row_max - 1, 0))
        c0 = np.minimum(np.floor(col).astype(np.intp), max(col_max - 1, 0))
        r1 = np.minimum(r0 + 1, row_max)
        c1 = np.minimum(c0 + 1, col_max)
        fr, fc = row - r0, col - c0
        z = self.heights_m
        top = z[r0, c0] * (1.0 - fc) + z[r0, c1] * fc
        bottom = z[r1, c0] * (1.0 - fc) + z[r1, c1] * fc
        return np.asarray(top * (1.0 - fr) + bottom * fr, dtype=np.float64)


# --------------------------------------------------------------------------- reading
def _sphere_radius(crs: dict[str, Any]) -> float:
    if "R" in crs:
        return float(crs["R"])
    if "a" in crs and ("b" in crs or "rf" in crs):
        a = float(crs["a"])
        if ("b" in crs and float(crs["b"]) == a) or ("rf" in crs and float(crs["rf"]) == 0.0):
            return a
    raise GeoreferencingError(f"cannot establish a spherical datum from CRS {crs}")


def grid_from_dataset(
    dataset: rasterio.DatasetReader, *, expected_radius_m: float = MOON_REFERENCE_RADIUS_M
) -> GridGeometry:
    """Check a raster's georeferencing and return its grid. Raises on any deviation."""
    if dataset.crs is None:
        raise GeoreferencingError("raster has no coordinate reference system")
    crs = dataset.crs.to_dict()
    problems: list[str] = []
    if crs.get("proj") != "stere":
        problems.append(f"projection is {crs.get('proj')!r}, expected 'stere'")
    for key, expected in (("lat_0", -90.0), ("lon_0", 0.0), ("x_0", 0.0), ("y_0", 0.0)):
        if float(crs.get(key, 0.0 if key != "lat_0" else math.nan)) != expected:
            problems.append(f"{key} is {crs.get(key)!r}, expected {expected}")
    if float(crs.get("lat_ts", -90.0)) != -90.0:
        problems.append(f"lat_ts is {crs.get('lat_ts')!r}, expected -90 (true scale at the pole)")
    for key in ("k", "k_0"):
        if float(crs.get(key, 1.0)) != 1.0:
            problems.append(f"{key} is {crs.get(key)!r}, expected 1")
    if crs.get("units", "m") != "m":
        problems.append(f"units are {crs.get('units')!r}, expected metres")
    try:
        radius = _sphere_radius(crs)
    except GeoreferencingError as error:
        problems.append(str(error))
    else:
        if radius != expected_radius_m:
            problems.append(f"sphere radius is {radius} m, expected {expected_radius_m} m")

    t = dataset.transform
    if t.b != 0.0 or t.d != 0.0:
        problems.append("grid is rotated or sheared")
    if not (t.a > 0.0 and t.e < 0.0 and abs(t.a + t.e) <= _GEOREF_TOLERANCE * t.a):
        problems.append(f"pixels are not square and north-up (dx={t.a}, dy={t.e})")
    if dataset.count != 1:
        problems.append(f"raster has {dataset.count} bands, expected 1")
    if problems:
        raise GeoreferencingError(f"{dataset.name}: " + "; ".join(problems))

    return GridGeometry(
        x_min=float(t.c),
        y_max=float(t.f),
        pixel_size_m=float(t.a),
        width=int(dataset.width),
        height=int(dataset.height),
        projection=SouthPolarStereographic(expected_radius_m),
    )


def read_raster(
    path: Path,
    bounds: Bounds | None = None,
    *,
    allow_nodata: bool = False,
) -> tuple[NDArray[np.float32], GridGeometry]:
    """Read a single-band raster (or the pixels covering ``bounds``) with strict checks.

    Used for DEMs and for co-registered auxiliary layers (height error, slope error,
    published illumination maps).
    """
    with rasterio.open(path) as dataset:
        full = grid_from_dataset(dataset)
        if bounds is None:
            row0, row1, col0, col1 = 0, full.height, 0, full.width
        else:
            row0, row1, col0, col1 = full.window(bounds)
        data = dataset.read(1, window=Window(col0, row0, col1 - col0, row1 - row0))
        nodata = dataset.nodata
        scale, offset = float(dataset.scales[0]), float(dataset.offsets[0])
    values = np.asarray(data, dtype=np.float32)
    if scale != 1.0 or offset != 0.0:  # packed integers (e.g. the PGDA illumination maps)
        values = (np.asarray(data, dtype=np.float64) * scale + offset).astype(np.float32)
    if nodata is not None and not math.isnan(nodata):
        values = np.where(data == nodata, np.float32(np.nan), values)
    if not allow_nodata and bool(np.isnan(values).any()):
        raise ValueError(f"{path.name}: window contains no-data pixels")
    return values, full.subgrid(row0, row1, col0, col1)


def read_dem(
    path: Path,
    bounds: Bounds | None = None,
    *,
    source_id: str = "",
    source_sha256: str = "",
) -> Dem:
    """Read a DEM in metres above the reference sphere.

    Beyond the georeferencing checks, heights must be free of no-data and inside the
    range of lunar topography. That bound catches gross errors only; a DEM stored in
    kilometres would pass it, so the unit of every catalog DEM is additionally pinned by
    a regression test on its height statistics.
    """
    heights, grid = read_raster(path, bounds)
    low, high = float(heights.min()), float(heights.max())
    if low < _PLAUSIBLE_HEIGHT_RANGE_M[0] or high > _PLAUSIBLE_HEIGHT_RANGE_M[1]:
        raise ValueError(
            f"{path.name}: heights span [{low}, {high}] m, outside the plausible lunar range; "
            "check the units of the source"
        )
    return Dem(heights, grid, source_id, source_sha256)


def load_dem(store: DataStore, entry_id: str, bounds: Bounds | None = None) -> Dem:
    """Read a catalog DEM from the data cache, carrying its checksum as provenance."""
    record = store.record(entry_id)
    return read_dem(store.path(entry_id), bounds, source_id=record.id, source_sha256=record.sha256)


# --------------------------------------------------------------------------- writing
def write_dem(dem: Dem, path: Path) -> None:
    """Write a DEM as a compressed GeoTIFF with the georeferencing ``read_dem`` expects."""
    grid = dem.grid
    crs = rasterio.crs.CRS.from_proj4(
        f"+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +x_0=0 +y_0=0 "
        f"+R={grid.projection.radius_m:.3f} +units=m +no_defs"
    )
    transform = rasterio.transform.from_origin(
        grid.x_min, grid.y_max, grid.pixel_size_m, grid.pixel_size_m
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=grid.width,
        height=grid.height,
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=float("nan"),
        compress="deflate",
        predictor=3,
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dataset:
        dataset.write(np.ascontiguousarray(dem.heights_m, dtype=np.float32), 1)


def ingest_crop(
    store: DataStore, entry_id: str, bounds: Bounds, *, keep_raw: bool = True
) -> DerivedRecord:
    """Crop a cached DEM, store the crop as a derived product and optionally drop the source.

    This is how large ensembles (the 100 error clones per site) fit on a laptop: each
    file is downloaded, verified against its pin, cropped to the window that the
    near-field horizon needs, and deleted. The derived record keeps the parent checksum,
    so provenance survives the deletion.
    """
    dem = load_dem(store, entry_id, bounds)
    b = dem.grid.bounds
    tag = f"x{b.x_min:.0f}_{b.x_max:.0f}_y{b.y_min:.0f}_{b.y_max:.0f}"
    derived_id = f"{entry_id}.crop.{tag}"
    relative = store.catalog.get(entry_id).relative_path
    out = store.new_derived_path(str(relative.with_name(f"{relative.stem}__{tag}.tif")))
    write_dem(dem, out)
    record = store.register_derived(
        derived_id,
        out,
        parent_id=entry_id,
        operation="crop",
        parameters={"x_min": b.x_min, "y_min": b.y_min, "x_max": b.x_max, "y_max": b.y_max},
        created_utc=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )
    if not keep_raw:
        store.discard_raw(entry_id)
    return record
