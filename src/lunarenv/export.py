"""Neutral export of a run: boundary conditions for an asset's own thermal model.

Writes what a thermal analyst feeds to another tool (Thermal Desktop, ESATAN-TMS, a
script): the direction and visible fraction of the Sun over time, the terrain horizon,
the irradiance on each asset surface and, on request, the terrain as a triangle mesh.
Plain CSV and Wavefront OBJ, with a ``manifest.json`` that states the frame, the units
of every column, the checksum of every file and the validation status of the engine.

Local frame of every vector and mesh (the grid-aligned frame of
:meth:`lunarenv.geometry.projection.SouthPolarStereographic.local_basis`): origin on the
ground under the observer, ``+x`` along map ``+x``, ``+y`` along map ``+y``, ``+z`` up
(radially outward). Azimuths are grid azimuths, from ``+y`` towards ``+x``.

Ground temperatures are not exported: the polar thermal results are outside the
validated range (VALIDATION.md, TH-V3 to TH-V3g). They stay in ``results.nc``.

Assumptions: the horizon, and so the visible fraction of the Sun and the horizon ring,
hold at the observer point only; the mesh is the DEM surface at its own resolution, with
no relief below it.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike, NDArray

from lunarenv.data.dem import Bounds, Dem, load_dem
from lunarenv.data.store import DataStore
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.sites import TERRAIN_SETS

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

EXPORT_SCHEMA_VERSION = 1

#: Variables of ``results.nc`` written to the export, whatever the file.
_EXPORTED = (
    "horizon_elevation",
    "horizon_distance",
    "sun_elevation",
    "sun_azimuth",
    "sun_visible_fraction",
    "sun_flux_fraction",
    "eclipse_factor",
    "solar_irradiance_normal",
    "earth_elevation",
    "earth_visible_fraction",
    "oriented_direct_irradiance",
    "oriented_diffuse_irradiance",
)

_TIME_COLUMNS = (
    ("time_utc", "ISO 8601", "start of the time step, UTC"),
    ("elapsed_s", "s", "seconds since the first time step"),
)

# column, units, description, variable of results.nc (None: derived here)
_SUN_COLUMNS: tuple[tuple[str, str, str, str | None], ...] = (
    (
        "sun_azimuth_deg",
        "degree",
        "grid azimuth of the Sun's centre, from +y towards +x",
        "sun_azimuth",
    ),
    (
        "sun_elevation_deg",
        "degree",
        "apparent elevation of the Sun's centre above the local horizontal",
        "sun_elevation",
    ),
    ("sun_x", "1", "unit vector from the observer to the Sun's centre, local +x (map +x)", None),
    ("sun_y", "1", "unit vector from the observer to the Sun's centre, local +y (map +y)", None),
    ("sun_z", "1", "unit vector from the observer to the Sun's centre, local +z (up)", None),
    (
        "sun_visible_fraction",
        "1",
        "fraction of the solar disc area above the terrain horizon",
        "sun_visible_fraction",
    ),
    (
        "sun_flux_fraction",
        "1",
        "limb-darkening-weighted fraction of the solar disc above the terrain horizon",
        "sun_flux_fraction",
    ),
    (
        "eclipse_factor",
        "1",
        "fraction of the solar disc not covered by the Earth",
        "eclipse_factor",
    ),
    (
        "solar_irradiance_normal_w_m2",
        "W m-2",
        "direct solar irradiance on a surface facing the Sun: solar constant scaled to the Sun's distance, times sun_flux_fraction and eclipse_factor",
        "solar_irradiance_normal",
    ),
    (
        "earth_elevation_deg",
        "degree",
        "apparent elevation of the Earth's centre above the local horizontal",
        "earth_elevation",
    ),
    (
        "earth_visible_fraction",
        "1",
        "fraction of the Earth's disc above the terrain horizon",
        "earth_visible_fraction",
    ),
)

_HORIZON_COLUMNS = (
    ("azimuth_deg", "degree", "grid azimuth, from +y towards +x", "azimuth"),
    (
        "horizon_elevation_deg",
        "degree",
        "terrain horizon elevation seen from the observer height",
        "horizon_elevation",
    ),
    (
        "horizon_distance_m",
        "m",
        "ground distance to the terrain point that defines the horizon",
        "horizon_distance",
    ),
)

_DIFFUSE_DESCRIPTIONS = {
    "terrain_scattered": "sunlight scattered once by the terrain (Lambertian, scenario albedo law)",
    "terrain_infrared": "infrared emitted by the terrain (its own 1-D temperature)",
    "earthshine_reflected": "sunlight reflected by the Earth",
    "earthshine_infrared": "infrared emitted by the Earth",
}


# --------------------------------------------------------------------------- geometry
def direction_vectors(azimuth_deg: ArrayLike, elevation_deg: ArrayLike) -> FloatArray:
    """Unit vectors in the local frame, shape ``(..., 3)``, from grid azimuth and elevation."""
    azimuth = np.radians(np.asarray(azimuth_deg, dtype=np.float64))
    elevation = np.radians(np.asarray(elevation_deg, dtype=np.float64))
    cosine = np.cos(elevation)
    return np.stack(
        [cosine * np.sin(azimuth), cosine * np.cos(azimuth), np.sin(elevation)], axis=-1
    )


def terrain_mesh(
    dem: Dem, x_m: float, y_m: float, *, half_width_m: float, stride: int = 1
) -> tuple[FloatArray, IntArray]:
    """Triangle mesh of the DEM in a square around ``(x_m, y_m)``, in the local frame.

    One vertex per DEM pixel centre (every ``stride``-th pixel), two triangles per cell,
    all facing up. Vertices are true positions: the sphere's curvature is kept, so the
    mesh drops below the tangent plane by about ``d^2 / 2R`` at distance ``d``. The
    origin is the ground at ``(x_m, y_m)``. Returns ``(vertices (n, 3) in metres, faces
    (m, 3) as zero-based vertex indices)``. Raises if the square leaves the DEM.
    """
    if stride < 1:
        raise ValueError("stride must be at least 1")
    if not half_width_m > 0.0:
        raise ValueError("half_width_m must be positive")
    window = dem.crop(Bounds.around(x_m, y_m, half_width_m))
    heights = np.asarray(window.heights_m[::stride, ::stride], dtype=np.float64)
    xx, yy = np.meshgrid(window.grid.x_centers()[::stride], window.grid.y_centers()[::stride])
    rows, cols = heights.shape
    if rows < 2 or cols < 2:
        raise ValueError("the mesh needs at least two vertices in each direction")

    projection = dem.grid.projection
    origin = projection.to_cartesian(x_m, y_m, float(dem.sample(x_m, y_m)))
    e_x, e_y, up = projection.local_basis(x_m, y_m)
    offset = projection.to_cartesian(xx, yy, heights).reshape(-1, 3) - origin
    vertices = np.column_stack([offset @ e_x, offset @ e_y, offset @ up])

    # rows run towards -y, columns towards +x: (r, c) -> (r + 1, c) -> (r, c + 1) is
    # counter-clockwise seen from above
    index = np.arange(rows * cols, dtype=np.int64).reshape(rows, cols)
    a, b = index[:-1, :-1].ravel(), index[:-1, 1:].ravel()
    c, d = index[1:, :-1].ravel(), index[1:, 1:].ravel()
    faces = np.concatenate([np.column_stack([a, c, b]), np.column_stack([b, c, d])])
    return vertices, faces


def horizon_ring(
    azimuth_deg: ArrayLike,
    elevation_deg: ArrayLike,
    *,
    observer_height_m: float,
    radius_m: float,
) -> tuple[FloatArray, IntArray]:
    """A wall around the observer whose top edge is the terrain horizon.

    Seen from the observer, at ``(0, 0, observer_height_m)`` in the local frame, the top
    of the wall at each azimuth has the given elevation, so the wall hides the Sun
    exactly when the terrain does. It stands in for all terrain beyond ``radius_m``
    (far-field terrain reduced to a horizon, HLS-UG-001 section 9.4.3) and is exact at
    the observer point only. Vertices: the top edge in azimuth order, then the bottom
    edge; faces look inwards.
    """
    azimuth = np.radians(np.asarray(azimuth_deg, dtype=np.float64))
    elevation = np.radians(np.asarray(elevation_deg, dtype=np.float64))
    if azimuth.ndim != 1 or azimuth.shape != elevation.shape or azimuth.size < 3:
        raise ValueError("need one elevation per azimuth, at least three of them")
    if not radius_m > 0.0:
        raise ValueError("radius_m must be positive")
    n = azimuth.size
    x, y = radius_m * np.sin(azimuth), radius_m * np.cos(azimuth)
    top = observer_height_m + radius_m * np.tan(elevation)
    bottom = np.full(n, top.min() - max(0.05 * radius_m, 1.0))
    vertices = np.concatenate([np.column_stack([x, y, top]), np.column_stack([x, y, bottom])])

    i = np.arange(n, dtype=np.int64)
    j = (i + 1) % n
    faces = np.concatenate([np.column_stack([i, n + i, j]), np.column_stack([j, n + i, n + j])])
    return vertices, faces


# --------------------------------------------------------------------------- writers
def _number(value: float) -> str:
    return format(float(value), ".10g")


def _write_csv(path: Path, header: list[str], columns: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(zip(*columns, strict=True))


def write_obj(path: Path, vertices: FloatArray, faces: IntArray, comments: list[str]) -> None:
    """Wavefront OBJ: comment lines, ``v x y z`` in metres (0.1 mm), ``f`` one-based."""
    with path.open("w", encoding="utf-8") as handle:
        handle.writelines(f"# {line}\n" for line in comments)
        handle.writelines(f"v {x:.4f} {y:.4f} {z:.4f}\n" for x, y, z in vertices)
        handle.writelines(f"f {a} {b} {c}\n" for a, b, c in faces + 1)


def _file_entry(path: Path, **fields: Any) -> dict[str, Any]:
    content = path.read_bytes()
    return {
        "name": path.name,
        "sha256": hashlib.sha256(content).hexdigest(),
        "bytes": len(content),
        **fields,
    }


def _columns(
    spec: tuple[tuple[str, str, str], ...] | list[tuple[str, str, str]],
) -> list[dict[str, str]]:
    return [{"name": name, "units": units, "description": text} for name, units, text in spec]


def _time_columns(dataset: xr.Dataset) -> tuple[list[str], list[str], dict[str, Any]]:
    time = dataset["time"].values.astype("datetime64[s]")
    elapsed = (time - time[0]).astype(np.float64)
    stamps = [f"{stamp}Z" for stamp in np.datetime_as_string(time, unit="s")]
    step = float(elapsed[1] - elapsed[0]) if time.size > 1 else 0.0
    summary = {
        "scale": "UTC",
        "start_utc": stamps[0],
        "end_utc": stamps[-1],
        "step_s": step,
        "steps": int(time.size),
    }
    return stamps, [_number(v) for v in elapsed], summary


# --------------------------------------------------------------------------- export
def export_run(
    run: str | Path,
    output: str | Path | None = None,
    *,
    mesh_half_width_m: float | None = None,
    mesh_stride: int = 1,
    ring_radius_m: float | None = None,
    store: DataStore | None = None,
) -> dict[str, Any]:
    """Write the neutral export of a run and return its manifest.

    ``run`` is a run directory (``results.nc`` and ``metadata.json``); ``output``
    defaults to ``<run>/export``. Always written: ``sun_earth.csv``, ``horizon.csv``,
    one ``surface_<name>.csv`` per asset surface of the scenario, ``manifest.json``.
    With ``ring_radius_m``: ``horizon_ring.obj`` (see :func:`horizon_ring`). With
    ``mesh_half_width_m``: ``terrain.obj``, the site DEM in a square of that half width
    (see :func:`terrain_mesh`); this reads the site terrain from the data cache.

    The same run and options give byte-identical files. The returned manifest has one
    key more than ``manifest.json``: ``directory``.
    """
    from lunarenv.runner import load_run

    run_dir = Path(run)
    if not (run_dir / "results.nc").is_file() or not (run_dir / "metadata.json").is_file():
        raise FileNotFoundError(f"{run_dir} is not a run directory (results.nc, metadata.json)")
    dataset, metadata = load_run(run_dir)
    out = Path(output) if output is not None else run_dir / "export"
    out.mkdir(parents=True, exist_ok=True)

    site = metadata["site"]
    x_m, y_m = float(site["x_m"]), float(site["y_m"])
    observer_height = float(metadata["scenario"]["site"]["height_above_surface_m"])
    stamps, elapsed, time_summary = _time_columns(dataset)
    files: list[dict[str, Any]] = []

    # ------------------------------------------------------------------ Sun and Earth
    vectors = direction_vectors(dataset["sun_azimuth"].values, dataset["sun_elevation"].values)
    derived = {"sun_x": vectors[:, 0], "sun_y": vectors[:, 1], "sun_z": vectors[:, 2]}
    spec = [
        (name, units, text, variable)
        for name, units, text, variable in _SUN_COLUMNS
        if variable is None or variable in dataset
    ]
    values = [
        derived[name] if variable is None else dataset[variable].values
        for name, _, _, variable in spec
    ]
    path = out / "sun_earth.csv"
    _write_csv(
        path,
        [c[0] for c in _TIME_COLUMNS] + [c[0] for c in spec],
        [stamps, elapsed] + [[_number(v) for v in column] for column in values],
    )
    files.append(
        _file_entry(
            path,
            rows=len(stamps),
            columns=_columns(list(_TIME_COLUMNS) + [c[:3] for c in spec]),
        )
    )

    # ------------------------------------------------------------------------ horizon
    path = out / "horizon.csv"
    _write_csv(
        path,
        [c[0] for c in _HORIZON_COLUMNS],
        [[_number(v) for v in dataset[variable].values] for *_, variable in _HORIZON_COLUMNS],
    )
    files.append(
        _file_entry(
            path,
            rows=int(dataset.sizes["azimuth"]),
            columns=_columns([c[:3] for c in _HORIZON_COLUMNS]),
        )
    )

    # ------------------------------------------------------------------ asset surfaces
    surfaces: dict[str, Any] = {}
    if "oriented_direct_irradiance" in dataset:
        described = metadata["summary"].get("oriented_surfaces", {}).get("surfaces", {})
        diffuse = dataset.get("oriented_diffuse_irradiance")
        components = [] if diffuse is None else [str(c) for c in dataset["component"].values]
        for index, name in enumerate(str(s) for s in dataset["surface"].values):
            direct = dataset["oriented_direct_irradiance"].values[index]
            parts = [] if diffuse is None else list(diffuse.values[index])
            spec3 = [
                (
                    "direct_w_m2",
                    "W m-2",
                    "direct solar irradiance on the surface: solar_irradiance_normal times the cosine of the incidence, zero from behind",
                )
            ]
            spec3 += [
                (
                    f"{c}_w_m2",
                    "W m-2",
                    f"irradiance on the surface: {_DIFFUSE_DESCRIPTIONS.get(c, c)}",
                )
                for c in components
            ]
            spec3.append(("total_w_m2", "W m-2", "sum of the irradiance columns of this file"))
            columns = [direct, *parts, direct + np.sum(parts, axis=0) if parts else direct]
            path = out / f"surface_{name}.csv"
            _write_csv(
                path,
                [c[0] for c in _TIME_COLUMNS] + [c[0] for c in spec3],
                [stamps, elapsed] + [[_number(v) for v in column] for column in columns],
            )
            files.append(
                _file_entry(path, rows=len(stamps), columns=_columns(list(_TIME_COLUMNS) + spec3))
            )
            info = described.get(name, {})
            entry: dict[str, Any] = {"file": path.name}
            if "tilt_deg" in info:
                tilt, facing = float(info["tilt_deg"]), float(info["facing_grid_azimuth_deg"])
                entry |= {
                    "tilt_deg": tilt,
                    "facing_grid_azimuth_deg": facing,
                    "normal_local": direction_vectors(facing, 90.0 - tilt).tolist(),
                }
            entry |= {k: info[k] for k in ("sky_view_factor", "ground_view_factor") if k in info}
            surfaces[name] = entry

    # ------------------------------------------------------------------------- meshes
    if ring_radius_m is not None:
        vertices, faces = horizon_ring(
            dataset["azimuth"].values,
            dataset["horizon_elevation"].values,
            observer_height_m=observer_height,
            radius_m=ring_radius_m,
        )
        path = out / "horizon_ring.obj"
        write_obj(
            path,
            vertices,
            faces,
            [
                f"lunarenv run {metadata['run_id']}: terrain horizon as a wall of radius {_number(ring_radius_m)} m",
                f"exact for an observer at (0, 0, {_number(observer_height)}) m only; metres, local frame (manifest.json)",
            ],
        )
        files.append(
            _file_entry(
                path,
                vertices=len(vertices),
                faces=len(faces),
                units="m",
                radius_m=float(ring_radius_m),
                description="wall whose top edge is the terrain horizon seen from the observer; stands in for the terrain beyond it",
            )
        )
    if mesh_half_width_m is not None:
        dem = load_dem(store or DataStore(), TERRAIN_SETS[site["terrain_set"]].near_id)
        vertices, faces = terrain_mesh(
            dem, x_m, y_m, half_width_m=mesh_half_width_m, stride=mesh_stride
        )
        spacing = dem.grid.pixel_size_m * mesh_stride
        path = out / "terrain.obj"
        write_obj(
            path,
            vertices,
            faces,
            [
                f"lunarenv run {metadata['run_id']}: site terrain {dem.source_id}, +/-{_number(mesh_half_width_m)} m, vertex spacing {_number(spacing)} m",
                "metres, local frame (manifest.json); origin on the ground under the observer",
            ],
        )
        files.append(
            _file_entry(
                path,
                vertices=len(vertices),
                faces=len(faces),
                units="m",
                half_width_m=float(mesh_half_width_m),
                vertex_spacing_m=float(spacing),
                source_id=dem.source_id,
                source_sha256=dem.source_sha256,
                description="site DEM surface, one vertex per DEM pixel centre, faces up; no relief below the DEM resolution. vertex_spacing_m is in map metres; true distances are larger by (R + height) / R, about 0.1 % at 1.7 km",
            )
        )

    # ----------------------------------------------------------------------- manifest
    projection = SouthPolarStereographic(float(site["reference_radius_m"]))
    e_x, e_y, up = projection.local_basis(x_m, y_m)
    validation = metadata.get("validation") or {}
    manifest: dict[str, Any] = {
        "export_schema_version": EXPORT_SCHEMA_VERSION,
        "run_id": metadata["run_id"],
        "result_checksum_sha256": metadata["result_checksum_sha256"],
        "lunarenv_version": metadata["lunarenv_version"],
        "scenario_name": metadata["scenario"].get("name"),
        "site": site,
        "frame": {
            "description": "Grid-aligned local frame. Origin on the ground under the observer; +x along map +x and +y along map +y of the south polar stereographic projection of the terrain products; +z up (radially outward). Right-handed. Azimuths are grid azimuths, from +y towards +x; true azimuth = grid azimuth - longitude.",
            "observer_local_m": [0.0, 0.0, observer_height],
            "body_fixed_frame": site["body_fixed_frame"],
            "local_to_body_fixed_rows": [e_x.tolist(), e_y.tolist(), up.tolist()],
            "local_to_body_fixed_note": "rows are the local +x, +y, +z unit vectors in body-fixed axes: body = x * row0 + y * row1 + z * row2",
        },
        "time": time_summary,
        "files": files,
        "surfaces": surfaces,
        "validity_warnings": list(metadata.get("validity", {}).get("warnings", [])),
        "validation": {
            "generated_utc": validation.get("generated_utc"),
            "counts": validation.get("counts"),
            "note": "Status of the engine's validation cases when the run was made (VALIDATION.md). The illumination cases bear on these files; the horizon and the visible fraction hold at the observer point and height.",
        },
        "not_exported": {
            "variables": sorted(str(v) for v in dataset.data_vars if v not in _EXPORTED),
            "reason": "Ground temperatures and their bounds stay in results.nc: the polar thermal results are outside the validated range. Use them as bounds for a ground boundary condition, not as predictions.",
        },
        "data": [
            {k: record.get(k) for k in ("id", "url", "version", "sha256", "citation")}
            for record in metadata.get("data", [])
        ],
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {**manifest, "directory": str(out)}


__all__ = [
    "EXPORT_SCHEMA_VERSION",
    "direction_vectors",
    "export_run",
    "horizon_ring",
    "terrain_mesh",
    "write_obj",
]
