"""Neutral export of a run: tables, meshes and manifest, on a synthetic run."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from lunarenv.cli import main as cli_main
from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.export import (
    direction_vectors,
    export_run,
    horizon_ring,
    terrain_mesh,
    write_obj,
)
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack
from lunarenv.thermal.surface import plane_normal
from tests.helpers import PROJ, analytic_dem

SITE_X, SITE_Y, HEIGHT = -2800.0, -9000.0, 2.0
COMPONENTS = (
    "terrain_scattered",
    "terrain_infrared",
    "earthshine_reflected",
    "earthshine_infrared",
)
SURFACES = {"roof": (0.0, 0.0), "wall-90": (90.0, 90.0)}


def _write_run(directory: Path, *, thermal: bool = True) -> xr.Dataset:
    rng = np.random.default_rng(7)
    n_time, n_az = 48, 72
    time = np.datetime64("2026-01-01T00:00:00") + np.arange(n_time) * np.timedelta64(30, "m")
    azimuth = np.arange(n_az) * 360.0 / n_az
    variables: dict[str, tuple] = {
        "horizon_elevation": (("azimuth",), rng.uniform(-1.0, 3.0, n_az), {"units": "degree"}),
        "horizon_distance": (("azimuth",), rng.uniform(50.0, 9e4, n_az), {"units": "m"}),
        "sun_elevation": (("time",), rng.uniform(-1.5, 1.8, n_time), {"units": "degree"}),
        "sun_azimuth": (("time",), rng.uniform(0.0, 360.0, n_time), {"units": "degree"}),
        "sun_visible_fraction": (("time",), rng.uniform(0.0, 1.0, n_time), {"units": "1"}),
        "sun_flux_fraction": (("time",), rng.uniform(0.0, 1.0, n_time), {"units": "1"}),
        "eclipse_factor": (("time",), np.ones(n_time), {"units": "1"}),
        "solar_irradiance_normal": (
            ("time",),
            rng.uniform(0.0, 1400.0, n_time),
            {"units": "W m-2"},
        ),
        "earth_elevation": (("time",), rng.uniform(-6.0, 6.0, n_time), {"units": "degree"}),
        "earth_visible_fraction": (("time",), rng.uniform(0.0, 1.0, n_time), {"units": "1"}),
    }
    coords: dict[str, tuple] = {"time": (("time",), time), "azimuth": (("azimuth",), azimuth)}
    summary: dict = {"illumination": {}}
    if thermal:
        coords["surface"] = (("surface",), np.array(list(SURFACES)))
        coords["component"] = (("component",), np.array(COMPONENTS))
        variables["oriented_direct_irradiance"] = (
            ("surface", "time"),
            rng.uniform(0.0, 1400.0, (len(SURFACES), n_time)),
            {"units": "W m-2"},
        )
        variables["oriented_diffuse_irradiance"] = (
            ("surface", "component", "time"),
            rng.uniform(0.0, 6.0, (len(SURFACES), len(COMPONENTS), n_time)),
            {"units": "W m-2"},
        )
        variables["surface_temperature"] = (
            ("time",),
            rng.uniform(40.0, 120.0, n_time),
            {"units": "K"},
        )
        summary["oriented_surfaces"] = {
            "surfaces": {
                name: {
                    "tilt_deg": tilt,
                    "facing_grid_azimuth_deg": facing,
                    "sky_view_factor": 0.75,
                    "ground_view_factor": 0.25,
                }
                for name, (tilt, facing) in SURFACES.items()
            }
        }
    dataset = xr.Dataset(
        {k: xr.Variable(d, v, a) for k, (d, v, a) in variables.items()},
        coords={k: xr.Variable(d, v) for k, (d, v) in coords.items()},
    )
    metadata = {
        "run_id": "0123456789abcdef",
        "result_checksum_sha256": "ab" * 32,
        "lunarenv_version": "0.0.test",
        "scenario": {
            "name": "synthetic",
            "site": {"terrain": "site04", "height_above_surface_m": HEIGHT},
        },
        "site": {
            "x_m": SITE_X,
            "y_m": SITE_Y,
            "longitude_deg": 197.3,
            "latitude_deg": -89.69,
            "terrain_height_m": 1730.5,
            "terrain_set": "site04",
            "body_fixed_frame": "MOON_ME_DE440_ME421",
            "reference_radius_m": R,
        },
        "summary": summary,
        "validity": {"warnings": ["thermal: outside the validated range"]},
        "validation": {
            "generated_utc": "2026-11-16T21:32:48Z",
            "counts": {"PASS": 19, "FAIL": 6},
            "cases": [{"id": "ILL-A1", "status": "PASS", "summary": "ok"}],
        },
        "data": [
            {
                "id": "pgda78.site04.surf",
                "url": "https://example.invalid/site04.tif",
                "version": "v1",
                "sha256": "cd" * 32,
                "citation": "Barker2021",
                "downloaded_utc": "2026-11-12T15:30:36Z",
            }
        ],
    }
    directory.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(directory / "results.nc", engine="h5netcdf")
    (directory / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    return dataset


def _read(path: Path) -> dict[str, list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    return {name: [row[i] for row in rows[1:]] for i, name in enumerate(rows[0])}


def _floats(table: dict[str, list[str]], name: str) -> np.ndarray:
    return np.array([float(v) for v in table[name]])


def test_direction_vectors_follow_the_grid_azimuth_convention():
    vectors = direction_vectors([0.0, 90.0, 180.0, 33.0], [0.0, 0.0, 0.0, 90.0])
    np.testing.assert_allclose(vectors[0], [0.0, 1.0, 0.0], atol=1e-15)  # azimuth 0: map +y
    np.testing.assert_allclose(vectors[1], [1.0, 0.0, 0.0], atol=1e-15)  # azimuth 90: map +x
    np.testing.assert_allclose(vectors[2], [0.0, -1.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(vectors[3], [0.0, 0.0, 1.0], atol=1e-15)
    # the same convention as the asset surface normals of the engine
    np.testing.assert_allclose(
        direction_vectors([40.0], [0.0])[0], plane_normal(np.pi / 2, np.radians(40.0)), atol=1e-15
    )


def test_tables_round_trip(tmp_path: Path):
    dataset = _write_run(tmp_path / "run")
    manifest = export_run(tmp_path / "run")
    out = tmp_path / "run" / "export"
    assert manifest["directory"] == str(out)

    sun = _read(out / "sun_earth.csv")
    assert sun["time_utc"][0] == "2026-01-01T00:00:00Z"
    assert sun["time_utc"][1] == "2026-01-01T00:30:00Z"
    np.testing.assert_array_equal(_floats(sun, "elapsed_s"), np.arange(48) * 1800.0)
    for column, variable in (
        ("sun_azimuth_deg", "sun_azimuth"),
        ("sun_elevation_deg", "sun_elevation"),
        ("sun_visible_fraction", "sun_visible_fraction"),
        ("sun_flux_fraction", "sun_flux_fraction"),
        ("eclipse_factor", "eclipse_factor"),
        ("solar_irradiance_normal_w_m2", "solar_irradiance_normal"),
        ("earth_elevation_deg", "earth_elevation"),
        ("earth_visible_fraction", "earth_visible_fraction"),
    ):
        np.testing.assert_allclose(_floats(sun, column), dataset[variable].values, rtol=1e-9)
    vectors = np.column_stack([_floats(sun, f"sun_{axis}") for axis in "xyz"])
    np.testing.assert_allclose(np.linalg.norm(vectors, axis=1), 1.0, rtol=1e-9)
    np.testing.assert_allclose(
        np.degrees(np.arcsin(vectors[:, 2])), dataset["sun_elevation"].values, atol=1e-7
    )
    np.testing.assert_allclose(
        np.degrees(np.arctan2(vectors[:, 0], vectors[:, 1])) % 360.0,
        dataset["sun_azimuth"].values,
        atol=1e-5,
    )

    horizon = _read(out / "horizon.csv")
    np.testing.assert_allclose(_floats(horizon, "azimuth_deg"), dataset["azimuth"].values)
    np.testing.assert_allclose(
        _floats(horizon, "horizon_elevation_deg"), dataset["horizon_elevation"].values, rtol=1e-9
    )
    np.testing.assert_allclose(
        _floats(horizon, "horizon_distance_m"), dataset["horizon_distance"].values, rtol=1e-9
    )

    for index, name in enumerate(SURFACES):
        table = _read(out / f"surface_{name}.csv")
        direct = dataset["oriented_direct_irradiance"].values[index]
        diffuse = dataset["oriented_diffuse_irradiance"].values[index]
        np.testing.assert_allclose(_floats(table, "direct_w_m2"), direct, rtol=1e-9)
        for k, component in enumerate(COMPONENTS):
            np.testing.assert_allclose(_floats(table, f"{component}_w_m2"), diffuse[k], rtol=1e-9)
        np.testing.assert_allclose(
            _floats(table, "total_w_m2"), direct + diffuse.sum(axis=0), rtol=1e-9
        )


def test_manifest_describes_every_file(tmp_path: Path):
    _write_run(tmp_path / "run")
    manifest = export_run(tmp_path / "run", output=tmp_path / "elsewhere")
    out = tmp_path / "elsewhere"
    on_disk = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert on_disk == {k: v for k, v in manifest.items() if k != "directory"}

    assert on_disk["run_id"] == "0123456789abcdef"
    assert on_disk["result_checksum_sha256"] == "ab" * 32
    names = {entry["name"] for entry in on_disk["files"]}
    assert names == {"sun_earth.csv", "horizon.csv", "surface_roof.csv", "surface_wall-90.csv"}
    assert names | {"manifest.json"} == {p.name for p in out.iterdir()}
    for entry in on_disk["files"]:
        content = (out / entry["name"]).read_bytes()
        assert entry["sha256"] == hashlib.sha256(content).hexdigest()
        assert entry["bytes"] == len(content)
        header = content.decode("utf-8").splitlines()[0].split(",")
        assert [c["name"] for c in entry["columns"]] == header
        assert all(c["units"] and c["description"] for c in entry["columns"])
        assert entry["rows"] == len(content.decode("utf-8").splitlines()) - 1

    frame = on_disk["frame"]
    e_x, e_y, up = PROJ.local_basis(SITE_X, SITE_Y)
    np.testing.assert_allclose(frame["local_to_body_fixed_rows"], [e_x, e_y, up], atol=1e-15)
    assert frame["body_fixed_frame"] == "MOON_ME_DE440_ME421"
    assert frame["observer_local_m"] == [0.0, 0.0, HEIGHT]

    wall = on_disk["surfaces"]["wall-90"]
    np.testing.assert_allclose(wall["normal_local"], [1.0, 0.0, 0.0], atol=1e-15)
    assert wall["sky_view_factor"] == 0.75
    assert wall["file"] == "surface_wall-90.csv"

    assert on_disk["validation"]["counts"] == {"PASS": 19, "FAIL": 6}
    assert on_disk["validity_warnings"] == ["thermal: outside the validated range"]
    assert "surface_temperature" in on_disk["not_exported"]["variables"]
    assert on_disk["data"] == [
        {
            "id": "pgda78.site04.surf",
            "url": "https://example.invalid/site04.tif",
            "version": "v1",
            "sha256": "cd" * 32,
            "citation": "Barker2021",
        }
    ]


def test_export_is_reproducible(tmp_path: Path):
    _write_run(tmp_path / "run")
    export_run(tmp_path / "run", output=tmp_path / "a")
    export_run(tmp_path / "run", output=tmp_path / "b")
    names = sorted(p.name for p in (tmp_path / "a").iterdir())
    assert names == sorted(p.name for p in (tmp_path / "b").iterdir())
    for name in names:
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes()


def test_run_without_asset_surfaces(tmp_path: Path):
    _write_run(tmp_path / "run", thermal=False)
    manifest = export_run(tmp_path / "run")
    assert {entry["name"] for entry in manifest["files"]} == {"sun_earth.csv", "horizon.csv"}
    assert manifest["surfaces"] == {}


def test_terrain_mesh_is_local_and_follows_the_sphere():
    dem = analytic_dem(
        lambda x, y: np.zeros_like(x), center=(SITE_X, SITE_Y), half_width_m=3000.0, pixel_m=20.0
    )
    vertices, faces = terrain_mesh(dem, SITE_X, SITE_Y, half_width_m=2000.0)
    side = math.isqrt(len(vertices))
    assert side * side == len(vertices) and side >= 200
    assert len(faces) == 2 * (side - 1) ** 2
    assert faces.min() == 0 and faces.max() == len(vertices) - 1

    distance = np.hypot(vertices[:, 0], vertices[:, 1])
    assert distance.max() < 2000.0 * np.sqrt(2.0) * 1.01
    # a sphere drops below the tangent plane by d^2 / 2R (9 mm at 180 m, 1.15 m at 2 km)
    np.testing.assert_allclose(vertices[:, 2], -(distance**2) / (2.0 * R), atol=2e-3)
    assert np.abs(vertices[np.argmin(distance)]).max() < 20.0

    a, b, c = (vertices[faces[:, k]] for k in range(3))
    assert (np.cross(b - a, c - a)[:, 2] > 0.0).all()  # every face looks up


def test_terrain_mesh_heights_and_stride():
    def hill(x, y):
        return 30.0 * np.exp(-((x - SITE_X - 100.0) ** 2 + (y - SITE_Y) ** 2) / 150.0**2)

    dem = analytic_dem(hill, center=(SITE_X, SITE_Y), half_width_m=1000.0, pixel_m=5.0)
    full, _ = terrain_mesh(dem, SITE_X, SITE_Y, half_width_m=400.0)
    coarse, faces = terrain_mesh(dem, SITE_X, SITE_Y, half_width_m=400.0, stride=4)
    assert len(coarse) < len(full) / 12
    assert faces.max() == len(coarse) - 1
    # heights are relative to the ground under the observer
    origin_height = float(dem.sample(SITE_X, SITE_Y))
    top = full[np.argmax(full[:, 2])]
    np.testing.assert_allclose(top[:2], [100.0, 0.0], atol=5.0)
    np.testing.assert_allclose(top[2], 30.0 - origin_height, atol=0.05)

    with pytest.raises(ValueError, match="beyond"):
        terrain_mesh(dem, SITE_X, SITE_Y, half_width_m=5000.0)
    with pytest.raises(ValueError, match="stride"):
        terrain_mesh(dem, SITE_X, SITE_Y, half_width_m=400.0, stride=0)


def test_horizon_ring_reproduces_the_horizon_from_the_observer():
    azimuth = np.arange(0.0, 360.0, 5.0)
    elevation = 2.0 * np.sin(np.radians(3.0 * azimuth)) - 0.5
    vertices, faces = horizon_ring(azimuth, elevation, observer_height_m=HEIGHT, radius_m=800.0)
    n = len(azimuth)
    assert vertices.shape == (2 * n, 3) and faces.shape == (2 * n, 3)
    top, bottom = vertices[:n], vertices[n:]
    seen = top - np.array([0.0, 0.0, HEIGHT])
    np.testing.assert_allclose(np.hypot(seen[:, 0], seen[:, 1]), 800.0, rtol=1e-12)
    np.testing.assert_allclose(
        np.degrees(np.arctan2(seen[:, 2], np.hypot(seen[:, 0], seen[:, 1]))), elevation, atol=1e-10
    )
    np.testing.assert_allclose(
        np.degrees(np.arctan2(seen[:, 0], seen[:, 1])) % 360.0, azimuth, atol=1e-9
    )
    assert (bottom[:, 2] < top[:, 2].min()).all()
    # faces look inwards, at the observer
    a, b, c = (vertices[faces[:, k]] for k in range(3))
    normal = np.cross(b - a, c - a)
    centre = (a + b + c) / 3.0
    assert (np.einsum("ij,ij->i", normal[:, :2], -centre[:, :2]) > 0.0).all()


def test_write_obj(tmp_path: Path):
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.5], [0.0, 1.0, -0.25]])
    faces = np.array([[0, 1, 2]])
    path = tmp_path / "t.obj"
    write_obj(path, vertices, faces, ["first line", "second line"])
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[:2] == ["# first line", "# second line"]
    assert [line for line in lines if line.startswith("v ")] == [
        "v 0.0000 0.0000 0.0000",
        "v 1.0000 0.0000 0.5000",
        "v 0.0000 1.0000 -0.2500",
    ]
    assert lines[-1] == "f 1 2 3"


def test_ring_file_and_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    _write_run(tmp_path / "run")
    code = cli_main(
        [
            "export",
            str(tmp_path / "run"),
            "-o",
            str(tmp_path / "out"),
            "--ring-radius",
            "500",
            "--json",
        ]
    )
    assert code == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["directory"] == str(tmp_path / "out")
    ring = next(entry for entry in printed["files"] if entry["name"] == "horizon_ring.obj")
    assert ring["vertices"] == 144 and ring["faces"] == 144
    text = (tmp_path / "out" / "horizon_ring.obj").read_text(encoding="utf-8")
    assert text.count("\nv ") == 144 and text.count("\nf ") == 144

    assert cli_main(["export", str(tmp_path / "missing")]) == 2
    assert "error:" in capsys.readouterr().err


def test_terrain_mesh_has_the_horizon_the_engine_traces():
    """Frame check: seen from the observer, the mesh rises where the tracer's horizon does."""

    def hills(x, y):
        dx, dy = x - SITE_X, y - SITE_Y
        return 40.0 * np.exp(-((dx - 300.0) ** 2 + (dy - 100.0) ** 2) / 80.0**2) + 15.0 * np.exp(
            -((dx + 150.0) ** 2 + (dy + 350.0) ** 2) / 60.0**2
        )

    dem = analytic_dem(hills, center=(SITE_X, SITE_Y), half_width_m=1500.0, pixel_m=5.0)
    stack = TerrainStack((TerrainLevel(dem, 0.0, 900.0, 2.5),))
    traced = compute_horizons(stack, SITE_X, SITE_Y, height_above_surface_m=HEIGHT, n_azimuth=360)
    horizon = np.degrees(traced.elevation_rad[0])
    azimuth = np.degrees(traced.azimuth_rad)

    vertices, _ = terrain_mesh(dem, SITE_X, SITE_Y, half_width_m=1000.0)
    distance = np.hypot(vertices[:, 0], vertices[:, 1])
    inside = (distance > 5.0) & (distance < 900.0)
    seen_azimuth = np.degrees(np.arctan2(vertices[inside, 0], vertices[inside, 1])) % 360.0
    seen = np.degrees(np.arctan2(vertices[inside, 2] - HEIGHT, distance[inside]))

    # no vertex stands above the traced horizon
    limit = np.interp(seen_azimuth, azimuth, horizon, period=360.0)
    assert (seen - limit).max() < 0.05
    # and towards each hill the mesh reaches it: the highest vertex within half a degree
    sector = np.round(seen_azimuth).astype(int) % 360
    mesh_top = np.full(360, -90.0)
    np.maximum.at(mesh_top, sector, seen)
    hilly = horizon > 1.0
    assert hilly.sum() >= 40
    assert (mesh_top - horizon)[hilly].min() > -0.05
    assert abs(azimuth[np.argmax(horizon)] - np.degrees(np.arctan2(300.0, 100.0))) < 2.0
