"""End-to-end runs on real data: outputs, provenance, reproducibility, command line."""

from __future__ import annotations

import json

import numpy as np
import pytest
import xarray as xr

from lunarenv.api import load_run, run_scenario
from lunarenv.cli import main
from lunarenv.geometry.ephemeris import DE440
from lunarenv.scenario import Scenario

pytestmark = pytest.mark.data

REQUIRED = ("pgda78.site04.surf", "pgda90.ldem_80s_80mpp", *DE440.kernel_ids)
SCENARIO = {
    "name": "end-to-end",
    "site": {"x_m": -2822.5, "y_m": -9062.5, "height_above_surface_m": 2.0},
    "period": {
        "start_utc": "2026-03-01T00:00:00",
        "end_utc": "2026-03-15T00:00:00",
        "step_minutes": 120,
    },
    "uncertainty": {"members": 4, "terrain_clones": 2, "seed": 7},
    "thermal": {
        "slope_grid": {"slopes_deg": [0.0, 5.0], "aspects": 2},
        "oriented_surfaces": [
            {"name": "roof", "tilt_deg": 0.0},
            {"name": "wall-north", "tilt_deg": 90.0, "facing_azimuth_deg": 0.0},
        ],
    },
    "illumination": {
        "settings": {"n_azimuth": 360},
        "map": {"half_width_m": 10.0, "n_azimuth": 90},
    },
}


@pytest.fixture(scope="module")
def first(store, require_data, tmp_path_factory):
    return run_scenario(SCENARIO, output_dir=tmp_path_factory.mktemp("runs"), store=store)


@pytest.fixture(scope="module")
def require_data(store):
    missing = [i for i in REQUIRED if not store.is_available(i)]
    if missing:
        pytest.skip(f"not in the data cache: {', '.join(missing)}")


def test_every_variable_has_units_and_the_expected_shape(first):
    dataset = first.dataset
    n = 14 * 12 + 1
    assert dataset.sizes["time"] == n and dataset.sizes["azimuth"] == 360
    for name in dataset.variables:
        if name != "time":
            assert dataset[name].attrs.get("units"), f"{name} has no units"
            assert dataset[name].attrs.get("long_name"), f"{name} has no long_name"
    expected = {
        "sun_elevation",
        "sun_visible_fraction",
        "sun_flux_fraction",
        "solar_irradiance_normal",
        "eclipse_factor",
        "earth_visible_fraction",
        "horizon_elevation",
        "absorbed_solar_flux",
        "surface_temperature",
        "temperature",
        "surface_temperature_constant_albedo",
        "surface_temperature_quantiles",
        "map_mean_visible_fraction",
        "surface_temperature_albedo_envelope",
        "surface_temperature_envelope_quantiles",
        "slope_grid_surface_temperature",
        "oriented_direct_irradiance",
    }
    assert expected <= set(dataset.data_vars)
    assert dataset["temperature"].dims == ("time", "depth")
    # +/-10 m around a pixel centre, rounded outwards to whole 5 m pixels
    assert dataset["map_mean_visible_fraction"].shape == (5, 5)
    q = dataset["surface_temperature_quantiles"].values
    assert np.all(q[0] <= q[1]) and np.all(q[1] <= q[2])
    assert (
        float(dataset["sun_visible_fraction"].min()) >= 0.0
        and float(dataset["sun_visible_fraction"].max()) <= 1.0
    )
    # the 2026-03-03 lunar eclipse is inside the period
    assert float(dataset["eclipse_factor"].min()) < 1.0


def test_metadata_carries_full_provenance(first, store):
    meta = first.metadata
    assert meta["run_id"] == first.run_id and len(first.run_id) == 16
    assert meta["scenario"] == Scenario.model_validate(SCENARIO).resolved()
    recorded = {record["id"]: record["sha256"] for record in meta["data"]}
    for entry in REQUIRED:
        assert recorded[entry] == store.record(entry).sha256
    assert [m["name"] for m in meta["models"]] == ["illumination", "thermal"]
    assert {"python", "numpy", "numba", "spice_toolkit", "platform"} <= set(meta["environment"])
    assert meta["numerical_error"]["thermal"]["surface_temperature_total_k"] == 0.2
    assert meta["validation"]["counts"] and any(
        c["id"] == "ILL-V2" for c in meta["validation"]["cases"]
    )
    # polar thermal results are outside the validated range, and say so
    assert meta["validity"]["within_validated_range"] is False
    assert any("extrapolation" in w for w in meta["validity"]["warnings"])
    assert meta["model_form_sensitivity"]["constant_albedo"]["difference_surface_mean_k"] > 5.0
    assert meta["uncertainty"]["thermal"]["members"] == 4
    assert set(meta["variables"]) == set(first.dataset.variables)
    json.dumps(meta)


def test_same_scenario_and_data_give_the_same_run_id_and_bits(first, store, tmp_path):
    again = run_scenario(SCENARIO, output_dir=tmp_path, store=store)
    assert again.run_id == first.run_id
    assert again.result_checksum == first.result_checksum
    xr.testing.assert_identical(again.dataset, first.dataset)
    different = run_scenario(
        SCENARIO | {"uncertainty": {"members": 4, "terrain_clones": 2, "seed": 8}}, store=store
    )
    assert different.run_id != first.run_id and different.result_checksum != first.result_checksum
    # the nominal results do not depend on the uncertainty seed
    assert np.array_equal(
        different.dataset["surface_temperature"], first.dataset["surface_temperature"]
    )


def test_written_run_reads_back_identically(first):
    assert first.directory is not None and first.directory.name == first.run_id
    assert {p.name for p in first.directory.iterdir()} == {
        "results.nc",
        "metadata.json",
        "scenario.json",
    }
    dataset, metadata = load_run(first.directory)
    assert metadata["result_checksum_sha256"] == first.result_checksum
    assert dataset.attrs["run_id"] == first.run_id
    for name in first.dataset.data_vars:
        assert np.array_equal(dataset[name].values, first.dataset[name].values, equal_nan=True), (
            name
        )
    assert np.array_equal(dataset["time"].values, first.dataset["time"].values)


def test_longitude_latitude_site_and_illumination_only(store, require_data):
    lon, lat = first_site_lonlat()
    geo = run_scenario(
        {
            "name": "geo",
            "site": {"longitude_deg": lon, "latitude_deg": lat, "height_above_surface_m": 2.0},
            "period": SCENARIO["period"],
            "thermal": {"enabled": False},
            "uncertainty": {"members": 0, "terrain_clones": 0},
            "illumination": {"settings": {"n_azimuth": 360}},
        },
        store=store,
    )
    assert "surface_temperature" not in geo.dataset and "thermal" not in geo.metadata["summary"]
    assert geo.metadata["site"]["x_m"] == pytest.approx(-2822.5, abs=0.01)
    assert geo.metadata["validity"]["within_validated_range"] is True


def first_site_lonlat() -> tuple[float, float]:
    from lunarenv.geometry.projection import SouthPolarStereographic

    lon, lat = SouthPolarStereographic().inverse(-2822.5, -9062.5)
    return float(np.degrees(lon)), float(np.degrees(lat))


def test_site_outside_the_tile_is_refused(store, require_data):
    with pytest.raises(ValueError, match="outside"):
        run_scenario(SCENARIO | {"site": {"x_m": 50_000.0, "y_m": 0.0}}, store=store)


def test_command_line_run(store, require_data, tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("LUNARENV_DATA_DIR", str(store.root))
    path = tmp_path / "scenario.json"
    path.write_text(json.dumps(SCENARIO))
    assert main(["run", str(path), "-o", str(tmp_path / "out"), "--json"]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert (
        len(printed["run_id"]) == 16
        and (tmp_path / "out" / printed["run_id"] / "results.nc").exists()
    )


def test_albedo_envelope_brackets_the_nominal_law(first):
    envelope = first.dataset["surface_temperature_albedo_envelope"]
    cold, nominal, hot = (envelope.sel(albedo_case=c).values for c in ("cold", "nominal", "hot"))
    assert np.array_equal(nominal, first.dataset["surface_temperature"].values)
    # at grazing polar incidence the Bürger law reflects most and the Foote law least
    assert cold.max() <= nominal.max() <= hot.max()
    bounds = first.metadata["uncertainty"]["thermal_envelope"]["bounds"]["surface_max_k"]
    assert bounds[0] < bounds[1]
    grid = first.dataset["slope_grid_surface_temperature"]
    assert grid.sizes["orientation"] == 1 + 2 and grid.sizes["albedo_case"] == 3


def test_oriented_surfaces_irradiance_and_view_factors(first):
    surfaces = first.metadata["summary"]["oriented_surfaces"]["surfaces"]
    roof, wall = surfaces["roof"], surfaces["wall-north"]
    assert roof["sky_view_factor"] + roof["ground_view_factor"] == pytest.approx(1.0)
    assert roof["sky_view_factor"] > 0.9 and 0.4 < wall["sky_view_factor"] < 0.6
    irradiance = first.dataset["oriented_direct_irradiance"]
    # with the Sun within 2 degrees of the horizon, a wall facing it gets far more than the roof
    assert float(irradiance.sel(surface="roof").max()) < 60.0
    assert np.all(irradiance.values >= 0.0)
