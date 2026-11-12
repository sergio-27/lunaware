"""Scenario schema, run identity, model descriptions and command line (no data needed)."""

from __future__ import annotations

import json

import numpy as np
import pytest
import yaml
from pydantic import ValidationError

import lunarenv
from lunarenv.cli import build_parser, main
from lunarenv.core.model import EnvironmentModel, ModelDescription, registered_models
from lunarenv.core.provenance import arrays_checksum, canonical_json, compute_run_id
from lunarenv.scenario import Scenario, json_schema, load_scenario

BASE = {
    "name": "t",
    "site": {"x_m": -2822.5, "y_m": -9062.5},
    "period": {"start_utc": "2026-01-01T00:00:00", "end_utc": "2026-01-03T00:00:00"},
}


def scenario(**changes) -> Scenario:
    return Scenario.model_validate(BASE | changes)


# ---------------------------------------------------------------------------- schema
def test_defaults_resolve_and_spelling_out_a_default_does_not_change_identity():
    implicit = scenario().resolved()
    explicit = scenario(uncertainty={"members": 64}, thermal={"enabled": True}).resolved()
    assert implicit == explicit
    assert implicit["period"]["start_utc"] == "2026-01-01T00:00:00Z"  # normalized to UTC
    assert implicit["thermal"]["properties"]["albedo_normal"] == 0.12  # defaults are spelled out
    assert implicit["illumination"]["settings"]["n_azimuth"] == 1440


@pytest.mark.parametrize(
    "bad",
    [
        {"site": {"x_m": 1.0, "y_m": 2.0, "elevation": 3.0}},  # unknown field
        {"site": {"x_m": 1.0}},  # half a coordinate pair
        {"site": {"x_m": 1.0, "y_m": 2.0, "longitude_deg": 10.0, "latitude_deg": -89.0}},  # both
        {"site": {"longitude_deg": 10.0, "latitude_deg": -45.0}},  # outside the polar region
        {"site": {"x_m": 1.0, "y_m": 2.0, "terrain": "tycho"}},  # unknown terrain set
        {"period": {"start_utc": "2026-01-01", "end_utc": "2026-01-01"}},  # empty
        {
            "period": {"start_utc": "2026-01-01T00:00", "end_utc": "2026-01-01T02:30"}
        },  # not whole steps
        {
            "period": {"start_utc": "2026-01-01", "end_utc": "2026-01-05", "step_minutes": 7}
        },  # does not divide a day
        {"thermal": {"propertys": {}}},  # misspelt block
        {"thermal": {"properties": {"emissivity": 1.4}}},
        {"illumination": {"enabled": False}},  # thermal needs illumination
        {"schema_version": 2},
        {"name": "../escape"},
        {"uncertainty": {"members": -1}},
    ],
)
def test_invalid_scenarios_are_rejected(bad):
    with pytest.raises(ValidationError):
        scenario(**bad)


def test_illumination_only_scenario_is_valid():
    only = scenario(thermal={"enabled": False}, illumination={"enabled": True})
    assert not only.thermal.enabled


def test_loading_from_yaml_json_and_dict(tmp_path):
    as_yaml = tmp_path / "s.yaml"
    as_yaml.write_text(yaml.safe_dump(BASE))
    as_json = tmp_path / "s.json"
    as_json.write_text(json.dumps(BASE))
    loaded = [
        load_scenario(as_yaml),
        load_scenario(as_json),
        load_scenario(BASE),
        load_scenario(scenario()),
    ]
    assert all(s == loaded[0] for s in loaded)
    (tmp_path / "list.yaml").write_text("- 1\n- 2\n")
    with pytest.raises(ValueError, match="mapping"):
        load_scenario(tmp_path / "list.yaml")


def test_json_schema_is_exportable_and_strict():
    schema = json_schema()
    assert schema["type"] == "object" and schema["additionalProperties"] is False
    assert {"name", "site", "period"} <= set(schema["required"])
    assert schema["$defs"]["Site"]["additionalProperties"] is False
    json.dumps(schema)  # serializable


# -------------------------------------------------------------------------- identity
def test_run_id_depends_on_every_input_and_on_nothing_else():
    resolved = scenario().resolved()
    data = {"a": "1" * 64, "b": "2" * 64}
    models = {"illumination": "0.1", "thermal": "0.1"}
    reference = compute_run_id(resolved, data, models, "0.1.0")
    assert len(reference) == 16
    assert (
        compute_run_id(resolved, dict(reversed(data.items())), models, "0.1.0") == reference
    )  # order-free
    assert compute_run_id(scenario(name="other").resolved(), data, models, "0.1.0") != reference
    assert compute_run_id(resolved, data | {"a": "3" * 64}, models, "0.1.0") != reference
    assert compute_run_id(resolved, data, models | {"thermal": "0.2"}, "0.1.0") != reference
    assert compute_run_id(resolved, data, models, "0.1.1") != reference


def test_array_checksum_is_sensitive_to_bits_shape_and_names():
    a = {"t": np.arange(6.0).reshape(2, 3), "n": np.array([1, 2], dtype=np.int64)}
    reference = arrays_checksum(a)
    assert arrays_checksum(dict(reversed(a.items()))) == reference
    assert arrays_checksum({**a, "t": a["t"].reshape(3, 2)}) != reference
    assert arrays_checksum({**a, "t": a["t"].astype(np.float32)}) != reference
    flipped = a["t"].copy()
    flipped[0, 0] = np.nextafter(flipped[0, 0], 1.0)
    assert arrays_checksum({**a, "t": flipped}) != reference
    assert arrays_checksum({"u": a["t"], "n": a["n"]}) != reference
    with pytest.raises(ValueError):
        canonical_json({"x": float("nan")})


# --------------------------------------------------------------------- models and API
def test_models_describe_themselves_completely():
    models = registered_models()
    assert list(models) == ["illumination", "thermal"]  # dependency order
    for model in models.values():
        assert isinstance(model, EnvironmentModel)
        description = model.describe()
        assert isinstance(description, ModelDescription)
        for field in (
            "assumptions",
            "simplifications",
            "validity_range",
            "uncertainty_sources",
            "validation_cases",
        ):
            assert getattr(description, field), (model.name, field)
        assert description.numerical_error and description.documentation.endswith(".md")
    described = lunarenv.describe_models()
    assert [m["name"] for m in described] == ["illumination", "thermal"]
    thermal_parameters = {p["name"]: p for p in described[1]["parameters"]}
    assert "0.047" in thermal_parameters["scale_height_m"]["uncertainty"]
    assert thermal_parameters["radiative_parameter"]["uncertainty"] == "not varied"
    json.dumps(described)


def test_public_api_surface():
    for name in ("run_scenario", "describe_models", "get_validation_status"):
        assert callable(getattr(lunarenv, name))
    with pytest.raises(AttributeError):
        lunarenv.no_such_thing  # noqa: B018
    status = lunarenv.get_validation_status()
    assert {"counts", "cases", "generated_utc", "lunarenv_version"} == set(status)


# ------------------------------------------------------------------------------- CLI
def test_cli_describe_schema_and_validate(capsys):
    assert main(["describe", "--json"]) == 0
    assert [m["name"] for m in json.loads(capsys.readouterr().out)] == ["illumination", "thermal"]
    assert main(["schema"]) == 0
    assert json.loads(capsys.readouterr().out)["title"] == "Scenario"
    assert main(["validate", "--json"]) == 0
    assert "cases" in json.loads(capsys.readouterr().out)
    assert main(["describe"]) == 0 and "assumptions" in capsys.readouterr().out


def test_cli_reports_errors_without_a_traceback(tmp_path, capsys):
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(BASE | {"site": {"x_m": 1.0}}))
    assert main(["run", str(bad), "-o", str(tmp_path)]) == 2
    assert "error:" in capsys.readouterr().err
    assert main(["run", str(tmp_path / "missing.yaml")]) == 2
    with pytest.raises(SystemExit):
        build_parser().parse_args(["frobnicate"])
