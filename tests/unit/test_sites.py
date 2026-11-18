"""Terrain sets: every PGDA 78 site is declared, consistently with the catalog and schema."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from lunarenv.cli import main as cli_main
from lunarenv.data.catalog import load_catalog, load_pins
from lunarenv.data.store import DataStore
from lunarenv.scenario import Scenario
from lunarenv.scenario.schema import json_schema
from lunarenv.sites import DEFAULT_TERRAIN, PGDA78_SITES, TERRAIN_SETS, site_status

PERIOD = {"start_utc": "2026-01-01T00:00:00", "end_utc": "2026-01-03T00:00:00"}


def test_all_sites_of_the_product_are_declared():
    assert len(PGDA78_SITES) == len(TERRAIN_SETS) == 27
    assert DEFAULT_TERRAIN in TERRAIN_SETS
    for name, product, _ in PGDA78_SITES:
        assert name == product.lower()
        assert TERRAIN_SETS[name].name == name


def test_every_set_refers_to_catalog_entries():
    catalog = load_catalog()
    for name, product, description in PGDA78_SITES:
        terrain = TERRAIN_SETS[name]
        entry = catalog.get(terrain.near_id)
        assert entry.url == (
            f"https://pgda.gsfc.nasa.gov/data/LOLA_5mpp/{product}/{product}_final_adj_5mpp_surf.tif"
        )
        assert name in entry.groups  # `lunarenv data fetch <terrain>` fetches the site DEM
        assert terrain.regional_id in catalog and terrain.regional_error_id in catalog
        assert description in terrain.description
        if terrain.clone_group is not None:
            assert len(catalog.group(terrain.clone_group)) == 100


def test_every_site_dem_is_pinned():
    pins = load_pins()
    assert [t.name for t in TERRAIN_SETS.values() if t.near_id not in pins] == []


@pytest.mark.parametrize("name", sorted(TERRAIN_SETS))
def test_scenario_accepts_every_set(name):
    scenario = Scenario.model_validate(
        {"name": "s", "site": {"terrain": name, "x_m": 0.0, "y_m": 0.0}, "period": PERIOD}
    )
    assert scenario.site.terrain == name


def test_unknown_set_is_refused_with_the_known_names():
    with pytest.raises(ValidationError, match=r"unknown terrain set 'tycho'.*haworth.*site04"):
        Scenario.model_validate(
            {"name": "s", "site": {"terrain": "tycho", "x_m": 0.0, "y_m": 0.0}, "period": PERIOD}
        )


def test_default_and_schema_enum():
    scenario = Scenario.model_validate(
        {"name": "s", "site": {"x_m": 0.0, "y_m": 0.0}, "period": PERIOD}
    )
    assert scenario.site.terrain == "site04"
    terrain = json_schema()["$defs"]["Site"]["properties"]["terrain"]
    assert terrain["enum"] == sorted(TERRAIN_SETS)
    assert terrain["default"] == "site04"


def test_site_status_on_an_empty_cache(tmp_path, capsys, monkeypatch):
    rows = site_status(DataStore(tmp_path))
    assert [row["name"] for row in rows] == list(TERRAIN_SETS)
    by_name = {row["name"]: row for row in rows}
    assert not any(row["cached"] for row in rows)
    assert by_name["site04"]["pinned"] and by_name["site04"]["clones_cached"] == 0
    assert by_name["site23"]["pinned"] and by_name["site23"]["clones_cached"] is None

    monkeypatch.setenv("LUNARENV_DATA_DIR", str(tmp_path))
    assert cli_main(["sites", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == rows
    assert cli_main(["sites"]) == 0
    text = capsys.readouterr().out
    assert "site23" in text and "Malapert massif" in text and "pinned" in text
