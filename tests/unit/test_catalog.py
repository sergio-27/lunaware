"""Catalog and pin-file consistency."""

from __future__ import annotations

import pytest

from lunarenv.data.catalog import Pin, load_catalog, load_pins, parse_catalog, save_pin

ENTRY = """
[[file]]
id = "a.b"
url = "{scheme}://example.org/dir/file.tif"
kind = "geotiff"
groups = ["g"]
product = "p"
version = "v"
citation = "c"
description = "d"
"""


def test_shipped_catalog_is_well_formed():
    catalog = load_catalog()
    assert len(catalog.ids) == len(set(catalog.ids))
    assert len(catalog.group("site04-clones")) == 100
    assert catalog.get("pgda78.site04.clone.0100").url.endswith("_0100_err.tif")
    assert all(catalog.get(i).url.startswith("https://") for i in catalog.ids)


def test_every_pin_refers_to_a_catalog_entry():
    catalog = load_catalog()
    assert set(load_pins()) <= set(catalog.ids)


def test_relative_path_mirrors_url():
    entry = parse_catalog(ENTRY.format(scheme="https")).get("a.b")
    assert str(entry.relative_path) == "example.org/dir/file.tif"


def test_plain_http_is_rejected():
    with pytest.raises(ValueError, match="https"):
        parse_catalog(ENTRY.format(scheme="http"))


def test_duplicate_ids_are_rejected():
    with pytest.raises(ValueError, match="duplicate"):
        parse_catalog(ENTRY.format(scheme="https") * 2)


def test_unknown_group_and_id_raise():
    catalog = parse_catalog(ENTRY.format(scheme="https"))
    with pytest.raises(KeyError):
        catalog.get("nope")
    with pytest.raises(KeyError):
        catalog.group("nope")


def test_pins_are_append_only(tmp_path):
    path = tmp_path / "pins.json"
    first = Pin(sha256="a" * 64, bytes=10)
    save_pin("x", first, path)
    save_pin("x", first, path)  # idempotent
    assert load_pins(path) == {"x": first}
    with pytest.raises(ValueError, match="already pinned"):
        save_pin("x", Pin(sha256="b" * 64, bytes=10), path)
