"""Scenario schema.

A scenario says where, when and with which model parameters. Everything that affects a
result is in it, so that scenario + data + code version identify a run. Unknown fields
are rejected rather than ignored: a misspelt option must not silently fall back to a
default.

Units are part of the field names (``_m``, ``_deg``, ``_minutes``, ``_utc``).
"""

from __future__ import annotations

import json
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from lunarenv.illumination.settings import IlluminationSettings
from lunarenv.thermal.model import ThermalSettings
from lunarenv.thermal.properties import RegolithProperties

SCHEMA_VERSION = 1


class _Strict(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class Site(_Strict):
    """Location of the point of interest.

    Give either map coordinates (south polar stereographic metres, as in the LOLA
    products) or longitude and latitude. ``terrain`` names the terrain data set that
    covers the site.
    """

    terrain: Literal["site04", "site01"] = Field(
        "site04", description="Terrain set: 5 m site DEM plus regional DEM (see lunarenv.sites)."
    )
    x_m: float | None = Field(None, description="Polar stereographic x, metres.")
    y_m: float | None = Field(None, description="Polar stereographic y, metres.")
    longitude_deg: float | None = Field(None, ge=-360.0, le=360.0, description="East longitude.")
    latitude_deg: float | None = Field(
        None, ge=-90.0, lt=-60.0, description="Planetocentric latitude."
    )
    height_above_surface_m: float = Field(
        2.0,
        ge=0.0,
        le=1000.0,
        description="Height of the point above the terrain, for illumination and Earth visibility.",
    )

    @model_validator(mode="after")
    def _one_coordinate_pair(self) -> Site:
        has_map = self.x_m is not None and self.y_m is not None
        has_geo = self.longitude_deg is not None and self.latitude_deg is not None
        partial = (self.x_m is None) != (self.y_m is None) or (self.longitude_deg is None) != (
            self.latitude_deg is None
        )
        if partial or has_map == has_geo:
            raise ValueError("give exactly one of (x_m, y_m) or (longitude_deg, latitude_deg)")
        return self


class Period(_Strict):
    """Time span of the results, in UTC, sampled uniformly."""

    start_utc: datetime
    end_utc: datetime
    step_minutes: int = Field(60, ge=1, le=360, description="Sampling step of all time series.")

    @model_validator(mode="after")
    def _check(self) -> Period:
        start, end = _as_utc(self.start_utc), _as_utc(self.end_utc)
        if end - start < timedelta(minutes=2 * self.step_minutes):
            raise ValueError("the period must span at least two steps")
        if (end - start) % timedelta(minutes=self.step_minutes):
            raise ValueError("the period must be a whole number of steps")
        if (24 * 60) % self.step_minutes:
            raise ValueError(
                "step_minutes must divide a day, so that the spin-up year is a whole number of steps"
            )
        return self


def _as_utc(value: datetime) -> datetime:
    return value.replace(tzinfo=UTC) if value.tzinfo is None else value.astimezone(UTC)


class MapBlock(_Strict):
    """Optional map of illumination statistics over a square window centred on the site."""

    half_width_m: float = Field(250.0, gt=0.0, le=2000.0)
    n_azimuth: int = Field(
        720, ge=8, description="Horizon directions per pixel (cost is linear in it)."
    )


class IlluminationBlock(_Strict):
    enabled: bool = True
    settings: IlluminationSettings = Field(default_factory=IlluminationSettings)
    earth_visibility: bool = Field(True, description="Also compute visibility of the Earth's disc.")
    map: MapBlock | None = None


class SurfaceOrientation(_Strict):
    """Orientation of the surface element whose temperature is computed."""

    source: Literal["dem", "explicit"] = Field(
        "dem",
        description="'dem': slope of the terrain model at the site; 'explicit': values below.",
    )
    slope_deg: float = Field(0.0, ge=0.0, lt=90.0)
    facing_azimuth_deg: float = Field(
        0.0,
        ge=0.0,
        lt=360.0,
        description="Grid azimuth (from map +y towards +x) the surface faces.",
    )
    dem_baseline_m: float | None = Field(
        None,
        gt=0.0,
        description="Half-length over which the DEM slope is measured (default: one pixel).",
    )


class ThermalBlock(_Strict):
    enabled: bool = True
    properties: RegolithProperties = Field(default_factory=RegolithProperties)
    settings: ThermalSettings = Field(default_factory=ThermalSettings)
    surface: SurfaceOrientation = Field(default_factory=SurfaceOrientation)
    constant_albedo_alternative: bool = Field(
        True,
        description=(
            "Also compute the result with a constant albedo. The incidence-dependent law was fitted on "
            "level equatorial ground; at grazing polar incidence the two differ by tens of kelvin. "
            "The first polar comparison with Diviner favours the constant albedo, on thin evidence "
            "(VALIDATION.md, TH-V3 and TH-V3b)."
        ),
    )


class UncertaintyBlock(_Strict):
    """Monte Carlo propagation of terrain and regolith-parameter uncertainty."""

    members: int = Field(64, ge=0, le=4096, description="Regolith-parameter samples (0 disables).")
    terrain_clones: int = Field(
        16,
        ge=0,
        le=100,
        description="Terrain error clones to use, taken in catalog order from the cache.",
    )
    seed: int = Field(20261112, ge=0)


class Scenario(_Strict):
    """A complete, self-contained request."""

    schema_version: Literal[1] = 1
    name: str = Field(min_length=1, max_length=120, pattern=r"^[A-Za-z0-9][A-Za-z0-9 _.\-]*$")
    description: str = ""
    site: Site
    period: Period
    illumination: IlluminationBlock = Field(default_factory=IlluminationBlock)
    thermal: ThermalBlock = Field(default_factory=ThermalBlock)
    uncertainty: UncertaintyBlock = Field(default_factory=UncertaintyBlock)

    @model_validator(mode="after")
    def _dependencies(self) -> Scenario:
        if self.thermal.enabled and not self.illumination.enabled:
            raise ValueError("the thermal model needs the illumination model")
        return self

    def resolved(self) -> dict[str, Any]:
        """Every field with its value, defaults included, as plain JSON types.

        This is what identifies the run: two scenarios that differ only in spelling out
        a default resolve to the same dictionary.
        """
        data = self.model_dump(mode="json")
        for key in ("start_utc", "end_utc"):
            moment = _as_utc(getattr(self.period, key))
            data["period"][key] = moment.strftime("%Y-%m-%dT%H:%M:%SZ")
        return data


def load_scenario(source: Scenario | dict[str, Any] | str | Path) -> Scenario:
    """Scenario from an object, a dictionary, or a YAML/JSON file."""
    if isinstance(source, Scenario):
        return source
    if isinstance(source, dict):
        return Scenario.model_validate(source)
    path = Path(source)
    text = path.read_text(encoding="utf-8")
    data = json.loads(text) if path.suffix.lower() == ".json" else yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: a scenario must be a mapping")
    return Scenario.model_validate(data)


def json_schema() -> dict[str, Any]:
    """JSON Schema of a scenario (draft 2020-12, as produced by pydantic)."""
    return Scenario.model_json_schema()
