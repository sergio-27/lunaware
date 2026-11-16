"""Execution of a scenario: illumination, thermal, uncertainty, results and provenance.

The runner contains no physics. It turns a scenario into calls of the model kernels,
collects their outputs in one self-describing dataset and writes, next to it, a metadata
document with everything needed to audit or repeat the run: resolved scenario, data
checksums, model descriptions, software environment, uncertainty, numerical error,
validity warnings and validation status.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from lunarenv import __version__
from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4
from lunarenv.core.model import ModelDescription, ParameterInfo, register, registered_models
from lunarenv.core.provenance import arrays_checksum, compute_run_id, environment
from lunarenv.data.dem import Bounds, load_dem
from lunarenv.data.store import DataNotAvailableError, DataStore
from lunarenv.geometry.ephemeris import DE440, spice_session
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.illumination.eclipse import earth_occultation_factor
from lunarenv.illumination.horizon import compute_horizons, sun_relevant_incomplete
from lunarenv.illumination.products import illumination_series, illumination_statistics, longest_run
from lunarenv.illumination.settings import NUMERICAL_ERROR_ESTIMATE as ILLUMINATION_ERROR
from lunarenv.illumination.settings import IlluminationSettings
from lunarenv.illumination.terrain import TerrainStack, build_stack
from lunarenv.scenario.schema import Scenario, load_scenario
from lunarenv.sites import TERRAIN_SETS
from lunarenv.thermal import terrain_radiation as tr
from lunarenv.thermal.albedo_laws import ENVELOPE
from lunarenv.thermal.model import NUMERICAL_ERROR_ESTIMATE as THERMAL_ERROR
from lunarenv.thermal.model import run_thermal, summarize
from lunarenv.thermal.properties import RegolithProperties, albedo
from lunarenv.thermal.surface import (
    absorbed_solar_flux,
    direct_irradiance,
    neglected_terrain_flux_bound,
    normal_from_slope,
    plane_normal,
    slope_and_facing,
    surface_normal,
    view_factors,
)
from lunarenv.thermal.uq import DEFAULT_RANGES, sample_properties
from lunarenv.validation.core import get_validation_status

PROJ = SouthPolarStereographic()
QUANTILES = (0.05, 0.5, 0.95)
_SPINUP_DAYS = 365
#: Warn when the neglected terrain flux bound exceeds this share of the emission at the
#: surface minimum (scientific review 2026-11-13, finding F06).
_NEGLECTED_FLUX_WARNING_RATIO = 0.05
#: Warn when the other conductivity law moves a mean at depth by more than this.
_CONDUCTIVITY_WARNING_K = 1.0
#: Diffuse irradiance components, in the order stored in the results.
DIFFUSE_COMPONENTS = (
    "terrain_scattered",
    "terrain_infrared",
    "earthshine_reflected",
    "earthshine_infrared",
)


# ------------------------------------------------------------------ model descriptions
class IlluminationModel:
    name = "illumination"
    version = "0.1"
    requires: tuple[str, ...] = ()

    def describe(self) -> ModelDescription:
        return ModelDescription(
            name=self.name,
            version=self.version,
            summary=(
                "Terrain horizon from nested LOLA elevation models by great-circle ray marching on the "
                "sphere; apparent Sun and Earth from SPICE (DE440); visible fraction of the solar disc "
                "by strip integration with limb darkening; Earth occultation of the Sun."
            ),
            outputs=(
                "horizon_elevation",
                "sun_elevation",
                "sun_azimuth",
                "sun_visible_fraction",
                "sun_flux_fraction",
                "solar_irradiance_normal",
                "eclipse_factor",
                "earth_elevation",
                "earth_visible_fraction",
            ),
            assumptions=(
                "The terrain is the DEM surface, bilinear between pixel centres.",
                "The Moon is a sphere of radius 1737.4 km; the local vertical is radial.",
                "Illumination is direct sunlight only.",
            ),
            simplifications=(
                "No relief below the DEM resolution (boulders, small craters).",
                "Sky treated as flat across the solar disc; horizon linear between azimuth samples.",
                "Earth occultation as two uniform discs with a sharp limb (no atmosphere).",
                "'Lit' means any visible fraction of the disc above zero.",
            ),
            validity_range=(
                "Observers inside the site tile of the chosen terrain set.",
                "Epochs 1849-2150 (DE440 short file); compared with JPL Horizons for 2026-2028.",
                "Observer height 0 m and above; results at 0 m depend on relief the DEM does not contain.",
            ),
            parameters=tuple(
                ParameterInfo(name=name, default=field.default, description=field.description or "")
                for name, field in IlluminationSettings.model_fields.items()
            ),
            numerical_error=dict(ILLUMINATION_ERROR),
            uncertainty_sources=(
                "Near-field terrain height: PGDA error clones of the site DEM (Barker2021).",
                "Not sampled: regional DEM errors; relief below the DEM resolution.",
            ),
            validation_cases=(
                "ILL-A1",
                "ILL-A2",
                "ILL-A3",
                "ILL-X1",
                "ILL-V1",
                "ILL-V2",
                "ILL-V3",
                "ILL-V4",
            ),
            references=("Mazarico2011", "Barker2021", "Barker2023", "Park2021"),
            documentation="docs/models/illumination.md",
        )


class ThermalModel:
    name = "thermal"
    version = "0.4"  # 0.4: low-temperature conductivity alternative; 0.3: single-bounce terrain radiation and Earthshine; 0.2: albedo envelope, slope grid, asset surfaces
    requires: tuple[str, ...] = ("illumination",)

    def describe(self) -> ModelDescription:
        nominal = RegolithProperties()
        intervals = {r.field: f"[{r.low:g}, {r.high:g}] - {r.basis}" for r in DEFAULT_RANGES}
        return ModelDescription(
            name=self.name,
            version=self.version,
            summary=(
                "One-dimensional heat conduction in regolith with depth-dependent density and "
                "temperature-dependent conductivity and heat capacity (Hayne2017); radiating surface "
                "driven by absorbed direct sunlight on the local slope, single-bounce terrain-scattered "
                "sunlight and terrain infrared, and Earthshine; interior heat flow at the base."
            ),
            requires=self.requires,
            outputs=(
                "absorbed_solar_flux",
                "absorbed_diffuse_flux",
                "diffuse_irradiance",
                "surface_temperature",
                "temperature",
            ),
            assumptions=(
                "Heat flows vertically only; properties depend on depth and temperature.",
                "The surface element has the slope of the DEM at the site unless given explicitly.",
            ),
            simplifications=(
                "Terrain radiation is single bounce: terrain is lit by the Sun only, scatters "
                "Lambertian with its albedo law and emits with its own 1-D temperature, grouped in "
                "patches; multiple scattering and sub-DEM roughness are not included.",
                "Earthshine from the HLS-UG-001 peak values scaled by the illuminated fraction of "
                "the Earth and the Earth-Moon distance.",
                "Ensemble members and the constant-albedo alternative reuse the nominal terrain's "
                "diffuse irradiance.",
                "No rocks, no lateral conduction, no volatiles.",
                "Incidence angle of the disc centre for partially visible discs.",
            ),
            validity_range=(
                "Regolith properties fitted to Diviner data within 70 degrees of the equator; polar use "
                "is an extrapolation.",
                "Heat-capacity fit stated valid from below 90 K to above 400 K.",
                "Albedo law fitted on level equatorial ground; at grazing polar incidence it is the "
                "largest model-form uncertainty; results are bounded by the published laws "
                "(cold: Bürger 2024, hot: Foote 2020) and a constant-albedo alternative.",
            ),
            parameters=tuple(
                ParameterInfo(
                    name=name,
                    default=getattr(nominal, name)
                    if not isinstance(getattr(nominal, name), tuple)
                    else str(getattr(nominal, name)),
                    description=field.description or "",
                    uncertainty=intervals.get(name, "not varied"),
                )
                for name, field in RegolithProperties.model_fields.items()
            ),
            numerical_error=dict(THERMAL_ERROR),
            uncertainty_sources=(
                "Regolith parameters: Latin hypercube over stated intervals (two of them assumptions).",
                "Terrain: slope and illumination of each DEM error clone.",
                "Not sampled: the albedo law (bounded by the envelope instead), chi, the "
                "heat-capacity law, the low-temperature conductivity law (the other law is run "
                "separately and its change reported instead).",
            ),
            validation_cases=(
                "TH-A1",
                "TH-A2",
                "TH-P1",
                "TH-V1",
                "TH-V2",
                "TH-V2b",
                "TH-X1",
                "TH-V3",
                "TH-A3",
                "TH-C1",
                "TH-V2c",
                "TH-A2b",
                "TH-C2",
            ),
            references=(
                "Hayne2017",
                "WoodsRobinson2019",
                "Vasavada2012",
                "Keihm1984",
                "heat1d",
                "HLS-UG-001",
                "Ingersoll1992",
            ),
            documentation="docs/models/thermal.md",
        )


register(IlluminationModel())
register(ThermalModel())


def describe_models() -> list[dict[str, Any]]:
    """Descriptions of every registered model, in dependency order."""
    return [m.describe().model_dump(mode="json") for m in registered_models().values()]


# ------------------------------------------------------------------------------ run
@dataclass(frozen=True, slots=True)
class _Ephemeris:
    """Apparent positions (Moon-centred, body-fixed, metres) over the whole time axis."""

    et: np.ndarray
    sun: np.ndarray
    earth: np.ndarray
    r_sun: float
    r_earth: float


@dataclass(frozen=True, slots=True)
class RunResult:
    """Outcome of :func:`run_scenario`."""

    run_id: str
    result_checksum: str
    dataset: xr.Dataset
    metadata: dict[str, Any]
    directory: Path | None


def _variable(values: np.ndarray, dims: tuple[str, ...], units: str, long_name: str) -> xr.Variable:
    return xr.Variable(dims, values, {"units": units, "long_name": long_name})


#: Albedo cases of the envelope as stored in the results: the scenario's own properties
#: are the middle case (the Hayne2017 law by default).
ALBEDO_CASES = ("cold", "nominal", "hot")


def _envelope_properties(properties: RegolithProperties) -> dict[str, RegolithProperties]:
    return {
        "cold": ENVELOPE["cold"].apply(properties),
        "nominal": properties,
        "hot": ENVELOPE["hot"].apply(properties),
    }


def _surface_stats(result: Any) -> dict[str, float]:
    surface = result.surface_temperature_k
    return {
        "surface_min_k": float(surface.min()),
        "surface_max_k": float(surface.max()),
        "surface_mean_k": float(surface.mean()),
        "mean_at_1m_k": float(result.at_depth(min(1.0, float(result.depth_m[-1]))).mean()),
    }


def _diffuse_absorbed(
    parts: dict[str, np.ndarray],
    earth_incidence_rad: np.ndarray | None,
    properties: RegolithProperties,
) -> np.ndarray:
    """Absorbed diffuse flux: terrain light with the hemispherical albedo, Earthshine with
    the albedo at its incidence, infrared with the emissivity."""
    eps = properties.emissivity
    out = (1.0 - tr.diffuse_albedo(properties)) * parts["terrain_scattered"] + eps * (
        parts["terrain_infrared"] + parts["earthshine_infrared"]
    )
    if earth_incidence_rad is not None:
        out = out + (1.0 - albedo(earth_incidence_rad, properties)) * parts["earthshine_reflected"]
    return np.asarray(out)


def _utc(moment: np.datetime64) -> str:
    return str(np.datetime_as_string(moment, unit="s")) + "Z"


def _terrain_products(
    stack: TerrainStack,
    xy: tuple[float, float],
    scenario: Scenario,
    eph: _Ephemeris,
    spin: int,
) -> dict[str, Any]:
    """Everything derived from one terrain realization at the site."""
    settings = scenario.illumination.settings
    height = scenario.site.height_above_surface_m
    structure = compute_horizons(
        stack, *xy, height_above_surface_m=height, n_azimuth=settings.n_azimuth
    )
    ground = (
        structure if height == 0.0 else compute_horizons(stack, *xy, n_azimuth=settings.n_azimuth)
    )
    period = slice(spin, None)
    sun = illumination_series(
        structure,
        eph.et[period],
        eph.sun[period],
        eph.r_sun,
        n_strips=settings.n_strips,
        limb_u=settings.limb_darkening,
        mode=settings.disc_mode_code,
    )
    out: dict[str, Any] = {"structure": structure, "ground": ground, "sun": sun}
    position = PROJ.to_cartesian(*xy, float(ground.observer_height_m[0]))
    out["eclipse"] = earth_occultation_factor(eph.sun, eph.earth, position, eph.r_sun, eph.r_earth)
    if scenario.thermal.enabled:
        out["sun_ground"] = illumination_series(
            ground,
            eph.et,
            eph.sun,
            eph.r_sun,
            n_strips=settings.n_strips,
            limb_u=settings.limb_darkening,
            mode=settings.disc_mode_code,
        )
        surface = scenario.thermal.surface
        if surface.source == "dem":
            out["normal"] = surface_normal(
                stack.levels[0].dem, *xy, baseline_m=surface.dem_baseline_m
            )
        else:
            out["normal"] = normal_from_slope(
                np.radians(surface.slope_deg), np.radians(surface.facing_azimuth_deg)
            )
    return out


def run_scenario(
    source: Scenario | dict[str, Any] | str | Path,
    *,
    output_dir: str | Path | None = None,
    store: DataStore | None = None,
) -> RunResult:
    """Run a scenario and return (and optionally write) its results.

    With ``output_dir`` the run is written to ``<output_dir>/<run_id>/`` as
    ``results.nc`` (NetCDF-4, CF attributes), ``metadata.json`` and ``scenario.json``.
    The same scenario with the same data and code always maps to the same ``run_id``.
    """
    started = time.perf_counter()
    scenario = load_scenario(source)
    store = store or DataStore()
    terrain = TERRAIN_SETS[scenario.site.terrain]
    required = (*terrain.required_ids, *DE440.kernel_ids)
    missing = [i for i in required if not store.is_available(i)]
    if missing:
        raise DataNotAvailableError(
            f"scenario needs data that is not in the cache at {store.root}: {', '.join(missing)}; "
            "fetch it with `lunarenv data fetch <id>`"
        )
    cached_clones = [
        e.id for e in store.catalog.group(terrain.clone_group) if store.is_available(e.id)
    ]
    clone_ids = cached_clones[: scenario.uncertainty.terrain_clones]
    used_ids = (*required, *clone_ids)
    checksums = {i: store.record(i).sha256 for i in used_ids}
    models = registered_models()
    versions = {name: model.version for name, model in models.items()}
    run_id = compute_run_id(scenario.resolved(), checksums, versions, __version__)
    warnings: list[str] = []
    timing: dict[str, float] = {}

    # ----------------------------------------------------------------- site and time
    near = load_dem(store, terrain.near_id)
    stack = build_stack(
        near,
        load_dem(store, terrain.regional_id),
        near_radius_m=scenario.illumination.settings.near_radius_m,
        mid_radius_m=scenario.illumination.settings.mid_radius_m,
        step_fraction=scenario.illumination.settings.step_fraction,
        far_coarsening=scenario.illumination.settings.far_coarsening,
    )
    site = scenario.site
    if site.x_m is not None and site.y_m is not None:
        xy = (float(site.x_m), float(site.y_m))
    else:
        assert site.longitude_deg is not None and site.latitude_deg is not None  # schema guarantee
        gx, gy = PROJ.forward(np.radians(site.longitude_deg), np.radians(site.latitude_deg))
        xy = (float(gx), float(gy))
    try:
        surface_height = float(near.sample(*xy))
    except ValueError as error:
        raise ValueError(
            f"site ({xy[0]:.1f}, {xy[1]:.1f}) m lies outside the '{terrain.name}' terrain tile {near.grid.bounds}"
        ) from error
    longitude, latitude = (float(v) for v in PROJ.inverse(*xy))

    step_min = scenario.period.step_minutes
    step_s = 60.0 * step_min
    start = np.datetime64(scenario.resolved()["period"]["start_utc"][:-1], "s")
    end = np.datetime64(scenario.resolved()["period"]["end_utc"][:-1], "s")
    n_period = int((end - start) / np.timedelta64(step_min, "m"))
    thermal_on = scenario.thermal.enabled
    spin = _SPINUP_DAYS * 24 * 60 // step_min if thermal_on else 0
    times = (
        start
        - spin * np.timedelta64(step_min, "m")
        + np.arange(spin + n_period + 1) * np.timedelta64(step_min, "m")
    )
    period_times = times[spin:]
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(times)
        eph = _Ephemeris(
            et=et,
            sun=ephemeris.position_from_center("SUN", et),
            earth=ephemeris.position_from_center("EARTH", et),
            r_sun=ephemeris.body_radius_m("SUN"),
            r_earth=ephemeris.body_radius_m("EARTH"),
        )

    # ------------------------------------------------------------------ illumination
    tick = time.perf_counter()
    nominal = _terrain_products(stack, xy, scenario, eph, spin)
    structure, sun = nominal["structure"], nominal["sun"]
    eclipse_period = nominal["eclipse"][spin:]
    variables: dict[str, xr.Variable] = {
        "horizon_elevation": _variable(
            np.degrees(structure.elevation_rad[0].astype(np.float64)),
            ("azimuth",),
            "degree",
            "terrain horizon elevation at the observer height",
        ),
        "horizon_distance": _variable(
            structure.distance_m[0].astype(np.float64),
            ("azimuth",),
            "m",
            "ground distance to the terrain point defining the horizon",
        ),
        "sun_elevation": _variable(
            np.degrees(sun.elevation_rad),
            ("time",),
            "degree",
            "apparent elevation of the Sun's centre above the local horizontal",
        ),
        "sun_azimuth": _variable(
            np.degrees(sun.azimuth_rad),
            ("time",),
            "degree",
            "grid azimuth of the Sun's centre (from map +y towards +x)",
        ),
        "sun_visible_fraction": _variable(
            sun.visible_fraction,
            ("time",),
            "1",
            "fraction of the solar disc area above the terrain horizon",
        ),
        "sun_flux_fraction": _variable(
            sun.flux_fraction,
            ("time",),
            "1",
            "limb-darkening-weighted fraction of the solar disc above the terrain horizon",
        ),
        "eclipse_factor": _variable(
            eclipse_period, ("time",), "1", "fraction of the solar disc not covered by the Earth"
        ),
        "solar_irradiance_normal": _variable(
            sun.normal_irradiance_w_m2 * eclipse_period,
            ("time",),
            "W m-2",
            "direct solar irradiance on a surface facing the Sun",
        ),
    }
    lit = sun.visible_fraction > 0.0
    shadow, light = sun.longest_shadow(), sun.longest_illumination()
    summary: dict[str, Any] = {
        "illumination": {
            "mean_visible_fraction": sun.mean_visible_fraction,
            "mean_flux_fraction": float(sun.flux_fraction.mean()),
            "fraction_of_time_lit": float(lit.mean()),
            "longest_shadow_hours": shadow.duration_s / 3600.0,
            "longest_shadow_start_utc": _utc(period_times[shadow.start_index])
            if shadow.count
            else None,
            "longest_illumination_hours": light.duration_s / 3600.0,
            "longest_illumination_start_utc": _utc(period_times[light.start_index])
            if light.count
            else None,
            "eclipse_hours": float(np.sum(eclipse_period < 1.0) * step_s / 3600.0),
        }
    }
    if bool(sun_relevant_incomplete(structure, stack).any()):
        warnings.append(
            "illumination: some horizon rays left the terrain model while terrain beyond it could still "
            "hide the Sun; the regional DEM is too small for this observer"
        )
    if scenario.illumination.earth_visibility:
        settings = scenario.illumination.settings
        earth = illumination_series(
            structure,
            eph.et[spin:],
            eph.earth[spin:],
            eph.r_earth,
            n_strips=settings.n_strips,
            limb_u=0.0,
            mode=settings.disc_mode_code,
        )
        variables["earth_elevation"] = _variable(
            np.degrees(earth.elevation_rad),
            ("time",),
            "degree",
            "apparent elevation of the Earth's centre above the local horizontal",
        )
        variables["earth_visible_fraction"] = _variable(
            earth.visible_fraction,
            ("time",),
            "1",
            "fraction of the Earth's disc above the terrain horizon",
        )
        summary["illumination"]["earth_fraction_of_time_visible"] = float(
            (earth.visible_fraction > 0.0).mean()
        )
        summary["illumination"]["earth_mean_visible_fraction"] = earth.mean_visible_fraction
    timing["illumination_s"] = time.perf_counter() - tick

    # ----------------------------------------------------------------------- thermal
    coords: dict[str, Any] = {
        "time": ("time", period_times.astype("datetime64[ns]")),
        "azimuth": (
            "azimuth",
            np.degrees(structure.azimuth_rad),
            {"units": "degree", "long_name": "grid azimuth (from map +y towards +x)"},
        ),
    }
    model_form: dict[str, Any] = {}
    thermal_nominal = None
    if thermal_on:
        tick = time.perf_counter()
        properties = scenario.thermal.properties
        laws = _envelope_properties(properties)
        grid = scenario.thermal.slope_grid
        orientations = (
            [
                (s, 360.0 * k / grid.aspects)
                for s in grid.slopes_deg
                for k in range(1 if s == 0.0 else grid.aspects)
            ]
            if grid is not None
            else []
        )
        ground_normals = [nominal["normal"]] + [
            normal_from_slope(np.radians(a), np.radians(b)) for a, b in orientations
        ]
        asset_normals = [
            plane_normal(np.radians(o.tilt_deg), np.radians(o.facing_azimuth_deg))
            for o in scenario.thermal.oriented_surfaces
        ]
        diffuse_cases = ALBEDO_CASES if scenario.thermal.albedo_envelope else ("nominal",)
        n_total = eph.et.size
        radiation: tr.SiteRadiation | None = None
        asset_view = "ground"
        tick_diffuse = time.perf_counter()
        if scenario.thermal.terrain_radiation:
            heights = {"ground": 0.0}
            if asset_normals and scenario.site.height_above_surface_m > 0.0:
                heights["asset"] = scenario.site.height_above_surface_m
                asset_view = "asset"
            settings_ill = scenario.illumination.settings
            radiation = tr.site_radiation(
                stack,
                *xy,
                heights,
                eph.et,
                eph.sun,
                eph.r_sun,
                eclipse_factor=nominal["eclipse"],
                refinement=scenario.thermal.terrain_patch_refinement,
                n_strips=settings_ill.n_strips,
                limb_u=settings_ill.limb_darkening,
                mode=settings_ill.disc_mode_code,
            )
            ground_mask = radiation.relevant([("ground", n) for n in ground_normals])
            asset_mask = radiation.relevant([(asset_view, n) for n in asset_normals])
            for case in diffuse_cases:
                radiation.compute(
                    case,
                    laws[case],
                    ground_mask | asset_mask if case == "nominal" else ground_mask,
                    step_s=step_s,
                    spinup_samples=spin,
                    settings=scenario.thermal.settings,
                )
        earth_ground = (
            tr.earthshine(
                nominal["ground"], eph.et, eph.earth, eph.sun, eph.r_earth, projection=PROJ
            )
            if scenario.thermal.earthshine
            else None
        )
        earth_asset = (
            tr.earthshine(structure, eph.et, eph.earth, eph.sun, eph.r_earth, projection=PROJ)
            if scenario.thermal.earthshine and asset_normals
            else None
        )
        timing["diffuse_s"] = time.perf_counter() - tick_diffuse

        def diffuse_parts(
            case: str, normal: np.ndarray, view: str = "ground"
        ) -> tuple[dict[str, np.ndarray], np.ndarray | None]:
            """Diffuse irradiance components on a receiver and the Earth incidence angle."""
            zeros = np.zeros(n_total)
            parts = dict.fromkeys(DIFFUSE_COMPONENTS, zeros)
            if radiation is not None:
                parts["terrain_scattered"], parts["terrain_infrared"] = radiation.irradiance(
                    view, case, normal
                )
            earth = earth_ground if view == "ground" else earth_asset
            incidence = None
            if earth is not None:
                parts["earthshine_reflected"], parts["earthshine_infrared"], incidence = earth.on(
                    normal
                )
            return parts, incidence

        def diffuse_flux(
            case: str, normal: np.ndarray, props: RegolithProperties, *, earth: bool = True
        ) -> np.ndarray:
            parts, incidence = diffuse_parts(case, normal)
            if not earth:
                parts = {
                    **parts,
                    "earthshine_reflected": np.zeros(n_total),
                    "earthshine_infrared": np.zeros(n_total),
                }
            return _diffuse_absorbed(parts, incidence, props)

        diffuse_on = radiation is not None or earth_ground is not None
        direct = absorbed_solar_flux(
            nominal["sun_ground"], nominal["normal"], properties, eclipse_factor=nominal["eclipse"]
        )
        flux = direct + diffuse_flux("nominal", nominal["normal"], properties)
        in_shadow = nominal["sun_ground"].visible_fraction[spin:] <= 0.0
        thermal_nominal = run_thermal(
            flux,
            step_s,
            spin,
            properties,
            scenario.thermal.settings,
            latitude_rad=latitude,
            shadow=in_shadow,
            diffuse_included=radiation is not None,
        )
        warnings.extend(f"thermal: {w}" for w in thermal_nominal.warnings)
        stats = summarize(thermal_nominal)
        slope, facing = slope_and_facing(nominal["normal"])
        bound = neglected_terrain_flux_bound(nominal["ground"].elevation_rad[0], properties)
        # compare the neglected flux with what the surface radiates at its coldest: the
        # night-time balance, not the mean absorbed flux, is what the minimum depends on
        emission_at_min = properties.emissivity * STEFAN_BOLTZMANN_W_M2_K4 * stats.surface_min_k**4
        bound_ratio = bound.flux_w_m2 / emission_at_min
        if radiation is None and bound_ratio > _NEGLECTED_FLUX_WARNING_RATIO:
            warnings.append(
                f"thermal: the neglected terrain flux may reach {bound.flux_w_m2:.3g} W m-2, "
                f"{100.0 * bound_ratio:.0f} % of the {emission_at_min:.3g} W m-2 the surface emits "
                f"at its minimum of {stats.surface_min_k:.1f} K; the minimum is a lower bound and "
                "may be too cold by several kelvin"
            )
        coords["depth"] = (
            "depth",
            thermal_nominal.depth_m,
            {
                "units": "m",
                "long_name": "depth of cell centre below the surface",
                "positive": "down",
            },
        )
        variables["absorbed_solar_flux"] = _variable(
            direct[spin:],
            ("time",),
            "W m-2",
            "direct solar flux absorbed by the surface element",
        )
        if diffuse_on:
            nominal_parts, _ = diffuse_parts("nominal", nominal["normal"])
            coords["component"] = (
                "component",
                np.array(DIFFUSE_COMPONENTS),
                {"units": "1", "long_name": "diffuse irradiance component"},
            )
            variables["diffuse_irradiance"] = _variable(
                np.array([nominal_parts[c][spin:] for c in DIFFUSE_COMPONENTS]),
                ("component", "time"),
                "W m-2",
                "irradiance on the surface element from terrain (single bounce) and Earthshine",
            )
            variables["absorbed_diffuse_flux"] = _variable(
                (flux - direct)[spin:],
                ("time",),
                "W m-2",
                "terrain and Earthshine flux absorbed by the surface element",
            )
        variables["surface_temperature"] = _variable(
            thermal_nominal.surface_temperature_k,
            ("time",),
            "K",
            "temperature of the regolith surface",
        )
        variables["temperature"] = _variable(
            thermal_nominal.temperature_k, ("time", "depth"), "K", "regolith temperature"
        )
        summary["thermal"] = {
            "surface_slope_deg": float(np.degrees(slope)),
            "surface_facing_grid_azimuth_deg": float(np.degrees(facing)),
            "surface_min_k": stats.surface_min_k,
            "surface_max_k": stats.surface_max_k,
            "surface_mean_k": stats.surface_mean_k,
            "lunation_amplitude_mean_k": stats.lunation_amplitude_mean_k,
            "lunation_amplitude_max_k": stats.lunation_amplitude_max_k,
            "max_heating_rate_k_per_h": stats.max_heating_rate_k_per_h,
            "max_cooling_rate_k_per_h": stats.max_cooling_rate_k_per_h,
            "at_depth": {
                f"{d:g} m": {"min_k": float(lo), "max_k": float(hi), "mean_k": float(mean)}
                for d, lo, hi, mean in zip(
                    stats.depth_m,
                    stats.depth_min_k,
                    stats.depth_max_k,
                    stats.depth_mean_k,
                    strict=True,
                )
            },
            "mean_absorbed_flux_w_m2": float(thermal_nominal.absorbed_flux_w_m2.mean()),
            "mean_absorbed_direct_flux_w_m2": float(direct[spin:].mean()),
            "spinup_cycles": thermal_nominal.spinup.cycles,
            "spinup_residual_k": thermal_nominal.spinup.residual_k,
            "neglected_terrain_flux_upper_bound_w_m2": bound.flux_w_m2,
            "shadow_temperature_floor_from_that_bound_k": bound.shadow_temperature_floor_k,
            "emission_at_surface_min_w_m2": emission_at_min,
            "neglected_flux_bound_over_emission_at_min": bound_ratio,
        }
        if scenario.thermal.constant_albedo_alternative:
            flat = properties.model_copy(update={"albedo_a": 0.0, "albedo_b": 0.0})
            flux_flat = absorbed_solar_flux(
                nominal["sun_ground"], nominal["normal"], flat, eclipse_factor=nominal["eclipse"]
            ) + diffuse_flux("nominal", nominal["normal"], flat)
            alternative = run_thermal(flux_flat, step_s, spin, flat, scenario.thermal.settings)
            variables["surface_temperature_constant_albedo"] = _variable(
                alternative.surface_temperature_k,
                ("time",),
                "K",
                "surface temperature with a constant albedo (model-form alternative)",
            )
            alt = summarize(alternative)
            model_form["constant_albedo"] = {
                "description": "same run with albedo independent of incidence angle (albedo_a = albedo_b = 0)",
                "surface_min_k": alt.surface_min_k,
                "surface_max_k": alt.surface_max_k,
                "surface_mean_k": alt.surface_mean_k,
                "difference_surface_mean_k": alt.surface_mean_k - stats.surface_mean_k,
                "difference_surface_max_k": alt.surface_max_k - stats.surface_max_k,
                "status": "not settled: a constant 0.12 lies outside the published albedo "
                "envelope (hot side); against Diviner at the site (8 sunlit summer bins) it is "
                "+5.7 K (TH-V3b) and the envelope's hot bound -0.4 K (TH-V3c)",
            }
        if scenario.thermal.conductivity_alternative:
            other = (
                "hayne2017"
                if properties.conductivity_model == "woods_robinson2019"
                else "woods_robinson2019"
            )
            alt_k = run_thermal(
                flux,
                step_s,
                spin,
                properties.model_copy(update={"conductivity_model": other}),
                scenario.thermal.settings,
            )
            k_stats = summarize(alt_k)
            depth_change = {
                f"{d:g} m": float(alt_mean - mean)
                for d, alt_mean, mean in zip(
                    k_stats.depth_m, k_stats.depth_mean_k, stats.depth_mean_k, strict=True
                )
            }
            model_form["conductivity"] = {
                "description": f"same run with conductivity_model = {other!r} (Woods-Robinson "
                "2019 low-temperature contact term below 150 K, or the Hayne2017 constant one; "
                "review finding F09, TH-C2)",
                "conductivity_model": other,
                "surface_min_k": k_stats.surface_min_k,
                "surface_max_k": k_stats.surface_max_k,
                "surface_mean_k": k_stats.surface_mean_k,
                "difference_surface_min_k": k_stats.surface_min_k - stats.surface_min_k,
                "difference_surface_mean_k": k_stats.surface_mean_k - stats.surface_mean_k,
                "depth_mean_k": {
                    f"{d:g} m": float(m)
                    for d, m in zip(k_stats.depth_m, k_stats.depth_mean_k, strict=True)
                },
                "difference_depth_mean_k": depth_change,
                "status": "not settled: no measurement of lunar regolith conductivity below "
                "100 K is known; the option is a sensitivity, not a refitted model",
            }
            largest = max((abs(v) for v in depth_change.values()), default=0.0)
            if largest > _CONDUCTIVITY_WARNING_K:
                warnings.append(
                    "thermal: mean temperatures at depth change by up to "
                    f"{largest:.1f} K with the other low-temperature conductivity law "
                    f"({other}); treat depth results as bounded by both "
                    "(model_form_sensitivity.conductivity)"
                )
        if scenario.thermal.albedo_envelope:
            envelope_results = {"nominal": thermal_nominal}
            for case, case_properties in _envelope_properties(properties).items():
                if case == "nominal":
                    continue
                case_flux = absorbed_solar_flux(
                    nominal["sun_ground"],
                    nominal["normal"],
                    case_properties,
                    eclipse_factor=nominal["eclipse"],
                ) + diffuse_flux(case, nominal["normal"], case_properties)
                envelope_results[case] = run_thermal(
                    case_flux, step_s, spin, case_properties, scenario.thermal.settings
                )
            coords["albedo_case"] = (
                "albedo_case",
                np.array(ALBEDO_CASES),
                {"units": "1", "long_name": "albedo law: cold bound, scenario law, hot bound"},
            )
            variables["surface_temperature_albedo_envelope"] = _variable(
                np.array([envelope_results[c].surface_temperature_k for c in ALBEDO_CASES]),
                ("albedo_case", "time"),
                "K",
                "surface temperature under the bounds of the published albedo laws",
            )
            summary["thermal"]["albedo_envelope"] = {
                "cold": {
                    "law": ENVELOPE["cold"].name,
                    "source": ENVELOPE["cold"].source,
                    "note": ENVELOPE["cold"].note,
                    **_surface_stats(envelope_results["cold"]),
                },
                "nominal": {
                    "law": "scenario properties (Hayne2017 by default)",
                    **_surface_stats(envelope_results["nominal"]),
                },
                "hot": {
                    "law": ENVELOPE["hot"].name,
                    "source": ENVELOPE["hot"].source,
                    "note": ENVELOPE["hot"].note,
                    **_surface_stats(envelope_results["hot"]),
                },
                "validation": "TH-V3c: Diviner summer temperatures on a 720 m block at the site "
                "against both bounds (VALIDATION.md)",
            }
        if diffuse_on:
            full = {"nominal": thermal_nominal}
            if scenario.thermal.albedo_envelope:
                full.update({c: envelope_results[c] for c in ("cold", "hot")})
            effects: dict[str, Any] = {}
            for case, result in full.items():
                case_props = laws[case]
                direct_case = absorbed_solar_flux(
                    nominal["sun_ground"],
                    nominal["normal"],
                    case_props,
                    eclipse_factor=nominal["eclipse"],
                )
                bare = run_thermal(direct_case, step_s, spin, case_props, scenario.thermal.settings)
                terrain_only = run_thermal(
                    direct_case + diffuse_flux(case, nominal["normal"], case_props, earth=False),
                    step_s,
                    spin,
                    case_props,
                    scenario.thermal.settings,
                )
                at_min = int(np.argmin(result.surface_temperature_k))
                parts, incidence = diffuse_parts(case, nominal["normal"])
                earth_parts = {
                    **dict.fromkeys(("terrain_scattered", "terrain_infrared"), np.zeros(n_total)),
                    "earthshine_reflected": parts["earthshine_reflected"],
                    "earthshine_infrared": parts["earthshine_infrared"],
                }
                earth_absorbed = _diffuse_absorbed(earth_parts, incidence, case_props)
                bare_s, terrain_s, full_s = (
                    _surface_stats(r) for r in (bare, terrain_only, result)
                )
                effects[case] = {
                    "time_of_minimum_utc": _utc(period_times[at_min]),
                    "surface_min_k": full_s["surface_min_k"],
                    "terrain_scattered_at_minimum_w_m2": float(
                        parts["terrain_scattered"][spin + at_min]
                    ),
                    "terrain_infrared_at_minimum_w_m2": float(
                        parts["terrain_infrared"][spin + at_min]
                    ),
                    "terrain_irradiance_at_minimum_w_m2": float(
                        parts["terrain_scattered"][spin + at_min]
                        + parts["terrain_infrared"][spin + at_min]
                    ),
                    "earthshine_absorbed_at_minimum_w_m2": float(earth_absorbed[spin + at_min]),
                    "mean_absorbed_diffuse_w_m2": float(
                        diffuse_flux(case, nominal["normal"], case_props)[spin:].mean()
                    ),
                    "terrain_effect_k": {
                        "surface_min": terrain_s["surface_min_k"] - bare_s["surface_min_k"],
                        "surface_mean": terrain_s["surface_mean_k"] - bare_s["surface_mean_k"],
                    },
                    "earthshine_effect_k": {
                        "surface_min": full_s["surface_min_k"] - terrain_s["surface_min_k"],
                        "surface_mean": full_s["surface_mean_k"] - terrain_s["surface_mean_k"],
                    },
                    "left_out_patches_upper_bound_w_m2": (
                        radiation.left_out_bound_w_m2("ground", case, nominal["normal"])
                        if radiation is not None
                        else None
                    ),
                }
            diffuse_min = float(
                (flux - direct)[spin + int(np.argmin(thermal_nominal.surface_temperature_k))]
            )
            summary["thermal"]["diffuse"] = {
                "method": "single-bounce terrain radiation (Lambertian patches, each with its own "
                "horizon, sunlight and 1-D temperature) and Earthshine (HLS-UG-001 section 11.1); "
                "multiple scattering and sub-DEM roughness not included (docs/models/thermal.md)",
                "terrain_radiation": radiation is not None,
                "earthshine": earth_ground is not None,
                "patches": radiation.patches.size if radiation is not None else 0,
                "patch_points": radiation.patches.n_points if radiation is not None else 0,
                "patches_computed": {
                    c: int(m.sum()) for c, m in (radiation.computed.items() if radiation else [])
                },
                "patch_refinement": scenario.thermal.terrain_patch_refinement,
                "terrain_view_weight": float(radiation.weights("ground", nominal["normal"]).sum())
                if radiation is not None
                else None,
                "absorbed_diffuse_at_minimum_w_m2": diffuse_min,
                "absorbed_diffuse_at_minimum_over_emission": diffuse_min / emission_at_min,
                "cases": effects,
                "validation": "TH-A3 (analytic bowl crater) and TH-C1 (flux at this minimum "
                "against the bound), VALIDATION.md",
            }
            if diffuse_min / emission_at_min > _NEGLECTED_FLUX_WARNING_RATIO:
                warnings.append(
                    f"thermal: terrain radiation and Earthshine supply {diffuse_min:.3g} W m-2 at the "
                    f"surface minimum, {100.0 * diffuse_min / emission_at_min:.0f} % of what the surface "
                    "emits there; they are single-bounce estimates (multiple scattering and roughness "
                    "not included), so the minimum carries their uncertainty"
                )
        if grid is not None:
            cases = ALBEDO_CASES if scenario.thermal.albedo_envelope else ("nominal",)
            table = np.empty((len(orientations), len(cases), 3))
            rows = []
            for i, (slope_deg, facing_deg) in enumerate(orientations):
                normal = normal_from_slope(np.radians(slope_deg), np.radians(facing_deg))
                row: dict[str, Any] = {
                    "slope_deg": slope_deg,
                    "facing_grid_azimuth_deg": facing_deg,
                }
                for j, case in enumerate(cases):
                    case_flux = absorbed_solar_flux(
                        nominal["sun_ground"], normal, laws[case], eclipse_factor=nominal["eclipse"]
                    ) + diffuse_flux(case, normal, laws[case])
                    surface = run_thermal(
                        case_flux, step_s, spin, laws[case], scenario.thermal.settings
                    ).surface_temperature_k
                    table[i, j] = (surface.min(), surface.max(), surface.mean())
                    row[case] = dict(
                        zip(
                            ("surface_min_k", "surface_max_k", "surface_mean_k"),
                            table[i, j].tolist(),
                            strict=True,
                        )
                    )
                rows.append(row)
            coords["orientation"] = (
                "orientation",
                np.arange(len(orientations)),
                {"units": "1", "long_name": "index of the slope grid orientation"},
            )
            coords["orientation_slope"] = (
                "orientation",
                np.array([o[0] for o in orientations]),
                {"units": "degree", "long_name": "ground slope"},
            )
            coords["orientation_facing"] = (
                "orientation",
                np.array([o[1] for o in orientations]),
                {"units": "degree", "long_name": "grid azimuth the slope faces"},
            )
            coords["statistic"] = (
                "statistic",
                np.array(["min", "max", "mean"]),
                {"units": "1", "long_name": "statistic over the period"},
            )
            if "albedo_case" not in coords:
                coords["albedo_case"] = (
                    "albedo_case",
                    np.array(cases),
                    {"units": "1", "long_name": "albedo law"},
                )
            variables["slope_grid_surface_temperature"] = _variable(
                table,
                ("orientation", "albedo_case", "statistic"),
                "K",
                "surface temperature statistics over the period for ground slopes and aspects",
            )
            summary["thermal"]["slope_grid"] = {
                "note": "nominal terrain horizon; slope of the ground under the observer only. The "
                "5 m DEM slope is uncertain by about 2 degrees (Barker2021).",
                "range_surface_max_k": [float(table[:, :, 1].min()), float(table[:, :, 1].max())],
                "range_surface_min_k": [float(table[:, :, 0].min()), float(table[:, :, 0].max())],
                "orientations": rows,
            }
        timing["thermal_s"] = time.perf_counter() - tick

    # --------------------------------------------------------------- asset surfaces
    if scenario.thermal.oriented_surfaces:
        surfaces_out: dict[str, Any] = {}
        irradiance_rows = []
        diffuse_rows: list[list[np.ndarray]] = []
        horizon = nominal["structure"].elevation_rad[0]
        for item in scenario.thermal.oriented_surfaces:
            normal = plane_normal(np.radians(item.tilt_deg), np.radians(item.facing_azimuth_deg))
            irradiance = direct_irradiance(
                sun, normal, eclipse_factor=np.asarray(nominal["eclipse"])[spin:]
            )
            irradiance_rows.append(irradiance)
            factors = view_factors(horizon, structure.azimuth_rad, normal)
            dark = longest_run(irradiance <= 0.0, step_s)
            diffuse_entry: dict[str, Any] = {}
            if thermal_on and diffuse_on:
                parts, _ = diffuse_parts("nominal", normal, asset_view)
                diffuse_rows.append([parts[c][spin:] for c in DIFFUSE_COMPONENTS])
                diffuse_entry = {
                    f"mean_{c}_w_m2": float(parts[c][spin:].mean()) for c in DIFFUSE_COMPONENTS
                } | {f"max_{c}_w_m2": float(parts[c][spin:].max()) for c in DIFFUSE_COMPONENTS}
                if radiation is not None:
                    diffuse_entry["terrain_view_weight"] = float(
                        radiation.weights(asset_view, normal).sum()
                    )
                    diffuse_entry["left_out_patches_upper_bound_w_m2"] = (
                        radiation.left_out_bound_w_m2(asset_view, "nominal", normal)
                    )
            surfaces_out[item.name] = {
                "tilt_deg": item.tilt_deg,
                "facing_grid_azimuth_deg": item.facing_azimuth_deg,
                "mean_direct_irradiance_w_m2": float(irradiance.mean()),
                "max_direct_irradiance_w_m2": float(irradiance.max()),
                "fraction_of_time_with_direct_sun": float((irradiance > 0.0).mean()),
                "longest_without_direct_sun_hours": dark.duration_s / 3600.0,
                "sky_view_factor": factors.sky,
                "ground_view_factor": factors.ground,
                **diffuse_entry,
            }
        coords["surface"] = (
            "surface",
            np.array([s.name for s in scenario.thermal.oriented_surfaces]),
            {"units": "1", "long_name": "asset surface name"},
        )
        variables["oriented_direct_irradiance"] = _variable(
            np.array(irradiance_rows),
            ("surface", "time"),
            "W m-2",
            "direct solar irradiance on the asset surface (no terrain-scattered light)",
        )
        if diffuse_rows:
            if "component" not in coords:
                coords["component"] = (
                    "component",
                    np.array(DIFFUSE_COMPONENTS),
                    {"units": "1", "long_name": "diffuse irradiance component"},
                )
            variables["oriented_diffuse_irradiance"] = _variable(
                np.array(diffuse_rows),
                ("surface", "component", "time"),
                "W m-2",
                "irradiance on the asset surface from terrain (single bounce, scenario albedo "
                "law) and Earthshine",
            )
        summary["oriented_surfaces"] = {
            "note": "At the observer height of the site. Direct sunlight, and (when computed) "
            "single-bounce terrain-scattered sunlight, terrain infrared and Earthshine with the "
            "scenario albedo law for the terrain; the asset does not shadow the ground. "
            "Temperatures of asset surfaces belong in the asset's thermal model "
            "(HLS-UG-001 section 9)."
            if diffuse_rows
            else "Direct sunlight only, at the observer height of the site; terrain-scattered "
            "sunlight, terrain infrared and Earthshine are not included (the ground view factor "
            "scales them). Temperatures of asset surfaces belong in the asset's thermal model "
            "(HLS-UG-001 section 9).",
            "surfaces": surfaces_out,
        }

    # ------------------------------------------------------------------- uncertainty
    uncertainty: dict[str, Any] = {
        "quantile_levels": list(QUANTILES),
        "terrain_clones_used": len(clone_ids),
        "terrain_clones_requested": scenario.uncertainty.terrain_clones,
    }
    tick = time.perf_counter()
    members = [
        _terrain_products(stack.with_near(load_dem(store, clone)), xy, scenario, eph, spin)
        for clone in clone_ids
    ]
    if len(clone_ids) < scenario.uncertainty.terrain_clones:
        warnings.append(
            f"uncertainty: {len(clone_ids)} terrain clones in the data cache, {scenario.uncertainty.terrain_clones} requested"
        )
    if members:
        means = np.array([m["sun"].mean_visible_fraction for m in members])
        shadows = np.array([m["sun"].longest_shadow().duration_s / 3600.0 for m in members])
        uncertainty["illumination"] = {
            "mean_visible_fraction": {
                "nominal": sun.mean_visible_fraction,
                "ensemble_mean": float(means.mean()),
                "ensemble_std": float(means.std(ddof=1)) if means.size > 1 else 0.0,
                "quantiles": np.quantile(means, QUANTILES).tolist(),
                "min": float(means.min()),
                "max": float(means.max()),
            },
            "longest_shadow_hours": {
                "nominal": shadow.duration_s / 3600.0,
                "quantiles": np.quantile(shadows, QUANTILES).tolist(),
            },
            "source": "PGDA terrain error clones of the site DEM; regional terrain not varied",
        }
    n_members = scenario.uncertainty.members if thermal_on else 0
    if n_members:
        property_sets, _ = sample_properties(
            n_members, scenario.uncertainty.seed, base=scenario.thermal.properties
        )
        rng = np.random.default_rng(scenario.uncertainty.seed)
        pick = rng.integers(0, len(members), size=n_members) if members else None
        surfaces = np.empty((n_members, n_period + 1))
        scalars = np.empty((n_members, 4))
        bounds_on = scenario.thermal.albedo_envelope
        bound_cases = ("cold", "hot") if bounds_on else ()
        bound_surfaces = {c: np.empty((n_members, n_period + 1)) for c in bound_cases}
        bound_scalars = {c: np.empty((n_members, 4)) for c in bound_cases}
        for i, member_properties in enumerate(property_sets):
            terrain_member = members[int(pick[i])] if pick is not None else nominal
            member_flux = absorbed_solar_flux(
                terrain_member["sun_ground"],
                terrain_member["normal"],
                member_properties,
                eclipse_factor=terrain_member["eclipse"],
            ) + diffuse_flux("nominal", nominal["normal"], member_properties)
            result = run_thermal(
                member_flux, step_s, spin, member_properties, scenario.thermal.settings
            )
            surfaces[i] = result.surface_temperature_k
            scalars[i] = (
                result.surface_temperature_k.min(),
                result.surface_temperature_k.max(),
                result.surface_temperature_k.mean(),
                result.at_depth(min(1.0, float(result.depth_m[-1]))).mean(),
            )
            # the same member under each bound of the albedo envelope (the law replaces
            # every albedo parameter, the sampled A0 included)
            for case in bound_cases:
                law_properties = _envelope_properties(member_properties)[case]
                law_flux = absorbed_solar_flux(
                    terrain_member["sun_ground"],
                    terrain_member["normal"],
                    law_properties,
                    eclipse_factor=terrain_member["eclipse"],
                ) + diffuse_flux(case, nominal["normal"], law_properties)
                law_result = run_thermal(
                    law_flux, step_s, spin, law_properties, scenario.thermal.settings
                )
                bound_surfaces[case][i] = law_result.surface_temperature_k
                bound_scalars[case][i] = tuple(_surface_stats(law_result).values())
        coords["quantile"] = (
            "quantile",
            np.array(QUANTILES),
            {"units": "1", "long_name": "quantile level over the uncertainty ensemble"},
        )
        variables["surface_temperature_quantiles"] = _variable(
            np.quantile(surfaces, QUANTILES, axis=0),
            ("quantile", "time"),
            "K",
            "quantiles of the surface temperature over the uncertainty ensemble",
        )
        names = ("surface_min_k", "surface_max_k", "surface_mean_k", "mean_at_1m_k")
        if bounds_on:
            per_case = {
                "cold": bound_scalars["cold"],
                "nominal": scalars,
                "hot": bound_scalars["hot"],
            }
            variables["surface_temperature_envelope_quantiles"] = _variable(
                np.array(
                    [
                        np.quantile(
                            surfaces if c == "nominal" else bound_surfaces[c], QUANTILES, axis=0
                        )
                        for c in ALBEDO_CASES
                    ]
                ),
                ("albedo_case", "quantile", "time"),
                "K",
                "quantiles of the surface temperature over the ensemble, per albedo case",
            )
            uncertainty["thermal_envelope"] = {
                "presentation": "Report polar temperatures as these bounds: the lower end is the "
                "5 % quantile of the coldest albedo case, the upper end the 95 % quantile of the "
                "hottest. The ensemble quantiles of one case are the spread within that "
                "scenario, not the uncertainty of the result.",
                "bounds": {
                    name: [
                        float(min(np.quantile(v[:, j], QUANTILES[0]) for v in per_case.values())),
                        float(max(np.quantile(v[:, j], QUANTILES[-1]) for v in per_case.values())),
                    ]
                    for j, name in enumerate(names)
                },
                "cases": {
                    case: {
                        name: {"quantiles": np.quantile(v[:, j], QUANTILES).tolist()}
                        for j, name in enumerate(names)
                    }
                    for case, v in per_case.items()
                },
                "not_included": "multiple scattering, the low-temperature conductivity law "
                "(its effect is reported in model_form_sensitivity.conductivity), sub-DEM "
                "roughness",
            }
        uncertainty["thermal"] = {
            "members": n_members,
            **{
                name: {"quantiles": np.quantile(scalars[:, j], QUANTILES).tolist()}
                for j, name in enumerate(names)
            },
            "source": "Latin hypercube over regolith parameter intervals, each member on a randomly chosen terrain clone"
            if members
            else "Latin hypercube over regolith parameter intervals; terrain not varied (no clones in the cache)",
            "parameter_intervals": {
                r.field: {"low": r.low, "high": r.high, "basis": r.basis} for r in DEFAULT_RANGES
            },
        }
    timing["uncertainty_s"] = time.perf_counter() - tick

    # --------------------------------------------------------------------------- map
    if scenario.illumination.map is not None:
        tick = time.perf_counter()
        block = scenario.illumination.map
        window = near.crop(Bounds.around(xy[0], xy[1], block.half_width_m))
        xx, yy = np.meshgrid(window.grid.x_centers(), window.grid.y_centers())
        horizons = compute_horizons(
            stack,
            xx.ravel(),
            yy.ravel(),
            height_above_surface_m=site.height_above_surface_m,
            n_azimuth=block.n_azimuth,
            diagnostics=False,
        )
        settings = scenario.illumination.settings
        stats_map = illumination_statistics(
            horizons,
            eph.et[spin:],
            eph.sun[spin:],
            eph.r_sun,
            n_strips=settings.n_strips,
            limb_u=settings.limb_darkening,
            mode=settings.disc_mode_code,
        )
        coords["x"] = (
            "x",
            window.grid.x_centers(),
            {"units": "m", "long_name": "south polar stereographic x"},
        )
        coords["y"] = (
            "y",
            window.grid.y_centers(),
            {"units": "m", "long_name": "south polar stereographic y"},
        )
        variables["map_terrain_height"] = _variable(
            window.heights_m.astype(np.float64),
            ("y", "x"),
            "m",
            "terrain height above the 1737.4 km sphere",
        )
        variables["map_mean_visible_fraction"] = _variable(
            stats_map.mean_visible_fraction.reshape(xx.shape),
            ("y", "x"),
            "1",
            "mean visible fraction of the solar disc over the period",
        )
        variables["map_longest_shadow"] = _variable(
            stats_map.longest_shadow_s.reshape(xx.shape) / 3600.0,
            ("y", "x"),
            "h",
            "longest continuous period with no part of the Sun visible",
        )
        variables["map_longest_illumination"] = _variable(
            stats_map.longest_illumination_s.reshape(xx.shape) / 3600.0,
            ("y", "x"),
            "h",
            "longest continuous period with any part of the Sun visible",
        )
        timing["map_s"] = time.perf_counter() - tick

    # ------------------------------------------------------------ dataset and metadata
    dataset = xr.Dataset(variables, coords=coords)
    arrays = {
        str(name): np.asarray(dataset[name].values) for name in dataset.variables if name != "time"
    }
    arrays["time"] = dataset["time"].values.astype("datetime64[s]").astype(np.int64)
    result_checksum = arrays_checksum(arrays)
    dataset.attrs = {
        "title": f"lunarenv run {run_id}: {scenario.name}",
        "Conventions": "CF-1.10",
        "run_id": run_id,
        "result_checksum_sha256": result_checksum,
        "lunarenv_version": __version__,
        "site_x_m": xy[0],
        "site_y_m": xy[1],
        "site_longitude_deg": float(np.degrees(longitude) % 360.0),
        "site_latitude_deg": float(np.degrees(latitude)),
        "observer_height_above_surface_m": site.height_above_surface_m,
        "frame": DE440.body_fixed_frame,
        "comment": "Full provenance, uncertainty and validity information are in metadata.json next to this file.",
    }
    timing["total_s"] = time.perf_counter() - started
    validation = get_validation_status()
    environment_info = environment()
    metadata: dict[str, Any] = {
        "run_id": run_id,
        "result_checksum_sha256": result_checksum,
        "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lunarenv_version": __version__,
        "scenario": scenario.resolved(),
        "site": {
            "x_m": xy[0],
            "y_m": xy[1],
            "longitude_deg": float(np.degrees(longitude) % 360.0),
            "latitude_deg": float(np.degrees(latitude)),
            "terrain_height_m": surface_height,
            "terrain_set": terrain.name,
            "terrain_description": terrain.description,
            "body_fixed_frame": DE440.body_fixed_frame,
            "reference_radius_m": PROJ.radius_m,
        },
        "summary": summary,
        "uncertainty": uncertainty,
        "numerical_error": {
            "illumination": dict(ILLUMINATION_ERROR),
            **({"thermal": dict(THERMAL_ERROR)} if thermal_on else {}),
        },
        "model_form_sensitivity": model_form,
        "validity": {"warnings": warnings, "within_validated_range": not warnings},
        "validation": {
            "generated_utc": validation["generated_utc"],
            "counts": validation["counts"],
            "cases": [
                {"id": c["id"], "status": c["status"], "summary": c["summary"]}
                for c in validation["cases"]
            ],
        },
        "models": [m.describe().model_dump(mode="json") for m in models.values()],
        "data": [store.record(i).model_dump(mode="json") for i in sorted(used_ids)],
        "variables": {
            name: {
                "dims": list(dataset[name].dims),
                "units": dataset[name].attrs.get("units", ""),
                "long_name": dataset[name].attrs.get("long_name", ""),
            }
            for name in dataset.variables
        },
        "environment": environment_info,
        "timing_s": timing,
    }

    directory: Path | None = None
    if output_dir is not None:
        directory = Path(output_dir) / run_id
        directory.mkdir(parents=True, exist_ok=True)
        dataset.to_netcdf(directory / "results.nc", engine="h5netcdf")
        (directory / "metadata.json").write_text(
            json.dumps(metadata, indent=1) + "\n", encoding="utf-8"
        )
        (directory / "scenario.json").write_text(
            json.dumps(scenario.resolved(), indent=1) + "\n", encoding="utf-8"
        )
    return RunResult(run_id, result_checksum, dataset, metadata, directory)


def load_run(directory: str | Path) -> tuple[xr.Dataset, dict[str, Any]]:
    """Dataset and metadata of a run written by :func:`run_scenario`."""
    path = Path(directory)
    metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    return xr.load_dataset(path / "results.nc", engine="h5netcdf"), metadata


__all__ = ["RunResult", "describe_models", "load_run", "run_scenario"]
