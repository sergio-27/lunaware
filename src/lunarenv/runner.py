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
from lunarenv.core.model import ModelDescription, ParameterInfo, register, registered_models
from lunarenv.core.provenance import arrays_checksum, compute_run_id, environment
from lunarenv.data.dem import Bounds, load_dem
from lunarenv.data.store import DataNotAvailableError, DataStore
from lunarenv.geometry.ephemeris import DE440, spice_session
from lunarenv.geometry.projection import SouthPolarStereographic
from lunarenv.illumination.eclipse import earth_occultation_factor
from lunarenv.illumination.horizon import compute_horizons, sun_relevant_incomplete
from lunarenv.illumination.products import illumination_series, illumination_statistics
from lunarenv.illumination.settings import NUMERICAL_ERROR_ESTIMATE as ILLUMINATION_ERROR
from lunarenv.illumination.settings import IlluminationSettings
from lunarenv.illumination.terrain import TerrainStack, build_stack
from lunarenv.scenario.schema import Scenario, load_scenario
from lunarenv.sites import TERRAIN_SETS
from lunarenv.thermal.model import NUMERICAL_ERROR_ESTIMATE as THERMAL_ERROR
from lunarenv.thermal.model import run_thermal, summarize
from lunarenv.thermal.properties import RegolithProperties
from lunarenv.thermal.surface import (
    absorbed_solar_flux,
    neglected_terrain_flux_bound,
    normal_from_slope,
    slope_and_facing,
    surface_normal,
)
from lunarenv.thermal.uq import DEFAULT_RANGES, sample_properties
from lunarenv.validation.core import get_validation_status

PROJ = SouthPolarStereographic()
QUANTILES = (0.05, 0.5, 0.95)
_SPINUP_DAYS = 365


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
    version = "0.1"
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
                "driven by absorbed direct sunlight on the local slope; interior heat flow at the base."
            ),
            requires=self.requires,
            outputs=("absorbed_solar_flux", "surface_temperature", "temperature"),
            assumptions=(
                "Heat flows vertically only; properties depend on depth and temperature.",
                "The surface element has the slope of the DEM at the site unless given explicitly.",
            ),
            simplifications=(
                "No sunlight scattered or heat radiated by surrounding terrain: temperatures in shadow "
                "are lower bounds; an upper bound of the neglected flux is reported with each run.",
                "No rocks, no lateral conduction, no volatiles.",
                "Incidence angle of the disc centre for partially visible discs.",
            ),
            validity_range=(
                "Regolith properties fitted to Diviner data within 70 degrees of the equator; polar use "
                "is an extrapolation.",
                "Heat-capacity fit stated valid from below 90 K to above 400 K.",
                "Albedo law fitted on level equatorial ground; at grazing polar incidence it is the "
                "largest model-form uncertainty (see constant-albedo alternative).",
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
                "Not sampled: the albedo law (reported separately), chi, the heat-capacity polynomial.",
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
            ),
            references=("Hayne2017", "Vasavada2012", "Keihm1984", "heat1d"),
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
        flux = absorbed_solar_flux(
            nominal["sun_ground"], nominal["normal"], properties, eclipse_factor=nominal["eclipse"]
        )
        in_shadow = nominal["sun_ground"].visible_fraction[spin:] <= 0.0
        thermal_nominal = run_thermal(
            flux,
            step_s,
            spin,
            properties,
            scenario.thermal.settings,
            latitude_rad=latitude,
            shadow=in_shadow,
        )
        warnings.extend(f"thermal: {w}" for w in thermal_nominal.warnings)
        stats = summarize(thermal_nominal)
        slope, facing = slope_and_facing(nominal["normal"])
        bound = neglected_terrain_flux_bound(nominal["ground"].elevation_rad[0], properties)
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
            thermal_nominal.absorbed_flux_w_m2,
            ("time",),
            "W m-2",
            "direct solar flux absorbed by the surface element",
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
            "spinup_cycles": thermal_nominal.spinup.cycles,
            "spinup_residual_k": thermal_nominal.spinup.residual_k,
            "neglected_terrain_flux_upper_bound_w_m2": bound.flux_w_m2,
            "shadow_temperature_floor_from_that_bound_k": bound.shadow_temperature_floor_k,
        }
        if scenario.thermal.constant_albedo_alternative:
            flat = properties.model_copy(update={"albedo_a": 0.0, "albedo_b": 0.0})
            flux_flat = absorbed_solar_flux(
                nominal["sun_ground"], nominal["normal"], flat, eclipse_factor=nominal["eclipse"]
            )
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
                "status": "validation has not settled which albedo applies at grazing incidence (TH-V2, TH-V2b)",
            }
        timing["thermal_s"] = time.perf_counter() - tick

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
        for i, member_properties in enumerate(property_sets):
            terrain_member = members[int(pick[i])] if pick is not None else nominal
            member_flux = absorbed_solar_flux(
                terrain_member["sun_ground"],
                terrain_member["normal"],
                member_properties,
                eclipse_factor=terrain_member["eclipse"],
            )
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
