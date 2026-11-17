"""Validation cases of the thermal model."""

from __future__ import annotations

import json
from importlib import resources

import numpy as np

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4 as SIGMA
from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import DE440
from lunarenv.sites import TERRAIN_SETS
from lunarenv.thermal.grid import geometric_grid, uniform_grid
from lunarenv.thermal.model import SYNODIC_MONTH_S, ThermalResult, ThermalSettings, run_thermal
from lunarenv.thermal.properties import RegolithProperties, kirchhoff
from lunarenv.thermal.solver import SchemeName, solve_conduction
from lunarenv.thermal.surface import analytic_absorbed_flux
from lunarenv.validation.core import Case, Check, Outcome, Status

_REFERENCE_REQUIRES = (*TERRAIN_SETS["site04"].required_ids, *DE440.kernel_ids)

_PER_CYCLE = 720


# ------------------------------------------------------------------------- analytic
def _constant_medium() -> tuple[RegolithProperties, float]:
    k, rho, cp = 1.0e-3, 1500.0, 600.0
    properties = RegolithProperties(
        surface_conductivity_w_m_k=k,
        deep_conductivity_w_m_k=k,
        surface_density_kg_m3=rho,
        deep_density_kg_m3=rho,
        radiative_parameter=0.0,
        heat_capacity_coefficients=(cp, 0.0, 0.0, 0.0, 0.0),
        heat_flow_w_m2=0.0,
    )
    return properties, k / (rho * cp)


def _wave_error(cells_per_skin: int, steps: int, scheme: SchemeName = "bdf2") -> float:
    properties, kappa = _constant_medium()
    omega = 2 * np.pi / SYNODIC_MONTH_S
    skin = np.sqrt(2 * kappa / omega)

    def exact(z: np.ndarray, t: np.ndarray) -> np.ndarray:
        return np.asarray(250.0 + 100.0 * np.exp(-z / skin) * np.cos(omega * t - z / skin))

    grid = uniform_grid(skin / cells_per_skin, 12 * skin)
    times = np.arange(steps + 1) * SYNODIC_MONTH_S / steps
    result = solve_conduction(
        grid,
        properties,
        exact(grid.centers_m, np.zeros(1)),
        exact(np.zeros(1), times),
        SYNODIC_MONTH_S / steps,
        scheme=scheme,
        substeps=1,
        boundary="temperature",
    )
    reference = exact(grid.centers_m[None, :], times[:, None])
    return float(np.abs(result.temperature_k - reference).max())


def _th_a1(_: DataStore) -> Outcome:
    error = _wave_error(20, 2000)
    space = [_wave_error(c, 20000) for c in (5, 10, 20)]
    order = float(np.log2(space[0] / space[2]) / 2.0)
    checks = [
        Check(
            quantity="max error, 20 cells per skin depth, 2000 steps per cycle",
            reference=0.0,
            measured=error,
            tolerance=0.25,
            unit="K",
        ),
        Check(
            quantity="observed order of accuracy in space",
            reference=2.0,
            measured=order,
            tolerance=0.3,
        ),
    ]
    return Outcome(checks, f"max error {error:.3f} K of a 100 K wave; spatial order {order:.2f}")


def _th_a2(_: DataStore) -> Outcome:
    q_abs, heat_flow, k = 60.0, 2.0, 1.0e-3
    properties = RegolithProperties(
        surface_conductivity_w_m_k=k,
        deep_conductivity_w_m_k=k,
        surface_density_kg_m3=1500.0,
        deep_density_kg_m3=1500.0,
        radiative_parameter=2.7,
        heat_flow_w_m2=heat_flow,
    )
    grid = geometric_grid(0.004, 1.3, 0.4)
    result = solve_conduction(
        grid,
        properties,
        np.full(grid.n_cells, 200.0),
        np.full(400, q_abs),
        1.0e9,
        scheme="backward_euler",
    )
    surface = ((q_abs + heat_flow) / (properties.emissivity * SIGMA)) ** 0.25
    potential = kirchhoff(result.temperature_k[-1], k, 2.7)
    expected = kirchhoff(surface, k, 2.7) + heat_flow * grid.centers_m
    error = float(np.max(np.abs(potential - expected) / expected))
    checks = [
        Check(
            quantity="surface temperature",
            reference=float(surface),
            measured=float(result.surface_temperature_k[-1]),
            tolerance=1e-6,
            unit="K",
        ),
        Check(
            quantity="relative error of the Kirchhoff potential profile",
            reference=0.0,
            measured=error,
            tolerance=1e-9,
        ),
    ]
    contrast = float(result.temperature_k[-1, -1] - surface)
    return Outcome(
        checks,
        f"profile error {error:.1e} (relative) on {grid.n_cells} cells with {contrast:.0f} K across the domain",
    )


def _th_p1(_: DataStore) -> Outcome:
    properties = RegolithProperties()
    grid = geometric_grid(0.003, 1.15, 1.2)
    phase = 2 * np.pi * np.arange(3 * 240 + 1) / 240
    flux = 1100.0 * np.clip(np.cos(phase), 0.0, None)
    step = SYNODIC_MONTH_S / 240
    checks = []
    worst = 0.0
    schemes: tuple[SchemeName, ...] = ("explicit", "backward_euler", "crank_nicolson", "bdf2")
    for scheme in schemes:
        result = solve_conduction(
            grid, properties, np.full(grid.n_cells, 240.0), flux, step, scheme=scheme
        )
        relative = abs(result.energy_residual_j_m2) / float(flux.sum() * step)
        worst = max(worst, relative)
        checks.append(
            Check(
                quantity=f"relative energy residual, {scheme}",
                reference=0.0,
                measured=relative,
                tolerance=1e-9,
            )
        )
    return Outcome(checks, f"largest relative residual {worst:.1e} over three lunations")


# ------------------------------------------------------------------------ published
def _cycle(
    latitude_deg: float,
    properties: RegolithProperties,
    *,
    distance_au: float = 1.0,
    declination_rad: float = 0.0,
    period_s: float = SYNODIC_MONTH_S,
) -> ThermalResult:
    hour = 2 * np.pi * np.arange(2 * _PER_CYCLE + 1) / _PER_CYCLE
    flux = analytic_absorbed_flux(
        hour,
        np.radians(latitude_deg),
        properties,
        declination_rad=declination_rad,
        sun_distance_au=distance_au,
    )
    return run_thermal(flux, period_s / _PER_CYCLE, _PER_CYCLE, properties, ThermalSettings())


def _equator_checks(properties: RegolithProperties) -> tuple[list[Check], tuple[float, ...]]:
    result = _cycle(0.0, properties)
    surface = result.surface_temperature_k
    values = (float(surface.max()), float(surface[_PER_CYCLE // 2]), float(surface.min()))
    checks = [
        Check(
            quantity="peak noon temperature",
            reference=385.0,
            measured=values[0],
            tolerance=2.0,
            unit="K",
        ),
        Check(
            quantity="midnight temperature",
            reference=101.0,
            measured=values[1],
            tolerance=2.0,
            unit="K",
        ),
        Check(
            quantity="minimum night-time temperature",
            reference=95.0,
            measured=values[2],
            tolerance=2.0,
            unit="K",
        ),
    ]
    return checks, values


def _th_v1(_: DataStore) -> Outcome:
    checks, values = _equator_checks(RegolithProperties())
    notes = (
        "Reference values are those quoted in Hayne2017 Table A2 from Vasavada2012 (Diviner, equator, "
        "A0 = 0.12), described there as typical and rounded to 1 K. Vasavada2012 itself was not accessed.",
        "Level ground at the equator, Sun at 1 au, declination 0, standard parameters (H = 0.06 m).",
        "Disclosure: the 2 K tolerance was written down after exploratory runs had already shown "
        "385.2 / 100.3 / 93.7 K. It is anchored to the residual stated by the source, not to the result.",
        "This case does not discriminate between albedo laws: with a constant albedo of 0.12 the "
        "same checks also pass (TH-V1b, pre-registered). It is therefore no evidence for the "
        "incidence-dependent albedo used at the pole.",
    )
    return Outcome(
        checks, "{:.1f} / {:.1f} / {:.1f} K against 385 / 101 / 95 K".format(*values), notes
    )


def _th_v1b(_: DataStore) -> Outcome:
    checks, values = _equator_checks(RegolithProperties(albedo_a=0.0, albedo_b=0.0))
    notes = (
        "Same as TH-V1 with a constant albedo of 0.12 (albedo_a = albedo_b = 0). Pre-registered "
        "in docs/review/preregistration-2026-11-13.md (finding F04 of the 2026-11-13 scientific "
        "review) before it was implemented or run; reference and tolerances are those of TH-V1.",
        "Reading fixed beforehand: a PASS means TH-V1 does not discriminate between the "
        "incidence-dependent and the constant albedo, so it is not evidence for the albedo law "
        "used at the pole; a miss by more than 2 K means it does.",
    )
    return Outcome(
        checks,
        "{:.1f} / {:.1f} / {:.1f} K against 385 / 101 / 95 K [constant albedo]".format(*values),
        notes,
    )


_APOLLO = (
    ("Apollo 17 (20 N)", 20.0, 216.0, 1.30, 256.0),
    ("Apollo 15 (26 N)", 26.0, 211.0, 0.83, 252.0),
)


def _apollo(properties: RegolithProperties) -> tuple[list[Check], str]:
    checks = []
    parts = []
    for name, latitude, surface_ref, depth, deep_ref in _APOLLO:
        result = _cycle(latitude, properties, distance_au=1.017)
        surface = float(result.surface_temperature_k[:-1].mean())
        deep = float(result.at_depth(depth)[:-1].mean())
        checks.append(
            Check(
                quantity=f"{name}: diurnal mean surface temperature",
                reference=surface_ref,
                measured=surface,
                tolerance=5.0,
                unit="K",
            )
        )
        checks.append(
            Check(
                quantity=f"{name}: diurnal mean temperature at {depth:.2f} m",
                reference=deep_ref,
                measured=deep,
                tolerance=5.0,
                unit="K",
            )
        )
        parts.append(f"{surface:.1f} / {deep:.1f} K")
    return checks, "; ".join(parts) + " against 216 / 256 and 211 / 252 K (+/-5 K)"


_APOLLO_NOTES = (
    "Reference: Apollo 15 and 17 heat-flow experiment values (Keihm1973) as quoted, with their +/-5 K, "
    "in Hayne2017 Table A2. The original papers were not accessed, so what the 'mean surface "
    "temperature' was measured on (regolith surface or probe hardware) is not established here.",
    "Sun at 1.017 au, as stated in the caption of Hayne2017 Figure A1.",
)


def _th_v2(_: DataStore) -> Outcome:
    checks, summary = _apollo(RegolithProperties())
    notes = (
        *_APOLLO_NOTES,
        "With the published parameter set, including the incidence-dependent albedo of Hayne2017 "
        "Table A1, the model is 8-10 K below all four values while it reproduces the equatorial "
        "noon, midnight and minimum temperatures (TH-V1) and agrees with the reference code at the "
        "surface (TH-X1). The cause is not resolved; see TH-V2b. No tolerance was changed.",
    )
    return Outcome(checks, summary, notes)


def _th_v2b(_: DataStore) -> Outcome:
    checks, summary = _apollo(RegolithProperties(albedo_a=0.0, albedo_b=0.0))
    notes = (
        *_APOLLO_NOTES,
        "Same model with a constant albedo of 0.12. It then also matches the model curves of "
        "Hayne2017 Figure A1 (surface ~215 K, 1 m ~256 K at 20 N). An earlier note read this as "
        "the figure having been computed with a constant albedo; the caption says 'All model "
        "parameters are from Table A1', which includes the incidence-dependent law, so that "
        "reading is not supported (checked in the arXiv version, 2026-11-13). It does not show "
        "which albedo is right: daytime Diviner temperatures are what the incidence-dependent law "
        "was fitted to.",
        "The Apollo comparison is also not like-for-like: per the 2026-11-13 scientific review, "
        "the 216 K of Apollo 17 is a brightness temperature (emissivity 1) from cable "
        "thermocouples above the surface, and both sites are dark mare. Not yet checked here.",
    )
    return Outcome(checks, summary + " [constant albedo]", notes)


def _heat1d_checks(
    file: str, grid: list[int] | None = None
) -> tuple[list[Check], list[float], bool, dict]:  # type: ignore[type-arg]
    text = resources.files("lunarenv.validation").joinpath(f"reference/{file}").read_text()
    reference = json.loads(text)
    planet = reference["planet"]
    properties = RegolithProperties(scale_height_m=planet["H"])
    checks: list[Check] = []
    surface_ok = True
    deep_difference = []
    for case in reference["cases"]:
        if case["equilibration_years"] != 20:
            continue
        if grid is not None and case.get("grid_m_n") != grid:
            continue
        result = _cycle(
            case["latitude_deg"],
            properties,
            distance_au=case["sun_distance_au"],
            declination_rad=case["declination_rad"],
            period_s=planet["day_s"],
        )
        surface = result.surface_temperature_k[:-1]
        label = f"latitude {case['latitude_deg']:.0f}"
        for name, measured in (
            ("mean", surface.mean()),
            ("max", surface.max()),
            ("min", surface.min()),
        ):
            check = Check(
                quantity=f"{label}: surface {name}",
                reference=case[f"surface_{name}_k"],
                measured=float(measured),
                tolerance=1.0,
                unit="K",
            )
            surface_ok = surface_ok and bool(check.passed)
            checks.append(check)
        deep = float(result.at_depth(case["bottom_depth_m"])[:-1].mean())
        deep_difference.append(deep - case["bottom_mean_k"])
        checks.append(
            Check(
                quantity=f"{label}: mean at {case['bottom_depth_m']:.2f} m",
                reference=case["bottom_mean_k"],
                measured=deep,
                tolerance=1.0,
                unit="K",
            )
        )
    return checks, deep_difference, surface_ok, reference


def _heat1d_status(checks: list[Check], surface_ok: bool) -> Status | None:
    if surface_ok and not all(bool(c.passed) for c in checks):
        return Status.PARTIAL
    return None


def _th_x1(_: DataStore) -> Outcome:
    checks, deep_difference, surface_ok, reference = _heat1d_checks("heat1d_reference.json")
    one_year = next(
        c for c in reference["cases"] if c["equilibration_years"] == 1 and c["latitude_deg"] == 20.0
    )
    notes = (
        "Reference: heat1d 0.3.2 (the public code of Hayne2017), run on 2026-11-12 with its lunar "
        "defaults except the equilibration time, raised from 1 to 20 years; inputs and outputs are "
        "stored in validation/reference/heat1d_reference.json with the script that produced them.",
        "lunarenv is run with heat1d's parameters (H = 0.07 m, its Sun distance and declination at "
        "output time). Tolerance 1 K: the level to which Hayne2017 states its grids are converged. "
        "Disclosure: it was written down after exploratory runs of both codes had been compared.",
        f"At depth lunarenv is {np.mean(deep_difference):+.1f} K from heat1d. lunarenv's deep "
        "temperature changes by 0.013 K between grids with first cells from 3.6 mm to 0.25 mm; heat1d "
        "returns NaN on every grid finer than its default, so its value cannot be refined. The gap is "
        "attributed to heat1d's averaging of conductivities across the near-surface gradient "
        "(a hypothesis consistent with the evidence, not demonstrated).",
        f"With its default one-year equilibration heat1d's deep temperature at 20 N is "
        f"{one_year['bottom_mean_k']:.1f} K, 9 K above its own equilibrated value.",
        "Superseded by TH-X1b: the heat1d README says the PyPI release used here is obsolete.",
    )
    return Outcome(
        checks,
        f"surface within 1 K at three latitudes; at depth {np.mean(deep_difference):+.1f} K",
        notes,
        _heat1d_status(checks, surface_ok),
    )


def _th_x1b(_: DataStore) -> Outcome:
    checks, deep_difference, surface_ok, reference = _heat1d_checks("heat1d_reference_current.json")
    auto = {c["latitude_deg"]: c for c in reference["cases"] if c["equilibration_years"] == 0}
    notes = (
        f"Reference: heat1d {reference['version']} from GitHub (commit {reference['commit']}), "
        "the current code of Hayne2017, run on 2026-11-15 with its lunar defaults, default grid "
        f"and solver, equilibration 20 years; stored in reference/heat1d_reference_current.json, "
        "produced by benchmarks/make_heat1d_reference.py in an isolated environment.",
        "Pre-registered in docs/review/preregistration-th-x1b-2026-11-15.md (finding F11 of the "
        "2026-11-13 scientific review) before the reference was generated. That file left the "
        "equilibration choice ambiguous (this version has both an automatic and a fixed "
        "setting); 20 years was chosen before the run, as TH-X1, and the automatic setting is "
        "stored for information: at 20 N its deep mean is "
        f"{auto[20.0]['bottom_mean_k']:.2f} K against {next(c for c in reference['cases'] if c['equilibration_years'] == 20 and c['latitude_deg'] == 20.0)['bottom_mean_k']:.2f} K.",
        "lunarenv is run with heat1d's parameters (scale height, Sun distance and declination "
        "at output time, day length), as in TH-X1; same checks and 1 K tolerance.",
        "Disclosure: the first run (2026-11-15) failed with surface means 7-12 K above "
        "lunarenv. Investigation showed heat1d 0.4.1 writes its output at non-uniform "
        "local-time steps (0.012 to 0.024 h, denser by day), and the reference script took a "
        "plain mean of the samples, biased warm. The reference means are now time-weighted "
        "(the plain sample means are kept in the file). This corrects how the reference is "
        "measured; no tolerance or check was changed.",
        "Result as pre-registered: FAIL. Surface maxima agree within 0.6 K and the deep means "
        "within 1 K (the 2.3 K depth gap of TH-X1 came from the obsolete release), but "
        "lunarenv's night minima are 1.5-1.6 K and its means 1.1-1.6 K warmer than heat1d's "
        "default run.",
        "Investigation (post hoc, after the result; not part of the pre-registered checks), "
        "equator, 20 years: heat1d's solver settings move the minimum by at most 0.3 K "
        "(explicit, Crank-Nicolson, implicit, fixed step), but its grid moves it by about "
        "2 K. Its default grid (m=10 layers per skin depth, n=4, i.e. thickness growth 1.25; "
        "about 20 layers) is not converged: with m=n=10, 20, 40 (33, 64, 125 layers) heat1d "
        "gives minimum 93.36, 93.25, 93.19 K, mean 211.53, 211.45, 211.42 K and maximum "
        "388.51, 388.50, 388.49 K, against lunarenv's 93.14, 211.75 and 388.05 K. Against "
        "the converged heat1d the differences are within 0.5 K. The failure is read as "
        "discretisation error of heat1d's default grid, not as a difference in physics; a "
        "pre-registered comparison against a converged heat1d grid at all three latitudes "
        "would be needed to replace this case.",
        "Follow-up: TH-X1c (pre-registered 2026-11-16) makes that comparison and passes; "
        "TH-X1b stays recorded as FAIL.",
    )
    return Outcome(
        checks,
        f"surface max |diff| {max(abs(c.measured - c.reference) for c in checks if 'surface' in c.quantity):.2f} K "  # type: ignore[operator]
        f"at three latitudes; at depth {np.mean(deep_difference):+.1f} K",
        notes,
        _heat1d_status(checks, surface_ok),
    )


def _th_x1c(_: DataStore) -> Outcome:
    checks, deep_difference, surface_ok, reference = _heat1d_checks(
        "heat1d_reference_converged.json", grid=[40, 40]
    )
    by_grid: dict[tuple[float, int], dict] = {  # type: ignore[type-arg]
        (c["latitude_deg"], c["grid_m_n"][0]): c for c in reference["cases"]
    }
    for lat in (0.0, 20.0, 26.0):
        coarse, fine = by_grid[(lat, 20)], by_grid[(lat, 40)]
        change = max(
            abs(fine[f"surface_{q}_k"] - coarse[f"surface_{q}_k"]) for q in ("mean", "max", "min")
        )
        checks.append(
            Check(
                quantity=f"latitude {lat:.0f}: heat1d surface change, grid m=n=20 to 40",
                reference=0.0,
                measured=change,
                tolerance=0.25,
                unit="K",
                informational=True,
            )
        )
    notes = (
        f"Reference: heat1d {reference['version']} from GitHub (commit {reference['commit']}) "
        "on refined grids m=n=20 and m=n=40 (the reference), default solver, 20-year "
        "equilibration, run on 2026-11-16; stored in reference/heat1d_reference_converged.json, "
        "produced by benchmarks/make_heat1d_reference.py in an isolated environment.",
        "Pre-registered in docs/review/preregistration-th-x1c-2026-11-16.md after TH-X1b "
        "failed and before the refined reference at 20 and 26 deg was generated. Not blind at "
        "the equator: heat1d's refined values there were known from the TH-X1b investigation.",
        "Same checks, parameters and 1 K tolerance as TH-X1b. Informational: the change in "
        "heat1d's surface values from m=n=20 to m=n=40; above 0.25 K the reference would be "
        "reported as not converged.",
    )
    return Outcome(
        checks,
        f"surface max |diff| {max(abs(c.measured - c.reference) for c in checks if 'surface ' in c.quantity and 'heat1d' not in c.quantity):.2f} K "  # type: ignore[operator]
        f"at three latitudes; at depth {np.mean(deep_difference):+.1f} K",
        notes,
        _heat1d_status([c for c in checks if not c.informational], surface_ok),
    )


# ------------------------------------------------- published: Diviner polar temperatures
_TH_V3_CACHE: dict[str, dict] = {}  # type: ignore[type-arg]
_TH_V3_TOLERANCE_K = 10.0


def _th_v3_data(store: DataStore) -> dict:  # type: ignore[type-arg]
    """Forward model and Diviner data, computed once per data store for TH-V3 and TH-V3b."""
    from lunarenv.validation import diviner_pcp as dp

    key = str(store.root)
    if key in _TH_V3_CACHE:
        return _TH_V3_CACHE[key]
    reference = dp.load_reference()
    keys = [c["key"] for c in reference["cells"]]
    centres = [(float(c["x_m"]), float(c["y_m"])) for c in reference["cells"]]
    cells = dp.forward(store, centres, 20.0)
    tx, ty = reference["test_point_m"]
    pixel = float(reference["pixel_m"])
    middle = keys.index(f"{round(tx / pixel)},{round(ty / pixel)}")  # cell holding the test point
    fine = dp.forward(store, [centres[middle]], 10.0)[0]
    diviner = np.full((len(keys), dp.N_BINS), np.nan)
    for b, entry in enumerate(reference["bins"]):
        for c, k in enumerate(keys):
            if k in entry["tbol_k"]:
                diviner[c, b] = entry["tbol_k"][k]
    _TH_V3_CACHE[key] = {
        "keys": keys,
        "cells": cells,
        "middle": middle,
        "fine": fine,
        "diviner": diviner,
        "reference": reference,
    }
    return _TH_V3_CACHE[key]


def _tbol(radiance_sum: np.ndarray, count: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.asarray((radiance_sum / count) ** 0.25)


def _th_v3_law(store: DataStore, law: str) -> Outcome:
    d = _th_v3_data(store)
    cells, diviner = d["cells"], d["diviner"]
    model_cell = np.array([_tbol(c.sums[law], c.sums["count"]) for c in cells])  # (cell, bin)
    visible = np.mean([c.sums["visible"] / c.sums["count"] for c in cells], axis=0)
    model_block = np.mean(model_cell**4, axis=0) ** 0.25
    complete = np.all(np.isfinite(diviner), axis=0)
    diviner_block = np.where(
        complete, np.mean(np.where(complete, diviner, 0.0) ** 4, axis=0) ** 0.25, np.nan
    )
    sunlit = visible >= 0.5
    used = sunlit & complete
    every = complete
    diff = model_block - diviner_block
    bias = float(np.mean(diff[used])) if used.any() else float("nan")
    rms = float(np.sqrt(np.mean(diff[used] ** 2))) if used.any() else float("nan")
    checks = [
        Check(
            quantity="bias of block Tbol over sunlit summer bins (model - Diviner)",
            reference=0.0,
            measured=bias,
            tolerance=_TH_V3_TOLERANCE_K,
            unit="K",
        ),
        Check(
            quantity="RMS over sunlit bins",
            reference=0.0,
            measured=rms,
            unit="K",
            informational=True,
        ),
        Check(
            quantity="bias over all summer bins with data",
            reference=0.0,
            measured=float(np.mean(diff[every])) if every.any() else None,
            unit="K",
            informational=True,
        ),
        Check(
            quantity="sunlit bins used / bins with all nine cells",
            reference=float(every.sum()),
            measured=float(used.sum()),
            informational=True,
        ),
    ]
    if every.any():
        hottest = int(np.nanargmax(np.where(every, diviner_block, -np.inf)))
        checks.append(
            Check(
                quantity=f"bin with the highest Diviner block Tbol (sunlon {3.75 * hottest:.2f} deg): Diviner / model",
                reference=float(diviner_block[hottest]),
                measured=float(model_block[hottest]),
                unit="K",
                informational=True,
            )
        )
    for c, k in enumerate(d["keys"]):
        ok = sunlit & np.isfinite(diviner[c])
        checks.append(
            Check(
                quantity=f"cell {k} (240 m units): bias over sunlit bins",
                reference=0.0,
                measured=float(np.mean(model_cell[c, ok] - diviner[c, ok])) if ok.any() else None,
                unit="K",
                informational=True,
            )
        )
    # spread between individual summers of the model's sunlit-bin block Tbol
    n_seasons = cells[0].sums["season_count"].size // model_cell.shape[1]
    season_rad = np.mean(
        [(c.sums[law + "_season"] / c.sums["season_count"]).reshape(n_seasons, -1) for c in cells],
        axis=0,
    )
    with np.errstate(invalid="ignore"):
        season_tbol = season_rad**0.25
        spread = np.nanstd(season_tbol[:, sunlit], axis=0)
    checks.append(
        Check(
            quantity="spread between summers of model sunlit-bin block Tbol (median of per-bin std)",
            measured=float(np.nanmedian(spread)) if sunlit.any() else None,
            unit="K",
            informational=True,
        )
    )
    fine = d["fine"]
    coarse = cells[d["middle"]]
    fine_t = _tbol(fine.sums[law], fine.sums["count"])
    coarse_t = _tbol(coarse.sums[law], coarse.sums["count"])
    checks.append(
        Check(
            quantity="centre cell sunlit-bin mean Tbol, pixels every 20 m / every 10 m",
            reference=float(np.mean(coarse_t[sunlit])) if sunlit.any() else None,
            measured=float(np.mean(fine_t[sunlit])) if sunlit.any() else None,
            unit="K",
            informational=True,
        )
    )
    reference = d["reference"]
    notes = (
        "Reference: Diviner Polar Cumulative Product, diurnal, summer, south pole "
        f"({reference['source']}), 3 x 3 cells of 240 m centred on the cell containing the test "
        "point, extracted on "
        + reference["retrieved_utc"]
        + " by benchmarks/extract_diviner_pcp.py. "
        "Radiances are binned before conversion to bolometric temperature (Diviner RDR SIS, "
        "section 9.1).",
        "Forward model: every 20 m inside each cell (144 pixels per cell) on the 5 m site DEM, "
        "ground-level horizon, DEM facet normal, Earth occultation, default regolith "
        "parameters (emissivity 0.95), hourly 2009-07-05 to 2019-02-17 after a one-year "
        "spin-up; summer samples (subsolar latitude < 0) binned by subsolar longitude; "
        "Tbol = (mean eps T^4)^(1/4) over pixels and samples; block = radiance mean of the "
        "nine cells. Sunlit bins: model block mean visible solar fraction >= 0.5.",
        "Pre-registered in docs/review/preregistration-th-v3-2026-11-15.md (finding F01 of the "
        "2026-11-13 scientific review) before implementation and before any Diviner value at "
        "the site was looked at. The 10 K tolerance is not a Diviner error budget (none was "
        "found in the sources read); it is set so the two albedo laws, about 26 K apart in the "
        "yearly mean, can be told apart.",
        "Not modelled: terrain-scattered sunlight and terrain emission, sub-DEM roughness, the "
        "timing of Diviner observations inside each bin (the model averages every summer hour), "
        "the Diviner footprint beyond the 240 m cell, spectral emissivity.",
    )
    status = None if used.any() else Status.PARTIAL
    summary = (
        f"bias {bias:+.1f} K, RMS {rms:.1f} K over {int(used.sum())} sunlit summer bins"
        if used.any()
        else "no sunlit bin with data in all nine cells"
    )
    return Outcome(checks, summary, notes, status)


def _th_v3(store: DataStore) -> Outcome:
    return _th_v3_law(store, "hayne")


def _th_v3b(store: DataStore) -> Outcome:
    outcome = _th_v3_law(store, "constant")
    return Outcome(
        outcome.checks, outcome.summary + " [constant albedo]", outcome.notes, outcome.status
    )


def _block_bias(store: DataStore, law: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Model minus Diviner block Tbol per bin, the sunlit-bin mask and the complete mask."""
    d = _th_v3_data(store)
    cells, diviner = d["cells"], d["diviner"]
    model_cell = np.array([_tbol(c.sums[law], c.sums["count"]) for c in cells])
    visible = np.mean([c.sums["visible"] / c.sums["count"] for c in cells], axis=0)
    complete = np.all(np.isfinite(diviner), axis=0)
    diviner_block = np.where(
        complete, np.mean(np.where(complete, diviner, 0.0) ** 4, axis=0) ** 0.25, np.nan
    )
    model_block = np.mean(model_cell**4, axis=0) ** 0.25
    return model_block - diviner_block, (visible >= 0.5) & complete, complete


def _block_from_radiance(
    cell_rad: np.ndarray, diviner: np.ndarray, visible: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Model block Tbol, model - Diviner per bin, sunlit-and-complete and complete masks."""
    with np.errstate(invalid="ignore"):
        model_cell = cell_rad**0.25
    complete = np.all(np.isfinite(diviner), axis=0)
    diviner_block = np.where(
        complete, np.mean(np.where(complete, diviner, 0.0) ** 4, axis=0) ** 0.25, np.nan
    )
    model_block = np.mean(model_cell**4, axis=0) ** 0.25
    return model_block, model_block - diviner_block, (visible >= 0.5) & complete, complete


def _th_v3d(store: DataStore) -> Outcome:
    from lunarenv.validation import diviner_pcp as dp

    d = _th_v3_data(store)
    cells, diviner = d["cells"], d["diviner"]
    centres = [(float(c["x_m"]), float(c["y_m"])) for c in d["reference"]["cells"]]
    receivers = dp.forward_diffuse(store, centres)
    visible = np.mean([c.sums["visible"] / c.sums["count"] for c in cells], axis=0)

    def change(key: str) -> np.ndarray:  # (cell, bin) mean radiance change over receivers
        return np.array([np.mean([r[key] / r["count"] for r in rs], axis=0) for rs in receivers])

    def base(law: str) -> np.ndarray:
        return np.array([c.sums[law] / c.sums["count"] for c in cells])

    results = {}
    for law in dp.DIFFUSE_LAWS:
        before = _block_from_radiance(base(law), diviner, visible)
        after = _block_from_radiance(base(law) + change(law + "_change"), diviner, visible)
        results[law] = (before, after)
    hayne_with = base("hayne") + change("hayne_change")
    low_t = _block_from_radiance(hayne_with + change("hayne_lowt_change"), diviner, visible)
    (_, diff0, used, every), (block, diff, _, _) = results["hayne"]
    (_, hot0, _, _), (_, hot, _, _) = results["foote2020"]

    def mean(x: np.ndarray, mask: np.ndarray) -> float:
        return float(np.mean(x[mask])) if mask.any() else float("nan")

    bias = mean(diff, used)
    checks = [
        Check(
            quantity="bias of block Tbol over sunlit summer bins (model - Diviner), "
            "Hayne2017 law with terrain radiation and Earthshine",
            reference=0.0,
            measured=bias,
            tolerance=_TH_V3_TOLERANCE_K,
            unit="K",
        ),
        Check(
            quantity="RMS over sunlit bins",
            reference=0.0,
            measured=float(np.sqrt(np.mean(diff[used] ** 2))) if used.any() else None,
            unit="K",
            informational=True,
        ),
        Check(
            quantity="bias over all summer bins with data",
            reference=0.0,
            measured=mean(diff, every),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="change of the sunlit-bin bias from TH-V3 (diffuse terms)",
            reference=0.0,
            measured=bias - mean(diff0, used),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="change of the all-bin bias from TH-V3 (diffuse terms)",
            reference=0.0,
            measured=mean(diff, every) - mean(diff0, every),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="hot law (Foote 2020) with the diffuse terms: sunlit-bin bias",
            reference=0.0,
            measured=mean(hot, used),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="hot law with the diffuse terms: all-bin bias",
            reference=0.0,
            measured=mean(hot, every),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="hot law: change of the sunlit-bin / all-bin bias from TH-V3c",
            reference=mean(hot, used) - mean(hot0, used),
            measured=mean(hot, every) - mean(hot0, every),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="low-temperature conductivity (Hayne2017, with diffuse): change of the "
            "sunlit-bin / all-bin block Tbol",
            reference=mean(low_t[0] - block, used),
            measured=mean(low_t[0] - block, every),
            unit="K",
            informational=True,
        ),
    ]
    flat = [r for rs in receivers for r in rs]
    for law, label in (("hayne", "Hayne2017"), ("foote2020", "hot law")):
        for part, name in (("terrain", "terrain radiation"), ("earth", "Earthshine")):
            values = np.array([float(r[f"{law}_{part}_abs"][0]) for r in flat])
            checks.append(
                Check(
                    quantity=f"{label}: absorbed {name}, summer mean over the {len(flat)} "
                    "receivers / largest receiver",
                    reference=float(values.mean()),
                    measured=float(values.max()),
                    unit="W m-2",
                    informational=True,
                )
            )
    checks.append(
        Check(
            quantity="largest bound of the flux from patches not computed (Hayne2017)",
            measured=max(float(r["left_out_bound"][0]) for r in flat),
            unit="W m-2",
            informational=True,
        )
    )
    with np.errstate(invalid="ignore"):
        before_cell = base("hayne") ** 0.25
        after_cell = hayne_with**0.25
    for c, k in enumerate(d["keys"]):
        ok = np.isfinite(diviner[c])
        delta = after_cell[c] - before_cell[c]
        checks.append(
            Check(
                quantity=f"cell {k}: change of mean Tbol, sunlit bins / all bins with data",
                reference=mean(delta, (visible >= 0.5) & ok),
                measured=mean(delta, ok),
                unit="K",
                informational=True,
            )
        )
    notes = (
        "Pre-registered in docs/review/preregistration-th-v3d-2026-11-16.md before it was "
        "implemented or run; TH-V3, TH-V3b and TH-V3c stay as recorded.",
        "Base: the TH-V3 forward model (144 pixels per cell). Added: the change of the summer "
        "radiance caused by single-bounce terrain radiation (engine method, default patch "
        "settings, patch temperatures with the same law) and Earthshine (HLS-UG-001), "
        "computed at 4 receivers per cell (+/-60 m from the centre) and averaged per cell. "
        "Sunlit bins, data and tolerance as TH-V3.",
        "The subsample estimates each cell's change, not its temperature. Not modelled: "
        "multiple scattering, sub-DEM roughness, the Diviner footprint and observation "
        "timing, spectral emissivity.",
    )
    summary = (
        f"bias {bias:+.1f} K over {int(used.sum())} sunlit summer bins "
        f"({bias - mean(diff0, used):+.2f} K from TH-V3); all summer bins "
        f"{mean(diff, every):+.1f} K ({mean(diff, every) - mean(diff0, every):+.2f} K) "
        "[with terrain radiation and Earthshine]"
    )
    return Outcome(checks, summary, notes, None if used.any() else Status.PARTIAL)


def _th_v3c(store: DataStore) -> Outcome:
    cold, used, complete = _block_bias(store, "burger2024")
    hot, _, _ = _block_bias(store, "foote2020")
    middle, _, _ = _block_bias(store, "hayne")
    inside = (cold <= 0.0) & (hot >= 0.0)
    cold_bias, hot_bias = float(np.mean(cold[used])), float(np.mean(hot[used]))
    checks = [
        Check(
            quantity="Diviner colder than the cold law (Bürger 2024), sunlit-bin mean; 0 if not",
            reference=0.0,
            measured=max(cold_bias, 0.0),
            tolerance=_TH_V3_TOLERANCE_K,
            unit="K",
        ),
        Check(
            quantity="Diviner warmer than the hot law (Foote 2020), sunlit-bin mean; 0 if not",
            reference=0.0,
            measured=max(-hot_bias, 0.0),
            tolerance=_TH_V3_TOLERANCE_K,
            unit="K",
        ),
        Check(
            quantity="bias of the cold law over sunlit bins (model - Diviner)",
            measured=cold_bias,
            unit="K",
            informational=True,
        ),
        Check(
            quantity="bias of the hot law over sunlit bins (model - Diviner)",
            measured=hot_bias,
            unit="K",
            informational=True,
        ),
        Check(
            quantity="bias of Hayne 2017 (middle) over sunlit bins",
            measured=float(np.mean(middle[used])),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="sunlit bins with Diviner between the two laws / sunlit bins",
            reference=float(used.sum()),
            measured=float(inside[used].sum()),
            informational=True,
        ),
        Check(
            quantity="summer bins with Diviner between the two laws / bins with all nine cells",
            reference=float(complete.sum()),
            measured=float(inside[complete].sum()),
            informational=True,
        ),
        Check(
            quantity="bias over all summer bins with data, cold law",
            measured=float(np.mean(cold[complete])),
            unit="K",
            informational=True,
        ),
        Check(
            quantity="bias over all summer bins with data, hot law",
            measured=float(np.mean(hot[complete])),
            unit="K",
            informational=True,
        ),
    ]
    notes = (
        "Pre-registered in docs/review/preregistration-th-v3c-2026-11-16.md (change 2 of the "
        "2026-11-13 scientific review) before either envelope law was run at the site.",
        "Same data, cells, forward model and sunlit-bin selection as TH-V3; the two envelope "
        "laws were added to the forward model's law list. Cold bound: Bürger2024 Eq. 15, A0 "
        "0.09, a 0.10, b 0.02 (0.860 at 88.2 deg). Hot bound: Foote2020 Eq. 9, Apollo 16 "
        "highland soil, laboratory hemispheric albedo measured to 75 deg and extrapolated "
        "beyond (0.238 at 88.2 deg). The review attributes the same bounds to DSNE Rev J, "
        "which could not be obtained.",
        "PASS means Diviner lies inside the envelope within the 10 K tolerance of TH-V3: cold "
        "bias at most +10 K and hot bias at least -10 K. It does not select a law. Caveats of "
        "TH-V3 apply: 8 sunlit bins, one 720 m block, no terrain radiation, per-cell "
        "patterns not reproduced.",
    )
    summary = (
        f"cold law {cold_bias:+.1f} K, hot law {hot_bias:+.1f} K over "
        f"{int(used.sum())} sunlit summer bins"
    )
    return Outcome(checks, summary, notes)


# ------------------------------------------------------------- terrain radiation (change 3)
_TH_A3_TOLERANCE = 0.05  # relative, pre-registered 2026-11-16
_TH_A3_FINE = {
    "n_sectors": 180,
    "band_edges_m": (*(10.0 * k for k in range(120)), float("inf")),
    "n_rays": 720,
    "patch_azimuths": 360,
}
_TH_C1_THRESHOLD_W_M2 = 0.05  # the review's falsification threshold


def _th_a3(_: DataStore) -> Outcome:
    from lunarenv.thermal import terrain_radiation as tr
    from lunarenv.validation import bowl

    fine = bowl.run_bowl(**_TH_A3_FINE)  # type: ignore[arg-type]
    if not all(fine.in_shadow):
        raise RuntimeError("a TH-A3 receiver is sunlit; the analytic reference does not apply")
    sectors, edges = tr.patch_layout()
    default = bowl.run_bowl(n_sectors=sectors, band_edges_m=edges, n_rays=360, patch_azimuths=360)
    single = bowl.run_bowl(**_TH_A3_FINE, points_per_patch=1)  # type: ignore[arg-type]
    checks: list[Check] = []
    for i, offset in enumerate(fine.offsets_m):
        where = "bowl centre" if offset == 0.0 else f"floor, {offset:+.0f} m across the Sun"
        checks += [
            Check(
                quantity=f"scattered sunlight, {where}",
                reference=fine.reference_visible_w_m2,
                measured=fine.visible_w_m2[i],
                tolerance=_TH_A3_TOLERANCE * fine.reference_visible_w_m2,
                unit="W m-2",
            ),
            Check(
                quantity=f"terrain infrared, {where}",
                reference=fine.reference_infrared_w_m2,
                measured=fine.infrared_w_m2[i],
                tolerance=_TH_A3_TOLERANCE * fine.reference_infrared_w_m2,
                unit="W m-2",
            ),
        ]
    for label, result in (("default patch layout", default), ("fine, one point per patch", single)):
        for i, offset in enumerate(result.offsets_m):
            checks.append(
                Check(
                    quantity=f"scattered sunlight, {label}, offset {offset:+.0f} m",
                    reference=result.reference_visible_w_m2,
                    measured=result.visible_w_m2[i],
                    unit="W m-2",
                    informational=True,
                )
            )
    worst = max(abs(v / fine.reference_visible_w_m2 - 1.0) for v in (*fine.visible_w_m2,))
    worst_default = max(abs(v / default.reference_visible_w_m2 - 1.0) for v in default.visible_w_m2)
    notes = (
        "Pre-registered in docs/review/preregistration-th-a3-2026-11-16.md (change 3 of the "
        "2026-11-13 scientific review), amended before recording: patches are represented by "
        "3 points stratified by solid angle, not by the one sample of largest solid angle.",
        f"Spherical bowl of rim radius {bowl.RIM_RADIUS_M:g} m and depth {bowl.DEPTH_M:g} m "
        f"(sphere radius {bowl.sphere_radius_m():g} m) in level ground, synthetic DEM of "
        f"{bowl.PIXEL_M:g} m pixels; Sun fixed at {bowl.SUN_ELEVATION_DEG:g} deg elevation, "
        f"constant albedo {bowl.ALBEDO:g}, patches in instantaneous radiative equilibrium.",
        "Reference: inside a sphere the view factor between elements is dA / (4 pi R^2), so "
        "a shadowed point receives A F sin(e) r^2 / (4 R^2) of scattered sunlight and "
        "(1 - A) F sin(e) r^2 / (4 R^2) of infrared, single bounce (Ingersoll1992).",
        f"Fine layout (checks): 180 sectors, 10 m bands, 720 rays, {fine.n_patches} patches. "
        f"Default scenario layout (informational): {sectors} sectors, bands growing by sqrt(2) "
        f"from 25 m, largest error {100 * worst_default:.1f} %.",
        "The ratio of infrared to scattered light is (1 - A) / A by construction; the two "
        "checks per receiver test the bookkeeping, not independent physics.",
    )
    summary = (
        f"scattered {', '.join(f'{v:.3f}' for v in fine.visible_w_m2)} W m-2 against "
        f"{fine.reference_visible_w_m2:.3f} (largest error {100 * worst:.1f} %); default "
        f"layout within {100 * worst_default:.1f} %"
    )
    return Outcome(checks, summary, notes)


def _reference_scenario(**thermal: object) -> dict:  # type: ignore[type-arg]
    from pathlib import Path

    import yaml

    path = Path(__file__).resolve().parents[3] / "scenarios" / "shackleton_rim_2026.yaml"
    if not path.is_file():
        raise FileNotFoundError(f"TH-C1 needs the reference scenario at {path} (source checkout)")
    scenario: dict = yaml.safe_load(path.read_text())  # type: ignore[type-arg]
    scenario["uncertainty"] = {"members": 0, "terrain_clones": 0, "seed": 1}
    scenario["thermal"].update(
        {"slope_grid": None, "constant_albedo_alternative": False, **thermal}
    )
    return scenario


def _th_c1(store: DataStore) -> Outcome:
    from lunarenv.runner import run_scenario

    run = run_scenario(_reference_scenario(), store=store)
    thermal = run.metadata["summary"]["thermal"]
    diffuse = thermal["diffuse"]
    bound = float(thermal["neglected_terrain_flux_upper_bound_w_m2"])
    checks: list[Check] = []
    labels = {
        "cold": "cold law (Bürger 2024)",
        "nominal": "scenario law (Hayne 2017)",
        "hot": "hot law (Foote 2020)",
    }
    for case, label in labels.items():
        flux = float(diffuse["cases"][case]["terrain_irradiance_at_minimum_w_m2"])
        checks += [
            Check(
                quantity=f"terrain irradiance at the minimum, {label}: within the earlier bound",
                reference=0.0,
                measured=flux,
                tolerance=bound,
                unit="W m-2",
            ),
            Check(
                quantity=f"terrain irradiance at the minimum, {label}: below the review threshold",
                reference=0.0,
                measured=flux,
                tolerance=_TH_C1_THRESHOLD_W_M2,
                unit="W m-2",
            ),
        ]
    for case, label in labels.items():
        entry = diffuse["cases"][case]
        checks += [
            Check(
                quantity=f"Earthshine absorbed at the minimum, {label}",
                measured=float(entry["earthshine_absorbed_at_minimum_w_m2"]),
                unit="W m-2",
                informational=True,
            ),
            Check(
                quantity=f"change of the minimum from terrain radiation, {label}",
                measured=float(entry["terrain_effect_k"]["surface_min"]),
                unit="K",
                informational=True,
            ),
            Check(
                quantity=f"change of the minimum from Earthshine, {label}",
                measured=float(entry["earthshine_effect_k"]["surface_min"]),
                unit="K",
                informational=True,
            ),
        ]
    # patch convergence on the ground receiver alone (asset surfaces off, so that the
    # doubled layout stays affordable)
    ground = {}
    for refinement in (1, 2):
        r = run_scenario(
            _reference_scenario(oriented_surfaces=[], terrain_patch_refinement=refinement),
            store=store,
        )
        ground[refinement] = float(
            r.metadata["summary"]["thermal"]["diffuse"]["cases"]["nominal"][
                "terrain_irradiance_at_minimum_w_m2"
            ]
        )
    checks.append(
        Check(
            quantity="terrain irradiance at the minimum, scenario law, twice the sectors and bands",
            reference=ground[1],
            measured=ground[2],
            unit="W m-2",
            informational=True,
        )
    )
    nominal = diffuse["cases"]["nominal"]
    notes = (
        "Pre-registered in docs/review/preregistration-th-a3-2026-11-16.md (change 3 of the "
        "2026-11-13 scientific review, criteria (b) and its falsification threshold).",
        "Run of scenarios/shackleton_rim_2026.yaml with the uncertainty ensemble, the slope "
        "grid and the constant-albedo alternative switched off (none of them changes the "
        "nominal surface element). Ground receiver at the site with the DEM normal.",
        f"The earlier bound ({bound:.4f} W m-2) assumed level ground and full sunlight on all "
        "terrain above the horizon; the receiver here has the DEM slope, whose plane hides "
        f"almost all terrain (terrain view weight {diffuse['terrain_view_weight']:.2e}).",
        f"Scenario law: minimum {nominal['surface_min_k']:.1f} K at {nominal['time_of_minimum_utc']}; "
        f"absorbed diffuse flux there {diffuse['absorbed_diffuse_at_minimum_w_m2']:.4f} W m-2 "
        f"({100 * diffuse['absorbed_diffuse_at_minimum_over_emission']:.1f} % of the emission).",
        "Convergence check on the ground receiver alone (no asset surfaces), default and "
        "doubled patch layout.",
    )
    fluxes = [float(diffuse["cases"][c]["terrain_irradiance_at_minimum_w_m2"]) for c in labels]
    summary = (
        f"terrain irradiance at the minimum {max(fluxes):.2e} W m-2 at most (bound {bound:.4f}, "
        f"threshold {_TH_C1_THRESHOLD_W_M2}); Earthshine absorbed "
        f"{nominal['earthshine_absorbed_at_minimum_w_m2']:.4f} W m-2"
    )
    return Outcome(checks, summary, notes)


# ------------------------------------------------------------ TH-V2c (like for like)
_A17_LATITUDE_DEG = 20.0
_A17_HEAT_FLOW = 0.016  # W m-2, Langseth final report Table 1 (Taurus Littrow)
_A15_HEAT_FLOW = 0.021  # W m-2, same table (Hadley Rille)
_SITE_ALBEDO = 0.08  # Keihm & Langseth 1973, reduction equation


def _apollo_like(
    latitude_deg: float, properties: RegolithProperties, depth_m: float, distance_au: float = 1.0
) -> tuple[float, float]:
    """Diurnal mean brightness temperature (e_m = 1) at the surface and mean temperature
    at depth."""
    result = _cycle(latitude_deg, properties, distance_au=distance_au)
    brightness = properties.emissivity**0.25 * result.surface_temperature_k[:-1]
    return float(brightness.mean()), float(result.at_depth(depth_m)[:-1].mean())


def _th_v2c(_: DataStore) -> Outcome:
    a17 = RegolithProperties(albedo_normal=_SITE_ALBEDO, heat_flow_w_m2=_A17_HEAT_FLOW)
    surface, deep = _apollo_like(_A17_LATITUDE_DEG, a17, 1.30)
    checks = [
        Check(
            quantity="Apollo 17 (20 N): diurnal mean surface brightness temperature (e_m = 1)",
            reference=216.0,
            measured=surface,
            tolerance=5.0,
            unit="K",
        ),
        Check(
            quantity="Apollo 17 (20 N): diurnal mean temperature at 1.30 m",
            reference=256.0,
            measured=deep,
            tolerance=5.0,
            unit="K",
        ),
    ]
    constant = a17.model_copy(update={"albedo_a": 0.0, "albedo_b": 0.0})
    c_surface, c_deep = _apollo_like(_A17_LATITUDE_DEG, constant, 1.30)
    f_surface, f_deep = _apollo_like(_A17_LATITUDE_DEG, a17, 1.30, distance_au=1.017)
    a15 = RegolithProperties(albedo_normal=_SITE_ALBEDO, heat_flow_w_m2=_A15_HEAT_FLOW)
    s15, d15_hayne = _apollo_like(26.0, a15, 0.83)
    d15_report = _apollo_like(26.0, a15, 1.38)[1]
    info = [
        ("Apollo 17, constant albedo 0.08: surface brightness mean", 216.0, c_surface),
        ("Apollo 17, constant albedo 0.08: mean at 1.30 m", 256.0, c_deep),
        ("Apollo 17, Sun at 1.017 au: surface brightness mean", 216.0, f_surface),
        ("Apollo 17, Sun at 1.017 au: mean at 1.30 m", 256.0, f_deep),
        ("Apollo 15 (26 N): surface brightness mean, against Hayne2017's 211 K", 211.0, s15),
        ("Apollo 15: mean at 0.83 m, against Hayne2017's 252 K", 252.0, d15_hayne),
        ("Apollo 15: surface brightness mean, against the final report's 207 K", 207.0, s15),
        ("Apollo 15: mean at 1.38 m, against the final report's 253.0 K", 253.0, d15_report),
    ]
    checks += [
        Check(quantity=q, reference=r, measured=m, unit="K", informational=True) for q, r, m in info
    ]
    notes = (
        "Pre-registered in docs/review/preregistration-th-v2c-2026-11-16.md (change 3 (d) and "
        "finding F05 of the 2026-11-13 scientific review) before it was implemented or run.",
        "Reference: Keihm & Langseth (1973), Proc. Lunar Sci. Conf. 4th, 2503-2513, read from "
        "the ADS scan: thermocouples about 15 cm above the surface, reduced with surface "
        "emissivity 1.0 (blackbody brightness temperatures), mean solar constant 0.1352 W cm-2 "
        "and lunar albedo 0.08; mean surface temperature 216 +/- 5 K and 256 K at 130 cm.",
        "Model: Hayne2017 parameters with the site's albedo A0 = 0.08 (incidence dependence "
        "kept) and the measured Taurus Littrow heat flow 0.016 W m-2 (Langseth et al. final "
        "report, NTRS 19780010035, Table 1); Sun at 1 au, the lunation of the source's data "
        "being unstated. The surface value compared is the mean of eps^(1/4) T_s, eps = 0.95.",
        "Apollo 15 is informational: the final report gives 207 +/- 7 K (not the 211 K quoted "
        "by Hayne2017, whose source, Keihm et al. 1973 EPSL, was not read), computed from "
        "night data supplemented by daytime models, and lowered mainly by early-morning "
        "shadowing from the local topography, which a 1-D level model does not have.",
    )
    summary = (
        f"Apollo 17 brightness mean {surface:.1f} K against 216 K, {deep:.1f} K at 1.30 m "
        f"against 256 K (+/-5 K) [site albedo 0.08, e_m = 1]"
    )
    return Outcome(checks, summary, notes)


# ------------------------------------------------- F09: low-temperature conductivity
_LOW_T = "woods_robinson2019"


def _th_a2b(_: DataStore) -> Outcome:
    heat_flow, k = 0.005, 1.0e-5
    properties = RegolithProperties(
        surface_conductivity_w_m_k=k,
        deep_conductivity_w_m_k=k,
        surface_density_kg_m3=1500.0,
        deep_density_kg_m3=1500.0,
        radiative_parameter=2.7,
        heat_flow_w_m2=heat_flow,
        conductivity_model=_LOW_T,
    )
    grid = geometric_grid(0.004, 1.3, 0.4)
    result = solve_conduction(
        grid,
        properties,
        np.full(grid.n_cells, 200.0),
        np.full(400, 0.0),
        1.0e9,
        scheme="backward_euler",
    )
    surface = (heat_flow / (properties.emissivity * SIGMA)) ** 0.25
    potential = kirchhoff(result.temperature_k[-1], k, 2.7, 1.0)
    expected = kirchhoff(surface, k, 2.7, 1.0) + heat_flow * grid.centers_m
    error = float(np.max(np.abs(potential - expected) / expected))
    checks = [
        Check(
            quantity="surface temperature",
            reference=float(surface),
            measured=float(result.surface_temperature_k[-1]),
            tolerance=1e-6,
            unit="K",
        ),
        Check(
            quantity="relative error of the Kirchhoff potential profile",
            reference=0.0,
            measured=error,
            tolerance=1e-9,
        ),
    ]
    deepest = float(result.temperature_k[-1, -1])
    notes = (
        "Pre-registered in docs/review/preregistration-f09-2026-11-16.md before the option was "
        "implemented or run.",
        "conductivity_model = woods_robinson2019: K = K_c [g(T) + chi (T/350)^3], g = 1 above "
        "150 K, k_am(T)/k_am(150 K) between 20 and 150 K (Woods-Robinson et al. 2019, sec. 4.3, "
        "G_am with the sign inferred from their Figure 6c) and linear to zero below 20 K. The "
        "profile runs from the surface through all three branches.",
    )
    return Outcome(
        checks,
        f"profile error {error:.1e} (relative) from {surface:.1f} K at the surface to "
        f"{deepest:.0f} K at the base [low-T conductivity]",
        notes,
    )


def _th_c2(store: DataStore) -> Outcome:
    from lunarenv.runner import run_scenario

    run = run_scenario(_reference_scenario(terrain_radiation=False, earthshine=False), store=store)
    alt = run.metadata["model_form_sensitivity"]["conductivity"]
    nominal = run.metadata["summary"]["thermal"]
    info: list[tuple[str, float, float]] = [
        ("Rim, surface minimum", nominal["surface_min_k"], alt["surface_min_k"]),
        ("Rim, surface mean", nominal["surface_mean_k"], alt["surface_mean_k"]),
    ]
    for depth in ("0.3 m", "1 m"):
        info.append(
            (
                f"Rim, mean at {depth}",
                float(nominal["at_depth"][depth]["mean_k"]),
                float(alt["depth_mean_k"][depth]),
            )
        )
    low_t = RegolithProperties(conductivity_model=_LOW_T)
    _, equator = _equator_checks(low_t)
    info += [
        ("Equator (TH-V1 setting), noon peak, against 385 K", 385.0, equator[0]),
        ("Equator, midnight, against 101 K", 101.0, equator[1]),
        ("Equator, night minimum, against 95 K", 95.0, equator[2]),
    ]
    a17 = RegolithProperties(
        albedo_normal=_SITE_ALBEDO, heat_flow_w_m2=_A17_HEAT_FLOW, conductivity_model=_LOW_T
    )
    surface, deep = _apollo_like(_A17_LATITUDE_DEG, a17, 1.30)
    info += [
        ("Apollo 17 (TH-V2c setting), surface brightness mean, against 216 K", 216.0, surface),
        ("Apollo 17, mean at 1.30 m, against 256 K", 256.0, deep),
    ]
    checks = [
        Check(quantity=q, reference=r, measured=m, unit="K", informational=True) for q, r, m in info
    ]
    change = alt["difference_depth_mean_k"]
    notes = (
        "Pre-registered in docs/review/preregistration-f09-2026-11-16.md (finding F09 of the "
        "2026-11-13 scientific review) before the option was implemented or run. Sensitivity "
        "case: no pass/fail.",
        "Rim rows: reference = the run with the Hayne2017 conductivity, measured = the same run "
        "with the Woods-Robinson 2019 low-temperature contact term (scenarios/"
        "shackleton_rim_2026.yaml, nominal albedo law, no ensemble, no slope grid, terrain "
        "radiation and Earthshine off: TH-C1 shows they move the rim minimum by less than "
        "0.3 K).",
        "Equator and Apollo rows: the TH-V1 and TH-V2c settings with the option; TH-V1's 2 K and "
        "TH-V2c's 5 K are not applied here.",
        "The option keeps Hayne2017's K_s, K_d and chi, which were fitted with a "
        "temperature-independent contact term; it is a bound on the conductivity form, not a "
        "refitted model. No measurement of lunar regolith conductivity below 100 K is known "
        "to this project.",
    )
    summary = (
        f"rim: minimum {alt['difference_surface_min_k']:+.1f} K, mean "
        f"{alt['difference_surface_mean_k']:+.1f} K, 0.3 m {change['0.3 m']:+.1f} K, 1 m "
        f"{change['1 m']:+.1f} K; equator night minimum {equator[2]:.1f} K (95 K); Apollo 17 "
        f"{surface:.1f} / {deep:.1f} K (216 / 256 K)"
    )
    return Outcome(checks, summary, notes)


# --------------------------------------------- TH-V3e: four more sites (2026-11-17)
_TH_V3E_MIN_SUNLIT = 4
_TH_V3E_LAWS = ("foote2020", "burger2024", "hayne", "constant")
_TH_V3E_NAMES = {
    "foote2020": "hot bound (Foote 2020)",
    "burger2024": "cold bound (Burger 2024)",
    "hayne": "Hayne2017",
    "constant": "constant 0.12",
}


def _th_v3e(store: DataStore) -> Outcome:
    import json
    from importlib import resources

    from lunarenv.validation import diviner_pcp as dp

    path = resources.files("lunarenv.validation") / "reference" / "diviner_pcp_sites.json"
    reference = json.loads(path.read_text(encoding="utf-8"))
    checks: list[Check] = []
    parts: list[str] = []
    qualifying: list[bool] = []
    passed: list[bool] = []
    cold_side: list[bool] = []
    for label, near_id in dp.TH_V3E_SITES:
        selection = reference["selection"][label]
        data = reference["sites"][label]
        keys = [c["key"] for c in data["cells"]]
        centres = [(float(c["x_m"]), float(c["y_m"])) for c in data["cells"]]
        diviner = np.full((len(keys), dp.N_BINS), np.nan)
        for b, entry in enumerate(data["bins"]):
            for c, k in enumerate(keys):
                if k in entry["tbol_k"]:
                    diviner[c, b] = entry["tbol_k"][k]
        cells = dp.forward(store, centres, 20.0, near_id=near_id)
        visible = np.mean([c.sums["visible"] / c.sums["count"] for c in cells], axis=0)
        blocks = {}
        for law in _TH_V3E_LAWS:
            radiance = np.array([c.sums[law] / c.sums["count"] for c in cells])
            blocks[law] = _block_from_radiance(radiance, diviner, visible)
        _, _, used, complete = blocks["hayne"]
        n_used = int(used.sum())
        ok = n_used >= _TH_V3E_MIN_SUNLIT
        qualifying.append(ok)
        bias = {
            law: float(np.mean(b[1][used])) if n_used else float("nan") for law, b in blocks.items()
        }
        hot = bias["foote2020"]
        if ok:
            passed.append(abs(hot) <= _TH_V3_TOLERANCE_K)
            cold_side.append(bias["burger2024"] < -_TH_V3_TOLERANCE_K)
        where = (
            f"{label} ({selection['latitude_deg']:.2f} deg, block centre "
            f"{selection['centre_m'][0]:.0f}, {selection['centre_m'][1]:.0f} m)"
        )
        checks.append(
            Check(
                quantity=f"{where}: hot bound, bias of block Tbol over sunlit summer bins"
                + ("" if ok else f" (fewer than {_TH_V3E_MIN_SUNLIT} sunlit bins: not in verdict)"),
                reference=0.0,
                measured=hot if n_used else None,
                tolerance=_TH_V3_TOLERANCE_K if ok else None,
                unit="K",
                informational=not ok,
            )
        )
        for law in _TH_V3E_LAWS[1:]:
            checks.append(
                Check(
                    quantity=f"{label}: {_TH_V3E_NAMES[law]}, bias over sunlit bins",
                    reference=0.0,
                    measured=bias[law] if n_used else None,
                    unit="K",
                    informational=True,
                )
            )
        diviner_block = blocks["hayne"][0] - blocks["hayne"][1]
        between = (blocks["burger2024"][0] <= diviner_block) & (
            diviner_block <= blocks["foote2020"][0]
        )
        checks += [
            Check(
                quantity=f"{label}: RMS over sunlit bins, hot bound / Hayne2017",
                reference=float(np.sqrt(np.mean(blocks["foote2020"][1][used] ** 2)))
                if n_used
                else None,
                measured=float(np.sqrt(np.mean(blocks["hayne"][1][used] ** 2))) if n_used else None,
                unit="K",
                informational=True,
            ),
            Check(
                quantity=f"{label}: bias over all summer bins with data, hot bound / Hayne2017",
                reference=float(np.mean(blocks["foote2020"][1][complete]))
                if complete.any()
                else None,
                measured=float(np.mean(blocks["hayne"][1][complete])) if complete.any() else None,
                unit="K",
                informational=True,
            ),
            Check(
                quantity=f"{label}: sunlit bins used / bins with all nine cells",
                reference=float(complete.sum()),
                measured=float(n_used),
                informational=True,
            ),
            Check(
                quantity=f"{label}: sunlit bins with Diviner between the cold and hot bounds",
                reference=float(n_used),
                measured=float(between[used].sum()),
                informational=True,
            ),
            Check(
                quantity=f"{label}: block mean visible fraction, selection (cell centres) / "
                "summer (all pixels)",
                reference=float(selection["block_mean_visible_fraction"]),
                measured=float(np.mean(visible)),
                informational=True,
            ),
        ]
        parts.append(
            f"{label} hot {hot:+.1f} / cold {bias['burger2024']:+.1f} / Hayne "
            f"{bias['hayne']:+.1f} K ({n_used} bins)"
        )
    n_ok = sum(qualifying)
    if n_ok < 3:
        status: Status | None = Status.PARTIAL
    else:
        status = None  # PASS if every toleranced (qualifying) check passes
    reading = (
        "hot bound within 10 K and cold bound more than 10 K too cold at every qualifying site"
        if passed and all(passed) and all(cold_side)
        else "the Shackleton pattern does not hold at every qualifying site"
    )
    notes = (
        "Pre-registered in docs/review/preregistration-th-v3e-2026-11-17.md before it was "
        "implemented, before the site DEMs were downloaded and before any Diviner value "
        "outside the TH-V3 block was read.",
        "Blocks chosen from model geometry only (highest summer mean visible disc at the nine "
        "cell centres, 1 km inside the tile, TH-V3 block excluded) and recorded in "
        "diviner_pcp_sites.json before the Diviner values were extracted.",
        "Forward model, bins and data as TH-V3 (no terrain radiation or Earthshine: under 1 K "
        "at Shackleton, TH-V3d). A site needs at least 4 sunlit bins with all nine cells to "
        f"enter the verdict; {n_ok} of {len(qualifying)} sites qualify.",
        f"Reading: {reading}.",
    )
    return Outcome(checks, "; ".join(parts), notes, status)


def _th_v3f_equator() -> tuple[float, float]:
    """Level ground at the equator, equinox, Hayne2017: rough (20 deg) minus smooth Tbol
    at noon and at midnight (K)."""
    from lunarenv.core.constants import SOLAR_CONSTANT_W_M2
    from lunarenv.thermal.roughness import rough_radiance, sun_vectors

    n = 709  # hourly samples per synodic month
    hour_angle = 2.0 * np.pi * np.arange(2 * n + 1) / n
    elevation = np.arcsin(np.clip(np.cos(hour_angle), -1.0, 1.0))
    azimuth = np.where(np.sin(hour_angle) < 0.0, np.pi / 2.0, -np.pi / 2.0)  # east a.m.
    irradiance = np.where(elevation > 0.0, SOLAR_CONSTANT_W_M2, 0.0)
    result = rough_radiance(
        irradiance,
        sun_vectors(elevation, azimuth),
        [0.0, 0.0, 1.0],
        RegolithProperties(),
        (20.0,),
        SYNODIC_MONTH_S / n,
        n,
        ThermalSettings(),
    )
    smooth = result.smooth**0.25
    rough = result.rough[20.0] ** 0.25
    return float(rough[0] - smooth[0]), float(rough[n // 2] - smooth[n // 2])


_TH_V3F_SLOPES = (0.0, 20.0, 35.0)


def _th_v3f(store: DataStore) -> Outcome:
    from lunarenv.illumination.settings import SITE04_NEAR
    from lunarenv.validation import diviner_pcp as dp

    path = resources.files("lunarenv.validation") / "reference" / "diviner_pcp_sites.json"
    sites_ref = json.loads(path.read_text(encoding="utf-8"))
    shackleton = dp.load_reference()
    sites: list[tuple[str, str, dict]] = [("Shackleton", SITE04_NEAR, shackleton)]  # type: ignore[type-arg]
    sites += [(label, near, sites_ref["sites"][label]) for label, near in dp.TH_V3E_SITES]
    checks: list[Check] = []
    parts: list[str] = []
    qualifying: list[bool] = []
    passed: list[bool] = []
    improvement: list[float] = []
    hot_too_warm: list[str] = []
    fine_check: Check | None = None
    for s_i, (label, near_id, data) in enumerate(sites):
        keys = [c["key"] for c in data["cells"]]
        centres = [(float(c["x_m"]), float(c["y_m"])) for c in data["cells"]]
        diviner = np.full((len(keys), dp.N_BINS), np.nan)
        for b, entry in enumerate(data["bins"]):
            for c, k in enumerate(keys):
                if k in entry["tbol_k"]:
                    diviner[c, b] = entry["tbol_k"][k]
        cells = dp.forward(store, centres, 20.0, near_id=near_id, laws=dp.ROUGH_LAWS)
        middle = len(centres) // 2
        px, py = dp.cell_pixels(*centres[middle], 120.0)
        point = (float(px[0]), float(py[0]))
        receivers, fine = dp.forward_rough(
            store, centres, near_id=near_id, convergence=point if s_i == 0 else None
        )
        visible = np.mean([c.sums["visible"] / c.sums["count"] for c in cells], axis=0)

        def change(key: str, law: str, rs_all: list = receivers) -> np.ndarray:  # type: ignore[type-arg]
            return np.array(
                [np.mean([(r[key] - r[law]) / r["count"] for r in rs], axis=0) for rs in rs_all]
            )

        blocks = {}
        for law in dp.ROUGH_LAWS:
            base = np.array([c.sums[law] / c.sums["count"] for c in cells])
            for rms in _TH_V3F_SLOPES:
                rad = base if rms == 0.0 else base + change(dp.rough_key(law, rms), law)
                blocks[(law, rms)] = _block_from_radiance(rad, diviner, visible)
        _, _, used, complete = blocks[("hayne", 0.0)]
        n_used = int(used.sum())
        ok = n_used >= _TH_V3E_MIN_SUNLIT
        qualifying.append(ok)

        def bias(
            law: str,
            rms: float,
            mask: np.ndarray = used,
            blocks: dict = blocks,  # type: ignore[type-arg]
        ) -> float:
            return float(np.mean(blocks[(law, rms)][1][mask])) if mask.any() else float("nan")

        main = bias("hayne", 20.0)
        if ok:
            passed.append(abs(main) <= _TH_V3_TOLERANCE_K)
            improvement.append(main - bias("hayne", 0.0))
            if abs(bias("foote2020", 0.0)) <= _TH_V3_TOLERANCE_K and (
                bias("foote2020", 20.0) > _TH_V3_TOLERANCE_K
            ):
                hot_too_warm.append(label)
        checks.append(
            Check(
                quantity=f"{label}: Hayne2017 with 20 deg RMS roughness, bias of block Tbol "
                "over sunlit summer bins"
                + ("" if ok else f" (fewer than {_TH_V3E_MIN_SUNLIT} sunlit bins: not in verdict)"),
                reference=0.0,
                measured=main if n_used else None,
                tolerance=_TH_V3_TOLERANCE_K if ok else None,
                unit="K",
                informational=not ok,
            )
        )
        for law, rms, name in (
            ("hayne", 35.0, "Hayne2017, 35 deg"),
            ("hayne", 0.0, "Hayne2017, smooth (as TH-V3 / TH-V3e)"),
            ("foote2020", 0.0, "hot law, smooth (as TH-V3c / TH-V3e)"),
            ("foote2020", 20.0, "hot law, 20 deg"),
            ("foote2020", 35.0, "hot law, 35 deg"),
        ):
            checks.append(
                Check(
                    quantity=f"{label}: {name}, bias over sunlit bins",
                    reference=0.0,
                    measured=bias(law, rms) if n_used else None,
                    unit="K",
                    informational=True,
                )
            )
        hayne20 = blocks[("hayne", 20.0)]
        checks += [
            Check(
                quantity=f"{label}: Hayne2017 20 deg, RMS over sunlit bins / bias over all "
                "summer bins with data",
                reference=float(np.sqrt(np.mean(hayne20[1][used] ** 2))) if n_used else None,
                measured=bias("hayne", 20.0, complete),
                unit="K",
                informational=True,
            ),
            Check(
                quantity=f"{label}: change of block Tbol from roughness (Hayne2017, sunlit "
                "bins), 20 / 35 deg",
                reference=float(np.mean(hayne20[0][used] - blocks[("hayne", 0.0)][0][used]))
                if n_used
                else None,
                measured=float(
                    np.mean(blocks[("hayne", 35.0)][0][used] - blocks[("hayne", 0.0)][0][used])
                )
                if n_used
                else None,
                unit="K",
                informational=True,
            ),
            Check(
                quantity=f"{label}: sunlit bins used / bins with all nine cells",
                reference=float(complete.sum()),
                measured=float(n_used),
                informational=True,
            ),
        ]
        if fine is not None:
            coarse_r = next(r for r in receivers[middle] if np.allclose(r["xy"], point))
            key = dp.rough_key("hayne", 20.0)
            sunlit_c = used & (coarse_r["count"] > 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                t_coarse = (coarse_r[key] / coarse_r["count"]) ** 0.25
                t_fine = (fine[key] / fine["count"]) ** 0.25
            fine_check = Check(
                quantity="quadrature: Shackleton centre-cell receiver, Hayne2017 20 deg, "
                "sunlit-bin mean Tbol, 3 nodes x 12 azimuths / 5 x 24",
                reference=float(np.mean(t_coarse[sunlit_c])) if sunlit_c.any() else None,
                measured=float(np.mean(t_fine[sunlit_c])) if sunlit_c.any() else None,
                unit="K",
                informational=True,
            )
        parts.append(
            f"{label} Hayne 20 deg {main:+.1f} K (smooth {bias('hayne', 0.0):+.1f}), hot 20 deg "
            f"{bias('foote2020', 20.0):+.1f} K ({n_used} bins)"
        )
    noon, midnight = _th_v3f_equator()
    checks.append(
        Check(
            quantity="equator, equinox, level, Hayne2017: Tbol change from 20 deg roughness, "
            "noon / midnight",
            reference=noon,
            measured=midnight,
            unit="K",
            informational=True,
        )
    )
    if fine_check is not None:
        checks.append(fine_check)
    n_ok = sum(qualifying)
    status: Status | None = Status.PARTIAL if n_ok < 3 else None
    if passed and all(passed):
        reading = (
            "roughness of the size measured by Diviner at the equator explains the cold bias of "
            "the default law at every qualifying block"
        )
    elif improvement and all(d > 10.0 for d in improvement):
        reading = "roughness is a major part of the gap (over 10 K at every site), not all of it"
    elif improvement and all(abs(d) < 5.0 for d in improvement):
        reading = "roughness of this kind does not explain the gap (under 5 K at every site)"
    else:
        reading = "roughness closes part of the gap at some sites; see the per-site values"
    extra = []
    if hot_too_warm:
        extra.append(
            "hot law with 20 deg roughness more than 10 K too warm at "
            + ", ".join(hot_too_warm)
            + ", where it passed smooth"
        )
    if abs(noon) > 5.0:
        extra.append(
            f"equatorial noon changes by {noon:+.1f} K: the default law and the roughness "
            "model are not independent"
        )
    notes = (
        "Pre-registered in docs/review/preregistration-th-v3f-2026-11-17.md before the "
        "roughness model was implemented or run.",
        "Roughness (lunarenv.thermal.roughness, after Bandfield et al. 2015, section 2): "
        "Gaussian slopes, tan(theta_rms)^2 = <|grad h|^2>, 3 Gauss-Laguerre nodes x 12 "
        "azimuths about the DEM normal; Smith (1967) cast shadows, lit and shadowed "
        "populations; adjacent-surface sunlight and infrared over (1 - cos beta)/2 of the "
        "hemisphere from a plane with the smooth facet's sunlight and temperature (radiative "
        "equilibrium over the spin-up year); nadir view along the DEM normal; each facet runs "
        "the 1-D model.",
        "Change (rough - smooth radiance) at 4 receivers per cell, averaged per cell and "
        "added to the smooth forward model at every 20 m (as TH-V3d). Sunlit bins from the "
        "smooth model geometry. Not included: lateral conduction, shadow-orientation "
        "correlation, multiple scattering, oblique views.",
        f"{n_ok} of {len(qualifying)} sites have at least {_TH_V3E_MIN_SUNLIT} sunlit bins.",
        f"Reading: {reading}." + ("" if not extra else " Also: " + "; ".join(extra) + "."),
    )
    return Outcome(checks, "; ".join(parts), notes, status)


CASES: tuple[Case, ...] = (
    Case(
        "TH-A1",
        "Damped temperature wave under a periodic surface temperature",
        "analytic",
        "CarslawJaeger1959",
        "max |T - exact|; observed spatial order",
        "0.25 K (0.25 % of amplitude); order 2 +/- 0.3",
        run=_th_a1,
    ),
    Case(
        "TH-A2",
        "Steady conduction with interior heat flow, K depending on T",
        "analytic",
        "closed form (Kirchhoff transformation)",
        "relative error of the potential profile",
        "1e-9",
        run=_th_a2,
    ),
    Case(
        "TH-A2b",
        "Steady conduction with the low-temperature conductivity option",
        "analytic",
        "closed form (Kirchhoff transformation)",
        "relative error of the potential profile, 17 K surface to above 150 K at the base",
        "1e-9 (pre-registered 2026-11-16)",
        run=_th_a2b,
    ),
    Case(
        "TH-P1",
        "Energy conservation of the four time integrators",
        "property",
        "energy identity of each scheme",
        "relative residual of stored energy vs boundary flux",
        "1e-9",
        run=_th_p1,
    ),
    Case(
        "TH-V1",
        "Equatorial diurnal temperatures",
        "published",
        "Hayne2017 Table A2 (Vasavada2012, Diviner)",
        "noon, midnight and minimum surface temperature",
        "2 K (twice the < 1 K RMS residual Hayne2017 reports for its fits)",
        run=_th_v1,
    ),
    Case(
        "TH-V1b",
        "Equatorial diurnal temperatures, constant albedo (discrimination check)",
        "published",
        "Hayne2017 Table A2 (Vasavada2012, Diviner)",
        "noon, midnight and minimum surface temperature",
        "2 K (as TH-V1; pre-registered 2026-11-13)",
        run=_th_v1b,
    ),
    Case(
        "TH-V2",
        "Apollo 15 and 17 mean temperatures, published parameter set",
        "published",
        "Hayne2017 Table A2 (Keihm1973)",
        "diurnal mean at the surface and at depth",
        "5 K (stated by the source)",
        run=_th_v2,
    ),
    Case(
        "TH-V2b",
        "Apollo 15 and 17 mean temperatures, constant albedo",
        "published",
        "Hayne2017 Table A2 (Keihm1973)",
        "diurnal mean at the surface and at depth",
        "5 K (stated by the source)",
        run=_th_v2b,
    ),
    Case(
        "TH-X1",
        "Code-to-code comparison with the reference implementation",
        "cross-check",
        "heat1d 0.3.2",
        "surface mean, max, min and deep mean at 0, 20, 26 deg latitude",
        "1 K (convergence level stated in Hayne2017)",
        run=_th_x1,
    ),
    Case(
        "TH-X1b",
        "Code-to-code comparison with the current reference implementation",
        "cross-check",
        "heat1d 0.4.1 (GitHub 5b73a6c)",
        "surface mean, max, min and deep mean at 0, 20, 26 deg latitude",
        "1 K (as TH-X1; pre-registered 2026-11-15)",
        run=_th_x1b,
    ),
    Case(
        "TH-X1c",
        "Code-to-code comparison with the reference implementation on a converged grid",
        "cross-check",
        "heat1d 0.4.1 (GitHub 5b73a6c), grid m=n=40",
        "surface mean, max, min and deep mean at 0, 20, 26 deg latitude",
        "1 K (as TH-X1b; pre-registered 2026-11-16)",
        run=_th_x1c,
    ),
    Case(
        "TH-V3",
        "Polar summer temperatures on the Shackleton rim against Diviner",
        "published",
        "Diviner PCP (diurnal, summer), PDS urn:nasa:pds:lro_diviner_derived1",
        "bias of 720 m block bolometric temperature over sunlit summer bins",
        "10 K (pre-registered 2026-11-15)",
        requires=("pgda78.site04.surf", "pgda90.ldem_80s_80mpp"),
        slow=True,
        run=_th_v3,
    ),
    Case(
        "TH-V3b",
        "Polar summer temperatures against Diviner, constant albedo",
        "published",
        "Diviner PCP (diurnal, summer), PDS urn:nasa:pds:lro_diviner_derived1",
        "bias of 720 m block bolometric temperature over sunlit summer bins",
        "10 K (pre-registered 2026-11-15)",
        requires=("pgda78.site04.surf", "pgda90.ldem_80s_80mpp"),
        slow=True,
        run=_th_v3b,
    ),
    Case(
        "TH-V3c",
        "Polar summer temperatures against Diviner, inside the albedo envelope",
        "published",
        "Diviner PCP (diurnal, summer), PDS urn:nasa:pds:lro_diviner_derived1",
        "bias of 720 m block bolometric temperature over sunlit summer bins, cold and hot law",
        "cold <= +10 K and hot >= -10 K (pre-registered 2026-11-16)",
        requires=("pgda78.site04.surf", "pgda90.ldem_80s_80mpp"),
        slow=True,
        run=_th_v3c,
    ),
    Case(
        "TH-V3d",
        "Polar summer temperatures with terrain radiation and Earthshine",
        "published",
        "Diviner PCP (diurnal, summer), PDS urn:nasa:pds:lro_diviner_derived1 (as TH-V3)",
        "bias of block Tbol over sunlit summer bins, Hayne2017 law with the diffuse terms",
        "10 K (pre-registered 2026-11-16)",
        requires=("pgda78.site04.surf", "pgda90.ldem_80s_80mpp"),
        slow=True,
        run=_th_v3d,
    ),
    Case(
        "TH-V3e",
        "Polar summer temperatures at four more sites",
        "published",
        "Diviner PCP (diurnal, summer), PDS urn:nasa:pds:lro_diviner_derived1 (as TH-V3); "
        "PGDA 78 Sites 1, 7, 11, 23",
        "bias of block Tbol over sunlit summer bins, hot bound of the albedo envelope, per site",
        "10 K at every site with at least 4 sunlit bins (pre-registered 2026-11-17)",
        requires=(
            "pgda78.site01.surf",
            "pgda78.site07.surf",
            "pgda78.site11.surf",
            "pgda78.site23.surf",
            "pgda90.ldem_80s_80mpp",
        ),
        slow=True,
        run=_th_v3e,
    ),
    Case(
        "TH-V3f",
        "Sub-DEM roughness against Diviner at five polar sites",
        "published",
        "Diviner PCP (diurnal, summer), as TH-V3 and TH-V3e; roughness after Bandfield et al. "
        "2015 (Icarus 248, 357), Smith 1967 shadowing",
        "bias of block Tbol over sunlit summer bins, Hayne2017 law with 20 deg RMS roughness, "
        "per site",
        "10 K at each of the five sites (pre-registered 2026-11-17)",
        requires=(
            "pgda78.site04.surf",
            "pgda78.site01.surf",
            "pgda78.site07.surf",
            "pgda78.site11.surf",
            "pgda78.site23.surf",
            "pgda90.ldem_80s_80mpp",
        ),
        slow=True,
        run=_th_v3f,
    ),
    Case(
        "TH-V2c",
        "Apollo 17 mean temperatures, like for like",
        "published",
        "Keihm & Langseth 1973 (Proc. LSC 4th); Langseth et al. final report (NTRS 19780010035)",
        "diurnal mean surface brightness temperature (e_m = 1) and mean at 1.30 m",
        "5 K (stated by the source; pre-registered 2026-11-16)",
        run=_th_v2c,
    ),
    Case(
        "TH-A3",
        "Shadowed floor of a spherical bowl crater, single bounce",
        "analytic",
        "Ingersoll1992 (sphere view factor, single-scattering term)",
        "scattered sunlight and terrain infrared at three shadowed floor points",
        "5 % (pre-registered 2026-11-16)",
        slow=True,
        run=_th_a3,
    ),
    Case(
        "TH-C1",
        "Terrain irradiance at the rim minimum against the earlier bound",
        "cross-check",
        "neglected-flux bound of the same run; review threshold 0.05 W m-2",
        "single-bounce terrain irradiance at the time of the minimum, per albedo case",
        "at most the bound and at most 0.05 W m-2 (pre-registered 2026-11-16)",
        requires=(*_REFERENCE_REQUIRES,),
        slow=True,
        run=_th_c1,
    ),
    Case(
        "TH-C2",
        "Effect of a low-temperature conductivity (Woods-Robinson 2019 shape)",
        "cross-check",
        "Woods-Robinson et al. 2019 (sec. 4.3); review finding F09",
        "change of the rim minimum, mean and means at 0.3 m and 1 m; equator and Apollo 17 "
        "with the option",
        "none: sensitivity (pre-registered 2026-11-16)",
        requires=(*_REFERENCE_REQUIRES,),
        slow=True,
        run=_th_c2,
    ),
)
