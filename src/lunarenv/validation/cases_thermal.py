"""Validation cases of the thermal model."""

from __future__ import annotations

import json
from importlib import resources

import numpy as np

from lunarenv.core.constants import STEFAN_BOLTZMANN_W_M2_K4 as SIGMA
from lunarenv.data.store import DataStore
from lunarenv.thermal.grid import geometric_grid, uniform_grid
from lunarenv.thermal.model import SYNODIC_MONTH_S, ThermalResult, ThermalSettings, run_thermal
from lunarenv.thermal.properties import RegolithProperties, kirchhoff
from lunarenv.thermal.solver import SchemeName, solve_conduction
from lunarenv.thermal.surface import analytic_absorbed_flux
from lunarenv.validation.core import Case, Check, Outcome, Status

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


def _th_x1(_: DataStore) -> Outcome:
    text = (
        resources.files("lunarenv.validation")
        .joinpath("reference/heat1d_reference.json")
        .read_text()
    )
    reference = json.loads(text)
    planet = reference["planet"]
    properties = RegolithProperties(scale_height_m=planet["H"])
    checks: list[Check] = []
    surface_ok = True
    deep_difference = []
    for case in reference["cases"]:
        if case["equilibration_years"] != 20:
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
    )
    status = None
    if surface_ok and not all(bool(c.passed) for c in checks):
        status = Status.PARTIAL
    return Outcome(
        checks,
        f"surface within 1 K at three latitudes; at depth {np.mean(deep_difference):+.1f} K",
        notes,
        status,
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
)
