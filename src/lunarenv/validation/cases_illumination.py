"""Validation cases of the illumination model."""

from __future__ import annotations

import os
from collections.abc import Callable
from importlib import resources

import numpy as np

from lunarenv.core.constants import MOON_REFERENCE_RADIUS_M as R
from lunarenv.data.dem import Bounds, load_dem, read_raster
from lunarenv.data.store import DataStore
from lunarenv.geometry.ephemeris import DE440, spice_session
from lunarenv.geometry.projection import SouthPolarStereographic, elevation_azimuth
from lunarenv.geometry.visibility import horizon_dip
from lunarenv.illumination.disc import CHORD, PROFILE, disc_fractions
from lunarenv.illumination.ensemble import regional_error_stacks
from lunarenv.illumination.horizon import compute_horizons
from lunarenv.illumination.products import illumination_statistics
from lunarenv.illumination.settings import IlluminationSettings, load_site_stack
from lunarenv.illumination.terrain import TerrainLevel, TerrainStack, build_stack, coarsen
from lunarenv.sites import TERRAIN_SETS
from lunarenv.validation.analytic import analytic_dem, reference_horizon
from lunarenv.validation.core import Case, Check, Outcome, Status

PROJ = SouthPolarStereographic()
_REGIONAL = "pgda90.ldem_80s_80mpp"
_REGIONAL_ERR = "pgda90.ldem_80s_80mpp_err"
_MONTHS = {
    m: i + 1 for i, m in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split())
}


# ------------------------------------------------------------------------- analytic
def _flat(height: float) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    return lambda x, y: np.full(np.shape(x), height)


def _ill_a1(_: DataStore) -> Outcome:
    worst = 0.0
    for terrain, above in ((0.0, 2.0), (1500.0, 10.0), (-2500.0, 50.0)):
        center = (4000.0, -6000.0)
        dem = analytic_dem(
            _flat(terrain),
            center=center,
            half_width_m=20_000.0,
            pixel_m=20.0,
        )
        stack = TerrainStack((TerrainLevel(dem, 0.0, 19_000.0, 10.0),))
        horizons = compute_horizons(stack, *center, height_above_surface_m=above, n_azimuth=36)
        expected = float(horizon_dip(above, R + terrain))
        worst = max(worst, float(np.abs(horizons.elevation_rad[0] - expected).max()))
    check = Check(
        quantity="max horizon error over three observer heights",
        reference=0.0,
        measured=worst,
        tolerance=1e-7,
        unit="rad",
    )
    return Outcome([check], f"max error {worst:.1e} rad")


def _hill(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.asarray(400.0 * np.exp(-((x - 3800.0) ** 2 + (y + 4500.0) ** 2) / (2.0 * 600.0**2)))


def _ill_a2(_: DataStore) -> Outcome:
    observer = (3000.0, -5000.0)
    dem = analytic_dem(_hill, center=(3500.0, -4500.0), half_width_m=6_000.0, pixel_m=5.0)
    stack = TerrainStack((TerrainLevel(dem, 0.0, 5_000.0, 2.5),))
    horizons = compute_horizons(stack, *observer, height_above_surface_m=2.0, n_azimuth=72)
    reference = reference_horizon(
        _hill,
        observer,
        float(horizons.observer_height_m[0]),
        horizons.azimuth_rad,
        max_distance_m=5_000.0,
        step_m=0.5,
    )
    error = float(np.abs(horizons.elevation_rad[0] - reference).max())
    check = Check(
        quantity="max horizon error over 72 azimuths",
        reference=0.0,
        measured=error,
        tolerance=2e-5,
        unit="rad",
    )
    return Outcome(
        [check],
        f"max error {error:.1e} rad",
        (
            "Tolerance is the bound of bilinear interpolation of a 5 m raster of this hill "
            "(3.5 mm in height at about 1 km).",
        ),
    )


def _ill_a3(_: DataStore) -> Outcome:
    rho, elevation = 0.00465, 0.02
    worst_area = worst_flux = 0.0
    for y0 in (-0.95, -0.5, 0.0, 0.3, 0.8, 0.99):
        horizon = np.full(1440, elevation + y0 * rho)
        geometric, flux = disc_fractions(elevation, 1.0, rho, horizon, 512, 0.6, PROFILE)
        root = np.sqrt(1.0 - y0 * y0)
        segment = (np.arccos(y0) - y0 * root) / np.pi
        uniform = 2.0 * (np.pi / 4.0 - 0.5 * (y0 * root + np.arcsin(y0)))
        darkened = 0.5 * np.pi * (2.0 / 3.0 - (y0 - y0**3 / 3.0))
        flux_exact = (0.4 * uniform + 0.6 * darkened) / (np.pi * 0.8)
        worst_area = max(worst_area, abs(geometric - segment))
        worst_flux = max(worst_flux, abs(flux - flux_exact))
    checks = [
        Check(
            quantity="visible area fraction, max error",
            reference=0.0,
            measured=worst_area,
            tolerance=2e-6,
        ),
        Check(
            quantity="limb-darkened flux fraction, max error",
            reference=0.0,
            measured=worst_flux,
            tolerance=2e-6,
        ),
    ]
    return Outcome(checks, f"area {worst_area:.1e}, flux {worst_flux:.1e} (512 strips)")


# ---------------------------------------------------------------------- JPL Horizons
def read_horizons_table(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """UTC stamps, azimuth and elevation (radians) from a stored Horizons CSV observer table."""
    text = resources.files("lunarenv.validation").joinpath(f"reference/{name}").read_text()
    lines = text.splitlines()
    start, stop = lines.index("$$SOE"), lines.index("$$EOE")
    stamps, azimuth, elevation = [], [], []
    for line in lines[start + 1 : stop]:
        fields = [f.strip() for f in line.split(",")]
        date, clock = fields[0].split()
        year, month, day = date.split("-")
        stamps.append(f"{year}-{_MONTHS[month]:02d}-{day}T{clock}")
        azimuth.append(float(fields[3]))
        elevation.append(float(fields[4]))
    return np.array(stamps), np.radians(azimuth), np.radians(elevation)


HORIZONS_SITE = (197.461177, -89.687048, 1731.7)  # E longitude, latitude (deg), height (m)


def horizons_separation(store: DataStore, body: str, fixture: str) -> tuple[np.ndarray, np.ndarray]:
    """Angular separation (rad) between lunarenv and Horizons, exact and fast observer paths."""
    stamps, azimuth_ref, elevation_ref = read_horizons_table(fixture)
    x, y = PROJ.forward(np.radians(HORIZONS_SITE[0]), np.radians(HORIZONS_SITE[1]))
    site = PROJ.to_cartesian(x, y, HORIZONS_SITE[2])
    basis = PROJ.local_basis(x, y)
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(stamps)
        exact = ephemeris.position_from_site(body, et, site)  # type: ignore[arg-type]
        approximate = ephemeris.position_from_center(body, et) - site  # type: ignore[arg-type]

    def separation(direction: np.ndarray) -> np.ndarray:
        elevation, grid_azimuth = elevation_azimuth(direction, *basis)
        d_azimuth = np.angle(
            np.exp(1j * (grid_azimuth - PROJ.grid_convergence(x, y) - azimuth_ref))
        )
        return np.asarray(np.hypot(elevation - elevation_ref, d_azimuth * np.cos(elevation_ref)))

    return separation(exact), separation(approximate)


def _ill_x1(store: DataStore) -> Outcome:
    checks = []
    for body, fixture in (
        ("SUN", "horizons_10_site04_rim.txt"),
        ("EARTH", "horizons_399_site04_rim.txt"),
    ):
        exact, _ = horizons_separation(store, body, fixture)
        checks.append(
            Check(
                quantity=f"{body.title()}: max angular difference, 288 epochs 2026-2028",
                reference=0.0,
                measured=float(exact.max()),
                tolerance=2e-6,
                unit="rad",
            )
        )
    notes = (
        "Reference: JPL Horizons API responses of 2026-11-12 (DE441, its 'MEAN_ME (high precision)' "
        "lunar frame), stored verbatim in validation/reference/. Horizons shares the JPL ephemerides "
        "with lunarenv, so this validates time conversion, frames, projection and topocentric "
        "reduction, not the ephemeris itself.",
    )
    worst = [float(c.measured or 0.0) for c in checks]
    return Outcome(checks, f"Sun {worst[0]:.1e} rad, Earth {worst[1]:.1e} rad", notes)


# --------------------------------------------------- published: PGDA average illumination
def _nodal_cycle(store: DataStore, start: str, hours: int):  # type: ignore[no-untyped-def]
    """Sun positions over one 18.6-year nodal cycle sampled every ``hours``."""
    count = int(18.6 * 365.25 * 24 / hours)
    utc = np.datetime64(start) + np.arange(count) * np.timedelta64(hours, "h")
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(utc)
        return et, ephemeris.position_from_center("SUN", et), ephemeris.body_radius_m("SUN")


def _pgda69_window(store: DataStore):  # type: ignore[no-untyped-def]
    """Published map and lunarenv statistics on the ILL-V2 window (shared by ILL-V2 and ILL-V2b)."""
    window = Bounds(-6000.0, -12000.0, 0.0, -6000.0)  # 100 x 100 cells of the published 60 m grid
    published, grid = read_raster(store.path("pgda69.avgvisib_85s_60m"), window, allow_nodata=True)
    # site DEM cropped to the published lattice, then block-averaged from 5 m to 60 m
    near = coarsen(
        load_dem(store, "pgda78.site04.surf", Bounds(-9000.0, -15000.0, 6960.0, 960.0)), 12
    )
    stack = build_stack(
        near, load_dem(store, _REGIONAL), near_radius_m=8000.0, mid_radius_m=100_000.0
    )
    xx, yy = np.meshgrid(grid.x_centers(), grid.y_centers())
    horizons = compute_horizons(stack, xx.ravel(), yy.ravel(), n_azimuth=720, diagnostics=False)
    et, sun, radius = _nodal_cycle(store, "2010-01-01T00:00", 2)
    stats = illumination_statistics(horizons, et, sun, radius, mode=CHORD)
    return 100.0 * published.astype(np.float64), stats, xx.shape


def _map_checks(ours: np.ndarray, theirs: np.ndarray) -> tuple[list[Check], str]:
    valid = np.isfinite(theirs)
    difference = (ours - theirs)[valid]
    correlation = float(np.corrcoef(ours[valid], theirs[valid])[0, 1])
    both_dark = float(np.mean((ours[valid] < 0.5) & (theirs[valid] < 0.5)))
    either_dark = float(np.mean((ours[valid] < 0.5) | (theirs[valid] < 0.5)))
    lit = valid & (theirs > 50.0)
    checks = [
        Check(
            quantity="mean difference (lunarenv - published)",
            reference=0.0,
            measured=float(difference.mean()),
            tolerance=2.0,
            unit="pp",
        ),
        Check(
            quantity="RMS difference",
            reference=0.0,
            measured=float(np.sqrt(np.mean(difference**2))),
            tolerance=6.0,
            unit="pp",
        ),
        Check(
            quantity="correlation coefficient", reference=1.0, measured=correlation, tolerance=0.05
        ),
        Check(
            quantity="published maximum in the window",
            reference=float(theirs[valid].max()),
            measured=float(ours[valid].max()),
            unit="%",
            informational=True,
        ),
        Check(
            quantity="mean over cells published above 50 %",
            reference=float(theirs[lit].mean()),
            measured=float(ours[lit].mean()),
            unit="%",
            informational=True,
        ),
        Check(
            quantity="cells below 0.5 % in both / in either",
            reference=either_dark,
            measured=both_dark,
            informational=True,
        ),
    ]
    bias = float(difference.mean())
    rms = float(np.sqrt(np.mean(difference**2)))
    summary = f"bias {bias:+.2f} pp, RMS {rms:.2f} pp, r = {correlation:.3f} over {int(valid.sum())} cells"
    return checks, summary


def _ill_v2(store: DataStore) -> Outcome:
    theirs, stats, shape = _pgda69_window(store)
    checks, summary = _map_checks(100.0 * stats.mean_visible_fraction.reshape(shape), theirs)
    notes = (
        "Reference: PGDA average solar illumination map AVGVISIB_85S_060M_201608 (method of "
        "Mazarico2011: hourly over one 18.6-year cycle on LOLA 60 m topography, 2016 processing).",
        "lunarenv was run like-for-like: Site 4 DEM block-averaged from 5 m to 60 m on the published "
        "lattice, observers at ground level, single-chord solar disc, 720 azimuths, 100 x 100 cells "
        "(6 x 6 km on the Shackleton rim), 2010-01-01 plus 18.6 years every 2 h (the 2 h step changes "
        "the mean by 0.004 pp; the published start date is not known and does not matter for a full cycle).",
        "Known differences that the tolerances were set to cover: DEM generation (2020 track-adjusted "
        "site model vs 2016 gridded product), far-field DEM version, and the observer height of the "
        "published map, which its documentation does not state.",
        "Tolerances were fixed before the comparison from the terrain sensitivity measured in Phase 2 "
        "(coarsening the near DEM twofold moves a rim point by 0.3-1.5 pp; terrain clones by 1 pp).",
        "The case passes, but not uniformly: on cells the published map puts above 50 %, lunarenv at "
        "ground level is about 6 pp darker. A follow-up run (2026-11-12, 6 h sampling) with the "
        "observer 2 m above the ground gave overall bias -0.7 pp, +1.2 pp on those cells and r = 0.950; "
        "at 10 m, +10 pp and r = 0.897. The difference therefore equals about 1.5 m of observer "
        "height. A likely cause is the treatment of terrain next to the observer (lunarenv samples the "
        "slope of the observer's own 60 m cell from 15 m away; a horizon built from neighbouring grid "
        "nodes does not), but the published method could not be checked, so this is not established. "
        "The registered comparison stays at ground level.",
        "Update 2026-11-13: this case compares unlike metrics. The PDS label of the published map "
        "defines its values as the percentage of time steps with any part of the disc visible, "
        "while this case uses the mean visible fraction. ILL-V2b (pre-registered) repeats the "
        "comparison with the label's metric: the difference on cells above 50 % falls from "
        "-5.9 pp to -1.8 pp and the overall bias becomes +0.5 pp. The metric therefore explains "
        "most of the deficit, and the observer-height explanation above is dropped. This case is "
        "kept as recorded; ILL-V2b is the like-for-like comparison.",
    )
    return Outcome(checks, summary, notes)


def _ill_v2b(store: DataStore) -> Outcome:
    theirs, stats, shape = _pgda69_window(store)
    checks, summary = _map_checks(100.0 * stats.sunlit_time_fraction.reshape(shape), theirs)
    lit = np.isfinite(theirs) & (theirs > 50.0)
    lit_mean = next(c for c in checks if c.quantity.startswith("mean over cells published above"))
    checks.append(
        Check(
            quantity="mean difference over cells published above 50 % (discriminating check)",
            reference=0.0,
            measured=float(lit_mean.measured - lit_mean.reference),  # type: ignore[operator]
            tolerance=2.0,
            unit="pp",
        )
    )
    notes = (
        "Same set-up as ILL-V2; the only change is the lunarenv statistic: the percentage of time "
        "steps in which any part of the solar disc is above the horizon, instead of the mean "
        "visible fraction of the disc.",
        "Why: the PDS label of AVGVISIB_85S_060M_201608 (read 2026-11-13 with TLS verification on) "
        "says the values 'indicate the percentage of timesteps where each pixel was sunlit (by any "
        "fraction of the solar disc)', while describing the image object as 'average solar disc "
        "visibility'. The label is ambiguous; this case decides between the two readings.",
        "Pre-registered in docs/review/preregistration-2026-11-13.md (finding F08 of the "
        "2026-11-13 scientific review) before it was implemented or run. The three ILL-V2 "
        "tolerances are unchanged; the added check (mean difference on cells published above "
        "50 % within 2 pp) is the reviewer's stated falsification criterion.",
        f"Cells published above 50 %: {int(lit.sum())}.",
    )
    return Outcome(checks, summary + " [any part of the disc visible]", notes)


# ------------------------------------------------------- published: Barker2021 Table 2
#: Barker2021 Table 2, rows "A" (statistic averaged over the region for each clone, then
#: percentile over 100 clones). Columns: region, X and Y of the centroid (km, as printed),
#: area (m2), then for dz = 1 m and 5 m: 1st percentile of average illumination (%),
#: 1st percentile of the longest continuous illumination period (days), 99th percentile of the
#: longest continuous shadow period (days).
BARKER_TABLE2 = (
    (1, 11.5825, 12.4625, 8625, (68.03, 63.35, 16.55), (89.26, 215.79, 5.29)),
    (2, 11.4875, 12.5575, 7550, (70.32, 113.80, 17.84), (88.51, 216.28, 7.31)),
    (3, 11.3525, 12.5025, 3825, (67.51, 117.86, 22.23), (84.53, 216.89, 8.41)),
    (4, 11.2875, 12.1675, 13300, (69.61, 136.49, 13.27), (88.12, 231.20, 4.00)),
    (5, 11.1025, 11.9425, 10075, (70.38, 142.17, 15.24), (85.71, 218.98, 4.52)),
    (6, 10.8025, 11.7375, 10000, (71.15, 81.01, 12.74), (87.35, 220.98, 4.16)),
    (7, 10.2275, 11.5525, 7075, (67.36, 100.23, 14.54), (82.79, 119.07, 5.42)),
)
_Z99 = 2.326  # normal quantile of the 1st / 99th percentile
_EMPIRICAL_MIN_MEMBERS = 50  # from here on percentiles are taken from the ensemble itself


def _region(mean: np.ndarray, seed: tuple[int, int], threshold: float) -> np.ndarray:
    """Pixels connected (8-neighbourhood) to ``seed`` with ``mean`` above ``threshold``."""
    mask = np.zeros(mean.shape, dtype=bool)
    if mean[seed] <= threshold:
        return mask
    stack = [seed]
    mask[seed] = True
    rows, cols = mean.shape
    while stack:
        r, c = stack.pop()
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                rr, cc = r + dr, c + dc
                if (
                    0 <= rr < rows
                    and 0 <= cc < cols
                    and not mask[rr, cc]
                    and mean[rr, cc] > threshold
                ):
                    mask[rr, cc] = True
                    stack.append((rr, cc))
    return mask


def _ill_v1(store: DataStore) -> Outcome:
    near = load_dem(store, "pgda78.site01.surf")
    stack = build_stack(
        near, load_dem(store, _REGIONAL), near_radius_m=8000.0, mid_radius_m=100_000.0
    )
    utc = np.arange(
        np.datetime64("2024-01-01T00:00"), np.datetime64("2026-01-01T00:00"), np.timedelta64(1, "h")
    )
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(utc)
        sun = ephemeris.position_from_center("SUN", et)
        radius = ephemeris.body_radius_m("SUN")

    def statistics(terrain: TerrainStack, x: np.ndarray, y: np.ndarray, height: float):  # type: ignore[no-untyped-def]
        horizons = compute_horizons(
            terrain, x, y, height_above_surface_m=height, n_azimuth=720, diagnostics=False
        )
        return illumination_statistics(horizons, et, sun, radius, mode=CHORD)

    # 1. reconstruct each region: pixels connected to the published centroid with > 70 % at 1 m
    regions = []
    notes = []
    for number, x_km, y_km, area, _, _ in BARKER_TABLE2:
        cx, cy = -1000.0 * x_km, -1000.0 * y_km
        half = 100.0
        while True:
            window = near.crop(Bounds.around(cx, cy, half))
            xx, yy = np.meshgrid(window.grid.x_centers(), window.grid.y_centers())
            stats = statistics(stack, xx.ravel(), yy.ravel(), 1.0)
            mean = stats.mean_visible_fraction.reshape(xx.shape)
            row, col = window.grid.fractional_index(cx, cy)
            seed = (round(float(row)), round(float(col)))
            mask = _region(mean, seed, 0.70)
            on_edge = mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any()
            if not on_edge or half >= 400.0:
                break
            half *= 2.0
        if not mask.any():
            notes.append(
                f"Region {number}: the published centroid pixel has {100 * mean[seed]:.1f} % in the "
                "nominal DEM, below the 70 % that defines the region; compared at the centroid pixel."
            )
            mask[seed] = True
        regions.append(
            {
                "number": number,
                "x": xx[mask],
                "y": yy[mask],
                "area": area,
                "our_area": 25.0 * int(mask.sum()),
                "edge": bool(on_edge),
            }
        )

    def region_means(terrain: TerrainStack, height: float) -> np.ndarray:
        x = np.concatenate([r["x"] for r in regions])
        y = np.concatenate([r["y"] for r in regions])
        s = statistics(terrain, x, y, height)
        out = []
        start = 0
        for r in regions:
            stop = start + r["x"].size
            out.append(
                [
                    100 * s.mean_visible_fraction[start:stop].mean(),
                    s.longest_illumination_s[start:stop].mean() / 86400,
                    s.longest_shadow_s[start:stop].mean() / 86400,
                ]
            )
            start = stop
        return np.array(out)  # (region, statistic)

    limit = int(os.environ.get("LUNARENV_VALIDATION_MAX_CLONES", "100"))
    clone_ids = [e.id for e in store.catalog.group("site01-clones") if store.is_available(e.id)][
        :limit
    ]
    checks: list[Check] = []
    for column, height in ((4, 1.0), (5, 5.0)):
        nominal = region_means(stack, height)
        members = np.array(
            [region_means(stack.with_near(load_dem(store, c)), height) for c in clone_ids]
        )
        for i, entry in enumerate(BARKER_TABLE2):
            published = entry[4] if column == 4 else entry[5]
            label = f"region {entry[0]}, dz = {height:g} m"
            checks.append(
                Check(
                    quantity=f"{label}: nominal DEM average illumination",
                    reference=published[0],
                    measured=float(nominal[i, 0]),
                    unit="%",
                    informational=True,
                )
            )
            if len(clone_ids) >= 8:
                if len(clone_ids) >= _EMPIRICAL_MIN_MEMBERS:
                    low = np.percentile(members[:, i], 1.0, axis=0)
                    high = np.percentile(members[:, i], 99.0, axis=0)
                else:  # too few members for a 1 % tail: normal approximation
                    mean, std = members[:, i].mean(axis=0), members[:, i].std(axis=0, ddof=1)
                    low, high = mean - _Z99 * std, mean + _Z99 * std
                checks.append(
                    Check(
                        quantity=f"{label}: 1st percentile of average illumination",
                        reference=published[0],
                        measured=float(low[0]),
                        tolerance=3.0,
                        unit="%",
                    )
                )
                checks.append(
                    Check(
                        quantity=f"{label}: clone minimum of average illumination",
                        reference=published[0],
                        measured=float(members[:, i, 0].min()),
                        unit="%",
                        informational=True,
                    )
                )
                checks.append(
                    Check(
                        quantity=f"{label}: 1st percentile of longest illumination period",
                        reference=published[1],
                        measured=float(low[1]),
                        unit="d",
                        informational=True,
                    )
                )
                checks.append(
                    Check(
                        quantity=f"{label}: 99th percentile of longest shadow period",
                        reference=published[2],
                        measured=float(high[2]),
                        unit="d",
                        informational=True,
                    )
                )
            if height == 1.0:
                checks.append(
                    Check(
                        quantity=f"region {entry[0]}: area above 70 % at 1 m",
                        reference=float(entry[3]),
                        measured=float(regions[i]["our_area"]),
                        unit="m2",
                        informational=True,
                    )
                )
    n = len(clone_ids)
    notes = [
        "Reference: Barker2021 Table 2, rows A, Site 1 (connecting ridge), 2024-01-01 to 2026-01-01, "
        "hourly, on the same 5 m LOLA site DEM and its error clones.",
        "The table prints centroid coordinates as positive kilometres; the Site 1 tile lies at negative "
        "x and y, so the centroids are taken as (-X, -Y). They fall on pixel centres of the tile.",
        "Regions are rebuilt here as the pixels connected to the published centroid with average "
        "illumination above 70 % at 1 m in the nominal DEM (the published definition); their areas are "
        "compared for information.",
        f"Published values are 1st/99th percentiles over 100 clones. Here {n} clone(s) were used. "
        f"With {_EMPIRICAL_MIN_MEMBERS} or more the percentiles are taken from the ensemble itself; with 8 to "
        f"{_EMPIRICAL_MIN_MEMBERS - 1} they are estimated as mean -/+ 2.326 standard deviations of the region means, "
        "which for 12 to 16 clones carries a sampling error of about half a standard deviation.",
        "Method differences: far-field terrain is the 2023 PGDA 80 m product (the paper used the "
        "earlier 80 m and 240 m LDEMs); single-chord disc with the Eddington limb-darkening law (the "
        "paper: limb darkening at 550 nm); longest-period statistics are discontinuous functionals and "
        "are reported without a tolerance.",
        "Tolerance for average illumination, fixed beforehand: 3 percentage points.",
        *notes,
    ]
    area_error = max(abs(r["our_area"] / r["area"] - 1.0) for r in regions)
    area_text = f"region areas within {100 * area_error:.1f} % of the published ones"
    toleranced = [c for c in checks if c.passed is not None]
    if not toleranced:
        status = Status.PARTIAL
        summary = (
            f"nominal DEM only ({n} clones available): {area_text}; no percentile comparison yet"
        )
    else:
        worst = max(abs(c.measured - c.reference) for c in toleranced)  # type: ignore[operator]
        failed = [c for c in toleranced if not c.passed]
        status = Status.FAIL if failed else (Status.PASS if n >= 100 else Status.PARTIAL)
        summary = (
            f"1st percentile of average illumination: largest difference {worst:.2f} pp over "
            f"{len(toleranced)} region-height pairs, from {n} of 100 clones; {area_text}"
        )
    return Outcome(checks, summary, tuple(notes), status)


# ------------------------------------------------------ published: Glaser2018 named points
#: Glaser2018 Tables 1 (surface) and 2 (2 m), points inside the cached 5 m site DEMs.
#: Columns: label, site DEM, observer height (m), latitude, east longitude (deg), average
#: illumination (%), solar visibility (%), maximum time in shadow (h).
GLASER2018 = (
    ("C1-2m", "pgda78.site01.surf", 2.0, -89.4398, 222.83, 88.0, 92.1, 112.0),
    ("S-2m-a", "pgda78.site04.surf", 2.0, -89.6866, 197.19, 85.5, 89.4, 66.0),
    ("S-2m-b", "pgda78.site04.surf", 2.0, -89.7815, 204.70, 82.1, 89.8, 65.0),
    ("C1-0m", "pgda78.site01.surf", 0.0, -89.4505, 222.84, 76.0, 81.3, 268.0),
    ("S-0m", "pgda78.site04.surf", 0.0, -89.7849, 203.97, 81.0, 85.5, 221.0),
)
_NEIGHBOURHOOD_M = 40.0  # two pixels of the source's 20 m DTM


def _ill_v5(store: DataStore) -> Outcome:
    settings = IlluminationSettings()
    utc = np.arange(
        np.datetime64("2017-01-01T00:00"), np.datetime64("2037-01-01T00:00"), np.timedelta64(1, "h")
    )
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(utc)
        sun = ephemeris.position_from_center("SUN", et)
        radius = ephemeris.body_radius_m("SUN")
    offsets = np.arange(-_NEIGHBOURHOOD_M, _NEIGHBOURHOOD_M + 1.0, 5.0)
    ox, oy = np.meshgrid(offsets, offsets)
    disc = np.hypot(ox, oy) <= _NEIGHBOURHOOD_M
    ox, oy = ox[disc], oy[disc]  # the point itself is offset (0, 0)
    centre = int(np.flatnonzero((ox == 0.0) & (oy == 0.0))[0])
    checks: list[Check] = []
    parts = []
    stacks: dict[str, TerrainStack] = {}
    for label, near_id, height, lat, lon, avg_ref, vis_ref, shadow_ref in GLASER2018:
        if near_id not in stacks:
            stacks[near_id] = load_site_stack(store, settings, near_id=near_id)
        x0, y0 = PROJ.forward(np.radians(lon), np.radians(lat))
        horizons = compute_horizons(
            stacks[near_id],
            float(x0) + ox,
            float(y0) + oy,
            height_above_surface_m=height,
            n_azimuth=settings.n_azimuth,
            diagnostics=False,
        )
        stats = illumination_statistics(
            horizons,
            et,
            sun,
            radius,
            n_strips=settings.n_strips,
            limb_u=settings.limb_darkening,
            mode=settings.disc_mode_code,
        )
        average = 100.0 * stats.mean_visible_fraction
        visibility = 100.0 * float(stats.sunlit_time_fraction[centre])
        shadow_h = float(stats.longest_shadow_s[centre]) / 3600.0
        checks += [
            Check(
                quantity=f"{label}: average illumination",
                reference=avg_ref,
                measured=float(average[centre]),
                tolerance=3.0,
                unit="%",
            ),
            Check(
                quantity=f"{label}: solar visibility",
                reference=vis_ref,
                measured=visibility,
                tolerance=3.0,
                unit="%",
            ),
            Check(
                quantity=f"{label}: maximum time in shadow",
                reference=shadow_ref,
                measured=shadow_h,
                unit="h",
                informational=True,
            ),
            Check(
                quantity=f"{label}: average illumination, max / min within 40 m",
                reference=float(average.max()),
                measured=float(average.min()),
                unit="%",
                informational=True,
            ),
        ]
        parts.append(f"{label} {average[centre]:.1f}/{visibility:.1f} %")
    notes = (
        "Reference: Glaser2018 Tables 1 (surface) and 2 (2 m above ground), read in full from the "
        "accepted manuscript on NTRS (20170007365); the final journal version was not checked. "
        "Source method: 20 m/pixel LOLA DTM, horizon method of Mazarico2011, 2017-01-01 to "
        "2037-01-01 sampled hourly. 'Average illumination' is the mean visible fraction of the "
        "disc (their section 2); 'solar visibility' is the share of time with any part visible.",
        "lunarenv: default illumination settings (1440 azimuths, per-strip disc profile, nested "
        "5 m / 80 m / coarsened regional terrain), nominal site DEM, the same period and step.",
        "Connecting Ridge C2 (both tables) lies outside the cached site DEMs and is not run.",
        "The published points are local maxima picked on a 20 m DTM; the source says such "
        "locations are 20-40 m wide. Picking a maximum biases the published value upwards "
        "relative to a point value on an independent DEM, so a lower lunarenv value is the "
        "expected direction of a miss. The 'max / min within 40 m' rows show the sensitivity to "
        "registration; they do not change the status.",
        "Pre-registered in docs/review/preregistration-2026-11-13.md (finding F13 of the "
        "2026-11-13 scientific review) before it was implemented or run. Tolerance 3 pp as "
        "ILL-V1; longest-shadow statistics carry no tolerance in this project.",
        "Exploratory, not pre-registered, does not change the status (2026-11-13): the same "
        "points with the site DEM block-averaged to 10 m and 20 m. At 2 m the values move by "
        "under 3 pp. At ground level they move a lot: S-0m 44.3 / 45.9 / 74.2 % and C1-0m "
        "37.1 / 44.3 / 49.2 % average illumination at 5 / 10 / 20 m, against 81.0 and 76.0 % "
        "published. Ground-level results depend strongly on terrain resolution, and the "
        "published values come from a different 20 m DTM, so the surface-level misses are not "
        "attributed to the engine or to the source; they show that 0 m outputs are not robust.",
    )
    return Outcome(checks, "; ".join(parts) + " (average / visibility)", notes)


# ------------------------------------------- sensitivity: regional DEM error (F14)
#: Observers of ILL-C1: label, site DEM, height (m), map x, y (m). The reference
#: scenario's point at 2 m and 0 m, then the Glaser2018 points of ILL-V5.
def _ill_c1_points() -> list[tuple[str, str, float, float, float]]:
    points = [
        ("test point 2 m", "pgda78.site04.surf", 2.0, -2822.5, -9062.5),
        ("test point 0 m", "pgda78.site04.surf", 0.0, -2822.5, -9062.5),
    ]
    for label, near_id, height, lat, lon, *_ in GLASER2018:
        x, y = PROJ.forward(np.radians(lon), np.radians(lat))
        points.append((label, near_id, height, float(x), float(y)))
    return points


def _ill_c1(store: DataStore) -> Outcome:
    settings = IlluminationSettings()
    utc = np.arange(
        np.datetime64("2026-01-01T00:00"), np.datetime64("2027-01-01T00:00"), np.timedelta64(1, "h")
    )
    with spice_session(store, DE440) as ephemeris:
        et = ephemeris.utc_to_et(utc)
        sun = ephemeris.position_from_center("SUN", et)
        radius = ephemeris.body_radius_m("SUN")
    regional = load_dem(store, _REGIONAL)
    error = load_dem(store, _REGIONAL_ERR)
    points = _ill_c1_points()
    # per point: variant -> (mean visible %, longest shadow h, horizon elevation rad)
    results: dict[str, dict[str, tuple[float, float, np.ndarray]]] = {p[0]: {} for p in points}
    clones: dict[str, list[tuple[float, float]]] = {p[0]: [] for p in points}

    def evaluate(
        stack: TerrainStack, members: list[tuple[str, str, float, float, float]]
    ) -> dict[str, tuple[float, float, np.ndarray]]:
        out: dict[str, tuple[float, float, np.ndarray]] = {}
        for height in sorted({m[2] for m in members}):
            group = [m for m in members if m[2] == height]
            horizons = compute_horizons(
                stack,
                np.array([m[3] for m in group]),
                np.array([m[4] for m in group]),
                height_above_surface_m=height,
                n_azimuth=settings.n_azimuth,
                diagnostics=False,
            )
            stats = illumination_statistics(
                horizons,
                et,
                sun,
                radius,
                n_strips=settings.n_strips,
                limb_u=settings.limb_darkening,
                mode=settings.disc_mode_code,
            )
            for i, m in enumerate(group):
                out[m[0]] = (
                    100.0 * float(stats.mean_visible_fraction[i]),
                    float(stats.longest_shadow_s[i]) / 3600.0,
                    np.asarray(horizons.elevation_rad[i], dtype=np.float64),
                )
        return out

    n_clones: dict[str, int] = {}
    for near_id in dict.fromkeys(p[1] for p in points):
        members = [p for p in points if p[1] == near_id]
        near = load_dem(store, near_id)
        nominal = build_stack(
            near,
            regional,
            near_radius_m=settings.near_radius_m,
            mid_radius_m=settings.mid_radius_m,
            step_fraction=settings.step_fraction,
            far_coarsening=settings.far_coarsening,
        )
        variants = {"nominal": nominal}
        for k in (1.0, 2.0):
            low, high = regional_error_stacks(near, regional, error, k, settings)
            variants[f"lowered {k:g}"], variants[f"raised {k:g}"] = low, high
        for name, stack in variants.items():
            for label, value in evaluate(stack, members).items():
                results[label][name] = value
        group = next(t.clone_group for t in TERRAIN_SETS.values() if t.near_id == near_id)
        clone_ids = [e.id for e in store.catalog.group(group) if store.is_available(e.id)]
        n_clones[near_id] = len(clone_ids)
        for clone_id in clone_ids:
            for label, value in evaluate(
                nominal.with_near(load_dem(store, clone_id)), members
            ).items():
                clones[label].append(value[:2])

    checks: list[Check] = []
    widest = (0.0, "")
    dominated = []
    for label, *_ in points:
        r = results[label]
        nom, low2, high2, low1, high1 = (
            r[v] for v in ("nominal", "lowered 2", "raised 2", "lowered 1", "raised 1")
        )
        clone = np.array(clones[label]) if clones[label] else np.empty((0, 2))
        width2 = low2[0] - high2[0]
        clone_width = (
            float(np.quantile(clone[:, 0], 0.95) - np.quantile(clone[:, 0], 0.05))
            if clone.size
            else float("nan")
        )
        if width2 > widest[0]:
            widest = (width2, label)
        if clone.size and width2 > clone_width:
            dominated.append(label)
        checks += [
            Check(
                quantity=f"{label}: mean visible disc, k = 2 raised / lowered regional terrain",
                reference=high2[0],
                measured=low2[0],
                unit="%",
                informational=True,
            ),
            Check(
                quantity=f"{label}: change of the mean visible disc from the nominal, k = 2 raised / lowered",
                reference=high2[0] - nom[0],
                measured=low2[0] - nom[0],
                unit="pp",
                informational=True,
            ),
            Check(
                quantity=f"{label}: width of the k = 2 bound / of the clone 5-95 % range, mean visible disc",
                reference=width2,
                measured=clone_width,
                unit="pp",
                informational=True,
            ),
            Check(
                quantity=f"{label}: width of the k = 1 bound, mean visible disc",
                measured=low1[0] - high1[0],
                unit="pp",
                informational=True,
            ),
            Check(
                quantity=f"{label}: longest shadow, k = 2 lowered / raised (nominal {nom[1]:.0f} h)",
                reference=low2[1],
                measured=high2[1],
                unit="h",
                informational=True,
            ),
            Check(
                quantity=f"{label}: longest shadow, clone 5 / 95 %",
                reference=float(np.quantile(clone[:, 1], 0.05)) if clone.size else None,
                measured=float(np.quantile(clone[:, 1], 0.95)) if clone.size else None,
                unit="h",
                informational=True,
            ),
            Check(
                quantity=f"{label}: largest horizon change over azimuth, k = 2 raised / lowered",
                reference=float(np.degrees(np.max(high2[2] - nom[2]))),
                measured=float(np.degrees(np.max(nom[2] - low2[2]))),
                unit="deg",
                informational=True,
            ),
        ]
    shadow_step = max(
        max(
            abs(results[p[0]]["raised 2"][1] - results[p[0]]["nominal"][1]),
            abs(results[p[0]]["lowered 2"][1] - results[p[0]]["nominal"][1]),
        )
        for p in points
    )
    summary = (
        f"k = 2 bound on the mean visible disc at most {widest[0]:.2f} pp wide ({widest[1]}); "
        f"longest shadow changes by at most {shadow_step:.0f} h; bound wider than the clone "
        f"5-95 % range at {len(dominated)} of {len(points)} points"
    )
    notes = (
        "Pre-registered in docs/review/preregistration-f14-2026-11-16.md (finding F14 of the "
        "2026-11-13 scientific review) before it was implemented or run.",
        "Error map: PGDA 90 LDEM_80S_80MPP_ADJ_ERR, rms height error of the 80 m regional DEM per "
        "pixel (Barker2023 section 3). Only a per-pixel rms is published, so the regional DEM is "
        "raised and lowered everywhere by k errors: a bound for any regional error within k "
        "errors, whatever its correlation, not a confidence interval.",
        "Observers: the reference scenario's point (2 m and 0 m) and the five Glaser2018 points of "
        "ILL-V5, 2026 hourly, default settings, nominal site DEMs. Clone ranges: "
        + ", ".join(f"{n} clones of {i}" for i, n in n_clones.items())
        + " (near field only).",
        "Points where the bound is wider than the clone range: "
        + (", ".join(dominated) or "none")
        + ".",
    )
    return Outcome(checks, summary, notes, Status.INFO)


_PAYWALL = "The paper is behind a paywall and was not accessed; only its abstract was read."

CASES: tuple[Case, ...] = (
    Case(
        "ILL-A1",
        "Horizon of a smooth sphere (dip)",
        "analytic",
        "closed form",
        "max |horizon - dip| over all azimuths",
        "1e-7 rad",
        run=_ill_a1,
    ),
    Case(
        "ILL-A2",
        "Horizon behind an analytic hill, off-pole observer",
        "analytic",
        "independent geodesic computation from true azimuths",
        "max horizon error",
        "2e-5 rad",
        run=_ill_a2,
    ),
    Case(
        "ILL-A3",
        "Visible fraction of a disc cut by a straight horizon",
        "analytic",
        "circular segment; limb-darkened segment in closed form",
        "max error of area and flux fractions",
        "2e-6",
        run=_ill_a3,
    ),
    Case(
        "ILL-X1",
        "Apparent direction of the Sun and the Earth at the rim site",
        "cross-check",
        "JPL Horizons (HORIZONS)",
        "max angular separation",
        "2e-6 rad",
        requires=DE440.kernel_ids,
        run=_ill_x1,
    ),
    Case(
        "ILL-V1",
        "Illumination of the Site 1 regions of interest",
        "published",
        "Barker2021 Table 2",
        "1st percentile over terrain clones of the region-average illumination, dz = 1 and 5 m",
        "3 percentage points",
        requires=("pgda78.site01.surf", _REGIONAL, *DE440.kernel_ids),
        slow=True,
        run=_ill_v1,
    ),
    Case(
        "ILL-V2",
        "18.6-year average illumination map, Shackleton rim",
        "published",
        "PGDA-69 (method of Mazarico2011)",
        "bias, RMS and correlation over 100 x 100 cells of 60 m",
        "bias 2 pp, RMS 6 pp, r >= 0.95",
        requires=("pgda69.avgvisib_85s_60m", "pgda78.site04.surf", _REGIONAL, *DE440.kernel_ids),
        slow=True,
        run=_ill_v2,
    ),
    Case(
        "ILL-V2b",
        "18.6-year average illumination map, metric of the published label",
        "published",
        "PGDA-69 (method of Mazarico2011)",
        "as ILL-V2, with lunarenv's percentage of time steps with any part of the disc visible",
        "bias 2 pp, RMS 6 pp, r >= 0.95; cells above 50 %: 2 pp (pre-registered 2026-11-13)",
        requires=("pgda69.avgvisib_85s_60m", "pgda78.site04.surf", _REGIONAL, *DE440.kernel_ids),
        slow=True,
        run=_ill_v2b,
    ),
    Case(
        "ILL-V3",
        "Continuous sunlight and longest darkness near the Shackleton rim",
        "published",
        "Mazarico2011",
        "days of continuous sunlight; longest total darkness",
        "not set",
        pending_reason=(
            _PAYWALL
            + " The abstract gives about 240 days of continuous sunlight and 1.5 days of longest "
            "darkness for a point near the rim, without its coordinates or observer height. Needs: Table 3 "
            "of the paper."
        ),
    ),
    Case(
        "ILL-V4",
        "Best-lit point of the connecting ridge over 19 years",
        "published",
        "Glaser2014",
        "average illumination at 2 m and 10 m",
        "not set",
        pending_reason=(
            _PAYWALL
            + " The abstract gives 92.27 % at 2 m and 95.65 % at 10 m on a 20 m terrain model, "
            "without the coordinates. Needs: the location from the paper."
        ),
    ),
    Case(
        "ILL-V5",
        "Named best-lit points at the south pole over 20 years",
        "published",
        "Glaser2018 Tables 1-2",
        "average illumination and solar visibility at the point, surface and 2 m",
        "3 percentage points (as ILL-V1; pre-registered 2026-11-13)",
        requires=(
            "pgda78.site01.surf",
            "pgda78.site04.surf",
            _REGIONAL,
            *DE440.kernel_ids,
        ),
        slow=True,
        run=_ill_v5,
    ),
    Case(
        "ILL-C1",
        "Effect of the regional DEM error on illumination (bound)",
        "cross-check",
        "PGDA 90 rms height error map (Barker2023); terrain error clones for comparison",
        "mean visible disc and longest shadow with the regional DEM raised and lowered by 1 and 2 "
        "errors, against the clone 5-95 % range",
        "none: sensitivity (pre-registered 2026-11-16)",
        requires=(
            "pgda78.site01.surf",
            "pgda78.site04.surf",
            _REGIONAL,
            _REGIONAL_ERR,
            *DE440.kernel_ids,
        ),
        slow=True,
        run=_ill_c1,
    ),
)

__all__ = ["BARKER_TABLE2", "CASES", "horizons_separation", "read_horizons_table"]
