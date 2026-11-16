# Illumination model

Version: `illumination 0.1` (Phase 2). Reference keys resolve in
[../references.md](../references.md). Numbers marked *measured* come from the scripts in
`benchmarks/` and the test-suite, run on 2026-11-12 on an Apple M1 Pro (8 performance +
2 efficiency cores).

## What it computes

For an observer at a map position and a height above the terrain:

- the **horizon**: elevation of the terrain skyline in `n_azimuth` directions;
- at each epoch, the apparent position of the Sun (and the Earth), its angular radius,
  and the **visible fraction** of its disc, geometric and flux-weighted;
- over a period: mean illumination, longest continuous shadow and illumination
  periods, for one point or for every pixel of a window;
- the distribution of those results over the terrain error ensemble.

Definitions of the statistics follow Barker2021 (section 5).

## Method

1. **Terrain stack.** Three DEM levels: the 5 m site product near the observer, the
   80 m regional product in the mid field, the same product block-averaged to 240 m in
   the far field. A ray that leaves a level early falls through to the next one.
2. **Horizon.** For each azimuth a ray follows the great circle through the observer.
   At each sample the terrain height is interpolated bilinearly and its elevation
   computed from exact spherical geometry,
   `tan e = ((R + h) cos psi - (R + h_obs)) / ((R + h) sin psi)`. The horizon is the
   maximum. Curvature is part of the geometry, not a correction.
3. **Stopping rule.** A ray ends when terrain as high as the highest point of the stack
   could no longer rise above the horizon found so far. Rays that run out of data
   before that are flagged; a second test says whether a flagged ray could matter for
   the Sun (it must rise above the lowest elevation the solar disc ever reaches at the
   site).
4. **Body direction.** Apparent position from SPICE (light-time and stellar
   aberration), DE440, in the frame `MOON_ME_DE440_ME421`, minus the observer's
   position. Angular radius from the actual distance and the IAU nominal solar radius.
5. **Visible fraction.** The disc is cut into vertical strips; in each strip the
   terrain hides what lies below the horizon interpolated at that strip's azimuth
   ("profile" mode). Limb darkening follows the linear law with `u = 0.6` (Eddington
   grey atmosphere); the flux-weighted fraction is what the thermal model uses. A
   "chord" mode reproduces the single horizontal cut of Mazarico2011/Barker2021.

## Alternatives considered

| Choice | Alternative | Why this one |
|---|---|---|
| Precomputed horizons | Ray tracing a mesh at every epoch | One-off cost, then any epoch is a lookup: a year, a nodal cycle or 100 terrain clones cost seconds per point. Same method as the published references, so validation is like-for-like. |
| Great-circle rays, 3-D elevation | Straight map lines plus curvature term | Exact on the sphere at no extra cost; no small-angle argument to defend. |
| Nested DEM levels | One resolution everywhere | 5 m out to 300 km is 1e5 samples per ray; nesting brings it to about 1e4. Measured penalty below. |
| Strip integration across the disc | One chord at the centre azimuth | Resolves skyline structure across the 0.53 deg solar width. Measured difference in yearly mean illumination: about 0.001 percentage points; larger at individual epochs. |
| Linear limb darkening, closed form | Uniform disc; multi-term laws | Changes the flux during partial occultation; closed-form strip integral, so no extra cost. Higher-order laws would refine a correction that is already small (see table). |

## Assumptions and simplifications

| Assumption | Expected effect |
|---|---|
| Terrain is the DEM, bilinear between pixel centres. About 90 % of the 5 m pixels are interpolated between LOLA tracks (Barker2021). | Dominant uncertainty. Handled statistically with the clone ensemble, not removed. |
| No relief below the DEM resolution (boulders, small craters). | Real skylines are rougher; shadows from unresolved relief are missing. Matters most for observers at zero height. |
| The Moon is a sphere; the vertical is radial. | Consistent with the LOLA products. |
| Sky around the disc treated as flat; horizon linear between azimuth samples. | Order `rho^2` = 2e-5 on the fraction; azimuth sampling error measured below. |
| Terrain beyond the regional DEM (more than 287 km away) ignored. | Proven irrelevant for the Sun for every observer above -1426 m in the tile ([data.md](../data.md)); rays are checked individually. |
| Illumination means direct sunlight only. Light scattered by terrain is not included. | Irrelevant for the visible fraction; matters for the energy budget in shadow (thermal model). |
| No diffraction, no solar corona, no eclipse by the Earth. | Earth eclipses of the Sun (lunar eclipses) are rare, hours long, and not modelled: the model reports sunlight during them. Flagged as a known limitation. |
| "Lit" means any visible fraction above zero. | Longest-period statistics depend on this threshold; a sliver counts. |

## Numerical error and uncertainty

Full study: [../convergence/illumination.md](../convergence/illumination.md).

- **Numerical error at default settings** (estimate, not a bound): yearly mean
  illumination up to +0.12 percentage points, biased towards more light; horizon
  30-65 arcsec RMS (the Sun's radius is 959 arcsec); period lengths resolved to one time
  step at each end. The largest term is the 80 m resolution of the mid-field terrain,
  which converges at first order.
- **Terrain uncertainty** is propagated by running the model on the PGDA error clones.
  On the rim crest at 2 m it is about +/-1 percentage point on the yearly mean and the
  nominal DEM sits at the optimistic end of the ensemble. At 10 m it falls to 0.01
  points, below the numerical estimate.
- **Regional terrain error** (review finding F14; ILL-C1, pre-registered 2026-11-16) is
  not available as clones, only as the PGDA 90 per-pixel rms height error map of the
  80 m regional DEM (Barker2023 section 3; median 0.29 m, largest 9.0 m). Runs bound it:
  the regional DEM is raised and lowered everywhere by `k` errors (`regional_error_k`,
  default 2), which bounds the horizon, and with it the visible fraction and the shadow
  and illumination periods, for any regional error within `k` errors whatever its
  correlation. It is a bound, not a confidence interval. At the reference point and the
  five Gläser 2018 points (2026) the k = 2 bound moves the yearly mean by at most
  0.013 pp and the horizon by at most 0.003 degrees, and changes no longest shadow; the
  site clones' 5-95 % range is at least 30 times wider. Regional error is negligible at
  these points.
- **Not represented**: relief below the DEM resolution; errors of terrain beyond the
  80 S product.
- Results carry both numbers; neither is folded into the other.

## Performance

Measured in CPU time, which does not depend on what else the machine is running (the
benchmarks ran while other workloads kept the load average above 25, so wall-clock
figures from this session are not representative). Scripts:
`benchmarks/illumination_timing.py`; raw numbers in `../benchmarks/`.

| Operation | Cost (CPU) | Scaling |
|---|---|---|
| Horizon, one observer, 1440 azimuths, 0.5 pixel step | 0.10-0.13 s (72-88 microseconds per ray, about 7000 terrain samples each) | linear in azimuths; 48 microseconds per ray at 1 pixel step, 153 at 0.25 |
| Visibility statistics | 0.95 microseconds per observer and epoch | linear in epochs: 8 ms for a year hourly, 0.16 s for 18.6 years |
| Ephemeris, 163 056 epochs (18.6 years hourly) | 4 s, single thread | linear |
| One point, one terrain member, one year | 1.2 s including reading the 41 MB DEM | 100 clones: about 2 CPU-minutes |
| Map, 2 x 2 km at 5 m (160 000 observers), 1440 azimuths | about 4.7 CPU-hours (extrapolated from 1000-4000 observers; cost is linear in observers) | about 35 minutes on eight idle performance cores; half with 720 azimuths |

Parallelism is over observers. Their cost differs by an order of magnitude (rays end
early inside craters), so observers are handed to threads in small blocks; with one
block per thread most cores sat idle.

## Range of validity

- South polar region covered by the stack: observers inside the Site 4 tile.
- Epochs 1849-2150 (DE440 short file); validated against JPL Horizons for 2026-2028.
- Observer height above the surface: 0 m and up. Results at 0 m depend on relief the
  DEM does not contain; heights of a metre or more are better supported. Evidence
  (ILL-V5): at 2 m three Gläser2018 points agree within 2.1 pp; at 0 m two points are
  37-39 pp low and move by up to 30 pp between 5 m and 20 m terrain. 0 m results are
  not validated.
- Outside these, the code runs but results are not validated.

## Verification (see VALIDATION.md for the full matrix)

| Case | Reference | Result |
|---|---|---|
| ILL-A1 smooth sphere | horizon dip, closed form | max error < 1e-7 rad at three heights |
| ILL-A2 Gaussian hill, off-pole observer | independent geodesic computation from true azimuths | max error < 2e-5 rad (bound set by bilinear interpolation of a 5 m raster) |
| ILL-A3 straight horizon | circular segment; limb-darkened segment in closed form | error < 2e-6 on the fraction with 512 strips |
| ILL-X1 Sun and Earth direction at the rim site | JPL Horizons (DE441, its own lunar frame and code), 288 epochs 2026-2028 | Sun 9.3e-10 rad, Earth 9.4e-9 rad (tolerance 2e-6 rad, fixed beforehand) |
| Properties | fraction in [0, 1]; raising terrain never adds light; mirror symmetry; fast and exact body direction agree | hold over randomized inputs |

Comparisons with published illumination results (Barker2021 Table 2 in ILL-V1, the
PGDA average illumination map in ILL-V2 and ILL-V2b, Glaser2018 named points in
ILL-V5) are in VALIDATION.md. The PGDA map counts time steps with any part of the disc
visible (its PDS label); ILL-V2b compares like for like.
