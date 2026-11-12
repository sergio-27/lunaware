# Illumination model: convergence, sensitivity and numerical error

Scripts: `benchmarks/illumination_convergence.py`, `benchmarks/illumination_ensemble.py`.
Raw output: [illumination.json](illumination.json),
[illumination_ensemble.json](illumination_ensemble.json); full tables:
[illumination_tables.md](illumination_tables.md). Run 2026-11-12.

Setup. Three observers in the Site 4 tile, 2 m above the surface, year 2026:
**A** rim crest, best lit (mean illumination 84.3 %); **B** rim crest, north (83.6 %);
**C** inner flank, poorly lit (8.9 %). Each study varies one parameter around the
baseline and compares with the finest value tried. Metrics: horizon profile difference
(arcsec; the Sun's radius is 959 arcsec), mean visible fraction of the solar disc
(percentage points, pp), longest continuous shadow period (hours).

## Results at a glance

Differences from the finest setting, worst of the three sites unless stated.

| Parameter | Default | Horizon RMS at default | Mean illumination at default | Observed behaviour |
|---|---|---|---|---|
| Azimuths per horizon | 1440 (0.25 deg) | 10 arcsec | 0.001 pp | RMS falls about 2.5x per doubling (180 to 5760). At 360 and fewer the longest shadow changes by up to 19 h. |
| Ray step | 0.5 pixel | 37 arcsec | +0.030 pp | Roughly first order. Coarser steps miss crests, so the horizon is biased low and illumination high. |
| Site DEM used out to | 8 km (tile limit) | reference | reference | 0.5 to 2 km all give the same error (40 arcsec, +0.05 pp); it drops from 4 km on. The 5 m data matter several kilometres out. |
| Mid-field DEM resolution | 80 m | 26 arcsec vs 40 m | +0.045 pp vs 40 m | First order in pixel size: 240, 160, 80 m differ from the next finer by 0.157, 0.078, 0.037 pp (site A). Extrapolated error of 80 m: about +0.07 to +0.09 pp. |
| Mid field extends to | 100 km | 3.6 arcsec | +0.009 pp | Nothing changes beyond 200 km. |
| Far-field block averaging | 3 (240 m) | included above | included above | 480 m: +0.03 pp. |
| Disc strips | 32 | - | < 0.001 pp | Converged at 16. |
| Time step | 60 min | - | 0.001 pp | Longest shadow within 0.4 h at 60 min. At 240 min it is wrong by up to 24 days at site C: short sunlit intervals are skipped and separate nights merge. |

Model-form sensitivities (not errors of discretization):

| Change | Effect on the 2026 mean |
|---|---|
| Single chord across the disc (Mazarico2011) instead of per-strip horizon | 0.001 pp |
| Uniform disc instead of limb darkening u = 0.6 | flux-weighted mean changes by up to 0.011 pp |
| Site DEM block-averaged to 10 m / 20 m / 40 m | +0.28 / +0.56 / +1.49 pp at site A; horizon RMS 205 / 308 / 702 arcsec |

The last row is the important one. Halving the resolution of the near-field terrain
changes the answer ten times more than all numerical parameters together. The model is
limited by knowledge of the terrain, not by its numerics.

## Numerical error estimate attached to results

Adding the measured contributions at the default settings (they share a sign, so they
are added linearly, not in quadrature):

| Quantity | Estimate | Dominant terms |
|---|---|---|
| Mean illumination over a year | up to +0.12 pp, biased towards more light | mid-field resolution (+0.07 to +0.09), ray step (+0.03), far field (+0.01) |
| Horizon profile | 30 to 65 arcsec RMS (3 to 7 % of the solar radius) | mid-field resolution, ray step |
| Longest shadow / illumination period | one time step at each end; up to 4 h more at sites where a short sunlit interval is marginal | time step, mid-field resolution |

This estimate is extrapolated from the first-order trend in mid-field resolution; it
is not a bound. Using the 40 m regional product for the mid field halves the largest
term at about 35 % more computing time.

## Why these defaults

- **1440 azimuths.** 720 would be enough for mean illumination (0.003 pp) but changed
  the longest shadow at one of three sites by 4 h. 1440 changed nothing.
- **0.5 pixel ray step.** 0.25 pixel reduces the bias from +0.030 to +0.009 pp and
  doubles the cost; both are below the mid-field term, so the cheaper one is kept.
- **80 m mid field.** The 40 m product is 716 MB instead of 189 MB and needs about
  1 GB of memory when loaded whole. The gain (0.04 pp) is small against the terrain
  uncertainty at low observer heights. It remains selectable.
- **60 min time step** for statistics, as in Barker2021. Period lengths are quantized
  to it.

## Terrain uncertainty (error clones)

All 100 Site 4 clones (regenerated on 2026-11-12 once the download completed; the first
version of this table used 13).

| Site, observer height | Nominal DEM | Ensemble mean +/- std | 1 % to 99 % | Longest shadow: nominal; ensemble |
|---|---|---|---|---|
| A, 2 m | 84.30 % | 83.17 +/- 0.86 % | 80.61 to 84.32 % | 2.67 d; 2.79 +/- 0.56 d |
| B, 2 m | 83.59 % | 83.53 +/- 0.18 % | 83.19 to 83.90 % | 6.42 d; 6.86 +/- 0.41 d |
| C, 2 m | 8.86 % | 8.86 +/- 0.00 % | 8.85 to 8.87 % | 176.5 d; 176.5 +/- 0.0 d |
| A, 10 m | 86.34 % | 86.32 +/- 0.01 % | 86.29 to 86.34 % | 2.62 d; 2.62 +/- 0.00 d |
| B, 10 m | 87.03 % | 87.03 +/- 0.01 % | 87.00 to 87.05 % | 2.54 d; 2.54 +/- 0.00 d |
| C, 10 m | 9.07 % | 9.06 +/- 0.00 % | 9.06 to 9.07 % | 176.5 d; 176.5 +/- 0.0 d |

Reading:

- At 2 m on the rim crest, terrain height error moves the yearly mean by about 1 pp and
  the nominal DEM sits at the optimistic end of the ensemble (site A: nominal 84.30 %,
  ensemble median 83.37 %, 99th percentile 84.32 %), as Barker2021 reports: one high
  pixel raises a horizon, lowering it needs many low ones.
- At 10 m the nearby terrain no longer forms the skyline and the ensemble spread falls
  to 0.01 pp. There the reported uncertainty is floored by the numerical estimate above
  (about 0.1 pp), and by errors of the regional DEM, which the clones do not sample.
- Where the skyline is a distant wall (site C), near-field errors do not matter.

Convergence with ensemble size (site A, 2 m; mean +/- std): 13 clones 83.30 +/- 0.97 %,
25 clones 83.19 +/- 0.87 %, 50 clones 83.07 +/- 0.91 %, 100 clones 83.17 +/- 0.86 %.
Mean and spread are stable to about 0.1 pp from 25 clones on; the default of 16 clones
in the reference scenario estimates the spread to roughly 20 %.
