# Pre-registration: regional DEM error bounds (F14), ILL-C1

Written before the option is implemented, before the error map has been downloaded and
before anything below has been run.

Origin: finding F14 of `docs/review/scientific-review-2026-11-13.md`: the nested-DEM
horizon is sound; recommended action "Add regional-DEM error (PGDA 90 `ADJ_ERR`) for
distant skylines, since the reported spread is a lower bound there." Today the terrain
uncertainty comes only from the error clones of the 5 m site tile; the 80 m regional DEM
used beyond `near_radius_m` (and its block-averaged far level) is the same in every
member.

What is already known (not blind): the reference run's clone spread at the test point
(2 m, 2026): mean visible disc 5/50/95 % = 85.35 / 85.54 / 85.69 %. Nothing about the
regional error has been computed with lunarenv.

## Source (read 2026-11-16)

- PGDA product 90 page (https://pgda.gsfc.nasa.gov/products/90): `LDEM_80S_80MPP_ADJ_ERR.TIF`
  is the "LDEM surface height error (in meters)" of `LDEM_80S_80MPP_ADJ.TIF`, the regional
  DEM lunarenv already uses (255 MB, last modified 2023-06-03).
- Barker et al. (2023), PSJ 4, 183, doi:10.3847/PSJ/acf3e1, section 3 (open access, PDF
  SHA-256 d1362725cd46d4dda62e13d148434458b5e6636377365db5bc5d6520c9b04025): the maps are
  rms height uncertainty per pixel, converted from the effective-resolution map with a
  spline fitted to clone ensembles of four representative tiles (spline residual MAD
  about 20 %), "the uncertainty in the LDEM gridded height ... relative to the true
  average height ... within each pixel". No error correlation length or covariance is
  published with the maps, and regional clones are not distributed.

## Method (fixed now)

Because only a per-pixel rms is available, the regional error is propagated as a
**bound**, not as a distribution:

- Two extra terrain stacks: the regional DEM with every pixel raised by `k * ERR` and with
  every pixel lowered by `k * ERR`, before the far level is block-averaged. The site DEM
  (near level) and the observer are unchanged.
- Why this is a bound: the horizon elevation along a ray is the maximum over samples of
  an increasing function of each sample's height. For any regional height error with
  `|dh| <= k * ERR` pixel by pixel, whatever its correlation, the horizon lies between
  the lowered-stack and raised-stack horizons, at every azimuth. Visible disc fraction,
  sunlit time and the longest shadow are monotone in the horizon, so they are bounded
  too. The bound is conservative (it assumes the worst case at every pixel at once).
- `k = 2` is reported (about 95 % per pixel if the per-pixel error is Gaussian); `k = 1`
  is shown in ILL-C1 for scale.
- Scenario option `uncertainty.regional_error_k` (default 2.0; 0 disables). The error map
  is a separate catalog file (`pgda90.ldem_80s_80mpp_err`, group `regional-error`); if it
  is not in the data cache the run warns, as it does for missing clones, and reports no
  bound.
- Run output: mean visible disc fraction, sunlit time fraction (where reported), longest
  shadow and longest illumination for the lowered and raised stacks, with their changes
  from the nominal terrain, labelled "regional DEM error bound (+/- k ERR)". The thermal
  model is not rerun on these stacks.
- The error map must be on the same grid as the regional DEM; if not, the option is not
  implemented as described and this note is amended before any run.

## ILL-C1: effect of the regional DEM error (sensitivity, no pass/fail)

Status INFO by construction.

- Observers: the reference scenario's point at 2 m and at 0 m (2026, hourly), and the five
  Gläser 2018 points of ILL-V5 at their heights (2026, hourly; one year rather than 20,
  for cost).
- Default illumination settings and terrain stacks; nominal site DEMs.
- Per observer, for k = 2 (and k = 1): mean visible disc fraction and longest shadow on
  the lowered and raised stacks, and their change from the nominal; largest change of the
  horizon elevation over azimuth (degrees). For the site 4 points, the 5-95 % range of the
  same quantities over the terrain clones in the data cache, for comparison.

## Reading (fixed now)

- If at every observer the k = 2 bound changes the mean visible disc fraction by less than
  0.5 pp and the longest shadow by less than one time step, regional DEM error is
  negligible at these points and the F14 concern is closed for them; the run still reports
  the bound.
- If the bound is wider than the clone 5-95 % range at an observer, the regional terrain
  dominates the terrain uncertainty there and the README's "lower bound" warning is
  replaced by the bound itself.
- Limits: the bound is not a confidence interval. Errors beyond the 80 S product (rays that
  leave it) and the effect on temperatures are not covered.
