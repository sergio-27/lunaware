# Pre-registration: TH-A3 (bowl crater, single bounce) and TH-C1 (terrain flux at the test-point minimum)

Written before the terrain-radiation code exists and before anything below has been run.

Origin: change 3 of `docs/review/scientific-review-2026-11-13.md` (finding F06). The
review asks for single-bounce terrain-scattered sunlight, terrain infrared and
Earthshine, and states how to tell whether that worked: "(a) the bowl-crater analytic
case (Ingersoll-type irradiance ...) is reproduced; (b) the reported flux at the
test-point minimum agrees with the 0.02 W m-2 bound", and how to tell it is wrong: "the
single-bounce flux at the test-point minimum exceeds about 0.05 W m-2".

## Method to be implemented (fixed now, so the cases test it rather than shape it)

- Visible terrain: rays from the receiver over the terrain stack, with the same
  geometry and sample plan as the horizon tracer. A sample is visible when its elevation
  exceeds every closer sample on the ray; it fills the elevation interval between that
  previous maximum and its own elevation (the first sample fills down from -90 deg).
- Patches: visible samples grouped by azimuth sector and distance band around the site.
  Each patch is represented by its visible sample of largest solid angle; that point gets
  its own horizon, DEM normal, direct sunlight (same disc model and Earth-eclipse factor
  as the site) and 1-D temperature (same thermal model and properties as the receiver).
- Patch radiance: Lambertian. Reflected sunlight A(theta) E cos(theta), emitted
  epsilon sigma T^4. Patches receive direct sunlight only (single bounce).
- Receiver irradiance: the radiance of each sample's patch times the cosine-weighted
  solid angle of the sample's elevation interval, clipped at the receiver plane
  (exact in elevation, midpoint in azimuth), for any receiver normal.
- Absorbed: (1 - A_diffuse) E_vis + epsilon E_IR, with A_diffuse the hemispherical
  average of the receiver's albedo law under isotropic illumination.
- Earthshine, HLS-UG-001 section 11.1 (NTRS 20210010030, read 2026-11-16): normal
  irradiance (0.100 gamma + 0.056) W m-2 at full Earth (gamma = 1); the guidebook's peak
  values, reflected and infrared. gamma = illuminated fraction of the Earth disc from the
  ephemeris (Sun-Earth-site phase angle). Scaled by (384400 km / d)^2; that the
  guidebook's values refer to the mean distance is an assumption, not stated there. Times
  the visible fraction of the Earth disc above the terrain horizon and the cosine of
  incidence on the receiver.

## TH-A3: spherical bowl, single bounce (analytic)

- Terrain: a spherical bowl of rim radius r = 500 m and depth H = 200 m (sphere radius
  R = (r^2 + H^2) / 2H = 725 m) cut into level ground, on a synthetic DEM of 2.5 m
  pixels, one terrain level.
- Sun: fixed direction, elevation 10 deg at the bowl centre, grid azimuth 0, 1 au,
  solar disc as in the illumination model. Constant albedo 0.12 (albedo_a = albedo_b =
  0). Patches in instantaneous radiative equilibrium (zero thermal inertia).
- Receivers: the bowl centre and the floor points 150 m either side of it along the
  map x axis (across the Sun direction), normal from the DEM. All three must be in
  shadow according to the model's own horizon; otherwise the case is an error.
- Reference: in a sphere the view factor between any two interior elements is
  dA / (4 pi R^2), so a shadowed interior point receives the power entering the opening,
  F sin(e) pi r^2, scattered or re-emitted once, divided by 4 pi R^2:
  E_vis = A F sin(e) r^2 / (4 R^2) and, in equilibrium, E_IR = (1 - A) F sin(e) r^2 /
  (4 R^2) (Ingersoll et al. 1992, Icarus 100, 40-47, single-scattering term).
- Settings for the check: 180 azimuth sectors, 10 m distance bands, 720 receiver rays,
  360 azimuths for patch horizons.
- Checks: E_vis and E_IR at each of the three receivers within 5 % of the reference (six
  checks). PASS if all pass.
- Informational: the same with the default patch settings used for scenarios.

## TH-C1: terrain flux at the minimum of the reference scenario

- Run `scenarios/shackleton_rim_2026.yaml` (nominal terrain; ground receiver at the
  site with the DEM normal) with terrain radiation, default patch settings.
- For each albedo case (cold, scenario law, hot): the single-bounce terrain irradiance
  E_vis + E_IR at the time of that case's minimum surface temperature.
- Checks (PASS/FAIL), each albedo case: (i) at most the run's neglected-flux upper bound
  (0.0217 W m-2 in the current reference run; the bound is recomputed by the same run);
  (ii) at most 0.05 W m-2, the review's falsification threshold. PASS if all six hold.
- Informational: Earthshine absorbed at the minimum; change of the minimum and mean
  surface temperature from terrain radiation and from Earthshine; the same irradiance
  with twice the sectors and bands (patch convergence).

## Reading

- TH-A3 PASS: the geometry, view factors and single-bounce bookkeeping reproduce the
  analytic bowl. A miss means the method or its implementation is wrong for concave
  terrain; it is reported as FAIL and the terrain term is not used to make claims.
- TH-C1 PASS: the earlier bound was consistent and the rim minimum moves by about 1 K or
  less from terrain radiation. A flux above the bound would mean the bound was wrong; above
  0.05 W m-2 the review's statement that the minimum was wrong by more than about 2 K
  holds.
- Not tested here: multiple scattering, sub-DEM roughness, non-Lambertian scattering,
  shadowing of the terrain by the asset itself.
