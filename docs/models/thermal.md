# Thermal model

Version: `thermal 0.1` (Phase 3). Reference keys resolve in
[../references.md](../references.md). Numbers marked *measured* come from
`benchmarks/thermal_site.py` and the test-suite, run on 2026-11-12.

## What it computes

Temperature of the regolith at a point, at the surface and at depth, over a period:
time series, extremes, mean, amplitude of the thermal cycles and fastest heating and
cooling rates, with the uncertainty that follows from the regolith parameters and from
the terrain.

It is the thermal environment of the *ground*. A structure has its own geometry,
optical properties and heat capacity; what this model provides for it is the solar
input, the ground temperature it sits on and radiates to, and how both vary. Asset
surfaces of any orientation get their direct solar irradiance and view factors (see
"Bounds" below), not a temperature.

## Physics

One-dimensional heat conduction, `rho(z) c_p(T) dT/dt = d/dz (K(T, rho) dT/dz)`, with

- density rising exponentially with depth, `rho_s = 1100` to `rho_d = 1800 kg m-3`,
  scale `H = 0.06 m`;
- conductivity `K = K_c(rho) [1 + chi (T/350)^3]`, the second term being radiative
  transfer between grains (alternative: `conductivity_model: woods_robinson2019`, which
  replaces the 1 by the WoodsRobinson2019 temperature shape `g(T)` below 150 K; see
  "Low-temperature conductivity" below);
- heat capacity as a quartic polynomial in temperature (alternative:
  `heat_capacity_model: woods_robinson2019`, the WoodsRobinson2019 estimate as tabulated
  by DSNE);
- surface: absorbed sunlight plus conduction from below balances thermal emission,
  `Q_abs + K dT/dz = eps sigma T_s^4`;
- base: interior heat flow `Q = 0.018 W m-2`.

All laws and default values are those of Hayne2017 (Appendix A, Table A1).

Absorbed sunlight: `Q_abs = (1 - A(theta)) E cos(theta) f_eclipse`, with `E` the direct
irradiance from the illumination model (distance, visible disc fraction weighted by
limb darkening), `theta` the incidence angle on the local slope taken from the DEM,
`A(theta)` the incidence-dependent albedo of Hayne2017, and `f_eclipse` the fraction of
the Sun not covered by the Earth.

Diffuse input (review change 3, on by default: `terrain_radiation`, `earthshine`):
`Q_diffuse = (1 - A_d) E_ts + eps E_ti + (1 - A(theta_E)) E_er + eps E_ei`, with `E_ts` and
`E_ti` the sunlight scattered and the infrared emitted by surrounding terrain (single
bounce), `E_er` and `E_ei` the reflected and infrared Earthshine, `A_d` the albedo law
averaged over isotropic illumination and `theta_E` the incidence angle of the Earth. See
"Terrain radiation and Earthshine" below.

## Numerical method and why

| Choice | Alternatives | Why |
|---|---|---|
| Conservative finite volumes, enthalpy as the stored variable | Node-centred finite differences on temperature (Hayne2017, heat1d) | With `c_p(T)`, advancing enthalpy is what makes stored energy change by exactly the boundary fluxes. *Measured*: relative energy residual 4e-15 over the test year. |
| Fluxes as differences of the Kirchhoff variable `psi(T)` | Average the conductivities of neighbouring cells | `K` is separable, `K_c(z) k(T)`, so `-K_c d(psi)/dz` is the flux exactly. Steady conduction is reproduced to round-off on any grid, however large the temperature contrast. The obsolete heat1d 0.3.2, which averages conductivities, is 2-3 K warmer at depth on its only stable grid; the current 0.4.1 agrees within 1 K (see below). |
| BDF2 in time, two steps per hourly forcing sample | Forward Euler; backward Euler; Crank-Nicolson | *Measured* at the test point against a 64-fold finer step: BDF2 0.019 K in 0.56 s; forward Euler 0.016 K in 1.5 s (27 sub-steps forced by stability); backward Euler 0.20 K; Crank-Nicolson 0.010 K but it overshoots the radiative-equilibrium bound by up to 15 K after an abrupt change, which no physical solution can do. BDF2 is second order and L-stable. |
| Newton iteration with the exact tridiagonal Jacobian | Lagged coefficients | Converges in two or three iterations with the T^4 and T^3 nonlinearities fully implicit. |
| Spin-up by imposing the exact mean profile | Integrate for decades | See next section. |

## Initial condition

The state at the start of the period is the periodic state under the forcing of the
preceding year. Deep layers respond over decades (conduction time across 2 m is about
25 years), so a plain spin-up is slow: heat1d's default of one year leaves its deep
temperature 9 K from the converged value (*measured*: 261.7 K after 1 year, 252.5 K
after 20 and after 80).

Here the deep profile is not integrated but imposed. In a periodic state the time-mean
of `psi` at any depth equals its time-mean at the surface plus the interior heat flow
times the thermal resistance in between; this is exact for the equations above. After
each cycle every cell is shifted to satisfy it. *Measured*: convergence to 1e-3 K in
8-10 cycles; starting from 60 K or from 320 K changes the result by 2e-4 K.

## Assumptions and simplifications

| Assumption | Effect |
|---|---|
| Regolith properties from Diviner fits at latitudes within 70 degrees | Polar use is an extrapolation; flagged in every result. |
| Heat-capacity fit stated valid from below 90 K | The test point reaches 43 K. Flagged. The WoodsRobinson2019 c_p adopted by DSNE (SLS-SPEC-159 Rev I, Table 3.4.6.2.2-1, 10-300 K, stated valid 10-400 K) is available as an option; it is 16 % lower at 43 K but changes the reference scenario by -0.13 K at the minimum, -0.03 K in the mean and at most -0.25 K at depth (nominal run, 2026-11-15). The default stays the Hayne2017 polynomial so that the code comparisons (TH-X1, TH-X1b) remain like-for-like. Low-temperature *conductivity* (DSNE: the Hayne relations "will likely over predict" it below 100 K; review finding F09) is available as an option and every run reports its effect (TH-C2: rim minimum -3.3 K, 1 m mean +7.8 K). |
| Albedo law fitted to level equatorial ground, applied at grazing polar incidence | **Largest model-form uncertainty.** At 88 degrees incidence the law gives an albedo of 0.79. With a constant albedo of 0.12 instead, the test point is 26 K warmer on average and 50 K warmer at its maximum. |
| Slope from the 5 m DEM | **Largest input uncertainty.** With the Sun never more than 1.8 degrees high, a tilt of 1 degree changes the absorbed flux severalfold. Slope over a 20 m baseline instead of 5 m raises the maximum by 7 K. |
| Terrain radiation single bounce; Earthshine from guidebook peak values | Before change 3 the neglected terrain flux was bounded at 0.02 W m-2 at the test point (11 % of the emission at the minimum). With single-bounce terrain radiation the ground at the test point receives 8e-5 W m-2 at its minimum (its 0.4 deg slope hides almost all terrain), and Earthshine 0.006 W m-2 absorbed: +0.01 K and +0.23 K on the minimum (TH-C1). Multiple scattering is not included: inside craters results are lower bounds. Runs warn when the diffuse input at the minimum exceeds 5 % of the emission. |
| One-dimensional, no rocks, no lateral conduction | Valid at scales above the skin depth (centimetres to decimetres). |
| Incidence angle of the disc centre for partial discs | Second-order effect during sunrise and sunset. |
| Lunar eclipses as overlap of two uniform discs, no atmosphere | Contact times within a minute or two of almanac values. Checked against the eclipse of 2026 March 3. |

## What the test point looks like (2026, nominal parameters)

Point (-2822.5, -9062.5) m on the Shackleton rim crest, latitude -89.69, slope 0.4
degrees. Direct sunlight reaches it 77 % of the year at ground level, but the ground is
almost level and the Sun is at most 1.8 degrees high, so it absorbs at most 12 W m-2
and nothing at all for the half of the year when the Sun is below the local horizontal.

| Quantity | Nominal | 5 % / median / 95 % of the ensemble |
|---|---|---|
| Surface minimum | 43 K | 40 / 45 / 50 K |
| Surface maximum | 119 K | 100 / 124 / 142 K |
| Surface mean | 70 K | 60 / 69 / 78 K |
| Mean amplitude of a lunation cycle | 20 K | 13 / 36 / 60 K |
| Mean at 0.3 m | 73 K | 61 / 71 / 80 K |
| Mean at 1 m | 77 K | 63 / 74 / 84 K |
| Fastest heating | 2.4 K/h | 1.5 / 4.3 / 11.4 K/h |

Figure: [../figures/thermal_site_2026.png](../figures/thermal_site_2026.png).

Lunar eclipses in the period (9 hourly samples) change nothing here: they occur when the
Sun is near the Moon's equatorial plane, which at the pole puts it below the local
horizontal, where level ground absorbs nothing anyway. A surface tilted towards the
Sun, or a structure, would feel them.

## Uncertainty

Ensemble of 128 members: eight regolith parameters sampled by Latin hypercube within
stated intervals (`lunarenv.thermal.uq.DEFAULT_RANGES`, each with its basis; two are
assumptions and say so), combined with terrain realizations (all 100 DEM clones; each brings its own horizon and
slope, 0.1 to 3.2 degrees).

Sensitivity ranking (rank correlation with the output):

| Output | First | Second | Third |
|---|---|---|---|
| Surface maximum | slope +0.83 | albedo -0.53 | heat flow +0.26 |
| Surface mean | albedo -0.68 | slope +0.68 | heat flow +0.27 |
| Surface minimum | slope +0.70 | albedo -0.52 | heat flow +0.29 |
| Cycle amplitude | slope +0.90 | albedo -0.33 | heat flow +0.23 |

Only the first two columns are robust. Each parameter sample is paired with a randomly
chosen terrain clone, and slope dominates, so a correlation of about 0.25 in the third
column is within what chance pairing produces with 128 members (it named a different
parameter when the ensemble was run on 28 clones).

Quantiles are stable with ensemble size: the 5-95 % range of the surface maximum is
100-142 K with 128 members, 100-141 K with 64, 106-142 K with 32. The median of the
ensemble (124 K) is above the nominal value (119 K) because the nominal DEM happens to
be flatter at this pixel (0.39 degrees) than most of its error clones.

Not sampled: the albedo law itself (bounded instead, see below), the radiative
parameter `chi` and the heat-capacity polynomial (no spread available), and relief below
the DEM resolution.

## Bounds: albedo envelope x ensemble (review change 2)

The ensemble above is the spread *within* one albedo law. At grazing polar incidence the
law itself is the larger unknown (review finding F03): published highland laws give 0.24
to 0.86 at the test point's 88.2 degree incidence. Results are therefore reported as
bounds (`albedo_envelope: true`, the default):

| Case | Law | Source | A at 88.2 deg |
|---|---|---|---|
| cold | Hayne form, A0 0.09, a 0.10, b 0.02 | Bürger2024 Eq. 15, section 4.2.2 (Diviner fit to 80 deg latitude) | 0.860 |
| nominal | scenario properties (Hayne2017 Table A1 by default) | Hayne2017 | 0.784 |
| hot | 1.3299e-5 i^2 - 1.0248e-4 i + 0.14382 (i in degrees) | Foote2020 Eq. 9, Apollo 16 highland soil, laboratory, measured to 75 deg | 0.238 |

The review attributes these two bounds to the DSNE Rev J highland envelope, which could
not be obtained; the laws are taken from the papers (`lunarenv.thermal.albedo_laws`).
Each ensemble member is run under every law (the law replaces all albedo parameters,
the sampled A0 included). The reported bounds are the 5 % quantile of the coldest case
and the 95 % quantile of the hottest. At the test point (2026): surface minimum 41-60 K,
maximum 103-200 K, mean 61-109 K, against 43.3 / 119.2 / 70.3 K for the nominal run.

Validation (TH-V3c, pre-registered): Diviner summer temperatures on the 720 m block lie
inside the envelope on average over the sunlit bins, but at its hot edge (hot law
-0.4 K, cold law -34.2 K), and above the hot law in 5 of the 8 sunlit bins and in 59
of 64 summer bins. The hot side of the envelope is therefore the realistic one at this
site. Terrain radiation and Earthshine, added in TH-V3d, warm the block by only 0.7 K
(sunlit bins) and 0.7 K (all summer bins), so they do not explain why Diviner is warmer
still; sub-DEM roughness, the Diviner footprint and observation timing remain. The constant albedo 0.12,
outside the envelope, is kept as a separate alternative.

**Slope grid** (`slope_grid`). The DEM slope at a point is uncertain by about 2 degrees
(Barker2021). The grid runs the nominal terrain horizon with ground slopes and aspects
under each law: at the test point, 0-10 degrees and 8 aspects give a surface maximum of
102 to 253 K.

**Asset surfaces** (`oriented_surfaces`). For plane surfaces of a structure at any tilt
up to vertical: direct solar irradiance (visible disc above the terrain horizon, Earth
occultation, incidence on the plane), time with direct sun, the longest period without
it, and the sky and ground view factors from the horizon profile (`sky + ground = 1`;
for a level plane `1 - mean(sin^2 e_h)`). Their temperature is not computed: it depends
on the asset and belongs in its thermal model (HLS-UG-001 section 9). At the test point
a vertical wall receives up to 1400 W m-2, level ground 45 W m-2. Terrain-scattered
sunlight, terrain infrared and Earthshine on these surfaces are computed with the
scenario's albedo law for the terrain (`oriented_diffuse_irradiance`; see below).

## Terrain radiation and Earthshine (review change 3)

Method (`lunarenv.thermal.terrain_radiation`; pre-registered in
docs/review/preregistration-th-a3-2026-11-16.md, with its amendment):

1. **Visible terrain.** Rays from the receiver (360 by default) over the terrain stack,
   with the horizon tracer's geometry and sample plan. A sample is visible when it rises
   above every closer sample on its ray; it fills the elevation interval from that
   previous maximum to its own elevation.
2. **Patches.** Visible samples are grouped in 36 azimuth sectors and distance bands
   growing by sqrt(2) from 25 m (`terrain_patch_refinement` multiplies both). Each patch
   is represented by 3 points, stratified by solid angle along distance; each point gets
   its own horizon, DEM normal, direct sunlight (same disc model, Earth-eclipse factor of
   the site) and 1-D temperature under the same albedo law and properties. The patch
   exitance is the mean over its points: scattered `A(theta) E cos(theta)` and emitted
   `eps sigma T^4`, both Lambertian.
3. **Receiver.** Each sample sends its patch's exitance times its cosine-weighted solid
   angle over pi, integrated exactly in elevation and clipped at the receiver plane, so
   any orientation works, including vertical walls. The weights of all terrain plus the
   sky add up to 1 (checked against the horizon view factors to 1e-3).
4. **Cost control.** Patches that cannot send more than 1e-4 W m-2 to any receiver of the
   run are not computed; the bound of what they could send is reported. Cold and hot
   albedo cases are computed for the ground receivers only (site normal, slope grid);
   asset surfaces get the scenario law.
5. **Earthshine** (HLS-UG-001 section 11.1): normal irradiance
   `(0.100 gamma + 0.056) (384400 km / d)^2` W m-2 (reflected, infrared), `gamma` the
   illuminated fraction of the Earth's disc from the Sun-Earth-site phase angle, times the
   visible fraction of the Earth's disc above the terrain horizon and the cosine of
   incidence. That the guidebook's peak values refer to the mean distance is an assumption.

Not included: multiple scattering between terrain elements, roughness below the DEM
resolution, non-Lambertian scattering and emission, shadowing of the ground by the asset,
and sunlight variation within a patch beyond its 3 points. Ensemble members and the
constant-albedo alternative reuse the nominal terrain's diffuse irradiance (absorbed with
their own albedo and emissivity). In the receiver's spin-up year, the terrain's own
temperatures come from a periodic spin-up of each patch point, then a run over both years.

Verification: TH-A3 (analytic spherical bowl, single bounce, Ingersoll1992) and TH-C1
(terrain irradiance at the rim minimum against the earlier bound), VALIDATION.md.

Results (2026-11-16): TH-A3 PASS, within 0.4 % at the pre-registered fine layout and
3.0 % at the default layout. TH-C1 PASS: at the rim minimum the ground receives at most
8.5e-5 W m-2 from terrain (bound 0.0217 W m-2, threshold 0.05 W m-2); Earthshine 0.006
W m-2 absorbed, raising the minimum by 0.23 K. Doubling sectors and bands changes the
terrain value from 8.0e-5 to 4.5e-5 W m-2 (ground receiver alone): at a point whose plane
hides almost all terrain, the few patches that graze it are not converged in relative
terms; the absolute change, 3.5e-5 W m-2, is far below anything that moves the
temperature. For receivers that see much terrain (walls, crater floors) the bowl case is
the relevant check.

## Numerical error at the default settings

Largest instantaneous difference of the surface temperature against the finest setting
tried, at the test point (full study: [../convergence/thermal.md](../convergence/thermal.md)):
time integration 0.02 K, hourly forcing 0.12 K, grid 0.05 K, spin-up 0.0002 K; about
0.2 K in total. Errors of extremes and means are below 0.01 K. This is three orders of
magnitude below the input uncertainty.

## Performance

One point, nine simulated years (eight spin-up cycles and the period) at hourly forcing,
43 cells: 0.3 CPU-seconds with one step per sample, 0.6 with the default two. Cost is
linear in steps and cells, about 90 ns per cell and step. The 128-member ensemble took
42 CPU-seconds.

Terrain radiation adds 1-2 CPU-seconds per patch point (spin-up, then two years). At the
test point the ground receivers (site slope and slope grid) needed 131 patches (393
points) per albedo case and the six asset surfaces 317 patches (934 points) under the
scenario law: about 530 s of the 932 s wall-clock reference run on 7 worker processes.

## Verification

| Case | Reference | Result |
|---|---|---|
| TH-A1 periodic surface temperature, constant properties | damped wave (CarslawJaeger1959) | second order in space; in time order 1 (backward Euler) and 2 (Crank-Nicolson, BDF2), measured |
| TH-A2 steady state with interior heat flow | closed form, constant and temperature-dependent conductivity | exact to 1e-9; second order with depth-dependent density |
| TH-A2 vanishing conductivity | radiative equilibrium | 1e-4 K |
| TH-P1 energy conservation | scheme identity | relative residual below 1e-9 for all four integrators |
| TH-P1 bounds and monotonicity | radiative-equilibrium bound; more flux never cools | hold for forward Euler, backward Euler, BDF2; Crank-Nicolson violates the bound (kept as a test) |
| Periodic state | mean-`psi` identity; independence of initial state | 2e-5 relative; 0.01 K |
| Equatorial cycle (pre-check of TH-V1) | Hayne2017 Table A2: 385 / 101 / 95 K at noon / midnight / minimum | 385.2 / 100.3 / 93.7 K |
| Same, constant albedo 0.12 (TH-V1b) | same | 385.2 / 101.8 / 95.1 K: also passes, so TH-V1 does not discriminate |
| Lunar eclipse of 2026-03-03 | NASA eclipse figure (Espenak) | site contacts inside the published windows |

## Findings to resolve in Phase 4

1. **Reference code** (validation cases TH-X1, TH-X1b, TH-X1c). Resolved. Against
   heat1d 0.4.1, the current GitHub code, with heat1d's own parameters:
   - default heat1d grid (TH-X1b, pre-registered): FAIL; lunarenv's night minima are
     1.5-1.6 K and its means 1.1-1.6 K warmer, maxima within 0.6 K, depth within 1 K;
   - heat1d grid m=n=40, 125 layers (TH-X1c, pre-registered after TH-X1b): PASS; all
     surface values within 0.8 K and the means at heat1d's deepest node within 0.8 K at
     0, 20 and 26 deg; heat1d itself changes by at most 0.06 K from m=n=20 to 40.
   The TH-X1b failure is therefore discretisation error of heat1d's default grid (about
   20 layers, thickness growth 1.25). The 2-3 K deep gap against the obsolete PyPI
   release 0.3.2 (TH-X1), previously attributed to heat1d's averaging of conductivities,
   does not appear with 0.4.1 and is attributed to that release.
2. **Apollo temperatures.** With the published parameter set the model gives a mean
   surface temperature of 205.6 K at 20 N and 203.0 K at 26 N, and 247.0 / 241.7 K at
   depth. The Apollo values quoted by Hayne2017 are 216 +/- 5 and 211 +/- 5 K at the
   surface, 256 +/- 5 and 252 +/- 5 K at depth: the model is 8-10 K low. With a
   constant albedo of 0.12 it gives 213.4 / 255.2 K and 211.2 / 250.7 K, inside the
   quoted uncertainties, and matches the curves of Hayne2017's Figure A1. An earlier
   version of this note concluded that the figure corresponds to a constant albedo. The
   figure caption says "All model parameters are from Table A1", which includes the
   incidence-dependent law, so that conclusion is not supported; the discrepancy is
   unresolved. Which albedo is right at grazing incidence decides the polar result
   (26 K on the mean at the test point).
   **Like for like (TH-V2c, pre-registered 2026-11-16): fail.** Read from the primary
   source (Keihm & Langseth 1973), the Apollo 17 surface mean of 216 +/- 5 K is a
   brightness temperature (reduced with e_m = 1, mean solar constant, albedo 0.08), not
   a physical temperature. Compared with the model's mean of eps^(1/4) T_s, with the
   site albedo A0 = 0.08 (incidence dependence of Hayne2017 kept) and the measured heat
   flow 0.016 W m-2, the model gives 207.0 K (-9.0 K, outside 5 K) at the surface and
   252.7 K at 1.30 m (-3.3 K, inside). Changing the quantity and the site albedo
   therefore does not close the gap: with the published incidence law the model stays
   about 9 K too cold at an equatorial mare site measured day and night. A constant
   albedo 0.08 gives 214.4 / 260.6 K (informational). Apollo 15 (informational, its mean
   partly modelled and lowered by topographic shadowing): 204.4 K against the final
   report's 207 +/- 7 K. TH-V2 stays FAIL in the ledger. This favours a weaker
   incidence dependence at low latitudes; it does not say which law holds at grazing
   polar incidence (TH-V3c).
3. **TH-V1 does not discriminate between albedo laws** (TH-V1b, pre-registered
   2026-11-13). With a constant albedo of 0.12 the equatorial case gives
   385.2 / 101.8 / 95.1 K and passes as well as the incidence-dependent law does
   (385.2 / 100.3 / 93.7 K). TH-V1 is therefore no evidence for the albedo law used at
   the pole.
4. **First polar comparison (TH-V3, TH-V3b, pre-registered 2026-11-15).** Diviner
   summer Polar Cumulative Product, 3 x 3 cells of 240 m around the test point,
   bolometric temperature of the 720 m block over the summer subsolar-longitude bins in
   which the model block sees at least half of the solar disc on average. Published
   albedo law: bias -23.8 K (FAIL at 10 K). Constant albedo 0.12: +5.7 K (PASS). Caveats
   that limit what this shows: only 8 of the 64 complete bins qualify as sunlit; per-cell
   biases range from -70 to +32 K, so the spatial pattern inside the block is not
   reproduced (Diviner's footprint is larger than a cell, and the model has no terrain
   radiation or sub-DEM roughness); the hottest Diviner bin (254 K) is 43-67 K above
   both laws; over all summer bins, sunlit or not, both laws are too cold (-39 / -13 K),
   as expected without terrain radiation. Pixel sampling (20 m vs 10 m) changes the
   centre cell by under 0.1 K; the model's summer-to-summer spread is about 1 K.
   Ground-level horizons computed on the site DEM block-averaged to 10 and 20 m instead
   of 5 m (surface normal kept at 5 m; `benchmarks/th_v3_horizon_sensitivity.py`,
   informational, after the result) raise the mean summer visible fraction of the block
   from 0.40 to 0.41 and 0.43 and change the block bias on the same 8 bins by at most
   0.5 K (published law -23.8 / -23.6 / -23.5 K; constant +5.7 / +6.0 / +6.2 K) and
   per-cell biases by at most 1.6 K. With 20 m horizons 10 bins qualify as sunlit and
   the biases become -24.2 / +4.5 K. The TH-V3 / TH-V3b outcome does not depend on the
   horizon resolution in this range (review concern on 0 m horizons, cf. ILL-V5).
5. **Low-temperature conductivity (F09; TH-A2b, TH-C2, pre-registered 2026-11-16).**
   Option `conductivity_model: woods_robinson2019`: `K = K_c(rho) [g(T) + chi (T/350)^3]`
   with `g = 1` above 150 K, `g = k_am(T) / k_am(150 K)` between 20 and 150 K (k_am: the
   WoodsRobinson2019 fit to amorphous SiO2, section 4.3) and linear to zero below 20 K,
   where the source calls its model nonphysical. `g` is 0.31 at 43 K and 0.69 at 95 K.
   The sign of the coefficient `G_am` is taken as negative: as printed (positive) the fit
   gives 6.4 W m-1 K-1 at 300 K, against 1-2 in the source's own Figure 6c; with the
   negative sign it gives 1.38. The full WoodsRobinson2019 model was not adopted: its
   constants come from basalt powder, it would replace the Diviner-fitted conductivity at
   all temperatures, and its density dependence is not separable. Martinez & Siegler
   (2021) was not accessible. Results: the steady analytic case with the option passes
   to 1e-15 (TH-A2b). At the rim the option lowers the minimum by 3.3 K and the surface
   mean by 1.6 K, and raises the mean at 0.3 m by 5.2 K and at 1 m by 7.8 K (TH-C2, terrain
   radiation off). At the equator it moves the night values by 0.6 K (93.1 K against
   95 K, still within TH-V1's 2 K); at Apollo 17 (TH-V2c setting) 206.7 / 254.1 K, so it
   does not close TH-V2c's surface gap. Reading fixed beforehand: the 1 m change exceeds
   5 K, so depth results are reported with both laws: every run computes the nominal case
   with the other law (`conductivity_alternative`, on by default) and warns when a mean at
   depth moves by more than 1 K. The option is a bound on the conductivity form, not a
   refitted model: Hayne2017's K_s, K_d and chi were fitted with a constant contact term,
   and no measurement of lunar regolith conductivity below 100 K is known to this project.
6. **Polar comparison with terrain radiation and Earthshine (TH-V3d, pre-registered
   2026-11-16).** TH-V3's forward model plus the change caused by single-bounce terrain
   radiation (engine method, default patch settings, patch temperatures with the same law)
   and Earthshine, computed at 4 receivers per 240 m cell (36 in all; the full 1296
   pixels would take over 100 h) and added per cell. Hayne2017 law: bias -23.2 K over the
   8 sunlit bins, still FAIL (+0.66 K from TH-V3); all summer bins -37.9 K (+0.75 K). Hot
   law: +0.1 K / -17.4 K (+0.49 / +0.56 K from TH-V3c). Absorbed terrain radiation over
   the summer samples is 0.77 W m-2 on average over the receivers and 4.0 W m-2 at the
   largest, Earthshine 0.002 W m-2; cell changes range from 0.1 to 3.7 K. Reading fixed
   beforehand: the change is far below 5 K, consistent with TH-C1, and the Hayne2017
   law's cold bias at this site is not due to the missing diffuse terms. Of TH-V3's
   -39 K over all summer bins they explain under 1 K; the remainder points to sub-DEM
   roughness, the 0 m horizon on a 5 m DEM, the Diviner footprint and observation timing,
   or the low-temperature properties (the low-temperature conductivity option changes the
   block by 0.05 K). Limits: 4 of 144 pixels per cell estimate the cell's change, not its
   temperature; multiple scattering is absent. Run time 3.0 h.
