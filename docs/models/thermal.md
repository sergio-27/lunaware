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
input, the ground temperature it sits on and radiates to, and how both vary.

## Physics

One-dimensional heat conduction, `rho(z) c_p(T) dT/dt = d/dz (K(T, rho) dT/dz)`, with

- density rising exponentially with depth, `rho_s = 1100` to `rho_d = 1800 kg m-3`,
  scale `H = 0.06 m`;
- conductivity `K = K_c(rho) [1 + chi (T/350)^3]`, the second term being radiative
  transfer between grains;
- heat capacity as a quartic polynomial in temperature;
- surface: absorbed sunlight plus conduction from below balances thermal emission,
  `Q_abs + K dT/dz = eps sigma T_s^4`;
- base: interior heat flow `Q = 0.018 W m-2`.

All laws and default values are those of Hayne2017 (Appendix A, Table A1).

Absorbed sunlight: `Q_abs = (1 - A(theta)) E cos(theta) f_eclipse`, with `E` the direct
irradiance from the illumination model (distance, visible disc fraction weighted by
limb darkening), `theta` the incidence angle on the local slope taken from the DEM,
`A(theta)` the incidence-dependent albedo of Hayne2017, and `f_eclipse` the fraction of
the Sun not covered by the Earth.

## Numerical method and why

| Choice | Alternatives | Why |
|---|---|---|
| Conservative finite volumes, enthalpy as the stored variable | Node-centred finite differences on temperature (Hayne2017, heat1d) | With `c_p(T)`, advancing enthalpy is what makes stored energy change by exactly the boundary fluxes. *Measured*: relative energy residual 4e-15 over the test year. |
| Fluxes as differences of the Kirchhoff variable `psi(T)` | Average the conductivities of neighbouring cells | `K` is separable, `K_c(z) k(T)`, so `-K_c d(psi)/dz` is the flux exactly. Steady conduction is reproduced to round-off on any grid, however large the temperature contrast. The reference code heat1d, which averages conductivities, is 2-3 K warmer at depth on its only stable grid (see below). |
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
| Heat-capacity fit stated valid from below 90 K | The test point reaches 43 K. Flagged. A low-temperature property model (WoodsRobinson2019) is the natural extension. |
| Albedo law fitted to level equatorial ground, applied at grazing polar incidence | **Largest model-form uncertainty.** At 88 degrees incidence the law gives an albedo of 0.79. With a constant albedo of 0.12 instead, the test point is 26 K warmer on average and 50 K warmer at its maximum. |
| Slope from the 5 m DEM | **Largest input uncertainty.** With the Sun never more than 1.8 degrees high, a tilt of 1 degree changes the absorbed flux severalfold. Slope over a 20 m baseline instead of 5 m raises the maximum by 7 K. |
| No light scattered or heat radiated by surrounding terrain | Upper bound from the horizon profile at the test point: 0.02 W m-2 against a mean absorbed flux of 2.5 W m-2, i.e. below 1 %. Inside craters it would dominate; the bound is reported with each result. In shadow the computed temperatures are lower bounds. |
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

Not sampled: the albedo law itself (reported above as a model-form sensitivity), the
radiative parameter `chi` and the heat-capacity polynomial (no spread available), and
relief below the DEM resolution.

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

1. **Reference code** (validation case TH-X1). With heat1d's own parameters (H = 0.07 m,
   its Sun distance and declination) the two codes agree at the surface within 0.75 K
   at latitudes 0, 20 and 26 (mean, maximum and minimum). At 0.68 m, heat1d's deepest
   node, lunarenv is 1.9 to 2.9 K colder than equilibrated heat1d. lunarenv's deep
   value does not move by more than 0.013 K between grids with first cells from 3.6 mm
   to 0.25 mm, including heat1d's own grid. heat1d returns NaN on every grid finer than
   its default, so its value cannot be refined. The difference is attributed to
   heat1d's averaging of conductivities across the steep near-surface gradient; this is
   a hypothesis consistent with the evidence, not a demonstration.
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
