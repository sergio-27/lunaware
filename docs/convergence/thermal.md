# Thermal model: convergence, sensitivity and numerical error

Script: `benchmarks/thermal_site.py`; raw output: [thermal.json](thermal.json). Run
2026-11-12. Test case: Shackleton rim point (-2822.5, -9062.5) m, year 2026 with 2025
as spin-up cycle, hourly forcing unless stated, nominal Hayne2017 parameters. Metric:
largest instantaneous difference of the surface temperature over the year against the
finest setting of each study; differences of the yearly extremes and mean are given
where they are not negligible.

## Time integration (hourly forcing; reference: BDF2 with 64 steps per hour)

| Scheme | Steps per hour | Max difference (K) | RMS (K) | CPU (s) |
|---|---|---|---|---|
| Forward Euler | 27 (stability limit) | 0.016 | 0.0008 | 1.53 |
| Backward Euler | 1 / 2 / 4 | 0.39 / 0.20 / 0.10 | 0.019 / 0.010 / 0.005 | 0.30 / 0.61 / 1.25 |
| Crank-Nicolson | 1 / 2 / 4 | 0.067 / 0.010 / 0.003 | 0.0022 / 0.0004 / 0.0001 | 0.31 / 0.58 / 1.17 |
| BDF2 | 1 / 2 / 4 | 0.092 / 0.019 / 0.0025 | 0.0024 / 0.0005 / 0.0001 | 0.30 / 0.56 / 1.12 |

Backward Euler halves its error per halving of the step (order 1); Crank-Nicolson and
BDF2 reduce it 4-7 fold (order 2). On this smooth forcing Crank-Nicolson is slightly
more accurate than BDF2, but it is not L-stable: in the unit test with an abrupt start
it exceeds the radiative-equilibrium temperature by 15 K in the first cell. Sunrise
over a ridge is abrupt, so BDF2 is the default, with two steps per hourly sample
(0.019 K, three times cheaper than forward Euler at the same accuracy).

## Forcing sampling (solver step fixed at 5 min; reference: 5 min illumination)

| Illumination sampled every | Max difference (K) | RMS (K) |
|---|---|---|
| 15 min | 0.007 | 0.0002 |
| 30 min | 0.027 | 0.0009 |
| 60 min | 0.122 | 0.0039 |

Hourly forcing is the default; the error is confined to the minutes around abrupt
sunrises and sunsets.

## Grid (reference: first cell 0.25 mm, growth 1.05, depth 4 m, 4 steps per hour)

| Parameter | Value | Cells | Max difference (K) | Difference of yearly mean (K) |
|---|---|---|---|---|
| First cell | 8 / 4 / **2** / 1 / 0.5 mm | 31 / 37 / 43 / 49 / 55 | 0.63 / 0.18 / **0.051** / 0.018 / 0.008 | 0.009 / 0.002 / 0.000 / 0.001 / 0.001 |
| Growth ratio | 1.3 / 1.2 / **1.12** / 1.06 | 22 / 30 / 43 / 71 | 0.080 / 0.062 / **0.051** / 0.046 | 0.008 / 0.003 / 0.000 / 0.001 |
| Depth | 0.3 / 0.5 / 1 / **2** / 4 m | 26 / 31 / 37 / 43 / 49 | 4.4 / 0.49 / 0.061 / **0.051** / 0.051 | 0.25 / 0.08 / 0.006 / 0.000 / 0.000 |

Defaults in bold. The error falls about 3-fold per halving of the first cell (it is
dominated by the steep gradient under the surface, where the temperature is not yet
smooth on the scale of a cell). The growth ratio hardly matters. The domain must reach
1 m; 2 m is used so that the yearly wave is damped with margin.

## Spin-up

Starting the spin-up from 60 K and from 320 K: 8 and 10 cycles to reach 1e-3 K; the two
resulting states differ by 1.9e-4 K at most, at any depth, and the surface temperature
in the period by 2.8e-5 K.

## Model-form sensitivities (largest difference of the surface temperature; yearly maximum; yearly mean)

| Change | Max difference | Yearly maximum | Yearly mean |
|---|---|---|---|
| Constant albedo 0.12 instead of the incidence-dependent law | 49.9 K | +49.9 K | +26.3 K |
| Surface slope over a 20 m baseline instead of 5 m | 19.6 K | +7.2 K | +0.9 K |
| Illumination computed 2 m above the surface instead of at it | 8.3 K | 0.0 K | 0.0 K |
| No lunar eclipses | 0.0 K | 0.0 K | 0.0 K |

## Numerical error attached to results

Surface temperature, largest instantaneous error: 0.02 K (time integration) + 0.12 K
(hourly forcing) + 0.05 K (grid) + 0.0002 K (spin-up), about 0.2 K. Yearly extremes and
means: below 0.01 K. Compare with the 5-95 % ensemble range of the yearly maximum,
100-142 K (128 parameter samples on 100 terrain clones).
