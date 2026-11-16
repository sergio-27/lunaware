# Pre-registration: low-temperature conductivity option (F09), TH-A2b and TH-C2

Written before the option is implemented and before anything below has been run.

Origin: finding F09 of `docs/review/scientific-review-2026-11-13.md`: "Low-temperature
conductivity is not represented and matters at depth." Recommended action: "Add a low-T
conductivity option (Woods-Robinson 2019 / Martinez & Siegler 2021) and report the 0.3 m
and 1 m results as unbounded until it exists." The heat-capacity half of F09 (Woods-Robinson
c_p) was done on 2026-11-15.

What is already known (not blind): the reviewer's crude test (conductivity times 0.5 / 0.2
below 80 K) moved the rim minimum by -2.7 / -5.6 K and the 1 m mean by +6.3 / +18.1 K, in
the reviewer's reduced model. Nothing below has been run with lunarenv.

## Sources (read 2026-11-16)

- Woods-Robinson, Siegler & Paige (2019), JGR Planets 124, 1989-2011,
  doi:10.1029/2019JE005955 (author copy https://luna1.diviner.ucla.edu/~dap/pubs/096.pdf,
  SHA-256 df61d7e12f922a3a87a170d5d1a3e3048d7e3238a8f96c5e6b91998ba922350d), section 4.3
  read in text and as page images. Conductivity model, eq. (30)-(32):
  k_eff = A T^3 + (B + C rho^2)(1 - p) k_am(T), with k_am(T) the conductivity of amorphous
  SiO2 fitted as A_am + B_am T^-4 + C_am T^-3 + D_am T^-2 + E_am T^-1 + F_am T + G_am T^2
  + H_am T^3 + I_am T^4. Statements used: "In the measured range 150-400 K, the differences
  between the old and new models are negligible"; "below 150 K the old and new fits begin
  to deviate significantly"; "Results at and below ~20 K in Figure 11a are nonphysical".
  The constants A, B, C were fitted to particulate basalt (Fountain & West 1970), not to
  Diviner data.
- Martinez & Siegler (2021), JGR Planets 126, doi:10.1029/2021JE006829: not read (closed
  access; no open copy found by Semantic Scholar or OpenAlex). Not used.

### A sign in the source's coefficient table

As printed, the table gives G_am = +2.80491e-5. With that sign k_am(300 K) = 6.43 W m-1 K-1,
k_am(150 K) = 2.21. The source's own Figure 6c (amorphous SiO2, Touloukian & Ho 1970) shows
about 0.1 W m-1 K-1 near 10 K and between 1 and 2 W m-1 K-1 at 300 K. With G_am =
-2.80491e-5 the fit gives 0.124 (10 K), 0.682 (100 K), 0.946 (150 K), 1.380 (300 K), which
matches the figure. The negative sign is therefore used; this is an inference by the
author of this note, not a statement of the source, and it is recorded in the code and the
references. Coefficients used:

| A_am | B_am | C_am | D_am | E_am | F_am | G_am | H_am | I_am |
|---|---|---|---|---|---|---|---|---|
| -2.03297e-1 | -11.472 | 22.5793 | -14.3084 | 3.41742 | 0.01101 | **-2.80491e-5** | 3.35837e-8 | -1.40021e-11 |

## Method (fixed now)

The Woods-Robinson model as a whole is not adopted: its constants come from basalt powder,
it would replace the Diviner-fitted Hayne2017 conductivity at all temperatures, and its
radiative and contact terms depend on density differently, which the solver's Kirchhoff
formulation does not allow. The option keeps Hayne2017 above 150 K and imposes the
Woods-Robinson temperature shape on the contact (phonon) term below it:

    K(T, z) = K_c(rho(z)) [ g(T) + chi (T / 350)^3 ]

    g(T) = 1                                   T >= 150 K
    g(T) = k_am(T) / k_am(150 K)               20 K <= T < 150 K
    g(T) = g(20 K) T / 20 K                    T < 20 K   (linear to zero; the source calls
                                                           its model nonphysical there)

- 150 K is where the source says the low-temperature behaviour starts to matter; above it
  the option is identical to the default.
- g(20 K) = 0.152, g(43 K) = 0.31 (about the rim minimum), g(95 K) = 0.69 (equatorial
  night).
- K stays separable in depth and temperature, so the Kirchhoff variable is the integral
  of g(T) + chi (T/350)^3 (analytic in each branch) and the solver is otherwise unchanged.
- Name: `RegolithProperties.conductivity_model = "woods_robinson2019"` (default
  `"hayne2017"`, unchanged results). Independent of `heat_capacity_model`.
- This is a sensitivity option, not a recalibrated model: Hayne2017's K_s, K_d and chi
  were fitted with a temperature-independent contact term, and are not refitted here.
- Runs report, for the scenario's nominal case, the same quantities with the other
  conductivity law (surface minimum and mean; mean at 0.3 m and 1 m), so that depth results
  are no longer unbounded with respect to this choice.

## TH-A2b: steady conduction with the low-temperature option (analytic)

- As TH-A2 (uniform density, steady state reached by implicit steps of 1e9 s, Kirchhoff
  closed form: the potential is linear in depth with slope Q), with
  `conductivity_model = "woods_robinson2019"`, K_c = 1e-5 W m-1 K-1, Q = 0.005 W m-2, no
  absorbed flux, grid 4 mm growing by 1.3 to 0.4 m. Surface (0.005 / eps sigma)^(1/4) =
  17.4 K, so the profile crosses all three branches of g(T).
- Checks: surface temperature within 1e-6 K; relative error of the potential profile at
  most 1e-9 (the TH-A2 tolerances). PASS if both pass.

## TH-C2: effect of the low-temperature conductivity (sensitivity, no pass/fail)

Status INFO by construction; numbers recorded, reading below.

- Rim: `scenarios/shackleton_rim_2026.yaml`, no ensemble, no slope grid, terrain radiation
  and Earthshine off (TH-C1: they move the rim minimum by less than 0.3 K), nominal albedo
  law; Hayne2017 against the option: change of the surface minimum and mean, and of the
  mean at 0.3 m and 1 m.
- Equator (TH-V1 setting): noon, midnight and 6 am temperatures with the option, against
  385 / 101 / 95 K, with TH-V1's 2 K shown for reference.
- Apollo 17 (TH-V2c setting): brightness mean and 1.30 m mean with the option, against
  216 / 256 K.

## Reading

- TH-A2b PASS: the option is implemented consistently with the solver's exact steady
  state. FAIL: the option is wrong and is not used.
- TH-C2, rim: the reviewer's crude test predicts a colder minimum (a few K) and a warmer
  1 m mean (several K or more). A change of the 1 m mean of more than 5 K means depth
  outputs must be reported with both laws; less than 1 K means the review's concern was
  overstated for this site.
- TH-C2, equator: if the option moves the night values by more than 2 K, it is not
  compatible with the equatorial fit of Hayne2017 without refitting, and is to be read as
  a bound on the conductivity form, not as a better model.
- Not tested: the true low-temperature conductivity of lunar regolith (no measurement of
  lunar samples below 100 K is known to this project).
