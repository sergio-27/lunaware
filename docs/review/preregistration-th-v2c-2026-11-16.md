# Pre-registration: TH-V2c, Apollo mean temperatures compared like for like

Written before the case is implemented or run.

Origin: change 3 (d) and finding F05 of `docs/review/scientific-review-2026-11-13.md`:
TH-V2 compares a physical model temperature with a brightness temperature, with the
generic highland albedo at dark mare sites. The review asks for "a like-for-like case,
registered before it is run: site albedo from the source, brightness-temperature output
with eps = 1, the source's solar distance. Tolerance unchanged (+/-5 K)."

What is already known (not blind): TH-V2 and TH-V2b results; the review's own model
values for related settings (A0 = 0.06 at 1 au: 210.9 / 208.2 K; constant 0.08:
215.5 / 213.4 K). The configuration below has not been run.

## Sources (read 2026-11-16)

- Keihm & Langseth (1973), Proc. Lunar Sci. Conf. 4th, 2503-2513 (ADS scan,
  1973LPSC....4.2503K, pages read as images). Thermocouples about 15 cm above the
  surface; reduction equation (1)-(2) with "e_m = infrared emissivity of the lunar
  surface = 1.0", "S = the mean solar constant = 0.1352 W/cm2", "A = the lunar albedo =
  0.08". "The choice of e_m = 1.0 implies that the surface temperatures deduced
  correspond to blackbody brightness temperatures." "From the data shown in Fig. 2, a mean
  surface temperature of 216 (+/-5) K was calculated for the Taurus Littrow site", against
  "the 256 K temperature at 130 cm measured by the heat flow probes". The lunation of
  Fig. 2 is not stated.
- Langseth, Keihm & Chute, Lunar heat-flow experiment, final report (NTRS 19780010035),
  Table 1 and section V-A: Hadley Rille (Apollo 15) mean surface temperature 207 K
  ("207 +/- 7 K"), 253.0 K at 1.38 m, heat flow 21 mW m-2; Taurus Littrow 216 K,
  256.7 K at 2.34 m, heat flow 16 mW m-2. At Apollo 15 the daytime thermocouple data
  could not be interpreted to better than 20 K, and "theoretical models taking into
  account latitude and previous albedo estimates were used to supplement the nighttime
  data"; the lower Apollo 15 mean "is due primarily to the early morning shadowing
  effect" of the local topography.
- Not read: Keihm et al. (1973), Earth Planet. Sci. Lett. 19, 337-351 (Apollo 15), the
  source of the 211 K / 252 K at 0.83 m that Hayne2017 Table A2 quotes.

## Case

- Apollo 17 only is checked: its mean comes from measured day and night temperatures on
  level ground. Apollo 15 is informational: its mean is partly a model, is lowered by
  topographic shadowing that a 1-D level model cannot have, and its reduction source was
  not read.
- Model: Hayne2017 parameter set with the site values the sources give: normal albedo
  A0 = 0.08 (the source's lunar albedo; the incidence dependence a, b of Hayne2017 kept)
  and interior heat flow 0.016 W m-2 (final report, Taurus Littrow). Latitude 20 deg,
  declination 0, Sun at 1 au (the source reduces with the mean solar constant; the
  measured lunation is not stated). Same numerical settings as TH-V2.
- Output compared with the surface value: the diurnal mean of the brightness temperature
  eps^(1/4) T_s with the model emissivity (0.95), i.e. what a sensor reduced with
  e_m = 1 would report. At 1.30 m: the diurnal mean of the physical temperature.
- Checks: surface brightness-temperature mean against 216 K, mean at 1.30 m against
  256 K, +/-5 K each (unchanged). PASS if both pass, FAIL otherwise.
- Informational: Apollo 17 with a constant albedo 0.08; Apollo 17 with the Sun at
  1.017 au (as TH-V2); Apollo 15 (26 deg) with the same treatment (A0 = 0.08, heat flow
  0.021 W m-2) against 211 K / 252 K at 0.83 m (Hayne2017's quote) and 207 K / 253 K at
  1.38 m (final report).

## Reading

- PASS: TH-V2's failure came from comparing unlike quantities with a highland albedo; the
  equatorial mare data do not contradict the incidence-dependent law. TH-V2 stays FAIL in
  the ledger.
- FAIL: even like for like, the model with the published law is too cold at Apollo 17;
  the albedo question at mare sites stays open.
- Neither outcome says which law holds at grazing polar incidence (TH-V3c does that, at
  its own resolution).
