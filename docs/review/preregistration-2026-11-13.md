# Pre-registration: follow-up cases from the scientific review of 2026-11-13

Written before any of these cases was implemented or run. It fixes, for each case, what is
compared, the tolerance and how each outcome will be read. Tolerances are taken from cases
already in `VALIDATION.md`; none was chosen after seeing a result. Existing cases and their
tolerances are not changed.

Origin: findings F04, F08 and F13 of `docs/review/scientific-review-2026-11-13.md`.

## TH-V1b: equatorial diurnal temperatures with a constant albedo (F04)

- **Question.** Does TH-V1 discriminate between the incidence-dependent albedo law and a
  constant albedo? If both pass, TH-V1 is not evidence for the albedo law used at the pole.
- **Set-up.** Identical to TH-V1 (equator, level ground, Sun at 1 au, declination 0,
  standard parameters) except `albedo_a = albedo_b = 0`, i.e. a constant albedo of 0.12.
- **Reference and tolerance.** Same as TH-V1: Hayne2017 Table A2, noon 385 K, midnight 101 K,
  minimum 95 K, each within 2 K.
- **Reading.**
  - PASS: TH-V1 does not discriminate between the two albedo laws. This is recorded in the
    notes of TH-V1 and in the thermal model documentation.
  - FAIL by more than 2 K on any quantity: TH-V1 does discriminate, and the reviewer's reduced
    model (which predicts 385.2 / 101.7 / 94.9 K) is wrong on this point.

## ILL-V2b: the same map comparison with the metric of the published label (F08)

- **Question.** Does the ~6 pp deficit on well-lit cells in ILL-V2 come from comparing
  different metrics?
- **Evidence for the hypothesis.** The PDS label of `AVGVISIB_85S_060M_201608` (read
  2026-11-13 from pds-geosciences.wustl.edu with TLS verification on) says the values
  "indicate the percentage of timesteps where each pixel was sunlit (by any fraction of the
  solar disc)". The same label describes the image object as "average solar disc
  visibility", so the label is ambiguous; this case decides between the two readings
  empirically.
- **Set-up.** Identical to ILL-V2 (same 100 x 100 cells of 60 m, Site 4 DEM block-averaged
  to 60 m, regional DEM, 720 azimuths, single-chord disc, observers at ground level,
  2010-01-01 plus 18.6 years every 2 h). The only change is the lunarenv statistic: the
  percentage of time steps in which any part of the solar disc is above the horizon, instead
  of the mean visible fraction of the disc.
- **Checks and tolerances.** The three ILL-V2 checks unchanged (bias within 2 pp, RMS within
  6 pp, correlation >= 0.95), plus one discriminating check: the mean difference over cells
  the published map puts above 50 % within 2 pp (the reviewer's stated falsification
  criterion).
- **Reading.**
  - Mean difference on cells above 50 % within +/-2 pp: the metric explains the deficit.
    ILL-V2 stays as recorded; its note is updated to say it compares unlike metrics and
    ILL-V2b is the like-for-like comparison. The 1.5 m observer-height explanation is
    dropped.
  - Mean difference on those cells still near -6 pp: the metric hypothesis is rejected and the
    observer-height or slope-sampling explanation remains open.
  - Anything in between is reported as it is, with no further case added on the strength of it.

## ILL-V5: named points of Gläser et al. (2018) (F13)

- **Source.** Gläser, Oberst, Neumann, Mazarico, Speyerer, Robinson (2018), *Illumination
  conditions at the lunar poles: Implications for future exploration*, PSS 162, 170-178.
  Read in full from the accepted manuscript on NTRS (20170007365); the final journal version
  was not checked. Method: 20 m/pixel LOLA DTM, horizon method of Mazarico2011, 2017-01-01
  to 2037-01-01 sampled hourly. "Average illumination" is the mean visible fraction of the
  disc (their section 2); "solar visibility" is reported separately.
- **Points** (only those inside the cached 5 m site DEMs):

  | Point | Table | Height | Latitude | Longitude E | Avg. illumination | Solar visibility | Max. shadow | DEM |
  |---|---|---|---|---|---|---|---|---|
  | C1-2m | 2 | 2 m | -89.4398 | 222.83 | 88.0 % | 92.1 % | 112 h | Site 1 |
  | S-2m-a | 2 | 2 m | -89.6866 | 197.19 | 85.5 % | 89.4 % | 66 h | Site 4 |
  | S-2m-b | 2 | 2 m | -89.7815 | 204.70 | 82.1 % | 89.8 % | 65 h | Site 4 |
  | C1-0m | 1 | 0 m | -89.4505 | 222.84 | 76.0 % | 81.3 % | 268 h | Site 1 |
  | S-0m | 1 | 0 m | -89.7849 | 203.97 | 81.0 % | 85.5 % | 221 h | Site 4 |

  Connecting Ridge C2 (both tables) lies outside both site DEMs and is not run.
- **Set-up.** Default illumination settings of the engine (1440 azimuths, per-strip disc
  profile, nested 5 m / 80 m / coarsened regional terrain), nominal DEM, Sun from DE440,
  2017-01-01T00:00 to 2037-01-01T00:00 every 1 h (as the source).
- **Checks and tolerances.**
  - Average illumination (mean visible fraction of the disc): 3 pp, the tolerance of ILL-V1.
  - Solar visibility (percentage of hourly steps with any part of the disc visible): 3 pp.
  - Maximum time in shadow: informational, no tolerance (longest-run statistics carry no
    tolerance anywhere in this project).
  - Informational: maximum and minimum average illumination within 40 m (two source pixels)
    of each point, to show the sensitivity to registration.
- **Expected sign of a discrepancy.** The published points are local maxima picked on a 20 m
  DTM, and the source says such locations are 20-40 m wide. Picking a maximum biases the
  published value upwards relative to a point value on an independent DEM, so a lunarenv
  value lower than the published one is the expected direction of any miss.
- **Reading.** PASS if all ten toleranced checks pass. A FAIL is recorded as a FAIL; the
  40 m neighbourhood is reported to help interpret it but does not change the status.
  ILL-V3 and ILL-V4 stay PENDING: ILL-V5 does not replace them, it adds an open-access
  comparison.
