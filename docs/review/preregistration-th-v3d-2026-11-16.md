# Pre-registration: TH-V3d, Diviner polar temperatures with terrain radiation and Earthshine

Written before the case is implemented and before any forward model with terrain
radiation has been compared with the Diviner data.

Origin: TH-V3 (pre-registered 2026-11-15) was run without terrain-scattered sunlight,
terrain infrared and Earthshine, which review change 3 added to the engine on 2026-11-16
(TH-A3, TH-C1). TH-V3 lists them under "Not modelled" and reads the cold bias over all
summer bins (-39 / -13 K) as "expected without terrain radiation". This case measures how
much they change the comparison. TH-V3, TH-V3b and TH-V3c stay recorded as they are.

What is already known (not blind): the TH-V3, TH-V3b and TH-V3c results (Hayne2017 law
-23.8 K over the 8 sunlit bins, FAIL; hot law -0.4 K; all-summer-bin biases -39 / -13 K);
TH-C1 (terrain irradiance at the test point at its minimum 8.5e-5 W m-2). A timing trial
(2026-11-16) at two block pixels gave a summer-mean terrain irradiance of 8.6e-4 W m-2 at
the test point and 0.23 W m-2 at (-2722, -9142) m, with the Hayne2017 law; no temperature
from it was compared with Diviner.

## Why a subsample

Per pixel, single-bounce terrain radiation over the 10.6-year axis costs about 35 s of
patch geometry and 25-80 s of patch thermal model per albedo law with 7 processes (timing
trial). For the 1296 pixels of TH-V3 that is over 100 hours. The case therefore computes
the change caused by the diffuse terms on a stratified subsample and adds it to the
full TH-V3 forward model.

## Method

- Everything of TH-V3 is kept: Diviner data, 3 x 3 cells of 240 m, 144 pixels per cell
  for the base model, hourly 2008-07-05 to 2019-02-18, summer samples binned by subsolar
  longitude, radiance averaging, block Tbol, sunlit bins (model block mean visible
  fraction >= 0.5, from the base model's geometry), all-nine-cells requirement.
- Receivers: 4 per cell, the pixel centres at +/-60 m from the cell centre in x and y
  (the 2 x 2 grid of 120 m spacing), 36 in all. Each gets its own ground-level horizon,
  DEM normal and Earth-occultation factor exactly as in TH-V3, plus:
  - single-bounce terrain radiation with the engine's method and default settings
    (`site_radiation`, ground view, 360 rays, patch layout refinement 1, 3 points per
    patch, patches computed when they can send more than 1e-4 W m-2 to the receiver),
    patch temperatures from the same 1-D model and the same albedo law as the receiver,
    spun up on the first year of the axis;
  - Earthshine (HLS-UG-001 section 11.1, as in the engine);
  - absorbed diffuse flux as in the engine: (1 - A_diffuse) E_ts + eps (E_ti + E_ei) +
    (1 - A(theta_E)) E_er.
- For each receiver and law, the 1-D model is run with and without the diffuse flux; the
  change of the summer radiance eps T^4 is summed per bin. Per cell, the radiance change
  is the mean over its 4 receivers. TH-V3d's model cell radiance per bin = TH-V3's
  (144 pixels) + that change.
- Laws: Hayne2017 (TH-V3) and the hot envelope law Foote2020 (TH-V3c's warm edge).
  The constant albedo and the cold law are not repeated (cost).

## Check

- **TH-V3d**: Hayne2017 law with the diffuse terms, bias of block Tbol over the sunlit
  summer bins (model - Diviner), |bias| <= 10 K, the TH-V3 tolerance. PASS / FAIL on this
  check only.
- Informational: RMS over sunlit bins; bias over all summer bins with data; the same two
  for the hot law with the diffuse terms; the change of each from its TH-V3 / TH-V3c value;
  per-cell change of the sunlit-bin and all-bin mean Tbol; mean absorbed diffuse flux over
  the summer samples of the 36 receivers (terrain and Earthshine separately) and its
  largest receiver value; on the same receivers, the change of the sunlit-bin and all-bin
  block Tbol from the low-temperature conductivity option (Hayne2017 law, with the diffuse
  terms), added in the same way.

## Reading

- From the timing trial and TH-C1 the expectation is a small change on the sunlit bins
  (well under 1 K) and a larger one on the bins without direct sunlight. If TH-V3d still
  fails by about the TH-V3 amount, the Hayne2017 law's cold bias at this site is not due
  to the missing terrain radiation or Earthshine.
- A change of the sunlit-bin bias of more than 5 K would contradict TH-C1 and the bowl
  case and would be investigated before any other use of the terrain term.
- On the all-summer-bin bias, the result says how much of TH-V3's -39 K is explained by
  diffuse light; the remainder points to sub-DEM roughness, the 0 m horizon on a 5 m DEM,
  the Diviner footprint and observation timing, or the low-temperature properties.
- Limits: the subsample (4 of 144 pixels per cell) estimates the cell's change, not its
  temperature; multiple scattering and sub-DEM roughness are still absent.
