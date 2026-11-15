# Pre-registration: TH-V3, polar temperatures against Diviner (finding F01)

Written before the case was implemented and before any Diviner temperature at the site
was looked at. The extraction script (`benchmarks/extract_diviner_pcp.py`) was run and
printed only bin counts. The model temperature at the test point is already known from
earlier work (yearly range 43-119 K with the incidence-dependent albedo, about 26 K
warmer on the mean with a constant albedo); that is why the data were not looked at
before this file was written.

Origin: finding F01 of `docs/review/scientific-review-2026-11-13.md`: no thermal
comparison exists at any latitude above 26 deg, while the polar result rests on
model-form choices worth 26-50 K.

## Reference data

- Diviner Polar Cumulative Product, diurnal, summer, south pole:
  `pcp_avg_tbol_pols_sum_slonNN_240.tab`, NN = 01..96 (PDS bundle
  `urn-nasa-pds-lro_diviner_derived1`, collection `data_derived_pcp`).
- Content, from the PDS3 labels and the Diviner RDR and Derived Products SIS, section 9.1
  (both read in full for the parts used): all nadir observations 2009-07-05 to
  2019-02-17 with subsolar latitude below the equator; radiances binned per channel
  (3-9) in 240 m polar stereographic cells (radius 1737.4 km) and 3.75 deg of subsolar
  longitude; the bolometric brightness temperature is derived from the binned
  radiances. Empty bins are omitted.
- Cells: the 3 x 3 cells centred on the cell that contains the test point
  (-2822.5, -9062.5) m, i.e. cell centres at x = -3120, -2880, -2640 m and
  y = -9360, -9120, -8880 m (a 720 x 720 m block on the Shackleton rim).
- The extracted values are stored in
  `src/lunarenv/validation/reference/diviner_pcp_shackleton.json` with their provenance.

## Forward model

- Pixels: the 5 m site DEM (`site04`), sampled every 20 m inside each cell (12 x 12
  pixels per cell, 1296 in all), each with its own horizon (ground level, default
  illumination settings: 1440 azimuths, per-strip disc profile, nested terrain), its own
  DEM normal (default baseline), and the Earth-occultation factor.
- Time: hourly from 2008-07-05 (one-year spin-up cycle) to 2019-02-18; the compared
  period is 2009-07-05 to 2019-02-17, as the data.
- Regolith: default parameters (Hayne2017 Table A1, emissivity 0.95).
- Bins: model samples with subsolar latitude < 0 (in the body-fixed frame of the
  engine), grouped by subsolar east longitude in 3.75 deg bins matching the files.
- Bolometric brightness temperature of the model: the radiance-equivalent mean,
  Tbol = (mean over pixels and samples of eps * T^4)^(1/4), consistent with radiances
  averaged before conversion. A gray emitter is assumed (no spectral emissivity).
- Block value: the nine cells combined with equal weight in radiance,
  Tbol_block = (mean over cells of Tbol^4)^(1/4), for both model and data. A bin enters
  the comparison only if all nine cells have data.
- Sunlit bins: bins in which the model's mean visible fraction of the solar disc over
  the block's pixels and the bin's samples is at least 0.5. Decided from the model
  geometry alone, before seeing the data.

## Cases, metric and tolerance

- **TH-V3**: incidence-dependent albedo of Hayne2017 (the engine default).
- **TH-V3b**: the same with a constant albedo of 0.12 (as TH-V2b).
- Metric: bias = mean over sunlit bins of (model - Diviner) block Tbol.
- Tolerance: |bias| <= 10 K.
- Why 10 K: no Diviner error budget for polar bolometric temperatures was found in the
  sources read, so this is not a calibration-derived tolerance. It is set so the test
  can separate the two albedo laws, which differ by about 26 K in the yearly mean at
  the test point: with a 10 K tolerance both cannot pass if their difference in the
  sunlit-bin mean exceeds 20 K. It also sets the level at which a law is usable for
  preliminary design (a 10 K bias).
- Informational, no tolerance: RMS over sunlit bins; bias over all summer bins; bias
  per cell over sunlit bins; bias in the bin with the highest Diviner value; number of
  bins used; the spread of the model's sunlit-bin block Tbol between individual summers
  (an indication of how much the unknown timing of the Diviner observations inside each
  bin can matter); the change in the centre cell's sunlit-bin mean when its pixels are
  sampled every 10 m instead of 20 m.

## Reading

- One law passes and the other fails: that law is supported at this site and season at
  the 10 K level; the other is not.
- Both pass: the comparison does not discriminate at this tolerance; the per-cell and
  maximum-bin values are reported to show whether a sharper test is possible.
- Both fail on the same side: the model form is wrong beyond the albedo law (sub-pixel
  roughness, terrain scattering and emission, the 0 m horizon on a 5 m DEM, which
  ILL-V5 showed to be resolution-sensitive, or the low-temperature properties).
- Both fail on opposite sides: the truth lies between the two laws; neither is usable as
  is.
- The status of TH-V3 is the outcome of its own check; no tolerance is changed after
  the run, and a FAIL is reported as a FAIL.
- Not covered: winter, shadowed or crater-floor cells, depth temperatures, other sites.
