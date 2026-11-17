# Pre-registration: TH-V3h, Diviner observation timing at five polar sites

Written before the model has been evaluated at the observation times.

Origin: TH-V3 to TH-V3g compare the model with the Diviner summer Polar Cumulative Product
(PCP) on 720 m blocks at five polar sites. The model's bin value is the mean over every
hour of the ten southern summers 2009-2019 whose subsolar longitude falls in the 3.75 deg
bin. Diviner's bin value is the mean over the observations actually made. Over the sunlit
summer bins the Hayne2017 law is 24-51 K too cold at every site. Neither sub-DEM roughness
(TH-V3f, 11-18 K, not enough) nor the footprint and binning (TH-V3g, -5 to +2 K) explain
it. The question here is whether the times at which Diviner observed each block change
the model's bin values by an amount that matters.

What is already known (not blind): every TH-V3 to TH-V3g number, including the per-bin
Diviner values at all five blocks. TH-V3 found a spread between individual summers of the
model's sunlit-bin block Tbol of about 1 K (median of the per-bin standard deviation).
No model value at the observation times has been computed.

## Observation times (extracted 2026-11-18, before this was written)

`benchmarks/extract_diviner_rdr_times.py` queried the PDS Geosciences Node Orbital Data
Explorer Diviner RDR point service (the service behind the "DIVINER RDR Query V2" page)
for channel 7 records in a latitude/longitude box around each block (the 3 x 3 cells
grown by 240 m), 2009-07-05 to 2019-02-17. Kept: activity flag 110 (on Moon, standard
nadir), subsolar latitude below 0, brightness temperature 10-450 K, the selection the SIS
states for the gridded products (sections 2.2.2 and 9.1). Stored per observation: UTC
Julian date, polar stereographic position of the centre, RDR subsolar longitude, orbit.
Radiances and temperatures were neither stored nor printed. Output:
`src/lunarenv/validation/reference/diviner_rdr_times.json.gz`.

Kept observations per year (2009 to 2018), counts only:

| Block | Kept | Per year 2009-2018 |
|---|---|---|
| Shackleton | 8301 | 3857, 4444, 0, ... |
| Site 1 | 11442 | 7104, 1460, 2878, 0, ... |
| Site 7 | 7194 | 2123, 1557, 2250, 1264, 0, ... |
| Site 11 | 5382 | 968, 866, 1353, 1849, 346, 0, ... |
| Site 23 | 12188 | 864, 816, 690, 726, 1031, 2076, 1722, 1733, 2530, 0 |

Standard nadir observations of the four blocks poleward of 88.5 deg S stop between 2011
and 2013 (the cause, presumably the orbit's inclination, was not checked). Coverage, from
positions and bins only: assigning each observation to every cell one of its nine IFOV
points falls in (SIS section 5.4, geometry as TH-V3g) gives observations in all nine cells
in exactly the bins where the PCP has all nine cells (64, 68, 56, 42, 63 bins). Assigning
by the centre only leaves 0-13 PCP cell-bins per site without an observation. This
supports the nine-point binning for the PCP, assumed in TH-V3g. Median 4-6 observations
per cell and bin.

## Model at the observation times (fixed now)

- TH-V3's forward model at each block: every 20 m (144 pixels per 240 m cell), own
  horizon and DEM normal, 1-D model hourly 2009-07-05 to 2019-02-17 after a one-year
  spin-up. Laws: Hayne2017 (default) and the hot bound Foote 2020.
- Radiance eps T^4 of each pixel interpolated linearly between the hourly samples to each
  observation time; the cell's value at that time is the mean over its 144 pixels.
- An observation is binned by its RDR subsolar longitude (3.75 deg bins) and counted in
  every cell that one of its nine IFOV points falls in, once per point (TH-V3g geometry:
  335 x 170 m at 50 km, in track along the meridian). The cell's bin value is the mean of
  the cell's values at the times of the observations counted in it. Block Tbol from the
  nine cells as before.
- The observation's own footprint average is not applied (TH-V3g covers it); the value
  counted is the cell's mean at the observation time. This isolates timing.
- Sunlit bins: unchanged (model block visible fraction >= 0.5 from every summer hour, and
  Diviner in all nine cells). Compared: sunlit bins that also have observations in all
  nine cells (by the coverage above, all of them).

## Checks

- **TH-V3h**: Hayne2017 at the observation times, bias of block Tbol over the sunlit
  summer bins, |bias| <= 10 K (the TH-V3 tolerance) at each of the five sites. PASS if it
  holds at every site with at least 4 such bins, FAIL otherwise (as TH-V3e).
- Informational, per site:
  - hot law at the observation times;
  - both laws over every summer hour on the same bins, and over all sunlit bins (must
    equal TH-V3 / TH-V3c / TH-V3e within 0.05 K);
  - change from timing (observation times minus every hour, same bins), both laws;
  - bound: both laws with the hottest summer hour of each pixel in each bin (no timing of
    the observations can give a warmer bin value), over all sunlit bins;
  - Hayne2017 at the observation times binned by the centre only, against nine-point, on
    the bins complete both ways;
  - number of observations binned in the block in the bins used, and how many before
    2012;
  - fraction of those observations whose RDR subsolar longitude falls in the bin of the
    model's subsolar longitude (SPICE) at the same time (a check of time and frame).

## Reading (fixed now)

- PASS: the observation timing, not the physics, explains the cold bias of the default
  law; TH-V3 to TH-V3g need re-reading at the observation times.
- Change negative at a site: timing makes the model colder there and cannot be the cause.
- |change| under 5 K at every site, both laws: timing does not explain the gap; with
  roughness, footprint, terrain radiation and timing each under the gap, the remaining
  candidates are the albedo law (and other surface properties) and the illumination of
  the blocks.
- Hottest-hour bound below -10 K at a site: no timing of the observations can bring the
  default law within 10 K there.
- Hot law at the observation times still below -10 K at Sites 1 and 7: timing does not
  explain why Diviner is above the hot bound there.
- If centre-only binning moves a site's bias by more than the timing change, the result
  there depends on the binning assumption and is reported as such.
- Expectation, stated before the run: the change is under 5 K at every site (the
  per-summer spread is about 1 K and each bin holds observations from several lunations);
  the hottest-hour bound is above the timed value but still below -10 K at Sites 1 and 7.
- Limits: channel 7 times stand for all channels; the RDR subsolar longitude defines the
  bin; cell mean instead of the footprint average; flat footprint geometry; 50 km
  altitude for the binning points; summer, sunlit bins, 720 m blocks.
- Run time estimate: about 1.5 h (1296 pixels per site, two laws).
