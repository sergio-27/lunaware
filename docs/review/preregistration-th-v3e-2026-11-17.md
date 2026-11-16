# Pre-registration: TH-V3e, Diviner polar temperatures at four more sites

Written before the case is implemented, before the new site DEMs have been downloaded and
before any Diviner value outside the TH-V3 block has been read.

Origin: TH-V3 to TH-V3d compare the model with Diviner on one 720 m block on the
Shackleton rim. There the Hayne2017 albedo law is 23-24 K too cold over the 8 sunlit
summer bins (FAIL), the hot bound of the albedo envelope (Foote 2020) is within 0.5 K
and the cold bound (Bürger 2024) 34 K too cold; terrain radiation and Earthshine change
this by under 1 K (TH-V3d). Every run reports temperatures under both envelope bounds.
The question here is whether the Shackleton result holds at other polar sites: if the
hot bound matches Diviner everywhere and the cold bound is too cold everywhere, the
cold bound is not realistic for sunlit summer conditions and the reported bounds could
be narrowed. One site with 8 bins cannot support that.

What is already known (not blind): the TH-V3 to TH-V3d results above; the model's
illumination at Site 1 from ILL-V1 and ILL-V5 (Gläser 2018 points C1). No Diviner value
outside the TH-V3 block, and no thermal result at any of the new sites, has been seen.

## Sites (fixed now)

PGDA 78 5 m site DEMs (Barker2021; https://pgda.gsfc.nasa.gov/products/78), chosen to
span latitude and setting, before their terrain has been looked at:

| Site | Name on the PGDA page | Catalog id |
|---|---|---|
| Site01 | Connecting ridge | `pgda78.site01.surf` (already cached) |
| Site07 | Peak near Shackleton | `pgda78.site07.surf` |
| Site11 | de Gerlache rim | `pgda78.site11.surf` |
| Site23 | Malapert massif | `pgda78.site23.surf` |

Regional terrain: the 80 m PGDA 90 DEM, as everywhere else. If a site DEM cannot be read
or its tile does not allow a block as defined below, the site is dropped and the drop is
reported; no replacement site is chosen.

## Block selection (fixed now, model geometry only)

- Candidate cells: the Diviner PCP 240 m cells (polar stereographic, radius 1737.4 km,
  cell centres at integer multiples of 240 m) whose centres lie at least 1 km inside the
  site tile. Candidate blocks: 3 x 3 such cells.
- For each candidate cell centre, the ground-level horizon (default illumination
  settings) and the mean visible fraction of the solar disc over the summer samples of
  the TH-V3 period (hourly, 2009-07-05 to 2019-02-18, subsolar latitude < 0).
- The block with the highest mean of its nine cell-centre values is used. Blocks that
  share a cell with the TH-V3 block are excluded. Ties: the first in row-major order
  of the block centre (largest y first, then smallest x).
- The selection is recorded (block centre, value) before the Diviner values are
  extracted. The extraction prints only cell counts.

## Forward model and comparison (as TH-V3)

- Pixels every 20 m inside each cell (144 per cell, 1296 per block), each with its own
  ground-level horizon on its site's nested terrain, its own DEM normal and the Earth
  occultation factor; hourly from 2008-07-05 (one-year spin-up) to 2019-02-18; summer
  samples binned by subsolar longitude (96 bins of 3.75 deg); radiance eps T^4 averaged
  over pixels and samples; cell Tbol; block Tbol as the radiance mean of the nine cells;
  a bin enters only if all nine cells have Diviner data.
- Sunlit bins: model block mean visible fraction >= 0.5 (model geometry only).
- Laws: Hayne2017 (default), constant 0.12, cold bound Bürger 2024, hot bound Foote 2020,
  all with the other Hayne2017 parameters, as in TH-V3 / TH-V3c.
- Not included: terrain radiation and Earthshine (TH-V3d: under 1 K at Shackleton),
  sub-DEM roughness, spectral emissivity, Diviner footprint and timing.

## Checks

- **TH-V3e**: per new site, the hot bound's bias of block Tbol over the sunlit summer
  bins (model - Diviner), |bias| <= 10 K, the TH-V3 tolerance. PASS if it passes at every
  site that qualifies; FAIL otherwise.
- A site qualifies if at least 4 sunlit bins have data in all nine cells. A site with
  fewer is reported (informational) and does not enter the verdict. If fewer than 3 of
  the 4 sites qualify, the status is PARTIAL.
- Informational, per site: bias for the other three laws; RMS over sunlit bins; bias over
  all summer bins; number of sunlit and complete bins; fraction of sunlit bins in which
  Diviner lies between the cold and hot bounds; block latitude and mean visible fraction.

## Reading (fixed now)

- Hot bound passes at every qualifying site and the cold bound is more than 10 K too cold
  at every qualifying site (and at Shackleton, TH-V3c): Diviner supports the hot side of
  the envelope in sunlit summer conditions at five sites. A change to the reported
  thermal bounds (for example, labelling the cold bound as excluded by Diviner for sunlit
  conditions) is then proposed to the project owner; this case does not change the
  engine.
- Hot bound fails at some sites: the Shackleton agreement is site-specific; the bounds
  stay as they are and the per-site pattern (latitude, setting) is reported.
- Hayne2017 within 10 K at some site: its cold bias is site-dependent, not a property of
  the law at grazing incidence.
- Limits: summer only, sunlit bins only (the night and shadow temperatures that set the
  design minimum are not tested); block-scale (720 m) temperatures, not pixel ones;
  sites are rims and massifs chosen for illumination, not a random sample of the pole.
