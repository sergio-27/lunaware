# Pre-registration: TH-V3g, Diviner footprint and binning at five polar sites

Written before the footprint model is implemented and before anything below has been run.

Origin: TH-V3 to TH-V3f compare the model with the Diviner summer Polar Cumulative Product
(PCP) on 720 m blocks at five polar sites. The model cell is the mean over the cell's own
240 m square. Over the sunlit summer bins the Hayne2017 law is 24-51 K too cold at every
site, and still 13-34 K too cold with 20 deg sub-DEM roughness (TH-V3f); the hot bound
(Foote 2020) is 13-22 K too cold at Sites 1 and 7. With the same constant albedo 0.12 the
sunlit-bin bias ranges from -15.7 K (Site 1) to +6.1 K (Site 11) (TH-V3b, TH-V3e), so part
of the gap depends on the site. The question here is whether the way Diviner measures and
bins radiance (a footprint larger than the cell in track, assigned to every bin it
touches) changes the comparison by an amount that matters.

What is already known (not blind): every TH-V3 to TH-V3f number, including the per-bin
Diviner values at all five blocks. No footprint result has been computed.

Expectation, stated before the run: the blocks were chosen as the best-lit 720 m blocks
in their areas (TH-V3e) or sit on the Shackleton rim, so neighbouring ground is on average
darker and the footprint should make the model colder, i.e. widen the gap. A result in
that direction rules the footprint out as a cause.

## Sources (read 2026-11-17)

- Diviner RDR and Derived Products SIS (UCLA), PDS
  `urn-nasa-pds-lro_diviner_derived2/document/Diviner_RDR_SIS.pdf`, SHA-256
  a2d2f02ec364021b0a3ff5d475395a8126192e32de94500e9a13d197d2974a46. Read [A]:
  - Table 1: detector geometric IFOV 6.7 mrad in track and 3.4 mrad cross track,
    "320 m on ground in track for 50 km altitude, 160 m ... cross track"; 0.128 s signal
    integration period; primarily nadir pushbroom mapping.
  - Section 2.2.2: for the gridded maps "the finite field of view of the Diviner
    footprints is also taken into account", footprints projected onto a LOLA DEM.
  - Section 5.4: the program "pfootprint" adds eight points to each observation (the
    +/- X and +/- Y edge midpoints and the four corners of the IFOV rectangle at the
    spacecraft altitude, rotated to the detector orientation); the points are "added to
    the pipes data stream and further processed along with the rest of the data",
    populating every bin a footprint touches.
  - Sections 2.2.7 and 9.1: PCP (diurnal) bins RDR radiances of all nadir observations
    2009-07-05 to 2019-02-17 in 240 m polar stereographic bins; they do not say whether
    the footprint points of section 5.4 are used for the PCP.
- Williams, J.-P., et al. (2023), LRO Diviner polar mapping data products, Workshop on
  Ices in the Solar System, LPI Contrib. 2809, abstract 4046, SHA-256
  fa4fba8f96721e4c91dfe02cd312c361674c583da233c8d5b0f37596382aa629. Read [A]: GDR data
  "binned by projecting the modeled effective field of view of each RDR observation onto
  a Digital Elevation Model"; the PCP "combine multiple years of observations" from the
  same RDR data (no separate binning method stated).
- Williams, J.-P., et al. (2019), Seasonal polar temperatures on the Moon, JGR Planets
  124, 2505, doi:10.1029/2019JE006028. Not read (publisher refused access). A search
  summary attributes to it that the south polar footprint stayed about 140 m cross track
  and 400 m in track at 50 km altitude because the altitude over the south pole varied
  little [S]. Used only to support the 50 km altitude, not as a number in the model.

## Footprint and binning model (fixed now)

- **Assumption (not stated in the sources):** the PCP use the section 5.4 binning. Each
  observation contributes its radiance once for each of its nine points (centre, four
  edge midpoints, four corners of the geometric IFOV rectangle) that falls in a bin; a
  bin's value is the mean over those contributions. This is the reading of "added to the
  pipes data stream" as separate data points. Observation centres are taken as uniformly
  distributed over the area (ten years of tracks).
- **Geometry:** altitude 50 km; rectangle 335 m in track x 170 m cross track
  (6.7 x 3.4 mrad x 50 km; the SIS rounds them to 320 x 160 m, the rounded values are not
  used).
- **Radiance of one observation:** the mean of the model radiance (eps T^4, per
  subsolar-longitude bin) over the footprint response: uniform cross track over 170 m; in
  track the 335 m rectangle convolved with the motion during the 0.128 s integration,
  0.128 s x ground speed. Ground speed at 50 km: sqrt(GM / r) R / r, GM = 4.9028e12 m^3
  s^-2 (standard lunar value [U]), R = 1737.4 km, r = R + 50 km: 1.61 km/s, 206 m. A
  trapezoid in track, full extent 541 m.
- **Orientation:** in track along the local meridian (a polar orbit), i.e. along the
  polar stereographic radius through the block centre; one orientation per block.
- **Model bin value of a cell:** (1/9) sum over the nine points p of the mean over
  observation centres c in (cell - d_p) of the observation radiance at c, with d_p the
  point's offset. Centres every 20 m. Block Tbol from the nine cells as before.
- **Model field:** TH-V3's forward model (every pixel its own horizon, normal, 1-D model;
  hourly 2009-07-05 to 2019-02-17), per pixel and bin: every 20 m inside the block (as
  TH-V3) and every 60 m in a 660 m wide ring around it, the ring values spread over the
  20 m grid by nearest neighbour. 660 m covers the farthest footprint response point:
  bin-point offset up to 188 m plus response half-diagonal 284 m (426 m for the 1.5 x
  sensitivity case below), 614 m in all.
- Sunlit and complete bins as before, from the smooth model geometry inside the block
  (unchanged, so the bins compared are TH-V3e's).
- Laws: Hayne2017 (default) and the hot bound Foote 2020.
- Not included: the detector's actual spatial response (blur beyond the geometric IFOV),
  altitude variation, pointing error, the DEM projection of the footprint (a flat
  footprint is used), the actual distribution of tracks and their timing in each bin,
  roughness (TH-V3f) together with the footprint.

## Checks

- **TH-V3g**: Hayne2017 law with the footprint and binning model, bias of block Tbol over
  the sunlit summer bins, |bias| <= 10 K (the TH-V3 tolerance) at each of the five sites.
  PASS if it holds at every site with at least 4 sunlit bins, FAIL otherwise (as
  TH-V3e).
- Informational, per site:
  - hot law with the footprint;
  - both laws without it (must equal TH-V3 / TH-V3e within 0.05 K);
  - change of block Tbol from the footprint, both laws;
  - Hayne2017 with the footprint response 1.5 times larger in both directions (binning
    points unchanged: they follow the geometric IFOV by definition);
  - Hayne2017 with the in-track direction averaged over 12 orientations instead of the
    meridian;
  - Hayne2017 with the centre point only (no eight extra points; footprint response
    kept), to separate binning from footprint averaging;
  - ring sampling: the block's Hayne2017 sunlit-bin Tbol from the inner pixels on the
    60 m sub-grid against the 20 m grid (the sampling error of the ring);
  - RMS over sunlit bins and bias over all summer bins with data.

## Reading (fixed now)

- PASS: the comparison method, not the physics, explains the cold bias of the default
  law; TH-V3 to TH-V3f would need re-reading with the footprint.
- Footprint change negative (colder) at a site: it widens the gap there and cannot be its
  cause at that site.
- |change| under 5 K at every site, both laws: the footprint and binning do not explain
  the gap at 720 m block scale; the remaining candidates are the albedo law, the
  observation timing within the bins and the illumination of the blocks.
- Hot law with the footprint still below -10 K at Sites 1 and 7: the footprint does not
  explain why Diviner is above the hot bound there.
- If the 1.5 x response or the orientation average moves a site's bias by more than the
  nominal footprint change, the result is sensitive to unknown instrument details and is
  reported as such.
- Limits: the PCP binning method is an assumption from the GDR description; flat
  footprints; summer, sunlit bins, 720 m blocks; ring at 60 m.
- Run time estimate: about 4 h (2300 pixels per site, two laws).
