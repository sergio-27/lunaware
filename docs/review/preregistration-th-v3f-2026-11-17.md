# Pre-registration: TH-V3f, sub-DEM roughness against Diviner at five polar sites

Written before the roughness model is implemented and before anything below has been run.

Origin: TH-V3 to TH-V3e compare the model with the Diviner summer Polar Cumulative Product
on 720 m blocks at five polar sites (the Shackleton rim block of TH-V3 and PGDA 78 Sites
1, 7, 11, 23). Over the sunlit summer bins the Hayne2017 albedo law is 23-51 K too cold at
every site; the hot bound of the albedo envelope (Foote 2020) is within 10 K at three
sites and 13-22 K too cold at two (Sites 1 and 7). Terrain radiation and Earthshine add
under 1 K (TH-V3d). The review (2026-11-13, Q3) names sub-DEM roughness as a candidate:
at grazing incidence, sunward facets below the 5 m DEM scale receive far more flux than
the DEM facet. The question here is whether roughness of the size measured by Diviner
elsewhere on the Moon closes the gap.

What is already known (not blind): every TH-V3 to TH-V3e number above, including the
per-bin Diviner values at all five blocks. No roughness result has been computed.

## Source (read 2026-11-17)

Bandfield, J.L., Hayne, P.O., Williams, J.-P., Greenhagen, B.T., Paige, D.A. (2015),
Lunar surface roughness derived from LRO Diviner Radiometer observations, Icarus 248,
357-372; author manuscript, https://gemelli.spacescience.org/jbandfield/publications/bandfield_lunar_surface_roughness.pdf,
SHA-256 f9b83878749b4164cef1ba820c055691219729c4c631c033727259f1ef40f22d. Read: abstract
and section 2 [A]. Nadir Diviner data (four equatorial regions) best match a 20 deg RMS
slope; multi-angle data 20-35 deg; roughness at millimetre scales (surfaces 0.5-5 mm
apart stay thermally isolated). Model: Gaussian slope distribution (Shepard 1995 Eq. 13,
theta0 = tangent of the RMS slope angle), random azimuths; each facet's radiance weighted
by its probability and by the dot product of its normal with the view direction; cast
shadows by Smith (1967), "any slope/azimuth surface that would otherwise be illuminated
has the same statistical chance of being within a shadow"; sloped facets receive
downwelling solar and infrared from adjacent surfaces in proportion to the fraction of
their hemisphere those fill, as if from a horizontal surface; shadowed facets at 100 K.
Their equations are images in the manuscript and were not transcribed; the forms below
are the standard ones for the stated assumptions, and are verified numerically in tests.
Their regions are equatorial, so 20 deg is not a polar measurement.

## Roughness model (fixed now)

For one DEM point (the receiver) with DEM normal n0:

- **Slopes.** Surface gradient relative to the plane normal to n0, two independent
  Gaussian components with equal variance, s^2 = <|grad h|^2> = tan^2(theta_rms) (the
  Shepard 1995 / Bandfield 2015 convention). Quadrature: |grad h|^2 / s^2 is
  exponentially distributed; 3 Gauss-Laguerre nodes in it times 12 azimuths (every 30
  deg), 36 facets, weights = node weight / 12.
- **Direct sunlight on a facet.** E(t) max(0, n_f . sun) (1 - A(i_f)), where E is the
  normal irradiance of the visible solar disc from the ground-level horizon on the nested
  DEM (as TH-V3), times the Earth occultation factor, and A the albedo law at the facet's
  own incidence i_f.
- **Cast shadows (Smith 1967, uncorrelated form).** At sun zenith theta_s relative to n0,
  a front-facing facet is lit with probability G1 = 1 / (1 + Lambda(nu)),
  nu = cot(theta_s) / (sqrt(2) sigma_c), sigma_c = s / sqrt(2) the RMS slope of one
  component, Lambda(nu) = (exp(-nu^2) / (nu sqrt(pi)) - erfc(nu)) / 2; G1 = 0 when the Sun
  is below the plane. The lit and shadowed parts of each facet are separate
  populations: radiance G1 L_lit + (1 - G1) L_shadow.
- **Adjacent surfaces (as Bandfield 2015).** A facet of slope beta sees the fraction
  F = (1 - cos beta) / 2 of its hemisphere filled by a surface parallel to the DEM plane
  that has the smooth DEM facet's incident sunlight and temperature: absorbed
  F [(1 - A_diffuse) A(i_0) E cos i_0 + eps * eps sigma T_0^4], A_diffuse the engine's
  hemispherical albedo, i_0 and T_0 the smooth facet's incidence and temperature (the
  smooth run below). Shadowed population: this flux alone (one run per slope node).
- **Temperatures.** Every lit facet and every shadow population runs the engine's 1-D
  model (Hayne2017 parameters, hourly, one-year spin-up), as TH-V3. Mean radiance per
  time: sum over facets of weight [G1 eps sigma T_lit^4 + (1 - G1) eps sigma T_shadow^4].
- **View.** Nadir along n0: facets weighted by their probability per unit area of the
  mean plane (the dot-product weighting of Bandfield 2015 for this view). No masking.
- Not included: lateral conduction between facets (Bandfield 2015: facets more than
  0.5-5 mm apart are isolated), correlation between shadowing and facet orientation,
  multiple scattering beyond the adjacent-surface term, emission angle of the actual
  Diviner observations (polar orbit, mostly nadir; not checked per bin).

## Forward model and comparison

- Smooth base: TH-V3's forward model at every 20 m in the nine cells of each block,
  unchanged (it must reproduce the TH-V3 / TH-V3c / TH-V3e sunlit-bin biases).
- Roughness change: at 4 receivers per cell (every 120 m, as TH-V3d), summer radiance of
  the rough surface minus that of the smooth DEM facet, per subsolar-longitude bin;
  mean over the 4 receivers of each cell, added to the cell's smooth radiance; block
  Tbol as before. Sunlit and complete bins as before (smooth model geometry).
- Laws: Hayne2017 (default) and the hot bound Foote 2020. RMS slope 20 deg (Bandfield
  2015 nadir best fit) and 35 deg (upper end of their multi-angle range).
- Sites: the TH-V3 Shackleton block and the four TH-V3e blocks (diviner_pcp_sites.json).

## Checks

- **TH-V3f**: Hayne2017 law with 20 deg RMS roughness, bias of block Tbol over the
  sunlit summer bins, |bias| <= 10 K (the TH-V3 tolerance) at each of the five sites.
  PASS if it holds at all five, FAIL otherwise.
- Informational, per site: Hayne2017 at 35 deg; hot law at 0 / 20 / 35 deg; smooth
  Hayne2017 (must equal TH-V3 / TH-V3e within 0.05 K); RMS over sunlit bins; bias over
  all summer bins; mean roughness change of block Tbol.
- Informational, equator: level ground at 0 deg latitude, equinox, Hayne2017 law, 20 deg
  roughness against smooth: change of the noon and midnight Tbol (the Hayne2017 law was
  fitted to equatorial Diviner data with a smooth model, so a large change there means
  the law and the roughness model double-count).
- Informational, quadrature: at one receiver (Shackleton centre cell), Hayne2017, 20 deg:
  sunlit-bin mean Tbol with 5 nodes x 24 azimuths against the default 3 x 12.

## Reading (fixed now)

- PASS: roughness of the size measured by Diviner at the equator explains the cold bias
  of the default law at all five polar blocks. Proposed to the project owner: a roughness
  option in the engine (not added by this case), with the equatorial check deciding
  whether the albedo law must be refitted with it.
- Hayne2017 improves by more than 10 K at every site but some site stays outside 10 K:
  roughness is a major part of the gap, not all of it; the 35 deg values are reported.
- Hayne2017 changes by less than 5 K at every site: roughness of this kind does not
  explain the gap; footprint, timing and the albedo law remain.
- Hot law with 20 deg roughness more than 10 K too warm at the sites where it passed
  smooth (Shackleton, Sites 11 and 23): the hot law and this roughness cannot both hold
  there.
- Equatorial noon change above 5 K at 20 deg: the default law and the roughness model
  are not independent; any engine option would need the law refitted with roughness.
- Limits: summer, sunlit bins, 720 m blocks; roughness parameter from equatorial
  regions; the change is estimated at 4 of 144 pixels per cell; Smith shadowing is an
  approximation that is poorest at grazing incidence, which is exactly the polar case.

## Amendment 1 (2026-11-17, before any polar run)

A smoke test of the implementation, run before the case, gave a 20 deg roughness change
of -3.2 K at equatorial noon and -33.6 K at equatorial midnight. The midnight value is an
artefact of the text above: with G1 = 0 when the Sun is below the plane, the whole night
radiance came from the cast-shadow population, which never sees the Sun and so stores no
daytime heat. Changed, before any polar receiver was computed: the lit / cast-shadow split
(G1 L_lit + (1 - G1) L_shadow) applies only while the facet faces the Sun and the Sun is
above the mean plane; otherwise the facet's radiance is that of its own run (which then
receives only the adjacent-surface flux). Nothing else changes. One receiver at the
Shackleton block was computed in the smoke test (all summer bins, mean Tbol per law); no
block bias or Diviner comparison was computed.
