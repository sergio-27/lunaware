# References

Status legend

- **[A]** Accessed during planning (2026-11-12); the facts quoted in PLAN.md were read from the source itself.
- **[S]** Only abstract, landing page or a secondary source was accessed. Numbers taken from it are flagged in PLAN.md.
- **[U]** Cited from prior knowledge; bibliographic details and any numbers must be verified against the source before the reference is used in code, docs or VALIDATION.md.

No number enters `VALIDATION.md` from a [U] reference. A case whose source stays [S]/[U] is reported as `PENDING`.

## 1. Data products

| Key | Reference | Status | Used for |
|---|---|---|---|
| PGDA-78 | NASA GSFC PGDA, *High-Resolution LOLA Topography for Lunar South Pole Sites*, https://pgda.gsfc.nasa.gov/products/78 ; data at https://pgda.gsfc.nasa.gov/data/LOLA_5mpp/ (README read). | A | 5 m/pix site DEMs (Site01, Site04, ...), `toterr`, `slperr`, 100 error clones per site |
| PGDA-81 | PGDA, *South Pole LOLA DEM Mosaic* (87–90°S, 5 m/pix, 3.3 GB), https://pgda.gsfc.nasa.gov/products/81 | A | Alternative near-field source |
| PGDA-90 | PGDA, *A New View of the Lunar South Pole from LOLA* (LDEM 60S/75S/80S/83S at 10–240 m/pix, `_ADJ`), https://pgda.gsfc.nasa.gov/products/90 ; data DOI 10.60903/gsfcpgda-lola-spole | A | Mid- and far-field DEMs for horizons |
| PGDA-69 | PGDA, *Lunar Polar Illumination* (`AVGVISIB_85S_060M_201608`, `_EARTH`, `LPSR_*`), https://pgda.gsfc.nasa.gov/products/69 | A | Published average-illumination and Earth-visibility maps (validation) |
| NAIF-GEN | NAIF generic kernels, https://naif.jpl.nasa.gov/pub/naif/generic_kernels/ : `lsk/naif0012.tls`, `spk/planets/de440s.bsp`, `pck/pck00011.tpc`, `pck/moon_pa_de440_200625.bpc`, `fk/satellites/moon_de440_250416.tf` (directory listings read). | A | Ephemerides, lunar orientation, frames |
| DIV-PCP | LRO Diviner Polar Cumulative Products, diurnal (summer/winter, 96 bins of 3.75 deg subsolar longitude, 2009-07-05 to 2019-02-17) and seasonal; 240 m/pix polar stereographic poleward of 80 deg. PDS bundle `urn-nasa-pds-lro_diviner_derived1`, collection `data_derived_pcp`, https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_pcp/ ; PDS3 labels and the Diviner RDR and Derived Products SIS (section 9) read 2026-11-15. | A | TH-V3, TH-V3b |
| HORIZONS | JPL Horizons system, https://ssd.jpl.nasa.gov/horizons/ ; API responses of 2026-11-12 stored in `src/lunarenv/validation/reference/`. | A | Independent check of topocentric Sun/Earth direction (ILL-X1) |

## 2. Topography and geodesy

| Key | Reference | Status |
|---|---|---|
| Barker2021 | Barker, M.K., et al. (2021). Improved LOLA elevation maps for south pole landing sites: Error estimates and their impact on illumination conditions. *Planet. Space Sci.* 203, 105119. doi:10.1016/j.pss.2020.105119 (NTRS manuscript read in full). | A |
| Barker2023 | Barker, M.K., et al. (2023). A new view of the lunar south pole from the Lunar Orbiter Laser Altimeter (LOLA). *Planet. Sci. J.* 4, 183. doi:10.3847/PSJ/acf3e1. Section 3 (uncertainty maps: rms height error per pixel from effective resolution, spline fits to clone ensembles) read 2026-11-16 from the open-access PDF. | A |
| Smith2010 | Smith, D.E., et al. (2010). The Lunar Orbiter Laser Altimeter investigation on the Lunar Reconnaissance Orbiter mission. *Space Sci. Rev.* 150, 209–241. | U |
| Smith2017 | Smith, D.E., et al. (2017). Summary of the results from the Lunar Orbiter Laser Altimeter after seven years in lunar orbit. *Icarus* 283, 70–91. | U |
| Zuber2012 | Zuber, M.T., et al. (2012). Constraints on the volatile distribution within Shackleton crater at the lunar south pole. *Nature* 486, 378–381. doi:10.1038/nature11216 (NTRS manuscript; Table 1 used for the georeferencing check in docs/data.md). | A |
| Park2021 | Park, R.S., et al. (2021). The JPL planetary and lunar ephemerides DE440 and DE441. *Astron. J.* 161, 105. doi:10.3847/1538-3881/abd414 | U |
| Acton1996 | Acton, C.H. (1996). Ancillary data services of NASA's Navigation and Ancillary Information Facility. *Planet. Space Sci.* 44, 65–70. | U |
| Acton2018 | Acton, C., et al. (2018). A look towards the future in the handling of space science mission geometry. *Planet. Space Sci.* 150, 9–12. | U |
| Annex2020 | Annex, A., et al. (2020). SpiceyPy: a Pythonic wrapper for the SPICE toolkit. *J. Open Source Softw.* 5(46), 2050. doi:10.21105/joss.02050 | U |
| Prsa2016 | Prša, A., et al. (2016). Nominal values for selected solar and planetary quantities: IAU 2015 Resolution B3. *Astron. J.* 152, 41. (nominal solar radius 695 700 km) | U |

## 3. Illumination

| Key | Reference | Status |
|---|---|---|
| Mazarico2011 | Mazarico, E., et al. (2011). Illumination conditions of the lunar polar regions using LOLA topography. *Icarus* 211, 1066–1081. doi:10.1016/j.icarus.2010.10.030. Abstract only; Table 3 (site list and per-site statistics) **not** accessed. | S |
| Mazarico2018 | Mazarico, E., et al. (2018). Advanced illumination modeling for data analysis and calibration. Application to the Moon. *Adv. Space Res.* 62, 3214–3228. (nested-grid horizon method cited by Barker2021) | U |
| Glaser2014 | Gläser, P., et al. (2014). Illumination conditions at the lunar south pole using high resolution Digital Terrain Models from LOLA. *Icarus* 243, 78–90. doi:10.1016/j.icarus.2014.08.013. Abstract only. | S |
| Glaser2018 | Gläser, P., Oberst, J., Neumann, G.A., Mazarico, E., Speyerer, E.J., Robinson, M.S. (2018). Illumination conditions at the lunar poles: Implications for future exploration. *Planet. Space Sci.* 162, 170–178. Accepted manuscript read in full from NTRS 20170007365 (2026-11-13); final journal version not checked. | A |
| Speyerer2013 | Speyerer, E.J., Robinson, M.S. (2013). Persistently illuminated regions at the lunar poles: Ideal sites for future exploration. *Icarus* 222, 122–136. doi:10.1016/j.icarus.2012.10.010 | S |
| Bussey2010 | Bussey, D.B.J., et al. (2010). Illumination conditions of the south pole of the Moon derived using Kaguya topography. *Icarus* 208, 558–564. doi:10.1016/j.icarus.2010.03.028 | S |
| Noda2008 | Noda, H., et al. (2008). Illumination conditions at the lunar polar regions by KAGUYA(SELENE) laser altimeter. *Geophys. Res. Lett.* 35, L24203. doi:10.1029/2008GL035692 | S |
| Espenak-LE2026 | Espenak, F. *Total Lunar Eclipse of 2026 Mar 03* (figure with contact times), NASA GSFC, https://eclipse.gsfc.nasa.gov/LEplot/LEplot2001/LE2026Mar03T.pdf | A |
| LimbDark | Solar limb-darkening law and coefficients (candidate sources: Neckel & Labs 1994, *Sol. Phys.* 153, 91; Hestroffer & Magnan 1998, *Astron. Astrophys.* 333, 338). Selection and verification in Phase 2. | U |

## 4. Thermal

| Key | Reference | Status |
|---|---|---|
| Hayne2017 | Hayne, P.O., et al. (2017). Global regolith thermophysical properties of the Moon from the Diviner Lunar Radiometer Experiment. *J. Geophys. Res. Planets* 122, 2371–2400. doi:10.1002/2017JE005387; arXiv:1711.00977 (Appendix A, Tables A1–A2 read in full). | A |
| heat1d | Hayne, P.O. *heat1d: thermal model for planetary science applications*. https://github.com/phayne/heat1d ; PyPI `heat1d` 0.3.2, installed and run on 2026-11-12 (TH-X1). GitHub version 0.4.1, commit `5b73a6c1ba5047f59f6210185b37539e46e68136`, run on 2026-11-15 (TH-X1b, default grid) and 2026-11-16 (TH-X1c, grids m=n=20 and 40). Only its outputs are stored; no heat1d code is included in lunarenv. | A |
| Vasavada2012 | Vasavada, A.R., et al. (2012). Lunar equatorial surface temperatures and regolith properties from the Diviner Lunar Radiometer Experiment. *J. Geophys. Res.* 117, E00H18. doi:10.1029/2011JE003987. Values used so far are those quoted in Hayne2017 Table A2. | S |
| Williams2017 | Williams, J.-P., et al. (2017). The global surface temperatures of the Moon as measured by the Diviner Lunar Radiometer Experiment. *Icarus* 283, 300–325. doi:10.1016/j.icarus.2016.08.012 | U |
| Williams2019 | Williams, J.-P., et al. (2019). Seasonal polar temperatures on the Moon. *J. Geophys. Res. Planets* 124, 2505–2521. doi:10.1029/2019JE006028 | S |
| Paige2010 | Paige, D.A., et al. (2010). Diviner Lunar Radiometer observations of cold traps in the Moon's south polar region. *Science* 330, 479–482. doi:10.1126/science.1187726 | U |
| WoodsRobinson2019 | Woods-Robinson, R., Siegler, M.A., Paige, D.A. (2019). A model for the thermophysical properties of lunar regolith at low temperatures. *J. Geophys. Res. Planets* 124, 1989–2011. (alternative c_p / K for T < 90 K). Its c_p is used through the DSNE table (DSNE-I); section 2.2 read 2026-11-15 from https://luna1.diviner.ucla.edu/~dap/pubs/096.pdf. Section 4.3 (k_am(T) fit to amorphous SiO2, eq. 30-32) read 2026-11-16 in text and as page images (PDF SHA-256 df61d7e12f922a3a87a170d5d1a3e3048d7e3238a8f96c5e6b91998ba922350d) for the low-temperature conductivity option; the sign of G_am is inferred from the source's Figure 6c, not printed (see docs/models/thermal.md). | A |
| DSNE-I | NASA SLS-SPEC-159 Rev I (27 Oct 2021), *Cross-Program Design Specification for Natural Environments (DSNE)*, §3.4.6.2.2, Table 3.4.6.2.2-1 (Woods-Robinson 2019 c_p, 10-300 K). NTRS 20210024522, PDF SHA-256 60f46c95ce4981ee48d58e516fe5c1150d184d3d88d3396ca87f94f094b031d4, read 2026-11-15. Rev J (M2M-30159, Dec 2025), cited by the review, was not found in NTRS; whether its table differs is not known. | A |
| KeihmLangseth1973 | Keihm, S.J., Langseth, M.G. (1973). Surface brightness temperatures at the Apollo 17 heat flow site: Thermal conductivity of the upper 15 cm of regolith. *Proc. Lunar Sci. Conf. 4th*, 2503–2513 (ADS 1973LPSC....4.2503K, scan read 2026-11-16). | A |
| LangsethHFE | Langseth, M.G., Keihm, S.J., Chute, J.L. Lunar heat-flow experiment, final report (date not recorded here). NTRS 19780010035 (read 2026-11-16, Table 1 and section V-A). | A |
| Keihm1973 | Keihm, S.J., et al. (1973a,b). Apollo 15 and 17 heat-flow probe surface and subsurface temperatures. As quoted in Hayne2017 Table A2. | S |
| Langseth1976 | Langseth, M.G., Keihm, S.J., Peters, K. (1976). Revised lunar heat-flow values. *Proc. Lunar Sci. Conf. 7th*, 3143–3171. As quoted in Hayne2017. | S |
| MartinezSiegler2021 | Martinez, A., Siegler, M.A. (2021). A global thermal conductivity model for lunar regolith at low temperatures. *J. Geophys. Res. Planets* 126, doi:10.1029/2021JE006829. Cited by the review; closed access, no open copy found (Semantic Scholar, OpenAlex, 2026-11-16); not read and not used. | U |
| Keihm1984 | Keihm, S.J. (1984). Interpretation of the lunar microwave brightness temperature spectrum: Feasibility of orbital heat flow mapping. *Icarus* 60, 568–589. (albedo vs incidence law) | S |
| KoppLean2011 | Kopp, G., Lean, J.L. (2011). A new, lower value of total solar irradiance. *Geophys. Res. Lett.* 38, L01706. doi:10.1029/2010GL045777 | U |
| Whipple1950, FountainWest1970, MitchelldePater1994, Ledlow1992, Hemingway1981, Carrier1973 | Sources of the conductivity, density and heat-capacity laws, as cited in Hayne2017 Appendix A. | S |
| AharonsonSchorghofer2006 | Aharonson, O., Schörghofer, N. (2006). Subsurface ice on Mars with rough topography. *J. Geophys. Res.* 111, E11007. (slope geometry formulas) | U |
| Schorghofer-PCC | Schörghofer, N. *Planetary-Code-Collection* (thermal and horizon codes, including Topo3D for terrain-scattered and emitted radiation), https://github.com/nschorgh/Planetary-Code-Collection ; no licence declared (GitHub API, 2026-11-15). Prior art, not used. | A |
| fluxpy | Potter, S. *fluxpy: fast thermal modeling and radiosity in Python*, https://github.com/sampotter/fluxpy ; no licence declared (GitHub API, 2026-11-15). Prior art, not used. | A |
| King2020 | King, O., Warren, T., Bowles, N., Sefton-Nash, E., Fisackerly, R., Trautner, R. (2020). The Oxford 3D thermophysical model with application to PROSPECT/Luna 27 study landing sites. *Planet. Space Sci.* 182, 104790. doi:10.1016/j.pss.2019.104790. Prior art (terrain shadowing and scattering), not used. | S |
| Burger2024 | Bürger, J., et al. (2024). A microphysical thermal model for the lunar regolith: investigating the latitudinal dependence of regolith properties. *J. Geophys. Res. Planets* 129, e2023JE008152. doi:10.1029/2023JE008152. Eq. 15 and section 4.2.2 (A0 0.09, a 0.10, b 0.02, highlands) read 2026-11-16 from the open copy at e-docs.geo-leo.de. Cold bound of the albedo envelope. | A |
| Foote2020 | Foote, E.J., Paige, D.A., Shepard, M.K., Johnson, J.R., Biggar, S. (2020). The bidirectional and directional hemispheric reflectance of Apollo 11 and 16 soils: laboratory and Diviner measurements. *Icarus* 336, 113456. doi:10.1016/j.icarus.2019.113456. Eq. 9 (Apollo 16 highland soil) read 2026-11-16 from https://luna1.diviner.ucla.edu/~dap/pubs/091.pdf. Hot bound of the albedo envelope; measured to 75 deg incidence. | A |
| HLS-UG-001 | NASA (2021). *Human Landing System Lunar Thermal Analysis Guidebook* (L-TAG), HLS-UG-001, baseline, 4 Jan 2021. NTRS 20210010030 (PDF SHA-256 e5885413...d50d). Section 11.1 (Earthshine: formula, 56 mW m-2 infrared and 100 mW m-2 reflected at peak, tilt) read 2026-11-16. Earthshine term of the thermal model; the Earth albedo, temperature and distance behind the peak values are not stated there. | A |
| Glenar2019 | Glenar, D.A., Stubbs, T.J., Schwieterman, E.W., Robinson, T.D., Livengood, T.A. (2019). Earthshine as an illumination source at the Moon. *Icarus* 321, 841–856. Cited by HLS-UG-001 for the Earthshine values; not read. | U |
| Ingersoll1992 | Ingersoll, A.P., Svitek, T., Murray, B.C. (1992). Stability of polar frosts in spherical bowl-shaped craters on the Moon, Mercury, and Mars. *Icarus* 100, 40–47. doi:10.1016/0019-1035(92)90016-Z. Bibliographic record checked (NTRS 19930034907, CaltechAUTHORS); full text not read. The single-bounce bowl expressions used in TH-A3 are derived in src/lunarenv/validation/bowl.py from the sphere view factor. | S |
| Bandfield2015 | Bandfield, J.L., Hayne, P.O., Williams, J.-P., Greenhagen, B.T., Paige, D.A. (2015). Lunar surface roughness derived from LRO Diviner Radiometer observations. *Icarus* 248, 357–372. Author manuscript https://gemelli.spacescience.org/jbandfield/publications/bandfield_lunar_surface_roughness.pdf (SHA-256 f9b83878...0f22d); abstract and section 2 read 2026-11-17. Roughness model and 20 deg RMS slope of TH-V3f. | A |
| Smith1967 | Smith, B.G. (1967). Geometrical shadowing of a random rough surface. *IEEE Trans. Antennas Propag.* 15, 668–671. Shadowing function as cited by Bandfield2015; the form used is checked numerically against its definition in tests. | S |
| DivinerSIS | UCLA Diviner team. *Diviner RDR and Derived Products Software Interface Specification*. PDS `urn-nasa-pds-lro_diviner_derived2/document/Diviner_RDR_SIS.pdf` (SHA-256 a2d2f02e...4a46). Table 1 (IFOV 6.7 x 3.4 mrad, 0.128 s integration) and sections 2.2.2, 2.2.7, 5.4 (nine-point footprint binning) and 9.1 read 2026-11-17. Footprint and binning model of TH-V3g. | A |
| Williams2023 | Williams, J.-P., et al. (2023). LRO Diviner polar mapping data products. Workshop on Ices in the Solar System, LPI Contrib. 2809, abstract 4046 (SHA-256 fa4fba8f...a629). Read 2026-11-17: GDR binned by projecting each observation's effective field of view; PCP from the same RDR data. TH-V3g. | A |
| CarslawJaeger1959 | Carslaw, H.S., Jaeger, J.C. (1959). *Conduction of Heat in Solids*, 2nd ed. Oxford. (analytic solutions) | U |

## 5. Numerical methods and verification practice

| Key | Reference | Status |
|---|---|---|
| Roache1994 | Roache, P.J. (1994). Perspective: a method for uniform reporting of grid refinement studies. *J. Fluids Eng.* 116, 405–413. (Richardson extrapolation, GCI) | U |
| McKay1979 | McKay, M.D., Beckman, R.J., Conover, W.J. (1979). A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. *Technometrics* 21, 239–245. (Latin hypercube) | U |
| CF | CF Metadata Conventions, https://cfconventions.org | U |
