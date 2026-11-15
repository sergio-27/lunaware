# Scientific review of lunarenv: is the approach sound enough to build on?

Reviewer: Claude Science (independent review, requested per `scientific-review-brief.md`).
Date: 2026-11-13. Scope: the scientific approach only (no code, style, performance or architecture).
Inputs: the lunarenv documents supplied (README, PLAN, VALIDATION, models/illumination, models/thermal, convergence/*, references), the literature and data listed in section 4, and a small independent thermal calculation by the reviewer (section 3, "Reviewer calculation", and Appendix A). Its numbers come from a reduced model, not from lunarenv, and are labelled as such wherever they appear.

Source marks: **[A]** the relevant part of the source was read in full text by this review. This includes reading done by the review's own research agents; values they extracted from data files are marked *provisional*. **[S]** only an abstract, landing page, search snippet or secondary source was read. **[L]** a lunarenv document supplied with the brief.

---

## 1. Verdict (Q6)

**Go with conditions.** The illumination half is sound: it uses the published method (Barker 2021), is validated against two independent published products, and agrees with a third that the author has not used yet (Gläser et al. 2018: 85.5 % at a point 21 m from the test point). It can be used now for Sun and Earth geometry and the visible solar fraction, with its clone spread. The thermal half is a well-verified solver applied to a physical model that has not been validated at the pole. Its central input, the albedo at 88° incidence, is not established in the literature, and the bound offered for it (constant albedo 0.12) lies outside the published envelope. It must not be presented as a prediction until it has been compared with polar Diviner data.

The three changes that would most increase scientific credibility, in order:

1. **Run TH-V3** against the Diviner Polar Cumulative Products, which are archived (path in section 3, Q4). Use summer, sunlit rim cells, with model output aggregated from 5 m to 240 m.
2. **Replace "nominal + constant albedo" with a bounded, orientation-resolved product.** Use the DSNE Rev J albedo envelope crossed with the existing ensemble.
3. **Add terrain-scattered sunlight, terrain infrared and fluxes on oriented surfaces** (including vertical faces). Fix the two non-like-for-like comparisons: ILL-V2 (metric) and TH-V2 (brightness temperature, site albedo).

**Minimum before an engineer uses a temperature in a design decision:**
- TH-V3 passes at the site against a tolerance registered before the run.
- The value comes from the bounded envelope, for the actual orientation of the design surface.
- It is used only as a boundary condition for a coupled asset model, per HLS-UG-001, not as a structure temperature.

---

## 2. Findings

Severity: **Critical**, the approach cannot support its stated use without this change. **Major**, results are presented with more confidence than the evidence supports. **Minor**, a gap in justification or documentation that does not change results.

| ID | Q | Finding | Severity | Evidence and source | Recommended action |
|---|---|---|---|---|---|
| F01 | Q4/Q5 | No thermal comparison exists at any latitude above 26°, yet the polar result rests on model-form choices worth 26-50 K. The data to test it exist and are archived. | **Critical** | All thermal cases are equatorial or at 20-26° N [L: VALIDATION]. Diviner PCP: binned bolometric brightness temperature at 240 m, 80° S to the pole, 6 Ls × 24 subsolar-longitude bins plus summer local-time bins, 2009-07-05 to 2020-03-15 [A: PDS4 label `pcp_avg_tbol_pols_ls030_slon01_240.xml`]. DSNE Rev J: "None at polar regions (relies on modeling)" for subsurface data [A: Moore 2026 slides]. | Run TH-V3 as specified in section 3, Q4, with its tolerance committed before the run. Until then, label every polar temperature "scenario, not validated". |
| F02 | Q1 | Level-ground regolith temperature at a point is a boundary condition, not the design environment of a structure. At 1.8° Sun elevation a vertical sun-facing panel receives 1360 W m⁻², against 42.8 W m⁻² on level ground (×32). | **Major** | DSNE Rev J §3.4.6.5: natural temperatures are "of limited accuracy" for a deployed asset; it defers to HLS-UG-001 [A]. HLS-UG-001 §9.0: "it is generally necessary to calculate surface temperatures as opposed to setting surface temperatures to boundary temperatures" [A]. The asset changes local regolith temperatures by more than 200 °C and nozzle temperatures by almost 50 °C (Gasbarre 2008) [A]. Terrain IR cuts radiator rejection by up to 30-31 % (horizontal) and 68-70 % (20° tilt) (Birmingham 2024) [A]. Flux ratio: reviewer arithmetic. | Output Sun and Earth vectors and direct flux on user-defined surface orientations, including vertical. Output sky and terrain view factors from the horizons already computed. State that the ground temperature is an input to a coupled asset model. |
| F03 | Q3/Q5 | The albedo at 88° incidence is not established. The constant albedo 0.12 used as the alternative lies **outside**, on the hot side, of the published design envelope. "Treat them as bounds" is therefore mis-specified. | **Major** | DSNE Rev J Table 3.4.6.1.1-1 highland envelope at θ = 88.2°: 0.238 (Foote 2020 min) to 0.860 (Bürger 2024 max). Hayne 2017 gives 0.784 [A; values evaluated by the reviewer from the DSNE formulas]. DSNE: constant albedo is hot-biased, +28 K at 84° S [A]. Reviewer model at the test point: envelope mean 62.9-91.9 K, max 101.5-156.7 K. Constant 0.12 gives 94.8 / 162.6 K. | Report the DSNE Rev J highland envelope as the thermal bounds. Keep constant 0.12 only as a labelled sensitivity. |
| F04 | Q3/Q4 | TH-V1 does not discriminate between the albedo laws, so it is not evidence for A(θ). | **Major** | Reviewer model, equator, constant 0.12: 385.2 / 101.7 / 94.9 K. A(θ): 385.2 / 100.1 / 93.5 K. Reference 385 / 101 / 95 K. Both pass a 2 K tolerance (Appendix A). | Run TH-V1 in lunarenv with the alternative law and report it. If it passes, state that TH-V1 does not discriminate. If it misses by more than 2 K, the reviewer's reduced model is wrong on this point. |
| F05 | Q3/Q4 | TH-V2 (Apollo) is not like-for-like. Its FAIL mixes three effects, so neither the FAIL nor TH-V2b's PASS settles the albedo question. | **Major** | Apollo 17's 216 K is a diurnal mean *brightness* temperature from cable thermocouples ~15 cm above the surface, reduced with ε = 1.0 and A = 0.08 (Keihm & Langseth 1973) [A, from an OCR'd scan]. The physical equivalent at ε = 0.95 is ~219 K (reviewer arithmetic). Both sites are dark mare: heat1d's validation uses A0 = 0.06 [A]. The Hayne 2017 Fig. A1 caption says "All model parameters are from Table A1", which includes A(θ) [A]; this contradicts the author's hypothesis that the figure used a constant albedo. Reviewer model, surface means A17 / A15: A(θ) with A0 = 0.06 at 1 AU gives 210.9 / 208.2 K; constant 0.08 gives 215.5 / 213.4 K. | Keep TH-V2 as recorded. Add a like-for-like case, registered before it is run: site albedo from the source, brightness-temperature output with ε = 1, the source's solar distance. Tolerance unchanged (±5 K). |
| F06 | Q3 | Neglected terrain radiation is compared with the wrong quantity. The 0.02 W m⁻² bound is "below 1 %" of the mean absorbed flux, but at the 43 K minimum the surface emits only 0.18 W m⁻². The bound is about 11 % of the night-time balance and equals the interior heat flow (0.018 W m⁻²) that the model does include. | **Major** | Emission εσT⁴ at 43 K = 0.184 W m⁻² (reviewer arithmetic). Reviewer model: an extra 0.02 / 0.1 / 0.3 / 1.0 W m⁻² raises the minimum by +1.0 / +4.4 / +11.0 / +24.8 K and the mean by +0.5 / +2.2 / +5.8 / +14.3 K. Earthshine can exceed 150 mW m⁻² locally (HLS-UG-001 §11.1) [A]. Scattering raises PSR maxima by more than 100 K, and 1-D errors exceed 100 K at high latitude in strong terrain (King 2020) [A]. | Report the bound relative to εσT_min⁴. Add single-bounce scattered sunlight, terrain IR and Earthshine (change 3). Until then, flag minima as unreliable whenever the bound exceeds about 5 % of εσT_min⁴. |
| F07 | Q3 | Temperature at one DEM pixel is ill-conditioned at the pole. Slope uncertainty of about 2° (the 5 m DEM's error) moves the maximum by tens of kelvin. Diviner-constrained facets around the test point are much warmer than the engine's level pixel. | **Major** | Reviewer model, fixed-aspect slope 1° / 2° / 5°: maximum +17.4 / +32.1 / +67.0 K (A(θ)). Barker 2021: median RMS slope error 1.90° [A]. PRP facets 0.14-0.38 km from the test point, slopes 2.9-6.7°: annual maximum 199-231 K [A, provisional extraction]. In the 89.67-89.71° S ring, facets with slope < 3° have a mean maximum of 119.9 K; 3-6° facets 182.0 K [A, provisional]. | Report temperature over a slope and aspect grid, or over the sub-pixel slope distribution, not at the nominal pixel only. Use the ensemble median rather than the nominal value (the nominal DEM is flatter than its clones). |
| F08 | Q4 | ILL-V2 probably compares different metrics. The PGDA map counts time steps with *any part* of the disc visible; a mean visible-disc fraction is lower by construction. This alone can explain the ~6 pp deficit on well-lit cells without a 1.5 m observer height. | **Major** (on the stated explanation); the ILL-V2 PASS is not affected | AVGVISIB label: values "indicate the percentage of timesteps where each pixel was sunlit (by any fraction of the solar disc)" [A, label read with TLS verification disabled]. Readmes describe illumination "on the lunar surface" [A]. Gläser 2018: solar visibility exceeds average illumination by 3.9-8.8 pp at the same points [A]. | Recompute ILL-V2 with the "any part visible" statistic at 0 m (see Q4). If the deficit on cells above 50 % stays near 6 pp, this hypothesis is wrong. |
| F09 | Q3/Q5 | Low-temperature conductivity is not represented and matters at depth. Heat capacity below 90 K does not matter. | **Major** for depth outputs; **Minor** for c_p | Reviewer model: Biele 2022 or Woods-Robinson c_p instead of the Hayne polynomial moves the minimum by -0.7 / -0.2 K and the mean by ≤ 0.3 K, although c_p differs by 51 % / 19 % at 43 K. A crude conductivity reduction below 80 K (×0.5 / ×0.2) moves the minimum by -2.7 / -5.6 K and the 1 m mean by +6.3 / +18.1 K. DSNE Rev J: below 100 K the Hayne relations "will likely over predict the thermal conductivity"; it specifies Woods-Robinson c_p for 10-400 K [A]. | Adopt the Woods-Robinson c_p (as DSNE does) and change the warning. Add a low-T conductivity option (Woods-Robinson 2019 / Martinez & Siegler 2021) and report the 0.3 m and 1 m results as unbounded until it exists. |
| F10 | Q2 | "Open, verified, traceable, with uncertainty" is mostly not new. The genuinely new elements are narrower than claimed. | Minor | DEM-clone uncertainty: Barker 2021 + open PGDA 78 clones [A]. Open thermal code with a Table A2 validation suite: heat1d (MIT) [A]. Terrain scattering and IR in open (unlicensed) codes: PCC Topo3D, fluxpy, O3DTM [A]. At least five 2026 open repositories with LOLA + SPICE illumination, some with Hayne thermal [A, README level]. | Rewrite the contribution claim (Q2). Cite prior art. Consider building terrain radiation on fluxpy or Topo3D rather than re-implementing it. |
| F11 | Q4 | TH-X1 used an obsolete reference. heat1d's README says the PyPI package (0.3.2, the version used) "is an obsolete release that is no longer updated". The findings "not converged at depth" and "NaN on finer grids" may not apply to the current code. | Minor | heat1d README, GitHub master [A, read directly by the reviewer]. | Re-run TH-X1 against the current GitHub heat1d, registered before the run. Keep the 0.3.2 result as history. |
| F12 | Q4 | Two tolerances (TH-V1, TH-X1) were written after exploratory runs. This is disclosed and anchored to source figures, but the 2 K tolerance is decisive for the TH-V1 minimum (1.26 K off; it would fail at 1 K). | Minor | [L: VALIDATION TH-V1, TH-X1]. Since TH-V1 does not discriminate (F04), little conclusion rests on it. | Commit tolerances for TH-V3, the ILL-V2 metric variant, Gläser 2018 and the like-for-like TH-V2 to git before running them. |
| F13 | Q4 | Pending illumination cases can be unblocked from an open source. Gläser 2018 (NTRS manuscript) gives coordinates and 2 m results for Shackleton and Connecting Ridge points; one is 21 m from the test point. | Minor (an opportunity) | Gläser 2018 Table 2 (2 m, 2017-2037): Shackleton S −89.6866, 197.19° E, 85.5 % average illumination, 89.4 % solar visibility, 66 h maximum shadow. Connecting Ridge C1 −89.4398, 222.83° E, 88.0 %, 112 h [A, read directly by the reviewer]. lunarenv test point = −89.6870, 197.30° E (reviewer conversion); 2026 mean 85.51 %, longest shadow 63 h [L: README]. | Add an ILL-V5 case (Gläser 2018 Tables 1-2, surface and 2 m) with a tolerance set beforehand that includes inter-annual variation (1 year against 20 years). |
| F14 | Q3 | The nested-DEM horizon is sound. | Minor (none needed) | Same recipe as Barker 2021 (5 m < 5 km, 80 m to 100 km, 240 m beyond) [A]. Convergence documented: numerical error ≤ +0.12 pp [L]. ILL-V1 within 2.2 pp [L]. | Keep. Add regional-DEM error (PGDA 90 `ADJ_ERR`) for distant skylines, since the reported spread is a lower bound there. |
| F15 | Q3 | Hayne 2017 property laws at 89° S: acceptable with stated limits. The latitude extrapolation itself matters less than the albedo law and slope. | Minor | Hayne 2017 fits to ±70°, with topographic artifacts above ~60° [A]. Bürger 2024: latitudinal gradient against the standard model up to ±80° [S]. Reviewer model: heat flow 0.012-0.030 W m⁻² changes the surface mean by ≤ 0.3 K and the 1 m mean by −2.0/+4.0 K. | Keep, with the warning. TH-V3 is the test. |
| F16 | Q3 | 1-D conduction, no rocks, no volatiles: acceptable for regolith at a sunlit rim. | Minor | Rock heating matters above ~3 % rock abundance; rocks larger than 10-100 cm cover < 1 % globally (Hayne 2017 App. B, citing Bandfield 2011) [A]. | Keep. Flag rocky pixels if a rock-abundance map is added. |
| F17 | Q5 | The ensemble quantiles omit the dominant terms: albedo model form, low-T conductivity, terrain radiation, sub-DEM slope. A reader of "5/50/95 %: 105.6 / 118.4 / 146.8 K" will take it as the uncertainty. | **Major** | [L: README, models/thermal "Not sampled"]. Reviewer model: the albedo envelope alone spans 101.5-156.7 K on the maximum at the nominal slope. | Present thermal output as scenario bounds (hot / cold envelope × ensemble 5-95 %), with the ensemble shown as the within-scenario spread. |
| F18 | Q4 | No high-latitude in-situ temperature measurement is comparable with an 89.7° S level-ground 1-D model. ChaSTE is the only peer-reviewed high-southern-latitude dataset. | Minor (establishes a limit) | ChaSTE: 69.373° S, daytime only, ~6° sunward slope, peak 355 K against Diviner 330 K (Durga Prasad 2025 CEE) [A]. Lander influence ~40 % (MNRAS 2025) [S]. Chang'e-7 slipped to 2027 [S]. IM-1 / IM-2: no regolith temperatures published [S]. | Use Diviner PCP for TH-V3. Treat ChaSTE at most as a moderate-incidence slope check. |

---

## 3. Answers to Q1-Q5

### Reviewer calculation (used in several answers)

To test the author's hypotheses without the source code, the reviewer wrote an independent reduced 1-D model from Hayne 2017 Appendix A and Table A1 [A]. It uses finite volumes on 43 cells (first cell 2 mm, growth 1.12, 2 m deep), backward Euler with Newton-linearized surface emission, and Table A1 parameters. Details and code are in Appendix A.

Checks before use:
- Equator: 385.2 / 100.1 / 93.5 K, against Hayne's 385 / 101 / 95 K and lunarenv's 385.2 / 100.3 / 93.7 K.
- Apollo surface means: 205.6 / 203.0 K, identical to lunarenv's TH-V2.
- Polar test-point forcing (synthetic ephemeris, flat horizon, level ground, no terrain shadowing): 42.7 / 113.6 / 69.0 K minimum / maximum / mean, against lunarenv's 43.3 / 119.2 / 70.3 K.

Its deep (≥ 0.8 m) Apollo means are ~2 K below lunarenv's and drift ~0.3 K with longer spin-up, so they are not used.

All reviewer-model numbers are sensitivities of a proxy, not lunarenv results. Each recommendation that relies on one says how lunarenv itself can confirm or refute it.

### Q1. Is the problem well posed for its users?

**Conclusion: partly.** Illumination (Sun and Earth vectors, visible disc fraction, shadow durations) is a quantity engineers need directly. Level-ground regolith temperature at a point is not. It is one boundary condition of a coupled asset problem. NASA's references present design environments as envelopes, and leave site-specific, asset-coupled, orientation-resolved analysis to the engineer.

**What the design references contain.**
- **DSNE.** The current revision is M2M-30159 Rev J (released 19 Dec 2025, 446 pp.), which supersedes SLS-SPEC-159 Rev I (27 Oct 2021) [A].
  - §3.4.6 says "Tabular values in this section are to be used for design. Figures are provided for illustration only." [A].
  - Surface temperatures are *extremes by latitude band* (Table 3.4.6.3-1, from Diviner zonal means of Williams 2017, ±1σ). Examples: 85° maximum 182 / 224 K (mean / 1σ) at local noon; 85° minimum 61 / 41 K; coldest PSR 18 K; hottest south-polar value 323 K, on an equator-facing slope at 84.6° S (new in Rev J) [A].
  - It has no site time series and nothing above 85° except the PSR and hottest-pole rows [A].
  - Rev J adds an incidence-dependent albedo envelope (Table 3.4.6.1.1-1). For 84-90° S it says emissivity and normal albedo are "Not available"; use equatorial highland values [A].
  - Illumination (§3.4.9) is descriptive only. It points to LROC WAC maps; it gives no Earthshine and no surface view factors [A].
  - §3.4.6.5 defers asset thermal analysis to HLS-UG-001 [A].
- **HLS-UG-001 (L-TAG, 2021)** prescribes a method rather than numbers [A]:
  - 1-D regolith conduction coupled to the asset (§9.2, the same model class as lunarenv);
  - ground-plane sizing by view factors (§9.3.1);
  - far-field terrain reduced to a horizon (§9.4.3);
  - a polar Earthshine model with tilt (§11.1, peaks above 150 mW m⁻²);
  - a cold-case ground boundary of 18-50 K, with crater floors near 50 K from wall-scattered light (§11.2).
- **ESA.** ECSS-E-ST-10-04C Rev.1 (2020) has no lunar-surface thermal or illumination content; the Moon is absent from its albedo/IR table (Table F-1) [A]. ESA's lunar-surface guidance appears to be an internal ESATAN-TMS technical note (ESA-TECMTV-TN-2022-003981), not public [S]. ESA engineers use a Hayne-2017-based engineering model (ICES-2025-202) [S]. For Luna-27/ProSPA the environment was specified as a single envelope flux, 707 W m⁻² including a 15° slope, plus a 20 W m⁻² IR margin (ICES-2020-291) [S]. No public Argonaut or EL3 environment specification was found.

**What an engineer needs first, and how lunarenv compares.**

| Need | Source | lunarenv |
|---|---|---|
| Design envelopes (hot/cold cases) | DSNE Table 3.4.6.3-1; L-TAG cold-case 18-50 K | Time series plus ensemble quantiles; no hot/cold case definition |
| Sun (and Earth) vectors at the site with terrain horizon | L-TAG §9.4.3 (horizon reduction) | **Yes**, validated (ILL-X1, ILL-V1, ILL-V2): this is its strongest contribution |
| Flux on arbitrarily oriented and vertical surfaces | L-TAG; hottest pole case on equator-facing slope (DSNE) | No (level or DEM-slope ground only) |
| Terrain IR and reflected sunlight; view factors | L-TAG §9.3, Birmingham 2024 | No (bound only) |
| Earthshine | L-TAG §11.1 | No |
| Ground temperature as coupled boundary | L-TAG §9.0 | Uncoupled ground temperature (stated in README) |
| Albedo and emissivity at the pole | DSNE: "Not available"; envelope from equatorial highlands | A single law (Hayne) near the cold end of the envelope (F03) |

So lunarenv fills the gap DSNE explicitly leaves: site-specific, DEM-based time series, for which DSNE points users to raw data. But its thermal output answers a narrower question than the one structure designers face (F02).

### Q2. Does it fill a real gap?

**Conclusion: a modest, real gap.** It is an integration and an auditable validation ledger, not new science. Most of the advertised ingredients exist.

**What exists.**
- **Illumination.**
  - PGDA 69 distributes open 18.6-year *average* maps of Sun, Earth and sky visibility at 60-240 m. It gives no time series, no stated observer height and no uncertainty layer, and the code is not released ("Contact us if you have specific needs") [A].
  - Barker 2021 published the exact recipe lunarenv adopts: Mazarico 2011 horizons, nested 5 / 80 / 240 m DEMs, limb-darkened disc, 100 DEM error clones and percentiles [A]. The clones are open (PGDA 78) [A].
  - Gläser 2018 gives 20 m results at surface and 2 m for named sites [A].
  - Mazarico 2011, Gläser 2014, Speyerer & Robinson 2013 (observed WAC illumination) and Bussey 2010 were read only as abstracts [S].
- **Thermal.**
  - heat1d (Hayne; MIT licence) implements exactly the Hayne 2017 laws, including A(θ), with an 8-check validation suite. It has no DEM horizon [A].
  - PCC Topo3D (Schorghofer), fluxpy (Potter) and O3DTM (King 2020) include terrain-scattered sunlight and terrain IR; none declares a licence [A].
  - O3DTM reduced its error against Diviner from ~10 K to ~2 K for mid-latitude craters [A].
  - The Paige 2010 polar model (ray-traced, with scattering) and Gläser & Gläser 2019 (scattering on a 180 m DTM) are not open [S].
  - No open polar code with a published site-specific Diviner polar validation was found.
- **Engineer-facing tools.** NASA uses DSNE + HLS-UG-001 with Thermal Desktop; ESA uses ESATAN-TMS toolboxes [A/S]. At least five open GitHub projects created in 2026 (e.g., ArtemisThermalBase, MIT; lunar-horizon, MIT) offer LOLA + SPICE illumination and Earth visibility, some with Hayne-based thermal and provenance sidecars. None reports a Diviner validation [A, README level; immature, 0-2 stars].

**What is new in lunarenv, specifically.**
1. Point time series at a user-defined height, with DEM-clone percentiles, from open code. PGDA offers this only on request; Barker 2021 published it for fixed regions.
2. Joint propagation of DEM clones and thermophysical parameters into point temperatures. Not found elsewhere.
3. A validation ledger with sources, metrics, tolerances, disclosed post-hoc tolerances and recorded FAIL and PENDING cases. Not found in any of the open codes examined.
4. Run hashing and full data provenance. Already claimed by two 2026 repositories.
5. An API and MCP surface.

Items 1-3 are a meaningful contribution *if* the validation ledger is extended to the polar thermal case. Without TH-V3, the thermal half repackages heat1d's physics with less terrain radiation than existing 3-D codes.

### Q3. Are the modelling choices scientifically defensible at the pole?

| # | Choice | Assessment | Key evidence |
|---|---|---|---|
| 1 | Nested-DEM horizon (5 / 80 / 240 m) | **Sound** | Same as Barker 2021 [A]; converged (≤ +0.12 pp) [L]; independent agreement with Gläser 2018 at the test point (F13) [A] |
| 2 | Hayne 2017 laws (fitted within ±70°) at 89° S | **Acceptable with stated limits** | Hayne 2017 §3.2 [A]; reviewer model: heat flow and c_p are minor; albedo and slope dominate |
| 3 | A(θ) at θ ≈ 88° | **Undermines the thermal results as presented**: neither A(θ) nor constant 0.12 is established; the alternative is outside the published envelope | See below |
| 4 | c_p fit below 90 K | **Sound in effect** (minor); low-T *conductivity* is the real gap | F09 |
| 5 | No terrain radiation | **Acceptable for maxima and means on a sunlit, near-level rim crest; undermines minima and any tilted or shadowed surface** | F06 |
| 6 | 1-D, no rocks or volatiles; slope from a 90 %-interpolated 5 m DEM | **1-D: acceptable. Point slope: undermines point predictions** unless presented over orientations | F07, F16 |

**3.3 Albedo at grazing incidence: what is established.**
- *Form and fit.* Hayne 2017 adopts A(θ) = A0 + a(θ/45°)³ + b(θ/90°)⁸ "following the empirical fits of Keihm et al. (1984) and Vasavada et al. (2012)", with a = 0.06, b = 0.25 fitted in "this study" [A]. The paper does not state the incidence range of the fit, and does not discuss polar use [A]. Vasavada 2012 used a = 0.045, b = 0.14, and noted that albedo depends on illumination, viewing geometry and roughness [S]. Keihm 1984: a = 0.03, b = 0.14 (as quoted by Foote 2020) [S].
- *Values at the test geometry.* At θ = 88.2° the three coefficient sets give A = 0.784 (Hayne), 0.578 (Vasavada form) and 0.465 (Keihm) (reviewer arithmetic). The absorbed fraction therefore varies 2.5-fold across the authors who built this law.
- *Laboratory evidence is split.* Birkebak (1974) integrating-sphere reflectances of Apollo fines, as reproduced in HLS-UG-001 Table 7.1-1, rise to 0.49-0.62 at 80° and 0.71-0.79 at 85°, close to Hayne's 0.555 and 0.683 [A via HLS-UG-001; Birkebak itself not read]. Modern goniometer data (Curtis 2024; Foote 2020) stop at 60° and are "considerably flatter than those assumed in previous thermal models" [A / S].
- *Interpretation (hypothesis, not established).* The Diviner-fitted law is an *apparent* albedo. Cast shadows in a rough nadir-viewed footprint lower the brightness temperature, which a smooth-surface fit absorbs as higher albedo (Nishiyama 2022 [A]; Vasavada 2012 [S]). Applied to a *resolved* 5 m facet with its own slope, it may double-count roughness. Real roughness at grazing incidence also adds sunward facets and multiple scattering, which raise absorbed flux (Rozitis & Green 2011) [A].
- *Polar studies are inconsistent.* Hayne 2021 (micro cold traps) uses a constant 0.12 [A]. King 2020 (O3DTM, south-polar sites) applies a Keihm-type A(θ) to direct and scattered light [A]. DSNE Rev J treats incidence dependence as the baseline, and calls constant albedo a hot bias: +28 K at 84° S against the Feng 2020 law [A]. The albedo choices of Paige 2010 and Williams 2019 could not be read [S].
- *Which is justified?* Neither is established at 88°. Constant 0.12 is less defensible than A(θ): it lies outside the DSNE envelope (0.238-0.860), and both lab reflectances (Birkebak) and Diviner fits show strong incidence dependence. A(θ) with Hayne's coefficients sits near the cold end.

  The defensible presentation is the DSNE Rev J highland envelope at the test point (reviewer model, level ground):
  - cold (Bürger max): minimum 40.8 K, maximum 101.5 K, mean 62.9 K;
  - Hayne: 42.7 / 113.6 / 69.0 K;
  - hot (Foote min): 49.2 / 156.7 / 91.9 K.

  Only a polar temperature comparison can narrow this (Q4).
- *The author's hypothesis about Hayne Fig. A1 is contradicted by its caption.* The author's claim that Hayne 2017 Fig. A1 "appears to correspond to a constant albedo" conflicts with the caption, "All model parameters are from Table A1, with the Moon at 1.017 AU" [A, read by the reviewer]. Either the figure's implementation differs from Eq. A8, or another input (site albedo) differs. heat1d's own Apollo check uses A0 = 0.06 [A] and documents 214.2 / 212.2 K; the reviewer's reduced model gives 210.9 / 208.2 K for that set-up at 1 AU. This discrepancy is unresolved, and the hypothesis is not established.

**3.4 Heat capacity.** Hayne 2017 states the fit is "valid for temperatures from < 90 K to >400 K, but is strictly invalid at temperatures < 1.3 K" [A]. Biele 2022 treats it as valid at 90-350 K [A].

At 43 K the polynomial gives 117.8 J kg⁻¹ K⁻¹. Alternatives give less:
- Biele 2022 Eq. 24: 78.1 J kg⁻¹ K⁻¹ (5-10 % uncertainty below 50 K) [A];
- Woods-Robinson 2019, as tabulated in DSNE Rev J: ~99 J kg⁻¹ K⁻¹ (log-interpolated between the 40 K and 50 K rows, 89 and 122) [A].

The reviewer model with either alternative changes the yearly minimum by -0.7 / -0.2 K and the mean by -0.3 / -0.1 K. The likely error from the extrapolation is therefore below 1 K at the surface. A better published set exists (Woods-Robinson 2019, adopted by DSNE Rev J for 10-400 K). The larger low-temperature issue is conductivity, which DSNE flags (F09).

**3.5 Terrain-scattered and emitted radiation.**
- *Established.* It dominates in shadow: "In large shadows … heating is dominated by radiation scattered and emitted by surrounding terrain" (Hayne 2021) [A]. In PSRs the effect exceeds 100 K (King 2020) [A]. Flat surfaces with no line of sight to other terrain are unaffected (King 2020) [A]. The best-lit Shackleton-rim and ridge sites receive essentially no extra flux (Gläser & Gläser 2019) [S].
- *Consequence for lunarenv.* For the sunlit, near-level rim test point, neglecting it is acceptable for the yearly maximum and mean (+0.1 K / +0.5 K at the reported 0.02 W m⁻² bound, reviewer model). It is not acceptable for the minimum, which is set by a night-time balance of 0.18 W m⁻². Nor is it acceptable for tilted surfaces or structures, which see the terrain directly. Crater rims are therefore not "meaningless", but the minima there are.

**3.6 Slope and 1-D.**
- *1-D and rocks.* 1-D conduction is standard for regolith (L-TAG §9.2 [A]). Rock heating matters only above ~3 % rock abundance [A].
- *Slope dominates.* The 5 m DEM's median RMS slope error is 1.90° (Barker 2021) [A]. In the reviewer model, 1° / 2° / 5° of fixed-aspect slope raise the maximum by +17 / +32 / +67 K.
- *Diviner-constrained facets near the point are much warmer.* PRP facets within 0.4 km of the test point, sloped 2.9-6.7°, have annual maxima of 199-231 K [A, provisional]. These are 0.5-0.7 km facets and model fits, so they neither confirm nor refute lunarenv's level-pixel 119 K. They do show that a single-pixel temperature is not representative of the ground a structure occupies.

### Q4. Is the validation strategy adequate to call the results "verified"?

**Conclusion.** Illumination: yes, it is verified and validated within its stated range. Thermal: verified as a solver (TH-A1, A2, P1), but not validated for the polar use; "verified" should not be applied to polar temperatures.

**Polar or high-latitude temperature data that exist.**
- **Diviner Polar Cumulative Products** (Williams; PDS Geosciences Node, bundle `urn-nasa-pds-lro_diviner_derived1`, collection `data_derived_pcp`).
  - Content: binned bolometric brightness temperature at 240 m, from 80° S to the pole, observations 5 Jul 2009 to 15 Mar 2020 [A].
  - Seasonal set: `seasonal/slon/pols/pcp_avg_tbol_pols_lsNNN_slonMM_240.tab`, 6 Ls bins × 24 subsolar-longitude bins. Diurnal summer set: `diurnal/ltim/pols/pcp_avg_tbol_pols_sum_ltimNN_240.tab` [A].
  - The label cites Schorghofer & Williams 2020 (PSJ 1, 54) [A]. Shackleton rim at 89.7° S is covered at 240 m [A].
  - Values at the test point were not extracted: the files are 190-220 MB each and unindexed. **This is the comparable dataset.**
- **Diviner Polar Resource Product** (Paige; `data_derived_prp/dlre_prp_south.tab`, LRO-L-DLRE-5-PRP-V2.0).
  - Content: model-fit annual mean temperature at 2 cm and annual surface maximum, on 0.5-0.7 km facets [A].
  - It is model output constrained by first-year Diviner data [S], so it is useful only as a cross-check.
- **Williams 2019** (seasonal polar maps): abstract only [S].
- **In situ.**
  - ChaSTE / Chandrayaan-3 (69.373° S, 32.319° E): peak surface 355 K on a ~6° sunward slope against Diviner 330 K; daytime only (9:42-15:00 LT); probe to 100 mm (Durga Prasad et al. 2025, Commun. Earth Environ. 6, 153) [A]. The abstract also states 335 K, an internal inconsistency [A]. In-situ conductivity 0.0115-0.0124 W m⁻¹ K⁻¹ at 80 mm (Mathew 2025) [A]. Lander influence estimated at ~40 % (MNRAS 2025) [S]. Not comparable as a direct validation: wrong latitude regime, no night, lander-perturbed, sloped.
  - Chang'e-7 (Shackleton rim target) slipped to 2027 [S]. IM-1 and IM-2 published no regolith temperatures [S]. Chang'e-4 (45° S, probes < 1 cm, 1-D fit) [A] and Blue Ghost LISTER (mid-latitude heat flow 13-14 mW m⁻²) [S] are not polar.
  - DSNE: "These results are not extensible to the lunar south pole where no measurements have been made." [A].

**Proposed TH-V3 (a check that can fail).**
- *Cells.* The 240 m PCP cells whose footprint contains the test point and its rim neighbours.
- *Bins.* Summer Ls bins and subsolar longitudes when the cells are sunlit.
- *Forward model.* Compute every 5 m pixel in each cell with its own slope and horizon, aggregate radiance (σT⁴ averaged, then inverted) into a Tbol, and do this for the hot and cold albedo bounds and for Hayne.
- *Why it discriminates.* At these cells scattering is small (Gläser 2019 [S]), and the albedo laws differ by about 55 K in the maximum (reviewer model). The comparison therefore discriminates by a wide margin.
- *Before running.* Fix the cell list, the bins and a tolerance from the Diviner calibration uncertainty plus the aggregation error.
- *How to read the outcome.* If Diviner falls outside both bounds, the model form is wrong beyond albedo (sub-pixel roughness, scattering), and the thermal model should be rethought. If it falls inside, the envelope narrows to the law that fits.

**Illumination agreement: good enough?** Yes, for preliminary design at ≥ 1-2 m height, with three caveats.
- *Mean illumination.* ILL-V1 agrees within 2.2 pp at 1 m and 0.33 pp at 5 m. Terrain-clone spread is ±1 pp at 2 m [L]. At 85 % illumination this is about 2.5 % relative, smaller than typical power margins.
- *Longest periods.* Longest-period statistics, which drive battery sizing, have no tolerance. They agree within 3.1 days (shadow) at 1 m, but the longest illumination period differs by up to 49 days [L]. They should be reported with their ensemble spread and not called validated.
- *The ILL-V2 deficit on well-lit cells has a simpler explanation* than the 1.5 m observer height (F08). Test it by recomputing lunarenv's fraction of time steps with *any part* of the disc visible, at 0 m, on the same 10 000 cells. Prediction: the bias on cells above 50 % falls from about −6 pp to within ±2 pp. If it stays near −6 pp, the metric is not the cause, and the observer-height or slope-sampling hypothesis remains open.

**Post-hoc tolerances.** This weakens the conclusions little, because both cases are anchored to source-stated figures and disclosed. TH-V1's minimum (−1.26 K) passes only because the tolerance is 2 K rather than 1 K. Since TH-V1 does not discriminate between the albedo laws (F04), nothing about the polar result depends on it.

**Pending cases: what was accessed and what it changes.**

| Item | Access | What it changes |
|---|---|---|
| Mazarico 2011 Table 3 | Not accessed (paywall). Its Shackleton and CR1 values are listed in Gläser 2018 Table 1: S −89.78, 204.27° E, 86.7 % average, 90.5 % visibility; CR1 −89.45, 222.69° E, 89.0 % [A for Gläser's listing, S for Mazarico] | Gläser 2018 says Mazarico 2011 overestimates illumination because its 240 m grid smooths close topography [A], so it is a weak benchmark. Use Gläser 2018 instead. Mazarico 2011's time step is stated as 6 h in a snippet and 1 h on the PGDA page [S]; unresolved. |
| Gläser 2014 best-lit ridge point | Not accessed. Snippet gives longitude 222.8445° E only [S]; Shackleton SR2 at 89.6871° S, 197.0983° E [S], 33 m from the test point | ILL-V4 can be replaced by Gläser 2018 Table 2 C1 (−89.4398, 222.83° E; 88.0 % at 2 m; 112 h) [A] |
| Vasavada 2012 | Not accessed (paywall) | TH-V1 tolerances remain unsourced; coefficients a = 0.045, b = 0.14 known only from snippets [S] |
| Keihm & Langseth 1973 (Apollo 17) | Accessed [A], via OCR of the ADS scan | 216 K is a blackbody brightness temperature from cable thermocouples ~15 cm above the surface, with A = 0.08 and ε = 1.0 assumed. Daytime values carry large errors. Redefine the comparison (F05). |
| Keihm et al. 1973 (Apollo 15, EPSL) | Not accessed | Apollo 15 is presumably measured the same way; not verified |
| Langseth et al. 1976 | Accessed [A], via OCR | Does not re-derive the mean surface temperature. Annual surface-mean variation ~3 K peak-to-peak [A] |
| Hayne 2017 Table A2 | Read [A] | The Apollo citation labels (1973a / 1973b) appear swapped relative to the text and reference list [A]; the values are correct for Apollo 17 |

### Q5. Is the uncertainty honest?

**Conclusion.** It is honest in what it says (it lists what is not sampled) but misleading in what it shows. The headline quantiles are not the uncertainty of the result. The thermal output should be presented as a **bounded range from explicit scenarios** (albedo envelope × slope / orientation), with the ensemble as within-scenario spread. It should not be presented as a prediction with quantiles.

**Evidence.**
1. The sampled ensemble gives a surface-maximum 5-95 % range of 100-142 K [L]. The unsampled albedo model form alone spans 101.5-156.7 K at the nominal slope (reviewer model, DSNE envelope).
2. Slope, which the ensemble does sample through the clones, contributes tens of kelvin (F07).
3. Low-T conductivity (unsampled) moves the 1 m mean by +6 to +18 K in a crude bound (F09).
4. Terrain flux (unsampled) moves the minimum by +1 to +25 K for 0.02-1 W m⁻² (F06).

By contrast, the two assumed ranges matter little at the surface. Heat flow over 0.012-0.030 W m⁻² changes the surface mean by ≤ 0.3 K (reviewer model). That is consistent with the author's observation that heat flow's third-rank correlation is chance pairing. LISTER's 13-14 mW m⁻² [S] lies inside the range.

The non-varied χ is also secondary to the albedo law and slope.

Suggested presentation, per output:
- the **hot scenario** (DSNE highland min albedo) and the **cold scenario** (DSNE highland max albedo), each with its ensemble 5-95 %;
- the Hayne law as a named middle case;
- terrain radiation shown as a one-sided correction for minima;
- separate entries for a level surface and for user-specified orientations.

**What would show this recommendation wrong.** TH-V3 placing Diviner Tbol consistently at one law with residuals smaller than the ensemble spread. In that case the scenario axis collapses, and quantiles become the honest presentation.

### Q6. Supporting detail for the verdict

**Change 1: TH-V3 against Diviner PCP** (protocol in Q4).
- *Shows it worked:* a pre-registered comparison over a set of sunlit rim cells and summer bins, with bias and RMS reported per albedo law, that selects or excludes laws.
- *Shows it wrong:* Diviner outside both bounds, which points to roughness or scattering, not albedo.

**Change 2: bounded, orientation-resolved thermal product** (DSNE Rev J highland envelope × terrain/parameter ensemble × surface orientation grid, including vertical).
- *Shows it worked:* the TH-V3 cells fall inside the envelope, and the hottest PRP / PCP facets near the site fall inside the envelope at their own slope.
- *Shows it wrong:* Diviner warmer than the hot bound on sunlit equator-facing slopes. DSNE's 323 K hottest value at 84.6° S is a reference point for this check.

**Change 3: terrain-scattered sunlight, terrain IR and Earthshine (single bounce is enough), fluxes on oriented surfaces, and correction of the two non-like-for-like comparisons.**
- *Shows it worked:* (a) the bowl-crater analytic case (Ingersoll-type irradiance, as used by Hayne 2021) is reproduced; (b) the reported flux at the test-point minimum agrees with the 0.02 W m⁻² bound; (c) ILL-V2 with the "any part visible" metric removes the well-lit deficit; (d) like-for-like TH-V2 is reported against brightness temperature.
- *Shows it wrong:* the single-bounce flux at the test-point minimum exceeds about 0.05 W m⁻² (> 25 % of εσT_min⁴). The rim minimum would then be wrong by more than about 2 K, and the 0.02 W m⁻² bound would be shown to underestimate the neglected flux.

**Minimum before an engineer uses a lunarenv temperature in a design decision:**
1. TH-V3 passed at or adjacent to the site against a pre-registered tolerance.
2. The temperature is taken from the bounded envelope, not from the nominal case.
3. It is taken for the actual orientation of the design surface.
4. For any shadowed, tilted or crater-facing location, terrain radiation is included.
5. It is used as a boundary condition in a coupled asset model, per HLS-UG-001 §9, never as the structure's temperature.

Illumination outputs (Sun and Earth vectors, visible fraction, shadow durations with their clone spread) can be used for preliminary design now.

---

## 4. Sources consulted

### Read for this review

| Source | Mark | Used for |
|---|---|---|
| lunarenv README, PLAN, VALIDATION, docs/models/illumination.md, docs/models/thermal.md, docs/convergence/*, docs/references.md | [L] | Everything attributed to the project |
| Hayne P.O. et al. (2017) JGR Planets 122, 2371, doi:10.1002/2017JE005387 (arXiv:1711.00977), §3.2-3.3, App. A, B, Tables A1-A2, Fig. A1 caption | [A] (read by the reviewer) | Laws, coefficients, c_p validity statement, Fig. A1 caption, latitude limits, rock heating |
| NASA M2M-30159 DSNE Rev J (19 Dec 2025), §3.3.9.1, §3.4.6-3.4.6.5, §3.4.9, Tables 3.4.6.1.1-1/-2, 3.4.6.2-1, 3.4.6.3-1; SLS-SPEC-159 Rev I (27 Oct 2021) | [A] | Design-spec content (Q1); albedo envelope (F03); c_p and conductivity statements (F09) |
| NASA HLS-UG-001, Lunar Thermal Analysis Guidebook (L-TAG), 4 Jan 2021, §4.3.1, §7.1 (Table 7.1-1), §9.0-9.4, §11.1-11.2 | [A] | Engineering method; Birkebak lab reflectance; Earthshine; cold case |
| Moore R., DSNE overview, SESWW June 2026 (NTRS 20260005205) | [A] | No polar subsurface measurements |
| ECSS-E-ST-10-04C Rev.1 (15 June 2020) | [A] | No lunar-surface thermal content |
| Gasbarre J.F. et al. (2008) Ground plane and near-surface thermal analysis, TFAWS (NTRS 20080034507) | [A] | Asset-ground coupling |
| Birmingham W., Schunk G., Evans B. (2024) Lunar latitude and terrain radiator sensitivity study, TFAWS (NTRS 20240010276) | [A] | Terrain IR effect on radiators |
| Burbridge E. et al. (2026) ICES-2026-490 (LEMS Artemis III TVAC) | [A] | Example of envelope-based polar design |
| Barker M.K. et al. (2021) PSS 203, 105119 (NTRS manuscript) | [A] | Method, slope error, clones |
| Gläser P. et al. (2018) PSS 162, 170 (NTRS 20170007365), Tables 1-2 | [A] (Table 2 checked directly by the reviewer) | F13; Mazarico 2011 values as listed; metric difference (F08) |
| PGDA products 69, 78, 81, 90, 104, 106 web pages | [A] | Prior art (Q2) |
| PDS labels: DLRE PRP (`dlre_prp_south.lbl` / `.xml`), PCP (`pcp_avg_tbol_pols_ls030_slon01_240.xml`), LOLA `avgvisib_85s_060m_201608.lbl` and illumination readmes | [A]; the review's research agent read the PDS server with TLS certificate verification disabled (the sandbox could not verify its chain), which the reviewer could not repeat with verification on | Diviner archive paths (Q4); AVGVISIB metric (F08) |
| DLRE_PRP_SOUTH.TAB (byte-range subset of 9,384 facets near the pole) and AVGVISIB IMG pixel reads | [A], *provisional* (own geometry; ±1 pixel registration sensitivity) | F07 context; not used for any pass/fail statement |
| heat1d GitHub README and `python/docs/validation.md` (master, 2026) | [A] (checked directly by the reviewer) | Obsolete PyPI release (F11); Apollo checks with A0 = 0.06 (F05) |
| Schorghofer Planetary-Code-Collection; Potter fluxpy; King O3DTM; 2026 GitHub repositories (READMEs, licences) | [A] | Prior art (Q2) |
| King O. et al. (2020) PSS 182, 104790 (accepted manuscript) | [A] | Polar 3-D thermal model, scattering magnitudes, A(θ) use |
| Hayne P.O., Aharonson O., Schorghofer N. (2021) Nat. Astron. 5, 169 (arXiv:2005.05369) | [A] | Constant 0.12 albedo in polar rough-surface model; shadow heating |
| Biele J. et al. (2022) Int. J. Thermophys. 43, 144, §5.2.2, Eq. 24 | [A] | Low-T c_p (F09) |
| Curtis R.J. et al. (2024) Meteorit. Planet. Sci., doi:10.1111/maps.14266 (accepted manuscript) | [A] | Lab hemispheric albedo to 60° |
| Rozitis B., Green S.F. (2011) MNRAS 415, 2042 (arXiv:1211.1844) | [A] | Roughness lowers effective Bond albedo |
| Nishiyama G. et al. (2022) Earth Planets Space 74, 105 | [A] | Diviner-fitted albedo includes roughness |
| Gomez C. et al. (2022) L-TAG short course, TFAWS (NTRS 20220013412) | [A] | Context on c_p fits |
| Keihm S.J., Langseth M.G. (1973) Proc. LSC 4th, 2503 (ADS scan, OCR) | [A] (OCR; possible character errors) | What the Apollo 17 mean surface temperature is (F05) |
| Langseth M.G., Keihm S.J., Peters K. (1976) Proc. LSC 7th, 3143 (ADS scan, OCR) | [A] (OCR) | Pending-reference status |
| Durga Prasad K. et al. (2025) Commun. Earth Environ. 6, 153, doi:10.1038/s43247-025-02114-6 | [A] | ChaSTE (Q4) |
| Mathew N. et al. (2025) Sci. Rep., doi:10.1038/s41598-025-91866-4 | [A] | ChaSTE conductivity |
| Xiao X. et al. (2022) Natl. Sci. Rev. 9, nwac175 | [A] | Chang'e-4 in-situ temperatures |
| Vasavada A.R. et al. (2012) JGR 117, E00H18 | [S] | Albedo coefficients (snippet) |
| Keihm S.J. (1984) Icarus 60, 568 | [S] (via Foote 2020 and HLS-UG-001) | Albedo-law origin |
| Foote E.J. et al. (2020) Icarus 336, 113456 | [S] | Lab albedo flatness; Keihm coefficients |
| Williams J.-P. et al. (2017) Icarus 283, 300; (2019) JGR Planets 124, 2505 | [S] | Polar maps context |
| Bürger J. et al. (2024) JGR Planets 129, e2023JE008152 | [S] | Latitudinal gradient of the standard model |
| Feng J. et al. (2020) JGR Planets 125, e2019JE006130 | [S] | Alternative albedo law |
| Woods-Robinson R. et al. (2019) JGR Planets 124, 1989; Martinez A., Siegler M.A. (2021) JGR Planets 126, e2021JE006829 | [S] (Woods-Robinson c_p values via DSNE Rev J [A]) | Low-T properties |
| Gläser P., Gläser D. (2019) A&A 627, A129 | [S] | Scattering negligible at best-lit rim sites |
| Rubanenko L., Aharonson O. (2017) Icarus 296, 99 | [S] | Rough-terrain cold traps |
| Mazarico E. et al. (2011) Icarus 211, 1066; Mazarico et al. (2018) ASR 62, 3214; Gläser et al. (2014) Icarus 243, 78; Speyerer & Robinson (2013); Bussey et al. (2010); Barker et al. (2023) PSJ 4, 183 | [S] | Prior art; pending cases |
| Durga Prasad K. et al. (2025) MNRAS 538, 2330 | [S] | ChaSTE lander influence |
| Hager & Vaughan ICES-2025-202; ICES-2020-291 (ProSPA); ESA-TECMTV-TN-2022-003981 (cited) | [S] | ESA practice |
| SpaceNews (Aug 2026) on Chang'e-7; NASA IM-2 releases; LISTER EGU26-6066 abstract | [S] | Mission status |

### Tried and not accessed

| Source | What was tried | Result |
|---|---|---|
| Mazarico et al. 2011, Table 3 | Full-text retrieval by DOI; NTRS record (no PDF) | Paywalled; not read |
| Gläser et al. 2014 | Full-text retrieval by DOI; repository DOI | Paywalled; coordinates of the best-lit ridge point not obtained |
| Vasavada et al. 2012 | Full-text retrieval; NTRS search | HTTP 403; stated uncertainties not obtained |
| Williams et al. 2019; Paige et al. 2010; Woods-Robinson et al. 2019; Martinez & Siegler 2021; Bürger et al. 2024; Bandfield et al. 2015; Gläser & Gläser 2019; Gläser et al. 2021 | Full-text retrieval; author-hosted copies on non-allowlisted hosts (not requested) | Abstracts and snippets only |
| Keihm et al. 1973 (Apollo 15, EPSL 19, 337) | Full-text retrieval | Closed; Apollo 15 sensor details not verified |
| Hayne 2017 Fig. A1 plotted values | — | Caption read; curve values not digitized |
| Diviner PCP values at the test point | Directory listing and labels read | Not extracted (190-220 MB unindexed files) |
| ESA ICES papers and ESA technical note | Institutional repository blocked; no public copy | Snippets only |
| heat1d current version, run locally | Installed the PyPI package (the obsolete one); its import requires a download from naif.jpl.nasa.gov, which the sandbox blocks | Not run. The reviewer wrote an independent reduced model instead (Appendix A) |

---

## Appendix A. Reviewer's reduced model and runs

**Purpose.** To test hypotheses quantitatively without lunarenv's code. These are not lunarenv results.

**Model.**
- *Equations and parameters.* 1-D conduction ρ c_p(T) ∂T/∂t = ∂/∂z (K ∂T/∂z), using the laws and constants of Hayne 2017 Table A1: H = 0.06 m, ρ 1100-1800 kg m⁻³, K_s 7.4 × 10⁻⁴, K_d 3.4 × 10⁻³ W m⁻¹ K⁻¹, χ = 2.7, ε = 0.95, Q = 0.018 W m⁻², c_p quartic.
- *Grid.* Finite volumes, first cell 2 mm, growth 1.12, 2 m deep (43 cells).
- *Time stepping.* Backward Euler, with two Newton iterations on the surface balance per step and lagged coefficients.
- *Forcing.* Equatorial and Apollo cases: Sun on the equatorial plane at 1 or 1.017 AU, as stated; spin-up 300-600 lunations at P/96, then 6 lunations at P/2880. Polar case: latitude −89.69°; subsolar latitude 1.543° sin(2πt / 346.62 d); synodic hour angle; uniform-disc segment above a flat 0° horizon; no terrain shadowing; level ground unless stated; 60 years of spin-up at 6 h, then one year at 30 min. A 120-year spin-up changed the results by ≤ 0.1 K.
- *Checks.* Equator 385.2 / 100.1 / 93.5 K. Apollo surface means 205.6 / 203.0 K, the same as lunarenv's TH-V2. Polar nominal 42.7 / 113.6 / 69.0 K, against lunarenv's 43.3 / 119.2 / 70.3 K. Deep Apollo means are not reliable to better than ~2 K and are not used.

**A1. Equator (TH-V1 discrimination).**

| Case | Noon | Midnight | Minimum | Mean |
|---|---|---|---|---|
| A(θ), A0 = 0.12 | 385.2 | 100.1 | 93.5 | 210.8 |
| constant 0.12 | 385.2 | 101.7 | 94.9 | 218.0 |
| Reference (Hayne 2017 Table A2) | 385 | 101 | 95 | — |

**A2. Apollo diurnal mean surface temperature (K).** Targets: 216 / 211 K brightness temperature; ~218.8 / 213.7 K physical-temperature equivalent at ε = 0.95.

| Case | A17 (20° N) | A15 (26° N) |
|---|---|---|
| A(θ), A0 = 0.12, 1.017 AU | 205.6 | 203.0 |
| A(θ), A0 = 0.12, 1 AU | 207.2 | 204.5 |
| A(θ), A0 = 0.12, H = 0.07 | 205.1 | 202.5 |
| A(θ), A0 = 0.10 / 0.08 / 0.07 | 206.9 / 208.1 / 208.7 | 204.2 / 205.5 / 206.1 |
| A(θ), A0 = 0.06, 1.017 AU | 209.3 | 206.7 |
| A(θ), A0 = 0.06, 1 AU | 210.9 | 208.2 |
| constant 0.12 | 213.4 | 211.3 |
| constant 0.08 (Keihm's assumed A) | 215.5 | 213.4 |
| constant 0.06 | 216.6 | 214.4 |

**A3. Polar test-point sensitivities (level ground unless stated; K).**

| Case | Min | Max | Mean | Mean at 1 m |
|---|---|---|---|---|
| Nominal (Hayne A(θ)) | 42.7 | 113.6 | 69.0 | 74.8 |
| constant albedo 0.12 | 49.9 | 162.6 | 94.8 | 101.5 |
| DSNE highland min albedo (Foote 2020), hot bound | 49.2 | 156.7 | 91.9 | 98.5 |
| DSNE highland max albedo (Bürger 2024), cold bound | 40.8 | 101.5 | 62.9 | 68.6 |
| + 0.02 W m⁻² isotropic | 43.7 | 113.6 | 69.5 | 75.3 |
| + 0.1 W m⁻² | 47.2 | 113.9 | 71.2 | 77.0 |
| + 0.3 W m⁻² | 53.8 | 114.6 | 74.8 | 80.6 |
| + 1.0 W m⁻² | 67.5 | 116.9 | 83.3 | 89.0 |
| heat flow 0.030 / 0.012 W m⁻² | 43.4 / 42.4 | 113.6 / 113.5 | 69.3 / 68.9 | 78.8 / 72.8 |
| slope 1° / 2° / 5° (fixed aspect), A(θ) | 43.3 / 44.3 / 47.5 | 130.9 / 145.7 / 180.6 | 70.1 / 73.0 / 83.0 | 76.2 / 79.6 / 91.1 |
| slope 1° / 2° / 5°, constant 0.12 | 50.0 / 50.5 / 53.0 | 180.9 / 195.0 / 239.9 | 93.5 / 94.3 / 102.2 | 101.4 / 103.2 / 113.7 |
| c_p Biele 2022 Eq. 24 | 42.1 | 113.6 | 68.8 | 74.6 |
| c_p Woods-Robinson (DSNE Rev J table) | 42.6 | 113.6 | 69.0 | 74.7 |
| K × 0.5 below 80 K (crude bound, not a published law) | 40.1 | 113.8 | 67.6 | 81.1 |
| K × 0.2 below 80 K (crude bound) | 37.1 | 114.1 | 65.9 | 92.9 |

Flux arithmetic at 1.8° Sun elevation:
- Level ground: 1361 sin 1.8° = 42.8 W m⁻² incident. Absorbed 9.2 W m⁻² with A(88.2°) = 0.784, or 37.6 W m⁻² with A = 0.12.
- Vertical sun-facing face: 1360 W m⁻².
- Emission at 43 K: 0.184 W m⁻².

**Test-point coordinates.** (x, y) = (−2822.5, −9062.5) m, south polar stereographic, R = 1737.4 km, longitude = atan2(x, y). This gives −89.6870°, 197.30° E. It lies 21 m from Gläser 2018's 2 m Shackleton point (−89.6866°, 197.19° E) and 33 m from Gläser 2014 SR2 (−89.6871°, 197.0983° E, snippet).

**Core code (Python / NumPy / SciPy):**

```python
import numpy as np
from scipy.linalg import solve_banded
SIG=5.670374e-8; S0=1361.0; P=2.55024e6
c=np.array([-3.6125,2.7431,2.3616e-3,-1.2340e-5,8.9093e-9])
cp_hayne=lambda T: c[0]+c[1]*T+c[2]*T**2+c[3]*T**3+c[4]*T**4
A_hayne=lambda th,A0=0.12: A0+0.06*(th/(np.pi/4))**3+0.25*(th/(np.pi/2))**8
class Model:
    def __init__(s,H=0.06,rs=1100.,rd=1800.,Ks=7.4e-4,Kd=3.4e-3,chi=2.7,eps=0.95,Q=0.018,
                 dz0=2e-3,g=1.12,depth=2.0,cp=cp_hayne):
        dz=[dz0]
        while sum(dz)<depth: dz.append(dz[-1]*g)
        s.dz=np.array(dz); s.z=np.cumsum(s.dz)-s.dz/2
        s.rho=rd-(rd-rs)*np.exp(-s.z/H); s.Kc=Kd-(Kd-Ks)*(rd-s.rho)/(rd-rs)
        s.chi,s.eps,s.Q,s.cp,s.n=chi,eps,Q,cp,len(s.dz)
    def K(s,T): return s.Kc*(1+s.chi*(T/350.)**3)
    def run(s,flux,t0,t1,dt,T,Ts):
        n=s.n; dist=np.concatenate([[s.dz[0]/2],(s.dz[:-1]+s.dz[1:])/2]); t=t0
        for _ in range(int(round((t1-t0)/dt))):
            t+=dt; F=flux(t); Kc=s.K(T)
            Ki=np.empty(n); Ki[0]=Kc[0]; Ki[1:]=0.5*(Kc[:-1]+Kc[1:]); G=Ki/dist
            C=s.rho*s.cp(T)*s.dz/dt; Gd=np.append(G[1:],0.0)
            for _ in range(2):                      # Newton on eps*sig*Ts^4
                ab=np.zeros((3,n+1)); r=np.zeros(n+1); d=4*s.eps*SIG*Ts**3
                ab[1,0]=d+G[0]; ab[0,1]=-G[0]; r[0]=F+3*s.eps*SIG*Ts**4
                ab[1,1:]=C+G+Gd; ab[2,:n]=-G; ab[0,2:]=-Gd[:-1]
                r[1:]=C*T; r[-1]+=s.Q
                x=solve_banded((1,1),ab,r); Ts=x[0]
            T=x[1:].copy()
        return T,Ts
```

The polar forcing, the albedo-envelope laws (DSNE Rev J Table 3.4.6.1.1-1, highland, θ > 60° and θ > 55° branches) and the c_p alternatives (Biele 2022 Eq. 24; DSNE Rev J Woods-Robinson table, log-interpolated) are as described above.
