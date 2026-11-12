# lunarenv — open, verified lunar environment engine: MVP plan

Status: plan approved 2026-11-12. **All seven phases complete; the MVP meets its definition of done** (section 13). Validation stands at 11 PASS, 1 FAIL, 1 PARTIAL, 3 PENDING after ILL-V1 was re-run with the full clone ensemble (addendum to Phase 4).
Reference keys (e.g. `Hayne2017`) resolve in [docs/references.md](docs/references.md), which also records which sources were actually read.

## 1. Goal and scope

Compute, deterministically and traceably, the illumination and thermal environment of a structure at a point on the lunar south pole over a period of time, with verification against published data. Test site: Shackleton crater rim.

In scope: data layer, illumination, 1-D thermal, verification suite, scenarios, Python API, CLI, MCP server, automatic report.
Out of scope (architecture must admit them later): structural/FEA, seismic, radiation, dust, ejecta, web UI, ML surrogates, multi-site campaigns.

## 2. Findings from research that shape the plan

1. **Best near-field DEM is the PGDA 5 m/pix site product** (`PGDA-78`, `Barker2021`). `Site04` is the Shackleton rim. Each site ships a height map, a formal height-error map, a slope-error map and **100 error "clones"** (statistically consistent DEM realizations). This gives a real, published input-uncertainty ensemble for Monte Carlo, not an invented one. Median RMS height error 0.30–0.50 m, slope error 1.5–2.5° (`Barker2021`).
2. **All LOLA polar DEMs are in the `MOON_ME` frame of DE421**, south polar stereographic, metres, pixel-registered (`PGDA-78` README). The reference radius is not stated in the README; it will be read from the GeoTIFF CRS and asserted (expected 1737.4 km, as stated for `PGDA-69`).
3. **`Barker2021` Table 2 is the strongest illumination benchmark available**, but it is for **Site 1 (connecting ridge), not Site 4**: 7 regions of interest with coordinates, average illumination, longest continuous illumination/shadow periods, at 1 m and 5 m above the surface, 2024-01-01 to 2026-01-01, hourly, as 1st/99th percentiles over the 100 clones. It uses the same DEM we can download, so differences are attributable to the algorithm rather than the terrain.
4. **`Mazarico2011` site table was not accessible** (paywalled). Only the abstract statement is available (a point near the Shackleton rim: ~240 days of continuous sunlight, ~1.5 days longest total darkness; observer height not confirmed). The published 18.6-year average-illumination map derived with that method **is** public (`PGDA-69`, 60 m/pix, 85–90°S) and can be compared pixel-wise.
5. **`Hayne2017` is fully specified** (Appendix A, Table A1 read): density profile, conductivity law, c_p polynomial, albedo law, boundary conditions, grid recipe and benchmark values. Two limits matter for a polar product:
   - the Diviner fit covers **latitudes within ±70°**; using it at 89.8°S is an extrapolation;
   - the c_p polynomial is stated valid from below 90 K to above 400 K; shadowed polar surfaces go colder.
   Both become explicit out-of-validated-range warnings in every result.
6. **Published Diviner-derived temperatures that a terrain-free 1-D model can reproduce are equatorial** (`Hayne2017` Table A2 quoting `Vasavada2012`: ~385 K noon, 101 K midnight, 95 K pre-dawn at the equator; Apollo 15/17 mean surface and subsurface temperatures ±5 K). Polar Diviner maps exist (`Williams2019`, `DIV-PCP`) but at 240 m/pix and strongly affected by terrain-scattered and re-emitted radiation, which the MVP neglects.
7. **Brute-force horizon maps for the full tile at 5 m are not a laptop job**: the Site 4 tile is 16 × 16 km (3200² ≈ 10 M pixels; the 20 km figure in the first draft was wrong) × 720 azimuths × ~1500 samples ≈ 10¹³ DEM samples (several hours on 8 cores, estimate). Maps must target a region of interest or a coarser grid; the far-field horizon varies slowly (`Barker2021`) and can be computed on a sparse observer grid and interpolated, with measured error.
8. Local machine has **24 GB free disk**. Data budget must stay small: crop on ingest, stream the clones.

## 3. Design decisions

Each row compares alternatives. "Verify" names the study that confirms or overturns the choice.

| # | Decision | Alternatives considered | Choice and reason | Verify |
|---|---|---|---|---|
| D1 | Shadowing method | (a) precomputed horizon map per observer; (b) per-epoch ray tracing against a mesh (BVH); (c) per-epoch shadow sweep on the grid | **(a)**. Cost is paid once; any epoch is then an O(1) lookup, which suits 1-year and 18.6-year series and UQ. It is the method of `Mazarico2011`/`Barker2021`, so validation is like-for-like. (b) is better for arbitrary receivers and multiple scattering; noted as the route for a future radiative module. | Convergence in azimuth and range step |
| D2 | Ray geometry | (a) straight lines in the stereographic plane with a curvature correction; (b) great-circle marching, elevation computed from 3-D body-fixed vectors | **(b)**. Exact on the sphere, no small-angle correction to justify; cost difference is negligible. (a) kept as a test oracle near the pole. | Agreement (a) vs (b); analytic sphere dip |
| D3 | DEM pyramid for horizons | (a) single resolution; (b) nested grids: 5 m near, coarser mid, coarsest far | **(b)**, as `Barker2021` (5 m < 5 km, 80 m to 100 km, 240 m beyond). Switch radii are parameters chosen by a convergence study, not copied. | Switch-radius study |
| D4 | DEM extent ("margin") | (a) fixed radius; (b) analytic visibility bound | **(b)**. A ray stops when a hypothetical point at the regional maximum height at that range could no longer rise above the horizon already found (exact spherical geometry). The largest stopping range over all observers is reported and must lie inside the far-field DEM (80°S product reaches ≥ 297 km from the site; the 60°S product is the fallback). | Recorded per run; test that extending the DEM changes nothing |
| D5 | Solar disc | (a) point Sun; (b) disc cut by a horizontal chord at the horizon elevation under the disc centre (`Mazarico2011`, `Barker2021`); (c) disc integrated in vertical strips against the horizon profile across its width | **(c)** as default, **(b)** as a selectable mode to reproduce published numbers. Both the geometric visible fraction and the limb-darkened flux fraction are output; the thermal model uses the flux fraction. | Strip-count convergence; (b) vs (c) difference reported |
| D6 | Ephemeris and frames | (a) DE421 SPK + `moon_pa_de421` (identical to LOLA processing, ends 2050); (b) DE440s + `moon_pa_de440` with the frame `MOON_ME_DE440_ME421` (DE440 orientation aligned to the DE421 mean-Earth frame) | **(b)**. Current ephemeris, coverage to 2150, and the frame was defined by NAIF for exactly this consistency problem. | Test: Sun direction (a) vs (b) at the site; cross-check against `HORIZONS` |
| D7 | Thermal physics | (a) `Hayne2017`; (b) `Vasavada2012` (two-layer predecessor); (c) `WoodsRobinson2019` low-temperature properties | **(a)** as baseline: current community standard, fully published, reference code exists. Property laws are pluggable so (c) can be added for cold surfaces; decision on shipping (c) in the MVP after reading the paper in Phase 3. | Benchmarks of `Hayne2017` Table A2 |
| D8 | Spatial discretization | (a) node-centred finite differences as in `Hayne2017`; (b) conservative finite volumes with interface conductivities | **(b)**. Discrete energy conservation holds to round-off and becomes a property test. (a) is implemented too, only as the reproduction of the published scheme. | Order of accuracy on analytic cases |
| D9 | Time integration | (a) explicit forward Euler (`Hayne2017`), Δt limited by the top cell; (b) implicit, L-stable, 2nd order (TR-BDF2 or BDF2) with Newton on the nonlinear surface balance; (c) Crank–Nicolson | Implement (a) and (b); default chosen by **measured error at equal cost**. Expectation: (b), because polar sunrise/sunset is abrupt and (c) is prone to oscillation with a stiff radiative boundary. | Temporal convergence, cost vs error plot |
| D10 | Uncertainty propagation | (a) local sensitivities (first-order); (b) Monte Carlo with Latin hypercube sampling plus the PGDA DEM clones | **(b)** for results (illumination is strongly nonlinear in terrain height: `Barker2021` shows errors bias illumination downward); **(a)** reported alongside as a sensitivity ranking. | Convergence of percentiles with ensemble size |
| D11 | Result format | (a) NetCDF-4 via xarray with CF attributes; (b) Parquet; (c) Zarr | **(a)**. N-dimensional (time, depth, y, x, quantile) data with units and metadata in one self-describing file; readable from every ecosystem engineers use. Parquet fits tables only. Zarr reconsidered if results outgrow a single file. | — |
| D12 | Units | (a) `pint` quantities end to end; (b) SI internally, unit-suffixed field names in the schema, CF `units` on every output variable | **(b)**. Numba kernels need plain arrays; the schema makes units unambiguous at the boundary. | Schema tests |
| D13 | Run identity | Hash of (canonical resolved scenario, data manifest checksums, model versions, package version) → `run_id`; a second hash over output arrays → `result_checksum` | Same inputs give the same `run_id` everywhere. Bitwise-identical output is guaranteed for the same platform and locked dependency set (`uv.lock`); across platforms equality is to a stated tolerance. No `fastmath`; parallel loops write per-index only. | Reproducibility test (run twice, compare checksums) |
| D14 | Interfaces | One service layer; Python API, CLI (`typer`) and MCP (`mcp` SDK) are thin adapters with no logic | Guarantees the three surfaces cannot drift. | Contract test: the three return identical payloads |

## 4. Architecture

```
src/lunarenv/
  api.py            run_scenario, describe_models, get_validation_status
  cli.py            lunarenv run | validate | report
  mcp_server.py     same three operations
  core/             model protocol + registry, provenance/hashing, result containers,
                    validity-range checks, cited constants
  data/             source catalog, fetch + sha256 + cache, manifest, DEM I/O and crop,
                    SPICE kernel set
  geometry/         projections and frames, ephemeris (SPICE), terrain normals
  illumination/     horizon kernels (numba), solar disc, products
  thermal/          property laws, grid, solvers, surface energy balance
  uq/               sampling, propagation, convergence tools (Richardson, GCI)
  scenario/         pydantic schema, YAML/JSON I/O, JSON Schema export
  validation/       case registry, runner, VALIDATION.md generator
  report/           HTML/Markdown report
tests/              unit, property (hypothesis), analytic, regression, validation
benchmarks/         timing and scaling scripts
docs/               references, model descriptions, convergence studies
scenarios/          shackleton_rim_1yr.yaml
```

**Extension contract.** A model implements `EnvironmentModel`: `name`, `version`, `requires` (outputs of other models), `describe()` (assumptions, simplifications, validity ranges, parameters with uncertainty, references) and `run(context, params) -> xarray.Dataset`. Models register through an entry-point group; the scenario holds one typed parameter block per model; the runner executes them as a dependency graph (illumination → thermal). The shared `context` carries site, terrain, ephemeris and time axis. Radiation, seismic or structural-load modules add a block and a class without touching the core.

**Every result carries**: units, uncertainty (quantiles over the ensemble plus a numerical-error estimate), model versions, data versions and checksums, `run_id`, `result_checksum`, and a list of validity warnings.

## 5. Data plan

| Item | Source | Size (approx.) | Handling |
|---|---|---|---|
| Site04 height, height error, slope error | `PGDA-78` | 41 MB each | Full files kept (downloaded, pinned) |
| Site04 clones (100) | `PGDA-78` | 4.1 GB raw | Streamed: download, crop to the near-field window, keep the crop, delete the raw file |
| Site01 height + clones | `PGDA-78` | as above | Same streaming; needed only for the `Barker2021` validation case |
| Mid/far-field LDEM, 80°S | `PGDA-90` | 189 MB (80 m), 716 MB (40 m), 2.7 GB (20 m) | 80 m downloaded and pinned; finer ones only if study D3 asks for them. No 120/240 m 80°S file exists; coarser levels are derived by block averaging |
| Average illumination and Earth visibility maps, 85°S 60 m | `PGDA-69` | < 0.1 GB | Validation only |
| SPICE kernels | `NAIF-GEN` | ~45 MB | `naif0012.tls`, `de440s.bsp`, `pck00011.tpc`, `moon_pa_de440_200625.bpc`, `moon_de440_250416.tf` |

Pinned checksums live in the package (`src/lunarenv/data/checksums.json`); the manifest (`<cache>/manifest.json`) records source, product version, URL, sha256, size, download time and citation of every cached file, and results embed the records they used. A mismatch against a pin is an error. Details, conventions and measured consistency checks: [docs/data.md](docs/data.md).

## 6. Illumination model

- Observer at height Δz above the DEM surface (parameter; structures are not at ground level).
- Horizon elevation per observer along N azimuths by great-circle marching over the nested DEMs, bilinear sampling, range step growing with distance, exact early termination (D4).
- Sun and Earth positions from SPICE with light-time and stellar aberration, converted to the local topocentric frame of the observer; solar angular radius from the actual distance.
- Visible fraction of the disc by strip integration against the interpolated horizon (D5); flux fraction with limb darkening.
- Products: time series at a point; mean illumination map over a period; longest continuous shadow and illumination periods (with their dates); optional Earth visibility with the same machinery.
- Stated simplifications: bare-surface DEM (no boulders below DEM resolution), no penumbra from diffraction, no light scattered from terrain (handled as a thermal-model caveat).
- Convergence studies: azimuth step, range step, switch radii, DEM resolution (5/10/20/40 m), time step (60/30/15/5 min), strip count. Reported metrics: mean illumination, longest shadow period, horizon RMS against the finest setting.

## 7. Thermal model

- 1-D heat conduction, `Hayne2017` property laws: ρ(z) exponential between 1100 and 1800 kg m⁻³ with scale H; K = K_c(ρ)[1 + χ(T/350)³], K_s = 7.4 × 10⁻⁴, K_d = 3.4 × 10⁻³ W m⁻¹ K⁻¹, χ = 2.7; c_p(T) quartic polynomial; ε = 0.95; albedo A(θ) = A₀ + a(θ/(π/4))³ + b(θ/(π/2))⁸, a = 0.06, b = 0.25; solar constant 1361 W m⁻² scaled by true distance; heat flow 0.018 W m⁻² at the base.
- Surface balance: absorbed direct solar flux on the local slope (normal from the DEM; incidence from the Sun vector; visible flux fraction from the illumination model) + conduction = εσT⁴.
- Grid: geometric, tied to the diurnal skin depth; domain deep enough for the **annual** wave (the polar seasonal signal is not negligible) — depth fixed by study.
- Spin-up: forced with the real illumination of the years preceding the scenario until the cycle-to-cycle change is below a stated threshold; the residual drift is reported.
- Outputs: T(t) at the surface and at depth, extremes, diurnal and seasonal amplitudes, maximum heating/cooling rates.
- Stated simplifications and expected effect:
  - **No terrain-scattered sunlight or terrain thermal emission.** Negligible when the surface is directly lit; dominant when it is in shadow, where the model will be too cold. An order-of-magnitude bound will be computed from the visible sunlit terrain and reported; results in shadow are flagged.
  - 1-D (no lateral conduction): valid when horizontal scales exceed the skin depth (centimetres to decimetres).
  - No rocks; regolith properties extrapolated from latitudes below 70° (finding 5).
  - Slope defined at the DEM baseline (5 m); sub-pixel roughness ignored.
- Convergence studies: first-cell size, growth factor, depth, time step, spin-up length; Richardson extrapolation and grid convergence index (`Roache1994`); observed versus theoretical order on analytic cases.

## 8. Uncertainty

| Input | Representation | Source |
|---|---|---|
| Near-field terrain height and slope | 100 PGDA clones | `Barker2021` |
| Far-field terrain | Treated as fixed; effect tested once | `Barker2021` states it is negligible |
| H, K_s, K_d, ρ_s, ρ_d, χ | Distributions bounded by published ranges (e.g. H from 0 to > 0.2 m, mean 0.06 m) | `Hayne2017`; ranges fixed and cited in Phase 3 |
| Albedo A₀, emissivity ε | Ranges from literature (ε: 0.95 adopted, 0.98 prior value) | `Hayne2017`, `Vasavada2012` |
| Heat flow | 0.018 W m⁻² with spread between Apollo sites | `Langseth1976` (to verify) |
| Ephemeris, solar constant | Shown negligible by a sensitivity test, then fixed | `Park2021`, `KoppLean2011` |
| Numerical error | Richardson estimate, reported separately from input uncertainty | Convergence studies |

Point results carry full Monte Carlo quantiles. Maps carry the nominal value plus Monte Carlo on a subsample of pixels unless the user requests the full ensemble (cost stated in the report). Ensemble size is set by a convergence check on the reported percentiles.

## 9. Verification and validation matrix

Tolerances are fixed **before** the case is first run and recorded in git; a failing case triggers investigation, not a tolerance change.

| ID | Case | Source | Metric | Status at planning |
|---|---|---|---|---|
| ILL-A1 | Smooth sphere: horizon dip versus observer height | Analytic | Max error in elevation | Ready |
| ILL-A2 | Synthetic terrain (wall, spherical-cap crater) with closed-form horizon | Analytic | Max error in elevation | Ready |
| ILL-A3 | Disc fraction for a straight horizon (circular segment) | Analytic | Max error in fraction | Ready |
| ILL-X1 | Topocentric Sun elevation/azimuth at the site | `HORIZONS` | Angular difference | Ready (network) |
| ILL-V1 | Site 1 regions of interest: mean illumination, longest illumination and shadow periods at 1 m and 5 m, 2024–2026, clone percentiles | `Barker2021` Table 2 | Difference per region and statistic | Ready (data public, table read) |
| ILL-V2 | 18.6-year mean illumination map, 60 m, around the site | `PGDA-69` (`Mazarico2011` method) | Bias, RMS, correlation; run on a DEM of matching resolution, with the 5 m result shown separately | Ready (data public); DEM-version mismatch is a known confounder |
| ILL-V3 | Shackleton rim point: continuous sunlight and longest darkness | `Mazarico2011` | Days | **Pending**: paper table not accessed |
| ILL-V4 | Connecting ridge: 92.27 % at 2 m, 95.65 % at 10 m over 19 years | `Glaser2014` (abstract) | Percentage points | **Pending**: location coordinates need the paper |
| TH-A1 | Constant properties, periodic surface temperature: damped wave | `CarslawJaeger1959` | Error and order of convergence | Ready |
| TH-A2 | Steady state with basal heat flow; zero-conductivity radiative equilibrium | Analytic | Error | Ready |
| TH-P1 | Energy conservation; T bounded by radiative-equilibrium limits; monotone response to albedo and flux | Property tests | Residual | Ready |
| TH-V1 | Equatorial diurnal curve: ~385 K noon, 101 K midnight, 95 K minimum | `Hayne2017` Table A2 (`Vasavada2012`) | K | Ready; original paper to be read for stated uncertainties |
| TH-V2 | Apollo 15 and 17 mean surface and subsurface temperatures (216/256 K, 211/252 K, ±5 K) | `Hayne2017` Table A2 (`Keihm1973`) | K, source-stated ±5 K | Ready |
| TH-X1 | Code-to-code against the reference implementation, identical parameters | `heat1d` | Max |ΔT| over a cycle | Ready (licence check first) |
| TH-V3 | Polar surface temperatures at the test site, sunlit pixels, seasonal maps | `Williams2019`, `DIV-PCP` | K at 240 m aggregation | **Informational**: data path unconfirmed; scattering neglected, so a bias is expected and will be explained, not tuned away |

Regression tests pin the outputs of the reference scenario by checksum.

## 10. Phases

Each phase ends with: tests green, convergence/timing results written to `docs/`, a sceptical self-review recorded in section 12, and a stop for the next go-ahead only where a physical decision is open.

| Phase | Content | Exit criteria |
|---|---|---|
| 1 Data | uv project, catalog, fetch/cache/manifest, DEM I/O with CRS assertions, crop, kernel set, frames and projection round-trip tests | Manifest complete; projection and frame tests pass; D6 check done |
| 2 Illumination | Horizon kernels, solar disc, products, Earth visibility, profiling, convergence studies, clone ensemble | ILL-A*, ILL-X1 pass; timings and scaling documented |
| 3 Thermal | Property laws, FV grid, both integrators, surface balance, spin-up, convergence, UQ | TH-A*, TH-P1 pass; default integrator chosen from data |
| 4 Verification | ILL-V*, TH-V*, TH-X1; VALIDATION.md generator | At least one illumination and one thermal published comparison reported with sources |
| 5 API/CLI | Scenario schema, runner, results + provenance, `run`/`validate`/`report` | One-year reference scenario runs; wall time documented; reproducibility test passes |
| 6 MCP | Thin server over the same service layer | Contract test passes; an agent session runs a scenario and reads validation status |
| 7 Report and docs | HTML/Markdown report, README, model docs, known limitations | Definition of done met |

## 11. Open decisions for the user

Defaults below are my recommendations; approving the plan without comment adopts them.

1. **Terrain-scattered and re-emitted radiation.** (A) Neglect, flag results in shadow, publish a bound — *recommended for the MVP*, consistent with the brief. (B) Add a single-bounce model at point receivers: makes shadowed temperatures and the polar Diviner comparison meaningful; costs roughly one extra phase.
2. **Second DEM tile for validation.** Download Site01 (connecting ridge) only to reproduce `Barker2021` Table 2 — *recommended*; it is the only like-for-like published benchmark found. The product's test site remains Site04.
3. **Map extent.** Default maps on a 2 × 2 km window at 5 m around the site, full tile available at coarser resolution — *recommended*. Full tile at 5 m is possible but is an overnight computation.
4. **Paywalled sources.** `Mazarico2011`, `Glaser2014`, `Vasavada2012`, `Williams2019`: if you can provide the PDFs, ILL-V3/V4 move from pending to ready and TH-V1 gets source-stated tolerances. Otherwise they stay pending, as required by principle 4.
5. **Licence.** Apache-2.0 (*recommended*: explicit patent grant, familiar to industry) or BSD-3-Clause (astropy ecosystem convention).

## 12. Self-review log

### Phase 0 — planning (2026-11-12)

Issues found while reviewing the first draft of this plan, and the resulting changes:

- The benchmark I intended for the Shackleton rim (`Barker2021` Table 2) is for Site 1. Validation now uses Site 1 for that case and says so, instead of comparing against numbers from a different location.
- The `Mazarico2011` site table could not be read (an NTRS search hit turned out to be `Zuber2012`). The case is marked pending rather than filled from memory.
- The Diviner-derived property model is fitted only within ±70° latitude and its c_p fit has a stated temperature range. Polar use is an extrapolation and must raise a validity warning; this was missing.
- A full-tile horizon map at 5 m was implicitly assumed feasible. The operation count says otherwise; maps are now scoped to a window, with sparse far-field computation as the measured optimization.
- "Reproducible" was stated without qualification. Bitwise identity across platforms cannot be promised with compiled numerical kernels; the guarantee is now stated precisely (D13).
- A fixed-radius DEM margin was arbitrary. Replaced by an analytic visibility bound with the realized maximum range recorded per run (D4).
- LOLA products are in the DE421 mean-Earth frame while the current ephemeris is DE440; mixing them silently would introduce a small unquantified rotation. The DE440 frame aligned to ME421 is used and the difference is tested (D6).
- Runtime figures in this document are estimates from operation counts. None has been measured; measured values replace them in Phases 2–5.
- Several references are cited from prior knowledge and are labelled [U] in the reference list; none of their numbers is used yet.

### Phase 1 — data (2026-11-12)

Delivered: uv project, data catalog with pinned checksums, verified download and cache with manifest, strict DEM reader/cropper, projection and local frames, SPICE session and ephemeris wrapper, line-of-sight bounds for the DEM extent. 78 tests pass (58 unit including property-based, 20 against real data); `ruff` and `mypy --strict` clean. Measured results are tabulated in [docs/data.md](docs/data.md).

Review as a sceptical reader, and what changed as a result:

- **A tolerance of mine was wrong, and the cause was physical.** I had asserted that computing Sun/Earth directions from the Moon's centre and subtracting the site vector is accurate to 1e-7 rad. The Earth case measured 4.7e-7 rad. Cause: stellar aberration (v/c ≈ 1e-4) is applied to the centre-to-body vector instead of the site-to-body vector, leaving (lunar radius / distance) × (v/c) = 5.1e-7 rad for the Earth. The tolerance is now that analytic bound, the docstring states it, and the measured value sits just under it (Sun: 1.1e-10 against 1.2e-9). Negligible for both bodies (Earth's angular radius is 1.7e-2 rad).
- **Server returns `200 OK` with an HTML page for missing files.** A status check would have cached garbage. Downloads now check the file signature, the length and the pinned checksum.
- **uv-managed Python ships no root certificates.** The first download failed TLS verification. Fixed by verifying against the OS trust store, not by relaxing verification.
- **Concurrent downloads lost a manifest record** (two of my own processes wrote the manifest at once). Manifest updates are now a locked read-modify-write, and a file already present that matches its pin is adopted instead of downloaded again.
- **My co-registration test assumed the tile sits on the 80 m lattice.** It is offset by 40 m. The test now compares whole 80 m cells; result 0.46 m RMS at zero shift against 21 m at one cell.
- **Common-mode risk.** Both PGDA products share one registration convention, so agreeing with each other does not prove either is placed correctly. Added an independent check against Shackleton's published centre, diameter, rim height and floor radius (`Zuber2012` Table 1): centre within 243 m, diameter 21.1 km.
- **The SPICE session cleared the global kernel pool on entry**, which would break a host application's own SPICE use. It now requires an empty pool and unloads only what it loaded.
- **A kilometre/metre heuristic in the DEM reader rejected flat synthetic terrain** needed for analytic tests. Removed; units are pinned per product by regression tests on height statistics instead.
- **The DEM-extent criterion had an unstated condition**: the bound at the DEM edge covers everything farther out only beyond the grazing point (~64 km). The function now enforces it.
- **Tile size**: 16 × 16 km, not 20 × 20 km; the pole lies inside the tile. Projection and frame code were already written to be regular at the pole; a test covers it.

Open items carried forward:

- The 80S product is proven sufficient for 98 % of the tile (all observers above −1426 m). The deep crater interior is not proven; Phase 2 applies the bound per ray and reports violations.
- 199 of 200 clones are not yet pinned or ingested. PGDA delivered 0.4–1 MB/s, so each site ensemble costs roughly one to three hours of download. Phase 2 must fix the near-field crop window before ingesting.
- Site 1 tile and the published illumination maps are catalogued but not downloaded; needed in Phases 2 and 4.
- No continuous integration yet; data-backed tests skip when the cache is absent.
- The terrain height cap (11 km) relies on a global maximum of about 10.8 km cited from memory (`Smith2010`, status U). It only needs to be an upper bound, but the source should be confirmed.

### Phase 2 — illumination (2026-11-12)

Delivered: nested terrain stack, great-circle horizon tracer with an exact stopping rule and per-ray completeness, strip-integrated solar disc with limb darkening, point time series and map statistics, terrain-clone ensemble, numerical settings with a documented error estimate, convergence, timing and ensemble scripts, two figures. 137 tests pass (115 unit, 22 against real data); `ruff` and `mypy --strict` clean. Documents: [docs/models/illumination.md](docs/models/illumination.md), [docs/convergence/illumination.md](docs/convergence/illumination.md).

Results that settle open decisions:

- **ILL-A1, A2, A3 pass** (smooth sphere < 1e-7 rad; analytic hill < 2e-5 rad against an independent geodesic computation; disc segment and limb-darkened segment < 2e-6).
- **ILL-X1 passes**: Sun direction at the rim site agrees with JPL Horizons to 9e-10 rad, Earth to 9e-9 rad (tolerance 2e-6, fixed beforehand). This validates time conversion, frames, projection and topocentric reduction. It does not validate the ephemeris itself, since both sides use JPL ephemerides.
- **D3 (nested levels)**: 5 m data matter out to at least 4 km; the 80 m mid field is the largest numerical error (+0.04 pp against the 40 m product, first-order convergence, about +0.08 pp extrapolated); the far field at 240 m costs 0.01 pp.
- **D4 (DEM extent)**: every ray from 25 600 observers spread over the whole tile was closed by the stopping rule inside the 80S product. The open item from Phase 1 (deep crater interior) is resolved for this tile.
- **D5 (disc)**: per-strip horizon and single chord differ by 0.001 pp on the yearly mean; limb darkening changes the flux-weighted mean by 0.01 pp. Both are kept because they matter at individual epochs, which is what the thermal model consumes.
- **Dominant uncertainty is the terrain**: block-averaging the site DEM to 10 m changes the yearly mean by +0.28 pp, ten times all numerical parameters together; the clone ensemble gives +/-1 pp at 2 m on the rim crest.
- **Candidate test point**: (-2847.5, -9052.5) m on the rim crest, 84.3 % mean illumination in 2026 at 2 m (ensemble 83.3 +/- 1.0 %), longest shadow 2.7 days.

Review as a sceptical reader, and what changed as a result:

- **My cost estimate for maps was wrong by a factor of three.** I extrapolated from three observers that each had a thread to themselves. The 2 x 2 km window at 5 m needs about 4.7 CPU-hours, not about 75 CPU-minutes. I also found the parallel loop was leaving most cores idle (one block of observers per thread, with very unequal cost) and fixed it with small blocks. Decision 3 of section 11 (default window) needs revisiting: see open items.
- **Wall-clock timings from this session are unreliable.** Other workloads on the machine kept the load average between 25 and 38. Timings are reported in CPU time; the "minutes on an idle machine" figures are derived, not measured.
- **A property test caught a real defect.** A fully visible limb-darkened disc returned 1 - 1e-6 because the flux was normalized analytically while integrated numerically. The normalization now uses the same quadrature; clear discs give exactly 1.
- **Three of my test expectations were wrong, not the code**: the tangent distance must be measured on the reference sphere; a cliff sampled at 0.25 deg is smeared on the scale of the solar disc; my brute-force flux reference had its own half-row error (replaced by a closed form integrated along rows, independent of the kernel's strips).
- **A time axis uniform in UTC is not uniform in TDB** at the millisecond level; my uniformity check rejected valid input. It now tolerates the relativistic periodic terms and a leap second.
- **A numba parallel-analysis failure** (tuples returned by inlined helpers inside a parallel loop) was worked around by writing the local basis in scalars; the scalar algebra is checked against the numpy implementation through the series test.
- **The largest numerical error is one I had not planned to measure directly.** Only the comparison against the 40 m regional product showed that the mid-field resolution converges at first order. The default stays at 80 m with the estimate attached to results, and the 40 m product is selectable.
- **The ensemble under-reports uncertainty when the skyline is distant.** Clones perturb only the site tile. At 10 m the spread is 0.01 pp, below the numerical estimate; the regional DEM's own error is not sampled. This is stated in the model document and both numbers are reported separately.
- **Ensemble size is 13 of 100.** The server limits each connection to about 100 kB/s and refused requests at eight parallel files. Downloads now use byte ranges, back off on 503/429, and continue in the background. Percentiles beyond the 5th and 95th are not meaningful until the set is complete.
- **Lunar eclipses are not modelled**: when the Earth covers the Sun the model still reports sunlight. Irrelevant for yearly means, relevant for the thermal model (hours of cooling). To be added in Phase 3 as a disc-overlap factor.
- **"Lit" means any sliver of the disc**, following Barker2021. Longest-period statistics are sensitive to that threshold; a configurable threshold belongs in the scenario schema (Phase 5).

Open items carried forward:

- Regenerate the ensemble table with all 100 Site 4 clones; download the Site 1 tile and clones for ILL-V1 (Phase 4). Budget several hours of download per site.
- The 1 x 1 km, 5 m window map and its figure were still computing when this was written.
- Default map extent: a 2 x 2 km window at 5 m and 1440 azimuths costs about 4.7 CPU-hours. Options: keep it and document the cost; default to 720 azimuths (half the cost, 0.003 pp); default to 1 x 1 km. Recommendation: 1 x 1 km at 1440 azimuths by default, larger windows on request.
- Comparison with published illumination (ILL-V1, ILL-V2) is the decisive validation and is Phase 4.

### Phase 3 — thermal (2026-11-12)

Delivered: Hayne2017 property laws, depth grid, conservative finite-volume solver with four time integrators, surface energy input on the DEM slope, lunar-eclipse factor, point model with exact-mean spin-up, parameter ranges with Latin hypercube sampling, site study script (results, convergence, model form, uncertainty), figure. Test suite: 188 pass (165 unit, 23 against real data); `ruff` and `mypy --strict` clean. Documents: [docs/models/thermal.md](docs/models/thermal.md), [docs/convergence/thermal.md](docs/convergence/thermal.md).

Decisions settled by measurement:

- **D8**: finite volumes on enthalpy conserve energy to 4e-15 (relative) over the test year.
- **D9**: BDF2, two steps per hourly sample: 0.019 K at a third of the cost of forward Euler; Crank-Nicolson rejected because it overshoots the physical bound by 15 K after an abrupt change (kept as a test).
- **Completed open item of Phase 2**: the 1 x 1 km window at 5 m (40 000 observers, 720 azimuths) gave the test point (-2822.5, -9062.5) m with 85.5 % mean illumination at 2 m; lunar eclipses are now modelled.

Review as a sceptical reader, and what changed as a result:

- **I replaced my own discretization mid-phase.** The first version averaged conductivities between cells and needed more than 300 spin-up cycles. Noticing that the conductivity is separable led to writing fluxes in the Kirchhoff variable: steady states became exact on any grid and the spin-up could impose the exact mean profile, converging in 8-10 cycles independently of the initial state (2e-4 K). A test that expected second-order convergence of the steady state had to be rewritten because the error became round-off.
- **Deviation from the plan**: I did not write a separate node-centred finite-difference solver "as in Hayne2017". The published reference code heat1d plays that role and was run instead.
- **The reference code is not converged in two ways, and I checked rather than assumed.** heat1d's default one-year equilibration leaves its deep temperature 9 K high (261.7 K against 252.5 K after 20 and 80 years). Once equilibrated it is still 2.8 K above lunarenv at depth while agreeing within 0.5 K at the surface. lunarenv's deep value moves by 0.013 K across a 14-fold range of grids; heat1d cannot be refined (NaN on every finer grid). I attribute the gap to its averaging of conductivities; that is a hypothesis, stated as such in the model document.
- **A discrepancy with published data that I could not resolve.** With the published parameters the model is 8-10 K below the Apollo mean temperatures quoted by Hayne2017 (outside their +/-5 K), while it matches the equatorial noon, midnight and minimum values of the same paper. A constant albedo reproduces both the Apollo values and the paper's Figure A1. No tolerance was touched; the case goes to VALIDATION.md as it stands (Phase 4).
- **The polar result hinges on that same question.** At the test point the incidence-dependent albedo law, fitted on level equatorial ground, is used at 88 degrees incidence where it gives 0.79. A constant albedo changes the yearly mean by +26 K and the maximum by +50 K. This is larger than every sampled uncertainty and is reported as a model-form sensitivity, not hidden inside the ensemble.
- **The headline number is not what a reader might expect.** The "best lit" point (77 % of the year at ground level) has a surface between 43 and 119 K, because the ground there is almost level and the Sun never rises above 1.8 degrees. Illumination fraction is not absorbed energy. Slope is the leading sensitivity (rank correlation 0.8-0.9), and the slope of a 5 m DEM that is 90 % interpolated is poorly known.
- **Temperatures go far below the stated range of the heat-capacity fit** (43 K against 90 K). Flagged in every result; a low-temperature property model is the obvious next addition.
- **A result that looked like a bug and was not**: removing lunar eclipses changed nothing. Eclipses occur with the Sun near the lunar equatorial plane, which at the pole is below the local horizontal, where level ground absorbs nothing.
- **Two parameter intervals are assumptions** (surface density, heat flow) and one parameter (chi) is not varied for lack of a published spread. Marked in code and documents.
- **The terrain ensemble has 28 of 100 clones.** Quantiles are stable from 32 to 128 members of the parameter sample, but the terrain part is still a partial sample.

Open items carried forward:

- VALIDATION.md entries for TH-V1 (equator: pass expected), TH-V2 (Apollo: fails with the published albedo law, passes with constant albedo; report both), TH-X1 (heat1d: surface agrees, depth differs by 2.8 K with an explanation).
- Obtain Vasavada2012 and Keihm1973/Langseth1976 to establish what the Apollo "mean surface temperature" measures and which albedo the published fits used.
- Low-temperature heat capacity (WoodsRobinson2019) as a selectable property model.
- Scenario runner (Phase 5) must expose: surface orientation override, albedo law choice, forcing sampling, and report both uncertainty and model-form sensitivities.

### Phase 4 — verification (2026-11-12)

Delivered: validation package (case registry with source, metric and tolerance in code; runner; results stored with the package; Markdown report), 16 cases, [VALIDATION.md](VALIDATION.md), reference fixtures (JPL Horizons tables, heat1d outputs and the script that produced them), `benchmarks/run_validation.py`. Test suite: 196 pass (173 unit, 23 against real data); `ruff` and `mypy --strict` clean.

Outcome: 10 PASS, 1 FAIL, 2 PARTIAL, 3 PENDING.

| Case | Result |
|---|---|
| ILL-V1, Barker2021 Table 2 (Site 1 regions) | 1st percentile of average illumination within 0.88 pp of the published values for all 14 region-height pairs (tolerance 3 pp); at 5 m within 0.2 pp. Region areas, rebuilt from the published definition, within 0.7-2.6 %. PARTIAL: 12 of 100 clones. |
| ILL-V2, PGDA 18.6-year illumination map | bias -1.60 pp, RMS 4.12 pp, r = 0.959 over 10 000 cells of 60 m (tolerances 2 pp, 6 pp, 0.95). PASS, with a systematic difference described below. |
| TH-V1, equatorial Diviner temperatures | 385.2 / 100.3 / 93.7 K against 385 / 101 / 95 K. PASS. |
| TH-V2, Apollo mean temperatures | 8-10 K low with the published parameters. FAIL, unresolved. Passes with a constant albedo (TH-V2b). |
| TH-X1, heat1d | surface within 0.75 K at three latitudes; 1.9-2.9 K colder at depth. PARTIAL. |
| ILL-V3, ILL-V4, TH-V3 | PENDING: references not obtained (two paywalled papers, one unconfirmed data path). |

Review as a sceptical reader, and what changed as a result:

- **Not every tolerance predates its result.** ILL-V1 (3 pp) and ILL-V2 (2 pp, 6 pp, 0.95) were fixed before their first run. TH-V1 (2 K) and TH-X1 (1 K) were written down after I had already seen exploratory numbers in Phase 3. Both are anchored to figures stated by the source, and both cases say so in their notes; the header of VALIDATION.md no longer claims otherwise.
- **ILL-V2 passes with thin margins and a structure in the error.** On cells the published map puts above 50 %, lunarenv at ground level is 6 pp darker. I investigated before accepting the pass: with the observer 2 m up the bias on those cells becomes +1.2 pp, at 10 m +10 pp. The difference is worth about 1.5 m of observer height; the published product's height is undocumented. Probable cause, not established: lunarenv samples the slope of the observer's own 60 m cell, a node-to-node horizon does not. The registered case was left at ground level, as defined before running.
- **ILL-V1 is strong but provisional.** The published numbers are percentiles over 100 clones; I estimate them from 12 with a normal approximation. The agreement (0.88 pp at worst) is better than the sampling error of that estimate (about half a standard deviation), so part of it is luck. It must be re-run with the full ensemble, using empirical percentiles.
- **An assumption about the reference had to be made.** Barker2021 prints the centroids as positive kilometres while the tile is at negative coordinates. I took them as (-X, -Y). Evidence that this is right: they fall on pixel centres, and the regions grown from them reproduce the published areas within 3 %.
- **Longest-period statistics are reported without a tolerance.** They are discontinuous functionals: one hour of grazing light splits a period. Most agree within a day or a few days; one differs by 44 days. I did not set a tolerance after seeing that.
- **The Apollo failure stands.** No parameter was adjusted. What is established: the model reproduces Diviner equatorial extremes and the reference code's surface temperatures with the published parameters, and reproduces the Apollo means and the paper's own Figure A1 only with a constant albedo. What is not: which albedo applies, and what the Apollo "mean surface temperature" measured. This matters directly for polar results (Phase 3).
- **heat1d is treated as a reference although it is demonstrably not converged at depth** (and unstable on finer grids). The case is PARTIAL rather than PASS or FAIL, with the evidence and the hypothesis stated; calling it PASS would overstate, calling it FAIL would blame the wrong code without proof.
- **Thermal validation at the pole does not exist yet.** Every thermal comparison is equatorial or mid-latitude and terrain-free. TH-V3 is pending for lack of data access, and would in any case be limited by the neglected terrain radiation.
- **A defect in my own report generator**: metrics containing a vertical bar broke the Markdown table. Fixed and tested.
- **Stored results can go stale.** `results.json` records the package version and data checksums but nothing forces a re-run when code changes. Phase 5 should record the source revision and warn when it differs.

Open items carried forward:

- Re-run ILL-V1 with all 100 Site 1 clones and empirical percentiles; finish the Site 4 clone download and regenerate the Phase 2 and 3 ensembles.
- Obtain Mazarico2011 Table 3, Glaser2014, Vasavada2012, Keihm1973/Langseth1976; locate the Diviner polar products.
- Decide the albedo question with the user before the report presents polar temperatures as predictions.

#### Addendum: ILL-V1 with all 100 clones (2026-11-12, later the same day)

Re-run with the complete Site 1 ensemble and percentiles taken from the ensemble itself (no normal approximation). Result: **PASS**. The 1st percentile of the region-average illumination differs from Barker2021 Table 2 by at most 2.22 pp at 1 m (region 1; the other six are within 1.6 pp) and by at most 0.33 pp at 5 m; tolerance 3 pp. The 99th percentile of the longest shadow period is within 3.1 days at 1 m and 0.7 days at 5 m. Region areas are within 2.6 %. The run took 87 minutes.

What this changes in the review above: the 0.88 pp obtained from 12 clones was indeed flattered by the normal approximation and the small sample, as suspected; the honest figure is 2.2 pp. At 1 m lunarenv's lower tail is slightly brighter than the published one in every region (+0.1 to +2.2 pp, mean +1.0). That is the same sign as the mid-field resolution bias measured in Phase 2 and as the far-field DEM difference, but it was not investigated further because it is inside the tolerance fixed beforehand. The longest illumination period at 1 m differs by up to 49 days in one region; it has no tolerance (discontinuous statistic) and is reported as found.

### Phase 5 — scenarios, API and command line (2026-11-12)

Delivered: strict scenario schema (pydantic, exported to [docs/scenario.schema.json](docs/scenario.schema.json)); model protocol and registry with self-descriptions; run identity and result checksum; the runner (illumination, thermal, constant-albedo alternative, terrain and parameter uncertainty, optional map); results as NetCDF-4 with CF attributes plus `metadata.json` and `scenario.json`; Python API (`run_scenario`, `describe_models`, `get_validation_status`); command line (`lunarenv run | validate | describe | schema | data`); reference scenario [scenarios/shackleton_rim_2026.yaml](scenarios/shackleton_rim_2026.yaml). Test suite: 227 pass (197 unit, 30 against real data); `ruff` and `mypy --strict` clean.

Reference run (`lunarenv run scenarios/shackleton_rim_2026.yaml`, run `f93b88484006a71f`): one year hourly, 16 terrain clones, 64 parameter samples, **109 s wall-clock** on the M1 Pro while other workloads held the load average near 40. Breakdown: nominal illumination 1 s, nominal thermal with its alternative 2 s, uncertainty 101 s. Output 4.1 MB.

What a result carries: units and a long name on every variable; quantiles over the ensemble; numerical error estimates; versions of the package and of each model; id, URL and SHA-256 of every input file; software environment and source revision; the list of validity warnings; the stored validation status; `run_id` and `result_checksum`.

Decision adopted without objection from the user (asked at the end of Phase 4): the published incidence-dependent albedo is the nominal result and the constant-albedo result is computed and shown alongside by default.

Review as a sceptical reader, and what changed as a result:

- **Tests found two defects before any user did.** The model registry was empty unless the runner module happened to be imported, so `describe_models` could return nothing; and one coordinate had no units, breaking the promise that everything does. Both fixed; both now have tests.
- **A test expectation of mine was wrong again** (a map window is rounded outwards to whole pixels: 5 x 5, not 4 x 4). The code was right.
- **Reproducibility is verified, with its limits stated.** Running the same scenario twice gives the same `run_id`, the same `result_checksum` and bit-identical arrays; changing the uncertainty seed changes the identity and the ensemble but not the nominal results. This is tested on one machine. `run_id` does not include the software environment, so across platforms the same id may come with different last bits; the environment is recorded to tell such cases apart.
- **"Within validated range" is blunt.** It is false whenever any warning exists, which for the thermal model at the pole is always (properties extrapolated beyond 70 degrees of latitude). That is honest, but a per-model flag would be more useful to a reader. Left as is for the MVP.
- **Deviations from the plan**: the command line uses `argparse`, not `typer` (one dependency fewer; the interface is small); `lunarenv report` is delivered with the report generator in Phase 7, not here.
- **Not done**: the configurable threshold for "lit" announced in Phase 2. Longest-period statistics still count any sliver of the disc.
- **Stale validation status** (raised in Phase 4): runs now record the source revision, but nothing yet compares it with the revision the validation results were produced with.
- **Cost is dominated by uncertainty**, which recomputes horizons and series for every terrain clone. It scales linearly with clones and parameter samples and is the place to optimize if runs become routine.

Open items carried forward:

- `lunarenv report <run>` (Phase 7); configurable lit threshold; per-model validity flags; comparison of run revision and validation revision.
- Re-run the reference scenario once the Site 4 clone ensemble is complete (62 of 100 at the time of writing).

### Phase 6 — MCP server (2026-11-12)

Delivered: `lunarenv-mcp` (stdio), three tools that call the API functions of the same name, the scenario schema as a resource, server instructions for agents, [docs/mcp.md](docs/mcp.md). One summary function (`api.run_summary`) now feeds the command line's `--json` output and the MCP tool, so the interfaces cannot drift. Test suite: 233 pass (202 unit, 31 against real data); `ruff` and `mypy --strict` clean.

Exit criteria met: a contract test runs one scenario through the Python API, the command line and an MCP session started as a subprocess over stdio, and requires identical payloads (wall-clock timings excepted); in the same session the client lists the tools and reads the validation status.

Review as a sceptical reader, and what changed as a result:

- **My knowledge of the SDK was out of date.** The installed `mcp` is 2.3, where the server class and the client API differ from what I would have written from memory. I inspected the installed package instead of guessing, and the tests exercise the real protocol path (subprocess, stdio), not only in-process calls.
- **The first version crashed on bad input.** An invalid scenario surfaced as a generic "error executing tool". Anticipated failures (invalid scenario, site outside the tile, data missing) are now returned as tool errors carrying the original message; this is the only code in the server beyond delegation, and it is tested.
- **"No logic of its own" has two small exceptions**, both deliberate: that error translation, and running the computation in a worker thread so the server keeps answering during a one-to-two-minute run.
- **What the server cannot guarantee.** Its instructions ask the agent to report uncertainty, warnings and the constant-albedo alternative with the nominal values. A client may ignore instructions; the safeguard is that those fields are in the same payload as the numbers, not behind another call.
- **Not tested**: a long run against real client timeouts; concurrent tool calls (SPICE's process-wide state makes a second simultaneous run fail with a clear error rather than interleave, but there is no test of that through MCP).
- **Scope kept narrow**: stdio only. The SDK's HTTP transports are not enabled; serving on a network would need authentication.

### Phase 7 — report and README (2026-11-12)

Delivered: report generator (`lunarenv report <run>`: one self-contained HTML file, or Markdown with PNG figures), its design system ([DESIGN.md](DESIGN.md)), and the [README](README.md) with installation, a complete example and known limitations. Test suite: 236 pass (202 unit, 34 against real data); `ruff` and `mypy --strict` clean.

The report reads only `results.nc` and `metadata.json`. Order of content: identity of the run, validity banner with every warning, nominal values next to their ensemble quantiles, the constant-albedo alternative, seven charts, uncertainty and numerical error, model assumptions, validation status, provenance.

Review as a sceptical reader, and what changed as a result:

- **I looked at the output instead of trusting it.** The report was opened in a browser at 1440 px and 375 px, in light and dark themes: no console errors, no horizontal scroll at phone width, charts follow the theme. Chart colours (two categorical series, a four-step ordered ramp) were run through a palette validator in both themes and pass.
- **A test caught a latent defect in my theming trick.** Charts are drawn with placeholder colours that are swapped for CSS variables; one placeholder had the same value as the real text colour, so any ink drawn in that colour would have been silently repainted. Placeholders now use a pattern that cannot collide, and an assertion enforces it.
- **The first report was 1.5 MB** because the Sun's track was 8760 separate marks. Drawn as lines it is 0.7 MB.
- **Design process shortcuts, stated**: the design system was written and applied without asking for approval first, and without studying reference products; it follows general rules (hierarchy before decoration, state as icon plus label, one accent) rather than a reviewed design. Charts are static: there is no hover read-out, which interactive charts would normally have; the numbers are in `results.nc`.
- **The report can only be as honest as the metadata.** It prints the stored validation status, which is not recomputed per run, and the warning list, which is blunt (Phase 5 review). It does put validity above the numbers and the albedo alternative beside the nominal value, which was the point.
- **README claims were checked against outputs** (commands run, numbers copied from the reference run and from VALIDATION.md), not written from memory.

## 13. Definition of done: status

| Criterion | Status | Evidence |
|---|---|---|
| `lunarenv run` executes a one-year scenario at the test site locally, with the time documented | Met | Run `f93b88484006a71f`: 109 s wall-clock on an M1 Pro under load (16 terrain clones, 64 parameter samples); README, Phase 5 |
| VALIDATION.md shows at least one illumination and one thermal case against published data, with sources | Met | ILL-V1 (Barker2021, pass), ILL-V2 (PGDA map, pass), TH-V1 (Diviner via Hayne2017, pass) |
| Each model has its convergence study and its uncertainty documented | Met | docs/convergence/illumination.md, docs/convergence/thermal.md, docs/models/ |
| The MCP server lets an agent run a scenario and query the validation status | Met | tests/data/test_mcp_contract.py: real stdio session, identical payloads across API, CLI and MCP |
| README with installation, complete example and known limitations | Met | README.md |

Non-negotiable principles, checked against the code:

1. *Deterministic physics, no LLM in the computation*: true by construction; the MCP server only delegates.
2. *Units, uncertainty, model and data versions with checksums, out-of-range warning on every result*: in `results.nc` attributes and `metadata.json`; tested.
3. *Reproducibility by hash*: `run_id` and `result_checksum`; tested bit for bit on one platform.
4. *No invented validation data*: every case cites its source; three cases are PENDING rather than filled from memory; two tolerances written after exploratory runs are disclosed.
5. *Assumptions, simplifications and validity range per module*: `describe_models`, module docstrings, docs/models/.

## 14. What remains open

Honest list, in order of importance for a user of the numbers:

1. **Albedo at grazing incidence** (26 K on the polar mean). Needs Vasavada2012 and the Apollo heat-flow papers, then a decision.
2. **No polar thermal validation**; terrain-scattered and terrain-emitted radiation are neglected (decision 1 of section 11).
3. **Low-temperature heat capacity**: results reach 43 K, the fit is stated from about 90 K.
4. ~~ILL-V1 with the full clone ensemble~~ done (PASS, see the Phase 4 addendum). Phase 2 and 3 ensemble tables regenerated with all 100 Site 4 clones.
5. **Pending references**: Mazarico2011 Table 3, Glaser2014, Diviner polar products.
6. Smaller items: configurable "lit" threshold; per-model validity flags; warning when validation results predate the code; hover values in report charts. (Done since: the repository has a commit history and continuous integration on Linux for lint, types, unit tests and kernel-only data tests; terrain-backed tests still run only locally.)
