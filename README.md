# lunarenv

Open, verified lunar surface environment engine. Given a point near the lunar south
pole and a period of time, it computes the **illumination** (terrain horizon, visible
fraction of the Sun, Earth visibility, eclipses) and the **regolith temperature** at the
surface and at depth, with uncertainty, and with every result traceable to its input
data, model version and validation status.

It is meant for engineers who design structures for the lunar surface and need
environmental conditions they can audit. The same three operations are available as a
Python API, a command line and an MCP server for AI agents. All physics is computed by
deterministic code.

**Status: MVP, pre-alpha (0.1.0.dev0).** One test site (Shackleton crater rim) and two
models. Read [Known limitations](#known-limitations) before using a number.

## Related work and what is new here

Most ingredients exist already; lunarenv is an integration with an auditable validation
ledger, not new science (finding F10 of the [scientific review](docs/review/)).

- **Illumination.** The method is that of Barker et al. 2021 (nested 5 / 80 / 240 m
  LOLA terrain, limb-darkened disc, 100 terrain error clones), whose clones are open
  (PGDA 78). PGDA 69 distributes open long-term average maps of Sun and Earth
  visibility; Gläser et al. 2018 give results for named polar sites.
- **Thermal.** [heat1d](https://github.com/phayne/heat1d) (Hayne, MIT licence) is the
  open implementation of the Hayne et al. 2017 regolith model that lunarenv's thermal
  model follows. Terrain-scattered sunlight and terrain infrared, which lunarenv models
  with a single bounce, are in open codes with multiple scattering: Topo3D in Schörghofer's
  [Planetary-Code-Collection](https://github.com/nschorgh/Planetary-Code-Collection),
  [fluxpy](https://github.com/sampotter/fluxpy) (Potter), and the Oxford 3D
  thermophysical model (King et al. 2020). The review also reports several open 2026
  repositories that combine LOLA terrain and SPICE geometry, some with Hayne-type
  thermal models (not examined here).

What lunarenv adds, as far as the review found: point time series at a chosen height
with terrain-clone percentiles from open code; joint propagation of terrain clones and
thermophysical parameters into point temperatures; a validation ledger with sources,
tolerances (pre-registered since the review) and recorded FAIL and pending cases,
including a comparison with Diviner polar temperatures; run hashing and data provenance;
and a Python, command-line and MCP interface.

## Install

Requires Python 3.12 and [uv](https://docs.astral.sh/uv/). Developed on macOS, Apple
Silicon; nothing depends on a GPU.

```bash
git clone https://github.com/sergio-27/lunaware.git
cd lunaware
uv sync
```

Fetch the data (about 360 MB: SPICE kernels, the 5 m site terrain, the 80 m regional
terrain). Every file is checked against a pinned SHA-256.

```bash
uv run lunarenv data fetch kernels site04 regional
```

Optional, for the regional terrain error bound: its error map (255 MB). Without it runs
warn and report no bound.

```bash
uv run lunarenv data fetch regional-error
```

Optional, for terrain uncertainty: the 100 error clones of the site terrain, 41 MB each.
The source server delivers about 100 kB/s per connection, so allow one to two hours.

```bash
uv run lunarenv data fetch site04-clones
```

The cache lives in `$LUNARENV_DATA_DIR` if set; otherwise in `.data/` at the root of the
repository, if that directory exists (the location the tests use); otherwise in the
per-user cache directory.

## A complete example

The reference scenario, [scenarios/shackleton_rim_2026.yaml](scenarios/shackleton_rim_2026.yaml):

```yaml
schema_version: 1
name: shackleton-rim-2026
site:
  terrain: site04
  x_m: -2822.5            # south polar stereographic metres, as in the LOLA products
  y_m: -9062.5
  height_above_surface_m: 2.0
period:
  start_utc: 2026-01-01T00:00:00
  end_utc: 2027-01-01T00:00:00
  step_minutes: 60
thermal:
  albedo_envelope: true   # cold (Bürger 2024) and hot (Foote 2020) albedo laws
  slope_grid: {slopes_deg: [0, 2, 5, 10], aspects: 8}
  oriented_surfaces:      # asset surfaces: direct sunlight and view factors
    - {name: roof, tilt_deg: 0}
    - {name: wall-0, tilt_deg: 90, facing_azimuth_deg: 0}
    # ... (three more walls and a 45 degree panel)
uncertainty:
  members: 64             # regolith-parameter samples
  terrain_clones: 16      # terrain error clones (uses those in the cache)
```

Run it and build the report:

```bash
uv run lunarenv run scenarios/shackleton_rim_2026.yaml
```

```bash
uv run lunarenv report <run_id>
```

Output of the run (16 clones in the cache and the regional error map; 941 s on an M1 Pro
with the test suite sharing the machine. Before the regional bound was added the same
scenario took 563 s on an otherwise idle machine, about 530 s of it for the terrain
patches seen by the asset walls when the machine was shared):

```text
run 8bdc60ddb18e5bd9  ->  runs/8bdc60ddb18e5bd9
illumination: mean visible solar disc 85.51 %, longest shadow 63 h, longest illumination 2803 h
              terrain ensemble 5/50/95 %: 85.35 / 85.54 / 85.69 %
              regional DEM error bound (+/-2 rms): 85.50 - 85.52 %, longest shadow 63 - 63 h
thermal bounds (albedo envelope x ensemble 5-95 %): surface minimum 41.3-60.4 K, maximum 102.9-200.1 K, mean 61.0-108.6 K
         nominal terrain: cold law max 106.7 K, scenario law 119.2 K, hot law 163.0 K
thermal, scenario law: surface 43.6 to 119.2 K (mean 70.4 K) on a 0.39 deg slope
         surface maximum, ensemble 5/50/95 %: 105.6 / 118.4 / 146.8 K
         with constant albedo: mean 96.7 K (+26.3 K); not settled by validation
         with woods_robinson2019 conductivity: minimum -3.3 K, mean -1.6 K; means at depth 0.05 m +2.5 K, 0.1 m +3.6 K, 0.3 m +5.2 K, 1 m +7.8 K
         slope grid 0-10 deg, all aspects and laws: surface maximum 101.9 to 252.9 K
         diffuse at the minimum: terrain 8.4e-05 W m-2, Earthshine absorbed 0.006 W m-2; minimum raised +0.01 K (terrain), +0.23 K (Earthshine)
surface roof (tilt 0, facing 0): direct sun mean 12 W m-2, max 45 W m-2, lit 50 %; view factor sky 1.000, ground 0.000
         diffuse mean: terrain 0.00049 W m-2 scattered, 0.00037 W m-2 infrared; Earthshine 0.0034 W m-2
surface wall-0 (tilt 90, facing 0): direct sun mean 373 W m-2, max 1401 W m-2, lit 46 %; view factor sky 0.517, ground 0.483
         diffuse mean: terrain 5.5 W m-2 scattered, 2.3 W m-2 infrared; Earthshine 0.053 W m-2
...
WARNING thermal: temperatures reach 43.6 K, below the 90 K down to which the heat-capacity fit is stated to be valid
WARNING thermal: regolith properties were derived from data within 70 degrees of the equator; their use at this latitude is an extrapolation
WARNING thermal: 17.7 % of the period has no direct sunlight; single-bounce terrain radiation and Earthshine are included, multiple scattering and sub-DEM roughness are not
WARNING thermal: mean temperatures at depth change by up to 7.8 K with the other low-temperature conductivity law (woods_robinson2019); treat depth results as bounded by both (model_form_sensitivity.conductivity)
```

The run directory holds `results.nc` (NetCDF-4, CF attributes, units on every
variable), `metadata.json` (provenance, uncertainty, validity, validation status),
`scenario.json` (the scenario with every default spelled out) and, after
`lunarenv report`, a self-contained `report.html`.

From Python:

```python
import lunarenv

result = lunarenv.run_scenario("scenarios/shackleton_rim_2026.yaml", output_dir="runs")
ds = result.dataset  # xarray.Dataset
ds["surface_temperature"].max()  # K
result.metadata["uncertainty"]["thermal"]["surface_max_k"]["quantiles"]  # 5, 50, 95 %
result.metadata["uncertainty"]["thermal_envelope"]["bounds"]["surface_max_k"]  # [lower, upper] K
result.metadata["validity"]["warnings"]

lunarenv.describe_models()  # assumptions, validity range, numerical error
lunarenv.get_validation_status()  # outcome of every validation case
```

Other commands: `lunarenv validate` (stored validation status with the date each case was
computed; `--run` to recompute, `--slow` to include the cases that take minutes; cases
not recomputed are listed with the reason),
`lunarenv describe`, `lunarenv schema` (JSON Schema of a scenario, also in
[docs/scenario.schema.json](docs/scenario.schema.json)), `lunarenv data list|verify`.
For agents: `lunarenv-mcp`, see [docs/mcp.md](docs/mcp.md).

## What a result carries

- units and a description on every variable;
- nominal values with 5 / 50 / 95 % quantiles over the uncertainty ensemble, and a
  separate estimate of the numerical error;
- the result of the alternative albedo model next to the nominal one;
- versions of the package and of each model; URL, version and SHA-256 of every input
  file; the software environment;
- a list of warnings whenever an input is outside the range the models were validated for;
- `run_id`, a hash of scenario, data checksums and code versions: the same inputs give
  the same id and, on the same platform, bit-identical results (tested).

## Validation

[VALIDATION.md](VALIDATION.md) lists every case with its source, metric, tolerance and
result. Summary (33 cases): 19 pass, 9 fail, 1 partial, 2 informational, 2 pending. Eighteen
cases (TH-V1b, ILL-V2b, ILL-V5, TH-V3, TH-V3b, TH-V3c, TH-X1b, TH-X1c, TH-A3, TH-C1,
TH-V2c, TH-A2b, TH-C2, TH-V3d, ILL-C1, TH-V3e, TH-V3f, TH-V3g) follow an independent scientific review and were
pre-registered in [docs/review/](docs/review/) before they were run.

| Compared with | Result |
|---|---|
| Closed-form solutions (horizon of a sphere and of an analytic hill, disc segment, damped thermal wave, steady conduction, energy conservation) | pass |
| JPL Horizons, direction of the Sun and the Earth at the site | pass: 9e-10 and 9e-9 rad |
| Barker et al. 2021, Table 2: illumination of seven regions at 1 m and 5 m, 1st percentile over 100 terrain clones | pass: within 2.2 percentage points at 1 m and 0.3 at 5 m (tolerance 3); region areas within 2.6 % |
| PGDA 18.6-year illumination map (method of Mazarico et al. 2011) | pass: bias -1.6 pp, RMS 4.1 pp, r = 0.96; about 6 pp darker on well-lit cells, mostly a metric mismatch: with the map's own metric (any part of the disc visible) bias +0.5 pp and -1.8 pp on well-lit cells |
| Gläser et al. 2018, named best-lit points over 20 years | 2 m: pass at all three points (within 2.1 pp; longest shadow 116 / 66 / 65 h against 112 / 66 / 65 h). Ground level: **fail** (37-44 % against 76-81 %), and strongly dependent on terrain resolution |
| Regional terrain error, PGDA 90 rms height error map, regional DEM raised and lowered by 2 errors (ILL-C1) | sensitivity: at most 0.013 pp on the yearly mean visible disc and 0.003 degrees on the horizon at seven points; no longest shadow changes. The site-terrain clones spread at least 30 times wider |
| Diviner equatorial temperatures quoted by Hayne et al. 2017 | pass: 385.2 / 100.3 / 93.7 K against 385 / 101 / 95 K. A constant albedo passes too (TH-V1b), so this case does not test the albedo law |
| Apollo 15 and 17 mean temperatures quoted by Hayne et al. 2017 | **fail** with the published parameters (8-10 K low); pass with a constant albedo |
| Apollo 17 like for like (TH-V2c): Keihm & Langseth 1973, brightness temperature (e_m = 1), site albedo 0.08, measured heat flow | **fail**: surface 207.0 K against 216 +/- 5 K (-9.0 K); 252.7 K at 1.30 m against 256 K (pass). With a constant albedo 0.08: 214.4 / 260.6 K, within 5 K |
| heat1d 0.4.1, the current reference code of Hayne et al. 2017 | converged grid (TH-X1c): **pass**, surface within 0.8 K and depth within 0.8 K at 0, 20 and 26 deg. Default grid (TH-X1b): **fail**, night minimum and mean 1.1-1.6 K warmer, because that grid is not converged. The obsolete 0.3.2 (TH-X1): partial, 2-3 K colder at depth |
| Diviner polar summer temperatures, 720 m block on the Shackleton rim (TH-V3, TH-V3b) | published albedo law: **fail**, 24 K too cold over the sunlit summer bins; constant albedo: pass, +5.7 K (tolerance 10 K). Only 8 bins qualify; single 240 m cells differ by up to 70 K, and both laws are 43-67 K below the hottest Diviner bin |
| Diviner, same block, against the albedo envelope (TH-V3c) | **pass**, at the hot edge: cold law (Bürger 2024) -34.2 K, hot law (Foote 2020) -0.4 K over the sunlit bins. Diviner lies between the two laws in only 3 of 8 sunlit bins and 5 of 64 summer bins; elsewhere it is warmer than the hot law |
| Diviner, same block, with terrain radiation and Earthshine added (TH-V3d) | published albedo law: still **fail**, -23.2 K over the sunlit bins (+0.7 K from TH-V3); hot law +0.1 K. Over all summer bins the diffuse terms explain +0.7 K of the -39 K. Absorbed terrain radiation 0.77 W m-2 on average over the summer (up to 4.0 W m-2 at one receiver); cells change by 0.1-3.7 K |
| Diviner at four more polar sites, 720 m blocks chosen from model geometry (TH-V3e: PGDA 78 Sites 1, 7, 11, 23) | **fail**: hot law -22.0 / -12.7 / -0.6 / -7.3 K over 28 / 26 / 14 / 28 sunlit bins (tolerance 10 K; passes at Sites 11 and 23). Cold law -39 to -63 K and published law -27 to -51 K at every site. Diviner lies between the two laws in 16 of 96 sunlit bins |
| Diviner at the five polar blocks with sub-DEM roughness added (TH-V3f: Gaussian slopes, Smith 1967 shadowing, after Bandfield et al. 2015) | **fail**: Hayne law with 20 deg RMS roughness -12.8 / -34.4 / -27.0 / -14.1 / -21.1 K (Shackleton, Sites 1, 7, 11, 23; tolerance 10 K), 11-18 K warmer than smooth at every site. 35 deg: -9.6 to -30.1 K. Hot law changes by under 1 K at 20 deg. Equatorial noon -3.2 K |
| Diviner at the five polar blocks with the Diviner footprint and binning applied to the model (TH-V3g: 335 x 170 m IFOV at 50 km, 206 m smear, nine-point binning of the RDR SIS) | **fail**: Hayne law -28.8 / -49.3 / -44.3 / -26.5 / -33.9 K (Shackleton, Sites 1, 7, 11, 23; tolerance 10 K). The footprint changes the block by -5.0 / +2.2 / +0.6 / +0.7 / -0.2 K (hot law -5.3 to +1.3 K): it does not explain the gap. Hot law with footprint -20.7 K at Site 1, -12.6 K at Site 7 |
| Spherical bowl crater, single-bounce terrain radiation (TH-A3, Ingersoll et al. 1992) | **pass**: scattered sunlight and terrain infrared within 0.4 % at three shadowed floor points (tolerance 5 %); default patch layout within 3 % |
| Terrain irradiance at the rim minimum against the earlier 0.0217 W m-2 bound and the review's 0.05 W m-2 threshold (TH-C1) | **pass**: 8.5e-5 W m-2 at most over the three albedo laws; Earthshine 0.006 W m-2 absorbed, +0.23 K on the minimum. With twice the patches the terrain value is 4.5e-5 W m-2: not converged in relative terms, negligible in absolute ones |
| Low-temperature conductivity, Woods-Robinson et al. 2019 shape below 150 K (TH-A2b, TH-C2) | analytic steady case with the option: **pass** (1e-15). Effect at the test point (sensitivity, no pass/fail): minimum -3.3 K, mean -1.6 K, mean at 0.3 m +5.2 K, at 1 m +7.8 K; equatorial night -0.6 K (still within TH-V1's 2 K) |
| Mazarico 2011 site table, Gläser 2014 | pending: references not obtained |

## Known limitations

Scope

- **One site.** Terrain sets exist for the Shackleton rim (`site04`) and the connecting
  ridge (`site01`, used for validation). Observers must lie inside the 16 km site tile.
- **Two models.** No structural, seismic, radiation, dust or ejecta model. The thermal
  model gives the temperature of the *ground*, not of a structure.

Illumination

- The terrain is the 5 m LOLA site model, about 90 % of whose pixels are interpolated
  between altimeter tracks. Relief below that scale (boulders, small craters) is absent.
  At ground level this matters; from a metre or more above the ground much less.
  Against Gläser et al. 2018, ground-level values at two ridge points are 37-39 pp
  below the published ones and move by up to 30 pp when the site DEM is coarsened from
  5 m to 20 m, while 2 m values agree within 2.1 pp. Do not use 0 m illumination for
  design.
- Terrain uncertainty comes from the published error clones of the site tile. The
  regional terrain has no clones; runs bound its effect with the published error map
  instead (regional DEM raised and lowered by twice its rms error). At the test point and
  the Gläser 2018 points this bound is at most 0.013 percentage points wide (ILL-C1).
- "Lit" means any visible fraction of the solar disc; longest-period statistics are
  sensitive to that definition and are resolved to one time step.
- Numerical error of the yearly mean illumination: up to +0.12 percentage points,
  biased towards more light.

Thermal

- **The polar thermal results are outside the validated range**, and say so in every
  run. The regolith properties were fitted to data within 70 degrees of the equator;
  the heat-capacity fit is stated valid from about 90 K, and polar ground gets colder
  (the low-temperature Woods-Robinson 2019 heat capacity used by NASA's DSNE is an
  option and moves the results by less than 0.3 K). Conductivity below 150 K is not
  known for lunar regolith; a Woods-Robinson 2019 low-temperature shape is an option, and
  every run reports its effect: at the test point the minimum drops by 3.3 K and the mean
  at 1 m rises by 7.8 K (TH-C2). Treat temperatures at depth as bounded by both laws.
- **The albedo law at grazing incidence is not known, so polar temperatures are given
  as bounds.** Published highland laws give 0.24 to 0.86 at the test point's 88 degree
  incidence. Every run reports the surface temperature under the cold (Bürger 2024) and
  hot (Foote 2020) laws, each with its ensemble, and the bounds over both: at the test
  point the yearly maximum lies between 103 and 200 K. Design with the bounds, not the
  nominal Hayne law. Against Diviner summer data on a 720 m block (TH-V3, TH-V3c) the
  Hayne law is 24 K too cold, the hot bound 0.4 K too cold and a constant albedo 0.12
  (outside the envelope) 6 K too warm; Diviner sits at the hot edge of the envelope,
  above it in most individual bins. The evidence is thin (8 sunlit bins, one site),
  single 240 m cells disagree by up to 70 K. Adding terrain radiation and Earthshine to
  the comparison (TH-V3d) moves it by only +0.7 K, so they do not explain the gap. At
  four more polar sites (TH-V3e) the Hayne law is 27-51 K too cold and the cold bound
  39-63 K too cold everywhere, and **Diviner is warmer than the hot bound** by 22 K
  (Site 1, connecting ridge) and 13 K (Site 7) over the sunlit summer bins. So the hot
  bound is not an upper bound for sunlit summer ground at 720 m scale. Sub-DEM roughness
  of the size Diviner measures at the equator (20 deg RMS, TH-V3f) warms the Hayne law
  by 11-18 K at all five blocks but leaves it 13-34 K too cold, and barely moves the hot
  bound (under 1 K), so it closes part of the Hayne gap and none of the hot-bound gap at
  Sites 1 and 7. The Diviner footprint and binning (TH-V3g) change the blocks by -5 to
  +2 K, so they explain neither gap; observation timing and the albedo law remain.
- **Slope dominates and is poorly known.** With the Sun never more than 1.8 degrees
  high, a tilt of one degree changes the absorbed flux severalfold. The slope of a 5 m
  terrain model is used unless the scenario gives one; a slope grid (`slope_grid`)
  reports the temperature over slopes and aspects under each albedo law (at the test
  point, 0-10 degrees: maximum 102 to 253 K).
- **Ground temperature is not the temperature of a structure.** For asset surfaces
  (`oriented_surfaces`, any tilt up to vertical) the engine gives the direct solar
  irradiance, the sky and ground view factors and the diffuse irradiance (terrain
  scattered, terrain infrared, Earthshine), the inputs of the asset's own thermal model
  (HLS-UG-001 section 9). At the test point a vertical wall receives up to 1400 W m-2 of
  direct sunlight and on average about 6 W m-2 scattered and 2 W m-2 infrared from the
  terrain; level ground 45 W m-2 direct. The terrain under the asset is the scenario
  albedo law, and the asset does not shadow it.
- **Terrain radiation is single bounce.** Terrain lit by the Sun scatters (Lambertian,
  its albedo law) and emits (its own 1-D temperature), grouped into patches; light
  bouncing more than once, roughness below the 5 m terrain model and non-Lambertian
  scattering are not modelled. The method reproduces the analytic bowl crater within
  0.4 % (3 % at the default patch layout, TH-A3). On the rim test point the ground's own
  slope hides almost all terrain: terrain raises the minimum by 0.01 K and Earthshine by
  0.23 K (TH-C1). Inside craters, where multiple scattering matters, results are lower
  bounds. Earthshine uses the peak values of HLS-UG-001 section 11.1.
- One-dimensional; no rocks; no volatiles.

Software

- Results are reproducible bit for bit on one platform with the locked dependencies;
  across platforms the last bits may differ. The environment is recorded with each run.
- Validation results are stored with the package and are not recomputed automatically
  when code changes.
- Maps cost about 0.1 CPU-second per pixel (a 2 x 2 km window at 5 m: about 4.7
  CPU-hours).
- Continuous integration (GitHub Actions, Linux) runs lint, type checks, the unit tests
  and the tests that only need SPICE kernels. Tests that need terrain data, including
  the end-to-end runs and the published-data comparisons, run only locally.

## Performance

Measured in CPU time on an Apple M1 Pro. Wall-clock figures were taken while other
workloads loaded the machine and are upper bounds.

| Operation | Cost |
|---|---|
| Horizon of one observer, 1440 directions | 0.10-0.13 CPU-s |
| One year of hourly illumination at a point | 8 ms |
| One year of hourly temperature at a point, with spin-up | 0.6 CPU-s |
| Reference scenario: one year, 16 terrain clones, 64 parameter samples, six asset surfaces | 563 s wall-clock (932 s with another job sharing the machine, about 530 s of it for terrain patches) |
| Terrain radiation at a ground point, one year hourly | 1-2 CPU-s per patch point (two-year run, 3 points per patch); at the test point the ground and its slope grid used 131 patches, the asset walls 317 |
| Map, per pixel, 1440 directions | about 0.1 CPU-s |

Details: [docs/models/illumination.md](docs/models/illumination.md),
[docs/models/thermal.md](docs/models/thermal.md).

## Documentation

| Document | Content |
|---|---|
| [PLAN.md](PLAN.md) | Architecture, design decisions with their alternatives, and a self-review of every phase |
| [VALIDATION.md](VALIDATION.md) | Every validation case with source, tolerance and result |
| [docs/data.md](docs/data.md) | Data sources, integrity model, conventions, how large the terrain model must be |
| [docs/models/illumination.md](docs/models/illumination.md), [docs/models/thermal.md](docs/models/thermal.md) | Method, assumptions, uncertainty, performance |
| [docs/convergence/](docs/convergence/) | Convergence and sensitivity studies behind every numerical default |
| [docs/mcp.md](docs/mcp.md) | MCP server and client configuration |
| [docs/references.md](docs/references.md) | References, with which ones were actually read |
| [DESIGN.md](DESIGN.md) | Design system of the HTML report |

## Development

```bash
uv run pytest
```

```bash
uv run ruff check . && uv run mypy
```

Tests marked `data` use the local data cache (`$LUNARENV_DATA_DIR`, or `.data/` in the
repository) and are skipped if it is empty. Studies and benchmarks are in `benchmarks/`.

## Data and citations

lunarenv downloads, and does not redistribute, products of NASA's Planetary Geodesy
Data Archive (LOLA elevation models and illumination maps) and of NAIF (SPICE kernels).
If you publish results, cite the data products listed in
[docs/references.md](docs/references.md), in particular Barker et al. (2021, 2023) for
the terrain and Hayne et al. (2017) for the regolith model.

## License

Apache-2.0. See [LICENSE](LICENSE).
