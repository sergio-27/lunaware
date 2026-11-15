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

Output of the run (16 clones in the cache, 109 s on an M1 Pro under load):

```text
run f93b88484006a71f  ->  runs/f93b88484006a71f
illumination: mean visible solar disc 85.51 %, longest shadow 63 h, longest illumination 2803 h
              terrain ensemble 5/50/95 %: 85.35 / 85.54 / 85.69 %
thermal: surface 43.3 to 119.2 K (mean 70.3 K) on a 0.39 deg slope
         surface maximum, ensemble 5/50/95 %: 105.6 / 118.4 / 146.8 K
         with constant albedo: mean 96.6 K (+26.3 K); not settled by validation
WARNING thermal: temperatures reach 43.3 K, below the 90 K down to which the heat-capacity fit is stated to be valid
WARNING thermal: regolith properties were derived from data within 70 degrees of the equator; their use at this latitude is an extrapolation
WARNING thermal: 17.7 % of the period has no direct sunlight; terrain-scattered and terrain-emitted radiation are neglected, so temperatures in shadow are lower bounds
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
result. Summary (19 cases): 13 pass, 2 fail, 1 partial, 3 pending. Three cases (TH-V1b,
ILL-V2b, ILL-V5) follow an independent scientific review and were pre-registered in
[docs/review/preregistration-2026-11-13.md](docs/review/preregistration-2026-11-13.md).

| Compared with | Result |
|---|---|
| Closed-form solutions (horizon of a sphere and of an analytic hill, disc segment, damped thermal wave, steady conduction, energy conservation) | pass |
| JPL Horizons, direction of the Sun and the Earth at the site | pass: 9e-10 and 9e-9 rad |
| Barker et al. 2021, Table 2: illumination of seven regions at 1 m and 5 m, 1st percentile over 100 terrain clones | pass: within 2.2 percentage points at 1 m and 0.3 at 5 m (tolerance 3); region areas within 2.6 % |
| PGDA 18.6-year illumination map (method of Mazarico et al. 2011) | pass: bias -1.6 pp, RMS 4.1 pp, r = 0.96; about 6 pp darker on well-lit cells, mostly a metric mismatch: with the map's own metric (any part of the disc visible) bias +0.5 pp and -1.8 pp on well-lit cells |
| Gläser et al. 2018, named best-lit points over 20 years | 2 m: pass at all three points (within 2.1 pp; longest shadow 116 / 66 / 65 h against 112 / 66 / 65 h). Ground level: **fail** (37-44 % against 76-81 %), and strongly dependent on terrain resolution |
| Diviner equatorial temperatures quoted by Hayne et al. 2017 | pass: 385.2 / 100.3 / 93.7 K against 385 / 101 / 95 K. A constant albedo passes too (TH-V1b), so this case does not test the albedo law |
| Apollo 15 and 17 mean temperatures quoted by Hayne et al. 2017 | **fail** with the published parameters (8-10 K low); pass with a constant albedo |
| heat1d, the reference code of Hayne et al. 2017 | **partial**: surface within 0.75 K; 2-3 K colder at depth |
| Mazarico 2011 site table, Gläser 2014, Diviner polar maps | pending: references not obtained (the Diviner polar products are archived; the comparison is not yet built) |

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
- Terrain uncertainty comes from the published error clones of the site tile only.
  Errors of the regional terrain are not sampled, so for observers whose skyline is
  distant the reported spread is a lower bound.
- "Lit" means any visible fraction of the solar disc; longest-period statistics are
  sensitive to that definition and are resolved to one time step.
- Numerical error of the yearly mean illumination: up to +0.12 percentage points,
  biased towards more light.

Thermal

- **The polar thermal results are outside the validated range**, and say so in every
  run. The regolith properties were fitted to data within 70 degrees of the equator;
  the heat-capacity fit is stated valid from about 90 K, and polar ground gets colder.
- **The albedo law is unresolved at grazing incidence.** The published
  incidence-dependent albedo and a constant albedo differ by 26 K in the yearly mean at
  the test point. Validation supports each in a different comparison. Both results are
  always reported; treat them as bounds of an open question.
- **Slope dominates and is poorly known.** With the Sun never more than 1.8 degrees
  high, a tilt of one degree changes the absorbed flux severalfold. The slope of a 5 m
  terrain model is used unless the scenario gives one.
- No sunlight scattered or heat radiated by surrounding terrain. Temperatures in shadow
  are lower bounds; an upper bound of the neglected flux is reported (0.02 W m-2 at the
  test point, where it is negligible; inside craters it would dominate).
- One-dimensional; no rocks; no volatiles.
- No thermal comparison with polar measurements exists yet.

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
| Reference scenario: one year, 16 terrain clones, 64 parameter samples | 109 s wall-clock |
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
