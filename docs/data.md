# Data layer

Scope: what lunarenv downloads, how it knows the bytes are right, the conventions it
assumes about them, and how large the terrain model must be. Reference keys resolve in
[references.md](references.md). All numbers in the "measured" tables come from the
test-suite (`pytest -s tests/data`) run on 2026-11-12.

## Sources

| Catalog id | Product | Grid | Used for |
|---|---|---|---|
| `pgda78.site04.surf` | PGDA-78 Site 4 (Shackleton rim) | 3200 x 3200 at 5 m, x in [-9, 7] km, y in [-15, 1] km | Near-field terrain of the test site |
| `pgda78.site04.toterr`, `.slperr` | same | same | Formal height and slope uncertainty |
| `pgda78.site04.clone.0001..0100` | same | same | Terrain ensemble for uncertainty propagation |
| `pgda78.site01.*` | PGDA-78 Site 1 (connecting ridge) | 5 m | Validation case ILL-V1 only |
| `pgda90.ldem_80s_80mpp` | PGDA-90 | 7600 x 7600 at 80 m, +/-304 km | Mid- and far-field terrain |
| `pgda90.ldem_60s_240mpp` | PGDA-90 | 240 m | Fallback for rays leaving the 80S product |
| `pgda69.avgvisib_85s_60m*` | PGDA-69 | 60 m | Published illumination and Earth-visibility maps (validation) |
| `naif.*` | NAIF generic kernels | - | Time, ephemerides, lunar orientation, frames |

The full list with URLs, versions and citations is `src/lunarenv/data/catalog.toml`.
The south-pole lies inside the Site 4 tile; the tile is 16 km on a side, not 20 km.

## Integrity and provenance

- **Pins.** `src/lunarenv/data/checksums.json` holds the SHA-256 and size of every
  accepted file. It is written once, when a maintainer first downloads a file with
  `--pin`, and is append-only afterwards. A user download that does not match its pin
  fails. Files without a pin are refused.
- **Three checks per download**: file signature, length against `Content-Length`,
  checksum against the pin. The signature check exists because pgda.gsfc.nasa.gov
  answers a missing file with `200 OK` and an HTML page.
- **Manifest.** `<cache>/manifest.json` records, for each cached file: id, URL, SHA-256,
  size, product, version, citation, download time and the server's `Last-Modified`.
  Results embed the manifest records of the files they used.
- **Derived files** (crops) are registered with the checksum of their parent, so the
  parent can be deleted without losing provenance. This is how the clone ensembles fit
  on a laptop.
- **TLS** is verified against the operating-system trust store.

Commands: `python -m lunarenv.data fetch <group|id>`, `... verify`, `... list`.
Cache location: `$LUNARENV_DATA_DIR`; else `.data/` at the root of the source checkout the
package runs from, if that directory exists (where the tests look); else the per-user
cache directory. `lunarenv validate --run` prints the location it used.

## Conventions, each one asserted in code

| Convention | Where enforced |
|---|---|
| South polar stereographic, sphere of radius 1 737 400 m, true scale at the pole, central meridian 0, metres | `grid_from_dataset` refuses anything else |
| North-up, square pixels, single band | same |
| Pixel registration (edges on the grid, centres half a pixel inside) | GMT `node_offset = 1` tag checked; co-registration test below |
| Heights in metres above the sphere, no no-data | `read_dem`; height statistics pinned per product |
| Longitude 0 along map +y, 90 E along map +x | Agreement with PROJ to 1e-6 m |
| Body-fixed frame `MOON_ME_DE440_ME421` | `KernelSet`; D6 check below |

## Measured consistency

| Check | Result | Criterion |
|---|---|---|
| Projection vs PROJ, 500 random points poleward of 60S | < 1e-6 m | 1e-6 m |
| Site 4 (5 m, 2020) block-averaged vs regional (80 m, 2023) | RMS 0.46 m at zero shift; 21.1 m when shifted one 80 m cell | zero shift is the minimum |
| Error clone 0001 minus nominal Site 4 DEM | mean < 1 cm, RMS 0.47 m | consistent with 0.30-0.50 m median error (Barker2021) |
| Shackleton rim-circle centre vs Zuber2012 Table 1 | 243 m offset | < 500 m |
| Shackleton rim diameter | 21.1 km | 21 km +/- 1 km |
| Shackleton floor height at published centre | within 210 m of -2770 m | published floor relief |
| Shackleton mean rim-crest height | 1249 m | 1.3 km +/- 150 m |
| Shackleton depth, rim to floor | 4018 m | 4.1 km; see note |
| D6: Sun and Earth direction at the site, DE440/ME421 chain vs pure DE421, 2024-2034 | max 3.0e-7 rad | 1e-5 rad, fixed before the run |
| Moon-centre observer approximation vs surface observer | Sun 1.1e-10 rad, Earth 4.7e-7 rad | analytic bounds 1.2e-9 and 5.1e-7 rad |

Note on the depth: Zuber2012 quotes 4.1 +/- 0.05 km. Our estimator (mean of the highest
point per azimuth between 8 and 13 km, minus the median floor height near the centre)
differs from theirs, and gives 82 m less. The check is a georeferencing sanity test,
not a reproduction of their measurement.

## How large must the terrain model be? (decision D4)

Terrain can hide part of the Sun only if it rises above the lowest elevation the solar
disc ever reaches at the site:

    e_low = -(max sub-solar latitude + colatitude + max solar radius + parallax)

With bounds verified against DE440 over 2024-2044 (measured maxima in brackets):
sub-solar latitude 1.60 deg (1.588), solar angular radius 0.2720 deg (0.2717), parallax
1.2e-5 rad (1.18e-5). For the rim point used in the tests, `e_low` is about -2.1 deg.

A terrain point at ground distance `d` must have height

    h_needed(d) = (R + h_obs) cos(e_low) / cos(d/R + e_low) - R

to reach that elevation. Past the grazing point (`d > R |e_low|`, about 64 km) this
grows with distance, so if `h_needed` at the DEM edge exceeds every lunar height (cap:
11 km), nothing beyond the edge can matter.

Applied to the Site 4 tile with the 80S product (edge at least 287 km away):

- provably sufficient for **98.0 %** of the tile, including every observer above
  -1426 m, hence the whole rim and its flanks;
- not proven for the deepest part of Shackleton's interior, where an observer looks up
  at distant terrain. Those points are in permanent shadow in any case; the horizon
  tracer (Phase 2) applies the same bound per ray and reports any ray that leaves the
  DEM while still able to reach the Sun, and the 60S product is the fallback.

This criterion is for the Sun. The Earth ranges over about +/-7 deg of elevation at the
pole, so for Earth visibility the per-ray rule is the only bound.

## Known limitations

- Roughly 90 % of the 5 m pixels of the site products are interpolated between LOLA
  tracks (Barker2021); the clones represent that uncertainty statistically, not the
  true small-scale relief. Boulders and craters below the effective resolution are
  absent.
- The Moon is a sphere here, as in the LOLA products.
- Download speed. pgda.gsfc.nasa.gov delivers roughly 100 kB/s per connection and
  answered `503` when eight files were requested at once. Large files are therefore
  fetched in parallel byte ranges (`--connections`), and `503`/`429` answers are retried
  after a pause. A 100-clone ensemble (4.1 GB) takes hours; clones are pinned as they
  arrive. Status on 2026-11-12: 13 Site 4 clones cached, the rest downloading; Site 1
  not started.
- Manifest locking is advisory and POSIX-only.
