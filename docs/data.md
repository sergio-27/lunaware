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
| `pgda78.site07.surf`, `.site11.surf`, `.site23.surf` | PGDA-78 Sites 7, 11, 23 | 5 m | Terrain sets `site07`, `site11`, `site23`; validation case TH-V3e. No clones in the catalog |
| `pgda78.<name>.surf` for the other 22 sites | PGDA-78 | 5 m; 3200 to 6000 pixels on a side | Terrain sets of the same name (`lunarenv sites`). Pinned 2026-11-18. No clones in the catalog |
| `pgda90.ldem_80s_80mpp` | PGDA-90 | 7600 x 7600 at 80 m, +/-304 km | Mid- and far-field terrain |
| `pgda90.ldem_80s_80mpp_err` | PGDA-90 | same grid | RMS height error of the 80 m DEM (m); regional error bound (optional, group `regional-error`) |
| `pgda90.ldem_60s_240mpp` | PGDA-90 | 7760 x 7760 at 240 m, +/-931 km | Terrain beyond the 80S product, for observers it is not provably large enough for (group `regional-wide`) |
| `pgda69.avgvisib_85s_60m*` | PGDA-69 | 60 m | Published illumination and Earth-visibility maps (validation) |
| `naif.*` | NAIF generic kernels | - | Time, ephemerides, lunar orientation, frames |

The full list with URLs, versions and citations is `src/lunarenv/data/catalog.toml`.

### Site tiles (PGDA 78)

Read from the site DEMs on 2026-11-18 (all at 5 m per pixel, no missing heights). Ranges
are the outer edges of the tile in south polar stereographic kilometres; heights are above
the 1737.4 km sphere.

| Terrain set | Site of the product | Tile | x range (km) | y range (km) | Centre (lat, east lon) | Height range (m) |
|---|---|---|---|---|---|---|
| `site01` | Site01: Connecting ridge | 3200 x 3200 | -19 to -3 | -20 to -4 | -89.46, 222.5 | -523 to 1959 |
| `site04` | Site04: Shackleton rim | 3200 x 3200 | -9 to 7 | -15 to 1 | -89.77, 188.1 | -2848 to 1805 |
| `site06` | Site06: Nobile rim 1 | 4000 x 4000 | 74 to 94 | 100 to 120 | -85.44, 37.4 | 860 to 6800 |
| `site07` | Site07: Peak near Shackleton | 3200 x 3200 | 22 to 38 | -28 to -12 | -88.81, 123.7 | -419 to 1649 |
| `site11` | Site11: de Gerlache rim | 3200 x 3200 | -45 to -29 | 7 to 23 | -88.68, 292.1 | -805 to 1816 |
| `site20` | Site20: Leibnitz beta plateau | 3200 x 3200 | 65 to 81 | 110 to 126 | -85.43, 31.7 | 2796 to 6464 |
| `site20v2` | Site20v2: Leibnitz beta plateau, extended boundaries | 5200 x 5200 | 60 to 86 | 105 to 131 | -85.43, 31.7 | 1330 to 6788 |
| `site23` | Site23: Malapert massif | 4200 x 4200 | -11 to 10 | 111 to 132 | -85.99, 359.8 | 902 to 5174 |
| `site42` | Site42: de Gerlache-Kocher massif | 4000 x 4000 | -123.4 to -103.4 | -66.1 to -46.1 | -85.83, 243.7 | -1093 to 2998 |
| `haworth` | Haworth | 5960 x 5960 | -52.9 to -23.1 | 75.6 to 105.4 | -86.76, 337.2 | -1643 to 2842 |
| `shoemaker` | Shoemaker | 4000 x 4000 | 66 to 86 | 29 to 49 | -87.18, 62.8 | -1075 to 714 |
| `dm1` | DM1: Amundsen rim | 4000 x 4000 | 154 to 174 | 51.5 to 71.5 | -84.23, 69.4 | -1839 to 1020 |
| `dm2` | DM2: Nobile rim 2 | 4000 x 4000 | 147 to 167 | 85 to 105 | -83.95, 58.8 | -1430 to 1818 |
| `sl2` | SL2: de Gerlache rim | 4000 x 4000 | -58.5 to -38.5 | 13 to 33 | -88.23, 295.4 | -1860 to 1765 |
| `sl3` | SL3: Connecting ridge extension | 4000 x 4000 | -38.7 to -18.7 | -16.1 to 3.9 | -89.03, 258.0 | -1170 to 1693 |
| `npa` | NPA: Cabeus exterior wall 1 | 4000 x 4000 | -98.5 to -78.5 | 19.5 to 39.5 | -86.92, 288.4 | -2612 to -1347 |
| `npb` | NPB: Amundsen 1 | 4000 x 4000 | 175 to 195 | -8.5 to 11.5 | -83.90, 89.5 | -3891 to -2859 |
| `npc` | NPC: Idel'son L crater 1 | 4000 x 4000 | 122.5 to 142.5 | -109.5 to -89.5 | -84.54, 126.9 | -3397 to -1948 |
| `npd` | NPD: Malapert crater 1 | 4000 x 4000 | 13 to 33 | 146 to 166 | -84.80, 8.4 | -2463 to -1337 |
| `lm1` | LM1: Shackleton Rim B | 4000 x 4000 | 2 to 22 | -2 to 18 | -89.52, 56.3 | -2084 to 1444 |
| `lm2` | LM2: Shoemaker Rim A | 4000 x 4000 | 2 to 22 | 47 to 67 | -88.08, 11.9 | -3907 to -153 |
| `lm3` | LM3: Shoemaker Rim B | 4000 x 4000 | 4 to 24 | 37 to 57 | -88.38, 16.6 | -2884 to 71 |
| `lm4` | LM4: Shoemaker Rim C | 4000 x 4000 | 9 to 29 | 15 to 35 | -88.96, 37.2 | -3792 to -48 |
| `lm5` | LM5: Shoemaker Rim D | 4000 x 4000 | 30 to 50 | 1 to 21 | -88.63, 74.6 | -2985 to -173 |
| `lm6` | LM6: Shoemaker Rim E | 4000 x 4000 | 49 to 69 | 10 to 30 | -87.95, 71.3 | -3721 to -123 |
| `lm7` | LM7: Faustini Rim A | 4000 x 4000 | 54 to 74 | -10 to 10 | -87.89, 90.0 | -2776 to 579 |
| `lm8` | LM8: Shoemaker Rim F | 4000 x 4000 | 60 to 80 | 30 to 50 | -87.34, 60.3 | -2498 to 714 |

The regional DEM of every set is the 80 m product poleward of 80 S (+/-304 km); beyond
it the 240 m product poleward of 60 S is used where needed (next section). A short
illumination run at the centre of each tile completes at all 27.
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
| Wide DEM (60S, 240 m) vs regional (80S, 80 m) block-averaged to 240 m, on the whole 80S area | mean 0.00 m, RMS 2.4 m at zero shift; RMS 33 m when shifted one 240 m cell in any direction | zero shift is the minimum |
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


### Beyond the 80S product: the wide DEM (2026-11-18)

The same test is applied to the observer of every run, with the shortest distance to the
edge of the 80S product (a lower bound: the map distance divided by the largest scale
factor inside the product, 1.0153) and the observer's height. Where it fails, the terrain
stack gets a fourth level, the 240 m product poleward of 60 S (`pgda90.ldem_60s_240mpp`,
+/-931 km), and a ray that leaves the 80S product continues on it. Inside the 80S product
the terrain sampled is the same as before; only the height cap of the stopping rule can
change (the highest point of the stack), which makes rays stop later, never earlier. If
the file is not in the cache the run keeps the three levels and warns.

Evaluated 2 m above the ground on every eighth pixel of the 27 site tiles:

| Where the 80S product is not provably sufficient | Tiles |
|---|---|
| Nowhere at the reference point of `site04` (rim crest); 4.3 % of that tile (the crater floor) | `site04` |
| 12 %, 45 %, 97 % of the tile | `site01`, `lm1`, `sl3` |
| The whole tile | the other 23 |

The test assumes terrain 11 km high right at the edge; the highest point of the 60S
product is 7.0 km, so it asks for the wide DEM far more often than terrain there turns
out to matter. Measured at the centres of seven tiles (`site04`, `site06`, `site20`,
`site23`, `site42`, `npb`, `dm2`; 1440 azimuths, 2 m): the horizon is identical with and
without the wide level at every azimuth. At `site06` and `site42`, 121 rays each ended at
the edge of the 80S product while terrain beyond could still have hidden the Sun; with
the wide level none does, and the 2026 mean visible disc is unchanged (30.1857 % and
79.2568 %, longest shadow 512 h and 214 h). The wide level therefore changed no number
that was checked; it replaces "upper bound" by a result at those points.

Not checked: observers other than tile centres, heights other than 2 m, and the 240 m
resolution of the wide level against a finer product (there is none beyond 80 S in the
catalog). Unit test: a ridge placed beyond the regional DEM is found within 4e-5 rad of
an independent geodesic computation, and azimuths away from it are unchanged to the bit.

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
