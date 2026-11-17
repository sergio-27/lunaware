# Export of a run for another thermal tool

`lunarenv export` writes a finished run as plain files that a thermal analyst can load
into another tool (Thermal Desktop, ESATAN-TMS, a script): where the Sun is and how much
of it the terrain leaves visible, the terrain horizon, the irradiance on each asset
surface and, on request, the terrain as a triangle mesh.

**Status.** The files are neutral (CSV, Wavefront OBJ, JSON). No importer for Thermal
Desktop or ESATAN-TMS exists yet, and the files have not been loaded into either tool.

```bash
uv run lunarenv export <run_id>
```

```bash
uv run lunarenv export <run_id> --ring-radius 1000 --mesh-half-width 500 --mesh-stride 2
```

From Python: `lunarenv.export_run(run_directory, output=None, mesh_half_width_m=None,
mesh_stride=1, ring_radius_m=None)` returns the manifest.

The output goes to `<run>/export/` unless `-o` names another directory. The same run and
options give byte-identical files.

## Frame

Every vector and mesh is in the grid-aligned local frame:

- origin on the ground under the observer;
- `+x` along map `+x` and `+y` along map `+y` of the south polar stereographic projection
  of the LOLA products; `+z` up (radially outward); right-handed;
- azimuths are grid azimuths, from `+y` towards `+x`. True azimuth = grid azimuth minus
  the longitude of the site. Grid axes are used because east and north are not defined
  at the pole;
- the observer is at `(0, 0, height_above_surface_m)`.

`manifest.json` gives the three local unit vectors in the body-fixed frame
(`frame.local_to_body_fixed_rows`), so the frame can be related to any other.

## Files

| File | Rows | Content |
|---|---|---|
| `sun_earth.csv` | one per time step | Sun azimuth and elevation, unit vector to the Sun, visible and flux fractions of the disc, Earth eclipse factor, direct normal irradiance, Earth elevation and visible fraction |
| `horizon.csv` | one per azimuth | terrain horizon elevation and the distance to the terrain that defines it |
| `surface_<name>.csv` | one per time step | for each asset surface of the scenario (`thermal.oriented_surfaces`): direct irradiance, the four diffuse components when computed, and their sum |
| `horizon_ring.obj` | `--ring-radius` | a wall around the observer whose top edge is the horizon |
| `terrain.obj` | `--mesh-half-width` | the site DEM as triangles |
| `manifest.json` | | frame, time axis, units and description of every column, SHA-256 of every file, surface normals and view factors, validation status, warnings, input data with checksums |

Time columns: `time_utc` (ISO 8601, UTC) and `elapsed_s` (seconds since the first step).
Numbers have 10 significant digits. The header of each CSV is its first line; units and
descriptions are in the manifest, not in the CSV.

`solar_irradiance_normal_w_m2` already includes the distance to the Sun, the
limb-darkened fraction of the disc above the horizon and the Earth eclipse factor. Do not
multiply it again by `sun_flux_fraction`.

`direct_w_m2` of a surface is that irradiance times the cosine of the incidence on the
surface. The diffuse columns are single bounce (see
[models/thermal.md](models/thermal.md)).

## Meshes

- **`terrain.obj`**: one vertex per pixel centre of the site DEM (every N-th with
  `--mesh-stride N`), two triangles per cell, facing up, in a square of the given half
  width around the site. Vertices are true positions: the curvature of the Moon is kept,
  so the mesh drops below the tangent plane by about d²/2R (0.07 m at 500 m, 18 m at
  8 km), and neighbouring vertices are (R + height) / R times the DEM pixel apart
  (5.004 m for a 5 m pixel at 1.35 km height). Relief below the DEM resolution is absent. Needs the site terrain in the data
  cache.
- **`horizon_ring.obj`**: at each azimuth of the run the top of the wall, seen from the
  observer, has the horizon elevation, so the wall hides the Sun exactly when the terrain
  does. It stands for all terrain beyond its radius (far-field terrain reduced to a
  horizon, HLS-UG-001 section 9.4.3). Faces look inwards.

Use the ring **or** rely on `sun_flux_fraction`, not both: a tool that ray-traces the
ring and is also given the already reduced irradiance counts the terrain shadow twice.

## Limits

- The horizon, the visible fractions and the ring hold at the observer point and height
  only. A surface a few metres away or at another height sees a different horizon; at the
  test site the illumination at ground level differs from that at 2 m by tens of
  percentage points (ILL-V5).
- The asset does not shadow the ground or itself; the diffuse columns do not include
  light scattered by the asset.
- The direction of the Earth is exported as an elevation only (the run does not store
  its azimuth).
- Ground temperatures are **not** exported. The polar thermal results are outside the
  validated range (TH-V3 to TH-V3g in [VALIDATION.md](../VALIDATION.md)); they stay in
  `results.nc` and are listed under `not_exported` in the manifest.
- Verification: unit tests on a synthetic run (tables read back against the run, mesh
  geometry against a sphere, ring against the horizon, and the horizon of the mesh
  against the one the engine traces on the same synthetic terrain). On the reference run
  the horizon ray-marched on `terrain.obj` (+/-1000 m) agrees with `horizon.csv` within
  0.05 degrees (RMS 0.011) at the 136 azimuths whose horizon lies inside the mesh; with
  x and y swapped or one axis flipped the RMS is 2.7 to 5.3 degrees. No comparison with
  another tool's import.
