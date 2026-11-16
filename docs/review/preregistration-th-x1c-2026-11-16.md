# Pre-registration: TH-X1c, code-to-code comparison with a converged heat1d grid

Written before the reference values at 20 and 26 deg are generated.

Origin: TH-X1b (pre-registered 2026-11-15) failed: lunarenv's night minima and surface
means were 1.1-1.6 K warmer than heat1d 0.4.1 on its default grid. A post-hoc check at the
equator (`benchmarks/heat1d_grid_check.py`) showed that heat1d's default grid (m=10,
n=4, about 20 layers) is not converged and that heat1d on finer grids is within 0.5 K of
lunarenv. TH-X1b stays as recorded; this case tests the explanation at all three
latitudes.

What is already known (not blind): heat1d's values at the equator on the grids m=n=10,
20 and 40, and lunarenv's values at all three latitudes. heat1d's values at 20 and 26 deg
on a refined grid have not been computed.

## Reference

- heat1d 0.4.1, GitHub `phayne/heat1d`, commit
  `5b73a6c1ba5047f59f6210185b37539e46e68136`, isolated environment, as TH-X1b.
- Lunar defaults, default solver, NYEARSEQ = 20, latitudes 0, 20 and 26 deg.
- Grid m = 40, n = 40 (125 layers at the equator). Chosen because at the equator the
  change from m=n=20 to m=n=40 is below 0.1 K at the surface. The m=n=20 grid is also run
  and stored, so convergence can be checked at every latitude.
- Means time-weighted over local time, as TH-X1b. The script
  `benchmarks/make_heat1d_reference.py` writes `heat1d_reference_converged.json`; the
  earlier reference files are not changed.

## Comparison

- lunarenv with heat1d's parameters (scale height, Sun distance and declination at output
  time, day length) and its default numerical settings, as TH-X1 and TH-X1b.
- Checks and tolerance as TH-X1b: surface mean, maximum and minimum, and the mean at
  heat1d's deepest node, at each latitude; 1 K each.
- Status: PASS if all checks pass; PARTIAL if all surface checks pass and a depth check
  does not; FAIL otherwise.
- Informational: the change in heat1d between m=n=20 and m=n=40 at each latitude. If it
  exceeds 0.25 K at the surface anywhere, the reference is reported as not converged
  there, whatever the status.

## Reading

- PASS: the TH-X1b failure is attributed to heat1d's default grid; TH-X1c becomes the
  current code-to-code comparison. TH-X1b stays in the ledger as FAIL.
- Surface checks failing at 20 or 26 deg: the grid explanation is incomplete at those
  latitudes and is reported as such.
