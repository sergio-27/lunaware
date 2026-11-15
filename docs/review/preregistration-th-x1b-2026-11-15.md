# Pre-registration: TH-X1b, code-to-code comparison with the current heat1d (finding F11)

Written before the reference values were generated. lunarenv's own values for these
cases are already known (TH-X1), so the heat1d output is not looked at before this file
is committed.

Origin: finding F11 of `docs/review/scientific-review-2026-11-13.md`. TH-X1 used heat1d
0.3.2 from PyPI. The heat1d README (GitHub, read 2026-11-15) says the PyPI package "is
an obsolete release that is no longer updated" and points to the GitHub version.

## Reference

- heat1d 0.4.1, GitHub `phayne/heat1d`, commit
  `5b73a6c1ba5047f59f6210185b37539e46e68136` (2026-10-26), Python package in
  `python/`, installed in an isolated environment.
- Lunar defaults of that version (`planets.Moon`), its default grid and solver.
  Latitudes 0, 20 and 26 deg, as TH-X1. Equilibration: the version's own setting if it
  has one; if it has a configurable equilibration time, 20 years, as TH-X1. The script
  `benchmarks/make_heat1d_reference.py` is extended to write a second reference file;
  the 0.3.2 file is kept unchanged as history.

## Comparison

- lunarenv run with heat1d's parameters (scale height, Sun distance, declination and day
  length as reported by that heat1d run), exactly as TH-X1 does.
- Checks and tolerance as TH-X1: surface mean, maximum and minimum, and the mean at
  heat1d's deepest node, at each latitude; 1 K each.
- Status: PASS if all checks pass; PARTIAL if all surface checks pass and a depth check
  does not (as TH-X1); FAIL otherwise.

## Reading

- Depth within 1 K: the TH-X1 depth gap came from the obsolete release; TH-X1 is kept as
  history and TH-X1b becomes the current code-to-code comparison.
- Depth gap of about 2-3 K as before: the gap does not come from the release, and the
  existing hypothesis (conductivity averaging in heat1d) stays open.
- Any surface check failing is reported as such and investigated before anything else.
