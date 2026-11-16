# Pre-registration: TH-V3c, Diviner polar temperatures inside the albedo envelope

Written before any forward model with the two envelope laws has been run at the site.

Origin: change 2 of `docs/review/scientific-review-2026-11-13.md` (findings F03, F17).
The review proposes to present polar temperatures as a bounded envelope over the albedo
law rather than as one nominal law, and states how to tell whether that worked: "the
TH-V3 cells fall inside the envelope".

## Envelope laws (read from the primary sources on 2026-11-16)

- Cold bound, Bürger et al. (2024), JGR Planets 129, e2023JE008152, Eq. 15 with the
  best fit of section 4.2.2 for the highlands high-conductivity model: A(theta) = A0 +
  a (theta / 45 deg)^3 + b (theta / 90 deg)^8, A0 = 0.09, a = 0.10, b = 0.02. A = 0.860
  at 88.2 deg.
- Hot bound, Foote et al. (2020), Icarus 336, 113456, Eq. 9 (Apollo 16 highland soil,
  solar-spectrum hemispheric albedo, laboratory): A(i) = 1.3299e-5 i^2 - 1.0248e-4 i +
  0.14382, i in degrees. A = 0.238 at 88.2 deg. The measurements reach 75 deg incidence;
  beyond that the fit is extrapolated.
- These are the values the review attributes to the DSNE Rev J highland envelope
  (0.238-0.860 at 88.2 deg). DSNE Rev J itself could not be obtained; the laws are taken
  from the two papers.

## Comparison

- Exactly the TH-V3 protocol (`docs/review/preregistration-th-v3-2026-11-15.md`): same
  Diviner PCP data, same 3 x 3 cells and 20 m pixel sampling, same forward model, same
  sunlit-bin selection (model block visible fraction >= 0.5, all nine cells with data),
  block Tbol as the radiance mean of the cells. The two laws are added to the forward
  model's law list; nothing else changes.
- Check (PASS/FAIL): over the sunlit bins, the mean of (model block Tbol - Diviner block
  Tbol) is at most +10 K for the cold law and at least -10 K for the hot law, i.e.
  Diviner lies inside the envelope within the TH-V3 tolerance of 10 K. PASS if both
  hold; FAIL otherwise.
- Informational: the bias of each law; the number of sunlit bins in which Diviner lies
  between the two laws; the same over all summer bins with data.

## Reading

- PASS: the envelope brackets the only polar data available at the site; it can be
  presented as the bounds of the polar thermal product, with the caveats of TH-V3
  (8 bins, one block, no terrain radiation).
- Diviner warmer than the hot bound by more than 10 K: the envelope is too cold; the
  hot side needs roughness or terrain radiation, not a lower albedo alone.
- Diviner colder than the cold bound by more than 10 K: the envelope is too warm.
