# Scientific review brief: is the approach of lunarenv sound enough to build on?

This brief is for an independent reviewer (Claude Science). It asks for a review of the
**scientific approach**, not of the code. The question to answer is: does the way this
project frames and solves the problem have a sound scientific basis, and should it go
ahead as it is, go ahead with changes, or be rethought?

## 1. What the project claims to be

lunarenv is an open engine that, for a point near the lunar south pole and a period of
time, computes:

- **illumination**: terrain horizon from LOLA elevation models, visible fraction of the
  solar disc, Earth visibility, eclipses;
- **regolith temperature** at the surface and at depth, with a one-dimensional heat
  conduction model using the property laws of Hayne et al. (2017);
- with every number carrying units, uncertainty, model and data versions and a warning
  when inputs are outside the validated range.

Intended users: **engineers designing structures for the lunar surface** who need
environmental conditions they can audit. The stated ambition is "open and verified".

Current state (MVP, one site: Shackleton crater rim): 16 validation cases, of which
11 pass, 1 fails, 1 is partial and 3 are pending for lack of references. At the test
point the yearly mean visible solar disc is 85.5 % and the surface temperature ranges
from 43 to 119 K.

## 2. Read first

All paths are relative to the repository root.

| File | Why |
|---|---|
| `README.md`, section "Known limitations" | What the author already admits |
| `PLAN.md` sections 1-3, 6-9 | Scope, research findings, design decisions with alternatives |
| `PLAN.md` sections 12-14 | The author's own sceptical review of each phase and the open list |
| `VALIDATION.md` | Every validation case: source, metric, tolerance, result |
| `docs/models/illumination.md`, `docs/models/thermal.md` | Method, assumptions, validity range |
| `docs/convergence/*.md` | Numerical evidence behind defaults (read only if a question needs it) |
| `docs/references.md` | Sources, marked [A] read in full or [S] abstract/secondary only |

Do **not** review the source code, style, performance or software architecture.

## 3. Questions

Answer each one with a short conclusion first, then the evidence.

### Q1. Is the problem well posed for its users?

- Is "ground illumination and regolith temperature at a point" the quantity an engineer
  designing a lunar surface structure needs? Or do they need something else first:
  temperature of the structure itself, the full radiative environment (direct, scattered
  and infrared from terrain, Earthshine), sky view factors, extremes or design envelopes
  rather than time series?
- How do existing design references for the lunar environment present these
  quantities? Check what NASA and ESA publish for natural-environment design
  specifications; one candidate to verify is NASA's *Cross-Program Design Specification
  for Natural Environments* (DSNE). Say what they contain and how lunarenv compares.

### Q2. Does it fill a real gap?

- What already exists, open or published, that computes polar lunar illumination and
  regolith temperatures (for example PGDA illumination products, the methods of
  Mazarico et al. 2011, Gläser et al. 2014, Barker et al. 2021, polar thermal models
  built on Diviner data)?
- Is "open, verified, traceable, with uncertainty" a meaningful contribution over those,
  or is it repackaging? Be specific about what is new, if anything.

### Q3. Are the modelling choices scientifically defensible at the pole?

Assess each choice against the literature, and say which ones are **sound**, which are
**acceptable with stated limits**, and which **undermine the results**:

1. Horizon method with nested DEMs (5 m near, 80 m and 240 m far), as in Barker 2021.
2. Using the Hayne 2017 regolith property laws, fitted on Diviner data within about
   70 degrees of the equator, at 89 degrees south.
3. **Incidence-dependent albedo law A(θ) at grazing incidence (θ ≈ 88°).** This is the
   largest open problem: it gives A = 0.79, while a constant albedo raises the yearly
   mean temperature by 26 K and the maximum by 50 K. With the published law the model
   matches Diviner equatorial extremes (TH-V1) but is 8-10 K below the Apollo mean
   temperatures (TH-V2, fail); with a constant albedo it matches Apollo (TH-V2b). Which
   is physically justified at grazing incidence? What do polar thermal studies use?
4. Heat-capacity fit stated valid from about 90 K, used down to 43 K. How large is the
   likely error? Is there a better published low-temperature property set?
5. Neglect of terrain-scattered sunlight and terrain thermal emission (results in shadow
   are flagged as lower bounds; a bound of 0.02 W m⁻² is reported at the test point).
   Is neglecting them acceptable for an MVP, or does it make polar temperatures
   meaningless near crater rims?
6. One-dimensional conduction, no rocks, no volatiles; slope from a 5 m DEM of which
   about 90 % of pixels are interpolated, while slope is the leading sensitivity.

### Q4. Is the validation strategy adequate to call the results "verified"?

- All thermal comparisons are equatorial or mid-latitude and terrain-free. Is there
  **any** polar or high-latitude temperature measurement that could validate the thermal
  model? Candidates to check: Diviner polar products (Williams et al. 2019; archive path
  not yet found by the author), in-situ regolith temperature measurements by landers at
  high southern latitudes (for example ChaSTE on Chandrayaan-3). Give what exists, where,
  and whether it is comparable.
- Illumination agrees with Barker 2021 within 2.2 percentage points and with the PGDA
  18.6-year map with a bias of -1.6 pp, but 6 pp darker on well-lit cells (the author
  attributes it to observer height or slope sampling, not established). Is that
  agreement good enough for the stated use?
- Two tolerances (TH-V1, TH-X1) were written after exploratory numbers were seen. Does
  that weaken the conclusions?
- Can the three pending cases be unblocked? Try to obtain Mazarico 2011 Table 3,
  Gläser 2014 (location of its best-lit point), Vasavada 2012, Keihm 1973 / Langseth
  1976 (what the Apollo "mean surface temperature" actually measures). Report what you
  could access and what it changes.

### Q5. Is the uncertainty honest?

- Uncertainty combines 100 terrain error clones and a Latin hypercube over regolith
  parameters; model-form alternatives (albedo law) are reported separately, not inside
  the ensemble. Two parameter ranges (surface density, heat flow) are assumptions; one
  parameter (χ) is not varied.
- Given that the albedo choice (26 K) and slope dominate everything sampled, should the
  thermal output be presented as a prediction at all, or as a bounded range or a
  scenario comparison?

### Q6. Verdict

- **Go / go with conditions / rethink**, with the reason in two or three sentences.
- The **three changes** that would most increase scientific credibility, in order, each
  with the evidence that would show it worked.
- What the minimum is before an engineer should use a temperature from this engine in a
  design decision.

## 4. Rules

1. **Do not invent data, values or references.** Every factual claim names its source.
   Mark each source **[A]** if you read the relevant part yourself, **[S]** if only the
   abstract or a secondary source. If you cannot access something, say so; do not
   reconstruct it.
2. Separate **established** facts from **hypotheses**, including the author's: the
   self-review in PLAN.md contains several stated hypotheses; test them, do not adopt
   them.
3. Be quantitative where you can: "changes the mean by about X K" beats "matters".
4. Prefer checks that can fail: when you recommend something, say what result would
   show the recommendation was wrong.
5. Do not change any file other than the deliverable, and do not propose loosening a
   tolerance to make a case pass.

## 5. Deliverable

One Markdown file, `docs/review/scientific-review-YYYY-MM-DD.md`, in English:

1. **Verdict** (Q6), at most 10 lines.
2. **Findings table**: ID, question (Q1-Q5), finding, severity, evidence and source
   [A]/[S], recommended action.
   Severity:
   - **Critical**: the approach cannot support its stated use without this change.
   - **Major**: results are presented with more confidence than the evidence supports.
   - **Minor**: a gap in justification or documentation that does not change results.
3. **Answers to Q1-Q5**, short conclusion first.
4. **Sources consulted**, with [A]/[S] and what each was used for, plus a list of
   sources you tried and could not access.
