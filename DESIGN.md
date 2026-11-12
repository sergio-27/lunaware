# Design system: lunarenv run report

Design reading: a technical audit report, dense and scannable, in the manner of
scientific documentation. Monochrome with one data accent; system typeface; compact
tables; light and dark.

## Principles

- The first thing a reader sees is what the run is, whether it is inside the validated
  range, and the warnings. Numbers come second, charts third, provenance last.
- Every nominal value sits next to its uncertainty. No hero numbers without a range.
- Colour is for data and for state, never decoration. State is always icon + label.
- One column, left-aligned, 1040 px at most. No cards, no gradients, no shadows.

## Tokens

| Token | Light | Dark | Use |
|---|---|---|---|
| `--surface` | `#fcfcfb` | `#1a1a19` | page and chart background |
| `--surface-2` | `#f3f3f1` | `#232322` | table header, code, banner |
| `--border` | `#deddd8` | `#3a3a37` | rules, table lines, chart grid |
| `--text` | `#0b0b0b` | `#f0efec` | body, values |
| `--text-2` | `#52514e` | `#c3c2b7` | captions, units, axis labels |
| `--series-1` | `#2a78d6` | `#3987e5` | first data series; links |
| `--series-2` | `#eb6834` | `#d95926` | second data series (the alternative model) |
| `--ramp-1..4` | `#86b6ef #3987e5 #1c5cab #0d366b` | `#9ec5f4 #5598e7 #2a78d6 #184f95` | ordered series (depths) |
| `--neutral` | `#c9c8c2` | `#4a4a46` | terrain silhouette, hidden-Sun marks |
| `--good / --warning / --serious / --critical` | `#0ca30c #fab219 #ec835a #d03b3b` | same | validation state, with icon and label |

Chart colours were checked with the palette validator (categorical pair and ordinal
ramp, light and dark): all pass.

## Type, space, shape

- Typeface: system UI stack; numbers and identifiers in the system monospace.
- Sizes: 28 / 20 / 16 / 14 / 13 px (title, section, body, table, caption). Line height 1.5.
- Spacing scale: 4, 8, 12, 16, 24, 32, 48 px.
- Radius: 4 px on badges, banners and code; none elsewhere.

## Components

- **Identity block**: scenario name, run id, checksum, site, period.
- **Validity banner**: state badge plus the list of warnings; always present, also when
  there are none ("No warnings").
- **Key-values table**: quantity, nominal, 5 %, median, 95 %, unit.
- **Chart**: title, one-sentence caption stating what to read from it, inline SVG that
  follows the theme, legend when there are two or more series. Single axis only.
- **State badge**: glyph + label (PASS, FAIL, PARTIAL, PENDING, NOT RUN).
- **Provenance table**: data id, version, checksum prefix, citation.
- **Disclosure** (`details`): scenario JSON, model descriptions.

## States

- No thermal block: thermal sections are omitted and the identity block says so.
- No uncertainty ensemble: the range columns read "not computed" and a note explains why.
- No warnings: the banner says "No warnings" with the good state.

## Known limits

- Charts are static (no hover values); the numbers are in `results.nc`.
- Dark mode follows the reader's system setting; there is no toggle.
