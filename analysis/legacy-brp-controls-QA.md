# Archived validation record — superseded BRP-controls concept

Checked on 20 October 2026. This is a computational and presentation check,
not biological validation or an endorsement by either study's authors.

## Numerical checks

- Eleven Python tests pass across the two analyses.
- A separate audit parsed the original XLSX XML without importing either
  production analysis. All 14 BRP conditions × two cell assays and six feeding
  times matched in observations, missingness, means and derived normalizations.
- That audit also matched all 15 cancer-study groups, every endpoint observation,
  sample size, mean, sample SD, median and quartile, plus all 280 longitudinal
  body-mass observations.
- Source SHA-256 values are pinned in the manifests. Rebuilding the inputs
  produces the same outputs; raw workbooks and upstream code are not published.
- Frontend tests cover 300 peptide/threshold/time combinations and 12 cancer
  contexts, including unavailable values, negative effects, CSV serialization,
  selected-peptide sequence mapping, and rounded chart bounds.

## Checked examples

- L8A residual responses at 30 minutes: NS1 24.1931%, Neuro2a 1.0688%,
  feeding inhibition 14.2857%.
- At the illustrative 25% ceiling, L8A is the only complete-data candidate at
  30 minutes and 1 hour. No candidate remains at 1.5, 2 or 3 hours.
- At a 20% ceiling, no candidate remains at any selectable time.
- R3A is marked for a missing NS1 measurement, not treated as a complete-data
  candidate or imputed as zero.
- The external semaglutide group has a mean within-mouse body-mass change of
  −7.8268% from source week 16 to 22; the comparison with WD is a difference
  of −12.3959 percentage points, not a ratio.
- The primary paper reports Fig. 3 n as individual mice. The immune gate is
  TRP-2-specific CD8+ TIL, not all CD8+ cells. Source IDs cannot match panels.

## Browser checks

- Desktop first viewport compared with the Ramani reference at 1280px.
- Fixed editorial-shell validation passed.
- A 390 × 844 embedded browser viewport checked responsive controls and text,
  with document width and scroll width both 390px. This is a responsive-layout
  test, not a claim of testing on a physical phone.
- Changing peptide, tolerance and feeding time updates the result and shortlist.
  Keyboard selection of a plot point updates the source-sequence mapping.
- Changing cancer endpoint, ICB and comparator updates the observations,
  group contrast and proposed follow-up question. ICB is disabled and reset
  for the body-mass endpoint.
- Both CSV controls produced downloaded files. Console inspection found no
  application errors during the tested interactions.

## Interpretation guardrails

The threshold is an illustrative descriptive choice, not an inactivity or
equivalence test. Reference-mean uncertainty is not propagated. No pairing,
batch assignment, significance test, receptor mechanism, cross-study treatment
effect or clinical recommendation is inferred. BRP and semaglutide remain
separate studies. Proposed experiments are hypotheses requiring validation.
