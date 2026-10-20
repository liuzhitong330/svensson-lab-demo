# Validation record — precursor evidence explorer

Checked 20 October 2026. These are computational and presentation checks,
not biological validation or endorsement by either study's authors.

## Data and arithmetic

- Raw matrix: 54,309 genes × 18 deposited specimens, nonnegative integer counts.
- Late-time-point analysis: 11 independent male mice, LFD 4 / HFD 4 /
  HFD + semaglutide 3. Earlier specimens are not pooled or paired.
- 333 precursors have supplied fragments with actual sequence length 5–25 aa;
  255 pass the full-protein sequence mapping audit. Failures remain in
  `mapping-audit.tsv`.
- 182 of the 255 precursors meet the coverage rule. Of these, 49 have at least
  one identical human/mouse peptide with predicted processing boundaries in
  both species. The default effect-size and leave-one-out rules retain 25.
- Nine analysis regression/invariant tests pass. A separate implementation
  checks normalized sample values, group means, log ratios, omission ranges
  and all 420 projected peptide sequences across the 255 candidates.
- Independent cached-input reproduction is byte-identical for the main JSON
  and both TSV outputs.
- Input SHA-256 values are recorded in `peptide-cancer-provenance.json`.
  The 292 cached inputs were unchanged at final analysis. Annotation endpoints
  are retrieval snapshots; retain the source cache for exact reproduction.

## Checked examples and useful uncertainty

- Coch: descriptive log2 SEMA/HFD contrast +3.2273; leave-one-animal-out range
  +1.6756 to +4.6026.
- Ntn1: +2.1472; omission range +1.7431 to +2.3949.
- Nrcam: −2.0905; omission range −2.4938 to −1.6311.
- These ranges describe single-animal influence, not confidence intervals.
  A narrower range does not establish biological reproducibility.
- Aggf1, Nucb2 and Wnt4 reverse contrast direction under the qualitative CPM
  sensitivity check. Requiring CPM agreement changes the default count from
  25 to 22. The +1 pseudocount has different units in the two normalizations.
- BRINP2 has low precursor RNA in this dataset; the human BRP sequence
  THRILRRLFNLC and mouse MHRIVRRLFNLC differ. It is an audit reference, not
  a nominated tumor-response hit.
- A protein's established function is not assigned to its predicted fragment.
  The supplied library includes membrane and other unconventional precursors;
  the displayed entries are not 255 validated secreted hormones.

## Interface checks

- Desktop first viewport compared with the Ramani reference at 1280 × 720.
- Editorial shell validation passed.
- Responsive browser override tested at 390 × 844; document and scroll width
  both 390 px. This is not testing on a physical phone.
- Candidate, peptide and evidence-rule changes update the visible evidence.
- Three-candidate comparison renders individual omission ranges and prevents
  adding a fourth candidate.
- The actual CSV download contains all qualifying candidates, their sequence
  coordinates, filters, normalization sensitivity and next measurement.
- Frontend regression checks cover all 255 candidate views and 6,120
  candidate/rule combinations, including finite SVG coordinates.

## Interpretation guardrails

RNA abundance does not measure peptide concentration, processing, secretion,
cellular source or anticancer activity. Semaglutide exposure began before
engraftment of pancreatic neoplastic lesions. Low-fat animals also bear these
lesions; they are not healthy or pair-fed controls. Bulk cellular composition
and host physiology remain alternative explanations. No p-values, activity
score, inferred receptor or causal mediation claim is supplied.

The previous concept's QA record is preserved in
`analysis/legacy-brp-controls-QA.md` and does not describe this interface.
