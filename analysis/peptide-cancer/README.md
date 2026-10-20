# Peptide-precursor evidence in a semaglutide-treated pancreatic-lesion dataset

This exploratory analysis connects the Svensson Lab's **published predicted peptide fragments** to an independently deposited mouse RNA-seq dataset. It does not measure peptide abundance, identify a receptor, infer cancer efficacy, or discover a new hormone.

## Reproduce

Use Python 3.12 or later with `numpy`, `pandas`, `openpyxl`, and `regex`. Create an isolated environment if these are unavailable. Run from the repository root:

```sh
python -m venv /path/to/private/analysis-venv
/path/to/private/analysis-venv/bin/pip install numpy pandas openpyxl regex
/path/to/private/analysis-venv/bin/python analysis/peptide-cancer/fetch_inputs.py --cache /path/to/private/source-cache
/path/to/private/analysis-venv/bin/python analysis/peptide-cancer/analyze.py --cache /path/to/private/source-cache --output peptide-cancer-data.json
/path/to/private/analysis-venv/bin/python analysis/peptide-cancer/test_analysis.py --data peptide-cancer-data.json
```

The cache contains downloaded workbooks, FASTA, upstream code and orthology alignments. **Do not publish the cache.** The Nature supplementary tables and the upstream repository are not licensed here for wholesale redistribution. The public outputs are our derived evidence rows, small sequence facts and reproducible analysis code.

Pinned upstream commit: `c478c45405adce4c62baa0719891c3c89f6edcda`. `FastaParser.py` and `secreted.fasta` must match fixed SHA256 values before execution. The adapter executes only the inspected cleavage-pattern expression and pure `get_cleavage_sites` method; it does not execute the full module, web-scraping logic, constructor or FASTA writer. Default upstream minimum spacing is **5 aa**, not a silent substitution for the paper's screening setting.

Ensembl annotation is a retrieval snapshot, not a pinned release API. The provenance output records all input hashes, including every cached protein alignment. Preserve the cache for exact reproduction. A new download may change mappings; compare hashes and mapping-failure counts rather than assuming identical annotation.

## Inputs and unit of analysis

- Coassolo et al., *Nature* 2025, DOI `10.1038/s41586-025-08683-y`, Supplementary Table 2. Its 2,672 peptide–protein rows cover 373 precursors. The manuscript reports 2,683 peptides; we use the actual supplied rows, without inventing the difference.
- The pinned Svensson `pro-hormone-predictor` repository's 2,082 human secreted-protein sequences.
- Ensembl BioMart exact reviewed UniProt accession→human gene mapping, one-to-one mouse orthology, then REST full aligned human/mouse proteins.
- GSE266899 integer raw gene counts and individual GSM metadata; Cencioni et al., DOI `10.1186/s13046-024-03263-w`. This is **external pancreatic-neoplastic-lesion data**, not a Svensson-lab experiment.

The raw matrix contains 54,309 genes and 18 specimens. Analysis uses only the 11 T2 (90-day) male mouse samples: low-fat diet B1/B2/B4/B5; high-fat diet D1/D2/D3/D5; high-fat diet plus semaglutide F2/F4/F5. Each is an independent mouse lesion. T1 is excluded, not paired to T2. The protocol mentions a planned four mice per T2 group, but only three semaglutide samples are deposited. No absent specimen is imputed.

## Predeclared selection rules

1. Start with every supplied fragment whose **actual displayed sequence** is 5–25 amino acids. Reported length differs from displayed sequence length in many rows (mostly by two); the pre-/post-cleavage convention is unresolved. We retain the original length and use verified coordinates, not the length column, for projection.
2. Require a unique exact UniProt→human-gene mapping, a unique one-to-one mouse ortholog and presence in the count matrix.
3. Require the entire ungapped human protein in the current Ensembl alignment to equal the pinned precursor FASTA. A canonical-isoform mismatch is excluded rather than patched by a guessed sequence.
4. Project each uniquely located peptide through that actual protein alignment. Repeated motifs are flagged ambiguous. Keep aligned sequences and insertion/deletion flags.
5. Assess exact local peptide identity separately from cleavage conservation. Both flanking dibasic sites must be recognized in both full species-specific proteins by the pinned upstream cleavage rule. Whole-protein identity is not substituted for local peptide identity.
6. A precursor is measurable if at least three of the seven main-comparison mice have raw count ≥10. The default shortlist requires at least one peptide with exact identity **and** conserved boundaries; same-direction leave-one-out effects; and absolute descriptive log2 ratio ≥0.5. The effect threshold is adjustable and is **not significance**.
7. Sort eligible candidates before audit entries, stable before unstable, then by absolute descriptive effect; break ties with mouse Ensembl ID. This is transparent experimental-priority ordering, not an anticancer or bioactivity score.

The complete sequence-verified precursor universe remains in JSON and TSV, including low-expression and sequence-mismatch audit entries. `mapping-audit.tsv` records failures, not just successful joins. Empty selections are valid results.

## Numerical analysis

Median-ratio size factors use the geometric mean of each gene's counts across all 11 T2 mice, restricted to genes positive in all 11. Each sample factor is the median of its gene/geometric-mean ratios. No additional recentering or length correction is performed.

Normalized count = raw count / sample size factor.

Descriptive effect = `log2[(mean normalized semaglutide count + 1)/(mean normalized high-fat count + 1)]`.

For each of the seven HFD/SEMA mice, omit that mouse, recompute size factors on the remaining ten T2 mice, then recompute the group effect. A stable direction means every omitted-mouse effect has the full-data effect's sign. The minimum/maximum is an **influence range**, not a confidence interval or rank-stability statistic.

The low-fat comparison is contextual. It is not a pair-fed, body-weight-matched or semaglutide-free weight-loss control. A CPM effect is also supplied as a qualitative normalization-sensitivity check; its +1 is in CPM units, so it is not numerically equivalent at low abundance.

## Interpretation and next experiment

These observations can help choose which precursor transcripts, species-matched sequences and processing assumptions to validate first. They cannot show that the peptide exists extracellularly, that the transcript originates in cancer cells, or that the peptide mediates a semaglutide response. Cell-composition changes, stromal remodeling, small group sizes, metabolic effects and bulk-RNA confounding remain.

Validate the candidate's cellular source and precursor expression first; then confirm the actual processed peptide with an appropriate targeted protein/peptide assay. Only after confirming processing and extracellular presence should a species-matched synthetic peptide be prioritized for a controlled cellular bioactivity assay. A positive result would justify mechanism experiments, not establish anticancer activity.

BRINP2/BRP is an explicit audit reference: its human and mouse short sequences differ, and its precursor counts are low in this dataset. Its prior anti-obesity biology does not make it a tumor hit.
