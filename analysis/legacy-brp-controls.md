# Archived first concept: BRP controls and separate cancer endpoints

This describes the superseded first demo, not the current peptide-precursor
explorer. Its original source data and scripts are retained for traceability.
Relative links and commands below refer to the original repository root.

A small, reproducible research demo by Cathy Liu, focused on the cellular peptide-bioactivity and in-vivo cancer-study analysis duties of Stanford research position 201427. It asks how control selection changes when a peptide's cellular response and a host-level outcome do not tell the same story.

The two analyses remain separate:

1. **Target-lab BRP control audit.** Reanalyse the Svensson Lab's published alanine-substitution experiments to shortlist low-residual-response control candidates across NS1 and Neuro2a FOS assays and a selected feeding time. Test how the shortlist changes with the descriptive threshold and time point.
2. **External semaglutide study.** Compare body-mass change, day-16 B16 tumor volume, and antigen-specific T-cell function within the treatment contexts of Piening et al. Changing the endpoint changes the next control question; it does not establish a mechanism.

BRP is a non-incretin peptide. Its measurements and mutants are **not** semaglutide measurements or controls. This demo does not establish receptor specificity, reproduce the full discovery pipeline, or provide clinical recommendations.

## Run the site

The checked-in JSON files are sufficient to view the site. From the repository root:

```sh
python3 -m http.server 8000
```

Open `http://localhost:8000`. Use a local HTTP server rather than opening `index.html` directly. No server-side application, account, API key, or build step is required.

## Rebuild the data

Requirements: Python 3.10 or later, Git, and curl. The analysis dependencies are pinned in `analysis/requirements.txt`.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r analysis/requirements.txt
```

Download inputs to a temporary directory outside the repository, or set `DEMO_INPUTS` to your own empty input directory. The raw workbooks and upstream repository are not distributed with this demo.

```sh
DEMO_INPUTS=$(mktemp -d)

git clone https://github.com/Svensson-Lab/pro-hormone-predictor.git "$DEMO_INPUTS/pro-hormone-predictor"
git -C "$DEMO_INPUTS/pro-hormone-predictor" checkout --detach c478c45405adce4c62baa0719891c3c89f6edcda

curl --fail --location \
  'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-08683-y/MediaObjects/41586_2025_8683_MOESM19_ESM.xlsx' \
  --output "$DEMO_INPUTS/brp-source-edfig10.xlsx"

curl --fail --location \
  'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-47359-5/MediaObjects/41467_2024_47359_MOESM4_ESM.xlsx' \
  --output "$DEMO_INPUTS/piening2024-source.xlsx"

python analysis/brp_analysis.py --input-dir "$DEMO_INPUTS" --output brp-data.json
python analysis/cancer_analysis.py --source "$DEMO_INPUTS/piening2024-source.xlsx" --output cancer-data.json
```

These commands intentionally regenerate the two derived JSON files; they never edit the inputs. Both scripts refuse source files whose SHA-256 hashes differ from the inspected versions. Full source URLs, hashes, ranges, and exclusions are recorded in [the BRP manifest](analysis/source_manifest_brp.json) and [the cancer-study manifest](analysis/source_manifest_cancer.json).

Alternatively, `python analysis/cancer_analysis.py` downloads its source into an operating-system temporary directory, verifies it, and writes `cancer-data.json`.

### Exactly what runs from the lab's GitHub repository

The BRP script checks the pinned `FastaParser.py` and `secreted.fasta` hashes, extracts the inspected `get_cleavage_sites` method and its cleavage-pattern expression using Python's AST, and executes that bounded component through an adapter. It checks BRINP2 and GCG at minimum cleavage spacings of 3 and 5.

It does **not** run upstream `main.py`, the tissue-expression network workflow, or the upstream FASTA-writing routine. The result is a sequence-identity and motif-context check, not a fresh discovery of BRP or proof of secretion, processing, or receptor engagement. Review the linked pinned upstream source before running a reproduction: the component executes code fetched from that repository.

## Check the outputs

From the repository root:

```sh
python analysis/test_brp.py
python analysis/test_cancer.py
node scripts/test_app.cjs
```

The optional frontend checks require Node.js 18 or later and no npm packages. They exercise all 300 peptide-condition/threshold/time combinations, cancer comparison contexts, CSV exports, and unavailable-value guards. Python checks cover source counts, independent arithmetic, missingness, sequence identity, and source-reported experimental units.

These tests verify the implemented computations and selected source facts, not the biological validity of a proposed experiment.

## Data and interpretation

### BRP assays

Source: Coassolo et al., [*Prohormone cleavage prediction uncovers a non-incretin anti-obesity peptide*](https://www.nature.com/articles/s41586-025-08683-y), Nature (2025), Extended Data 10a–c.

- Cell assays measure vehicle-relative FOS expression after 1 hour at 100 µg/ml. There are three biological replicates per condition except NS1 R3A, where one technical-error entry is missing. Missing observations remain missing; R3A is not eligible for the complete-data shortlist.
- Feeding observations are **three independent experiment means, each averaging three mice per group**, not individual-mouse measurements. The recorded outcome is cumulative intake in grams, 0–3 hours after 5 mg/kg intraperitoneal dosing.
- Cell residual response is `(condition mean − vehicle mean) / (BRP mean − vehicle mean)`. Feeding residual response is `(vehicle mean − condition mean) / (vehicle mean − BRP mean)` at each time. The interface displays these ratios as percentages.
- A candidate passes only when the absolute residual is within the selected ceiling in both cell lines and the selected feeding time, with complete source observations. This is an exploratory triage rule, not an equivalence test, proof of inactivity, or a paper-specified cutoff.
- Values are not clipped to 0–100%. Time-zero feeding normalization is undefined and excluded. Reference means are shared, not paired with replicate columns; their uncertainty is not propagated. SD describes observed spread, not a confidence interval. Vehicle values normalized to 1 do not imply noiseless assays.
- Calculations use explicit observations in 10a–c, not the 10d summary sheet. The manifest records a summary-sheet discrepancy without assigning an unverified cause.

### External cancer experiment

Source: Piening et al., [*Obesity-related T cell dysfunction impairs immunosurveillance and increases cancer risk*](https://www.nature.com/articles/s41467-024-47359-5), Nature Communications 15, 2835 (2024), Fig. 3B/G/I. This is **not Svensson Lab data**.

- Body mass: within each source mouse row, compute `100 × (week 22 / week 16 − 1)`, then summarize. There are 10 NC, 10 WD, and 20 WD+Sema rows. Week 16 is the first recorded source time used here, not a verified pre-treatment baseline. This panel has no ICB subdivision.
- Tumor endpoint: day-16 volume in mm³; n = 18/7/19/14/35/35 for NC, NC+ICB, WD, WD+ICB, WD+Sema, and WD+Sema+ICB.
- Immune endpoint: percentage of **TRP-2-specific CD8+ tumor-infiltrating lymphocytes co-expressing GzmB and Perforin**, not the percentage among all CD8+ cells; n = 9/3/9/7/17/19 in the same group order.
- The Fig. 3 caption reports each n as an individual mouse. Source labels are panel-local: they do not establish cross-panel matching. Tumor and immune data pool two experiments without batch identifiers; body mass is a representative experiment. Do not join the three panels by row number or infer individual-animal correlations between them.
- Summaries are means, sample SDs, medians, and linearly interpolated quartiles. Blank endpoints are omitted without imputation, while recorded numeric zero is retained. Differences and ratios of group means are descriptive only. Body-mass contrasts use percentage-point differences, not ratios.
- No new P values, confidence intervals, adjusted treatment effects, or human-outcome predictions are computed. The smallest immune group contains only three mice. A weight response cannot substitute for an antitumor or immune response; these data do not resolve direct tumor-cell action versus host-mediated effects.

## Reuse and attribution

The research utility is the auditable comparison and control-review workflow, not a claim of new experimental findings. CSV exports retain source information and expose the observations behind each displayed comparison.

See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for source attribution, the separate source-code restrictions, and the external study's CC BY 4.0 terms. Public availability of a paper or GitHub repository is not a blanket permission to redistribute it. This repository contains original analysis/interface work and selected numerical facts or transformations, not upstream code, source workbooks, or copied paper figures.
