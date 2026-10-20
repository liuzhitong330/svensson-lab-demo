# Which peptide precursors deserve the next experiment?

An exploratory research demo prepared for Cathy Liu, inspired by the Svensson
Lab's peptide-discovery programme and the GLP-1/cancer duties of Stanford
Life Science Research Professional 1 position 201427.

**Live demo:** <https://liuzhitong330.github.io/svensson-lab-demo/>

The lab's published peptide library and an inspected component of its
[`pro-hormone-predictor`](https://github.com/Svensson-Lab/pro-hormone-predictor)
are connected to an external pancreatic-neoplastic-lesion RNA-sequencing study,
[GSE266899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266899).
The proposed contribution is a short, inspectable list of precursor–sequence
pairs for peptide detection and subsequent cellular follow-up. It is not a
list of validated anticancer peptides or an assertion about an undisclosed
Svensson Lab project.

## What the interaction supports

- Inspect the precursor RNA measurements in individual mice.
- Examine the predicted human peptide, the corresponding mouse sequence and
  the evidence for predicted cleavage in both species separately.
- Change explicit evidence filters and examine sensitivity to individual
  animals before choosing a candidate.
- Export the resulting shortlist with its assumptions.
- Identify the next missing measurement: an endogenous peptide identity,
  producing-cell context, processing evidence or a controlled functional test.

RNA, sequence conservation and a predicted cleavage motif are distinct
evidence layers. None establishes secretion, mature-peptide abundance,
receptor binding or cancer activity.

## Data and scope

The original count matrix contains **54,309 gene rows and 18 samples**. The
main analysis uses the **11 T2 samples only**: four low-fat-diet (LFD), four
high-fat-diet (HFD), and three HFD plus semaglutide mice. LFD animals also bear
pancreatic neoplastic lesions; they are not healthy controls. The main treatment contrast is
HFD plus semaglutide versus HFD. Earlier samples remain separate.

This is an external study by Cencioni and colleagues, not Svensson-lab data.
Semaglutide began before lesion engraftment. Bulk cellular composition,
systemic metabolism and treatment timing are alternative explanations for
RNA differences. No weight-independent effect or mechanism of established
tumour treatment is inferred.

The inspected published peptide table contains 2,672 peptide rows from
373 precursor proteins. We use the actual sequence strings and coordinates,
not the reported length field when it differs. The source table count is
reported as inspected rather than forced to equal the paper's headline total.

Of 333 precursors with supplied 5–25-aa fragments, 255 pass the sequence
mapping audit. There are 25 default exploratory candidates, or 22 when also
requiring consistent effect direction under the CPM sensitivity check.
Ntn1 illustrates a less single-animal-sensitive contrast than the larger
Coch contrast; neither establishes a biologically active fragment.

## Run locally

The checked-in derived data are sufficient. No account, key or build step is
needed:

```sh
python3 -m http.server 8769
```

Open `http://localhost:8769`. The browser loads the data over HTTP; do not
open the HTML using `file://`.

## Reproduce and audit

See [the analysis guide](analysis/peptide-cancer/README.md) for the exact
commands, dependencies, input URLs, hashes, mapping coverage, selection
criteria and checks. [QA.md](QA.md) records release checks.

The lab repository is pinned at
`c478c45405adce4c62baa0719891c3c89f6edcda`. An inspected, bounded cleavage
component is used, not the full upstream discovery pipeline. Reproduction
fetches the upstream source; it is not vendored or relicensed here.

## Sustained use

The same approach can accept a new gene-by-sample count table with biological
replicate identifiers, species, condition and treatment timing. Recheck
annotation coverage, orthology, library composition and normalization before
selection. Endogenous peptide measurements would be a new evidence layer,
not something to infer from RNA. Any functional claim requires new controlled
experiments and independent validation.

## Attribution and ownership

This is an assistant-supported exploratory project prepared for Cathy Liu.
Her documented PDAC transcriptomic project, cancer-model work and analytical
validation experience motivate the approach; this demo does not establish
independent authorship of every implementation or prior hands-on peptide
bioactivity/peptidomics experience.

See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for sources and reuse
boundaries. The original interface and analysis code are MIT-licensed; that
license does not cover third-party publications, inputs or upstream software.
No CV, transcript, reference contacts, application records or credentials are
included.

The earlier BRP-control/separate-cancer-endpoint concept is superseded.
[Its documentation](analysis/legacy-brp-controls.md) and original derived
files remain only for historical traceability; the current page does not use
them.
