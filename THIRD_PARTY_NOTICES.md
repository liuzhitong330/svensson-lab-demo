# Third-party sources and reuse boundaries

This demo distinguishes original analysis and interface code from source publications, source data, and upstream software. Attribution is not an endorsement by any cited lab, author, institution, or publisher.

## Original demo work

The comparison rules, data transformations, interface, and explanatory text are original demo work prepared for Cathy Liu. They are not the authors' published analysis software. The original website and analysis code are released under the repository's [MIT license](LICENSE); that license does not apply to third-party sources.

Source numerical facts and newly computed summaries are not assigned the upstream software's license. Conversely, a source's open-access status does not license this repository's unrelated original code. Check the applicable terms before reusing either.

## Svensson Lab BRP publication and numerical observations

**Citation:** Coassolo et al. (2025). *Prohormone cleavage prediction uncovers a non-incretin anti-obesity peptide*. Nature 641, 192–201. [DOI: 10.1038/s41586-025-08683-y](https://doi.org/10.1038/s41586-025-08683-y).

Used: selected numerical observations from the source workbook for Extended Data 10a–c, minimal peptide identity information, and independently calculated descriptive transformations.

The article is published under publisher-exclusive rights; this project does not assert that it or its workbook is CC BY. The original workbook, article text, and figures are not redistributed here. Reproduction downloads the public source from the publisher and checks its recorded hash. Access and any reuse of the original source remain subject to the source's applicable terms.

Changes made: explicit-replicate summaries; vehicle/BRP-relative residual-response calculations; adjustable cross-assay shortlist rules; and missingness and aggregation checks. These are new exploratory analyses, not author-validated results or a reproduction of every analysis in the paper.

Exact source ranges, units, SHA-256, and limitations: [source_manifest_brp.json](analysis/source_manifest_brp.json).

## Svensson Lab pro-hormone-predictor software

Repository: [Svensson-Lab/pro-hormone-predictor](https://github.com/Svensson-Lab/pro-hormone-predictor). Inspected commit: [`c478c45405adce4c62baa0719891c3c89f6edcda`](https://github.com/Svensson-Lab/pro-hormone-predictor/tree/c478c45405adce4c62baa0719891c3c89f6edcda).

No `LICENSE` file was present at this commit. Public visibility does not supply a blanket license to redistribute the upstream program. No upstream Python source or full FASTA is copied into this demo repository.

To reproduce the bounded identity check, the user fetches the original repository. The demo verifies `FastaParser.py` and `secreted.fasta`, then runs only the inspected cleavage-pattern method and expression through an adapter. This is not a full upstream pipeline rerun. This notice neither grants upstream rights nor represents permission from its authors.

## External Piening cancer study

**Citation:** Piening et al. (2024). *Obesity-related T cell dysfunction impairs immunosurveillance and increases cancer risk*. Nature Communications 15, 2835. [DOI: 10.1038/s41467-024-47359-5](https://doi.org/10.1038/s41467-024-47359-5).

The article is licensed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). Its notice includes article material under that license unless separately credited otherwise. The selected source data are attributed here; no contrary credit was identified for these selected Fig. 3 source tables.

Used: numeric source observations for Fig. 3B, 3G, and 3I. This external study is not represented as work from the Svensson Lab.

Changes made: within-row week-16-to-week-22 body-mass percentage changes, consistent endpoint summaries, descriptive group contrasts, source-cell provenance, and interactive rendering. These changes do not imply endorsement or a new causal or clinical result. No raw workbook, article text, or paper figure is republished.

Exact source URL, SHA-256, ranges, units, and exclusions: [source_manifest_cancer.json](analysis/source_manifest_cancer.json).

## Dependencies

Python reproduction uses `openpyxl` and `regex`; the frontend and its tests require no third-party JavaScript packages. Dependencies are installed separately rather than vendored. Their respective distribution licenses continue to apply. The dependency list is in [analysis/requirements.txt](analysis/requirements.txt).
