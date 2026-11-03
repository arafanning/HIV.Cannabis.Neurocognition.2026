# Neurocognition syntax recovery: version reconciliation

**Repository:** `Nonlinear.Associations.2026` (current GitHub name pending user-initiated rename to `HIV.Cannabis.Neurocognition.2026`)

**Review status:** A proposed reconciliation, **not** a verified final analysis. Nothing in this PR was executed using participant-level data. All original files remain untouched on `main`.

## What was recovered

The existing repository has six R files (`01_setup.R`, `02_data_preparation.R`, `03_disease_severity_factor.R`, `04_descriptive_statistics.R`, `05_primary_models.R`, and `08_tables_figures.R`). The second file is a large working record containing candidacy analyses, figures, exploratory models, repeated transformations, and later statistical checks. `05_primary_models.R` is explicitly labeled as primary **candidacy** models and uses a simpler severity-by-cannabis specification. Do **not** treat those files as one sequential, finalized pipeline.

An uploaded **October 4, 2026** consolidated document, *Through the Haze: Cannabis Use Level and Neurocognition Among People With HIV*, contains approximately 11,588 text lines. It explicitly calls **HIV-negative/cannabis-negative control-referenced scoring primary** and within-PWH z-scoring a sensitivity analysis. It also contains the original within-PWH test-score transformations and domain definitions. This was an intermediate specification, not the last analytic decision.

Subsequent analysis discussions in **October 2026** reverted to the original **within-PWH (N=113) standardization** as primary. They preserved the hierarchical specification: Step 1 HIV severity; Step 2 cannabis; Step 3 demographic covariates; Step 4 cannabis-by-severity interactions. Inference uses **HC3** robust covariance estimates and planned contrasts. Race/ethnicity removal and observed-use spacing were discussed as checks.

## New file

[`09_recovered_2026_primary_analysis.R`](09_recovered_2026_primary_analysis.R) is a **new, self-contained reconstruction**, not an unmodified historic script. It consolidates:

- Ten neuropsychological tests, reverse-oriented where needed, standardized **within the PWH analytic sample**, with the original five domain composites and 9-of-10 test minimum for global.
- Three-indicator latent HIV severity from negative square-root current and nadir CD4 and log1p viral load, retaining the original 0.001 fixed residual in the CFA for the documented Heywood problem.
- Past-year cannabis groups None, Low and High, categorized at the **median among non-zero users**. Previously recorded groups: 62, 27 and 24.
- Disease severity, cannabis, demographics, and interaction as **four hierarchical steps on identical complete cases for each outcome**.
- HC3 coefficient inference and joint 2-df Wald F tests for cannabis groups; HC3 interaction tests; robust estimated marginal means; Tukey pairwise comparisons; planned equal-spaced and log-median-spaced linear/quadratic group contrasts.
- A model **without race/ethnicity** as a robustness check, without selecting the preferred specification by p value.

**Previously reported use-group medians** were 0, 36 and 720 grams; group scores were 0, 3.6109 and 6.5806 on a log1p scale. These numbers are comments / reference checks, not hardcoded as individual participant values.

## Items requiring validation before the script could be authoritative

1. The original source dataset `data/HIV+113.sav` is **restricted and unavailable here**. This script has not been run. Do not present p values or model results as recovered from the new code.
2. Confirm the current population, data dictionary, and coding of `phq_7_degree` as numeric education. Confirm the dichotomous `race_eth_binary` mapping and overlap across cannabis groups; only one White non-Hispanic participant was previously reported in the High group.
3. Refit/inspect the HIV severity CFA and its fixed residual, verify factor-score orientation and participant alignment, and document the Heywood constraint. A one-degree-of-freedom CFA does not by itself establish measurement validity.
4. The archived **October 4 control-referenced** analysis is not reproduced by this script; it remains an **earlier specification**. Neither scoring approach represents published population norms.
5. Validate all six cognitive composites, missing-data rules, treatment of participants whose domains rely on one test, and the common-sample strategy within each outcome.
6. Decide how to report planned contrasts, Tukey comparisons, and multiplicity across six outcomes before manuscript submission. The working script reports HC3 contrasts without an automatic across-domain correction.
7. Recover the user-only continuous amount model and detailed influence, alcohol, and alternative HIV-index checks separately. These are **not** silently included in the primary file.
8. Inspect the previously observed discrepancies between the lower-use range (0.05 versus 0.5 grams) in manuscript versions; do not alter category membership or text labels without raw data verification.
9. Verify all code execution and manuscript-consistent numerical estimates before merging or reporting.

## What this PR does not do

It does not change any existing analysis file, manuscript, results, PDF, dataset, or collaborator's work. It does not rename the repository or its URLs. It does not imply that the current script reconstructs the entire final analysis. It provides an inspectable starting point for recovering and validating the current primary analysis **one repository at a time**.

**Recommended next action:** compare the newly assembled file with the user's local final script and run it only in an authorized environment containing the original data. After confirming model outputs, migrate the validated code into the primary documented analysis workflow. Then address the ADL repository.
