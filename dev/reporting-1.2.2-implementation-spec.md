# Implementation Specification
## cifmodeling 1.2.2 — Integration with gtsummary

### 1. Objective

Improve interoperability between `cifmodeling` and `gtsummary` so that users can produce publication-ready tables directly from the results of survival and competing risks analyses.

The primary goal is to establish a simple, reproducible workflow:

**Analysis with cifmodeling → Reporting with gtsummary**

Implement the necessary compatibility improvements, automated tests, documentation updates, and two new vignettes.

Target version: **1.2.2**

Repository: https://github.com/gestimation/cifmodeling

Reference repositories:

- https://github.com/ddsjoberg/gtsummary
- https://github.com/larmarange/broom.helpers
- https://github.com/pharmaverse/cardx

### 2. General requirements

Follow these principles throughout the implementation.

1. Preserve all existing statistical estimation algorithms. Use the current working source, including the implemented IF-t and Fay inference changes, as the numerical baseline.
2. Do not modify estimating equations, variance estimators, confidence interval calculations, or optimization algorithms.
3. Preserve existing numerical results.
4. Maintain compatibility with `modelsummary`, `broom`, and existing user code wherever possible.
5. Prefer established S3 methods and standard R interfaces over package-specific workarounds.
6. Do not introduce new public table-generation functions such as `tbl_cif()` or `tbl_polyreg()`.
7. Keep `gtsummary` and other reporting packages as optional dependencies.
8. Use the existing package coding conventions and documentation style.
9. Do not commit generated HTML files or unnecessary build artifacts.
10. Document any unavoidable changes to existing interfaces.

Do not modify the estimation algorithm merely to satisfy the requirements of `broom.helpers`.

Numerical preservation applies to the stored coefficients, covariance matrices, probability estimates, influence functions, and confidence limits on their original scale. Correcting a reporting adapter that rescales a standard error incorrectly or selects the wrong event is permitted and must be tested against the stored values. Exponentiation is a reporting transformation and must not mutate the fitted object.

This release covers reporting compatibility. It does not add continuous-exposure effect estimation, propensity-score estimation, noninferiority testing, or new small-sample inference methods. Preserve the implemented `error = "fay"`, `conf.dist`, and `df.method` interfaces. Add classes or extraction methods only when a demonstrated compatibility need requires them.

Keep the core package's declared R requirement and the user's default package library unchanged. Test optional reporting dependencies in an isolated library with recorded versions. Record the initial working-tree diff and source hashes before implementation; preserve unrelated work and identify the exact scope of the final PR.

#### 2.1 User interface and output simplicity

The ordinary workflow must stay short: fit the model, call the existing gtsummary function, and optionally customize or export the resulting table. The compatibility and validation requirements belong in package internals and tests, not in additional setup steps for users.

- The standard regression route must work with `tbl_regression(fit, exponentiate = TRUE)`. Selecting the other event adds only `event = "event2"`. Do not require a custom `tidy_fun`, manual coefficient renaming, metadata injection, or a reporting configuration object for these routes.
- The standard CIF route must work with `tbl_survfit(fit, times = ..., type = "risk")`. Do not require users to remove fields, set `logse`, repair strata, or transform confidence limits themselves.
- Keep the direct `tidy()` interface backward compatible. For new documentation, teach the usual broom meanings of logical `conf.int` and numeric `conf.level`; keep the legacy numeric `conf.int` handling in the compatibility layer and a short compatibility note. Do not introduce another confidence-level alias or require users to learn both conventions for an ordinary table.
- Use existing gtsummary arguments for labels, formatting, merging, and export. Do not add cifmodeling-specific formatting arguments, table builders, export wrappers, or a second event-selection mechanism.
- Reuse `estimand` and existing fit information. Keep any necessary additional coefficient mapping in one internal reporting component or attributes, rather than duplicating it across public output fields. Add S3 methods only when the tested reporting route requires them.
- Preserve existing default tidy columns, readable labels, ordering, and modelsummary behavior except for a demonstrated defect. Supply mapping information internally rather than replacing readable labels with opaque identifiers in user-facing output.
- A default regression table should show exposure labels, the correctly named effect measure, its confidence interval, and p-values. Event and evaluation-time information belong in a concise heading or caption. Do not display nuisance parameters by default.
- A default CIF table should show groups and the requested-time probability with its confidence interval. A short footnote should identify non-default inference methods. Degrees of freedom, the Fay cap, trace calculations, engines, and other diagnostics are not default visible columns.
- A ciftest reporting example should show the test name and p-value, with the temporal weighting identified where needed. Statistic, degrees of freedom, and diagnostics may be shown using ordinary column selection when requested. Do not remove existing machine-readable fields to simplify the displayed table.
- Keep calculation failures and unavailable limits visible through the established missing-value representation and an informative note where needed. Do not add repeated diagnostic messages to successful table creation or hide an actual failure to keep the display clean.
- Do not add empty result columns, public arguments, or classes solely for future noninferiority or propensity-score support. Future extensions can use the established result conventions when their statistical methods are implemented.

Each vignette should open with a complete copyable example using the target API before introducing optional settings. Put unusual confidence levels, alternative variance selections, diagnostic extraction, and version troubleshooting in later sections. Avoid documenting internal dependency machinery as a normal user task.

### 3. Phase 1 — Improve polyreg interoperability

#### 3.1 Target API

The following workflow should work:

```r
fit <- polyreg(...)

gtsummary::tbl_regression(
  fit,
  exponentiate = TRUE
)
```

For competing risks:

```r
tbl1 <- gtsummary::tbl_regression(
  fit,
  exponentiate = TRUE,
  event = "event1"
)

tbl2 <- gtsummary::tbl_regression(
  fit,
  exponentiate = TRUE,
  event = "event2"
)

gtsummary::tbl_merge(
  list(tbl1, tbl2),
  tab_spanner = c("**Event 1**", "**Event 2**")
)
```

The default event should remain `"event1"`.

The direct `tbl_regression()` workflow reports one event at a time. `event = "both"` remains supported by `tidy.polyreg()`, but a combined regression table must not silently select one event or mix effect measures under one heading. If `tbl_regression(event = "both")` cannot represent the result correctly, fail with an instruction to create separate event tables and merge them.

If direct forwarding of `event` through `tbl_regression()` and `broom.helpers` is incompatible with the current implementation, investigate the cause and implement the smallest maintainable solution.

Do not silently substitute a different event.

#### 3.2 Improve tidy.polyreg()

Inspect:

- `R/polyreg-methods.R`
- `R/polyreg-outputs.R`
- `R/polyreg.R`

Add explicit support for:

```r
tidy.polyreg(
  x,
  event = c("event1", "event2", "both"),
  conf.int = x$conf.int,
  variance = NULL,
  exponentiate = FALSE,
  conf.level = NULL,
  ...
)
```

Required behavior:

Preserve the positions and meaning of existing positional arguments. Add the new arguments without shifting `event`, the legacy numeric `conf.int`, or `variance`.

Support the broom convention as well as the existing interface:

- A numeric `conf.int` in `(0, 1)` retains its existing meaning as the requested confidence level.
- Logical `conf.int = TRUE` requests confidence limits, using `conf.level` when supplied and the fitted level otherwise.
- Logical `conf.int = FALSE` omits the confidence-limit columns without altering the fit or its stored limits.
- `conf.level` is a numeric confidence level in `(0, 1)`. Conflicting explicitly supplied numeric levels must produce an informative error.
- With unchanged inference settings, extract the stored limits, including stored BCa limits. A supported explicit change of level or variance may use the existing reconstruction path; an unsupported change must fail instead of switching interval methods or refitting the model.
- A table's confidence-level heading must match the level of the returned limits. For a fit stored at a level other than the `tbl_regression()` default, document `conf.level = fit$conf.int` when the intention is to retain those fitted limits.

With `exponentiate = FALSE`, return the existing log-scale coefficient estimates and confidence limits.

With `exponentiate = TRUE`, exponentiate:

- `estimate`
- `conf.low`
- `conf.high`

Keep `std.error` on the original log scale.

Do not exponentiate p-values or test statistics.

Preserve the `event` column.

For `event = "both"`, preserve event identifiers and row ordering.

Ensure that confidence intervals remain consistent with the original variance estimator and confidence interval method.

Do not recompute confidence intervals unnecessarily.

#### 3.3 Investigate broom.helpers compatibility

Inspect the following functions:

- `broom.helpers::tidy_plus_plus()`
- `broom.helpers::tidy_and_attach()`
- `broom.helpers::tidy_identify_variables()`
- `broom.helpers::model_identify_variables()`
- `broom.helpers::model_get_model_frame()`
- `broom.helpers::model_get_terms()`

Determine exactly which model metadata are required for successful execution of `tbl_regression()`.

The current `polyreg` object stores:

- `coef`
- `vcov`
- `summary`
- `diagnostics`
- `outcome.type`
- `exposure`
- `estimand`
- `call`

However, the reported coefficients do not necessarily correspond to the nuisance-model coefficients.

This distinction must be preserved.

Do not implement misleading `model.frame.polyreg()`, `terms.polyreg()`, or `model.matrix.polyreg()` methods that describe the nuisance model as though it were the target effect model.

Add standard S3 methods only when their statistical meaning is correct.

If necessary, introduce minimal reporting metadata in the `polyreg` object to support coefficient identification and variable labeling.

#### 3.4 Coefficient names and labels

The existing tidy output uses descriptive terms such as:

```text
treatment, 1 vs 0
```

These are useful for reporting but may not correspond to standard model-matrix column names.

Investigate whether `broom.helpers` can identify the variables correctly.

If changes are required, distinguish:

- Machine-readable coefficient identifiers
- Original variables
- Human-readable reporting labels
- Event identifiers
- Effect measures

Preserve the existing readable labels where possible.

Keep an explicit mapping from coefficient indices to event, exposure variable, exposure contrast, reference level, readable label, and effect measure. A display-only reference row must not be treated as an estimated coefficient with a standard error or p-value.

Reporting must work from an RDS round-trip of the fitted object after the original data and formula symbols have been removed from the calling environment. Store only the reporting metadata needed to make that work. Do not re-run the original fit or evaluate its call to recover labels.

Avoid changing the default `tidy()` output unless necessary for interoperability.

#### 3.5 Effect-measure labeling

`polyreg()` supports RR, OR, and SHR.

The selected effect measure must be correctly identified for each event.

When event 1 and event 2 use different effect measures, the table headings must reflect this distinction.

Use the existing `estimand` metadata rather than inferring effect measures from coefficient values.

Event selection must apply consistently to the numeric rows, reference labels, effect-measure heading, and evaluation time or over-time interpretation. An event-2 table must not inherit an event-1 heading. Keep the original readable terms available for existing `modelsummary` users.

### 4. Phase 2 — Automated tests

Create focused tests under `tests/testthat/`.

#### 4.1 tidy.polyreg()

Verify:

- Default log-scale output remains unchanged.
- `exponentiate = TRUE` transforms estimates correctly.
- Confidence limits are transformed correctly.
- Standard errors remain on the log scale.
- p-values remain unchanged.
- Event 1 and event 2 are selected correctly.
- `event = "both"` preserves both events.
- Existing `modelsummary` functionality is unaffected.

#### 4.2 tbl_regression()

Test:

```r
gtsummary::tbl_regression(
  fit,
  exponentiate = TRUE
)
```

and:

```r
gtsummary::tbl_regression(
  fit,
  exponentiate = TRUE,
  event = "event2"
)
```

Verify numerical values against the original `polyreg()` results.

Test binary and multilevel categorical exposures represented by factors and numeric category codes, and models with continuous nuisance covariates. The current estimator treats exposure as categorical; do not describe distinct values of a continuous exposure as a per-unit effect or add a continuous-exposure estimator in this reporting change.

Verify that reference categories are correctly identified or explicitly documented when unavailable.

Verify that confidence intervals are not exponentiated twice.

Also test logical `conf.int`, legacy numeric `conf.int`, a non-default fitted level, `conf.level`, conflicting explicit levels, stored BCa limits, and unsupported interval requests. Test extraction without mutation, RDS round-trips without the original data environment, non-first reference levels, and event-specific headings when the effect measures differ.

Test the target calls without user-supplied tidiers or metadata repairs. Check that ordinary tables do not expose internal mapping fields or diagnostic columns, and that optional formatting uses existing gtsummary arguments. Verify the unchanged default tidy schema and readable terms against the baseline.

Test different combinations of RR, OR, and SHR.

Check that no nuisance parameters are incorrectly labeled as treatment effects.

#### 4.3 Existing regression tests

Run the complete existing test suite.

Compare numerical results before and after the changes.

The changes must not alter:

- Estimated coefficients
- Variance-covariance matrices
- Confidence intervals on the original scale
- p-values
- Predicted CIFs
- Influence-function calculations

### 5. Phase 3 — CIF reporting compatibility

Verify the existing workflow:

```r
fit <- cifcurve(...)

gtsummary::tbl_survfit(
  fit,
  times = c(12, 24, 36),
  type = "risk"
)
```

The `cifcurve()` object stores:

\[
S^*(t)=1-F_1(t)
\]

for competing-risk outcomes.

Confirm that `cardx::ard_survival_survfit()` correctly produces:

\[
F_1(t)=1-S^*(t)
\]

and transforms confidence intervals as:

\[
[L_F(t),U_F(t)]
=
[1-U_{S^*}(t),1-L_{S^*}(t)].
\]

Verify compatibility with:

- `survival::survfit0()`
- `summary.survfit()`
- `cardx::ard_survival_survfit()`
- `gtsummary::tbl_survfit()`

Test stratified and unstratified models.

Check estimates at specified time points against the original `cifcurve()` output.

Also verify that event counts, risk-set counts, and confidence limits are interpreted correctly.

Do not change `cifcurve()` unless an actual compatibility defect is identified.

Repair the demonstrated compatibility defects while preserving the stored probability estimates, standard errors on their declared scale, and confidence-limit calculations:

- Set `logse` consistently with `report.survfit.std.err`, so `summary.survfit()` returns the probability-scale SE once rather than multiplying it by the probability incorrectly. Cover both settings and boundary values.
- Represent an unstratified fit correctly. Removing an empty `strata` field alone is insufficient if `strata.levels` can then be matched by R's partial `$strata` lookup. Check `survfit0()` and summary extraction together.
- Retain the effective confidence level, reference distribution, transformation, variance choice, df allocation, and fixed Fay cap. The fixed `fay.bound` value is output metadata, not a new tuning argument.
- Keep time-specific `df`, `df.reason`, `critical.value`, and `conf.status` aligned with strata and times. If a downstream table cannot retain these fields, demonstrate how to report the method and obtain the diagnostics from the fit without claiming that they were retained automatically.
- Use right-continuous evaluation for the curve and its stored inference. Determine the risk set immediately before the requested time separately; it need not use the same step index as the probability.
- Do not carry a valid interval past the last positive-weight follow-up, including a requested time between that endpoint and a later zero-weight row. Test all-zero-weight strata and exact structural boundaries. Preserve the established missing-limit and reason policy; do not substitute normal intervals for unavailable t intervals.
- Keep weighted risk-set sums, actual row counts, and effective sample sizes distinct. Explain the existing rounding behavior. For competing risks, explicitly identify whether a reported event count includes all causes or the event of interest; do not label an all-cause count as cause 1.
- Cover IF-normal, IF-t with both df allocations, FG-normal, FG-t, the existing default inference paths, and the supported transformations. The optional one-sided p-value remains a separate normal-approximation result; `conf.dist = "t"` does not change that test.

The required reporting route remains `gtsummary::tbl_survfit()`. A new dedicated class or `tidy.cifcurve()` is not required merely for this specification; add one only if the validated standard route needs it. Do not replace the required route with a new public table function.

### 6. Phase 4 — ciftest reporting

Use the existing methods:

```r
generics::tidy(fit)
generics::glance(fit)
```

Create a reproducible workflow for reporting Gray and time-directed score tests.

Use existing `gtsummary` functionality wherever possible.

The reporting workflow should preserve:

- Test name
- Test statistic
- p-value
- Outcome definition
- Temporal weighting, when relevant

Distinguish tests comparing entire CIF curves from estimates at individual time points.

Do not create a new public `cifmodeling` function solely for formatting p-values.

If a custom transformation is needed, include a short, reusable example in the vignette.

Use one result row per test, including one row for a multiple-direction omnibus test. Preserve the actual method, variance source, analysis horizon, analysis and nuisance strata, and the temporal weighting used. Distinguish the original requested selector from the resolved test family. Do not interpret the test's omnibus degrees of freedom as a small-sample t degrees of freedom.

Do not use default `add_p()` in a way that silently substitutes an ordinary unweighted log-rank test for the selected Gray, score, augmented, or weighted test. Prefer a saved-result workflow with an explicit match between the displayed outcome/groups and the fitted test. If a custom test re-fits a model, document that behavior and propagate the same weights, exclusions, strata, horizon, event coding, and test specification. A built-in Gray example must respect its integer frequency-weight restriction.

Do not add a public function solely to extract or format p-values. Keep future noninferiority fields extensible, but do not implement noninferiority testing in this change.

### 7. Phase 5 — Create two new vignettes

Follow the existing vignette structure and visual conventions.

Current vignettes:

```text
v1_overview.Rmd
v2_polyreg.Rmd
v3_formulas.Rmd
v4_examples.Rmd
v5_gallery.Rmd
v6_information.Rmd
v7_arguments.Rmd
```

Add:

```text
v8_reporting_cif.Rmd
v9_reporting_polyreg.Rmd
```

Use the existing R Markdown format and package website configuration.

Do not redesign the website.

#### Vignette 8

Title:

**Reporting Cumulative Incidence with gtsummary**

Focus on `cifcurve()` and `ciftest()`.

Suggested structure:

1. Introduction
2. Preparing the data
3. Estimating cumulative incidence
4. Creating tables with `tbl_survfit()`
5. Adding confidence intervals
6. Reporting Gray and time-directed score tests
7. Formatting tables for publication
8. Exporting tables to Word and HTML

Use existing package datasets whenever possible.

Prefer the same dataset used in the package overview so that readers can follow the workflow without learning a new data structure.

Use `type = "risk"` when reporting CIF estimates from `cifcurve()`.

Explain why this transformation is required.

Keep statistical methodology explanations concise and refer to the existing vignettes for details.

#### Vignette 9

Title:

**Reporting CIF Regression with gtsummary**

Focus on `polyreg()`.

Suggested structure:

1. Introduction
2. Fitting a direct CIF regression model
3. Creating a regression table
4. Reporting RR, OR, and SHR
5. Reporting competing events separately
6. Combining event-specific tables
7. Customizing labels and confidence intervals
8. Exporting tables

The central workflow should be:

```r
fit <- polyreg(...)

gtsummary::tbl_regression(
  fit,
  exponentiate = TRUE
)
```

Demonstrate the event-specific workflow and `tbl_merge()`.

Explain how effect measures are selected and how their interpretation differs.

Do not reproduce the full mathematical explanation of the polytomous log-odds product model from the existing `v2_polyreg.Rmd`.

#### Vignette writing style

Follow these requirements:

- Use concise, publication-oriented English.
- Prefer short explanatory paragraphs.
- Place executable code immediately after its explanation.
- Avoid lengthy descriptions of basic R syntax.
- Avoid duplicating existing statistical methodology sections.
- Use actual computed results rather than manually entered numerical examples.
- Explain the statistical meaning of each table.
- Keep the code sufficiently short for readers to reuse.
- Ensure every evaluated code chunk executes successfully.
- Declare every directly used optional reporting or export package in `Suggests`, with supported versions where required. Analysis and ordinary extraction must not load a reporting package implicitly.
- In the full reporting validation environment, execute all substantive analysis, table, and export examples. Missing optional packages may cause documented conditional evaluation in a minimal installation, but an entirely skipped vignette does not satisfy the rendering acceptance criteria.
- Use temporary or build-output paths for generated Word and HTML examples. Do not write generated exports into tracked source directories.

### 8. Phase 6 — Documentation and website

Update:

- `DESCRIPTION`: version 1.2.2
- `NEWS.md`
- Relevant roxygen documentation
- Generated `.Rd` files
- Relevant examples in `polyreg()` documentation
- Package website configuration, if necessary

Preserve the existing website navigation and vignette ordering.

Ensure that both new vignettes appear under the Articles section.

Keep reporting dependencies optional.

### 9. Phase 7 — Validation

Run:

```r
devtools::document()
devtools::test()
devtools::build_vignettes()
devtools::check()
```

Also build the pkgdown website.

Verify that both new vignettes render successfully and that all tables contain correct numerical results.

Check for:

- Missing packages
- Broken links
- Deprecated functions
- Incorrect effect-measure labels
- Incorrect confidence intervals
- Incorrect event selection
- Unexpected changes in existing tests

Do not suppress warnings merely to obtain a successful build.

Before changing code, record a baseline test/check result and compare it with the final result. Separate pre-existing failures and external dependency failures from changes introduced by this PR; resolve required acceptance failures and list any remaining baseline issues explicitly. Record intentional skips rather than counting them as passed tests.

Use an isolated, mutually compatible set of `gtsummary`, `broom.helpers`, `cards`, and `cardx` versions. Reproduce the reporting routes on a supported released combination; source inspection of a development branch is not execution validation. Record package versions and session information. The installed `gtsummary 2.1.0`/`cardx 0.2.4` combination has a reproduced table-construction failure even for an ordinary `survfit`, so it must not be the only acceptance environment.

The available RStudio Pandoc can be configured with the process-local `RSTUDIO_PANDOC` path. Check the TeX toolchain required for a PDF manual separately. Document failures caused by missing build tools instead of disabling the requested checks without explanation.

Build pkgdown output in a temporary destination or exclude its generated files from the PR. Review the generated pages, numerical table bodies, both event-specific headings, and Word/HTML exports. Preserve the existing website navigation and add the two articles without redesigning it.

### 10. Deliverables

Produce one Pull Request containing:

1. The compatibility improvements
2. Automated tests
3. Two new vignettes
4. Updated documentation
5. Version 1.2.2
6. A concise implementation report

The implementation report should describe:

- Files changed
- New or modified S3 methods
- Compatibility issues discovered
- How each issue was resolved
- Tests performed and their results
- Any remaining limitations

### 11. Acceptance criteria

The implementation is complete only when:

**A.** `tbl_regression()` works directly with `polyreg()` for both events.

**B.** Estimates and confidence intervals are correctly exponentiated.

**C.** RR, OR, and SHR are labeled correctly.

**D.** Existing statistical results remain unchanged.

**E.** `cifcurve()` results can be reported through `tbl_survfit()` with numerically correct CIF estimates and confidence intervals.

**F.** Gray and time-directed score test results can be reported using the documented workflow.

**G.** Both new vignettes render successfully.

**H.** The package passes its existing test suite and package checks.

**I.** Reported standard errors have the correct scale, requested-time CIF limits respect positive-weight follow-up, and event/risk-set count labels match their definitions.

**J.** Event selection, effect-measure headings, reference categories, and confidence-level labels agree with the numeric results. Stored BCa intervals are not silently replaced.

**K.** Reporting does not re-fit the model or mutate it, and event-specific regression reporting works after saving and reloading the fit without its original data environment.

**L.** The acceptance environment has recorded compatible reporting dependencies, and optional dependencies remain optional for estimation and ordinary extraction.

**M.** The primary examples work without package-specific reporting setup, custom tidiers, manual fit-object repairs, or manual exponentiation. Default tables show the requested statistical results without internal diagnostics or nuisance coefficients; existing machine-readable extraction remains available.

### 12. Implementation strategy

Work in the following order:

1. Inspect the current repository, preserve the initial working-tree state, record numerical and check baselines, and identify the minimal changes required.
2. Implement and test `polyreg()` interoperability.
3. Verify `cifcurve()` interoperability.
4. Establish the `ciftest()` reporting workflow.
5. Write the two vignettes using the verified APIs.
6. Update documentation and package version.
7. Run all tests and builds.
8. Review the final diff for unnecessary changes.

Prefer minimal, maintainable modifications.

If a requested API cannot be implemented safely using the existing S3 structure, explain the incompatibility and propose the smallest alternative rather than introducing misleading model metadata or silently changing statistical behavior.

**The primary success criterion is a simple, statistically correct, and reproducible reporting workflow—not the number of new functions implemented.**
