# cifmodeling 1.2.2 reporting implementation

The ordinary reporting routes are now `gtsummary::tbl_regression(fit,
exponentiate = TRUE)` and `gtsummary::tbl_survfit(curve, times = ...,
type = "risk")`. Selecting the competing regression event adds only
`event = "event2"`. Tables use saved inference and ordinary gtsummary
formatting, merging and export methods.

## Scope and numerical baseline

The baseline was the initial working source, including the previously
implemented IF/Fay inference and shared output/input helpers. Before editing,
the working diff, source hashes, complete tests and a limited package check
were recorded under `tmp/reporting-1.2.2/baseline*`.

[PR #28](https://github.com/gestimation/cifmodeling/pull/28) included those
necessary existing source changes as well as this reporting implementation.
It was merged into `master` on 2026-10-31 as
[`6de460e`](https://github.com/gestimation/cifmodeling/commit/6de460e7be6b45741400fead55d13681483e1201).
Its final code commit was
[`5ee612f`](https://github.com/gestimation/cifmodeling/commit/5ee612f06f59d3aca9e48066bc3db2ebb2766133);
the merge commit has the same tracked-file tree. The earlier C++ IF correction is
documented separately in NEWS; it was already present in the numerical baseline
for this task. Reporting changes do not modify that recursion, estimating
equations, variance estimators, confidence-limit formulas or optimization.
Research files, simulation outputs, generated site/article files, validation
libraries, and the unrelated pre-existing README development note are excluded.
No propensity-score estimation or noninferiority algorithm is added.

## Implementation

- `R/polyreg-methods.R`: `tidy.polyreg()` accepts `exponentiate`, logical
  `conf.int` and numeric `conf.level`, retaining legacy numeric `conf.int`
  and existing positional arguments. Exponentiation changes only the estimate
  and confidence limits. Stored asymmetric/BCa limits survive unchanged when
  inference settings are unchanged; unsupported BCa reconstruction fails.
- `R/polyreg-reporting.R` and `R/polyreg.R`: one internal contrast mapping
  records the exposure, actual reference order, readable terms, coefficient
  indices and event-specific effect measures. `tbl_regression.polyreg()`
  reports exposure effects from a local copy, selects one event, and sets
  RR/OR/SHR headings and an event/horizon footnote. Four optional broom.helpers
  methods identify variables, contrasts and levels. No nuisance-model
  `terms()`, `model.frame()` or `model.matrix()` method is manufactured.
  A variable's label attribute is saved before analysis-data preprocessing
  can discard it. Real fitted-model tests cover factors and numeric category
  codes, a non-first reference category, and continuous nuisance covariates.
- `R/cifcurve-reporting.R` and `R/cifcurve.R`: `cifcurve` is a subclass of
  `survfit`. Its summary/subsetting methods retain saved inference through
  `survfit0()` and requested-time evaluation. They correct SE scale handling,
  prevent an empty unstratified fit from partially matching `strata.levels`,
  and mask interval extrapolation past positive-weight follow-up. An internal
  native grid keeps limits and time-varying diagnostics aligned; unavailable
  t limits are not replaced by normal limits.
- Existing `tidy.ciftest()`/`glance.ciftest()` supply the Gray and time-directed
  score-test examples. A short `as_gtsummary()` recipe reports their saved
  methods, temporal weights and p-values without substituting a log-rank test.
- `DESCRIPTION`/`NAMESPACE`, three new reporting test files, generated Rd,
  NEWS, and two reporting vignettes complete the package changes. The website
  retains its existing order, adds both articles, and registers the existing
  output helpers in the Reference list.

The new reporting S3 methods are `tbl_regression.polyreg`,
`model_identify_variables.polyreg`, `model_list_variables.polyreg`,
`model_list_contrasts.polyreg`, `model_list_terms_levels.polyreg`,
`summary.cifcurve` and `[.cifcurve`. Optional namespace registration does not
load reporting packages during estimation or ordinary extraction.

## Verification

### Final merged source

The final-code checks below are separate from the historical implementation
checks. In particular, the earlier count of 2,628 passing expectations is not
a test count for the final commit.

| Check | Code and environment | Result |
|---|---|---|
| [Pull-request CI](https://github.com/gestimation/cifmodeling/actions/runs/41199803111) | `5ee612f`; Ubuntu, R 4.6.2; `--no-manual --ignore-vignettes --as-cran` | Passed: 0 errors, 0 warnings, 2 notes |
| [Push CI](https://github.com/gestimation/cifmodeling/actions/runs/41199799956) | `5ee612f`; same workflow | Passed |
| Merged-source audit | `6de460e`; 40 existing R/C++ files compared with the initial working-source snapshot, normalizing CRLF/LF line endings | 5 existing R files changed; the other 35 existing sources, including all existing C++ sources, matched |
| Event-count scope check | `6de460e`; five observations with status codes `1, 2, 0, 1, 0` | C++ engine: 3 total events; legacy AJ path: 2 cause-1 events |
| Final-source vignettes | Clean Git export of `6de460e`, freshly installed in an isolated library; Windows, R 4.4.2 | All 9 vignettes rebuilt, including the substantive reporting-table and export examples |
| Final-source website | Separate clean export of `6de460e`, using the same fresh installation | Complete pkgdown site rebuilt, including all 9 articles and the Reference pages |
| Final-site audit | Generated reporting articles | 3 curve tables and 5 regression tables; article-local links resolved; saved variable label and effect headings verified; both article layouts visually inspected |

The five changed existing source files are `R/cifcurve.R`, `R/cifplot.R`,
`R/output-interop.R`, `R/polyreg-methods.R`, and `R/polyreg.R`. The initial
four-file audit preceded the later weight-forwarding correction in
`R/cifplot.R`. This audit concerns files already present in the initial working
source; it does not count newly added files. Its baseline already contained
the C++ IF correction, so it is not a claim that released v1.2.1 results are
unchanged.

The final Ubuntu CI's two notes concern the absent prebuilt vignette index
under the no-build-vignettes workflow and the existing static-analysis
findings below. They differ from the notes in the historical Windows checks.
That CI did not rebuild vignettes or the PDF manual.

The final local rebuilds used the Windows dependency versions listed below,
recorded again in `package-versions.csv`. Source provenance, installation and
build logs, session information, the event-count check, HTML audit and page
screenshots are under `tmp/reporting-1.2.2/final-report-validation/`, outside
the PR. These builds used only tracked merged source, without the unrelated
working-tree README edit or research files. The generated site was not
deployed. The PDF manual was not rebuilt in this final documentation audit.

### Initial implementation checks (historical)

These results describe the implementation checkpoints before the later
decimal-time, CI-free BCa, weight-forwarding and CI graphics-device corrections.

| Check | Result |
|---|---|
| Initial working-source tests | Passed; 7 existing intentional skips and 1 existing Rcpp build-version warning |
| Complete implementation tests | Passed with the same skips and Rcpp warning |
| Targeted reporting/inference/output tests | Passed, including requested-time risk/event counts |
| Curve numerical preservation | 150 configurations matched the initial baseline after normalizing intentional metadata changes |
| Regression numerical preservation | 6 real-data refits matched coefficients, covariance, stored intervals/p-values, potential CIFs and IFs exactly |
| Event-specific regression tables | 12 real-data tables matched extracted estimates, limits, SEs and p-values exactly |
| Frozen inference replay | 64 comparisons; maximum absolute difference 1.78e-14 |
| Initial computational-source audit | 40 existing R/C++ files checked; 4 existing reporting/wrapper files changed at this checkpoint, before the later plotting correction |
| Optional dependency loading | `library(cifmodeling)` did not load gtsummary, broom.helpers, gt, flextable or officer |
| Vignettes and exports | All 9 vignettes built; substantive chunks in both new articles and Word/HTML exports executed |
| Website | Complete pkgdown site built in a temporary destination |
| Export/page audit | Word ZIP/XML, OR/SHR headings, HTML tables and article-local links checked; representative tables and the regression article visually inspected |
| Full `devtools::check()` | PDF manual, all tests/examples and all 9 vignettes: 0 errors, 0 warnings, 3 notes |
| Initial Windows package recheck | All tests/examples and PDF manual: 0 errors, 0 warnings, 2 notes; 2,628 passing expectations, 7 existing skips and the existing Rcpp build-version warning |

The historical nine-vignette check preceded the variable-label fix.
The initial Windows package recheck used `--ignore-vignettes`; both affected articles
were additionally rendered from the package installed at that checkpoint, including their
substantive table/export examples. Their generated pages passed local-link
and table-content checks, including the saved variable label. This avoids
attributing an older installed package's output to the corrected source at
that checkpoint; the merged-source rebuild is recorded separately above.

The historical Windows full check's three notes concern installed size (7.4 MB, including 4.8 MB
of vignette documentation), inability to verify network time, and the existing
static-analysis findings below. The initial Windows check without built vignette
documentation has the latter two notes. An intermediate missing direct
`cardx` dependency declaration was corrected; the rechecks have no dependency
warnings. Validation logs, session information, numerical snapshots, output
audit and source hashes remain under `tmp/reporting-1.2.2/`, outside the PR.

Historical local validation used R 4.4.2 on Windows, gtsummary 2.6.1, broom.helpers 1.20.0,
cards 0.9.0, cardx 0.3.4, gt 0.11.1, flextable 0.10.3, officer 0.7.7,
rmarkdown 2.29, pkgdown 2.2.0 and roxygen2 7.3.3. Reporting dependencies were
installed in an isolated validation library; the user's default library and
the package's core R requirement were unchanged. Existing RStudio Pandoc 3.10
and TinyTeX were used.

## Compatibility and limits

- The declared direct regression-table route requires gtsummary >= 2.6.1
  and broom.helpers >= 1.20.0. The old installed gtsummary 2.1.0/cardx 0.2.4
  combination also failed for an ordinary survival fit and was not used as
  the acceptance environment. Estimation does not require these packages.
- `event = "both"` remains available for tidy extraction. For regression
  tables it gives an instruction to build separate event tables and merge;
  it never silently selects one event.
- For ordinary survfit helpers, use `summary()` to extract requested-time
  diagnostics. Upstream `survfit0()` does not expand arbitrary added raw
  vectors when inserting a time-zero row; the internal native grid allows
  the cifcurve summary method to return aligned diagnostics. Table conversion
  does not claim to retain all those diagnostics.
- An unstratified fit now omits the empty `strata` and `strata.levels` fields;
  the class adds `cifcurve` ahead of `survfit`. Non-default confidence levels
  are reflected in the regression table heading. No-CI curve summaries use
  missing limits, while the original fit still has no limit vectors.
- For competing-risk fits, `n.event` from the default `calculateAJ_Rcpp`
  engine counts all event causes. Its `n.risk.weighted`, `n.risk.unweighted`
  and `n.risk.ess` fields distinguish weight sums, actual row counts and ESS;
  the legacy display fields retain upward rounding.
- The alternative `engine = "calculateKM"` competing-risk path reports
  `n.event` from `n.event1`, counting only the event of interest on its own
  event-time grid. Do not describe that count as all-cause events or assume
  the C++ engine's additional risk-set fields are present. For status codes
  `1, 2, 0, 1, 0`, the respective event totals are 3 and 2. This report
  clarifies the existing engine-specific meanings; it does not change them.
- Boundary singletons and missing intervals retain the documented inference
  policy. Weights are fixed for inference; an external estimated IPTW is not
  converted into a full propensity-score IF by reporting.
- Initial check findings were not hidden: regenerating documentation fixed
  the existing missing `tie.order` argument entry. Existing static-analysis
  notes about `.constructor`, `.permitted`, `head` and `tail` remain outside
  this reporting change. Pandoc emits its deprecated `--highlight-style`
  tooling notice; devtools emits a `build_vignettes()` deprecation notice.
  Neither was suppressed.

## Review corrections (2026-10-31)

`summary.cifcurve()` now uses native selected times to align inference when
`scale != 1`, instead of multiplying scaled times back into native units.
The latter could round an event time to its left and return a pre-event SE,
interval and degrees of freedom alongside the post-event probability.
The standard survival method still supplies the scaled output and aggregate
table. Regression tests cover decimal event times, neighboring requested
times, strata, normal/Student intervals, positional scale arguments,
data-frame output, and summaries without explicitly requested times.

`tidy.polyreg()` now separates covariance/SE selection from CI reconstruction.
With `conf.int = FALSE`, changing only `conf.level` preserves the saved
estimates, SEs and p-values, including for legacy fits without coefficient
indices. An explicit covariance choice can update SEs without attempting to
reconstruct BCa limits. Unsupported requests for actual BCa reconstruction
still fail. Tests cover both events, combined extraction and exponentiation.

The targeted reporting/inference/output suite and the complete
`devtools::test()` suite passed after these fixes, with the same seven
pre-existing skips and the existing Rcpp build-version warning. Logs are
`tmp/reporting-1.2.2/review-fixes-tests.log` and
`tmp/reporting-1.2.2/review-fixes-full-tests.log`. The package-check and website
results in the historical table describe the initial implementation. They
were not repeated at the time of this follow-up; the final merged-source
verification is recorded separately above.

## Absent-weight forwarding correction (2026-10-31)

Unweighted `cifplot()` calls forwarded a local `weights = NULL` argument to
`cifcurve()`. The latter classified weights from the argument expression,
incorrectly selecting IF inference and falling back from explicit Aalen or
delta requests. Weight selection now retains bare data-column lookup while
resolving other supplied values in their original promise environment.
Plotting forwards the resolved value so an unused column named `weights`
cannot select weights accidentally. Variance-estimator formulas are unchanged.

The reported diabetes-cohort calls now select delta by default, Aalen when
requested, and delta when requested. Their estimates, SEs and confidence
limits match direct `cifcurve()` calls exactly; the three estimate-table plots
render successfully. The existing axis-range warning remains because maximum
follow-up is 11.0034 and the requested display ends at 8.

Regression tests cover forwarded/computed NULL, local vectors, bare/quoted
column names, actual supplied weights including unit weights, the unused
`weights` column, and the complete reported plotting calls. Targeted tests and
the complete test suite passed, retaining the seven existing skips and Rcpp
build-version warning. Rendering tests write to temporary PDF files. Logs are
`tmp/reporting-1.2.2/weight-fallback-tests.log`,
`tmp/reporting-1.2.2/weight-fallback-full-tests.log`, and
`tmp/reporting-1.2.2/weight-fallback-after.log`.

## CI graphics-device correction and merge (2026-10-31)

The [earlier PR run](https://github.com/gestimation/cifmodeling/actions/runs/41199304319)
tested `f96cc88` and failed one rendering test on Ubuntu/R 4.6.2: the base PDF
device could not encode the default risk-table square (`U+25A0`). This was a
graphics-device failure in the test, not evidence of a failed variance
calculation. The rendering helper now uses Cairo PDF when available, retaining
the default symbol and all numerical assertions rather than skipping the test.

The subsequent `5ee612f` pull-request and push runs both passed, as linked in
the final verification table. PR #28 was then merged as `6de460e`. The earlier
failed run remains in the Actions history and must not be mistaken for the
result of the merged code.
