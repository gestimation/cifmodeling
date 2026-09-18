# cifmodeling 1.2.0

* Added `ciftest()` for standard and weighted log-rank/Gray tests and
  closed-form augmented competing-risk score tests with separately specified
  censoring and competing-risk nuisance strata.

* Redesigned the public `ciftest()` selector around `test = "auto"`,
  `"logrank"`, `"gray"`, `"score"`, or `"augmented"`; `"early"` and `"late"`
  select outcome-specific Fleming--Harrington presets. Survival outcomes now
  default to log-rank and competing-risk outcomes to the augmented score.
  Added `"L"`, `"LR"`, and `"log-rank"` aliases for log-rank, `"G"` for Gray,
  and `"A"`, `"aug"`, and `"augmentation"` for the augmented score test.
  Added a null robust score option, individual score-contribution matrices,
  two-sided multigroup omnibus tests, and `survdiff` compatibility for the
  ordinary unweighted log-rank result.

* Renamed the public analysis-stratification argument from `strata.event` to
  `strata`. `strata.censor` is formally available for the competing-risk null
  score, while `strata.competing.risk` remains specific to augmentation.

* Added `prob.truncation` for derivative-aware probability lower truncation
  of censoring KM and working AJ survival denominators. Raw probabilities
  continue to determine
  structural positivity failures; raw and used minima and truncation counts
  are returned in diagnostics.

* Added `ciftest_mdir()` to combine user-selected unweighted, early, late, or
  custom Fleming--Harrington directions in a Ditzhaus--Friedrich quadratic
  score test with rank-adaptive asymptotic chi-square calibration. Classical
  log-rank directions use a joint finite-population-corrected hypergeometric
  covariance, while classical Gray directions use the bilinear extension of
  Gray's event, competing-event, and censoring covariance. Consequently,
  linearly dependent FH directions do not change the quadratic test.
  `ciftest(test = "multiple")` and its `"multi"` and `"m"` aliases use the
  linearly independent default directions `(2, 0)`, `(0, 2)`, and `(0, 0)`;
  the single-direction early and late presets remain `(1, 0)` and `(0, 1)`.

* Added analysis stratification to `ciftest()` independently of
  `strata.censor` and `strata.competing.risk`. Log-rank and
  Gray scores/covariances, Fine--Gray risk sets, and FH weight processes are
  constructed within event strata and then summed. Closed-form augmentation
  uses event-by-exposure-by-competing-risk-by-censoring computational cells;
  positive fixed-point iterations remain restricted to pooled models. All
  three strata arguments accept one or more column names; multiple columns
  define their observed interaction and share one complete-case filter.

* Added finite Scheike fixed-point refinement to `ciftest()`. The single
  `iteration` argument selects the closed-form augmented score at zero and a
  fixed number of iterated time-weighted Fine--Gray updates at positive
  integers, with subject-level scores and residual diagnostics returned.
  Appendix-E iteration now evaluates the working competing-risk ratio with
  the coherent Aalen--Johansen survival and enforces positivity only where a
  future cause-1 score contribution is active; censoring times beyond the last
  cause-1 event are excluded from the correction map and reported in
  `diagnostics$iteration.support`.

* Anchored positive finite iterations and the internal direct fixed-point
  reference at the closed-form augmented score using
  `closed-form + AIPWCC(value) - AIPWCC(seed)`. This preserves the zeroth-step
  identity in each observed sample while retaining the raw `AIPWCC(seed)` as
  a batch-only diagnostic. Anchoring discrepancies are returned under
  `diagnostics$iteration.anchor`.

* Added an internal batched `ciftest` engine for simulation studies. It reuses
  censoring Kaplan--Meier and working Aalen--Johansen nuisance fits within one
  data set, preserves scalar `ciftest()` result objects, and moves the nested
  Fine--Gray score-IID and closed-form augmentation loops to Rcpp kernels.
  Fixed Fleming--Harrington weights with a common nuisance specification are
  evaluated together in a multiweight traversal.

* Replaced the default Fine--Gray score-IID traversal with an internal
  sort-plus-prefix/reverse-suffix C++ kernel. Risk moments, the base IID,
  censoring derivative, and censoring-martingale IID are evaluated without the
  former event-by-censoring-time-by-subject loop. The independent legacy
  kernel remains available through
  `options(cifmodeling.ciftest.fg.engine = "legacy")`; `"check"` evaluates
  both kernels and stops if their numerical outputs differ beyond tolerance.

* Added a cell-aggregated prefix C++ kernel for closed-form augmentation.
  Reverse at-risk sums replace the censoring-time-by-subject centering loop,
  and cell-specific censoring-hazard prefixes replace the subject-level
  martingale traversal. The legacy implementation remains available through
  `options(cifmodeling.ciftest.augmentation.engine = "legacy")`; `"check"`
  evaluates both implementations and verifies all returned IID and diagnostic
  quantities.

* Added a batch-only direct Appendix-E fixed-point reference solver. It builds
  the affine iteration operator, reuses one QR factorization across FH weights,
  and reports residual, spectral-radius, and conditioning diagnostics. Finite
  iteration and the public scalar `ciftest()` interface remain unchanged.
  A direct solution is now marked converged only when its algebraic residual
  meets tolerance and the empirical operator has spectral radius below one;
  solving a non-contractive linear system is retained as a diagnostic result.

* Added a batch-only Fine--Gray score-IID reference for simulation comparisons.
  It uses the Fine--Gray influence-function cross-product as its test variance
  and is kept separate from the public standard Gray test.

* Added `diagnostics$fh.weight.process` to `ciftest()` results so the exact
  event-time FH process and its pooled-AJ or native-Gray construction can be
  audited and compared with independent implementations.

# cifmodeling 1.0.0

* Third CRAN release.

* Added `cifflowchart()` for lightweight flowchart summaries of exclusions, groups, and outcome status.
* Added one-sided normal approximation tests in `cifcurve()` via `time.point` and `null.hypothesis`.
* Improved handling of risk-set display and event coding options.
* Updated documentation, examples, vignettes, and pkgdown site.
* Fixed bugs and expanded test coverage.

# cifmodeling 0.9.8

* Second CRAN release.

# cifmodeling 0.9.7

* Added `n.risk.type` to `cifcurve()` to optionally display weighted, unweighted,
  or Kish ESS risk set sizes without changing estimates or SEs.
* Modified the risk table of `cifplot()` to display integer values.
* Simplified the special-case handling for `cloglog` by unifying axis-application
  logic, improving maintainability without changing user-facing results.
* Added/strengthened robust tests for axis limits and breaks behavior.

# cifmodeling 0.9.6

* Documentation and tests were polished for CRAN submission.
* Fixed GCC r-devel warning (`-Wstringop-overread`) triggered in `calculateAJ_Rcpp` (no change to user-facing results).

# cifmodeling 0.9.5

* First CRAN release.
* Documentation and tests were polished for CRAN submission.

# cifmodeling 0.9.4

* Sixth CRAN submission.

Invalid DESCRIPTION file, software names were single quoted. Rebuilt and resubmitted.

# cifmodeling 0.9.3

* Fifth CRAN submission.

Addressed CRAN feedback regarding file-system and environment management.
Removed development-time scripts under tools/ that performed bulk renaming
and wrote to the package directory or the current working directory.
These files are no longer included in the source tarball.
Eliminated the remaining use of .GlobalEnv in internal helpers.
panel_as_formula_global() now evaluates formulas in the calling environment
(parent.frame()) instead of .GlobalEnv, fully complying with the
“no modification of the .GlobalEnv” CRAN policy.
No user-visible changes; only internal cleanup for CRAN compliance.

# cifmodeling 0.9.2

* Fourth CRAN submission.
* Invalid file URLs  
  (articles/gallery.html, articles/polyreg.html, articles/formulas.html)
  were replaced by URLs in the web and verified using urlchecker::url_check(). 
  Rebuilt and resubmitted.

# cifmodeling 0.9.1

* Third CRAN submission.
* The DOI reported by CRAN (`10.21105/joss.00510`) does not appear to be present
  in the package source. Verified using recursive search and a freshly created
  tarball. Rebuilt and resubmitted.

# cifmodeling 0.9.0

* Second CRAN submission.

# cifmodeling 0.8.2

* Bug fixes
* README.Rmd, Vignettes and site were updated. 
* Cumulative hazard and log-log plots were implemented

# cifmodeling 0.8.1

* README.Rmd, Vignettes and site were updated. 
* polyreg S3 methods were implemented

# cifmodeling 0.8.0

* Breaking change: `cifplot()` and `cifpanel()` now return structured objects
  with explicit `plot`/`list.plot` and `patchwork` elements, metadata, and
  printing handled via new S3 methods. Automatic printing is limited to
  interactive sessions, and panel outputs store engine metadata in
  `print.info$engine`.

* Standardized outcome.type canonical values to lowercase; legacy uppercase and short aliases remain accepted. Internals now use tolower().

# cifmodeling 0.7.0

* Bug fixes
* README.Rmd, Vignettes and site were updated. 

# cifmodeling 0.6.0

* Initial CRAN submission.

# cifmodeling 0.5.0

* Bug fixes

# cifmodeling 0.4.0

* Bug fixes

# cifmodeling 0.3.0

* `cifpanel()` and panel mode of * `cifplot()` were developed.

# cifmodeling 0.2.0

* Bug fixes

# cifmodeling 0.1.0

* `cifcurve()`, `cifplot()` and `polyreg()` were developed.
