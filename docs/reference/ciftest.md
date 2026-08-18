# Tests for survival and cumulative incidence curves

Provides a formula-based interface for log-rank, Fleming-Harrington,
Gray, and augmented score tests. `test = "auto"` selects log-rank for
survival outcomes and the augmented score test for competing-risk
outcomes. The `"early"` and `"late"` presets select that same
outcome-specific test with Fleming–Harrington parameters
`(rho, gamma) = (1, 0)` and `(0, 1)`, respectively. Analysis strata can
be specified independently with `strata`. For augmented tests, working
Aalen-Johansen distributions are estimated within
exposure-by-`strata.competing.risk` cells, censoring distributions are
estimated within `strata.censor`, and the closed-form augmentation is
used. Setting `iteration` to a positive integer applies that many
Scheike fixed-point refinements. Finite iterates are anchored as the
closed-form score plus the change in the AIPWCC image from the seed to
the refined direction, preserving the zeroth-step identity in the
observed sample.

## Usage

``` r
ciftest(
  formula,
  data,
  weights = NULL,
  subset.condition = NULL,
  na.action = na.omit,
  outcome.type = "auto",
  test = "auto",
  code.event1 = 1,
  code.event2 = 2,
  code.censoring = 0,
  rho = NULL,
  gamma = NULL,
  iteration = 0L,
  tolerance = NULL,
  strata = NULL,
  strata.censor = NULL,
  strata.competing.risk = NULL,
  tau = NULL,
  prob.bound = 1e-07,
  prob.truncation = NULL,
  ...
)
```

## Arguments

- formula:

  A two-sided formula with an
  [`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
  or `Surv()` response and one grouping variable on the right-hand side.

- data:

  A data frame containing variables in the formula.

- weights:

  Optional numeric case weights or the name of a weight column. The
  standard Gray branch currently accepts integer frequency weights.

- subset.condition:

  Optional character string giving a logical condition to subset `data`
  (default `NULL`).

- na.action:

  A function specifying the action to take on missing values (default
  `na.omit`).

- outcome.type:

  One of `"auto"`, `"competing-risk"`, or `"survival"`. The aliases
  `"C"` and `"S"` are accepted.

- test:

  One of `"auto"`, `"early"`, `"late"`, `"logrank"`, `"gray"`,
  `"score"`, `"augmented"`, or `"multiple"`. The aliases `"L"`, `"LR"`,
  and `"log-rank"` select log-rank; `"G"` selects Gray; and `"A"`,
  `"aug"`, and `"augmentation"` select the augmented score test. The
  early and late choices are outcome-specific single-direction presets.
  `"multiple"` (aliases `"multi"` and `"m"`) returns
  [`ciftest_mdir()`](https://gestimation.github.io/cifmodeling/reference/ciftest_mdir.md)
  using the default directions `(rho, gamma) = (2, 0), (0, 2), (0, 0)`.
  `"score"` uses the null score-IID variance; for competing-risk
  outcomes it may use censoring Kaplan–Meier strata.

- code.event1:

  Integer code of the event of interest (default `1`).

- code.event2:

  Integer code of the competing risk (default `2`).

- code.censoring:

  Integer code of censoring (default `0`).

- rho, gamma:

  Optional non-negative Fleming-Harrington weight parameters. Omitted
  values default to zero. They cannot be combined with the fixed
  `"early"` and `"late"` presets.

- iteration:

  Non-negative integer giving the number of Scheike fixed-point
  refinements. `0` returns the closed-form augmented score; positive
  values use the closed-form-anchored AIPWCC difference.

- tolerance:

  Optional positive convergence tolerance. `NULL` performs exactly the
  requested number of finite iterations.

- strata:

  Optional character vector of one or more column names defining
  analysis strata. Multiple columns define their observed interaction.
  The grouping variable cannot be included. Event risk sets, null event
  distributions, and Fleming–Harrington weight processes are constructed
  separately within these strata, and their score vectors are summed.
  This is distinct from the two nuisance-model strata.

- strata.censor:

  Optional character vector of one or more column names defining
  censoring Kaplan-Meier strata. The grouping variable may be included
  explicitly when group-specific censoring distributions are required.

- strata.competing.risk:

  Optional character vector of one or more column names defining the
  exposure-by-stratum working Aalen-Johansen models used for the
  competing-risk nuisance distribution. The grouping variable must not
  be included because it is crossed with these strata automatically.

- tau:

  Optional finite non-negative analysis horizon.

- prob.bound:

  Strictly positive numerical/positivity bound. Required nuisance
  probabilities at or below this value produce an error rather than
  being silently replaced.

- prob.truncation:

  Optional probability lower truncation strictly above `prob.bound`. It
  is applied only to positive censoring and working-survival
  denominators used by score and augmented tests.

- ...:

  Reserved for future extensions.

## Value

An object inheriting from `"ciftest"` and `"htest"`, with a
test-specific subclass. An ordinary unweighted log-rank result also
inherits from `"survdiff"` and contains a fully compatible object in
`survdiff`. Every score-based branch returns an analysis-row by contrast
matrix of individual null-score contributions in `score.iid`. Augmented
results include `score.iid.base`, `score.iid.censor`, and
`score.iid.augment` matrices and use their summed empirical
cross-product as `vcov.score`. Standard Gray results retain the
classical Gray covariance. If the optional standard-Gray score-IID
diagnostic cannot be constructed, the Gray test is still returned; its
score-IID matrices contain `NA` and the reason is recorded in
`diagnostics$score.iid.error`. Positive `iteration` values return the
requested finite-iteration score, its subject-level score matrix,
anchoring diagnostics, and fixed-point diagnostics. The event-time FH
process actually used by the fit is returned in
`diagnostics$fh.weight.process`, including whether it came from the
pooled or event-stratified AJ left limit or the native Gray
subdistribution construction.
