# Multiple-direction quadratic score test

Combines several Fleming–Harrington score directions in the quadratic
form proposed by Ditzhaus and Friedrich. Each component is fitted with
the same outcome-specific
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
branch. Classical log-rank components use the joint,
finite-population-corrected hypergeometric covariance, with
cross-direction blocks weighted by `w_r(t) * w_s(t)`. Classical Gray
components analogously use the joint Gray (1988) covariance, including
the target-event, competing-event, and censoring contributions. Score
and augmented components use the cross-product of the stacked individual
score contributions. The default chi-square reference uses the numerical
rank of the joint covariance. This function does not perform the
permutation calibration from the original two-sample survival proposal.

## Usage

``` r
ciftest_mdir(
  formula,
  data,
  directions = c("early", "late", "unweighted"),
  weights = NULL,
  subset.condition = NULL,
  outcome.type = "auto",
  test = "auto",
  code.event1 = 1,
  code.event2 = 2,
  code.censoring = 0,
  iteration = 0L,
  tolerance = NULL,
  strata = NULL,
  strata.censor = NULL,
  strata.competing.risk = NULL,
  tau = NULL,
  prob.bound = 1e-07,
  prob.truncation = NULL,
  na.action = stats::na.omit
)
```

## Arguments

- formula:

  A two-sided formula with an
  [`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
  or `Surv()` response and one grouping variable on the right-hand side.

- data:

  A data frame containing variables in the formula.

- directions:

  Fleming–Harrington directions. A character vector may contain
  `"unweighted"`, `"early"`, and `"late"`. Alternatively supply a
  two-column numeric matrix/data frame with `rho` and `gamma`, or a
  named list of numeric length-two vectors. In this multiple-direction
  interface, the default character presets correspond to `(2, 0)`,
  `(0, 2)`, and `(0, 0)`. The stronger early/late directions keep the
  default set linearly independent.

- weights:

  Optional numeric case weights or the name of a weight column. The
  standard Gray branch currently accepts integer frequency weights.

- subset.condition:

  Optional character string giving a logical condition to subset `data`
  (default `NULL`).

- outcome.type:

  One of `"auto"`, `"competing-risk"`, or `"survival"`. The aliases
  `"C"` and `"S"` are accepted.

- test:

  Underlying component test: `"auto"`, `"logrank"`, `"gray"`, `"score"`,
  or `"augmented"`. The corresponding aliases documented for
  [`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
  are accepted. Put early and late choices in `directions`.

- code.event1:

  Integer code of the event of interest (default `1`).

- code.event2:

  Integer code of the competing risk (default `2`).

- code.censoring:

  Integer code of censoring (default `0`).

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

- na.action:

  A function specifying the action to take on missing values (default
  `na.omit`).

## Value

An object inheriting from `"ciftest_mdir"`, `"ciftest"`, and `"htest"`.
Component fits are available in `components`; the unmodified individual
contributions are in `score.iid`, and the covariance actually used by
the quadratic statistic is in `vcov.score`.

## References

Ditzhaus M, Friedrich S. More powerful logrank permutation tests for
two-sample survival data. arXiv preprint (2018).
