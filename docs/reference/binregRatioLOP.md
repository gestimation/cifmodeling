# Fit direct binomial regression for restricted mean time lost using a log-odds product parameterization

Fits an inverse probability of censoring weighted (IPCW) ratio
regression model for the conditional proportion of outcome attributable
to a specific cause in competing risks data. The target quantity is
\$\$E\\Z_1(\tau)\mid X\\ / E\\Z(\tau)\mid X\\,\$\$ where \\Z_1(\tau)\\
is the cause-specific restricted time lost up to the time horizon
\\\tau\\ for cause 1 and \\Z(\tau)\\ is the corresponding total
restricted time lost up to \\\tau\\. The model is parameterized through
a nuisance log-odds product model together with a log percentage ratio
parameter for a binary exposure. The final column of the design matrix
is assumed to be a binary exposure coded as 0/1, and the remaining
columns are treated as adjustment covariates.

## Usage

``` r
binregRatioLOP(
  formula,
  data,
  cause = 1,
  time = NULL,
  beta = NULL,
  type = c("II", "I"),
  offset = NULL,
  weights = NULL,
  cens.weights = NULL,
  cens.model = ~+1,
  se = TRUE,
  kaplan.meier = TRUE,
  cens.code = 0,
  no.opt = FALSE,
  augmentation = NULL,
  outcome = "rmtl",
  model = "log-odds",
  Ydirect = NULL,
  ...
)
```

## Arguments

- formula:

  A model formula with an
  [`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
  response on the left-hand side. The right-hand side defines the
  regression design matrix. The final column of the resulting design
  matrix must correspond to a binary exposure coded as 0/1.

- data:

  A data frame containing the variables in `formula`.

- cause:

  Integer or vector of integers specifying the event code(s) treated as
  the primary cause of interest.

- time:

  Numeric scalar giving the time horizon \\\tau\\.

- beta:

  Optional numeric vector of starting values. Its length must equal the
  number of columns in the design matrix. By default, nuisance
  coefficients are initialized at `0.1` and the exposure log percentage
  ratio parameter at `log(1.2)`.

- type:

  Character string specifying the estimating equation. `"I"` gives the
  basic IPCW estimator, and `"II"` gives the augmented estimator that
  updates the estimating equation using censoring-related augmentation
  terms.

- offset:

  Optional numeric vector of offsets. Defaults to a vector of zeros.

- weights:

  Optional observation weights for the estimating equation. Defaults to
  1 for all observations.

- cens.weights:

  Optional censoring survival probabilities used for IPCW. When
  supplied, internal fitting of the censoring model is skipped.

- cens.model:

  A right-hand side formula for the censoring model used when
  `cens.weights` is `NULL`. The default `~ + 1` gives marginal censoring
  weights.

- se:

  Logical; if `TRUE`, compute robust influence-function-based standard
  errors.

- kaplan.meier:

  Logical; if `TRUE`, allow Kaplan-Meier-type baseline estimation in the
  censoring model fit when applicable.

- cens.code:

  Integer code used for censoring in the event variable. Defaults to
  `0`.

- no.opt:

  Logical; if `TRUE`, skip numerical optimization and evaluate the
  estimating equations at the supplied `beta`.

- augmentation:

  Optional numeric vector used to augment the estimating equation.
  Defaults to a vector of zeros.

- outcome:

  Character string specifying the outcome scale. Currently only
  `"rmtl"`, restricted mean time lost decomposition up to `time`, is
  implemented.

- model:

  Character string identifying the model family. Currently only
  `"log-odds"` is implemented.

- Ydirect:

  Optional user-supplied outcome matrix replacing the internally
  constructed IPCW outcome.

- ...:

  Additional arguments passed through model-frame construction.

## Value

An object of class `"binreg"` containing at least the following
components:

- coef:

  Estimated regression coefficients. The last coefficient is the log
  percentage ratio parameter for the binary exposure.

- se.robust:

  Robust standard errors based on the estimated influence functions.

- robvar:

  Robust variance-covariance matrix.

- iid:

  Estimated influence functions used for robust variance estimation.

- p0, p1:

  Estimated conditional percentages under exposure levels 0 and 1.

- p:

  Estimated conditional percentage at the observed exposure level.

- converged:

  Logical indicating whether the numerical solver reported convergence.

## Details

Right censoring is handled through IPCW. When `cens.weights` is not
supplied, censoring weights are estimated internally from `cens.model`
using
[`mets::phreg()`](http://kkholst.github.io/mets/reference/phreg.md) and
predicted censoring survival probabilities.

Let \\A\\ denote the binary exposure in the final column of the design
matrix and let \\L\\ denote the remaining covariates. The function
models the conditional percentage for the primary cause under a log-odds
product parameterization, returning fitted percentages under `A = 0` and
`A = 1`.

The returned object includes coefficient estimates, naive and robust
variance estimates, estimated influence functions, fitted percentages,
and quantities related to the censoring model fit.

## See also

[`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md),
[`calculatePercentageLOP()`](https://gestimation.github.io/cifmodeling/reference/calculatePercentageLOP.md),
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md),
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)

## Examples

``` r
## event: 0 = censoring, 1 = primary cause, 2 = competing cause
data(diabetes.complications)

fit <- binregRatioLOP(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  time = 8,
  cause = 1,
  type = "I"
)
#> Error in binregRatioLOP(Event(t, epsilon) ~ fruitq1, data = diabetes.complications,     time = 8, cause = 1, type = "I"): could not find function "binregRatioLOP"

fit$coef
#> Error: object 'fit' not found
fit$se.robust
#> Error: object 'fit' not found

## Estimated percentages under A = 0 and A = 1 when there are no adjustment covariates
calculatePercentageLOP(
  fit$coef,
  X_L = matrix(1, nrow = 1, ncol = 1),
  offset = 0
)
#> Error in calculatePercentageLOP(fit$coef, X_L = matrix(1, nrow = 1, ncol = 1),     offset = 0): could not find function "calculatePercentageLOP"
```
