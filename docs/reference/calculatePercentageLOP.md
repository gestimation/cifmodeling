# Calculate fitted percentages from log-odds product regression coefficients

Converts regression coefficients from
[`binregRatioLOP()`](https://gestimation.github.io/cifmodeling/reference/binregRatioLOP.md)
into fitted percentages under exposure levels 0 and 1.

## Usage

``` r
calculatePercentageLOP(beta, X_L, offset, tol = 1e-08, eps = 1e-10)
```

## Arguments

- beta:

  Numeric vector of regression coefficients. Its length must equal
  `ncol(X_L) + 1`, where the final element is the exposure log
  percentage ratio parameter.

- X_L:

  Numeric matrix of adjustment covariates excluding the binary exposure.
  A vector is coerced to a one-column matrix.

- offset:

  Numeric vector of offsets with one value per row of `X_L`.

- tol:

  Numeric tolerance used to switch to a numerically stable expression
  when the nuisance linear predictor is close to zero.

- eps:

  Small positive constant used to truncate fitted percentages away from
  0 and 1.

## Value

A numeric matrix with two columns:

- p_0:

  Estimated percentage under exposure level 0.

- p_1:

  Estimated percentage under exposure level 1.

## Details

The coefficient vector is assumed to consist of nuisance coefficients
for the covariates `X_L`, followed by a final coefficient representing
the log percentage ratio parameter for the binary exposure.

For each observation, the function returns two fitted percentages: one
for exposure level 0 and one for exposure level 1. These are obtained
from the nuisance linear predictor and the log percentage ratio
parameter under the log-odds product parameterization.

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

## Estimated percentages of restricted mean time lost under A = 0 and A = 1 when there are no adjustment covariates
calculatePercentageLOP(
  fit$coef,
  X_L = matrix(1, nrow = 1, ncol = 1),
  offset = 0
)
#> Error in calculatePercentageLOP(fit$coef, X_L = matrix(1, nrow = 1, ncol = 1),     offset = 0): could not find function "calculatePercentageLOP"
```
