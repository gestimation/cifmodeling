# Methods for `polyreg` model objects

S3 methods to extract coefficients, variance–covariance matrices, sample
size, and formatted summaries from objects returned by
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md).

## Usage

``` r
# S3 method for class 'polyreg'
coef(object, ...)

# S3 method for class 'polyreg'
vcov(object, type = c("default", "sandwich", "bootstrap"), ...)

# S3 method for class 'polyreg'
nobs(object, ...)

# S3 method for class 'polyreg'
summary(object, ...)

# S3 method for class 'summary.polyreg'
print(summary, digits = 3, ...)
```

## Arguments

- object:

  A `"polyreg"` object returned by
  [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md).

- ...:

  Further arguments passed to or from methods.

- type:

  Character string; one of `"default"`, `"sandwich"`, or `"bootstrap"`.
  When `"default"`, the function chooses between sandwich and bootstrap
  variance based on the original
  [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
  settings, using `outcome.type`, `report.sandwich.conf`, and
  `report.boot.conf`. (Used only by `vcov.polyreg()`.)

- digits:

  Number of digits to print for parameter estimates. (Used only by
  `print.summary.polyreg()`.)

- x:

  A `"summary.polyreg"` object, as returned by `summary.polyreg()`.

## Value

- `coef.polyreg()` returns a numeric vector of regression coefficients.

- `vcov.polyreg()` returns a variance–covariance matrix.

- `nobs.polyreg()` returns the number of observations.

- `summary.polyreg()` returns a list of tidy and glance summaries by
  event.

- `print.summary.polyreg()` is called for its side effect of printing a
  formatted, `modelsummary`-like table to the console and returns `x`
  invisibly.

## See also

[`polyreg`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
