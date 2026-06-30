# Methods for polyreg objects

S3 methods to extract coefficients, variance-covariance matrix, sample
size, formatted summaries, and tidy/glance/augment from objects returned
by
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
print(x, digits = 3, ...)

effect_label.polyreg(
  x,
  event = c("event1", "event2"),
  add.time.point = TRUE,
  add.outcome = TRUE,
  add.exposure.levels = TRUE,
  add.conf = TRUE,
  add.p = TRUE,
  value.time = NULL,
  unit.time = NULL,
  digits = 2,
  p_digits = 2,
  p_cut = 0.05,
  ...
)

# S3 method for class 'polyreg'
tidy(x, event = c("event1", "event2", "both"), ...)

# S3 method for class 'polyreg'
glance(x, event = c("event1", "event2"), ...)

# S3 method for class 'polyreg'
augment(x, ...)
```

## Arguments

- object:

  A polyreg object returned by
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

- x:

  Object to be printed or summarised. Typically a `"summary.polyreg"`
  object for `print.summary.polyreg()`, or a `"polyreg"` object for
  `tidy.polyreg()`, `glance.polyreg()`, `augment.polyreg()`, and
  `effect_label.polyreg()`.

- digits:

  Number of digits to print for parameter estimates or effect measures.
  Used by `print.summary.polyreg()` and `effect_label.polyreg()`.

- event:

  Character string indicating which event to extract. For
  `effect_label.polyreg()` and `glance.polyreg()` this is one of
  `"event1"` or `"event2"`. For `tidy.polyreg()` it can also be `"both"`
  to return rows for all events.

- add.time.point:

  Logical; if `TRUE`, `effect_label.polyreg()` appends the time point to
  the label (e.g., “at 5 years”).

- add.outcome:

  Logical; if `TRUE`, `effect_label.polyreg()` appends the outcome/event
  description (e.g., “of event 1”).

- add.exposure.levels:

  Logical; if `TRUE`, `effect_label.polyreg()` includes the exposure
  level in the label (e.g., treatment group).

- add.conf:

  Logical; if `TRUE`, `effect_label.polyreg()` includes a confidence
  interval in the label.

- add.p:

  Logical; if `TRUE`, `effect_label.polyreg()` includes a p-value or
  thresholded p-value (e.g. p \< 0.05).

- value.time:

  Optional numeric value overriding the time point stored in the
  `"polyreg"` object when constructing labels in
  `effect_label.polyreg()`.

- unit.time:

  Optional character string giving the time unit to display in labels
  constructed by `effect_label.polyreg()`, such as `"years"` or
  `"months"`.

- p_digits:

  Integer; number of digits used to format p-values in
  `effect_label.polyreg()`.

- p_cut:

  Numeric threshold used by `effect_label.polyreg()` to decide between
  printing `p < p_cut` and an exact p-value.

## Value

- `coef.polyreg()` returns a numeric vector of regression coefficients.

- `vcov.polyreg()` returns a variance-covariance matrix.

- `nobs.polyreg()` returns the number of observations.

- `summary.polyreg()` returns a list of tidy and glance summaries by
  event.

- `print.summary.polyreg()` is called for its side effect of printing a
  formatted, modelsummary-like table to the console and returns `x`
  invisibly.

- `tidy.polyreg()` returns a data frame of tidy coefficients by event.

- `glance.polyreg()` returns a data frame of model-level summaries by
  event.

- `augment.polyreg()` returns an augmented data frame containing
  diagnostics, weights, and predicted CIFs.

## See also

[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
for log odds product modeling of CIFs
