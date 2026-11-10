# Extract prediction data in the common long format

Extract prediction data in the common long format

## Usage

``` r
as_prediction_data(x, ...)

# S3 method for class 'polyreg'
as_prediction_data(
  x,
  exposure.level = x$estimand$code.exposure.ref,
  event = c("both", "event1", "event2"),
  ...
)
```

## Arguments

- x:

  Prediction data or a supported fitted model.

- ...:

  Model-specific arguments. GAM/CIF methods take newdata and inference
  settings. polyreg extracts stored training-cohort potential CIFs.

- exposure.level:

  For polyreg, one fitted exposure label; defaults to the reference.
  Values are counterfactual predictions for this level.

- event:

  For polyreg, event1, event2, or both (the default).

## Value

A data frame with .row, estimate, optional category, and optional
conf.low/conf.high columns. No newdata prediction is added to polyreg.
