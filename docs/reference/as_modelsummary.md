# Convert regression results to modelsummary input

Convert regression results to modelsummary input

## Usage

``` r
as_modelsummary(x, ...)

# S3 method for class 'polyreg'
as_modelsummary(
  x,
  event = c("both", "event1", "event2"),
  conf.int = x$conf.int,
  variance = NULL,
  ...
)
```

## Arguments

- x:

  A supported fitted model.

- ...:

  Model-specific inference and coefficient-selection arguments.

- event:

  For polyreg, event1, event2, or both (the default).

- conf.int:

  Interval level. For polyreg, defaults to its fitted level.

- variance:

  Covariance estimator: if or bootstrap for polyreg. Omitting it
  preserves the stored inference.

## Value

A modelsummary_list, or a named list of these for multiple outcomes.
Each contains tidy and glance data frames. Extraction does not require
modelsummary. Inference settings are recorded in glance and metadata.
