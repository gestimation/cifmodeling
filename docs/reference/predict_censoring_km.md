# Evaluate a censoring Kaplan-Meier distribution

Evaluate a censoring Kaplan-Meier distribution

## Usage

``` r
predict_censoring_km(
  object,
  time,
  strata = rep.int("pooled", length(time)),
  side = c("left", "right"),
  use = c("used", "raw")
)
```
