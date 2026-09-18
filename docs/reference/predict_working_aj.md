# Evaluate a working Aalen-Johansen distribution

Evaluate a working Aalen-Johansen distribution

## Usage

``` r
predict_working_aj(
  object,
  time,
  exposure,
  strata = rep.int("pooled", length(time)),
  side = c("left", "right"),
  use = c("used", "raw")
)
```
