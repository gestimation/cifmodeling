# Internal Kaplan-Meier interface

The numerical engine always retains fractional weighted counts. The
public package interface defaults to integer reporting for plotting
compatibility; internal estimation code can request the unrounded
weighted totals.

## Usage

``` r
calculateKM(
  t,
  d,
  w = numeric(),
  strata = integer(),
  error = "greenwood",
  count.type = c("integer", "numeric")
)
```

## Arguments

- t:

  Follow-up times.

- d:

  Event indicator.

- w:

  Optional non-negative case weights.

- strata:

  Optional integer stratum codes.

- error:

  Variance calculation method.

- count.type:

  Return weighted `n.risk`, `n.event`, and `n.censor` as upward-rounded
  integers (`"integer"`) or exact numeric totals (`"numeric"`).
