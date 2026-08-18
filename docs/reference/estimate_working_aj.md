# Estimate working Aalen-Johansen distributions by exposure and stratum

The package Aalen-Johansen engine is evaluated separately within every
observed exposure-by-stratum cell. A second call with the two event
codes exchanged supplies the competing-event cumulative incidence.

## Usage

``` r
estimate_working_aj(
  t,
  epsilon,
  exposure,
  strata = rep.int("pooled", length(t)),
  weights = rep.int(1, length(t)),
  code.event1 = 1L,
  code.event2 = 2L,
  code.censoring = 0L,
  prob.bound = 1e-07,
  prob.truncation = NULL
)
```
