# Compute Gray's K-sample test

Implements the predictable-weight estimating equations and covariance
from Gray (1988). The pooled weight is
`(1 - F_1(t-))^rho * F_1(t-)^gamma`; `gamma = 0` recovers the family
used by
[`cmprsk::cuminc()`](https://rdrr.io/pkg/cmprsk/man/cuminc.html).

## Usage

``` r
calculate_gray(
  t,
  epsilon,
  exposure,
  code.exposure.ref = NULL,
  prefix = "a",
  weights = rep.int(1, length(t)),
  strata = rep.int(1L, length(t)),
  data,
  rho = 0,
  gamma = 0,
  prob.bound = 1e-10
)
```

## Details

`weights` are frequency weights. `strata` defines independent analysis
strata; stratum-specific scores and covariance matrices are summed.
