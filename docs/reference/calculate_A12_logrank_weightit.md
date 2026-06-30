# A12 for log-rank-type score using directional finite differences with dw/dB from WeightIt

Returns A12 on the "mean score" scale: (1/n) \* dU_total/dB^T

## Usage

``` r
calculate_A12_logrank_weightit(
  t,
  epsilon,
  strata,
  data,
  exposure,
  weightit,
  code.exposure.ref = NULL,
  prefix = "a",
  rho = 0,
  gamma = 0,
  prob.bound = 1e-07,
  fd_rel_step = 1e-06
)
```
