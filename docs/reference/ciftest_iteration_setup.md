# Build a finite-iteration time-directed Fine-Gray score

`iterations = 0` is intentionally handled by the closed-form caller.
This function computes one or more empirical Appendix-E refinements and
maps the resulting full-data direction to observed-data score
contributions.

## Usage

``` r
ciftest_iteration_setup(
  base,
  working,
  t,
  epsilon,
  x,
  exposure,
  weights = rep.int(1, length(t)),
  code.event1 = 1L,
  code.event2 = 2L,
  code.censoring = 0L,
  prob.bound = 1e-07,
  prob.truncation = NULL
)
```
