# Build a finite-iteration time-directed Fine-Gray score

`iterations = 0` is intentionally handled by the closed-form caller.
This function computes one or more empirical Appendix-E refinements and
maps the resulting full-data direction to observed-data score
contributions.

## Usage

``` r
build_iterated_score_iid(
  base,
  working,
  t,
  epsilon,
  x,
  exposure,
  weights = rep.int(1, length(t)),
  iterations = 1L,
  code.event1 = 1L,
  code.event2 = 2L,
  code.censoring = 0L,
  prob.bound = 1e-07,
  tolerance = NULL
)
```
