# Build the closed-form competing-risk augmentation

This computes the working-model process \$\$H_a(t,X)=\int_t^\tau
a(u)(X-e_X(u))G_C(u)d\Lambda_1(u\|X)\$\$ in each
exposure-by-competing-risk-by-censoring nuisance cell, then evaluates
the centered censoring-martingale projection. The returned matrix is on
the score scale; its column sums are the closed-form augmentation score.

## Usage

``` r
build_closed_form_augmentation(
  base,
  working,
  t,
  epsilon,
  x,
  exposure,
  strata.censor = rep.int("pooled", length(t)),
  strata.competing.risk = rep.int("pooled", length(t)),
  weights = rep.int(1, length(t)),
  code.event1 = 1L,
  code.event2 = 2L,
  code.censoring = 0L,
  prob.bound = 1e-07,
  strata.event = rep.int("pooled", length(t))
)
```
