# Estimate a stratified censoring Kaplan-Meier distribution

Uses the package Kaplan-Meier engine with censoring as the event. Both
left and right limits are retained because Fine-Gray weights must be
predictable at tied event/censoring times.

## Usage

``` r
estimate_censoring_km(
  t,
  epsilon,
  code.censoring = 0L,
  strata = rep.int("pooled", length(t)),
  weights = rep.int(1, length(t)),
  prob.bound = 1e-07,
  prob.truncation = NULL,
  censoring.event = epsilon == code.censoring
)
```
