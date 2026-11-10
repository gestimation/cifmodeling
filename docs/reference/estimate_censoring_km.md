# Estimate a stratified censoring Kaplan-Meier distribution

Uses the package Kaplan-Meier counts, retaining both left and right
limits. The default keeps the existing censoring-first convention for
ciftest. With tie.order="event_first", actual failures leave before the
censoring jump, consistently with observation of failures when T \<= C.

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
  censoring.event = epsilon == code.censoring,
  tie.order = c("censoring_first", "event_first")
)
```
