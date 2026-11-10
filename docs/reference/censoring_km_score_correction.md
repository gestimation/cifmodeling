# Censoring-estimation correction to left-limit IPCW regression scores

Returns q_i such that IF_beta_i = n H^-1 (s_i - mean(s) + q_i). Uses the
same source cohort as KM, with unit sampling weights. The supplied
scores already contain IPCW, including zero scores for unknown outcomes.
This product-limit derivative includes the (1-h) factor at tied jumps.
No n-by-n influence matrix is constructed: time-grouped score sums and
cumulative risk contributions give O(n log n + n*p) work and O(n*p)
storage.

## Usage

``` r
censoring_km_score_correction(
  object,
  t,
  epsilon,
  strata,
  evaluation.time,
  score
)
```
