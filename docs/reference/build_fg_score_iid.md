# Build individual Fine-Gray score-scale contributions

The base component is the known-censoring-distribution martingale score.
The censoring component is the delta-method contribution from estimating
the stratified censoring hazards. Their column sums equal the plug-in
score and zero, respectively. The matrices are diagnostic infrastructure
for the augmented branch; they do not replace Gray's classical
covariance.

## Usage

``` r
build_fg_score_iid(
  t,
  epsilon,
  x,
  code.event1 = 1L,
  code.event2 = 2L,
  code.censoring = 0L,
  strata = rep.int("pooled", length(t)),
  weights = rep.int(1, length(t)),
  rho = 0,
  gamma = 0,
  fh.weight = NULL,
  censoring = NULL,
  prob.bound = 1e-07,
  prob.truncation = NULL,
  strata.event = rep.int("pooled", length(t))
)
```
