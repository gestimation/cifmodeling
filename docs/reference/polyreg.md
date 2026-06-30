# Fit coherent regression models of CIFs using polytomous log odds products

`polyreg()` fits regression models of CIFs, targeting familiar effect
measures (risk ratios, odds ratios and subdistribution hazard ratios).
Modeling the nuisance structure using polytomous log odds products
ensures that the sum of cause-specific CIFs does not exceed one, and
enables coherent modelling of the multiplicative effects.

This function follows a familiar formula–data workflow: the outcome and
covariates other than the exposure are specified through a formula in
`nuisance.model` (with
[`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
or [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) on
the LHS), and the exposure of interest is given by a separate variable
name in `exposure`. The fitted object contains tidy summaries of
exposure effects (point estimates, SEs, CIs, and p-values) and can be
summarised with
[`summary.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
or formatted with external tools such as
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html).

## Usage

``` r
polyreg(
  nuisance.model,
  exposure,
  strata = NULL,
  data,
  subset.condition = NULL,
  na.action = na.omit,
  code.event1 = 1,
  code.event2 = 2,
  code.censoring = 0,
  code.exposure.ref = 0,
  effect.measure1 = "RR",
  effect.measure2 = "RR",
  time.point = NULL,
  outcome.type = "competing-risk",
  conf.int = 0.95,
  report.nuisance.parameter = FALSE,
  report.optim.convergence = FALSE,
  report.sandwich.conf = TRUE,
  report.boot.conf = NULL,
  boot.bca = FALSE,
  boot.multiplier = "rademacher",
  boot.replications = 200,
  boot.seed = 46,
  nleqslv.method = "Newton",
  optim.parameter1 = 1e-06,
  optim.parameter2 = 1e-06,
  optim.parameter3 = 100,
  optim.parameter4 = 50,
  optim.parameter5 = 50,
  optim.parameter6 = 50,
  optim.parameter7 = 1e-10,
  optim.parameter8 = 1e-06,
  optim.parameter9 = 0.001,
  optim.parameter10 = 40,
  optim.parameter11 = 1e-06,
  optim.parameter12 = 2,
  optim.parameter13 = 0.5,
  data.initial.values = NULL,
  normalize.covariate = TRUE,
  terminate.time.point = TRUE,
  prob.bound = 1e-07
)
```

## Arguments

- nuisance.model:

  A `formula` describing the outcome and nuisance covariates, excluding
  the exposure of interest. The LHS must be `Event(time, status)` or
  `survival::Surv(time, status)`.

- exposure:

  A character string giving the name of the categorical exposure
  variable in `data`.

- strata:

  Optional character string with the name of the stratification variable
  used to adjust for dependent censoring (default `NULL`).

- data:

  A data frame containing the outcome, exposure and nuisance covariates
  referenced by `nuisance.model`.

- subset.condition:

  Optional character string giving a logical condition to subset `data`
  (default `NULL`).

- na.action:

  A function specifying the action to take on missing values (default
  `na.omit`).

- code.event1:

  Integer code of the event of interest (default `1`).

- code.event2:

  Integer code of the competing event (default `2`).

- code.censoring:

  Integer code of censoring (default `0`).

- code.exposure.ref:

  Integer code identifying the reference exposure category (default
  `0`).

- effect.measure1:

  Character string specifying the effect measure for the primary event.
  Supported values are `"RR"`, `"OR"` and `"SHR"`.

- effect.measure2:

  Character string specifying the effect measure for the competing
  event. Supported values are `"RR"`, `"OR"` and `"SHR"`.

- time.point:

  Numeric time point at which the exposure effect is evaluated for
  time-point models. Required for `"competing-risk"` and `"survival"`
  outcomes.

- outcome.type:

  Character string selecting the outcome type. Valid values are
  `"competing-risk"`, `"survival"`, `"binomial"`,
  `"proportional-survival"`, and `"proportional-competing-risk"`. The
  default is `"competing-risk"`. If explicitly set to `NULL`,
  `polyreg()` attempts to infer the outcome type from the data: if the
  event variable has more than two distinct levels, `"competing-risk"`
  is assumed; otherwise, `"survival"` is assumed. Abbreviations such as
  `"S"` or `"C"` are accepted; mixed or ambiguous inputs trigger
  automatic detection from the event coding in `data`.

- conf.int:

  Numeric two-sided level of CIs (default `0.95`).

- report.nuisance.parameter:

  Logical; if `TRUE`, the returned object includes estimates of the
  nuisance model parameters (default `FALSE`).

- report.optim.convergence:

  Logical; if `TRUE`, optimization convergence summaries are returned
  (default `FALSE`).

- report.sandwich.conf:

  Logical or `NULL`. When `TRUE`, CIs based on sandwich variance are
  computed. When `FALSE`, they are omitted (default `TRUE`). This CI is
  default for time-point models (`"outcome.type=competing-risk"`,
  `"survival"` or `"binomial"`) and is not available otherwise.

- report.boot.conf:

  Logical or `NULL`. When `TRUE`, bootstrap CIs are computed. When
  `FALSE`, they are omitted. If `NULL`, the function chooses based on
  `outcome.type` (default `NULL`). This CI is default for proportional
  models (`outcome.type="proportional-competing-risk"` or
  `"proportional-survival"`).

- boot.bca:

  Logical indicating the bootstrap CI method. Use `TRUE` for
  bias-corrected and accelerated intervals or `FALSE` for the normal
  approximation (default `FALSE`).

- boot.multiplier:

  Character string specifying the wild bootstrap weight distribution.
  One of `"rademacher"`, `"mammen"` or `"gaussian"` (default
  `"rademacher"`).

- boot.replications:

  Integer giving the number of bootstrap replications (default `200`).

- boot.seed:

  Numeric seed used for resampling of bootstrap.

- nleqslv.method:

  Character string specifying the solver used in nleqslv(). Available
  choices are `"Broyden"` and `"Newton"`.

- optim.parameter1:

  Numeric tolerance for convergence of the outer loop (default `1e-6`).

- optim.parameter2:

  Numeric tolerance for convergence of the inner loop (default `1e-6`).

- optim.parameter3:

  Numeric constraint on the absolute value of parameters (default
  `100`).

- optim.parameter4:

  Integer maximum number of outer loop iterations (default `50`).

- optim.parameter5:

  Integer maximum number of `nleqslv` iterations per outer iteration
  (default `50`).

- optim.parameter6:

  Integer maximum number of iterations for the Levenberg-Marquardt
  routine (default `50`).

- optim.parameter7:

  Numeric convergence tolerance for the Levenberg-Marquardt routine
  (default `1e-10`).

- optim.parameter8:

  Numeric tolerance for updating the Hessian in the Levenberg-Marquardt
  routine (default `1e-6`).

- optim.parameter9:

  Numeric starting value for the Levenberg-Marquardt damping parameter
  lambda (default `1e-3`).

- optim.parameter10:

  Numeric upper bound for lambda in the Levenberg-Marquardt routine
  (default `40`).

- optim.parameter11:

  Numeric lower bound for lambda in the Levenberg-Marquardt routine
  (default `1e-6`).

- optim.parameter12:

  Numeric multiplicative decrement applied to lambda when the
  Levenberg-Marquardt step is unsuccessful (default `2`).

- optim.parameter13:

  Numeric multiplicative increment applied to lambda when the
  Levenberg-Marquardt step is successful (default `0.5`).

- data.initial.values:

  Optional data frame providing starting values for the optimization
  (default `NULL`).

- normalize.covariate:

  Logical indicating whether covariates should be centered and scaled
  prior to optimization (default `TRUE`).

- terminate.time.point:

  Logical indicating whether time points that contribute estimation are
  terminated by min of max follow-up times of each exposure level
  (default `TRUE`).

- prob.bound:

  Numeric lower bound used to internally truncate probabilities away
  from 0 and 1 (default `1e-7`).

## Value

A list of class `"polyreg"` containing the fitted exposure effects and
supporting results. Key components and methods include:

- `coef`: regression coefficients on the chosen effect-measure scale

- `vcov`: variance–covariance matrix of the regression coefficients

- `diagnostic.statistics`: a data frame with inverse probability
  weights, influence function contributions, and predicted potential
  outcomes

- `summary`: event-wise tidy/glance summaries used by
  [`summary.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md)
  or
  [`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html)

- additional elements storing convergence information and internal
  tuning parameters.

Standard S3 methods are available:
[`coef.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md),
[`vcov.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md),
[`nobs.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md),
and
[`summary.polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg-methods.md).

## Details

### Overview

`polyreg()` implements **log odds product modeling** for CIFs at
user-specified time points, focusing on multiplicative effects of a
categorical exposure, or constant effects over time like Cox regression
and Fine-Gray models. It estimates multiplicative effects such as **risk
ratios**, **odds ratios**, or **subdistribution hazard ratios**, while
ensuring that the probabilities across competing events sum to one. This
is achieved through **reparameterization using polytomous log odds
products**, which fits so-called effect-measure models and nuisance
models on multiple competing events simultaneously. Additionally,
`polyreg()` supports direct binomial regression for survival outcomes
and the Richardson model for binomial outcomes, both of which use log
odds products.

### Key arguments

- `nuisance.model`: a formula with
  [`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
  or `survivai::Surv()` describing the outcome and nuisance covariates,
  excluding the exposure of interest.

- `exposure`: name of the categorical exposure variable

- `effect.measure1` and `effect.measure2`: the effect measures for
  event1 and event2 (`"RR"`, `"OR"` or `"SHR"`).

- `outcome.type`: type of the outcome variable (`"competing-risk"`,
  `"survival"`, `"binomial"`, `"proportional-survival"` or
  `"proportional-competing-risk"`).

- `time.point`: time point(s) at which the exposure effect is evaluated.
  Required for `"competing-risk"` and `"survival"` outcomes.

- `strata`: name of the stratification variable used for IPCW adjustment
  for dependent censoring.

### Outcome type and event status coding

The `outcome.type` argument must be set to:

- Effects on cumulative incidence probabilities at a specific time:
  `"competing-risk"`.

- Effects on a risk at a specific time: `"survival"`.

- Common effects on cumulative incidence probabilities over time:
  `"proportional-competing-risk"`.

- Common effects on a risk over time: `"proportional-survival"`.

- Effects on a risk of a binomial outcome: `"binomial"`.

|  |  |  |
|----|----|----|
| Setting | Codes | Meaning |
| competing-risk | `code.event1`, `code.event2`, `code.censoring` | event of interest / competing event / censoring |
| competing-risk (default) | `code.event1 = 1`, `code.event2 = 2`, `code.censoring = 0` | event of interest / competing event / censoring |
| survival | `code.event1`, `code.censoring` | event / censoring |
| survival (default) | `code.event1 = 1`, `code.censoring = 0` | event / censoring |
| survival (ADaM-ADTTE) | `code.event1 = 0`, `code.censoring = 1` | set to match ADaM convention |
| proportional-survival | `code.event1`, `code.censoring` | event / censoring |
| proportional-survival (default) | `code.event1 = 1`, `code.censoring = 0` | event / censoring |
| proportional-survival (ADaM) | `code.event1 = 0`, `code.censoring = 1` | set to match ADaM convention |
| proportional-competing-risk | `code.event1`, `code.event2`, `code.censoring` | event of interest / competing event / censoring |
| proportional-competing-risk (default) | `code.event1 = 1`, `code.event2 = 2`, `code.censoring = 0` | event of interest / competing event / censoring |

### Effect measures for categorical exposure

Choose the effect scale for event 1 and (optionally) event 2:

|                   |                   |                         |         |
|-------------------|-------------------|-------------------------|---------|
| Argument          | Applies to        | Choices                 | Default |
| `effect.measure1` | event of interest | `"RR"`, `"OR"`, `"SHR"` | `"RR"`  |
| `effect.measure2` | competing event   | `"RR"`, `"OR"`, `"SHR"` | `"RR"`  |

- `RR`: risk ratio at `time.point` or common over time.

- `OR`: odds ratio at `time.point` or common over time.

- `SHR`: subdistribution hazard ratio or common over time.

### Inference and intervals (advanced)

|  |  |  |
|----|----|----|
| Argument | Meaning | Default |
| `conf.int` | Wald-type CI level | `0.95` |
| `report.sandwich.conf` | Sandwich variance CIs | `TRUE` |
| `report.boot.conf` | Bootstrap CIs (used by `"proportional-*"` types) | `NULL` |
| `boot.bca` | Use BCa intervals (else normal approximation) | `FALSE` |
| `boot.multiplier` | Method for wild bootstrap | `"rademacher"` |
| `boot.replications` | Bootstrap replications | `200` |
| `boot.seed` | Seed for resampling | `46` |

### Optimization & solver controls (advanced)

`polyreg()` solves estimating equations with optional inner routines.

|  |  |  |
|----|----|----|
| Argument | Role | Default |
| `nleqslv.method` | Root solver | `"Newton"` |
| `optim.parameter1`, `optim.parameter2` | Outer / inner convergence tolerances | `1e-6`, `1e-6` |
| `optim.parameter3` | Parameter absolute bound | `100` |
| `optim.parameter4` | Max outer iterations | `50` |
| `optim.parameter5` | Max `nleqslv` iterations per outer | `50` |
| `optim.parameter6:13` | Levenberg–Marquardt controls (iterations, tolerances, lambda) | see defaults |

### Data handling and stability

|  |  |  |
|----|----|----|
| Argument | Meaning | Default |
| `subset.condition` | Expression (as character) to subset `data` | `NULL` |
| `na.action` | NA handling function | [`stats::na.omit`](https://rdrr.io/r/stats/na.fail.html) |
| `normalize.covariate` | Center/scale nuisance covariates | `TRUE` |
| `terminate.time.point` | Truncate support by exposure-wise follow-up maxima | `TRUE` |
| `prob.bound` | Truncate probabilities away from 0/1 (numerical guard) | `1e-5` |
| `data.initial.values` | Optional starting values data frame | `NULL` |

### Downstream use

`polyreg()` returns an object of class `"polyreg"` that contains
regression coefficients (`coef`), variance-covariance matrix (`vcov`)
and a list of event-wise *tidy* and *glance* tables (`summary`). Users
should typically access results via the S3 methods:

- [`coef()`](https://rdrr.io/r/stats/coef.html) — extract regression
  coefficients.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) — extract the
  variance–covariance matrix (sandwich or bootstrap, depending on
  `outcome.type` and the `report.*` arguments).

- [`nobs()`](https://rdrr.io/r/stats/nobs.html) — number of observations
  used in the fit.

- [`summary()`](https://rdrr.io/r/base/summary.html) — print an
  event-wise, modelsummary-like table of estimates, CIs and p-values,
  and return the underlying list of tidy/glance tables invisibly.

For backward compatibility, components named `coefficient` and `cov` may
also be present and mirror `coef` and `vcov`, respectively. The
`summary` component can be passed to external functions such as
[`modelsummary()`](https://modelsummary.com/man/modelsummary.html) for
further formatting, if desired.

### Reproducibility and conventions

- If convergence warnings appear, relax/tighten tolerances or cap the
  parameter bound (`optim.parameter1`–`3`) and inspect the output with
  `report.optim.convergence = TRUE`.

- If necessary, modify other `optim.parameter`, provide user-specified
  initial values, or reduce the number of nuisance parameters (e.g.,
  provide a small set of time points contributing to estimation when
  using `"proportional-survival"` or `"proportional-competing-risk"`).

- Set `boot.seed` for reproducible bootstrap results.

- Match CDISC ADaM conventions via `code.event1 = 0`,
  `code.censoring = 1` (and, if applicable, `code.event2` for competing
  events).

## Lifecycle

**\[experimental\]**

## See also

[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
for KM/AJ estimators;
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
for display of a CIF;
[`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
for display of multiple CIFs;
[ggsurvfit::ggsurvfit](http://www.danieldsjoberg.com/ggsurvfit/reference/ggsurvfit.md),
[patchwork::patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
and
[modelsummary::modelsummary](https://modelsummary.com/man/modelsummary.html)
for display helpers.

## Examples

``` r
data(diabetes.complications)
output <- polyreg(
  nuisance.model = Event(t, epsilon) ~ +1,
  exposure = "fruitq1",
  data = diabetes.complications,
  effect.measure1 = "RR",
  effect.measure2 = "RR",
  time.point = 8,
  outcome.type = "competing-risk"
)

coef(output)
#> [1] -1.38313159  0.30043899 -3.99147265  0.07582589
vcov(output)
#>              [,1]         [,2]         [,3]         [,4]
#> [1,]  0.007132823 -0.004524224  0.002772872 -0.002210804
#> [2,] -0.004524224  0.009639840 -0.001178880  0.004588230
#> [3,]  0.002772872 -0.001178880  0.018638168 -0.016918563
#> [4,] -0.002210804  0.004588230 -0.016918563  0.054838178
nobs(output)
#> [1] 978
summary(output)
#> 
#>                       event1        event2      
#> ---------------------------------------------- 
#> fruitq1, 1 vs 0      
#>                       1.350         1.079       
#>                       [1.114, 1.637]  [0.682, 1.707]
#>                       (p=0.002)     (p=0.746)   
#> 
#> ---------------------------------------------- 
#> 
#> effect.measure        RR at 8       RR at 8     
#> n.events              279 in N = 978  79 in N = 978
#> median.follow.up      8             -           
#> range.follow.up       [0.05, 11.00]  -           
#> n.parameters          4             -           
#> converged.by          Converged in objective function  -           
#> nleqslv.message       Function criterion near zero  -           
#> 
```
