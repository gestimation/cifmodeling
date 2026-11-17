# Data conventions and implementation details

This vignette explains how `cifmodeling` interprets time-to-event data
and processes them in its main functions. Use it to check input
conventions, review computational steps, or investigate unexpected
results. For analysis examples, see the [worked
examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md);
for statistical formulas, see [Computational formulas in
`cifcurve()`](https://gestimation.github.io/cifmodeling/articles/v4_formulas.md).

## 1. Formulas and outcome types

The main functions use a formula-and-data interface, but the formula’s
right-hand side serves different purposes in curve estimation and
regression.

| Function | Formula | Role of other arguments |
|:---|:---|:---|
| [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md) | `Event(time, status) ~ group` | Estimates separate curves by `group`; optional case weights are supplied through `weights`. |
| [`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md) | `Event(time, status) ~ group` | Compares groups; `test` selects the test and optional case weights are supplied separately. |
| [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md) | `nuisance.model = Event(time, status) ~ age + sex` | Models the non-exposure part of the outcome distribution; `exposure = "arm"` names the exposure whose effect is estimated, while `strata` specifies censoring-estimation strata. |

The right-hand-side `group` defines **curve groups** in
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
and **comparison groups** in
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md).
Within
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md),
`strata` defines **analysis strata**, while `strata.censor` and
`strata.competing.risk` specify **censoring strata** and **working-AJ
strata**, respectively.

The response may be specified with
[`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
or a supported right-censored
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object.
Binomial models use a standard response formula. The event coding must
agree with the chosen `code.event1`, `code.event2`, and `code.censoring`
arguments. In particular,
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
does **not** put the exposure in `nuisance.model`; it is supplied
separately. The nuisance component describes part of the outcome
probability structure; its covariates may also be used for confounding
adjustment. It is not a propensity-score model.

The supported outcome types also differ by function:

| Function | Outcome types | Time-point requirement |
|:---|:---|:---|
| [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md) | `"survival"`, `"competing-risk"` | Curves are estimated over follow-up; no fixed evaluation time is required. |
| [`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md) | `"survival"`, `"competing-risk"` | Tests use follow-up information; an analysis horizon may be set using `tau`. |
| [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md) | `"survival"`, `"competing-risk"`, `"binomial"`, `"proportional-survival"`, `"proportional-competing-risk"` | `time.point` is required for fixed-time survival and competing-risk models. Models with time-constant effect parameters can use observed event times or a supplied grid. |

The following calls estimate a curve by fruit-intake group and an age-
and sex-adjusted risk ratio at eight years:

``` r

library(cifmodeling)
data(diabetes.complications)

curve <- cifcurve(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk"
)

fit <- polyreg(
  nuisance.model = Event(t, epsilon) ~ age + sex,
  exposure = "fruitq1",
  data = diabetes.complications,
  outcome.type = "competing-risk",
  effect.measure1 = "RR",
  effect.measure2 = "RR",
  time.point = 8
)
```

In
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md),
`outcome.type = NULL` determines the outcome type from the observed
event codes. For descriptive CIF estimation, specify
`outcome.type = "competing-risk"` explicitly if an event category is
absent from an analysis subset. In contrast,
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
defaults to `"competing-risk"`, and its fixed-time competing-risk model
**requires both event types to have been observed by `time.point`**.
Explicitly selecting the outcome type does not remove this requirement.
See [Direct polytomous
regression](https://gestimation.github.io/cifmodeling/articles/v2_polyreg.md)
for the model and its effect measures.

## 2. Event codes and CDISC ADaM ADTTE

Event and censoring codes differ across datasets. Check the meaning of
the status variable before fitting a model.

| Data convention            | Event | Censoring |
|:---------------------------|:------|:----------|
| Common `Surv()` convention | `1`   | `0`       |
| Binary ADaM ADTTE `CNSR`   | `0`   | `1`       |

For an ADTTE dataset, `AVAL` is typically the analysis time and
`CNSR = 0` indicates an observed event. If `CNSR` is coded 0/1, use
[`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
with explicit event and censoring codes without reversing the values:

``` r

# Illustrative ADTTE dataset: AVAL, CNSR, and ARM are its columns.
curve_adtte <- cifcurve(
  Event(AVAL, CNSR) ~ ARM,
  data = adtte,
  outcome.type = "survival",
  code.event1 = 0,
  code.censoring = 1
)

fit_adtte <- polyreg(
  nuisance.model = Event(AVAL, CNSR) ~ 1,
  exposure = "ARM",
  data = adtte,
  outcome.type = "survival",
  time.point = 12,                 # Adapt to the endpoint's time unit.
  code.event1 = 0,
  code.censoring = 1,
  code.exposure.ref = "Placebo"    # Replace with an observed ARM level.
)
```

Do not pass an unrecoded ADTTE `CNSR` variable directly to
`Surv(AVAL, CNSR)`: `Surv()` ordinarily interprets 0 as censoring. Some
ADTTE datasets use multiple positive `CNSR` codes for different
censoring reasons. Check the dataset’s metadata and, if necessary,
recode them to a single censoring status.

Competing-risk analyses require a status variable that **distinguishes
the primary event, competing event, and censoring**. A binary `CNSR`
flag alone cannot identify which type of event occurred. For example,
the usual `cifmodeling` codes are 1 (primary event), 2 (competing
event), and 0 (censoring), but other distinct codes can be specified
explicitly.

## 3. How `cifcurve()` computes estimates

The
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
workflow has four stages:

1.  **Read the data.** `util_read_surv()` interprets the formula, checks
    time and event variables, and collects group and weight information.
2.  **Estimate the curves.** The C++ engine computes Kaplan–Meier
    survival or Aalen–Johansen cumulative incidence estimates and the
    quantities needed for the selected variance method.
3.  **Construct inference.** The `error` argument selects the variance
    method, `conf.type` the confidence-interval transformation, and
    `conf.dist` the normal or Student reference distribution. For
    Student intervals, `df.method` controls the working degrees of
    freedom.
4.  **Prepare the output.** The result is a `survfit`-compatible object
    containing curves, uncertainty estimates, and available inference
    diagnostics.

For competing risks, one output convention is particularly important:
`$surv` contains **`1 - CIF` for the selected event**, not the CIF
itself. The plotting and reporting methods account for this convention.
Also, `n.risk.type` controls how risk-set counts are *displayed*; it
does not change the curve estimate or its standard error.

The default engine is `"calculateAJ_Rcpp"`. The older `"calculateKM"`
engine supports a narrower set of normal-reference variance methods. For
IF-based and Fay–Graubard-type intervals, see the [formulas
vignette](https://gestimation.github.io/cifmodeling/articles/v4_formulas.md).

**Important limitation:** weights supplied to
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
are treated as fixed. If they were estimated externally (for example,
propensity-score weights), the reported variance does not include
uncertainty from estimating the weight model.

## 4. How `polyreg()` fits a model

[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
estimates exposure effects on the selected effect-measure scale and uses
a nuisance model to maintain coherent cause-specific CIFs. Its
computation has five stages:

1.  **Prepare inputs.** The function validates event and exposure codes,
    builds the analysis dataset, sets the exposure reference, and
    usually centers and scales nuisance covariates for numerical
    stability.
2.  **Set up censoring adjustment.** For right-censored outcomes, it
    constructs inverse-probability-of-censoring weights (IPCW). The
    optional `strata` argument defines censoring-estimation strata; it
    is different from the covariates in `nuisance.model`.
3.  **Fit effect and nuisance components.** Starting values are
    calculated or provided through `data.initial.values`. Numerical
    solvers then iteratively solve the estimating equations.
4.  **Estimate uncertainty.** Fixed-time models use sandwich intervals
    by default; models with time-constant effect parameters use
    bootstrap intervals by default. The `report.*.conf` arguments
    control optional intervals.
5.  **Return results and diagnostics.** Use
    [`coef()`](https://rdrr.io/r/stats/coef.html),
    [`vcov()`](https://rdrr.io/r/stats/vcov.html),
    [`summary()`](https://rdrr.io/r/base/summary.html), and
    [`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
    to extract results. Additional diagnostic information includes
    fitted quantities, weights, and influence-function contributions
    when available.

The underlying workflow uses helpers including
`createAnalysisDataset()`, `calculateIPCW()`, `getInitialValues()` (or
`getInitialValuesProportional()`), and `solveEstimatingEquation()`.
These are internal implementation steps; users do not normally call the
helpers directly.

## 5. Convergence and reproducibility

If a model produces a warning or unexpected result, check the inputs and
model specification before changing optimization controls:

- **Verify the endpoint:** confirm the event codes, time units,
  evaluation time, and whether the intended event was observed in the
  analysis population.
- **Check the exposure:** verify that the intended reference group is
  present; set `code.exposure.ref` explicitly when necessary.
- **Inspect convergence:** set `report.optim.convergence = TRUE` in
  [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md).
  For an unstable fit, consider a simpler nuisance model, fewer
  evaluation times for models with time-constant effect parameters, or
  alternative starting values.
- **Adjust numerical controls only if needed:** `optim.parameter1` and
  `optim.parameter2` control convergence tolerances; other
  `optim.parameter*` arguments control parameter bounds and solver
  behavior. Record any changes when reporting an analysis.
- **Make bootstrap results reproducible:** use `boot.seed` and report
  the number of replications (`boot.replications`).

A convergence flag does not, by itself, establish that the model is
suitable for the data. Also check for sparse events, limited follow-up,
and extreme censoring weights.

## Where to go next

For practical analyses, see [Worked
examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md).
For the mathematical basis of the curve estimators and confidence
intervals, see [Computational formulas in
`cifcurve()`](https://gestimation.github.io/cifmodeling/articles/v4_formulas.md).
For tables and file exports, see [Reporting
results](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.md).
Complete function arguments are collected in the [function
reference](https://gestimation.github.io/cifmodeling/articles/v9_arguments.md).
