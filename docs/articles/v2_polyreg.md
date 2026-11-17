# Direct polytomous regression

## Why model cumulative incidence directly?

In competing-risk studies, researchers often want to compare event
probabilities between groups. For example, is an exposure associated
with a higher eight-year risk of diabetic retinopathy, and how does it
relate to the risk of macrovascular complications?

The cumulative incidence function (CIF) is the probability that a given
event occurs first by time $`t`$. Fine–Gray regression models its
*subdistribution hazard*. By contrast, direct polytomous regression can
express exposure effects directly as risk ratios (RRs) or odds ratios
(ORs) of cumulative incidence. It also supports the subdistribution
hazard ratio (SHR) scale.

At a given time, the cause-specific CIFs must sum to **no more than
one**. Separate models for each event need not satisfy this constraint.
The approach of Tanaka and Scheike (2026), implemented in
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md),
models the events **jointly** to produce coherent probabilities.

## How does the model work?

### One set of probabilities for all events

Let $`A=0`$ denote the reference exposure group and $`A=a`$ another
exposure group. Let $`L`$ denote other covariates. For event $`j`$,
write

``` math
p_{ja}(t,l)=P(T\leq t,\,\epsilon=j\mid A=a,L=l).
```

The probability of remaining event-free at time $`t`$ is
$`p_{0a}(t,l)=1-\sum_{j=1}^{J}p_{ja}(t,l)`$. For each exposure level,
the $`J`$ event probabilities and the event-free probability therefore
sum to one.

For example, suppose the reference-group probabilities of two events are
0.20 and 0.30. Independently specifying RRs of 2 and 3 would give
probabilities of 0.40 and 0.90 in the exposed group. Their sum is 1.30,
which is impossible. Joint modeling prevents such incompatible
probabilities.

### Exposure effects and the nuisance component

[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
combines two model components:

1.  **Exposure-effect components** specify how each event probability is
    compared between exposure groups. For example, the RR for event
    $`j`$ is

    ``` math
    \mathrm{RR}_{ja}(t,l)=\frac{p_{ja}(t,l)}{p_{j0}(t,l)}.
    ```

2.  **Nuisance components** specify the remaining probability structure.
    Rather than modeling each event’s baseline risk independently, the
    method uses *polytomous log-odds products*:

    ``` math
    \phi_j(t,l)=\log\left\{
    \frac{\prod_{a=0}^{K}p_{ja}(t,l)}
         {\prod_{a=0}^{K}p_{0a}(t,l)}
    \right\},\qquad j=1,\ldots,J.
    ```

Here, exposure $`A`$ has categories $`0,\ldots,K`$. Importantly,
$`p_{0a}`$ is the probability of **remaining free of all event types**,
not just event $`j`$. With one event type and binary exposure, the
expression reduces to the usual log-odds product for a binary outcome.

Under the conditions established by Tanaka and Scheike (2026), the
exposure-effect parameters and nuisance components determine a unique
set of valid probabilities for the RR, OR, and SHR parameterizations
they study. This *variation independence* allows exposure effects to be
specified without violating the constraints on CIFs. See the original
paper for the proof.

The nuisance component describes part of the **joint outcome
distribution**; it is not simply a model for confounders. Covariates can
be included to adjust exposure comparisons, but adjustment alone does
not establish causality.

### Which effect measure should I choose?

| Measure | What it compares at a chosen time | [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md) setting |
|:---|:---|:---|
| RR | Cumulative incidence probabilities | `"RR"` |
| OR | Odds of the event occurring by that time | `"OR"` |
| SHR | Cumulative subdistribution hazards, $`-\log(1-\mathrm{CIF})`$ | `"SHR"` |

In a **fixed-time** model, SHR denotes a ratio of *cumulative*
subdistribution hazards at the chosen time, not a ratio of instantaneous
hazards. Under a proportional subdistribution hazards model, the ratio
is constant over time and corresponds to the common SHR parameter.

## Specifying a model with `polyreg()`

Specify a
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
model in three steps:

- **Outcome and covariates:** In `nuisance.model`, place
  `Event(time, status)` (or `survival::Surv(time, status)`) on the left
  and covariates other than the exposure on the right. `~ 1` specifies
  an intercept-only nuisance model.
- **Exposure and effect measures:** Specify the categorical exposure
  through `exposure`, its reference level through `code.exposure.ref`,
  and the cause-specific effect measures through `effect.measure1` and
  `effect.measure2`. For fixed-time models, set `time.point`.
- **Censoring:**
  [`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
  uses inverse probability of censoring weighting (IPCW). Specify
  `strata` to estimate the censoring distribution within
  censoring-estimation strata; otherwise, the function uses a pooled
  Kaplan–Meier estimate. This choice should reflect the censoring
  assumptions.

Choose the analysis through `outcome.type`:

| `outcome.type` | Analysis |
|:---|:---|
| `"competing-risk"` | Cause-specific CIF effects at a chosen time |
| `"survival"` | Effect on the risk of a single event at a chosen time |
| `"proportional-competing-risk"` | Time-constant effect parameters (RR, OR, or SHR) for competing events |
| `"proportional-survival"` | Time-constant effect parameters for a survival outcome |
| `"binomial"` | Effect on a binary outcome |

The default competing-risk event codes are 1 and 2, with 0 for
censoring. Change them through `code.event1`, `code.event2`, and
`code.censoring` if needed. A fixed-time competing-risk model requires
both event types to be observed by `time.point`, with sufficient
follow-up to estimate their risks.

The `"proportional-"` names refer to models with an exposure-effect
parameter that is **constant over time on the chosen effect-measure
scale**. They do not generally assume that the CIFs themselves are
proportional: a constant RR gives proportional CIFs, whereas a constant
OR or SHR applies to the respective odds or cumulative
subdistribution-hazard scale.

Set `outcome.type = "proportional-competing-risk"` (or
`"proportional-survival"` for a survival outcome) to fit such a model.
You may omit `time.point` to use the available event times. See
[`?polyreg`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
for time selection and bootstrap inference.

## Example: eight-year risks in a diabetes cohort

We use `diabetes.complications`, with follow-up time `t` in years and
`epsilon` coding diabetic retinopathy (1), macrovascular complications
(2), and censoring (0). We compare low fruit intake (the lowest
quartile, `fruitq == "Q1"`) with the other quartiles.

First, create a binary exposure indicator, using 0 as the reference
level:

``` r

analysis_data <- diabetes.complications
analysis_data$low_fruit <- as.integer(analysis_data$fruitq == "Q1")
```

### Unadjusted (crude) RRs for both events

``` r

fit_unadjusted <- polyreg(
  nuisance.model = Event(t, epsilon) ~ 1,
  exposure = "low_fruit",
  data = analysis_data,
  code.exposure.ref = 0,
  outcome.type = "competing-risk",
  time.point = 8,
  effect.measure1 = "RR",
  effect.measure2 = "RR"
)
summary(fit_unadjusted)
#> 
#>                       event1        event2      
#> ---------------------------------------------- 
#> low_fruit, 1 vs 0    
#>                       1.350         1.079       
#>                       [1.099, 1.659]  [0.674, 1.728]
#>                       (p=0.004)     (p=0.752)   
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
```

This fit estimates an eight-year RR for **each** competing event,
comparing low fruit intake with the other quartiles. An RR above 1
indicates higher cumulative incidence, and an RR below 1 indicates lower
cumulative incidence. These are association measures unless the
assumptions for causal interpretation are met.

### Adjusting for baseline covariates

To adjust for age and sex, add them to `nuisance.model` while leaving
the exposure in `exposure`. This example is displayed but not run during
vignette building to limit computation time.

``` r

fit_adjusted <- polyreg(
  nuisance.model = Event(t, epsilon) ~ age + sex,
  exposure = "low_fruit",
  data = analysis_data,
  code.exposure.ref = 0,
  outcome.type = "competing-risk",
  time.point = 8,
  effect.measure1 = "RR",
  effect.measure2 = "RR"
)
summary(fit_adjusted)
```

This model estimates exposure RRs conditional on the specified
covariates. A causal interpretation additionally requires assumptions
about confounding, positivity, and censoring.

### Using another effect scale

Change `effect.measure1` or `effect.measure2` to `"OR"` or `"SHR"` to
use another effect scale. The two events can also be modeled on
different scales:

``` r

fit_mixed <- polyreg(
  nuisance.model = Event(t, epsilon) ~ age + sex,
  exposure = "low_fruit",
  data = analysis_data,
  code.exposure.ref = 0,
  outcome.type = "competing-risk",
  time.point = 8,
  effect.measure1 = "RR",
  effect.measure2 = "SHR"
)
summary(fit_mixed)
```

Here, event 1 is summarized by an RR and event 2 by a fixed-time
cumulative subdistribution-hazard ratio. These measures are on different
scales and should not be compared numerically. The coherence theorem of
Tanaka and Scheike (2026) assumes a common effect-measure family across
causes. This mixed-scale example illustrates a software option, not a
separate result established by that theorem.

## Estimation, confidence intervals, and results

For right-censored outcomes,
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
uses IPCW-based estimation, with optional censoring strata. Tanaka and
Scheike (2026) establish an asymptotic minimum-variance property
**within a specified class of augmented IPCW estimators**, under their
assumptions; this is not a general claim of semiparametric efficiency.

For fixed-time models, confidence intervals default to a sandwich
variance estimate. For models with time-constant effect parameters, they
default to a refitting bootstrap with a normal approximation. See
[`?polyreg`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
for other bootstrap options. These intervals concern regression effects,
not the curve-based intervals in [Computational
formulas](https://gestimation.github.io/cifmodeling/articles/v4_formulas.md).

A fitted model has class `"polyreg"`, with standard accessors:

``` r

coef(fit_unadjusted)
#> [1] -1.38313159  0.30043899 -3.99147265  0.07582589
vcov(fit_unadjusted)
#>              [,1]          [,2]         [,3]          [,4]
#> [1,] 2.539857e-02  5.688786e-03 9.114685e-03  3.927105e-05
#> [2,] 5.688786e-03  1.104372e-02 2.138178e-05 -5.263603e-03
#> [3,] 9.114685e-03  2.138178e-05 7.213912e-02  2.795468e-02
#> [4,] 3.927105e-05 -5.263603e-03 2.795468e-02  5.773891e-02
nobs(fit_unadjusted)
#> [1] 978
```

`summary(fit_unadjusted)` displays event-specific estimates, confidence
intervals, and p-values. The fitted object also contains weights,
observation-level diagnostics, and summary tables. See [Reporting
results](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.md)
to format regression tables.

## Further reading

- [Overview](https://gestimation.github.io/cifmodeling/articles/v1_overview.md)
  for the full analysis workflow.
- [Worked
  examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md)
  for nonparametric estimation, regression, and plotting.
- [Reporting
  results](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.md)
  for formatted regression tables.
- [Function
  reference](https://gestimation.github.io/cifmodeling/articles/v9_arguments.md)
  for all arguments and defaults.

## References

- Tanaka S, Scheike TH. Estimating effects on cumulative incidence
  probabilities by direct polytomous regression and polytomous log-odds
  product. *Biometrics*. 2026;82(3):ujag160.
  [doi:10.1093/biomtc/ujag160](https://doi.org/10.1093/biomtc/ujag160)
