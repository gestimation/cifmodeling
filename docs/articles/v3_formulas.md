# Computational formulas in cifcurve()

## Formulas for the Kaplan–Meier estimator

Let $`t_1 < t_2 < \cdots < t_D`$ represent the distinct event times. For
each event time $`j=1,\ldots,D`$, let $`Y_j`$ be the size of the risk
set (the number of surviving observations) just prior to $`t_j`$.

Survival data is assumed when `outcome.type = "survival"`. Let $`d_j`$
be the number of failures at $`t_j`$. The Kaplan–Meier estimator of the
survival function $`S(t)`$ is as below.

**Kaplan–Meier estimator**

``` math
\hat{S}(t) = \prod_{t_j \le t} \left(1 - \frac{d_j}{Y_j} \right)
```
Note that the estimator is defined to be right-continuous, so the events
at $`t_j`$ are included in the estimate of $`\hat{S}(t_j)`$.

The variance or standard error of the Kaplan–Meier estimator is often
calculated with the Greenwood formula. This formula is derived from a
binomial argument, so extension to the weighted case is ad hoc.
Alternatively, Tsiatis (1981) proposes a slightly different formula
based on a counting process argument which includes the weighted case.

**Greenwood variance**
``` math
\textrm{Var}\!\left\{\hat{S}(t)\right \}
= \hat{S}(t)^2\sum_{t_j \le t} \frac{d_j}{Y_j\,(Y_j - d_j)}
```

**Tsiatis variance**
``` math
\textrm{Var}\!\left\{\hat{S}(t)\right \}
= \hat{S}(t)^2\sum_{t_j \le t} \frac{d_j}{Y_j^2}
```

Xie and Liu (2005) considered Greenwood-type variance for the weighted
Kaplan–Meier estimator. Suppose that the survival data consists of
$`\{T_i, d_i, w_i\}`$ ($`i=1,...,N`$), an independent sample of
right-censored survival data with weights. Let $`d_j^w`$ and $`Y_j^w`$
be the weighted number of failures and the weighted number at risk,
respectively, at time $`t_j`$. The weighted Kaplan–Meier estimator of
the survival function $`S(t)`$ is

``` math
\hat{S}(t) = \prod_{t_j \le t} \left( 1 - \frac{d_j^w}{Y_j^w} \right). 
```
The variance proposed by Xie and Liu (2005) is expressed as an adjusted
Greenwood formula,

**Greenwood variance for weighted Kaplan–Meier**
``` math
\textrm{Var}\!\left\{\hat{S}(t)\right \}
= \hat{S}(t)^2\sum_{t_j \le t} \frac{d_j^w}{M_j (Y_j^w - d_j^w)}, 
```
where $`M_j`$ is an adjustment factor defined as
``` math
M_j = \frac{\left(\sum_{t_i \geq t_j} w_i \right)^2}{\sum_{t_i \geq t_j} w_i^2}. 
```
A Tsiatis-type variance, obtained in the same spirit by replacing the
Greenwood denominator with its first-order counterpart, is calculated as
follows. This extension is not the variance formula proposed by Xie and
Liu (2005).

**Tsiatis variance for weighted Kaplan–Meier**
``` math
\textrm{Var}\!\left\{\hat{S}(t)\right \}
= \hat{S}(t)^2\sum_{t_j \le t} \frac{d_j^w}{M_j Y_j^w}
```

These weighted Greenwood/Tsiatis expressions provide background. In the
current
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
API, supplied weights select IF inference by default, and the supported
weighted variance choices are `"if"` and `"fay"`. Requesting another
weighted variance method produces a warning and falls back to IF.

## Formulas for the Aalen-Johansen estimator

Competing risks (`outcome.type = "competing-risk"`) arise in studies in
which individuals are exposed to two or more mutually exclusive failure
events. When a failure occurs, we observe the time to event $`T`$ and
the cause of failure $`\epsilon`$. Suppose that $`\epsilon=1`$ and
$`\epsilon=2`$ represent the event of interest and the competing risk,
respectively. Let $`d_{jk}`$ be the number of failures of cause $`k`$ at
time $`t_j`$, and now the total number of failures at $`t_j`$ is
$`d_j = d_{j1} + d_{j2}`$.

The Aalen-Johansen estimator of CIF for cause $`\epsilon=k`$ is as
below.

**Aalen-Johansen estimator**

``` math
\hat{F}_k(t) = \sum_{t_j \le t} \frac{d_{jk}}{Y_j}\,\hat{S}(t_{j-1}).
```

where $`\hat{S}(t)`$ is the overall survival function.

Two variance estimators of the Aalen-Johansen estimator are commonly
used: one based on counting process theory (Aalen, 1978) and the other
based on the delta method.

**Aalen variance**
``` math
\begin{aligned}
\textrm{Var}\!\left\{\hat{F}_k(t)\right \}
&= \sum_{t_j \le t}
\bigl[\hat{F}_k(t) - \hat{F}_k(t_j)\bigr]^2
\frac{d_j}{(Y_j-1)(Y_j - d_j)}
\\[2pt]
&\quad + \sum_{t_j \le t}
\hat{S}^2(t_{j-1})
\frac{d_{jk}\,(Y_j - d_{jk})}{Y_j^2\,(Y_j-1)}
\\[2pt]
&\quad - 2 \sum_{t_j \le t}
\bigl[\hat{F}_k(t) - \hat{F}_k(t_j)\bigr]\,
\hat{S}(t_{j-1})
\frac{d_{jk}\,(Y_j - d_{jk})}{Y_j\,(Y_j - d_j)\,(Y_j-1)}.
\end{aligned}
```

**Delta method variance**
``` math
\begin{aligned}
\textrm{Var}\!\left\{\hat{F}_k(t)\right \}
&= \sum_{t_j \le t}
\bigl[\hat{F}_k(t) - \hat{F}_k(t_j)\bigr]^2
\frac{d_j}{Y_j\,(Y_j - d_j)}
\\[2pt]
&\quad + \sum_{t_j \le t}
\hat{S}^2(t_{j-1})
\frac{d_{jk}\,(Y_j - d_{jk})}{Y_j^3}
\\[2pt]
&\quad - 2 \sum_{t_j \le t}
\bigl[\hat{F}_k(t) - \hat{F}_k(t_j)\bigr]\,
\hat{S}(t_{j-1})
\frac{d_{jk}}{Y_j^2}.
\end{aligned}
```

**Variance based on influence functions**

For the unweighted Aalen-Johansen estimator, Zhang and Fine (2008) give
the asymptotic linear expansion
``` math
  n^{1/2}\{\hat F_k(t) - F_k(t)\}
  = n^{-1/2} \sum_{i=1}^n IF_{ik}(t) + o_p(1),
```
and the process $`n^{1/2}\{\hat F_k(t) - F_k(t)\}`$ converges weakly to
a tight Gaussian process. Let $`Y(t)`$ denote the sample at-risk process
and let
``` math
  y(t) = \lim_{n\to\infty} n^{-1}Y(t).
```
If $`M_i(t)`$ and $`M_{ik}(t)`$ are the martingales for the total event
count and the cause-$`k`$ event count of observation $`i`$,
respectively, and the cumulative hazards are continuous, the influence
function may be written as
``` math
\begin{aligned}
  IF_{ik}(t)
  &= \int_0^t
       S(u^-)
       \frac{dM_{ik}(u)}{y(u)} \\
  &\quad - \int_0^t
       S(u^-)
       \left\{
         \int_0^{u^-} \frac{dM_i(v)}{y(v)}
       \right\}
       d\Lambda_k(u).
\end{aligned}
```
Equivalently, after changing the order of integration,
``` math
  IF_{ik}(t)
  = \int_0^t \frac{S(u^-)}{y(u)}\,dM_{ik}(u)
    - \int_0^t
      \frac{F_k(t)-F_k(u^-)}{y(u)}\,dM_i(u).
```
A consistent estimator of the asymptotic variance of
$`n^{1/2}\{\hat F_k(t)-F_k(t)\}`$ is
``` math
  \frac{1}{n}\sum_{i=1}^n \{\widehat{IF}_{ik}(t)\}^2,
```
so the variance of $`\hat F_k(t)`$ itself is estimated by
``` math
  \widehat{\textrm{Var}}\!\left\{\hat F_k(t)\right\}
  = \frac{1}{n^2}\sum_{i=1}^n \{\widehat{IF}_{ik}(t)\}^2.
```

For weighted data,
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
differentiates the weighted empirical product limit, including tied
event times, as described below. Internally, the C++ routine stores the
empirical contribution $`\widehat{IF}_{ik}(t)/n`$; consequently, its
probability-scale standard error is calculated as the square root of the
sum of squared empirical contributions. The Zhang and Fine (2008) result
above directly concerns the unweighted estimator. Supplied weights are
held fixed in the implemented weighted calculation, including weights
estimated externally from a propensity-score model; uncertainty from
estimating those weights is not included.

## Fixed-weight IF and Fay–Graubard-type variances

The following calculations are performed separately within each stratum.
Only observations with positive weights enter the variance and
degrees-of-freedom system; write their number as $`m`$. With no supplied
weights, set $`w_i=1`$. Define
``` math
Y_j^w=\sum_i w_iY_{ij},\qquad
\hat h_{jk}=\frac{\sum_i w_i\Delta N_{ijk}}{Y_j^w},\qquad
\hat h_j=\sum_k\hat h_{jk},
```
where $`Y_{ij}=1(T_i\ge t_j)`$ and
$`\Delta N_{ijk}=1(T_i=t_j,\epsilon_i=k)`$. For survival data there is
one event cause. Set
``` math
r_{ijk}=w_i(\Delta N_{ijk}-Y_{ij}\hat h_{jk}),\qquad
r_{ij}=\sum_k r_{ijk}.
```
Starting from $`a^S_{i0}=a^F_{ik0}=0`$ and $`\hat S_0=1`$, the
uncorrected empirical IF contributions obey
``` math
\begin{aligned}
a^F_{ikj}&=a^F_{ik,j-1}+\hat h_{jk}a^S_{i,j-1}
                  +\hat S_{j-1}r_{ijk}/Y_j^w,\\
a^S_{ij}&=(1-\hat h_j)a^S_{i,j-1}
                  -\hat S_{j-1}r_{ij}/Y_j^w.
\end{aligned}
```
Thus the IF variance (`error = "if"`) is $`\sum_i(a^S_{ij})^2`$ for KM
survival or $`\sum_i(a^F_{ikj})^2`$ for the cause-$`k`$ CIF. These are
contributions on the estimator scale, already divided by the sample size
in the usual IF convention. The recursion differentiates the actual
product limit rather than a continuous-hazard approximation.

For `error = "fay"`, replace each residual $`r_{ijk}`$ and $`r_{ij}`$ in
this recursion by $`H_{ij}r_{ijk}`$ and $`H_{ij}r_{ij}`$, respectively,
where
``` math
\ell_{ij}=\frac{w_iY_{ij}}{Y_j^w},\qquad
H_{ij}=\{1-\min(0.75,\ell_{ij})\}^{-1/2}.
```
The FG variance is the sum of squares of these corrected contributions.
The estimate itself is unchanged. The leverage cap is fixed at 0.75: the
returned `fay.bound` records this value and is not a user-supplied
argument. `"fg"` and `"fay-graubard"` are aliases for `"fay"`.

The returned influence-function matrices, when requested with
`report.influence.function = TRUE`, always contain the **uncorrected**
contributions for the event probability `1 - fit$surv`, even when the
reported SE uses the FG correction. Zero-weight observations have zero
contributions.

## Degrees of freedom for IF-t and FG-t

Both Student constructions use a trace ratio, with a choice of
covariance allocation for IF-t. To specify the calculation at time
$`t`$, stack the individual weighted hazard scores
$`U_{h_j i}=w_i\{\Delta\mathbf N_i(t_j)-Y_{ij}\mathbf h_j\}`$ over event
times $`t_j\le t`$, using hazard logits with the no-event category as
baseline. There are no additional time- or cause-specific equation
multipliers. Only positive-probability coordinates are retained; at
complete exit a positive event cause is used as baseline instead. At the
fitted parameter $`\hat\vartheta`$, define
``` math
\Omega_i=-\frac{\partial U_i}{\partial\vartheta^\top},\qquad
B=\sum_i\Omega_i,\qquad c=\nabla_\vartheta\hat\theta(t),
```
where $`\hat\theta(t)`$ is the KM survival estimate or the AJ CIF. This
production calculation contains hazard scores only, with no
propensity-score estimating equations.

**Deletion allocation (FDF).** The default `df.method = "deletion"` uses
``` math
\rho_i=c^\top\{(B-\Omega_i)^{-1}-B^{-1}\}c,\qquad
\alpha_i=\rho_i/\sum_j\rho_j.
```
These are derivative-based deletions at the original fit, not delete-one
refits. They require invertible relevant deleted blocks and a
nonnegative, nonzero allocation. The allocation can depend on the
normalization of the estimating equations, so the score convention above
is part of the method definition.

**Equal allocation (EDF).** With `error = "if"`, `df.method = "equal"`
instead sets $`\alpha_i=1/m`$. For this uncorrected IF construction, the
resulting degrees of freedom are invariant to fixed invertible
transformations of the same score system, with the corresponding
transformation of derivatives. This statement holds for a fixed active
set and does not cover arbitrary parameter-dependent changes of
estimating equations. EDF is generally not $`m-1`$, and its
finite-sample coverage need not equal FDF coverage. Equal allocation is
not available for FG-t.

For either allocation, let $`H_i=I`$ for IF-t, and use the diagonal
matrix containing the factors $`H_{ij}`$ above for FG-t. Form
``` math
\begin{aligned}
T_i^H&=H_iU_i,& M_H&=\sum_iT_i^H(T_i^H)^\top,\\
L^H_{ij}&=c^\top B^{-1}H_i\{1(i=j)I-\Omega_iB^{-1}\},\\
Q_H&=L^H\operatorname{blockdiag}(\alpha_1M_H,\ldots,\alpha_mM_H)(L^H)^\top,\\
\hat\nu_H(t)&=\frac{\{\operatorname{tr}(Q_H)\}^2}{\operatorname{tr}(Q_H^2)}.
\end{aligned}
```
Here $`L^H`$ is the block matrix with blocks $`L^H_{ij}`$. FG-t uses the
original, uncorrected derivatives $`B,\Omega_i`$ and deletion
allocation; $`H_i`$ is not differentiated. Both $`M_H`$ and $`L^H`$
include the correction at the displayed positions. Consequently, IF-t
(deletion) and FG-t generally have different degrees of freedom as well
as different variances.

The trace ratio is a Satterthwaite moment match for a specified Gaussian
working quadratic form. It does **not** equate the actual sampling
moments of the estimated variance or establish an exact Student
distribution. FDF denotes the FG-type deletion allocation; EDF denotes
the equal-allocation alternative. Neither name implies a uniform
improvement in small-sample coverage.

## Confidence interval options

**Defaults and method selection.** With `weights = NULL` and
`conf.dist = "normal"`, the default is Greenwood for KM and delta for
AJ. With supplied weights, the default is IF. Specifying
`conf.dist = "t"` without `error` selects IF. The variance, reference
distribution, and transformation are separate choices:

| Construction | `error` | `conf.dist` | `df.method`            |
|:-------------|:--------|:------------|:-----------------------|
| IF-normal    | `"if"`  | `"normal"`  | Not used               |
| IF-t (FDF)   | `"if"`  | `"t"`       | `"deletion"` (default) |
| IF-t (EDF)   | `"if"`  | `"t"`       | `"equal"`              |
| FG-normal    | `"fay"` | `"normal"`  | Not used               |
| FG-t         | `"fay"` | `"t"`       | `"deletion"`           |

Student intervals require `error = "if"` or `"fay"`. Student and FG
calculations require `engine = "calculateAJ_Rcpp"` (the default) or
`"auto"`; they are unavailable with `engine = "calculateKM"`. Degrees of
freedom are calculated separately for each time and stratum; this
calculation can be expensive for large risk sets.

**Transformations and probability scale.** Let $`p=\hat S(t)`$ for KM
and $`p=1-\hat F_1(t)`$ for AJ, matching the returned `fit$surv`. At an
interior estimate with positive variance, the transformed interval is
obtained from
``` math
g(p)\ \pm\ q\,|g'(p)|\,\widehat{\operatorname{SE}}(p),\qquad
q=\begin{cases}
z_{(1+\mathrm{conf.int})/2},&\text{normal},\\
t_{\hat\nu(t),(1+\mathrm{conf.int})/2},&\text{Student}.
\end{cases}
```
The endpoints are mapped back using $`g^{-1}`$ and ordered on the
probability scale. Bounds are restricted to the transformation range
before inversion where necessary, and to $`[0,1]`$. The supported
transformations are:

| `conf.type` | $`g(p)`$ |
|:---|:---|
| `"plain"` | $`p`$ |
| `"log"` | $`\log p`$ |
| `"log-log"` | $`\log(-\log p)`$ |
| `"logit"` | $`\log\{p/(1-p)\}`$ |
| `"arcsin"` (default; aliases `"arcsine-square root"`, `"a"`) | $`\arcsin\sqrt p`$ |

For AJ, `fit$lower` and `fit$upper` bound $`1-\hat F_1(t)`$; the CIF
limits are `1 - fit$upper` and `1 - fit$lower`. Thus `"log-log"` uses
$`\log[-\log\{1-\hat F_1(t)\}]`$ for AJ, not $`\log[-\log\hat F_1(t)]`$.
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
displays the CIF limits when plotting risk. Passing `conf.type = "none"`
(or `"n"`) skips intervals. By default, `std.err` is on the probability
scale; `report.survfit.std.err = TRUE` requests the log scale of the
stored `surv` component instead.

**Boundary and failure rules.** On the C++ IF/FG path, a valid estimate
at 0 or 1 gives the corresponding singleton interval. This is an
endpoint convention, not a guarantee of boundary coverage. Interior
estimates with zero or invalid variance, unavailable Student degrees of
freedom, or evaluation times beyond the positive-weight follow-up
support have missing limits. Failed Student intervals are not replaced
by normal intervals. The fit records time- and stratum-specific `df`,
`df.reason`, `critical.value`, and `conf.status`;
`summary(fit, times = ..., extend = TRUE)` extracts the saved inference
results and preserves missing intervals outside support.

These are pointwise confidence intervals, not simultaneous confidence
bands. The optional one-sided test requested through `time.point` and
`null.hypothesis` retains its normal reference distribution even when
`conf.dist = "t"` is used for confidence intervals.

## Examples

The following small dataset is illustrative, not clinical data. The same
inference options apply to KM and AJ, with or without fixed weights.

``` r

library(cifmodeling)
example_data <- data.frame(
  time = 1:12,
  status = c(1, 0, 2, 1, 2, 0, 1, 2, 0, 1, 0, 0),
  weight = c(1, 1.5, 0.8, 1.2, 1, 0.7, 1.3, 1, 0.9, 1.1, 1, 1)
)

# AJ with IF-t and equal-allocation degrees of freedom.
fit_edf <- cifcurve(
  Event(time, status) ~ 1, data = example_data, weights = "weight",
  outcome.type = "competing-risk", error = "if",
  conf.dist = "t", df.method = "equal", conf.type = "log-log"
)
at_five <- summary(fit_edf, times = 5)
data.frame(
  time = at_five$time, CIF = 1 - at_five$surv,
  SE = at_five$std.err, lower = 1 - at_five$upper,
  upper = 1 - at_five$lower, df = at_five$df,
  status = at_five$conf.status
)
```

    ##   time    CIF        SE      lower     upper       df status
    ## 1    5 0.1904 0.1220426 0.04291352 0.6383596 10.94327     ok

``` r

# KM for the composite of either event, with FG-t.
fit_fg <- cifcurve(
  Event(time, status > 0) ~ 1, data = example_data,
  outcome.type = "survival", error = "fay",
  conf.dist = "t", conf.type = "arcsin"
)
```

For plotting Student intervals, calculate them first and pass the fitted
object to
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md).
Its formula interface does not forward `conf.dist` or `df.method` to
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md).

``` r

cifplot(fit_edf, outcome.type = "competing-risk", add.risktable = FALSE)
```

See [Full lists of
arguments](https://gestimation.github.io/cifmodeling/articles/v7_arguments.md)
for the API options and [Reporting cumulative
incidence](https://gestimation.github.io/cifmodeling/articles/v8_reporting_cif.md)
for extracting estimates and constructing tables.

## References

- Xie J, Liu C. Adjusted Kaplan–Meier estimator and log-rank test with
  inverse probability of treatment weighting for survival data.
  *Statistics in Medicine*. 2005, 24(20), 3089–3110.
  [doi:10.1002/sim.2174](https://doi.org/10.1002/sim.2174)
