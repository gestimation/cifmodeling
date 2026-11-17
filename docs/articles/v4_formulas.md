# Computational formulas in cifcurve()

## How to use this vignette

[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
computes Kaplan–Meier (KM) survival curves and Aalen–Johansen (AJ)
cumulative incidence functions (CIFs), along with pointwise standard
errors and confidence intervals. This vignette documents **how
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
computes its estimates, variances, and intervals**. For analysis
workflows and publication tables, see [Worked
examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md)
and [Reporting
results](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.md).

The main arguments control distinct parts of the calculation:

| Argument | What it controls | Options or default |
|:---|:---|:---|
| `outcome.type` | Estimated probability | `"survival"` (KM) or `"competing-risk"` (AJ) |
| `weights` | Case weights | Optional, treated as fixed |
| `error` | Standard-error calculation | Greenwood/Tsiatis (KM); delta/Aalen (AJ); `"if"` or `"fay"` (both) |
| `conf.dist` | CI reference distribution | `"normal"` (default) or `"t"` |
| `conf.type` | CI transformation | `"arcsine-square root"` (default), `"plain"`, `"log"`, `"log-log"`, `"logit"` |
| `df.method` | Working df for Student CIs | `"deletion"` (FDF, default) or `"equal"` (EDF, IF-t only) |
| `report.influence.function` | Return empirical IF contributions | `FALSE` (default) |
| `engine` | Computation back end | `"calculateAJ_Rcpp"` (default) |
| `n.risk.type` | Displayed risk-set counts | `"weighted"` (default), `"unweighted"`, or `"ess"`; estimates unchanged |

Here, **IF** denotes influence function and **FG** denotes Fay–Graubard
(not Fine–Gray regression). **FDF** denotes Fay-type working degrees of
freedom, and **EDF** denotes equal-allocation working degrees of
freedom. FDF and EDF describe working-df allocations, not variance
corrections.

With **unweighted** data, `error = NULL` selects Greenwood variance for
KM and delta-method variance for AJ. When case weights are supplied, the
default is `error = "if"`; the available weighted variance methods are
`"if"` and `"fay"`. Setting `conf.dist = "t"` without specifying `error`
also selects `"if"` and computes working degrees of freedom. None of
`error`, `conf.dist`, `conf.type`, or `df.method` changes the point
estimate; they affect its uncertainty estimates or confidence limits.

For competing-risk outcomes, `fit$surv` contains **one minus the event-1
CIF**, not overall event-free survival. To obtain CIF confidence limits,
complement the stored limits and reverse their order. All supplied
weights, including externally estimated IPTW, are treated as **fixed**;
uncertainty from propensity-score estimation is not included. The
default C++ engine handles tied event times and implements the IF and FG
calculations described below.

The sections below cover KM/AJ formulas, the implemented IF and
Fay–Graubard calculations, working degrees of freedom, and
confidence-interval transformations. Weighted Greenwood/Tsiatis formulas
are included for comparison; they are **not** available as weighted
`error` options in
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md).

## Kaplan–Meier estimation and classical variances

For `outcome.type = "survival"`, `fit$surv` contains the KM survival
estimate. With unweighted data, `error = "greenwood"` is the default;
`error = "tsiatis"` selects the alternative classical standard error.
Let $`t_1 < t_2 < \cdots < t_D`$ denote the distinct event times,
$`Y_j`$ the number at risk just before $`t_j`$, and $`d_j`$ the number
of events at $`t_j`$. The KM estimator of $`S(t)`$ is

**Kaplan–Meier estimator**

``` math
\hat{S}(t) = \prod_{t_j \le t} \left(1 - \frac{d_j}{Y_j} \right)
```
The estimator is right-continuous, so events at $`t_j`$ contribute to
$`\hat S(t_j)`$.

For unweighted KM estimates,
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
offers Greenwood and Tsiatis variance formulas (Tsiatis, 1981). The
formulas differ in their denominators:

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

For weighted data, let $`w_i`$ be the supplied case weight, and let
$`d_j^w`$ and $`Y_j^w`$ be the weighted event count and risk-set size at
$`t_j`$. The weighted KM estimator is

``` math
\hat{S}(t) = \prod_{t_j \le t} \left( 1 - \frac{d_j^w}{Y_j^w} \right). 
```
Xie and Liu (2005) proposed the following adjusted Greenwood expression:

**Greenwood variance for weighted Kaplan–Meier**
``` math
\textrm{Var}\!\left\{\hat{S}(t)\right \}
= \hat{S}(t)^2\sum_{t_j \le t} \frac{d_j^w}{M_j (Y_j^w - d_j^w)}, 
```
where $`M_j`$ is an adjustment factor defined as
``` math
M_j = \frac{\left(\sum_{t_i \geq t_j} w_i \right)^2}{\sum_{t_i \geq t_j} w_i^2}. 
```
For comparison, the analogous Tsiatis-type expression uses the
corresponding first-order denominator. It is not the variance formula
proposed by Xie and Liu (2005).

**Tsiatis variance for weighted Kaplan–Meier**
``` math
\textrm{Var}\!\left\{\hat{S}(t)\right \}
= \hat{S}(t)^2\sum_{t_j \le t} \frac{d_j^w}{M_j Y_j^w}
```

Neither weighted formula is selectable through `error` in the current
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
implementation. With weights, use `error = "if"` or `"fay"`. Requesting
another variance method issues a warning and falls back to IF.

## Aalen–Johansen estimation and classical variances

With `outcome.type = "competing-risk"`,
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
estimates the event-1 CIF using AJ. For unweighted data,
`error = "delta"` is the default, while `error = "aalen"` selects the
alternative classical variance. In the formulas, cause 1 is the event of
interest, cause 2 is the competing event, $`d_{jk}`$ counts cause-$`k`$
events at $`t_j`$, and $`d_j=d_{j1}+d_{j2}`$. AJ weights each
cause-specific event increment by the probability of remaining
event-free just before that time.

For cause $`k`$, the estimator is

**Aalen-Johansen estimator**

``` math
\hat{F}_k(t) = \sum_{t_j \le t} \frac{d_{jk}}{Y_j}\,\hat{S}(t_{j-1}).
```

Here $`\hat S(t_{j-1})=\hat S(t_j-)`$ is overall event-free survival
immediately before $`t_j`$, with $`\hat S(t_0)=1`$. It differs from the
`$surv` component of a competing-risk fit, which stores
$`1-\hat F_1(t)`$.

For **unweighted** AJ estimates, `error` selects either the Aalen (1978)
or delta-method variance. The implemented formulas are:

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

## Influence-function representation

Set `error = "if"` to estimate variance from individual
influence-function (IF) contributions. Add
`report.influence.function = TRUE` to return the uncorrected empirical
contributions. The classical AJ representation below provides
background; the next section describes the **finite-jump, fixed-weight
calculation implemented in
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)**.
For unweighted AJ estimates, Zhang and Fine (2008) give the asymptotic
linear expansion
``` math
  n^{1/2}\{\hat F_k(t) - F_k(t)\}
  = n^{-1/2} \sum_{i=1}^n IF_{ik}(t) + o_p(1),
```
The process $`n^{1/2}\{\hat F_k(t) - F_k(t)\}`$ also converges weakly to
a tight Gaussian process. Let $`Y(t)`$ be the sample at-risk process,
and define
``` math
  y(t) = \lim_{n\to\infty} n^{-1}Y(t).
```
For continuous cumulative hazards, let $`M_i(t)`$ and $`M_{ik}(t)`$ be
the martingales for all events and cause-$`k`$ events, respectively, for
observation $`i`$. The influence function is then
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

The preceding representation applies to unweighted AJ estimation under
the stated assumptions, including continuous cumulative hazards. By
contrast,
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
differentiates the **empirical weighted product-limit estimator**,
allowing tied event times. The C++ routine stores contributions on the
estimator scale (conventionally $`\widehat{IF}_{ik}(t)/n`$); their
squared sum estimates the variance of $`\hat F_k(t)`$. Externally
supplied weights remain fixed, so propensity-score estimation is not
included. Deng and Wang (2025) derive influence functions that account
for estimating propensity scores.

## Fixed-weight IF and Fay–Graubard variances

This section details the calculations behind `error = "if"` and
`error = "fay"`. Both retain the **same KM/AJ point estimate** but use
different contributions to estimate variance. The C++ engine
differentiates the finite-jump product limit, with or without weights,
within each stratum. Only positive-weight observations enter the
variance and degrees-of-freedom calculations; let $`m`$ be their number.
Without supplied weights, set $`w_i=1`$. Define
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
With `error = "if"`, the variance is $`\sum_i(a^S_{ij})^2`$ for KM
survival or $`\sum_i(a^F_{ikj})^2`$ for the cause-$`k`$ CIF. The
contributions are already on the estimator scale, with the usual $`1/n`$
factor incorporated. This recursion differentiates the finite-jump
product limit used by the software, rather than relying on a
continuous-hazard approximation.

With `error = "fay"`, the Fay–Graubard correction (Fay and Graubard
2001) applies the same recursion after replacing $`r_{ijk}`$ and
$`r_{ij}`$ by $`H_{ij}r_{ijk}`$ and $`H_{ij}r_{ij}`$, respectively,
where
``` math
\ell_{ij}=\frac{w_iY_{ij}}{Y_j^w},\qquad
H_{ij}=\{1-\min(0.75,\ell_{ij})\}^{-1/2}.
```
The Fay–Graubard (FG) variance sums the squared corrected contributions;
**the KM/AJ point estimate is unchanged**. The leverage cap is fixed at
0.75 and reported as `fay.bound`; users cannot set it directly. Both
`"fg"` and `"fay-graubard"` are aliases for `"fay"`.

Even when `error = "fay"`, `report.influence.function = TRUE` returns
**uncorrected** contributions for the event probability `1 - fit$surv`.
Observations with zero weight contribute zero.

## Working degrees of freedom for IF-t and FG-t

With `conf.dist = "t"`,
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
uses a time- and stratum-specific t critical value rather than a normal
one. For `error = "if"`, `df.method` selects deletion allocation
(`"deletion"`, FDF) or equal allocation (`"equal"`, EDF). For
`error = "fay"`, only deletion allocation is available (FG-t). The trace
ratio below provides **working**, not exact finite-sample Student,
degrees of freedom.

The calculation starts with individual weighted hazard scores
$`U_{h_j i}=w_i\{\Delta\mathbf N_i(t_j)-Y_{ij}\mathbf h_j\}`$, stacked
over event times $`t_j\le t`$. Hazard logits use the no-event category
as the baseline, except at complete exit, when a positive cause becomes
the baseline; only positive-probability coordinates are retained. No
additional time- or cause-specific equation multipliers are used. At the
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
These quantities use derivatives at the original fit; they do not
require refitting after deleting each observation. The relevant deleted
blocks must be invertible, and the allocation must be nonnegative and
nonzero. Because the allocation can depend on how estimating equations
are normalized, the score convention above is part of the method
definition.

**Equal allocation (EDF).** For `error = "if"`, `df.method = "equal"`
sets $`\alpha_i=1/m`$. In the uncorrected IF calculation, the resulting
degrees of freedom are invariant to fixed invertible score
transformations when the derivatives are transformed consistently. This
holds for a fixed active set, not for arbitrary parameter-dependent
transformations. EDF degrees of freedom generally **do not equal**
$`m-1`$, and EDF is unavailable with FG-t.

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
Here $`L^H`$ has blocks $`L^H_{ij}`$. FG-t uses the original,
uncorrected derivatives $`B,\Omega_i`$ and deletion allocation; it does
not differentiate $`H_i`$. The correction enters $`M_H`$ and $`L^H`$ at
the indicated positions. Thus IF-t (deletion) and FG-t can differ in
both variance and degrees of freedom.

The trace ratio matches moments of a Gaussian working quadratic form,
following the Satterthwaite principle. It does **not** imply an exact t
distribution or guarantee improved small-sample coverage. FDF uses
Fay-type deletion allocation; EDF uses equal allocation. Even with the
same IF variance, different allocations may produce different critical
values.

## Confidence interval options

### Method combinations

The following argument combinations select the available IF/FG interval
methods. Setting `conf.dist = "t"` without `error` selects IF, even for
unweighted data. `df.method` affects IF-t, but not normal-reference
intervals.

| Construction | `error` | `conf.dist` | `df.method`            |
|:-------------|:--------|:------------|:-----------------------|
| IF-normal    | `"if"`  | `"normal"`  | Not used               |
| IF-t (FDF)   | `"if"`  | `"t"`       | `"deletion"` (default) |
| IF-t (EDF)   | `"if"`  | `"t"`       | `"equal"`              |
| FG-normal    | `"fay"` | `"normal"`  | Not used               |
| FG-t         | `"fay"` | `"t"`       | `"deletion"`           |

Student intervals require `error = "if"` or `"fay"` and the C++ engine
(`engine = "calculateAJ_Rcpp"`, the default, or `"auto"`). The legacy
`"calculateKM"` engine does not support them. Degrees of freedom are
calculated separately for each time and stratum and may be costly for
large risk sets.

### Transformations and probability scale

The transformation is applied to the **stored** probability:
$`p=\hat S(t)`$ for KM, but $`p=1-\hat F_1(t)`$ for AJ (`fit$surv`). At
an interior estimate with positive variance, the transformed limits are
obtained from
``` math
g(p)\ \pm\ q\,|g'(p)|\,\widehat{\operatorname{SE}}(p),\qquad
q=\begin{cases}
z_{(1+\mathrm{conf.int})/2},&\text{normal},\\
t_{\hat\nu(t),(1+\mathrm{conf.int})/2},&\text{Student}.
\end{cases}
```
The limits are mapped back through $`g^{-1}`$ and ordered on the
probability scale. Where needed, they are restricted to the
transformation’s range before inversion and to $`[0,1]`$ afterward. The
supported transformations are:

| `conf.type` | $`g(p)`$ |
|:---|:---|
| `"plain"` | $`p`$ |
| `"log"` | $`\log p`$ |
| `"log-log"` | $`\log(-\log p)`$ |
| `"logit"` | $`\log\{p/(1-p)\}`$ |
| `"arcsin"` (default; aliases `"arcsine-square root"`, `"a"`) | $`\arcsin\sqrt p`$ |

The default is the **arcsine–square-root CI transformation**. This
option is distinct from `type.y`, which controls the vertical plotting
scale.

For AJ, `fit$lower` and `fit$upper` are limits for $`1-\hat F_1(t)`$, so
the CIF limits are `1 - fit$upper` and `1 - fit$lower`. Accordingly,
`"log-log"` applies $`\log[-\log\{1-\hat F_1(t)\}]`$, not
$`\log[-\log\hat F_1(t)]`$. When displaying risk,
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
uses the CIF limits. `conf.type = "none"` (or `"n"`) omits intervals. By
default, `std.err` is on the probability scale;
`report.survfit.std.err = TRUE` instead uses the log scale of the stored
`surv` component.

### Boundaries, failures, and one-sided tests

For the C++ IF/FG calculations, a valid estimate of 0 or 1 receives a
singleton interval. This is a boundary convention, **not** a coverage
guarantee. Limits are missing for interior estimates with zero or
invalid variance, unavailable Student degrees of freedom, or times
beyond positive-weight follow-up. A failed Student interval does **not**
fall back to a normal interval. The fit records `df`, `df.reason`,
`critical.value`, and `conf.status` by time and stratum;
`summary(fit, times = ..., extend = TRUE)` preserves these diagnostics,
including missing limits outside follow-up support.

These intervals are **pointwise**, not simultaneous confidence bands.
The optional one-sided fixed-time probability test using `time.point`
and `null.hypothesis` retains a normal reference even when
`conf.dist = "t"` is selected for intervals. Its default alternative is
higher survival for KM and lower event-1 cumulative incidence for AJ.

## Worked examples

These short examples demonstrate
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
options rather than a full analysis workflow. With this small synthetic
dataset, working degrees of freedom or CIs may be unavailable at some
times; `conf.status` indicates why. For complete analyses and tables,
see [Worked
examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md)
and [Reporting
results](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.md).

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
  critical.value = at_five$critical.value,
  df.reason = at_five$df.reason,
  status = at_five$conf.status
)
```

    ##   time    CIF        SE      lower     upper       df critical.value df.reason
    ## 1    5 0.1904 0.1220426 0.04291352 0.6383596 10.94327       2.202378        ok
    ##   status
    ## 1     ok

``` r

# KM for the composite of either event, with FG-t.
example_data$composite_event <- as.integer(example_data$status > 0)
fit_fg <- cifcurve(
  Event(time, composite_event) ~ 1, data = example_data,
  outcome.type = "survival", error = "fay",
  conf.dist = "t", conf.type = "arcsin"
)
```

The first output table shows the event-1 CIF and its CI on the
probability scale, along with degrees of freedom, the critical value,
and inference status. The same diagnostics are stored as `fit_edf$df`,
`fit_edf$df.reason`, `fit_edf$critical.value`, and `fit_edf$conf.status`
on the original event-time grid; use `summary(fit_edf, times = ...)` to
extract them at a chosen time. In the second fit, `error = "fay"`
changes the variance, while `conf.dist = "t"` changes the reference
critical value.

To plot Student intervals, first fit them with
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
and pass the result to
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md).
The
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
**formula interface** does not forward `conf.dist` or `df.method` to
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md).

``` r

cifplot(fit_edf, outcome.type = "competing-risk", add.risktable = FALSE)
```

See [Function
reference](https://gestimation.github.io/cifmodeling/articles/v9_arguments.md)
for the full set of arguments and [Reporting cumulative
incidence](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.html#reporting-cumulative-incidence)
for extracting estimates and building tables.

## References

- Aalen OO. Nonparametric inference for a family of counting processes.
  *Ann Statist*. 1978;6(4):701–726.
  [doi:10.1214/aos/1176344247](https://doi.org/10.1214/aos/1176344247)
- Tsiatis AA. A large sample study of Cox’s regression model. *Ann
  Statist*. 1981;9(1):93–108.
  [doi:10.1214/aos/1176345335](https://doi.org/10.1214/aos/1176345335)
- Zhang M-J, Fine JP. Summarizing differences in cumulative incidence
  functions. *Stat Med*. 2008;27(24):4939–4949.
  [doi:10.1002/sim.3339](https://doi.org/10.1002/sim.3339)
- Xie J, Liu C. Adjusted Kaplan–Meier estimator and log-rank test with
  inverse probability of treatment weighting for survival data. *Stat
  Med* 2005;24(20):3089–110.
  [doi:10.1002/sim.2174](https://doi.org/10.1002/sim.2174)
- Deng Y, Wang R. Adjusted Nelson–Aalen estimators by inverse treatment
  probability weighting with an estimated propensity score. *Statistics
  in Medicine*. 2025;44(10–12):e70085.
  [doi:10.1002/sim.70085](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.70085)
- Fay MP, Graubard BI. Small-sample adjustments for Wald-type tests
  using sandwich estimators. *Biometrics*. 2001;57(4):1198–1206.
  [doi:10.1111/j.0006-341X.2001.01198.x](https://onlinelibrary.wiley.com/doi/10.1111/j.0006-341X.2001.01198.x)
