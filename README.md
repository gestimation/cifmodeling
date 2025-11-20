
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

# Visualization and Modeling of Survival and Competing Risks Data — cifmodeling

## Quick start

This package is a compact, high-level extension of the existing
`survival` ecosystem. It provides a unified interface for Kaplan-Meier
and Aalen-Johansen curves, modern visualization, and direct polytomous
regression for survival and competing risks data.

``` r
library(cifmodeling)
data(diabetes.complications)
cifplot(Event(t,epsilon)~fruitq, data=diabetes.complications, 
        outcome.type="competing-risk", panel.per.event=TRUE)
```

<div class="figure">

<img src="man/figures/README-example01-1-1.png" alt="Aalen-Johansen cumulative incidence curves from cifplot()" width="70%" />
<p class="caption">

Aalen-Johansen cumulative incidence curves from cifplot()
</p>

</div>

In competing risks data, censoring is often coded as 0, the event of
interest as 1, and competing risks as 2. In the `diabetes.complications`
data frame, `epsilon` follows this convention. With
`panel.per.event = TRUE`, `cifplot()` visualizes the cumulative
incidence functions (CIFs), with the CIF of diabetic retinopathy
(`epsilon = 1`) shown on the left and the CIF of macrovascular
complications (`epsilon = 2`) on the right.

## Why cifmodeling?

- **Unified interface** for Kaplan–Meier and Aalen–Johansen curves, with
<<<<<<< HEAD
  survival and competing risks handled by the same syntax.

- **Publication-ready graphics** built on `ggsurvfit` and `ggplot2`,
  including risk/estimate tables,
  censoring/competing-risks/intercurrent-events marks, and multi-panel
  layouts. \<\<\<\<\<\<\< HEAD

- # **Coherent regression models** for CIFs, targeting familiar effect measures such as risk ratios (RR), odds ratios (OR), and subdistribution hazard ratios (SHR).

- **Tidy summaries and reporting**: regression results from `polyreg()`
  support `generics::tidy()`, `glance()`, and `augment()`, which
  integrate smoothly with `modelsummary` and other broom-style tools.

- **Coherent regression models** of CIFs, targeting familiar effect
  measures (risk ratios, odds ratios, and subdistribution hazard
  ratios). Modeling the nuisance structure using polytomous log odds
  products ensures that the sum of cause-specific CIFs does not exceed
  one. \>\>\>\>\>\>\> 52e2a3008aaa2cc10374d4a220ad7b6c757aba48
=======
  survival and competing risks handled by the same `Event()` + formula +
  data syntax.
- **Effects on the CIF scale**: while Fine-Gray models subdistribution
  hazards, `polyreg()` directly targets ratios of CIFs (risk ratios,
  odds ratios, subdistribution hazard ratios), so parameters align
  closely with differences seen in CIF curves.
- **Coherent, joint modeling of all competing events**: `polyreg()`
  models all cause-specific CIFs together, parameterizing the nuisance
  structure with polytomous log odds products and enforcing that their
  CIFs sum to at most one.
- **Tidy summaries and reporting**: support for `generics::tidy()`,
  `glance()`, and `augment()`, which integrate `polyreg()` smoothly with
  `modelsummary` and other broom-style tools.
- **Publication-ready graphics** built on `ggsurvfit` and `ggplot2`,
  including at-risk/CIF+CI tables,
  censoring/competing-risk/intercurrent-event marks, and multi-panel
  layouts.
>>>>>>> ee5005b4869f1b0b293ed0d404c0115d2926b6ed

## Tools for survival and competing risks analysis

In clinical and epidemiological research, analysts often need to handle
censoring, competing risks, and intercurrent events (e.g. treatment
switching), but existing R packages typically separate these tasks
across different interfaces. `cifmodeling` provides a unified,
publication-ready toolkit that integrates nonparametric estimation,
regression modelling, and visualization for survival and competing risks
data. The package is centered around three tightly connected functions:

- `cifplot()` typically generates a survival or CIF curve with marks
  that represent censoring, competing risks and intercurrent events.
  **Multiple standard error (SE) estimators and confidence interval (CI)
  methods valid for unweighted and weighted data are supported.** The
  visualization is built on top of `ggsurvfit` and `ggplot2`.

- `cifpanel()` creates a multi-panel figure for survival/CIF curves,
  arranged either **in a grid layout or as an inset overlay**.

- `polyreg()` fits **coherent regression models** of CIFs using
  polytomous log odds products.

These functions adopt a formula + data syntax, return tidy,
publication-ready outputs, and integrate seamlessly with `ggsurvfit` and
`modelsummary` for visualization and reporting.

## Position in the survival ecosystem

Several excellent R packages exist for survival and competing risks
analysis. The **survival** package provides the canonical API for
survival data. In combination with the **ggsurvfit** package,
`survival::survfit()` can produce publication-ready survival plots. For
CIF plots, however, integration in the general ecosystem is less
streamlined. `cifmodeling` fills this gap by offering `cifplot()` for
survival/CIF plots and multi-panel figures via a single, unified
interface.

Beyond providing a unified interface, `cifcurve()` also extends
`survfit()` in a few targeted ways. For unweighted survival data, it
reproduces the standard Kaplan-Meier estimator with **Greenwood and
Tsiatis SEs** and a unified set of CI transformations. For competing
risks data, it computes Aalen-Johansen CIFs with both **Aalen-type and
delta-method SEs**. For weighted survival or competing risks data
(e.g. inverse probability weighting), it implements **influence-function
based SEs** (Deng and Wang 2025) as well as **modified Greenwood- and
Tsiatis-type SEs** (Xie and Liu 2005), which are valid under general
positive weights.

If you need very fine-grained plot customization, you can compute the
estimator and keep a survfit-compatible object with `cifcurve()` (or
supply your own survfit object) and then style it using
`ggsurvfit/ggplot2` layers. In other words:

- use `cifcurve()` for estimation,
- use `cifplot()` / `cifpanel()` for quick, high-quality figures, and
- fall back to the `ggplot` ecosystem when you want full artistic
  control.

The **mets** package is a more specialised toolkit that provides
advanced methods for competing risks analysis. `cifmodeling::polyreg()`
focuses on coherent modelling of all CIFs simultaneously to estimate the
exposure effects in terms of RR/OR/SHR. This coherence can come with
longer runtimes for large problems. If you prefer fitting separate
regression models for each competing event or specifically need the
Fine-Gray models (Fine and Gray 1999) and the direct binomial models
(Scheike, Zhang and Gerds 2008), `mets::cifreg()` and `mets::binreg()`
are excellent choices.

## Installation

The package is implemented in R and relies on `Rcpp`, `nleqslv` and
`boot` for its numerical back-end. The examples in this document also
use `ggplot2`, `ggsurvfit`, `patchwork` and `modelsummary` for
tabulation and plotting. Install the core package and these companion
packages with:

``` r
# Install cifmodeling from GitHub
devtools::install_github("gestimation/cifmodeling")

# Core dependencies
install.packages(c("Rcpp", "nleqslv", "boot"))

# Recommended packages for plotting and tabulation in this README
install.packages(c("ggplot2", "ggsurvfit", "patchwork", "modelsummary"))
```

## Quality control

`cifmodeling` includes an extensive test suite built with **testthat**,
which checks the numerical accuracy and graphical consistency of all
core functions (`cifcurve()`, `cifplot()`, `cifpanel()`, and
`polyreg()`). The estimators are routinely compared against related
functions in **survival**, **cmprsk** and **mets** packages to ensure
consistency. The package is continuously tested on GitHub Actions
(Windows, macOS, Linux) to maintain reproducibility and CRAN-level
compliance.

## An example of competing risks analysis

For the initial illustration, we focus on unadjusted cumulative
incidence of diabetic retinopathy (event 1) and macrovascular
complications (event 2) at 8 years of follow-up. To visualize each
covariate separately when multiple strata are supplied, set
`panel.per.variable = TRUE`. Each variable on the right-hand side is
plotted in its own panel, and the layout can be controlled with
`rows.columns.panel`. The figure below contrasts the CIFs of diabetic
retinopathy for the quartiles `fruitq` and a binary exposure `fruitq1`,
low (Q1) and high (Q2 to 4) intake of fruit, generated by `cifplot()`.
The `add.conf=TRUE` argument adds CIs to the plot. This helps visualize
the statistical uncertainty of estimated probabilities across exposure
levels. When using numeric variables for stratification, discretize them
beforehand with `cut()` or `factor()`. The labels of x-axis (Time) and
y-axis (Cumulative incidence) in these panels are default labels.

``` r
data(diabetes.complications)
diabetes.complications$fruitq1 <- ifelse(
  diabetes.complications$fruitq == "Q1","Q1","Q2 to Q4"
)
cifplot(Event(t,epsilon)~fruitq+fruitq1, data=diabetes.complications, 
        outcome.type="competing-risk",
        add.conf=TRUE, add.censor.mark=FALSE, 
        add.competing.risk.mark=FALSE, panel.per.variable=TRUE)
```

<div class="figure">

<img src="man/figures/README-example04-1-1-1.png" alt="Cumulative incidence curves per each stratification variable" width="70%" />
<p class="caption">

Cumulative incidence curves per each stratification variable
</p>

</div>

In the second figure, **competing-risk marks** are added
(`add.competing.risk.mark = TRUE`) to indicate individuals who
experienced the competing event (macrovascular complications) before
diabetic retinopathy. Here we show a workflow slightly different from
the previous code. First, we compute a survfit-compatible object
`output1` using `cifcurve()` with `outcome.type="competing-risk"` by
calculating Aalen–Johansen estimator stratified by `fruitq1`. The time
points at which the macrovascular complications occurred were obtained
as `output2` for each strata using a helper function
`extract_time_to_event()`. Then, `cifplot()` is used to generate the
figure. These marks help distinguish between events due to the primary
cause and those attributable to competing causes. Note that the names of
`competing.risk.time` and `intercurrent.event.time` must match the
strata labels used in the plot if supplied by the user. The `label.y`,
`label.x` and `limit.x` arguments are also used to customize **the axis
labels and limits**.

``` r
output1 <- cifcurve(Event(t,epsilon)~fruitq1, data=diabetes.complications, 
                    outcome.type="competing-risk")
output2 <- extract_time_to_event(Event(t,epsilon)~fruitq1, 
                                 data=diabetes.complications, which.event="event2")
cifplot(output1, add.conf=FALSE, add.risktable=FALSE, 
        add.censor.mark=FALSE, add.competing.risk.mark=TRUE, competing.risk.time=output2, 
        label.y="CIF of diabetic retinopathy", label.x="Years from registration",
        limits.x=c(0,8))
```

<div class="figure">

<img src="man/figures/README-example04-1-2-1.png" alt="Cumulative incidence curves with competing risk marks" width="70%" />
<p class="caption">

Cumulative incidence curves with competing risk marks
</p>

</div>

The `label.strata` is another argument for customizing labels, but when
inputting a survfit object, it becomes invalid because it does not
contain stratum information. Therefore, the following code inputs the
formula and data. `label.strata` is used by combining `level.strata` and
`order.strata`. The `level.strata` specifies **the levels of the
stratification variable corresponding to each label** in `label.strata`.
The levels specified in `level.strata` are then displayed in the figure
**in the order defined by `order.strata`**. A figure enclosed in a
square was generated, which is due to `style="framed"` specification.

``` r
cifplot(Event(t,epsilon)~fruitq1, data=diabetes.complications, 
        outcome.type="competing-risk", add.conf=FALSE, add.risktable=FALSE, 
        add.estimate.table=TRUE, add.censor.mark=FALSE, add.competing.risk.mark=TRUE, 
        competing.risk.time=output2, label.y="CIF of diabetic retinopathy", 
        label.x="Years from registration", limits.x=c(0,8),
        label.strata=c("High intake","Low intake"), level.strata=c("Q2 to Q4","Q1"), 
        order.strata=c("Q1", "Q2 to Q4"), style="framed")
```

<div class="figure">

<img src="man/figures/README-example04-1-3-1.png" alt="Cumulative incidence curves with strata labels and framed style" width="70%" />
<p class="caption">

Cumulative incidence curves with strata labels and framed style
</p>

</div>

By specifying `add.estimate.table = TRUE`, **the risks of diabetic
retinopathy (estimates for CIFs) along with their CIs** are shown in the
table at the bottom of the figure. The risk ratios at a specific time
point (e.g. 8 years) for competing events can be jointly and coherently
estimated using `polyreg()` with `outcome.type = "competing-risk"`. In
the code of `polyreg()` below, no covariates are included in the
nuisance model (`~1` specifies intercept only). The effect of low fruit
intake `fruitq1` is estimated as **an unadjusted risk ratio**
(`effect.measure1="RR"`) for diabetic retinopathy (event 1) and
macrovascular complications (event 2) at 8 years (`time.point=8`).

``` r
output3 <- polyreg(nuisance.model=Event(t,epsilon)~1, exposure="fruitq1", 
          data=diabetes.complications, effect.measure1="RR", effect.measure2="RR", 
          time.point=8, outcome.type="competing-risk")
coef(output3)
#> [1] -1.38313159 -0.30043942 -3.99147406 -0.07582595
vcov(output3)
#>              [,1]         [,2]         [,3]         [,4]
#> [1,]  0.017018160 -0.012351309  0.009609321 -0.008372500
#> [2,] -0.012351309  0.012789187 -0.006012254  0.006540183
#> [3,]  0.009609321 -0.006012254  0.048161715 -0.044070501
#> [4,] -0.008372500  0.006540183 -0.044070501  0.055992232
summary(output3)
#> 
#>                       event1        event2      
#> ---------------------------------------------- 
#> fruitq1, Q2 to Q4 vs Q1 
#>                       0.740         0.927       
#>                       [0.593, 0.924]  [0.583, 1.474]
#>                       (p=0.008)     (p=0.749)   
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

The `summary()` method prints an event-wise table of point estimates,
CIs, and p-values. Internally, a `"polyreg"` object also supports the
**generics** API:

- `tidy()`: coefficient-level summaries (one row per term and per
  event),
- `glance()`: model-level summaries (follow-up, convergence, number of
  events),
- `augment()`: observation-level diagnostics (weights for IPCW,
  predicted CIFs, and influence-functions).

This means that `polyreg()` fits integrate naturally with the broader
`broom/modelsummary` ecosystem. For publication-ready tables, you can
pass `polyreg` objects directly to `modelsummary::msummary()` and
`modelsummary::modelplot()`, including exponentiated summaries (risk
ratios, odds ratios, subdistribution hazard ratios) via the
`exponentiate = TRUE` option.

``` r
output4 <- polyreg(nuisance.model=Event(t,epsilon)~1, exposure="fruitq1", 
                   data=diabetes.complications, effect.measure1="RR", effect.measure2="RR", 
                   time.point=2, outcome.type="competing-risk")
output5 <- polyreg(nuisance.model=Event(t,epsilon)~1, exposure="fruitq1", 
                   data=diabetes.complications, effect.measure1="RR", effect.measure2="RR", 
                   time.point=4, outcome.type="competing-risk")
output6 <- polyreg(nuisance.model=Event(t,epsilon)~1, exposure="fruitq1", 
                   data=diabetes.complications, effect.measure1="RR", effect.measure2="RR", 
                   time.point=6, outcome.type="competing-risk")
summary <- list(
  "RR of diabetes retinopathy at 2 years" = output4$summary$event1, 
  "RR of diabetes retinopathy at 4 years" = output5$summary$event1, 
  "RR of diabetes retinopathy at 6 years" = output6$summary$event1, 
  "RR of diabetes retinopathy at 8 years" = output3$summary$event1,
  "RR of marcovascular complications at 2 years" = output4$summary$event2, 
  "RR of marcovascular complications at 4 years" = output5$summary$event2, 
  "RR of marcovascular complications at 6 years" = output6$summary$event2, 
  "RR of marcovascular complications at 8 years" = output3$summary$event2
)
library(modelsummary)
modelplot(summary, coef_rename="", exponentiate = TRUE)
```

<div class="figure">

<img src="man/figures/README-example04-1-5-1.png" alt="Unadjusted risk ratios at 2, 4, 6 and 8 years by polyreg()" width="70%" />
<p class="caption">

Unadjusted risk ratios at 2, 4, 6 and 8 years by polyreg()
</p>

</div>
