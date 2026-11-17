# Tests of cumulative incidence functions

[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
compares survival or cumulative incidence curves between groups using a
formula interface. For competing-risk outcomes, it provides the
classical Gray test, score and augmented score tests, and tests
combining multiple Fleming–Harrington (FH) temporal-weight directions.
Use the Ditzhaus test for a quadratic statistic or MaxCombo for a
maximum standardized score.

This vignette focuses on **test selection, function calls, and
interpreting
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
output**. For a complete analysis with CIF estimation and risk ratios,
see [Worked
examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md).

## Run a single CIF comparison with `ciftest()`

In `diabetes.complications`, `epsilon` denotes diabetic retinopathy (1),
a competing macrovascular event (2), or censoring (0). The formula’s
right-hand side identifies the groups. The variable `fruitq1` compares
the lowest fruit-intake quartile with the other quartiles. The event of
interest defaults to `code.event1 = 1`.

``` r

library(cifmodeling)
data(diabetes.complications)
```

Request the classical Gray test with `test = "gray"`:

``` r

fit_gray <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "gray"
)
fit_gray
```

    ## 
    ## Gray's test
    ## 
    ## Outcome: competing-risk
    ## Test: gray
    ## FH weights: rho = 0, gamma = 0
    ## Analysis strata: pooled
    ## Censoring strata: not used by this test
    ## Working AJ strata: not used by this test
    ## Variance: gray
    ## 
    ## Chi-squared = 8.8786, df = 1, p-value = 0.002885

``` r

generics::tidy(fit_gray)
```

    ##   statistic     p.value parameter      method   outcome.type test     weight
    ## 1  8.878578 0.002885361         1 Gray's test competing-risk gray unweighted
    ##   iteration fixed.point.residual   n
    ## 1         0                   NA 978

To fit the augmented score test, change `test` to `"augmented"`. This
procedure uses working models for event and censoring distributions when
constructing the score and its variance:

``` r

fit_aug <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "augmented"
)
fit_aug
```

    ## 
    ## Closed-form augmented Fine-Gray score test
    ## 
    ## Outcome: competing-risk
    ## Test: augmented
    ## FH weights: rho = 0, gamma = 0
    ## Analysis strata: pooled
    ## Censoring strata: pooled
    ## Working AJ strata: fruitq1 only (2 fitted cells)
    ## Iteration: 0 (closed-form augmentation)
    ## Variance: score-iid
    ## 
    ## Chi-squared = 7.931, df = 1, p-value = 0.004859

``` r

generics::tidy(fit_aug)
```

    ##   statistic     p.value parameter                                     method
    ## 1   7.93102 0.004859441         1 Closed-form augmented Fine-Gray score test
    ##     outcome.type      test     weight iteration fixed.point.residual   n
    ## 1 competing-risk augmented unweighted         0                   NA 978

Both fitted objects contain `statistic`, `parameter` (degrees of
freedom), `p.value`, and `method`. Use
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) to
collect these fields in a data frame. The p-value concerns differences
between CIFs over follow-up, **not** a risk ratio at a fixed time.
Because Gray and augmented tests use different scores and covariance
estimators, their p-values may differ; neither test is uniformly more
powerful.

By default, `test = "auto"` selects the augmented score test for
competing risks and the log-rank test for survival outcomes. Use
`test = "gray"` when you specifically want the classical Gray test.
These tests compare groups but do not, by themselves, establish causal
effects.

## Select early or late weighting

To emphasize differences early or late in follow-up,
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
accepts FH temporal weights of the form

``` math
w(t)=\{1-F_1(t-)\}^{\rho}F_1(t-)^{\gamma},
```

Here, $`F_1(t-)`$ is the left-limit CIF used to construct the weight.
Larger $`\rho`$ emphasizes earlier differences, while larger $`\gamma`$
emphasizes later differences. The **estimated weight process depends on
the component test**: Gray and augmented score tests need not use the
same estimated CIF. Inspect the process used in a fit through
`fit$diagnostics$fh.weight.process`.

For competing-risk outcomes, `test = "early"` and `test = "late"` select
augmented score tests with FH(1,0) and FH(0,1), respectively:

``` r

fit_early <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "early"
)
fit_late <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "late"
)
fit_early
```

    ## 
    ## Closed-form augmented Fine-Gray score test
    ## 
    ## Outcome: competing-risk
    ## Test: augmented
    ## FH weights: rho = 1, gamma = 0
    ## Analysis strata: pooled
    ## Censoring strata: pooled
    ## Working AJ strata: fruitq1 only (2 fitted cells)
    ## Iteration: 0 (closed-form augmentation)
    ## Variance: score-iid
    ## 
    ## Chi-squared = 8.0774, df = 1, p-value = 0.004482

``` r

fit_late
```

    ## 
    ## Closed-form augmented Fine-Gray score test
    ## 
    ## Outcome: competing-risk
    ## Test: augmented
    ## FH weights: rho = 0, gamma = 1
    ## Analysis strata: pooled
    ## Censoring strata: pooled
    ## Working AJ strata: fruitq1 only (2 fitted cells)
    ## Iteration: 0 (closed-form augmentation)
    ## Variance: score-iid
    ## 
    ## Chi-squared = 4.9028, df = 1, p-value = 0.02681

To use a different direction, specify `rho` and `gamma` alongside the
test, for example, `test = "augmented", rho = 2, gamma = 0`. FH(0,0) is
called `"unweighted"` because no additional temporal weight is applied;
this does **not** mean that censoring adjustment is omitted. The terms
`"early"` and `"late"` describe *when* differences are emphasized, not
*which group* has the higher CIF.

## Combine directions with the Ditzhaus test

`test = "ditzhaus"` (or `"d"`) runs a quadratic-form test. For
competing-risk outcomes, the default combines augmented scores with
FH(0,0), FH(0,1), and FH(2,0) weights to give one omnibus p-value.

``` r

fit_ditzhaus <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "ditzhaus"
)
fit_ditzhaus
```

    ## 
    ## Multiple-direction Closed-form augmented Fine-Gray score test (asymptotic chi-square)
    ## 
    ## Outcome: competing-risk
    ## Component test: augmented
    ## Directions: unweighted (0, 0), FH01 (0, 1), early (2, 0)
    ## 
    ## Chi-squared = 8.5265, df = 3, p-value = 0.0363

``` r

generics::tidy(fit_ditzhaus)
```

    ##   statistic   p.value parameter
    ## 1  8.526544 0.0362953         3
    ##                                                                                  method
    ## 1 Multiple-direction Closed-form augmented Fine-Gray score test (asymptotic chi-square)
    ##     outcome.type      test directions   n
    ## 1 competing-risk augmented          3 978

``` r

fit_ditzhaus$directions
```

    ##    direction rho gamma
    ## 1 unweighted   0     0
    ## 2       FH01   0     1
    ## 3      early   2     0

[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
combines the scores into a quadratic-form statistic:

``` math
Q=U^\top V^+U.
```

The fitted object stores the score vector $`U`$ in `fit_ditzhaus$score`
and its covariance matrix $`V`$ in `fit_ditzhaus$vcov.score`. The test
uses the generalized inverse $`V^+`$ and a chi-squared reference
distribution with degrees of freedom equal to the **numerical rank** of
$`V`$. Consequently, three directions need not yield three degrees of
freedom. `fit_ditzhaus$directions` records the FH weights, and
`fit_ditzhaus$components` contains the individual fits. This omnibus
statistic does not estimate a risk difference or ratio.

The `"early"` and `"late"` presets within `directions` use **FH(2,0)**
and **FH(0,2)**. These differ from the single-direction shortcuts
`test = "early"` and `test = "late"`. For custom weights, use
`directions` and, if necessary, `component.test`:

``` r

fit_custom <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  component.test = "augmented",
  directions = rbind(
    early = c(rho = 1, gamma = 0),
    late  = c(rho = 0, gamma = 1)
  )
)
fit_custom
```

    ## 
    ## Multiple-direction Closed-form augmented Fine-Gray score test (asymptotic chi-square)
    ## 
    ## Outcome: competing-risk
    ## Component test: augmented
    ## Directions: early (1, 0), late (0, 1)
    ## 
    ## Chi-squared = 8.0954, df = 2, p-value = 0.01746

Supplying `directions` alone selects the Ditzhaus test. For competing
risks, the default component test is `"augmented"`; set
`component.test = "gray"` or `"score"` to change the score and
covariance calculation. For survival outcomes, you can use
`component.test = "logrank"`.

## Run MaxCombo with `test = "maxcombo"`

For two groups with competing risks, `test = "maxcombo"` (or `"m"`)
combines augmented scores using FH(0,0), FH(2,0), and FH(0,1) by
default. Supply `directions` to change the weights. The default
two-sided test uses the largest absolute standardized score and a joint
normal reference.

``` r

fit_maxcombo <- ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "m",
  code.exposure.ref = 0
)
fit_maxcombo
```

    ## 
    ## MaxCombo augmented score test (joint normal reference)
    ## 
    ## Outcome: competing-risk
    ## Component test: augmented
    ## Directions: unweighted (0, 0), early (2, 0), FH01 (0, 1)
    ## 
    ## Comparison: 1 vs 0, alternative = two.sided
    ## Standardized scores: unweighted = 2.8162, early = 2.8371, FH01 = 2.2142
    ## Max Z = 2.8371, p-value = 0.008381

For a prespecified lower CIF in the comparison group, set
`alternative = "less"`. Use `"greater"` for a higher CIF.
`code.exposure.ref` selects the reference group by its code or factor
level. The default is the first factor level present in the analysis
data. The returned comparison identifies the two groups.

``` r

ciftest(
  Event(t, epsilon) ~ fruitq1,
  data = diabetes.complications,
  outcome.type = "competing-risk",
  test = "maxcombo",
  alternative = "less",
  code.exposure.ref = 0
)
```

    ## 
    ## MaxCombo augmented score test (joint normal reference)
    ## 
    ## Outcome: competing-risk
    ## Component test: augmented
    ## Directions: unweighted (0, 0), early (2, 0), FH01 (0, 1)
    ## 
    ## Comparison: 1 vs 0, alternative = less
    ## Standardized scores: unweighted = 2.8162, early = 2.8371, FH01 = 2.2142
    ## Max Z = -2.2142, p-value = 0.9955

## Choosing a test

Use this table to select a
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)
option:

| Goal | Call or result | What you obtain |
|:---|:---|:---|
| Standard survival-curve comparison | `ciftest(..., outcome.type = "survival", test = "auto")` | Log-rank test |
| Classical CIF comparison | `ciftest(..., test = "gray")` | Gray test statistic and p-value |
| Augmented CIF comparison (default) | `ciftest(..., test = "augmented")` | Augmented score test and p-value |
| Unaugmented score test | `ciftest(..., test = "score")` | Score test using the null-score IID variance |
| Focus on an earlier or later difference | `ciftest(..., test = "early")` or `"late"` | Single-direction FH-weighted test |
| Combine directions in a quadratic form | `ciftest(..., test = "ditzhaus")` or `directions = ...` | Omnibus quadratic statistic, rank-based df, p-value |
| Use the largest signed or absolute standardized score | `ciftest(..., test = "maxcombo", alternative = ...)` | MaxCombo statistic and global p-value |

These options test related but distinct alternatives. For confirmatory
analyses, prespecify the test and weight directions; choosing the
smallest p-value after inspecting the data does not preserve the nominal
significance level. The built-in Ditzhaus test uses an asymptotic
chi-squared reference. Ditzhaus and Friedrich (2020) studied permutation
calibration for weighted log-rank combinations. See [Worked
examples](https://gestimation.github.io/cifmodeling/articles/v5_examples.md)
for a simpler comparison using this dataset and [Function
reference](https://gestimation.github.io/cifmodeling/articles/v9_arguments.md)
for the full argument list.

## References

- Gray RJ (1988). A class of K-sample tests for comparing the cumulative
  incidence of a competing risk. *The Annals of Statistics*, **16**,
  1141–1154.
- Ditzhaus M, Friedrich S (2020). More powerful logrank permutation
  tests for two-sample survival data. *Journal of Statistical
  Computation and Simulation*, **90**, 2209–2227.
  [doi:10.1080/00949655.2020.1773463](https://doi.org/10.1080/00949655.2020.1773463).
