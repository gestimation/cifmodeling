# Create a survival or competing-risks response

A lightweight response constructor used in
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
and
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
to pass survival and competing-risks data via a model formula.

## Usage

``` r
Event(time, event, allowed = getOption("cifmodeling.allowed", c(0, 1, 2)))
```

## Arguments

- time:

  Numeric vector of follow-up times (non-negative).

- event:

  Integer (0=censor, 1,2,...) or a character/factor vector whose levels
  are numeric codes "0","1","2",... for competing events.

- allowed:

  Numeric vector of acceptable event codes.

## Value

An object of class `"Event"` (a 2-column matrix) with columns `time`,
`event`.

## Lifecycle

**\[stable\]**

## See also

[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
for log-odds product modeling of CIFs;
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
## event: 0=censor, 1=primary, 2=competing
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
```
