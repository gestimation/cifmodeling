# Extract per-stratum event times from a formula and data

Creates a list of event times that can be passed to downstream
visualization or analysis functions such as `competing.risk.time` or
`intercurrent.event.time` in
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
and
[`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md).
Event types are specified by event 1, event 2, censoring, or
user-specified codes.

## Usage

``` r
extract_time_to_event(
  formula,
  data,
  subset.condition = NULL,
  na.action = na.omit,
  which.event = c("event2", "event1", "censor", "censoring", "user_specified"),
  code.event1 = 1,
  code.event2 = 2,
  code.censoring = 0,
  code.user.specified = NULL,
  read.unique.time = TRUE,
  drop.empty = TRUE
)
```

## Arguments

- formula:

  A model formula specifying the outcome and (optionally) `strata()`.

- data:

  A data frame containing variables in `formula`.

- subset.condition:

  Optional expression (as a character string) defining a subset of
  `data` to analyse. Defaults to `NULL`.

- na.action:

  Function to handle missing values (default: `na.omit` in stats).

- which.event:

  One of `"event1"`, `"event2"`, `"censor"`, `"censoring"`, or
  `"user_specified"`, indicating which event type to extract times for.

- code.event1, code.event2, code.censoring:

  Integer codes representing the event and censoring categories.
  Defaults are `1`, `2`, and `0`, respectively.

- code.user.specified:

  When `which.event = "user_specified"`, the integer event code to
  extract (e.g., 3 for an intercurrent event).

- read.unique.time:

  Logical if `TRUE`, only unique and sorted time points are returned for
  each stratum.

- drop.empty:

  Logical if `TRUE` (default), strata with no events are dropped from
  the returned list. Set to `FALSE` to retain empty strata as
  `numeric(0)` vectors (useful for diagnostics or consistent list
  length).

## Value

A named list of numeric vectors, where each element corresponds to a
stratum and contains the event times of the selected type.

## Details

This function is typically used internally by plotting and model
functions, but can also be called directly to inspect the per-stratum
event-time structure of a data frame.

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
data(diabetes.complications)
output <- extract_time_to_event(Event(t,epsilon) ~ fruitq,
                                data = diabetes.complications,
                                which.event = "event2")
cifplot(Event(t,epsilon) ~ fruitq,
        data = diabetes.complications,
        outcome.type="competing-risk",
        add.conf=FALSE,
        add.risktable=FALSE,
        add.censor.mark=FALSE,
        add.competing.risk.mark=TRUE,
        competing.risk.time=output,
        label.y="CIF of diabetic retinopathy",
        label.x="Years from registration")
```
