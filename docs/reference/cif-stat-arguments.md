# These arguments are shared by [`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md), [`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md), and [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md).

These arguments are shared by
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md),
[`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md),
and
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md).

## Arguments

- data:

  A data frame containing variables in the formula.

- weights:

  Optional name of the weight variable in `data`. Weights must be
  nonnegative.

- n.risk.type:

  Character string; one of `"weighted"`, `"unweighted"`, or `"ess"`.
  Controls which risk set size is returned in `$n.risk` without
  affecting estimates or SEs (default `"weighted"`).

- subset.condition:

  Optional character string giving a logical condition to subset `data`
  (default `NULL`).

- na.action:

  A function specifying the action to take on missing values (default
  `na.omit`).

- outcome.type:

  Character string specifying the type of time-to-event outcome. One of
  `"survival"` (Kaplan-Meier) or `"competing-risk"` (Aalen-Johansen). If
  `NULL` (default), the function automatically infers the outcome type
  from the data: if the event variable has more than two unique levels,
  `"competing-risk"` is assumed; otherwise, `"survival"` is used. You
  can also use abbreviations such as `"S"` or `"C"`. Mixed or ambiguous
  inputs (e.g., `c("S", "C")`) trigger automatic detection based on the
  event coding.

- code.event1:

  Integer code of the event of interest (default `1`).

- code.event2:

  Integer code of the competing risk (default `2`).

- code.censoring:

  Integer code of censoring (default `0`).

- error:

  Character string specifying the method for SEs and CIs used
  internally. For `"survival"` without weights, choose one of
  `"greenwood"` (default), `"tsiatis"`, or `"if"`. For
  `"competing-risk"` without weights, choose one of `"delta"` (default),
  `"aalen"`, or `"if"`. SEs and CIs based on influence functions
  (`"if"`) is recommended for weighted analysis.

- conf.type:

  Character specifying the method of transformation for CIs used
  internally (default `arcsine-square root`).

- conf.int:

  Numeric two-sided level of CIs (default `0.95`).
