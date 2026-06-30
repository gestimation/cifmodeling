# Generate a survival/CIF curve with marks that represent censoring, competing risks and intercurrent events

This function generates a survival or CIF curve from a unified
formula–data interface or from an existing survfit object. When a
formula is supplied, the LHS is typically
[`Event()`](https://gestimation.github.io/cifmodeling/reference/Event.md)
or [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html), and
the RHS specifies an optional stratification variable. In addition to
the curves themselves, `cifplot()` can add numbers-at-risk tables,
tables of point estimates and CIs, censoring marks, competing-risk
marks, and intercurrent-event marks.

For usual single-panel mode, the function returns an object whose `plot`
component is a regular ggplot object that can be further modified
(compatible with `+` and `%+%`). For more complex multi-panel displays,
`cifplot()` can internally call
[`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
via several “panel modes” (per event, per variable, or
censoring-focused).

## Usage

``` r
cifplot(
  formula_or_fit,
  data = NULL,
  weights = NULL,
  subset.condition = NULL,
  na.action = na.omit,
  outcome.type = c("competing-risk", "survival"),
  code.event1 = 1,
  code.event2 = 2,
  code.censoring = 0,
  code.events = NULL,
  error = NULL,
  conf.type = "arcsine-square root",
  conf.int = 0.95,
  n.risk.type = c("weighted", "unweighted", "ess"),
  type.y = NULL,
  label.x = "Time",
  label.y = NULL,
  label.strata = NULL,
  level.strata = NULL,
  order.strata = NULL,
  limits.x = NULL,
  limits.y = NULL,
  breaks.x = NULL,
  breaks.y = NULL,
  use.coord.cartesian = FALSE,
  add.conf = TRUE,
  add.risktable = TRUE,
  add.estimate.table = FALSE,
  symbol.risk.table = "square",
  font.size.risk.table = 3,
  add.censor.mark = TRUE,
  shape.censor.mark = 3,
  size.censor.mark = 2,
  add.competing.risk.mark = FALSE,
  competing.risk.time = list(),
  shape.competing.risk.mark = 16,
  size.competing.risk.mark = 2,
  add.intercurrent.event.mark = FALSE,
  intercurrent.event.time = list(),
  shape.intercurrent.event.mark = 1,
  size.intercurrent.event.mark = 2,
  add.quantile = FALSE,
  level.quantile = 0.5,
  panel.per.event = FALSE,
  panel.censoring = FALSE,
  panel.per.variable = FALSE,
  panel.mode = "auto",
  rows.columns.panel = NULL,
  style = "classic",
  palette = NULL,
  linewidth = 0.8,
  linetype = FALSE,
  font.family = "sans",
  font.size = 12,
  legend.position = "top",
  print.panel = FALSE,
  filename.ggsave = NULL,
  width.ggsave = 6,
  height.ggsave = 6,
  dpi.ggsave = 300,
  survfit.info = NULL,
  axis.info = NULL,
  visual.info = NULL,
  panel.info = NULL,
  style.info = NULL,
  inset.info = NULL,
  print.info = NULL,
  ggsave.info = NULL,
  ...
)
```

## Arguments

- formula_or_fit:

  Either a model formula or a survfit object. When a formula is
  supplied, the LHS must be `Event(time, status)` or
  `Surv(time, status)`. The RHS specifies an optional stratification
  variable.

- data:

  A data frame containing variables in the formula.

- weights:

  Optional name of the weight variable in `data`. Weights must be
  nonnegative.

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

- code.events:

  Optional specification of event/censoring codes. For single-panel
  calls, supply a numeric vector. For competing-risk outcomes, use
  `c(event1, event2, censoring)`. For survival outcomes, a length-2 or
  length-3 vector is allowed: `c(event, censoring)` or
  `c(event, *, censoring)`, where any middle element is ignored. When
  supplied, this argument overrides `code.event1`, `code.event2`, and
  `code.censoring` for the purpose of estimation. For panel displays
  (e.g.
  [`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
  or when `panel.per.event = TRUE` or `panel.censoring = TRUE`),
  `code.events` may also be a list of such numeric vectors, one per
  panel.

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

- n.risk.type:

  Character string; one of `"weighted"`, `"unweighted"`, or `"ess"`.
  Controls which risk set size is returned in `$n.risk` without
  affecting estimates or SEs (default `"weighted"`).

- type.y:

  Character string specifying the y-scale. For survival/CIF curves,
  `"surv"` implies survival probabilities and `"risk"` implies CIF
  (1-survival in simple survival settings). Specify `"cumhaz"` to plot
  cumulative hazard or `"cloglog"` to generate a complementary log-log
  plot. If `NULL`, a default is chosen from `outcome.type` or the
  survfit object.

- label.x:

  Character x-axis label (default `"Time"`).

- label.y:

  Character y-axis label (default is chosen automatically from
  `outcome.type` and `type.y`, e.g. "Survival", "Cumulative incidence"
  or "Cumulative hazard").

- label.strata:

  Character vector or named character vector specifying labels for
  strata. Names (if present) must match the (re-ordered) underlying
  strata levels. **Note:** when any of the panel modes is active
  (`panel.per.variable = TRUE`, `panel.per.event = TRUE`,
  `panel.censoring = TRUE`, or `panel.mode = "auto"` and it actually
  dispatches to a panel), strata labels are suppressed to avoid
  duplicated legends across sub-plots.

- level.strata:

  Optional character vector giving the full set of expected strata
  levels. When provided, both `order.strata` and `label.strata` are
  validated against it before application.

- order.strata:

  Optional character vector specifying the display order of strata in
  the legend/number-at-risk table. Specify the levels of strata. Levels
  not listed are dropped.

- limits.x:

  Numeric length-2 vector specifying x-axis limits. If `NULL`, it is set
  from the fitted object (typically `c(0, max(time))`).

- limits.y:

  Numeric length-2 vector specifying y-axis limits. If `NULL`, it is set
  to `c(0, 1)` for probability-type outcomes.

- breaks.x:

  Numeric vector of x-axis breaks (default `NULL`).

- breaks.y:

  Numeric vector of y-axis breaks (default `NULL`).

- use.coord.cartesian:

  Logical; if `TRUE`, uses
  [`ggplot2::coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
  for zooming instead of changing the scale limits (default `FALSE`).

- add.conf:

  Logical; if `TRUE`, adds a CI ribbon (via
  [`ggsurvfit::add_confidence_interval()`](http://www.danieldsjoberg.com/ggsurvfit/reference/add_confidence_interval.md)).
  Default `TRUE`.

- add.risktable:

  Logical; if `TRUE`, adds a numbers-at-risk table under the plot.
  Default `TRUE`. **Note:** when a panel mode is active, tables are
  suppressed.

- add.estimate.table:

  Logical; if `TRUE`, adds a table of estimates and CIs. Default
  `FALSE`. **Note:** when a panel mode is active, tables are suppressed.

- symbol.risk.table:

  Character specifying the symbol used in the risk table to denote
  strata: `"square"`, `"circle"`, or `"triangle"` (default `"square"`).

- font.size.risk.table:

  Numeric font size for texts in risk / estimate tables (default `3`).

- add.censor.mark:

  Logical; if `TRUE`, draws censoring marks on each curve (via
  [`ggsurvfit::add_censor_mark()`](http://www.danieldsjoberg.com/ggsurvfit/reference/add_censor_mark.md)).
  Default `TRUE`.

- shape.censor.mark:

  Integer point shape used for censoring marks (default `3`).

- size.censor.mark:

  Numeric point size used for censoring marks (default `2`).

- add.competing.risk.mark:

  Logical; if `TRUE`, draws time marks for the competing event (event
  2). If no times are supplied via `competing.risk.time`, the function
  tries to extract them automatically from the data. Default `FALSE`.

- competing.risk.time:

  A **named list** of numeric vectors. Each name must correspond to a
  strata label, and its numeric vector gives the times at which the
  competing event occurred in that stratum. Typically left as
  [`list()`](https://rdrr.io/r/base/list.html) and filled internally.

- shape.competing.risk.mark:

  Integer point shape for competing-risk marks (default `16`).

- size.competing.risk.mark:

  Numeric point size for competing-risk marks (default `2`).

- add.intercurrent.event.mark:

  Logical; if `TRUE`, overlays user-specified intercurrent-event times
  per stratum. Default `FALSE`.

- intercurrent.event.time:

  A **named list** of numeric vectors for intercurrent events (names
  must match strata labels).

- shape.intercurrent.event.mark:

  Integer point shape for intercurrent-event marks (default `1`).

- size.intercurrent.event.mark:

  Numeric point size for intercurrent-event marks (default `2`).

- add.quantile:

  Logical; if `TRUE`, adds a quantile reference line (via
  [`ggsurvfit::add_quantile()`](http://www.danieldsjoberg.com/ggsurvfit/reference/add_quantile.md)).
  Default `FALSE`.

- level.quantile:

  Numeric quantile level to be shown (default `0.5` for the median).

- panel.per.event:

  Logical. **Explicit panel mode.** If `TRUE` and
  `outcome.type == "competing-risk"`, `cifplot()` internally calls
  [`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
  to display two event-specific CIFs side-by-side (event 1 and event 2)
  using reversed `code.events`. Ignored for non-competing-risk outcomes.

- panel.censoring:

  Logical. **Explicit panel mode.** If `TRUE` and
  `outcome.type == "survival"`, `cifplot()` internally calls
  [`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
  to display KM curves for `(event, censor)` and `(censor, event)` so
  that censoring patterns can be inspected.

- panel.per.variable:

  Logical. **Explicit panel mode.** If `TRUE` and the RHS of the formula
  has multiple covariates (e.g. `~ a + b + c`), the function produces a
  panel where each variable in the RHS is used once as the
  stratification factor.

- panel.mode:

  Character specifying **Automatic panel mode.** If `"auto"` and none of
  `panel.per.variable`, `panel.per.event`, `panel.censoring` has been
  set to `TRUE`, the function chooses a suitable panel mode
  automatically: (i) if the formula RHS has 2+ variables, it behaves
  like `panel.per.variable = TRUE`; (ii) otherwise, if
  `outcome.type == "competing-risk"`, it behaves like
  `panel.per.event = TRUE`; (iii) otherwise, if
  `outcome.type == "survival"`, it behaves like
  `panel.censoring = TRUE`. If a panel mode is explicitly specified,
  `panel.mode` is ignored.

- rows.columns.panel:

  Optional integer vector `c(nrow, ncol)` controlling the layout of the
  panel returned by the panel modes. If `NULL`, an automatic layout is
  determined from the number of subplots.

- style:

  Character choosing the base plot style: `"classic"`, `"bold"`,
  `"framed"`, `"grid"`, `"gray"` or `"ggsurvfit"` (default `"classic"`).
  Abbreviations such as `"C"`, `"B"`, `"F"`, or `"G"` are also accepted.

- palette:

  Optional character vector specifying the color palette to use across
  strata.

- linewidth:

  Optional numeric specifying the line width of curve (default `0.8`).

- linetype:

  Optional logical using different line types of curve (default
  `FALSE`).

- font.family:

  Character specifying the font family: `"sans"`, `"serif"`, or `"mono"`
  (default `"sans"`).

- font.size:

  Integer specifying the base font size (default `12`).

- legend.position:

  Character specifying the legend position: `"top"`, `"right"`,
  `"bottom"`, `"left"`, or `"none"` (default `"top"`).

- print.panel:

  Logical. When `TRUE`, panel displays created internally are printed
  automatically in interactive sessions; otherwise they are returned
  invisibly for further modification (default `FALSE`).

- filename.ggsave:

  Character; if non-`NULL`, save the plot to this file.

- width.ggsave:

  Numeric width passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  (default `6`).

- height.ggsave:

  Numeric height passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  (default `6`).

- dpi.ggsave:

  Numeric DPI passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  (default `300`).

- survfit.info, axis.info, visual.info, panel.info, style.info,
  inset.info, print.info, ggsave.info:

  Internal lists used for programmatic control. Not intended for direct
  user input.

- ...:

  Additional arguments passed to internal helper functions.

## Value

A `"cifplot"` object (a list) with at least the following elements:

- `plot`: a ggplot object containing the main plot

- `patchwork`: reserved for compatibility with panel displays (typically
  `NULL` for single-panel plots)

- `survfit.info`, `axis.info`, `visual.info`, `panel.info`,
  `style.info`, `inset.info`, `print.info`, `ggsave.info`: internal
  lists storing the fitted curves and display settings

- `version`: a character string giving the cifmodeling version used

- `call`: the original function call

When a panel mode is active and `print.panel = TRUE`, the panel is also
printed in interactive sessions.

## Details

### Typical use cases

- Draw one survival/CIF curve set by exposure groups (e.g., treatment vs
  control).

- Call
  [`cifpanel()`](https://gestimation.github.io/cifmodeling/reference/cifpanel.md)
  with a simplified code to create a panel displaying plots of multiple
  stratified survival/CIF curves or CIF curves for each event type.

- Add CIs and censor/competing-risk/intercurrent-event marks.

- Add number-at-risk table to display the number at risk or the
  estimated survival probabilities or CIFs and CIs at each point in
  time.

### Key arguments shared with cifcurve()

- **Outcome type and estimator**

  - `outcome.type = "survival"`: Kaplan-Meier estimator

  - `outcome.type = "competing-risk"`: Aalen-Johansen estimator

- **Confidence intervals**

  - `conf.int` sets the two-sided level (default 0.95)

  - `conf.type` chooses the transformation (`"arcsine-square root"`,
    `"plain"`, `"log"`, `"log-log"`, `"logit"`, or `"none"`)

  - `error` chooses the estimator for SE (`"greenwood"`, `"tsiatis"` or
    `"if"` for survival curves and `"delta"`, `"aalen"` or `"if"` for
    CIFs)

- **Risk sets used in tables**

  - `n.risk.type` controls whether `$n.risk` reflects weighted,
    unweighted, or effective sample size counts when building risk
    tables (forwarded to
    [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md));
    when a fitted `survfit` object is supplied, existing risk sets are
    used as-is.

### Key arguments for cifplot()

- **Data visualization**

  - `add.conf` adds CIs on the ggplot2-based plot

  - `add.competing.risk.mark` and `add.intercurrent.event.mark` adds
    symbols to describe competing risks or intercurrent events in
    addition to conventional censoring marks with `add.censor.mark`

  - `add.risktable` adds numbers at risk

  - `add.estimate.table` adds time-by-time estimates and CIs

  - `add.quantile` adds a reference line at a chosen quantile level

- **Plot customization**

  - `type.y` chooses y-axis (`"surv"` for survival and `"risk"` for
    1-survival/CIF)

  - `limits.x`, `limits.y`, `breaks.x`, `breaks.y`: numeric vectors for
    axis control

  - `style` specifies the appearance of plot (`"classic"`, `"bold"`,
    `"framed"`, `"grid"`, `"gray"` or `"ggsurvfit"`)

  - `palette` specifies color of each curve (e.g.
    `palette=c("blue1", "cyan3", "navy", "deepskyblue3"))`)

- **Panel display**

  - `panel.per.variable` produces multiple survival/CIF curves per
    stratification variable specified in the formula

  - `panel.per.event` produces CIF curves for each event type

  - `panel.censoring` produces the Kaplan–Meier curves for (event,
    censor) and (censor, event) so that censoring patterns can be
    inspected

  - `panel.mode` uses automatic panel mode

When `panel.per.event = TRUE`, two panels are created with
`code.events = list(c(e1, e2, c), c(e2, e1, c))`, where
`code.events = c(e1, e2, c)` is the input coding for event1, event2, and
censoring. Common legend is collected by default
(`legend.collect = TRUE`).

Numeric stratification variables are normalized automatically. Columns
with fewer than nine distinct numeric values are coerced to factors;
columns with nine or more distinct numeric values are split at the
median into “Below median” and “Above median” strata.

### Advanced control not required for typical use

The arguments below fine-tune internal estimation and figure appearance.
**Most users do not need to change these defaults.**

#### Graphical layers

|  |  |  |
|----|----|----|
| Argument | Description | Default |
| `add.conf` | Add confidence interval ribbon. | `TRUE` |
| `add.risktable` | Add numbers-at-risk table below the plot. | `TRUE` |
| `add.estimate.table` | Add estimates and confidence intervals table. | `FALSE` |
| `symbol.risk.table` | Symbol for strata in risk / estimate tables | `"square"` |
| `font.size.risk.table` | Font size for texts in risk / estimate tables | `3` |
| `add.censor.mark` | Add censoring marks. | `TRUE` |
| `add.competing.risk.mark` | Add marks for event2 of "competing-risk" outcome. | `FALSE` |
| `add.intercurrent.event.mark` | Add intercurrent event marks at user-specified times. | `FALSE` |
| `add.quantile` | Add quantile reference lines. | `FALSE` |
| `level.quantile` | Quantile level for `add.quantile`. | `0.5` |

#### Time for marks

|  |  |
|----|----|
| Argument | Description |
| `competing.risk.time` | **Named list** of numeric vectors that contains times of competing risks. Names must match strata labels. Typically created internally |
| `intercurrent.event.time` | **Named list** of numeric vectors that contains times of intercurrent events. Names must match strata labels. Typically created by [`extract_time_to_event()`](https://gestimation.github.io/cifmodeling/reference/extract_time_to_event.md). |

#### Appearance of marks

|                                 |                      |                      |
|---------------------------------|----------------------|----------------------|
| Argument                        | Applies to           | Default              |
| `shape.censor.mark`             | Censoring marks      | `3` (cross)          |
| `size.censor.mark`              | Censoring marks      | `2`                  |
| `shape.competing.risk.mark`     | Competing-risk marks | `16` (filled circle) |
| `size.competing.risk.mark`      | Competing-risk marks | `2`                  |
| `shape.intercurrent.event.mark` | Intercurrent marks   | `1` (circle)         |
| `size.intercurrent.event.mark`  | Intercurrent marks   | `2`                  |

#### Panel display

|  |  |
|----|----|
| Argument | Description |
| `panel.per.variable` | One panel per stratification variable |
| `panel.per.event` | For `"competing-risk"`, show CIFs of event 1 and event 2 |
| `panel.censoring` | For survival, show (event, censor) vs (censor, event) |
| `panel.mode` with 2+ stratification variables | Behave like `panel.per.variable` |
| `panel.mode` with outcome.type = "competing-risk" | Behave like `panel.per.event` |
| `panel.mode` with outcome.type = "survival" | Behave like `panel.censoring` |

#### Axes and legend

|  |  |  |
|----|----|----|
| Argument | Description | Default |
| `limits.x`, `limits.y` | Axis limits (`c(min, max)`) | Auto |
| `breaks.x`, `breaks.y` | Tick breaks for x and y axes | Auto |
| `use.coord.cartesian` | For zooming use `coord_cartesian()` | `FALSE` |
| `legend.position` | `"top"`, `"right"`, `"bottom"`, `"left"`, `"none"` | `"top"` |

#### Export

|                   |                                               |         |
|-------------------|-----------------------------------------------|---------|
| Argument          | Description                                   | Default |
| `filename.ggsave` | If non-`NULL`, save the plot using `ggsave()` | `NULL`  |
| `width.ggsave`    | Size passed to `ggsave()`                     | `6`     |
| `height.ggsave`   | Size passed to `ggsave()`                     | `6`     |
| `dpi.ggsave`      | DPI passed to `ggsave()`                      | `300`   |

**Notes**

- For CIF displays, set `type.y = "risk"`. For survival scale, use
  `type.y = NULL` or `= "surv"`. For a cumulative hazard plot, use
  `type.y = "cumhaz"`. To generate a log-log plot, use
  `type.y = "cloglog"`.

- Event coding can be controlled via `code.event1`, `code.event2`,
  `code.censoring`. For ADaM-style data, use `code.event1 = 0`,
  `code.censoring = 1`.

- Per-stratum time lists should have names identical to plotted strata
  labels.

## Lifecycle

**\[stable\]**

## See also

[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md)
for log-odds product modeling of CIFs;
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
for KM/AJ estimators;
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
cifplot(Event(t,epsilon) ~ fruitq,
        data = diabetes.complications,
        outcome.type="competing-risk",
        add.risktable = FALSE,
        label.y='CIF of diabetic retinopathy',
        label.x='Years from registration')

```
