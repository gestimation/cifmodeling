# Create a survival and competing-risk report

Create a survival and competing-risk report

## Usage

``` r
cifreport(
  data,
  outcomes,
  exposure,
  characteristics = character(),
  time.point = NULL,
  subset.condition = NULL,
  filename.cifreport = NULL,
  title = NULL,
  overall = FALSE,
  conf.int = 0.95,
  table1.args = list(),
  code.exposure.ref = NULL,
  data.name = NULL,
  subgroup = NULL,
  label.strata = NULL,
  order.strata = NULL,
  level.strata = NULL,
  outcome.type = NULL,
  cifplot.args = list(),
  output.file = NULL
)
```

## Arguments

- data:

  A data frame with one row per subject.

- outcomes:

  A named list of outcome definitions. Each is a list with column names
  `time` and `event`, and optionally `outcome.type` to override the
  shared outcome type, `code.event1` (default 1), `code.event2` (default
  2, for competing risks), `code.censoring` (default 0), and
  `time.point` to override the shared time points. Use
  `estimator = "KM"` or `"1-KM"` for survival and `"AJ"` or `"1-AJ"` for
  competing risks. Defaults are `"KM"` and `"AJ"`. The name of each
  definition is its outcome label.

- exposure:

  Character vector of grouping column names, or a named list whose
  entries are column names or lists with `variable` and optional
  `code.exposure.ref`, `level.strata`, `order.strata`, and
  `label.strata`. Vector or list names identify exposure definitions.

- characteristics:

  Character vector of baseline columns for Table 1.

- time.point:

  Shared numeric evaluation times. May be omitted when each outcome
  supplies its own times.

- subset.condition:

  One condition defining the analysis population. Use a logical vector,
  character condition, one-sided formula, or quoted expression. Missing
  condition values are excluded. `NULL` includes all rows.

- filename.cifreport:

  One or more paths ending in `.docx`, `.pdf`, `.html`, or `.htm`.
  Relative paths are resolved from the working directory. Parent
  directories must exist. Existing files are overwritten after
  rendering. If omitted, saves
  `yyyymmdd_<data.name>_analysis_report.docx` in the working directory.
  Supply `NULL` to return the results without writing any files.
  Multiple filenames share one analysis. Word output uses editable
  tables with 9-point text and requires the `xml2` and `zip` packages.
  PDF requires Chrome, Chromium, or Microsoft Edge. Set the
  `CIFREPORT_BROWSER` environment variable to the browser executable if
  it is not found automatically.

- title:

  Report title. The default is
  `Statistical Analysis Report of <data.name> on <yyyymmdd>`.

- overall:

  Include an overall column in Table 1 when there are multiple groups.

- conf.int:

  Pointwise confidence level for the curve and risk ratio.

- table1.args:

  Named list of additional arguments to
  [`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html),
  such as `label` and `type`. Column `label` attributes are also
  supported. Characteristics with a named `labels` attribute are
  summarized as categorical using those value labels. When comparison is
  possible, Table 1 includes p-values from the default tests of
  [`gtsummary::add_p()`](https://www.danieldsjoberg.com/gtsummary/reference/add_p.html)
  across the observed exposure groups.

- code.exposure.ref:

  Reference level for the risk ratios. Supply one level for all
  exposures or a named vector/list keyed by exposure section labels. A
  reference in an exposure definition takes precedence. By default the
  first observed factor level is the reference.

- data.name:

  Dataset name used in the default filename and title. By default it
  uses the supplied data object's name, or `"data"` for an inline
  expression.

- subgroup:

  Name of one categorical column. In addition to the analysis population
  as a whole, the report includes each observed level of this column
  within that population. Missing values contribute only to the analysis
  population as a whole. Numeric category codes are supported.

- label.strata:

  Display labels, either named by the original exposure codes or given
  in the order of `level.strata`. Supply any strata argument as a
  vector, or as a named list keyed by exposure definition for separate
  settings. Exposure-specific settings in `exposure` take precedence.
  The same labels and order are used in tables and figures. Use original
  exposure codes for `code.exposure.ref`.

- order.strata:

  Exposure codes in display order. Codes omitted from this vector follow
  in their default order.

- level.strata:

  Expected exposure codes. Unnamed labels are matched to these codes by
  position. If omitted, the levels in the data are used.

- outcome.type:

  Shared outcome type, `"survival"` or `"competing-risk"`. The default
  `NULL` (also `"auto"`) uses the same rule as
  [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md):
  more than two observed event codes indicate competing risks. Otherwise
  the outcome is treated as survival. Detection uses the population
  selected by `subset.condition` before splitting into subgroups.
  Specify `"competing-risk"` explicitly if an expected event code is
  absent.

- cifplot.args:

  Named list of
  [`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md)
  display arguments shared by all figures. Use the individual arguments
  for axis labels and scales, legends, confidence bands, risk and
  estimate tables, event marks, quantiles, and styling. For example,
  `list(label.x = "Months", add.risktable = FALSE, legend.position = "bottom")`.
  Set stratum labels and order through the report's strata arguments.
  Use `add.conf` to show or hide the bands and `conf.int` at report
  level to set their confidence level. By default the y-axis follows the
  outcome's estimator. `type.y` changes the figure's scale only.
  Competing-event marks use the observed event times for each
  competing-risk outcome. Settings apply to single-curve figures, with
  multiple exposure groups overlaid when present.

- output.file:

  Older name for `filename.cifreport`. Supply only one.

## Value

Invisibly, a list with `output.file` (absolute paths in the order
requested, or `NULL` if no file was saved), `specifications`, `tables`,
`figures`, `fits`, `status` and `versions`. Tables, figures and fits are
nested by population and grouping definition. Table 1 is a gtsummary
object. Tables 2–5 are numeric result data frames with outcome and group
identifiers. Group codes and display labels are available in Tables 2–4.
Figures and fits are indexed by outcome, with risk-ratio fits indexed by
time.

## Details

Each population and grouping variable produces a section containing a
baseline table, full-follow-up event counts, cumulative counts and
probabilities at requested times, crude risk ratios, curve comparison
tests, and outcome-specific plots. Observed cumulative counts include
events at exactly the requested time. Table 3 and, by default, the
figure use the selected estimator scale. Percentages in the tables are
displayed to one decimal place by default. KM is survival, 1-KM is event
probability, AJ is primary-event cumulative incidence, and 1-AJ is its
complement. Probabilities and confidence limits beyond the last observed
follow-up in a group are omitted. Risk ratios compare event
probabilities (1 minus survival for survival outcomes), with one
unadjusted fit per time. Competing-risk fits jointly model both causes.
Log-rank and Gray tests compare curves over follow-up. Table titles use
the exposure's variable label (`label` attribute) when available and
include its displayed levels when it has value labels. Subgroup titles
use the variable and value labels when available. For exposure levels,
`label.strata` takes precedence over the data's value labels. Unlabelled
codes are displayed as given. Within a subgroup level, an outcome with
zero observed events over the entire follow-up stops after Table 2. For
competing risks, this requires zero primary and zero competing events
across the exposure groups.

Missing baseline values are summarized by default. Missing group values
skip the population/group combination. Missing outcome values skip that
outcome within the population. Single-level groups retain descriptive
results. A missing reference level affects risk ratios only. The
analysis log records invalid, skipped, unsupported and failed analyses
and warnings.

The report supports unweighted, right-censored data with one row per
subject, one primary event code and, for competing risks, one competing
event code. File output requires optional reporting packages and Pandoc.
See
[`vignette("v7_reporting_results")`](https://gestimation.github.io/cifmodeling/articles/v7_reporting_results.md)
for examples and supported configurations.

## See also

[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md),
[`cifplot()`](https://gestimation.github.io/cifmodeling/reference/cifplot.md),
[`polyreg()`](https://gestimation.github.io/cifmodeling/reference/polyreg.md),
[`ciftest()`](https://gestimation.github.io/cifmodeling/reference/ciftest.md)

## Examples

``` r
if (requireNamespace("gtsummary", quietly = TRUE) &&
    utils::packageVersion("gtsummary") >= "2.6.1") {
  data(diabetes.complications)
  report <- cifreport(
    diabetes.complications,
    outcomes = list(Retinopathy = list(time = "t", event = "epsilon")),
    exposure = "fruitq1", characteristics = c("age", "sex"),
    time.point = 8, filename.cifreport = NULL
  )
}
```
