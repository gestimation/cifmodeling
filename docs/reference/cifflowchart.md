# Create a flowchart of exclusions, treatment groups, and outcome status

`cifflowchart()` creates a lightweight flowchart of exclusions,
treatment groups, and outcome status with `DiagrammeR`. It can be used
as a simple CONSORT-like flowchart helper, but it is not intended to
fully automate a CONSORT-compliant diagram.

## Usage

``` r
cifflowchart(
  formula,
  data,
  time.point = NULL,
  withdraw.consent = NULL,
  ineligible = NULL,
  pre.exclude = NULL,
  post.exclude = NULL,
  subset.condition = NULL,
  na.action = na.pass,
  outcome.type = NULL,
  code.event1 = 1,
  code.event2 = 2,
  code.censoring = 0,
  label.strata = NULL,
  order.strata = NULL,
  label.events = NULL,
  title = NULL,
  percent = TRUE,
  digits = 1,
  ...
)
```

## Arguments

- formula:

  A formula. Supported forms are `response ~ arm`, `response ~ 1`, and
  `Event(time, status) ~ arm`.

- data:

  A data frame.

- time.point:

  Optional time point for `Event(time, status)` formulas.

- withdraw.consent, ineligible, pre.exclude, post.exclude:

  Optional exclusion variables. Logical vectors use `TRUE` as excluded.
  Character and factor vectors use non-missing, non-empty values as
  exclusion reasons.

- subset.condition:

  Optional expression evaluated in `data` before counting.

- na.action:

  Included for consistency with other formula interfaces. The default is
  `na.pass` because missing outcomes and missing groups are displayed in
  the flowchart.

- outcome.type:

  Reserved for future use.

- code.event1, code.event2, code.censoring:

  Status codes for event-history formulas.

- label.strata:

  Optional labels for treatment groups.

- order.strata:

  Optional order for treatment groups.

- label.events:

  Optional labels for outcome states. A named vector may be used to map
  existing labels to display labels.

- title:

  Optional graph title.

- percent:

  Logical; show percentages in nodes.

- digits:

  Number of digits for percentages.

- ...:

  Reserved for future use.

## Value

A
[`DiagrammeR::grViz()`](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html)
htmlwidget.

## Lifecycle

**\[experimental\]**

## Examples

``` r
if (requireNamespace("DiagrammeR", quietly = TRUE)) {
  data(diabetes.complications)

  cifflowchart(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    time.point = 8,
    outcome.type = "competing-risk"
  )
}

{"x":{"diagram":"digraph cifflowchart {\ngraph [rankdir = TB];\nnode [shape = box, style = rounded, fontname = Helvetica];\nedge [fontname = Helvetica];\nall [label = \"All patients\\nn = 978 (100.0%)\"];\ngroup1 [label = \"0\\nn = 720 (73.6%)\"];\nall -> group1;\nout1_1 [label = \"Event 1 by 8\\nn = 188 (26.1%)\"];\ngroup1 -> out1_1;\nout1_2 [label = \"Event 2 by 8\\nn = 57 (7.9%)\"];\ngroup1 -> out1_2;\nout1_3 [label = \"Censored before 8\\nn = 225 (31.2%)\"];\ngroup1 -> out1_3;\nout1_4 [label = \"Event-free at 8\\nn = 250 (34.7%)\"];\ngroup1 -> out1_4;\ngroup2 [label = \"1\\nn = 258 (26.4%)\"];\nall -> group2;\nout2_1 [label = \"Event 1 by 8\\nn = 91 (35.3%)\"];\ngroup2 -> out2_1;\nout2_2 [label = \"Event 2 by 8\\nn = 22 (8.5%)\"];\ngroup2 -> out2_2;\nout2_3 [label = \"Censored before 8\\nn = 69 (26.7%)\"];\ngroup2 -> out2_3;\nout2_4 [label = \"Event-free at 8\\nn = 76 (29.5%)\"];\ngroup2 -> out2_4;\n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}
```
