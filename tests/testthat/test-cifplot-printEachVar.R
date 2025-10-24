test_that("printEachVar returns one plot per RHS var", {
  data(diabetes.complications)
  diabetes.complications$sex_ <- as.factor(diabetes.complications$sex)
  plt <- cifplot(
    Event(t, epsilon) ~ sex_ + fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    printEachVar = TRUE,
    rows.columns.panel = c(1, 2)
  )
  expect_true(inherits(plt, "patchwork") || inherits(plt, "ggplot"))
  pa <- attr(plt, "plots")
  expect_false(is.null(pa))
  expect_length(pa, 2L)
  expect_true(all(vapply(pa, function(p) inherits(p, "ggplot"), logical(1))))
})

test_that("order.strata and label.strata (positional) are respected", {
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    printEachVar = TRUE,
    rows.columns.panel = c(1, 1),
    order.strata = list(fruitq = c("Q1", "Q2", "Q3", "Q4")),
    label.strata = list(fruitq = c("Q1", "Q2", "Q3", "Q4"))
  )
  expect_true(inherits(plt, "patchwork") || inherits(plt, "ggplot"))
})

test_that("label.strata named mapping works with order.strata", {
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    printEachVar = TRUE,
    rows.columns.panel = c(1, 1),
    order.strata = list(fruitq1 = c("Q1", "Q2", "Q3", "Q4")),
    label.strata = list(fruitq1 = c(Q1 = "Q1", Q2 = "Q2", Q3 = "Q3", Q4 = "Q4"))
  )
  expect_true(inherits(plt, "patchwork") || inherits(plt, "ggplot"))
})

test_that("numeric RHS var errors under printEachVar", {
  data(diabetes.complications)
  expect_error(
    cifplot(
      Event(t, epsilon) ~ age,
      data = diabetes.complications,
      outcome.type = "COMPETING-RISK",
      code.event1 = 1, code.event2 = 2, code.censoring = 0,
      printEachVar = TRUE
    ),
    "numeric.*discretize",
    ignore.case = TRUE
  )
})
