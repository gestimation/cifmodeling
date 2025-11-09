skip_on_cran()

test_that("new dot.case arguments are recognized", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    t = c(1, 2, 3, 4),
    d = c(0, 1, 2, 0),
    x = factor(c("A", "A", "B", "B"))
  )

  plt <- cifplot(
    Event(t, d) ~ x,
    data = df,
    outcome.type = "competing-risk",
    add.ci = TRUE,
    add.risktable = TRUE,
    add.censor.mark = TRUE,
    add.quantile = TRUE,
    level.quantile = c(0.5),
    use.coord.cartesian = TRUE,
    print.panel = FALSE
  )

  expect_true(inherits(plt, "gg"))
})

test_that("panel flags work", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t = c(1, 2, 3, 4, 5, 6),
    d = c(0, 1, 2, 0, 1, 2),
    x = factor(c("A", "A", "B", "B", "A", "B"))
  )

  pw <- cifpanel(
    formulas = list(Event(t, d) ~ x, Event(t, d) ~ 1),
    data = df,
    code.events = list(c(1, 2, 0), c(1, 2, 0)),
    panel.per.event = TRUE,
    panel.censoring = FALSE,
    panel.per.variable = TRUE,
    print.panel = FALSE
  )

  expect_true(inherits(pw, "patchwork") || inherits(pw, "gg"))
})
