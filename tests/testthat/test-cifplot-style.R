test_that("unnamed palette with patterns assigns color and linetype in order", {
  example_data <- data.frame(
    t = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5, 1.2, 2.2, 3.2, 4.2),
    epsilon = c(1, 0, 2, 0, 0, 1, 0, 2, 2, 0, 1, 0),
    strata = factor(rep(c("Low", "Medium", "High"), each = 4), levels = c("Low", "Medium", "High"))
  )

  p <- cifplot(
    Event(t, epsilon) ~ strata,
    data = example_data,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    addConfidenceInterval = FALSE,
    addRiskTable = FALSE,
    addCensorMark = FALSE,
    addCompetingRiskMark = FALSE,
    addIntercurrentEventMark = FALSE,
    addQuantileLine = FALSE,
    palette = c("red", "blue-dot", "green-broken")
  )

  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))

  sc_col <- p$scales$get_scales("colour")
  sc_lin <- p$scales$get_scales("linetype")
  expect_s3_class(sc_col, "ScaleDiscrete")
  expect_s3_class(sc_lin, "ScaleDiscrete")
  expect_identical(sc_col$scale_name, "manual")
  expect_identical(sc_lin$scale_name, "manual")
  expect_equal(sc_col$palette(3), c("red", "blue", "green"))
  expect_equal(sc_lin$palette(3), c("solid", "dotted", "dashed"))
})

test_that("named palette maps by final labels", {
  example_data <- data.frame(
    t = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5, 1.2, 2.2, 3.2, 4.2),
    epsilon = c(1, 0, 2, 0, 0, 1, 0, 2, 2, 0, 1, 0),
    strata = factor(rep(c("Low", "Medium", "High"), each = 4), levels = c("Low", "Medium", "High"))
  )

  p <- cifplot(
    Event(t, epsilon) ~ strata,
    data = example_data,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    addConfidenceInterval = FALSE,
    addRiskTable = FALSE,
    addCensorMark = FALSE,
    addCompetingRiskMark = FALSE,
    addIntercurrentEventMark = FALSE,
    addQuantileLine = FALSE,
    label.strata = c(Low = "Group A", Medium = "Group B", High = "Group C"),
    palette = c("Group A" = "red-dot", "Group B" = "blue", "Group C" = "green-broken")
  )

  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))

  sc_col <- p$scales$get_scales("colour")
  sc_lin <- p$scales$get_scales("linetype")
  expect_identical(sc_col$palette(3), c("red", "blue", "green"))
  expect_equal(sc_lin$palette(3), c("dotted", "solid", "dashed"))
})

test_that("'dot' keyword yields black dotted line", {
  example_data <- data.frame(
    t = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5, 1.2, 2.2, 3.2, 4.2),
    epsilon = c(1, 0, 2, 0, 0, 1, 0, 2, 2, 0, 1, 0),
    strata = factor(rep(c("Low", "Medium", "High"), each = 4), levels = c("Low", "Medium", "High"))
  )

  p <- cifplot(
    Event(t, epsilon) ~ strata,
    data = example_data,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    addConfidenceInterval = FALSE,
    addRiskTable = FALSE,
    addCensorMark = FALSE,
    addCompetingRiskMark = FALSE,
    addIntercurrentEventMark = FALSE,
    addQuantileLine = FALSE,
    palette = c("dot", "red", "blue")
  )

  sc_col <- p$scales$get_scales("colour")
  sc_lin <- p$scales$get_scales("linetype")
  expect_equal(sc_col$palette(3)[1], "black")
  expect_equal(sc_lin$palette(3)[1], "dotted")
})
