library(ggplot2)

test_that("named palette maps to strata labels accurately", {
  data(diabetes.complications)
  labs <- c(Q1 = "Quartile 1", Q2 = "Quartile 2", Q3 = "Quartile 3", Q4 = "Quartile 4")
  pal <- c(
    "Quartile 1" = "red",
    "Quartile 2" = "blue",
    "Quartile 3" = "green",
    "Quartile 4" = "purple"
  )

  p <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    label.strata = labs,
    palette = pal,
    addRiskTable = FALSE
  )

  expect_true(inherits(p, c("gg", "ggplot")) || inherits(p, "patchwork"))
  sc <- p$scales$get_scales("colour")
  expect_true(inherits(sc, "ScaleDiscrete"))
  labels_vec <- sc$labels
  if (is.null(labels_vec)) labels_vec <- character()
  expect_true(all(unname(labs) %in% labels_vec))
  values_vec <- sc$palette.cache
  if (is.null(values_vec)) values_vec <- sc$range$range
  expect_true(all(unname(pal) %in% values_vec))
})


test_that("unnamed palette uses order of strata", {
  data(diabetes.complications)
  p <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    order.strata = c("Q4", "Q2", "Q1", "Q3"),
    palette = c("red", "blue", "green", "orange"),
    addRiskTable = FALSE
  )
  expect_true(inherits(p, c("gg", "ggplot")) || inherits(p, "patchwork"))
})


test_that("preset palette key expands", {
  data(diabetes.complications)
  p <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    palette = "basic",
    addRiskTable = FALSE
  )
  expect_true(inherits(p, c("gg", "ggplot")) || inherits(p, "patchwork"))
})


test_that("single-stratum color/fill works", {
  data(diabetes.complications)
  p <- cifplot(
    Event(t, epsilon) ~ 1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    color = "blue",
    fill = "skyblue",
    addRiskTable = FALSE
  )
  expect_true(inherits(p, c("gg", "ggplot")) || inherits(p, "patchwork"))
  sc_col <- p$scales$get_scales("colour")
  sc_fil <- p$scales$get_scales("fill")
  if (!is.null(sc_col)) {
    vals <- sc_col$palette.cache
    if (is.null(vals)) vals <- sc_col$range$range
    expect_true("blue" %in% vals)
  }
  if (!is.null(sc_fil)) {
    vals <- sc_fil$palette.cache
    if (is.null(vals)) vals <- sc_fil$range$range
    expect_true("skyblue" %in% vals || "blue" %in% vals)
  }
})
