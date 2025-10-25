test_that("default keeps ggsurvfit scales (no manual scales added)", {
  skip_if_not_installed("ggsurvfit")
  data(diabetes.complications, package = "ggsurvfit")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0)
  )
  expect_true(inherits(p, c("gg", "ggplot")) || inherits(p, "patchwork"))
  sc_col <- try(p$scales$get_scales("colour"), silent = TRUE)
  sc_lin <- try(p$scales$get_scales("linetype"), silent = TRUE)
  expect_true(!inherits(sc_col, "try-error"))
  expect_true(!inherits(sc_lin, "try-error"))
  if (!inherits(sc_col, "try-error")) {
    expect_false(inherits(sc_col, "ScaleDiscreteManual"))
  }
  if (!inherits(sc_lin, "try-error")) {
    expect_false(inherits(sc_lin, "ScaleDiscreteManual"))
  }
})

test_that("when all colors identical, linetype varies", {
  skip_if_not_installed("ggsurvfit")
  data(diabetes.complications, package = "ggsurvfit")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    palette = rep("red", 3)
  )
  expect_true(inherits(p, c("gg", "ggplot")) || inherits(p, "patchwork"))
  sc_lin <- p$scales$get_scales("linetype")
  lts <- sc_lin$get_limits()
  expect_gt(length(unique(lts)), 1)
})

test_that("when colors differ, all linetypes are solid", {
  skip_if_not_installed("ggsurvfit")
  data(diabetes.complications, package = "ggsurvfit")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    palette = c("red", "blue", "green")
  )
  sc_lin <- p$scales$get_scales("linetype")
  lts <- sc_lin$get_limits()
  expect_true(all(lts == "solid"))
})
