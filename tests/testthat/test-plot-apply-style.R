test_that("label.strata only adjusts labels and suppresses fill legend", {
  data(diabetes.complications)
  lbls <- c("0" = "Low intake", "1" = "High intake")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    label.strata = lbls,
    addRiskTable = FALSE
  )

  expect_s3_class(p, "ggplot")
  sc_col <- p$scales$get_scales("colour")
  sc_lin <- p$scales$get_scales("linetype")
  sc_fill <- p$scales$get_scales("fill")
  sc_shape <- p$scales$get_scales("shape")

  expect_equal(sc_col$get_labels(), unname(lbls))
  expect_equal(sc_lin$get_labels(), unname(lbls))
  expect_identical(sc_fill$guide, "none")
  expect_identical(sc_shape$guide, "none")
})

test_that("palette only uses manual color scale", {
  data(diabetes.complications)
  pal <- c("#FF0000", "#0000FF")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    palette = pal,
    addRiskTable = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  expect_s3_class(sc_col, "ScaleDiscreteManual")
  expect_identical(sc_col$scale_name, "manual")
  expect_equal(sc_col$palette(seq_along(pal)), pal)

  sc_fill <- p$scales$get_scales("fill")
  expect_identical(sc_fill$guide, "none")
})

test_that("label.strata overrides palette labels", {
  data(diabetes.complications)
  pal <- c("#FF0000", "#0000FF")
  lbls <- c("0" = "Group A", "1" = "Group B")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    palette = pal,
    label.strata = lbls,
    addRiskTable = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  expect_equal(sc_col$get_labels(), unname(lbls))
  expect_equal(sc_col$palette(seq_along(pal)), pal)

  sc_fill <- p$scales$get_scales("fill")
  expect_identical(sc_fill$guide, "none")
})

test_that("order.strata sets scale limits", {
  data(diabetes.complications)
  lbls <- c("0" = "Low", "1" = "High")
  ord <- c("1", "0")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    label.strata = lbls,
    order.strata = ord,
    addRiskTable = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  expect_identical(sc_col$get_limits(), ord)
  expect_equal(sc_col$get_labels(), unname(lbls[ord]))
})
