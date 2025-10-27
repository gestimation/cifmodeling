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

test_that("numeric vectors with <9 unique values are converted to factor", {
  x <- c(1, 2, 2, 3, NA_real_)
  out <- cifmodeling:::plot_normalize_strata(x)
  expect_true(is.factor(out$values))
  expect_identical(levels(out$values), c("1", "2", "3"))
  expect_identical(as.character(out$values[1:4]), c("1", "2", "2", "3"))
  expect_true(is.na(out$values[5]))
  expect_identical(out$strategy, "factor")
})

test_that("numeric vectors with >=9 unique values are split at the median", {
  x <- 1:10
  out <- cifmodeling:::plot_normalize_strata(x)
  expect_true(is.factor(out$values))
  expect_identical(levels(out$values),
                   c(cifmodeling:::plot_default_labels$strata$below_median,
                     cifmodeling:::plot_default_labels$strata$above_median))
  med <- stats::median(x)
  expect_true(all(out$values[x <= med] == cifmodeling:::plot_default_labels$strata$below_median))
  expect_true(all(out$values[x >  med] == cifmodeling:::plot_default_labels$strata$above_median))
  expect_identical(out$strategy, "median")
})

test_that("formula normalization converts numeric RHS strata", {
  df <- data.frame(time = 1:10, status = c(0, 1, rep(0, 8)), z = seq_len(10))
  fml <- Event(time, status) ~ z
  out <- cifmodeling:::plot_normalize_formula_data(fml, df)
  expect_true(is.factor(out$data$z))
  expect_identical(levels(out$data$z),
                   c(cifmodeling:::plot_default_labels$strata$below_median,
                     cifmodeling:::plot_default_labels$strata$above_median))
})
