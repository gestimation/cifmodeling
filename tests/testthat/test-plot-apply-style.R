setup_plot_data <- function() {
  skip_if_not_installed("ggsurvfit")
  data(diabetes.complications, package = "ggsurvfit")
  diabetes.complications
}

test_that("label.strata only adjusts labels and suppresses fill legend", {
  df <- setup_plot_data()
  lbls <- c("0" = "Low intake", "1" = "High intake")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = df,
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
  df <- setup_plot_data()
  pal <- c("#FF0000", "#0000FF")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = df,
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
  df <- setup_plot_data()
  pal <- c("#FF0000", "#0000FF")
  lbls <- c("0" = "Group A", "1" = "Group B")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = df,
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

test_that("MONOCHROME style forces black color and linetype scale", {
  df <- setup_plot_data()
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = df,
    outcome.type = "COMPETING-RISK",
    code.events  = c(1, 2, 0),
    style = "MONOCHROME",
    addRiskTable = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  sc_lin <- p$scales$get_scales("linetype")
  sc_fill <- p$scales$get_scales("fill")

  expect_s3_class(sc_lin, "ScaleDiscreteManual")
  breaks <- sc_col$get_breaks()
  n_breaks <- if (length(breaks)) length(breaks) else 1L
  expect_equal(sc_col$palette(seq_len(n_breaks)), rep("black", n_breaks))
  exp_fill <- grDevices::gray(seq(0.85, 0.30, length.out = n_breaks))
  expect_equal(sc_fill$palette(seq_len(n_breaks)), exp_fill)
  expect_identical(sc_fill$guide, "none")
})

test_that("palette is ignored when MONOCHROME style is requested", {
  df <- setup_plot_data()
  pal <- c("#FF0000", "#0000FF")
  lbls <- c("0" = "Low", "1" = "High")
  expect_warning({
    p <- cifplot(
      Event(t, epsilon) ~ fruitq1,
      data = df,
      outcome.type = "COMPETING-RISK",
      code.events  = c(1, 2, 0),
      palette = pal,
      style = "MONOCHROME",
      label.strata = lbls,
      addRiskTable = FALSE
    )
  }, regexp = NA)

  sc_col <- p$scales$get_scales("colour")
  sc_lin <- p$scales$get_scales("linetype")

  breaks <- sc_col$get_breaks()
  n_breaks <- if (length(breaks)) length(breaks) else 1L
  expect_equal(sc_col$palette(seq_len(n_breaks)), rep("black", n_breaks))
  expect_s3_class(sc_lin, "ScaleDiscreteManual")
  expect_equal(sc_lin$get_labels(), unname(lbls))
})

test_that("order.strata sets scale limits", {
  df <- setup_plot_data()
  lbls <- c("0" = "Low", "1" = "High")
  ord <- c("1", "0")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = df,
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
