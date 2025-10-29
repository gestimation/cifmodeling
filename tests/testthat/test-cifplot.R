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

test_that("order & labels must match as sets; otherwise order is ignored", {
  cur_full  <- c("grp=A","grp=B","grp=C")
  cur_short <- c("A","B","C")

  lbl <- c("Alpha","Bravo","Charlie")
  names(lbl) <- cur_full

  ord <- c("grp=B","grp=C","grp=A")
  out <- plot_reconcile_order_and_labels(cur_full, cur_short, ord, lbl)
  expect_true(out$used_order)
  expect_identical(out$limits_arg, ord)

  ord_bad <- c("grp=A","grp=B")
  expect_warning(
    out2 <- plot_reconcile_order_and_labels(cur_full, cur_short, ord_bad, lbl),
    "`order.strata` and `label.strata` must contain the same set of levels"
  )
  expect_false(out2$used_order)
  expect_identical(out2$limits_arg, names(lbl))

  ord_none <- c("grp=X","grp=Y")
  expect_warning(
    out3 <- plot_reconcile_order_and_labels(cur_full, cur_short, ord_none, NULL),
    "`order.strata` has no overlap"
  )
  expect_null(out3$limits_arg)
  expect_true(out3$forbid_limits_due_to_order)

  out4 <- plot_reconcile_order_and_labels(cur_full, cur_short, NULL, lbl)
  expect_identical(out4$limits_arg, names(lbl))
})
