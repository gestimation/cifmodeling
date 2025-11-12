get_plot_from_cifplot <- function(...) {
  res <- cifplot(...)
  testthat::expect_s3_class(res, "cifplot")
  testthat::expect_true(inherits(res$plot, "ggplot"))
  res$plot
}

get_panel_from_cifpanel <- function(...) {
  res <- cifpanel(...)
  testthat::expect_s3_class(res, "cifpanel")
  testthat::expect_true(is.list(res$list.plot))
  testthat::expect_true(all(vapply(res$list.plot, inherits, logical(1), what = "ggplot")))
  testthat::expect_true(inherits(res$patchwork, "patchwork"))
  res
}

expect_no_print <- function(expr) {
  out <- capture.output(res <- force(expr))
  testthat::expect_equal(length(out), 0L)
  invisible(res)
}

test_that("label.strata only adjusts labels and suppresses fill legend", {
  skip()
  data(diabetes.complications)
  label <- c("Low intake", "High intake")
  level <- c(0, 1)
  p <- get_plot_from_cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    label.strata = label,
    level.strata = level,
    add.risktable = FALSE,
    print.panel = FALSE  # ← 明示
  )

  sc_col   <- p$scales$get_scales("colour")
  sc_lin   <- p$scales$get_scales("linetype")
  sc_fill  <- p$scales$get_scales("fill")
  sc_shape <- p$scales$get_scales("shape")

  expect_setequal(sc_col$get_labels(), unname(label))
  expect_setequal(sc_lin$get_labels(), unname(label))
  expect_identical(sc_fill$guide, "none")
  expect_identical(sc_shape$guide, "none")
})

test_that("label.strata is reflected in color legend NEW", {
  data(diabetes.complications)
  p <- get_plot_from_cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    label.strata = c("Low intake", "High intake"),
    level.strata = c(0, 1),
    add.risktable = FALSE,
    print.panel = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  br     <- sc_col$get_breaks()

  get_labels_safe <- function(sc, br) {
    lab <- sc$labels
    if (is.function(lab)) {
      lab(br)
    } else if (is.null(lab)) {
      br
    } else {
      lab
    }
  }

  lbl <- get_labels_safe(sc_col, br)
  expect_equal(unname(lbl), c("Low intake", "High intake"))
})

test_that("label.strata overrides palette labels", {
  skip()
  data(diabetes.complications)
  pal  <- c("#FF0000", "#0000FF")
  lbls <- c("0" = "Group A", "1" = "Group B")
  p <- get_plot_from_cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    palette = pal,
    label.strata = lbls,
    add.risktable = FALSE,
    print.panel = FALSE
  )

  sc_col  <- p$scales$get_scales("colour")
  sc_fill <- p$scales$get_scales("fill")

  expect_setequal(sc_col$get_labels(), unname(lbls))
  expect_equal(sc_col$palette(seq_along(pal)), pal)
  expect_identical(sc_fill$guide, "none")
})

test_that("palette only uses manual color scale", {
  data(diabetes.complications)
  pal <- c("#FF0000", "#0000FF")
  p <- get_plot_from_cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    palette = pal,
    add.risktable = FALSE,
    print.panel = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  expect_s3_class(sc_col, "ScaleDiscreteManual")
  expect_identical(sc_col$scale_name, "manual")
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


test_that("cifplot does not print by default", {
  data(diabetes.complications)
  expect_no_print(
    cifplot(Event(t, epsilon) ~ fruitq1,
            data = diabetes.complications,
            outcome.type = "competing-risk",
            code.events = c(1,2,0))
  )
})



test_that("cifpanel returns list.plot and patchwork (no out_patchwork)", {
  data(diabetes.complications)
  res <- get_panel_from_cifpanel(
    formulas = list(
      Event(t, epsilon) ~ fruitq1,
      Event(t, epsilon) ~ 1
    ),
    data = diabetes.complications,
    code.events = list(c(1,2,0), c(1,2,0)),
    rows.columns.panel = c(1,2),
    print.panel = FALSE
  )
  expect_true(is.null(res$plot))
  expect_true(is.list(res$list.plot))
  expect_true(inherits(res$patchwork, "patchwork"))
  expect_false("out_patchwork" %in% names(res))
})


test_that("survfit.info includes meta fields and data.name", {
  data(diabetes.complications)
  res <- cifplot(Event(t, epsilon) ~ fruitq1,
                 data = diabetes.complications,
                 outcome.type = "competing-risk",
                 code.events = c(1,2,0))
  s <- res$survfit.info
  expect_true(is.list(s))
  expect_true(all(c("formula_or_fit","outcome.type","code.event1",
                    "code.event2","code.censoring","code.events","data.name") %in% names(s)))
  expect_type(s$data.name, "character")
})
