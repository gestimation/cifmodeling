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





get_ggplot <- function(x) {
  if (inherits(x, "cifplot")) return(x$plot)
  x
}

extract_x_range <- function(x) {
  p <- get_ggplot(x)
  b <- ggplot2::ggplot_build(p)

  # ggplot2 のバージョン差分吸収
  panel_params <- b$layout$panel_params
  if (is.null(panel_params)) {
    stop("Cannot find panel_params in ggplot_build output.")
  }

  rngs <- lapply(panel_params, function(pp) {
    if (!is.null(pp$x.range)) return(pp$x.range)
    if (!is.null(pp$x$range)) return(pp$x$range)
    NULL
  })
  rngs <- Filter(Negate(is.null), rngs)

  if (length(rngs) == 0L) stop("Cannot extract x range from panel_params.")

  rng_mat <- do.call(rbind, rngs)
  c(min(rng_mat[, 1], na.rm = TRUE), max(rng_mat[, 2], na.rm = TRUE))
}

test_that("cifplot respects limits.x even when breaks.x is supplied", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")

  set.seed(1)
  n <- 300
  df <- data.frame(
    time   = pmin(rexp(n, rate = 1/80), 240),  # tmax > 120 を保証
    status = rbinom(n, 1, 0.6),
    group  = factor(sample(c("A", "B"), n, TRUE))
  )

  old <- getOption("warn")
  options(warn = 0)
  on.exit(options(warn = old), add = TRUE)

  obj <- cifplot(
    survival::Surv(time, status) ~ group,
    data         = df,
    outcome.type = "survival",
    limits.x     = c(0, 120),
    breaks.x     = seq(0, 120, 12),
    add.conf        = FALSE,
    add.risktable   = FALSE,
    add.censor.mark = FALSE
  )
  xr <- extract_x_range(obj)
  expect_lte(xr[2], 120 + 1e-8)
})




test_that("cifplot respects limits.x with breaks.x when coord_cartesian is used", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")

  set.seed(2)
  n <- 300
  df <- data.frame(
    time   = pmin(rexp(n, rate = 1/80), 240),
    status = rbinom(n, 1, 0.6),
    group  = factor(sample(c("A", "B"), n, TRUE))
  )

  old <- getOption("warn")
  options(warn = 0)
  on.exit(options(warn = old), add = TRUE)

  obj <- cifplot(
    survival::Surv(time, status) ~ group,
    data         = df,
    outcome.type = "survival",
    limits.x     = c(0, 120),
    breaks.x     = seq(0, 120, 12),
    use.coord.cartesian = TRUE,
    add.conf        = FALSE,
    add.risktable   = FALSE,
    add.censor.mark = FALSE
  )

  xr <- extract_x_range(obj)
  expect_lte(xr[2], 120 + 1e-8)
})


testthat::test_that("cifplot returns a cifplot object whose $plot is ggplot", {
  testthat::skip_if_not_installed("ggplot2")

  data(diabetes.complications, package = "cifmodeling")

  out <- cifmodeling::cifplot(
    cifmodeling::Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  testthat::expect_s3_class(out, "cifplot")
  testthat::expect_true(is.list(out))
  testthat::expect_true(!is.null(out$plot))
  testthat::expect_s3_class(out$plot, "ggplot")
})

# ---- helpers ---------------------------------------------------------------

.extract_panel_range <- function(p, axis = c("x", "y")) {
  axis <- match.arg(axis)
  b <- ggplot2::ggplot_build(p)

  # ggplot2 version differences: panel_params structure varies.
  pp <- b$layout$panel_params[[1]]

  if (axis == "x") {
    if (!is.null(pp$x.range)) return(pp$x.range)
    if (!is.null(pp$x$range$range)) return(pp$x$range$range)
    if (!is.null(pp$x$range)) return(pp$x$range)
    if (!is.null(b$layout$panel_scales_x[[1]]$range$range)) return(b$layout$panel_scales_x[[1]]$range$range)
  } else {
    if (!is.null(pp$y.range)) return(pp$y.range)
    if (!is.null(pp$y$range$range)) return(pp$y$range$range)
    if (!is.null(pp$y$range)) return(pp$y$range)
    if (!is.null(b$layout$panel_scales_y[[1]]$range$range)) return(b$layout$panel_scales_y[[1]]$range$range)
  }

  stop("Could not extract panel range from ggplot_build().")
}

.extract_x_breaks <- function(p) {
  b <- ggplot2::ggplot_build(p)
  pp <- b$layout$panel_params[[1]]

  # Try common fields first
  if (!is.null(pp$x$breaks)) return(pp$x$breaks)
  if (!is.null(pp$x.major))  return(pp$x.major)

  # Fallback: scale object after build
  sc <- p$scales$get_scales("x")
  if (!is.null(sc)) {
    br <- sc$get_breaks()
    return(br)
  }

  NULL
}

# ---- axis regression tests -------------------------------------------------

testthat::test_that("cifplot respects limits.x even when breaks.x is supplied (scale_x_continuous path)", {
  testthat::skip_if_not_installed("ggplot2")

  data(diabetes.complications, package = "cifmodeling")

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  out <- cifmodeling::cifplot(
    cifmodeling::Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    limits.x = c(0, 120),
    breaks.x = seq(0, 120, 12),
    use.coord.cartesian = FALSE,
    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  xr <- .extract_panel_range(out$plot, "x")
  testthat::expect_true(is.numeric(xr) && length(xr) == 2L)

  tol <- 1e-8
  # 左端は pad により負になり得るため主張しない
  testthat::expect_lte(xr[2], 120 + tol)
  testthat::expect_gte(xr[2], 120 - tol)  # 任意だが安定なら有益

  br <- .extract_x_breaks(out$plot)
  testthat::expect_true(is.numeric(br))

  # breaks は "含まれていること" を確認（厳密一致にしない）
  expected <- seq(0, 120, 12)
  testthat::expect_true(all(expected %in% br | is.na(expected)))
})


testthat::test_that("cifplot respects limits.x with breaks.x when coord_cartesian is used", {
  testthat::skip_if_not_installed("ggplot2")

  data(diabetes.complications, package = "cifmodeling")

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  out <- cifmodeling::cifplot(
    cifmodeling::Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    limits.x = c(0, 120),
    breaks.x = seq(0, 120, 12),
    use.coord.cartesian = TRUE,
    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  xr <- .extract_panel_range(out$plot, "x")
  testthat::expect_true(is.numeric(xr) && length(xr) == 2L)

  tol <- 1e-8

  # 左端は ggsurvfit 互換の pad により負になり得るため主張しない
  testthat::expect_lte(xr[2], 120 + tol)

  # 任意：左端が極端に負に飛んでいないことの弱い担保（安定なら残す）
  # 0..120 の表示で、-2.4 程度は許容、-50 などは異常として検知したい、という意図。
  testthat::expect_gt(xr[1], -60)
})


testthat::test_that("cifplot respects limits.y when breaks.y is supplied (no warnings; simple survival risk)", {
  testthat::skip_if_not_installed("ggplot2")

  df0 <- data.frame(
    time  = c(1, 2, 3, 4, 5, 6),
    status = c(0, 0, 0, 0, 0, 0),
    grp = rep(c("A", "B"), each = 3)
  )

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  out <- cifmodeling::cifplot(
    survival::Surv(time, status) ~ grp,
    data = df0,
    outcome.type = "survival",
    type.y = "risk",
    limits.y = c(0, 0.5),
    breaks.y = seq(0, 0.5, 0.1),
    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  yr <- .extract_panel_range(out$plot, "y")
  testthat::expect_true(is.numeric(yr) && length(yr) == 2L)
  testthat::expect_gte(yr[1], 0 - 1e-8)
  testthat::expect_lte(yr[2], 0.5 + 1e-8)
})

testthat::test_that("cifplot warns when breaks.x are outside limits.x (optional behavior check)", {
  testthat::skip_if_not_installed("ggplot2")

  data(diabetes.complications, package = "cifmodeling")

  old <- getOption("warn")
  options(warn = 1)
  on.exit(options(warn = old), add = TRUE)

  testthat::expect_warning(
    cifmodeling::cifplot(
      cifmodeling::Event(t, epsilon) ~ fruitq,
      data = diabetes.complications,
      outcome.type = "competing-risk",
      limits.x = c(0, 120),
      breaks.x = c(0, 60, 120, 240),  # 240 が outside
      use.coord.cartesian = FALSE,
      add.risktable = FALSE,
      add.conf = FALSE,
      add.censor.mark = FALSE
    ),
    regexp = "breaks\\.x|outside plotting range|plotting range",
    fixed = FALSE
  )
})


testthat::test_that("cifplot: type.y='cumhaz' matches -log(KM survival)", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,3,4,5,6),
    status = c(1,1,1,1,1,0)
  )

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  sf <- survival::survfit(survival::Surv(time, status) ~ 1, data = df)
  surv5 <- base::summary(sf, times = 5)$surv
  testthat::expect_true(is.numeric(surv5) && length(surv5) == 1L)
  testthat::expect_gt(surv5, 0)

  expected <- -log(surv5)

  out <- cifmodeling::cifplot(
    survival::Surv(time, status) ~ 1,
    data = df,
    outcome.type = "survival",
    code.event1 = 1,
    code.censoring = 0,
    type.y = "cumhaz",
    add.conf = FALSE,
    add.risktable = FALSE,
    add.censor.mark = FALSE
  )

  p <- out$plot %||% out
  testthat::expect_s3_class(p, "ggplot")

  b <- ggplot2::ggplot_build(p)

  # Find the curve layer (x/y numeric, no ribbons since add.conf=FALSE)
  get_xy_layer <- function(build_obj) {
    for (i in seq_along(build_obj$data)) {
      d <- build_obj$data[[i]]
      if (is.data.frame(d) && nrow(d) > 0L &&
          "x" %in% names(d) && "y" %in% names(d) &&
          is.numeric(d$x) && is.numeric(d$y)) {
        return(d)
      }
    }
    stop("No layer with numeric x/y found.")
  }

  d1 <- get_xy_layer(b)
  y_max <- max(d1$y, na.rm = TRUE)

  testthat::expect_true(is.finite(y_max))
  testthat::expect_equal(y_max, expected, tolerance = 1e-6)

  # Cumhaz should be non-decreasing over time (step curve)
  o <- order(d1$x)
  yy <- d1$y[o]
  testthat::expect_true(all(diff(yy) >= -1e-10))
})


testthat::test_that("cifplot: type.y='cloglog' matches log(-log(KM survival))", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,3,4,5,6),
    status = c(1,1,1,1,1,0)
  )

  sf <- survival::survfit(survival::Surv(time, status) ~ 1, data = df)
  surv5 <- base::summary(sf, times = 5)$surv
  testthat::expect_gt(surv5, 0)

  expected <- log(-log(surv5))

  out <- cifmodeling::cifplot(
    survival::Surv(time, status) ~ 1,
    data = df,
    outcome.type = "survival",
    code.event1 = 1,
    code.censoring = 0,
    type.y = "cloglog",
    add.conf = FALSE,
    add.risktable = FALSE,
    add.censor.mark = FALSE
  )

  p <- out$plot %||% out
  testthat::expect_s3_class(p, "ggplot")

  b <- ggplot2::ggplot_build(p)

  get_xy_layer <- function(build_obj) {
    for (i in seq_along(build_obj$data)) {
      d <- build_obj$data[[i]]
      if (is.data.frame(d) && nrow(d) > 0L &&
          "x" %in% names(d) && "y" %in% names(d) &&
          is.numeric(d$x) && is.numeric(d$y)) {
        return(d)
      }
    }
    stop("No layer with numeric x/y found.")
  }

  d1 <- get_xy_layer(b)
  y_max <- max(d1$y, na.rm = TRUE)

  testthat::expect_true(is.finite(y_max))
  testthat::expect_equal(y_max, expected, tolerance = 1e-6)

  # cloglog(-) should also be non-decreasing over time (because -log(S) increases)
  o <- order(d1$x)
  yy <- d1$y[o]
  testthat::expect_true(all(diff(yy) >= -1e-10))
})


testthat::test_that("cifpanel: per-panel type.y supports cumhaz and cloglog", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,3,4,5,6),
    status = c(1,1,1,1,1,0)
  )

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  sf <- survival::survfit(survival::Surv(time, status) ~ 1, data = df)
  surv5 <- base::summary(sf, times = 5)$surv
  expected_cumhaz  <- -log(surv5)
  expected_cloglog <- log(-log(surv5))

  res <- cifmodeling::cifpanel(
    formulas = list(
      survival::Surv(time, status) ~ 1,
      survival::Surv(time, status) ~ 1
    ),
    data = df,
    outcome.type = "survival",
    code.events = list(c(1,0), c(1,0)),
    rows.columns.panel = c(1, 2),

    # per-panel y scale
    type.y = list("cumhaz", "cloglog"),

    add.conf = FALSE,
    add.risktable = FALSE,
    add.censor.mark = FALSE,

    # Avoid warn=2 surprises from auto-limits in edge environments
    limits.y = list(c(0, 5), c(-5, 5))
  )

  testthat::expect_true(is.list(res))
  testthat::expect_true(!is.null(res$list.plot))
  testthat::expect_gte(length(res$list.plot), 2L)

  plots <- res$list.plot[1:2]
  for (p in plots) testthat::expect_s3_class(p, "ggplot")

  b1 <- ggplot2::ggplot_build(plots[[1]])
  b2 <- ggplot2::ggplot_build(plots[[2]])

  get_xy_layer <- function(build_obj) {
    for (i in seq_along(build_obj$data)) {
      d <- build_obj$data[[i]]
      if (is.data.frame(d) && nrow(d) > 0L &&
          "x" %in% names(d) && "y" %in% names(d) &&
          is.numeric(d$x) && is.numeric(d$y)) {
        return(d)
      }
    }
    stop("No layer with numeric x/y found.")
  }

  y1 <- max(get_xy_layer(b1)$y, na.rm = TRUE)
  y2 <- max(get_xy_layer(b2)$y, na.rm = TRUE)

  testthat::expect_equal(y1, expected_cumhaz,  tolerance = 1e-6)
  testthat::expect_equal(y2, expected_cloglog, tolerance = 1e-6)
})

testthat::test_that("cifplot: cumhaz/cloglog respect limits.y and breaks.y (robust)", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("survival")

  data(diabetes.complications, package = "cifmodeling")
  df <- diabetes.complications
  df$status1 <- ifelse(df$epsilon == 1, 1L, 0L)

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  # panel_params から y.range と breaks を取り出す（ggplot2差分吸収）
  get_panel_yrange <- function(p) {
    b <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]
    pp$y.range %||% pp$y$range$range
  }
  get_panel_ybreaks <- function(p) {
    b <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]
    y <- pp$y
    out <- NULL
    if (!is.null(y$breaks)) out <- y$breaks
    if (is.null(out) && is.function(y$get_breaks)) out <- y$get_breaks()
    if (is.null(out) && is.function(y$break_positions)) out <- y$break_positions()
    out
  }

  # cloglog の -Inf を避けるため、最初のイベント以降を x 範囲にする
  first_event_time <- suppressWarnings(min(df$t[df$status1 == 1], na.rm = TRUE))
  if (!is.finite(first_event_time)) testthat::skip("No events found for KM in this dataset.")
  x_start <- first_event_time + 1e-6
  x_end   <- min(max(df$t, na.rm = TRUE), x_start + 120)

  cases <- list(
    list(type.y = "cumhaz",  limits.y = c(0, 2),  breaks.y = seq(0, 2, 0.5)),
    list(type.y = "cloglog", limits.y = c(-4, 1), breaks.y = seq(-4, 1, 1))
  )

  for (use_cc in c(FALSE, TRUE)) {
    for (cs in cases) {
      res <- cifmodeling::cifplot(
        formula       = survival::Surv(t, status1) ~ fruitq,
        data          = df,
        outcome.type  = "survival",
        code.events   = c(1, 0),
        type.y        = cs$type.y,

        # ここが重要：cloglogの -Inf(=t0, S=1) を検査領域から外す
        limits.x      = c(x_start, x_end),

        limits.y      = cs$limits.y,
        breaks.y      = cs$breaks.y,

        use.coord.cartesian = use_cc,  # ← axis.info ではなく ... で渡す
        add.conf      = FALSE,
        add.risktable = FALSE,
        add.censor.mark = FALSE
      )

      testthat::expect_true(is.list(res))
      testthat::expect_s3_class(res$plot, "ggplot")
      p <- res$plot

      yr <- get_panel_yrange(p)
      testthat::expect_true(is.numeric(yr) && length(yr) == 2L)
      testthat::expect_gte(yr[1], cs$limits.y[1] - 1e-6)
      testthat::expect_lte(yr[2], cs$limits.y[2] + 1e-6)

      # breaks は「指定した breaks のうち limits 内のものが、実際の breaks に含まれる」ことを確認
      yb <- get_panel_ybreaks(p)
      testthat::expect_true(is.numeric(yb) && length(yb) >= 2L)

      exp_in <- cs$breaks.y[cs$breaks.y >= cs$limits.y[1] & cs$breaks.y <= cs$limits.y[2]]
      matched <- vapply(exp_in, function(e) any(abs(yb - e) < 1e-6), logical(1))
      testthat::expect_true(all(matched))
    }
  }
})

testthat::test_that("cifpanel: per-panel cumhaz/cloglog respect limits.y and breaks.y (robust)", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("survival")

  data(diabetes.complications, package = "cifmodeling")
  df <- diabetes.complications
  df$status1 <- ifelse(df$epsilon == 1, 1L, 0L)

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  get_panel_yrange <- function(p) {
    b <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]
    pp$y.range %||% pp$y$range$range
  }
  get_panel_ybreaks <- function(p) {
    b <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]
    y <- pp$y
    out <- NULL
    if (!is.null(y$breaks)) out <- y$breaks
    if (is.null(out) && is.function(y$get_breaks)) out <- y$get_breaks()
    if (is.null(out) && is.function(y$break_positions)) out <- y$break_positions()
    out
  }

  first_event_time <- suppressWarnings(min(df$t[df$status1 == 1], na.rm = TRUE))
  if (!is.finite(first_event_time)) testthat::skip("No events found for KM in this dataset.")
  x_start <- first_event_time + 1e-6
  x_end   <- min(max(df$t, na.rm = TRUE), x_start + 120)

  typey <- list("cumhaz", "cloglog")
  limsy <- list(c(0, 2), c(-4, 1))
  brksy <- list(seq(0, 2, 0.5), seq(-4, 1, 1))

  for (use_cc in c(FALSE, TRUE)) {
    res <- cifmodeling::cifpanel(
      formulas = list(
        survival::Surv(t, status1) ~ fruitq,
        survival::Surv(t, status1) ~ sex
      ),
      data = df,
      outcome.type = "survival",
      code.events  = list(c(1, 0), c(1, 0)),
      rows.columns.panel = c(1, 2),

      type.y   = typey,
      limits.x = c(x_start, x_end),
      limits.y = limsy,
      breaks.y = brksy,

      use.coord.cartesian = use_cc, # ← ... で渡す
      add.conf = FALSE,
      add.risktable = FALSE,
      add.censor.mark = FALSE
    )

    testthat::expect_true(is.list(res))
    testthat::expect_gte(length(res$list.plot), 2L)

    plots <- res$list.plot[1:2]
    for (i in 1:2) {
      p <- plots[[i]]
      testthat::expect_s3_class(p, "ggplot")

      yr <- get_panel_yrange(p)
      testthat::expect_gte(yr[1], limsy[[i]][1] - 1e-6)
      testthat::expect_lte(yr[2], limsy[[i]][2] + 1e-6)

      yb <- get_panel_ybreaks(p)
      testthat::expect_true(is.numeric(yb) && length(yb) >= 2L)

      exp_in <- brksy[[i]][brksy[[i]] >= limsy[[i]][1] & brksy[[i]] <= limsy[[i]][2]]
      matched <- vapply(exp_in, function(e) any(abs(yb - e) < 1e-6), logical(1))
      testthat::expect_true(all(matched))
    }
  }
})

