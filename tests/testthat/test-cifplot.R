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

.get_x_range <- function(p) {
  b <- ggplot2::ggplot_build(p)
  pp <- b$layout$panel_params[[1]]

  xr <- pp$x.range
  if (is.null(xr)) xr <- pp$x$range$range
  xr
}

.expect_xrange_includes_limits_with_padding <- function(xr, lim, pad_frac_max = 0.25, tol = 1e-8) {
  testthat::expect_true(is.numeric(xr) && length(xr) == 2L)

  span <- lim[2] - lim[1]

  # 重要：指定した limits 区間を「含む」こと
  testthat::expect_lte(xr[1], lim[1] + tol)
  testthat::expect_gte(xr[2], lim[2] - tol)

  # 余白は「小さい」こと（パッチで 2%〜程度を想定、上限は緩めに 25%）
  testthat::expect_gte(xr[1], lim[1] - pad_frac_max * span - tol)
  testthat::expect_lte(xr[2], lim[2] + pad_frac_max * span + tol)

  # 今回の仕様：lim[1]=0 のときは左が負にずれる（原点が 0 固定にならない）
  if (abs(lim[1]) < 1e-12) testthat::expect_lt(xr[1], 0 + tol)

  # 右マージンも入れる仕様なら、右端は lim[2] より大きいはず
  testthat::expect_gt(xr[2], lim[2] - tol)
}

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

  lim <- c(0, 120)
  xr <- .get_x_range(out$plot)
  .expect_xrange_includes_limits_with_padding(xr, lim, pad_frac_max = 0.25)

  br <- .extract_x_breaks(out$plot)
  testthat::expect_true(is.numeric(br))
  testthat::expect_true(all(seq(0, 120, 12) %in% br))
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

  lim <- c(0, 120)
  xr <- .get_x_range(out$plot)
  .expect_xrange_includes_limits_with_padding(xr, lim, pad_frac_max = 0.25)

  br <- .extract_x_breaks(out$plot)
  testthat::expect_true(is.numeric(br))
  testthat::expect_true(all(seq(0, 120, 12) %in% br))
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

test_that("cifplot palette is not lost when level.strata/order.strata are supplied", {
  skip_if_not_installed("ggsurvfit")
  skip_if_not_installed("patchwork")

  data("diabetes.complications", package = "cifmodeling")

  get_main_plot <- function(x) {
    if (inherits(x, "patchwork")) {
      if (!is.null(x$plots) && length(x$plots) >= 1L) return(x$plots[[1L]])
      if (!is.null(x$patches$plots) && length(x$patches$plots) >= 1L) return(x$patches$plots[[1L]])
    }
    x
  }

  # ---- 4 strata example ----
  lv4  <- levels(diabetes.complications$fruitq)
  pal4 <- c("#E7B800", "#2E9FDF", "#E7B800", "#2E9FDF")

  obj4 <- cifplot(
    Event(t,epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    add.conf = FALSE,
    add.censor.mark = FALSE,
    add.risktable = TRUE,
    add.estimate.table = TRUE,
    style = "classic",
    limits.x = c(0, 12),
    breaks.x = seq(0, 12, 1),
    limits.y = c(0, 1),
    level.strata = lv4,
    order.strata = lv4,
    palette = pal4
  )

  main4 <- get_main_plot(obj4$plot)
  sc4 <- main4$scales$get_scales("colour")
  expect_true(inherits(sc4, "ScaleDiscreteManual"))

  lims4 <- sc4$get_limits()
  mapped4 <- sc4$map(lims4)
  expect_equal(unname(mapped4), pal4[seq_along(lims4)])

  # ---- 2 strata example (fruitq1) ----
  # fruitq1 の levels が無い場合にも備える
  lv2 <- levels(diabetes.complications$fruitq1)
  if (is.null(lv2)) lv2 <- sort(unique(as.character(diabetes.complications$fruitq1)))
  pal2 <- c("#E7B800", "#2E9FDF")

  obj2 <- cifplot(
    Event(t,epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    add.conf = FALSE,
    add.censor.mark = FALSE,
    add.risktable = TRUE,
    add.estimate.table = TRUE,
    style = "classic",
    limits.x = c(0, 12),
    breaks.x = seq(0, 12, 1),
    limits.y = c(0, 1),
    level.strata = lv2,
    order.strata = lv2,
    palette = pal2
  )

  main2 <- get_main_plot(obj2$plot)
  sc2 <- main2$scales$get_scales("colour")
  expect_true(inherits(sc2, "ScaleDiscreteManual"))

  lims2 <- sc2$get_limits()
  mapped2 <- sc2$map(lims2)
  expect_equal(unname(mapped2), pal2[seq_along(lims2)])
})


test_that("apply_strata_to_plots preserves manual colour/linetype palettes", {
  skip_if_not_installed("ggplot2")

  brks <- c("C", "A", "B")
  label_map <- c(A = "Alpha", B = "Beta", C = "Gamma")

  df <- data.frame(
    time   = rep(1:3, times = 3),
    y      = c(1, 2, 3,  1, 1.5, 2,  2, 2.5, 3),
    strata = factor(rep(c("A","B","C"), each = 3), levels = c("A","B","C"))
  )

  p0 <- ggplot2::ggplot(df, ggplot2::aes(time, y, colour = strata, linetype = strata)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::scale_color_manual(values = c(A = "#000000", B = "#111111", C = "#222222")) +
    ggplot2::scale_linetype_manual(values = c(A = "solid", B = "dashed", C = "dotted"))

  build0 <- ggplot2::ggplot_build(p0)$data[[1]][, c("colour","linetype")]

  p1 <- cifmodeling:::apply_strata_to_plots(list(p0),
                                            order_data = brks, label_map = label_map, touch_colour = TRUE
  )[[1]]

  build1 <- ggplot2::ggplot_build(p1)$data[[1]][, c("colour","linetype")]

  # manual palette が壊れていない（= colour/linetype の割当が同一）
  expect_identical(build1$colour,   build0$colour)
  expect_identical(build1$linetype, build0$linetype)

  # breaks/labels が意図通りに更新されている（scale を差し替えずに）
  sc_col <- p1$scales$get_scales("colour")
  sc_lt  <- p1$scales$get_scales("linetype")
  expect_identical(sc_col$breaks, brks)
  expect_identical(sc_lt$breaks,  brks)
  expect_identical(sc_col$labels, unname(label_map[brks]))
  expect_identical(sc_lt$labels,  unname(label_map[brks]))
})

test_that("apply_strata_to_plots preserves manual fill palette (e.g., CI ribbon)", {
  skip_if_not_installed("ggplot2")

  brks <- c("B", "C", "A")
  label_map <- c(A = "Alpha", B = "Beta", C = "Gamma")

  df <- data.frame(
    time   = rep(1:3, times = 3),
    y      = c(1, 2, 3,  1, 1.5, 2,  2, 2.5, 3),
    strata = factor(rep(c("A","B","C"), each = 3), levels = c("A","B","C"))
  )
  df$ymin <- df$y - 0.1
  df$ymax <- df$y + 0.1

  p0 <- ggplot2::ggplot(df, ggplot2::aes(time, y, group = strata)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax, fill = strata),
                         alpha = 0.3, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(colour = strata), linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c(A = "#f0f0f0", B = "#d0d0d0", C = "#b0b0b0")) +
    ggplot2::scale_color_manual(values = c(A = "#000000", B = "#111111", C = "#222222"))

  # ribbon data は data[[1]] に入ることが多い（geom_ribbon が先）
  build0 <- ggplot2::ggplot_build(p0)
  fill0  <- build0$data[[1]]$fill

  p1 <- cifmodeling:::apply_strata_to_plots(list(p0),
                                            order_data = brks, label_map = label_map, touch_colour = TRUE
  )[[1]]

  build1 <- ggplot2::ggplot_build(p1)
  fill1  <- build1$data[[1]]$fill

  # manual fill が壊れていない（= fill の割当が同一）
  expect_identical(fill1, fill0)

  sc_fill <- p1$scales$get_scales("fill")
  expect_identical(sc_fill$breaks, brks)
  expect_identical(sc_fill$labels, unname(label_map[brks]))
})

test_that("apply_strata_to_plots is idempotent (does not accumulate scales)", {
  skip_if_not_installed("ggplot2")

  brks <- c("C", "A", "B")
  label_map <- c(A = "Alpha", B = "Beta", C = "Gamma")

  df <- data.frame(
    time   = rep(1:3, times = 3),
    y      = c(1, 2, 3,  1, 1.5, 2,  2, 2.5, 3),
    strata = factor(rep(c("A","B","C"), each = 3), levels = c("A","B","C"))
  )

  p0 <- ggplot2::ggplot(df, ggplot2::aes(time, y, colour = strata)) +
    ggplot2::geom_line()

  p1 <- cifmodeling:::apply_strata_to_plots(list(p0),
                                            order_data = brks, label_map = label_map, touch_colour = TRUE
  )[[1]]
  n_scales_1 <- length(p1$scales$scales)
  build1 <- ggplot2::ggplot_build(p1)$data[[1]]

  p2 <- cifmodeling:::apply_strata_to_plots(list(p1),
                                            order_data = brks, label_map = label_map, touch_colour = TRUE
  )[[1]]
  n_scales_2 <- length(p2$scales$scales)
  build2 <- ggplot2::ggplot_build(p2)$data[[1]]

  # 2回目で scale が増殖しない
  expect_identical(n_scales_2, n_scales_1)

  # 描画結果も変わらない（少なくとも data レベルで）
  expect_identical(build2$colour, build1$colour)
  expect_identical(build2$group,  build1$group)
})

test_that("touch_colour=FALSE does not inject an explicit colour scale", {
  skip_if_not_installed("ggplot2")

  brks <- c("C", "A", "B")
  label_map <- c(A = "Alpha", B = "Beta", C = "Gamma")

  df <- data.frame(
    time   = rep(1:3, times = 3),
    y      = c(1, 2, 3,  1, 1.5, 2,  2, 2.5, 3),
    strata = factor(rep(c("A","B","C"), each = 3), levels = c("A","B","C"))
  )

  p0 <- ggplot2::ggplot(df, ggplot2::aes(time, y, colour = strata)) +
    ggplot2::geom_line()

  # もともと明示的な colour scale は無い（default は build 時に入る）
  expect_null(p0$scales$get_scales("colour"))

  p_no_touch <- cifmodeling:::apply_strata_to_plots(list(p0),
                                                    order_data = brks, label_map = label_map, touch_colour = FALSE
  )[[1]]
  expect_null(p_no_touch$scales$get_scales("colour"))

  p_touch <- cifmodeling:::apply_strata_to_plots(list(p0),
                                                 order_data = brks, label_map = label_map, touch_colour = TRUE
  )[[1]]
  expect_false(is.null(p_touch$scales$get_scales("colour")))
})

test_that("apply_strata_to_plots is a no-op when order_data/label_map is NULL or plots empty", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:3, y = 1:3, strata = factor(c("A","B","C")))
  p0 <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = strata)) + ggplot2::geom_line()

  # order_data NULL -> 早期 return（同一オブジェクト）
  out1 <- cifmodeling:::apply_strata_to_plots(list(p0), order_data = NULL, label_map = c(A="a"), touch_colour = TRUE)
  expect_identical(out1[[1]], p0)

  # label_map NULL -> 早期 return
  out2 <- cifmodeling:::apply_strata_to_plots(list(p0), order_data = c("A","B","C"), label_map = NULL, touch_colour = TRUE)
  expect_identical(out2[[1]], p0)

  # plots empty -> そのまま
  out3 <- cifmodeling:::apply_strata_to_plots(list(), order_data = c("A"), label_map = c(A="a"), touch_colour = TRUE)
  expect_identical(out3, list())
})


testthat::test_that("normalize_strata_info() handles unnamed and named label.strata", {
  norm <- cifmodeling:::normalize_strata_info(
    level.strata = c("A","B"),
    order.strata = NULL,
    label.strata = c("Arm A","Arm B")   # unnamed, length matches
  )
  testthat::expect_identical(norm$level, c("A","B"))
  testthat::expect_identical(norm$order_data, c("A","B"))
  testthat::expect_identical(unname(norm$label_map), c("Arm A","Arm B"))
  testthat::expect_identical(names(norm$label_map), c("A","B"))

  norm2 <- cifmodeling:::normalize_strata_info(
    level.strata = c("A","B"),
    order.strata = c("B","A"),
    label.strata = c(A="Arm A", C="EXTRA")  # B missing, C extra
  )
  testthat::expect_identical(norm2$order_data, c("B","A"))
  # B should be filled with itself (default)
  testthat::expect_identical(unname(norm2$label_map[c("A","B")]), c("Arm A","B"))
})

testthat::test_that("normalize_strata_info() returns NULL order/label when order.strata invalid", {
  norm <- cifmodeling:::normalize_strata_info(
    level.strata = c("A","B"),
    order.strata = c("B","Z"),          # invalid
    label.strata = c(A="Arm A", B="Arm B")
  )
  testthat::expect_identical(norm$level, c("A","B"))
  testthat::expect_null(norm$order_data)
  testthat::expect_null(norm$label_map)
})

testthat::test_that("normalize_strata_info() returns all NULL when level.strata missing/empty", {
  norm1 <- cifmodeling:::normalize_strata_info(level.strata = NULL, order.strata = NULL, label.strata = NULL)
  testthat::expect_null(norm1$level)
  testthat::expect_null(norm1$order_data)
  testthat::expect_null(norm1$label_map)

  norm2 <- cifmodeling:::normalize_strata_info(level.strata = character(0), order.strata = NULL, label.strata = NULL)
  testthat::expect_null(norm2$level)
  testthat::expect_null(norm2$order_data)
  testthat::expect_null(norm2$label_map)
})

testthat::test_that("apply_strata_to_plots() is no-op when order/label missing", {
  testthat::skip_if_not_installed("ggplot2")
  df <- data.frame(x = 1:2, y = 1:2, strata = factor(c("A","B")))
  p  <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = strata)) + ggplot2::geom_point()

  out1 <- cifmodeling:::apply_strata_to_plots(list(p), order_data = NULL, label_map = c(A="a", B="b"))
  out2 <- cifmodeling:::apply_strata_to_plots(list(p), order_data = c("B","A"), label_map = NULL)

  testthat::expect_true(inherits(out1[[1]], "ggplot"))
  testthat::expect_true(inherits(out2[[1]], "ggplot"))
})

testthat::test_that("apply_strata_to_plots() updates breaks/labels and preserves manual colour mapping", {
  testthat::skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:4, y = 1:4, strata = factor(c("A","B","A","B")))
  p0 <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = strata)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c(A = "#111111", B = "#999999"))

  out <- cifmodeling:::apply_strata_to_plots(
    plots         = list(p0),
    order_data    = c("B","A"),
    label_map     = c(A = "Arm A", B = "Arm B"),
    touch_colour  = TRUE
  )[[1]]

  sc1 <- out$scales$get_scales("colour")
  if (is.null(sc1)) sc1 <- out$scales$get_scales("color")

  # breaks/labels は指定順に更新されていること
  testthat::expect_identical(sc1$breaks, c("B","A"))
  testthat::expect_identical(sc1$labels, c("Arm B","Arm A"))

  # 実際の色が manual のまま（touch_colour=TRUE が manual palette を潰していない）
  gb <- ggplot2::ggplot_build(out)
  cols <- gb$data[[1]]$colour

  testthat::expect_equal(length(cols), nrow(df))
  testthat::expect_true(all(cols[df$strata == "A"] == "#111111"))
  testthat::expect_true(all(cols[df$strata == "B"] == "#999999"))
})

testthat::test_that("apply_strata_to_plots() does not add colour changes when touch_colour = FALSE", {
  testthat::skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:4, y = 1:4, strata = factor(c("A","B","A","B")))
  p0 <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = strata)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c(A = "#111111", B = "#999999"),
                                breaks = c("A","B"),
                                labels = c("A","B"))

  out <- cifmodeling:::apply_strata_to_plots(
    plots         = list(p0),
    order_data    = c("B","A"),
    label_map     = c(A = "Arm A", B = "Arm B"),
    touch_colour  = FALSE
  )[[1]]

  sc1 <- out$scales$get_scales("colour")
  if (is.null(sc1)) sc1 <- out$scales$get_scales("color")

  # touch_colour=FALSE では colour scale の breaks/labels を触らない
  testthat::expect_identical(sc1$breaks, c("A","B"))
  testthat::expect_identical(sc1$labels, c("A","B"))

  # manual palette が実際に使われている（build後の colour を検証）
  gb <- ggplot2::ggplot_build(out)
  cols <- gb$data[[1]]$colour

  testthat::expect_equal(length(cols), nrow(df))
  testthat::expect_true(all(cols[df$strata == "A"] == "#111111"))
  testthat::expect_true(all(cols[df$strata == "B"] == "#999999"))
})

testthat::test_that("apply_strata_to_plots() tolerates plot lists with mixed aesthetics", {
  testthat::skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:4, y = 1:4, strata = factor(c("A","B","A","B")))

  p_strata <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = strata)) +
    ggplot2::geom_point()

  p_nostrata <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  out <- cifmodeling:::apply_strata_to_plots(
    plots       = list(p_strata, p_nostrata),
    order_data  = c("B","A"),
    label_map   = c(A = "Arm A", B = "Arm B"),
    touch_colour = TRUE
  )

  testthat::expect_equal(length(out), 2L)
  testthat::expect_true(inherits(out[[1]], "ggplot"))
  testthat::expect_true(inherits(out[[2]], "ggplot"))

  # second plot should still build cleanly (no accidental breakage)
  testthat::expect_silent(ggplot2::ggplot_build(out[[2]]))
})

testthat::test_that("cifpanel(): order/label are applied to returned plots (engine=cifplot)", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("survival")

  dat <- data.frame(
    t       = c(1,2,3,4,5,6,7,8),
    epsilon = c(1,0,1,0,1,0,1,0),
    arm     = factor(c("A","A","B","B","A","B","A","B"))
  )

  res <- cifmodeling::cifpanel(
    formula      = cifmodeling::Event(t, epsilon) ~ arm,
    data         = dat,
    code.events  = list(c(1, 0)),  # survival setting
    level.strata = c("A","B"),
    order.strata = c("B","A"),
    label.strata = c(A="Arm A", B="Arm B"),
    rows.columns.panel = c(1,1),
    legend.position = "top"
  )

  p1 <- res$list.plot[[1]]
  sc <- p1$scales$get_scales("colour")
  if (is.null(sc)) sc <- p1$scales$get_scales("color")

  # if your apply_strata_to_plots() implementation only edits existing scales,
  # this verifies cifpanel produced a colour scale and it was re-labeled.
  testthat::expect_false(is.null(sc))
  testthat::expect_identical(sc$breaks, c("B","A"))
  testthat::expect_identical(sc$labels, c("Arm B","Arm A"))
})

testthat::test_that("cifpanel(): invalid order.strata triggers warning and disables relabeling", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("survival")

  dat <- data.frame(
    t       = c(1,2,3,4,5,6),
    epsilon = c(1,0,1,0,1,0),
    arm     = factor(c("A","A","B","B","A","B"))
  )

  res <- testthat::expect_warning(
    cifmodeling::cifpanel(
      formula      = cifmodeling::Event(t, epsilon) ~ arm,
      data         = dat,
      code.events  = list(c(1, 0)),
      level.strata = c("A","B"),
      order.strata = c("B","Z"),           # invalid
      label.strata = c(A="Arm A", B="Arm B"),
      rows.columns.panel = c(1,1)
    ),
    regexp = "order\\.strata has unknown levels"
  )

  testthat::expect_null(res$axis.info$order.strata)
  testthat::expect_null(res$axis.info$label.strata)
})
