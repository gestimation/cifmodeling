test_that("cifpanel() accepts pre-built plots (grid mode)", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(disp, qsec)) + ggplot2::geom_point()

  out <- cifpanel(
    plots = list(p1, p2),
    rows.columns.panel = c(1, 2),
    legend.collect = TRUE,
    legend.position = "bottom",
    print.panel = FALSE
  )

  expect_true(is.list(out))
  expect_true(all(c("list.plot", "patchwork") %in% names(out)))
  expect_equal(length(out$list.plot), 2L)
})

test_that("cifpanel() accepts pre-built plots (inset mode)", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(disp, qsec)) + ggplot2::geom_point()

  out <- cifpanel(
    plots = list(p1, p2),
    inset.panel = TRUE,
    title.plot = c("Base", "Inset"),
    print.panel = FALSE
  )

  expect_true(is.list(out))
  expect_true(all(c("list.plot", "patchwork") %in% names(out)))
})


make_min_data <- function() {
  data.frame(
    time = c(1, 2, 3),
    status = c(0, 1, 2),
    A = c(0, 1, 0)
  )
}

make_min_inputs <- function() {
  list(
    K = 1L,
    formulas = list(Event(time, status) ~ A),
    data = make_min_data(),
    code.events = list(c(1, 2, 0)),
    outcome.flags = c("C"),
    outcome.list = list("competing-risk"),
    typey.list = list("risk"),
    labely.list = list("CIF"),
    labelx.list = list("Time"),
    limsx.list = NULL,
    limsy.list = NULL,
    breakx.list = NULL,
    breaky.list = NULL,
    addCI.list = list(TRUE),
    addCen.list = list(TRUE),
    addCR.list = list(FALSE),
    addIC.list = list(FALSE),
    addQ.list = list(FALSE),
    strata.list = NULL,
    legend.position = "top",
    dots = list(style = "classic", font.family = NULL, font.size = NULL)
  )
}

make_inputs_2panel <- function(df) {
  K <- 2L
  list(
    K = K,
    formulas = rep(list(Event(t, epsilon) ~ fruitq), K),
    data = df,
    code.events = list(c(1, 2, 0), c(2, 1, 0)),
    outcome.flags = c("C", "C"),
    outcome.list = list("competing-risk", "competing-risk"),
    typey.list = list("risk", "risk"),
    labely.list = list("Diabetic retinopathy", "Macrovascular complications"),
    labelx.list = list("Years from registration", "Years from registration"),
    limsx.list = NULL,
    limsy.list = NULL,
    breakx.list = NULL,
    breaky.list = NULL,
    addCI.list = list(TRUE, TRUE),
    addCen.list = list(TRUE, TRUE),
    addCR.list = list(FALSE, FALSE),
    addIC.list = list(FALSE, FALSE),
    addQ.list = list(FALSE, FALSE),
    strata.list = NULL,
    legend.position = "top",
    dots = list(style = "classic", font.family = NULL, font.size = NULL)
  )
}

test_that("cifpanel() produces expected outputs with competing risks data (no plotting)", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggsurvfit")

  data(diabetes.complications)

  out <- cifpanel(
    title.panel      = "A comparison of cumulative incidence of competing events",
    rows.columns.panel = c(1, 2),
    formula          = Event(t, epsilon) ~ fruitq,
    data             = diabetes.complications,
    outcome.type     = "competing-risk",
    code.events      = list(c(1, 2, 0), c(2, 1, 0)),
    label.y          = c("Diabetic retinopathy", "Macrovascular complications"),
    label.x          = "Years from registration",
    subtitle.panel   = "Stratified by fruit intake",
    caption.panel    = "Data: diabetes.complications",
    title.plot       = c("Diabetic retinopathy", "Macrovascular complications"),
    legend.position  = "bottom",
    legend.collect   = TRUE,
    print.panel      = FALSE
  )
  expect_true(is.list(out))
  expect_true(all(c("list.plot", "patchwork") %in% names(out)))
})

test_that("panel_prepare() returns curves and plot arguments", {
  skip_on_cran()
  skip_if_not_installed("survival")
  inputs <- make_min_inputs()
  prep <- do.call(cifmodeling:::panel_prepare, inputs)

  expect_true(is.list(prep$curves))
  expect_equal(length(prep$curves), inputs$K)
  expect_true(all(vapply(prep$curves, inherits, logical(1), what = "survfit")))

  expect_true(is.list(prep$plot_args))
  expect_equal(length(prep$plot_args), inputs$K)
  expect_true(all(vapply(prep$plot_args, is.list, logical(1))))
  expect_true(all(vapply(prep$plot_args, function(x) identical(x$legend.position, inputs$legend.position), logical(1))))

  expect_identical(prep$K, inputs$K)
})

test_that("panel_prepare() returns per-panel objects with expected structure", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggsurvfit")

  data(diabetes.complications)

  inputs <- make_inputs_2panel(diabetes.complications)
  prep <- do.call(cifmodeling:::panel_prepare, inputs)

  expect_true(is.list(prep))
  expect_identical(prep$K, inputs$K)

  expect_true(is.list(prep$curves))
  expect_equal(length(prep$curves), inputs$K)
  expect_true(all(vapply(prep$curves, inherits, logical(1), what = "survfit")))

  expect_true(is.list(prep$plot_args))
  expect_equal(length(prep$plot_args), inputs$K)
  expect_true(all(vapply(prep$plot_args, is.list, logical(1))))

  expect_true(all(vapply(prep$plot_args, function(x) identical(x$legend.position, inputs$legend.position), logical(1))))
  expect_true(identical(prep$plot_args[[1]]$label.x, inputs$labelx.list[[1]]))
  expect_true(identical(prep$plot_args[[1]]$label.y, inputs$labely.list[[1]]))
  expect_true(identical(prep$plot_args[[2]]$label.x, inputs$labelx.list[[2]]))
  expect_true(identical(prep$plot_args[[2]]$label.y, inputs$labely.list[[2]]))
})

test_that("cifpanel() with two formulas shows per-plot titles", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(5, 7, 9, 10, 3, 4, 6, 8),
    epsilon = c(1, 0, 2, 0, 1, 2, 0, 1),
    x1      = factor(c("A","B","A","B","A","B","A","B")),
    x2      = c(0,1,0,1,1,0,1,0)
  )

  res <- cifpanel(
    formulas     = list(Event(t, epsilon) ~ x1,
                        Event(t, epsilon) ~ x2),
    data         = df,
    outcome.type = "competing-risk",
    code.events  = list(c(1, 2, 0), c(1, 2, 0)),
    title.plot   = c("Plot-x1", "Plot-x2"),
    print.panel  = FALSE
  )

  expect_type(res, "list")
  expect_length(res$list.plot, 2)
  expect_equal(res$list.plot[[1]]$labels$title, "Plot-x1")
  expect_equal(res$list.plot[[2]]$labels$title, "Plot-x2")
})

test_that("panel flags work", {
  skip()
  skip_if_not_installed("survival")
  skip_if_not_installed("ggsurvfit")
  skip_if_not_installed("patchwork")

  data(diabetes.complications)

  f1 <- Event(t, epsilon) ~ fruitq
  f2 <- Event(t, epsilon) ~ sex

  res <- cifpanel(
    formulas           = list(f1, f2),
    data               = diabetes.complications,
    outcome.type       = "competing-risk",
    panel.per.event    = TRUE,
    panel.censoring    = FALSE,
    panel.per.variable = TRUE,
    print.panel        = FALSE
  )

  expect_true(is.list(res))
  expect_true(
    inherits(res$patchwork, "patchwork") || inherits(res$patchwork, "gg")
  )
})

test_that("cifpanel() applies panel-level tag_levels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(5, 7, 9, 10),
    epsilon = c(1, 0, 2, 0),
    x       = c(0, 1, 0, 1)
  )

  res <- cifpanel(
    formulas     = list(Event(t, epsilon) ~ x,
                        Event(t, epsilon) ~ 1),
    data         = df,
    outcome.type = "competing-risk",
    code.events  = list(c(1, 2, 0), c(1, 2, 0)),
    tag.panel = "A",
    print.panel  = FALSE
  )

  ann <- res$patchwork$patches$annotation
  expect_false(is.null(ann))
  expect_equal(ann$tag_levels, "A")
})

test_that("cifplot(panel.per.event=TRUE) creates 2-event panel", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(2, 4, 6, 8, 10, 12),
    epsilon = c(1, 0, 2, 0, 1, 2),
    trt     = c(0, 1, 0, 1, 0, 1)
  )

  p <- cifplot(
    Event(t, epsilon) ~ trt,
    data             = df,
    outcome.type     = "competing-risk",
    code.events      = c(1, 2, 0),
    panel.per.event   = TRUE,
    add.risktable     = FALSE,
    label.y          = c("Cumulative incidence of interest", "Cumulative incidence of competing risk"),
    print.panel      = FALSE
  )

  expect_s3_class(p, "cifpanel")
  expect_true(inherits(p$patchwork, "patchwork"))
  expect_equal(length(p$list.plot), 2)

  expect_equal(p$list.plot[[1]]$labels$y, "Cumulative incidence of interest")
  expect_equal(p$list.plot[[2]]$labels$y, "Cumulative incidence of competing risk")
})

test_that("cifpanel() can take per-panel label.x/label.y", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(1, 2, 3, 4, 5, 6),
    epsilon = c(1, 0, 2, 0, 1, 2),
    g       = c(0, 0, 1, 1, 0, 1)
  )

  res <- cifpanel(
    formulas     = list(Event(t, epsilon) ~ g,
                        Event(t, epsilon) ~ 1),
    data         = df,
    outcome.type = "competing-risk",
    code.events  = list(c(1, 2, 0), c(1, 2, 0)),
    label.y      = list("CIF by group", "Overall CIF"),
    label.x      = list("Time (days)", "Time (days)"),
    print.panel  = FALSE
  )

  expect_equal(res$list.plot[[1]]$labels$y, "CIF by group")
  expect_equal(res$list.plot[[2]]$labels$y, "Overall CIF")
})

test_that("panel_recycle_to() enforces 1-or-n length rule", {
  rec <- cifmodeling:::panel_recycle_to

  x1 <- list("A")
  out1 <- rec(x1, n = 3L)
  expect_equal(length(out1), 3L)
  expect_equal(out1, list("A", "A", "A"))

  x2 <- list("A", "B")
  out2 <- rec(x2, n = 2L)
  expect_identical(out2, x2)

  x_bad <- list("A", "B", "C")
  expect_error(
    rec(x_bad, n = 2L),
    regexp = "Length mismatch"
  )
})

test_that("cifpanel() recycles scalar panel-wise arguments to all panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(5, 7, 9, 10),
    epsilon = c(1, 0, 2, 0),
    x       = c(0, 1, 0, 1)
  )

  res <- cifpanel(
    formulas     = list(Event(t, epsilon) ~ x,
                        Event(t, epsilon) ~ 1),
    data         = df,
    outcome.type = "competing-risk",
    code.events  = list(c(1, 2, 0), c(1, 2, 0)),
    label.y      = "Shared Y label",   # 長さ1 → リサイクル
    print.panel  = FALSE
  )

  expect_length(res$list.plot, 2L)
  expect_equal(res$list.plot[[1]]$labels$y, "Shared Y label")
  expect_equal(res$list.plot[[2]]$labels$y, "Shared Y label")
})

test_that("cifpanel() errors when panel-wise argument length is not 1 or K", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(5, 7, 9, 10),
    epsilon = c(1, 0, 2, 0),
    x       = c(0, 1, 0, 1)
  )

  # formulas が 2 なので K = 2
  expect_error(
    cifpanel(
      formulas     = list(Event(t, epsilon) ~ x,
                          Event(t, epsilon) ~ 1),
      data         = df,
      outcome.type = "competing-risk",
      code.events  = list(c(1, 2, 0), c(1, 2, 0)),
      label.y      = c("Y1", "Y2", "Y3"),  # 長さ 3 → ダメ
      print.panel  = FALSE
    ),
    regexp = "Length mismatch|length .*1 or 2",
    ignore.case = TRUE
  )
})


test_that("cifpanel() recycles scalar panel-wise arguments to all panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(5, 7, 9, 10),
    epsilon = c(1, 0, 2, 0),
    x       = c(0, 1, 0, 1)
  )

  # formulas が 2 個 → K = 2
  res <- cifpanel(
    formulas     = list(Event(t, epsilon) ~ x,
                        Event(t, epsilon) ~ 1),
    data         = df,
    outcome.type = "competing-risk",
    code.events  = list(c(1, 2, 0), c(1, 2, 0)),
    label.y      = "Shared Y label",   # 長さ1（リサイクル期待）
    label.x      = "Time",             # 長さ1（リサイクル期待）
    print.panel  = FALSE
  )

  expect_equal(length(res$list.plot), 2L)
  expect_equal(res$list.plot[[1]]$labels$y, "Shared Y label")
  expect_equal(res$list.plot[[2]]$labels$y, "Shared Y label")
  expect_equal(res$list.plot[[1]]$labels$x, "Time")
  expect_equal(res$list.plot[[2]]$labels$x, "Time")
})


test_that("cifpanel() errors when panel-wise argument length is not 1 or K", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(5, 7, 9, 10),
    epsilon = c(1, 0, 2, 0),
    x       = c(0, 1, 0, 1)
  )

  expect_error(
    cifpanel(
      formulas          = list(Event(t, epsilon) ~ x,
                               Event(t, epsilon) ~ 1),
      data              = df,
      outcome.type      = "competing-risk",
      code.events       = list(c(1, 2, 0), c(1, 2, 0)),
      rows.columns.panel = c(1, 2),
      label.y           = c("Y1", "Y2", "Y3"),
      print.panel       = FALSE
    ),
    regexp = "Length mismatch"
  )
})

test_that("cifpanel() enforces 1-or-K rule for panel-wise listable arguments", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  # K = 2 パネルのベース設定
  df <- data.frame(
    t       = c(5, 7, 9, 10),
    epsilon = c(1, 0, 2, 0),
    g       = c(0, 1, 0, 1)
  )

  base_args <- list(
    formulas          = list(Event(t, epsilon) ~ g,
                             Event(t, epsilon) ~ 1),
    data              = df,
    outcome.type      = "competing-risk",
    code.events       = list(c(1, 2, 0), c(1, 2, 0)),
    rows.columns.panel = c(1, 2),
    print.panel       = FALSE
  )

  check_bad_length <- function(arg_name, bad_value) {
    args2 <- base_args
    args2[[arg_name]] <- bad_value

    expect_error(
      do.call(cifpanel, args2),
      regexp = "Length mismatch",
      info   = paste("Argument:", arg_name)
    )
  }

  ## ここから 1-or-K ルールを期待したい panel-wise 引数たち

  # 文字列系
  check_bad_length("outcome.type",   list("competing-risk", "survival", "competing-risk"))
  check_bad_length("type.y",         list("risk", "risk", "risk"))
  check_bad_length("label.x",        list("Time1", "Time2", "Time3"))
  check_bad_length("label.y",        list("Y1", "Y2", "Y3"))

  # 数値レンジ系
  check_bad_length("limits.x",       list(c(0, 5), c(0, 5), c(0, 5)))
  check_bad_length("limits.y",       list(c(0, 1), c(0, 1), c(0, 1)))
  check_bad_length("breaks.x",       list(c(0, 5, 10), c(0, 5, 10), c(0, 5, 10)))
  check_bad_length("breaks.y",       list(c(0, 0.5, 1), c(0, 0.5, 1), c(0, 0.5, 1)))

  # logical トグル系
  check_bad_length("add.conf",                 list(TRUE, FALSE, TRUE))
  check_bad_length("add.censor.mark",          list(TRUE, TRUE, FALSE))
  check_bad_length("add.competing.risk.mark",  list(FALSE, TRUE, FALSE))
  check_bad_length("add.intercurrent.event.mark", list(FALSE, FALSE, TRUE))
  check_bad_length("add.quantile",             list(FALSE, TRUE, TRUE))
})


test_that("normalize_strata_info() respects level/labels consistency", {
  lvl <- c("0", "1")
  lab <- c("A", "B")
  out <- cifmodeling:::normalize_strata_info(
    level.strata = lvl,
    label.strata = lab
  )
  expect_equal(out$level, lvl)
  expect_equal(unname(out$label_map), lab)
})


test_that("cifpanel() accepts panel-wise list arguments (happy path)", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  df <- data.frame(
    t       = c(2, 4, 6, 8, 10, 12),
    epsilon = c(1, 0, 2, 0, 1, 2),
    trt     = c(0, 1, 0, 1, 0, 1)
  )

  ## 2 パネル想定
  args <- list(
    formulas     = list(
      Event(t, epsilon) ~ trt,
      Event(t, epsilon) ~ 1
    ),
    data         = df,
    outcome.type = list("competing-risk", "competing-risk"),
    code.events  = list(c(1, 2, 0), c(1, 2, 0)),

    ## per-panel ラベル
    label.x      = list("Time (years)", "Time (years)"),
    label.y      = list("CIF by treatment", "Overall CIF"),

    ## per-panel limits / breaks（実際の値はそこまで重要でない）
    limits.y     = list(c(0, 1), c(0, 1)),
    breaks.y     = list(seq(0, 1, 0.25), seq(0, 0.5, 0.1)),

    ## per-panel トグル類
    add.conf                 = list(TRUE, FALSE),
    add.censor.mark          = list(TRUE, FALSE),
    add.competing.risk.mark  = list(FALSE, FALSE),
    add.intercurrent.event.mark = list(FALSE, FALSE),
    add.quantile             = list(FALSE, TRUE),

    ## レイアウト
    rows.columns.panel = c(1, 2),
    print.panel        = FALSE
  )

  res <- do.call(cifpanel, args)

  ## cifpanel オブジェクトとして正常に生成されるか
  expect_s3_class(res, "cifpanel")
  expect_true(inherits(res$patchwork, "patchwork"))
  expect_length(res$list.plot, 2L)

  ## per-panel の label.x / label.y がちゃんと効いているか
  expect_equal(res$list.plot[[1]]$labels$y, "CIF by treatment")
  expect_equal(res$list.plot[[2]]$labels$y, "Overall CIF")
  expect_equal(res$list.plot[[1]]$labels$x, "Time (years)")
  expect_equal(res$list.plot[[2]]$labels$x, "Time (years)")

  ## ここまで来ていれば、list で渡した各種フラグや limits/breaks が
  ## 少なくともエラーなく処理されていることは確認できる
})


testthat::test_that("cifpanel returns a patchwork object for multi-panel inputs", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")

  data(diabetes.complications, package = "cifmodeling")

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  res <- cifmodeling::cifpanel(
    formulas = list(
      cifmodeling::Event(t, epsilon) ~ fruitq,
      cifmodeling::Event(t, epsilon) ~ sex
    ),
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events = list(c(1, 2, 0), c(1, 2, 0)),
    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  testthat::expect_true(is.list(res))
  testthat::expect_true(!is.null(res$patchwork))
  testthat::expect_s3_class(res$patchwork, "patchwork")
})

# ---- helpers ---------------------------------------------------------------

.extract_panel_range <- function(p, axis = c("x", "y")) {
  axis <- match.arg(axis)
  b <- ggplot2::ggplot_build(p)
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

.extract_cifpanel_plots <- function(res) {
  pw <- NULL
  if (is.list(res) && !is.null(res$patchwork)) pw <- res$patchwork
  if (inherits(res, "patchwork")) pw <- res
  if (is.null(pw) || !inherits(pw, "patchwork")) stop("No patchwork found in cifpanel result.")

  plots <- pw$patches$plots
  plots <- Filter(Negate(is.null), plots)
  plots
}

# ---- axis regression tests -------------------------------------------------

testthat::test_that("cifpanel respects limits.x even when breaks.x is supplied (scale_x_continuous path)", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")

  data(diabetes.complications, package = "cifmodeling")

  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  res <- cifmodeling::cifpanel(
    formulas = list(
      cifmodeling::Event(t, epsilon) ~ fruitq,
      cifmodeling::Event(t, epsilon) ~ sex
    ),
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events = list(c(1, 2, 0), c(1, 2, 0)),
    rows.columns.panel = c(1, 2),
    axis.info = list(use.coord.cartesian = FALSE),

    limits.x = list(c(0, 120), c(0, 120)),
    breaks.x = list(seq(0, 120, 12), seq(0, 120, 12)),

    limits.y = list(c(0, 1), c(0, 1)),
    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  testthat::expect_true(is.list(res))
  testthat::expect_true(!is.null(res$list.plot))
  testthat::expect_gte(length(res$list.plot), 2L)

  plots <- res$list.plot[1:2]
  tol <- 1e-8

  for (p in plots) {
    testthat::expect_s3_class(p, "ggplot")

    b <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]

    xr <- NULL
    if (!is.null(pp$x.range)) {
      xr <- pp$x.range
    } else if (!is.null(pp$x) && !is.null(pp$x$range) && !is.null(pp$x$range$range)) {
      xr <- pp$x$range$range
    }

    testthat::expect_true(is.numeric(xr) && length(xr) == 2L)

    # 左端は ggsurvfit 互換の pad により負になり得るため主張しない
    testthat::expect_lte(xr[2], 120 + tol)

    # 任意：右端が意図せず縮み過ぎないことの弱い担保（安定なら残す）
    testthat::expect_gte(xr[2], 120 - tol)
  }
})


testthat::test_that("cifpanel respects limits.x with breaks.x when coord_cartesian is used", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")

  data(diabetes.complications, package = "cifmodeling")
  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  res <- cifmodeling::cifpanel(
    formulas = list(
      cifmodeling::Event(t, epsilon) ~ fruitq,
      cifmodeling::Event(t, epsilon) ~ sex
    ),
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events = list(c(1, 2, 0), c(1, 2, 0)),
    rows.columns.panel = c(1, 2),
    axis.info = list(use.coord.cartesian = TRUE),

    limits.x = list(c(0, 120), c(0, 120)),
    breaks.x = list(seq(0, 120, 12), seq(0, 120, 12)),
    limits.y = list(c(0, 1), c(0, 1)),

    add.risktable = FALSE,
    add.conf = FALSE,
    add.censor.mark = FALSE
  )

  testthat::expect_gte(length(res$list.plot), 2L)
  plots <- res$list.plot[1:2]

  for (p in plots) {
    b <- ggplot2::ggplot_build(p)
    pp <- b$layout$panel_params[[1]]
    xr <- pp$x.range %||% pp$x$range$range
    testthat::expect_lte(xr[2], 120 + 1e-8)
  }
})

testthat::test_that("cifpanel errors if a panel-wise list argument length is neither 1 nor K", {

  data(diabetes.complications, package = "cifmodeling")
  old <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old), add = TRUE)

  testthat::expect_error(
    cifmodeling::cifpanel(
      formulas = list(
        cifmodeling::Event(t, epsilon) ~ fruitq,
        cifmodeling::Event(t, epsilon) ~ sex
      ),
      data = diabetes.complications,
      outcome.type = "competing-risk",
      code.events = list(c(1, 2, 0), c(1, 2, 0)),
      rows.columns.panel = c(1, 2),

      # K=2 に対して長さ3（1でも2でもない）→ ここで落ちるのが期待
      limits.x = list(c(0, 120), c(0, 120), c(0, 120)),
      breaks.x = list(seq(0, 120, 12), seq(0, 120, 12)),
      limits.y = list(c(0, 1), c(0, 1)),

      add.risktable = FALSE,
      add.conf = FALSE,
      add.censor.mark = FALSE
    )
  )
})

