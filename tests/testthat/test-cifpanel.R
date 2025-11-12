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
    dots = list(style = "classsic", font.family = NULL, font.size = NULL)
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
    dots = list(style = "classsic", font.family = NULL, font.size = NULL)
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
