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
    outcome.list = list("COMPETING-RISK"),
    typey.list = list("risk"),
    labely.list = list("CIF"),
    typex.list = NULL,
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
    dots = list(style = "CLASSIC", font.family = NULL, font.size = NULL)
  )
}


test_that("panel_prepare returns curves and plot args", {
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
