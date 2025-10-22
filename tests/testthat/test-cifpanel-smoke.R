test_that("cifpanel returns non-null structure (smoke)", {
  skip_on_cran()

  d <- data.frame(
    time = c(1, 2, 3),
    status = c(0, 1, 2),
    A = c(0, 1, 0)
  )

  out <- cifpanel(
    formula = Event(time, status) ~ A,
    data = d,
    outcome.type = "COMPETING-RISK",
    code.events = list(c(1, 2, 0)),
    print.panel = FALSE
  )

  expect_false(is.null(out))
  expect_true(is.list(out))
  expect_true(all(c("plots", "out_patchwork") %in% names(out)))
})
