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
  expect_true(any(c("plots", "out_patchwork") %in% names(out)))
  expect_equal(length(out$plots), 2L)
})

test_that("cifpanel() accepts pre-built plots (inset mode)", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(disp, qsec)) + ggplot2::geom_point()

  out <- cifpanel(
    plots = list(p1, p2),
    use_inset_element = TRUE,
    title.plot = c("Base", "Inset"),
    print.panel = FALSE
  )

  expect_true(is.list(out))
  expect_true(any(c("plots", "out_patchwork") %in% names(out)))
})
