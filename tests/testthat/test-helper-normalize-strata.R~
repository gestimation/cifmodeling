test_that("numeric vectors with <9 unique values are converted to factor", {
  x <- c(1, 2, 2, 3, NA_real_)
  out <- cifmodeling:::cifplot_normalize_strata_var(x)
  expect_true(is.factor(out$values))
  expect_identical(levels(out$values), c("1", "2", "3"))
  expect_identical(as.character(out$values[1:4]), c("1", "2", "2", "3"))
  expect_true(is.na(out$values[5]))
  expect_identical(out$strategy, "factor")
})

test_that("numeric vectors with >=9 unique values are split at the median", {
  x <- 1:10
  out <- cifmodeling:::cifplot_normalize_strata_var(x)
  expect_true(is.factor(out$values))
  expect_identical(levels(out$values),
                   c(cifmodeling:::.cifplot_labels$strata$below_median,
                     cifmodeling:::.cifplot_labels$strata$above_median))
  med <- stats::median(x)
  expect_true(all(out$values[x <= med] == cifmodeling:::.cifplot_labels$strata$below_median))
  expect_true(all(out$values[x >  med] == cifmodeling:::.cifplot_labels$strata$above_median))
  expect_identical(out$strategy, "median")
})

test_that("formula normalization converts numeric RHS strata", {
  df <- data.frame(time = 1:10, status = c(0, 1, rep(0, 8)), z = seq_len(10))
  fml <- Event(time, status) ~ z
  out <- cifmodeling:::cifplot_normalize_formula_data(fml, df)
  expect_true(is.factor(out$data$z))
  expect_identical(levels(out$data$z),
                   c(cifmodeling:::.cifplot_labels$strata$below_median,
                     cifmodeling:::.cifplot_labels$strata$above_median))
})
