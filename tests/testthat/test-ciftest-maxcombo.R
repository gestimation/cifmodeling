maxcombo_test_data <- function() {
  data.frame(
    time = rep(1:12, 2),
    status = c(1, 1, 2, 1, 1, 0, 1, 2, 1, 0, 0, 0,
               0, 2, 0, 2, 1, 2, 0, 1, 2, 0, 0, 0),
    group = factor(rep(c("Control", "Treatment"), each = 12))
  )
}

testthat::test_that("MaxCombo defaults to 00+20+01 and reverses one-sided tests with group order", {
  dat <- maxcombo_test_data()
  form <- Event(time, status) ~ group
  two <- ciftest(form, dat, test = "maxcombo")
  less <- ciftest(form, dat, test = "maxcombo", alternative = "less")
  greater <- ciftest(form, dat, test = "maxcombo", alternative = "greater")
  reversed <- transform(dat, group = factor(group, levels = rev(levels(group))))
  rev_two <- ciftest(form, reversed, test = "maxcombo")
  rev_less <- ciftest(form, reversed, test = "maxcombo", alternative = "less")
  rev_greater <- ciftest(form, reversed, test = "maxcombo", alternative = "greater")

  testthat::expect_identical(two$test, "augmented")
  testthat::expect_identical(two$alternative, "two.sided")
  for (fit in list(two, less, greater)) {
    testthat::expect_identical(fit$directions$rho, c(0, 2, 0))
    testthat::expect_identical(fit$directions$gamma, c(0, 0, 1))
  }
  testthat::expect_identical(less$comparison, c(reference = "Control", comparison = "Treatment"))
  testthat::expect_identical(rev_less$comparison, c(reference = "Treatment", comparison = "Control"))
  testthat::expect_true(all(less$z < 0))
  testthat::expect_lt(less$p.value, greater$p.value)
  testthat::expect_equal(unname(less$z), -unname(rev_less$z), tolerance = 1e-12)
  testthat::expect_equal(less$p.value, rev_greater$p.value, tolerance = 1e-12)
  testthat::expect_equal(greater$p.value, rev_less$p.value, tolerance = 1e-12)
  testthat::expect_equal(two$p.value, rev_two$p.value, tolerance = 1e-12)
  testthat::expect_output(print(less), "Treatment vs Control, alternative = less")
  ref_less <- ciftest(form, dat, test = "m", alternative = "less", code.exposure.ref = "Treatment")
  testthat::expect_equal(ref_less$p.value, greater$p.value, tolerance = 1e-12)
  testthat::expect_equal(generics::tidy(less)$p.value, less$p.value)
})

testthat::test_that("m selects MaxCombo and code.exposure.ref controls the contrast", {
  dat <- maxcombo_test_data()
  dat$group <- as.integer(dat$group) - 1L
  attr(dat$group, "labels") <- c(Control = 0L, Treatment = 1L)
  form <- Event(time, status) ~ group
  for (alternative in c("two.sided", "less", "greater")) {
    fit <- ciftest(form, dat, test = "m", alternative = alternative, code.exposure.ref = 1L)
    full <- ciftest(form, dat, test = "maxcombo", alternative = alternative, code.exposure.ref = "1")
    opposite <- switch(alternative, two.sided = "two.sided", less = "greater", greater = "less")
    default <- ciftest(form, dat, test = "maxcombo", alternative = opposite)
    testthat::expect_identical(fit$test.requested, "maxcombo")
    testthat::expect_identical(fit$directions, full$directions)
    testthat::expect_identical(fit$p.value, full$p.value)
    testthat::expect_identical(fit$comparison, c(reference = "1", comparison = "0"))
    testthat::expect_equal(fit$p.value, default$p.value, tolerance = 1e-12)
    testthat::expect_equal(unname(fit$z), -unname(default$z), tolerance = 1e-12)
    testthat::expect_output(print(fit), "Control vs Treatment")
    testthat::expect_identical(fit$data, dat)
  }
  reversed <- transform(dat, group = factor(group, levels = c(1, 0)))
  fit <- ciftest(form, reversed, test = "m", alternative = "less", code.exposure.ref = 0)
  default <- ciftest(form, dat, test = "maxcombo", alternative = "less")
  testthat::expect_equal(fit$p.value, default$p.value, tolerance = 1e-12)
  testthat::expect_error(ciftest(form, dat, test = "m", code.exposure.ref = 2),
                         "must match a group in the analysis data")
  testthat::expect_error(ciftest(form, dat, test = "m", code.exposure.ref = c(0, 1)),
                         "must be one non-missing exposure code or level")
  testthat::expect_error(ciftest(form, dat, test = "ditzhaus", code.exposure.ref = 0),
                         'requires `test = "maxcombo"`')
})

testthat::test_that("MaxCombo accepts custom directions and requires augmented binary scores", {
  dat <- maxcombo_test_data()
  form <- Event(time, status) ~ group
  single <- ciftest(form, dat, test = "maxcombo", directions = "unweighted",
                    alternative = "less")
  testthat::expect_equal(single$p.value, unname(pnorm(single$z)))
  testthat::expect_error(ciftest(form, dat, test = "maxcombo", component.test = "gray"),
                         "Omit `component.test`")
  testthat::expect_error(ciftest(form, dat, test = "ditzhaus", alternative = "less"),
                         "requires `test")
  three <- transform(dat, group = factor(rep(c("A", "B", "C"), each = 8)))
  testthat::expect_error(ciftest(form, three, test = "maxcombo"), "exactly two groups")
  survival <- transform(dat, status = ifelse(status == 2, 0, status))
  testthat::expect_error(ciftest(form, survival, test = "maxcombo"), 'outcome.type = "survival"')
})

testthat::test_that("MaxCombo uses joint normal tail probabilities", {
  z <- c(-0.4, 1.25)
  for (alternative in c("two.sided", "less", "greater")) {
    fit <- ciftest_maxcombo_statistic(z, diag(2), alternative)
    marginal <- if (alternative == "two.sided") 2 * pnorm(-fit$statistic) else pnorm(-fit$statistic)
    testthat::expect_equal(fit$p.value, 1 - (1 - marginal)^2, tolerance = 1e-12)
    perfectly_correlated <- ciftest_maxcombo_statistic(rep(1, 2), matrix(1, 2, 2), alternative)
    single <- ciftest_maxcombo_statistic(1, matrix(1), alternative)
    testthat::expect_equal(perfectly_correlated$p.value, single$p.value, tolerance = 1e-12)
    fit3 <- ciftest_maxcombo_statistic(c(z, 1), diag(3), alternative)
    marginal3 <- if (alternative == "two.sided") 2 * pnorm(-fit3$statistic) else pnorm(-fit3$statistic)
    testthat::expect_equal(fit3$p.value, 1 - (1 - marginal3)^3, tolerance = 1e-9)
    perfectly_correlated3 <- ciftest_maxcombo_statistic(rep(1, 3), matrix(1, 3, 3), alternative)
    testthat::expect_equal(perfectly_correlated3$p.value, single$p.value, tolerance = 1e-9)
  }
  set.seed(1)
  state <- .Random.seed
  fit <- ciftest_maxcombo_statistic(c(z, 1, 0), diag(4), "greater")
  testthat::expect_equal(fit$p.value, 1 - pnorm(1.25)^4, tolerance = 1e-6)
  testthat::expect_identical(.Random.seed, state)
})

testthat::test_that("custom MaxCombo directions retain their dependent constraints", {
  dat <- maxcombo_test_data()
  directions <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  # Values from a high-accuracy joint normal reference calculation on this fit.
  expected <- c(two.sided = 0.1555895875424, less = 0.07779494432402,
                greater = 0.9088315669356)
  for (alternative in names(expected)) {
    fit <- ciftest(Event(time, status) ~ group, dat, test = "maxcombo",
                   directions = directions, alternative = alternative)
    testthat::expect_equal(fit$diagnostics$covariance.rank, 3)
    testthat::expect_lt(abs(fit$p.value - expected[[alternative]]), 5e-6)
  }
})

testthat::test_that("internal normal integration agrees with a conditional normal reference", {
  rho <- 0.7
  correlation <- matrix(rho, 4L, 4L) + diag(1 - rho, 4L)
  score <- c(-1.6, 0.2, 1.3, 2.2)
  set.seed(17)
  state <- .Random.seed
  for (alternative in c("two.sided", "less", "greater")) {
    fit <- ciftest_maxcombo_statistic(score, correlation, alternative)
    threshold <- fit$statistic
    two_sided <- alternative == "two.sided"
    # Here the coordinates are independent conditional on a common normal.
    reference <- stats::integrate(function(common) {
      mean <- sqrt(rho) * common
      upper <- pnorm((threshold - mean) / sqrt(1 - rho), lower.tail = FALSE)
      lower <- if (two_sided) pnorm((-threshold - mean) / sqrt(1 - rho)) else 0
      -expm1(4 * log1p(-pmin(1, lower + upper))) * dnorm(common)
    }, -Inf, Inf, abs.tol = 1e-12)$value
    testthat::expect_lt(abs(fit$p.value - reference), 5e-6)
  }
  testthat::expect_identical(.Random.seed, state)
})

testthat::test_that("very small MaxCombo tail probabilities retain relative precision", {
  # Four directions encode two independent normal variables and their negatives.
  correlation <- crossprod(cbind(diag(2), -diag(2)))
  t <- 8
  tail <- stats::pnorm(t, lower.tail = FALSE)
  expected <- -expm1(2 * log1p(-2 * tail))
  for (alternative in c("two.sided", "less", "greater")) {
    fit <- ciftest_maxcombo_statistic(c(t, 1, -t, -1), correlation, alternative)
    testthat::expect_lt(abs(fit$p.value / expected - 1), 1e-7)
  }
})
