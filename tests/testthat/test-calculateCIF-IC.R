test_that("the CIF bounds respect open left endpoints, exact events and the tail", {
  fit <- calculateCIF_IC(
    left = c(3, 6, 5, 9), right = c(6, 7, 7, Inf),
    epsilon = c(1, 1, 2, 0), times = c(3, 5, 6, 7, 9, 10)
  )
  expect_true(fit$diagnostics$converged)
  expect_equal(fit$cif.lower[, 1], c(0, 0, .25, .5, .5, .5))
  expect_equal(fit$cif.upper[, 1], c(0, .25, .25, .5, .5, .75))
  expect_equal(fit$cif.lower[, 2], c(0, 0, 0, .25, .25, .25))
  expect_equal(fit$cif.upper[, 2], c(0, 0, .25, .25, .25, .5))
  expect_equal(fit$cif[, 1], c(0, NA, .25, .5, .5, NA))
  expect_equal(fit$cif[, 2], c(0, 0, NA, .25, .25, NA))

  # An exact event at 2 and an interval starting at 2 cannot share its mass.
  mixed <- calculateCIF_IC(c(2, 2, 1), c(2, 3, Inf), c(1, 1, 0),
                           times = c(1.5, 2, 2.5, 3))
  expect_equal(mixed$cif.lower[, 1], c(0, .5, .5, 1))
  expect_equal(mixed$cif.upper[, 1], c(0, .5, 1, 1))
  expect_equal(mixed$cif[, 2], rep(0, 4))
  expect_equal(mixed$diagnostics$loglik, 2 * log(.5))
})

test_that("unresolved right-censored mass remains beyond observed follow-up", {
  fit <- calculateCIF_IC(c(0, 4), c(Inf, Inf), c(0, 0), times = c(0, 4, 5))
  expect_equal(fit$cif.lower, matrix(0, 3, 2, dimnames = list(NULL, c("1", "2"))))
  expect_equal(unname(fit$cif.upper[3, ]), c(1, 1))
  expect_true(all(is.na(fit$cif[3, ])))
  expect_equal(fit$diagnostics$loglik, 0)
  expect_error(calculateCIF_IC(1, Inf, 1), "codes")
  expect_error(calculateCIF_IC(0, 1, 1, w = 0), "positive weight")
  expect_warning(
    calculateCIF_IC(c(0, 0, 2), c(1, 2, Inf), c(1, 1, 0), max.iter = 1),
    "did not converge"
  )
})

test_that("exact times with right censoring recover the current AJ curves", {
  time <- c(4, 0, 2, 2, 3, 2, 1, 3, 5, 5, 1, 1, 4)
  epsilon <- c(0, 1, 1, 2, 0, 0, 0, 2, 1, 0, 2, 2, 0)
  strata <- factor(rep(c("B", "A", "C"), c(6, 5, 2)), levels = c("A", "B", "C"))
  w <- c(.5, 1, 2, 1.2, .7, 1, 1.5, 2, .75, 3, 0, 1, 0)
  for (weight in list(NULL, w)) {
    fit <- calculateCIF_IC(time, ifelse(epsilon == 0, Inf, time), epsilon,
                           w = weight, strata = strata)
    reference <- lapply(1:2, function(cause) {
      keep <- if (is.null(weight)) rep(TRUE, length(time)) else weight > 0
      event <- ifelse(epsilon == 0, 0L, ifelse(epsilon == cause, 1L, 2L))
      calculateAJ_Rcpp(time[keep], event[keep], w = weight[keep], strata = strata[keep],
                       conf_type = "none", return_if = FALSE, competing_risk = TRUE)
    })
    expect_equal(fit$time, reference[[1]]$time)
    expect_equal(fit$cif.lower, fit$cif.upper)
    expect_equal(unname(fit$cif), do.call(cbind, lapply(reference, `[[`, "aj")),
                 tolerance = 1e-14)
  }
})
