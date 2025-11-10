test_that("getPotentialRisk uses closed form for binary RR", {
  n <- 4
  x_l <- cbind(1, c(-0.2, 0.1, 0.3, -0.4))
  exposure <- factor(c(0, 1, 0, 1))
  design <- cifmodeling:::reg_read_exposure_design(
    data.frame(a = exposure),
    exposure = "a",
    code.exposure.ref = 0
  )
  x_a <- design$x_a
  iv <- cifmodeling:::reg_index_for_parameter(NA, x_l, x_a, length.time.point = 1L)
  estimand <- list(
    exposure.levels = design$exposure.levels,
    effect.measure1 = "RR",
    index.vector = iv
  )
  alpha_beta <- c(0.1, -0.3, 0.2)
  offset <- rep(0, n)
  optim.method <- list()
  prob.bound <- 1e-5

  risk <- cifmodeling:::getPotentialRisk(alpha_beta, x_a, x_l, offset, estimand, optim.method, prob.bound)
  expect_identical(attr(risk, "solver"), "closed-form")
  expect_equal(dim(risk), c(n, 2))
})

test_that("LM path handles multi-level exposure", {
  set.seed(123)
  n <- 6
  df <- data.frame(a = factor(rep(c("a0", "a1", "a2"), length.out = n)))
  design <- cifmodeling:::reg_read_exposure_design(df, exposure = "a", code.exposure.ref = "a0")
  x_a <- design$x_a
  x_l <- cbind(1, stats::rnorm(n))
  iv <- cifmodeling:::reg_index_for_parameter(NA, x_l, x_a, length.time.point = 1L)
  estimand <- list(
    exposure.levels = design$exposure.levels,
    effect.measure1 = "RR",
    index.vector = iv
  )
  alpha_beta <- c(-0.2, 0.15, log(1.1), log(0.9))
  offset <- rep(0, n)
  optim.method <- list()
  prob.bound <- 1e-5

  risk <- cifmodeling:::getPotentialRisk(alpha_beta, x_a, x_l, offset, estimand, optim.method, prob.bound)
  expect_identical(attr(risk, "solver"), "LM")
  expect_equal(dim(risk), c(n, 3))
  ey <- cifmodeling:::calculateEY1(risk, x_a)
  expect_true(all(is.finite(ey)))
  expect_true(all(ey > 0 & ey < 1))
})

test_that("PROPORTIONAL-SURVIVAL style loop yields finite EY", {
  n <- 5
  exposure <- factor(c(0, 1, 0, 1, 0))
  design <- cifmodeling:::reg_read_exposure_design(
    data.frame(a = exposure),
    exposure = "a",
    code.exposure.ref = 0
  )
  x_a <- design$x_a
  x_l <- matrix(1, nrow = n, ncol = 1)
  time.point <- c(1, 2)
  iv <- cifmodeling:::reg_index_for_parameter(NA, x_l, x_a, length.time.point = length(time.point))
  estimand <- list(
    exposure.levels = design$exposure.levels,
    effect.measure1 = "SHR",
    index.vector = iv
  )
  alpha_beta <- c(-0.1, -0.05, log(1.2))
  offset <- rep(0, n)
  optim.method <- list()
  prob.bound <- 1e-5

  alpha_beta_i <- rep(NA_real_, iv[7])
  for (i_time in seq_along(time.point)) {
    i_para <- iv[1] * (i_time - 1L) + 1L
    idx_alpha_src <- seq.int(i_para, i_para + iv[1] - 1L)
    alpha_beta_i[seq_len(iv[1])] <- alpha_beta[idx_alpha_src]
    idx_beta <- seq.int(iv[2], iv[3])
    alpha_beta_i[idx_beta] <- alpha_beta[iv[8] / 2]

    risk <- cifmodeling:::getPotentialRisk(alpha_beta_i, x_a, x_l, offset, estimand, optim.method, prob.bound)
    ey <- cifmodeling:::calculateEY1(risk, x_a)
    expect_true(all(is.finite(ey)))
    expect_true(all(ey > 0 & ey < 1))
  }
})

test_that("calculateCovSurvival names exposures when multiple levels", {
  n <- 4
  x_l <- cbind(1, stats::rnorm(n))
  exposure <- factor(c(0, 1, 2, 1))
  design <- cifmodeling:::reg_read_exposure_design(
    data.frame(a = exposure),
    exposure = "a",
    code.exposure.ref = 0
  )
  x_a <- design$x_a
  iv <- cifmodeling:::reg_index_for_parameter(NA, x_l, x_a, length.time.point = 1L)
  estimand <- list(
    index.vector = iv,
    exposure.levels = design$exposure.levels,
    effect.measure1 = "RR"
  )
  set.seed(321)
  potential <- matrix(stats::runif(n * design$exposure.levels, 0.1, 0.4), nrow = n)
  obj <- list(
    score = matrix(stats::runif(n * iv[3], 0.01, 0.2), nrow = n, ncol = iv[3]),
    ey_1 = rep(0.2, n),
    w11 = rep(5, n),
    t = seq_len(n),
    y_0_ = rep(c(1, 0), length.out = n),
    y_1 = rep(0, n),
    x_a = x_a,
    x_l = x_l,
    potential.CIFs = potential
  )
  boot.method <- list(report.sandwich.conf = FALSE, report.boot.conf = FALSE)
  out <- cifmodeling:::calculateCovSurvival(obj, estimand, boot.method, prob.bound = 1e-5)
  expect_equal(
    tail(colnames(out$influence.function), design$exposure.levels - 1L),
    paste0("exposure", seq_len(design$exposure.levels - 1L))
  )
})
