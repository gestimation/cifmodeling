test_that("cifsimulate returns analysis-ready competing-risks data", {
  set.seed(20260815)
  dat <- cifsimulate(500, effect = "early", censor.rate = 0.2)

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 500L)
  expect_true(all(c(
    "id", "time", "status", "A", "L", "event_time", "event_status",
    "censor_time", "random_censored", "p1", "p2", "p1_analysis",
    "p2_analysis"
  ) %in% names(dat)))
  expect_setequal(unique(dat$status), intersect(0:2, unique(dat$status)))
  expect_identical(levels(dat$A), c("0", "1"))
  expect_identical(levels(dat$L), c("0", "1"))
  expect_true(all(dat$p1 + dat$p2 <= 1 + 1e-12))
  expect_true(all(dat$time <= 6))

  truth <- attr(dat, "truth")
  expect_named(truth, c("curves", "marginal", "effect_process", "summary"))
  expect_true(min(truth$curves$survival) >= -1e-12)
  expect_equal(truth$summary$effect, "early")
})

test_that("all time-varying effects remain on one side of the null", {
  effects <- c("constant", "early", "late", "fh-early", "fh-late")
  for (effect in effects) {
    dat <- cifsimulate(20, effect = effect, shr = 0.65)
    process <- attr(dat, "truth")$effect_process
    expect_true(all(process$shr_A1_vs_A0 >= 0.65 - 1e-12))
    expect_true(all(process$shr_A1_vs_A0 <= 1 + 1e-12))
  }

  early <- attr(cifsimulate(20, effect = "early"), "truth")
  early_change <- early$summary$change.time
  expect_equal(
    early$effect_process$shr_A1_vs_A0[
      early$effect_process$start < early_change
    ],
    rep(0.65, sum(early$effect_process$start < early_change)),
    tolerance = 1e-12
  )
  expect_equal(
    early$effect_process$shr_A1_vs_A0[
      early$effect_process$start >= early_change
    ],
    rep(1, sum(early$effect_process$start >= early_change)),
    tolerance = 1e-12
  )
})

test_that("null and constant truths agree with their closed forms", {
  null <- attr(cifsimulate(20, shr = 1), "truth")$marginal
  expect_equal(
    null$cif1[null$A == 0], null$cif1[null$A == 1],
    tolerance = 1e-12, ignore_attr = TRUE
  )

  constant <- attr(cifsimulate(20, effect = "constant"), "truth")$curves
  cell <- subset(constant, A == 1 & L == 0)
  expected <- 1 - exp(-0.2 * (1 - exp(-cell$time)) * 0.65)
  expect_equal(cell$cif1, expected, tolerance = 1e-12)
})

test_that("supplied random inputs give exact common-random-number runs", {
  set.seed(20260816)
  n <- 200L
  A <- stats::rbinom(n, 1, 0.5)
  L <- stats::rbinom(n, 1, 0.5)
  uniforms <- list(
    cause1 = stats::runif(n),
    cause2 = stats::runif(n),
    event = stats::runif(n),
    censor = stats::runif(n)
  )
  first <- cifsimulate(
    n, effect = "late", A = A, L = L, uniforms = uniforms,
    censor.rate = 0.2
  )
  set.seed(1)
  second <- cifsimulate(
    n, effect = "late", A = A, L = L, uniforms = uniforms,
    censor.rate = 0.2
  )
  expect_identical(first, second)
})

test_that("empirical CIFs reproduce the returned population truth", {
  set.seed(20260817)
  dat <- cifsimulate(12000, effect = "early")
  endpoint <- attr(dat, "truth")$summary$endpoint
  empirical1 <- tapply(
    dat$event_status == 1L & dat$event_time <= 1, dat$A, mean
  )
  empirical2 <- tapply(
    dat$event_status == 2L & dat$event_time <= 1, dat$A, mean
  )
  expect_lt(max(abs(empirical1 - endpoint$cif1)), 0.025)
  expect_lt(max(abs(empirical2 - endpoint$cif2)), 0.025)
})

test_that("cifsimulate output can be analyzed directly by ciftest", {
  set.seed(20260818)
  dat <- cifsimulate(400, effect = "constant")
  fit <- ciftest(
    Event(time, status) ~ A,
    data = dat,
    augmentation = FALSE,
    tau = 1
  )
  expect_s3_class(fit, "ciftest")
  expect_true(is.finite(fit$p.value))
})

test_that("cifsimulate rejects invalid inputs", {
  expect_error(cifsimulate(0), "positive integer")
  expect_error(cifsimulate(10.5), "positive integer")
  expect_error(cifsimulate(10, shr = 0), "positive finite")
  expect_error(cifsimulate(10, analysis.tau = 7), "analysis.tau")
  expect_error(cifsimulate(10, A = rep(2, 10)), "binary 0/1")
  expect_error(
    cifsimulate(10, uniforms = list(bad = stats::runif(10))),
    "may contain"
  )
})
