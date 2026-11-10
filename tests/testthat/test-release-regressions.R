test_that("multilevel binomial polyreg reports the fitted risk ratios", {
  group <- rep(LETTERS[1:4], each = 100L)
  cases <- c(10L, 20L, 30L, 40L)
  status <- as.integer(rep(seq_len(100L), 4L) <= rep(cases, each = 100L))
  d <- data.frame(time = 1, status, group)

  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    outcome.type = "binomial", code.exposure.ref = "A", effect.measure1 = "RR")
  reported <- generics::tidy(fit, exponentiate = TRUE)
  expected <- cases[-1L] / cases[1L]
  terms <- paste0("group, ", LETTERS[2:4], " vs A")
  risks <- vapply(LETTERS[1:4], function(level) {
    as_prediction_data(fit, exposure.level = level)$estimate[1L]
  }, numeric(1))

  expect_equal(unname(risks), cases / 100, tolerance = 1e-5)
  expect_equal(reported$estimate[match(terms, reported$term)], expected, tolerance = 1e-5)
  probability <- cases / 100
  log_risk_variance <- (1 - probability) / cases
  expect_equal(unname(vcov(fit)[2:4, 2:4]),
    diag(log_risk_variance[-1]) + log_risk_variance[1], tolerance = 1e-5)
})

test_that("multilevel binomial polyreg supports categorical covariates and a different reference", {
  d <- expand.grid(row = 1:50, group = c("arm_A", "arm_B", "arm_C"),
    site = c("north", "south"))
  rates <- c(.1, .2, .3) * 50
  d$status <- as.integer(d$row <= rates[match(d$group, levels(d$group))] *
    ifelse(d$site == "north", 1, 2))
  d$time <- 1
  fit <- polyreg(Event(time, status) ~ site, exposure = "group", data = d,
    outcome.type = "binomial", code.exposure.ref = "arm_B")
  result <- generics::tidy(fit, exponentiate = TRUE)
  expect_equal(result$term, c("group, arm_A vs arm_B", "group, arm_C vs arm_B"))
  expect_equal(result$estimate, c(.5, 1.5), tolerance = 1e-5)
  expect_true(all(is.finite(result$std.error) & result$std.error > 0))
})

test_that("polyreg IPCW uses the specified censoring code", {
  d <- data.frame(time = rep(1:8, 2L), group = rep(0:1, each = 8L),
    status = c(1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0))
  d$reversed <- 1 - d$status
  d$custom <- ifelse(d$status == 1, 7, 9)
  original_allowed <- getOption("cifmodeling.allowed")
  default <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    outcome.type = "survival", time.point = 4.5)
  reversed <- polyreg(Event(time, reversed) ~ 1, exposure = "group", data = d,
    outcome.type = "survival", time.point = 4.5, code.event1 = 0, code.censoring = 1)
  custom <- polyreg(cifmodeling::Event(time, custom) ~ 1, exposure = "group", data = d,
    outcome.type = "survival", time.point = 4.5, code.event1 = 7, code.censoring = 9)

  # Censoring KM jumps at 1, 2 and 4 with risk sets of size 16, 14 and 10.
  g_after_two <- (15 / 16) * (12 / 14)
  g_at_horizon <- g_after_two * (9 / 10)
  expected <- ifelse(d$time > 4.5, 1 / g_at_horizon,
    ifelse(d$status == 0, 0, ifelse(d$time == 1, 1, 1 / g_after_two)))
  expect_equal(as.numeric(default$diagnostics$ip.weight), expected)
  expect_equal(as.numeric(reversed$diagnostics$ip.weight), expected)
  expect_equal(as.numeric(custom$diagnostics$ip.weight), expected)
  expect_equal(generics::tidy(reversed), generics::tidy(default), tolerance = 1e-6)
  expect_equal(generics::tidy(custom), generics::tidy(default), tolerance = 1e-6)
  expect_identical(getOption("cifmodeling.allowed"), original_allowed)
})

test_that("polyreg IPCW preserves character censoring strata", {
  d <- data.frame(time = rep(1:8, 2L), group = rep(0:1, each = 8L),
    status = c(1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0))
  d$stratum <- rep(c("A", "B"), each = 8L)
  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    strata = "stratum", outcome.type = "survival", time.point = 4.5)
  factor_fit <- polyreg(Event(time, status) ~ 1, exposure = "group",
    data = transform(d, stratum = factor(stratum)), strata = "stratum",
    outcome.type = "survival", time.point = 4.5)

  # At the horizon, G_A = (6/7)(4/5) and G_B = (7/8)(6/7).
  expected <- c(1, 0, 7 / 6, 0, rep(35 / 24, 4), 0, 0, rep(4 / 3, 6))
  expect_equal(as.numeric(fit$diagnostics$ip.weight), expected)
  expect_equal(generics::tidy(fit), generics::tidy(factor_fit))
})

test_that("joint Gray covariance is invariant to event recoding", {
  d <- data.frame(time = rep(1:8, 3L), group = rep(LETTERS[1:3], each = 8L),
    status = c(1, 0, 2, 1, 0, 2, 0, 1, 0, 1, 2, 0, 1, 0, 2, 2,
      2, 0, 1, 2, 1, 0, 1, 0))
  d$reversed <- ifelse(d$status == 2, 2, 1 - d$status)
  directions <- c("unweighted", "early", "late")
  default <- ciftest_mdir(Event(time, status) ~ group, d,
    outcome.type = "competing-risk", test = "gray", directions = directions)
  reversed <- ciftest_mdir(Event(time, reversed) ~ group, d,
    outcome.type = "competing-risk", test = "gray", directions = directions,
    code.event1 = 0, code.censoring = 1)

  reference <- gray_reference_slow(d$time, d$status, d$group)
  expect_equal(unname(default$vcov.score[1:2, 1:2]), reference$var)
  expect_equal(reversed$vcov.score, default$vcov.score)
  expect_equal(reversed$p.value, default$p.value)
})

test_that("CIF standard errors depend only on events through the reported time", {
  d <- data.frame(time = 1:6, status = c(0, 1, 0, 2, 1, 0))
  recoded <- transform(d, status = ifelse(time == 5, 0, status))
  # At time 2 one of five at risk fails. Through time 4, CIF stays at 1/5.
  # The documented variance formulas reduce to their second term there.
  expected_se <- c(delta = sqrt(4 / 125), aalen = sqrt(4 / (25 * 4)))
  # At time 5, the earlier jumps enter through F(5) - F(2) = F(5) - F(4).
  difference <- 4 / 15
  survival_left <- 8 / 15
  terminal_var <- c(
    delta = difference^2 * (1 / 20 + 1 / 6) + 4 / 125 + survival_left^2 / 8 -
      2 * difference / 25,
    aalen = difference^2 * (1 / 16 + 1 / 4) + 1 / 25 + survival_left^2 / 4 -
      2 * difference / 20
  )
  for (method in names(expected_se)) {
    full <- cifcurve(Event(time, status) ~ 1, d, outcome.type = "competing-risk",
      conf.type = "plain", error = method)
    without_later_event <- cifcurve(Event(time, status) ~ 1, recoded,
      outcome.type = "competing-risk", conf.type = "plain", error = method)
    expect_equal(1 - full$surv[1:4], c(0, rep(1 / 5, 3)), info = method)
    expect_equal(full$std.err[1:4], c(0, rep(expected_se[[method]], 3)), info = method)
    expect_equal(full$std.err[5:6], rep(sqrt(terminal_var[[method]]), 2), info = method)
    expect_equal(full$lower[1:4], without_later_event$lower[1:4], info = method)
    expect_equal(full$upper[1:4], without_later_event$upper[1:4], info = method)
  }
})

test_that("legacy competing-risk curves respect constant case weights", {
  d <- data.frame(time = 1:6, status = c(0, 1, 0, 2, 1, 0), weight = 3)
  fit <- cifcurve(Event(time, status) ~ 1, d, outcome.type = "competing-risk",
    weights = "weight", engine = "calculateKM", n.risk.type = "unweighted", conf.type = "plain")
  # The AJ jumps are 1/5 at time 2 and (4/5) * (2/3) * (1/2) at time 5.
  expect_equal(1 - fit$surv, c(0, rep(1 / 5, 3), rep(7 / 15, 2)))
  expect_equal(fit$n.risk, 6:1)
})

test_that("cifcurve rejects Surv responses that are not right-censored", {
  skip_if_not_installed("survival")
  d <- data.frame(start = 0, stop = 1, status = c(1, 0, 1, 0))
  expect_error(cifcurve(survival::Surv(start, stop, status) ~ 1, d,
    outcome.type = "survival"), "right.censored")
  expect_error(cifcurve(survival::Surv(stop, status, type = "left") ~ 1, d,
    outcome.type = "survival"), "right.censored")
})

test_that("proportional polyreg normal bootstrap honors the confidence level", {
  d <- data.frame(group = rep(0:1, each = 24L),
    status = c(rep(1L, 8L), rep(0L, 16L), rep(1L, 14L), rep(0L, 10L)))
  d$time <- ifelse(d$status == 1L, 1, 2)
  for (level in c(.80, .90, .95)) {
    fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
      outcome.type = "proportional-survival", time.point = 1, conf.int = level,
      report.boot.conf = TRUE, boot.replications = 20, boot.seed = 46)
    result <- generics::tidy(fit)
    half_width <- (result$conf.high - result$conf.low) / 2
    expect_true(is.finite(result$std.error) && result$std.error > 0)
    expect_equal(half_width, stats::qnorm((1 + level) / 2) * result$std.error,
      info = paste("confidence level", level))
  }

  changed <- generics::tidy(fit, conf.level = .90)
  expected_limits <- result$estimate + c(-1, 1) * stats::qnorm(.95) * result$std.error
  expect_equal(c(changed$conf.low, changed$conf.high), expected_limits)
})

test_that("proportional polyreg BCa selects the exposure coefficient", {
  d <- data.frame(group = rep(0:1, each = 24L),
    status = c(rep(1L, 8L), rep(0L, 16L), rep(1L, 14L), rep(0L, 10L)))
  d$time <- ifelse(d$status == 1L, 1, 2)
  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    outcome.type = "proportional-survival", time.point = 1, conf.int = .90,
    report.boot.conf = TRUE, boot.bca = TRUE, boot.replications = 199, boot.seed = 46)
  result <- generics::tidy(fit)
  expect_true(all(is.finite(c(result$conf.low, result$conf.high))))
  expect_lt(result$conf.low, result$conf.high)
  # With one time point and no earlier censoring, the target is log((14/24)/(8/24)).
  expect_lt(result$conf.low, log(14 / 8))
  expect_gt(result$conf.high, log(14 / 8))
})

test_that("proportional bootstrap follows the coefficient layout across times and causes", {
  d <- do.call(rbind, lapply(0:1, function(group) {
    n1 <- c(6, 9)[group + 1]
    n2 <- c(8, 6)[group + 1]
    counts <- c(n1, n2, n1, n2, 60 - 2 * (n1 + n2))
    data.frame(group = group, time = rep(c(1, 1, 2, 2, 3), counts),
      status = rep(c(1, 2, 1, 2, 0), counts))
  }))
  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    outcome.type = "proportional-competing-risk", time.point = 1:2,
    report.nuisance.parameter = TRUE, conf.int = .80, boot.replications = 20, boot.seed = 49)
  rows <- c(3, 6)
  expect_equal(coef(fit)[rows], log(c(1.5, .75)), tolerance = 1e-5)
  result <- generics::tidy(fit, event = "both")
  expect_identical(result$term, rep(c("Intercept at 1", "Intercept at 2", "group, 1 vs 0"), 2))
  expect_identical(result$event[rows], c("event1", "event2"))
  expect_equal(result$std.error, sqrt(diag(vcov(fit, variance = "bootstrap"))))
  expect_equal((result$conf.high - result$conf.low) / 2,
    qnorm(.90) * result$std.error)
})
