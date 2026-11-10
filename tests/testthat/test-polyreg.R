test_that("polyreg estimates competing-risk effects for categorical and stratified binary exposure", {
  skip_on_cran()
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', code.exposure.ref='Q1', time.point=8, outcome.type='C')
  expect_equal(round(unname(coef(output)), 3),
    c(-1.383, -0.276, -0.067, -0.555, -3.991, -0.104, 0.045, -0.169))
  expect_true(all(is.finite(vcov(output))))

  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, strata = 'strata', effect.measure1='SHR', effect.measure2='SHR', time.point=8, outcome.type='C')
  expect_equal(round(unname(coef(output)), 3), c(-1.383, 0.363, -3.988, 0.081))
})

test_that("fixed-time polyreg covariance matches perturbations of the censoring KM", {
  d <- data.frame(
    time = c(1, 1, 2, 3, 4, 4, 5, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 6),
    status = c(0, 1, 0, 2, 1, 0, 2, 1, 0, 0, 2, 1, 0, 1, 2, 0, 1, 0),
    group = rep(0:1, each = 9)
  )
  d$stratum <- factor(d$group)
  horizon <- 4
  n <- nrow(d)
  weight <- rep(1, n)

  # In a saturated fit the risk in each group is the mean of the weighted
  # event indicators. Perturb the case weights in that mean and in the KM.
  risk <- function(w, cause, status, stratum) {
    vapply(0:1, function(group) {
      events <- which(d$group == group & status == cause & d$time <= horizon)
      weighted_events <- vapply(events, function(i) {
        same_stratum <- stratum == stratum[i]
        censored <- same_stratum & status == 0
        jumps <- unique(d$time[censored & d$time < d$time[i]])
        survival <- prod(vapply(jumps, function(u) {
          1 - sum(w[censored & d$time == u]) / sum(w[same_stratum & d$time >= u])
        }, numeric(1)))
        w[i] / survival
      }, numeric(1))
      sum(weighted_events) / sum(w[d$group == group])
    }, numeric(1))
  }
  influence <- function(estimand) {
    vapply(seq_len(n), function(i) {
      delta <- numeric(n)
      delta[i] <- 1e-5
      (estimand(weight + delta) - estimand(weight - delta)) / (2e-5)
    }, numeric(length(estimand(weight))))
  }

  d$failure <- as.integer(d$status != 0)
  log_rr <- function(w) diff(log(risk(w, 1, d$failure, d$stratum)))
  fit <- polyreg(Event(time, failure) ~ 1, "group", data = d, strata = "stratum",
    time.point = horizon, outcome.type = "survival")
  expect_equal(unname(coef(fit)[2]), log_rr(weight), tolerance = 1e-6)
  expect_equal(unname(vcov(fit)[2, 2]), sum(influence(log_rr)^2), tolerance = 1e-5)

  for (stratified in c(TRUE, FALSE)) {
    stratum <- if (stratified) d$stratum else rep(1, n)
    log_effects <- function(w) c(
      diff(log(risk(w, 1, d$status, stratum))),
      diff(log(-log1p(-risk(w, 2, d$status, stratum))))
    )
    fit <- polyreg(Event(time, status) ~ 1, "group", data = d,
      strata = if (stratified) "stratum" else NULL, time.point = horizon,
      effect.measure1 = "RR", effect.measure2 = "SHR")
    contrasts <- c(2, 4)
    expect_equal(unname(coef(fit)[contrasts]), log_effects(weight), tolerance = 1e-5)
    expected_if <- influence(log_effects)
    expect_equal(unname(vcov(fit)[contrasts, contrasts]), tcrossprod(expected_if),
      tolerance = 1e-4)
    expect_equal(unname(fit$diagnostics$influence.function[, contrasts]) / n,
      t(expected_if), tolerance = 1e-4)
  }
})

test_that("polyreg() produced expected coefficients and variance covariance matrix from survival data with missing variables", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  expected_df <- diabetes.complications[-(1:10), ]
  expected_output <- polyreg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = expected_df, effect.measure1='RR', time.point=8, outcome.type='S')
  expected <- round(expected_output$coef,digits=3)

  diabetes.complications$t[1:2] <- NA
  diabetes.complications$d[3:4] <- NA
  diabetes.complications$sex[5:6] <- NA
  diabetes.complications$fruitq1[7:8] <- NA
  diabetes.complications$strata[9:10] <- NA
  output <- polyreg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested <- round(output$coef,digits=3)
  expect_equal(tested, expected)
})

test_that("polyreg() produced expected common effects from competing risks data", {
  skip_on_cran()
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=1:8, outcome.type='PC', report.boot.conf=FALSE)
  tested <- round(output$coef,digits=3)
  expected <- c(-7.225, -4.051, -3.067, -2.534, -2.114, -1.803, -1.572, -1.392, 0.296, -9.319, -7.712, -6.977, -6.087, -5.156, -4.731, -4.366, -4.040, -0.022)
  expect_equal(tested, expected)
})

test_that("polyreg() produced expected SE for common effects from survival data", {
  df_test <- createTestData1(100, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  df_test$d <- as.numeric(df_test$epsilon>0)
  df_test$a <- as.numeric(df_test$strata)
  output <- polyreg(nuisance.model = Event(t,d) ~ +1, exposure = 'a', data = df_test,
                    effect.measure1='RR', time.point=20, outcome.type='PS', report.boot.conf=TRUE, boot.replications=10, boot.seed=46)
  tested <- round(output$summary$`event 1 (no competing risk`$tidy$std.error,digits=3)
  expected <- c(0.194)
  expect_equal(tested, expected)
})

#test_that("polyreg() produced expected SE for common effects from competing-risks data", {
#  df_test <- createTestData1(100, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
#  df_test$a <- as.numeric(df_test$strata)
#  output <- polyreg(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'a', data = df_test,
#                    effect.measure1='RR', effect.measure2='RR', time.point=20, outcome.type='PROPORTIONAL-COMPETING-RISK', report.boot.conf=TRUE, `boot.replications`=10, boot.seed=46)
#  tested <- round(output$summary$event2$tidy$std.error,digits=3)
#  expected <- c(4.886)
#  expect_equal(tested, expected)
#})

test_that("tidy/glance methods for polyreg are registered", {
  # generics パッケージの generic を対象に S3 メソッド一覧を取る
  ms_tidy   <- methods(generics::tidy)
  ms_glance <- methods(generics::glance)

  expect_true("tidy.polyreg"   %in% ms_tidy)
  expect_true("glance.polyreg" %in% ms_glance)

  # getS3method: generic は generics 名前空間で探す
  tidy_m <- utils::getS3method(
    f     = "tidy",
    class = "polyreg",
    envir = asNamespace("generics"),
    optional = TRUE
  )
  glance_m <- utils::getS3method(
    f     = "glance",
    class = "polyreg",
    envir = asNamespace("generics"),
    optional = TRUE
  )

  expect_true(is.function(tidy_m))
  expect_true(is.function(glance_m))
})

test_that("tidy.polyreg() and glance.polyreg() work for survival models", {
  skip_on_cran()

  data(diabetes.complications)
  diabetes.complications$d <- as.integer(diabetes.complications$epsilon>0)
  diabetes.complications$fruitq1 <- ifelse(diabetes.complications$fruitq == "Q1",1,0)

  fit <- polyreg(
    nuisance.model = Event(t, d) ~ 1,
    exposure       = "fruitq1",
    data           = diabetes.complications,
    effect.measure1 = "RR",
    time.point      = 8,
    outcome.type    = "survival",
    report.optim.convergence = FALSE
  )

  expect_s3_class(fit, "polyreg")

  td <- generics::tidy(fit)
  expect_true(is.data.frame(td))
  expect_gt(nrow(td), 0L)

  expect_true(all(c("term", "estimate", "std.error",
                    "p.value", "conf.low", "conf.high") %in% names(td)))
  expect_true("event" %in% names(td))
  expect_setequal(unique(td$event), "event1")

  gl <- generics::glance(fit)
  expect_true(is.data.frame(gl))
  expect_gte(nrow(gl), 1L)
  expect_true(all(c("effect.measure", "n.events") %in% names(gl)))
})

test_that("tidy/glance work for competing-risk polyreg", {
  skip_on_cran()

  data(diabetes.complications)
  diabetes.complications$fruitq1 <- ifelse(diabetes.complications$fruitq == "Q1",1,0)

  fit <- polyreg(
    nuisance.model = Event(t, epsilon) ~ 1,
    exposure       = "fruitq1",
    data           = diabetes.complications,
    effect.measure1 = "RR",
    effect.measure2 = "RR",
    time.point      = 8,
    outcome.type    = "competing-risk",
    report.optim.convergence = FALSE
  )

  td <- generics::tidy(fit)
  expect_true(is.data.frame(td))
  expect_gt(nrow(td), 0L)
  expect_true(all(c("term", "estimate", "std.error",
                    "p.value", "conf.low", "conf.high") %in% names(td)))
  expect_true("event" %in% names(td))
  expect_true("event1" %in% unique(td$event))  # まずはここだけ保証

  gl <- generics::glance(fit)
  expect_true(is.data.frame(gl))
  expect_gt(nrow(gl), 0L)
  expect_true(all(c("effect.measure", "n.events") %in% names(gl)))
})


test_that("augment.polyreg() returns row-wise diagnostics", {
  skip_on_cran()
  data(diabetes.complications)
  diabetes.complications$fruitq1 <- ifelse(
    diabetes.complications$fruitq == "Q1",1,0
  )

  fit <- polyreg(
    nuisance.model = Event(t, epsilon) ~ 1,
    exposure       = "fruitq1",
    data           = diabetes.complications,
    effect.measure1 = "RR",
    time.point      = 8,
    outcome.type    = "competing-risk",
    report.optim.convergence = FALSE
  )

  au <- generics::augment(fit)
  expect_true(is.data.frame(au))
  expect_equal(nrow(au), nrow(diabetes.complications))
  expect_true(".weights" %in% names(au))
  expect_true(any(grepl("^\\.fitted_", names(au))))
})

test_that("coef(), vcov(), nobs() work for survival polyreg models", {
  skip_on_cran()

  data(diabetes.complications)
  diabetes.complications$d <- as.integer(diabetes.complications$epsilon > 0)
  diabetes.complications$fruitq1 <- ifelse(diabetes.complications$fruitq == "Q1",1,0)

  fit <- polyreg(
    nuisance.model = Event(t, d) ~ 1,
    exposure       = "fruitq1",
    data           = diabetes.complications,
    effect.measure1 = "RR",
    time.point      = 8,
    outcome.type    = "survival",
    report.optim.convergence = FALSE
  )

  co <- coef(fit)
  expect_true(is.numeric(co))
  expect_gt(length(co), 0L)

  expect_equal(nobs(fit), nrow(diabetes.complications))

  V_def  <- vcov(fit)
  V_sand <- vcov(fit, type = "sandwich")

  expect_true(is.matrix(V_def))
  expect_equal(dim(V_def), c(length(co), length(co)))
  expect_identical(V_def, V_sand)

  expect_warning(
    V_boot <- vcov(fit, type = "bootstrap"),
    "Bootstrap variance is not available; returning sandwich variance instead.",
    fixed = TRUE
  )
  expect_identical(V_boot, V_sand)
})

test_that("vcov.polyreg() chooses bootstrap variance when available", {

  fake <- list(
    coef           = c(beta1 = 0.1, beta2 = -0.2),
    vcov           = NULL,
    vcov_bootstrap = matrix(c(0.04, 0.01,
                              0.01, 0.09), nrow = 2),
    outcome.type   = "competing-risk",
    boot.method    = list(
      report.sandwich.conf = FALSE,
      report.boot.conf     = TRUE
    )
  )
  class(fake) <- "polyreg"

  V_def <- vcov(fake, type = "default")
  expect_identical(V_def, fake$vcov_bootstrap)

  expect_warning(
    V_sand <- vcov(fake, type = "sandwich"),
    "Sandwich variance is not available; returning bootstrap variance instead.",
    fixed = TRUE
  )
  expect_identical(V_sand, fake$vcov_bootstrap)

  V_boot <- vcov(fake, type = "bootstrap")
  expect_identical(V_boot, fake$vcov_bootstrap)
})
