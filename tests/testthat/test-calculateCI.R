test_that("Greenwood standard error of cifcurve() yields the same outputs as separate analysis when strata is present", {
  testdata <- createTestData1(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata1 <- subset(testdata, strata==FALSE)
  testdata2 <- subset(testdata, strata==TRUE)
  t <- cifcurve(Event(t, d)~strata, testdata, outcome.type="S", weight="w", error = "greenwood")
  e1 <- cifcurve(Event(t, d)~1, testdata1, outcome.type="S", weight="w", error = "greenwood")
  e2 <- cifcurve(Event(t, d)~1, testdata2, outcome.type="S", weight="w", error = "greenwood")
  expected <- c(e1$std.err, e2$std.err)
  tested <- t$std.err
  expect_equal(expected, tested)
})

test_that("cifcurve() produced the same outputs as survfit() in survival in survival data", {
  testthat::skip_if_not_installed("survival")
  testdata <- createTestData1(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survival::survfit(survival::Surv(t, d)~strata, testdata, weight=w, conf.type = "none")
  t <- cifcurve(Event(t, d) ~ strata, data = testdata, weight="w", conf.type = "none", outcome.type = "survival", engine="calculateKM")
  #  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  #  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, e$lower, e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, t$lower, t$strata))
  expect_equal(expected, tested)
})

test_that("cifcurve() produced the same estimates as cif() in mets in competing risks data", {
  testthat::skip_if_not_installed("mets")
  library(mets)
  data(diabetes.complications)
  cif_fit <- mets::cif(Event(t,epsilon) ~ +1, data=diabetes.complications, cause=1)
  surv_fit <- cifcurve(Event(t, epsilon) ~ +1, data = diabetes.complications, conf.type = "none", outcome.type = "C", error="delta")
  index<-c(1,3,5,6,8,10,12,13,15,16,17,18)
  expected <- round(1-surv_fit$surv[index],digit=5)
  tested <- round(cif_fit$mu[1:12],digit=5)
  expect_equal(expected, tested)
})

test_that("cifcurve yields the same outputs as cif of mets", {
  skip_on_cran()
  skip_if_not_installed("mets")
  skip_if_not_installed("Rcpp")
  testdata <- createTestData1(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  library(mets)
  library(Rcpp)
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e_surv <- 1-e$mu
  e_time <- e$times
  expected <- as.numeric(c(e_surv[c(1,2,4,6,8,10,12,14)], e_time[c(1,2,4,6,8,10,12,14)]))
  t <- cifcurve(Event(t, epsilon)~1, data=testdata, conf.type = "none", outcome.type = "competing-risk")
  tested <- as.numeric(c(t$surv[c(2:9)], t$time[c(2:9)]))
  expect_equal(tested, expected)
})

test_that("cifcurve by strata yields the same outputs as subsetting", {
  skip_on_cran()
  skip_if_not_installed("mets")
  skip_if_not_installed("Rcpp")
  testdata <- createTestData1(20, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata_f <- testdata[testdata$strata==FALSE,]
  library(Rcpp)
  e <- cif(Event(t,epsilon)~1, data=testdata, cause=1)
  e <- cifcurve(Event(t, epsilon)~1, data=testdata_f, conf.type = "none", outcome.type = "competing-risk")
  expected <- c(e$surv, e$time)
  t <- cifcurve(Event(t, epsilon)~strata, data=testdata, conf.type = "none", outcome.type = "competing-risk")
  tested <- c(t$surv[c(1:5)], t$time[c(1:5)])
  expect_equal(tested, expected)
})

test_that("cifcurve() yields the same outputs as survfit() with log-log transformation", {
  skip_on_cran()
  skip_if_not_installed("survival")

  df_test <- createTestData1(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survival::survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "log-log")
  t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log-log", report.survfit.std.err = TRUE, outcome.type = "survival", engine="calculateKM")
  #t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log-log", report.survfit.std.err = TRUE, outcome.type = "survival")
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.infinite(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.infinite(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.infinite(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x<=0, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.infinite(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.infinite(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.infinite(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("cifcurve() yields the same outputs as survfit() with log transformation", {
  skip_on_cran()
  skip_if_not_installed("survival")

  df_test <- createTestData1(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survival::survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "log")
  #t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log", report.survfit.std.err = TRUE, outcome.type = "survival")
  t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "log", report.survfit.std.err = TRUE, outcome.type = "survival", engine="calculateKM")
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.infinite(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.infinite(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.infinite(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x<=0, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.infinite(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.infinite(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.infinite(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("cifcurve() yields the same outputs as survfit() with arcsine transformation", {
  skip_on_cran()
  skip_if_not_installed("survival")

  df_test <- createTestData1(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survival::survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "a")
  #t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "a", report.survfit.std.err = TRUE, outcome.type = "survival")
  t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "a", report.survfit.std.err = TRUE, outcome.type = "survival", engine="calculateKM")
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.infinite(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.infinite(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.infinite(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x<=0, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.infinite(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.infinite(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.infinite(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("cifcurve() yields the same outputs as survfit()", {
  skip_on_cran()
  skip_if_not_installed("survival")

  #df_test <- createTestData1(200, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  #e <- survival::survfit(Surv(t, d)~strata, df_test, weight=w, conf.type = "plain")
  #t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "plain", report.survfit.std.err = TRUE, outcome.type = "survival", engine="calculateAJ_Rcpp", error="greenwood")

  df_test <- createTestData1(200, 1, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE)
  e <- survival::survfit(Surv(t, d)~strata, df_test, conf.type = "plain")
  t <- cifcurve(Surv(t, d)~strata, df_test, conf.type = "plain", report.survfit.std.err = TRUE, outcome.type = "survival", error="greenwood")
  #t <- cifcurve(Surv(t, d)~strata, df_test, weight="w", conf.type = "plain", report.survfit.std.err = TRUE, outcome.type = "survival", engine="calculateKM")
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.nan(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.nan(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.nan(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(is.infinite(x), NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(is.infinite(x), NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(is.infinite(x), NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x>=1, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x>=1, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x>=1, NA, x))
  e$std.err <- sapply(e$std.err, function(x) ifelse(x<=0, NA, x))
  e$lower <- sapply(e$lower, function(x) ifelse(x<=0, NA, x))
  e$upper <- sapply(e$upper, function(x) ifelse(x<=0, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.nan(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.nan(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.nan(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(is.infinite(x), NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(is.infinite(x), NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(is.infinite(x), NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x>=1, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x>=1, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x>=1, NA, x))
  t$std.err <- sapply(t$std.err, function(x) ifelse(x<=0, NA, x))
  t$lower <- sapply(t$lower, function(x) ifelse(x<=0, NA, x))
  t$upper <- sapply(t$upper, function(x) ifelse(x<=0, NA, x))
  expected <- as.numeric(c(round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
  tested <- as.numeric(c(round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
#  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), round(e$lower,digit=5), round(e$upper,digit=5), e$strata))
#  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), round(t$lower,digit=5), round(t$upper,digit=5), t$strata))
  expect_equal(expected, tested)
})

test_that("Greenwood standard error of cifcurve() yields the same outputs as separate analysis when strata is present", {
  testdata <- createTestData1(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata1 <- subset(testdata, strata==FALSE)
  testdata2 <- subset(testdata, strata==TRUE)
  t <- cifcurve(Surv(t, d)~strata, testdata, weight="w", error = "greenwood", outcome.type = "survival")
  e1 <- cifcurve(Surv(t, d)~1, testdata1, weight="w", error = "greenwood", outcome.type = "survival")
  e2 <- cifcurve(Surv(t, d)~1, testdata2, weight="w", error = "greenwood", outcome.type = "survival")
  expected <- c(e1$std.err, e2$std.err)
  tested <- t$std.err
  expect_equal(expected, tested)
})

test_that("AJ influence function roughly matches jackknife influence", {
  skip_on_cran()

  set.seed(123)
  n  <- 40L
  t  <- rexp(n, rate = 0.1)
  epsilon <- sample(c(0L, 1L, 2L), n, replace = TRUE,
                    prob = c(0.3, 0.4, 0.3))

  fit <- calculateAJ_Rcpp(
    t        = t,
    epsilon  = epsilon,
    error    = "if",
    conf_type= "none",
    return_if= TRUE
  )

  time0 <- fit$time
  theta0 <- fit$aj
  IF_cpp <- fit$influence.function[[1L]]

  skip_if(
    is.null(IF_cpp) || length(IF_cpp) == 0L || ncol(IF_cpp) == 0L,
    "Influence function matrix is empty (return_if=FALSE or error_if=FALSE)."
  )

  expect_true(is.matrix(IF_cpp))
  expect_equal(dim(IF_cpp)[1L], n)
  K <- length(time0)
  expect_equal(dim(IF_cpp)[2L], K)

  theta_minus <- matrix(NA_real_, nrow = n, ncol = K)

  for (i in seq_len(n)) {
    fit_i <- calculateAJ_Rcpp(
      t        = t[-i],
      epsilon  = epsilon[-i],
      error    = "if",
      conf_type= "none",
      return_if= FALSE
    )
    f_step <- stats::stepfun(fit_i$time, c(0, fit_i$aj))
    theta_minus[i, ] <- f_step(time0)
  }

  theta_mat <- matrix(theta0, nrow = n, ncol = K, byrow = TRUE)
  pseudo <- n * theta_mat - (n - 1) * theta_minus
  IF_jk  <- pseudo - theta_mat
  IF_jk  <- sweep(IF_jk, 2, colMeans(IF_jk), "-")

  mid <- which(theta0 > 0.05 & theta0 < 0.95)
  skip_if(length(mid) < 3L, "Too few interior time points for correlation check.")

  cors <- vapply(mid, function(j) {
    cor(IF_cpp[, j], IF_jk[, j])
  }, numeric(1))

  expect_true(all(cors > 0.9),
              info = paste("column-wise correlations:",
                           paste(round(cors, 2), collapse = ", ")))

  se_cpp <- fit$std.err.aj
  se_jk  <- sqrt(colSums(IF_jk^2) / (n * n))

  rel_diff <- abs(se_cpp[mid] - se_jk[mid]) / pmax(se_cpp[mid], 1e-8)
  expect_lt(max(rel_diff), 0.2)
})

test_that("curve_check_error fallback matches spec", {
  expect_identical(curve_check_error(NULL, "survival", weights = NULL), "greenwood")
  expect_identical(curve_check_error(NULL, "competing-risk", weights = NULL), "delta")

  expect_identical(curve_check_error(NULL, "survival", weights = numeric(0)), "greenwood")
  expect_identical(curve_check_error(NULL, "competing-risk", weights = numeric(0)), "delta")

  expect_identical(curve_check_error(NULL, "survival", weights = 1, has_weights = TRUE), "if")
  expect_identical(curve_check_error(NULL, "competing-risk", weights = 1, has_weights = TRUE), "if")
})


# tests/testthat/test-calculateAJ_Rcpp-if.R

test_that("SE/CI are broadly consistent across error methods (outcome='s' and 'c')", {

  calc <- function(t, eps, w = NULL, error, conf_type = "plain", conf_int = 0.95, return_if = TRUE) {
    if (is.null(w)) {
      calculateAJ_Rcpp(
        t = t, epsilon = eps,
        error = error, conf_type = conf_type, conf_int = conf_int,
        return_if = return_if
      )
    } else {
      calculateAJ_Rcpp(
        t = t, epsilon = eps, w = w,
        error = error, conf_type = conf_type, conf_int = conf_int,
        return_if = return_if
      )
    }
  }

  extract_core <- function(out) {
    list(
      time = out$time,
      est  = out$surv,
      se   = out$std.err,
      low  = out$low,
      high = out$high
    )
  }

  core_mask <- function(core, min_est = 0.02, max_est = 0.98) {
    is.finite(core$est) &
      is.finite(core$se)  & core$se > 0 &
      is.finite(core$low) & is.finite(core$high) &
      core$est > min_est & core$est < max_est
  }

  expect_ci_valid <- function(core) {
    ok <- is.finite(core$est) & is.finite(core$low) & is.finite(core$high)
    expect_true(all(core$low[ok] <= core$high[ok]))
    expect_true(all(core$low[ok] >= 0 & core$high[ok] <= 1))
    expect_true(all(core$low[ok] <= core$est[ok] & core$est[ok] <= core$high[ok]))
  }

  expect_methods_close <- function(core_if, core_ref,
                                   se_median_ratio_bounds = c(0.85, 1.15),
                                   se_q05_q95_bounds      = c(0.60, 1.60),
                                   ci_max_abs_diff        = 0.06,
                                   label = "") {

    expect_equal(core_if$time, core_ref$time, info = paste0(label, " time grid differs"))
    expect_equal(core_if$est, core_ref$est, tolerance = 1e-12,
                 info = paste0(label, " point estimates differ"))
    expect_ci_valid(core_if)
    expect_ci_valid(core_ref)

    m <- core_mask(core_ref)
    expect_true(sum(m) >= 10, info = paste0(label, " too few comparable points"))

    # SE 比
    ratio <- core_if$se[m] / core_ref$se[m]
    med   <- median(ratio)
    q05   <- as.numeric(stats::quantile(ratio, 0.05, names = FALSE))
    q95   <- as.numeric(stats::quantile(ratio, 0.95, names = FALSE))

    expect_true(med >= se_median_ratio_bounds[1] && med <= se_median_ratio_bounds[2],
                info = paste0(label, " median(SE_if/SE_ref) = ", signif(med, 3)))
    expect_true(q05 >= se_q05_q95_bounds[1] && q95 <= se_q05_q95_bounds[2],
                info = paste0(label, " q05/q95(SE_if/SE_ref) = ",
                              signif(q05, 3), "/", signif(q95, 3)))

    d_low  <- max(abs(core_if$low[m]  - core_ref$low[m]))
    d_high <- max(abs(core_if$high[m] - core_ref$high[m]))
    expect_true(
      d_low < ci_max_abs_diff,
      info = paste0(label, " max|low diff| = ", signif(d_low, 3),
                    " (threshold=", ci_max_abs_diff, ")")
    )
    expect_true(
      d_high < ci_max_abs_diff,
      info = paste0(label, " max|high diff| = ", signif(d_high, 3),
                    " (threshold=", ci_max_abs_diff, ")")
    )
  }

  expect_if_dims <- function(out, n) {
    expect_true(is.list(out$influence.function))
    expect_true(length(out$influence.function) >= 1)
    mat <- out$influence.function[[1]]
    expect_true(is.matrix(mat))
    expect_equal(nrow(mat), n)
    expect_equal(ncol(mat), length(out$time))
  }

  set.seed(1)
  n  <- 400
  te <- rexp(n, rate = 0.08)
  tc <- rexp(n, rate = 0.05)
  t_s   <- pmin(te, tc)
  eps_s <- as.integer(te <= tc)

  out_if_s <- calc(t_s, eps_s, error = "if",        conf_type = "plain", return_if = TRUE)
  out_gw_s <- calc(t_s, eps_s, error = "greenwood", conf_type = "plain", return_if = TRUE)

  core_if_s <- extract_core(out_if_s)
  core_gw_s <- extract_core(out_gw_s)

  expect_equal(out_if_s$error, "if")
  expect_equal(out_gw_s$error, "greenwood")

  expect_if_dims(out_if_s, n)
  expect_if_dims(out_gw_s, n)

  expect_methods_close(
    core_if_s, core_gw_s,
    se_median_ratio_bounds = c(0.85, 1.15),
    se_q05_q95_bounds      = c(0.60, 1.60),
    ci_max_abs_diff        = 0.05,
    label = "unweighted outcome='s' (if vs greenwood):"
  )

  set.seed(2)
  w_s <- runif(n, min = 0.5, max = 2.0)

  out_if_sw <- calc(t_s, eps_s, w = w_s, error = "if",        conf_type = "plain", return_if = TRUE)
  out_gw_sw <- calc(t_s, eps_s, w = w_s, error = "greenwood", conf_type = "plain", return_if = TRUE)

  core_if_sw <- extract_core(out_if_sw)
  core_gw_sw <- extract_core(out_gw_sw)

  expect_equal(out_if_sw$error, "if")
  expect_equal(out_gw_sw$error, "greenwood")

  expect_if_dims(out_if_sw, n)
  expect_if_dims(out_gw_sw, n)

  expect_methods_close(
    core_if_sw, core_gw_sw,
    se_median_ratio_bounds = c(0.80, 1.20),
    se_q05_q95_bounds      = c(0.50, 1.80),
    ci_max_abs_diff        = 0.08,
    label = "weighted outcome='s' (if vs greenwood):"
  )

  set.seed(3)
  n2 <- 350
  t1 <- rexp(n2, rate = 0.06)
  t2 <- rexp(n2, rate = 0.04)
  tc2 <- rexp(n2, rate = 0.03)
  mat <- cbind(t1, t2, tc2)
  ix  <- max.col(-mat, ties.method = "first")
  t_c   <- mat[cbind(seq_len(n2), ix)]
  eps_c <- c(1L, 2L, 0L)[ix]

  out_if_c    <- calc(t_c, eps_c, error = "if",    conf_type = "plain", return_if = TRUE)
  out_delta_c <- calc(t_c, eps_c, error = "delta", conf_type = "plain", return_if = TRUE)

  core_if_c    <- extract_core(out_if_c)
  core_delta_c <- extract_core(out_delta_c)

  expect_equal(out_if_c$error, "if")
  expect_equal(out_delta_c$error, "delta")

  expect_if_dims(out_if_c, n2)
  expect_if_dims(out_delta_c, n2)

  expect_true(length(out_if_c$aj) > 0)
  expect_equal(out_if_c$aj + out_if_c$surv, rep(1, length(out_if_c$time)), tolerance = 1e-10)

  expect_methods_close(
    core_if_c, core_delta_c,
    se_median_ratio_bounds = c(0.80, 1.20),
    se_q05_q95_bounds      = c(0.50, 1.80),
    ci_max_abs_diff        = 0.08,
    label = "unweighted outcome='c' (if vs delta):"
  )
})
