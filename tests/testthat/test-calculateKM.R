test_that("calculateKM() and util_get_surv() produced expected KM for survival data", {
  skip_on_cran()
  testthat::skip_if_not_installed("mets")
  testthat::skip_if_not_installed("Rcpp")
  data("diabetes.complications")
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon==0)
  resC <- suppressWarnings(mets::phreg(Surv(t, d) ~ 1, data=diabetes.complications))
  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = diabetes.complications, type = "survival", times = diabetes.complications$t, individual.time = TRUE, se = FALSE, km = TRUE, tminus = TRUE))
  expected <- out_predict1$surv
  expected <- round(expected[1:10,],digit=3)
  out_km <- calculateKM(diabetes.complications$t, diabetes.complications$d)
  tested <- as.matrix(util_get_surv(diabetes.complications$t, out_km$surv, out_km$time))
  tested <- round(tested[1:10,], digit=3)
  expect_equal(tested, expected)
})

test_that("calculateKM() and util_get_surv() produced expected stratified KM for survival data", {
  skip_on_cran()
  testthat::skip_if_not_installed("mets")
  testthat::skip_if_not_installed("Rcpp")
  library(mets)
  data("diabetes.complications")
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon==0)
  resC <- suppressWarnings(mets::phreg(Surv(t, d) ~ strata(strata), data=diabetes.complications))
  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = diabetes.complications, type = "survival", times = diabetes.complications$t, individual.time = TRUE, se = FALSE, km = TRUE, tminus = TRUE))
  expected <- out_predict1$surv
  expected <- round(expected[1:10,],digit=3)
  out_km <- calculateKM(t=diabetes.complications$t, d=diabetes.complications$d, strata=diabetes.complications$strata)
  tested <- as.matrix(util_get_surv(diabetes.complications$t, out_km$surv, out_km$time, diabetes.complications$strata, out_km$strata, out_km$strata.levels))
  tested <- round(tested[1:10,], digit=3)
  expect_equal(tested, expected)
})

test_that("calculateAJ_Rcpp() and util_get_surv() produced expected KM for survival data", {
  skip_on_cran()
  testthat::skip_if_not_installed("mets")
  testthat::skip_if_not_installed("Rcpp")
  data("diabetes.complications")
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon==0)
  resC <- suppressWarnings(mets::phreg(Surv(t, d) ~ 1, data=diabetes.complications))
  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = diabetes.complications, type = "survival", times = diabetes.complications$t, individual.time = TRUE, se = FALSE, km = TRUE, tminus = TRUE))
  expected <- out_predict1$surv
  expected <- round(expected[1:10,],digit=3)
  out_km <- calculateAJ_Rcpp(t=diabetes.complications$t, epsilon=diabetes.complications$d)
  tested <- as.matrix(util_get_surv(diabetes.complications$t, out_km$surv, out_km$time))
  tested <- round(tested[1:10,], digit=3)
  expect_equal(tested, expected)
})

test_that("calculateAJ_Rcpp() and util_get_surv() produced expected stratified KM for survival data", {
  skip_on_cran()
  testthat::skip_if_not_installed("mets")
  testthat::skip_if_not_installed("Rcpp")
  library(mets)
  data("diabetes.complications")
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon==0)
  resC <- suppressWarnings(mets::phreg(Surv(t, d) ~ strata(strata), data=diabetes.complications))
  out_predict1 <- suppressWarnings(stats::predict(resC, newdata = diabetes.complications, type = "survival", times = diabetes.complications$t, individual.time = TRUE, se = FALSE, km = TRUE, tminus = TRUE))
  expected <- out_predict1$surv
  expected <- round(expected[1:10,],digit=3)
  out_km <- calculateAJ_Rcpp(t=diabetes.complications$t, epsilon=diabetes.complications$d, strata=diabetes.complications$strata)
  tested <- as.matrix(util_get_surv(diabetes.complications$t, out_km$surv, out_km$time, diabetes.complications$strata, out_km$strata, out_km$strata.levels))
  tested <- round(tested[1:10,], digit=3)
  expect_equal(tested, expected)
})

test_that("drop-in replacement: no strata, no weights (greenwood / tsiatis)", {
  skip_on_cran()
  set.seed(123)

  n  <- 120
  t  <- rexp(n, rate = 0.2)
  ep <- rbinom(n, 1, 0.6)  # {0,1}
  stopifnot(any(ep == 1))  # at least one event

  call_calculateKM <- function(t, eps, w = NULL, strata = NULL, error = "greenwood") {
    fmls <- names(formals(calculateKM))
    args <- list()
    if ("t" %in% fmls) args$t <- t else if ("time" %in% fmls) args$time <- t else stop("calculateKM: time arg not found")
    if ("d" %in% fmls) args$d <- eps else if ("epsilon" %in% fmls) args$epsilon <- eps else stop("calculateKM: d/epsilon arg not found")
    if (!is.null(w)      && "w"      %in% fmls) args$w      <- w
    if (!is.null(strata) && "strata" %in% fmls) args$strata <- strata
    if ("error" %in% fmls) args$error <- error
    do.call(calculateKM, args)
  }

  for (err in c("greenwood", "tsiatis")) {
    old <- call_calculateKM(t, ep, error = err)
    new <- calculateAJ_Rcpp(t, ep, error = err, return_if = FALSE)
    old$std.err <- old$surv*old$std.err

    # Same time grid
    expect_equal(unname(old$time),   unname(new$time),   tolerance = 1e-12)

    # Core KM outputs (numeric)
    expect_equal(unname(old$surv),    unname(new$surv),    tolerance = 1e-10)
    expect_equal(unname(old$`n.risk`),unname(new$`n.risk`),tolerance = 1e-10)
    expect_equal(unname(old$`n.event`),unname(new$`n.event`),tolerance = 1e-10)
    expect_equal(unname(old$`n.censor`),unname(new$`n.censor`),tolerance = 1e-10)
    expect_equal(unname(old$`std.err`),unname(new$`std.err`),tolerance = 1e-8)

    # Meta fields (character / integers)
#   expect_identical(old$conf.type, new$conf.type)
    expect_identical(old$type,      new$type)
    expect_identical(old$method,    new$method)
    expect_identical(unname(old$n), unname(new$n))
    expect_identical(unname(old$strata),        unname(new$strata))
    expect_identical(unname(old$strata.levels), unname(new$strata.levels))

    # New fields exist (not compared to old)
    expect_true(!is.null(new$`std.err.aj`))
  }
})

test_that("drop-in replacement: with strata - weights", {
  skip_on_cran()
  set.seed(456)

  n  <- 150
  t  <- rexp(n, rate = 0.15)
  ep <- rbinom(n, 1, 0.5)        # {0,1}
  g  <- sample(1:3, n, replace = TRUE)

  for (err in c("greenwood", "tsiatis")) {
    old <- calculateKM(t, ep, strata = g, error = err)
    new <- calculateAJ_Rcpp(t, ep, strata = g, error = err, return_if = FALSE)
    old$std.err <- old$surv*old$std.err

    expect_equal(unname(old$time),    unname(new$time),    tolerance = 1e-12)
    expect_equal(unname(old$surv),    unname(new$surv),    tolerance = 1e-8)
    expect_equal(unname(old$`n.risk`),unname(new$`n.risk`),tolerance = 1e-8)
    expect_equal(unname(old$`n.event`),unname(new$`n.event`),tolerance = 1e-8)
    expect_equal(unname(old$`n.censor`),unname(new$`n.censor`),tolerance = 1e-8)
    expect_equal(unname(old$`std.err`),unname(new$`std.err`),tolerance = 1e-6)

#    expect_identical(old$conf.type, new$conf.type)
    expect_identical(old$type,      new$type)
    expect_identical(old$method,    new$method)
    expect_identical(unname(old$n), unname(new$n))
    expect_identical(unname(old$strata),        unname(new$strata))
    expect_identical(unname(old$strata.levels), unname(new$strata.levels))
  }
})
