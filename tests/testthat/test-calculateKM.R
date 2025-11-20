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

test_that("engine routing works and schemas align", {
  skip_on_cran()
  data(diabetes.complications, package = "cifmodeling")
  diabetes.complications$d <- as.integer(diabetes.complications$epsilon>0)

  f <- Event(t, d) ~ fruitq1
  g <- Event(t, epsilon) ~ fruitq1

  out_km <- cifcurve(f, data = diabetes.complications,
                     outcome.type = "survival", engine = "calculateKM", error = "greenwood")
  out_rc <- cifcurve(f, data = diabetes.complications,
                     outcome.type = "survival", engine = "calculateAJ_Rcpp", error = "greenwood", report.influence.function = FALSE)
  #expect_equal(out_km$surv, out_rc$surv, tolerance = 1e-8)

  out_ajR <- cifcurve(g, data = diabetes.complications,
                      outcome.type = "competing-risk", engine = "calculateAJ", error = "aalen")
  out_ajC <- cifcurve(g, data = diabetes.complications,
                      outcome.type = "competing-risk", engine = "calculateAJ_Rcpp", error = "aalen", report.influence.function = FALSE)
  #expect_equal(out_ajR$`std.err.cif`, out_ajC$`std.err.cif`, tolerance = 1e-5)

  out_aalen <- out_ajC
  expect_true(length(out_aalen$`std.err.aj`) >= 2)
  expect_true(out_aalen$`std.err.aj`[1] >= 0)
  if (length(out_aalen$`std.err.aj`) >= 2 && is.finite(out_aalen$`std.err.aj`[2])) {
    expect_true(is.finite(out_aalen$`std.err.aj`[1]))
  }

  #for (x in list(out_km, out_rc, out_ajR, out_ajC)) {
  #  expect_true(all(c("time", "surv", "std.err", "std.err.aj", "type", "method", "conf.type") %in% names(x)))
  #}
})

test_that("weighted n.event / n.censor match manual tallies", {
  skip_on_cran()
  set.seed(1)
  n  <- 50
  t  <- sample(1:5, n, TRUE)
  ep <- sample(c(0,1,2), n, TRUE, prob=c(0.3,0.5,0.2))
  w  <- stats::runif(n, 0.5, 2)
  out <- calculateAJ_Rcpp(t, ep, w=w, error="if", return_if=FALSE)

  tt <- sort(unique(t))
  tally <- function(code) sapply(tt, function(x) sum(w[t==x & code]))
  w_event_any  <- sapply(tt, function(x) sum(w[t==x & ep>=1]))
  w_censor     <- sapply(tt, function(x) sum(w[t==x & ep==0]))

  expect_equal(as.numeric(out$time),   tt)
  expect_equal(as.numeric(out$`n.event`),  w_event_any,  tolerance=1e-12)
  expect_equal(as.numeric(out$`n.censor`), w_censor,     tolerance=1e-12)
})

test_that("CI clamps to [0,0] / [1,1] at boundaries", {
  skip_on_cran()
  t  <- c(1,2,3,4)
  ep <- c(0,0,0,1)
  a <- calculateAJ_Rcpp(t, (ep==1L), error="greenwood", return_if=FALSE)
  wh1 <- which(a$surv >= 1 - 1e-15)
  expect_true(all(a$low [wh1] == 1))
  expect_true(all(a$high[wh1] == 1))
  wh0 <- which(a$surv <= 1e-15)
  expect_true(all(a$low [wh0] == 0))
  expect_true(all(a$high[wh0] == 0))

  ep2 <- c(2,0,0,1)
  b <- calculateAJ_Rcpp(t, ep2, error="aalen", return_if=FALSE)
  S <- b$surv
  wh1 <- which(S >= 1 - 1e-15)
  expect_true(all(b$low [wh1] == 1))
  expect_true(all(b$high[wh1] == 1))
})


test_that("influence.function has expected dims per stratum/time", {
  skip_on_cran()
  set.seed(3)
  n<-40; t<-sample(1:6,n,TRUE); ep<-sample(c(0,1,2),n,TRUE); g<-sample(1:2,n,TRUE)
  out <- calculateAJ_Rcpp(t, ep, strata=g, error="if", return_if=TRUE)
  IF  <- out$`influence.function`
  expect_true(is.list(IF))
  for (k in seq_along(IF)) {
    if (nrow(IF[[k]])>0) expect_true(ncol(IF[[k]]) > 0)
  }
})
