mkdf <- function(n = 20) {
  set.seed(1)
  df <- data.frame(
    t = rexp(n, rate = 0.2),
    d = sample(c(0L, 1L), n, replace = TRUE),
    sex = factor(sample(c("F","M"), n, replace = TRUE)),
    fruitq1 = factor(sample(c(0,1,2), n, replace = TRUE)),
    strata = factor(sample(letters[1:2], n, replace = TRUE)),
    w = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )
  df
}

test_that("Event() allowed NA and retained length; use na.omit as necessary", {
  df <- mkdf()
  df$t[1] <- NA
  df$d[2] <- NA

  ev <- with(df, Event(t, d))
  expect_s3_class(ev, "Event")
  expect_equal(nrow(ev), nrow(df))
  expect_true(any(is.na(ev[,1])))
  expect_true(any(is.na(ev[,2])))

  mf <- model.frame(Event(t, d) ~ sex, data = df, na.action = na.omit)
  y <- model.extract(mf, "response")
  expect_true(nrow(y) < nrow(df))
  expect_false(any(is.na(y)))
})

test_that("normalize_time_event() allowed only expected codes for events", {
  df <- mkdf(6)
  res <- normalize_time_event(df$t, c(0L,1L,NA,0L,1L,2L), allowed = c(0,1,2))
  expect_equal(length(res$event), 6L)

  expect_error(
    normalize_time_event(df$t, c(0,-1,1,0,1,0)),
    class = "cifmodeling_ev_codes"
  )

  expect_error(
    normalize_time_event(df$t, c(0, 1.2, 1, 0, 1, 0)),
    class = "cifmodeling_ev_codes"
  )
})

test_that("Surv() and Event() in readSurv() yield the same event status", {
  df <- mkdf()
  df$t[1:2] <- NA; df$d[3] <- NA; df$sex[4] <- NA

  a <- readSurv(Surv(t, d) ~ sex, data = df, na.action = na.omit)
  b <- readSurv(Event(t, d) ~ sex, data = df, na.action = na.omit)
  expect_equal(a$d, b$d)
  expect_equal(a$t, b$t)
})

test_that("readSurv() handled strata/weights as expected", {
  df <- mkdf()
  df$t[1] <- NA
  df$d[2] <- NA
  df$sex[3] <- NA
  df$w[4] <- NA

  out <- readSurv(Event(t, d) ~ sex, data = df, weights = "w", na.action = na.omit)
  k <- length(out$t)
  expect_equal(unname(lengths(out[c("epsilon","d","d0","d1","d2","w","strata")])), rep(k, 7))
  expect_false(any(is.na(out$t)))
  expect_false(any(is.na(out$epsilon)))
})

test_that("readExposureDesign() allowed only expected codes for exposure", {
  df <- mkdf()
  de1 <- readExposureDesign(df, exposure = "fruitq1", code.exposure.ref = 1)
  expect_true(is.matrix(de1$x_a))
  expect_equal(nrow(de1$x_a), nrow(df))
  expect_true(all(colnames(de1$x_a) != "a__2"))

  df$fruitq1 <- factor(0)
  expect_error(readExposureDesign(df, "fruitq1"), "only one level")
})

test_that("check_outcome.type() / check_effect.measure() / check_error()", {
  expect_equal(check_outcome.type("s"), "SURVIVAL")
  expect_equal(check_outcome.type("competing risk"), "COMPETING-RISK")

  em <- check_effect.measure("rr", "Or")
  expect_equal(em$effect.measure1, "RR")
  expect_equal(em$effect.measure2, "OR")

  expect_equal(check_error(NULL, "SURVIVAL"), "greenwood")
  expect_equal(check_error(NULL, "COMPETING-RISK"), "delta")

  expect_warning(
    expect_equal(check_error("bad", "SURVIVAL"), "greenwood"),
    "SURVIVAL"
  )
  expect_warning(
    expect_equal(check_error("bad", "COMPETING-RISK"), "delta"),
    "COMPETING-RISK"
  )
})

test_that("read_time.point() yields expected outputs according to outcome.type", {
  df <- mkdf()
  expect_error(read_time.point(Event(t,d)~1, df, matrix(1, nrow(df), 1),
                               "SURVIVAL", code.censoring = 0,
                               should.terminate.time.point = TRUE,
                               time.point = NULL),
               "time.point is required")

  tp <- read_time.point(Event(t,d)~1, df, matrix(1, nrow(df), 1),
                        "BINOMIAL", code.censoring = 0,
                        should.terminate.time.point = TRUE,
                        time.point = NULL)
  expect_true(is.infinite(tp))

  x_a <- model.matrix(~ fruitq1, data = df)[, -1, drop = FALSE]
  tp2 <- read_time.point(Event(t,d)~sex, df, x_a,
                         "PROPORTIONAL", code.censoring = 0,
                         should.terminate.time.point = TRUE,
                         time.point = NULL)
  expect_true(all(tp2 >= 0))
  expect_true(is.unsorted(tp2, strictly = FALSE) == FALSE)
})

test_that("checkInput() handles NA as expected", {
  df <- mkdf()
  df$t[1] <- NA; df$d[2] <- NA; df$sex[3] <- NA; df$fruitq1[4] <- NA
  out <- checkInput(
    data = df,
    formula = Event(t,d) ~ sex,
    exposure = "fruitq1",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    code.exposure.ref = 0,
    outcome.type = "SURVIVAL",
    conf.level = 0.95,
    report.sandwich.conf = TRUE,
    report.boot.conf = NULL,
    nleqslv.method = "nleqslv",
    should.normalize.covariate = TRUE,
    strata = "strata",
    subset.condition = NULL,
    na.action = na.omit
  )
  expect_true(is.matrix(out$x_a))
  expect_true(ncol(out$x_l) >= 1)
})
