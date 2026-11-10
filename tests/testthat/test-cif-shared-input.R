test_that("shared event preparation retains transformed-missing row provenance", {
  d <- data.frame(time=1:6,status=c(0,1,2,1,2,0),x=c(1,-1,2,3,4,5),
                  group=c("a","a",NA,"b","b","b"))
  prep <- suppressWarnings(cifmodeling:::cif_prepare_event(
    Event(time,status)~log(x), d, other.variables.analyzed="group",na.action=na.exclude))
  expect_identical(prep$row.index,c(1L,4L,5L,6L))
  expect_identical(which(!prep$included),c(2L,3L))
  expect_equal(prep$exclusion.reason[c(2,3)],c("missing","missing"))
  expect_s3_class(stats::model.response(prep$frame),"Event")
  expect_s3_class(prep$na.action,"exclude")
  expect_equal(prep$t,c(1,4,5,6))
})

test_that("shared response decoding supports arbitrary causes without global option changes", {
  d <- data.frame(time=1:4,status=c(9,0,7,8))
  prior <- getOption("cifmodeling.allowed")
  prep <- cifmodeling:::cif_prepare_event(Event(time,status)~1,d)
  expect_equal(prep$epsilon,d$status)
  expect_identical(getOption("cifmodeling.allowed"),prior)
  expect_error(cifmodeling:::cif_prepare_event(
    Event(time,status,allowed=c(0,1,2))~1,d),"codes")
})

test_that("legacy CIF preparation retains its RHS grouping and row meanings", {
  d <- data.frame(time=1:6,status=c(0,1,2,1,2,0),group=c("a","a",NA,"b","b","b"))
  prep <- cifmodeling:::cif_prepare_input(Event(time,status)~group,d,
                                        outcome.type="competing-risk")
  expect_equal(prep$t,c(1,2,4,5,6))
  expect_identical(prep$row.index,c(1L,2L,4L,5L,6L))
  expect_equal(as.character(prep$strata),c("a","a","b","b","b"))
  expect_equal(prep$strata_name,"group")
  expect_equal(prep$d0,as.integer(prep$epsilon==0))
})

test_that("explicit competing-risk curves retain AJ with absent status codes", {
  d <- data.frame(time = 1:4, status = c(0, 2, 0, 2))
  fit <- cifcurve(Event(time, status) ~ 1, d, outcome.type = "competing-risk")
  expect_identical(fit$type, "aalen-johansen")
  expect_equal(fit$aj, rep(0, 4))
  expect_error(cifcurve(Event(time, status) ~ 1, d, outcome.type = "survival"),
    "Observed status codes are incompatible")

  d$status <- c(0, 1, 0, 1)
  fit <- cifcurve(Event(time, status) ~ 1, d, outcome.type = "competing-risk",
    error = "aalen")
  expect_identical(fit$type, "aalen-johansen")
  expect_equal(fit$aj, c(0, 1 / 3, 1 / 3, 1))

  d$status <- c(2, 1, 2, 1)
  fit <- cifcurve(Event(time, status) ~ 1, d, outcome.type = "competing-risk")
  expect_equal(fit$aj, c(0, .25, .25, .5))
})

test_that("explicit Gray tests accept data without observed competing events", {
  d <- data.frame(time = rep(1:6, 2), group = rep(c("a", "b"), each = 6),
    status = c(1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1))
  reference <- gray_reference_slow(d$time, d$status, d$group)
  scalar <- ciftest(Event(time, status) ~ group, d,
    outcome.type = "competing-risk", test = "gray")
  joint <- ciftest_mdir(Event(time, status) ~ group, d,
    outcome.type = "competing-risk", test = "gray", directions = "unweighted")
  expect_equal(scalar$p.value, reference$p.value)
  expect_equal(joint$p.value, reference$p.value)
})

test_that("time-to-event models validate Surv response types before decoding", {
  skip_if_not_installed("survival")
  d <- data.frame(start = 0, stop = 1:4, status = c(0, 2, 0, 2),
    group = c(0, 0, 1, 1))
  counting <- survival::Surv(start, stop, status != 0) ~ group
  expect_error(ciftest(counting, d), "right.censored")
  expect_error(polyreg(survival::Surv(start, stop, status != 0) ~ 1,
    exposure = "group", data = d, outcome.type = "survival", time.point = 3),
    "right.censored")

  d$state <- factor(d$status, levels = 0:2)
  fit <- cifcurve(survival::Surv(stop, state) ~ 1, d, outcome.type = "competing-risk")
  expect_equal(fit$aj, rep(0, 4))
})

test_that("competing-risk polyreg allows two observed causes without censoring", {
  d <- data.frame(time = rep(rep(1:2, each = 4), 2),
    status = c(1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
    group = rep(0:1, each = 8))
  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    outcome.type = "competing-risk", time.point = 1)
  expect_identical(fit$outcome.type, "competing-risk")
  expect_equal(unname(coef(fit)[c(2, 4)]), log(c(2, 2 / 3)), tolerance = 1e-6)
  expect_error(polyreg(Event(time, status) ~ 1, exposure = "group",
    data = transform(d, status = ifelse(status == 2, 0, status)),
    outcome.type = "competing-risk", time.point = 1), "both event types")
})

test_that("fixed-time competing-risk polyreg requires both causes by the analysis time", {
  d <- data.frame(time = rep(1:4, 2), status = rep(c(0, 1, 2, 0), 2),
    group = rep(0:1, each = 4))
  expect_error(polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    time.point = 2), "time.point = 2.*code.event2 = 2")
  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d, time.point = 3)
  expect_true(all(is.finite(coef(fit))))
})

test_that("polyreg validates active event codes before evaluating the response", {
  d <- data.frame(time = rep(1:4, 2), status = rep(c(0, 1, 1, 0), 2),
    group = rep(0:1, each = 4))
  formula <- Event(stop("response evaluated"), status) ~ 1
  expect_error(polyreg(formula, exposure = "group", data = d,
    code.event1 = 0), "distinct non-negative integers")
  expect_error(polyreg(formula, exposure = "group", data = d,
    outcome.type = "proportional-competing-risk", code.event2 = 1),
    "distinct non-negative integers")
  expect_error(polyreg(formula, exposure = "group", data = d,
    outcome.type = "survival", code.censoring = 1), "distinct non-negative integers")
  expect_error(polyreg(formula, exposure = "group", data = d,
    outcome.type = NULL, code.event1 = c(0, 1)), "distinct non-negative integers")

  fit <- polyreg(Event(time, status) ~ 1, exposure = "group", data = d,
    outcome.type = "survival", code.event2 = NULL, time.point = 3)
  expect_true(all(is.finite(coef(fit))))

  recoded <- transform(d, status = ifelse(status == 1, 8L, 9L))
  custom <- polyreg(Event(time, status) ~ 1, exposure = "group", data = recoded,
    outcome.type = "survival", code.event1 = 8L, code.censoring = 9L, time.point = 3)
  expect_equal(coef(custom), coef(fit))
})
