make_toy_df <- function() {
  data.frame(
    t       = c(5,  7,  7,  9, 10, 12, NA,  4),
    status  = c(1,  0,  2,  0,  1,  2,  1,  0),
    x       = c(0,  1,  1,  0,  0,  1,  1,  0),
    strata  = factor(c("A","A","B","B","A","B","A","B")),
    w       = c(1,  2,  1,  1,  0.5, 1,  1,  1),
    stringsAsFactors = FALSE
  )
}

testthat::skip_on_cran()

test_that("extract_time_to_event() with basic flows for event1/event2/censor works", {
  dat <- createTestData2()
  f <- Event(t, status) ~ 1

  te1 <- extract_time_to_event(
    formula = f, data = dat, which.event = "event1",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    read.unique.time = TRUE, drop.empty = TRUE
  )
  expect_true(is.numeric(te1))
  expect_true(length(te1) >= 0)
  expect_true(length(unique(te1)) == length(te1)) # readUniqueTime

  te2 <- extract_time_to_event(
    formula = f, data = dat, which.event = "event2",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    read.unique.time = FALSE, drop.empty = FALSE
  )
  testthat::expect_true(is.numeric(te2))

  tc <- extract_time_to_event(
    formula = f, data = dat, which.event = "censor",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    read.unique.time = TRUE, drop.empty = TRUE
  )
  expect_true(is.numeric(tc))
})

test_that("In extract_time_to_event() subset.condition and na.action work", {
  dat <- createTestData2()
  f   <- Event(t, status) ~ x

  te <- extract_time_to_event(
    formula = f, data = dat, which.event = "event1",
    subset.condition = "x == 1",
    na.action = stats::na.omit,
    read.unique.time = TRUE, drop.empty = TRUE
  )
  manual <- unique(na.omit(dat$t[dat$status == 1 & dat$x == 1]))
  testthat::expect_equal(sort(te), sort(manual))
})

testthat::test_that("In extract_time_to_event() remapped codes are respected", {
  dat <- createTestData2()
  dat2 <- within(dat, { status <- ifelse(status == 1, 10, ifelse(status == 2, 20, 0)) })

  f <- Event(t, status) ~ 1
  te <- extract_time_to_event(
    formula = f, data = dat2,
    which.event = "event2",
    code.event1 = 10, code.event2 = 20, code.censoring = 0,
    read.unique.time = TRUE, drop.empty = TRUE
  )
  manual <- unique(na.omit(dat2$t[dat2$status == 20]))
  testthat::expect_equal(sort(te), sort(manual))
})

testthat::test_that("extract_time_to_event() with user_specified requires code", {
  dat <- createTestData2()
  f   <- Event(t, status) ~ 1

  testthat::expect_error(
    extract_time_to_event(formula = f, data = dat, which.event = "user_specified"),
    regexp = "user_specified", ignore.case = TRUE
  )

  te <- extract_time_to_event(
    formula = f, data = dat, which.event = "user_specified",
    code.user.specified = 2,
    read.unique.time = TRUE, drop.empty = TRUE
  )
  manual <- unique(na.omit(dat$t[dat$status == 2]))
  testthat::expect_equal(sort(te), sort(manual))
})

testthat::test_that("In extract_time_to_event() invalid which_event fails clearly", {
  dat <- createTestData2()
  f   <- Event(t, status) ~ 1
  testthat::expect_error(
    extract_time_to_event(formula = f, data = dat, which.event = "ev3"),
    regexp = "should be one of", ignore.case = TRUE
  )
})
