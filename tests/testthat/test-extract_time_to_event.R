# helper-data.R
# Small deterministic toy dataset for tests

make_toy_df <- function() {
  # status codes: 0=censor, 1=event1, 2=event2
  data.frame(
    t       = c(5,  7,  7,  9, 10, 12, NA,  4),
    status  = c(1,  0,  2,  0,  1,  2,  1,  0),
    x       = c(0,  1,  1,  0,  0,  1,  1,  0),
    strata  = factor(c("A","A","B","B","A","B","A","B")),
    w       = c(1,  2,  1,  1,  0.5, 1,  1,  1),
    stringsAsFactors = FALSE
  )
}

# test-extract_time_to_event.R
testthat::skip_on_cran()

.local_has_dataset <- function() {
  # If package ships 'diabetes.complications', use it; otherwise, use toy data.
  ok <- FALSE
  try({
    data("diabetes.complications", package = utils::packageName(), envir = environment())
    ok <- exists("diabetes.complications", inherits = FALSE)
  }, silent = TRUE)
  ok
}

get_data_for_tests <- function() {
  if (.local_has_dataset()) {
    get("diabetes.complications", envir = environment())
  } else {
    make_toy_df()
  }
}

testthat::test_that("extract_time_to_event(): basic flows for event1/event2/censor", {
  dat <- get_data_for_tests()

  # fall back to toy data if named columns differ
  has_expected <- all(c("t","status") %in% names(dat))
  if (!has_expected) {
    dat <- make_toy_df()
  }

  f <- Event(t, status) ~ 1

  # event1
  te1 <- extract_time_to_event(
    formula = f, data = dat, which_event = "event1",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    readUniqueTime = TRUE, dropEmpty = TRUE
  )
  testthat::expect_true(is.numeric(te1))
  testthat::expect_true(length(te1) >= 0)
  testthat::expect_true(length(unique(te1)) == length(te1)) # readUniqueTime

  # event2
  te2 <- extract_time_to_event(
    formula = f, data = dat, which_event = "event2",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    readUniqueTime = FALSE, dropEmpty = FALSE
  )
  testthat::expect_true(is.numeric(te2))
  # Not enforcing uniqueness here

  # censor
  tc <- extract_time_to_event(
    formula = f, data = dat, which_event = "censor",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    readUniqueTime = TRUE, dropEmpty = TRUE
  )
  testthat::expect_true(is.numeric(tc))
})

testthat::test_that("extract_time_to_event(): subset.condition and na.action work", {
  dat <- make_toy_df()
  f   <- Event(t, status) ~ x

  # keep only x==1
  te <- extract_time_to_event(
    formula = f, data = dat, which_event = "event1",
    subset.condition = "x == 1",
    na.action = stats::na.omit,
    readUniqueTime = TRUE, dropEmpty = TRUE
  )
  # Compare to manual filter
  manual <- unique(na.omit(dat$t[dat$status == 1 & dat$x == 1]))
  testthat::expect_equal(sort(te), sort(manual))
})

testthat::test_that("extract_time_to_event(): remapped codes are respected", {
  dat <- make_toy_df()
  # Remap: event1=10, event2=20, censor=0
  dat2 <- within(dat, { status <- ifelse(status == 1, 10, ifelse(status == 2, 20, 0)) })

  f <- Event(t, status) ~ 1
  te <- extract_time_to_event(
    formula = f, data = dat2,
    which_event = "event2",
    code.event1 = 10, code.event2 = 20, code.censoring = 0,
    readUniqueTime = TRUE, dropEmpty = TRUE
  )
  manual <- unique(na.omit(dat2$t[dat2$status == 20]))
  testthat::expect_equal(sort(te), sort(manual))
})

testthat::test_that("extract_time_to_event(): user_specified requires code", {
  dat <- make_toy_df()
  f   <- Event(t, status) ~ 1

  testthat::expect_error(
    extract_time_to_event(formula = f, data = dat, which_event = "user_specified"),
    regexp = "user_specified", ignore.case = TRUE
  )

  te <- extract_time_to_event(
    formula = f, data = dat, which_event = "user_specified",
    user_specified_code = 2,
    readUniqueTime = TRUE, dropEmpty = TRUE
  )
  manual <- unique(na.omit(dat$t[dat$status == 2]))
  testthat::expect_equal(sort(te), sort(manual))
})

testthat::test_that("extract_time_to_event(): invalid which_event fails clearly", {
  dat <- make_toy_df()
  f   <- Event(t, status) ~ 1
  testthat::expect_error(
    extract_time_to_event(formula = f, data = dat, which_event = "ev3"),
    regexp = "should be one of", ignore.case = TRUE
  )
})


