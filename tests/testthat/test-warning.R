# tests/testthat/test-warning.R

get_warn_fns <- function() {
  fn   <- cifmodeling:::.warn
  msgs <- cifmodeling:::.msg
  list(fn = fn, msgs = msgs)
}

test_that(".warn: emits simple warning", {
  w <- get_warn_fns()
  res <- NULL
  expect_warning(
    { res <- w$fn("panel_disables_tables", .messages = w$msgs) },
    regexp = "^add.risktable/add.estimate.table are ignored in panel mode",
    fixed  = FALSE
  )
  expect_true(isTRUE(res) || is.null(res))
})

test_that(".warn: placeholder replacement (limits)", {
  w <- get_warn_fns()
  res <- NULL
  expect_warning(
    { res <- w$fn("est_outside_limits_y", arg="limits.y", a=0, b=1, .messages=w$msgs) },
    regexp = "fall outside `limits.y` = \\[0, 1\\]\\.",
    fixed  = FALSE
  )
  expect_true(isTRUE(res) || is.null(res))
})

test_that(".warn: placeholder replacement (counts, two-sentence message)", {
  w <- get_warn_fns()
  res <- NULL
  expect_warning(
    { res <- w$fn("plots_extra_dropped", n_plots=5, n_slots=3, .messages=w$msgs) },
    # 1文目だけにマッチさせ、末尾アンカーは使わない
    regexp = "^There are 5 plots but grid holds 3\\.",
    fixed  = FALSE
  )
  expect_true(isTRUE(res) || is.null(res))
})

test_that(".warn: unknown key -> fallback warning and FALSE", {
  w <- get_warn_fns()
  res <- NULL
  expect_warning(
    { res <- w$fn("this_key_does_not_exist", .messages = w$msgs) },
    regexp = "^Unknown warning key: this_key_does_not_exist$",
    fixed  = FALSE
  )
  expect_identical(isTRUE(res), FALSE)  # FALSE 期待
})

test_that(".warn: multiple replacements/types", {
  w <- get_warn_fns()
  res <- NULL
  expect_warning(
    { res <- w$fn("upper_outside_limits_y", arg="limits.y", a=0, b=100, .messages=w$msgs) },
    regexp = "upper CI values fall outside `limits.y` = \\[0, 100\\]\\.",
    fixed  = FALSE
  )
  expect_true(isTRUE(res) || is.null(res))
})
