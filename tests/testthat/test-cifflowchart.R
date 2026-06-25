test_that("categorical flowchart supports grouped and overall formulas", {
  dat <- data.frame(response = c("CR", "PD", NA, "CR"), arm = c("B", "A", "A", "B"))
  x <- cifmodeling:::.flowchart_prepare_data(response ~ arm, dat)
  expect_equal(x$total_n, 4)
  expect_true("Missing outcome" %in% names(x$outcome_counts[[1]]) || "Missing outcome" %in% names(x$outcome_counts[[2]]))
  y <- cifmodeling:::.flowchart_prepare_data(response ~ 1, dat)
  expect_false(y$has_group)
  expect_equal(y$analysis_n, 4)
})

test_that("Event(time, status) formulas classify final and tau status", {
  dat <- data.frame(time = c(10, 20, 30, 40), status = c(1, 2, 0, 0), arm = c("A", "A", "B", "B"))
  x <- cifmodeling:::.flowchart_prepare_data(Event(time, status) ~ arm, dat)
  expect_true("Event 1" %in% names(x$outcome_counts[[1]]))
  y <- cifmodeling:::.flowchart_prepare_data(Event(time, status) ~ arm, dat, time.point = 30)
  all_out <- unlist(lapply(y$outcome_counts, names), use.names = FALSE)
  expect_true("Event 1 by 30" %in% all_out)
  expect_true("Event 2 by 30" %in% all_out)
  expect_true("Event-free at 30" %in% all_out)
})

test_that("unsupported formula shapes error", {
  dat <- data.frame(response = 1:2, arm = 1:2, sex = 1:2)
  expect_error(cifmodeling:::.flowchart_prepare_data(response ~ arm, dat, time.point = 1), "time.point")
  expect_error(cifmodeling:::.flowchart_prepare_data(response ~ arm + sex, dat), "right-hand side")
  expect_error(cifmodeling:::.flowchart_prepare_data(Event(status) ~ arm, data.frame(status = 1:2, arm = 1:2)), "Event")
})

test_that("logical exclusions are counted and pre exclusion takes precedence", {
  dat <- data.frame(response = c("A", "B", "C", "D"), arm = c("x", "x", "y", "y"), pre = c(TRUE, FALSE, TRUE, FALSE), post = c(TRUE, TRUE, FALSE, TRUE))
  expect_silent(
    x <- cifmodeling:::.flowchart_prepare_data(
      response ~ arm,
      dat,
      pre.exclude = quote(pre),
      post.exclude = quote(post)
    )
  )
  expect_equal(sum(x$pre_counts), 2)
  expect_equal(sum(vapply(x$post_counts, sum, integer(1))), 2)
  expect_error(cifmodeling:::.flowchart_prepare_data(response ~ arm, transform(dat, pre = c(TRUE, NA, FALSE, FALSE)), pre.exclude = quote(pre)), "must not contain NA")
})

test_that("character and factor exclusions are counted by reason", {
  dat <- data.frame(response = c("A", "B", "C"), arm = c("x", "x", "y"), pre = factor(c("reason 1", "", NA)), post = c("", "reason 2", "reason 2"))
  x <- cifmodeling:::.flowchart_prepare_data(response ~ arm, dat, pre.exclude = quote(pre), post.exclude = quote(post))
  expect_equal(names(x$pre_counts), "reason 1")
  expect_true(all(vapply(x$post_counts, function(z) "reason 2" %in% names(z) || length(z) == 0, logical(1))))
  expect_error(cifmodeling:::.flowchart_prepare_data(response ~ arm, transform(dat, pre = 1:3), pre.exclude = quote(pre)), "Numeric")
})

test_that("missing group warns and labels/orders are reflected", {
  dat <- data.frame(response = factor(c("PD", "CR", "PD"), levels = c("CR", "PD")), arm = factor(c("B", "A", NA), levels = c("A", "B")))
  expect_warning(x <- cifmodeling:::.flowchart_prepare_data(response ~ arm, dat, order.strata = c("B", "A"), label.strata = c(B = "Arm B", A = "Arm A"), label.events = c(CR = "Complete", PD = "Progression")), "Missing group")
  expect_equal(x$groups[1:2], c("B", "A"))
  expect_equal(x$group_labels[1:2], c("Arm B", "Arm A"))
  expect_true(any(unlist(lapply(x$outcome_counts, names)) %in% c("Complete", "Progression")))
  expect_true("Missing group" %in% x$groups)
})

test_that("events exactly at time.point are not overwritten as event-free", {
  dat <- data.frame(
    time = c(30, 30, 30, 40),
    status = c(1, 2, 0, 0),
    arm = c("A", "A", "A", "A")
  )

  x <- cifmodeling:::.flowchart_prepare_data(
    Event(time, status) ~ arm,
    dat,
    time.point = 30
  )

  out <- names(x$outcome_counts[[1]])

  expect_true("Event 1 by 30" %in% out)
  expect_true("Event 2 by 30" %in% out)
  expect_true("Event-free at 30" %in% out)
})

test_that("factor outcome levels are preserved", {
  dat <- data.frame(
    response = factor(c("PD", "CR", "PD"), levels = c("CR", "PD")),
    arm = "A"
  )

  x <- cifmodeling:::.flowchart_prepare_data(response ~ arm, dat)

  expect_equal(names(x$outcome_counts[[1]]), c("CR", "PD"))
})

test_that("named label.events maps factor outcome levels", {
  dat <- data.frame(
    response = factor(c("PD", "CR", "PD"), levels = c("CR", "PD")),
    arm = "A"
  )

  x <- cifmodeling:::.flowchart_prepare_data(
    response ~ arm,
    dat,
    label.events = c(CR = "Complete", PD = "Progression")
  )

  expect_equal(names(x$outcome_counts[[1]]), c("Complete", "Progression"))
})

test_that("overall formula does not create a group node in dot output", {
  dat <- data.frame(response = c("CR", "PD", "CR"))

  x <- cifmodeling:::.flowchart_prepare_data(response ~ 1, dat)
  dot <- cifmodeling:::.flowchart_make_dot(x)

  expect_false(grepl("group1 \\[label", dot))
})

test_that("censoring before time.point is classified separately", {
  dat <- data.frame(
    time = c(10, 30, 40),
    status = c(0, 0, 0),
    arm = c("A", "A", "A")
  )

  x <- cifmodeling:::.flowchart_prepare_data(
    Event(time, status) ~ arm,
    dat,
    time.point = 30
  )

  expect_equal(unname(x$outcome_counts[[1]]["Censored before 30"]), 1)
  expect_equal(unname(x$outcome_counts[[1]]["Event-free at 30"]), 2)
})

test_that("withdraw.consent and ineligible are counted before pre.exclude", {
  dat <- data.frame(
    response = c("A", "B", "C", "D", "E"),
    arm = c("x", "x", "y", "y", "y"),
    wd = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    inelig = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    pre = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    post = c(TRUE, TRUE, TRUE, TRUE, FALSE)
  )

  expect_silent(
    x <- cifmodeling:::.flowchart_prepare_data(
      response ~ arm,
      dat,
      withdraw.consent = quote(wd),
      ineligible = quote(inelig),
      pre.exclude = quote(pre),
      post.exclude = quote(post)
    )
  )

  expect_equal(sum(x$withdraw_counts), 1)
  expect_equal(sum(x$ineligible_counts), 1)
  expect_equal(sum(x$pre_counts), 1)
  expect_equal(sum(vapply(x$post_counts, sum, integer(1))), 1)
  expect_equal(sum(x$analysis_n), 1)
})
