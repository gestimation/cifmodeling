labelled_trial <- function() {
  d <- data.frame(time = rep(1:10, 2), status = rep(c(1, 2, 1, 2, rep(0, 6)), 2),
    arm = rep(0:1, each = 10))
  attr(d$time, "label") <- "Follow-up (months)"
  attr(d$arm, "label") <- "Assigned treatment"
  attr(d$arm, "labels") <- c(Placebo = 0, DES = 1)
  d
}

restore_labelled_fit <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(x, path)
  readRDS(path)
}

test_that("curve labels survive serialization and match formula and panel plots", {
  d <- labelled_trial()
  f <- Event(time, status) ~ arm
  fit <- cifcurve(f, d, outcome.type = "competing-risk")
  original <- fit$strata
  fit <- restore_labelled_fit(fit)
  expect_identical(fit$strata, original)
  expect_equal(fit[2]$reporting$display$levels, c("arm=1" = "DES"))
  expect_identical(as.character(summary(fit, times = 4)$strata), names(original))

  draw <- function(x, ...) cifplot(x, ..., outcome.type = "competing-risk", add.risktable = FALSE)$plot
  labels <- function(p) unname(p$scales$get_scales("colour")$labels)
  p <- cifplot(fit, add.risktable = FALSE)$plot
  expect_equal(labels(p), c("Placebo", "DES"))
  expect_equal(labels(draw(f, data = d)), labels(p))
  expect_equal(p$labels$x, "Follow-up (months)")
  expect_equal(p$labels$colour, "Assigned treatment")
  overrides <- c("1" = "Active")
  expect_equal(labels(draw(fit, label.strata = overrides)), c("Placebo", "Active"))
  expect_equal(draw(fit, label.x = "Months", label.legend = "Arm")$labels$x, "Months")
  expect_equal(draw(fit, label.legend = "Arm")$labels$colour, "Arm")

  panel <- cifpanel(plots = list(draw(fit, label.strata = overrides)), rows.columns.panel = c(1, 1))
  panel <- restore_labelled_fit(panel)
  expect_equal(labels(panel$list.plot[[1L]]), c("Placebo", "Active"))
  reordered <- cifpanel(plots = list(draw(fit)), rows.columns.panel = c(1, 1),
    level.strata = names(original), order.strata = "arm=1", label.strata = overrides)
  expect_equal(labels(reordered$list.plot[[1L]]), "Active")
  relabelled <- cifpanel(plots = list(draw(fit)), rows.columns.panel = c(1, 1),
    label.strata = c("1" = "Active"))
  relabelled <- restore_labelled_fit(relabelled)
  scale <- relabelled$list.plot[[1L]]$scales$get_scales("colour")
  expect_equal(scale$limits, names(original))
  expect_equal(scale$labels(names(original)), c("Placebo", "Active"))
  from_formula <- cifpanel(formula = f, data = d, outcome.type = "competing-risk",
    label.strata = overrides, add.risktable = FALSE, rows.columns.panel = c(1, 1))
  expect_equal(labels(from_formula$list.plot[[1L]]), labels(panel$list.plot[[1L]]))

  limited <- draw(fit, level.strata = 0:1, order.strata = 1, label.strata = overrides)
  expect_equal(limited$scales$get_scales("colour")$limits, "1")
  expect_equal(labels(limited), "Active")

  survival <- cifcurve(Event(time, status != 0) ~ arm, d, outcome.type = "survival")
  survival <- restore_labelled_fit(survival)
  expect_equal(cifplot(survival, add.risktable = FALSE)$plot$labels$y, "Survival")
})

test_that("event and flowchart labels preserve analytical codes and counts", {
  d <- labelled_trial()
  attr(d$status, "labels") <- c(Censored = 0, Failure = 1, Other = 2)
  ev <- Event(d$time, d$status)
  expect_equal(attr(ev, "inputAttributes")$event$labels, attr(d$status, "labels"))
  expect_equal(unname(ev[, "event"]), as.integer(d$status))
  x <- .flowchart_prepare_data(Event(time, status) ~ arm, d,
    label.strata = c("1" = "Active"), label.events = c("Event 2" = "Competing"))
  expect_identical(x$groups, c("0", "1"))
  expect_identical(x$group_labels, c("Placebo", "Active"))
  expect_identical(names(x$outcome_counts[[1L]]), c("Failure", "Competing", "Censored"))
  expect_equal(as.integer(x$outcome_counts[[1L]]), c(2, 2, 6))
})

test_that("joint strata are labelled without parsing or replacing their keys", {
  d <- labelled_trial()
  d$region <- rep(c("E", "W"), 10)
  attr(d$region, "labels") <- c(East = "E")
  fit <- cifcurve(Event(time, status) ~ arm + region, d, outcome.type = "competing-risk")
  expect_equal(names(fit$reporting$display$levels), names(fit$strata))
  expect_setequal(unname(fit$reporting$display$levels),
                  c("Placebo / East", "DES / East", "Placebo / W", "DES / W"))
})

test_that("modelsummary displays stored labels while fitted terms keep their keys", {
  skip_if_not_installed("modelsummary")
  d <- labelled_trial()
  raw <- polyreg(Event(time, status) ~ 1, "arm", data = d, time.point = 4)
  fit <- polyreg(Event(time, status) ~ 1, "arm", data = d, time.point = 4,
    label.exposure = "Treatment", label.strata = c("1" = "Active"))
  fit <- restore_labelled_fit(fit)
  expect_equal(coef(fit), coef(raw))
  expect_equal(generics::tidy(fit), generics::tidy(raw))
  expect_identical(fit$estimand$exposure.labels, c("0", "1"))
  expect_match(effect_label.polyreg(fit, add.exposure.levels = TRUE), "Active vs Placebo", fixed = TRUE)
  for (event in c("event1", "event2")) {
    shown <- as_modelsummary(fit, event = event)
    plain <- as_modelsummary(fit, event = event, display.labels = FALSE)
    expect_equal(shown$tidy$estimate, plain$tidy$estimate)
    expect_equal(shown$tidy$conf.low, plain$tidy$conf.low)
    expect_equal(shown$tidy$conf.high, plain$tidy$conf.high)
    tab <- modelsummary::msummary(shown, output = "data.frame", exponentiate = TRUE)
    expect_equal(unique(tab$term[tab$part == "estimates"]), "Treatment (Active vs Placebo)")
    expect_equal(plain$tidy$term, "arm, 1 vs 0")
  }
  both <- modelsummary::msummary(as_modelsummary(fit), output = "data.frame")
  expect_equal(unique(both$term[both$part == "estimates"]), "Treatment (Active vs Placebo)")
  if (requireNamespace("gtsummary", quietly = TRUE) &&
      requireNamespace("broom.helpers", quietly = TRUE)) {
    tab <- gtsummary::tbl_regression(fit, label = list(arm ~ "Explicit heading"))
    expect_equal(unname(tab$table_body$label), c("Explicit heading", "Placebo", "Active"))
  }
})
