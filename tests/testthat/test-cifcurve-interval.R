test_that("cifcurve interpolates within finite support and keeps exact jumps", {
  d <- data.frame(
    start = c(0, 4, 0), stop = c(2, 4, 6), code = c(5, 5, 7),
    w = c(1, 2, 1), include = c(TRUE, TRUE, TRUE)
  )
  fit <- cifcurve(Interval(l = start, r = stop, epsilon = code) ~ 1, d, weights = w,
                  subset.condition = "include", code.event1 = 5, code.event2 = 7)
  expect_s3_class(fit, "cifcurve")
  expect_s3_class(fit, "survfit")
  expect_identical(fit$reporting$outcome.type, "competing-risk")
  expected <- c(0, 1 / 16, 3 / 16, 1 / 4, 1 / 4, 3 / 4)
  times <- c(0, .5, 1.5, 2, 3.99, 4)
  s <- summary(fit, times = times, scale = 2, data.frame = TRUE)
  expect_equal(s$time, times / 2)
  expect_equal(s$cif, expected)
  expect_equal(s$surv, 1 - expected)
  expect_equal(s$estimate.status,
               c("estimated", "interpolated", "interpolated", rep("estimated", 3)))
  expect_equal(s$n.risk.confirmed, c(3L, rep(1L, 5)))
  expect_true(all(is.na(s[c("n.risk", "n.event", "std.err", "lower", "upper")])))
  expect_equal(summary(fit)$surv, fit$surv)
  expect_equal(summary(fit)$n.risk.confirmed, fit$n.risk.confirmed)
  expect_equal(summary(fit)$estimate.status, fit$estimate.status)
  expect_s3_class(summary(fit, times), "summary.survfit")
  expect_error(cifcurve(Interval(start, stop, code) ~ 1, d,
                       code.event1 = 5, code.event2 = 7, null.hypothesis = .1),
               "point estimation only")
})

test_that("unresolved tails are tracked separately in each stratum", {
  d <- data.frame(
    left = c(0, 3, 8, 0, 3), right = c(2, 3, Inf, 2, 3),
    event = c(1, 1, 0, 1, 1), group = c("A", "A", "A", "B", "B")
  )
  f <- Interval(left, right, event) ~ group
  fit <- cifcurve(f, d)
  expect_named(fit$strata, c("group=A", "group=B"))
  s <- summary(unserialize(serialize(fit, NULL)), c(.5, 8, 9), extend = TRUE)
  expect_equal(s$cif, c(1 / 12, 2 / 3, NA, 1 / 8, 1, 1))
  expect_equal(s$estimate.status, c("interpolated", "estimated", "unresolved_tail",
                                     "interpolated", "estimated", "estimated"))
  expect_equal(s$n.risk.confirmed, c(2L, 1L, 0L, 1L, 0L, 0L))
  expect_equal(summary(fit["group=B"], 9, extend = TRUE)$cif, 1)
  expect_equal(summary(fit[c(2, 1)], .5)$cif, c(1 / 8, 1 / 12))
  reordered <- fit[c(2, 1)]
  expect_equal(reordered$n.risk.confirmed, summary(reordered)$n.risk.confirmed)
  expect_equal(reordered$estimate.status, summary(reordered)$estimate.status)
  expect_length(summary(fit, 9)$cif, 0)
  expect_error(plot(fit), "Plotting Interval", fixed = TRUE)
  expect_error(cifplot(f, d), "Plotting Interval", fixed = TRUE)
  expect_error(cifpanel(formula = f, data = d), "Plotting Interval", fixed = TRUE)
})

test_that("cifplot draws interpolated segments and exact jumps from a saved interval fit", {
  d <- data.frame(
    l = c(0, 4, 0, 0, 4, 0), r = c(2, 4, 6, 2, 4, 6),
    epsilon = rep(c(1, 1, 2), 2), group = rep(c("A", "B"), each = 3)
  )
  attr(d$group, "labels") <- c(Control = "A", Active = "B")
  fit <- cifcurve(Interval(l, r, epsilon) ~ group, d)
  fit <- unserialize(serialize(fit, NULL))
  result <- cifplot(fit)
  expect_s3_class(result, "cifplot")
  expect_identical(result$survfit.info$survfit, fit)
  expect_equal(unname(attr(result$plot, "cifmodeling_display_labels")[c("group=A", "group=B")]),
               c("Control", "Active"))
  layer <- ggplot2::ggplot_build(result$plot)$data[[1L]]
  first <- layer[layer$group == layer$group[1L], ]
  expect_equal(first$x, c(0, 2, 4, 4, 6))
  expect_equal(first$y, c(0, 1 / 3, 1 / 3, 2 / 3, 2 / 3))
  expect_equal(mean(first$y[1:2]), summary(fit[1], 1)$cif)
  surv <- cifplot(fit, type.y = "surv", limits.x = c(.5, 2.5))$plot
  expect_equal(surv$data$estimate, 1 - result$plot$data$estimate)
  expect_s3_class(surv$coordinates, "CoordCartesian")
  expect_warning(cifplot(fit, add.conf = TRUE), "Options unavailable")
  expect_error(cifplot(fit, type.y = "cloglog"), "support type.y")
})

test_that("interval plots align confirmed counts and point estimates to the axis breaks", {
  d <- data.frame(
    l = c(0, 4, 0, 8, 0, 2, 0, 7), r = c(2, 4, 6, Inf, 2, 2, 4, Inf),
    epsilon = rep(c(1, 1, 2, 0), 2), group = rep(c("A", "B"), each = 4)
  )
  fit <- cifcurve(Interval(l, r, epsilon) ~ group, d)
  out <- cifplot(fit, style = "ggsurvfit", add.estimate.table = TRUE,
                 breaks.x = c(0, 1, 4, 7, 8), label.strata = c(A = "Control", B = "Active"))
  marks <- out$plot$layers[[2L]]$data
  expect_equal(marks$time, c(8, 7))
  expect_equal(marks$estimate, c(summary(fit[1], 8)$cif, summary(fit[2], 7)$cif))
  expect_s3_class(out$patchwork, "patchwork")
  as_text <- function(g) {
    c(g$label, unlist(lapply(c(g$grobs, g$children), as_text), use.names = FALSE))
  }
  g <- patchwork::patchworkGrob(out$patchwork)
  panels <- g$layout[grepl("^panel_patch-[123]$", g$layout$name), ]
  expect_equal(length(unique(panels$l)), 1L)
  expect_true(all(diff(panels$t) > 0))
  text <- as_text(g)
  expect_true(all(c("Confirmed at risk", "Estimate", "\u25A0 Control", "\u25A0 Active",
                    "0.12*", "\u2014",
                    "* Interpolated assuming uniform allocation within the support interval.") %in% text))
  endpoints <- cifplot(fit, add.estimate.table = TRUE, breaks.x = c(0, 2, 4, 8))
  endpoint.text <- as_text(patchwork::patchworkGrob(endpoints$patchwork))
  expect_false(any(startsWith(endpoint.text, "* Interpolated")))
  expect_null(cifplot(fit, add.risktable = FALSE, add.censor.mark = FALSE)$patchwork)
  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  cifplot(fit, filename.ggsave = path)
  expect_gt(file.info(path)$size, 0)
})

test_that("exact Interval observations with right censoring agree with AJ", {
  d <- data.frame(
    time = c(0, 2, 2, 3, 1, 2, 3, 4), event = c(1, 0, 2, 0, 1, 2, 0, 0),
    group = rep(c("A", "B"), each = 4), w = c(2, 1, 3, 1, 1, 1, 2, 2)
  )
  d$right <- ifelse(d$event == 0, Inf, d$time)
  interval <- cifcurve(Interval(time, right, event) ~ group, d, weights = "w")
  aj <- cifcurve(Event(time, event) ~ group, d, weights = "w", n.risk.type = "unweighted",
                 outcome.type = "competing-risk")
  times <- c(0, 1.5, 2, 3)
  expect_equal(summary(interval, times)$surv, summary(aj, times)$surv)
  expect_equal(interval$surv, aj$surv)
  expect_equal(interval$n.risk.confirmed, aj$n.risk.unweighted)
  expect_equal(summary(interval, times)$n.risk.confirmed,
               summary(aj, times)$n.risk)
})

test_that("interval responses retain aligned weights and strata after row selection", {
  d <- data.frame(
    left = c(0, 0, 3, 3, 4), right = c(2, 2, 3, 3, Inf),
    event = c(5, 5, 5, 5, 9), w = c(1, 100, 3, 0, 50),
    group = c("A", "unused", "A", "zero", "A"),
    include = c(TRUE, TRUE, TRUE, TRUE, FALSE)
  )
  keep_group <- function(x) replace(x, x == "unused", NA_character_)
  fit <- cifcurve(Interval(left, right, event) ~ keep_group(group), d,
                  weights = w, subset.condition = "include",
                  code.event1 = 5, code.censoring = 9)
  expect_equal(fit$n, 2L)
  expect_equal(summary(fit, .5)$cif, 1 / 16)
  expect_equal(summary(fit, c(0, 2, 3, 4), extend = TRUE)$n.risk.confirmed,
               c(2L, 1L, 1L, 0L))
  expect_identical(names(fit$strata), "keep_group(group)=A")
})
