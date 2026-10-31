test_that("single-plot extraction preserves the ggplot without rebuilding it", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_line()
  expect_identical(as_ggplot(p), p)
  expect_identical(as_ggplot(structure(list(plot = p), class = "cifplot")), p)
  expect_error(as_ggplot(structure(list(plot = NULL), class = "cifplot")), "single plot")
  expect_error(as_ggplot(list(plot = p)), "Supply")
  expect_error(as_ggplot(p, typo = TRUE), "Unsupported")
  pred <- data.frame(.row = 1:3, estimate = c(.1, .2, .3))
  expect_identical(as_prediction_data(pred), pred)
  expect_error(as_prediction_data(pred, newdata = pred), "Unsupported")
})

test_that("polyreg summaries keep stored inference and explicitly recompute requested inference", {
  table <- data.frame(term = "a", estimate = .3, std.error = .2,
    p.value = .13, conf.low = -.09, conf.high = .69)
  attr(table, "coefficient.index") <- 2L
  part <- structure(list(tidy = table, glance = data.frame(effect.measure = "RR")),
                    class = "modelsummary_list")
  x <- structure(list(coef = c(0, .3), vcov = diag(c(.1, .04)),
    vcov_bootstrap = diag(c(.2, .09)), summary = list(event1 = part, event2 = part),
    conf.int = .95, outcome.type = "competing-risk",
    boot.method = list(report.sandwich.conf = TRUE, report.boot.conf = FALSE)), class = "polyreg")
  expect_equal(generics::tidy(x)$conf.low, table$conf.low)
  td <- generics::tidy(x, variance = "if", conf.int = .9)
  expect_equal(td$std.error, .2)
  expect_equal(td$conf.low, .3 - qnorm(.95) * .2)
  expect_equal(generics::tidy(x, variance = "bootstrap")$std.error, .3)
  expect_equal(vcov(x, variance = "if"), x$vcov)
  expect_equal(vcov(x, type = "sandwich"), x$vcov)
  expect_error(vcov(x, type = "bootstrap", variance = "if"), "different")
  expect_error(generics::tidy(x, varaince = "if"), "Unsupported")
  expect_error(generics::tidy(x, conf.int = NULL), "conf.int")
  expect_error(generics::glance(x, typo = TRUE), "Unsupported")
  expect_error(generics::augment(x, newdata = data.frame(a = 1)), "Unsupported")
  out <- as_modelsummary(x, variance = "if", conf.int = .9)
  expect_named(out, c("event1", "event2"))
  expect_s3_class(out$event1, "modelsummary_list")
  expect_equal(out$event1$glance$variance, "if")
  expect_equal(out$event1$glance$conf.level, .9)
  expect_equal(out$event1$glance$coefficient.scale, "log effect")
  expect_match(attr(out$event1, "notes"), "confidence")
  expect_s3_class(as_modelsummary(x, event = "event1"), "modelsummary_list")
  x$boot.method$boot.bca <- TRUE
  expect_equal(generics::tidy(x, variance = "bootstrap")$std.error, .3)
  x$outcome.type <- "proportional-competing-risk"
  expect_error(generics::tidy(x, variance = "bootstrap"), "BCa")
  x$outcome.type <- "competing-risk"
  x$vcov <- NULL
  expect_error(generics::tidy(x, variance = "if"), "unavailable")
})

test_that("the polyreg fitter stores coefficient indices for later extraction", {
  skip_on_cran()
  data(diabetes.complications)
  d <- diabetes.complications
  d$a <- as.integer(d$fruitq == "Q1")
  fit <- polyreg(Event(t, epsilon) ~ 1, exposure = "a", data = d,
    outcome.type = "competing-risk", time.point = 8, code.exposure.ref = 1,
    report.optim.convergence = FALSE)
  expect_identical(fit$estimand$exposure.labels, c("1", "0"))
  for (label in c("event1", "event2")) {
    td <- generics::tidy(fit, event = label, variance = "if", conf.int = .9)
    ix <- attr(fit$summary[[label]]$tidy, "coefficient.index")
    expect_equal(td$estimate, unname(coef(fit)[ix]))
    expect_equal(td$std.error, unname(sqrt(diag(vcov(fit, variance = "if"))[ix])))
  }
  predictions <- as_prediction_data(fit, exposure.level = "0")
  k <- fit$estimand$exposure.levels
  expect_equal(nrow(predictions), 2L * nobs(fit))
  expect_equal(predictions$estimate[predictions$category == "event1"],
    as.numeric(fit$diagnostics$potential.CIFs[, 2L]))
  expect_equal(predictions$estimate[predictions$category == "event2"],
    as.numeric(fit$diagnostics$potential.CIFs[, k + 2L]))
  expect_true(all(as.character(predictions$a) == "0"))
  expect_error(as_prediction_data(fit, newdata = d), "Unsupported")
  # Non-first references must map to the first potential-outcome column.
  expect_equal(as_prediction_data(fit, event = "event1")$estimate,
    as.numeric(fit$diagnostics$potential.CIFs[, 1L]))
  fit$estimand$exposure.labels <- NULL
  expect_error(as_prediction_data(fit), "legacy")
  skip_if_not_installed("modelsummary")
  tab <- modelsummary::modelsummary(as_modelsummary(fit), output = "data.frame")
  expect_true(is.data.frame(tab))
  expect_gt(nrow(tab), 0L)
})
