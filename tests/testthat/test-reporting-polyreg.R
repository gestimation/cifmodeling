reporting_polyreg_fixture <- function(levels = c("B", "A")) {
  k <- length(levels) - 1L
  part <- function(event, offset) {
    beta <- if (event == "event1") seq(.2, by = .1, length.out = k) else seq(-.1, by = -.1, length.out = k)
    td <- data.frame(term = paste0("group, ", levels[-1], " vs ", levels[1]),
      estimate = beta, std.error = .2, p.value = .123,
      conf.low = beta - .37, conf.high = beta + .48)
    attr(td, "coefficient.index") <- offset + seq_len(k)
    list(tidy = td, glance = data.frame(effect.measure = event))
  }
  a <- part("event1", 1L); b <- part("event2", k + 2L)
  d <- data.frame(group = factor(rep(levels, each = 3), levels = sort(levels)), age = seq_len(3 * length(levels)))
  attr(d$group, "label") <- "Treatment"
  structure(list(coef = c(-1, a$tidy$estimate, -2, b$tidy$estimate),
    vcov = diag(.04, 2L * (k + 1L)), summary = list(event1 = a, event2 = b),
    diagnostics = d, exposure = "group", conf.int = .9,
    outcome.type = "competing-risk",
    estimand = list(exposure.labels = levels, effect.measure1 = "RR", effect.measure2 = "SHR", time.point = 8),
    boot.method = list(report.sandwich.conf = TRUE, report.boot.conf = FALSE),
    call = quote(polyreg(unavailable_formula, data = unavailable_data))), class = "polyreg")
}

test_that("polyreg exponentiation preserves stored inference and the original fit", {
  f <- reporting_polyreg_fixture(c("B", "A", "C")); original <- serialize(f, NULL)
  log <- generics::tidy(f, event = "both")
  ratio <- generics::tidy(f, event = "both", exponentiate = TRUE)
  for (nm in c("estimate", "conf.low", "conf.high")) expect_equal(ratio[[nm]], exp(log[[nm]]), tolerance = 0)
  for (nm in c("std.error", "statistic", "p.value", "event", "term")) expect_identical(ratio[[nm]], log[[nm]])
  expect_identical(ratio$event, rep(c("event1", "event2"), each = 2))
  expect_identical(serialize(f, NULL), original)
  expect_equal(as_modelsummary(f)$event1$tidy$estimate, log$estimate[log$event == "event1"])
})

test_that("polyreg accepts broom confidence arguments without replacing stored limits", {
  f <- reporting_polyreg_fixture()
  stored <- generics::tidy(f)
  expect_identical(generics::tidy(f, conf.int = TRUE), stored)
  expect_identical(generics::tidy(f, conf.int = TRUE, conf.level = .9), stored)
  expect_identical(generics::tidy(f, conf.level = .9), stored)
  no.ci <- generics::tidy(f, conf.int = FALSE, exponentiate = TRUE)
  expect_false(any(c("conf.low", "conf.high") %in% names(no.ci)))
  expect_equal(no.ci$estimate, exp(stored$estimate), tolerance = 0)
  expect_error(generics::tidy(f, conf.int = .8, conf.level = .9), "different")
  expect_error(generics::tidy(f, conf.int = NA), "conf.int")
  expect_error(generics::tidy(f, exponentiate = 1), "exponentiate")
  expect_error(generics::tidy(f, conf.level = 1), "conf.level")
  f$outcome.type <- "proportional-competing-risk"
  f$boot.method <- list(report.sandwich.conf = FALSE, report.boot.conf = TRUE, boot.bca = TRUE)
  f$vcov_bootstrap <- f$vcov
  expect_equal(generics::tidy(f, conf.int = TRUE)$conf.low, stored$conf.low, tolerance = 0)
  expect_error(generics::tidy(f, conf.level = .8), "BCa")
})

test_that("tbl_regression selects the event, contrast and effect label automatically", {
  skip_if_not_installed("gtsummary", "2.6.1")
  skip_if_not_installed("broom.helpers", "1.20.0")
  for (levels in list(c("B", "A"), c("B", "A", "C"))) {
    f <- reporting_polyreg_fixture(levels)
    f <- unserialize(serialize(f, NULL))
    original <- serialize(f, NULL)
    for (event in c("event1", "event2")) {
      expect_silent(tbl <- gtsummary::tbl_regression(f, exponentiate = TRUE, event = event))
      td <- generics::tidy(f, event = event, exponentiate = TRUE)
      body <- tbl$table_body
      rows <- which(!is.na(body$reference_row) & !body$reference_row)
      reference <- which(!is.na(body$reference_row) & body$reference_row)
      for (nm in c("estimate", "conf.low", "conf.high", "std.error", "p.value"))
        expect_equal(body[[nm]][rows], td[[nm]], tolerance = 0)
      expect_true(all(body$variable == "group"))
      expect_true(all(body$var_label == "Treatment"))
      expect_equal(unname(body$label[reference]), "B")
      expect_true(all(is.na(body$std.error[reference])))
      expect_true(all(is.na(body$p.value[reference])))
      expected <- if (event == "event1") "**RR**" else "**SHR**"
      expect_identical(tbl$table_styling$header$label[tbl$table_styling$header$column == "estimate"], expected)
      ci.header <- tbl$table_styling$header$label[tbl$table_styling$header$column == "conf.low"]
      expect_match(ci.header, "90%")
    }
    expect_identical(serialize(f, NULL), original)
    default <- gtsummary::tbl_regression(f, exponentiate = TRUE)
    explicit <- gtsummary::tbl_regression(f, exponentiate = TRUE, event = "event1")
    expect_identical(default$table_body, explicit$table_body)
    expect_error(gtsummary::tbl_regression(f, event = "both"), "one event")
  }
})

test_that("fitted reporting retains variable labels and categorical numeric exposure", {
  skip_if_not_installed("gtsummary", "2.6.1")
  skip_if_not_installed("broom.helpers", "1.20.0")
  data(diabetes.complications)
  for (numeric.codes in c(FALSE, TRUE)) {
    d <- diabetes.complications
    if (numeric.codes) d$fruitq <- as.integer(d$fruitq)
    attr(d$fruitq, "label") <- "Fruit intake"
    f <- polyreg(Event(t,epsilon)~age+sex, exposure="fruitq", data=d,
      code.exposure.ref=if (numeric.codes) 2L else "Q2", time.point=8,
      effect.measure1="OR", effect.measure2="SHR")
    f <- unserialize(serialize(f,NULL))
    expect_identical(f$reporting$label, "Fruit intake")
    for (event in c("event1", "event2")) {
      tab <- gtsummary::tbl_regression(f, exponentiate=TRUE, event=event)
      td <- generics::tidy(f, event=event, exponentiate=TRUE)
      rows <- which(!is.na(tab$table_body$reference_row) & !tab$table_body$reference_row)
      expect_true(all(tab$table_body$var_label == "Fruit intake"))
      for (nm in c("estimate", "conf.low", "conf.high", "std.error", "p.value"))
        expect_equal(tab$table_body[[nm]][rows], td[[nm]], tolerance=0)
    }
  }
})
