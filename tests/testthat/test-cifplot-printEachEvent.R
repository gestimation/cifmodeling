test_that("panel.per.event ignored for non-CR outcomes with warning", {
  skip_on_cran()
  data(diabetes.complications)
  diabetes.complications$status <- as.integer(diabetes.complications$epsilon > 0)
  expect_warning(
    {
      plt <- cifplot(
        survival::Surv(t, status) ~ fruitq1,
        data           = diabetes.complications,
        outcome.type   = "survival",
        add.risktable   = FALSE,
        panel.per.event = TRUE
      )
      expect_s3_class(plt, "cifplot")
      expect_s3_class(plt$plot, "ggplot")
    },
    "panel.per.event=TRUE is only for COMPETING-RISK"
  )
})

test_that("cifplot(panel.per.event=TRUE) returns a patchwork object", {
  skip_on_cran()
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data           = diabetes.complications,
    outcome.type   = "competing-risk",
    code.events    = c(1, 2, 0),
    add.risktable   = FALSE,
    panel.per.event = TRUE
  )
  expect_s3_class(plt, "cifpanel")
  expect_true(inherits(plt$patchwork, "patchwork"))
  expect_equal(length(plt$list.plot), 2L)
  expect_true(all(vapply(plt$list.plot, function(p) inherits(p, "ggplot"), logical(1))))
  expect_identical(
    vapply(plt$list.plot, function(p) p$labels$y %||% NA_character_, character(1)),
    c(
      "Cumulative incidence of interest",
      "Cumulative incidence of competing risk"
    )
  )
})

test_that("cifplot(panel.per.event=TRUE) returns two panels and passes y labels", {
  skip()
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events = c(1, 2, 0),
    panel.per.event = TRUE,
    label.y = "Left axis"
  )

  expect_s3_class(plt, "cifpanel")
  expect_length(plt$list.plot, 2L)

  titles <- vapply(plt$list.plot, function(p) p$labels$title %||% NA_character_, character(1))
  expect_true(all(is.na(titles)))

  ylabels <- vapply(plt$list.plot, function(p) p$labels$y %||% NA_character_, character(1))
  expect_equal(ylabels, c("Left axis", "Left axis"))
})
