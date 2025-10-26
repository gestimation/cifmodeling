test_that("printEachEvent ignored for non-CR outcomes with warning", {
  data(diabetes.complications)
  diabetes.complications$status <- as.integer(diabetes.complications$epsilon > 0)
  expect_warning(
    {
      plt <- cifplot(
        survival::Surv(t, status) ~ fruitq1,
        data = diabetes.complications,
        outcome.type = "SURVIVAL",
        printEachEvent = TRUE
      )
      expect_s3_class(plt, "ggplot")
    },
    "printEachEvent=TRUE is only for COMPETING-RISK"
  )
})

test_that("cifplot(printEachEvent=TRUE) returns a patchwork object", {
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    printEachEvent = TRUE
  )
  expect_true(
    inherits(plt, c("gg", "ggplot")) ||
      inherits(plt, "patchwork") ||
      inherits(plt, "gtable")
  )
  plots_attr <- attr(plt, "plots")
  expect_false(is.null(plots_attr))
  expect_equal(length(plots_attr), 2L)
  expect_true(all(vapply(plots_attr, function(p) inherits(p, "ggplot"), logical(1))))
  expect_identical(
    vapply(plots_attr, function(p) p$labels$y %||% NA_character_, character(1)),
    c(
      "Cumulative incidence of interest",
      "Cumulative incidence of competing risk"
    )
  )
})

test_that("cifplot(printEachEvent=TRUE) returns two panels and passes y labels", {
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.events = c(1, 2, 0),
    printEachEvent = TRUE,
    label.y = "Left axis"
  )

  plots_attr <- attr(plt, "plots")
  expect_true(!is.null(plots_attr))
  expect_length(plots_attr, 2L)

  titles <- vapply(plots_attr, function(p) p$labels$title %||% NA_character_, character(1))
  expect_true(all(is.na(titles)))

  ylabels <- vapply(plots_attr, function(p) p$labels$y %||% NA_character_, character(1))
  expect_equal(ylabels, c("Left axis", "Left axis"))
})
