test_that("panel.per.variable returns one plot per RHS var", {
  skip_on_cran()
  data(diabetes.complications)
  diabetes.complications$sex_ <- as.factor(diabetes.complications$sex)
  plt <- cifplot(
    Event(t, epsilon) ~ sex_ + fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    panel.per.variable = TRUE,
    add.risktable = FALSE,
    rows.columns.panel = c(1, 2)
  )
  expect_s3_class(plt, "cifpanel")
  expect_true(inherits(plt$patchwork, "patchwork"))
  expect_length(plt$list.plot, 2L)
  expect_true(all(vapply(plt$list.plot, function(p) inherits(p, "ggplot"), logical(1))))
})

test_that("order.strata and label.strata (positional) are respected", {
  skip()
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    panel.per.variable = TRUE,
    rows.columns.panel = c(1, 1),
    order.strata = list(fruitq = c("Q1", "Q2", "Q3", "Q4")),
    label.strata = list(fruitq = c("Q1", "Q2", "Q3", "Q4"))
  )
  expect_s3_class(plt, "cifpanel")
})

test_that("label.strata named mapping works with order.strata", {
  skip()
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    panel.per.variable = TRUE,
    rows.columns.panel = c(1, 1),
    order.strata = list(fruitq1 = c("Q1", "Q2", "Q3", "Q4")),
    label.strata = list(fruitq1 = c(Q1 = "Q1", Q2 = "Q2", Q3 = "Q3", Q4 = "Q4"))
  )
  expect_s3_class(plt, "cifpanel")
})

test_that("order.strata works per-variable when panel.per.variable = TRUE", {
  skip()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  data(diabetes.complications)
  df <- diabetes.complications

  set.seed(1)
  df$z2 <- factor(sample(c("Female","Male"), size = nrow(df), replace = TRUE))

  ord_list <- list(
    fruitq = c("Q4","Q2","Q1"),
    z2     = c("Female","Male")
  )

  patch <- cifplot(
    Event(t, epsilon) ~ fruitq + z2,
    data = df,
    outcome.type = "competing-risk",
    type.y = "risk",
    panel.per.variable = TRUE,
    order.strata = ord_list,
    add.risktable = FALSE
  )

  expect_s3_class(patch, "cifpanel")
  expect_true(inherits(patch$patchwork, "patchwork"))
  expect_equal(length(patch$list.plot), 2L)
  expect_true(all(vapply(patch$list.plot, function(p) inherits(p, "ggplot"), logical(1))))

  sc1_col <- patch$list.plot[[1]]$scales$get_scales("colour")
  sc1_fil <- patch$list.plot[[1]]$scales$get_scales("fill")
  expect_identical(sc1_col$limits, ord_list$fruitq)
  expect_identical(sc1_fil$limits, ord_list$fruitq)

  sc2_col <- patch$list.plot[[2]]$scales$get_scales("colour")
  sc2_fil <- patch$list.plot[[2]]$scales$get_scales("fill")
  expect_identical(sc2_col$limits, ord_list$z2)
  expect_identical(sc2_fil$limits, ord_list$z2)
})

