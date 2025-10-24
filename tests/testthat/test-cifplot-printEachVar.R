test_that("printEachVar returns one plot per RHS var", {
  data(diabetes.complications)
  diabetes.complications$sex_ <- as.factor(diabetes.complications$sex)
  plt <- cifplot(
    Event(t, epsilon) ~ sex_ + fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    printEachVar = TRUE,
    rows.columns.panel = c(1, 2)
  )
  expect_true(inherits(plt, "patchwork") || inherits(plt, "ggplot"))
  pa <- attr(plt, "plots")
  expect_false(is.null(pa))
  expect_length(pa, 2L)
  expect_true(all(vapply(pa, function(p) inherits(p, "ggplot"), logical(1))))
})

test_that("order.strata and label.strata (positional) are respected", {
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    printEachVar = TRUE,
    rows.columns.panel = c(1, 1),
    order.strata = list(fruitq = c("Q1", "Q2", "Q3", "Q4")),
    label.strata = list(fruitq = c("Q1", "Q2", "Q3", "Q4"))
  )
  expect_true(inherits(plt, "patchwork") || inherits(plt, "ggplot"))
})

test_that("label.strata named mapping works with order.strata", {
  data(diabetes.complications)
  plt <- cifplot(
    Event(t, epsilon) ~ fruitq,
    data = diabetes.complications,
    outcome.type = "COMPETING-RISK",
    code.event1 = 1, code.event2 = 2, code.censoring = 0,
    printEachVar = TRUE,
    rows.columns.panel = c(1, 1),
    order.strata = list(fruitq1 = c("Q1", "Q2", "Q3", "Q4")),
    label.strata = list(fruitq1 = c(Q1 = "Q1", Q2 = "Q2", Q3 = "Q3", Q4 = "Q4"))
  )
  expect_true(inherits(plt, "patchwork") || inherits(plt, "ggplot"))
})

test_that("numeric RHS var errors under printEachVar", {
  data(diabetes.complications)
  expect_error(
    cifplot(
      Event(t, epsilon) ~ age,
      data = diabetes.complications,
      outcome.type = "COMPETING-RISK",
      code.event1 = 1, code.event2 = 2, code.censoring = 0,
      printEachVar = TRUE
    ),
    "numeric.*discretize",
    ignore.case = TRUE
  )
})

test_that("order.strata works per-variable when printEachVar = TRUE", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  data(diabetes.complications)
  df <- diabetes.complications

  # ダミーの第2因子を追加（確実に存在させる）
  set.seed(1)
  df$z2 <- factor(sample(c("Female","Male"), size = nrow(df), replace = TRUE))

  ord_list <- list(
    fruitq = c("Q4","Q2","Q1"),     # 一部だけ載せる
    z2     = c("Female","Male")     # 全部載せる
  )

  patch <- cifplot(
    Event(t, epsilon) ~ fruitq + z2,
    data = df,
    outcome.type = "COMPETING-RISK",
    type.y = "risk",
    printEachVar = TRUE,
    order.strata = ord_list,
    addRiskTable = FALSE
  )

  expect_true(
    inherits(patch, c("gg", "ggplot")) ||
      inherits(patch, "patchwork") ||
      inherits(patch, "gtable")
  )
  plots_attr <- attr(patch, "plots")
  expect_false(is.null(plots_attr))
  expect_equal(length(plots_attr), 2L)
  expect_true(all(vapply(plots_attr, function(p) inherits(p, "ggplot"), logical(1))))

  # 1枚目（fruitq）
  sc1_col <- plots_attr[[1]]$scales$get_scales("colour")
  sc1_fil <- plots_attr[[1]]$scales$get_scales("fill")
  expect_identical(sc1_col$limits, ord_list$fruitq)
  expect_identical(sc1_fil$limits, ord_list$fruitq)

  # 2枚目（z2）
  sc2_col <- plots_attr[[2]]$scales$get_scales("colour")
  sc2_fil <- plots_attr[[2]]$scales$get_scales("fill")
  expect_identical(sc2_col$limits, ord_list$z2)
  expect_identical(sc2_fil$limits, ord_list$z2)
})

test_that("order.strata with no overlap issues a warning and is ignored", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  data(diabetes.complications)

  ord_bad <- c("AAA","BBB")

  expect_warning(
    p <- cifplot(
      Event(t, epsilon) ~ fruitq,
      data = diabetes.complications,
      outcome.type = "COMPETING-RISK",
      type.y = "risk",
      order.strata = ord_bad,
      addRiskTable = FALSE
    ),
    regexp = "no overlap", fixed = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  sc_fil <- p$scales$get_scales("fill")

  # オーバーラップ無し → order は無視され、明示的な scale が付かない場合がある
  # よって scale オブジェクトが NULL でも合格にする
  if (!is.null(sc_col)) {
    expect_true(is.null(sc_col$limits) || identical(sc_col$limits, character()))
  }
  if (!is.null(sc_fil)) {
    expect_true(is.null(sc_fil$limits) || identical(sc_fil$limits, character()))
  }
})
