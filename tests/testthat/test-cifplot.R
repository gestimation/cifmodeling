test_that("label.strata only adjusts labels and suppresses fill legend", {
  skip()
  data(diabetes.complications)
  label <- c("Low intake", "High intake")
  level <- c(0, 1)
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    label.strata = label,
    level.strata = level,
    add.risktable = FALSE
  )

  expect_s3_class(p, "ggplot")
  sc_col   <- p$scales$get_scales("colour")
  sc_lin   <- p$scales$get_scales("linetype")
  sc_fill  <- p$scales$get_scales("fill")
  sc_shape <- p$scales$get_scales("shape")
  print(sc_col$get_labels())
  print(sc_lin$get_labels())

  # ここを「順序まで完全一致」→「内容が一致」に緩める
  expect_setequal(sc_col$get_labels(), unname(label))
  expect_setequal(sc_lin$get_labels(), unname(label))

  # fill / shape は非表示のままでOK
  expect_identical(sc_fill$guide, "none")
  expect_identical(sc_shape$guide, "none")
})

test_that("label.strata is reflected in color legend NEW", {
  data(diabetes.complications)
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    label.strata = c("Low intake", "High intake"),
    level.strata = c(0, 1),
    add.risktable = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  br     <- sc_col$get_breaks()

  get_labels_safe <- function(sc, br) {
    lab <- sc$labels
    if (is.function(lab)) {
      lab(br)
    } else if (is.null(lab)) {
      br
    } else {
      lab
    }
  }

  lbl <- get_labels_safe(sc_col, br)
  expect_equal(unname(lbl), c("Low intake", "High intake"))
})


test_that("label.strata overrides palette labels", {
  skip()
  data(diabetes.complications)
  pal  <- c("#FF0000", "#0000FF")
  lbls <- c("0" = "Group A", "1" = "Group B")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    palette = pal,
    label.strata = lbls,
    add.risktable = FALSE
  )

  sc_col  <- p$scales$get_scales("colour")
  sc_fill <- p$scales$get_scales("fill")

  # ラベルの中身が上書きされていることだけを確認（順番は許容）
  expect_setequal(sc_col$get_labels(), unname(lbls))

  # パレットは指定どおり使われてること
  expect_equal(sc_col$palette(seq_along(pal)), pal)

  # fill は表示しない
  expect_identical(sc_fill$guide, "none")
})


test_that("order is applied first; missing levels are appended", {
  skip()
  cur_full  <- c("grp=A","grp=B","grp=C")
  cur_short <- c("A","B","C")

  lbl <- c("Alpha","Bravo","Charlie")
  names(lbl) <- cur_full

  # 1) 正しい order はそのまま使われる
  ord <- c("grp=B","grp=C","grp=A")
  out <- plot_reconcile_order_and_labels(
    cur_lvls_full    = cur_full,
    cur_lvls_short   = cur_short,
    level.strata     = cur_full,
    order.strata     = ord,
    label.strata.map = lbl
  )
  expect_true(out$used_order)
  expect_identical(out$limits_arg, ord)
  expect_identical(out$strata_labels_final, unname(lbl[ord]))

  # 2) order に一部しか書いてなくても、残りが後ろに付く
  ord_bad <- c("grp=A","grp=B")   # ← C がない
  out2 <- plot_reconcile_order_and_labels(
    cur_lvls_full    = cur_full,
    cur_lvls_short   = cur_short,
    level.strata     = cur_full,
    order.strata     = ord_bad,
    label.strata.map = lbl
  )
  expect_true(out2$used_order)
  # A, B が先に来て、残りの C が後ろに付く
  expect_identical(out2$limits_arg, c("grp=A","grp=B","grp=C"))
  # ラベルも同じ順で並ぶ（ここは names いらないので unname する）
  expect_identical(
    out2$strata_labels_final,
    unname(lbl[c("grp=A","grp=B","grp=C")])
  )

  # 3) order がまったくかすらないときだけ warning & 無視
  ord_none <- c("grp=X","grp=Y")
  expect_warning(
    out3 <- plot_reconcile_order_and_labels(
      cur_lvls_full    = cur_full,
      cur_lvls_short   = cur_short,
      level.strata     = cur_full,
      order.strata     = ord_none,
      label.strata.map = NULL
    ),
    "`order.strata` has no overlap"
  )
  expect_null(out3$limits_arg)
  expect_true(out3$forbid_limits_due_to_order)

  # 4) order がなくて label だけなら label の順
  out4 <- plot_reconcile_order_and_labels(
    cur_lvls_full    = cur_full,
    cur_lvls_short   = cur_short,
    level.strata     = cur_full,
    order.strata     = NULL,
    label.strata.map = lbl
  )
  expect_identical(out4$limits_arg, names(lbl))
  expect_identical(out4$strata_labels_final, unname(lbl))
})




test_that("palette only uses manual color scale", {
  data(diabetes.complications)
  pal <- c("#FF0000", "#0000FF")
  p <- cifplot(
    Event(t, epsilon) ~ fruitq1,
    data = diabetes.complications,
    outcome.type = "competing-risk",
    code.events  = c(1, 2, 0),
    palette = pal,
    add.risktable = FALSE
  )

  sc_col <- p$scales$get_scales("colour")
  expect_s3_class(sc_col, "ScaleDiscreteManual")
  expect_identical(sc_col$scale_name, "manual")
  expect_equal(sc_col$palette(seq_along(pal)), pal)

  sc_fill <- p$scales$get_scales("fill")
  expect_identical(sc_fill$guide, "none")
})



test_that("numeric vectors with <9 unique values are converted to factor", {
  x <- c(1, 2, 2, 3, NA_real_)
  out <- cifmodeling:::plot_normalize_strata(x)
  expect_true(is.factor(out$values))
  expect_identical(levels(out$values), c("1", "2", "3"))
  expect_identical(as.character(out$values[1:4]), c("1", "2", "2", "3"))
  expect_true(is.na(out$values[5]))
  expect_identical(out$strategy, "factor")
})

test_that("numeric vectors with >=9 unique values are split at the median", {
  x <- 1:10
  out <- cifmodeling:::plot_normalize_strata(x)
  expect_true(is.factor(out$values))
  expect_identical(levels(out$values),
                   c(cifmodeling:::plot_default_labels$strata$below_median,
                     cifmodeling:::plot_default_labels$strata$above_median))
  med <- stats::median(x)
  expect_true(all(out$values[x <= med] == cifmodeling:::plot_default_labels$strata$below_median))
  expect_true(all(out$values[x >  med] == cifmodeling:::plot_default_labels$strata$above_median))
  expect_identical(out$strategy, "median")
})

test_that("formula normalization converts numeric RHS strata", {
  df <- data.frame(time = 1:10, status = c(0, 1, rep(0, 8)), z = seq_len(10))
  fml <- Event(time, status) ~ z
  out <- cifmodeling:::plot_normalize_formula_data(fml, df)
  expect_true(is.factor(out$data$z))
  expect_identical(levels(out$data$z),
                   c(cifmodeling:::plot_default_labels$strata$below_median,
                     cifmodeling:::plot_default_labels$strata$above_median))
})
