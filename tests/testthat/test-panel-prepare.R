make_min_data <- function() {
  data.frame(
    time = c(1, 2, 3),
    status = c(0, 1, 2),
    A = c(0, 1, 0)
  )
}

make_min_inputs <- function() {
  list(
    K = 1L,
    formulas = list(Event(time, status) ~ A),
    data = make_min_data(),
    code.events = list(c(1, 2, 0)),
    outcome.flags = c("C"),
    outcome.list = list("COMPETING-RISK"),
    typey.list = list("risk"),
    labely.list = list("CIF"),
    typex.list = NULL,
    labelx.list = list("Time"),
    limsx.list = NULL,
    limsy.list = NULL,
    breakx.list = NULL,
    breaky.list = NULL,
    addCI.list = list(TRUE),
    addCen.list = list(TRUE),
    addCR.list = list(FALSE),
    addIC.list = list(FALSE),
    addQ.list = list(FALSE),
    strata.list = NULL,
    legend.position = "top",
    dots = list(style = "CLASSIC", font.family = NULL, font.size = NULL)
  )
}


test_that("panel_prepare returns curves and plot args", {
  skip_on_cran()
  skip_if_not_installed("survival")
  inputs <- make_min_inputs()
  prep <- do.call(cifmodeling:::panel_prepare, inputs)

  expect_true(is.list(prep$curves))
  expect_equal(length(prep$curves), inputs$K)
  expect_true(all(vapply(prep$curves, inherits, logical(1), what = "survfit")))

  expect_true(is.list(prep$plot_args))
  expect_equal(length(prep$plot_args), inputs$K)
  expect_true(all(vapply(prep$plot_args, is.list, logical(1))))
  expect_true(all(vapply(prep$plot_args, function(x) identical(x$legend.position, inputs$legend.position), logical(1))))

  expect_identical(prep$K, inputs$K)
})


# tests/testthat/test-panel-prepare.R
testthat::test_that("panel_prepare() returns per-panel objects with expected structure", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("ggsurvfit")  # cifplot/cifcurve の内部で利用している場合に備えて

  data(diabetes.complications, package = "cifmodeling")

  # 2パネル（event1 / event2）を作る想定
  code.events <- list(c(1, 2, 0), c(2, 1, 0))

  # 内部ヘルパを直接呼ぶ
  out <- cifmodeling:::panel_prepare(
    formula              = Event(t, epsilon) ~ fruitq,
    data                 = diabetes.complications,
    outcome.type         = "COMPETING-RISK",
    code.events          = code.events,
    # 以下は panel_prepare が引き継ぐはずの代表的な引数（存在すれば pass-through を確認）
    label.y              = c("Diabetic retinopathy", "Macrovascular complications"),
    label.x              = "Years from registration",
    title.plot           = c("Diabetic retinopathy", "Macrovascular complications")
  )

  # --- 形の検証（実装が少し変わっても壊れにくいチェック） ---
  testthat::expect_true(is.list(out))
  testthat::expect_equal(length(out), length(code.events))  # パネル数

  # 各パネル要素は list であること
  lapply(out, function(p) testthat::expect_true(is.list(p)))

  # 各パネルに「推定カーブ（survfit 由来 or ggsurvfit 由来）っぽいオブジェクト」が最低1つ含まれる
  has_curve <- vapply(
    out,
    function(p) {
      any(vapply(p, function(x) inherits(x, c("survfit", "ggsurvfit", "tbl_survfit")), logical(1)))
    },
    logical(1)
  )
  testthat::expect_true(all(has_curve))

  # パネルごとのメタ情報（ラベル等）が pass-through されているなら存在をゆるく確認
  # （実装により格納先キーが異なる可能性があるため、キー名は緩めにチェック）
  maybe_label_keys <- c("label.y", "label_x", "label.x", "title.plot", "plot_args", "args")
  testthat::expect_true(any(names(out[[1]]) %in% maybe_label_keys))
  testthat::expect_true(any(names(out[[2]]) %in% maybe_label_keys))
})

testthat::test_that("cifpanel() smoke test with competing risks (no plotting)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("ggsurvfit")

  data(diabetes.complications, package = "cifmodeling")

  out <- cifpanel(
    title.panel      = "A comparison of cumulative incidence of competing events",
    rows.columns.panel = c(1, 2),
    formula          = Event(t, epsilon) ~ fruitq,
    data             = diabetes.complications,
    outcome.type     = "COMPETING-RISK",
    code.events      = list(c(1, 2, 0), c(2, 1, 0)),
    label.y          = c("Diabetic retinopathy", "Macrovascular complications"),
    label.x          = "Years from registration",
    subtitle.panel   = "Stratified by fruit intake",
    caption.panel    = "Data: diabetes.complications",
    title.plot       = c("Diabetic retinopathy", "Macrovascular complications"),
    legend.position  = "bottom",
    legend.collect   = TRUE,
    print.panel      = FALSE   # 重要：描画せずにオブジェクトだけ返す
  )

  testthat::expect_true(is.list(out))
  # 以前の実装メモに合わせて、代表キーの存在を穏当チェック
  testthat::expect_true(any(c("plots", "out_patchwork") %in% names(out)))
})
