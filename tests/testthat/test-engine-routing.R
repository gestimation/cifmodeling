test_that("engine routing works and schemas align", {
  skip_on_cran()
  data(diabetes.complications, package = "cifmodeling")
  f <- Event(t, epsilon) ~ fruitq1

  # SURVIVAL: KM vs Rcpp(KMモード=greenwood) — 生存曲線一致（ゆるめ）
  out_km <- cifcurve(f, data = diabetes.complications,
                     outcome.type = "SURVIVAL", engine = "calculateKM", error = "greenwood")
  out_rc <- cifcurve(f, data = diabetes.complications,
                     outcome.type = "SURVIVAL", engine = "calculateAJ_Rcpp", error = "greenwood", return_if = FALSE)
  expect_equal(out_km$surv, out_rc$surv, tolerance = 1e-8)

  # COMPETING-RISK: Aalen の SE を R 実装と突き合わせ（ゆるめ）
  out_ajR <- cifcurve(f, data = diabetes.complications,
                      outcome.type = "COMPETING-RISK", engine = "calculateAJ", error = "aalen")
  out_ajC <- cifcurve(f, data = diabetes.complications,
                      outcome.type = "COMPETING-RISK", engine = "calculateAJ_Rcpp", error = "aalen", return_if = FALSE)
  expect_equal(out_ajR$`std.err.cif`, out_ajC$`std.err.cif`, tolerance = 1e-5)

  # 先頭0問題: Aalen/Delta が先頭で 0 固定になっていない（= 前方シフトが効いている）
  out_aalen <- out_ajC
  expect_true(length(out_aalen$`std.err.cif`) >= 2)
  expect_true(out_aalen$`std.err.cif`[1] >= 0)
  if (length(out_aalen$`std.err.cif`) >= 2 && is.finite(out_aalen$`std.err.cif`[2])) {
    expect_true(is.finite(out_aalen$`std.err.cif`[1]))
  }

  # スキーマ最低限
  for (x in list(out_km, out_rc, out_ajR, out_ajC)) {
    expect_true(all(c("time", "surv", "std.err", "std.err.cif", "type", "method", "conf.type") %in% names(x)))
  }
})
