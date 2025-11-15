test_that("type.y passes through engines", {
  skip_if_not_installed("survival")
  sf <- survival::survfit(survival::Surv(c(1, 2, 3, 4), c(1, 0, 1, 0)) ~ 1)
  types <- c("survival", "risk", "cumhaz", "cloglog")
  for (ty in types) {
    p <- cifmodeling::cifplot(sf,
                              type.y = ty,
                              engine = "ggsurvfit",
                              add.conf = FALSE,
                              add.risktable = FALSE,
                              add.estimate.table = FALSE,
                              add.censor.mark = FALSE,
                              add.competing.risk.mark = FALSE,
                              add.intercurrent.event.mark = FALSE,
                              add.quantile = FALSE)
    expect_s3_class(p, "ggplot")
  }
})

test_that("ggplot engine handles cumhaz and cloglog", {
  skip_if_not_installed("survival")
  sf <- survival::survfit(survival::Surv(c(1, 2, 3, 4), c(1, 0, 1, 0)) ~ 1)
  for (ty in c("cumhaz", "cloglog")) {
    p <- cifmodeling::cifplot(sf,
                              type.y = ty,
                              engine = "ggplot",
                              add.conf = FALSE,
                              add.risktable = FALSE,
                              add.estimate.table = FALSE,
                              add.censor.mark = FALSE,
                              add.competing.risk.mark = FALSE,
                              add.intercurrent.event.mark = FALSE,
                              add.quantile = FALSE)
    expect_s3_class(p, "ggplot")
  }
})
