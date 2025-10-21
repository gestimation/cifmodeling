test_that("createAnalysisDataset() produced expected variables with missing data", {
  data(diabetes.complications)
  nuisance.model <- t~fruitq1

  other.variables.analyzed <- NULL
  all_vars <- c(all.vars(nuisance.model), other.variables.analyzed)
  expected <- diabetes.complications[, all_vars, drop = FALSE]
  expected <- expected[-(1), ]
  expected <- expected$t
  diabetes.complications$t[1] <- NA
  tested <- createAnalysisDataset(formula=nuisance.model, data=diabetes.complications, other.variables.analyzed=NULL, subset.condition=NULL, na.action=na.omit)
  tested <- tested$t
  expect_equal(tested, expected)
})

test_that("createAnalysisDataset() produced expected a subset dataset of men", {
  data(diabetes.complications)
  nuisance.model <- t~fruitq1

  other.variables.analyzed <- "sex"
  all_vars <- c(all.vars(nuisance.model), other.variables.analyzed)
  expected <- diabetes.complications[, all_vars, drop = FALSE]
  subset.condition="(diabetes.complications$sex == 1)"
  expected <- subset(expected, eval(parse(text = subset.condition)))
  tested <- createAnalysisDataset(formula=nuisance.model, data=diabetes.complications, other.variables.analyzed="sex", subset.condition="(diabetes.complications$sex == 1)", na.action=na.omit)
  expect_equal(tested, expected)
})

test_that("Greenwood standard error of cifcurve() yields the same outputs as separate analysis when strata is present", {
  testdata <- createTestData(20, 1, first_zero=TRUE, last_zero=FALSE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  testdata1 <- subset(testdata, strata==FALSE)
  testdata2 <- subset(testdata, strata==TRUE)
  t <- cifcurve(Event(t, d)~strata, testdata, outcome.type="S", weight="w", error = "greenwood")
  e1 <- cifcurve(Event(t, d)~1, testdata1, outcome.type="S", weight="w", error = "greenwood")
  e2 <- cifcurve(Event(t, d)~1, testdata2, outcome.type="S", weight="w", error = "greenwood")
  expected <- c(e1$std.err, e2$std.err)
  tested <- t$std.err
  expect_equal(expected, tested)
})

test_that("cifcurve() produced the same outputs as survfit() in survival in survival data", {
  testthat::skip_if_not_installed("survival")
  testdata <- createTestData(20, 1, first_zero=FALSE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  e <- survival::survfit(Surv(t, d)~strata, testdata, weight=w, conf.type = "none")
  t <- cifcurve(Event(t, d) ~ strata, data = testdata, weight="w", conf.type = "none", outcome.type = "SURVIVAL")
  #  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, round(e$std.err,digit=5), e$strata))
  #  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, round(t$std.err,digit=5), t$strata))
  expected <- as.numeric(c(e$time, round(e$surv,digit=5), e$n, e$n.risk, e$n.event, e$n.censor, e$lower, e$strata))
  tested <- as.numeric(c(t$time, round(t$surv,digit=5), t$n, t$n.risk, t$n.event, t$n.censor, t$lower, t$strata))
  expect_equal(expected, tested)
})

test_that("cifcurve() produced the same estimates as cif() in mets in competing risks data", {
  testthat::skip_if_not_installed("mets")
  library(mets)
  data(diabetes.complications)
  cif_fit <- mets::cif(Event(t,epsilon) ~ +1, data=diabetes.complications, cause=1)
  surv_fit <- cifcurve(Event(t, epsilon) ~ +1, data = diabetes.complications, conf.type = "none", outcome.type = "C", error="delta")
  index<-c(1,3,5,6,8,10,12,13,15,16,17,18)
  expected <- round(1-surv_fit$surv[index],digit=5)
  tested <- round(cif_fit$mu[1:12],digit=5)
  expect_equal(expected, tested)
})
