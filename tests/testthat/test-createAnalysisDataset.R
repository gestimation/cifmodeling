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
  t <- cifcurve(Event(t, d)~strata, testdata, weight="w", error = "greenwood", report.ggsurvfit = FALSE)
  e1 <- cifcurve(Event(t, d)~1, testdata1, weight="w", error = "greenwood", report.ggsurvfit = FALSE)
  e2 <- cifcurve(Event(t, d)~1, testdata2, weight="w", error = "greenwood", report.ggsurvfit = FALSE)
  expected <- c(e1$std.err, e2$std.err)
  tested <- t$std.err
  expect_equal(expected, tested)
})
