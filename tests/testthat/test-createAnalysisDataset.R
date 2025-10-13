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
