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

test_that("In createAnalysisDataset() basic subsetting works for logical vector", {
  set.seed(1)
  df <- data.frame(x = 1:6, y = letters[1:6], z = c(1, NA, 3, 4, 5, NA))
  keep <- c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)

  res <- createAnalysisDataset(
    y ~ x,
    data = df,
    subset.condition = keep,
    na.action = na.pass
  )

  expect_s3_class(res, "data.frame")
  expect_identical(colnames(res), c("y", "x"))
  expect_identical(nrow(res), sum(keep, na.rm = TRUE))
})

test_that("In createAnalysisDataset() subset.condition accepts one-sided formula, character, expression, and language", {
  df <- data.frame(x = 1:6, y = letters[1:6], g = c("a","a","b","b","b","a"), stringsAsFactors = FALSE)

  # one-sided formula
  r1 <- createAnalysisDataset(y ~ x, df, subset.condition = ~ g == "b")
  expect_true(all(r1$g == "b"))

  # character
  r2 <- createAnalysisDataset(y ~ x, df, subset.condition = 'g == "a"')
  expect_true(all(r2$g == "a"))

  # expression
  r3 <- createAnalysisDataset(y ~ x, df, subset.condition = expression(x >= 3))
  expect_true(all(r3$x >= 3))

  # language (call)
  r4 <- createAnalysisDataset(y ~ x, df, subset.condition = quote(x <= 2 | g == "b"))
  expect_true(all(r4$x <= 2 | r4$g == "b"))
})

test_that("In createAnalysisDataset() subset.condition validation errors", {
  df <- data.frame(x = 1:4, y = 1:4)

  # logical length mismatch
  expect_error(
    createAnalysisDataset(y ~ x, df, subset.condition = c(TRUE, FALSE)),
    "logical length must equal"
  )

  # non-logical evaluated result
  expect_error(
    createAnalysisDataset(y ~ x, df, subset.condition = ~ x + 1),
    "not logical"
  )

  # unsupported type
  expect_error(
    createAnalysisDataset(y ~ x, df, subset.condition = 123),
    "Unsupported `subset.condition` type"
  )
})

test_that("In createAnalysisDataset() missing columns -> stop or fill with NA depending on fill_missing", {
  df <- data.frame(x = 1:3, y = 2:4)

  # stop when missing
  expect_error(
    createAnalysisDataset(y ~ x + z, df, fill_missing = FALSE),
    "Undefined columns selected"
  )

  # fill when requested
  expect_warning(
    r <- createAnalysisDataset(y ~ x + z, df, fill_missing = TRUE),
    "will be filled with NA"
  )
  expect_true("z" %in% names(r))
  expect_true(all(is.na(r$z)))
})

test_that("In createAnalysisDataset() column order follows formula + other.variables.analyzed order", {
  df <- data.frame(a = 1:3, b = 2:4, c = 5:7, d = 8:10)

  res <- createAnalysisDataset(
    a + c ~ b,
    df,
    other.variables.analyzed = c("d")
  )
  expect_identical(colnames(res), unique(c(all.vars(a + c ~ b), "d")))
})

test_that("In createAnalysisDataset() na.action is applied once at the end", {
  df <- data.frame(x = c(1, 2, NA, 4), y = c("a","b","c","d"))

  r_pass <- createAnalysisDataset(y ~ x, df, na.action = na.pass)
  expect_identical(nrow(r_pass), nrow(df))        # no row dropped

  r_omit <- createAnalysisDataset(y ~ x, df, na.action = na.omit)
  expect_true(nrow(r_omit) < nrow(df))            # rows with NA in used cols dropped
  expect_false(anyNA(r_omit$x))
})
