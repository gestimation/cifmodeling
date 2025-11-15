test_that("polyreg() produced expected coefficients and variance covariance matrix for categorical exposure and binary exposure coded other than the default and stratified in IPCW", {
  skip_on_cran()
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', code.exposure.ref='Q1', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coef,digits=3)
  tested_cov <- round(output$vcov[1,],digits=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.383, -0.276, -0.067, -0.555, -3.991, -0.104, 0.045, -0.169, 0.017, -0.012, -0.012, -0.013, 0.010, -0.008, -0.008, -0.008)
  expect_equal(tested, expected)

  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', data = diabetes.complications, strata = 'strata', effect.measure1='SHR', effect.measure2='SHR', time.point=8, outcome.type='C')
  tested_coefficient <- round(output$coef,digits=3)
  tested_cov <- round(output$vcov[1,],digits=3)
  tested <- as.vector(cbind(tested_coefficient,tested_cov))
  expected <- c(-1.383, 0.363, -3.988, 0.081, 0.009, -0.006, 0.003, -0.003)
  expect_equal(tested, expected)
})

test_that("polyreg() produced expected coefficients and variance covariance matrix from survival data with missing variables", {
  data(diabetes.complications)
  diabetes.complications$d <- as.numeric(diabetes.complications$epsilon>0)
  expected_df <- diabetes.complications[-(1:10), ]
  expected_output <- polyreg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = expected_df, effect.measure1='RR', time.point=8, outcome.type='S')
  expected <- round(expected_output$coef,digits=3)

  diabetes.complications$t[1:2] <- NA
  diabetes.complications$d[3:4] <- NA
  diabetes.complications$sex[5:6] <- NA
  diabetes.complications$fruitq1[7:8] <- NA
  diabetes.complications$strata[9:10] <- NA
  output <- polyreg(nuisance.model = Event(t,d)~sex, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='S')
  tested <- round(output$coef,digits=3)
  expect_equal(tested, expected)
})

test_that("polyreg() produced expected common effects from competing risks data", {
  skip_on_cran()
  data(diabetes.complications)
  output <- polyreg(nuisance.model = Event(t,epsilon)~+1, exposure = 'fruitq1', strata = 'strata', data = diabetes.complications, effect.measure1='RR', effect.measure2='RR', time.point=1:8, outcome.type='PC', report.boot.conf=FALSE)
  tested <- round(output$coef,digits=3)
  expected <- c(-7.225, -4.051, -3.067, -2.534, -2.114, -1.803, -1.572, -1.392, 0.296, -9.319, -7.712, -6.977, -6.087, -5.156, -4.731, -4.366, -4.040, -0.022)
  expect_equal(tested, expected)
})

test_that("polyreg() produced expected SE for common effects from survival data", {
  df_test <- createTestData1(100, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
  df_test$d <- as.numeric(df_test$epsilon>0)
  df_test$a <- as.numeric(df_test$strata)
  output <- polyreg(nuisance.model = Event(t,d) ~ +1, exposure = 'a', data = df_test,
                    effect.measure1='RR', time.point=20, outcome.type='PS', report.boot.conf=TRUE, boot.replications=10, boot.seed=46)
  tested <- round(output$summary$`event 1 (no competing risk`$tidy$std.error,digits=3)
  expected <- c(0.194)
  expect_equal(tested, expected)
})

#test_that("polyreg() produced expected SE for common effects from competing-risks data", {
#  df_test <- createTestData1(100, 2, first_zero=TRUE, last_zero=TRUE, subset_present=FALSE, logical_strata=TRUE, na_strata=FALSE)
#  df_test$a <- as.numeric(df_test$strata)
#  output <- polyreg(nuisance.model = Event(t,epsilon) ~ +1, exposure = 'a', data = df_test,
#                    effect.measure1='RR', effect.measure2='RR', time.point=20, outcome.type='PROPORTIONAL-COMPETING-RISK', report.boot.conf=TRUE, `boot.replications`=10, boot.seed=46)
#  tested <- round(output$summary$event2$tidy$std.error,digits=3)
#  expected <- c(4.886)
#  expect_equal(tested, expected)
#})
