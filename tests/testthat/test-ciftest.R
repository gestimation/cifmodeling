# tests/testthat/test-ciftest-logrank-components.R
testthat::local_edition(3)

testthat::test_that("the public test resolver follows outcome defaults and presets", {
  survival_data <- data.frame(
    time = 1:8,
    status = c(1, 0, 1, 0, 1, 1, 0, 1),
    group = factor(rep(c("A", "B"), 4))
  )
  competing_data <- transform(
    survival_data,
    status = c(1, 0, 2, 0, 1, 2, 0, 1)
  )

  survival_fit <- ciftest(Event(time, status) ~ group, survival_data)
  competing_fit <- ciftest(Event(time, status) ~ group, competing_data)
  early_fit <- ciftest(
    Event(time, status) ~ group, competing_data, test = "early"
  )
  late_fit <- ciftest(
    Event(time, status) ~ group, competing_data, test = "late"
  )
  multiple_fit <- ciftest(
    Event(time, status) ~ group, survival_data, test = "m"
  )

  testthat::expect_identical(survival_fit$test, "logrank")
  testthat::expect_identical(competing_fit$test, "augmented")
  for (alias in c("L", "LR", "log-rank")) {
    spec <- ciftest_resolve_test_spec("survival", test = alias)
    testthat::expect_identical(spec$test, "logrank")
  }
  testthat::expect_identical(
    ciftest_resolve_test_spec("competing-risk", test = "G")$test,
    "gray"
  )
  for (alias in c("A", "aug", "augmentation")) {
    spec <- ciftest_resolve_test_spec("competing-risk", test = alias)
    testthat::expect_identical(spec$test, "augmented")
  }
  logrank_alias_fit <- ciftest(
    Event(time, status) ~ group, survival_data, test = "LR"
  )
  testthat::expect_equal(
    logrank_alias_fit$statistic, survival_fit$statistic,
    tolerance = 1e-12
  )
  testthat::expect_identical(c(early_fit$rho, early_fit$gamma), c(1, 0))
  testthat::expect_identical(c(late_fit$rho, late_fit$gamma), c(0, 1))
  testthat::expect_s3_class(multiple_fit, "ciftest_mdir")
  testthat::expect_identical(
    unname(as.matrix(multiple_fit$directions[, c("rho", "gamma")])),
    matrix(c(2, 0, 0, 0, 2, 0), ncol = 2)
  )
  testthat::expect_error(
    ciftest(Event(time, status) ~ group, competing_data,
            test = "early", rho = 0.5),
    "fixed preset"
  )
})

testthat::test_that("survival score uses individual robust score covariance", {
  df <- data.frame(
    time = 1:10,
    status = c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1),
    group = factor(rep(c("A", "B"), 5)),
    block = factor(rep(c("X", "Y"), each = 5))
  )
  fit <- ciftest(
    Event(time, status) ~ group, df,
    outcome.type = "S", test = "score", strata = "block"
  )
  testthat::expect_s3_class(fit, "ciftest_score")
  testthat::expect_equal(fit$vcov.score, crossprod(fit$score.iid))
  testthat::expect_equal(colSums(fit$score.iid), fit$score)
  testthat::expect_identical(fit$variance.method, "score-iid")
})

testthat::test_that("competing-risk score permits censoring strata", {
  df <- data.frame(
    time = 1:12,
    status = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    group = factor(rep(c("A", "B"), 6)),
    censor = factor(rep(c("L", "H"), each = 6))
  )
  fit <- ciftest(
    Event(time, status) ~ group, df,
    test = "score", strata.censor = "censor"
  )
  testthat::expect_s3_class(fit, "ciftest_score")
  testthat::expect_identical(fit$strata.censor.info$columns, "censor")
  testthat::expect_length(fit$strata.competing.risk.info$columns, 0L)
})

testthat::test_that("probability truncation reports raw and used nuisance support", {
  df <- data.frame(
    time = 1:10,
    status = c(0, 0, 2, 2, 1, 1, 2, 1, 1, 1),
    group = factor(rep(c("A", "B"), 5))
  )
  fit <- ciftest(
    Event(time, status) ~ group, df,
    test = "score", prob.truncation = 0.8
  )
  diagnostic <- fit$diagnostics$truncation
  testthat::expect_true(diagnostic$requested)
  testthat::expect_true(diagnostic$applied)
  testthat::expect_gt(diagnostic$censoring.count, 0L)
  testthat::expect_lt(diagnostic$censoring.minimum.raw,
                      diagnostic$censoring.minimum.used)
  testthat::expect_identical(fit$diagnostics$score.engine, "R")
  testthat::expect_error(
    ciftest(Event(time, status) ~ group, df,
            test = "gray", prob.truncation = 0.8),
    "only for competing-risk score and augmented"
  )
})

testthat::test_that("multigroup augmented score is an omnibus test", {
  df <- data.frame(
    time = rep(1:6, 3),
    status = rep(c(1, 2, 0, 1, 2, 0), 3),
    group = factor(rep(c("A", "B", "C"), each = 6))
  )
  fit <- ciftest(Event(time, status) ~ group, df)
  testthat::expect_s3_class(fit, "ciftest_augmented")
  testthat::expect_identical(ncol(fit$score.iid), 2L)
  testthat::expect_true(unname(fit$parameter) %in% 1:2)
  testthat::expect_true(is.finite(fit$p.value))
})

testthat::test_that("multiple-direction wrapper preserves one-direction tests", {
  survival_data <- data.frame(
    time = 1:10,
    status = c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1),
    group = factor(rep(c("A", "B"), 5))
  )
  scalar <- ciftest(
    Event(time, status) ~ group, survival_data, test = "logrank"
  )
  combined <- ciftest_mdir(
    Event(time, status) ~ group, survival_data,
    directions = "unweighted", test = "logrank"
  )
  testthat::expect_s3_class(combined, "ciftest_mdir")
  testthat::expect_equal(combined$statistic, scalar$statistic,
                         tolerance = 1e-10)
  testthat::expect_equal(combined$vcov.score, scalar$vcov.score,
                         ignore_attr = TRUE, tolerance = 1e-10)

  multiple <- ciftest_mdir(
    Event(time, status) ~ group, survival_data,
    directions = c("early", "late", "unweighted")
  )
  testthat::expect_identical(ncol(multiple$score.iid), 3L)
  testthat::expect_equal(
    unname(as.matrix(multiple$directions[, c("rho", "gamma")])),
    matrix(c(2, 0, 0, 0, 2, 0), ncol = 2)
  )
  testthat::expect_true(unname(multiple$parameter) %in% 1:3)
  testthat::expect_true(is.finite(multiple$p.value))
})

testthat::test_that("multiple-direction classical covariances preserve FH spans", {
  survival_data <- data.frame(
    time = rep(1:6, 3),
    status = c(
      1, 1, 0, 1, 0, 1,
      1, 0, 1, 1, 1, 0,
      0, 1, 1, 0, 1, 1
    ),
    group = factor(rep(c("A", "B", "C"), each = 6))
  )
  competing_data <- transform(
    survival_data,
    status = c(
      1, 2, 0, 1, 2, 1,
      1, 0, 2, 1, 1, 2,
      2, 1, 1, 0, 2, 1
    )
  )

  check_span <- function(data, test) {
    two_directions <- rbind(
      early1 = c(rho = 1, gamma = 0),
      late1 = c(rho = 0, gamma = 1)
    )
    three_directions <- rbind(
      unweighted = c(rho = 0, gamma = 0),
      early1 = c(rho = 1, gamma = 0),
      late1 = c(rho = 0, gamma = 1)
    )
    two <- ciftest_mdir(
      Event(time, status) ~ group, data,
      directions = two_directions, test = test
    )
    three <- ciftest_mdir(
      Event(time, status) ~ group, data,
      directions = three_directions, test = test
    )
    scalar <- ciftest(
      Event(time, status) ~ group, data,
      test = test, rho = 1, gamma = 0
    )

    testthat::expect_equal(three$statistic, two$statistic,
                           tolerance = 1e-10)
    testthat::expect_equal(three$p.value, two$p.value,
                           tolerance = 1e-10)
    testthat::expect_identical(unname(three$parameter), 4L)
    testthat::expect_equal(
      two$vcov.score[1:2, 1:2, drop = FALSE],
      scalar$vcov.score,
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
    two
  }

  logrank <- check_span(survival_data, "logrank")
  gray <- check_span(competing_data, "gray")
  testthat::expect_identical(
    logrank$diagnostics$covariance.source,
    "joint hypergeometric"
  )
  testthat::expect_identical(
    gray$diagnostics$covariance.source,
    "joint Gray"
  )
})

testthat::test_that("ciftest survival UI returns an htest-compatible result", {
  testthat::skip_if_not_installed("survival")
  df <- data.frame(
    time = c(1, 2, 2, 3, 4, 5, 6, 7),
    status = c(1, 0, 1, 0, 1, 1, 0, 1),
    group = factor(c("A", "A", "B", "B", "A", "B", "A", "B"))
  )

  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    outcome.type = "survival"
  )
  ref <- survival::survdiff(survival::Surv(time, status) ~ group, data = df)

  testthat::expect_s3_class(fit, "ciftest")
  testthat::expect_s3_class(fit, "htest")
  testthat::expect_s3_class(fit, "survdiff")
  testthat::expect_s3_class(fit$survdiff, "survdiff")
  testthat::expect_equal(unname(fit$statistic), unname(ref$chisq), tolerance = 1e-10)
  testthat::expect_identical(fit$variance.method, "hypergeometric")
  testthat::expect_true(is.matrix(fit$score.iid))
  testthat::expect_identical(dim(fit$score.iid), c(nrow(df), 1L))
  testthat::expect_identical(nobs(fit), nrow(df))
})

testthat::test_that("ciftest applies subset and missingness before automatic outcome detection", {
  df <- data.frame(
    time = 1:7,
    status = c(0, 1, 2, 0, 1, 2, 0),
    group = factor(c("A", "B", "A", "B", "A", "B", "A")),
    keep = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
    weight = c(1, 1, 1, NA, 1, 1, 1)
  )

  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    weights = "weight",
    subset.condition = ~ keep,
    outcome.type = NULL
  )

  testthat::expect_identical(fit$outcome.type, "survival")
  testthat::expect_identical(fit$n, 4L)
  testthat::expect_identical(fit$diagnostics$analysis.row.index, c(1L, 2L, 5L, 7L))

  augmented <- generics::augment(fit)
  testthat::expect_identical(augmented$.analysis_included,
                             c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE))
  testthat::expect_true(is.matrix(augmented$.score_iid))
  testthat::expect_identical(dim(augmented$.score_iid), c(nrow(df), 1L))
})

testthat::test_that("ciftest defaults to augmentation and retains standard Gray", {
  df <- data.frame(
    time = 1:6,
    status = c(0, 1, 2, 0, 1, 2),
    group = factor(c("A", "B", "A", "B", "A", "B"))
  )

  augmented <- ciftest(Event(time, status) ~ group, data = df)
  testthat::expect_s3_class(augmented, "ciftest")
  testthat::expect_true(augmented$augmentation)
  testthat::expect_identical(
    augmented$method,
    "Closed-form augmented Fine-Gray score test"
  )
  testthat::expect_identical(augmented$iteration, 0L)
  testthat::expect_identical(augmented$variance.method, "score-iid")
  fit <- ciftest(Event(time, status) ~ group, data = df, augmentation = FALSE)
  testthat::expect_s3_class(fit, "ciftest")
  testthat::expect_identical(fit$method, "Gray's test")
  testthat::expect_identical(fit$variance.method, "gray")
  testthat::expect_true(is.finite(unname(fit$statistic)))
  testthat::expect_identical(unname(fit$parameter), 1L)
  testthat::expect_error(
    ciftest(Event(time, status) ~ group, data = df,
            augmentation = FALSE, iteration = 1),
    "requires `test = \"augmented\"`"
  )
})

testthat::test_that("standard Gray matches frozen cmprsk 2.2.12 fixtures", {
  fixture_environment <- new.env(parent = baseenv())
  sys.source(
    testthat::test_path("fixtures", "gray_cmprsk_fixtures.R"),
    envir = fixture_environment
  )

  for (fixture in fixture_environment$gray_cmprsk_fixtures) {
    df <- data.frame(
      time = fixture$time,
      status = fixture$status,
      group = factor(fixture$group, levels = unique(fixture$group))
    )
    fit <- ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = FALSE,
      rho = fixture$rho
    )

    testthat::expect_equal(
      unname(fit$statistic), fixture$statistic,
      tolerance = 2e-12,
      info = fixture$id
    )
    testthat::expect_equal(
      fit$p.value, fixture$p.value,
      tolerance = 2e-12,
      info = fixture$id
    )
    testthat::expect_identical(unname(fit$parameter), fixture$df)
  }
})

testthat::test_that("Gray score and covariance match an independent slow reference", {
  df <- data.frame(
    time = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 9, 10),
    status = c(1L, 2L, 0L, 1L, 2L, 1L, 0L, 2L, 1L, 0L,
               2L, 1L, 0L, 1L, 2L, 1L, 0L, 2L),
    group = factor(
      c("A", "B", "C", "A", "B", "C", "A", "B", "C",
        "A", "B", "C", "A", "B", "C", "A", "B", "C"),
      levels = c("A", "B", "C")
    )
  )

  for (parameters in list(c(0, 0), c(0.5, 0), c(0.25, 0.75))) {
    production <- cifmodeling:::calculate_gray(
      t = df$time,
      epsilon = df$status,
      exposure = "group",
      weights = rep(1, nrow(df)),
      strata = rep(1L, nrow(df)),
      data = df,
      rho = parameters[1],
      gamma = parameters[2]
    )
    reference <- gray_reference_slow(
      time = df$time,
      status = df$status,
      group = df$group,
      rho = parameters[1],
      gamma = parameters[2]
    )

    testthat::expect_equal(
      unname(production$score), reference$score,
      tolerance = 2e-12,
      info = paste("rho/gamma", paste(parameters, collapse = "/"))
    )
    testthat::expect_equal(
      unname(production$var), unname(reference$var),
      tolerance = 2e-12,
      info = paste("rho/gamma", paste(parameters, collapse = "/"))
    )
  }
})

testthat::test_that("Gray integer frequency weights equal row replication", {
  df <- data.frame(
    time = c(1, 2, 2, 3, 4, 5, 5, 6, 7),
    status = c(1L, 0L, 2L, 1L, 2L, 0L, 1L, 1L, 0L),
    group = factor(c("A", "A", "B", "B", "A", "B", "A", "B", "A"))
  )
  frequency <- c(1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L)

  weighted <- cifmodeling:::calculate_gray(
    df$time, df$status, "group",
    weights = frequency,
    strata = rep(1L, nrow(df)),
    data = df
  )
  replicated <- df[rep(seq_len(nrow(df)), frequency), , drop = FALSE]
  expanded <- cifmodeling:::calculate_gray(
    replicated$time, replicated$status, "group",
    weights = rep(1, nrow(replicated)),
    strata = rep(1L, nrow(replicated)),
    data = replicated
  )

  testthat::expect_equal(weighted$score, expanded$score, tolerance = 2e-12)
  testthat::expect_equal(weighted$var, expanded$var, tolerance = 2e-12)
})

testthat::test_that("standard Gray reports unsupported weighting boundaries", {
  df <- data.frame(
    time = 1:6,
    status = c(0L, 1L, 2L, 0L, 1L, 2L),
    group = factor(rep(c("A", "B"), 3)),
    censor_stratum = factor(rep(c("X", "Y"), 3)),
    analysis_stratum = factor(c("X", "X", "Y", "Y", "X", "X"))
  )

  testthat::expect_error(
    ciftest(Event(time, status) ~ group, df,
            weights = rep(0.5, nrow(df)), augmentation = FALSE),
    "integer frequency weights"
  )
  testthat::expect_error(
    ciftest(Event(time, status) ~ group, df,
            augmentation = FALSE, strata.censor = "censor_stratum"),
    "only for competing-risk score and augmented"
  )
  stratified <- ciftest(
    Event(time, status) ~ group, df,
    test = "gray", strata = "analysis_stratum"
  )
  testthat::expect_identical(stratified$strata.info$columns,
                             "analysis_stratum")
})

testthat::test_that("standard Gray survives an unavailable diagnostic score iid", {
  df <- data.frame(
    time = c(0.2, 0.5, 1, 0.3, 2, 3),
    status = c(1L, 2L, 0L, 1L, 1L, 0L),
    group = factor(
      c("A", "A", "A", "B", "B", "B"),
      levels = c("A", "B")
    )
  )

  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = FALSE
  )

  testthat::expect_s3_class(fit, "ciftest")
  testthat::expect_true(is.finite(unname(fit$statistic)))
  testthat::expect_false(fit$diagnostics$score.iid.available)
  testthat::expect_match(
    fit$diagnostics$score.iid.error,
    "Censoring positivity"
  )
  testthat::expect_true(all(is.na(fit$score.iid)))
  testthat::expect_identical(fit$variance.method, "gray")
})

testthat::test_that("ciftest validates formula and control arguments", {
  df <- data.frame(
    time = 1:6,
    status = c(0, 1, 0, 1, 0, 1),
    g1 = factor(rep(c("A", "B"), 3)),
    g2 = factor(rep(c("X", "Y", "X"), 2))
  )

  testthat::expect_error(
    ciftest(Event(time, status) ~ g1 + g2, df, outcome.type = "survival"),
    "one untransformed grouping variable"
  )
  testthat::expect_error(
    ciftest(Event(time, status) ~ g1, df, outcome.type = "survival", rho = -1),
    "`rho`"
  )
  testthat::expect_error(
    ciftest(Event(time, status) ~ g1, df,
            outcome.type = "survival", iteration = -1),
    "non-negative integer"
  )
  testthat::expect_error(
    ciftest(Event(time, status) ~ g1, df,
            outcome.type = "survival", iteration = 1.5),
    "non-negative integer"
  )
  testthat::expect_error(
    ciftest(Event(time, status) ~ g1, df,
            outcome.type = "survival", iteration = TRUE),
    "non-negative integer"
  )
})

testthat::test_that("calculate_log_rank() matches survdiff for unweighted log-rank (rho=0,gamma=0)", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,4,5,6,6,7),
    status = c(1,0,1,2,1,0,1,1,0,0), # 1=event1, 2=other (treated as censor here), 0=censor
    trt    = factor(c("A","A","B","B","A","C","B","C","C","A"), levels = c("A","B","C")),
    str    = factor(c("s1","s1","s1","s1","s2","s2","s2","s2","s1","s2"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)   # treat status==2 as censor for cause1
  w   <- rep(1, nrow(df))

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = 0, gamma = 0
  )

  # survdiff comparison (gamma not supported; gamma=0 here)
  ev1 <- as.integer(df$status == 1)
  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt, data = df, rho = 0)

  # chisq from our score/var (df=K-1)
  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)

  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-8)

  # also check O-E matches our score components for non-ref groups
  oe <- as.numeric(sdiff$obs - sdiff$exp)
  names(oe) <- names(sdiff$obs)
  # our score corresponds to groups except reference (first level)
  testthat::expect_equal(unname(comp$score), unname(oe[-1]), tolerance = 1e-8)
})

testthat::test_that("matches survdiff for FH rho>0 (gamma=0)", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,4,5,6,6,7),
    status = c(1,0,1,2,1,0,1,1,0,0),
    trt    = factor(c("A","A","B","B","A","C","B","C","C","A"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))
  ev1 <- as.integer(df$status == 1)

  rho <- 1
  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = rho, gamma = 0
  )
  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt, data = df, rho = rho)

  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)
  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-6)
})

testthat::test_that("stratified version matches survdiff with strata()", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,4,5,6,6,7),
    status = c(1,0,1,2,1,0,1,1,0,0),
    trt    = factor(c("A","A","B","B","A","C","B","C","C","A"), levels = c("A","B","C")),
    str    = factor(c("s1","s1","s1","s1","s2","s2","s2","s2","s1","s2"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))
  ev1 <- as.integer(df$status == 1)

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = df$str,
    data = df, exposure = "trt",
    rho = 0, gamma = 0
  )

  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt + survival::strata(str), data = df, rho = 0)
  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)

  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-8)
})

testthat::test_that("reference level via code.exposure.ref is respected", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,4,5,6,6,7),
    status = c(1,0,1,1,1,1,0,0),
    trt    = factor(c("A","A","B","A","C","B","C","A"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))
  ev1 <- as.integer(df$status == 1)

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    code.exposure.ref = "B",
    rho = 0, gamma = 0,
    prefix = "a"
  )

  testthat::expect_equal(comp$ref, "B")
  testthat::expect_true(all(grepl("^a_", names(comp$score))))

  # Compare with survdiff on relevelled trt
  df2 <- df
  df2$trt <- stats::relevel(df2$trt, ref = "B")
  sdiff <- survival::survdiff(survival::Surv(time, ev1) ~ trt, data = df2, rho = 0)
  chisq_ours <- as.numeric(t(comp$score) %*% solve(comp$var) %*% comp$score)

  testthat::expect_equal(chisq_ours, as.numeric(sdiff$chisq), tolerance = 1e-8)
})

testthat::test_that("integer weights equal replication (frequency-weight sanity check)", {
  testthat::skip_if_not_installed("survival")

  df <- data.frame(
    time   = c(1,2,2,3,4,5,6,6),
    status = c(1,0,1,2,1,1,1,0),
    trt    = factor(c("A","A","B","B","A","C","B","C"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)

  w_int <- c(1,2,1,1,2,1,1,2)   # integer weights
  comp_w <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w_int,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = 0, gamma = 0
  )

  df_rep <- df[rep(seq_len(nrow(df)), w_int), , drop = FALSE]
  eps_rep <- ifelse(df_rep$status == 1, 1L, 0L)
  comp_rep <- cifmodeling:::calculate_log_rank(
    t = df_rep$time, epsilon = eps_rep, weights = rep(1, nrow(df_rep)),
    strata = rep.int(1L, nrow(df_rep)),
    data = df_rep, exposure = "trt",
    rho = 0, gamma = 0
  )

  testthat::expect_equal(comp_w$score, comp_rep$score, tolerance = 1e-10)
  testthat::expect_equal(comp_w$var,   comp_rep$var,   tolerance = 1e-10)
})

testthat::test_that("gamma>0 runs and returns finite outputs (smoke test)", {
  df <- data.frame(
    time   = c(1,2,2,4,5,6,7,7),
    status = c(1,0,1,1,0,1,1,0),
    trt    = factor(c("A","A","B","A","C","B","C","A"), levels = c("A","B","C"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w   <- rep(1, nrow(df))

  comp <- cifmodeling:::calculate_log_rank(
    t = df$time, epsilon = eps, weights = w,
    strata = rep.int(1L, nrow(df)),
    data = df, exposure = "trt",
    rho = 0, gamma = 1
  )

  testthat::expect_true(all(is.finite(comp$score)))
  testthat::expect_true(all(is.finite(comp$var)))
})

testthat::test_that("errors on one-level exposure", {
  df <- data.frame(
    time = c(1,2,3),
    status = c(1,0,1),
    trt = factor(c("A","A","A"))
  )
  eps <- ifelse(df$status == 1, 1L, 0L)
  w <- rep(1, nrow(df))

  testthat::expect_error(
    cifmodeling:::calculate_log_rank(
      t = df$time, epsilon = eps, weights = w,
      strata = rep.int(1L, nrow(df)),
      data = df, exposure = "trt",
      rho = 0, gamma = 0
    ),
    "only one level|no valid levels|Effect estimation is not possible",
    ignore.case = TRUE
  )
})

testthat::test_that("A12 from calculate_A12_logrank_weightit matches direct finite-difference under linear weights", {
  # Skip gracefully if not yet available
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found in cifmodeling namespace.")
  }
  if (!exists("calculate_log_rank", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_log_rank not found in cifmodeling namespace.")
  }

  set.seed(1)

  n <- 60
  df <- data.frame(
    time = sample(1:8, n, replace = TRUE),
    status_raw = sample(c(0L, 1L, 2L), n, replace = TRUE, prob = c(0.55, 0.35, 0.10)),
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trt = factor(sample(c("A","B","C"), n, replace = TRUE), levels = c("A","B","C"))
  )

  # epsilon normalized for log-rank: 1=event1, others treated as censor (0)
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)

  # Ensure at least one event in each stratum (otherwise score in that stratum is zero and derivative can be degenerate)
  # If not, tweak a couple of records deterministically
  if (sum(epsilon[df$str == "s1"]) == 0L) epsilon[which(df$str == "s1")[1L]] <- 1L
  if (sum(epsilon[df$str == "s2"]) == 0L) epsilon[which(df$str == "s2")[1L]] <- 1L

  t <- df$time
  strata <- df$str

  # ----- Known (linear) weight mechanism -----
  # w(beta) = w0 + X %*% beta, with beta0=0 => w0
  # then dw/dB = X exactly
  q <- 2
  Xw <- matrix(runif(n*q, min = -0.2, max = 0.2), nrow = n, ncol = q)
  colnames(Xw) <- paste0("B", 1:q)

  w0 <- rep(1, n)  # positive baseline weights

  # Mock weightit object with required fields
  weightit_mock <- list(
    weights = w0,
    Mparts = list(dw_dBtreat = Xw)
  )

  rho <- 0
  gamma <- 0
  prob.bound <- 1e-7

  # ----- Run A12 function (under test) -----
  out <- cifmodeling:::calculate_A12_logrank_weightit(
    t = t,
    epsilon = as.integer(epsilon),
    strata = strata,
    data = df,
    exposure = "trt",
    weightit = weightit_mock,
    rho = rho,
    gamma = gamma,
    prob.bound = prob.bound,
    fd_rel_step = 1e-4   # make the step not too tiny for numerical stability
  )

  A12_hat <- out$A12
  steps <- out$fd_steps

  # ----- "Truth" via direct finite difference w.r.t beta using the known mapping w(beta)=w0+Xw beta -----
  # Use the SAME step sizes as calculate_A12... to avoid step-size confounding.
  base <- cifmodeling:::calculate_log_rank(
    t = t,
    epsilon = as.integer(epsilon),
    weights = w0,
    strata = strata,
    data = df,
    exposure = "trt",
    rho = rho,
    gamma = gamma,
    prob.bound = prob.bound
  )

  p <- length(base$score)   # p = K-1
  A12_true <- matrix(0, nrow = p, ncol = q)
  rownames(A12_true) <- names(base$score)
  colnames(A12_true) <- colnames(Xw)

  for (j in seq_len(q)) {
    h <- steps[j]
    if (h == 0) next

    # beta0 +/- h e_j => w +/- h * Xw[,j]
    w_plus  <- w0 + h * Xw[, j]
    w_minus <- w0 - h * Xw[, j]

    s_plus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_plus,
      strata = strata, data = df, exposure = "trt",
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    s_minus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_minus,
      strata = strata, data = df, exposure = "trt",
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    dU_dBj <- (as.numeric(s_plus) - as.numeric(s_minus)) / (2 * h)

    # A12 is on mean-score scale (1/n) dU_total/dB^T
    A12_true[, j] <- dU_dBj / n
  }

  # ----- Assertions -----
  testthat::expect_equal(dim(A12_hat), dim(A12_true))
  testthat::expect_equal(rownames(A12_hat), rownames(A12_true))
  testthat::expect_equal(colnames(A12_hat), colnames(A12_true))

  # Linear weights => should match essentially to machine precision
  testthat::expect_equal(A12_hat, A12_true, tolerance = 1e-10)
})

testthat::test_that("A12 is zero when dw/dB is zero", {
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found in cifmodeling namespace.")
  }

  set.seed(2)
  n <- 40
  df <- data.frame(
    time = sample(1:6, n, replace = TRUE),
    status_raw = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.6, 0.4)),
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trt = factor(sample(c("A","B","C"), n, replace = TRUE), levels = c("A","B","C"))
  )
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)

  w0 <- rep(1, n)
  Xw0 <- matrix(0, nrow = n, ncol = 3)
  colnames(Xw0) <- paste0("B", 1:3)

  weightit_mock <- list(weights = w0, Mparts = list(dw_dBtreat = Xw0))

  out <- cifmodeling:::calculate_A12_logrank_weightit(
    t = df$time,
    epsilon = as.integer(epsilon),
    strata = df$str,
    data = df,
    exposure = "trt",
    weightit = weightit_mock,
    rho = 0,
    gamma = 0,
    fd_rel_step = 1e-4
  )

  testthat::expect_true(all(abs(out$A12) < 1e-14))
})


# --- helper: "true" A12 via direct finite-difference in beta (exact w(beta±h e_j)) ---
A12_true_beta_fd <- function(t, epsilon, strata, data, exposure,
                             beta0, w_map, h_vec,
                             rho = 0, gamma = 0, prob.bound = 1e-7) {
  base <- cifmodeling:::calculate_log_rank(
    t = t, epsilon = as.integer(epsilon), weights = w_map(beta0),
    strata = strata, data = data, exposure = exposure,
    rho = rho, gamma = gamma, prob.bound = prob.bound
  )
  score0 <- as.numeric(base$score)
  p <- length(score0)
  q <- length(beta0)
  n <- length(t)

  A12 <- matrix(0, nrow = p, ncol = q)
  rownames(A12) <- names(base$score) %||% paste0("score", seq_len(p))
  colnames(A12) <- paste0("B", seq_len(q))

  for (j in seq_len(q)) {
    h <- h_vec[j]
    if (!is.finite(h) || h == 0) next

    bp <- beta0; bp[j] <- bp[j] + h
    bm <- beta0; bm[j] <- bm[j] - h

    s_plus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_map(bp),
      strata = strata, data = data, exposure = exposure,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    s_minus <- cifmodeling:::calculate_log_rank(
      t = t, epsilon = as.integer(epsilon), weights = w_map(bm),
      strata = strata, data = data, exposure = exposure,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    dU <- (as.numeric(s_plus) - as.numeric(s_minus)) / (2 * h)
    A12[, j] <- dU / n  # mean-score scale
  }
  A12
}

# --- helper: run A12 function under test and compute error vs truth ---
run_and_error <- function(t, epsilon, strata, data, exposure,
                          weightit_mock, beta0, w_map,
                          fd_rel_step, rho = 0, gamma = 0, prob.bound = 1e-7) {

  out <- cifmodeling:::calculate_A12_logrank_weightit(
    t = t,
    epsilon = as.integer(epsilon),
    strata = strata,
    data = data,
    exposure = exposure,
    weightit = weightit_mock,
    rho = rho, gamma = gamma, prob.bound = prob.bound,
    fd_rel_step = fd_rel_step
  )
  A12_hat <- out$A12
  h_vec <- out$fd_steps

  A12_true <- A12_true_beta_fd(
    t = t, epsilon = epsilon, strata = strata, data = data, exposure = exposure,
    beta0 = beta0, w_map = w_map, h_vec = h_vec,
    rho = rho, gamma = gamma, prob.bound = prob.bound
  )

  err <- max(abs(A12_hat - A12_true))
  step_scale <- max(h_vec^2)

  list(err = err, step2 = step_scale, out = out, A12_true = A12_true)
}

testthat::test_that("Nonlinear exp weights: A12 approximation error scales like O(h^2)", {
  # skip if functions are not available yet
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found.")
  }
  if (!exists("calculate_log_rank", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_log_rank not found.")
  }

  set.seed(10)
  n <- 120
  df <- data.frame(
    time = sample(1:10, n, replace = TRUE),
    status_raw = rbinom(n, 1, 0.35),          # 1=event1, 0=censor
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trt = factor(sample(c("A","B","C"), n, replace = TRUE), levels = c("A","B","C"))
  )
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)
  # ensure at least one event per stratum
  if (sum(epsilon[df$str == "s1"]) == 0L) epsilon[which(df$str == "s1")[1L]] <- 1L
  if (sum(epsilon[df$str == "s2"]) == 0L) epsilon[which(df$str == "s2")[1L]] <- 1L

  t <- df$time
  strata <- df$str
  exposure <- "trt"

  # nonlinear exp weights: w(beta)=exp(X beta)
  q <- 2
  X <- matrix(runif(n*q, min = -0.08, max = 0.08), nrow = n, ncol = q)
  colnames(X) <- paste0("B", 1:q)
  beta0 <- c(0.15, -0.10)

  w_map <- function(beta) as.numeric(exp(X %*% beta))

  w0 <- w_map(beta0)
  dw <- w0 * X  # exact dw/dB at beta0

  weightit_mock <- list(
    weights = w0,
    Mparts = list(dw_dBtreat = dw)
  )

  # two step sizes (halve rel step)
  r1 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 1e-3, rho = 0, gamma = 0)
  r2 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 5e-4, rho = 0, gamma = 0)

  testthat::expect_true(is.finite(r1$err) && is.finite(r2$err))
  testthat::expect_true(r2$err > 0)  # should not be exactly 0 for nonlinear mapping

  # If error is O(h^2), err_ratio ~ step2_ratio
  err_ratio  <- r1$err / r2$err
  step_ratio <- r1$step2 / r2$step2

  # allow slack (numerical noise / piecewise behavior)
  testthat::expect_true(err_ratio / step_ratio > 0.4 && err_ratio / step_ratio < 2.5)

  # and error should decrease when step decreases
  testthat::expect_true(r2$err < r1$err)
})

testthat::test_that("Nonlinear logit-IPW weights (binary): A12 approximation error scales like O(h^2)", {
  if (!exists("calculate_A12_logrank_weightit", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_A12_logrank_weightit not found.")
  }
  if (!exists("calculate_log_rank", envir = asNamespace("cifmodeling"), inherits = FALSE)) {
    testthat::skip("calculate_log_rank not found.")
  }

  set.seed(11)
  n <- 150
  df <- data.frame(
    time = sample(1:12, n, replace = TRUE),
    status_raw = rbinom(n, 1, 0.30),
    str = factor(sample(c("s1", "s2"), n, replace = TRUE)),
    trtbin = factor(rbinom(n, 1, 0.5), levels = c(0,1))  # binary exposure
  )
  epsilon <- ifelse(df$status_raw == 1L, 1L, 0L)
  if (sum(epsilon[df$str == "s1"]) == 0L) epsilon[which(df$str == "s1")[1L]] <- 1L
  if (sum(epsilon[df$str == "s2"]) == 0L) epsilon[which(df$str == "s2")[1L]] <- 1L

  t <- df$time
  strata <- df$str
  exposure <- "trtbin"
  A <- as.integer(df$trtbin == 1L)

  q <- 2
  X <- matrix(runif(n*q, min = -0.12, max = 0.12), nrow = n, ncol = q)
  colnames(X) <- paste0("B", 1:q)
  beta0 <- c(0.10, -0.05)

  # define propensity and IPW weights
  p_map <- function(beta) as.numeric(stats::plogis(X %*% beta))

  # keep away from 0/1 but choose beta0, X so it should not bind
  clip <- function(p, lo = 1e-3, hi = 1 - 1e-3) pmin(pmax(p, lo), hi)

  w_map <- function(beta) {
    p <- clip(p_map(beta))
    A / p + (1 - A) / (1 - p)
  }

  p0 <- clip(p_map(beta0))
  w0 <- w_map(beta0)

  # analytic dw/dB at beta0 (no clipping active assumed)
  # dp/dB = p(1-p) X
  # dw/dB = X * [ -(A)(1-p)/p + (1-A)p/(1-p) ]
  fac <- (-(A) * (1 - p0) / p0) + ((1 - A) * p0 / (1 - p0))
  dw <- X * fac  # rowwise multiply

  weightit_mock <- list(
    weights = w0,
    Mparts = list(dw_dBtreat = dw)
  )

  r1 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 1e-3, rho = 0, gamma = 0)
  r2 <- run_and_error(t, epsilon, strata, df, exposure,
                      weightit_mock, beta0, w_map,
                      fd_rel_step = 5e-4, rho = 0, gamma = 0)

  testthat::expect_true(is.finite(r1$err) && is.finite(r2$err))
  testthat::expect_true(r2$err > 0)

  err_ratio  <- r1$err / r2$err
  step_ratio <- r1$step2 / r2$step2

  testthat::expect_true(err_ratio / step_ratio > 0.4 && err_ratio / step_ratio < 2.5)
  testthat::expect_true(r2$err < r1$err)
})
