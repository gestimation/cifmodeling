event_strata_example <- function() {
  n_per_center <- 60L
  within <- rep(seq_len(n_per_center), 2L)
  center <- factor(rep(c("C1", "C2"), each = n_per_center))
  data.frame(
    time = within + ifelse(center == "C2", 0.05, 0),
    status = rep(
      c(1L, 1L, 2L, 2L, 0L, 0L, 1L, 2L, 0L, 1L, 2L, 0L),
      length.out = 2L * n_per_center
    ),
    group = factor(rep(c("A", "B"), length.out = 2L * n_per_center)),
    center = center
  )
}

testthat::test_that("strata.event stratifies standard Gray and log-rank tests", {
  df <- event_strata_example()
  gray <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = FALSE,
    strata.event = "center"
  )
  gray_reference <- cifmodeling:::calculate_gray(
    t = df$time,
    epsilon = df$status,
    exposure = "group",
    strata = df$center,
    data = df
  )
  testthat::expect_equal(gray$score, gray_reference$score, tolerance = 1e-12)
  testthat::expect_equal(gray$vcov.score, gray_reference$var,
                         tolerance = 1e-12)
  testthat::expect_identical(gray$strata.event.info$name, "center")
  testthat::expect_false(gray$diagnostics$score.iid.available)
  testthat::expect_match(
    gray$diagnostics$score.iid.error,
    "stratified standard-Gray"
  )
  testthat::expect_setequal(
    unique(gray$diagnostics$fh.weight.process$stratum),
    levels(df$center)
  )
  testthat::expect_identical(
    unique(gray$diagnostics$fh.weight.process$source),
    "gray-event-stratified-subdistribution-left"
  )

  survival_df <- df
  survival_df$status <- as.integer(survival_df$status == 1L)
  logrank <- ciftest(
    Event(time, status) ~ group,
    data = survival_df,
    outcome.type = "survival",
    augmentation = FALSE,
    strata.event = "center"
  )
  logrank_reference <- cifmodeling:::calculate_log_rank(
    t = survival_df$time,
    epsilon = survival_df$status,
    exposure = "group",
    weights = rep.int(1, nrow(survival_df)),
    strata = survival_df$center,
    data = survival_df
  )
  testthat::expect_equal(logrank$score, logrank_reference$score,
                         tolerance = 1e-12)
  testthat::expect_equal(logrank$vcov.score, logrank_reference$var,
                         tolerance = 1e-12)
})

testthat::test_that("event-stratified closed form equals stacked stratum fits", {
  df <- event_strata_example()
  combined <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = "center",
    strata.censor = "center",
    strata.competing.risk = "center"
  )
  separate <- lapply(levels(df$center), function(level) {
    ciftest(
      Event(time, status) ~ group,
      data = df[df$center == level, , drop = FALSE],
      augmentation = TRUE
    )
  })
  stacked_iid <- do.call(rbind, lapply(separate, `[[`, "score.iid"))

  testthat::expect_equal(
    combined$score,
    Reduce(`+`, lapply(separate, `[[`, "score")),
    tolerance = 1e-10
  )
  testthat::expect_equal(combined$score.iid, stacked_iid, tolerance = 1e-10)
  testthat::expect_equal(
    combined$vcov.score,
    crossprod(stacked_iid),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    colSums(combined$score.iid), combined$score,
    tolerance = 1e-10
  )
  testthat::expect_identical(
    combined$diagnostics$score.engine, "R-event-stratified"
  )
  testthat::expect_setequal(
    combined$diagnostics$augmentation.cells$event.stratum,
    levels(df$center)
  )
})

testthat::test_that("event strata do not implicitly stratify nuisance models", {
  df <- event_strata_example()
  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = "center"
  )

  testthat::expect_identical(fit$strata.censor.info$name, NULL)
  testthat::expect_identical(fit$strata.competing.risk.info$name, NULL)
  testthat::expect_identical(
    fit$diagnostics$censoring.km$strata.levels, "pooled"
  )
  testthat::expect_identical(
    unique(fit$diagnostics$working.aj$cells$stratum), "pooled"
  )
  testthat::expect_identical(
    unique(fit$diagnostics$augmentation.cells$censor.stratum), "pooled"
  )
  testthat::expect_identical(
    unique(fit$diagnostics$augmentation.cells$competing.risk.stratum),
    "pooled"
  )
  testthat::expect_setequal(
    fit$diagnostics$augmentation.cells$event.stratum,
    levels(df$center)
  )
})

testthat::test_that("batch event strata agree with scalar fits", {
  df <- event_strata_example()
  methods <- data.frame(
    method_id = c("gray", "closed"),
    augmentation = c(FALSE, TRUE),
    iteration = 0L,
    strata.event = "center",
    strata.censor = c(NA_character_, "center"),
    strata.competing.risk = c(NA_character_, "center"),
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )
  scalar_gray <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = FALSE,
    strata.event = "center"
  )
  scalar_closed <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = "center",
    strata.censor = "center",
    strata.competing.risk = "center"
  )

  testthat::expect_equal(batch$gray$p.value, scalar_gray$p.value,
                         tolerance = 1e-12)
  testthat::expect_equal(batch$closed$score.iid, scalar_closed$score.iid,
                         tolerance = 1e-10)
  testthat::expect_equal(batch$closed$p.value, scalar_closed$p.value,
                         tolerance = 1e-12)
})

testthat::test_that("event strata use common missing-data handling and guard iteration", {
  df <- event_strata_example()
  df$center[3L] <- NA
  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = "center"
  )
  augmented <- augment(fit)
  testthat::expect_equal(nobs(fit), nrow(df) - 1L)
  testthat::expect_true(all(is.na(augmented$.score_iid[[3L]])))

  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group,
      data = event_strata_example(),
      augmentation = TRUE,
      iteration = 1L,
      strata.event = "center"
    ),
    "Analysis and nuisance-model strata"
  )
})

testthat::test_that("multiple strata columns equal explicit interaction columns", {
  df <- event_strata_example()
  df$region <- factor(rep(rep(c("R1", "R2"), each = 6L),
                          length.out = nrow(df)))
  df$event.cell <- interaction(df[c("center", "region")], drop = TRUE)
  df$censor.cell <- interaction(df[c("center", "region")], drop = TRUE)
  df$working.cell <- interaction(df[c("center", "region")], drop = TRUE)

  multiple <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = c("center", "region"),
    strata.censor = c("center", "region"),
    strata.competing.risk = c("center", "region"),
    tau = 46
  )
  explicit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = "event.cell",
    strata.censor = "censor.cell",
    strata.competing.risk = "working.cell",
    tau = 46
  )

  testthat::expect_equal(multiple$score, explicit$score, tolerance = 1e-10)
  testthat::expect_equal(
    multiple$score.iid, explicit$score.iid, tolerance = 1e-10
  )
  testthat::expect_equal(
    multiple$vcov.score, explicit$vcov.score, tolerance = 1e-10
  )
  testthat::expect_identical(
    multiple$strata.event.info$columns, c("center", "region")
  )
  testthat::expect_identical(
    multiple$strata.event.info$name, "center:region"
  )
  testthat::expect_identical(multiple$strata.event.info$role, "event")
  testthat::expect_equal(
    sum(multiple$strata.event.info$counts$n), nrow(df)
  )
  testthat::expect_equal(
    sum(multiple$strata.event.info$counts$weight), nrow(df)
  )
  testthat::expect_setequal(
    names(multiple$strata.event.info$mapping),
    c("stratum", "center", "region")
  )
})

testthat::test_that("strata roles enforce exposure and branch semantics", {
  df <- event_strata_example()

  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = FALSE,
      strata.event = "group"
    ),
    "must not include the grouping variable"
  )
  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = TRUE,
      strata.competing.risk = c("center", "group")
    ),
    "already fitted within exposure-by-stratum"
  )
  censor_by_group <- cifmodeling:::ciftest_prepare(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.censor = "group"
  )
  testthat::expect_identical(
    censor_by_group$strata.censor.info$columns, "group"
  )

  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = FALSE,
      strata.censor = c("center")
    ),
    "only for competing-risk score and augmented"
  )
  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = FALSE,
      strata.competing.risk = c("center")
    ),
    "requires `test = \"augmented\"`"
  )
})

testthat::test_that("multiple strata columns share missing-data handling", {
  df <- event_strata_example()
  df$region <- factor(rep(rep(c("R1", "R2"), each = 6L),
                          length.out = nrow(df)))
  df$region[3L] <- NA
  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = FALSE,
    strata.event = c("center", "region")
  )
  augmented <- augment(fit)
  testthat::expect_equal(nobs(fit), nrow(df) - 1L)
  testthat::expect_false(fit$diagnostics$analysis.included[3L])
  testthat::expect_true(all(is.na(augmented$.score_iid[[3L]])))
})

testthat::test_that("invalid multiple strata specifications fail early", {
  df <- event_strata_example()
  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group, df,
      strata.event = character()
    ),
    "non-empty character vector"
  )
  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group, df,
      strata.event = c("center", "center")
    ),
    "duplicated column names"
  )
  testthat::expect_error(
    ciftest(
      Event(time, status) ~ group, df,
      strata.event = c("center", "unknown")
    ),
    "not found in data"
  )
})

testthat::test_that("batch list-column strata agree with scalar fits", {
  df <- event_strata_example()
  df$region <- factor(rep(rep(c("R1", "R2"), each = 6L),
                          length.out = nrow(df)))
  methods <- data.frame(
    method_id = "closed-multiple",
    augmentation = TRUE,
    iteration = 0L,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  methods$strata.event <- I(list(c("center", "region")))
  methods$strata.censor <- I(list(c("center", "region")))
  methods$strata.competing.risk <- I(list(c("center", "region")))

  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods,
    tau = 46
  )
  scalar <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    strata.event = c("center", "region"),
    strata.censor = c("center", "region"),
    strata.competing.risk = c("center", "region"),
    tau = 46
  )
  testthat::expect_equal(
    batch[["closed-multiple"]]$score.iid,
    scalar$score.iid,
    tolerance = 1e-10
  )
  testthat::expect_identical(
    batch[["closed-multiple"]]$strata.censor.info$columns,
    c("center", "region")
  )
})

testthat::test_that("batch list-column strata use the scalar validator", {
  df <- event_strata_example()
  methods <- data.frame(
    method_id = "invalid-empty",
    augmentation = TRUE,
    iteration = 0L,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  methods$strata.event <- I(list(character()))
  methods$strata.censor <- I(list(NULL))
  methods$strata.competing.risk <- I(list(NULL))

  testthat::expect_error(
    cifmodeling:::ciftest_batch_internal(
      Event(time, status) ~ group,
      data = df,
      methods = methods
    ),
    "non-empty character vector"
  )
})
