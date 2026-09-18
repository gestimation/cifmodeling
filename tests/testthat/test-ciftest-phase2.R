testthat::local_edition(3)

perturb_censoring_hazard <- function(object, hazard_row, amount) {
  out <- object
  selected <- object$hazard.table[hazard_row, , drop = FALSE]
  table_row <- which(
    out$table$stratum == selected$stratum &
      out$table$time == selected$time
  )
  stopifnot(length(table_row) == 1L)
  out$table$hazard[table_row] <- out$table$hazard[table_row] + amount

  stratum_rows <- which(out$table$stratum == selected$stratum)
  right <- cumprod(1 - out$table$hazard[stratum_rows])
  out$table$survival.right[stratum_rows] <- right
  out$table$survival.left[stratum_rows] <- c(1, head(right, -1L))
  if ("survival.right.raw" %in% names(out$table)) {
    out$table$survival.right.raw[stratum_rows] <- right
    out$table$survival.left.raw[stratum_rows] <- c(1, head(right, -1L))
    truncation <- out$prob.truncation
    out$table$survival.right.used[stratum_rows] <- if (is.null(truncation)) {
      right
    } else {
      pmax(right, truncation)
    }
    out$table$survival.left.used[stratum_rows] <- if (is.null(truncation)) {
      c(1, head(right, -1L))
    } else {
      pmax(c(1, head(right, -1L)), truncation)
    }
  }
  out$hazard.table <- out$table[out$table$n.censor > 0, , drop = FALSE]
  out$minimum.survival <- min(out$table$survival.right)
  out
}

closed_form_example <- function() {
  make_cell <- function(group, stratum, status) {
    data.frame(
      time = seq_along(status),
      status = status,
      group = factor(rep(group, length(status)), levels = c("A", "B")),
      stratum = stratum
    )
  }
  do.call(rbind, list(
    make_cell("A", "X", c(2L, 0L, 1L, 0L, 1L, 2L, 0L, 1L)),
    make_cell("A", "Y", c(0L, 2L, 1L, 0L, 2L, 1L, 0L, 1L)),
    make_cell("B", "X", c(2L, 0L, 0L, 1L, 2L, 1L, 0L, 1L)),
    make_cell("B", "Y", c(0L, 2L, 0L, 1L, 1L, 0L, 2L, 1L))
  ))
}

# Direct transcription of Scheike et al. (2023), Section 6.1 and p. 2487.
# It deliberately reconstructs Lambda1 from the working CIF and does not use
# build_closed_form_augmentation() or its stored H process.
paper_augmentation_reference <- function(
    base, working, time, status, x, exposure,
    strata.censor, strata.competing.risk, weights
) {
  event_time <- base$event.time
  working_cell <- cifmodeling:::ciftest_working_cell(
    exposure,
    strata.competing.risk
  )
  cell <- cifmodeling:::ciftest_augmentation_cell(
    exposure,
    strata.competing.risk,
    strata.censor
  )
  cells <- unique(data.frame(
    cell = cell,
    working.cell = working_cell,
    censor.stratum = as.character(strata.censor),
    stringsAsFactors = FALSE
  ))
  censor_label <- as.character(strata.censor)
  p <- ncol(x)
  h_process <- vector("list", nrow(cells))
  names(h_process) <- cells$cell

  for (k in seq_len(nrow(cells))) {
    level <- cells$cell[k]
    member <- which(cell == level & weights > 0)
    tab <- working$table[
      working$table$cell == cells$working.cell[k],
      ,
      drop = FALSE
    ]
    matched <- match(event_time, tab$time)
    present <- !is.na(matched)
    d_lambda1 <- numeric(length(event_time))
    d_lambda1[present] <-
      log1p(-tab$cif1.left[matched[present]]) -
      log1p(-tab$cif1.right[matched[present]])
    g_left <- cifmodeling:::predict_censoring_km(
      base$censoring,
      event_time,
      rep.int(cells$censor.stratum[k], length(event_time)),
      side = "left"
    )
    x_cell <- x[member[1L], ]
    integrand <- sweep(base$xbar, 2L, x_cell,
                       function(mean, value) value - mean)
    integrand <- integrand * (base$fh.weight * g_left * d_lambda1)
    h_process[[level]] <- apply(
      integrand,
      2L,
      function(column) rev(cumsum(rev(column)))
    ) |>
      matrix(nrow = length(event_time), ncol = p)
  }

  iid <- matrix(0, nrow(x), p)
  for (k in seq_len(nrow(base$censoring$hazard.table))) {
    hazard <- base$censoring$hazard.table[k, ]
    in_stratum <- censor_label == hazard$stratum
    at_risk <- in_stratum & time >= hazard$time
    censored <- in_stratum & time == hazard$time & status == 0L
    martingale <- weights * (
      as.numeric(censored) - hazard$hazard * as.numeric(at_risk)
    )
    hstar <- matrix(0, nrow(x), p)
    relevant <- which(in_stratum)
    working_left <- cifmodeling:::predict_working_aj(
      working,
      rep.int(hazard$time, length(relevant)),
      exposure[relevant],
      strata.competing.risk[relevant],
      side = "left"
    )
    g_left <- cifmodeling:::predict_censoring_km(
      base$censoring, hazard$time, hazard$stratum, side = "left"
    )
    h_index <- which(event_time >= hazard$time)[1L]
    if (!is.na(h_index)) {
      for (level in unique(cell[relevant])) {
        target <- relevant[cell[relevant] == level]
        local <- match(target, relevant)
        ratio <- working_left[local, "cif2"] /
          (working_left[local, "survival"] * g_left)
        ratio[working_left[local, "cif2"] == 0] <- 0
        hstar[target, ] <- outer(ratio, h_process[[level]][h_index, ])
      }
    }
    hbar <- colSums(hstar[at_risk, , drop = FALSE] * weights[at_risk]) /
      sum(weights[at_risk])
    iid <- iid + sweep(hstar, 2L, hbar, "-") * martingale
  }
  colnames(iid) <- colnames(x)
  list(score.iid.augment = iid, h.process = h_process)
}

testthat::test_that("censoring KM retains strata, fractional weights, and limits", {
  time <- c(1, 2, 3, 1, 2, 4)
  status <- c(2L, 0L, 1L, 0L, 2L, 1L)
  stratum <- factor(c("A", "A", "A", "B", "B", "B"))
  weight <- c(0.5, 1.5, 2, 1, 2, 1)

  fit <- cifmodeling:::estimate_censoring_km(
    time, status,
    strata = stratum,
    weights = weight
  )

  testthat::expect_s3_class(fit, "ciftest_censoring_km")
  testthat::expect_equal(
    cifmodeling:::predict_censoring_km(fit, c(2, 1), c("A", "B"), "left"),
    c(1, 1),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    cifmodeling:::predict_censoring_km(fit, c(2, 1), c("A", "B"), "right"),
    c(4 / 7, 3 / 4),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    cifmodeling:::predict_censoring_km(fit, c(3, 4), c("A", "B"), "left"),
    c(4 / 7, 3 / 4),
    tolerance = 1e-12
  )
})

testthat::test_that("administrative horizon completion is not a censoring jump", {
  time <- c(1, 2, 3, 4, 6, 6)
  status <- c(1L, 2L, 0L, 1L, 0L, 0L)
  censoring_event <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
  fit <- cifmodeling:::estimate_censoring_km(
    time,
    status,
    censoring.event = censoring_event
  )

  testthat::expect_equal(fit$hazard.table$time, 3)
  testthat::expect_false(any(fit$hazard.table$time == 6))
  testthat::expect_identical(fit$censoring.event, censoring_event)
})

testthat::test_that("AIPWCC includes event-free completion through tau", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- data.frame(
    time = c(1, 2, 3, 4, 8, 9, 10, 12),
    status = c(1L, 2L, 1L, 2L, 0L, 0L, 0L, 0L),
    group = factor(c("A", "A", "B", "B", "A", "B", "A", "B"))
  )
  methods <- data.frame(
    method_id = "seed_map",
    augmentation = TRUE,
    iteration = 0L,
    fixed.point.solver = "seed-map",
    strata.censor = NA_character_,
    strata.competing.risk = NA_character_,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  fit <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods,
    tau = 6
  )[[1L]]
  components <- fit$diagnostics$iteration.components

  testthat::expect_true(any(abs(components$horizon.completion) > 0))
  testthat::expect_equal(
    components$event + components$censor.past +
      components$censor.working.aj,
    fit$score.iid,
    tolerance = 1e-12
  )
  testthat::expect_false(any(
    fit$diagnostics$censoring.km$hazard.table$time == 6
  ))
})

testthat::test_that("working AJ retains A-by-strata cells and predictable limits", {
  single <- cifmodeling:::estimate_working_aj(
    t = 1:4,
    epsilon = c(1L, 2L, 0L, 1L),
    exposure = rep("A", 4),
    strata = rep("X", 4)
  )
  testthat::expect_equal(
    single$table$survival.right,
    c(3 / 4, 1 / 2, 1 / 2, 0),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    single$table$d.cif2,
    pmax(0, single$table$cif2.right - single$table$cif2.left),
    tolerance = 0
  )
  testthat::expect_equal(single$table$cif1.right,
                         c(1 / 4, 1 / 4, 1 / 4, 3 / 4),
                         tolerance = 1e-12)
  testthat::expect_equal(single$table$cif2.right,
                         c(0, 1 / 4, 1 / 4, 1 / 4),
                         tolerance = 1e-12)
  testthat::expect_equal(
    single$table$d.lambda1[c(1, 4)],
    c(-log(3 / 4), log(3)),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    single$table$d.lambda1,
    diff(c(0, -log1p(-single$table$cif1.right))),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    cifmodeling:::predict_working_aj(single, 2, "A", "X", "left"),
    matrix(c(3 / 4, 1 / 4, 0), nrow = 1,
           dimnames = list(NULL, c("survival", "cif1", "cif2"))),
    tolerance = 1e-12
  )

  df <- closed_form_example()
  stratified <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$group, df$stratum
  )
  testthat::expect_identical(nrow(stratified$cells), 4L)
  testthat::expect_setequal(stratified$cells$exposure, c("A", "B"))
  testthat::expect_setequal(stratified$cells$stratum, c("X", "Y"))
  testthat::expect_lt(
    max(abs(with(stratified$table,
                 survival.right + cif1.right + cif2.right - 1))),
    1e-12
  )
})

testthat::test_that("Fine-Gray working hazard rejects a terminal unit CIF", {
  testthat::expect_error(
    cifmodeling:::estimate_working_aj(
      t = 1:3,
      epsilon = c(0L, 0L, 1L),
      exposure = rep("A", 3)
    ),
    "subdistribution-hazard positivity"
  )
})

testthat::test_that("FH weights remain aligned for mixed-width decimal times", {
  df <- data.frame(
    time = c(0.25, 2.5, 12.75, 103.125, 0.5, 4.75, 15.5, 120.25),
    status = c(2L, 0L, 1L, 1L, 0L, 2L, 1L, 1L),
    group = factor(rep(c("A", "B"), each = 4))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  event_time <- sort(unique(df$time[df$status == 1L]))
  supplied <- seq_along(event_time) / 10

  fit <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    fh.weight = supplied
  )
  automatic <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    rho = 0.5,
    gamma = 0.5
  )

  testthat::expect_equal(fit$event.time, event_time)
  testthat::expect_equal(fit$fh.weight, supplied)
  testthat::expect_length(automatic$fh.weight, length(event_time))
  testthat::expect_true(all(is.finite(automatic$fh.weight)))
})

testthat::test_that("Fine-Gray score iid has an exact matrix decomposition", {
  df <- data.frame(
    time = 1:12,
    status = c(2L, 0L, 1L, 1L, 2L, 0L, 1L, 2L, 1L, 0L, 1L, 2L),
    group = factor(rep(c("A", "B"), 6))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  colnames(x) <- "groupB"

  fit <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    rho = 0.25,
    gamma = 0.5
  )

  testthat::expect_identical(dim(fit$score.iid), c(nrow(df), 1L))
  testthat::expect_equal(colSums(fit$score.iid.base), fit$score, tolerance = 1e-12)
  testthat::expect_equal(fit$score, fit$event.score, tolerance = 1e-12)
  testthat::expect_equal(
    colSums(fit$score.iid.censor),
    stats::setNames(0, colnames(fit$score.iid.censor)),
    tolerance = 1e-12
  )
  testthat::expect_true(any(abs(fit$score.iid.censor) > 1e-8))
  testthat::expect_lt(fit$diagnostics$score.decomposition.error, 1e-12)
  testthat::expect_lt(fit$diagnostics$censor.centering.error, 1e-12)
})

testthat::test_that("censoring-score derivative matches finite differences", {
  df <- data.frame(
    time = 1:12,
    status = c(2L, 0L, 1L, 1L, 2L, 0L, 1L, 2L, 1L, 0L, 1L, 2L),
    group = factor(rep(c("A", "B"), 6))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  base <- cifmodeling:::build_fg_score_iid(df$time, df$status, x)
  step <- 1e-6

  for (k in seq_len(nrow(base$censoring$hazard.table))) {
    plus <- cifmodeling:::build_fg_score_iid(
      df$time, df$status, x,
      censoring = perturb_censoring_hazard(base$censoring, k, step)
    )
    minus <- cifmodeling:::build_fg_score_iid(
      df$time, df$status, x,
      censoring = perturb_censoring_hazard(base$censoring, k, -step)
    )
    numerical <- (plus$event.score - minus$event.score) / (2 * step)
    testthat::expect_equal(
      unname(base$censor.derivative[k, ]),
      unname(numerical),
      tolerance = 2e-7,
      info = paste("censoring hazard row", k)
    )
  }
})

testthat::test_that("censoring contribution is zero without competing events", {
  df <- data.frame(
    time = 1:10,
    status = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    group = factor(rep(c("A", "B"), 5))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  fit <- cifmodeling:::build_fg_score_iid(df$time, df$status, x)

  testthat::expect_equal(fit$score.iid.censor, fit$score.iid.censor * 0)
  testthat::expect_equal(fit$censor.derivative, fit$censor.derivative * 0)
})

testthat::test_that("standard Gray exposes diagnostic iid without changing its variance", {
  df <- data.frame(
    time = c(1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9),
    status = c(1L, 0L, 2L, 1L, 2L, 0L, 1L, 1L, 0L, 2L, 1L, 0L),
    group = factor(c("A", "A", "B", "B", "A", "B",
                     "A", "B", "A", "B", "A", "B"))
  )
  fit <- ciftest(
    Event(time, status) ~ group,
    df,
    augmentation = FALSE,
    rho = 1
  )

  testthat::expect_true(fit$diagnostics$score.iid.available)
  testthat::expect_equal(colSums(fit$score.iid.base), fit$score, tolerance = 1e-12)
  testthat::expect_equal(
    colSums(fit$score.iid.censor),
    stats::setNames(0, colnames(fit$score.iid.censor)),
    tolerance = 1e-12
  )
  testthat::expect_equal(fit$score.iid.augment, fit$score.iid.augment * 0)
  testthat::expect_identical(fit$variance.method, "gray")
  testthat::expect_match(fit$diagnostics$score.iid.variance.role, "diagnostic only")
  process <- fit$diagnostics$fh.weight.process
  testthat::expect_s3_class(process, "data.frame")
  testthat::expect_identical(
    names(process), c("time", "weight", "stratum", "source")
  )
  testthat::expect_identical(
    unique(process$source), "gray-pooled-subdistribution-left"
  )
  testthat::expect_true(all(process$time %in% df$time[df$status == 1L]))
})

testthat::test_that("truncated censoring-score derivative matches finite differences", {
  df <- data.frame(
    time = 1:12,
    status = c(2L, 0L, 1L, 1L, 2L, 0L, 1L, 2L, 1L, 0L, 1L, 2L),
    group = factor(rep(c("A", "B"), 6))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  censoring <- cifmodeling:::estimate_censoring_km(
    df$time, df$status, prob.truncation = 0.8
  )
  base <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    censoring = censoring, prob.truncation = 0.8
  )
  step <- 1e-6

  for (k in seq_len(nrow(base$censoring$hazard.table))) {
    plus <- cifmodeling:::build_fg_score_iid(
      df$time, df$status, x,
      censoring = perturb_censoring_hazard(base$censoring, k, step),
      prob.truncation = 0.8
    )
    minus <- cifmodeling:::build_fg_score_iid(
      df$time, df$status, x,
      censoring = perturb_censoring_hazard(base$censoring, k, -step),
      prob.truncation = 0.8
    )
    numerical <- (plus$event.score - minus$event.score) / (2 * step)
    testthat::expect_equal(
      unname(base$censor.derivative[k, ]),
      unname(numerical),
      tolerance = 2e-7,
      info = paste("truncated censoring hazard row", k)
    )
  }
})

testthat::test_that("closed-form augmentation has a score-scale iid decomposition", {
  df <- closed_form_example()
  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    strata.censor = "stratum",
    strata.competing.risk = "stratum",
    augmentation = TRUE
  )

  testthat::expect_identical(fit$method,
                             "Closed-form augmented Fine-Gray score test")
  testthat::expect_identical(fit$variance.method, "score-iid")
  testthat::expect_identical(dim(fit$score.iid.augment), c(nrow(df), 1L))
  testthat::expect_true(any(abs(fit$score.iid.augment) > 1e-8))
  testthat::expect_equal(
    fit$score,
    colSums(fit$score.iid.base) + colSums(fit$score.iid.censor) +
      colSums(fit$score.iid.augment),
    tolerance = 1e-12
  )
  testthat::expect_equal(fit$score.iid,
                         fit$score.iid.base + fit$score.iid.censor +
                           fit$score.iid.augment,
                         tolerance = 1e-12)
  testthat::expect_equal(fit$vcov.score, crossprod(fit$score.iid),
                         tolerance = 1e-12)
  testthat::expect_lt(fit$diagnostics$augmentation.centering.error, 1e-12)
  testthat::expect_s3_class(fit$diagnostics$working.aj,
                            "ciftest_working_aj")
  process <- fit$diagnostics$fh.weight.process
  expected <- cifmodeling:::ciftest_null_fh_weight(
    df$time, df$status, 1L, 2L, rep.int(1, nrow(df)), 0, 0
  )
  testthat::expect_equal(process$time, sort(unique(df$time[df$status == 1L])))
  testthat::expect_equal(process$weight, expected)
  testthat::expect_identical(unique(process$source), "pooled-aj-left")

  # Independent raw martingale sum. Risk-set centering leaves its column sum
  # unchanged because every fitted censoring-martingale increment sums to zero.
  working <- fit$diagnostics$working.aj
  censoring <- fit$diagnostics$censoring.km
  event_time <- sort(unique(df$time[df$status == 1L]))
  cell <- cifmodeling:::ciftest_augmentation_cell(
    df$group,
    df$stratum,
    df$stratum
  )
  raw <- 0
  for (k in seq_len(nrow(censoring$hazard.table))) {
    row <- censoring$hazard.table[k, ]
    in_stratum <- df$stratum == row$stratum
    at_risk <- in_stratum & df$time >= row$time
    censored <- in_stratum & df$time == row$time & df$status == 0L
    martingale <- as.numeric(censored) - row$hazard * as.numeric(at_risk)
    prediction <- cifmodeling:::predict_working_aj(
      working,
      rep(row$time, nrow(df)),
      df$group,
      df$stratum,
      "left"
    )
    event_index <- which(event_time >= row$time)[1L]
    if (!is.na(event_index)) {
      h <- vapply(cell, function(level) {
        fit$diagnostics$h.process[[level]][event_index, 1L]
      }, numeric(1))
      g <- cifmodeling:::predict_censoring_km(
        censoring, row$time, row$stratum, "left"
      )
      raw <- raw + sum(
        martingale * h * prediction[, "cif2"] /
          (prediction[, "survival"] * g),
        na.rm = TRUE
      )
    }
  }
  testthat::expect_equal(unname(colSums(fit$score.iid.augment)), raw,
                         tolerance = 1e-12)
})

testthat::test_that("closed-form augmentation matches the paper equation", {
  df <- closed_form_example()
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  colnames(x) <- "groupB"
  weights <- rep.int(1, nrow(df))
  base <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    strata = df$stratum,
    weights = weights,
    rho = 0.5,
    gamma = 0.25
  )
  working <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$group, df$stratum, weights
  )
  actual <- cifmodeling:::build_closed_form_augmentation(
    base, working, df$time, df$status, x, df$group,
    df$stratum, df$stratum, weights
  )
  reference <- paper_augmentation_reference(
    base, working, df$time, df$status, x, df$group,
    df$stratum, df$stratum, weights
  )

  testthat::expect_equal(
    actual$score.iid.augment,
    reference$score.iid.augment,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    unname(actual$score),
    unname(colSums(reference$score.iid.augment)),
    tolerance = 1e-12
  )
  for (level in names(actual$h.process)) {
    testthat::expect_equal(
      unname(actual$h.process[[level]]),
      unname(reference$h.process[[level]]),
      tolerance = 1e-12,
      info = level
    )
  }
})

testthat::test_that("censoring and competing-risk nuisance strata are independent", {
  df <- closed_form_example()
  df$censor.stratum <- factor(rep(c("C1", "C2"), each = 4L, times = 4L))
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  colnames(x) <- "groupB"
  weights <- rep.int(1, nrow(df))

  base <- cifmodeling:::build_fg_score_iid(
    df$time,
    df$status,
    x,
    strata = df$censor.stratum,
    weights = weights,
    rho = 0.5,
    gamma = 0.25
  )
  working <- cifmodeling:::estimate_working_aj(
    df$time,
    df$status,
    df$group,
    df$stratum,
    weights
  )
  actual <- cifmodeling:::build_closed_form_augmentation(
    base,
    working,
    df$time,
    df$status,
    x,
    df$group,
    df$censor.stratum,
    df$stratum,
    weights
  )
  reference <- paper_augmentation_reference(
    base,
    working,
    df$time,
    df$status,
    x,
    df$group,
    df$censor.stratum,
    df$stratum,
    weights
  )

  testthat::expect_equal(
    actual$score.iid.augment,
    reference$score.iid.augment,
    tolerance = 1e-12
  )
  testthat::expect_setequal(
    actual$augmentation.cells$censor.stratum,
    levels(df$censor.stratum)
  )
  testthat::expect_setequal(
    actual$augmentation.cells$competing.risk.stratum,
    unique(df$stratum)
  )

  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    strata.censor = "censor.stratum",
    strata.competing.risk = "stratum",
    rho = 0.5,
    gamma = 0.25,
    augmentation = TRUE
  )
  testthat::expect_identical(fit$strata.censor.info$name, "censor.stratum")
  testthat::expect_identical(fit$strata.competing.risk.info$name, "stratum")
  testthat::expect_equal(
    fit$diagnostics$augmentation.cells,
    actual$augmentation.cells
  )
  testthat::expect_equal(unname(fit$score.iid.augment),
                         unname(actual$score.iid.augment),
                         tolerance = 1e-12)
})

testthat::test_that("both nuisance strata participate in complete-case selection", {
  df <- closed_form_example()
  df$censor.stratum <- factor(rep(c("C1", "C2"), each = 4L, times = 4L))
  df$competing.stratum <- factor(df$stratum)
  df$censor.stratum[1L] <- NA
  df$competing.stratum[10L] <- NA

  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    strata.censor = "censor.stratum",
    strata.competing.risk = "competing.stratum",
    augmentation = TRUE
  )

  testthat::expect_identical(fit$n, nrow(df) - 2L)
  testthat::expect_identical(
    fit$diagnostics$analysis.row.index,
    setdiff(seq_len(nrow(df)), c(1L, 10L))
  )
  testthat::expect_false(anyNA(fit$strata.censor.info$values))
  testthat::expect_false(anyNA(fit$strata.competing.risk.info$values))
})

testthat::test_that("closed-form augmentation vanishes without competing events", {
  df <- data.frame(
    time = rep(1:6, 2),
    status = rep(c(1L, 0L, 1L, 0L, 1L, 0L), 2),
    group = factor(rep(c("A", "B"), each = 6))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  base <- cifmodeling:::build_fg_score_iid(df$time, df$status, x)
  working <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$group
  )
  augment <- cifmodeling:::build_closed_form_augmentation(
    base, working, df$time, df$status, x, df$group
  )
  testthat::expect_equal(augment$score.iid.augment,
                         augment$score.iid.augment * 0)
  testthat::expect_equal(unname(augment$score), 0)
})

testthat::test_that("augmented ciftest records its nontrivial pooled-AJ FH process", {
  df <- closed_form_example()
  fit <- ciftest(
    Event(time, status) ~ group,
    data = df,
    augmentation = TRUE,
    rho = 0.5,
    gamma = 0.5
  )
  process <- fit$diagnostics$fh.weight.process
  expected <- cifmodeling:::ciftest_null_fh_weight(
    df$time, df$status, 1L, 2L, rep.int(1, nrow(df)), 0.5, 0.5
  )

  testthat::expect_equal(process$time, sort(unique(df$time[df$status == 1L])))
  testthat::expect_equal(process$weight, expected, tolerance = 1e-14)
  testthat::expect_identical(unique(process$stratum), "pooled")
  testthat::expect_identical(unique(process$source), "pooled-aj-left")
})

testthat::test_that("integer iteration refines the closed-form seed", {
  df <- closed_form_example()
  augmented <- ciftest(
    Event(time, status) ~ group,
    data = df,
    iteration = 0
  )
  iterated1 <- ciftest(
    Event(time, status) ~ group,
    data = df,
    iteration = 1
  )
  iterated2 <- ciftest(
    Event(time, status) ~ group,
    data = df,
    iteration = 2
  )

  testthat::expect_identical(
    augmented$method,
    "Closed-form augmented Fine-Gray score test"
  )
  testthat::expect_identical(
    iterated1$method,
    "Iterated time-weighted Fine-Gray score test"
  )
  testthat::expect_identical(iterated1$iteration, 1L)
  testthat::expect_identical(iterated2$iterations, 2L)
  testthat::expect_identical(dim(iterated2$score.iid), c(nrow(df), 1L))
  testthat::expect_equal(
    iterated2$vcov.score,
    crossprod(iterated2$score.iid),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    iterated2$score,
    colSums(iterated2$score.iid),
    tolerance = 1e-12
  )
  testthat::expect_lt(
    iterated2$fixed.point.residual,
    iterated1$fixed.point.residual
  )
  testthat::expect_equal(
    iterated2$score.iid.iterated,
    iterated2$score.iid,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    iterated2$score.iid.augment,
    augmented$score.iid.augment,
    tolerance = 1e-12
  )
  testthat::expect_true(all(is.na(augmented$score.iid.iterated)))
  printed <- paste(utils::capture.output(print(iterated2)), collapse = "\n")
  testthat::expect_match(
    printed, "Iteration: 2 of 2 requested refinements", fixed = TRUE
  )
  testthat::expect_match(printed, "Fixed-point residual", fixed = TRUE)
})

testthat::test_that("Appendix-E correction vanishes without censoring", {
  n <- 40L
  df <- data.frame(
    time = c(seq_len(n), 100),
    status = c(rep(c(1L, 2L, 2L, 1L), n / 4L), 0L),
    group = factor(c(rep(c("A", "A", "B", "B"), n / 4L), "A")),
    weight = c(rep.int(1, n), 0)
  )
  augmented <- ciftest(
    Event(time, status) ~ group,
    data = df,
    weights = "weight",
    iteration = 0
  )
  iterated <- ciftest(
    Event(time, status) ~ group,
    data = df,
    weights = "weight",
    iteration = 1
  )

  testthat::expect_equal(iterated$last.increment, 0, tolerance = 1e-12)
  testthat::expect_equal(iterated$fixed.point.residual, 0,
                         tolerance = 1e-12)
  testthat::expect_equal(iterated$score, augmented$score,
                         tolerance = 1e-12)
  testthat::expect_equal(iterated$vcov.score, augmented$vcov.score,
                         tolerance = 1e-12)
})

testthat::test_that("Appendix-E positivity is limited to its active support", {
  set.seed(1668883259L)
  df <- cifsimulate(
    n = 200,
    effect = "constant",
    shr = 0.3,
    rho1 = 0.2,
    rho2 = 3,
    rate = c(1, 1),
    beta1.L = -0.1,
    beta2 = c(A = 0, L = 0.3),
    tau = 6,
    analysis.tau = 6,
    effect.rho = 1,
    effect.gamma = 1,
    grid.width = 0.01,
    prob.A = 0.5,
    prob.L = 0.5,
    censor.rate = 0.719323755452953,
    censor.log.hazard = c(A = 0, L = 0)
  )
  testthat::expect_identical(
    as.integer(table(factor(df$status, levels = 0:2))),
    c(60L, 17L, 123L)
  )

  fit <- ciftest(Event(time, status) ~ A, data = df, iteration = 1)
  support <- fit$diagnostics$iteration.support
  testthat::expect_true(is.finite(fit$p.value))
  testthat::expect_gt(support$ignored.tail.censoring.times, 0L)
  testthat::expect_gt(support$skipped.zero.future.evaluations, 0L)
  testthat::expect_gt(support$minimum.active.working.survival, 0)

  working <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$A
  )
  late_time <- max(df$time[df$status == 1L]) + 0.1
  l_process <- cifmodeling:::ciftest_appendix_e_l_process(
    working,
    late_time,
    exposure = "1",
    active = TRUE
  )
  testthat::expect_equal(
    l_process$value,
    unname(l_process$prediction[, "cif2"] /
      l_process$prediction[, "survival"]),
    tolerance = 1e-12
  )
})

testthat::test_that("weighted iteration uses the same FH direction", {
  df <- closed_form_example()
  early <- ciftest(
    Event(time, status) ~ group,
    data = df,
    rho = 1,
    gamma = 0,
    iteration = 2
  )
  late <- ciftest(
    Event(time, status) ~ group,
    data = df,
    rho = 0,
    gamma = 1,
    iteration = 2
  )

  testthat::expect_true(is.finite(early$p.value))
  testthat::expect_true(is.finite(late$p.value))
  testthat::expect_false(isTRUE(all.equal(early$score, late$score)))
  testthat::expect_identical(nrow(early$diagnostics$iteration.history), 2L)
})

testthat::test_that("iterated score is invariant to analysis-row order", {
  df <- closed_form_example()
  original <- ciftest(
    Event(time, status) ~ group,
    data = df,
    rho = 1,
    iteration = 2
  )
  permutation <- c(32, 1, 24, 9, 17, 4, 28, 13, 20, 7, 31, 2,
                   25, 10, 18, 5, 29, 14, 21, 8, 30, 3, 26, 11,
                   19, 6, 27, 12, 16, 15, 23, 22)
  permuted <- ciftest(
    Event(time, status) ~ group,
    data = df[permutation, , drop = FALSE],
    rho = 1,
    iteration = 2
  )

  testthat::expect_equal(permuted$score, original$score, tolerance = 1e-12)
  testthat::expect_equal(permuted$vcov.score, original$vcov.score,
                         tolerance = 1e-12)
  testthat::expect_equal(
    permuted$score.iid[order(permutation), , drop = FALSE],
    original$score.iid,
    tolerance = 1e-12
  )
})

testthat::test_that("Fine-Gray score is invariant to analysis-row order", {
  df <- data.frame(
    time = 1:12,
    status = c(2L, 0L, 1L, 1L, 2L, 0L, 1L, 2L, 1L, 0L, 1L, 2L),
    group = factor(rep(c("A", "B"), 6))
  )
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  original <- cifmodeling:::build_fg_score_iid(df$time, df$status, x)
  order <- c(12, 1, 8, 3, 10, 5, 2, 11, 7, 4, 9, 6)
  permuted <- cifmodeling:::build_fg_score_iid(
    df$time[order], df$status[order], x[order, , drop = FALSE]
  )

  testthat::expect_equal(permuted$score, original$score, tolerance = 1e-12)
  testthat::expect_equal(
    permuted$score.iid[order(order), , drop = FALSE],
    original$score.iid,
    tolerance = 1e-12
  )
})

testthat::test_that("batch engine reuses nuisance fits and matches scalar ciftest", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  methods <- data.frame(
    method_id = c(
      "gray__unweighted", "gray__early",
      "pooled__unweighted", "pooled__early",
      "L__unweighted", "L__early"
    ),
    augmentation = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    strata.censor = c(NA, NA, NA, NA, "stratum", "stratum"),
    strata.competing.risk = c(NA, NA, NA, NA, "stratum", "stratum"),
    rho = c(0, 1, 0, 1, 0, 1),
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )
  testthat::expect_s3_class(batch, "ciftest_batch")
  testthat::expect_identical(names(batch), methods$method_id)
  timing <- attr(batch, "timing")
  batch_timing <- attr(batch, "batch.timing")
  testthat::expect_identical(timing$method_id, methods$method_id)
  testthat::expect_true(all(timing$elapsed_seconds >= 0))
  testthat::expect_equal(batch_timing$cache_misses, 5L)
  testthat::expect_equal(batch_timing$cache_hits, 5L)

  for (i in seq_len(nrow(methods))) {
    scalar <- ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = methods$augmentation[i],
      strata.censor = if (is.na(methods$strata.censor[i])) {
        NULL
      } else methods$strata.censor[i],
      strata.competing.risk = if (
        is.na(methods$strata.competing.risk[i])
      ) {
        NULL
      } else methods$strata.competing.risk[i],
      rho = methods$rho[i],
      gamma = methods$gamma[i]
    )
    actual <- batch[[methods$method_id[i]]]
    testthat::expect_s3_class(actual, "ciftest")
    testthat::expect_equal(actual$score, scalar$score, tolerance = 1e-12)
    testthat::expect_equal(actual$vcov.score, scalar$vcov.score,
                           tolerance = 1e-12)
    testthat::expect_equal(actual$score.iid, scalar$score.iid,
                           tolerance = 1e-12)
    testthat::expect_equal(actual$statistic, scalar$statistic,
                           tolerance = 1e-12)
    testthat::expect_equal(actual$p.value, scalar$p.value,
                           tolerance = 1e-12)
  }
})

testthat::test_that("batch engine carries finite iteration counts", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  methods <- data.frame(
    method_id = c("augmented", "iterated_1", "iterated_2"),
    augmentation = TRUE,
    iteration = 0:2,
    strata.censor = NA_character_,
    strata.competing.risk = NA_character_,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )

  for (i in seq_len(nrow(methods))) {
    scalar <- ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = TRUE,
      iteration = methods$iteration[i]
    )
    actual <- batch[[methods$method_id[i]]]
    testthat::expect_equal(actual$score.iid, scalar$score.iid,
                           tolerance = 1e-12)
    testthat::expect_equal(actual$p.value, scalar$p.value,
                           tolerance = 1e-12)
    testthat::expect_identical(actual$iterations,
                               as.integer(methods$iteration[i]))
  }
})

testthat::test_that("batch Fine-Gray score-IID method uses its own variance", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  methods <- data.frame(
    method_id = c("fine_gray", "closed_form", "iterated_1"),
    augmentation = c(FALSE, TRUE, TRUE),
    iteration = c(0L, 0L, 1L),
    score.construction = c("fine-gray", "standard", "standard"),
    strata.censor = NA_character_,
    strata.competing.risk = NA_character_,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )
  fit <- batch$fine_gray
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  colnames(x) <- "a_B"
  reference <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x
  )

  testthat::expect_s3_class(fit, "ciftest")
  testthat::expect_identical(fit$method, "Fine-Gray score test")
  testthat::expect_identical(fit$variance.method, "score-iid")
  testthat::expect_false(fit$augmentation)
  testthat::expect_equal(fit$score, reference$score, tolerance = 1e-12)
  testthat::expect_equal(fit$score.iid, reference$score.iid,
                         tolerance = 1e-12)
  testthat::expect_equal(fit$vcov.score, crossprod(reference$score.iid),
                         tolerance = 1e-12)
  testthat::expect_false(isTRUE(all.equal(
    fit$vcov.score, batch$closed_form$vcov.score
  )))
})

testthat::test_that("direct fixed-point batch agrees with converged iteration", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  methods <- data.frame(
    method_id = c("finite_10", "direct_0", "direct_early"),
    augmentation = TRUE,
    iteration = c(10L, 0L, 0L),
    fixed.point.solver = c("finite", "direct", "direct"),
    strata.censor = NA_character_,
    strata.competing.risk = NA_character_,
    rho = c(0, 0, 5),
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )
  finite <- batch$finite_10
  direct <- batch$direct_0
  early <- batch$direct_early

  testthat::expect_s3_class(direct, "ciftest")
  testthat::expect_identical(
    direct$method,
    "Direct fixed-point time-weighted Fine-Gray score test"
  )
  testthat::expect_identical(direct$iterations, NA_integer_)
  testthat::expect_true(direct$converged)
  testthat::expect_lt(direct$fixed.point.residual, 1e-10)
  testthat::expect_equal(direct$score.iid, finite$score.iid,
                         tolerance = 1e-7)
  testthat::expect_equal(direct$p.value, finite$p.value,
                         tolerance = 1e-8)
  testthat::expect_true(is.finite(
    direct$diagnostics$iteration.support$operator.spectral.radius
  ))
  testthat::expect_true(
    direct$diagnostics$iteration.support$operator.contractive
  )
  testthat::expect_true(is.finite(
    direct$diagnostics$iteration.support$system.reciprocal.condition
  ))
  testthat::expect_true(is.finite(early$p.value))
  testthat::expect_identical(
    early$diagnostics$iteration.support$state.dimension,
    direct$diagnostics$iteration.support$state.dimension
  )
})

testthat::test_that("batch exposes the Appendix-E seed-map diagnostic", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  n <- 40L
  df <- data.frame(
    time = c(seq_len(n), 100),
    status = c(rep(c(1L, 2L, 2L, 1L), n / 4L), 0L),
    group = factor(c(rep(c("A", "A", "B", "B"), n / 4L), "A")),
    weight = c(rep.int(1, n), 0)
  )
  methods <- data.frame(
    method_id = c("closed", "seed_map"),
    augmentation = TRUE,
    iteration = 0L,
    fixed.point.solver = c("finite", "seed-map"),
    strata.censor = NA_character_,
    strata.competing.risk = NA_character_,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    weights = "weight",
    methods = methods
  )
  seed_map <- batch$seed_map

  testthat::expect_s3_class(seed_map, "ciftest")
  testthat::expect_identical(
    seed_map$method,
    "AIPWCC seed-map diagnostic score test"
  )
  testthat::expect_identical(seed_map$iterations, 0L)
  testthat::expect_identical(
    seed_map$diagnostics$iteration.support$solver,
    "seed-map"
  )
  testthat::expect_equal(
    seed_map$diagnostics$fixed.point$value,
    seed_map$diagnostics$fixed.point$seed,
    tolerance = 0
  )
  testthat::expect_equal(seed_map$score, batch$closed$score,
                         tolerance = 1e-12)
  testthat::expect_equal(seed_map$vcov.score, batch$closed$vcov.score,
                         tolerance = 1e-12)
  components <- seed_map$diagnostics$iteration.components
  testthat::expect_named(
    components,
    c(
      "event", "censor.past", "censor.working.aj", "censor.total",
      "horizon.completion"
    )
  )
  testthat::expect_equal(
    components$event + components$censor.past +
      components$censor.working.aj,
    seed_map$score.iid,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    components$censor.past + components$censor.working.aj,
    components$censor.total,
    tolerance = 0
  )
})

testthat::test_that("finite and direct iterates are anchored at closed form", {
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  methods <- data.frame(
    method_id = c("closed", "seed_map", "iterated_1", "direct"),
    augmentation = TRUE,
    iteration = c(0L, 0L, 1L, 0L),
    fixed.point.solver = c("finite", "seed-map", "finite", "direct"),
    strata.censor = NA_character_,
    strata.competing.risk = NA_character_,
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )

  testthat::expect_false(
    batch$seed_map$diagnostics$iteration.anchor$applied
  )
  for (method in c("iterated_1", "direct")) {
    fit <- batch[[method]]
    anchor <- fit$diagnostics$iteration.anchor
    testthat::expect_true(anchor$applied)
    testthat::expect_identical(
      anchor$construction,
      "closed-form + AIPWCC(value) - AIPWCC(seed)"
    )
    testthat::expect_equal(
      anchor$seed.aipwcc.score,
      batch$seed_map$score,
      tolerance = 1e-12
    )
    testthat::expect_equal(
      anchor$score.adjustment,
      batch$closed$score - batch$seed_map$score,
      tolerance = 1e-12
    )
    testthat::expect_equal(
      fit$score,
      anchor$raw.aipwcc.score + anchor$score.adjustment,
      tolerance = 1e-12
    )
    testthat::expect_equal(anchor$identity.error, 0, tolerance = 0)
  }
})

testthat::test_that("AIPWCC anchoring preserves the zeroth step", {
  closed <- matrix(c(1, 2, 3, 4), ncol = 2L)
  seed <- matrix(c(4, 3, 2, 1), ncol = 2L)
  value <- seed + matrix(c(0.1, -0.2, 0.3, -0.4), ncol = 2L)
  anchored <- cifmodeling:::ciftest_anchor_aipwcc_iid(
    closed, value, seed
  )
  zeroth <- cifmodeling:::ciftest_anchor_aipwcc_iid(
    closed, seed, seed
  )

  testthat::expect_equal(zeroth$score.iid, closed, tolerance = 0)
  testthat::expect_equal(
    anchored$score.iid,
    closed + value - seed,
    tolerance = 1e-15
  )
  testthat::expect_equal(anchored$identity.error, 0, tolerance = 0)
  testthat::expect_error(
    cifmodeling:::ciftest_anchor_aipwcc_iid(
      closed, value[-1L, , drop = FALSE], seed
    ),
    "must conform"
  )
})

testthat::test_that("direct convergence requires a contraction", {
  testthat::expect_true(cifmodeling:::ciftest_direct_converged(
    residual = 1e-12, tolerance = 1e-8, spectral.radius = 0.9
  ))
  testthat::expect_false(cifmodeling:::ciftest_direct_converged(
    residual = 1e-12, tolerance = 1e-8, spectral.radius = 1
  ))
  testthat::expect_false(cifmodeling:::ciftest_direct_converged(
    residual = 1e-12, tolerance = 1e-8, spectral.radius = 1.2
  ))
  testthat::expect_false(cifmodeling:::ciftest_direct_converged(
    residual = 1e-4, tolerance = 1e-8, spectral.radius = 0.5
  ))
})

testthat::test_that("Appendix-E operator is linear", {
  df <- closed_form_example()
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  base <- cifmodeling:::build_fg_score_iid(df$time, df$status, x)
  working <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$group
  )
  setup <- cifmodeling:::ciftest_iteration_setup(
    base, working, df$time, df$status, x, df$group
  )
  operator <- cifmodeling:::ciftest_appendix_e_operator(
    setup, base, working
  )
  state_dimension <- operator$state.dimension
  v <- array(seq_len(state_dimension) / state_dimension,
             dim = c(length(setup$working.cell), setup$event.count, 1L))
  zero <- array(0, dim = dim(v))
  mapped <- cifmodeling:::ciftest_appendix_e_map(
    v, zero, base, working,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    setup$risk.weight
  )
  testthat::expect_equal(
    as.numeric(mapped),
    as.numeric(operator$K %*% as.numeric(v)),
    tolerance = 1e-12
  )
})

testthat::test_that("batch engine applies common nuisance complete cases", {
  df <- closed_form_example()
  df$stratum[1L] <- NA
  methods <- data.frame(
    method_id = "L",
    augmentation = TRUE,
    strata.censor = "stratum",
    strata.competing.risk = "stratum",
    rho = 0,
    gamma = 0,
    stringsAsFactors = FALSE
  )
  fit <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group,
    data = df,
    methods = methods
  )
  testthat::expect_equal(nobs(fit$L), nrow(df) - 1L)
  testthat::expect_false(fit$L$diagnostics$analysis.included[1L])
})

testthat::test_that("Rcpp score and augmentation kernels match R reference", {
  namespace <- asNamespace("cifmodeling")
  testthat::skip_if_not(exists(
    "_cifmodeling_ciftest_fg_iid_kernel_cpp",
    envir = namespace,
    inherits = FALSE
  ))
  testthat::skip_if_not(exists(
    "_cifmodeling_ciftest_augmentation_iid_kernel_cpp",
    envir = namespace,
    inherits = FALSE
  ))
  old <- options(cifmodeling.ciftest.engine = "R")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  x <- stats::model.matrix(~ group, df)[, -1L, drop = FALSE]
  colnames(x) <- "groupB"
  base_r <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    strata = df$stratum,
    rho = 0.5,
    gamma = 0.25
  )
  working <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$group, df$stratum
  )
  augment_r <- cifmodeling:::build_closed_form_augmentation(
    base_r, working, df$time, df$status, x, df$group,
    df$stratum, df$stratum
  )

  options(cifmodeling.ciftest.engine = "Rcpp")
  base_cpp <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    strata = df$stratum,
    rho = 0.5,
    gamma = 0.25
  )
  augment_cpp <- cifmodeling:::build_closed_form_augmentation(
    base_cpp, working, df$time, df$status, x, df$group,
    df$stratum, df$stratum
  )
  fields <- c(
    "score", "event.score", "score.iid", "score.iid.base",
    "score.iid.censor", "xbar", "risk.total", "event.count",
    "censor.derivative"
  )
  for (field in fields) {
    testthat::expect_equal(
      base_cpp[[field]], base_r[[field]], tolerance = 1e-10,
      info = field
    )
  }
  testthat::expect_equal(
    augment_cpp$score.iid.augment,
    augment_r$score.iid.augment,
    tolerance = 1e-10
  )
  testthat::expect_equal(augment_cpp$score, augment_r$score,
                         tolerance = 1e-10)
  testthat::expect_identical(base_cpp$diagnostics$engine, "Rcpp-prefix")
  testthat::expect_identical(augment_cpp$diagnostics$engine, "Rcpp-prefix")
})

testthat::test_that("multiweight Rcpp batch matches scalar Rcpp fits", {
  namespace <- asNamespace("cifmodeling")
  symbols <- c(
    "_cifmodeling_ciftest_fg_iid_multi_kernel_cpp",
    "_cifmodeling_ciftest_augmentation_iid_multi_kernel_cpp"
  )
  for (symbol in symbols) {
    testthat::skip_if_not(exists(
      symbol, envir = namespace, inherits = FALSE
    ))
  }
  old <- options(cifmodeling.ciftest.engine = "Rcpp")
  on.exit(options(old), add = TRUE)
  df <- closed_form_example()
  methods <- data.frame(
    method_id = c("pooled_0", "pooled_1", "L_0", "L_1"),
    augmentation = TRUE,
    strata.censor = c(NA, NA, "stratum", "stratum"),
    strata.competing.risk = c(NA, NA, "stratum", "stratum"),
    rho = c(0, 1, 0, 1),
    gamma = 0,
    stringsAsFactors = FALSE
  )
  batch <- cifmodeling:::ciftest_batch_internal(
    Event(time, status) ~ group, df, methods
  )
  for (i in seq_len(nrow(methods))) {
    scalar <- ciftest(
      Event(time, status) ~ group,
      data = df,
      augmentation = TRUE,
      strata.censor = if (is.na(methods$strata.censor[i])) {
        NULL
      } else methods$strata.censor[i],
      strata.competing.risk = if (
        is.na(methods$strata.competing.risk[i])
      ) {
        NULL
      } else methods$strata.competing.risk[i],
      rho = methods$rho[i],
      gamma = methods$gamma[i]
    )
    actual <- batch[[methods$method_id[i]]]
    testthat::expect_equal(actual$score.iid, scalar$score.iid,
                           tolerance = 1e-10)
    testthat::expect_equal(actual$vcov.score, scalar$vcov.score,
                           tolerance = 1e-8)
    testthat::expect_equal(actual$p.value, scalar$p.value,
                           tolerance = 1e-8)
    testthat::expect_identical(
      actual$diagnostics$score.engine, "Rcpp-prefix-multi"
    )
    testthat::expect_identical(
      actual$diagnostics$augmentation.engine, "Rcpp-prefix-multi"
    )
  }
  timing <- attr(batch, "batch.timing")
  testthat::expect_true(timing$elapsed_multiweight_precompute_seconds >= 0)
  testthat::expect_length(timing$multiweight_errors, 0L)
})

testthat::test_that("Fine-Gray prefix kernels match legacy with ties and weights", {
  namespace <- asNamespace("cifmodeling")
  symbols <- c(
    "_cifmodeling_ciftest_fg_iid_prefix_kernel_cpp",
    "_cifmodeling_ciftest_fg_iid_prefix_multi_kernel_cpp"
  )
  for (symbol in symbols) {
    testthat::skip_if_not(exists(
      symbol, envir = namespace, inherits = FALSE
    ))
  }
  old <- options(
    cifmodeling.ciftest.engine = "Rcpp",
    cifmodeling.ciftest.fg.engine = "check"
  )
  on.exit(options(old), add = TRUE)

  time <- rep(seq_len(8L), each = 3L)
  status <- rep(c(11L, 22L, 0L), 8L)
  status[24L] <- 11L
  strata <- rep(c("S1", "S2"), length.out = length(time))
  weights <- rep(c(1, 0.5, 1.25, 0.75), length.out = length(time))
  weights[c(5L, 18L)] <- 0
  x <- cbind(
    group = rep(c(0, 1), length.out = length(time)),
    z = (seq_along(time) - mean(seq_along(time))) / length(time)
  )

  scalar <- cifmodeling:::build_fg_score_iid(
    time, status, x,
    code.event1 = 11L,
    code.event2 = 22L,
    code.censoring = 0L,
    strata = strata,
    weights = weights,
    rho = 0.5,
    gamma = 0.25
  )
  testthat::expect_identical(scalar$diagnostics$engine, "Rcpp-check")
  testthat::expect_lt(
    max(abs(colSums(scalar$score.iid.censor))), 1e-10
  )

  multi <- cifmodeling:::build_fg_score_iid_multi(
    time, status, x,
    rho = c(0, 0.5, 0),
    gamma = c(0, 0.25, 1),
    code.event1 = 11L,
    code.event2 = 22L,
    code.censoring = 0L,
    strata = strata,
    weights = weights
  )
  testthat::expect_length(multi, 3L)
  testthat::expect_true(all(vapply(
    multi,
    function(value) identical(value$diagnostics$engine,
                               "Rcpp-check-multi"),
    logical(1L)
  )))
})

testthat::test_that("augmentation prefix kernels match legacy by nuisance cell", {
  namespace <- asNamespace("cifmodeling")
  symbols <- c(
    "_cifmodeling_ciftest_augmentation_iid_prefix_kernel_cpp",
    "_cifmodeling_ciftest_augmentation_iid_prefix_multi_kernel_cpp"
  )
  for (symbol in symbols) {
    testthat::skip_if_not(exists(
      symbol, envir = namespace, inherits = FALSE
    ))
  }
  old <- options(
    cifmodeling.ciftest.engine = "Rcpp",
    cifmodeling.ciftest.fg.engine = "prefix",
    cifmodeling.ciftest.augmentation.engine = "check"
  )
  on.exit(options(old), add = TRUE)

  df <- closed_form_example()
  group_b <- as.numeric(df$group == levels(df$group)[2L])
  x <- cbind(groupB = group_b, groupB2 = 2 * group_b)
  weights <- rep(c(1, 0.5, 1.25, 0.75), length.out = nrow(df))
  weights[2L] <- 0
  base <- cifmodeling:::build_fg_score_iid(
    df$time, df$status, x,
    strata = df$stratum,
    weights = weights,
    rho = 0.5,
    gamma = 0.25
  )
  working <- cifmodeling:::estimate_working_aj(
    df$time, df$status, df$group, df$stratum,
    weights = weights
  )
  scalar <- cifmodeling:::build_closed_form_augmentation(
    base, working, df$time, df$status, x, df$group,
    df$stratum, df$stratum, weights = weights
  )
  testthat::expect_identical(scalar$diagnostics$engine, "Rcpp-check")

  bases <- cifmodeling:::build_fg_score_iid_multi(
    df$time, df$status, x,
    rho = c(0, 1), gamma = c(0, 0),
    strata = df$stratum, weights = weights
  )
  multi <- cifmodeling:::build_closed_form_augmentation_multi(
    bases, working, df$time, df$status, x, df$group,
    df$stratum, df$stratum, weights = weights
  )
  testthat::expect_length(multi, 2L)
  testthat::expect_true(all(vapply(
    multi,
    function(value) identical(value$diagnostics$engine,
                               "Rcpp-check-multi"),
    logical(1L)
  )))
})
