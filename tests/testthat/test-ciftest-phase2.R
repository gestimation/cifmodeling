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
