#' Estimate a stratified censoring Kaplan-Meier distribution
#'
#' Uses the package Kaplan-Meier counts, retaining both left and right limits.
#' The default keeps the existing censoring-first convention for ciftest.
#' With tie.order="event_first", actual failures leave before the censoring
#' jump, consistently with observation of failures when T <= C.
#'
#' @keywords internal
estimate_censoring_km <- function(
    t,
    epsilon,
    code.censoring = 0L,
    strata = rep.int("pooled", length(t)),
    weights = rep.int(1, length(t)),
    prob.bound = 1e-7,
    prob.truncation = NULL,
    censoring.event = epsilon == code.censoring,
    tie.order = c("censoring_first", "event_first")
) {
  tie.order <- match.arg(tie.order)
  n <- length(t)
  if (length(epsilon) != n || length(strata) != n || length(weights) != n ||
      length(censoring.event) != n) {
    stop("Censoring KM inputs must have the same length.", call. = FALSE)
  }
  if (!n || any(!is.finite(t)) || any(t < 0)) {
    stop("Censoring KM requires finite non-negative follow-up times.", call. = FALSE)
  }
  if (anyNA(epsilon) || anyNA(strata)) {
    stop("Censoring KM inputs must not contain missing status or strata.", call. = FALSE)
  }
  if (!is.logical(censoring.event) || anyNA(censoring.event)) {
    stop("`censoring.event` must be a non-missing logical vector.",
         call. = FALSE)
  }
  if (any(censoring.event & epsilon != code.censoring)) {
    stop("A censoring event must have the censoring status code.",
         call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Censoring KM weights must be finite and non-negative.", call. = FALSE)
  }
  if (!is.numeric(prob.bound) || length(prob.bound) != 1L ||
      !is.finite(prob.bound) || prob.bound <= 0 || prob.bound >= 1) {
    stop("`prob.bound` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.null(prob.truncation) &&
      (!is.numeric(prob.truncation) || length(prob.truncation) != 1L ||
       !is.finite(prob.truncation) || prob.truncation <= prob.bound ||
       prob.truncation >= 1)) {
    stop("`prob.truncation` must be NULL or strictly between `prob.bound` and one.",
         call. = FALSE)
  }

  strata_factor <- droplevels(factor(strata))
  pieces <- vector("list", nlevels(strata_factor))
  names(pieces) <- levels(strata_factor)

  for (level in levels(strata_factor)) {
    index <- which(strata_factor == level & weights > 0)
    if (!length(index)) {
      stop("Every censoring stratum must have positive total weight.", call. = FALSE)
    }
    km <- calculateKM(
      t = t[index],
      d = as.integer(censoring.event[index]),
      w = as.numeric(weights[index]),
      strata = rep.int(1L, length(index)),
      error = "none",
      count.type = "numeric"
    )
    # Actual failures leave before the censoring jump under T <= C. Do not
    # confuse administrative horizon completion with an observed failure.
    event_count <- numeric(length(km$time))
    counts <- rowsum(weights[index] * (epsilon[index] != code.censoring),
                    match(t[index], km$time), reorder = FALSE)
    event_count[as.integer(rownames(counts))] <- counts[, 1L]
    censor_risk <- as.numeric(km$n.risk)
    if (tie.order == "event_first") censor_risk <- censor_risk - event_count
    hazard <- numeric(length(censor_risk))
    positive <- km$n.event > 0
    if (any(censor_risk[positive] <= 0))
      stop("Censoring events require a positive censoring risk set.", call. = FALSE)
    hazard[positive] <- km$n.event[positive] / censor_risk[positive]
    survival_right <- if (tie.order == "censoring_first") as.numeric(km$surv) else
      cumprod(1 - hazard)
    survival_left <- c(1, utils::head(survival_right, -1L))
    pieces[[level]] <- data.frame(
      stratum = level,
      time = as.numeric(km$time),
      n.risk = as.numeric(km$n.risk),
      n.risk.censoring = censor_risk,
      n.event = event_count,
      n.censor = as.numeric(km$n.event),
      hazard = as.numeric(hazard),
      survival.left = survival_left,
      survival.right = survival_right,
      stringsAsFactors = FALSE
    )
  }

  table <- do.call(rbind, pieces)
  rownames(table) <- NULL
  table$survival.left.raw <- table$survival.left
  table$survival.right.raw <- table$survival.right
  table$survival.left.used <- if (is.null(prob.truncation)) {
    table$survival.left.raw
  } else {
    pmax(table$survival.left.raw, prob.truncation)
  }
  table$survival.right.used <- if (is.null(prob.truncation)) {
    table$survival.right.raw
  } else {
    pmax(table$survival.right.raw, prob.truncation)
  }
  table$truncated.left <- table$survival.left.raw < table$survival.left.used
  table$truncated.right <- table$survival.right.raw < table$survival.right.used
  low <- table$survival.left <= prob.bound | table$survival.right <= prob.bound

  structure(
    list(
      table = table,
      hazard.table = table[table$n.censor > 0, , drop = FALSE],
      strata.levels = levels(strata_factor),
      code.censoring = code.censoring,
      prob.bound = prob.bound,
      prob.truncation = prob.truncation,
      positivity.warning = any(low),
      minimum.survival = min(table$survival.right),
      minimum.survival.used = min(table$survival.right.used),
      truncation.count = sum(table$truncated.left | table$truncated.right),
      n = n,
      weights = as.numeric(weights),
      censoring.event = censoring.event,
      tie.order = tie.order
    ),
    class = c("ciftest_censoring_km", "cif_censoring_km")
  )
}

#' Evaluate a censoring Kaplan-Meier distribution
#'
#' @keywords internal
predict_censoring_km <- function(
    object,
    time,
    strata = rep.int("pooled", length(time)),
    side = c("left", "right"),
    use = c("used", "raw")
) {
  if (!inherits(object, "ciftest_censoring_km")) {
    stop("`object` must be a censoring KM fit.", call. = FALSE)
  }
  side <- match.arg(side)
  use <- match.arg(use)
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0) || anyNA(strata))
    stop("Censoring prediction requires finite non-negative times and nonmissing strata.", call. = FALSE)
  if (length(strata) != length(time)) {
    stop("`time` and `strata` must have the same length.", call. = FALSE)
  }
  strata <- as.character(strata)
  if (any(!strata %in% object$strata.levels)) {
    stop("Prediction contains an unknown censoring stratum.", call. = FALSE)
  }

  answer <- numeric(length(time))
  for (level in unique(strata)) {
    target <- which(strata == level)
    tab <- object$table[object$table$stratum == level, , drop = FALSE]
    position <- findInterval(time[target], tab$time)
    if (identical(side, "left")) {
      tied <- position > 0L
      tied[tied] <- tab$time[position[tied]] == time[target][tied]
      position[tied] <- position[tied] - 1L
    }
    value <- rep.int(1, length(target))
    has_history <- position > 0L
    survival_column <- if (is.null(object$prob.truncation)) {
      "survival.right"
    } else {
      paste0("survival.right.", use)
    }
    value[has_history] <- tab[[survival_column]][position[has_history]]
    answer[target] <- value
  }
  answer
}

#' Censoring-estimation correction to left-limit IPCW regression scores
#'
#' Returns q_i such that IF_beta_i = n H^-1 (s_i - mean(s) + q_i).
#' Uses the same source cohort as KM, with unit sampling weights. The supplied
#' scores already contain IPCW, including zero scores for unknown outcomes.
#' This product-limit derivative includes the (1-h) factor at tied jumps.
#' No n-by-n influence matrix is constructed: time-grouped score sums and
#' cumulative risk contributions give O(n log n + n*p) work and O(n*p) storage.
#' @keywords internal
censoring_km_score_correction <- function(object, t, epsilon, strata, evaluation.time, score) {
  if (!inherits(object, "cif_censoring_km")) stop("object must be a censoring KM fit")
  n <- length(t)
  if (n != object$n || length(epsilon) != n || length(strata) != n ||
      length(evaluation.time) != n || !is.matrix(score) || nrow(score) != n ||
      !ncol(score) || !n || anyNA(epsilon) || anyNA(strata) ||
      any(!is.finite(t) | t < 0) || any(!is.finite(evaluation.time) | evaluation.time < 0) ||
      any(!is.finite(score))) stop("KM score correction requires aligned finite analysis rows")
  if (any(object$weights != 1)) stop("KM score correction requires unit sampling weights")
  if (!is.null(object$prob.truncation)) stop("KM score correction does not support survival truncation")
  if (!object$tie.order %in% c("event_first", "censoring_first")) stop("Unknown KM tie order")
  strata <- as.character(strata)
  if (!setequal(unique(strata), object$strata.levels)) stop("KM score correction strata do not match")
  if (any(object$censoring.event & epsilon != object$code.censoring))
    stop("KM score correction censoring indicators do not match")
  correction <- matrix(0, n, ncol(score), dimnames = dimnames(score))
  for (z in object$strata.levels) {
    ix <- which(strata == z)
    jumps <- object$hazard.table
    jumps <- jumps[jumps$stratum == z & jumps$time < max(evaluation.time[ix]), , drop = FALSE]
    if (!nrow(jumps)) next
    if (any(jumps$n.risk.censoring <= 0 | jumps$hazard <= 0 | jumps$hazard >= 1))
      stop("KM score correction requires positive survival at all evaluated left limits")
    times <- sort(unique(evaluation.time[ix]))
    grouped <- rowsum(score[ix, , drop = FALSE], match(evaluation.time[ix], times))
    cumulative <- matrix(apply(grouped, 2L, cumsum), nrow = length(times))
    prefix <- rbind(rep(0, ncol(score)), cumulative)
    # A(u) = sum_j s_j I(v_j > u, z_j=z), with a strict inequality for G(v-).
    A <- matrix(colSums(score[ix, , drop = FALSE]), nrow(jumps), ncol(score), byrow = TRUE) -
      prefix[findInterval(jumps$time, times) + 1L, , drop = FALSE]
    B <- A / (jumps$n.risk.censoring * (1 - jumps$hazard))
    risk.cumulative <- matrix(apply(B * jumps$hazard, 2L, cumsum), nrow = nrow(jumps))
    risk.prefix <- rbind(rep(0, ncol(score)), risk.cumulative)
    position <- findInterval(t[ix], jumps$time)
    if (object$tie.order == "event_first") {
      # Actual failures leave before censoring at a shared time; administrative
      # censoring-status completion stays at risk, even without a censor jump.
      tied <- position > 0L
      tied[tied] <- jumps$time[position[tied]] == t[ix][tied]
      position[tied & epsilon[ix] != object$code.censoring] <-
        position[tied & epsilon[ix] != object$code.censoring] - 1L
    }
    correction[ix, ] <- -risk.prefix[position + 1L, , drop = FALSE]
    at.jump <- match(t[ix], jumps$time)
    censored <- object$censoring.event[ix] & !is.na(at.jump)
    correction[ix[censored], ] <- correction[ix[censored], , drop = FALSE] +
      B[at.jump[censored], , drop = FALSE]
  }
  correction
}

