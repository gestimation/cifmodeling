# Internal censoring and Fine-Gray score infrastructure --------------------

#' Estimate a stratified censoring Kaplan-Meier distribution
#'
#' Uses the package Kaplan-Meier engine with censoring as the event. Both left
#' and right limits are retained because Fine-Gray weights must be predictable
#' at tied event/censoring times.
#'
#' @keywords internal
estimate_censoring_km <- function(
    t,
    epsilon,
    code.censoring = 0L,
    strata = rep.int("pooled", length(t)),
    weights = rep.int(1, length(t)),
    prob.bound = 1e-7
) {
  n <- length(t)
  if (length(epsilon) != n || length(strata) != n || length(weights) != n) {
    stop("Censoring KM inputs must have the same length.", call. = FALSE)
  }
  if (!n || any(!is.finite(t)) || any(t < 0)) {
    stop("Censoring KM requires finite non-negative follow-up times.", call. = FALSE)
  }
  if (anyNA(epsilon) || anyNA(strata)) {
    stop("Censoring KM inputs must not contain missing status or strata.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Censoring KM weights must be finite and non-negative.", call. = FALSE)
  }
  if (!is.numeric(prob.bound) || length(prob.bound) != 1L ||
      !is.finite(prob.bound) || prob.bound <= 0 || prob.bound >= 1) {
    stop("`prob.bound` must be one number strictly between zero and one.",
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
      d = as.integer(epsilon[index] == code.censoring),
      w = as.numeric(weights[index]),
      strata = rep.int(1L, length(index)),
      error = "none",
      count.type = "numeric"
    )
    hazard <- km$n.event / km$n.risk
    survival_right <- as.numeric(km$surv)
    survival_left <- c(1, utils::head(survival_right, -1L))
    pieces[[level]] <- data.frame(
      stratum = level,
      time = as.numeric(km$time),
      n.risk = as.numeric(km$n.risk),
      n.censor = as.numeric(km$n.event),
      hazard = as.numeric(hazard),
      survival.left = survival_left,
      survival.right = survival_right,
      stringsAsFactors = FALSE
    )
  }

  table <- do.call(rbind, pieces)
  rownames(table) <- NULL
  low <- table$survival.left <= prob.bound | table$survival.right <= prob.bound

  structure(
    list(
      table = table,
      hazard.table = table[table$n.censor > 0, , drop = FALSE],
      strata.levels = levels(strata_factor),
      code.censoring = code.censoring,
      prob.bound = prob.bound,
      positivity.warning = any(low),
      minimum.survival = min(table$survival.right),
      n = n,
      weights = as.numeric(weights)
    ),
    class = "ciftest_censoring_km"
  )
}

#' Evaluate a censoring Kaplan-Meier distribution
#'
#' @keywords internal
predict_censoring_km <- function(
    object,
    time,
    strata = rep.int("pooled", length(time)),
    side = c("left", "right")
) {
  if (!inherits(object, "ciftest_censoring_km")) {
    stop("`object` must be a censoring KM fit.", call. = FALSE)
  }
  side <- match.arg(side)
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
    value[has_history] <- tab$survival.right[position[has_history]]
    answer[target] <- value
  }
  answer
}

# Working Aalen-Johansen distributions used by the closed-form augmentation.
ciftest_working_cell <- function(exposure, strata) {
  paste(as.character(exposure), as.character(strata), sep = "\r")
}

ciftest_augmentation_cell <- function(
    exposure,
    strata.competing.risk,
    strata.censor
) {
  paste(
    as.character(exposure),
    as.character(strata.competing.risk),
    as.character(strata.censor),
    sep = "\r"
  )
}

#' Estimate working Aalen-Johansen distributions by exposure and stratum
#'
#' The package Aalen-Johansen engine is evaluated separately within every
#' observed exposure-by-stratum cell. A second call with the two event codes
#' exchanged supplies the competing-event cumulative incidence.
#'
#' @keywords internal
estimate_working_aj <- function(
    t,
    epsilon,
    exposure,
    strata = rep.int("pooled", length(t)),
    weights = rep.int(1, length(t)),
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    prob.bound = 1e-7
) {
  n <- length(t)
  if (length(epsilon) != n || length(exposure) != n ||
      length(strata) != n || length(weights) != n) {
    stop("Working AJ inputs must have the same length.", call. = FALSE)
  }
  if (!n || any(!is.finite(t)) || any(t < 0) || anyNA(epsilon) ||
      anyNA(exposure) || anyNA(strata)) {
    stop("Working AJ inputs must be finite, non-missing, and have non-negative times.",
         call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Working AJ weights must be finite and non-negative.", call. = FALSE)
  }
  allowed <- c(code.censoring, code.event1, code.event2)
  if (length(unique(allowed)) != 3L || any(!epsilon %in% allowed)) {
    stop("Working AJ event codes must be distinct and cover every status.",
         call. = FALSE)
  }
  if (!is.numeric(prob.bound) || length(prob.bound) != 1L ||
      !is.finite(prob.bound) || prob.bound <= 0 || prob.bound >= 1) {
    stop("`prob.bound` must be one number strictly between zero and one.",
         call. = FALSE)
  }

  exposure_label <- as.character(droplevels(factor(exposure)))
  strata_label <- as.character(droplevels(factor(strata)))
  cell <- ciftest_working_cell(exposure_label, strata_label)
  cells <- unique(data.frame(
    cell = cell,
    exposure = exposure_label,
    stratum = strata_label,
    stringsAsFactors = FALSE
  ))
  pieces <- vector("list", nrow(cells))
  names(pieces) <- cells$cell

  normalized <- ifelse(
    epsilon == code.event1, 1L,
    ifelse(epsilon == code.event2, 2L, 0L)
  )
  exchanged <- ifelse(normalized == 1L, 2L,
                      ifelse(normalized == 2L, 1L, 0L))

  for (k in seq_len(nrow(cells))) {
    index <- which(cell == cells$cell[k] & weights > 0)
    if (!length(index)) {
      stop("Every working AJ cell must have positive total weight.",
           call. = FALSE)
    }
    fit1 <- calculateAJ_Rcpp(
      t = t[index],
      epsilon = normalized[index],
      w = as.numeric(weights[index]),
      strata = rep.int(1L, length(index)),
      error = "none",
      return_if = FALSE
    )
    fit2 <- calculateAJ_Rcpp(
      t = t[index],
      epsilon = exchanged[index],
      w = as.numeric(weights[index]),
      strata = rep.int(1L, length(index)),
      error = "none",
      return_if = FALSE
    )
    fit_survival <- calculateKM(
      t = t[index],
      d = as.integer(normalized[index] != 0L),
      w = as.numeric(weights[index]),
      strata = rep.int(1L, length(index)),
      error = "none",
      count.type = "numeric"
    )
    if (length(fit1$time) != length(fit2$time) ||
        length(fit1$time) != length(fit_survival$time) ||
        any(fit1$time != fit2$time) ||
        any(fit1$time != fit_survival$time)) {
      stop("Working AJ event-specific fits have incompatible time grids.",
           call. = FALSE)
    }
    survival_right <- as.numeric(fit_survival$surv)
    cif1_right <- if (length(fit1$aj) == length(fit1$time)) {
      as.numeric(fit1$aj)
    } else if (any(normalized[index] == 1L)) {
      1 - as.numeric(fit1$surv)
    } else {
      numeric(length(fit1$time))
    }
    cif2_right <- if (length(fit2$aj) == length(fit2$time)) {
      as.numeric(fit2$aj)
    } else if (any(exchanged[index] == 1L)) {
      1 - as.numeric(fit2$surv)
    } else {
      numeric(length(fit2$time))
    }
    survival_left <- c(1, utils::head(survival_right, -1L))
    cif1_left <- c(0, utils::head(cif1_right, -1L))
    cif2_left <- c(0, utils::head(cif2_right, -1L))
    d_cif1 <- pmax(0, cif1_right - cif1_left)
    fg_survival_left <- 1 - cif1_left
    fg_survival_right <- 1 - cif1_right
    bad <- d_cif1 > sqrt(.Machine$double.eps) &
      (fg_survival_left <= prob.bound | fg_survival_right <= prob.bound)
    if (any(bad)) {
      stop("Working AJ subdistribution-hazard positivity is violated.",
           call. = FALSE)
    }
    d_lambda1 <- numeric(length(d_cif1))
    usable <- d_cif1 > 0
    d_lambda1[usable] <-
      log1p(-cif1_left[usable]) - log1p(-cif1_right[usable])
    coherence <- abs(survival_right + cif1_right + cif2_right - 1)
    if (max(coherence) > 1e-8) {
      stop("Working AJ probabilities do not sum to one.", call. = FALSE)
    }
    pieces[[k]] <- data.frame(
      cell = cells$cell[k],
      exposure = cells$exposure[k],
      stratum = cells$stratum[k],
      time = as.numeric(fit1$time),
      survival.left = survival_left,
      survival.right = survival_right,
      cif1.left = cif1_left,
      cif1.right = cif1_right,
      cif2.left = cif2_left,
      cif2.right = cif2_right,
      d.cif1 = d_cif1,
      d.lambda1 = d_lambda1,
      stringsAsFactors = FALSE
    )
  }

  table <- do.call(rbind, pieces)
  rownames(table) <- NULL
  structure(
    list(
      table = table,
      cells = cells,
      code.event1 = code.event1,
      code.event2 = code.event2,
      code.censoring = code.censoring,
      prob.bound = prob.bound,
      minimum.survival = min(table$survival.right),
      positivity.warning = any(table$survival.left <= prob.bound),
      n = n,
      weights = as.numeric(weights)
    ),
    class = "ciftest_working_aj"
  )
}

#' Evaluate a working Aalen-Johansen distribution
#'
#' @keywords internal
predict_working_aj <- function(
    object,
    time,
    exposure,
    strata = rep.int("pooled", length(time)),
    side = c("left", "right")
) {
  if (!inherits(object, "ciftest_working_aj")) {
    stop("`object` must be a working AJ fit.", call. = FALSE)
  }
  side <- match.arg(side)
  n <- length(time)
  if (length(exposure) != n || length(strata) != n) {
    stop("Working AJ prediction inputs must have the same length.",
         call. = FALSE)
  }
  cell <- ciftest_working_cell(exposure, strata)
  if (any(!cell %in% object$cells$cell)) {
    stop("Prediction contains an unknown working AJ cell.", call. = FALSE)
  }

  answer <- matrix(0, n, 3L,
                   dimnames = list(NULL, c("survival", "cif1", "cif2")))
  answer[, "survival"] <- 1
  for (level in unique(cell)) {
    target <- which(cell == level)
    tab <- object$table[object$table$cell == level, , drop = FALSE]
    position <- findInterval(time[target], tab$time)
    has_history <- position > 0L
    if (any(has_history)) {
      tied <- rep.int(FALSE, length(target))
      tied[has_history] <-
        tab$time[position[has_history]] == time[target][has_history]
      use_left <- identical(side, "left") & tied
      use_right <- has_history & !use_left
      if (any(use_left)) {
        answer[target[use_left], "survival"] <-
          tab$survival.left[position[use_left]]
        answer[target[use_left], "cif1"] <- tab$cif1.left[position[use_left]]
        answer[target[use_left], "cif2"] <- tab$cif2.left[position[use_left]]
      }
      if (any(use_right)) {
        answer[target[use_right], "survival"] <-
          tab$survival.right[position[use_right]]
        answer[target[use_right], "cif1"] <- tab$cif1.right[position[use_right]]
        answer[target[use_right], "cif2"] <- tab$cif2.right[position[use_right]]
      }
    }
  }
  answer
}

# Pooled Aalen-Johansen left limits used for the common FH weight process.
ciftest_null_fh_weight <- function(
    t,
    epsilon,
    code.event1,
    code.event2,
    weights,
    rho,
    gamma
) {
  target_times <- sort(unique(t[epsilon == code.event1 & weights > 0]))
  if (!length(target_times)) {
    stop("At least one event of interest is required.", call. = FALSE)
  }
  all_times <- sort(unique(t))
  survival <- 1
  cif1 <- 0
  result <- numeric(length(target_times))

  for (current in all_times) {
    target_index <- match(current, target_times)
    if (!is.na(target_index)) {
      result[target_index] <- (1 - cif1)^rho * cif1^gamma
    }
    risk <- sum(weights[t >= current])
    if (risk <= 0) next
    d1 <- sum(weights[t == current & epsilon == code.event1])
    d2 <- sum(weights[t == current & epsilon == code.event2])
    cif1 <- cif1 + survival * d1 / risk
    survival <- survival * (1 - (d1 + d2) / risk)
  }
  result
}

#' Build individual Fine-Gray score-scale contributions
#'
#' The base component is the known-censoring-distribution martingale score.
#' The censoring component is the delta-method contribution from estimating
#' the stratified censoring hazards. Their column sums equal the plug-in score
#' and zero, respectively. The matrices are diagnostic infrastructure for the
#' augmented branch; they do not replace Gray's classical covariance.
#'
#' @keywords internal
build_fg_score_iid <- function(
    t,
    epsilon,
    x,
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    strata = rep.int("pooled", length(t)),
    weights = rep.int(1, length(t)),
    rho = 0,
    gamma = 0,
    fh.weight = NULL,
    censoring = NULL,
    prob.bound = 1e-7
) {
  n <- length(t)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) != n || ncol(x) < 1L || length(epsilon) != n ||
      length(strata) != n || length(weights) != n) {
    stop("Fine-Gray score inputs have incompatible dimensions.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("Fine-Gray score design values must be finite and non-missing.",
         call. = FALSE)
  }
  allowed <- c(code.censoring, code.event1, code.event2)
  if (anyNA(epsilon) || any(!epsilon %in% allowed)) {
    stop("Fine-Gray score status contains an unsupported event code.",
         call. = FALSE)
  }
  if (length(unique(allowed)) != 3L) {
    stop("Fine-Gray event codes must be distinct.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Fine-Gray score weights must be finite and non-negative.",
         call. = FALSE)
  }
  if (any(!is.finite(t)) || any(t < 0)) {
    stop("Fine-Gray score times must be finite and non-negative.", call. = FALSE)
  }
  if (anyNA(strata)) {
    stop("Fine-Gray censoring strata must not be missing.", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0 ||
      !is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("Fine-Gray `rho` and `gamma` must be finite non-negative numbers.",
         call. = FALSE)
  }
  if (!is.numeric(prob.bound) || length(prob.bound) != 1L ||
      !is.finite(prob.bound) || prob.bound <= 0 || prob.bound >= 1) {
    stop("`prob.bound` must be one number strictly between zero and one.",
         call. = FALSE)
  }

  strata_factor <- droplevels(factor(strata))
  strata_label <- as.character(strata_factor)
  if (is.null(censoring)) {
    censoring <- estimate_censoring_km(
      t = t,
      epsilon = epsilon,
      code.censoring = code.censoring,
      strata = strata_factor,
      weights = weights,
      prob.bound = prob.bound
    )
  }

  event_times <- sort(unique(t[epsilon == code.event1 & weights > 0]))
  if (!length(event_times)) {
    stop("At least one event of interest is required.", call. = FALSE)
  }
  if (is.null(fh.weight)) {
    fh_weight <- ciftest_null_fh_weight(
      t, epsilon, code.event1, code.event2, weights, rho, gamma
    )
  } else {
    if (!is.numeric(fh.weight) || length(fh.weight) != length(event_times) ||
        any(!is.finite(fh.weight)) || any(fh.weight < 0)) {
      stop("`fh.weight` must contain one finite non-negative value per event time.",
           call. = FALSE)
    }
    fh_weight <- as.numeric(fh.weight)
  }
  p <- ncol(x)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))
  base_iid <- matrix(0, n, p, dimnames = list(NULL, score_names))
  censor_iid <- matrix(0, n, p, dimnames = list(NULL, score_names))
  event_score <- numeric(p)
  xbar <- matrix(NA_real_, length(event_times), p,
                 dimnames = list(NULL, score_names))
  risk_total <- event_count <- numeric(length(event_times))

  censor_hazard <- censoring$hazard.table
  censor_derivative <- matrix(0, nrow(censor_hazard), p,
                              dimnames = list(NULL, score_names))

  g_at_competing <- predict_censoring_km(
    censoring,
    time = t,
    strata = strata_label,
    side = "left"
  )

  for (j in seq_along(event_times)) {
    current <- event_times[j]
    is_event <- epsilon == code.event1 & t == current
    is_competing_before <- epsilon == code.event2 & t < current
    subrisk <- t >= current | is_competing_before
    censor_weight <- rep.int(1, n)

    if (any(is_competing_before)) {
      g_current <- predict_censoring_km(
        censoring,
        time = rep.int(current, sum(is_competing_before)),
        strata = strata_label[is_competing_before],
        side = "left"
      )
      denominator <- g_at_competing[is_competing_before]
      if (any(denominator <= prob.bound) || any(g_current <= prob.bound)) {
        stop("Censoring positivity is violated in the Fine-Gray risk set.",
             call. = FALSE)
      }
      censor_weight[is_competing_before] <- g_current / denominator
    }

    weighted_risk <- weights * subrisk * censor_weight
    s0 <- sum(weighted_risk)
    if (!is.finite(s0) || s0 <= prob.bound) {
      stop("Fine-Gray subdistribution risk set is empty.", call. = FALSE)
    }
    mean_x <- colSums(x * weighted_risk) / s0
    centered_x <- sweep(x, 2L, mean_x, "-")
    d1 <- sum(weights[is_event])
    if (d1 <= 0) next
    hazard1 <- d1 / s0
    a <- fh_weight[j]

    coefficient <- weights * (
      as.numeric(is_event) - as.numeric(subrisk) * censor_weight * hazard1
    )
    base_iid <- base_iid + centered_x * (a * coefficient)
    event_score <- event_score +
      a * colSums(x[is_event, , drop = FALSE] * weights[is_event]) -
      a * d1 * mean_x
    xbar[j, ] <- mean_x
    risk_total[j] <- s0
    event_count[j] <- d1

    if (nrow(censor_hazard) && any(is_competing_before)) {
      eligible_hazard <- which(censor_hazard$time < current)
      for (k in eligible_hazard) {
        affected <- is_competing_before &
          strata_label == censor_hazard$stratum[k] &
          t <= censor_hazard$time[k]
        if (!any(affected)) next
        one_minus_hazard <- 1 - censor_hazard$hazard[k]
        if (one_minus_hazard <= prob.bound) {
          stop("Censoring positivity is violated after a competing event.",
               call. = FALSE)
        }
        derivative_weight <- -censor_weight[affected] / one_minus_hazard
        derivative_mean <- colSums(
          centered_x[affected, , drop = FALSE] *
            (weights[affected] * derivative_weight)
        ) / s0
        censor_derivative[k, ] <- censor_derivative[k, ] -
          a * d1 * derivative_mean
      }
    }
  }

  if (nrow(censor_hazard)) {
    for (k in seq_len(nrow(censor_hazard))) {
      in_stratum <- strata_label == censor_hazard$stratum[k]
      at_risk <- in_stratum & t >= censor_hazard$time[k]
      censored <- in_stratum & t == censor_hazard$time[k] &
        epsilon == code.censoring
      martingale <- weights * (
        as.numeric(censored) - censor_hazard$hazard[k] * as.numeric(at_risk)
      )
      censor_iid <- censor_iid + outer(
        martingale / censor_hazard$n.risk[k],
        censor_derivative[k, ]
      )
    }
  }

  score <- colSums(base_iid)
  names(score) <- score_names
  score_iid <- base_iid + censor_iid
  list(
    score = score,
    event.score = stats::setNames(event_score, score_names),
    score.iid = score_iid,
    score.iid.base = base_iid,
    score.iid.censor = censor_iid,
    censoring = censoring,
    event.time = event_times,
    fh.weight = unname(fh_weight),
    xbar = xbar,
    risk.total = risk_total,
    event.count = event_count,
    censor.derivative = censor_derivative,
    diagnostics = list(
      score.decomposition.error = max(abs(score - event_score)),
      censor.centering.error = max(abs(colSums(censor_iid))),
      minimum.censoring.survival = censoring$minimum.survival,
      positivity.warning = censoring$positivity.warning
    )
  )
}

# Reverse cumulative sum applied separately to every matrix column.
ciftest_reverse_cumsum <- function(x) {
  if (!nrow(x)) return(x)
  apply(x, 2L, function(column) rev(cumsum(rev(column)))) |>
    matrix(nrow = nrow(x), ncol = ncol(x),
           dimnames = dimnames(x))
}

#' Build the closed-form competing-risk augmentation
#'
#' This computes the working-model process
#' deqn{H_a(t,X)=\int_t^\tau a(u)(X-e_X(u))G_C(u)d\Lambda_1(u|X)}
#' in each exposure-by-competing-risk-by-censoring nuisance cell, then
#' evaluates the centered censoring-martingale projection. The returned matrix
#' is on the score scale; its column sums are the closed-form augmentation
#' score.
#'
#' @keywords internal
build_closed_form_augmentation <- function(
    base,
    working,
    t,
    epsilon,
    x,
    exposure,
    strata.censor = rep.int("pooled", length(t)),
    strata.competing.risk = rep.int("pooled", length(t)),
    weights = rep.int(1, length(t)),
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    prob.bound = 1e-7
) {
  if (!is.list(base) || is.null(base$censoring) ||
      !inherits(working, "ciftest_working_aj")) {
    stop("Closed-form augmentation requires Fine-Gray and working AJ fits.",
         call. = FALSE)
  }
  n <- length(t)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) != n || length(epsilon) != n || length(exposure) != n ||
      length(strata.censor) != n || length(strata.competing.risk) != n ||
      length(weights) != n) {
    stop("Closed-form augmentation inputs have incompatible dimensions.",
         call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x)) || anyNA(exposure) ||
      anyNA(strata.censor) || anyNA(strata.competing.risk) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("Closed-form augmentation inputs must be finite and non-missing.",
         call. = FALSE)
  }
  if (!is.numeric(prob.bound) || length(prob.bound) != 1L ||
      !is.finite(prob.bound) || prob.bound <= 0 || prob.bound >= 1) {
    stop("`prob.bound` must be one number strictly between zero and one.",
         call. = FALSE)
  }

  event_time <- as.numeric(base$event.time)
  p <- ncol(x)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))
  if (length(base$fh.weight) != length(event_time) ||
      !all(dim(base$xbar) == c(length(event_time), p))) {
    stop("The Fine-Gray fit has an incompatible event-time process.",
         call. = FALSE)
  }
  # H_a is constant within an exposure-by-competing-risk-by-censoring cell
  # because the exposure design row and both nuisance distributions are then
  # cell-specific.
  exposure_label <- as.character(exposure)
  censor_label <- as.character(strata.censor)
  competing_risk_label <- as.character(strata.competing.risk)
  working_cell <- ciftest_working_cell(exposure_label, competing_risk_label)
  augmentation_cell <- ciftest_augmentation_cell(
    exposure_label,
    competing_risk_label,
    censor_label
  )
  augmentation_cells <- unique(data.frame(
    cell = augmentation_cell,
    working.cell = working_cell,
    exposure = exposure_label,
    competing.risk.stratum = competing_risk_label,
    censor.stratum = censor_label,
    stringsAsFactors = FALSE
  ))
  if (any(!augmentation_cells$working.cell %in% working$cells$cell)) {
    stop("Closed-form augmentation contains an unknown working AJ cell.",
         call. = FALSE)
  }

  h_process <- vector("list", nrow(augmentation_cells))
  names(h_process) <- augmentation_cells$cell
  for (k in seq_len(nrow(augmentation_cells))) {
    level <- augmentation_cells$cell[k]
    member <- which(augmentation_cell == level & weights > 0)
    if (!length(member)) {
      stop("Every augmentation nuisance cell must have positive total weight.",
           call. = FALSE)
    }
    x_cell <- x[member[1L], ]
    if (length(member) > 1L &&
        max(abs(sweep(x[member, , drop = FALSE], 2L, x_cell, "-"))) > 1e-10) {
      stop("The exposure design must be constant within each working AJ cell.",
           call. = FALSE)
    }
    tab <- working$table[
      working$table$cell == augmentation_cells$working.cell[k],
      ,
      drop = FALSE
    ]
    d_lambda1 <- numeric(length(event_time))
    matched <- match(event_time, tab$time)
    present <- !is.na(matched)
    d_lambda1[present] <- tab$d.lambda1[matched[present]]
    g_left <- predict_censoring_km(
      base$censoring,
      time = event_time,
      strata = rep.int(
        augmentation_cells$censor.stratum[k],
        length(event_time)
      ),
      side = "left"
    )
    if (any(g_left[d_lambda1 > 0] <= prob.bound)) {
      stop("Censoring positivity is violated while constructing H_a.",
           call. = FALSE)
    }
    increment <- sweep(base$xbar, 2L, x_cell, function(mean, value) value - mean)
    increment <- increment * (base$fh.weight * g_left * d_lambda1)
    h_process[[k]] <- ciftest_reverse_cumsum(increment)
    colnames(h_process[[k]]) <- score_names
  }

  augment_iid <- matrix(0, n, p,
                        dimnames = list(NULL, score_names))
  censor_hazard <- base$censoring$hazard.table
  max_raw_center_error <- 0
  min_working_survival <- 1
  min_censoring_survival <- 1

  if (nrow(censor_hazard)) {
    for (k in seq_len(nrow(censor_hazard))) {
      current <- censor_hazard$time[k]
      stratum <- censor_hazard$stratum[k]
      in_stratum <- censor_label == stratum
      at_risk <- in_stratum & t >= current
      censored <- in_stratum & t == current & epsilon == code.censoring
      risk_weight <- sum(weights[at_risk])
      if (risk_weight <= prob.bound) next
      martingale <- weights * (
        as.numeric(censored) - censor_hazard$hazard[k] * as.numeric(at_risk)
      )

      hstar <- matrix(0, n, p, dimnames = list(NULL, score_names))
      relevant <- which(in_stratum)
      working_left <- predict_working_aj(
        working,
        time = rep.int(current, length(relevant)),
        exposure = exposure_label[relevant],
        strata = competing_risk_label[relevant],
        side = "left"
      )
      g_left <- predict_censoring_km(
        base$censoring,
        time = current,
        strata = stratum,
        side = "left"
      )
      min_working_survival <- min(min_working_survival,
                                  working_left[, "survival"])
      min_censoring_survival <- min(min_censoring_survival, g_left)

      h_index <- which(event_time >= current)[1L]
      if (!is.na(h_index)) {
        for (level in unique(augmentation_cell[relevant])) {
          target <- relevant[augmentation_cell[relevant] == level]
          h_value <- h_process[[level]][h_index, ]
          local <- match(target, relevant)
          numerator_active <- working_left[local, "cif2"] > 0 &
            any(abs(h_value) > sqrt(.Machine$double.eps))
          if (any(numerator_active) &&
              (any(working_left[local[numerator_active], "survival"] <= prob.bound) ||
               g_left <= prob.bound)) {
            stop("Working-model positivity is violated in the closed-form augmentation.",
                 call. = FALSE)
          }
          ratio <- numeric(length(target))
          usable <- rep.int(
            working_left[local[1L], "cif2"] > 0 &&
              any(abs(h_value) > sqrt(.Machine$double.eps)),
            length(target)
          )
          ratio[usable] <- working_left[local[usable], "cif2"] /
            (working_left[local[usable], "survival"] * g_left)
          hstar[target, ] <- outer(ratio, h_value)
        }
      }

      hbar <- colSums(hstar[at_risk, , drop = FALSE] * weights[at_risk]) /
        risk_weight
      centered <- sweep(hstar, 2L, hbar, "-")
      raw_increment <- colSums(hstar * martingale)
      centered_increment <- colSums(centered * martingale)
      max_raw_center_error <- max(
        max_raw_center_error,
        max(abs(raw_increment - centered_increment))
      )
      augment_iid <- augment_iid + centered * martingale
    }
  }

  score <- colSums(augment_iid)
  names(score) <- score_names
  list(
    score = score,
    score.iid.augment = augment_iid,
    working.aj = working,
    h.process = h_process,
    augmentation.cells = augmentation_cells,
    diagnostics = list(
      augmentation.centering.error = max_raw_center_error,
      minimum.working.survival = min_working_survival,
      minimum.censoring.survival = min_censoring_survival,
      positivity.warning = working$positivity.warning ||
        base$censoring$positivity.warning
    )
  )
}
