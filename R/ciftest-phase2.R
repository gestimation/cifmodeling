# Internal censoring and Fine-Gray score infrastructure --------------------

ciftest_use_cpp_kernel <- function(symbol) {
  engine <- getOption("cifmodeling.ciftest.engine", "auto")
  if (!engine %in% c("auto", "R", "Rcpp")) {
    stop(
      "Option `cifmodeling.ciftest.engine` must be 'auto', 'R', or 'Rcpp'.",
      call. = FALSE
    )
  }
  available <- exists(symbol, inherits = TRUE)
  if (identical(engine, "Rcpp") && !available) {
    stop("Requested ciftest Rcpp kernel is not available: ", symbol,
         call. = FALSE)
  }
  !identical(engine, "R") && available
}

# Select the independently switchable Fine-Gray C++ implementation.  The
# package-wide option above still controls R versus C++; this option only
# chooses between the two C++ algorithms and makes rollback immediate.
ciftest_fg_cpp_variant <- function(multi = FALSE) {
  variant <- getOption("cifmodeling.ciftest.fg.engine", "prefix")
  if (!is.character(variant) || length(variant) != 1L ||
      !variant %in% c("legacy", "prefix", "check")) {
    stop(
      paste0(
        "Option `cifmodeling.ciftest.fg.engine` must be ",
        "'legacy', 'prefix', or 'check'."
      ),
      call. = FALSE
    )
  }
  if (!identical(variant, "legacy")) {
    symbol <- if (multi) {
      "_cifmodeling_ciftest_fg_iid_prefix_multi_kernel_cpp"
    } else {
      "_cifmodeling_ciftest_fg_iid_prefix_kernel_cpp"
    }
    if (!exists(symbol, inherits = TRUE)) {
      stop("Requested ciftest prefix kernel is not available: ", symbol,
           call. = FALSE)
    }
  }
  variant
}

ciftest_compare_fg_kernels <- function(prefix, legacy, tolerance = 1e-9) {
  fields <- c(
    "score", "event.score", "score.iid.base", "score.iid.censor",
    "xbar", "risk.total", "event.count", "censor.derivative"
  )
  errors <- vapply(fields, function(field) {
    actual <- as.numeric(prefix[[field]])
    reference <- as.numeric(legacy[[field]])
    if (length(actual) != length(reference) ||
        !identical(is.na(actual), is.na(reference))) {
      return(Inf)
    }
    keep <- !is.na(reference)
    if (!any(keep)) return(0)
    max(abs(actual[keep] - reference[keep])) /
      max(1, max(abs(reference[keep])))
  }, numeric(1L))
  if (any(!is.finite(errors)) || max(errors) > tolerance) {
    stop(
      "Fine-Gray prefix/legacy kernel check failed; largest scaled error = ",
      format(max(errors), digits = 6L), ".",
      call. = FALSE
    )
  }
  invisible(errors)
}

ciftest_run_fg_cpp_kernel <- function(arguments, multi = FALSE) {
  variant <- ciftest_fg_cpp_variant(multi)
  legacy_fun <- if (multi) {
    ciftest_fg_iid_multi_kernel_cpp
  } else {
    ciftest_fg_iid_kernel_cpp
  }
  prefix_fun <- if (multi) {
    ciftest_fg_iid_prefix_multi_kernel_cpp
  } else {
    ciftest_fg_iid_prefix_kernel_cpp
  }
  if (identical(variant, "legacy")) {
    return(list(value = do.call(legacy_fun, arguments), engine = "legacy"))
  }
  prefix <- do.call(prefix_fun, arguments)
  if (identical(variant, "check")) {
    legacy <- do.call(legacy_fun, arguments)
    tolerance <- getOption("cifmodeling.ciftest.fg.check.tolerance", 1e-9)
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance <= 0) {
      stop(
        "Option `cifmodeling.ciftest.fg.check.tolerance` must be positive.",
        call. = FALSE
      )
    }
    ciftest_compare_fg_kernels(prefix, legacy, tolerance)
  }
  list(value = prefix, engine = variant)
}

ciftest_augmentation_cpp_variant <- function(multi = FALSE) {
  variant <- getOption(
    "cifmodeling.ciftest.augmentation.engine", "prefix"
  )
  if (!is.character(variant) || length(variant) != 1L ||
      !variant %in% c("legacy", "prefix", "check")) {
    stop(
      paste0(
        "Option `cifmodeling.ciftest.augmentation.engine` must be ",
        "'legacy', 'prefix', or 'check'."
      ),
      call. = FALSE
    )
  }
  if (!identical(variant, "legacy")) {
    symbol <- if (multi) {
      "_cifmodeling_ciftest_augmentation_iid_prefix_multi_kernel_cpp"
    } else {
      "_cifmodeling_ciftest_augmentation_iid_prefix_kernel_cpp"
    }
    if (!exists(symbol, inherits = TRUE)) {
      stop("Requested ciftest augmentation prefix kernel is unavailable: ",
           symbol, call. = FALSE)
    }
  }
  variant
}

ciftest_compare_augmentation_kernels <- function(
    prefix,
    legacy,
    tolerance = 1e-9
) {
  fields <- c(
    "score.iid.augment", "augmentation.centering.error",
    "minimum.working.survival", "minimum.censoring.survival"
  )
  errors <- vapply(fields, function(field) {
    actual <- as.numeric(prefix[[field]])
    reference <- as.numeric(legacy[[field]])
    if (length(actual) != length(reference) ||
        !identical(is.na(actual), is.na(reference))) {
      return(Inf)
    }
    keep <- !is.na(reference)
    if (!any(keep)) return(0)
    max(abs(actual[keep] - reference[keep])) /
      max(1, max(abs(reference[keep])))
  }, numeric(1L))
  if (any(!is.finite(errors)) || max(errors) > tolerance) {
    stop(
      paste0(
        "Augmentation prefix/legacy kernel check failed; ",
        "largest scaled error = "
      ),
      format(max(errors), digits = 6L), ".",
      call. = FALSE
    )
  }
  invisible(errors)
}

ciftest_run_augmentation_cpp_kernel <- function(arguments, multi = FALSE) {
  variant <- ciftest_augmentation_cpp_variant(multi)
  legacy_fun <- if (multi) {
    ciftest_augmentation_iid_multi_kernel_cpp
  } else {
    ciftest_augmentation_iid_kernel_cpp
  }
  prefix_fun <- if (multi) {
    ciftest_augmentation_iid_prefix_multi_kernel_cpp
  } else {
    ciftest_augmentation_iid_prefix_kernel_cpp
  }
  if (identical(variant, "legacy")) {
    return(list(value = do.call(legacy_fun, arguments), engine = "legacy"))
  }
  if (identical(variant, "prefix")) {
    return(list(value = do.call(prefix_fun, arguments), engine = "prefix"))
  }

  evaluate <- function(fun) {
    tryCatch(list(value = do.call(fun, arguments)), error = function(error) {
      list(error = error)
    })
  }
  prefix <- evaluate(prefix_fun)
  legacy <- evaluate(legacy_fun)
  if (!is.null(prefix$error) || !is.null(legacy$error)) {
    if (!is.null(prefix$error) && !is.null(legacy$error) &&
        identical(conditionMessage(prefix$error),
                  conditionMessage(legacy$error))) {
      stop(prefix$error)
    }
    stop(
      "Augmentation prefix/legacy kernels disagreed on whether to fail.",
      call. = FALSE
    )
  }
  tolerance <- getOption(
    "cifmodeling.ciftest.augmentation.check.tolerance", 1e-9
  )
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop(
      paste0(
        "Option `cifmodeling.ciftest.augmentation.check.tolerance` ",
        "must be positive."
      ),
      call. = FALSE
    )
  }
  ciftest_compare_augmentation_kernels(
    prefix$value, legacy$value, tolerance
  )
  list(value = prefix$value, engine = "check")
}

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
    prob.bound = 1e-7,
    prob.truncation = NULL,
    censoring.event = epsilon == code.censoring
) {
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
      censoring.event = censoring.event
    ),
    class = "ciftest_censoring_km"
  )
}

ciftest_censoring_event_indicator <- function(
    censoring,
    epsilon,
    code.censoring
) {
  indicator <- censoring$censoring.event
  if (is.null(indicator)) indicator <- epsilon == code.censoring
  if (!is.logical(indicator) || length(indicator) != length(epsilon) ||
      anyNA(indicator)) {
    stop("The censoring-event indicator is incompatible with the data.",
         call. = FALSE)
  }
  indicator
}

ciftest_kernel_status <- function(epsilon, censoring.event, code.censoring) {
  status <- as.integer(epsilon)
  administrative <- status == code.censoring & !censoring.event
  status[administrative] <- NA_integer_
  status
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

# Working Aalen-Johansen distributions used by the closed-form augmentation.
ciftest_working_cell <- function(exposure, strata) {
  paste(as.character(exposure), as.character(strata), sep = "\r")
}

ciftest_augmentation_cell <- function(
    exposure,
    strata.competing.risk,
    strata.censor,
    strata.event = NULL
) {
  values <- list(
    as.character(exposure),
    as.character(strata.competing.risk),
    as.character(strata.censor)
  )
  if (!is.null(strata.event)) {
    values <- c(list(as.character(strata.event)), values)
  }
  do.call(paste, c(values, sep = "\r"))
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
    prob.bound = 1e-7,
    prob.truncation = NULL
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

  if (!is.null(prob.truncation) &&
      (!is.numeric(prob.truncation) || length(prob.truncation) != 1L ||
       !is.finite(prob.truncation) || prob.truncation <= prob.bound ||
       prob.truncation >= 1)) {
    stop("`prob.truncation` must be NULL or strictly between `prob.bound` and one.",
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
    d_cif2 <- pmax(0, cif2_right - cif2_left)
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
      d.cif2 = d_cif2,
      d.lambda1 = d_lambda1,
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
  structure(
    list(
      table = table,
      cells = cells,
      code.event1 = code.event1,
      code.event2 = code.event2,
      code.censoring = code.censoring,
      prob.bound = prob.bound,
      prob.truncation = prob.truncation,
      minimum.survival = min(table$survival.right),
      minimum.survival.used = min(table$survival.right.used),
      truncation.count = sum(table$truncated.left | table$truncated.right),
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
    side = c("left", "right"),
    use = c("used", "raw")
) {
  if (!inherits(object, "ciftest_working_aj")) {
    stop("`object` must be a working AJ fit.", call. = FALSE)
  }
  side <- match.arg(side)
  use <- match.arg(use)
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
        survival_column <- if (is.null(object$prob.truncation)) {
          "survival.left"
        } else {
          paste0("survival.left.", use)
        }
        answer[target[use_left], "survival"] <-
          tab[[survival_column]][position[use_left]]
        answer[target[use_left], "cif1"] <- tab$cif1.left[position[use_left]]
        answer[target[use_left], "cif2"] <- tab$cif2.left[position[use_left]]
      }
      if (any(use_right)) {
        survival_column <- if (is.null(object$prob.truncation)) {
          "survival.right"
        } else {
          paste0("survival.right.", use)
        }
        answer[target[use_right], "survival"] <-
          tab[[survival_column]][position[use_right]]
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

# Evaluate the stratum-specific null FH process on a supplied common grid.
ciftest_null_fh_weight_at <- function(
    t,
    epsilon,
    code.event1,
    code.event2,
    weights,
    rho,
    gamma,
    evaluation.time
) {
  evaluation.time <- sort(unique(as.numeric(evaluation.time)))
  result <- numeric(length(evaluation.time))
  if (!length(evaluation.time)) return(result)
  all_times <- sort(unique(c(t, evaluation.time)))
  survival <- 1
  cif1 <- 0
  for (current in all_times) {
    target <- match(current, evaluation.time)
    if (!is.na(target)) {
      result[target] <- (1 - cif1)^rho * cif1^gamma
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

# Event-stratified Fine--Gray score with independently specified censoring
# strata. Event risk sets and FH processes are stratum-specific, while the
# censoring-hazard influence is accumulated over the shared censoring model.
build_fg_score_iid_event_stratified <- function(
    t,
    epsilon,
    x,
    strata.event,
    strata.censor,
    weights,
    code.event1,
    code.event2,
    code.censoring,
    rho,
    gamma,
    fh.weight = NULL,
    censoring = NULL,
    prob.bound = 1e-7,
    prob.truncation = NULL
) {
  n <- length(t)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  event_factor <- droplevels(factor(strata.event))
  censor_factor <- droplevels(factor(strata.censor))
  event_label <- as.character(event_factor)
  censor_label <- as.character(censor_factor)
  event_levels <- levels(event_factor)
  L <- length(event_levels)
  p <- ncol(x)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))

  if (nrow(x) != n || length(epsilon) != n ||
      length(strata.event) != n || length(strata.censor) != n ||
      length(weights) != n || anyNA(event_factor) || anyNA(censor_factor)) {
    stop("Event-stratified Fine-Gray inputs have incompatible dimensions or missing strata.",
         call. = FALSE)
  }
  if (is.null(censoring)) {
    censoring <- estimate_censoring_km(
      t = t,
      epsilon = epsilon,
      code.censoring = code.censoring,
      strata = censor_factor,
      weights = weights,
      prob.bound = prob.bound,
      prob.truncation = prob.truncation
    )
  }
  censoring_event <- ciftest_censoring_event_indicator(
    censoring, epsilon, code.censoring
  )
  event_times <- sort(unique(t[epsilon == code.event1 & weights > 0]))
  if (!length(event_times)) {
    stop("At least one event of interest is required.", call. = FALSE)
  }
  m <- length(event_times)
  if (is.null(fh.weight)) {
    fh_weight <- vapply(event_levels, function(level) {
      use <- event_label == level
      ciftest_null_fh_weight_at(
        t = t[use], epsilon = epsilon[use],
        code.event1 = code.event1, code.event2 = code.event2,
        weights = weights[use], rho = rho, gamma = gamma,
        evaluation.time = event_times
      )
    }, numeric(m))
  } else {
    fh_weight <- as.matrix(fh.weight)
    if (!all(dim(fh_weight) == c(m, L)) ||
        any(!is.finite(fh_weight)) || any(fh_weight < 0)) {
      stop("Stratified `fh.weight` must be an event-time by event-stratum matrix.",
           call. = FALSE)
    }
  }
  colnames(fh_weight) <- event_levels

  base_iid <- matrix(0, n, p, dimnames = list(NULL, score_names))
  censor_iid <- matrix(0, n, p, dimnames = list(NULL, score_names))
  event_score <- numeric(p)
  xbar <- array(
    NA_real_, dim = c(m, p, L),
    dimnames = list(NULL, score_names, event_levels)
  )
  risk_total <- event_count <- matrix(
    0, nrow = m, ncol = L,
    dimnames = list(NULL, event_levels)
  )
  censor_hazard <- censoring$hazard.table
  censor_derivative <- matrix(
    0, nrow = nrow(censor_hazard), ncol = p,
    dimnames = list(NULL, score_names)
  )
  g_at_competing <- predict_censoring_km(
    censoring, time = t, strata = censor_label, side = "left"
  )
  g_at_competing_raw <- predict_censoring_km(
    censoring, time = t, strata = censor_label, side = "left", use = "raw"
  )

  for (ell in seq_along(event_levels)) {
    in_event_stratum <- event_label == event_levels[ell]
    for (j in seq_along(event_times)) {
      current <- event_times[j]
      is_event <- in_event_stratum & epsilon == code.event1 & t == current
      is_competing_before <-
        in_event_stratum & epsilon == code.event2 & t < current
      subrisk <- in_event_stratum & (t >= current | is_competing_before)
      censor_weight <- rep.int(1, n)
      if (any(is_competing_before)) {
        g_current <- predict_censoring_km(
          censoring,
          time = rep.int(current, sum(is_competing_before)),
          strata = censor_label[is_competing_before],
          side = "left"
        )
        g_current_raw <- predict_censoring_km(
          censoring,
          time = rep.int(current, sum(is_competing_before)),
          strata = censor_label[is_competing_before],
          side = "left",
          use = "raw"
        )
        denominator <- g_at_competing[is_competing_before]
        denominator_raw <- g_at_competing_raw[is_competing_before]
        if (any(denominator_raw <= prob.bound) ||
            any(g_current_raw <= prob.bound)) {
          stop("Censoring positivity is violated in the Fine-Gray risk set.",
               call. = FALSE)
        }
        censor_weight[is_competing_before] <- g_current / denominator
      }
      weighted_risk <- weights * as.numeric(subrisk) * censor_weight
      s0 <- sum(weighted_risk)
      d1 <- sum(weights[is_event])
      event_count[j, ell] <- d1
      if (!is.finite(s0) || s0 <= prob.bound) {
        if (d1 > 0) {
          stop("Fine-Gray event-stratum risk set is empty at an event time.",
               call. = FALSE)
        }
        next
      }
      mean_x <- colSums(x * weighted_risk) / s0
      xbar[j, , ell] <- mean_x
      risk_total[j, ell] <- s0
      if (d1 <= 0) next
      centered_x <- sweep(x, 2L, mean_x, "-")
      hazard1 <- d1 / s0
      a <- fh_weight[j, ell]
      coefficient <- weights * (
        as.numeric(is_event) - as.numeric(subrisk) * censor_weight * hazard1
      )
      base_iid <- base_iid + centered_x * (a * coefficient)
      event_score <- event_score +
        a * colSums(x[is_event, , drop = FALSE] * weights[is_event]) -
        a * d1 * mean_x

      if (nrow(censor_hazard) && any(is_competing_before)) {
        eligible_hazard <- which(censor_hazard$time < current)
        for (k in eligible_hazard) {
          affected <- is_competing_before &
            censor_label == censor_hazard$stratum[k]
          if (!any(affected)) next
          one_minus_hazard <- 1 - censor_hazard$hazard[k]
          if (one_minus_hazard <= prob.bound) {
            stop("Censoring positivity is violated after a competing event.",
                 call. = FALSE)
          }
          current_active <- if (is.null(censoring$prob.truncation)) {
            rep.int(TRUE, sum(affected))
          } else {
            g_current_raw[
              match(which(affected), which(is_competing_before))
            ] > censoring$prob.truncation
          }
          denominator_active <- if (is.null(censoring$prob.truncation)) {
            rep.int(TRUE, sum(affected))
          } else {
            g_at_competing_raw[affected] > censoring$prob.truncation
          }
          dlog_ratio <- (
            -as.numeric(current_active) +
              as.numeric(denominator_active &
                censor_hazard$time[k] < t[affected])
          ) / one_minus_hazard
          derivative_weight <- censor_weight[affected] * dlog_ratio
          if (!any(derivative_weight != 0)) next
          derivative_mean <- colSums(
            centered_x[affected, , drop = FALSE] *
              (weights[affected] * derivative_weight)
          ) / s0
          censor_derivative[k, ] <- censor_derivative[k, ] -
            a * d1 * derivative_mean
        }
      }
    }
  }

  if (nrow(censor_hazard)) {
    for (k in seq_len(nrow(censor_hazard))) {
      in_censor_stratum <- censor_label == censor_hazard$stratum[k]
      at_risk <- in_censor_stratum & t >= censor_hazard$time[k]
      censored <- in_censor_stratum & t == censor_hazard$time[k] &
        censoring_event
      martingale <- weights * (
        as.numeric(censored) - censor_hazard$hazard[k] * as.numeric(at_risk)
      )
      censor_iid <- censor_iid + outer(
        martingale / censor_hazard$n.risk[k], censor_derivative[k, ]
      )
    }
  }
  score <- stats::setNames(colSums(base_iid), score_names)
  event_score <- stats::setNames(event_score, score_names)
  score_iid <- base_iid + censor_iid
  list(
    score = score,
    event.score = event_score,
    score.iid = score_iid,
    score.iid.base = base_iid,
    score.iid.censor = censor_iid,
    censoring = censoring,
    event.time = event_times,
    event.strata = event_levels,
    fh.weight = fh_weight,
    xbar = xbar,
    risk.total = risk_total,
    event.count = event_count,
    censor.derivative = censor_derivative,
    diagnostics = list(
      score.decomposition.error = max(abs(score - event_score)),
      censor.centering.error = max(abs(colSums(censor_iid))),
      minimum.censoring.survival = censoring$minimum.survival,
      positivity.warning = censoring$positivity.warning,
      engine = "R-event-stratified",
      event.stratified = TRUE
    )
  )
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
    prob.bound = 1e-7,
    prob.truncation = NULL,
    strata.event = rep.int("pooled", length(t))
) {
  n <- length(t)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) != n || ncol(x) < 1L || length(epsilon) != n ||
      length(strata.event) != n || length(strata) != n ||
      length(weights) != n) {
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
  if (anyNA(strata.event) || anyNA(strata)) {
    stop("Fine-Gray event and censoring strata must not be missing.",
         call. = FALSE)
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

  event_factor <- droplevels(factor(strata.event))
  event_levels <- levels(event_factor)
  event_stratified <- nlevels(event_factor) > 1L ||
    !identical(event_levels, "pooled")
  if (event_stratified) {
    return(build_fg_score_iid_event_stratified(
      t = t, epsilon = epsilon, x = x,
      strata.event = event_factor, strata.censor = strata,
      weights = weights, code.event1 = code.event1,
      code.event2 = code.event2, code.censoring = code.censoring,
      rho = rho, gamma = gamma, fh.weight = fh.weight,
      censoring = censoring, prob.bound = prob.bound,
      prob.truncation = prob.truncation
    ))
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
      prob.bound = prob.bound,
      prob.truncation = prob.truncation
    )
  }
  censoring_event <- ciftest_censoring_event_indicator(
    censoring, epsilon, code.censoring
  )
  kernel_epsilon <- ciftest_kernel_status(
    epsilon, censoring_event, code.censoring
  )

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
  g_at_competing_raw <- predict_censoring_km(
    censoring,
    time = t,
    strata = strata_label,
    side = "left",
    use = "raw"
  )

  if (is.null(censoring$prob.truncation) && ciftest_use_cpp_kernel(
    "_cifmodeling_ciftest_fg_iid_kernel_cpp"
  )) {
    strata_id <- as.integer(strata_factor)
    g_event_stratum <- vapply(
      levels(strata_factor),
      function(level) {
        predict_censoring_km(
          censoring,
          time = event_times,
          strata = rep.int(level, length(event_times)),
          side = "left"
        )
      },
      numeric(length(event_times))
    )
    g_event_stratum <- matrix(
      g_event_stratum,
      nrow = length(event_times),
      ncol = nlevels(strata_factor)
    )
    hazard_stratum <- match(
      as.character(censor_hazard$stratum),
      levels(strata_factor)
    )
    kernel_run <- ciftest_run_fg_cpp_kernel(list(
      t = as.numeric(t),
      epsilon = kernel_epsilon,
      x = x,
      weights = as.numeric(weights),
      strata_id = strata_id,
      code_event1 = as.integer(code.event1),
      code_event2 = as.integer(code.event2),
      code_censoring = as.integer(code.censoring),
      event_times = event_times,
      fh_weight = fh_weight,
      g_at_competing = g_at_competing,
      g_event_stratum = g_event_stratum,
      hazard_time = as.numeric(censor_hazard$time),
      hazard_stratum = as.integer(hazard_stratum),
      hazard = as.numeric(censor_hazard$hazard),
      hazard_n_risk = as.numeric(censor_hazard$n.risk),
      prob_bound = prob.bound
    ))
    kernel <- kernel_run$value
    base_iid <- kernel$score.iid.base
    censor_iid <- kernel$score.iid.censor
    xbar <- kernel$xbar
    censor_derivative <- kernel$censor.derivative
    dimnames(base_iid) <- dimnames(censor_iid) <- list(NULL, score_names)
    colnames(xbar) <- score_names
    colnames(censor_derivative) <- score_names
    score <- stats::setNames(as.numeric(kernel$score), score_names)
    event_score <- stats::setNames(
      as.numeric(kernel$event.score), score_names
    )
    score_iid <- base_iid + censor_iid
    return(list(
      score = score,
      event.score = event_score,
      score.iid = score_iid,
      score.iid.base = base_iid,
      score.iid.censor = censor_iid,
      censoring = censoring,
      event.time = event_times,
      fh.weight = unname(fh_weight),
      xbar = xbar,
      risk.total = as.numeric(kernel$risk.total),
      event.count = as.numeric(kernel$event.count),
      censor.derivative = censor_derivative,
      diagnostics = list(
        score.decomposition.error = max(abs(score - event_score)),
        censor.centering.error = max(abs(colSums(censor_iid))),
        minimum.censoring.survival = censoring$minimum.survival,
        positivity.warning = censoring$positivity.warning,
        engine = if (identical(kernel_run$engine, "legacy")) {
          "Rcpp"
        } else {
          paste0("Rcpp-", kernel_run$engine)
        }
      )
    ))
  }

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
      g_current_raw <- predict_censoring_km(
        censoring,
        time = rep.int(current, sum(is_competing_before)),
        strata = strata_label[is_competing_before],
        side = "left",
        use = "raw"
      )
      denominator <- g_at_competing[is_competing_before]
      denominator_raw <- g_at_competing_raw[is_competing_before]
      if (any(denominator_raw <= prob.bound) ||
          any(g_current_raw <= prob.bound)) {
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
          strata_label == censor_hazard$stratum[k]
        if (!any(affected)) next
        one_minus_hazard <- 1 - censor_hazard$hazard[k]
        if (one_minus_hazard <= prob.bound) {
          stop("Censoring positivity is violated after a competing event.",
               call. = FALSE)
        }
        current_active <- if (is.null(censoring$prob.truncation)) {
          rep.int(TRUE, sum(affected))
        } else {
          g_current_raw[
            match(which(affected), which(is_competing_before))
          ] > censoring$prob.truncation
        }
        denominator_active <- if (is.null(censoring$prob.truncation)) {
          rep.int(TRUE, sum(affected))
        } else {
          g_at_competing_raw[affected] > censoring$prob.truncation
        }
        dlog_ratio <- (
          -as.numeric(current_active) +
            as.numeric(denominator_active &
              censor_hazard$time[k] < t[affected])
        ) / one_minus_hazard
        derivative_weight <- censor_weight[affected] * dlog_ratio
        if (!any(derivative_weight != 0)) next
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
        censoring_event
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
      positivity.warning = censoring$positivity.warning,
      engine = "R"
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
#' \deqn{H_a(t,X)=\int_t^\tau a(u)(X-e_X(u))G_C(u)d\Lambda_1(u|X)}
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
    prob.bound = 1e-7,
    strata.event = rep.int("pooled", length(t))
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
      length(strata.event) != n || length(strata.censor) != n ||
      length(strata.competing.risk) != n || length(weights) != n) {
    stop("Closed-form augmentation inputs have incompatible dimensions.",
         call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x)) || anyNA(exposure) ||
      anyNA(strata.event) || anyNA(strata.censor) ||
      anyNA(strata.competing.risk) ||
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
  censoring_event <- ciftest_censoring_event_indicator(
    base$censoring, epsilon, code.censoring
  )
  kernel_epsilon <- ciftest_kernel_status(
    epsilon, censoring_event, code.censoring
  )
  p <- ncol(x)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))
  event_stratified <- isTRUE(base$diagnostics$event.stratified)
  event_levels <- if (event_stratified) base$event.strata else "pooled"
  if (!event_stratified &&
      (length(base$fh.weight) != length(event_time) ||
       !all(dim(base$xbar) == c(length(event_time), p)))) {
    stop("The Fine-Gray fit has an incompatible event-time process.",
         call. = FALSE)
  }
  if (event_stratified &&
      (!all(dim(base$fh.weight) == c(length(event_time), length(event_levels))) ||
       !all(dim(base$xbar) == c(length(event_time), p, length(event_levels))) ||
       !all(dim(base$risk.total) == c(length(event_time), length(event_levels))))) {
    stop("The stratified Fine-Gray fit has an incompatible event-time process.",
         call. = FALSE)
  }
  # H_a is constant within an exposure-by-competing-risk-by-censoring cell
  # because the exposure design row and both nuisance distributions are then
  # cell-specific.
  exposure_label <- as.character(exposure)
  event_label <- as.character(strata.event)
  censor_label <- as.character(strata.censor)
  competing_risk_label <- as.character(strata.competing.risk)
  working_cell <- ciftest_working_cell(exposure_label, competing_risk_label)
  augmentation_cell <- ciftest_augmentation_cell(
    exposure_label,
    competing_risk_label,
    censor_label,
    strata.event = if (event_stratified) event_label else NULL
  )
  augmentation_cells <- unique(data.frame(
    cell = augmentation_cell,
    working.cell = working_cell,
    exposure = exposure_label,
    event.stratum = if (event_stratified) event_label else "pooled",
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
    if (event_stratified) {
      event_index <- match(
        augmentation_cells$event.stratum[k], event_levels
      )
      if (is.na(event_index)) {
        stop("Closed-form augmentation contains an unknown event stratum.",
             call. = FALSE)
      }
      xbar_cell <- base$xbar[, , event_index, drop = FALSE]
      dim(xbar_cell) <- c(length(event_time), p)
      active_support <- base$risk.total[, event_index] > prob.bound &
        apply(is.finite(xbar_cell), 1L, all)
      d_lambda1[!active_support] <- 0
      fh_cell <- base$fh.weight[, event_index]
    } else {
      xbar_cell <- base$xbar
      fh_cell <- base$fh.weight
    }
    g_left <- predict_censoring_km(
      base$censoring,
      time = event_time,
      strata = rep.int(
        augmentation_cells$censor.stratum[k],
        length(event_time)
      ),
      side = "left"
    )
    g_left_raw <- predict_censoring_km(
      base$censoring,
      time = event_time,
      strata = rep.int(
        augmentation_cells$censor.stratum[k],
        length(event_time)
      ),
      side = "left",
      use = "raw"
    )
    if (any(g_left_raw[d_lambda1 > 0] <= prob.bound)) {
      stop("Censoring positivity is violated while constructing H_a.",
           call. = FALSE)
    }
    increment <- sweep(xbar_cell, 2L, x_cell, function(mean, value) value - mean)
    increment <- increment * (fh_cell * g_left * d_lambda1)
    h_process[[k]] <- ciftest_reverse_cumsum(increment)
    colnames(h_process[[k]]) <- score_names
  }

  if (is.null(base$censoring$prob.truncation) &&
      is.null(working$prob.truncation) && ciftest_use_cpp_kernel(
        "_cifmodeling_ciftest_augmentation_iid_kernel_cpp"
      )) {
    cell_levels <- names(h_process)
    working_levels <- working$cells$cell
    augmentation_cell_id <- match(augmentation_cell, cell_levels)
    working_cell_id <- match(working_cell, working_levels)
    h_array <- array(
      0,
      dim = c(length(cell_levels), length(event_time), p)
    )
    for (cell_index in seq_along(cell_levels)) {
      h_array[cell_index, , ] <- h_process[[cell_index]]
    }
    censor_hazard <- base$censoring$hazard.table
    working_survival <- working_cif2 <- matrix(
      numeric(nrow(censor_hazard) * length(working_levels)),
      nrow = nrow(censor_hazard),
      ncol = length(working_levels)
    )
    if (nrow(censor_hazard)) {
      for (working_index in seq_along(working_levels)) {
        predicted <- predict_working_aj(
          working,
          time = censor_hazard$time,
          exposure = rep.int(
            working$cells$exposure[working_index],
            nrow(censor_hazard)
          ),
          strata = rep.int(
            working$cells$stratum[working_index],
            nrow(censor_hazard)
          ),
          side = "left"
        )
        working_survival[, working_index] <- predicted[, "survival"]
        working_cif2[, working_index] <- predicted[, "cif2"]
      }
    }
    censor_levels <- base$censoring$strata.levels
    kernel_run <- ciftest_run_augmentation_cpp_kernel(list(
      t = as.numeric(t),
      epsilon = kernel_epsilon,
      weights = as.numeric(weights),
      censor_stratum_id = as.integer(match(censor_label, censor_levels)),
      augmentation_cell_id = as.integer(augmentation_cell_id),
      working_cell_id = as.integer(working_cell_id),
      code_censoring = as.integer(code.censoring),
      event_times = event_time,
      hazard_time = as.numeric(censor_hazard$time),
      hazard_stratum = as.integer(match(
        as.character(censor_hazard$stratum), censor_levels
      )),
      hazard = as.numeric(censor_hazard$hazard),
      hazard_n_risk = as.numeric(censor_hazard$n.risk),
      hazard_g_left = as.numeric(censor_hazard$survival.left),
      working_survival = working_survival,
      working_cif2 = working_cif2,
      h_process = h_array,
      prob_bound = prob.bound
    ))
    kernel <- kernel_run$value
    augment_iid <- kernel$score.iid.augment
    dimnames(augment_iid) <- list(NULL, score_names)
    score <- colSums(augment_iid)
    names(score) <- score_names
    return(list(
      score = score,
      score.iid.augment = augment_iid,
      working.aj = working,
      h.process = h_process,
      augmentation.cells = augmentation_cells,
      diagnostics = list(
        augmentation.centering.error =
          kernel$augmentation.centering.error,
        minimum.working.survival = kernel$minimum.working.survival,
        minimum.censoring.survival = kernel$minimum.censoring.survival,
        positivity.warning = working$positivity.warning ||
          base$censoring$positivity.warning,
        engine = if (identical(kernel_run$engine, "legacy")) {
          "Rcpp"
        } else {
          paste0("Rcpp-", kernel_run$engine)
        }
      )
    ))
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
      censored <- in_stratum & t == current & censoring_event
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
      working_left_raw <- predict_working_aj(
        working,
        time = rep.int(current, length(relevant)),
        exposure = exposure_label[relevant],
        strata = competing_risk_label[relevant],
        side = "left",
        use = "raw"
      )
      g_left <- predict_censoring_km(
        base$censoring,
        time = current,
        strata = stratum,
        side = "left"
      )
      g_left_raw <- predict_censoring_km(
        base$censoring,
        time = current,
        strata = stratum,
        side = "left",
        use = "raw"
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
              (any(working_left_raw[
                local[numerator_active], "survival"
              ] <= prob.bound) || g_left_raw <= prob.bound)) {
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
        base$censoring$positivity.warning,
      engine = "R"
    )
  )
}

# Scheike Appendix-E fixed-point infrastructure ---------------------------

# Apply a finite-sample anchor to an AIPWCC image. At the population nuisance
# law AIPWCC(seed) equals the closed-form augmented contribution. Replacing
# the empirical seed image by the empirical closed-form contribution preserves
# that zeroth-step identity while retaining only the fixed-point increment.
ciftest_anchor_aipwcc_iid <- function(
    closed.form.iid,
    value.aipwcc.iid,
    seed.aipwcc.iid
) {
  closed.form.iid <- as.matrix(closed.form.iid)
  value.aipwcc.iid <- as.matrix(value.aipwcc.iid)
  seed.aipwcc.iid <- as.matrix(seed.aipwcc.iid)
  if (!identical(dim(closed.form.iid), dim(value.aipwcc.iid)) ||
      !identical(dim(closed.form.iid), dim(seed.aipwcc.iid))) {
    stop("Closed-form, value, and seed score-IID matrices must conform.",
         call. = FALSE)
  }
  adjustment <- closed.form.iid - seed.aipwcc.iid
  aipwcc_increment <- value.aipwcc.iid - seed.aipwcc.iid
  score_iid <- closed.form.iid + aipwcc_increment
  list(
    score.iid = score_iid,
    score.adjustment = adjustment,
    aipwcc.increment = aipwcc_increment,
    seed.score.difference =
      colSums(seed.aipwcc.iid) - colSums(closed.form.iid),
    seed.iid.rms.difference = sqrt(mean(
      (seed.aipwcc.iid - closed.form.iid)^2
    )),
    identity.error = max(abs(
      score_iid - (closed.form.iid + aipwcc_increment)
    ))
  )
}

# A small algebraic residual only shows that the sampled linear system was
# solved. Interpreting that solution as the limit of Appendix-E iteration also
# requires a contractive empirical operator.
ciftest_direct_converged <- function(
    residual,
    tolerance,
    spectral.radius
) {
  is.numeric(residual) && length(residual) == 1L && is.finite(residual) &&
    is.numeric(tolerance) && length(tolerance) == 1L &&
    is.finite(tolerance) && tolerance >= 0 && residual <= tolerance &&
    is.numeric(spectral.radius) && length(spectral.radius) == 1L &&
    is.finite(spectral.radius) && spectral.radius < 1
}

# Evaluate L(t, X) = F2(t | X) / S(t | X) from one coherent working-AJ
# distribution. Positivity is required only where the caller's integrand
# actually uses L; this avoids rejecting irrelevant tail time points.
ciftest_appendix_e_l_process <- function(
    working,
    time,
    exposure,
    active = rep.int(TRUE, length(time)),
    prob.bound = 1e-7
) {
  time <- as.numeric(time)
  active <- as.logical(active)
  if (length(active) != length(time) || anyNA(active)) {
    stop("Appendix-E activity indicators are invalid.", call. = FALSE)
  }
  if (!length(time)) {
    return(list(
      value = numeric(),
      prediction = matrix(numeric(), 0L, 3L),
      minimum.active.survival = NA_real_,
      active.count = 0L
    ))
  }
  prediction <- predict_working_aj(
    working,
    time,
    rep.int(as.character(exposure), length(time)),
    rep.int("pooled", length(time)),
    side = "left"
  )
  prediction_raw <- predict_working_aj(
    working,
    time,
    rep.int(as.character(exposure), length(time)),
    rep.int("pooled", length(time)),
    side = "left",
    use = "raw"
  )
  required <- active &
    prediction[, "cif2"] > sqrt(.Machine$double.eps)
  if (any(required & prediction_raw[, "survival"] <= prob.bound)) {
    stop("Working-model positivity is violated on the active Appendix-E support.",
         call. = FALSE)
  }
  value <- numeric(length(time))
  value[required] <- prediction[required, "cif2"] /
    prediction[required, "survival"]
  list(
    value = value,
    prediction = prediction,
    minimum.active.survival = if (any(required)) {
      min(prediction[required, "survival"])
    } else {
      NA_real_
    },
    active.count = sum(required)
  )
}

# Evaluate the empirical Appendix-E correction J(v) in one working-AJ cell.
# The implementation is a discrete-time transcription of Supplement E.3,
# with the Fine-Gray baseline subdistribution hazard evaluated under H0.
ciftest_appendix_e_j_cell <- function(
    v,
    base,
    working,
    working.cell,
    exposure,
    prob.bound = 1e-7
) {
  event_time <- as.numeric(base$event.time)
  event_count <- as.numeric(base$event.count)
  risk_total <- as.numeric(base$risk.total)
  d_lambda1 <- event_count / risk_total
  if (any(!is.finite(d_lambda1)) || any(d_lambda1 < 0)) {
    stop("The null Fine-Gray baseline hazard is not finite.", call. = FALSE)
  }
  cumulative_lambda1 <- cumsum(d_lambda1)
  lambda1_left <- c(0, utils::head(cumulative_lambda1, -1L))
  d_f1 <- exp(-lambda1_left) * (-expm1(-d_lambda1))

  censoring <- base$censoring
  if (length(censoring$strata.levels) != 1L) {
    stop("Appendix-E iteration requires a marginal censoring distribution.",
         call. = FALSE)
  }
  censor_level <- censoring$strata.levels[[1L]]
  g_event <- predict_censoring_km(
    censoring,
    event_time,
    rep.int(censor_level, length(event_time)),
    side = "left"
  )
  g_event_raw <- predict_censoring_km(
    censoring,
    event_time,
    rep.int(censor_level, length(event_time)),
    side = "left",
    use = "raw"
  )
  if (any(g_event_raw <= prob.bound)) {
    stop("Censoring positivity is violated in Appendix-E iteration.",
         call. = FALSE)
  }

  censoring <- base$censoring
  censor_hazard_all <- censoring$hazard.table
  # J(v) is evaluated only at observed cause-1 event times. With the strict
  # left-limit convention in b_at(), censoring jumps at or after the final
  # cause-1 event cannot enter the fixed-point map.
  censor_hazard <- censor_hazard_all[
    censor_hazard_all$time < max(event_time),
    ,
    drop = FALSE
  ]
  ignored_tail <- nrow(censor_hazard_all) - nrow(censor_hazard)
  if (!nrow(censor_hazard) ||
      !any(censor_hazard$time < max(event_time))) {
    answer <- matrix(0, nrow(v), ncol(v), dimnames = dimnames(v))
    attr(answer, "appendix_e.diagnostics") <- list(
      ignored.tail.censoring.times = ignored_tail,
      minimum.active.working.survival = NA_real_,
      skipped.zero.future.evaluations = 0L
    )
    return(answer)
  }

  cell_table <- working$table[
    working$table$cell == working.cell,
    ,
    drop = FALSE
  ]
  d_f2_all <- if ("d.cif2" %in% names(cell_table)) {
    pmax(0, cell_table$d.cif2)
  } else {
    pmax(0, cell_table$cif2.right - cell_table$cif2.left)
  }
  f2_keep <- d_f2_all > sqrt(.Machine$double.eps)
  f2_time <- cell_table$time[f2_keep]
  d_f2 <- d_f2_all[f2_keep]

  censor_time <- as.numeric(censor_hazard$time)
  censor_increment <- as.numeric(censor_hazard$hazard)
  vd_lambda <- v * d_lambda1
  if (length(censor_time)) {
    future_active <- vapply(
      censor_time,
      function(current) {
        future <- event_time > current
        if (!any(future)) return(FALSE)
        any(abs(colSums(vd_lambda[future, , drop = FALSE])) >
              sqrt(.Machine$double.eps))
      },
      logical(1L)
    )
    l_process <- ciftest_appendix_e_l_process(
      working,
      censor_time,
      exposure,
      active = censor_increment > 0 & future_active,
      prob.bound = prob.bound
    )
    censor_integrand <- l_process$value * censor_increment
    skipped_zero_future <- sum(
      censor_increment > 0 & !future_active &
        l_process$prediction[, "cif2"] > sqrt(.Machine$double.eps)
    )
  } else {
    l_process <- list(minimum.active.survival = NA_real_)
    censor_integrand <- numeric()
    skipped_zero_future <- 0L
  }

  b_at <- function(query) {
    if (!length(censor_time)) return(numeric(length(query)))
    vapply(
      query,
      function(current) sum(censor_integrand[censor_time < current]),
      numeric(1L)
    )
  }

  b_event <- b_at(event_time)
  prefix_b_v <- apply(vd_lambda * b_event, 2L, cumsum)
  if (is.null(dim(prefix_b_v))) prefix_b_v <- matrix(prefix_b_v, ncol = 1L)
  total_v <- colSums(vd_lambda)

  a_at <- function(query) {
    answer <- matrix(0, length(query), ncol(v))
    b_query <- b_at(query)
    for (q in seq_along(query)) {
      through <- which(event_time <= query[q])
      prefix <- if (length(through)) {
        prefix_b_v[max(through), ]
      } else {
        rep.int(0, ncol(v))
      }
      used <- if (length(through)) {
        colSums(vd_lambda[through, , drop = FALSE])
      } else {
        rep.int(0, ncol(v))
      }
      answer[q, ] <- -prefix - b_query[q] * (total_v - used)
    }
    answer
  }

  a_event <- a_at(event_time)
  a_f2 <- if (length(f2_time)) a_at(f2_time) else
    matrix(0, 0L, ncol(v))
  competing_projection <- if (length(f2_time)) {
    colSums(a_f2 * d_f2)
  } else {
    rep.int(0, ncol(v))
  }
  future_f1 <- ciftest_reverse_cumsum(a_event * d_f1)
  v2 <- a_event - exp(cumulative_lambda1) *
    sweep(future_f1, 2L, competing_projection, "+")

  if (length(f2_time)) {
    g_f2 <- predict_censoring_km(
      censoring,
      f2_time,
      rep.int(censor_level, length(f2_time)),
      side = "left"
    )
    g_f2_raw <- predict_censoring_km(
      censoring,
      f2_time,
      rep.int(censor_level, length(f2_time)),
      side = "left",
      use = "raw"
    )
    if (any(g_f2_raw <= prob.bound)) {
      stop("Censoring positivity is violated at a competing-event jump.",
           call. = FALSE)
    }
    inner_f2 <- vapply(
      event_time,
      function(current) {
        eligible <- f2_time <= current
        if (!any(eligible)) return(0)
        g_current <- predict_censoring_km(
          censoring, current, censor_level, side = "left"
        )
        g_current_raw <- predict_censoring_km(
          censoring, current, censor_level, side = "left", use = "raw"
        )
        if (g_current_raw <= prob.bound) {
          stop("Censoring positivity is violated in Appendix-E iteration.",
               call. = FALSE)
        }
        sum((1 / g_current - 1 / g_f2[eligible]) * d_f2[eligible])
      },
      numeric(1L)
    )
    l_v <- colSums(vd_lambda * inner_f2)
  } else {
    l_v <- rep.int(0, ncol(v))
  }

  answer <- v2 - outer(exp(cumulative_lambda1), l_v)
  attr(answer, "appendix_e.diagnostics") <- list(
    ignored.tail.censoring.times = ignored_tail,
    minimum.active.working.survival =
      l_process$minimum.active.survival,
    skipped.zero.future.evaluations = skipped_zero_future
  )
  answer
}

ciftest_fg_risk_weights <- function(
    base,
    t,
    epsilon,
    weights,
    code.event2,
    prob.bound = 1e-7
) {
  event_time <- as.numeric(base$event.time)
  censoring <- base$censoring
  censor_level <- censoring$strata.levels[[1L]]
  result <- matrix(0, length(t), length(event_time))
  g_at_competing <- predict_censoring_km(
    censoring, t, rep.int(censor_level, length(t)), side = "left"
  )
  g_at_competing_raw <- predict_censoring_km(
    censoring, t, rep.int(censor_level, length(t)), side = "left", use = "raw"
  )
  for (j in seq_along(event_time)) {
    current <- event_time[j]
    competing_before <- epsilon == code.event2 & t < current
    subrisk <- t >= current | competing_before
    censor_weight <- rep.int(1, length(t))
    if (any(competing_before)) {
      g_current <- predict_censoring_km(
        censoring,
        rep.int(current, sum(competing_before)),
        rep.int(censor_level, sum(competing_before)),
        side = "left"
      )
      g_current_raw <- predict_censoring_km(
        censoring,
        rep.int(current, sum(competing_before)),
        rep.int(censor_level, sum(competing_before)),
        side = "left",
        use = "raw"
      )
      denominator <- g_at_competing[competing_before]
      denominator_raw <- g_at_competing_raw[competing_before]
      if (any(g_current_raw <= prob.bound) ||
          any(denominator_raw <= prob.bound)) {
        stop("Censoring positivity is violated in an iterated risk set.",
             call. = FALSE)
      }
      censor_weight[competing_before] <- g_current / denominator
    }
    result[, j] <- weights * subrisk * censor_weight
  }
  result
}

ciftest_appendix_e_map <- function(
    v,
    seed,
    base,
    working,
    working_cell,
    exposure_by_cell,
    subject_cell,
    risk_weight,
    prob.bound = 1e-7
) {
  cell_count <- dim(v)[1L]
  event_count <- dim(v)[2L]
  p <- dim(v)[3L]
  j_value <- array(0, dim(v), dimnames = dimnames(v))
  cell_diagnostics <- vector("list", cell_count)
  for (cell_index in seq_len(cell_count)) {
    cell_value <- ciftest_appendix_e_j_cell(
      v = matrix(v[cell_index, , ], nrow = event_count, ncol = p),
      base = base,
      working = working,
      working.cell = working_cell[cell_index],
      exposure = exposure_by_cell[cell_index],
      prob.bound = prob.bound
    )
    cell_diagnostics[[cell_index]] <-
      attr(cell_value, "appendix_e.diagnostics")
    j_value[cell_index, , ] <- cell_value
  }

  e_j <- matrix(0, event_count, p)
  for (event_index in seq_len(event_count)) {
    denominator <- sum(risk_weight[, event_index])
    if (denominator <= prob.bound) {
      stop("An iterated Fine-Gray risk set is empty.", call. = FALSE)
    }
    subject_j <- matrix(
      j_value[subject_cell, event_index, ],
      nrow = length(subject_cell), ncol = p
    )
    e_j[event_index, ] <- colSums(
      subject_j * risk_weight[, event_index]
    ) / denominator
  }

  g_event <- predict_censoring_km(
    base$censoring,
    base$event.time,
    rep.int(base$censoring$strata.levels[[1L]], length(base$event.time)),
    side = "left"
  )
  answer <- seed
  for (cell_index in seq_len(cell_count)) {
    centered_j <-
      matrix(j_value[cell_index, , ], nrow = event_count, ncol = p) - e_j
    answer[cell_index, , ] <- seed[cell_index, , ] - g_event * centered_j
  }
  active_survival <- vapply(
    cell_diagnostics,
    function(x) x$minimum.active.working.survival,
    numeric(1L)
  )
  attr(answer, "appendix_e.diagnostics") <- list(
    ignored.tail.censoring.times = max(vapply(
      cell_diagnostics,
      function(x) x$ignored.tail.censoring.times,
      integer(1L)
    )),
    minimum.active.working.survival = if (any(is.finite(active_survival))) {
      min(active_survival[is.finite(active_survival)])
    } else {
      NA_real_
    },
    skipped.zero.future.evaluations = sum(vapply(
      cell_diagnostics,
      function(x) x$skipped.zero.future.evaluations,
      integer(1L)
    ))
  )
  answer
}

# Map a full-data martingale direction v to its observed-data AIPWCC
# contribution, Scheike et al. (2023), equation (7) and Supplement E.3.
ciftest_iterated_aipwcc_iid <- function(
    v,
    base,
    working,
    t,
    epsilon,
    weights,
    working_cell,
    exposure_by_cell,
    subject_cell,
    code.event1,
    code.event2,
    code.censoring,
    prob.bound = 1e-7,
    return.components = FALSE
) {
  n <- length(t)
  event_time <- as.numeric(base$event.time)
  event_count <- length(event_time)
  p <- dim(v)[3L]
  d_lambda1 <- as.numeric(base$event.count / base$risk.total)
  censoring <- base$censoring
  censor_level <- censoring$strata.levels[[1L]]
  answer <- matrix(0, n, p, dimnames = list(NULL, dimnames(v)[[3L]]))
  if (!is.logical(return.components) || length(return.components) != 1L ||
      is.na(return.components)) {
    stop("`return.components` must be TRUE or FALSE.", call. = FALSE)
  }
  if (return.components) {
    event_answer <- censor_past_answer <- censor_working_answer <-
      matrix(0, n, p, dimnames = dimnames(answer))
    horizon_completion_answer <- matrix(
      0, n, p, dimnames = dimnames(answer)
    )
  }
  minimum_active_survival <- Inf
  skipped_zero_future <- 0L
  censoring_event <- censoring$censoring.event
  if (is.null(censoring_event)) {
    censoring_event <- epsilon == code.censoring
  }
  if (length(censoring_event) != n) {
    stop("The censoring-event indicator has an incompatible length.",
         call. = FALSE)
  }
  horizon_complete <- epsilon == code.censoring & !censoring_event

  for (cell_index in seq_len(dim(v)[1L])) {
    member <- which(subject_cell == cell_index)
    v_cell <- matrix(v[cell_index, , ], nrow = event_count, ncol = p)
    vd_lambda <- v_cell * d_lambda1
    cumulative_v <- apply(vd_lambda, 2L, cumsum)
    if (is.null(dim(cumulative_v))) cumulative_v <- matrix(cumulative_v, ncol = 1L)
    total_v <- colSums(vd_lambda)
    for (i in member) {
      if (epsilon[i] == code.event1) {
        jump <- match(t[i], event_time)
        if (is.na(jump)) {
          stop("A cause-1 event is absent from the iteration grid.",
               call. = FALSE)
        }
        b_full <- v_cell[jump, ] - cumulative_v[jump, ]
      } else if (epsilon[i] == code.event2) {
        b_full <- -total_v
      } else if (horizon_complete[i]) {
        b_full <- -total_v
      } else {
        b_full <- rep.int(0, p)
      }
      terminal_observed <- epsilon[i] != code.censoring || horizon_complete[i]
      if (terminal_observed) {
        g_terminal <- predict_censoring_km(
          censoring, t[i], censor_level, side = "left"
        )
        g_terminal_raw <- predict_censoring_km(
          censoring, t[i], censor_level, side = "left", use = "raw"
        )
        if (g_terminal_raw <= prob.bound) {
          stop("Censoring positivity is violated at an observed event.",
               call. = FALSE)
        }
        increment <- weights[i] * b_full / g_terminal
        answer[i, ] <- answer[i, ] + increment
        if (return.components) {
          event_answer[i, ] <- event_answer[i, ] + increment
          if (horizon_complete[i]) {
            horizon_completion_answer[i, ] <-
              horizon_completion_answer[i, ] + increment
          }
        }
      }
    }
  }

  censor_hazard <- censoring$hazard.table
  if (nrow(censor_hazard)) {
    for (hazard_index in seq_len(nrow(censor_hazard))) {
      current <- censor_hazard$time[hazard_index]
      at_risk <- t >= current
      censored <- t == current & censoring_event
      martingale <- weights * (
        as.numeric(censored) -
          censor_hazard$hazard[hazard_index] * as.numeric(at_risk)
      )
      if (!any(martingale != 0)) next
      g_left_raw <- censor_hazard$survival.left.raw[hazard_index]
      if (g_left_raw <= prob.bound) {
        stop("Censoring positivity is violated in the AIPWCC map.",
             call. = FALSE)
      }
      g_left <- censor_hazard$survival.left.used[hazard_index]
      for (cell_index in seq_len(dim(v)[1L])) {
        member <- which(subject_cell == cell_index)
        v_cell <- matrix(v[cell_index, , ], nrow = event_count, ncol = p)
        vd_lambda <- v_cell * d_lambda1
        through <- event_time <= current
        past <- if (any(through)) {
          colSums(vd_lambda[through, , drop = FALSE])
        } else {
          rep.int(0, p)
        }
        future <- if (any(!through)) {
          colSums(vd_lambda[!through, , drop = FALSE])
        } else {
          rep.int(0, p)
        }
        needs_l <- any(abs(future) > sqrt(.Machine$double.eps))
        l_process <- ciftest_appendix_e_l_process(
          working,
          current,
          exposure_by_cell[cell_index],
          active = needs_l,
          prob.bound = prob.bound
        )
        l_value <- l_process$value
        if (is.finite(l_process$minimum.active.survival)) {
          minimum_active_survival <- min(
            minimum_active_survival,
            l_process$minimum.active.survival
          )
        }
        if (!needs_l && l_process$prediction[, "cif2"] >
            sqrt(.Machine$double.eps)) {
          skipped_zero_future <- skipped_zero_future + 1L
        }
        past_increment <- outer(martingale[member] / g_left, -past)
        working_increment <- outer(
          martingale[member] / g_left, l_value * future
        )
        answer[member, ] <- answer[member, , drop = FALSE] +
          past_increment + working_increment
        if (return.components) {
          censor_past_answer[member, ] <-
            censor_past_answer[member, , drop = FALSE] + past_increment
          censor_working_answer[member, ] <-
            censor_working_answer[member, , drop = FALSE] + working_increment
        }
      }
    }
  }
  attr(answer, "appendix_e.diagnostics") <- list(
    minimum.active.working.survival = if (
      is.finite(minimum_active_survival)
    ) minimum_active_survival else NA_real_,
    skipped.zero.future.evaluations = skipped_zero_future
  )
  if (return.components) {
    attr(answer, "appendix_e.components") <- list(
      event = event_answer,
      censor.past = censor_past_answer,
      censor.working.aj = censor_working_answer,
      censor.total = censor_past_answer + censor_working_answer,
      horizon.completion = horizon_completion_answer
    )
  }
  answer
}

#' Build a finite-iteration time-directed Fine-Gray score
#'
#' `iterations = 0` is intentionally handled by the closed-form caller. This
#' function computes one or more empirical Appendix-E refinements and maps the
#' resulting full-data direction to observed-data score contributions.
#'
#' @keywords internal
ciftest_iteration_setup <- function(
    base,
    working,
    t,
    epsilon,
    x,
    exposure,
    weights = rep.int(1, length(t)),
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    prob.bound = 1e-7,
    prob.truncation = NULL
) {
  if (length(base$censoring$strata.levels) != 1L ||
      any(working$cells$stratum != "pooled")) {
    stop("Iteration currently requires pooled censoring and working models.",
         call. = FALSE)
  }
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- length(t)
  p <- ncol(x)
  event_time <- as.numeric(base$event.time)
  event_count <- length(event_time)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))

  exposure_label <- as.character(exposure)
  working_cell_subject <- ciftest_working_cell(
    exposure_label, rep.int("pooled", n)
  )
  working_cell <- unique(working_cell_subject)
  subject_cell <- match(working_cell_subject, working_cell)
  exposure_by_cell <- vapply(
    working_cell,
    function(cell) exposure_label[match(cell, working_cell_subject)],
    character(1L)
  )
  if (any(!working_cell %in% working$cells$cell)) {
    stop("Iteration contains an unknown working-AJ cell.", call. = FALSE)
  }

  g_event <- predict_censoring_km(
    base$censoring,
    event_time,
    rep.int(base$censoring$strata.levels[[1L]], event_count),
    side = "left"
  )
  seed <- array(
    0,
    dim = c(length(working_cell), event_count, p),
    dimnames = list(working_cell, NULL, score_names)
  )
  for (cell_index in seq_along(working_cell)) {
    member <- which(subject_cell == cell_index)[1L]
    centered <- sweep(base$xbar, 2L, x[member, ],
                      function(mean, value) value - mean)
    seed[cell_index, , ] <- centered * (g_event * base$fh.weight)
  }
  risk_weight <- ciftest_fg_risk_weights(
    base, t, epsilon, weights, code.event2, prob.bound
  )

  list(
    n = n,
    p = p,
    event.time = event_time,
    event.count = event_count,
    score.names = score_names,
    working.cell = working_cell,
    exposure.by.cell = exposure_by_cell,
    subject.cell = subject_cell,
    g.event = g_event,
    seed = seed,
    risk.weight = risk_weight
  )
}

# Construct the affine-linear Appendix-E operator T(v) = seed + K v.
# The factorization is kept with the operator so several FH seeds can reuse it.
ciftest_appendix_e_operator <- function(
    setup,
    base,
    working,
    prob.bound = 1e-7,
    max.states = 600L
) {
  state_dimension <- length(setup$working.cell) * setup$event.count
  if (!is.numeric(max.states) || length(max.states) != 1L ||
      !is.finite(max.states) || max.states < 1L) {
    stop("`max.states` must be one positive finite number.", call. = FALSE)
  }
  if (state_dimension > max.states) {
    stop(
      "The direct fixed-point state dimension is ", state_dimension,
      ", exceeding the safety limit of ", as.integer(max.states), ".",
      call. = FALSE
    )
  }
  started <- proc.time()[["elapsed"]]
  basis <- array(
    diag(state_dimension),
    dim = c(length(setup$working.cell), setup$event.count, state_dimension)
  )
  zero_seed <- array(0, dim = dim(basis))
  mapped <- ciftest_appendix_e_map(
    basis, zero_seed, base, working,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    setup$risk.weight, prob.bound
  )
  map_diagnostics <- attr(mapped, "appendix_e.diagnostics")
  K <- matrix(mapped, nrow = state_dimension, ncol = state_dimension)
  system <- diag(state_dimension) - K
  factorization <- qr(system, LAPACK = FALSE)
  if (factorization$rank < state_dimension) {
    stop("The direct fixed-point system is rank deficient.", call. = FALSE)
  }
  reciprocal_condition <- tryCatch(
    rcond(system),
    error = function(condition) NA_real_
  )
  eigenvalues <- tryCatch(
    eigen(K, only.values = TRUE)$values,
    error = function(condition) complex(real = NA_real_, imaginary = NA_real_)
  )
  spectral_radius <- if (all(is.finite(eigenvalues))) {
    max(Mod(eigenvalues))
  } else {
    NA_real_
  }
  list(
    K = K,
    system = system,
    qr = factorization,
    state.dimension = state_dimension,
    spectral.radius = spectral_radius,
    reciprocal.condition = reciprocal_condition,
    condition.number = if (is.finite(reciprocal_condition) &&
      reciprocal_condition > 0) 1 / reciprocal_condition else Inf,
    build.seconds = proc.time()[["elapsed"]] - started,
    diagnostics = map_diagnostics
  )
}

# Solve the Appendix-E fixed point directly. This is intentionally internal:
# it is a reference path for simulations, not part of the scalar public UI.
build_direct_fixed_point_score_iid <- function(
    base,
    working,
    t,
    epsilon,
    x,
    exposure,
    weights = rep.int(1, length(t)),
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    prob.bound = 1e-7,
    tolerance = 1e-8,
    max.states = 600L,
    setup = NULL,
    operator = NULL
) {
  if (is.null(setup)) {
    setup <- ciftest_iteration_setup(
      base, working, t, epsilon, x, exposure, weights,
      code.event1, code.event2, code.censoring, prob.bound
    )
  }
  if (is.null(operator)) {
    operator <- ciftest_appendix_e_operator(
      setup, base, working, prob.bound, max.states
    )
  }
  state_dimension <- operator$state.dimension
  rhs <- matrix(setup$seed, nrow = state_dimension, ncol = setup$p)
  solve_started <- proc.time()[["elapsed"]]
  value_matrix <- qr.coef(operator$qr, rhs)
  solve_seconds <- proc.time()[["elapsed"]] - solve_started
  if (any(!is.finite(value_matrix))) {
    stop("The direct fixed-point solution is non-finite.", call. = FALSE)
  }
  value <- array(value_matrix, dim = dim(setup$seed),
                 dimnames = dimnames(setup$seed))
  one_more <- ciftest_appendix_e_map(
    value, setup$seed, base, working,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    setup$risk.weight, prob.bound
  )
  map_diagnostics <- attr(one_more, "appendix_e.diagnostics")
  residual <- max(abs(one_more - value)) / max(1, max(abs(value)))

  score_iid <- ciftest_iterated_aipwcc_iid(
    value, base, working, t, epsilon, weights,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    code.event1, code.event2, code.censoring, prob.bound
  )
  aipwcc_diagnostics <- attr(score_iid, "appendix_e.diagnostics")
  attr(score_iid, "appendix_e.diagnostics") <- NULL
  colnames(score_iid) <- setup$score.names
  seed_score_iid <- ciftest_iterated_aipwcc_iid(
    setup$seed, base, working, t, epsilon, weights,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    code.event1, code.event2, code.censoring, prob.bound
  )
  seed_aipwcc_diagnostics <- attr(
    seed_score_iid, "appendix_e.diagnostics"
  )
  attr(seed_score_iid, "appendix_e.diagnostics") <- NULL
  colnames(seed_score_iid) <- setup$score.names
  active_survival <- c(
    map_diagnostics$minimum.active.working.survival,
    aipwcc_diagnostics$minimum.active.working.survival,
    seed_aipwcc_diagnostics$minimum.active.working.survival
  )
  list(
    score = stats::setNames(colSums(score_iid), setup$score.names),
    score.iid = score_iid,
    score.iid.seed = seed_score_iid,
    iterations = NA_integer_,
    converged = ciftest_direct_converged(
      residual, tolerance, operator$spectral.radius
    ),
    fixed.point.residual = residual,
    last.increment = NA_real_,
    contraction.ratio = NA_real_,
    diagnostics = list(
      ignored.tail.censoring.times =
        map_diagnostics$ignored.tail.censoring.times,
      minimum.active.working.survival = if (any(is.finite(active_survival))) {
        min(active_survival[is.finite(active_survival)])
      } else {
        NA_real_
      },
      skipped.zero.future.evaluations =
        map_diagnostics$skipped.zero.future.evaluations +
          aipwcc_diagnostics$skipped.zero.future.evaluations,
      seed.skipped.zero.future.evaluations =
        seed_aipwcc_diagnostics$skipped.zero.future.evaluations,
      solver = "direct",
      state.dimension = state_dimension,
      operator.spectral.radius = operator$spectral.radius,
      operator.contractive = is.finite(operator$spectral.radius) &&
        operator$spectral.radius < 1,
      system.reciprocal.condition = operator$reciprocal.condition,
      system.condition.number = operator$condition.number,
      operator.build.seconds = operator$build.seconds,
      solve.seconds = solve_seconds
    ),
    history = NULL,
    fixed.point = list(
      seed = setup$seed,
      value = value,
      event.time = setup$event.time,
      working.cell = setup$working.cell
    )
  )
}

# Map the unrefined Appendix-E seed through the observed-data AIPWCC map.
# This batch-only diagnostic isolates the empirical identity
# AIPWCC(seed) = closed-form augmentation from the subsequent fixed-point
# refinements. It is intentionally not exposed by the scalar ciftest() UI.
build_seed_aipwcc_score_iid <- function(
    base,
    working,
    t,
    epsilon,
    x,
    exposure,
    weights = rep.int(1, length(t)),
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    prob.bound = 1e-7,
    prob.truncation = NULL
) {
  setup <- ciftest_iteration_setup(
    base, working, t, epsilon, x, exposure, weights,
    code.event1, code.event2, code.censoring, prob.bound
  )
  value <- setup$seed
  one_more <- ciftest_appendix_e_map(
    value, setup$seed, base, working,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    setup$risk.weight, prob.bound
  )
  map_diagnostics <- attr(one_more, "appendix_e.diagnostics")
  residual <- max(abs(one_more - value)) / max(1, max(abs(value)))

  score_iid <- ciftest_iterated_aipwcc_iid(
    value, base, working, t, epsilon, weights,
    setup$working.cell, setup$exposure.by.cell, setup$subject.cell,
    code.event1, code.event2, code.censoring, prob.bound,
    return.components = TRUE
  )
  aipwcc_diagnostics <- attr(score_iid, "appendix_e.diagnostics")
  aipwcc_components <- attr(score_iid, "appendix_e.components")
  attr(score_iid, "appendix_e.diagnostics") <- NULL
  attr(score_iid, "appendix_e.components") <- NULL
  colnames(score_iid) <- setup$score.names
  active_survival <- c(
    map_diagnostics$minimum.active.working.survival,
    aipwcc_diagnostics$minimum.active.working.survival
  )
  list(
    score = stats::setNames(colSums(score_iid), setup$score.names),
    score.iid = score_iid,
    score.iid.seed = score_iid,
    components = aipwcc_components,
    iterations = 0L,
    converged = NA,
    fixed.point.residual = residual,
    last.increment = 0,
    contraction.ratio = NA_real_,
    diagnostics = list(
      ignored.tail.censoring.times =
        map_diagnostics$ignored.tail.censoring.times,
      minimum.active.working.survival = if (any(is.finite(active_survival))) {
        min(active_survival[is.finite(active_survival)])
      } else {
        NA_real_
      },
      skipped.zero.future.evaluations =
        map_diagnostics$skipped.zero.future.evaluations +
          aipwcc_diagnostics$skipped.zero.future.evaluations,
      solver = "seed-map",
      state.dimension = length(value),
      operator.spectral.radius = NA_real_,
      system.reciprocal.condition = NA_real_,
      system.condition.number = NA_real_,
      operator.build.seconds = 0,
      solve.seconds = 0
    ),
    history = NULL,
    fixed.point = list(
      seed = setup$seed,
      value = value,
      event.time = setup$event.time,
      working.cell = setup$working.cell
    )
  )
}

#' Build a finite-iteration time-directed Fine-Gray score
#'
#' `iterations = 0` is intentionally handled by the closed-form caller. This
#' function computes one or more empirical Appendix-E refinements and maps the
#' resulting full-data direction to observed-data score contributions.
#'
#' @keywords internal
build_iterated_score_iid <- function(
    base,
    working,
    t,
    epsilon,
    x,
    exposure,
    weights = rep.int(1, length(t)),
    iterations = 1L,
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    prob.bound = 1e-7,
    tolerance = NULL
) {
  if (!is.numeric(iterations) || is.logical(iterations) ||
      length(iterations) != 1L || !is.finite(iterations) ||
      iterations < 1 || iterations != floor(iterations)) {
    stop("`iterations` must be one positive integer.", call. = FALSE)
  }
  if (!is.null(tolerance) &&
      (!is.numeric(tolerance) || length(tolerance) != 1L ||
       !is.finite(tolerance) || tolerance <= 0)) {
    stop("`tolerance` must be NULL or one positive finite number.",
         call. = FALSE)
  }
  setup <- ciftest_iteration_setup(
    base, working, t, epsilon, x, exposure, weights,
    code.event1, code.event2, code.censoring, prob.bound
  )
  event_time <- setup$event.time
  event_count <- setup$event.count
  p <- setup$p
  score_names <- setup$score.names
  working_cell <- setup$working.cell
  exposure_by_cell <- setup$exposure.by.cell
  subject_cell <- setup$subject.cell
  seed <- setup$seed
  risk_weight <- setup$risk.weight

  value <- seed
  increments <- rep.int(NA_real_, iterations)
  completed <- 0L
  for (iteration_index in seq_len(iterations)) {
    updated <- ciftest_appendix_e_map(
      value, seed, base, working, working_cell, exposure_by_cell,
      subject_cell, risk_weight, prob.bound
    )
    scale <- max(1, max(abs(value)))
    increments[iteration_index] <- max(abs(updated - value)) / scale
    value <- updated
    completed <- iteration_index
    if (!is.null(tolerance) && increments[iteration_index] <= tolerance) {
      break
    }
  }
  increments <- increments[seq_len(completed)]
  one_more <- ciftest_appendix_e_map(
    value, seed, base, working, working_cell, exposure_by_cell,
    subject_cell, risk_weight, prob.bound
  )
  map_diagnostics <- attr(one_more, "appendix_e.diagnostics")
  residual <- max(abs(one_more - value)) / max(1, max(abs(value)))
  ratios <- c(NA_real_, increments[-1L] / increments[-length(increments)])
  ratios[!is.finite(ratios)] <- NA_real_

  score_iid <- ciftest_iterated_aipwcc_iid(
    value, base, working, t, epsilon, weights,
    working_cell, exposure_by_cell, subject_cell,
    code.event1, code.event2, code.censoring, prob.bound
  )
  aipwcc_diagnostics <- attr(score_iid, "appendix_e.diagnostics")
  attr(score_iid, "appendix_e.diagnostics") <- NULL
  colnames(score_iid) <- score_names
  seed_score_iid <- ciftest_iterated_aipwcc_iid(
    seed, base, working, t, epsilon, weights,
    working_cell, exposure_by_cell, subject_cell,
    code.event1, code.event2, code.censoring, prob.bound
  )
  seed_aipwcc_diagnostics <- attr(
    seed_score_iid, "appendix_e.diagnostics"
  )
  attr(seed_score_iid, "appendix_e.diagnostics") <- NULL
  colnames(seed_score_iid) <- score_names
  active_survival <- c(
    map_diagnostics$minimum.active.working.survival,
    aipwcc_diagnostics$minimum.active.working.survival,
    seed_aipwcc_diagnostics$minimum.active.working.survival
  )
  list(
    score = stats::setNames(colSums(score_iid), score_names),
    score.iid = score_iid,
    score.iid.seed = seed_score_iid,
    iterations = as.integer(completed),
    converged = if (is.null(tolerance)) NA else residual <= tolerance,
    fixed.point.residual = residual,
    last.increment = increments[[length(increments)]],
    contraction.ratio = if (length(increments) > 1L) {
      ratios[[length(ratios)]]
    } else {
      NA_real_
    },
    diagnostics = list(
      ignored.tail.censoring.times =
        map_diagnostics$ignored.tail.censoring.times,
      minimum.active.working.survival = if (any(is.finite(active_survival))) {
        min(active_survival[is.finite(active_survival)])
      } else {
        NA_real_
      },
      skipped.zero.future.evaluations =
        map_diagnostics$skipped.zero.future.evaluations +
          aipwcc_diagnostics$skipped.zero.future.evaluations,
      seed.skipped.zero.future.evaluations =
        seed_aipwcc_diagnostics$skipped.zero.future.evaluations
    ),
    history = data.frame(
      iteration = seq_len(completed),
      scaled.increment = increments,
      contraction.ratio = ratios
    ),
    fixed.point = list(
      seed = seed,
      value = value,
      event.time = event_time,
      working.cell = working_cell
    )
  )
}

# Batch-only multiweight Fine-Gray engine. Risk sets, censoring derivatives,
# and censoring martingales are shared across all fixed FH weight processes.
build_fg_score_iid_multi <- function(
    t,
    epsilon,
    x,
    rho,
    gamma,
    code.event1 = 1L,
    code.event2 = 2L,
    code.censoring = 0L,
    strata = rep.int("pooled", length(t)),
    weights = rep.int(1, length(t)),
    censoring = NULL,
    prob.bound = 1e-7
) {
  if (!ciftest_use_cpp_kernel(
    "_cifmodeling_ciftest_fg_iid_multi_kernel_cpp"
  )) {
    stop("The multiweight Fine-Gray Rcpp kernel is unavailable.", call. = FALSE)
  }
  rho <- as.numeric(rho)
  gamma <- as.numeric(gamma)
  if (!length(rho) || length(rho) != length(gamma) ||
      any(!is.finite(rho)) || any(rho < 0) ||
      any(!is.finite(gamma)) || any(gamma < 0)) {
    stop("Multiweight rho and gamma must be finite non-negative vectors.",
         call. = FALSE)
  }
  n <- length(t)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  strata_factor <- droplevels(factor(strata))
  strata_label <- as.character(strata_factor)
  event_times <- sort(unique(t[epsilon == code.event1 & weights > 0]))
  if (!length(event_times)) {
    stop("At least one event of interest is required.", call. = FALSE)
  }
  if (is.null(censoring)) {
    censoring <- estimate_censoring_km(
      t, epsilon, code.censoring, strata_factor, weights, prob.bound
    )
  }
  censoring_event <- ciftest_censoring_event_indicator(
    censoring, epsilon, code.censoring
  )
  kernel_epsilon <- ciftest_kernel_status(
    epsilon, censoring_event, code.censoring
  )
  fh_weight <- vapply(
    seq_along(rho),
    function(index) ciftest_null_fh_weight(
      t, epsilon, code.event1, code.event2, weights,
      rho[index], gamma[index]
    ),
    numeric(length(event_times))
  )
  fh_weight <- matrix(fh_weight, nrow = length(event_times))
  censor_hazard <- censoring$hazard.table
  g_at_competing <- predict_censoring_km(
    censoring, t, strata_label, side = "left"
  )
  g_event_stratum <- vapply(
    levels(strata_factor),
    function(level) predict_censoring_km(
      censoring,
      event_times,
      rep.int(level, length(event_times)),
      side = "left"
    ),
    numeric(length(event_times))
  )
  g_event_stratum <- matrix(
    g_event_stratum,
    nrow = length(event_times),
    ncol = nlevels(strata_factor)
  )
  kernel_run <- ciftest_run_fg_cpp_kernel(list(
    t = as.numeric(t),
    epsilon = kernel_epsilon,
    x = x,
    weights = as.numeric(weights),
    strata_id = as.integer(strata_factor),
    code_event1 = as.integer(code.event1),
    code_event2 = as.integer(code.event2),
    code_censoring = as.integer(code.censoring),
    event_times = event_times,
    fh_weight = fh_weight,
    g_at_competing = g_at_competing,
    g_event_stratum = g_event_stratum,
    hazard_time = as.numeric(censor_hazard$time),
    hazard_stratum = as.integer(match(
      as.character(censor_hazard$stratum), levels(strata_factor)
    )),
    hazard = as.numeric(censor_hazard$hazard),
    hazard_n_risk = as.numeric(censor_hazard$n.risk),
    prob_bound = prob.bound
  ), multi = TRUE)
  kernel <- kernel_run$value
  p <- ncol(x)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))
  lapply(seq_along(rho), function(index) {
    base_iid <- matrix(
      kernel$score.iid.base[, , index], nrow = n, ncol = p,
      dimnames = list(NULL, score_names)
    )
    censor_iid <- matrix(
      kernel$score.iid.censor[, , index], nrow = n, ncol = p,
      dimnames = list(NULL, score_names)
    )
    derivative <- matrix(
      kernel$censor.derivative[, , index],
      nrow = nrow(censor_hazard), ncol = p,
      dimnames = list(NULL, score_names)
    )
    score <- stats::setNames(
      as.numeric(kernel$score[, index]), score_names
    )
    event_score <- stats::setNames(
      as.numeric(kernel$event.score[, index]), score_names
    )
    list(
      score = score,
      event.score = event_score,
      score.iid = base_iid + censor_iid,
      score.iid.base = base_iid,
      score.iid.censor = censor_iid,
      censoring = censoring,
      event.time = event_times,
      fh.weight = unname(fh_weight[, index]),
      xbar = structure(kernel$xbar, dimnames = list(NULL, score_names)),
      risk.total = as.numeric(kernel$risk.total),
      event.count = as.numeric(kernel$event.count),
      censor.derivative = derivative,
      diagnostics = list(
        score.decomposition.error = max(abs(score - event_score)),
        censor.centering.error = max(abs(colSums(censor_iid))),
        minimum.censoring.survival = censoring$minimum.survival,
        positivity.warning = censoring$positivity.warning,
        engine = if (identical(kernel_run$engine, "legacy")) {
          "Rcpp-multi"
        } else {
          paste0("Rcpp-", kernel_run$engine, "-multi")
        }
      )
    )
  })
}

# Batch-only multiweight closed-form augmentation engine.
build_closed_form_augmentation_multi <- function(
    bases,
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
  if (!length(bases) || !ciftest_use_cpp_kernel(
    "_cifmodeling_ciftest_augmentation_iid_multi_kernel_cpp"
  )) {
    stop("The multiweight augmentation Rcpp kernel is unavailable.",
         call. = FALSE)
  }
  n <- length(t)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  p <- ncol(x)
  weight_count <- length(bases)
  event_time <- as.numeric(bases[[1L]]$event.time)
  score_names <- colnames(x)
  if (is.null(score_names)) score_names <- paste0("score", seq_len(p))
  exposure_label <- as.character(exposure)
  censor_label <- as.character(strata.censor)
  competing_label <- as.character(strata.competing.risk)
  working_cell <- ciftest_working_cell(exposure_label, competing_label)
  augmentation_cell <- ciftest_augmentation_cell(
    exposure_label, competing_label, censor_label
  )
  augmentation_cells <- unique(data.frame(
    cell = augmentation_cell,
    working.cell = working_cell,
    exposure = exposure_label,
    competing.risk.stratum = competing_label,
    censor.stratum = censor_label,
    stringsAsFactors = FALSE
  ))
  cell_count <- nrow(augmentation_cells)
  h_array <- array(
    0, dim = c(cell_count, length(event_time), p, weight_count)
  )
  h_processes <- lapply(seq_len(weight_count), function(index) {
    out <- vector("list", cell_count)
    names(out) <- augmentation_cells$cell
    out
  })
  for (cell_index in seq_len(cell_count)) {
    member <- which(
      augmentation_cell == augmentation_cells$cell[cell_index] & weights > 0
    )
    if (!length(member)) {
      stop("Every augmentation nuisance cell must have positive total weight.",
           call. = FALSE)
    }
    x_cell <- x[member[1L], ]
    tab <- working$table[
      working$table$cell == augmentation_cells$working.cell[cell_index],
      , drop = FALSE
    ]
    d_lambda1 <- numeric(length(event_time))
    matched <- match(event_time, tab$time)
    present <- !is.na(matched)
    d_lambda1[present] <- tab$d.lambda1[matched[present]]
    g_left <- predict_censoring_km(
      bases[[1L]]$censoring,
      event_time,
      rep.int(augmentation_cells$censor.stratum[cell_index],
              length(event_time)),
      side = "left"
    )
    if (any(g_left[d_lambda1 > 0] <= prob.bound)) {
      stop("Censoring positivity is violated while constructing H_a.",
           call. = FALSE)
    }
    for (weight_index in seq_len(weight_count)) {
      increment <- sweep(
        bases[[weight_index]]$xbar, 2L, x_cell,
        function(mean, value) value - mean
      )
      increment <- increment * (
        bases[[weight_index]]$fh.weight * g_left * d_lambda1
      )
      h_value <- ciftest_reverse_cumsum(increment)
      colnames(h_value) <- score_names
      h_processes[[weight_index]][[cell_index]] <- h_value
      h_array[cell_index, , , weight_index] <- h_value
    }
  }
  censor_hazard <- bases[[1L]]$censoring$hazard.table
  working_levels <- working$cells$cell
  working_survival <- working_cif2 <- matrix(
    numeric(nrow(censor_hazard) * length(working_levels)),
    nrow = nrow(censor_hazard), ncol = length(working_levels)
  )
  if (nrow(censor_hazard)) {
    for (working_index in seq_along(working_levels)) {
      predicted <- predict_working_aj(
        working,
        censor_hazard$time,
        rep.int(working$cells$exposure[working_index], nrow(censor_hazard)),
        rep.int(working$cells$stratum[working_index], nrow(censor_hazard)),
        side = "left"
      )
      working_survival[, working_index] <- predicted[, "survival"]
      working_cif2[, working_index] <- predicted[, "cif2"]
    }
  }
  censor_levels <- bases[[1L]]$censoring$strata.levels
  censoring_event <- ciftest_censoring_event_indicator(
    bases[[1L]]$censoring, epsilon, code.censoring
  )
  kernel_epsilon <- ciftest_kernel_status(
    epsilon, censoring_event, code.censoring
  )
  kernel_run <- ciftest_run_augmentation_cpp_kernel(list(
    t = as.numeric(t),
    epsilon = kernel_epsilon,
    weights = as.numeric(weights),
    censor_stratum_id = as.integer(match(censor_label, censor_levels)),
    augmentation_cell_id = as.integer(match(
      augmentation_cell, augmentation_cells$cell
    )),
    working_cell_id = as.integer(match(working_cell, working_levels)),
    code_censoring = as.integer(code.censoring),
    event_times = event_time,
    hazard_time = as.numeric(censor_hazard$time),
    hazard_stratum = as.integer(match(
      as.character(censor_hazard$stratum), censor_levels
    )),
    hazard = as.numeric(censor_hazard$hazard),
    hazard_n_risk = as.numeric(censor_hazard$n.risk),
    hazard_g_left = as.numeric(censor_hazard$survival.left),
    working_survival = working_survival,
    working_cif2 = working_cif2,
    h_process = h_array,
    prob_bound = prob.bound
  ), multi = TRUE)
  kernel <- kernel_run$value
  lapply(seq_len(weight_count), function(index) {
    augment_iid <- matrix(
      kernel$score.iid.augment[, , index], nrow = n, ncol = p,
      dimnames = list(NULL, score_names)
    )
    score <- stats::setNames(colSums(augment_iid), score_names)
    list(
      score = score,
      score.iid.augment = augment_iid,
      working.aj = working,
      h.process = h_processes[[index]],
      augmentation.cells = augmentation_cells,
      diagnostics = list(
        augmentation.centering.error =
          kernel$augmentation.centering.error[index],
        minimum.working.survival = kernel$minimum.working.survival,
        minimum.censoring.survival = kernel$minimum.censoring.survival,
        positivity.warning = working$positivity.warning ||
          bases[[index]]$censoring$positivity.warning,
        engine = if (identical(kernel_run$engine, "legacy")) {
          "Rcpp-multi"
        } else {
          paste0("Rcpp-", kernel_run$engine, "-multi")
        }
      )
    )
  })
}
