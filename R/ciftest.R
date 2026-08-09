#' Score tests for survival and cumulative incidence curves
#'
#' @description
#' Provides a formula-based interface for log-rank, Fleming-Harrington, Gray,
#' and augmented Fine-Gray score tests. Standard and Fleming-Harrington-type
#' Gray tests are selected with `augmentation = FALSE`. With
#' `augmentation = TRUE`, working Aalen-Johansen distributions are estimated
#' within exposure-by-`strata.competing.risk` cells, censoring distributions
#' are estimated within `strata.censor`, and the closed-form augmentation is
#' used.
#'
#' @param formula A two-sided formula with an `Event()` or `Surv()` response
#'   and one grouping variable on the right-hand side.
#' @param data A data frame.
#' @param weights Optional numeric case weights or the name of a weight column.
#'   The standard Gray branch currently accepts integer frequency weights.
#' @param subset.condition Optional logical subset specification accepted by
#'   the package input interface.
#' @param outcome.type One of `"competing-risk"` or `"survival"`. The default
#'   is `"competing-risk"`; use `NULL` for automatic detection.
#' @param code.event1,code.event2,code.censoring Distinct non-negative integer
#'   event codes.
#' @param rho,gamma Non-negative Fleming-Harrington weight parameters.
#' @param augmentation Logical or `NULL`. `NULL` selects the outcome-specific
#'   default (`FALSE` for survival and `TRUE` for competing risks).
#' @param iteration Use fixed-point refinement of the augmented score.
#' @param iter.max Maximum number of fixed-point iterations.
#' @param iter.tol Positive convergence tolerance, or `NULL` for the eventual
#'   simulation-calibrated default.
#' @param strata.censor Optional single column name defining the censoring
#'   Kaplan-Meier strata.
#' @param strata.competing.risk Optional single column name defining the
#'   exposure-by-stratum working Aalen-Johansen models used for the competing
#'   risk nuisance distribution.
#' @param tau Optional finite non-negative analysis horizon.
#' @param na.action Missing-data action.
#' @param ... Reserved for future extensions.
#'
#' @return An object inheriting from `"ciftest"` and `"htest"` for an
#'   implemented test branch. Augmented results include `score.iid.base`,
#'   `score.iid.censor`, and `score.iid.augment` matrices and use their summed
#'   empirical cross-product as `vcov.score`. Standard Gray results retain the
#'   classical Gray covariance. If the optional standard-Gray score-IID
#'   diagnostic cannot be constructed, the Gray test is still returned; its
#'   score-IID matrices contain `NA` and the reason is recorded in
#'   `diagnostics$score.iid.error`.
#' @export
ciftest <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    outcome.type = "competing-risk",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    rho = 0,
    gamma = 0,
    augmentation = NULL,
    iteration = FALSE,
    iter.max = 50L,
    iter.tol = NULL,
    strata.censor = NULL,
    strata.competing.risk = NULL,
    tau = NULL,
    na.action = stats::na.omit,
    ...
) {
  call <- match.call()
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    displayed <- if (is.null(dot_names)) {
      rep.int("<unnamed>", length(dots))
    } else {
      ifelse(nzchar(dot_names), dot_names, "<unnamed>")
    }
    stop(
      "Unused argument", if (length(dots) == 1L) "" else "s", ": ",
      paste(displayed, collapse = ", "),
      call. = FALSE
    )
  }
  weights.expr <- substitute(weights)
  weights.resolved <- util_resolve_weights(
    weights.expr = weights.expr,
    data = data,
    envir = parent.frame(),
    missing = missing(weights)
  )

  input <- ciftest_prepare(
    formula = formula,
    data = data,
    weights = weights.resolved,
    subset.condition = subset.condition,
    outcome.type = outcome.type,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring,
    rho = rho,
    gamma = gamma,
    augmentation = augmentation,
    iteration = iteration,
    iter.max = iter.max,
    iter.tol = iter.tol,
    strata.censor = strata.censor,
    strata.competing.risk = strata.competing.risk,
    tau = tau,
    na.action = na.action
  )
  score_parts <- NULL
  augmentation_parts <- NULL
  score_iid_error <- NULL

  if (identical(input$outcome.type, "competing-risk")) {
    if (isTRUE(input$iteration)) {
      stop(
        "Phase 3 fixed-point iteration is not implemented yet; ",
        "use `iteration = FALSE` for closed-form augmentation.",
        call. = FALSE
      )
    }
    exposure_design <- reg_read_exposure_design(
      input$data,
      exposure = input$exposure
    )

    if (isTRUE(input$augmentation)) {
      score_parts <- build_fg_score_iid(
        t = input$t,
        epsilon = input$epsilon,
        x = exposure_design$x_a,
        code.event1 = input$code.event1,
        code.event2 = input$code.event2,
        code.censoring = input$code.censoring,
        strata = input$strata.censor.info$values,
        weights = input$weights,
        rho = input$rho,
        gamma = input$gamma
      )
      working_aj <- estimate_working_aj(
        t = input$t,
        epsilon = input$epsilon,
        exposure = input$data[[input$exposure]],
        strata = input$strata.competing.risk.info$values,
        weights = input$weights,
        code.event1 = input$code.event1,
        code.event2 = input$code.event2,
        code.censoring = input$code.censoring
      )
      augmentation_parts <- build_closed_form_augmentation(
        base = score_parts,
        working = working_aj,
        t = input$t,
        epsilon = input$epsilon,
        x = exposure_design$x_a,
        exposure = input$data[[input$exposure]],
        strata.censor = input$strata.censor.info$values,
        strata.competing.risk = input$strata.competing.risk.info$values,
        weights = input$weights,
        code.event1 = input$code.event1,
        code.event2 = input$code.event2,
        code.censoring = input$code.censoring
      )
      total_iid <- score_parts$score.iid +
        augmentation_parts$score.iid.augment
      total_score <- colSums(total_iid)
      comp <- list(
        score = stats::setNames(total_score, colnames(exposure_design$x_a)),
        var = crossprod(total_iid),
        exposure.labels = exposure_design$exposure.labels
      )
      method <- "Closed-form augmented Fine-Gray score test"
      variance_method <- "score-iid"
    } else {
      if (!is.null(input$strata.censor.info$name) ||
          !is.null(input$strata.competing.risk.info$name)) {
        stop(
          "Nuisance-model strata are not available for the standard Gray test; ",
          "use `augmentation = TRUE`.",
          call. = FALSE
        )
      }
      comp <- calculate_gray(
        t = input$t,
        epsilon = as.integer(input$epsilon == input$code.event1) +
          2L * as.integer(input$epsilon == input$code.event2),
        exposure = input$exposure,
        weights = input$weights,
        strata = rep.int(1L, length(input$t)),
        data = input$data,
        rho = input$rho,
        gamma = input$gamma
      )
      gray_censoring_strata <- factor(
        input$data[[input$exposure]],
        levels = comp$exposure.labels
      )
      score_attempt <- tryCatch(
        build_fg_score_iid(
          t = input$t,
          epsilon = input$epsilon,
          x = exposure_design$x_a,
          code.event1 = input$code.event1,
          code.event2 = input$code.event2,
          code.censoring = input$code.censoring,
          strata = gray_censoring_strata,
          weights = input$weights,
          rho = input$rho,
          gamma = input$gamma,
          fh.weight = comp$fh.weight.process[[1L]]$weight
        ),
        error = identity
      )
      if (inherits(score_attempt, "error")) {
        score_iid_error <- conditionMessage(score_attempt)
      } else {
        decomposition_tolerance <- 1e-8 * max(1, max(abs(comp$score)))
        decomposition_error <- max(abs(score_attempt$score - comp$score))
        if (decomposition_error > decomposition_tolerance) {
          score_iid_error <- paste0(
            "Internal Gray score decomposition check failed: error = ",
            format(decomposition_error, digits = 6),
            ", tolerance = ", format(decomposition_tolerance, digits = 6),
            "."
          )
        } else {
          score_parts <- score_attempt
        }
      }
      method <- if (input$rho == 0 && input$gamma == 0) {
        "Gray's test"
      } else {
        "Fleming-Harrington-type weighted Gray test"
      }
      variance_method <- "gray"
    }
  } else {
    comp <- calculate_log_rank(
      t = input$t,
      epsilon = as.integer(input$epsilon == input$code.event1),
      exposure = input$exposure,
      weights = input$weights,
      strata = rep.int(1L, length(input$t)),
      data = input$data,
      rho = input$rho,
      gamma = input$gamma
    )
    method <- if (input$rho == 0 && input$gamma == 0) {
      "Log-rank score test"
    } else {
      "Fleming-Harrington weighted log-rank test"
    }
    variance_method <- "hypergeometric"
  }
  test_stat <- ciftest_quadratic_form(comp$score, comp$var)

  if (is.null(score_parts)) {
    score_iid <- matrix(
      NA_real_,
      nrow = length(input$t),
      ncol = length(comp$score),
      dimnames = list(NULL, names(comp$score))
    )
    score_iid_base <- score_iid_censor <- score_iid
    score_iid_augment <- score_iid_iterated <- score_iid
  } else {
    score_iid_augment <- if (is.null(augmentation_parts)) {
      matrix(
        0,
        nrow = nrow(score_parts$score.iid),
        ncol = ncol(score_parts$score.iid),
        dimnames = dimnames(score_parts$score.iid)
      )
    } else {
      augmentation_parts$score.iid.augment
    }
    score_iid <- score_parts$score.iid + score_iid_augment
    score_iid_base <- score_parts$score.iid.base
    score_iid_censor <- score_parts$score.iid.censor
    score_iid_iterated <- score_iid
  }

  out <- list(
    statistic = stats::setNames(test_stat$statistic, "X-squared"),
    parameter = stats::setNames(test_stat$rank, "df"),
    p.value = stats::pchisq(test_stat$statistic, df = test_stat$rank, lower.tail = FALSE),
    method = method,
    data.name = paste(deparse(formula), collapse = " "),
    call = call,
    outcome.type = input$outcome.type,
    code.event1 = input$code.event1,
    code.event2 = input$code.event2,
    code.censoring = input$code.censoring,
    rho = input$rho,
    gamma = input$gamma,
    augmentation = input$augmentation,
    iteration = input$iteration,
    tau = input$tau,
    variance.method = variance_method,
    score = comp$score,
    vcov.score = comp$var,
    score.iid = score_iid,
    score.iid.base = score_iid_base,
    score.iid.censor = score_iid_censor,
    score.iid.augment = score_iid_augment,
    score.iid.iterated = score_iid_iterated,
    iterations = 0L,
    converged = NA,
    fixed.point.residual = NA_real_,
    last.increment = NA_real_,
    contraction.ratio = NA_real_,
    n = length(input$t),
    n.events = table(factor(input$epsilon,
                            levels = c(input$code.censoring,
                                       input$code.event1,
                                       input$code.event2))),
    strata.censor.info = input$strata.censor.info,
    strata.competing.risk.info = input$strata.competing.risk.info,
    diagnostics = list(
      iter.max = input$iter.max,
      iter.tol = input$iter.tol,
      score.iid.available = !is.null(score_parts),
      score.iid.error = score_iid_error,
      score.iid.variance.role = if (is.null(score_parts)) {
        "not available"
      } else if (!is.null(augmentation_parts)) {
        "empirical score-iid cross-product used as the test variance"
      } else {
        "diagnostic only; Gray covariance remains the test variance"
      },
      censoring.km = if (is.null(score_parts)) NULL else score_parts$censoring,
      working.aj = if (is.null(augmentation_parts)) NULL else
        augmentation_parts$working.aj,
      augmentation.cells = if (is.null(augmentation_parts)) NULL else
        augmentation_parts$augmentation.cells,
      h.process = if (is.null(augmentation_parts)) NULL else
        augmentation_parts$h.process,
      augmentation.centering.error = if (is.null(augmentation_parts)) {
        NA_real_
      } else {
        augmentation_parts$diagnostics$augmentation.centering.error
      },
      score.decomposition.error = if (is.null(score_parts)) {
        NA_real_
      } else {
        score_parts$diagnostics$score.decomposition.error
      },
      analysis.row.index = input$row.index,
      analysis.included = input$included,
      exclusion.reason = input$exclusion.reason
    ),
    data = input$data.original
  )
  class(out) <- c("ciftest", "htest")
  out
}

ciftest_prepare <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    outcome.type = "competing-risk",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    rho = 0,
    gamma = 0,
    augmentation = NULL,
    iteration = FALSE,
    iter.max = 50L,
    iter.tol = NULL,
    strata.censor = NULL,
    strata.competing.risk = NULL,
    tau = NULL,
    na.action = stats::na.omit
) {
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0) {
    stop("`rho` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("`gamma` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.logical(iteration) || length(iteration) != 1L || is.na(iteration)) {
    stop("`iteration` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(iter.max) != 1L || is.na(iter.max) || !is.finite(iter.max) ||
      iter.max < 1 || iter.max != floor(iter.max)) {
    stop("`iter.max` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(iter.tol) &&
      (!is.numeric(iter.tol) || length(iter.tol) != 1L ||
       !is.finite(iter.tol) || iter.tol <= 0)) {
    stop("`iter.tol` must be NULL or one finite positive number.", call. = FALSE)
  }
  if (!is.null(tau) &&
      (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0)) {
    stop("`tau` must be NULL or one finite non-negative number.", call. = FALSE)
  }
  nuisance_strata <- list(
    strata.censor = strata.censor,
    strata.competing.risk = strata.competing.risk
  )
  for (argument in names(nuisance_strata)) {
    value <- nuisance_strata[[argument]]
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L || !nzchar(value))) {
      stop("`", argument, "` must be NULL or one column name.", call. = FALSE)
    }
    if (!is.null(value) && !value %in% names(data)) {
      stop(argument, " = '", value, "' is not found in data.", call. = FALSE)
    }
  }

  Terms0 <- stats::terms(formula, specials = c("strata", "offset", "cluster"), data = data)
  term_labels <- attr(Terms0, "term.labels")
  if (length(term_labels) != 1L || grepl("[:*(]", term_labels, fixed = FALSE)) {
    stop("The right-hand side of `formula` must be one untransformed grouping variable.", call. = FALSE)
  }
  exposure_vars <- all.vars(formula[[3L]])
  if (length(exposure_vars) != 1L || !exposure_vars %in% names(data)) {
    stop("The grouping variable must be one column in `data`.", call. = FALSE)
  }
  exposure <- exposure_vars[[1L]]

  prepared <- cif_prepare_input(
    formula = formula,
    data = data,
    weights = weights,
    other.variables.analyzed = unique(c(strata.censor, strata.competing.risk)),
    subset.condition = subset.condition,
    na.action = na.action,
    outcome.type = outcome.type,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring
  )
  reg_read_exposure_design(prepared$data, exposure = exposure)

  augmentation <- if (is.null(augmentation)) {
    identical(prepared$outcome.type, "competing-risk")
  } else {
    if (!is.logical(augmentation) || length(augmentation) != 1L || is.na(augmentation)) {
      stop("`augmentation` must be NULL, TRUE, or FALSE.", call. = FALSE)
    }
    augmentation
  }
  if (identical(prepared$outcome.type, "survival") && isTRUE(augmentation)) {
    stop("`augmentation = TRUE` is not available for survival outcomes.", call. = FALSE)
  }
  if (isTRUE(iteration) && !isTRUE(augmentation)) {
    stop("`iteration = TRUE` requires `augmentation = TRUE`.", call. = FALSE)
  }
  if (identical(prepared$outcome.type, "survival") && isTRUE(iteration)) {
    stop("`iteration = TRUE` is not available for survival outcomes.", call. = FALSE)
  }
  if (!isTRUE(augmentation) &&
      (!is.null(strata.censor) || !is.null(strata.competing.risk))) {
    stop(
      "`strata.censor` and `strata.competing.risk` require ",
      "`augmentation = TRUE`.",
      call. = FALSE
    )
  }
  if (isTRUE(iteration) &&
      (!is.null(strata.censor) || !is.null(strata.competing.risk))) {
    stop(
      "Nuisance-model strata are not available with iteration in the initial release.",
      call. = FALSE
    )
  }
  if (isTRUE(iteration) && rho + gamma > 0) {
    stop("Weighted iteration is release-gated pending simulation validation.", call. = FALSE)
  }

  t <- prepared$t
  epsilon <- prepared$epsilon
  tau.used <- if (is.null(tau)) max(t) else tau
  if (!is.null(tau)) {
    after_tau <- t > tau
    t[after_tau] <- tau
    epsilon[after_tau] <- prepared$code.censoring
  }

  strata_censor_values <- if (is.null(strata.censor)) {
    factor(rep.int("pooled", length(t)))
  } else {
    factor(prepared$data[[strata.censor]])
  }
  strata_competing_risk_values <- if (is.null(strata.competing.risk)) {
    factor(rep.int("pooled", length(t)))
  } else {
    factor(prepared$data[[strata.competing.risk]])
  }

  list(
    formula = formula,
    data = prepared$data,
    data.original = data,
    row.index = prepared$row.index,
    included = prepared$included,
    exclusion.reason = prepared$exclusion.reason,
    exposure = exposure,
    outcome.type = prepared$outcome.type,
    code.event1 = prepared$code.event1,
    code.event2 = prepared$code.event2,
    code.censoring = prepared$code.censoring,
    t = t,
    epsilon = epsilon,
    weights = prepared$w,
    rho = as.numeric(rho),
    gamma = as.numeric(gamma),
    augmentation = augmentation,
    iteration = iteration,
    iter.max = as.integer(iter.max),
    iter.tol = iter.tol,
    tau = tau.used,
    strata.censor.info = list(
      name = strata.censor,
      values = strata_censor_values
    ),
    strata.competing.risk.info = list(
      name = strata.competing.risk,
      values = strata_competing_risk_values
    )
  )
}

ciftest_quadratic_form <- function(score, variance) {
  variance <- (variance + t(variance)) / 2
  eig <- eigen(variance, symmetric = TRUE)
  tol <- max(1, max(abs(eig$values))) * sqrt(.Machine$double.eps)
  keep <- eig$values > tol
  rank <- sum(keep)
  if (rank < 1L) stop("The score covariance matrix has rank zero.", call. = FALSE)
  z <- crossprod(eig$vectors[, keep, drop = FALSE], as.numeric(score))
  statistic <- sum((z^2) / eig$values[keep])
  list(statistic = as.numeric(statistic), rank = as.integer(rank))
}

#' Compute Gray's K-sample test
#'
#' Implements the predictable-weight estimating equations and covariance from
#' Gray (1988). The pooled weight is
#' `(1 - F_1(t-))^rho * F_1(t-)^gamma`; `gamma = 0` recovers the family used
#' by `cmprsk::cuminc()`.
#'
#' `weights` are frequency weights. The public UI currently supplies one test
#' stratum; the internal `strata` argument is retained for reference checks and
#' a future explicitly stratified test interface.
#'
#' @keywords internal
calculate_gray <- function(
    t,
    epsilon,
    exposure,
    code.exposure.ref = NULL,
    prefix = "a",
    weights = rep.int(1, length(t)),
    strata = rep.int(1L, length(t)),
    data,
    rho = 0,
    gamma = 0,
    prob.bound = 1e-10
) {
  n <- length(t)
  stopifnot(
    length(epsilon) == n,
    length(weights) == n,
    length(strata) == n
  )
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!is.character(exposure) || length(exposure) != 1L ||
      !exposure %in% names(data)) {
    stop("`exposure` must name one column in `data`.", call. = FALSE)
  }
  if (any(!is.finite(t)) || any(t < 0)) {
    stop("Gray's test requires finite non-negative follow-up times.", call. = FALSE)
  }
  if (anyNA(epsilon) || any(!epsilon %in% 0:2)) {
    stop("`epsilon` must use 0 = censoring, 1 = event of interest, 2 = competing event.",
         call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0) ||
      any(abs(weights - round(weights)) > sqrt(.Machine$double.eps))) {
    stop("Standard Gray currently supports non-negative integer frequency weights only.",
         call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < 0 ||
      !is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("`rho` and `gamma` must be finite non-negative numbers.", call. = FALSE)
  }

  exp_info <- reg_read_exposure_design(
    data = data,
    exposure = exposure,
    code.exposure.ref = code.exposure.ref,
    prefix = prefix
  )
  K <- exp_info$exposure.levels
  if (K < 2L) stop("Exposure must have at least two levels.", call. = FALSE)

  group <- factor(data[[exposure]], levels = exp_info$exposure.labels)
  gid <- as.integer(group)
  if (anyNA(gid)) {
    stop("Missing or unknown exposure levels remain after preprocessing.", call. = FALSE)
  }

  q <- K - 1L
  score_base <- numeric(q)
  variance_base <- matrix(0, q, q)
  stratum <- droplevels(factor(strata))
  weight_process <- vector("list", nlevels(stratum))
  names(weight_process) <- levels(stratum)

  for (lev in levels(stratum)) {
    use <- which(stratum == lev & weights > 0)
    if (!length(use)) next
    part <- gray_stratum_components(
      t = t[use],
      epsilon = epsilon[use],
      gid = gid[use],
      weights = as.numeric(weights[use]),
      K = K,
      rho = rho,
      gamma = gamma,
      prob.bound = prob.bound
    )
    score_base <- score_base + part$score
    variance_base <- variance_base + part$var
    weight_process[[lev]] <- data.frame(
      time = part$event.time,
      weight = part$fh.weight
    )
  }

  # Gray's covariance is naturally parameterized by the first K-1 groups.
  # Transform it to the package convention: level 1 is the reference and the
  # returned coordinates correspond to levels 2,...,K.
  transform <- matrix(0, q, q)
  if (q > 1L) {
    transform[cbind(seq_len(q - 1L), 2:q)] <- 1
  }
  transform[q, ] <- -1
  score <- as.vector(transform %*% score_base)
  variance <- transform %*% variance_base %*% t(transform)

  score_names <- colnames(exp_info$x_a)
  names(score) <- score_names
  dimnames(variance) <- list(score_names, score_names)

  list(
    score = score,
    var = (variance + t(variance)) / 2,
    df = q,
    exposure.levels = K,
    exposure.labels = exp_info$exposure.labels,
    ref = exp_info$ref,
    fh.weight.process = weight_process
  )
}

# One-stratum Gray score and covariance in the first K-1 group coordinates.
gray_stratum_components <- function(
    t,
    epsilon,
    gid,
    weights,
    K,
    rho,
    gamma,
    prob.bound
) {
  q <- K - 1L
  times <- sort(unique(t))
  time_id <- match(t, times)
  m <- length(times)

  all_count <- event1_count <- event2_count <- matrix(0, m, K)
  for (i in seq_along(t)) {
    j <- time_id[i]
    k <- gid[i]
    all_count[j, k] <- all_count[j, k] + weights[i]
    if (epsilon[i] == 1L) {
      event1_count[j, k] <- event1_count[j, k] + weights[i]
    } else if (epsilon[i] == 2L) {
      event2_count[j, k] <- event2_count[j, k] + weights[i]
    }
  }

  risk <- colSums(all_count)
  survival_left <- rep.int(1, K)
  cif1_left <- numeric(K)
  pooled_cif_left <- 0
  score <- numeric(q)
  variance <- matrix(0, q, q)
  censor_link <- matrix(0, q, K)
  cross_term <- matrix(0, q, K)
  censor_variance <- numeric(K)
  event_time <- numeric()
  event_weight <- numeric()

  for (j in seq_len(m)) {
    d1 <- event1_count[j, ]
    d2 <- event2_count[j, ]
    d_all <- all_count[j, ]
    total_d1 <- sum(d1)
    total_d2 <- sum(d2)

    if (total_d1 + total_d2 > 0) {
      invalid <- risk > 0 & survival_left <= prob.bound
      if (any(invalid)) {
        stop("Gray covariance is undefined after a group survival estimate reaches zero.",
             call. = FALSE)
      }

      risk_over_survival <- numeric(K)
      active <- risk > 0
      risk_over_survival[active] <- risk[active] / survival_left[active]
      pooled_risk <- sum(risk_over_survival)
      if (!is.finite(pooled_risk) || pooled_risk <= 0) {
        stop("Gray's test encountered an empty pooled risk set.", call. = FALSE)
      }

      pooled_cif_right <- pooled_cif_left + total_d1 / pooled_risk
      gray_weight <- (1 - pooled_cif_left)^rho * pooled_cif_left^gamma

      influence_map <- -gray_weight *
        outer(risk_over_survival[seq_len(q)], risk_over_survival) / pooled_risk
      influence_map[cbind(seq_len(q), seq_len(q))] <-
        influence_map[cbind(seq_len(q), seq_len(q))] +
        gray_weight * risk_over_survival[seq_len(q)]

      if (total_d1 > 0) {
        event_time <- c(event_time, times[j])
        event_weight <- c(event_weight, gray_weight)
        cif_risk <- numeric(K)
        cif_risk[active] <- risk[active] * (1 - cif1_left[active]) /
          survival_left[active]
        cif_risk_total <- sum(cif_risk)
        if (!is.finite(cif_risk_total) || cif_risk_total <= 0) {
          stop("Gray's test encountered an empty subdistribution risk set.",
               call. = FALSE)
        }
        score <- score + gray_weight * (
          d1[seq_len(q)] - total_d1 * cif_risk[seq_len(q)] / cif_risk_total
        )

        remaining_cif <- 1 - pooled_cif_left
        if (remaining_cif <= prob.bound) {
          stop("Gray covariance is undefined when the pooled CIF reaches one.",
               call. = FALSE)
        }
        censor_link <- censor_link +
          influence_map * total_d1 / (pooled_risk * remaining_cif)
      }

      survival_right <- survival_left
      cif1_right <- cif1_left
      survival_right[active] <- survival_left[active] *
        (1 - (d1[active] + d2[active]) / risk[active])
      cif1_right[active] <- cif1_left[active] +
        survival_left[active] * d1[active] / risk[active]

      if (total_d1 > 0) {
        for (k in which(active)) {
          censor_ratio <- 1
          if (survival_right[k] > 0) {
            censor_ratio <- 1 - (1 - pooled_cif_right) / survival_right[k]
          }
          tie_adjust <- 1
          if (total_d1 > 1) {
            denominator <- pooled_risk * survival_left[k] - 1
            if (denominator <= 0) {
              stop("Gray target-event tie correction is undefined.", call. = FALSE)
            }
            tie_adjust <- 1 - (total_d1 - 1) / denominator
          }
          increment <- tie_adjust * survival_left[k] * total_d1 /
            (pooled_risk * risk[k])
          residual <- influence_map[, k] - censor_ratio * censor_link[, k]
          variance <- variance + tcrossprod(residual) * increment
          cross_term[, k] <- cross_term[, k] +
            residual * censor_ratio * increment
          censor_variance[k] <- censor_variance[k] +
            censor_ratio^2 * increment
        }
      }

      if (total_d2 > 0) {
        for (k in which(d2 > 0 & survival_right > 0)) {
          censor_ratio <- (1 - pooled_cif_right) / survival_right[k]
          tie_adjust <- 1
          if (d2[k] > 1) {
            if (risk[k] <= 1) {
              stop("Gray competing-event tie correction is undefined.", call. = FALSE)
            }
            tie_adjust <- 1 - (d2[k] - 1) / (risk[k] - 1)
          }
          increment <- tie_adjust * survival_left[k]^2 * d2[k] / risk[k]^2
          linked <- censor_ratio * censor_link[, k]
          variance <- variance + tcrossprod(linked) * increment
          cross_term[, k] <- cross_term[, k] -
            linked * censor_ratio * increment
          censor_variance[k] <- censor_variance[k] +
            censor_ratio^2 * increment
        }
      }

      pooled_cif_left <- pooled_cif_right
      survival_left <- survival_right
      cif1_left <- cif1_right
    }

    risk <- risk - d_all
    risk[abs(risk) < prob.bound] <- 0
  }

  variance <- variance +
    (censor_link * rep(censor_variance, each = q)) %*% t(censor_link) +
    censor_link %*% t(cross_term) + cross_term %*% t(censor_link)

  list(
    score = score,
    var = (variance + t(variance)) / 2,
    event.time = event_time,
    fh.weight = event_weight
  )
}

#' @export
print.ciftest <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Outcome: ", x$outcome.type, "\n", sep = "")
  cat("FH weights: rho = ", x$rho, ", gamma = ", x$gamma, "\n", sep = "")
  cat("Variance: ", x$variance.method, "\n\n", sep = "")
  cat(
    "Chi-squared = ", format(unname(x$statistic), digits = 5L),
    ", df = ", unname(x$parameter),
    ", p-value = ", format.pval(x$p.value, digits = 4L), "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.ciftest <- function(object, ...) object

#' @export
nobs.ciftest <- function(object, ...) object$n

#' @export
tidy.ciftest <- function(x, ...) {
  data.frame(
    statistic = unname(x$statistic),
    p.value = x$p.value,
    parameter = unname(x$parameter),
    method = x$method,
    outcome.type = x$outcome.type,
    n = x$n,
    stringsAsFactors = FALSE
  )
}

#' @export
glance.ciftest <- function(x, ...) tidy.ciftest(x, ...)

#' @export
augment.ciftest <- function(x, ...) {
  n0 <- nrow(x$data)
  p <- ncol(x$score.iid)
  restore <- function(value) {
    out <- matrix(NA_real_, nrow = n0, ncol = p, dimnames = list(NULL, colnames(value)))
    out[x$diagnostics$analysis.row.index, ] <- value
    I(out)
  }
  out <- x$data
  out$.score_iid <- restore(x$score.iid)
  out$.score_base <- restore(x$score.iid.base)
  out$.score_censoring <- restore(x$score.iid.censor)
  out$.score_augmentation <- restore(x$score.iid.augment)
  out$.score_iterated <- restore(x$score.iid.iterated)
  out$.censoring_weight <- NA_real_
  out$.analysis_included <- x$diagnostics$analysis.included
  out
}

#' @keywords internal
calculate_log_rank <- function(
    t,
    epsilon,
    exposure,
    code.exposure.ref = NULL,
    prefix = "a",
    weights,
    strata,
    data,
    rho = 0,
    gamma = 0,
    prob.bound = 1e-7
) {
  # --- basic checks ---
  n <- length(t)
  stopifnot(length(epsilon) == n, length(weights) == n, length(strata) == n)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(exposure) || length(exposure) != 1L) stop("`exposure` must be a single character string.")
  if (!exposure %in% names(data)) stop("exposure = '", exposure, "' is not found in data.")

  # --- exposure design (K levels -> K-1 dummies) ---
  exp_info <- reg_read_exposure_design(
    data = data,
    exposure = exposure,
    code.exposure.ref = code.exposure.ref,
    prefix = prefix
  )
  K <- exp_info$exposure.levels
  if (K < 2L) stop("Exposure must have >= 2 levels.")

  # Reconstruct the same factor coding/order as exp_info to get group IDs (1..K)
  a_ <- factor(data[[exposure]])
  a_ <- base::droplevels(a_)
  a_ <- factor(a_, levels = exp_info$exposure.labels) # ensures ref is level 1
  gid <- as.integer(a_)
  if (anyNA(gid)) stop("Missing values in exposure after preprocessing; handle via na.action before calling.")

  # event of interest indicator (epsilon==1)
  is_event1 <- (as.integer(epsilon) == 1L)

  # unique strata
  strata_fac <- factor(strata)
  L <- nlevels(strata_fac)

  p <- K - 1L
  score_total <- rep.int(0, p)
  var_total <- matrix(0, nrow = p, ncol = p)

  # helper: compute FH weights from pooled survival in a stratum
  fh_weights <- function(Yw, dNw, rho, gamma, prob.bound) {
    m <- length(Yw)
    Kt <- numeric(m)
    S <- 1
    for (j in seq_len(m)) {
      S_clip <- min(max(S, prob.bound), 1 - prob.bound)
      Kt[j] <- (S_clip^rho) * ((1 - S_clip)^gamma)
      if (is.finite(Yw[j]) && Yw[j] > 0) {
        S <- S * (1 - dNw[j] / Yw[j])
      }
    }
    Kt
  }

  # --- loop strata ---
  for (ll in seq_len(L)) {
    idx <- which(strata_fac == levels(strata_fac)[ll])
    if (length(idx) == 0L) next

    tt <- t[idx]
    ww <- weights[idx]
    ee <- is_event1[idx]
    gg <- gid[idx]

    # event times in this stratum
    times <- sort(unique(tt[ee]))
    M <- length(times)
    if (M == 0L) next

    # order by time for risk set suffix sums
    o <- order(tt)
    tt_o <- tt[o]
    ww_o <- ww[o]
    gg_o <- gg[o]

    # position of each event time within ordered times (first occurrence)
    pos <- match(times, tt_o)
    if (anyNA(pos)) stop("Internal error: could not match event times in ordered times.")

    # --- risk set sums by group: Yw_mat (M x K) ---
    Yw_mat <- matrix(0, nrow = M, ncol = K)
    for (k in seq_len(K)) {
      suf <- rev(cumsum(rev(ww_o * (gg_o == k))))
      Yw_mat[, k] <- suf[pos]
    }
    Yw <- rowSums(Yw_mat)

    # --- event sums by group at each event time: dNw_mat (M x K) ---
    dNw_mat <- matrix(0, nrow = M, ncol = K)
    ev_idx <- which(ee)
    if (length(ev_idx) > 0L) {
      m_idx <- match(tt[ev_idx], times)
      k_idx <- gg[ev_idx]
      for (j in seq_along(ev_idx)) {
        dNw_mat[m_idx[j], k_idx[j]] <- dNw_mat[m_idx[j], k_idx[j]] + ww[ev_idx[j]]
      }
    }
    dNw <- rowSums(dNw_mat)

    # --- FH weights K(t-) computed from pooled survival in this stratum ---
    Kt <- fh_weights(Yw = Yw, dNw = dNw, rho = rho, gamma = gamma, prob.bound = prob.bound)

    # --- score: U = sum_t K(t) (dN_g - Y_g/Y * dN) for non-ref groups ---
    # Avoid division by zero (shouldn't happen at event times, but guard anyway)
    P <- Yw_mat
    for (j in seq_len(M)) {
      if (Yw[j] > 0) P[j, ] <- Yw_mat[j, ] / Yw[j] else P[j, ] <- 0
    }
    U_full <- dNw_mat - P * dNw
    U_red  <- U_full[, -1, drop = FALSE]             # drop reference group
    score_l <- as.vector(crossprod(Kt, U_red))       # (K-1)-vector

    # --- variance-covariance: sum_t K(t)^2 * Cov( O-E at t ) ---
    # Multigroup log-rank covariance increment (survdiff-style):
    # c = d*(Y-d)/(Y*(Y-1));  Cov = c*( diag(Yg) - (Yg Yh)/Y )
    var_l <- matrix(0, nrow = p, ncol = p)
    for (j in seq_len(M)) {
      Y <- Yw[j]
      d <- dNw[j]
      if (!is.finite(Y) || !is.finite(d) || Y <= 1 || d <= 0) next
      if (Y - d < 0) next  # numerical guard for extreme weights
      cfac <- d * (Y - d) / (Y * (Y - 1))

      Yvec <- Yw_mat[j, ]
      C_full <- cfac * (diag(Yvec, nrow = K, ncol = K) - tcrossprod(Yvec) / Y)
      C_red  <- C_full[-1, -1, drop = FALSE]
      var_l  <- var_l + (Kt[j]^2) * C_red
    }

    score_total <- score_total + score_l
    var_total   <- var_total + var_l
  }

  colnames(var_total) <- rownames(var_total) <- colnames(exp_info$x_a)
  names(score_total) <- colnames(exp_info$x_a)

  list(
    score = score_total,
    var = var_total,
    df = length(score_total),
    exposure.levels = exp_info$exposure.levels,
    exposure.labels = exp_info$exposure.labels,
    ref = exp_info$ref
  )
}

#' Extract Mparts from a weightit object (robustly)
#' @keywords internal
get_weightit_mparts <- function(weightit) {
  if (!is.null(weightit$Mparts)) return(weightit$Mparts)
  mp <- attr(weightit, "Mparts", exact = TRUE)
  if (!is.null(mp)) return(mp)
  # NOTE: weightitMSM uses "Mparts.list" attribute (not for standard weightit)
  mp_list <- attr(weightit, "Mparts.list", exact = TRUE)
  if (!is.null(mp_list)) return(mp_list)
  stop("No Mparts found in `weightit`. Make sure the fit stored M-estimation parts.",
       call. = FALSE)
}

#' Choose a safe finite-difference step so w +/- h*delta stays >= 0
#' @keywords internal
fd_step_safe <- function(w, delta, rel = 1e-6, max_tries = 20L) {
  stopifnot(length(w) == length(delta))
  # target: max|h*delta| ≈ rel * max(1, median(w))
  scale_w <- max(1, stats::median(abs(w), na.rm = TRUE))
  max_abs_delta <- max(abs(delta), na.rm = TRUE)
  if (!is.finite(max_abs_delta) || max_abs_delta == 0) return(0)

  h <- rel * scale_w / max_abs_delta

  # ensure nonnegativity for both plus and minus
  idx <- which(abs(delta) > 0)
  if (length(idx)) {
    h_cap <- min(w[idx] / abs(delta[idx]), na.rm = TRUE)
    if (is.finite(h_cap)) h <- min(h, 0.49 * h_cap)
  }

  # backoff if still violates (due to numerical issues)
  for (k in seq_len(max_tries)) {
    if (all(w + h * delta >= 0) && all(w - h * delta >= 0)) return(h)
    h <- h * 0.5
  }
  stop("Could not find a safe finite-difference step (weights hit negative).",
       call. = FALSE)
}

#' A12 for log-rank-type score using directional finite differences with dw/dB from WeightIt
#'
#' Returns A12 on the "mean score" scale: (1/n) * dU_total/dB^T
#' @keywords internal
calculate_A12_logrank_weightit <- function(
    t,
    epsilon,
    strata,
    data,              # data_sync used by reg_read_exposure_design inside calculate_log_rank_components()
    exposure,
    weightit,
    code.exposure.ref = NULL,
    prefix = "a",
    rho = 0,
    gamma = 0,
    prob.bound = 1e-7,
    fd_rel_step = 1e-6
) {
  stopifnot(length(t) == length(epsilon),
            length(t) == length(strata))

  w0 <- weightit$weights
  if (is.null(w0)) stop("weightit$weights is NULL.", call. = FALSE)
  if (length(w0) != length(t)) {
    stop("Length mismatch: length(weightit$weights) != length(t).",
         call. = FALSE)
  }

  mp <- get_weightit_mparts(weightit)
  dw <- mp$dw_dBtreat
  if (is.null(dw)) stop("Mparts$dw_dBtreat is missing.", call. = FALSE)

  dw <- as.matrix(dw)
  if (nrow(dw) != length(t)) {
    stop("nrow(dw_dBtreat) != length(t). Check row alignment after subsetting/na.omit.",
         call. = FALSE)
  }

  # baseline score (total score, not averaged)
  base <- calculate_log_rank(
    t = t, epsilon = epsilon, weights = w0, strata = strata,
    data = data, exposure = exposure,
    code.exposure.ref = code.exposure.ref, prefix = prefix,
    rho = rho, gamma = gamma, prob.bound = prob.bound
  )
  if (is.null(base$score)) stop("calculate_log_rank() must return $score.", call. = FALSE)

  score0 <- as.numeric(base$score)
  p <- length(score0)
  n <- length(t)
  q <- ncol(dw)

  A12 <- matrix(NA_real_, nrow = p, ncol = q)
  rownames(A12) <- names(base$score) %||% paste0("score", seq_len(p))
  colnames(A12) <- colnames(dw) %||% paste0("Btreat", seq_len(q))

  steps <- numeric(q)

  for (j in seq_len(q)) {
    delta <- dw[, j]
    h <- fd_step_safe(w0, delta, rel = fd_rel_step)
    steps[j] <- h

    if (h == 0) {
      A12[, j] <- 0
      next
    }

    w_plus  <- w0 + h * delta
    w_minus <- w0 - h * delta

    # Compute scores at perturbed weights
    s_plus <- calculate_log_rank(
      t = t, epsilon = epsilon, weights = w_plus, strata = strata,
      data = data, exposure = exposure,
      code.exposure.ref = code.exposure.ref, prefix = prefix,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    s_minus <- calculate_log_rank(
      t = t, epsilon = epsilon, weights = w_minus, strata = strata,
      data = data, exposure = exposure,
      code.exposure.ref = code.exposure.ref, prefix = prefix,
      rho = rho, gamma = gamma, prob.bound = prob.bound
    )$score

    dU_dB_j <- (as.numeric(s_plus) - as.numeric(s_minus)) / (2 * h)

    # Put on mean-score scale: A12 = (1/n) * dU_total/dB^T
    A12[, j] <- dU_dB_j / n
  }

  list(A12 = A12, fd_steps = steps, score0 = base$score)
}
