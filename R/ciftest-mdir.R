#' Multiple-direction quadratic score test
#'
#' Combines several Fleming--Harrington score directions in the quadratic
#' form proposed by Ditzhaus and Friedrich. Each component is fitted with the
#' same outcome-specific `ciftest()` branch. Classical log-rank components use
#' the joint, finite-population-corrected hypergeometric covariance, with
#' cross-direction blocks weighted by `w_r(t) * w_s(t)`. Classical Gray
#' components analogously use the joint Gray (1988) covariance, including the
#' target-event, competing-event, and censoring contributions. Score and
#' augmented components use the cross-product of the stacked individual score
#' contributions. The default chi-square reference uses the numerical rank of
#' the joint covariance. This function does not perform the permutation
#' calibration from the original two-sample survival proposal.
#'
#' @inheritParams ciftest
#' @param test Underlying component test: `"auto"`, `"logrank"`, `"gray"`,
#'   `"score"`, or `"augmented"`. The corresponding aliases documented for
#'   [ciftest()] are accepted. Put early and late choices in `directions`.
#' @param directions Fleming--Harrington directions. A character vector may
#'   contain `"unweighted"`, `"early"`, and `"late"`. Alternatively supply a
#'   two-column numeric matrix/data frame with `rho` and `gamma`, or a named
#'   list of numeric length-two vectors. In this multiple-direction interface,
#'   the default character presets correspond to `(2, 0)`, `(0, 2)`, and
#'   `(0, 0)`. The stronger early/late directions keep the
#'   default set linearly independent.
#'
#' @return An object inheriting from `"ciftest_mdir"`, `"ciftest"`, and
#'   `"htest"`. Component fits are available in `components`; the unmodified
#'   individual contributions are in `score.iid`, and the covariance actually
#'   used by the quadratic statistic is in `vcov.score`.
#' @references Ditzhaus M, Friedrich S. More powerful logrank permutation
#'   tests for two-sample survival data. arXiv preprint (2018).
#' @export
ciftest_mdir <- function(
    formula,
    data,
    directions = c("early", "late", "unweighted"),
    weights = NULL,
    subset.condition = NULL,
    outcome.type = "auto",
    test = "auto",
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    iteration = 0L,
    tolerance = NULL,
    strata = NULL,
    strata.censor = NULL,
    strata.competing.risk = NULL,
    tau = NULL,
    prob.bound = 1e-7,
    prob.truncation = NULL,
    na.action = stats::na.omit
) {
  call <- match.call()
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("`test` must be one non-missing character value.", call. = FALSE)
  }
  if (tolower(test) %in% c("early", "late")) {
    stop("Specify early and late weights in `directions`, not in `test`.",
         call. = FALSE)
  }
  direction_table <- ciftest_mdir_directions(directions)

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
    test = test,
    code.event1 = code.event1,
    code.event2 = code.event2,
    code.censoring = code.censoring,
    rho = direction_table$rho[1L],
    gamma = direction_table$gamma[1L],
    iteration = iteration,
    tolerance = tolerance,
    strata = strata,
    strata.censor = strata.censor,
    strata.competing.risk = strata.competing.risk,
    tau = tau,
    prob.bound = prob.bound,
    prob.truncation = prob.truncation,
    na.action = na.action
  )
  cache <- new_ciftest_nuisance_cache()
  fits <- lapply(seq_len(nrow(direction_table)), function(index) {
    component_input <- input
    component_input$rho <- direction_table$rho[index]
    component_input$gamma <- direction_table$gamma[index]
    component_input$weight.label <- direction_table$direction[index]
    ciftest_fit_prepared(
      input = component_input,
      formula = formula,
      call = call,
      nuisance.cache = cache
    )
  })
  names(fits) <- direction_table$direction

  unavailable <- vapply(fits, function(fit) {
    is.null(fit$score.iid) || any(!is.finite(fit$score.iid))
  }, logical(1L))
  classical_covariance <- input$test %in% c("logrank", "gray")
  if (any(unavailable) && !classical_covariance) {
    stop(
      "Joint individual score contributions are unavailable for direction(s): ",
      paste(names(fits)[unavailable], collapse = ", "),
      ". Use `test = \"score\"` or `test = \"augmented\"`, or inspect the component diagnostic.",
      call. = FALSE
    )
  }

  score <- unlist(lapply(fits, function(fit) unname(fit$score)),
                  use.names = FALSE)
  contribution_blocks <- Map(function(fit, bad) {
    if (!bad) return(fit$score.iid)
    matrix(
      NA_real_, nrow = length(input$t), ncol = length(fit$score),
      dimnames = list(NULL, names(fit$score))
    )
  }, fits, unavailable)
  block_names <- unlist(Map(function(label, block, fit) {
    component_names <- colnames(block)
    if (is.null(component_names)) component_names <- names(fit$score)
    paste(label, component_names, sep = "::")
  }, direction_table$direction, contribution_blocks, fits),
  use.names = FALSE)
  score_iid <- do.call(cbind, contribution_blocks)
  colnames(score_iid) <- names(score) <- block_names

  raw_covariance <- if (any(unavailable)) {
    matrix(NA_real_, ncol(score_iid), ncol(score_iid))
  } else {
    crossprod(score_iid)
  }
  covariance_source <- "score-iid cross-product"
  if (identical(input$test, "logrank")) {
    covariance <- ciftest_mdir_logrank_covariance(input, direction_table)
    covariance_source <- "joint hypergeometric"
  } else if (identical(input$test, "gray")) {
    covariance <- ciftest_mdir_gray_covariance(input, direction_table)
    covariance_source <- "joint Gray"
  } else {
    covariance <- raw_covariance
  }
  calibration <- diag(ncol(covariance))
  dimnames(covariance) <- list(block_names, block_names)
  dimnames(calibration) <- list(block_names, block_names)
  statistic <- ciftest_quadratic_form(score, covariance)

  out <- list(
    statistic = stats::setNames(statistic$statistic, "X-squared"),
    parameter = stats::setNames(statistic$rank, "df"),
    p.value = stats::pchisq(
      statistic$statistic, statistic$rank, lower.tail = FALSE
    ),
    method = paste0(
      "Multiple-direction ", fits[[1L]]$method,
      " (asymptotic chi-square)"
    ),
    data.name = paste(deparse(formula), collapse = " "),
    call = call,
    outcome.type = input$outcome.type,
    test = input$test,
    directions = direction_table,
    score = score,
    vcov.score = covariance,
    vcov.score.raw = raw_covariance,
    covariance.calibration = calibration,
    score.iid = score_iid,
    components = fits,
    n = length(input$t),
    diagnostics = list(
      covariance.rank = statistic$rank,
      component.count = nrow(direction_table),
      covariance.source = covariance_source,
      nuisance.cache.hits = cache$hits,
      nuisance.cache.misses = cache$misses,
      permutation.calibration = FALSE,
      analysis.row.index = input$row.index,
      analysis.included = input$included,
      exclusion.reason = input$exclusion.reason
    ),
    data = input$data.original
  )
  class(out) <- c("ciftest_mdir", "ciftest", "htest")
  out
}

ciftest_mdir_directions <- function(directions) {
  presets <- rbind(
    unweighted = c(rho = 0, gamma = 0),
    early = c(rho = 2, gamma = 0),
    late = c(rho = 0, gamma = 2)
  )
  if (is.character(directions)) {
    key <- tolower(trimws(directions))
    if (!length(key) || anyNA(key) || any(!key %in% rownames(presets))) {
      stop("Character `directions` must use 'unweighted', 'early', or 'late'.",
           call. = FALSE)
    }
    values <- presets[key, , drop = FALSE]
    labels <- key
  } else if (is.data.frame(directions) || is.matrix(directions)) {
    values <- as.matrix(directions)
    if (!is.numeric(values) || ncol(values) != 2L || !nrow(values)) {
      stop("Numeric `directions` must have two columns and at least one row.",
           call. = FALSE)
    }
    if (!is.null(colnames(values)) &&
        all(c("rho", "gamma") %in% colnames(values))) {
      values <- values[, c("rho", "gamma"), drop = FALSE]
    } else {
      colnames(values) <- c("rho", "gamma")
    }
    labels <- rownames(values)
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- paste0("direction", seq_len(nrow(values)))
    }
  } else if (is.list(directions) && length(directions)) {
    valid <- vapply(directions, function(value) {
      is.numeric(value) && length(value) == 2L && all(is.finite(value))
    }, logical(1L))
    if (!all(valid)) {
      stop("Every list element in `directions` must be a finite numeric pair.",
           call. = FALSE)
    }
    values <- do.call(rbind, directions)
    colnames(values) <- c("rho", "gamma")
    labels <- names(directions)
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- paste0("direction", seq_len(nrow(values)))
    }
  } else {
    stop("`directions` has an unsupported format.", call. = FALSE)
  }
  storage.mode(values) <- "double"
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("Direction rho and gamma values must be finite and non-negative.",
         call. = FALSE)
  }
  if (anyDuplicated(data.frame(rho = values[, 1L], gamma = values[, 2L]))) {
    stop("`directions` must not contain duplicate rho/gamma pairs.",
         call. = FALSE)
  }
  labels <- make.unique(as.character(labels), sep = "_")
  data.frame(
    direction = labels,
    rho = values[, 1L],
    gamma = values[, 2L],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Joint survdiff-style hypergeometric covariance for several FH directions.
# Direction is the outer block index and exposure contrast is the inner index.
ciftest_mdir_logrank_covariance <- function(input, directions) {
  exp_info <- reg_read_exposure_design(
    data = input$data,
    exposure = input$exposure,
    code.exposure.ref = NULL,
    prefix = "a"
  )
  K <- exp_info$exposure.levels
  p <- K - 1L
  m_direction <- nrow(directions)
  group <- factor(input$data[[input$exposure]],
                  levels = exp_info$exposure.labels)
  gid <- as.integer(group)
  event <- as.integer(input$epsilon == input$code.event1)
  stratum <- factor(input$strata.event.info$values)
  covariance <- matrix(0, m_direction * p, m_direction * p)

  for (level in levels(stratum)) {
    use <- which(stratum == level)
    if (!length(use)) next
    time <- input$t[use]
    frequency <- input$weights[use]
    status <- event[use]
    group_id <- gid[use]
    event_time <- sort(unique(time[status == 1L]))
    if (!length(event_time)) next

    ordering <- order(time)
    ordered_time <- time[ordering]
    ordered_frequency <- frequency[ordering]
    ordered_group <- group_id[ordering]
    position <- match(event_time, ordered_time)

    risk <- matrix(0, length(event_time), K)
    for (k in seq_len(K)) {
      suffix <- rev(cumsum(rev(
        ordered_frequency * (ordered_group == k)
      )))
      risk[, k] <- suffix[position]
    }
    total_risk <- rowSums(risk)

    death <- matrix(0, length(event_time), K)
    event_index <- which(status == 1L)
    for (i in event_index) {
      j <- match(time[i], event_time)
      death[j, group_id[i]] <- death[j, group_id[i]] + frequency[i]
    }
    total_death <- rowSums(death)

    survival <- 1
    for (j in seq_along(event_time)) {
      Y <- total_risk[j]
      d <- total_death[j]
      survival_used <- min(max(survival, input$prob.bound),
                           1 - input$prob.bound)
      direction_weight <- survival_used^directions$rho *
        (1 - survival_used)^directions$gamma

      if (is.finite(Y) && is.finite(d) && Y > 1 && d > 0 && Y >= d) {
        finite_population <- d * (Y - d) / (Y * (Y - 1))
        group_risk <- risk[j, ]
        full_increment <- finite_population * (
          diag(group_risk, nrow = K, ncol = K) -
            tcrossprod(group_risk) / Y
        )
        reduced_increment <- full_increment[-1, -1, drop = FALSE]
        covariance <- covariance + kronecker(
          tcrossprod(direction_weight), reduced_increment
        )
      }
      if (is.finite(Y) && Y > 0) survival <- survival * (1 - d / Y)
    }
  }
  (covariance + t(covariance)) / 2
}

# Joint Gray (1988) covariance for several FH directions. This is the
# bilinear extension of gray_stratum_components(): stacking the direction-
# specific influence coefficients makes every diagonal block exactly the
# corresponding scalar Gray covariance and supplies the required cross blocks.
ciftest_mdir_gray_covariance <- function(input, directions) {
  exp_info <- reg_read_exposure_design(
    data = input$data,
    exposure = input$exposure,
    code.exposure.ref = NULL,
    prefix = "a"
  )
  K <- exp_info$exposure.levels
  q <- K - 1L
  m_direction <- nrow(directions)
  group <- factor(input$data[[input$exposure]],
                  levels = exp_info$exposure.labels)
  gid <- as.integer(group)
  stratum <- droplevels(factor(input$strata.event.info$values))
  covariance_base <- matrix(0, m_direction * q, m_direction * q)

  for (level in levels(stratum)) {
    use <- which(stratum == level & input$weights > 0)
    if (!length(use)) next
    covariance_base <- covariance_base + ciftest_mdir_gray_stratum_covariance(
      t = input$t[use],
      epsilon = input$epsilon[use],
      gid = gid[use],
      weights = as.numeric(input$weights[use]),
      K = K,
      rho = directions$rho,
      gamma = directions$gamma,
      prob.bound = input$prob.bound
    )
  }

  transform <- matrix(0, q, q)
  if (q > 1L) transform[cbind(seq_len(q - 1L), 2:q)] <- 1
  transform[q, ] <- -1
  joint_transform <- kronecker(diag(m_direction), transform)
  covariance <- joint_transform %*% covariance_base %*% t(joint_transform)
  (covariance + t(covariance)) / 2
}

ciftest_mdir_gray_stratum_covariance <- function(
    t, epsilon, gid, weights, K, rho, gamma, prob.bound
) {
  q <- K - 1L
  m_direction <- length(rho)
  dimension <- m_direction * q
  times <- sort(unique(t))
  time_id <- match(t, times)
  all_count <- event1_count <- event2_count <- matrix(0, length(times), K)
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
  variance <- matrix(0, dimension, dimension)
  censor_link <- matrix(0, dimension, K)
  cross_term <- matrix(0, dimension, K)
  censor_variance <- numeric(K)

  for (j in seq_along(times)) {
    d1 <- event1_count[j, ]
    d2 <- event2_count[j, ]
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
      direction_weight <- (1 - pooled_cif_left)^rho * pooled_cif_left^gamma

      influence_map <- matrix(0, dimension, K)
      for (direction in seq_len(m_direction)) {
        block <- (direction - 1L) * q + seq_len(q)
        derivative <- -direction_weight[direction] *
          outer(risk_over_survival[seq_len(q)], risk_over_survival) /
          pooled_risk
        derivative[cbind(seq_len(q), seq_len(q))] <-
          derivative[cbind(seq_len(q), seq_len(q))] +
          direction_weight[direction] * risk_over_survival[seq_len(q)]
        influence_map[block, ] <- derivative
      }

      if (total_d1 > 0) {
        remaining_cif <- 1 - pooled_cif_left
        if (remaining_cif <= prob.bound) {
          stop("Gray covariance is undefined when the pooled CIF reaches one.",
               call. = FALSE)
        }
        censor_link <- censor_link + influence_map * total_d1 /
          (pooled_risk * remaining_cif)
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
    risk <- risk - all_count[j, ]
    risk[abs(risk) < prob.bound] <- 0
  }

  variance <- variance +
    (censor_link * rep(censor_variance, each = dimension)) %*%
      t(censor_link) +
    censor_link %*% t(cross_term) + cross_term %*% t(censor_link)
  (variance + t(variance)) / 2
}

#' @export
print.ciftest_mdir <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Outcome: ", x$outcome.type, "\n", sep = "")
  cat("Component test: ", x$test, "\n", sep = "")
  cat(
    "Directions: ",
    paste0(
      x$directions$direction, " (", x$directions$rho, ", ",
      x$directions$gamma, ")", collapse = ", "
    ),
    "\n\n", sep = ""
  )
  cat(
    "Chi-squared = ", format(unname(x$statistic), digits = 5L),
    ", df = ", unname(x$parameter),
    ", p-value = ", format.pval(x$p.value, digits = 4L), "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
tidy.ciftest_mdir <- function(x, ...) {
  data.frame(
    statistic = unname(x$statistic),
    p.value = x$p.value,
    parameter = unname(x$parameter),
    method = x$method,
    outcome.type = x$outcome.type,
    test = x$test,
    directions = nrow(x$directions),
    n = x$n,
    stringsAsFactors = FALSE
  )
}

#' @export
glance.ciftest_mdir <- function(x, ...) tidy.ciftest_mdir(x, ...)

#' @export
augment.ciftest_mdir <- function(x, ...) {
  n0 <- nrow(x$data)
  restore <- matrix(
    NA_real_, nrow = n0, ncol = ncol(x$score.iid),
    dimnames = list(NULL, colnames(x$score.iid))
  )
  restore[x$diagnostics$analysis.row.index, ] <- x$score.iid
  out <- x$data
  out$.score_iid <- I(restore)
  out$.analysis_included <- x$diagnostics$analysis.included
  out
}
