estimating_equation_ipcw <- function(
    formula,
    data,
    outcome.type,
    exposure,
    ip.weight.matrix,
    alpha_beta,
    estimand,
    optim.method,
    prob.bound,
    initial.CIFs = NULL
){
  cl <- match.call()
  mf <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster", "offset")
  out_terms <- terms(formula, special, data = data)
  mf$formula <- out_terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  Y <- model.extract(mf, "response")
  t <- Y[, 1]
  epsilon <- Y[, 2]

  if (!is.null(attributes(out_terms)$specials$offset)) {
    ts <- untangle.specials(out_terms, "offset")
    if (length(ts$vars) > 0) {
      out_terms <- out_terms[-ts$terms]
      offset <- mf[[ts$vars]]
    } else {
      offset <- rep(0, nrow(mf))
    }
  } else {
    offset <- rep(0, nrow(mf))
  }
  if (is.null(estimand$time.point)) {
    time.point <- sort(unique(with(data, t[epsilon > 0])))
  } else {
    time.point <- as.numeric(estimand$time.point)
  }

  out_readExposureDesign <- reg_read_exposure_design(data, exposure, estimand$code.exposure.ref)
  x_a <- out_readExposureDesign$x_a
  x_l <- model.matrix(out_terms, mf)
  x_la <- cbind(x_l, x_a)

  n  <- nrow(x_l)
  ql <- ncol(x_l)
  qa <- ncol(x_a)
  q  <- ql + qa
  iv <- estimand$index.vector
  if (ql != iv[1L] || q != iv[3L]) {
    stop("The regression design does not match the fitted coefficient layout.", call. = FALSE)
  }

  if (outcome.type == "survival" || outcome.type == "binomial") {
    y_0 <- ifelse(epsilon == estimand$code.censoring | t > estimand$time.point, 1, 0)
    y_1 <- ifelse(epsilon == estimand$code.event1   & t <= estimand$time.point, 1, 0)
    y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
    y_1_ <- ifelse(epsilon == estimand$code.event1,   1, 0)

    potential.Risk <- getPotentialRisk(alpha_beta, x_a, x_l, offset, estimand, optim.method, prob.bound)
    ey_1 <- calculateEY1(potential.Risk, x_a)

    ip.weight <- if (is.null(ip.weight.matrix)) rep(1, nrow(x_l)) else as.vector(ip.weight.matrix)

    w11 <- 1 / (ey_1 * (1 - ey_1))
    wy_1 <- ip.weight * y_1
    wy_1ey_1 <- w11*(wy_1 - ey_1)

    d <- cbind(x_l, x_a)
    residual <- wy_1ey_1
    ret <- drop(crossprod(d, residual)) / nrow(x_l)

    n_col_d <- ncol(d)
    score.matrix <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
    for (j in seq_len(n_col_d)) score.matrix[,j] <- d[,j]*residual

    return(list(
      ret   = ret,
      score = score.matrix,
      ey_1  = ey_1,
      w11   = w11,
      t     = t,
      y_0   = y_0,
      y_1   = y_1,
      y_0_  = y_0_,
      y_1_  = y_1_,
      x_a   = x_a,
      x_l   = x_l,
      potential.CIFs = potential.Risk,
      ip.weight = ip.weight
    ))
  }

  if (outcome.type == "proportional-survival") {
    if (is.null(estimand$time.point)) {
      tp <- with(data, t[epsilon > 0]); time.point <- unique(tp)
    }

    y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
    y_1_ <- ifelse(epsilon == estimand$code.event1,   1, 0)

    layout <- reg_proportional_index(estimand)
    score <- numeric(length(alpha_beta))
    i_time <- 0L

    for (specific.time in time.point) {
      i_time <- i_time + 1L
      index <- c(layout$nuisance[, i_time], layout$exposure)
      alpha_beta_i <- alpha_beta[index]

      y_0 <- ifelse(epsilon == estimand$code.censoring | t > specific.time, 1, 0)
      y_1 <- ifelse(epsilon == estimand$code.event1   & t <= specific.time, 1, 0)

      potential.Risk <- getPotentialRisk(alpha_beta_i, x_a, x_l, offset, estimand, optim.method, prob.bound)
      ey_1 <- calculateEY1(potential.Risk, x_a)

      w11 <- 1 / (ey_1 * (1 - ey_1))
      wy_1 <- ip.weight.matrix[,i_time] * y_1
      wy_1ey_1 <- w11*(wy_1 - ey_1)

      d <- cbind(x_l, x_a)
      residual <- wy_1ey_1

      score[index] <- score[index] + drop(crossprod(d, residual)) / n
    }

    return(list(
      ret   = score,
      t     = t,
      y_0   = y_0,
      y_1   = y_1,
      y_0_  = y_0_,
      y_1_  = y_1_,
      x_a   = x_a,
      x_l   = x_l
    ))
  }

  zero <- matrix(0, nrow=n, ncol=q)
  tmp1 <- cbind(x_la, zero)
  tmp2 <- cbind(zero, x_la)
  d    <- rbind(tmp1, tmp2)

  if (length(time.point) == 1) {
    y_0 <- ifelse(epsilon == estimand$code.censoring | t > time.point[1], 1, 0)
    y_1 <- ifelse(epsilon == estimand$code.event1 & t <= time.point[1], 1, 0)
    y_2 <- ifelse(epsilon == estimand$code.event2 & t <= time.point[1], 1, 0)
    y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
    y_1_ <- ifelse(epsilon == estimand$code.event1, 1, 0)
    y_2_ <- ifelse(epsilon == estimand$code.event2, 1, 0)

    potential.CIFs <- calculatePotentialCIFs(
      alpha_beta, x_a, x_l, offset, y_1, y_2, estimand, optim.method, prob.bound, initial.CIFs
    )
    ey <- calculateEY(potential.CIFs, x_a)
    ey_1 <- ey$ey_1
    ey_2 <- ey$ey_2

    v11 <- ey_1 * (1 - ey_1)
    v12 <- -ey_1 * ey_2
    v22 <- ey_2 * (1 - ey_2)
    denom <- v11*v22 - v12*v12
    w11 <- v22 / denom
    w12 <- -v12 / denom
    w22 <- v11 / denom

    wy_1 <- as.vector(ip.weight.matrix) * y_1
    wy_2 <- as.vector(ip.weight.matrix) * y_2
    wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
    wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

    residual <- c(wy_1ey_1, wy_2ey_2)
    ret <- drop(crossprod(d, residual)) / nrow(x_l)

    n_col_d <- ncol(d)
    score.matrix <- matrix(NA, nrow=nrow(d), ncol=n_col_d)
    for (j in seq_len(n_col_d)) score.matrix[,j] <- d[,j]*residual

    return(list(
      ret   = ret,
      score = score.matrix,
      ey_1  = ey_1,
      ey_2  = ey_2,
      w11   = w11,
      w12   = w12,
      w22   = w22,
      t     = t,
      y_0   = y_0,
      y_1   = y_1,
      y_2   = y_2,
      y_0_  = y_0_,
      y_1_  = y_1_,
      y_2_  = y_2_,
      x_a   = x_a,
      x_l   = x_l,
      potential.CIFs = potential.CIFs,
      ip.weight.matrix = ip.weight.matrix
    ))
  }

  y_0_ <- ifelse(epsilon == estimand$code.censoring, 1, 0)
  y_1_ <- ifelse(epsilon == estimand$code.event1, 1, 0)
  y_2_ <- ifelse(epsilon == estimand$code.event2, 1, 0)

  layout <- reg_proportional_index(estimand)
  score <- numeric(length(alpha_beta))
  i_time <- 0L

  for (specific.time in time.point) {
    i_time <- i_time + 1L
    index <- c(layout$nuisance[, i_time], layout$exposure)
    index <- c(index, index + layout$block)
    alpha_beta_i <- alpha_beta[index]

    y_0 <- ifelse(epsilon == estimand$code.censoring | t > specific.time, 1, 0)
    y_1 <- ifelse(epsilon == estimand$code.event1 & t <= specific.time, 1, 0)
    y_2 <- ifelse(epsilon == estimand$code.event2 & t <= specific.time, 1, 0)

    potential.CIFs <- calculatePotentialCIFs(
      alpha_beta_i, x_a, x_l, offset, y_1_, y_2_, estimand, optim.method, prob.bound, initial.CIFs
    )
    ey <- calculateEY(potential.CIFs, x_a)
    ey_1 <- ey$ey_1
    ey_2 <- ey$ey_2

    v11 <- ey_1 * (1 - ey_1)
    v12 <- -ey_1 * ey_2
    v22 <- ey_2 * (1 - ey_2)
    denom <- v11*v22 - v12*v12

    w11 <- v22 / denom
    w12 <- -v12 / denom
    w22 <- v11 / denom

    wy_1 <- ip.weight.matrix[, i_time] * y_1
    wy_2 <- ip.weight.matrix[, i_time] * y_2

    wy_1ey_1 <- w11*(wy_1 - ey_1) + w12*(wy_2 - ey_2)
    wy_2ey_2 <- w12*(wy_1 - ey_1) + w22*(wy_2 - ey_2)

    residual <- c(wy_1ey_1, wy_2ey_2)
    score[index] <- score[index] + drop(crossprod(d, residual)) / n
  }

  return(list(
    ret   = score,
    t     = t,
    y_0   = y_0,
    y_1   = y_1,
    y_2   = y_2,
    y_0_  = y_0_ ,
    y_1_  = y_1_ ,
    y_2_  = y_2_ ,
    x_a   = x_a,
    x_l   = x_l
  ))
}

calculateEY <- function(potential.CIFs, x_a) {
  n <- nrow(x_a)
  K <- ncol(x_a) + 1
  if (nrow(potential.CIFs) != n) stop("Row dimension of potential.CIFs is more than rows of x_a")
  if (ncol(potential.CIFs) != 2*K) stop("Column dimension of potential.CIFs is not 2*K")
  index <- max.col(cbind(0, x_a), ties.method = "first")
  CIF1 <- potential.CIFs[, seq_len(K), drop = FALSE]
  CIF2 <- potential.CIFs[, K + seq_len(K), drop = FALSE]
  ey_1 <- CIF1[cbind(seq_len(n), index)]
  ey_2 <- CIF2[cbind(seq_len(n), index)]
  return(list(ey_1 = ey_1, ey_2 = ey_2))
}

calculateEY1 <- function(potential.Risk, x_a) {
  n <- nrow(x_a)
  K <- ncol(x_a) + 1
  if (nrow(potential.Risk) != n) stop("Row dimension of potential.Risk is inconsistent with x_a")
  if (ncol(potential.Risk) != K) stop("Column dimension of potential.Risk must equal exposure levels")
  index <- max.col(cbind(0, x_a), ties.method = "first")
  potential.Risk[cbind(seq_len(n), index)]
}

polyreg_score_sensitivity <- function(equation, coef) {
  # Use the fitted estimating equation, including the coupled potential risks.
  step <- 1e-4 * pmax(1, abs(coef))
  -vapply(seq_along(coef), function(j) {
    lower <- upper <- coef
    lower[j] <- coef[j] - step[j]
    upper[j] <- coef[j] + step[j]
    (equation(upper) - equation(lower)) / (2 * step[j])
  }, numeric(length(coef)))
}

calculateCov <- function(objget_results, estimand, boot.method, ipcw, sensitivity)
{
  score <- objget_results$score
  w11 <- objget_results$w11
  w12 <- objget_results$w12
  w22 <- objget_results$w22
  t <- objget_results$t
  y_1 <- objget_results$y_1
  y_2 <- objget_results$y_2
  x_a <- objget_results$x_a
  x_l <- objget_results$x_l
  iv <- estimand$index.vector
  n <- length(t)

  x_la <- cbind(x_l, x_a)
  AB1 <- score[1:n, 1:iv[3], drop = FALSE]
  AB2 <- score[(n + 1):(2 * n), iv[4]:iv[7], drop = FALSE]
  weight <- drop(ipcw$weights)
  weighted_score <- cbind(
    x_la * (weight * (w11 * y_1 + w12 * y_2)),
    x_la * (weight * (w12 * y_1 + w22 * y_2))
  )
  total_score <- cbind(AB1, AB2) + censoring_km_score_correction(
    ipcw$km, t, ipcw$epsilon, ipcw$strata, pmin(t, estimand$time.point), weighted_score
  )
  centered_score <- sweep(total_score, 2L, colMeans(total_score), `-`)
  influence.function <- t(solve(sensitivity, t(centered_score)))

  cov_estimated <- NULL
  cov_bootstrap <- NULL
  if (isTRUE(boot.method$report.sandwich.conf))
  {
    cov_estimated <- crossprod(influence.function) / n^2
  }
  if (isTRUE(boot.method$report.boot.conf))
  {
    cov_bootstrap <- cov_wild_bootstrap(influence.function=influence.function, multiplier=boot.method$boot.multiplier, R=boot.method$boot.replications, seed=boot.method$boot.seed)
  }

  if (ncol(influence.function)==4) {
    colnames(influence.function) <- c("intercept", "exposure", "intercept", "exposure")
  } else if (ncol(influence.function)>4) {
    colnames(influence.function) <- c("intercept", paste("covariate", 1:(ncol(influence.function)/2 - 2), sep = ""), "exposure",
                                      "intercept", paste("covariate", 1:(ncol(influence.function)/2 - 2), sep = ""), "exposure")
  }
  return(list(cov_estimated = cov_estimated, cov_bootstrap = cov_bootstrap, score.function = total_score, influence.function = influence.function))
}

calculateCovSurvival <- function(objget_results, estimand, boot.method, prob.bound,
                                 ipcw = NULL, sensitivity = NULL)
{
  score <- objget_results$score
  w11 <- objget_results$w11
  t <- objget_results$t
  y_1 <- objget_results$y_1
  x_a <- objget_results$x_a
  x_l <- objget_results$x_l
  potential.CIFs <- objget_results$potential.CIFs
  n <- length(t)

  x_la <- cbind(x_l, x_a)
  total_score <- score
  if (!is.null(ipcw)) {
    weighted_score <- x_la * (drop(ipcw$weights) * w11 * y_1)
    total_score <- total_score + censoring_km_score_correction(
      ipcw$km, t, ipcw$epsilon, ipcw$strata, pmin(t, estimand$time.point), weighted_score
    )
  }

  if (is.null(sensitivity)) {
    derivative <- calculateDSurvival(potential.CIFs, x_a, x_l, estimand, prob.bound)
    sensitivity <- crossprod(x_la, w11 * derivative) / n
  }
  centered_score <- sweep(total_score, 2L, colMeans(total_score), `-`)
  influence.function <- t(solve(sensitivity, t(centered_score)))

  cov_estimated <- NULL
  cov_bootstrap <- NULL
  if (isTRUE(boot.method$report.sandwich.conf)) {
    cov_estimated <- crossprod(influence.function) / n^2
  }
  if (isTRUE(boot.method$report.boot.conf)) {
    cov_bootstrap <- cov_wild_bootstrap(
      influence.function = influence.function,
      multiplier = boot.method$boot.multiplier,
      R = boot.method$boot.replications,
      seed = boot.method$boot.seed
    )
  }

  qa <- ncol(x_a); ql <- ncol(x_l)
  if (ncol(influence.function) > 0) {
    intercept_name <- if (ql >= 1L) "intercept" else character(0)
    covar_names    <- if (ql > 1L) paste0("covariate", seq_len(ql - 1L)) else character(0)
    exposure_names <- if (qa == 1L) "exposure" else paste0("exposure", seq_len(qa))
    colnames(influence.function) <- c(intercept_name, covar_names, exposure_names)
  }

  list(
    cov_estimated = cov_estimated,
    cov_bootstrap = cov_bootstrap,
    score.function = total_score,
    influence.function = influence.function
  )
}

calculateDSurvival <- function(potential.CIFs, x_a, x_l, estimand, prob.bound) {
  n <- nrow(x_l)
  K <- estimand$exposure.levels
  if (ncol(potential.CIFs) < K) stop("Column dimension of potential.CIFs is insufficient for the number of exposure levels.")

  seq_n <- seq_len(n)
  idx <- max.col(cbind(0, x_a), ties.method = "first")

  CIF_mat <- clampP(potential.CIFs[, seq_len(K), drop = FALSE], prob.bound)
  CIF_sel <- CIF_mat[cbind(seq_n, idx)]

  calculateA <- function(effect_measure, CIFs_selected) {
    CIFs_selected <- clampP(CIFs_selected, prob.bound)
    if (effect_measure == "RR") {
      return(1 / CIFs_selected)
    } else if (effect_measure == "OR") {
      denom <- clampP(1 - CIFs_selected, prob.bound)
      return(1 / CIFs_selected + 1 / denom)
    } else if (effect_measure == "SHR") {
      denom <- clampP(1 - CIFs_selected, prob.bound)
      tmp1 <- -1 / denom
      tmp2 <- log(denom)
      return(tmp1 / tmp2)
    } else {
      stop("Invalid effect_measure. Must be RR, OR or SHR.")
    }
  }

  a11 <- calculateA(estimand$effect.measure1, CIF_sel)
  d_ey_d_beta_11 <- 1 / a11
  d_ey_d_alpha_11 <- 1 / a11

  beta_11 <- if (ncol(x_a) == 0) NULL else sweep(x_a, 1, d_ey_d_beta_11, `*`)
  d_11 <- cbind(sweep(x_l, 1, d_ey_d_alpha_11, `*`), beta_11)
  return(d_11)
}

solveEstimatingEquation <- function(
    nuisance.model,
    exposure,
    strata,
    normalized_data,
    outcome.type,
    estimand,
    optim.method,
    prob.bound,
    alpha_beta_0
) {
  if (outcome.type == "binomial") {
    ip.weight.matrix <- matrix(1,nrow(normalized_data),1)
  } else {
    ip.weight.matrix <- calculateIPCW(nuisance.model, normalized_data, estimand$code.censoring, strata, estimand$time.point)
  }

  makeObjectiveFunction <- function() {
    out_ipcw <- list()
    initial.CIFs <- NULL

    call_and_capture <- function(fun, ...) {
      out_ipcw <<- do.call(fun, list(...))
      out_ipcw$ret
    }

    estimating_equation_i <- function(p) call_and_capture(
      estimating_equation_ipcw,
      formula = nuisance.model, data = normalized_data, exposure = exposure, outcome.type=outcome.type,
      ip.weight.matrix = ip.weight.matrix, alpha_beta = p, estimand = estimand,
      optim.method = optim.method, prob.bound = prob.bound,
      initial.CIFs = initial.CIFs
    )

    setInitialCIFs <- function(new.CIFs) initial.CIFs <<- new.CIFs
    getResults     <- function() out_ipcw
    list(
      estimating_equation_i  = estimating_equation_i,
      setInitialCIFs = setInitialCIFs,
      getResults     = getResults
    )
  }

  assessConvergence <- function(new_params, current_params, current_obj_value, optim.parameter1, optim.parameter2, optim.parameter3) {
    assessRelativeDifference <- function(new, old) {
      max(abs(new - old) / pmax(1, abs(old)))
    }
    is_stalled <- function(x, stall_patience = 3, stall_eps = 1e-3) {
      n <- length(x)
      if (n < stall_patience) return(FALSE)
      recent <- x[(n - stall_patience + 1L):n]
      rng <- range(recent)
      rel_diff <- (diff(rng) / max(1e-12, mean(recent)))
      rel_diff <= stall_eps
    }
    if (any(abs(new_params) > optim.parameter3)) {
      stop("Estimates are either too large or too small, and convergence might not be achieved.")
    }
    param_diff <- abs(new_params - current_params)
    max.absolute.difference <- max(param_diff)
    relative.difference <- assessRelativeDifference(new_params, current_params)
    obj_value <- drop(crossprod(obj$estimating_equation_i(new_params)))

    converged <- (relative.difference <= optim.parameter1) || (obj_value <= optim.parameter2) || is_stalled(c(current_obj_value, obj_value))
    criteria1 <- (relative.difference <= optim.parameter1)
    criteria2 <- (obj_value <= optim.parameter2)
    criteria3 <- is_stalled(c(current_obj_value, obj_value))
    converged  <- (criteria1 || criteria2 || criteria3)
    converged.by <- if (!converged) NA_character_
    else if (criteria1) "Converged in relative difference"
    else if (criteria2) "Converged in objective function"
    else "Stalled"

    list(converged = converged, converged.by=converged.by, relative.difference = relative.difference, max.absolute.difference = max.absolute.difference, obj_value = obj_value)
  }

  obj <- makeObjectiveFunction()
  iteration <- 0L
  max.absolute.difference <- Inf
  out_nleqslv <- NULL
  current_params <- alpha_beta_0
  current_obj_value <- numeric(0)
  trace_df  <- NULL

  computation.time0 <- proc.time()

  while ((iteration < optim.method$optim.parameter4) & (max.absolute.difference > optim.method$optim.parameter1)) {
    iteration <- iteration + 1L
    prev_params <- current_params

    out_nleqslv <- nleqslv(
      prev_params,
      obj$estimating_equation_i,
      method  = reg_choose_nleqslv_method(optim.method$nleqslv.method),
      control = list(maxit = optim.method$optim.parameter5, allowSingular = FALSE)
    )
    new_params <- out_nleqslv$x

    current_obj_value <- drop(crossprod(obj$estimating_equation_i(new_params)))
    obj$setInitialCIFs(obj$getResults()$potential.CIFs)

    ac <- assessConvergence(
      new_params, prev_params, current_obj_value,
      optim.method$optim.parameter1, optim.method$optim.parameter2, optim.method$optim.parameter3
    )
    current_params <- new_params
    current_obj_value <- ac$obj_value
    max.absolute.difference <- ac$max.absolute.difference
    if (ac$converged) break
  }
  if (!length(current_obj_value) || !is.finite(current_obj_value) ||
      current_obj_value > optim.method$optim.parameter2) {
    return(rep(NA_real_, length(current_params)))
  }
  return(current_params)
}

#' @name cov_wild_bootstrap
#' @noRd
#' @param influence.function  n x p numeric matrix of individual IFs (rows = units, cols = parameters).
#' @param R                   Integer, number of bootstrap replications (default 1000).
#' @param multiplier          "rademacher", "mammen", or "gaussian" (default "rademacher").
#' @param center              Logical, center IFs by column before bootstrapping (default FALSE).
#' @param seed                Optional integer seed for reproducibility (default NULL).
#' @return                    p x p covariance matrix estimated via wild bootstrap.
#'
#' @examples
#' # V <- cov_wild_bootstrap(IF, R = 1000, multiplier = "rademacher")
cov_wild_bootstrap <- function(
    influence.function,
    R = 1000,
    multiplier = c("rademacher", "mammen", "gaussian"),
    center = FALSE,
    seed = NULL
) {
  stopifnot(is.matrix(influence.function), is.numeric(influence.function))
  multiplier <- match.arg(multiplier)

  n <- nrow(influence.function)
  p <- ncol(influence.function)
  if (n <= 1L || p < 1L) stop("`influence.function` must be n x p with n>1 and p>=1.")
  if (!is.null(seed)) set.seed(seed)

  IF <- influence.function
  if (isTRUE(center)) {
    IF <- scale(IF, center = TRUE, scale = FALSE)
  }
  draw_multipliers <- function(n, kind) {
    if (kind == "rademacher") {
      sample(c(-1, 1), size = n, replace = TRUE)
    } else if (kind == "mammen") {
      sq5 <- sqrt(5)
      a <- (1 - sq5) / 2
      b <- (1 + sq5) / 2
      pa <- (sq5 + 1) / (2 * sq5)
      u <- stats::runif(n)
      ifelse(u < pa, a, b)
    } else if (kind == "gaussian") {
      stats::rnorm(n)
    } else {
      stop("Unknown `multiplier`.")
    }
  }
  IFt <- t(IF)
  deltas <- matrix(NA_real_, nrow = p, ncol = R)
  for (b in seq_len(R)) {
    g <- draw_multipliers(n, multiplier)
    deltas[, b] <- (IFt %*% g) / n
  }
  return(stats::cov(t(deltas)))
}
