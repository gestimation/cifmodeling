calculatePotentialCIFs <- function(alpha_beta_tmp, x_a, x_l, offset, epsilon, estimand, optim.method, prob.bound, initial.CIFs = NULL) {
  K <- estimand$exposure.levels
  n  <- length(epsilon)
  iv <- estimand$index.vector

  alpha_1    <- alpha_beta_tmp[seq_len(iv[1])]
  beta_tmp_1 <- alpha_beta_tmp[seq.int(iv[2], iv[3])]
  alpha_2    <- alpha_beta_tmp[seq.int(iv[4], iv[5])]
  beta_tmp_2 <- alpha_beta_tmp[seq.int(iv[6], iv[7])]
  alpha_tmp_1 <- as.numeric(x_l %*% matrix(alpha_1, ncol = 1) + offset)
  alpha_tmp_2 <- as.numeric(x_l %*% matrix(alpha_2, ncol = 1) + offset)

  tmp1 <- sum(epsilon == estimand$code.event1) / n + prob.bound
  CIF1 <- rep(tmp1, K)
  tmp2 <- sum(epsilon == estimand$code.event2) / n + prob.bound
  CIF2 <- rep(tmp2, K)
  p0 <- c(CIF1,CIF2)
  p0     <- clampP(p0, prob.bound)
  log_p0 <- log(p0)

  keys <- apply(x_l, 1, function(r) paste0(r, collapse = "\r"))
  uniq <- match(keys, unique(keys))
  cache_log_CIFs <- vector("list", length = max(uniq))

  potential.CIFs <- matrix(NA_real_, nrow = n, ncol = 2*K)
  for (i in seq_len(nrow(x_l))) {
    k <- uniq[i]

    if (!is.null(cache_log_CIFs[[k]])) {
      log_CIFs <- cache_log_CIFs[[k]]
    } else {
      CIFs0 <- if (!is.null(initial.CIFs)) as.numeric(initial.CIFs[i, seq_len(2*K), drop = FALSE]) else exp(log_p0)
      log_CIFs0 <- log(clampP(CIFs0, prob.bound))
      log_CIFs <- callLevenbergMarquardt(
        log_CIFs0  = log_CIFs0,
        alpha1    = alpha_tmp_1[i],
        beta1     = beta_tmp_1,
        alpha2    = alpha_tmp_2[i],
        beta2     = beta_tmp_2,
        estimand  = estimand,
        optim.method = optim.method,
        prob.bound = prob.bound
      )
      cache_log_CIFs[[k]] <- log_CIFs
    }
    potential.CIFs[i, ] <- clampP(exp(log_CIFs), prob.bound)
  }
  return(potential.CIFs)
}

calculatePotentialRisk <- function(alpha_beta, x_a, x_l, offset, estimand) {
  iv <- estimand$index.vector
  one <- rep(1, nrow(x_l))
  alpha_beta_ <- as.matrix(as.vector(alpha_beta))
  if (ncol(alpha_beta_) == 1L) alpha_beta_ <- t(alpha_beta_)
  phi <- x_l %*% alpha_beta_[, seq_len(iv[1])] + offset
  theta <- one * alpha_beta_[, iv[2]]

  if (estimand$effect.measure1 == "RR") {
    expphi <- exp(phi)
    exptheta <- exp(theta)
    tol <- 1e-8
    if (all(abs(phi) < tol)) {
      p_10 <- one / (one + exptheta)
      p_11 <- exptheta * p_10
    } else {
      denomi_1 <- -(exptheta + one) * expphi
      denomi_2 <- sqrt(exp(2 * phi) * (exptheta + one) * (exptheta + one) + 4 * exp(theta + phi) * (one - expphi))
      denomi <- denomi_1 + denomi_2
      numera <- 2 * exptheta * (one - expphi)
      p_10 <- denomi / numera
      p_11 <- exptheta * p_10
    }
  } else if (estimand$effect.measure1 == "OR") {
    sqrt1 <- sqrt(exp(-theta - phi))
    sqrt2 <- sqrt(exp(theta - phi))
    if (all(phi == theta)) {
      p_10 <- 0.5 * one
      p_11 <- one/(one + sqrt1)
    } else {
      p_10 <- one/(one + sqrt2)
      p_11 <- one/(one + sqrt1)
    }
  } else {
    stop("Invalid effect_measure. Must be RR or OR.")
  }
  return(cbind(p_10, p_11))
}

residuals_CIFs <- function(log_CIFs, alpha1, beta1, alpha2, beta2, estimand, prob.bound) {
  K <- estimand$exposure.levels
  clog_CIFs <- clampLogP(as.numeric(log_CIFs))
  if (length(clog_CIFs) != 2*K || length(clog_CIFs) %% 2 != 0) stop("log_CIFs length must be even (2K).")

  i1 <- 1:K                      # indices for event1 block
  i2 <- (K+1):(2*K)              # indices for event2 block
  p1 <- exp(clog_CIFs[i1])
  p2 <- exp(clog_CIFs[i2])
  p0 <- pmax(1 - p1 - p2, prob.bound)
  logp0 <- log(p0)

  if (K < 2) stop("exposure.levels must be at least 2.")

  r <- numeric(4 * (K - 1))
  for (k in 2:K) {
    idx <- k - 1
    r[idx] <- alpha1 - clog_CIFs[i1[1]] - clog_CIFs[i1[k]] + logp0[1] + logp0[k]
    r[(K - 1) + idx] <- alpha2 - clog_CIFs[i2[1]] - clog_CIFs[i2[k]] + logp0[1] + logp0[k]
  }
  calculateResidualBeta <- function(effect.measure, pref, pL, clogpref, clogpL, beta) {
    if (length(beta) != length(pL) || length(beta) != length(clogpL)) {
      stop("beta, pL, and clogpL must have the same length.")
    }
    if (effect.measure == "RR") {
      beta - clogpL + clogpref
    } else if (effect.measure == "OR") {
      beta - (log(pL) - log1p(-pL)) + (log(pref) - log1p(-pref))
    } else if (effect.measure == "SHR") {
      Bi <- log1p(-pref)
      Bj <- log1p(-pL)
      exp(beta) - (Bj / Bi)
    } else stop("effect.measure must be RR, OR, or SHR.")
  }
  beta_idx1 <- (2 * (K - 1) + 1):(3 * (K - 1))
  beta_idx2 <- (3 * (K - 1) + 1):(4 * (K - 1))
  r[beta_idx1] <- calculateResidualBeta(estimand$effect.measure1, p1[1], p1[2:K], clog_CIFs[i1[1]], clog_CIFs[i1[2:K]], beta1)
  r[beta_idx2] <- calculateResidualBeta(estimand$effect.measure2, p2[1], p2[2:K], clog_CIFs[i2[1]], clog_CIFs[i2[2:K]], beta2)
  return(r)
}

jacobian_CIFs <- function(log_CIFs, alpha1, beta1, alpha2, beta2, estimand, prob.bound) {
  K <- estimand$exposure.levels
  clog_CIFs <- clampLogP(as.numeric(log_CIFs))
  i1 <- 1:K
  i2 <- (K+1):(2*K)
  p1 <- exp(clog_CIFs[i1])
  p2 <- exp(clog_CIFs[i2])
  p0 <- pmax(1 - p1 - p2, prob.bound)

  dlogrem_dlp1 <- -p1/p0
  dlogrem_dlp2 <- -p2/p0

  if (K < 2) stop("exposure.levels must be at least 2.")

  J <- matrix(0.0, nrow = 4 * (K - 1), ncol = 2*K)
  for (k in 2:K) {
    r1 <- k - 1
    J[r1,   i1[1]] <- -1 + dlogrem_dlp1[1]
    J[r1,   i2[1]] <-      dlogrem_dlp2[1]
    J[r1,   i1[k]] <- -1 + dlogrem_dlp1[k]
    J[r1,   i2[k]] <-      dlogrem_dlp2[k]

    r2 <- (K - 1) + (k - 1)
    J[r2,   i2[1]] <- -1 + dlogrem_dlp2[1]
    J[r2,   i1[1]] <-      dlogrem_dlp1[1]
    J[r2,   i2[k]] <- -1 + dlogrem_dlp2[k]
    J[r2,   i1[k]] <-      dlogrem_dlp1[k]
  }

  row_beta1_start <- 2 * (K - 1) + 1
  row_beta2_start <- 3 * (K - 1) + 1
  idx_L1 <- i1[2:K]
  idx_L2 <- i2[2:K]
  p1_L <- p1[2:K]
  p2_L <- p2[2:K]

  for (k in seq_along(idx_L1)) {
    row <- row_beta1_start + k - 1
    idx_L <- idx_L1[k]
    if (estimand$effect.measure1 == "RR") {
      J[row, i1[1]] <- 1
      J[row, idx_L] <- -1
    } else if (estimand$effect.measure1 == "OR") {
      pl <- p1_L[k]
      J[row, i1[1]] <- 1 + p1[1] / (1 - p1[1])
      J[row, idx_L] <- -1 - pl / (1 - pl)
    } else if (estimand$effect.measure1 == "SHR") {
      Bi <- log1p(-p1[1])
      Bj <- log1p(-p1_L[k])
      pl <- p1_L[k]
      J[row, i1[1]] <- -Bj * p1[1] / ((1 - p1[1]) * Bi^2)
      J[row, idx_L] <- pl / ((1 - pl) * Bi)
    } else stop("effect.measure must be RR, OR, or SHR.")
  }

  for (k in seq_along(idx_L2)) {
    row <- row_beta2_start + k - 1
    idx_L <- idx_L2[k]
    if (estimand$effect.measure2 == "RR") {
      J[row, i2[1]] <- 1
      J[row, idx_L] <- -1
    } else if (estimand$effect.measure2 == "OR") {
      pl <- p2_L[k]
      J[row, i2[1]] <- 1 + p2[1] / (1 - p2[1])
      J[row, idx_L] <- -1 - pl / (1 - pl)
    } else if (estimand$effect.measure2 == "SHR") {
      Bi <- log1p(-p2[1])
      Bj <- log1p(-p2_L[k])
      pl <- p2_L[k]
      J[row, i2[1]] <- -Bj * p2[1] / ((1 - p2[1]) * Bi^2)
      J[row, idx_L] <- pl / ((1 - pl) * Bi)
    } else stop("effect.measure must be RR, OR, or SHR.")
  }
  return(J)
}

LevenbergMarquardt <- function(start,
                               res_fun,
                               jac_fun,
                               maxit = 30,
                               ftol  = 1e-10,
                               ptol  = 1e-10,
                               lambda0    = 1e-3,
                               lambda_max = 10,
                               lambda_min = 0.1,
                               lambda_increment    = 2.0,
                               lambda_decrement    = 0.5,
                               do_linesearch = TRUE,
                               ls_c      = 1e-4,
                               ls_shrink = 0.5,
                               verbose   = FALSE) {
  lp   <- as.numeric(start)
  r    <- res_fun(lp)
  f2   <- drop(crossprod(r))
  J    <- jac_fun(lp)
  g    <- drop(crossprod(J, r))
  lam  <- lambda0
  prev_f2 <- Inf

  if (isTRUE(verbose)) {
    hist <- list(it=integer(), f2=double(),
                 grad_inf=double(), lambda=double(),
                 step_norm=double(), rho=double())
  }

  for (k in seq_len(maxit)) {
    grad_inf <- max(abs(g))
    if (is.finite(grad_inf) && grad_inf < 1e-6) break
    A <- crossprod(J)
    diag(A) <- diag(A) + lam
    R <- try(chol(A), silent = TRUE)
    if (inherits(R, "try-error")) {
      lam <- min(lam * 10, lambda_max)
      next
    }
    delta <- -backsolve(R, forwardsolve(t(R), g))
    step_norm <- max(abs(delta))

    if (is.finite(step_norm) &&
        step_norm <= ptol * (ptol + max(1, max(abs(lp))))) break

    pred <- 0.5 * sum(delta * (lam * delta - g))
    if (!is.finite(pred) || pred <= 0) pred <- .Machine$double.eps

    lp_try <- lp + delta
    r_try  <- res_fun(lp_try)
    f2_try <- drop(crossprod(r_try))

    rho <- (f2 - f2_try) / pred
    accepted <- FALSE

    if (is.finite(rho) && rho > 0) {
      lp <- lp_try; r <- r_try; f2 <- f2_try
      lam <- max(lambda_min, lam * max(lambda_decrement, 1/(1 + rho)))
      accepted <- TRUE
    } else {
      lam <- min(lambda_max, lam * (lambda_increment * (1 + abs(ifelse(is.finite(rho), rho, 0)))))
      if (do_linesearch) {
        t <- 1.0; gTd <- sum(g * delta)
        while (t > 1e-6) {
          lp_ls <- lp + t * delta
          r_ls  <- res_fun(lp_ls)
          f2_ls <- drop(crossprod(r_ls))
          if (is.finite(f2_ls) &&
              f2_ls <= f2 + ls_c * t * gTd) {
            lp <- lp_ls; r <- r_ls; f2 <- f2_ls
            accepted <- TRUE
            break
          }
          t <- t * ls_shrink
        }
      }
    }
    J <- jac_fun(lp)
    g <- drop(crossprod(J, r))

    if (k > 1L &&
        is.finite(prev_f2) && is.finite(f2) &&
        abs(prev_f2 - f2) <= ftol * (abs(f2) + ftol)) break

    if (isTRUE(verbose)) { # ログの記録
      hist$it        <- c(hist$it, k)
      hist$f2        <- c(hist$f2, f2)
      hist$grad_inf  <- c(hist$grad_inf, grad_inf)
      hist$lambda    <- c(hist$lambda, lam)
      hist$step_norm <- c(hist$step_norm, step_norm)
      hist$rho       <- c(hist$rho, if (exists("rho")) rho else NA_real_)
    }
    prev_f2 <- f2
  }
  if (isTRUE(verbose)) attr(lp, "history") <- hist
  return(lp)
}

callLevenbergMarquardt <- function(log_CIFs0, alpha1, beta1, alpha2, beta2, optim.method, estimand, prob.bound)
{
  res_fun <- function(lp) residuals_CIFs(lp, alpha1, beta1, alpha2, beta2, estimand, prob.bound)
  jac_fun <- function(lp) jacobian_CIFs(lp, alpha1, beta1, alpha2, beta2, estimand, prob.bound)
  out_LevenbergMarquardt <- LevenbergMarquardt(
    start         = log_CIFs0,
    res_fun       = res_fun,
    jac_fun       = jac_fun,
    maxit         = optim.method$optim.parameter6   %||% 30,
    ftol          = optim.method$optim.parameter7   %||% 1e-10,
    ptol          = optim.method$optim.parameter8   %||% 1e-10,
    lambda0       = optim.method$optim.parameter9   %||% 1e-3,
    lambda_max    = optim.method$optim.parameter10  %||% 10,
    lambda_min    = optim.method$optim.parameter11  %||% 0.1,
    lambda_increment  = optim.method$optim.parameter12  %||% 2.0,
    lambda_decrement  = optim.method$optim.parameter13  %||% 0.5,
    do_linesearch = TRUE,
    ls_c          = 1e-4,
    ls_shrink     = 0.5,
    verbose       = FALSE
  )
  return(out_LevenbergMarquardt)
}

