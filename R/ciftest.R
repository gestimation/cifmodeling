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
