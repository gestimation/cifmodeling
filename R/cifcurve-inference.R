# Fixed-weight inference only. Scores use hazard logits with the no-event
# category as baseline, and no time/cause-specific equation rescaling.

# Validate the public inference combination and select the shared CI path.
.cifcurve_check_inference <- function(error, conf.dist, df.method, engine, level) {
  student <- conf.dist == "t"
  adjusted <- error == "fay" || student
  shared <- adjusted || error == "if"
  if (shared && (length(level) != 1L || !is.numeric(level) ||
                 !is.finite(level) || level <= 0 || level >= 1))
    stop("'conf.int' must be a finite number in (0, 1).", call. = FALSE)
  if (student && !error %in% c("if", "fay"))
    stop("Student intervals require error = 'if' or 'fay'.", call. = FALSE)
  if (student && error == "fay" && df.method == "equal")
    stop("Equal allocation is available for IF-t (error = 'if') only.", call. = FALSE)
  if (adjusted && identical(engine, "calculateKM"))
    stop("Student and Fay--Graubard intervals require engine = 'calculateAJ_Rcpp' or 'auto'.", call. = FALSE)
  shared
}

.cifcurve_df_system <- function(time, event, weight, times, corrected) {
  keep <- weight > 0
  time <- time[keep]; event <- event[keep]; weight <- weight[keep]
  n <- length(time)
  if (!n) return(list(n = 0L))
  weight <- weight / max(weight)
  group <- findInterval(time, times)
  causes <- if (any(event >= 2L)) 1:2 else 1L
  knots <- sort(unique(group[event > 0L]))
  blocks <- list(); pos <- 0L
  for (j in knots) {
    risk <- group >= j
    Y <- sum(weight[risk])
    h <- vapply(causes, function(k) sum(weight[group == j & event == k]) / Y, numeric(1))
    complete <- all(group[risk] == j & event[risk] > 0L)
    if (complete) h <- h / sum(h)
    cats <- which(h > 0)
    coords <- if (complete) cats[-length(cats)] else cats
    r <- length(coords); ix <- if (r) pos + seq_len(r) else integer()
    V <- diag(h[coords], r) - tcrossprod(h[coords])
    Vinv <- if (r && is.finite(rcond(V)) && rcond(V) > 1e-12) solve(V) else NULL
    U <- matrix(0, n, r); dh <- matrix(0, length(causes), r)
    if (r) for (k in seq_len(r)) {
      U[, k] <- weight * ((group == j & event == causes[coords[k]]) - risk * h[coords[k]])
      dh[, k] <- h * ((seq_along(causes) == coords[k]) - h[coords[k]])
    }
    blocks[[length(blocks) + 1L]] <- list(
      group = j, Y = Y, h = h, stay = if (complete) 0 else 1 - sum(h),
      complete = complete, ix = ix, V = V, Vinv = Vinv, U = U,
      leverage = weight * risk / Y, dh = dh
    )
    pos <- pos + r
  }
  U <- matrix(0, n, pos); leverage <- U
  for (b in blocks) if (length(b$ix)) {
    U[, b$ix] <- b$U
    leverage[, b$ix] <- b$leverage
  }
  H <- if (corrected) 1 / sqrt(1 - pmin(leverage, .75)) else matrix(1, n, pos)
  list(n = n, P = pos, causes = causes, blocks = blocks, U = U, H = H)
}

# The audited Gram trace calculation, specialized to a block-diagonal
# fixed-weight bread. No PS machinery and no delete-one refits are involved.
.cifcurve_df_coordinate <- function(sys, C, blocks, allocation) {
  fail <- function(reason) list(df = NA_real_, reason = reason)
  n <- sys$n; P <- length(C)
  if (!P || !any(C != 0)) return(fail("zero_contrast_information"))
  if (any(!is.finite(C))) return(fail("nonfinite_contrast"))
  C <- C / max(abs(C))
  v <- numeric(P); wg <- numeric(n)
  for (b in blocks) if (length(b$ix) && any(C[b$ix] != 0)) {
    if (is.null(b$Vinv)) return(fail("ill_conditioned_hazard_bread"))
    v[b$ix] <- drop(b$Vinv %*% C[b$ix]) / b$Y
    if (allocation == "deletion") {
      mass <- sum(C[b$ix] * v[b$ix])
      if (!is.finite(mass) || mass < 0) return(fail("invalid_deletion_weight"))
      if (mass > 0) {
        if (any(b$leverage >= 1)) return(fail("singular_hazard_deletion_bread"))
        wg <- wg + mass * b$leverage / (1 - b$leverage)
      }
    }
  }
  if (allocation == "equal") alpha <- rep(1 / n, n) else {
    if (any(!is.finite(wg)) || !any(wg > 0)) return(fail("invalid_allocation"))
    alpha <- wg / max(wg); alpha <- alpha / sum(alpha)
  }
  H <- sys$H[, seq_len(P), drop = FALSE]
  T <- sys$U[, seq_len(P), drop = FALSE] * H
  W <- sweep(H, 2L, v, `*`)
  B <- tcrossprod(W, T)
  E <- matrix(0, n, P); R <- matrix(0, P, n)
  for (b in blocks) if (length(b$ix) && any(v[b$ix] != 0)) {
    E[, b$ix] <- (W[, b$ix, drop = FALSE] %*% b$V) * (b$Y * b$leverage)
    R[b$ix, ] <- b$Vinv %*% t(T[, b$ix, drop = FALSE]) / b$Y
  }
  D <- E %*% R
  scale <- max(abs(B), abs(D))
  if (!is.finite(scale) || scale <= 0) return(fail("zero_working_covariance"))
  B <- B / scale; D <- D / scale
  AB <- B * alpha
  dense <- function() {
    Q <- tcrossprod(D) - tcrossprod(AB, D) - tcrossprod(D, AB)
    diag(Q) <- diag(Q) + alpha * rowSums(B^2)
    c(sum(diag(Q)), sum(Q^2))
  }
  di <- which(colSums(abs(E)) > 0 & rowSums(abs(R)) > 0)
  Rd <- R[di, , drop = FALSE]; Ed <- E[, di, drop = FALSE] / scale
  if (all(H == 1)) {
    Z <- cbind(Ed, alpha / scale)
    cross <- -Rd %*% drop(T %*% v)
  } else {
    ai <- which(colSums(abs(W)) > 0 & colSums(abs(T)) > 0)
    Z <- cbind(Ed, W[, ai, drop = FALSE] * alpha / scale)
    cross <- -Rd %*% T[, ai, drop = FALSE]
  }
  if (ncol(Z) >= n) tt <- dense() else {
    d <- length(di); ix <- seq_len(d); jx <- d + seq_len(ncol(Z) - d)
    M <- matrix(0, ncol(Z), ncol(Z))
    M[ix, ix] <- tcrossprod(Rd); M[ix, jx] <- cross; M[jx, ix] <- t(cross)
    largest <- apply(abs(Z), 2L, max); keep <- largest > 0
    Z <- sweep(Z[, keep, drop = FALSE], 2L, largest[keep], `/`)
    norm <- sqrt(colSums(Z^2)); Z <- sweep(Z, 2L, norm, `/`)
    M <- M[keep, keep, drop = FALSE] * tcrossprod(largest[keep] * norm)
    G <- crossprod(Z); MG <- M %*% G
    qdiag <- alpha * rowSums(B^2); adiag <- rowSums((Z %*% M) * Z)
    tt <- c(sum(M * G) + sum(qdiag),
            sum(MG * t(MG)) + 2 * sum(qdiag * adiag) + sum(qdiag^2))
    size <- c(sum(abs(M * G)) + sum(abs(qdiag)),
              sum(abs(MG * t(MG))) + 2 * sum(abs(qdiag * adiag)) + sum(qdiag^2))
    if (any(!is.finite(c(tt, size))) || any(abs(tt) <= 1e-8 * size)) tt <- dense()
  }
  if (any(!is.finite(tt)) || any(tt <= 0)) return(fail("invalid_trace_ratio"))
  nu <- tt[1]^2 / tt[2]
  if (!is.finite(nu) || nu < 1 - 1e-7 || nu > n + 1e-7)
    return(fail("invalid_effective_rank"))
  list(df = max(1, nu), reason = "ok")
}

.cifcurve_df_stratum <- function(time, event, weight, times, corrected, allocation) {
  sys <- .cifcurve_df_system(time, event, weight, times, corrected)
  df <- rep(NA_real_, length(times)); reason <- rep("no_positive_weights", length(times))
  if (!sys$n) return(list(df = df, reason = reason))
  G <- matrix(0, sys$P, length(sys$causes) + 1L)
  S <- 1; used <- 0L; last <- list(df = NA_real_, reason = "zero_contrast_information")
  for (j in seq_along(times)) {
    changed <- FALSE
    while (used < length(sys$blocks) && sys$blocks[[used + 1L]]$group <= j) {
      used <- used + 1L; b <- sys$blocks[[used]]; old <- G[, 1L]
      for (k in seq_along(sys$causes)) {
        G[, k + 1L] <- G[, k + 1L] + b$h[k] * old
        if (length(b$ix)) G[b$ix, k + 1L] <- G[b$ix, k + 1L] + S * b$dh[k, ]
      }
      G[, 1L] <- b$stay * old
      if (b$complete) G[, 1L] <- 0 else if (length(b$ix))
        G[b$ix, 1L] <- G[b$ix, 1L] - S * colSums(b$dh)
      S <- S * b$stay; changed <- TRUE
    }
    if (changed) {
      P <- sum(vapply(sys$blocks[seq_len(used)], function(b) length(b$ix), integer(1)))
      last <- .cifcurve_df_coordinate(sys, G[seq_len(P), 2L], sys$blocks[seq_len(used)], allocation)
    }
    df[j] <- last$df; reason[j] <- last$reason
  }
  list(df = df, reason = reason)
}

# Preserve the diagnostic precedence shared by normal and Student intervals.
.cifcurve_ci_validity <- function(p, se, critical, df.reason, supported) {
  status <- rep("ok", length(p))
  valid.p <- is.finite(p) & p >= 0 & p <= 1
  valid.v <- is.finite(se) & se >= 0
  boundary <- valid.p & (p == 0 | p == 1)
  status[!valid.p] <- "invalid_estimate"
  status[valid.p & !valid.v] <- "invalid_variance"
  status[valid.p & valid.v & !boundary & se == 0] <- "zero_interior_variance"
  bad.df <- valid.p & valid.v & !boundary & se > 0 & !is.finite(critical)
  status[bad.df] <- df.reason[bad.df]
  status[valid.p & valid.v & boundary] <- "endpoint_degenerate"
  status[!supported] <- "outside_observed_followup"
  list(status = status, boundary = boundary,
       valid = status %in% c("ok", "endpoint_degenerate"))
}

.cifcurve_adjusted_ci <- function(out, time, event, weight, strata,
                                 conf.type, level, conf.dist, df.method) {
  nr <- length(out$time); df <- rep(Inf, nr); df.reason <- rep("normal", nr)
  corrected <- out$error == "fay"
  supported <- logical(nr)
  counts <- if (length(out$strata)) as.integer(out$strata) else nr
  groups <- sort(unique(strata)); offset <- 0L
  for (g in seq_along(groups)) {
    ii <- which(strata == groups[g]); jj <- offset + seq_len(counts[g]); offset <- offset + counts[g]
    positive <- ii[weight[ii] > 0]
    if (length(positive)) supported[jj] <- out$time[jj] <= max(time[positive])
    if (conf.dist == "t") {
      d <- .cifcurve_df_stratum(time[ii], event[ii], weight[ii], out$time[jj], corrected, df.method)
      df[jj] <- d$df; df.reason[jj] <- d$reason
    }
  }
  critical <- rep(stats::qnorm((1 + level) / 2), nr)
  if (conf.dist == "t") {
    critical[] <- NA_real_; ok <- is.finite(df) & df > 0
    critical[ok] <- stats::qt((1 + level) / 2, df[ok])
  }
  p <- out$surv; se <- out$std.err
  validity <- .cifcurve_ci_validity(p, se, critical, df.reason, supported)
  none <- is.null(conf.type) || tolower(as.character(conf.type)[1L]) %in% c("none", "n")
  out$conf.type <- if (none) "none" else .cifcurve_normalize_conf_type(conf.type)
  if (none) { out$lower <- out$upper <- out$low <- out$high <- NULL } else {
    safe <- list(surv = ifelse(validity$valid, p, NA_real_),
                 std.err = ifelse(validity$valid, se, NA_real_))
    limits <- calculateCI(safe, level, out$conf.type, critical_value = critical)
    endpoint <- validity$boundary & validity$valid
    limits$lower[endpoint] <- limits$upper[endpoint] <- p[endpoint]
    out$lower <- out$low <- limits$lower; out$upper <- out$high <- limits$upper
  }
  out$conf.dist <- conf.dist; out$conf.int <- level
  out$df.method <- if (conf.dist == "t") df.method else NULL
  out$df <- df; out$df.reason <- df.reason; out$critical.value <- critical
  out$conf.status <- validity$status
  if (corrected) out$fay.bound <- .75
  out
}
