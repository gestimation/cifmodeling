# Hudgens, Satten and Longini (2001), doi:10.1111/j.0006-341X.2001.00074.x,
# without truncation. Inputs follow calculateAJ_Rcpp(): weights and strata
# are optional, and epsilon uses 0 for censoring and 1/2 for the two causes.
# Event observations are (left, right] or {left} when left == right.
# Right-censored observations use (left, Inf]. The visit/censoring process
# is assumed ignorable and the first event's cause is known.
# Zero-weight rows are excluded before constructing the intervals.
#
# Output rows are grouped by stratum and ordered by time. With times = NULL,
# each stratum uses its distinct finite observation endpoints. CIF matrices
# have columns "1" and "2". cif is NA when the likelihood leaves a range of
# values. cif.lower and cif.upper are pointwise NPMLE bounds, not confidence
# limits. support contains fitted masses (cause NA denotes the unresolved
# tail). strata gives the number of output rows per stratum, while n counts
# the positive-weight observations. diagnostics describes maximization.
calculateCIF_IC <- function(left, right, epsilon, w = NULL, strata = NULL,
                            times = NULL, tol = 1e-9, max.iter = 10000L) {
  n <- length(left)
  if (!n || !is.numeric(left) || any(!is.finite(left) | left < 0) ||
      !is.numeric(right) || length(right) != n ||
      anyNA(right) || any(right < left)) {
    stop("Use finite nonnegative left endpoints and right >= left.", call. = FALSE)
  }
  if ((!is.numeric(epsilon) && !is.logical(epsilon)) || length(epsilon) != n ||
      any(!(epsilon %in% 0:2)) || any((epsilon == 0) != (right == Inf))) {
    stop("Use codes 1 or 2 for finite events and 0 for right censoring (right = Inf).",
         call. = FALSE)
  }
  if (!length(w)) w <- rep(1, n)
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w) | w < 0) ||
      !any(w > 0)) {
    stop("w must contain one finite nonnegative weight per observation and some positive weight.",
         call. = FALSE)
  }
  if (!length(strata)) strata <- rep(1L, n)
  if ((!is.factor(strata) && !is.numeric(strata)) || length(strata) != n ||
      anyNA(strata) ||
      (is.numeric(strata) && any(!is.finite(strata) | strata < 1 | strata != floor(strata)))) {
    stop("strata must contain a factor or positive integer codes.", call. = FALSE)
  }
  if (!is.null(times) &&
      (!is.numeric(times) || !length(times) || any(!is.finite(times) | times < 0))) {
    stop("times must be finite and nonnegative.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      !is.numeric(max.iter) || length(max.iter) != 1L || !is.finite(max.iter) ||
      max.iter < 1 || max.iter > .Machine$integer.max || max.iter != floor(max.iter)) {
    stop("tol must be positive and max.iter must be a positive integer.", call. = FALSE)
  }

  keep <- w > 0
  groups <- split(which(keep), strata[keep], drop = TRUE)
  fits <- lapply(groups, function(i) {
    eval.time <- if (is.null(times)) {
      sort(unique(c(left[i], right[i][is.finite(right[i])])))
    } else {
      sort(unique(times))
    }
    fit <- .cif_ic_fit(left[i], right[i], epsilon[i], w[i], tol, max.iter)
    support <- fit$support
    exact <- support$left == support$right
    # A support interval can contribute to the upper bound at t if it has
    # any admissible time <= t. Its whole mass is forced in at its right end.
    can.occur <- outer(eval.time, support$left, `>`) |
      (outer(eval.time, support$left, `==`) & rep(exact, each = length(eval.time)))
    must.occur <- outer(eval.time, support$right, `>=`)
    cause <- outer(support$cause, 1:2, `==`)
    cause[is.na(cause)] <- FALSE
    cause.mass <- cause * support$mass
    possible.mass <- cause.mass + ifelse(is.na(support$cause), support$mass, 0)
    lower <- must.occur %*% cause.mass
    ambiguity <- (can.occur & !must.occur) %*% possible.mass
    lower <- pmin(lower, 1)
    upper <- pmin(lower + ambiguity, 1)
    colnames(lower) <- colnames(upper) <- c("1", "2")
    cif <- lower
    cif[ambiguity > 0] <- NA_real_
    list(time = eval.time, cif = cif, cif.lower = lower, cif.upper = upper,
         support = support, diagnostics = fit$diagnostics)
  })
  diagnostics <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))
  if (any(!diagnostics$converged)) {
    warning("CIF NPMLE did not converge in strata: ",
            paste(rownames(diagnostics)[!diagnostics$converged], collapse = ", "),
            call. = FALSE)
  }
  list(
    time = unlist(lapply(fits, `[[`, "time"), use.names = FALSE),
    cif = do.call(rbind, lapply(fits, `[[`, "cif")),
    cif.lower = do.call(rbind, lapply(fits, `[[`, "cif.lower")),
    cif.upper = do.call(rbind, lapply(fits, `[[`, "cif.upper")),
    n = lengths(groups),
    strata = lengths(lapply(fits, `[[`, "time")),
    support = lapply(fits, `[[`, "support"),
    diagnostics = diagnostics
  )
}

.cif_ic_fit <- function(left, right, epsilon, w, tol, max.iter) {
  censored <- epsilon == 0
  if (all(censored | left == right)) {
    # In the exact/right-censored case the likelihood factors over event
    # times. Estimate the conditional event probabilities by d_k(t) / Y(t).
    exits <- sort(unique(left))
    index <- match(left, exits)
    counts <- rowsum(cbind(w * censored, w * (epsilon == 1), w * (epsilon == 2)),
                     index)
    risk <- rev(cumsum(rev(rowSums(counts))))
    hazard <- counts[, 2:3, drop = FALSE] / risk
    survival <- cumprod((c(risk[-1L], 0) + counts[, 1L]) / risk)
    mass <- hazard * c(1, utils::head(survival, -1L))
    at <- which(mass > 0, arr.ind = TRUE)
    support <- data.frame(left = exits[at[, 1L]], right = exits[at[, 1L]],
                          cause = at[, 2L], mass = mass[at])
    remainder <- utils::tail(survival, 1L)
    if (remainder > 0) {
      support <- rbind(support, data.frame(left = max(left), right = Inf,
                                           cause = NA_integer_, mass = remainder))
    }
    probability <- survival[index]
    event <- which(!censored)
    probability[event] <- mass[cbind(index[event], epsilon[event])]
    return(list(support = support, diagnostics = data.frame(
      converged = TRUE, iterations = 0L, kkt = 0,
      loglik = sum(w * log(probability))
    )))
  }

  # For each cause, an event ending at each retained right endpoint covers
  # that cell and no later cell. These rows identify the finite cell masses.
  support <- do.call(rbind, lapply(1:2, function(cause) {
    event <- epsilon == cause
    hi <- sort(unique(right[event]))
    if (!length(hi)) return(NULL)
    exact <- hi %in% left[event & left == right]
    starts <- sort(unique(left[event | censored]))
    lo <- c(-Inf, starts)[findInterval(hi, starts, left.open = TRUE) + 1L]
    # At a non-exact right endpoint, every compatible observation contains
    # (lo, hi]. An earlier right endpoint strictly inside it belongs to one
    # additional event, so the later cell is dominated and can be omitted.
    maximal <- exact | c(-Inf, utils::head(hi, -1L)) <= lo
    data.frame(left = ifelse(exact, hi, lo)[maximal], right = hi[maximal],
               cause = cause)
  }))
  last.censor <- max(c(-Inf, left[censored]))
  if (last.censor >= max(right[!censored])) {
    support <- rbind(support, data.frame(left = last.censor, right = Inf,
                                         cause = NA_integer_))
  }
  # Aggregate identical observations before iterating.
  obs <- stats::aggregate(w, list(left = left, right = right, epsilon = epsilon), sum)
  compatible <- .cif_ic_compatible(obs$left, obs$right, obs$epsilon, support)
  weight <- obs$x / sum(obs$x)
  mass <- rep(1 / nrow(support), nrow(support))
  # Turnbull's self-consistency iteration, accelerated by a likelihood-
  # increasing exchange of mass. The simplex optimality test also checks
  # zero-mass intervals, which the exchange can bring back into the fit.
  for (iteration in seq_len(max.iter)) {
    probability <- drop(compatible %*% mass)
    score <- drop(crossprod(compatible, weight / probability))
    gap <- max(score) - 1
    if (gap <= tol || iteration == max.iter) break
    mass <- mass * score
    mass <- .cif_ic_exchange(mass, compatible, weight)
  }
  support$mass <- mass
  list(support = support, diagnostics = data.frame(
    converged = gap <= tol, iterations = iteration, kkt = max(0, gap),
    loglik = sum(obs$x * log(probability))
  ))
}

.cif_ic_exchange <- function(mass, compatible, weight) {
  probability <- drop(compatible %*% mass)
  score <- drop(crossprod(compatible, weight / probability))
  add <- which.max(score)
  active <- which(mass > 0)
  remove <- active[which.min(score[active])]
  if (score[add] <= score[remove]) return(mass)

  direction <- compatible[, add] - compatible[, remove]
  changed <- direction != 0
  weight <- weight[changed]
  probability <- probability[changed]
  direction <- direction[changed]
  limit <- mass[remove]
  slope <- function(step) {
    prob <- probability + step * direction
    if (any(prob <= 0)) return(-Inf)
    sum(weight * direction / prob)
  }
  step <- if (slope(limit) >= 0) {
    limit
  } else {
    stats::uniroot(slope, c(0, limit), tol = .Machine$double.eps)$root
  }
  mass[add] <- mass[add] + step
  mass[remove] <- mass[remove] - step
  mass
}

.cif_ic_compatible <- function(left, right, epsilon, support) {
  event <- outer(epsilon, support$cause, `==`)
  event[is.na(event)] <- FALSE
  exact <- support$left == support$right
  lo <- outer(left, support$left, `<=`) &
    !(outer(left, support$left, `==`) &
        outer(left != right, exact, `&`))
  hi <- outer(right, support$right, `>=`)
  compatible <- lo & hi & (event | epsilon == 0)
  storage.mode(compatible) <- "double"
  compatible
}
