# Plackett (1954), doi:10.1093/biomet/41.3-4.351. Integrate each
# correlation derivative along the straight line from independence. The
# angle substitution cancels the bivariate density's singular factor.
ciftest_normal_intersection <- function(threshold, correlation) {
  dimension <- nrow(correlation)
  tail <- stats::pnorm(threshold, lower.tail = FALSE)
  if (dimension == 2L) {
    rho <- correlation[1L, 2L]
    if (rho >= 1) return(list(probability = tail, error = 0))
    if (rho <= -1) return(list(probability = max(0, 2 * tail - 1), error = 0))
  }
  pairs <- which(upper.tri(correlation), arr.ind = TRUE)
  probability <- tail^dimension
  error <- 0
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i, ]
    rho <- max(-1, min(1, correlation[pair[1L], pair[2L]]))
    if (rho == 0) next
    other <- if (dimension == 3L) correlation[setdiff(seq_len(dimension), pair), pair] else NULL
    integrand <- function(angle) {
      r <- sin(angle)
      plus <- if (rho < 0) cos(angle)^2 / (1 - r) else 1 + r
      density <- exp(-threshold^2 / plus) / (2 * pi)
      if (dimension == 2L) return(density)
      minus <- if (rho > 0) cos(angle)^2 / (1 + r) else 1 - r
      path <- r / rho
      sum_cor <- sum(other) * path
      diff_cor <- diff(other) * path
      mean <- threshold * sum_cor / plus
      variance <- 1 - (sum_cor^2 / plus + diff_cor^2 / minus) / 2
      standardized <- (threshold - mean) / sqrt(pmax(variance, .Machine$double.eps))
      density * stats::pnorm(standardized, lower.tail = FALSE)
    }
    step <- stats::integrate(integrand, min(0, asin(rho)), max(0, asin(rho)),
                             abs.tol = 1e-13, rel.tol = 1e-9, subdivisions = 250L,
                             stop.on.error = FALSE)
    probability <- probability + sign(rho) * step$value
    error <- error + step$abs.error
  }
  list(probability = max(0, min(tail, probability)), error = error)
}

# Successive conditional normal probabilities transform the Gaussian box
# integral to a unit cube (Genz 1992, doi:10.1080/10618600.1992.10477010).
# For a singular covariance, intersect all constraints at the last latent
# variable they involve, so dependent directions still affect the event.
ciftest_normal_box <- function(lower, upper, correlation, normalization = 1) {
  decomposition <- suppressWarnings(chol(correlation, pivot = TRUE, tol = 1e-13))
  rank <- attr(decomposition, "rank")
  pivot <- attr(decomposition, "pivot")
  loadings <- t(decomposition[seq_len(rank), , drop = FALSE])
  dimension <- nrow(loadings)
  loadings[abs(loadings) < 1e-13] <- 0
  last <- max.col((loadings != 0) * rep(seq_len(rank), each = dimension), ties.method = "last")
  rows <- split(seq_len(dimension), factor(last, levels = seq_len(rank)))
  lower <- rep_len(lower, dimension)[pivot]
  upper <- rep_len(upper, dimension)[pivot]
  integrand <- function(points) {
    n <- nrow(points)
    latent <- matrix(0, n, rank - 1L)
    weight <- rep.int(1, n)
    for (k in seq_len(rank)) {
      pick <- rows[[k]]
      coefficient <- loadings[pick, k]
      predecessors <- seq_len(k - 1L)
      location <- latent[, predecessors, drop = FALSE] %*% t(loadings[pick, predecessors, drop = FALSE])
      lo <- rep.int(-Inf, n)
      hi <- rep.int(Inf, n)
      for (i in seq_along(pick)) {
        a <- (lower[pick[i]] - location[, i]) / coefficient[i]
        b <- (upper[pick[i]] - location[, i]) / coefficient[i]
        lo <- pmax(lo, pmin(a, b))
        hi <- pmin(hi, pmax(a, b))
      }
      # Work on the survival-probability scale for positive intervals.
      right <- lo > 0
      start <- stats::pnorm(lo)
      end <- stats::pnorm(hi)
      start[right] <- stats::pnorm(hi[right], lower.tail = FALSE)
      end[right] <- stats::pnorm(lo[right], lower.tail = FALSE)
      width <- pmax(0, end - start)
      weight <- weight * width
      if (k < rank) {
        probability <- start + points[, k] * width
        probability <- pmin(1 - .Machine$double.eps, pmax(.Machine$double.xmin, probability))
        latent[, k] <- stats::qnorm(probability)
        latent[right, k] <- -latent[right, k]
      }
    }
    weight / normalization
  }
  if (rank == 1L) {
    return(list(probability = normalization * integrand(matrix(numeric(), 1L, 0L)), error = 0))
  }
  if (rank == 2L) {
    result <- stats::integrate(function(u) integrand(matrix(u, ncol = 1L)), 0, 1,
                               abs.tol = 1e-10, rel.tol = 1e-9, subdivisions = 250L,
                               stop.on.error = FALSE)
    return(list(probability = normalization * result$value, error = normalization * result$abs.error))
  }
  result <- ciftest_normal_cube(integrand, rank - 1L)
  list(probability = normalization * result$probability, error = normalization * result$error)
}

# Evaluate very small tails directly by partitioning on the first direction
# outside the acceptance box. For a symmetric box the first positive and
# first negative exceedances have equal probability by central symmetry.
ciftest_normal_rare_tail <- function(threshold, correlation, two_sided, upper_tail) {
  dimension <- nrow(correlation)
  p <- upper_tail
  error <- 0
  lower <- if (two_sided) -threshold else -Inf
  for (j in seq.int(2L, dimension)) {
    order <- c(j, seq_len(j - 1L))
    box <- ciftest_normal_box(c(threshold, rep.int(lower, j - 1L)),
                             c(Inf, rep.int(threshold, j - 1L)),
                             correlation[order, order, drop = FALSE], upper_tail)
    p <- p + box$probability
    error <- error + box$error
  }
  multiplier <- if (two_sided) 2 else 1
  list(probability = multiplier * p, error = multiplier * error)
}

# Shifted Halton cubature. The fixed shift stream makes repeated evaluations
# reproducible without changing the caller's random number state. The error
# estimate combines between-shift variation and change with sample size.
ciftest_normal_cube <- function(integrand, dimension) {
  primes <- 2L
  candidate <- 3L
  while (length(primes) < dimension) {
    if (all(candidate %% primes != 0L)) primes <- c(primes, candidate)
    candidate <- candidate + 2L
  }
  radical_inverse <- function(index, base) {
    values <- numeric(length(index))
    place <- 1 / base
    while (any(index > 0L)) {
      values <- values + place * (index %% base)
      index <- index %/% base
      place <- place / base
    }
    values
  }
  shifts <- matrix(0, 12L, dimension)
  state <- 104729
  for (i in seq_along(shifts)) {
    state <- (16807 * state) %% 2147483647
    shifts[i] <- state / 2147483647
  }
  totals <- numeric(nrow(shifts))
  n <- 0L
  previous <- NA_real_
  for (power in 10:17) {
    size <- 2L^power
    index <- seq.int(n + 1L, size)
    points <- vapply(primes, function(base) radical_inverse(index, base), numeric(length(index)))
    for (i in seq_len(nrow(shifts))) {
      shifted <- sweep(points, 2L, shifts[i, ], `+`) %% 1
      totals[i] <- totals[i] + sum(integrand(abs(2 * shifted - 1)))
    }
    means <- totals / size
    probability <- mean(means)
    error <- max(3 * stats::sd(means) / sqrt(length(means)), abs(probability - previous), na.rm = TRUE)
    if (power >= 11L && error <= 5e-7) break
    previous <- probability
    n <- size
  }
  list(probability = probability, error = error)
}
