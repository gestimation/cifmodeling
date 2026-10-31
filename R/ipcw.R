# Shared observation weights. The caller specifies observability, evaluation
# times and the limit consistent with its estimating equation.
calculate_ipcw_weights <- function(censoring, time, observed,
                                   strata = rep.int("pooled", length(time)),
                                   side = "left") {
  if (!is.logical(observed) || anyNA(observed) || length(observed) != length(time))
    stop("`observed` must have one nonmissing logical value per evaluation time.", call. = FALSE)
  if (length(strata) != length(time))
    stop("IPCW times and strata must have the same length.", call. = FALSE)
  raw <- predict_censoring_km(censoring, time, strata, side = side, use = "raw")
  used <- predict_censoring_km(censoring, time, strata, side = side, use = "used")
  if (any(used[observed] <= 0))
    stop("IPCW positivity failed: an observed outcome has zero censoring survival.", call. = FALSE)
  weights <- numeric(length(time))
  weights[observed] <- 1 / used[observed]
  if (any(!is.finite(weights))) stop("IPCW weights must be finite.", call. = FALSE)
  n.positive <- sum(weights > 0)
  list(weights = weights, survival.raw = raw, survival.used = used,
    observed = observed, evaluation.time = time, side = side,
    diagnostics = list(n.positive = n.positive, n.zero = sum(!observed),
      maximum.weight = if (n.positive) max(weights) else NA_real_,
      effective.n = if (n.positive) sum(weights)^2 / sum(weights^2) else 0,
      minimum.survival = if (n.positive) min(raw[observed]) else NA_real_,
      truncated = sum(raw[observed] != used[observed]), normalized = FALSE))
}
