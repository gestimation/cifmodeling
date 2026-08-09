# Deliberately slow, test-only implementation of Gray's estimating equations.
#
# This routine scans every observation at every time and shares no numerical
# helpers with R/ciftest.R. It is intended as an independent executable version
# of the equations, not as package code.
gray_reference_slow <- function(time, status, group, rho = 0, gamma = 0) {
  stopifnot(length(time) == length(status), length(status) == length(group))
  group <- droplevels(factor(group))
  group_id <- as.integer(group)
  group_count <- nlevels(group)
  contrast_count <- group_count - 1L
  stopifnot(contrast_count >= 1L, all(status %in% 0:2))

  at_risk <- tabulate(group_id, nbins = group_count)
  overall_survival <- rep(1, group_count)
  cause1_cif <- numeric(group_count)
  pooled_cif <- 0

  u <- numeric(contrast_count)
  sigma <- matrix(0, contrast_count, contrast_count)
  history_coef <- matrix(0, contrast_count, group_count)
  history_cross <- matrix(0, contrast_count, group_count)
  history_var <- numeric(group_count)

  for (current_time in sort(unique(time))) {
    cause1_deaths <- cause2_deaths <- departures <- numeric(group_count)
    for (g in seq_len(group_count)) {
      in_cell <- group_id == g & time == current_time
      departures[g] <- sum(in_cell)
      cause1_deaths[g] <- sum(in_cell & status == 1L)
      cause2_deaths[g] <- sum(in_cell & status == 2L)
    }
    cause1_total <- sum(cause1_deaths)
    cause2_total <- sum(cause2_deaths)

    if (cause1_total + cause2_total > 0) {
      standardized_risk <- ifelse(
        at_risk > 0,
        at_risk / overall_survival,
        0
      )
      pooled_standardized_risk <- sum(standardized_risk)
      next_pooled_cif <- pooled_cif + cause1_total / pooled_standardized_risk
      w <- (1 - pooled_cif)^rho * pooled_cif^gamma

      derivative <- matrix(0, contrast_count, group_count)
      for (a in seq_len(contrast_count)) {
        for (b in seq_len(group_count)) {
          derivative[a, b] <- -w * standardized_risk[a] *
            standardized_risk[b] / pooled_standardized_risk
          if (a == b) {
            derivative[a, b] <- derivative[a, b] + w * standardized_risk[a]
          }
        }
      }

      if (cause1_total > 0) {
        subdistribution_risk <- ifelse(
          at_risk > 0,
          at_risk * (1 - cause1_cif) / overall_survival,
          0
        )
        subdistribution_risk_total <- sum(subdistribution_risk)
        for (a in seq_len(contrast_count)) {
          u[a] <- u[a] + w * (
            cause1_deaths[a] - cause1_total *
              subdistribution_risk[a] / subdistribution_risk_total
          )
        }
        history_coef <- history_coef + derivative * cause1_total /
          (pooled_standardized_risk * (1 - pooled_cif))
      }

      next_survival <- overall_survival
      next_cause1_cif <- cause1_cif
      occupied <- which(at_risk > 0)
      for (g in occupied) {
        next_survival[g] <- overall_survival[g] *
          (1 - (cause1_deaths[g] + cause2_deaths[g]) / at_risk[g])
        next_cause1_cif[g] <- cause1_cif[g] +
          overall_survival[g] * cause1_deaths[g] / at_risk[g]
      }

      if (cause1_total > 0) {
        for (g in occupied) {
          censoring_factor <- 1
          if (next_survival[g] > 0) {
            censoring_factor <- 1 - (1 - next_pooled_cif) / next_survival[g]
          }
          finite_population <- 1
          if (cause1_total > 1) {
            finite_population <- 1 - (cause1_total - 1) /
              (pooled_standardized_risk * overall_survival[g] - 1)
          }
          hazard_variance <- finite_population * overall_survival[g] *
            cause1_total / (pooled_standardized_risk * at_risk[g])
          innovation <- derivative[, g] -
            censoring_factor * history_coef[, g]
          sigma <- sigma + outer(innovation, innovation) * hazard_variance
          history_cross[, g] <- history_cross[, g] + innovation *
            censoring_factor * hazard_variance
          history_var[g] <- history_var[g] +
            censoring_factor^2 * hazard_variance
        }
      }

      for (g in which(cause2_deaths > 0 & next_survival > 0)) {
        censoring_factor <- (1 - next_pooled_cif) / next_survival[g]
        finite_population <- 1
        if (cause2_deaths[g] > 1) {
          finite_population <- 1 - (cause2_deaths[g] - 1) / (at_risk[g] - 1)
        }
        hazard_variance <- finite_population * overall_survival[g]^2 *
          cause2_deaths[g] / at_risk[g]^2
        linked_history <- censoring_factor * history_coef[, g]
        sigma <- sigma + outer(linked_history, linked_history) * hazard_variance
        history_cross[, g] <- history_cross[, g] - linked_history *
          censoring_factor * hazard_variance
        history_var[g] <- history_var[g] +
          censoring_factor^2 * hazard_variance
      }

      pooled_cif <- next_pooled_cif
      overall_survival <- next_survival
      cause1_cif <- next_cause1_cif
    }

    at_risk <- at_risk - departures
  }

  for (a in seq_len(contrast_count)) {
    for (b in seq_len(contrast_count)) {
      for (g in seq_len(group_count)) {
        sigma[a, b] <- sigma[a, b] +
          history_coef[a, g] * history_coef[b, g] * history_var[g] +
          history_coef[a, g] * history_cross[b, g] +
          history_coef[b, g] * history_cross[a, g]
      }
    }
  }

  # Convert first K-1 group coordinates to non-reference (groups 2,...,K).
  conversion <- matrix(0, contrast_count, contrast_count)
  if (contrast_count > 1L) {
    conversion[cbind(seq_len(contrast_count - 1L), 2:contrast_count)] <- 1
  }
  conversion[contrast_count, ] <- -1
  score <- as.vector(conversion %*% u)
  covariance <- conversion %*% sigma %*% t(conversion)
  statistic <- drop(crossprod(score, solve(covariance, score)))

  list(
    score = score,
    var = covariance,
    statistic = statistic,
    p.value = stats::pchisq(statistic, contrast_count, lower.tail = FALSE)
  )
}
