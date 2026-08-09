#' Simulate reusable competing-risks data
#'
#' @description
#' Generates two-cause competing-risks data from a coherent restricted-CIF
#' construction. Cause 1 follows a possibly time-varying proportional
#' subdistribution-hazard model,
#'
#' \deqn{\lambda_1^{SD}(t \mid A,L) = \lambda_{10}^{SD}(t)
#'       \exp\{\beta_{1L}L + \log(SHR)b(t)A\}.}
#'
#' The raw cause-2 CIF is multiplied by `1 - F1(tau | A, L)`, as in the
#' Scheike/mets restricted construction, so the two terminal cause
#' probabilities cannot sum above one. The result can be passed directly to
#' [cifcurve()], [cifplot()], [cifpanel()], [ciftest()], or [polyreg()].
#'
#' @param n Positive number of subjects.
#' @param effect Treatment-effect shape for cause 1: `"constant"`, `"early"`,
#'   `"late"`, `"fh-early"`, or `"fh-late"`.
#' @param shr Positive subdistribution hazard ratio for `A = 1` when the
#'   effect function equals one. Set `shr = 1` for a null scenario.
#' @param rho1,rho2 Positive scale parameters of the baseline cumulative
#'   subdistribution hazards.
#' @param rate Two positive time-scale parameters for causes 1 and 2.
#' @param beta1.L Log subdistribution hazard ratio for `L` on cause 1.
#' @param beta2 Two log subdistribution hazard coefficients for `A` and `L`
#'   on the raw cause-2 distribution. The default keeps raw competing risk
#'   equal between treatment groups.
#' @param tau Positive data-generating and administrative-censoring horizon.
#' @param analysis.tau Positive analysis horizon no greater than `tau`.
#' @param change.time Optional early/late change time. By default it is the
#'   time where the baseline cause-1 CIF reaches half its value at
#'   `analysis.tau`.
#' @param effect.rho,effect.gamma Non-negative shape parameters for the
#'   `"fh-early"` and `"fh-late"` alternatives.
#' @param grid.width Positive integration-grid width.
#' @param prob.A,prob.L Probabilities used to generate binary `A` and `L` when
#'   these vectors are not supplied. `prob.L` also defines the population
#'   weights in the marginal truth curves.
#' @param censor.rate Non-negative baseline exponential censoring rate.
#' @param censor.log.hazard Two log-hazard coefficients for `A` and `L` in
#'   the censoring distribution.
#' @param A,L Optional binary vectors of length `n`. Missing vectors are
#'   generated independently using `prob.A` and `prob.L`.
#' @param uniforms Optional named list containing any of `cause1`, `cause2`,
#'   `event`, and `censor`, each a length-`n` vector in `[0, 1]`. This supports
#'   common-random-number comparisons; ordinary use can call [set.seed()].
#'
#' @return A data frame with observed `time`, integer `status` (0 censoring,
#'   1 cause 1, 2 cause 2), factor `A` and `L`, latent event/censoring details,
#'   and subject-specific true probabilities. Attribute `"truth"` contains
#'   conditional and marginal CIF curves, the true time-varying SHR, and
#'   time-averaged CIF contrasts. Attribute `"simulation_parameters"` records
#'   the generating parameters.
#'
#' @details
#' For `"early"`, the SHR equals `shr` before `change.time` and one afterward;
#' `"late"` reverses those periods. These alternatives do not cross the null.
#' `"fh-early"` uses
#' \eqn{b(t)=\{1-F_{10}(t-)\}^{\rho}}, while `"fh-late"` uses
#' \eqn{b(t)=\min\{F_{10}(t-)/F_{10}(analysis.tau),1\}^{\gamma}}.
#'
#' The final cause-2 distribution can differ slightly by `A` even when
#' `beta2[1] = 0`, because its raw CIF is scaled by the remaining probability
#' `1 - F1(tau | A, L)` to maintain a valid joint competing-risks law.
#'
#' @examples
#' set.seed(2026)
#' dat <- cifsimulate(500, effect = "early", shr = 0.65)
#' head(dat[c("time", "status", "A", "L")])
#'
#' fit <- ciftest(
#'   Event(time, status) ~ A,
#'   data = dat,
#'   augmentation = FALSE,
#'   tau = 1
#' )
#'
#' truth <- attr(dat, "truth")
#' truth$summary$time_averaged_cif_difference
#'
#' @seealso [ciftest()], [cifcurve()], [polyreg()]
#' @export
cifsimulate <- function(
    n,
    effect = c("constant", "early", "late", "fh-early", "fh-late"),
    shr = 0.65,
    rho1 = 0.2,
    rho2 = 1,
    rate = c(1, 1),
    beta1.L = -0.1,
    beta2 = c(A = 0, L = 0.3),
    tau = 6,
    analysis.tau = 1,
    change.time = NULL,
    effect.rho = 1,
    effect.gamma = 1,
    grid.width = 0.01,
    prob.A = 0.5,
    prob.L = 0.5,
    censor.rate = 0,
    censor.log.hazard = c(A = 0, L = 0),
    A = NULL,
    L = NULL,
    uniforms = NULL) {
  effect <- match.arg(effect)
  if (length(n) != 1L || is.na(n) || !is.finite(n) || n <= 0 ||
      n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  scalar_finite <- function(x) {
    length(x) == 1L && !is.na(x) && is.finite(x)
  }
  if (!scalar_finite(shr) || shr <= 0) {
    stop("`shr` must be a positive finite scalar.", call. = FALSE)
  }
  if (!scalar_finite(rho1) || rho1 <= 0 ||
      !scalar_finite(rho2) || rho2 <= 0) {
    stop("`rho1` and `rho2` must be positive finite scalars.", call. = FALSE)
  }
  if (length(rate) != 2L || anyNA(rate) || any(!is.finite(rate)) ||
      any(rate <= 0)) {
    stop("`rate` must contain two positive finite scales.", call. = FALSE)
  }
  if (!scalar_finite(beta1.L)) {
    stop("`beta1.L` must be a finite scalar.", call. = FALSE)
  }
  beta2 <- stats::setNames(as.numeric(beta2), c("A", "L"))
  if (length(beta2) != 2L || anyNA(beta2) || any(!is.finite(beta2))) {
    stop("`beta2` must contain finite coefficients for A and L.",
         call. = FALSE)
  }
  if (!scalar_finite(tau) || tau <= 0 ||
      !scalar_finite(analysis.tau) || analysis.tau <= 0 ||
      analysis.tau > tau) {
    stop("Require 0 < `analysis.tau` <= `tau`.", call. = FALSE)
  }
  if (!scalar_finite(grid.width) || grid.width <= 0) {
    stop("`grid.width` must be a positive finite scalar.", call. = FALSE)
  }
  if (!scalar_finite(effect.rho) || effect.rho < 0 ||
      !scalar_finite(effect.gamma) || effect.gamma < 0) {
    stop("`effect.rho` and `effect.gamma` must be non-negative finite scalars.",
         call. = FALSE)
  }
  for (argument in c("prob.A", "prob.L")) {
    value <- get(argument, inherits = FALSE)
    if (!scalar_finite(value) || value < 0 || value > 1) {
      stop(sprintf("`%s` must be a finite scalar in [0, 1].", argument),
           call. = FALSE)
    }
  }
  if (!scalar_finite(censor.rate) || censor.rate < 0) {
    stop("`censor.rate` must be a non-negative finite scalar.", call. = FALSE)
  }
  censor.log.hazard <- stats::setNames(
    as.numeric(censor.log.hazard), c("A", "L")
  )
  if (length(censor.log.hazard) != 2L || anyNA(censor.log.hazard) ||
      any(!is.finite(censor.log.hazard))) {
    stop("`censor.log.hazard` must contain finite coefficients for A and L.",
         call. = FALSE)
  }

  binary_vector <- function(x, probability, name) {
    if (is.null(x)) return(stats::rbinom(n, 1, probability))
    if (is.factor(x) || is.character(x)) {
      x <- suppressWarnings(as.numeric(as.character(x)))
    } else {
      x <- as.numeric(x)
    }
    if (length(x) != n || anyNA(x) || any(!is.finite(x)) ||
        any(!x %in% c(0, 1))) {
      stop(sprintf("`%s` must be a binary 0/1 vector of length n.", name),
           call. = FALSE)
    }
    x
  }
  A <- binary_vector(A, prob.A, "A")
  L <- binary_vector(L, prob.L, "L")

  if (is.null(uniforms)) uniforms <- list()
  if (!is.list(uniforms) || is.null(names(uniforms)) && length(uniforms)) {
    stop("`uniforms` must be a named list.", call. = FALSE)
  }
  allowed_uniforms <- c("cause1", "cause2", "event", "censor")
  if (any(!names(uniforms) %in% allowed_uniforms) ||
      anyDuplicated(names(uniforms))) {
    stop("`uniforms` may contain cause1, cause2, event, and censor once each.",
         call. = FALSE)
  }
  get_uniform <- function(name) {
    x <- uniforms[[name]]
    if (is.null(x)) return(stats::runif(n))
    x <- as.numeric(x)
    if (length(x) != n || anyNA(x) || any(!is.finite(x)) ||
        any(x < 0 | x > 1)) {
      stop(sprintf("`uniforms$%s` must contain n values in [0, 1].", name),
           call. = FALSE)
    }
    x
  }
  u1 <- get_uniform("cause1")
  u2 <- get_uniform("cause2")
  ue <- get_uniform("event")
  uc <- if (censor.rate > 0) get_uniform("censor") else NULL

  lambda10 <- function(x) rho1 * (1 - exp(-x / rate[1L]))
  cif10 <- function(x) 1 - exp(-lambda10(x))
  lambda20 <- function(x) rho2 * (1 - exp(-x / rate[2L]))
  if (is.null(change.time)) {
    target <- 0.5 * cif10(analysis.tau)
    change.time <- stats::uniroot(
      function(x) cif10(x) - target,
      interval = c(0, analysis.tau), tol = 1e-12
    )$root
  }
  if (!scalar_finite(change.time) || change.time <= 0 ||
      change.time >= analysis.tau) {
    stop("`change.time` must be strictly between 0 and `analysis.tau`.",
         call. = FALSE)
  }

  time_grid <- sort(unique(c(
    seq(0, tau, by = grid.width), tau, analysis.tau, change.time
  )))
  time_grid <- time_grid[time_grid <= tau]
  interval_start <- head(time_grid, -1L)
  interval_end <- tail(time_grid, -1L)
  baseline_cif_left <- cif10(interval_start)
  effect_value <- switch(
    effect,
    constant = rep(1, length(interval_start)),
    early = as.numeric(interval_start < change.time),
    late = as.numeric(interval_start >= change.time),
    `fh-early` = (1 - baseline_cif_left)^effect.rho,
    `fh-late` = pmin(
      baseline_cif_left / cif10(analysis.tau), 1
    )^effect.gamma
  )
  theta <- log(shr)
  baseline_increment1 <- diff(lambda10(time_grid))

  inverse_curve <- function(cumulative, target) {
    target <- pmax(0, pmin(as.numeric(target), tail(cumulative, 1L)))
    keep <- !duplicated(cumulative)
    stats::approx(
      x = cumulative[keep], y = time_grid[keep], xout = target,
      method = "linear", rule = 2, ties = "ordered"
    )$y
  }
  cif1_array <- array(
    NA_real_, c(length(time_grid), 2L, 2L),
    dimnames = list(NULL, A = c("0", "1"), L = c("0", "1"))
  )
  cif2_array <- cif1_array
  cause1_time <- cause2_time <- numeric(n)
  p1 <- p2_raw <- p2 <- numeric(n)
  p1_analysis <- p2_analysis <- numeric(n)
  analysis_index <- match(analysis.tau, time_grid)

  for (a in 0:1) {
    for (ell in 0:1) {
      index <- A == a & L == ell
      cumulative1 <- c(
        0,
        cumsum(baseline_increment1 *
                 exp(beta1.L * ell + theta * a * effect_value))
      )
      cif1 <- 1 - exp(-cumulative1)
      rr2 <- exp(beta2[1L] * a + beta2[2L] * ell)
      cumulative2 <- lambda20(time_grid) * rr2
      raw_cif2 <- 1 - exp(-cumulative2)
      terminal1 <- tail(cif1, 1L)
      terminal2_raw <- tail(raw_cif2, 1L)
      cif2 <- raw_cif2 * (1 - terminal1)
      cif1_array[, a + 1L, ell + 1L] <- cif1
      cif2_array[, a + 1L, ell + 1L] <- cif2
      if (any(index)) {
        p1[index] <- terminal1
        p2_raw[index] <- terminal2_raw
        p2[index] <- terminal2_raw * (1 - terminal1)
        p1_analysis[index] <- cif1[analysis_index]
        p2_analysis[index] <- cif2[analysis_index]
        cause1_time[index] <- inverse_curve(
          cumulative1, -log1p(-u1[index] * terminal1)
        )
        cause2_time[index] <- inverse_curve(
          cumulative2, -log1p(-u2[index] * terminal2_raw)
        )
      }
    }
  }
  if (any(p1 + p2 > 1 + 1e-10)) {
    stop("Generated terminal cause probabilities exceed one.", call. = FALSE)
  }

  event_status <- integer(n)
  event_time <- rep(tau, n)
  cause1 <- ue <= p1
  cause2 <- !cause1 & ue <= p1 + p2
  event_status[cause1] <- 1L
  event_status[cause2] <- 2L
  event_time[cause1] <- cause1_time[cause1]
  event_time[cause2] <- cause2_time[cause2]

  if (censor.rate == 0) {
    censor_time <- rep(Inf, n)
  } else {
    censor_lp <- censor.log.hazard[1L] * A + censor.log.hazard[2L] * L
    censor_time <- -log(pmax(uc, .Machine$double.xmin)) /
      (censor.rate * exp(censor_lp))
  }
  pre_censor_time <- pmin(event_time, tau)
  random_censored <- censor_time < pre_censor_time
  observed_time <- pmin(pre_censor_time, censor_time)
  observed_status <- ifelse(
    !random_censored & event_time <= tau, event_status, 0L
  )

  out <- data.frame(
    id = seq_len(n),
    time = observed_time,
    status = as.integer(observed_status),
    A = factor(A, levels = c(0, 1)),
    L = factor(L, levels = c(0, 1)),
    event_time = event_time,
    event_status = event_status,
    censor_time = censor_time,
    random_censored = random_censored,
    admin_censored = !random_censored & observed_status == 0L,
    p1 = p1,
    p2 = p2,
    p2_raw = p2_raw,
    p1_analysis = p1_analysis,
    p2_analysis = p2_analysis,
    stringsAsFactors = FALSE
  )

  l_weights <- c(1 - prob.L, prob.L)
  marginal_cif1 <- sapply(1:2, function(a) {
    drop(cif1_array[, a, ] %*% l_weights)
  })
  marginal_cif2 <- sapply(1:2, function(a) {
    drop(cif2_array[, a, ] %*% l_weights)
  })
  colnames(marginal_cif1) <- colnames(marginal_cif2) <- c("A0", "A1")
  cif_difference <- marginal_cif1[, "A1"] - marginal_cif1[, "A0"]
  effect_at_grid <- c(effect_value, tail(effect_value, 1L))

  trapezoid_average <- function(value, lower, upper) {
    knots <- sort(unique(c(
      lower, time_grid[time_grid > lower & time_grid < upper], upper
    )))
    values <- stats::approx(time_grid, value, xout = knots, rule = 2)$y
    sum(diff(knots) * (head(values, -1L) + tail(values, -1L)) / 2) /
      (upper - lower)
  }
  analysis_rows <- which(time_grid <= analysis.tau)
  max_index <- analysis_rows[which.max(abs(cif_difference[analysis_rows]))]
  weight_average <- trapezoid_average(effect_at_grid, 0, analysis.tau)
  weighted_difference <- if (weight_average > 0) {
    trapezoid_average(cif_difference * effect_at_grid, 0, analysis.tau) /
      weight_average
  } else {
    NA_real_
  }
  truth_curves <- do.call(rbind, lapply(0:1, function(a) {
    do.call(rbind, lapply(0:1, function(ell) {
      data.frame(
        time = time_grid,
        A = factor(a, levels = c(0, 1)),
        L = factor(ell, levels = c(0, 1)),
        cif1 = cif1_array[, a + 1L, ell + 1L],
        cif2 = cif2_array[, a + 1L, ell + 1L],
        survival = 1 - cif1_array[, a + 1L, ell + 1L] -
          cif2_array[, a + 1L, ell + 1L],
        stringsAsFactors = FALSE
      )
    }))
  }))
  truth_summary <- list(
    effect = effect,
    shr = shr,
    analysis.tau = analysis.tau,
    change.time = change.time,
    time_averaged_cif_difference = trapezoid_average(
      cif_difference, 0, analysis.tau
    ),
    early_time_averaged_cif_difference = trapezoid_average(
      cif_difference, 0, change.time
    ),
    late_time_averaged_cif_difference = trapezoid_average(
      cif_difference, change.time, analysis.tau
    ),
    effect_weighted_cif_difference = weighted_difference,
    max_abs_cif_difference = abs(cif_difference[max_index]),
    cif_difference_at_max = cif_difference[max_index],
    time_of_max_abs_cif_difference = time_grid[max_index],
    endpoint = data.frame(
      A = factor(0:1, levels = c(0, 1)),
      cif1 = marginal_cif1[analysis_index, ],
      cif2 = marginal_cif2[analysis_index, ],
      stringsAsFactors = FALSE
    )
  )
  attr(out, "truth") <- list(
    curves = truth_curves,
    marginal = data.frame(
      time = rep(time_grid, 2L),
      A = factor(rep(0:1, each = length(time_grid)), levels = c(0, 1)),
      cif1 = c(marginal_cif1[, "A0"], marginal_cif1[, "A1"]),
      cif2 = c(marginal_cif2[, "A0"], marginal_cif2[, "A1"]),
      stringsAsFactors = FALSE
    ),
    effect_process = data.frame(
      start = interval_start,
      end = interval_end,
      baseline_cif_left = baseline_cif_left,
      effect_value = effect_value,
      shr_A1_vs_A0 = exp(theta * effect_value),
      stringsAsFactors = FALSE
    ),
    summary = truth_summary
  )
  attr(out, "simulation_parameters") <- list(
    effect = effect, shr = shr, rho1 = rho1, rho2 = rho2, rate = rate,
    beta1.L = beta1.L, beta2 = beta2, tau = tau,
    analysis.tau = analysis.tau, change.time = change.time,
    effect.rho = effect.rho, effect.gamma = effect.gamma,
    grid.width = grid.width, prob.A = prob.A, prob.L = prob.L,
    censor.rate = censor.rate, censor.log.hazard = censor.log.hazard
  )
  out
}
