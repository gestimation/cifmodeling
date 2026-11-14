#' @title Calculate the Kaplan-Meier estimator and the Aalen-Johansen estimator
#' @description
#' Core estimation routine that computes a survfit-compatible object
#' from a formula + data interface (`Event()` or `survival::Surv()` on
#' the LHS, and a stratification variable on the RHS if necessary).
#' The back-end  C++ routine supports both weighted and stratified data. Use this
#' when you want **numbers only** (e.g. estimates, SEs, CIs and influence functions)
#' and will plot it yourself.
#'
#' @inheritParams cif-stat-arguments
#'
#' @param formula A model formula specifying the time-to-event outcome on the LHS
#'   (typically `Event(time, status)` or `survival::Surv(time, status)`)
#'   and, optionally, a stratification variable on the RHS.
#'   `Surv()` responses must be right-censored.
#'   Unlike [cifplot()], this function does not accept a fitted survfit object.
#' @param time.point Optional finite, non-negative time for a one-sided normal
#'   approximation test when `null.hypothesis` is specified. The test returns
#'   a missing p-value for strata with no positive-weight follow-up at this time.
#' @param null.hypothesis Optional null value for the survival probability or
#'   cumulative incidence at `time.point`. For `outcome.type = "survival"`,
#'   the one-sided alternative is that the observed survival probability is
#'   higher than the null value. For `outcome.type = "competing-risk"`, the
#'   one-sided alternative is that the observed CIF for `code.event1` is lower
#'   than the null value. A scalar value is recycled across strata. A named
#'   numeric vector can be used to specify stratum-specific null values.
#' @param report.influence.function Logical. When `TRUE` and `engine = "calculateAJ_Rcpp"`,
#' the influence function is also computed and returned (default `FALSE`).
#' @param report.survfit.std.err Logical. If `TRUE`, report SE on the log-survival
#' scale (survfit's convention). Otherwise SE is on the probability scale.
#' @param engine Character. One of `"auto"`, `"calculateKM"`, or `"calculateAJ_Rcpp"` (default).
#'   `"calculateAJ"` is also accepted as an alias for `"calculateAJ_Rcpp"`.
#'   The legacy engine accepts no weights or positive weights constant within
#'   each stratum, and normal inference with Greenwood/Tsiatis for survival or
#'   delta/Aalen for competing risks. Its SEs use subject counts.
#' @param prob.bound Numeric lower bound used to internally truncate probabilities away from 0 and 1 (default `1e-7`).
#' @param conf.dist Reference distribution for pointwise CIs: `"normal"`
#'   (default) or `"t"`. Student intervals require `error = "if"` or `"fay"`;
#'   omitting `error` with `conf.dist = "t"` selects `"if"`.
#' @param df.method Allocation used for Student degrees of freedom:
#'   `"deletion"` (default, FG-type) or `"equal"` (IF-t only).
#'
#' @details
#'
#' ### Typical use cases
#' - When `outcome.type = "survival"`, this is a thin wrapper around
#' the KM estimator with the chosen variance / CI transformation.
#' - When `outcome.type = "competing-risk"`, this computes the AJ estimator of CIF for `code.event1`.
#'  The returned `$surv` is **1 - CIF**, i.e. in the format that ggsurvfit expects.
#' - Use [cifplot()] if you want to go straight to a figure; use [cifcurve()] if you only want the numbers.
#'
#' ### Risk set display
#' - Set `n.risk.type` to control whether `$n.risk` reflects weighted, unweighted,
#'   or Kish effective sample size (ESS) counts. This only affects the reported
#'   counts (e.g., for plotting or debugging) and leaves estimates and SEs unchanged.
#'
#' ### Standard error and confidence intervals
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `error` | SE for KM: `"greenwood"`, `"tsiatis"`, `"if"`, `"fay"`. For CIF: `"aalen"`, `"delta"`, `"if"`, `"fay"`. | `"greenwood"`, `"delta"` or `"if"` |
#' | `conf.type` | Transformation for CIs: `"plain"`, `"log"`, `"log-log"`, `"arcsin"`, `"logit"`, or `"none"`. | `"arcsin"` |
#' | `conf.int` | Two-sided CI level. | `0.95` |
#'
#' With the C++ engine, `error = "fay"` applies the Fay--Graubard-type
#' sandwich correction with leverage capped at 0.75. Combining `error = "if"`
#' or `"fay"` with `conf.dist = "t"` gives IF-t or FG-t, respectively.
#' `"fg"` and `"fay-graubard"` are accepted aliases for `"fay"`;
#' the returned `error` is always `"fay"` for this correction.
#' All supplied weights are held fixed, including externally estimated IPTW;
#' uncertainty from fitting a propensity-score model is not included.
#' Deletion degrees of freedom use unscaled hazard-score equations and
#' no-event-baseline hazard logits (a positive cause is the baseline at complete
#' exit). Equal allocation is invariant to fixed invertible score transformations
#' within this system; its small-sample coverage need not equal deletion allocation.
#' These are working quadratic-form approximations, not exact Student laws.
#'
#' With the C++ engine, IF and FG intervals use exact boundary singletons; interior zero
#' variance, unavailable degrees of freedom, and times beyond positive-weight
#' follow-up produce missing limits with reasons in `$conf.status`. No normal
#' fallback is used for a failed Student interval. Degrees of freedom vary by
#' stratum and time; calculation can be expensive for large risk sets. The
#' optional one-sided test retains its documented normal reference distribution.
#' Delta and Aalen variances use all event times up to the reported time.
#' If the Aalen formula gives a negative variance with tied events, its SE and
#' confidence limits are missing.
#' Returned influence functions are uncorrected case-weight derivatives (IF/n)
#' of the event probability `1 - $surv`,
#' even when the SE uses FG-adjusted contributions.
#'
#' ### One-sided normal approximation test
#'
#' When both `time.point` and `null.hypothesis` are specified, `cifcurve()`
#' compares the estimate at `time.point` with `null.hypothesis` using a
#' normal approximation. The test statistic is
#' \deqn{
#' Z = (g(\hat{p}) - g(p_0)) / \widehat{SE}(g(\hat{p})).
#' }
#' where `g()` is the transformation specified by `conf.type`. For survival
#' outcomes, \eqn{\hat p} is the Kaplan-Meier survival estimate and
#' the default one-sided alternative is \eqn{S(t) > S_0(t)}. For competing-risk
#' outcomes, the reported estimate is the Aalen-Johansen cumulative incidence
#' estimate \eqn{F_1(t)} for `code.event1`, but the transformation for the
#' one-sided test is applied to \eqn{1 - F_1(t)} and \eqn{1 - F_{10}(t)},
#' consistently with the survfit-compatible `$surv` component. Thus,
#' `conf.type = "log-log"` corresponds to
#' \eqn{\log[-\log(1 - F_1(t))]} for competing-risk outcomes. The default
#' one-sided alternative is \eqn{F_1(t) < F_{10}(t)}, equivalently
#' \eqn{1 - F_1(t) > 1 - F_{10}(t)}. When strata are present, the test is
#' calculated separately within each stratum.

#' @return
#' An object inheriting from `"cifcurve"` and `"survfit"`.
#' For `outcome.type="survival"`, `$surv` is the survival function.
#' For `outcome.type="competing-risk"`, `$surv` equals `1 - CIF` for `code.event1`.
#' SE and CIs are provided per  `error`, `conf.type` and `conf.int`.
#' This enables an independent use of standard methods for `survfit` such as:
#' -   `summary()`: time-by-time estimates with SEs and CIs
#' -   `plot()`: base R stepwise survival/CIF curves
#' -   `mean()`: restricted mean survival estimates with CIs
#' -   `quantile()`: quantile estimates with CIs
#'
#' If `null.hypothesis` is specified, the returned object additionally contains
#' a `one.sided.p` list with the one-sided normal approximation test results.
#' Its table marks strata beyond positive-weight follow-up with
#' `test.status = "outside_observed_followup"`.
#'
#' Note that `$n.risk`, `$n.event`, and `$n.censor` are rounded up to the nearest integer
#' regardless of whether the data is weighted or not.
#' Some methods (e.g. `residuals.survfit`) may not be supported.
#' IF/FG results also record `conf.dist`, `df`, `df.reason`,
#' `critical.value`, and `conf.status`; Student results record `df.method`.
#' These vectors follow the fitted curve's stratum/time ordering. FG results
#' record the fixed leverage cap in `fay.bound`. Influence matrices have one
#' row per original observation within each stratum, including zero-weight rows.
#' `summary(fit, times = ..., extend = TRUE)` extracts the saved intervals,
#' probability-scale SEs and time-specific inference diagnostics. Unsupported
#' times retain missing limits. Use this summary rather than indexing added
#' inference vectors on copies whose time grid was expanded by `survfit0()`.
#' `logse` records the scale of the stored `std.err`. The internal `reporting`
#' component preserves the native inference grid for this extraction.
#' Time and strata display labels are saved with the fit for [cifplot()].
#' @param object A fitted `cifcurve` object to summarize.
#' @param ... Arguments passed to `summary.survfit()`, including `times`.
#'
#' @examples
#' data(diabetes.complications)
#' output1 <- cifcurve(Event(t,epsilon) ~ fruitq,
#'                     data = diabetes.complications,
#'                     outcome.type="competing-risk")
#' cifplot(output1,
#'         outcome.type = "competing-risk",
#'         type.y = "risk",
#'         add.risktable = FALSE,
#'         label.y = "CIF of diabetic retinopathy",
#'         label.x = "Years from registration")
#'
#' output2 <- cifcurve(Event(t,epsilon) ~ fruitq,
#'                     data = diabetes.complications,
#'                     outcome.type = "competing-risk",
#'                     time.point = 8,
#'                     null.hypothesis = 0.30)
#' output2$one.sided.p
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats formula
#'
#' @name cifcurve
#' @section Lifecycle:
#' \lifecycle{stable}
#'
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifplot()] for display of a CIF; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifcurve <- function(
    formula,
    data,
    weights = NULL,
    n.risk.type = "weighted",
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = NULL,
    time.point = NULL,
    null.hypothesis = NULL,
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    report.influence.function = FALSE,
    report.survfit.std.err = FALSE,
    engine = "calculateAJ_Rcpp",
    prob.bound = 1e-7,
    conf.dist = c("normal", "t"),
    df.method = c("deletion", "equal")
) {
  conf.dist <- match.arg(conf.dist)
  df.method <- match.arg(df.method)
  engine <- match.arg(engine, c("calculateAJ_Rcpp", "auto", "calculateKM", "calculateAJ"))
  if (engine == "calculateAJ") engine <- "calculateAJ_Rcpp"
  n.risk.type <- util_check_n_risk_type(n.risk.type)
  substitute_weights <- substitute(weights)
  # Keep bare data-column lookup, but resolve other inputs in their original
  # promise environment. In particular, a forwarded NULL is not a weight.
  if (!is.name(substitute_weights) || !as.character(substitute_weights) %in% names(data))
    substitute_weights <- weights
  has_weights_user <- !is.null(substitute_weights)
  out_read_surv <- eval(substitute(
    util_read_surv(
      formula = formula, data = data, weights = arg_weights,
      code.event1 = code.event1, code.event2 = code.event2, code.censoring = code.censoring,
      subset.condition = subset.condition, na.action = na.action,
      outcome.type = outcome.type
    ),
    list(arg_weights = substitute_weights)
  ))
  outcome.type <- out_read_surv$outcome.type
  if (conf.dist == "t" && is.null(error)) error <- "if"
  if (engine == "calculateKM") {
    choices <- if (outcome.type == "survival") c("greenwood", "tsiatis") else c("delta", "aalen")
    error <- if (is.null(error)) choices[1L] else .cifcurve_normalize_error(error)
    if (conf.dist != "normal" || !error %in% choices || isTRUE(report.influence.function)) {
      stop("Legacy engine requires normal inference with ", paste(choices, collapse = " or "),
        ". Use engine = 'calculateAJ_Rcpp' for IF, Fay--Graubard or Student inference.",
        call. = FALSE)
    }
  } else {
    error <- curve_check_error(error, outcome.type, weights = out_read_surv$w,
      has_weights = has_weights_user)
  }
  shared_ci <- .cifcurve_check_inference(error, conf.dist, df.method, engine, conf.int)
  call <- match.call()

  strata_fac   <- as.factor(out_read_surv$strata)
  strata_lvls  <- levels(strata_fac)
  strata_var   <- out_read_surv$strata_name %||% NULL
  if (!is.null(strata_var)) {
    strata_fullnames <- paste0(strata_var, "=", strata_lvls)
  } else {
    strata_fullnames <- strata_lvls
  }

  epsilon_norm <- rep.int(0L, length(out_read_surv$epsilon))
  epsilon_norm[out_read_surv$epsilon == code.event1]    <- 1L
  epsilon_norm[out_read_surv$epsilon == code.event2]    <- 2L
  epsilon_norm[out_read_surv$epsilon == code.censoring] <- 0L

  if (engine == "calculateKM") {
    fit <- .cifcurve_legacy(out_read_surv, error, n.risk.type, strata_fullnames)
    fit$call <- call
    fit <- .cifcurve_adjusted_ci(fit, out_read_surv$t, epsilon_norm, out_read_surv$w,
      as.integer(out_read_surv$strata), conf.type, conf.int, conf.dist, df.method)
    fit <- .cifcurve_reporting(harmonize_engine_output(fit), out_read_surv, error,
      conf.type, conf.int, conf.dist, report.survfit.std.err)
    return(.cifcurve_add_one_sided_p(fit, time.point, null.hypothesis, outcome.type,
      conf.type, prob.bound, report.survfit.std.err))
  }

  out_cpp <- call_calculateAJ_Rcpp(
    t = out_read_surv$t,
    epsilon = as.integer(epsilon_norm),
    w = if (has_weights_user) out_read_surv$w else NULL,
    strata = as.integer(out_read_surv$strata),
    error = error,
    conf.type = if (shared_ci) "none" else conf.type %||% "none",
    conf.int = conf.int,
    report.influence.function = report.influence.function,
    prob.bound = prob.bound,
    n.risk.type = n.risk.type,
    competing_risk = outcome.type == "competing-risk"
  )
  if (length(strata_fullnames) && length(out_cpp$strata)) {
    names(out_cpp$strata) <- strata_fullnames
  }
  if (length(strata_lvls)) {
    out_cpp$`strata.levels` <- strata_lvls
  }
  out_cpp$call <- call
  out_cpp <- harmonize_engine_output(out_cpp)
  if (shared_ci) {
    out_cpp <- .cifcurve_adjusted_ci(
      out_cpp, out_read_surv$t, epsilon_norm, out_read_surv$w,
      as.integer(out_read_surv$strata), conf.type = conf.type, level = conf.int,
      conf.dist = conf.dist, df.method = df.method
    )
  }
  out_cpp$n.risk <- ceiling(out_cpp$n.risk)
  out_cpp$n.censor <- ceiling(out_cpp$n.censor)
  out_cpp$n.event <- ceiling(out_cpp$n.event)
  out_cpp$n.risk.type <- n.risk.type
  out_cpp <- .cifcurve_reporting(out_cpp, out_read_surv, error,
    conf.type, conf.int, conf.dist, report.survfit.std.err)
  out_cpp <- .cifcurve_add_one_sided_p(
    fit = out_cpp,
    time.point = time.point,
    null.hypothesis = null.hypothesis,
    outcome.type = outcome.type,
    conf.type = conf.type,
    prob.bound = prob.bound,
    report.survfit.std.err = report.survfit.std.err
  )
  return(out_cpp)
}

.cifcurve_legacy <- function(input, error, n.risk.type, labels) {
  competing <- input$outcome.type == "competing-risk"
  groups <- split(seq_along(input$t), as.integer(input$strata))
  pieces <- lapply(groups, function(i) {
    w <- input$w[i]
    if (w[1L] <= 0 || any(w != w[1L])) {
      stop("Legacy engine requires positive weights constant within each stratum. ",
        "Use engine = 'calculateAJ_Rcpp' for other weights.", call. = FALSE)
    }
    time <- sort(unique(input$t[i]))
    first <- input$d1[i]
    events <- if (competing) first + input$d2[i] else first
    counts <- rowsum(cbind(first, events, censor = 1 - events, exits = 1),
      match(input$t[i], time))
    Y <- rev(cumsum(rev(counts[, "exits"])))
    d <- counts[, "events"]
    survival <- cumprod(1 - d / Y)
    if (competing) {
      S <- c(1, utils::head(survival, -1L))
      increment <- S * counts[, "first"] / Y
      cif <- cumsum(increment)
      se <- .cifcurve_legacy_aj_se(Y, d, counts[, "first"], S, increment, error)
      if (utils::tail(survival, 1L) == 0 && all(counts[, "first"] == d)) {
        cif[length(cif)] <- 1
        se[length(se)] <- 0
      }
      survival <- 1 - cif
    } else {
      denominator <- if (error == "greenwood") Y - d else Y
      se <- survival * sqrt(cumsum(d / (Y * ifelse(denominator > 0, denominator, Inf))))
    }
    data.frame(time = time, surv = survival, std.err = se,
      n.risk.weighted = w[1L] * Y, n.risk.unweighted = Y, n.risk.ess = Y,
      n.event = ceiling(w[1L] * d), n.censor = ceiling(w[1L] * counts[, "censor"]))
  })
  out <- as.list(do.call(rbind, pieces))
  out$n.risk <- ceiling(out[[paste0("n.risk.", n.risk.type)]])
  out$n.risk.type <- n.risk.type
  out$error <- error
  out$n <- lengths(groups, use.names = FALSE)
  if (length(groups) > 1L) out$strata <- stats::setNames(vapply(pieces, nrow, integer(1)), labels)
  out$type <- if (competing) "aalen-johansen" else "kaplan-meier"
  out$method <- if (competing) "aalen-johansen" else "Kaplan-Meier"
  if (competing) out$std.err.cif <- out$std.err.aj <- out$std.err
  out
}

.cifcurve_legacy_aj_se <- function(Y, d, d1, S, increment, error) {
  aalen <- error == "aalen"
  quadratic <- linear <- variance <- scale <- 0
  se <- numeric(length(Y))
  for (j in seq_along(Y)) {
    x <- increment[j]
    divisor <- if (aalen) Y[j] - 1 else Y[j]
    direct <- if (divisor > 0) x * (S[j] - x) / divisor else 0
    change <- x * (2 * linear + x * quadratic)
    variance <- variance + change + direct
    scale <- scale + abs(change) + direct
    linear <- linear + x * quadratic
    remaining <- Y[j] - d[j]
    if (d[j] > 0 && remaining > 0 && divisor > 0) {
      quadratic <- quadratic + d[j] / (divisor * remaining)
      linear <- linear - if (aalen) x * (Y[j] - d1[j]) / (divisor * remaining) else x / Y[j]
    }
    se[j] <- if (variance < -64 * .Machine$double.eps * scale) NA_real_ else sqrt(max(0, variance))
  }
  se
}

.cifcurve_normalize_error <- function(x) {
  z <- tolower(trimws(as.character(x)))
  if (length(z) != 1L || is.na(z)) stop("'error' must be one method name.", call. = FALSE)
  if (z %in% c("g", "greenwood")) return("greenwood")
  if (z %in% c("t", "tsiatis")) return("tsiatis")
  if (z %in% c("if", "influence function", "influence_function",
               "influence curve", "ic")) return("if")
  if (z %in% c("fay", "fg", "fay-graubard")) return("fay")
  if (z %in% c("a", "aalen")) return("aalen")
  if (z %in% c("d", "delta")) return("delta")
  z
}

curve_check_error <- function(x, outcome.type, weights = NULL, has_weights = NULL) {
  ot <- util_check_outcome_type(x = outcome.type, auto_message = FALSE)

  choices <- switch(
    ot,
    "survival"       = c("greenwood", "tsiatis", "if", "fay"),
    "competing-risk" = c("aalen", "delta", "if", "fay"),
    stop(sprintf("Invalid outcome.type: %s", outcome.type), call. = FALSE)
  )

  if (is.null(has_weights)) {
    has_weights <- !is.null(weights) && length(weights) > 0
  }

  fallback <- if (has_weights) {
    "if"
  } else if (ot == "survival") {
    "greenwood"
  } else {
    "delta"
  }

  if (is.null(x)) return(fallback)

  x_norm <- .cifcurve_normalize_error(x)
  if (has_weights && !x_norm %in% c("if", "fay")) {
    warning(
      sprintf("%s with weights: error='%s' is not supported; falling back to 'if'.",
              ot, as.character(x)),
      call. = FALSE
    )
    return("if")
  }

  if (x_norm %in% choices) return(x_norm)

  warning(
    sprintf("%s: unsupported error='%s'; falling back to '%s'.",
            ot, as.character(x), fallback),
    call. = FALSE
  )
  fallback
}


call_calculateAJ_Rcpp <- function(t, epsilon, w = NULL, strata = NULL,
                                  error = "greenwood",
                                  conf.type = "arcsin",
                                  conf.int = 0.95,
                                  report.influence.function = FALSE,
                                  prob.bound = 1e-5,
                                  n.risk.type = "weighted",
                                  competing_risk = FALSE) {
  calculateAJ_Rcpp(
    t = t,
    epsilon = epsilon,
    w = w %||% numeric(),
    strata = strata %||% integer(),
    error = error,
    conf_type = conf.type,
    conf_int = conf.int,
    return_if = isTRUE(report.influence.function),
    prob_bound = prob.bound,
    n_risk_type = n.risk.type,
    competing_risk = competing_risk
  )
}

#' @keywords internal
harmonize_engine_output <- function(out) {
  n <- length(out$time %||% numeric())
  out$lower               <- out$lower               %||% out$low
  out$upper               <- out$upper               %||% out$high
  if (length(out$std.err.aj)) out$std.err.cif <- out$std.err.aj
  out$std.err.cif          <- out$std.err.cif        %||% rep(NA_real_, n)
  out$`influence.function` <- out$`influence.function` %||% list()
  out$strata               <- out$strata               %||% integer()
  out$`strata.levels`      <- out$`strata.levels`      %||% integer()
  out$type                 <- out$type                 %||% "kaplan-meier"
  out$method               <- out$method               %||% "Kaplan-Meier"
  out$`conf.type`          <- out$`conf.type`          %||% "arcsine-square root"
  out
}

#' @keywords internal
.cifcurve_add_one_sided_p <- function(fit,
                                      time.point,
                                      null.hypothesis,
                                      outcome.type,
                                      conf.type,
                                      prob.bound,
                                      report.survfit.std.err) {
  if (is.null(null.hypothesis)) {
    return(fit)
  }

  fit$one.sided.p <- .cifcurve_one_sided_p(
    fit = fit,
    time.point = time.point,
    null.hypothesis = null.hypothesis,
    outcome.type = outcome.type,
    conf.type = conf.type,
    prob.bound = prob.bound,
    report.survfit.std.err = report.survfit.std.err
  )

  .cifcurve_print_one_sided_p(fit$one.sided.p)

  fit
}

#' @keywords internal
.cifcurve_one_sided_p <- function(fit,
                                  time.point,
                                  null.hypothesis,
                                  outcome.type,
                                  conf.type,
                                  prob.bound = 1e-7,
                                  report.survfit.std.err = FALSE) {
  if (is.null(time.point)) {
    stop(
      "'time.point' must be specified when 'null.hypothesis' is specified.",
      call. = FALSE
    )
  }

  if (!is.numeric(time.point) || length(time.point) != 1L ||
      !is.finite(time.point) || time.point < 0) {
    stop("'time.point' must be a single finite, non-negative numeric value.", call. = FALSE)
  }

  if (!is.numeric(null.hypothesis)) {
    stop("'null.hypothesis' must be numeric.", call. = FALSE)
  }

  if (any(is.na(null.hypothesis))) {
    stop("'null.hypothesis' must not contain NA.", call. = FALSE)
  }

  if (any(null.hypothesis <= 0 | null.hypothesis >= 1)) {
    stop("'null.hypothesis' must be strictly between 0 and 1.", call. = FALSE)
  }

  if (isTRUE(report.survfit.std.err)) {
    stop(
      "'null.hypothesis' requires probability-scale standard errors; ",
      "set 'report.survfit.std.err = FALSE'.",
      call. = FALSE
    )
  }

  estimate_at_time <- .cifcurve_extract_at_time(
    fit = fit,
    time.point = time.point
  )

  strata <- estimate_at_time$strata
  surv_est <- estimate_at_time$surv
  se_prob <- estimate_at_time$std.err
  support <- vapply(fit$reporting$curves[strata], `[[`, numeric(1), "support")
  supported <- time.point <= unname(support)
  surv_est[!supported] <- se_prob[!supported] <- NA_real_

  null_value <- .cifcurve_match_null_hypothesis(
    null.hypothesis = null.hypothesis,
    strata = strata
  )

  if (any(!is.finite(se_prob[supported]))) {
    stop(
      "Standard errors are missing at 'time.point'; one-sided p-values cannot be computed.",
      call. = FALSE
    )
  }

  if (any(se_prob[supported] <= 0)) {
    stop(
      "Standard errors must be positive at 'time.point'; one-sided p-values cannot be computed.",
      call. = FALSE
    )
  }

  outcome.type <- .cifcurve_normalize_outcome_type(outcome.type)

  if (identical(outcome.type, "survival")) {
    estimate <- surv_est
    estimate.name <- "survival"

    test_estimate <- estimate
    test_null <- null_value
    test_scale <- "survival"

    alternative <- "observed survival probability is higher than null.hypothesis"
    alternative.direction <- "greater"
    test_alternative.direction <- "greater"

  } else if (identical(outcome.type, "competing-risk")) {
    estimate <- 1 - surv_est
    estimate.name <- "cumulative incidence"

    ## For competing-risk outcomes, the returned survfit object stores
    ## surv = 1 - CIF. Therefore transformations for the one-sided test
    ## are applied to 1 - CIF and 1 - null.hypothesis.
    test_estimate <- 1 - estimate
    test_null <- 1 - null_value
    test_scale <- "1 - cumulative incidence"

    alternative <- "observed cumulative incidence is lower than null.hypothesis"

    ## CIF < null is equivalent to 1 - CIF > 1 - null.
    test_alternative.direction <- "greater"

  } else {
    stop(
      "'null.hypothesis' one-sided testing is currently supported only for ",
      "'survival' and 'competing-risk' outcomes.",
      call. = FALSE
    )
  }

  tr_est <- .cifcurve_transform_probability(
    p = test_estimate,
    se = se_prob,
    conf.type = conf.type,
    prob.bound = prob.bound
  )

  tr_null <- .cifcurve_transform_probability(
    p = test_null,
    se = rep(NA_real_, length(test_null)),
    conf.type = conf.type,
    prob.bound = prob.bound
  )

  z <- (tr_est$value - tr_null$value) / tr_est$se

  p_value <- .cifcurve_one_sided_tail_pvalue(
    z = z,
    transform.sign = tr_est$sign,
    alternative.direction = test_alternative.direction
  )

  out <- data.frame(
    strata = strata,
    time.point = rep(time.point, length(estimate)),
    estimate.type = rep(estimate.name, length(estimate)),
    estimate = estimate,
    null.hypothesis = null_value,
    std.err = se_prob,
    transformation.scale = rep(test_scale, length(estimate)),
    estimate.for.test = test_estimate,
    null.for.test = test_null,
    estimate.transformed = tr_est$value,
    null.transformed = tr_null$value,
    std.err.transformed = tr_est$se,
    z = z,
    p.value = p_value,
    test.status = ifelse(supported, "ok", "outside_observed_followup"),
    stringsAsFactors = FALSE
  )

  list(
    method = "One-sided normal approximation test",
    outcome.type = outcome.type,
    alternative = alternative,
    conf.type = .cifcurve_normalize_conf_type(conf.type),
    time.point = time.point,
    table = out
  )
}

#' @keywords internal
.cifcurve_match_null_hypothesis <- function(null.hypothesis, strata) {
  if (length(null.hypothesis) == 1L) {
    return(rep(as.numeric(null.hypothesis), length(strata)))
  }

  if (!is.null(names(null.hypothesis))) {
    out <- as.numeric(null.hypothesis[match(strata, names(null.hypothesis))])

    if (any(is.na(out))) {
      stop(
        "When 'null.hypothesis' is named, its names must match all strata names.",
        call. = FALSE
      )
    }

    return(out)
  }

  if (length(null.hypothesis) != length(strata)) {
    stop(
      "'null.hypothesis' must be length 1, named by strata, or have one value per stratum.",
      call. = FALSE
    )
  }

  as.numeric(null.hypothesis)
}

#' @keywords internal
.cifcurve_normalize_outcome_type <- function(outcome.type) {
  if (is.null(outcome.type)) {
    stop(
      "'outcome.type' must be resolved before computing one-sided p-values.",
      call. = FALSE
    )
  }

  x <- tolower(as.character(outcome.type)[1])
  x <- gsub("_", "-", x, fixed = TRUE)

  if (x %in% c("s", "surv", "survival")) {
    return("survival")
  }

  if (x %in% c("c", "cr", "competing-risk", "competing risk")) {
    return("competing-risk")
  }

  x
}

#' @keywords internal
.cifcurve_transform_probability <- function(p,
                                            se,
                                            conf.type,
                                            prob.bound = 1e-7) {
  type <- .cifcurve_normalize_conf_type(conf.type)

  p0 <- pmin(pmax(as.numeric(p), prob.bound), 1 - prob.bound)

  if (type == "plain") {
    value <- p0
    deriv <- rep(1, length(p0))
  } else if (type == "arcsine-square root") {
    value <- asin(sqrt(p0))
    deriv <- 1 / (2 * sqrt(p0) * sqrt(1 - p0))
  } else if (type == "log") {
    value <- log(p0)
    deriv <- 1 / p0
  } else if (type == "log-log") {
    value <- log(-log(p0))
    deriv <- 1 / (p0 * log(p0))
  } else if (type == "logit") {
    value <- log(p0 / (1 - p0))
    deriv <- 1 / (p0 * (1 - p0))
  } else {
    stop("Unsupported 'conf.type'.", call. = FALSE)
  }

  list(
    value = value,
    se = abs(deriv) * se,
    sign = sign(deriv)
  )
}

#' @keywords internal
.cifcurve_normalize_conf_type <- function(conf.type) {
  if (is.null(conf.type)) {
    return("arcsine-square root")
  }

  x <- tolower(trimws(as.character(conf.type)[1]))
  x <- gsub("_", "-", x, fixed = TRUE)

  if (x %in% c("none", "n", "plain", "identity", "p", "linear")) {
    return("plain")
  }

  if (x %in% c("arcsine-square root", "arcsine", "arcsin", "asin", "a")) {
    return("arcsine-square root")
  }

  if (x %in% c("log")) {
    return("log")
  }

  if (x %in% c("log-log", "loglog")) {
    return("log-log")
  }

  if (x %in% c("logit")) {
    return("logit")
  }

  stop(
    "'conf.type' must be one of 'arcsine-square root', 'plain', 'log', 'log-log', or 'logit'.",
    call. = FALSE
  )
}

#' @keywords internal
.cifcurve_one_sided_tail_pvalue <- function(z,
                                            transform.sign,
                                            alternative.direction) {
  if (alternative.direction == "greater") {
    ifelse(
      transform.sign >= 0,
      stats::pnorm(z, lower.tail = FALSE),
      stats::pnorm(z, lower.tail = TRUE)
    )
  } else if (alternative.direction == "less") {
    ifelse(
      transform.sign >= 0,
      stats::pnorm(z, lower.tail = TRUE),
      stats::pnorm(z, lower.tail = FALSE)
    )
  } else {
    stop("Unknown alternative direction.", call. = FALSE)
  }
}

#' @keywords internal
.cifcurve_print_one_sided_p <- function(x) {
  cat("\n")
  cat(x$method, "\n", sep = "")
  cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("Time point: ", x$time.point, "\n", sep = "")
  cat("Transformation: ", x$conf.type, "\n", sep = "")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' @keywords internal
.cifcurve_extract_at_time <- function(fit, time.point) {
  n_time <- length(fit$time %||% numeric())

  if (n_time == 0L) {
    stop("No estimate is available in 'fit'.", call. = FALSE)
  }

  has_strata <- !is.null(fit$strata) &&
    length(fit$strata) > 0L &&
    sum(as.integer(fit$strata)) > 0L

  if (has_strata) {
    counts <- as.integer(fit$strata)

    if (sum(counts) != n_time) {
      stop(
        "The 'strata' component of the survfit object is inconsistent with 'time'.",
        call. = FALSE
      )
    }

    strata_names <- names(fit$strata)
    if (is.null(strata_names) || any(strata_names == "")) {
      strata_names <- paste0("strata=", seq_along(counts))
    }

    end <- cumsum(counts)
    start <- c(1L, utils::head(end, -1L) + 1L)

    index_list <- Map(seq.int, start, end)
  } else {
    strata_names <- "All"
    index_list <- list(seq_len(n_time))
  }

  out <- lapply(index_list, function(idx) {
    idx_at_time <- idx[fit$time[idx] <= time.point]

    if (!length(idx_at_time)) {
      return(list(
        surv = 1,
        std.err = 0
      ))
    }

    j <- utils::tail(idx_at_time, 1L)

    list(
      surv = as.numeric(fit$surv[j]),
      std.err = as.numeric(fit$std.err[j])
    )
  })

  data.frame(
    strata = strata_names,
    surv = vapply(out, `[[`, numeric(1), "surv"),
    std.err = vapply(out, `[[`, numeric(1), "std.err"),
    stringsAsFactors = FALSE
  )
}
