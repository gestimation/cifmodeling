#' @title Calculate the Kaplan–Meier estimator and the Aalen–Johansen estimator
#' @description
#' Core estimation routine that computes a \code{survfit}-compatible object
#' from a formula + data interface (\code{Event()} or \code{survival::Surv()} on
#' the LHS, and a stratification variable on the RHS if necessary).
#' Use this when you want **numbers only** (KM / CIF + SE + CI) and
#' you will plot it yourself (for example with \code{ggsurvfit} or \code{\link{cifplot}}).
#'
#' @inheritParams cif-stat-arguments
#'
#' @param formula A model formula specifying the time-to-event outcome on the LHS
#'   (typically \code{Event(time, status)} or \code{survival::Surv(time, status)})
#'   and, optionally, a stratification variable on the RHS.
#'   Unlike \code{\link{cifplot}}, this function does not accept a fitted
#'   \code{survfit} object.
#' @param report.influence.function Logical. When \code{TRUE} and \code{engine = "calculateAJ_Rcpp"},
#' the influence function is also computed and returned (default \code{FALSE}).
#' @param report.survfit.std.err Logical. If \code{TRUE}, report SE on the log-survival
#' scale (survfit's convention). Otherwise SE is on the probability scale.
#' @param engine Character. One of \code{"auto"}, \code{"calculateKM"}, or \code{"calculateAJ_Rcpp"} (default \code{"calculateAJ_Rcpp"}).
#' @param prob.bound Numeric lower bound used to internally truncate probabilities away from 0 and 1 (default \code{1e-7}).
#'
#' @details
#'
#' ### Typical use cases
#' - When \code{outcome.type = "survival"}, this is a thin wrapper around
#' the KM estimator with the chosen variance / CI transformation.
#' - When \code{outcome.type = "competing-risk"}, this computes the AJ
#'   estimator of CIF for \code{code.event1}. The returned \code{$surv} is
#'   \code{1 - CIF}, i.e. in the format that \pkg{ggsurvfit} expects.
#' - Use \code{\link{cifplot}} if you want to go straight to a figure; use
#'   \code{cifcurve()} if you only want the numbers.
#'
#' ### Standard error and confidence intervals
#'
#' | Argument | Description | Default |
#' |---|---|---|
#' | `error` | SE for KM: `"greenwood"`, `"tsiatis"`, `"if"`. For CIF: `"aalen"`, `"delta"`, `"if"`. | `"greenwood"` or `"delta"` |
#' | `conf.type` | Transformation for CIs: `"plain"`, `"log"`, `"log-log"`, `"arcsin"`, `"logit"`, or `"none"`. | `"arcsin"` |
#' | `conf.int` | Two-sided CI level. | `0.95` |
#'
#' @returns A \code{survfit} object. For \code{outcome.type="survival"}, \code{$surv} is the survival function.
#' For \code{outcome.type="competing-risk"}, \code{$surv} equals \code{1 - CIF} for \code{code.event1}.
#' SE and CIs are provided per  \code{error}, \code{conf.type} and \code{conf.int}.
#' This enables an independent use of standard methods for \code{survfit} such as:
#' -   `summary()` — time-by-time estimates with SEs and CIs
#' -   `plot()` — base R stepwise survival/CIF curves
#' -   `mean()` — restricted mean survival estimates with CIs
#' -   `quantile()` — quantile estimates with CIs
#'
#' Note that some methods (e.g., \code{residuals.survfit}) may not be supported.
#'
#' @examples
#' if (requireNamespace("cifmodeling", quietly = TRUE) &&
#'              requireNamespace("ggplot2", quietly = TRUE) &&
#'              requireNamespace("ggsurvfit", quietly = TRUE) &&
#'              requireNamespace("patchwork", quietly = TRUE)) {
#' data(diabetes.complications)
#' out_cifcurve <- cifcurve(Event(t,epsilon) ~ fruitq,
#'                          data = diabetes.complications,
#'                          outcome.type="competing-risk")
#' cifplot(out_cifcurve,
#'         outcome.type = "competing-risk",
#'         type.y = "risk",
#'         add.risktable = FALSE,
#'         label.y = "CIF of diabetic retinopathy",
#'         label.x = "Years from registration")
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats formula
#'
#' @name cifcurve
#' @seealso [polyreg()] for log-odds product modeling of CIFs; [cifplot()] for display of a CIF; [cifpanel()] for display of multiple CIFs; [ggsurvfit][ggsurvfit], [patchwork][patchwork] and [modelsummary][modelsummary] for display helpers.
#' @export
cifcurve <- function(
    formula,
    data,
    weights = NULL,
    subset.condition = NULL,
    na.action = na.omit,
    outcome.type = c("survival","competing-risk"),
    code.event1 = 1,
    code.event2 = 2,
    code.censoring = 0,
    error = NULL,
    conf.type = "arcsine-square root",
    conf.int = 0.95,
    report.influence.function = FALSE,
    report.survfit.std.err = FALSE,
    engine = "calculateAJ_Rcpp",
    prob.bound = 1e-7
) {
  outcome.type  <- util_check_outcome_type(outcome.type, formula = formula, data = data)
  out_readSurv  <- util_read_surv(formula, data, weights,
                                  code.event1, code.event2, code.censoring,
                                  subset.condition, na.action)
  error <- curve_check_error(error, outcome.type)
  call <- match.call()

  strata_fac   <- as.factor(out_readSurv$strata)
  strata_lvls  <- levels(strata_fac)
  strata_var   <- out_readSurv$strata_name %||% NULL
  if (!is.null(strata_var)) {
    strata_fullnames <- paste0(strata_var, "=", strata_lvls)
  } else {
    strata_fullnames <- strata_lvls
  }

  epsilon_norm <- rep.int(0L, length(out_readSurv$epsilon))
  epsilon_norm[out_readSurv$epsilon == code.event1]    <- 1L
  epsilon_norm[out_readSurv$epsilon == code.event2]    <- 2L
  epsilon_norm[out_readSurv$epsilon == code.censoring] <- 0L

  if (identical(outcome.type, "survival") && identical(engine, "calculateKM")) {
    out_km <- calculateKM(out_readSurv$t, out_readSurv$d,
                          out_readSurv$w, as.integer(out_readSurv$strata), error)
    out_km$std.err <- out_km$surv * out_km$std.err
    out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower = NULL)
    if (isTRUE(report.survfit.std.err))
      out_km$std.err <- out_km$std.err / out_km$surv

    survfit_object <- list(
      time      = out_km$time,
      surv      = out_km$surv,
      n         = out_km$n,
      n.risk    = out_km$n.risk,
      n.event   = out_km$n.event,
      n.censor  = out_km$n.censor,
      std.err   = out_km$std.err,
      std.err.cif = out_km$`std.err.cif`,
      upper     = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower     = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type = conf.type,
      call      = call,
      type      = "kaplan-meier",
      method    = "Kaplan-Meier"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      names(out_km$strata) <- strata_fullnames
      survfit_object$strata <- out_km$strata
    }
    survfit_object <- harmonize_engine_output(survfit_object)
    class(survfit_object) <- "survfit"
    return(survfit_object)

  } else if (identical(outcome.type, "competing-risk") && identical(engine, "calculateKM")) {
    out_aj <- calculateAJ(out_readSurv)
    names(out_aj$strata1) <- strata_fullnames

    if (any(as.integer(out_readSurv$strata) != 1)) {
      n <- table(as.integer(out_readSurv$strata))
      rep_list <- mapply(rep, n, out_aj$strata1, SIMPLIFY = FALSE)
      n.risk <- do.call(c, rep_list) -
        out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    } else {
      n <- length(out_readSurv$strata)
      n.risk <- n - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    }

    std_err_cif <- calculateAalenDeltaSE(
      out_aj$time1, out_aj$aj1,
      out_aj$n.event1, out_aj$n.event2,
      n.risk,
      out_aj$time0, out_aj$km0, out_aj$strata1, error
    )
    out_aj$std.err <- std_err_cif
    out_aj$surv <- 1 - out_aj$aj1
    out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower = NULL)
    if (isTRUE(report.survfit.std.err))
      out_aj$std.err <- out_aj$std.err / out_aj$surv

    survfit_object <- list(
      time        = out_aj$time1,
      surv        = out_aj$surv,
      n           = n,
      n.risk      = n.risk,
      n.event     = out_aj$n.event1,
      n.censor    = out_aj$n.censor,
      std.err     = out_aj$std.err,
      std.err.cif = std_err_cif,
      upper       = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$upper,
      lower       = if (is.null(conf.type) || conf.type %in% c("none","n")) NULL else out_ci$lower,
      conf.type   = conf.type,
      call        = call,
      type        = "aalen-johansen",
      method      = "aalen-johansen"
    )
    if (any(as.integer(out_readSurv$strata) != 1))
      survfit_object$strata <- out_aj$strata1

    survfit_object <- harmonize_engine_output(survfit_object)
    class(survfit_object) <- "survfit"
    return(survfit_object)
  }

  out_cpp <- call_calculateAJ_Rcpp(
    t = out_readSurv$t,
    epsilon = as.integer(epsilon_norm),
    w = out_readSurv$w,
    strata = as.integer(out_readSurv$strata),
    error = error,
    conf.type = conf.type,
    conf.int = conf.int,
    report.influence.function = report.influence.function,
    prob.bound = prob.bound
  )
  if (length(strata_fullnames) && length(out_cpp$strata)) {
    names(out_cpp$strata) <- strata_fullnames
  }
  if (length(strata_lvls)) {
    out_cpp$`strata.levels` <- strata_lvls
  }
  out_cpp$call <- call
  out_cpp <- harmonize_engine_output(out_cpp)
  if (isTRUE(report.survfit.std.err)) out_cpp$std.err <- out_cpp$std.err / out_cpp$surv
  class(out_cpp) <- "survfit"
  out_cpp
}

calculateAJ <- function(data) {
  out_km0 <- calculateKM(data$t, data$d0, data$w, as.integer(data$strata), "none")
  km0 <- util_get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata, out_km0$strata.levels)
  ip.weight <- (data$d0 == 0) * ifelse(km0 > 0, 1 / km0, 0)
  d1_ipw <- as.matrix(data$w * data$d1 * ip.weight)

  aj1 <- time1 <- integer(0)
  n.cum.event1 <- n.cum.event2 <- n.cum.censor <- numeric(0)
  n.event1 <- n.event2 <- n.censor <- numeric(0)
  strata1 <- integer(0)
  strata_vec <- as.integer(data$strata)

  for (level in sort(unique(strata_vec))) {
    idx <- which(strata_vec == level)

    sub_t  <- data$t[idx]
    sub_d0 <- data$d0[idx]
    sub_d1 <- data$d1[idx]
    sub_d2 <- data$d2[idx]
    sub_d1_ipw <- d1_ipw[idx, , drop = FALSE]

    o <- order(sub_t)
    sub_t  <- sub_t[o]
    sub_d0 <- sub_d0[o]
    sub_d1 <- sub_d1[o]
    sub_d2 <- sub_d2[o]
    sub_d1_ipw <- sub_d1_ipw[o, , drop = FALSE]

    not_atrisk <- outer(sub_t, sub_t, FUN = ">=")
    sub_aj1 <- as.vector(not_atrisk %*% sub_d1_ipw) / length(sub_t)
    sub_n.censor <- as.vector(not_atrisk %*% as.matrix(sub_d0))
    sub_n.event1 <- as.vector(not_atrisk %*% as.matrix(sub_d1))
    sub_n.event2 <- as.vector(not_atrisk %*% as.matrix(sub_d2))

    keep <- !duplicated(rev(sub_t))
    keep <- rev(keep)

    u_t  <- sub_t[keep]
    u_aj1 <- sub_aj1[keep]
    u_nc  <- sub_n.censor[keep]
    u_ne1 <- sub_n.event1[keep]
    u_ne2 <- sub_n.event2[keep]

    oo <- order(u_t)
    u_t  <- u_t[oo]
    u_aj1 <- u_aj1[oo]
    u_nc  <- u_nc[oo]
    u_ne1 <- u_ne1[oo]
    u_ne2 <- u_ne2[oo]

    inc_nc  <- c(u_nc[1],  diff(u_nc))
    inc_ne1 <- c(u_ne1[1], diff(u_ne1))
    inc_ne2 <- c(u_ne2[1], diff(u_ne2))

    time1 <- c(time1, u_t)
    aj1   <- c(aj1, u_aj1)

    n.cum.censor <- c(n.cum.censor, u_nc)
    n.cum.event1 <- c(n.cum.event1, u_ne1)
    n.cum.event2 <- c(n.cum.event2, u_ne2)

    n.censor <- c(n.censor, inc_nc)
    n.event1 <- c(n.event1, inc_ne1)
    n.event2 <- c(n.event2, inc_ne2)

    strata1 <- c(strata1, length(u_t))
  }

  list(
    time1        = time1,
    aj1          = aj1,
    n.event1     = n.event1,
    n.event2     = n.event2,
    n.censor     = n.censor,
    n.cum.event1 = n.cum.event1,
    n.cum.event2 = n.cum.event2,
    n.cum.censor = n.cum.censor,
    strata1      = strata1,
    time0        = out_km0$time,
    km0          = out_km0$surv
  )
}

curve_check_error <- function(x, outcome.type) {
  ot <- util_check_outcome_type(x = outcome.type, auto_message = FALSE)
  choices <- switch(ot,
                    "survival"          = c("greenwood", "tsiatis", "jackknife", "if"),
                    "competing-risk"    = c("aalen", "delta", "jackknife", "if"),
                    stop(sprintf("Invalid outcome.type: %s", outcome.type), call. = FALSE)
  )
  fallback <- if (ot == "survival") "greenwood" else "delta"

  if (is.null(x)) return(fallback)

  x_norm <- tolower(as.character(x))
  if (x_norm %in% choices) return(x_norm)

  warning(sprintf("%s: unsupported error='%s'; falling back to '%s'.", ot, x_norm, fallback), call. = FALSE)
  return(fallback)
}

call_calculateAJ_Rcpp <- function(t, epsilon, w = NULL, strata = NULL,
                                   error = "greenwood",
                                   conf.type = "arcsin",
                                   conf.int = 0.95,
                                   report.influence.function = FALSE,
                                   prob.bound = 1e-5) {
  calculateAJ_Rcpp(
    t = t,
    epsilon = epsilon,
    w = w %||% numeric(),
    strata = strata %||% integer(),
    error = error,
    conf_type = conf.type,
    conf_int = conf.int,
    return_if = isTRUE(report.influence.function),
    prob_bound = prob.bound
  )
}

#' @keywords internal
harmonize_engine_output <- function(out) {
  n <- length(out$time %||% numeric())
  out$lower               <- out$lower               %||% out$low
  out$upper               <- out$upper               %||% out$high
  out$`std.err.cif`        <- out$`std.err.cif`        %||% rep(NA_real_, n)
  out$`influence.function` <- out$`influence.function` %||% list()
  out$strata               <- out$strata               %||% integer()
  out$`strata.levels`      <- out$`strata.levels`      %||% integer()
  out$type                 <- out$type                 %||% "kaplan-meier"
  out$method               <- out$method               %||% "Kaplan-Meier"
  out$`conf.type`          <- out$`conf.type`          %||% "arcsine-square root"
  out
}

